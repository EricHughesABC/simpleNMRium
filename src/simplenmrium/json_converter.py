"""
json_converter.py

Converts a parsed NMRium/NOMAD dataset (nmrium_reader.NMRDataset) into a
simpleNMR submission payload via simplenmr_builder.SimpleNMRBuilder,
mirroring the structure of simpleNMRjeolTools3's json_converter.py and
simpleNMRbrukerTools's json_converter.py (same builder, same overall
class shape: parse once in __init__, then createJsonDict() assembles the
payload through the builder).

Design decisions confirmed with Eric (2026-09-06), NOT guesses:

- source="nmrium" — added to simplenmr_builder.constants.SOURCES,
  SOURCE_REQUIRED_EXTRA, SOURCE_FORBIDDEN (same profile as "jeol": no
  nmrAssignments, no exptIdentifiers, needs carbonCalcPositionsMethod).
- correlations block is deliberately IGNORED for v1. NMRium's own
  correlation-table atom labels ("C1", "C2", ... "H1", "H11" — see
  correlations.values[i].label.origin) do NOT correspond to the molfile's
  atom order (confirmed: the molfile's own SGROUP custom-label block is
  just sequential 1..14, unrelated to the correlation labels, which skip
  the carbonyl carbon and start renumbering from the HSQC-correlated
  carbons). Reconciling NMRium's correlation labels against RDKit atom
  indices — so that simpleNMR's *output* assignments could be written
  back into the NMRium JSON — is a real, separate piece of work Eric
  wants to tackle together later, once this converter's basic path is
  working. allAtomsInfo/carbonAtomsInfo here are built purely from RDKit
  atom order, exactly like the Bruker/JEOL converters.
- MNOVAcalcMethod = "NMRSHIFTDB2 Predict", carbonCalcPositionsMethod =
  "Calculated Positions" — copied verbatim from both real Bruker and
  real JEOL converters (confirmed identical in both
  simpleNMRbrukerTools/core/json_converter.py and
  simpleNMRjeolTools3/src/simplenmrjeol/json_converter.py). This value
  triggers prediction_from_nmrshiftdb2()==True server-side per the
  contract manifest, i.e. the server computes expected carbon shifts
  itself — matches Eric's own expectation. c13predictions IS still
  submitted, as a present-but-empty envelope (count=0), not omitted —
  confirmed by a real server run (2026-09-06): the server's simulated-
  annealing code path reads json_data["oldjsondata"]["c13predictions"]
  ["data"] unconditionally whenever simulatedAnnealing=True, and KeyErrors
  if the field is absent rather than empty. The real bruker_real_exam_cmcse
  .json golden file confirms this same present-but-empty shape.
- 2D zones: one row per zone using the zone-averaged position
  (Zone2D.average_x/average_y/average_intensity from nmrium_reader), not
  one row per raw sub-peak. Eric's instruction: use the average for now
  and revisit if the server doesn't like it once real submissions are
  tested.
"""

from __future__ import annotations

import platform
import re
import uuid
from pathlib import Path
from typing import Any, Optional, Union

from simplenmr_builder import SimpleNMRBuilder
from simplenmr_builder.constants import NMREXPERIMENTS, SKIP_TOKEN

try:
    from rdkit import Chem

    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

from .nmrium_reader import NMRDataset, Spectrum

# ---------------------------------------------------------------------------
# Experiment-type mapping: NMRium's `info.experiment` (+ nucleus/dimension)
# -> simplenmr_builder.constants.NMREXPERIMENTS token.
#
# Confirmed against the real uploaded file: NMRium uses a generic "1d"
# experiment name for BOTH 1H and 13C 1D spectra (distinguished only by
# `nucleus`), unlike Bruker/JEOL which already know their own experiment
# type from the pulse program name. 2D experiment names ("cosy", "hsqc",
# "hmbc") map straightforwardly. Tokens not seen in any real NMRium file
# yet (dept, hsqc-tocsy, etc.) are listed here as best-guess mappings, not
# confirmed — flagged in the docstring of map_experiment_to_token below.
# ---------------------------------------------------------------------------
_2D_EXPERIMENT_MAP = {
    "cosy": "COSY",
    "hsqc": "HSQC",
    "hmbc": "HMBC",
    "noesy": "NOESY",
    "roesy": "NOESY",  # best guess: ROESY is inert-passthrough like NOESY anyway
    "hsqctocsy": "HSQC_CLIPCOSY",  # unconfirmed guess, flagged below
}


def map_experiment_to_token(spectrum: Spectrum) -> Optional[str]:
    """
    Map an NMRium Spectrum's (experiment, nucleus, dimension) to a
    NMREXPERIMENTS token, or None if there's no known mapping (caller
    should skip the spectrum and print a warning rather than guess).
    """
    exp = spectrum.experiment.lower()

    if spectrum.is_1d:
        if exp in ("1d", "id"):  # "id" typo-guard: seen as a real NMRium quirk in some exports
            nucleus = spectrum.primary_nucleus.upper()
            if nucleus == "1H":
                return "H1_1D"
            if nucleus == "13C":
                return "C13_1D"
            return None  # unrecognized nucleus for a 1D spectrum — skip, don't guess
        if exp == "dept":
            return "DEPT135"  # unconfirmed guess — no real DEPT NMRium file seen yet
        return None

    return _2D_EXPERIMENT_MAP.get(exp)


# ---------------------------------------------------------------------------
# Molfile handling
# ---------------------------------------------------------------------------

# NMRium/OpenChemLib V3000 molfiles carry an SGROUP block (custom atom
# labels, NOSEARCH_OCL_CUSTOM_LABEL) that RDKit's MolFromMolBlock cannot
# parse ("BEGIN SGROUP found but Sgroups NOT expected"), confirmed against
# the real uploaded file. The SGROUP's own labels turned out to be plain
# sequential 1..N (not a match for the `correlations` block's C1/H1-style
# labels — see module docstring), so there's nothing lost by discarding
# it: strip it out before parsing, keep RDKit's own atom order as the
# atom_idx/atomNumber source of truth, exactly like Bruker/JEOL do.
_SGROUP_BLOCK_RE = re.compile(r"M {2}V30 BEGIN SGROUP\n.*?M {2}V30 END SGROUP\n", re.DOTALL)


def strip_sgroup_block(molblock: str) -> str:
    """Remove the V3000 SGROUP block RDKit can't parse. No-op if absent."""
    return _SGROUP_BLOCK_RE.sub("", molblock)


def parse_molecule(molblock: str):
    """Parse an NMRium molfile into an RDKit Mol, stripping SGROUP first."""
    if not RDKIT_AVAILABLE:
        raise RuntimeError("rdkit is required to convert NMRium files (pip install rdkit)")
    cleaned = strip_sgroup_block(molblock)
    mol = Chem.MolFromMolBlock(cleaned)
    if mol is None:
        raise ValueError(
            "RDKit could not parse the NMRium molfile even after stripping the "
            "SGROUP block. The molfile text has been left in the raised context "
            "for inspection."
        )
    return mol


def build_all_atoms_info_records(mol) -> list[dict]:
    """One record per atom, atom_idx/atomNumber from RDKit's own atom order."""
    records = []
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        records.append(
            {
                "atom_idx": idx,
                "id": idx,
                "atomNumber": str(idx + 1),
                "symbol": atom.GetSymbol(),
                "numProtons": atom.GetTotalNumHs(),
            }
        )
    return records


def build_carbon_atoms_info_records(mol) -> list[dict]:
    """Same as build_all_atoms_info_records, filtered to carbon atoms only."""
    return [rec for rec in build_all_atoms_info_records(mol) if rec["symbol"] == "C"]


def get_hostname() -> str:
    """
    MAC-based host identifier. Delegates to simplenmr_builder.gui.
    submission.machine_id() (the actual shared implementation Bruker/JEOL
    both submit and register against) rather than re-deriving it locally
    — falls back to a local hex(uuid.getnode()) only if the gui extra
    isn't installed, so this module still works for pure conversion
    (no Qt/requests) use.
    """
    try:
        from simplenmr_builder.gui.submission import machine_id

        return machine_id()
    except ImportError:
        return hex(uuid.getnode())


# ---------------------------------------------------------------------------
# Spectrum block construction
# ---------------------------------------------------------------------------

def _base_block_fields(spectrum: Spectrum, dim: str) -> dict:
    return {
        "datatype": "nmrspectrum",
        "type": dim,
        "experimenttype": spectrum.experiment,
        "pulsesequence": spectrum.pulse_sequence[0] if spectrum.pulse_sequence else "",
        "solvent": spectrum.solvent,
        "nucleus": spectrum.nucleus,
        "specfrequency": spectrum.base_frequency or spectrum.origin_frequency,
        "temperature": spectrum.temperature[0] if spectrum.temperature else None,
    }


def build_1d_block(spectrum: Spectrum) -> dict:
    """
    Build a H1_1D/C13_1D-shaped block from an NMRium 1D Spectrum.

    peaks: one row per Peak1D, delta1=shift (ppm), delta2=0.0 — confirmed
    shape from the real JEOL C13_1D_0 golden block.

    integrals: one row per Range1D, using a single-dimension analogue of
    the confirmed 2D integrals shape (delta1/rangeMin1/rangeMax1 + a
    representative intensity). NOT confirmed against a real populated 1D
    integrals block (none exists in any golden file — every real example
    has integrals count=0), so treat this as best-effort passthrough data
    rather than load-bearing; the manifest marks most of these sibling
    fields server_consumed=false anyway.

    multiplets: left empty, matching every real Bruker/JEOL/MNova
    converter (none of them populate this despite having multiplicity/J
    data available, e.g. via Signal1D.multiplicity/coupling_constants
    here) — don't diverge from working precedent without evidence the
    server reads it.
    """
    block = _base_block_fields(spectrum, "1D")

    peaks_data = {}
    for i, (shift, intensity) in enumerate(spectrum.get_all_peaks_1d()):
        peaks_data[str(i)] = {
            "intensity": intensity,
            "delta1": shift,
            "delta2": 0.0,
            "annotation": "",
            "type": 0,
        }
    block["peaks"] = {"datatype": "peaks", "count": len(peaks_data), "data": peaks_data}

    integrals_data = {}
    for i, range_obj in enumerate(spectrum.ranges):
        integrals_data[str(i)] = {
            "intensity": range_obj.integration,
            "rangeMin1": range_obj.from_ppm,
            "rangeMax1": range_obj.to_ppm,
            "delta1": range_obj.center,
            "delta2": 0.0,
            "annotation": "",
            "type": 0,
        }
    block["integrals"] = {
        "datatype": "integrals",
        "count": len(integrals_data),
        "normValue": 1,
        "data": integrals_data,
    }

    block["multiplets"] = {"datatype": "multiplets", "count": 0, "normValue": 1, "data": {}}

    return block


def build_2d_block(spectrum: Spectrum) -> dict:
    """
    Build a HSQC/HMBC/COSY-shaped block from an NMRium 2D Spectrum.

    peaks: one row per zone, using the zone-averaged position (Eric's
    instruction — revisit if the server doesn't like collapsed zones once
    tested for real). delta1=F1 (nucleus[1], e.g. 13C), delta2=F2
    (nucleus[0], e.g. 1H) — confirmed convention from the real Bruker
    HSQC_0 golden block.

    integrals: one row per zone too, using its full x/y range as
    rangeMin/rangeMax. IMPORTANT, confirmed by direct inspection of the
    real JEOL golden HSQC_0 block: integrals uses the OPPOSITE delta1/
    delta2 convention from peaks in the very same block — delta1=F2
    (1H)/rangeMin1,rangeMax1 bracket it, delta2=F1 (13C)/rangeMin2,
    rangeMax2 bracket it. This looks like a real quirk in the existing
    data rather than a mistake on this converter's part, so it's
    reproduced deliberately rather than "fixed" to be consistent with
    peaks.
    """
    block = _base_block_fields(spectrum, "2D")

    peaks_data = {}
    integrals_data = {}
    for i, zone in enumerate(spectrum.zones):
        avg_x, avg_y, avg_intensity, total_intensity = (
            zone.average_x,
            zone.average_y,
            zone.average_intensity,
            zone.total_intensity,
        )
        peaks_data[str(i)] = {
            "intensity": avg_intensity,
            "delta1": avg_y,  # F1 (e.g. 13C)
            "delta2": avg_x,  # F2 (e.g. 1H)
            "annotation": "",
            "type": 0,
        }
        integrals_data[str(i)] = {
            "intensity": total_intensity,
            "rangeMin1": zone.x_from,
            "rangeMax1": zone.x_to,
            "rangeMin2": zone.y_from,
            "rangeMax2": zone.y_to,
            "delta1": avg_x,  # F2 (e.g. 1H) — swapped vs. peaks, see docstring
            "delta2": avg_y,  # F1 (e.g. 13C)
            "annotation": "",
            "type": 0,
        }

    block["peaks"] = {"datatype": "peaks", "count": len(peaks_data), "data": peaks_data}
    block["integrals"] = {
        "datatype": "integrals",
        "count": len(integrals_data),
        "normValue": 1,
        "data": integrals_data,
    }
    block["multiplets"] = {"datatype": "multiplets", "count": 0, "normValue": 1, "data": {}}

    return block


def build_chosen_entry(spectrum: Spectrum, spectrum_id: str, token: str) -> str:
    """
    Human-readable chosenSpectra entry string. The exact wording isn't
    contract-parsed (SpectrumBuilder tracks the real block registration
    separately via add_block/add_chosen_candidate insertion order) — this
    follows the Bruker converter's format:
    "<nucleus> <dim>D <experiment> <spectrum_id> <token>".
    """
    if spectrum.is_1d:
        nucleus_str = spectrum.primary_nucleus
    else:
        nucleus_str = f"[{', '.join(spectrum.nucleus)}]"
    return f"{nucleus_str} {spectrum.dimension}D {spectrum.experiment} {spectrum_id} {token}"


def spectrum_display_name(spectrum: Spectrum, index: int) -> str:
    """
    Row label for the spectrum-selection dialog / HTML preview. Includes
    the peak/zone count so the person can tell an empty or near-empty
    spectrum apart from a real one at a glance, which the dialog's plain
    entry-name list otherwise can't show.
    """
    nucleus_str = "/".join(spectrum.nucleus)
    n_peaks = spectrum.total_peaks_1d if spectrum.is_1d else spectrum.total_peaks_2d
    unit = "peaks" if spectrum.is_1d else "zones"
    count = spectrum.n_ranges if spectrum.is_1d else spectrum.n_zones
    return f"[{index}] {spectrum.experiment.upper()} {nucleus_str} ({spectrum.dimension}D, {count} {unit}, {n_peaks} total peaks)"


# ---------------------------------------------------------------------------
# Top-level converter
# ---------------------------------------------------------------------------

class NMRiumData:
    """
    Loads an NMRium/NOMAD JSON file and builds a simpleNMR submission
    payload from it. Mirrors jeolData/brukerData's shape: heavy parsing
    happens once in __init__, createJsonDict() assembles the payload via
    SimpleNMRBuilder.

    Spectrum-to-token assignment is a two-step, dialog-friendly process,
    matching how jeolData.choosePeakPickedSpectaforSimpleNMR() works:

    1. __init__ computes .suggested_assignments — one
       {"experiment_name", "experiment_type", "index"} dict per parsed
       spectrum (same shape SpectrumAssignmentDialog.get_assignments()
       returns), pre-filled from map_experiment_to_token()'s automatic
       guess, or "SKIP" if nothing matched. Nothing is submitted yet.
    2. A caller (the CLI, showing simplenmr_builder.gui.
       spectrum_assignment_dialog.SpectrumAssignmentDialog pre-populated
       with these suggestions so the person can confirm, override, or
       SKIP each row) calls apply_assignments(...) with the confirmed
       list. This builds .blocks / .skipped_spectra from the CONFIRMED
       types, not the auto-detected ones.
    3. If apply_assignments() is never called (e.g. programmatic/test
       use with no dialog), createJsonDict() and has_hsqc() fall back to
       auto-applying .suggested_assignments themselves, so non-interactive
       callers still get a sensible default rather than an empty payload.
    """

    def __init__(self, file_path: Union[str, Path]):
        self.file_path = Path(file_path)
        self.dataset = NMRDataset.from_file(self.file_path)

        molblock = self.dataset.primary_molfile
        if molblock is None:
            raise ValueError(f"No molecule/molfile found in {self.file_path}")
        self.rdkit_mol = parse_molecule(molblock)

        # RDKit round-trips the molfile through its own writer so the
        # atom order recorded in allAtomsInfo/carbonAtomsInfo is
        # GUARANTEED to match the molfile actually submitted alongside
        # it — matching Bruker/JEOL, which both build molfile from the
        # same rdkit_mol object rather than passing NMRium's raw text
        # straight through untouched.
        self.molfile = Chem.MolToMolBlock(self.rdkit_mol)
        self.smiles = self.dataset.primary_smiles or Chem.MolToSmiles(self.rdkit_mol)

        self.all_atoms_info_records = build_all_atoms_info_records(self.rdkit_mol)
        self.carbon_atoms_info_records = build_carbon_atoms_info_records(self.rdkit_mol)

        self.hostname = get_hostname()

        self.suggested_assignments: list[dict] = [
            {
                "experiment_name": spectrum_display_name(spectrum, i),
                "experiment_type": map_experiment_to_token(spectrum) or SKIP_TOKEN,
                "index": i,
            }
            for i, spectrum in enumerate(self.dataset.spectra)
        ]

        # Populated by apply_assignments(); None means "not decided yet"
        # so has_hsqc()/createJsonDict() know whether to fall back to
        # auto-applying suggested_assignments.
        self.blocks: Optional[list[tuple[str, dict, str]]] = None
        self.skipped_spectra: list[Spectrum] = []

    def apply_assignments(self, assignments: list[dict]) -> None:
        """
        Build .blocks/.skipped_spectra from a CONFIRMED assignment list —
        the exact shape simplenmr_builder.gui.spectrum_assignment_dialog
        .SpectrumAssignmentDialog.get_assignments() returns:
        [{"experiment_name", "experiment_type", "index"}, ...], one entry
        per row in the same order .suggested_assignments was built in
        (so "index" maps straight back to self.dataset.spectra[index]).

        A row whose experiment_type is "SKIP" (or anything not in
        NMREXPERIMENTS — defensive, matches the server-rejecting-unknown-
        tokens behavior SpectrumBuilder.add_block already enforces) is
        excluded from submission and recorded in .skipped_spectra instead.
        """
        token_counts: dict[str, int] = {}
        blocks: list[tuple[str, dict, str]] = []
        skipped: list[Spectrum] = []

        for entry in assignments:
            spectrum = self.dataset.spectra[entry["index"]]
            token = entry["experiment_type"]

            if token == SKIP_TOKEN or token not in NMREXPERIMENTS:
                skipped.append(spectrum)
                continue

            count = token_counts.get(token, 0)
            token_counts[token] = count + 1
            spectrum_id = f"{token}_{count}"

            block = build_1d_block(spectrum) if spectrum.is_1d else build_2d_block(spectrum)
            chosen_entry = build_chosen_entry(spectrum, spectrum_id, token)
            blocks.append((token, block, chosen_entry))

        self.blocks = blocks
        self.skipped_spectra = skipped

    def _ensure_assignments_applied(self) -> None:
        if self.blocks is None:
            self.apply_assignments(self.suggested_assignments)

    def has_hsqc(self) -> bool:
        self._ensure_assignments_applied()
        return any(token == "HSQC" for token, _, _ in self.blocks)

    def createJsonDict(
        self,
        simulated_annealing: bool = True,
        ml_consent: bool = False,
        validator_path: Optional[str] = None,
        run_validation: bool = True,
    ) -> dict:
        """
        Build the simpleNMR submission payload. Raises HSQCMissingError /
        ContractError (from simplenmr_builder) if the payload would be
        invalid — matching the "structurally cannot produce an invalid
        submission" guarantee both other converters rely on.
        """
        self._ensure_assignments_applied()

        builder = SimpleNMRBuilder(source="nmrium")

        builder.set_scalar("smiles", self.smiles)
        builder.set_scalar("molfile", self.molfile)
        builder.set_scalar("hostname", self.hostname)
        builder.set_scalar("MNOVAcalcMethod", "NMRSHIFTDB2 Predict")
        builder.set_scalar("carbonCalcPositionsMethod", "Calculated Positions")
        builder.set_scalar("simulatedAnnealing", simulated_annealing)
        builder.set_scalar("ml_consent", ml_consent)
        # Real bug found 2026-09 (Windows only): a raw Windows path here
        # (backslash-separated) can contain a backslash immediately
        # followed by a digit purely by coincidence of the person's own
        # folder structure (e.g. "...\2025\python\..."). The server
        # embeds this string unescaped into a JS template literal
        # (backtick string) in the generated result HTML, and template
        # literals — unlike plain quoted strings — forbid legacy octal
        # escape sequences (\2, \02, etc.), which is exactly what
        # backslash-digit forms as valid-looking JS. Result: "Uncaught
        # SyntaxError: Octal escape sequences are not allowed in
        # template strings" when the results viewer opens the page.
        # Normalizing to forward slashes sidesteps this entirely — valid
        # on Windows, and can never accidentally form an escape
        # sequence. workingFilename is just a bare filename (no
        # separators), so it isn't at risk and is left as-is.
        builder.set_scalar("workingDirectory", str(self.file_path.parent).replace("\\", "/"))
        builder.set_scalar("workingFilename", self.file_path.name)

        builder.set_all_atoms_info(self.all_atoms_info_records)
        builder.set_carbon_atoms_info(self.carbon_atoms_info_records)
        # Confirmed by a real server run (2026-09-06): c13predictions must
        # be PRESENT in the payload whenever simulatedAnnealing=True, even
        # with no real prediction data — the server's simulated_annealing
        # code path does json_data["oldjsondata"]["c13predictions"]["data"]
        # unconditionally, which KeyErrors if the field is missing
        # entirely rather than present-but-empty. Confirmed against the
        # real bruker_real_exam_cmcse.json golden file: it submits
        # {"datatype": "c13predictions", "count": 0, "data": {}} — present
        # and empty, NOT omitted. MNOVAcalcMethod="NMRSHIFTDB2 Predict"
        # still means the server recomputes its own predictions rather
        # than trusting this (confirmed in the same server run's log:
        # "Using nmrshiftDB2 for C13 predictions" happens independently of
        # this field's contents) — so an empty list here is correct, not
        # a placeholder that should eventually be filled with real values.
        builder.set_c13predictions([])

        for token, block, chosen_entry in self.blocks:
            builder.spectra.add_block(token, block)
            builder.spectra.add_chosen_candidate(chosen_entry, skip=False)
            builder.spectra.add_spectrum_with_peaks(chosen_entry)

        return builder.build(validator_path=validator_path, run_validation=run_validation)
