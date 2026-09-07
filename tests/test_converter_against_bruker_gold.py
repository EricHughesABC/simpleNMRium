"""
Cross-checks the NMRium converter's output against a real, independently
produced submission of the SAME physical sample (exam_CMCse_1,
alpha-ionone) from simpleNMRbuilder's own Bruker golden test file.

This is deliberately NOT a byte-for-byte comparison — the two sources
were acquired/processed independently (different instrument software,
different peak-picking), so exact numeric equality is not expected or
desired. What IS checked: same molecule (via RDKit canonical SMILES),
same carbon count, and HSQC 13C shifts agreeing within a instrument-
variance tolerance. A regression here means the converter has drifted
from a result already confirmed sane by hand, not that every future
decimal must match Bruker's exactly.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from rdkit import Chem

from simplenmrium import NMRiumData

HERE = Path(__file__).parent
NMRIUM_FILE = HERE / "golden_files" / "exam_CMCse_1_NOMAD-pp.json"

# simplenmr_builder must be installed editable from a local checkout that
# sits next to this repo for this specific cross-check to find the real
# Bruker golden file; skip gracefully if it's not there (e.g. CI checking
# out simpleNMRium alone).
BRUKER_GOLD = (
    HERE.parent.parent / "simpleNMRbuilder" / "tests" / "golden_files" / "valid" / "bruker_real_exam_cmcse.json"
)


@pytest.fixture(scope="module")
def nmrium_payload():
    data = NMRiumData(NMRIUM_FILE)
    return data.createJsonDict()


def test_converter_produces_valid_payload(nmrium_payload):
    assert "HSQC_0" in nmrium_payload
    assert nmrium_payload["chosenSpectra"]["count"] == 5


def test_all_five_expected_blocks_present(nmrium_payload):
    for token in ("H1_1D_0", "C13_1D_0", "COSY_0", "HSQC_0", "HMBC_0"):
        assert token in nmrium_payload, f"missing expected block {token}"


def test_carbon_count_matches_molecule(nmrium_payload):
    assert nmrium_payload["carbonAtomsInfo"]["count"] == 13


def test_zone_counts_match_known_zone_averaging_csvs(nmrium_payload):
    # From ZONE_AVERAGING_GUIDE.md / the CSVs Eric already generated:
    # cosy 8 zones, hsqc 11 zones, hmbc 17 zones.
    assert nmrium_payload["COSY_0"]["peaks"]["count"] == 8
    assert nmrium_payload["HSQC_0"]["peaks"]["count"] == 11
    assert nmrium_payload["HMBC_0"]["peaks"]["count"] == 17


def test_c13predictions_present_but_empty(nmrium_payload):
    # Regression test for a real server crash (2026-09-06): the server's
    # simulated-annealing code path does
    # json_data["oldjsondata"]["c13predictions"]["data"] unconditionally
    # whenever simulatedAnnealing=True, and KeyErrors if the field is
    # absent from the payload rather than present-but-empty. Confirmed
    # against the real bruker_real_exam_cmcse.json golden file, which
    # submits {"datatype": "c13predictions", "count": 0, "data": {}} —
    # present and empty, not omitted.
    assert "c13predictions" in nmrium_payload
    assert nmrium_payload["c13predictions"]["count"] == 0
    assert nmrium_payload["c13predictions"]["data"] == {}


def test_suggested_assignments_auto_detect_all_five():
    data = NMRiumData(NMRIUM_FILE)
    types = {a["experiment_type"] for a in data.suggested_assignments}
    assert types == {"H1_1D", "C13_1D", "COSY", "HSQC", "HMBC"}


def test_apply_assignments_honors_user_skip():
    # Simulates the spectrum-selection dialog: person confirms everything
    # except manually marking COSY as SKIP.
    data = NMRiumData(NMRIUM_FILE)
    confirmed = []
    for a in data.suggested_assignments:
        a = dict(a)
        if a["experiment_type"] == "COSY":
            a["experiment_type"] = "SKIP"
        confirmed.append(a)

    data.apply_assignments(confirmed)
    tokens = [t for t, _, _ in data.blocks]
    assert "COSY" not in tokens
    assert "HSQC" in tokens
    assert len(data.skipped_spectra) == 1
    assert data.skipped_spectra[0].experiment == "cosy"

    payload = data.createJsonDict()
    assert "COSY_0" not in payload
    assert "HSQC_0" in payload


def test_working_directory_has_no_backslashes_even_on_windows_paths():
    # Regression test for a real bug found 2026-09 on a Windows machine:
    # a raw backslash-separated path submitted as workingDirectory can
    # contain a backslash immediately followed by a digit purely from
    # the person's own folder structure (e.g. "...\2025\python\...").
    # The server embeds this unescaped into a JS template literal in the
    # generated results HTML, and template literals forbid legacy octal
    # escape sequences (\2, \02, ...) — exactly what backslash-digit
    # looks like — causing "Uncaught SyntaxError: Octal escape sequences
    # are not allowed in template strings" when the results viewer opens
    # the page. Bruker's and JEOL's own converters already sanitize this
    # the same way; this was purely a gap in NMRiumData specifically.
    from pathlib import PureWindowsPath

    data = NMRiumData(NMRIUM_FILE)

    fake_windows_dir = PureWindowsPath(
        r"C:\Users\vsmw51\OneDrive - Durham University\projects\programming\2025\python\awh"
    )

    class FakeWindowsPath:
        parent = fake_windows_dir
        name = "exam_CMCse_1_NOMAD-pp.json"
        suffix = ".json"

    data.file_path = FakeWindowsPath()
    payload = data.createJsonDict()

    working_dir = payload["workingDirectory"]["data"]["0"]
    assert "\\" not in working_dir
    assert working_dir == "C:/Users/vsmw51/OneDrive - Durham University/projects/programming/2025/python/awh"


@pytest.mark.skipif(not BRUKER_GOLD.exists(), reason="bruker gold file not found next to this checkout")
def test_same_molecule_as_bruker_gold(nmrium_payload):
    with open(BRUKER_GOLD) as f:
        bruker = json.load(f)

    nmrium_smiles = Chem.CanonSmiles(nmrium_payload["smiles"]["data"]["0"])
    bruker_smiles = Chem.CanonSmiles(bruker["smiles"]["data"]["0"])
    assert nmrium_smiles == bruker_smiles


@pytest.mark.skipif(not BRUKER_GOLD.exists(), reason="bruker gold file not found next to this checkout")
def test_hsqc_carbon_shifts_agree_with_bruker_gold(nmrium_payload):
    with open(BRUKER_GOLD) as f:
        bruker = json.load(f)

    nmrium_shifts = sorted(v["delta1"] for v in nmrium_payload["HSQC_0"]["peaks"]["data"].values())
    bruker_shifts = sorted(v["delta1"] for v in bruker["HSQC_0"]["peaks"]["data"].values())

    assert len(nmrium_shifts) == len(bruker_shifts) == 11

    # Independent acquisitions of the same sample: allow generous
    # instrument/referencing variance, not exact match.
    TOLERANCE_PPM = 0.3
    for n, b in zip(nmrium_shifts, bruker_shifts):
        assert abs(n - b) < TOLERANCE_PPM, f"NMRium {n} vs Bruker {b} exceeds {TOLERANCE_PPM} ppm"
