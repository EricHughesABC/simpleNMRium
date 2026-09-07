# simpleNMRium

NMRium/NOMAD client for `simpleNMR` — converts an already peak-picked
NMRium export JSON file into the shared `simpleNMR` submission contract,
lets you confirm or override which detected spectra get submitted and
as what experiment type, previews exactly what will be sent, and submits
it to the simpleNMR server.

Sibling projects in the `simpleNMR` ecosystem:
[`simpleNMRbuilder`](https://github.com/EricHughesABC/simpleNMRbuilder) (shared contract/builder library — a required dependency, see below),
[`simpleNMRbrukerTools`](https://github.com/EricHughesABC/simpleNMRbrukerTools),
[`simpleNMRjeolTools3`](https://github.com/EricHughesABC/simpleNMRjeolTools3).

## What it does

1. Reads an NMRium/NOMAD export JSON file (molecule + peak-picked 1D/2D
   spectra) into a clean object model.
2. Auto-detects which `simpleNMR` experiment type each spectrum
   corresponds to (`H1_1D`, `C13_1D`, `COSY`, `HSQC`, `HMBC`, ...).
3. Shows a dialog — the same shared `SpectrumAssignmentDialog` the
   Bruker and JEOL clients use — so you can confirm, override, or SKIP
   each detected spectrum, and set the simulated-annealing/database-
   consent options, before anything is built.
4. Builds the `simpleNMR` submission payload via `simplenmr_builder
   .SimpleNMRBuilder`, which validates it against the real JSON schema
   and HSQC-presence rule before returning.
5. Shows an HTML preview of exactly what will be submitted, in a Qt
   window, with Submit/Cancel.
6. On Submit, sends it to the `simpleNMR` server and opens the results
   viewer.

## Installing

This package depends on `simpleNMRbuilder`, which must have `"nmrium"`
registered as a recognized submission source (`simplenmr_builder
.constants.SOURCES`) — this repo's own `pyproject.toml` already declares
that dependency, so a normal install pulls it in automatically:

```bash
git clone https://github.com/EricHughesABC/simpleNMRium.git
cd simpleNMRium
pip install -e ".[test]"
pytest   # should show all tests passing
```

### If you see `ValueError: unknown source 'nmrium', expected one of ('bruker', 'jeol')`

That means the `simpleNMRbuilder` copy pip installed doesn't have
`"nmrium"` registered yet — either the change hasn't been pushed to
`simpleNMRbuilder`'s `main` branch, or an older cached copy got
installed first. Fix it by installing `simpleNMRbuilder` explicitly,
editable, from a checkout that has the change, so it takes priority over
whatever pip fetched automatically:

```bash
git clone https://github.com/EricHughesABC/simpleNMRbuilder.git
cd simpleNMRbuilder
pip uninstall simpleNMRbuilder -y
pip install -e ".[gui,viewer]"
python -c "from simplenmr_builder.constants import SOURCES; print(SOURCES)"
# should print: ('bruker', 'jeol', 'nmrium')
```

Then re-run `pytest` in `simpleNMRium`.

### Qt binding note (macOS/conda environments especially)

This package forces `QT_API=pyside6` before any Qt import (see the top
of `simplenmrium_cli.py`) and runs the HTML preview window in a separate
OS process (`preview_viewer.py`). Both exist to work around a real,
confirmed crash (`QWidget: Must construct a QApplication before a
QWidget`) that happens if a conda/pip environment has more than one Qt
binding installed (e.g. an unrelated package pulling in PyQt5 alongside
this project's PySide6) — `qtpy`'s auto-detection isn't reliably
consistent across every call site once that happens. If you hit a
similar Qt crash somewhere new, the fix is almost always "make sure
PySide6 is the only thing constructing widgets in that process" —
either force `QT_API` earlier, or isolate the widget-constructing code
into its own subprocess like `preview_viewer.py` does.

## Usage

```python
from simplenmrium import NMRiumData

data = NMRiumData("export.json")
payload = data.createJsonDict()   # auto-detected spectrum types, no dialog
```

Or from the command line — opens a file picker if no path is given,
then the spectrum-selection dialog, then the HTML preview:

```bash
simplenmr-nmrium path/to/export.json
```

## Structure

```
src/simplenmrium/
    nmrium_reader.py       parses NMRium/NOMAD JSON into an object model
                             (Spectrum/Range1D/Zone2D/Signal/Peak)
    json_converter.py       NMRiumData — builds the simpleNMR payload via
                             simplenmr_builder.SimpleNMRBuilder; owns the
                             suggested-vs-confirmed spectrum assignment
                             split (see Design decisions below)
    html_report.py           builds the standalone HTML preview
    simplenmrium_cli.py       `simplenmr-nmrium` console script: file →
                             spectrum-selection dialog → convert →
                             preview → optional submit
    preview_viewer.py         the HTML preview window, run as a SEPARATE
                             process (see Qt binding note above)
tests/
    golden_files/             exam_CMCse_1_NOMAD-pp.json (real sample,
                             alpha-ionone, exam_CMCse_1)
    test_converter_against_bruker_gold.py
```

## Design decisions

- **Spectrum selection reuses the shared dialog.** `simplenmr_builder.gui
  .spectrum_assignment_dialog.SpectrumAssignmentDialog` — the exact same
  dialog Bruker/JEOL already use — lets you confirm, override, or SKIP
  each detected spectrum before anything is submitted.
  `json_converter.NMRiumData` computes automatic type suggestions
  (`.suggested_assignments`) from `map_experiment_to_token()`, but
  nothing is committed to the payload until `.apply_assignments(...)` is
  called with the CONFIRMED list — either from the dialog (interactive
  CLI use) or, for non-interactive/test use, a lazy fallback that
  auto-applies the suggestions if no dialog was ever shown.
- **`correlations` block is ignored for v1.** NMRium's own correlation
  labels (`"C1"`, `"C2"`, ... `"H1"`, `"H11"` in
  `correlations.values[i].label.origin`) do NOT correspond to the
  molfile's RDKit atom order — the molfile's own embedded custom-label
  SGROUP turned out to be plain sequential `1..N`, unrelated to the
  correlation labels (which skip the carbonyl carbon and renumber from
  the HSQC-correlated carbons only). `allAtomsInfo`/`carbonAtomsInfo`
  are built purely from RDKit atom order here, exactly like the
  Bruker/JEOL converters — no atom assignments are submitted
  (`nmrAssignments` is forbidden for this source, same as JEOL).
  **Writing simpleNMR's solved assignments back into the NMRium JSON via
  these correlation labels is real follow-up work**, deliberately
  deferred — needs its own reconciliation approach, not attempted here.
- **`MNOVAcalcMethod` = `"NMRSHIFTDB2 Predict"`, `carbonCalcPositionsMethod`
  = `"Calculated Positions"`** — copied verbatim from both the real
  Bruker and real JEOL converters. This value makes the server compute
  expected carbon shifts itself (`prediction_from_nmrshiftdb2()` per the
  contract manifest, confirmed by a real server log:
  `"Using nmrshiftDB2 for C13 predictions"`). `c13predictions` IS still
  submitted, as a **present-but-empty** envelope (`count: 0, data: {}`),
  not omitted — confirmed necessary by a real server crash: the server's
  simulated-annealing code path reads `json_data["oldjsondata"]
  ["c13predictions"]["data"]` unconditionally whenever
  `simulatedAnnealing=True`, and raises `KeyError` if the field is
  missing entirely rather than present-and-empty. Matches the real
  `bruker_real_exam_cmcse.json` golden file's shape exactly.
- **2D zones collapse to one averaged peak per zone** (`Zone2D
  .average_x/average_y/average_intensity`), not one row per raw
  sub-peak. Revisit if real server submissions show this doesn't solve
  well for multi-peak zones (e.g. the 6-peak COSY zone in the sample
  file).

## A real technical gotcha, fixed here: NMRium's molfile SGROUP block

NMRium/OpenChemLib V3000 molfiles carry a `SGROUP` block (custom atom
labels via `NOSEARCH_OCL_CUSTOM_LABEL`) that RDKit's `MolFromMolBlock`
cannot parse at all (`"BEGIN SGROUP found but Sgroups NOT expected"`,
confirmed against the real uploaded sample). `json_converter.
strip_sgroup_block()` removes it before parsing — confirmed the labels
it carries are just plain sequential `1..N` in the one real file tested,
so nothing is lost by discarding it. This will affect **any** NMRium
export, not just this one file, so it's handled generically rather than
patched around for this sample specifically.

## Cross-validated against real data

`tests/test_converter_against_bruker_gold.py` checks this converter's
output for `exam_CMCse_1_NOMAD-pp.json` (alpha-ionone) against
`simpleNMRbuilder`'s own `bruker_real_exam_cmcse.json` golden file — a
real, independently-produced submission of the exact same physical
sample from the Bruker/TopSpin side. Same molecule (confirmed via RDKit
canonical SMILES), same carbon count (13), and HSQC ¹³C shifts agree to
within ~0.3 ppm — real evidence the conversion logic is sound, not just
a schema-shaped guess.

**Also confirmed end-to-end against a real running `simpleNMR` server**
(2026-09-06): the full flow — spectrum-selection dialog → payload build
→ validation → submission → registration check → server-side NMRSHIFTDB2
prediction → HSQC indexing → symmetry reconciliation → CH assignment →
HMBC graph construction → simulated annealing — completed successfully
for the sample file.

## What this deliberately does not do (v1)

- No override for `simulatedAnnealing`/`ml_consent` beyond what the
  spectrum-selection dialog's own checkboxes provide (matches Bruker/JEOL
  — those are the only processing-option controls the shared dialog
  exposes).
- No handling of experiment types not seen in any real NMRium file yet
  (DEPT, HSQC-CLIPCOSY-style experiments) beyond a best-guess mapping in
  `json_converter.map_experiment_to_token` — flagged in its docstring as
  unconfirmed, not silently trusted. Only one real NMRium file has been
  tested against so far (`exam_CMCse_1_NOMAD-pp.json`) — more real
  samples, especially with different experiment types or multi-block
  spectra, would be the most valuable next input.
- Does not attempt the `correlations`-label reconciliation problem (see
  Design decisions above) — writing solved assignments back into the
  NMRium JSON is separate, deferred work.

## Testing

```bash
pytest
```

9 tests: schema/HSQC-presence validation, all expected blocks present,
carbon count, known zone counts (from the sample file's own
zone-averaging CSVs), `c13predictions` shape regression, auto-detected
assignment coverage, user-SKIP override behavior, and the two
cross-checks against the real Bruker golden file described above.
