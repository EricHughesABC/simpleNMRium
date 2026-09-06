"""
html_report.py

Builds a standalone, self-contained HTML summary of what was extracted
from an NMRium file and what would be submitted to simpleNMR — a review
page shown before submission, not the post-submission D3 result page
that simplenmr_builder.gui.html_viewer handles.

Deliberately NOT built on html_viewer.py's QWebChannel export-back
machinery: that module exists to intercept a specific "Export" button on
the SERVER's result page and route it back into Python. This is a
read-only local preview with no export interaction, so a plain static
HTML string loaded via QWebEngineView.setHtml() is enough — no bridge,
no injected JS, no PySide6-specific permission-API probing needed beyond
what a bare QWebEngineView already provides.
"""

from __future__ import annotations

import html
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .json_converter import NMRiumData

_CSS = """
body { font-family: -apple-system, Segoe UI, Helvetica, Arial, sans-serif; margin: 2em; color: #1a1a1a; }
h1 { font-size: 1.4em; margin-bottom: 0.1em; }
h2 { font-size: 1.1em; margin-top: 1.6em; border-bottom: 1px solid #ddd; padding-bottom: 0.2em; }
.subtitle { color: #666; margin-top: 0; }
table { border-collapse: collapse; width: 100%; margin-top: 0.5em; }
th, td { text-align: left; padding: 4px 10px; border-bottom: 1px solid #eee; font-size: 0.92em; }
th { background: #f5f5f5; }
.ok { color: #157a15; font-weight: 600; }
.warn { color: #a15c00; font-weight: 600; }
.err { color: #b30000; font-weight: 600; }
.pill { display: inline-block; padding: 1px 8px; border-radius: 10px; font-size: 0.8em; margin-left: 6px; }
.pill.ok { background: #e5f5e5; }
.pill.warn { background: #fdf0d5; }
.mono { font-family: Menlo, Consolas, monospace; font-size: 0.85em; background: #f7f7f7; padding: 6px 10px; border-radius: 4px; word-break: break-all; }
"""


def _esc(s) -> str:
    return html.escape(str(s))


def build_html_report(nmrium_data: "NMRiumData", validation_ok: bool = True, validation_messages=None) -> str:
    """
    Render an HTML summary from an already-constructed NMRiumData
    instance (i.e. after __init__ has parsed the file and built
    .blocks, but independent of whether createJsonDict() has been
    called/succeeded yet — validation_ok/validation_messages are passed
    in separately so the caller can show a failed-validation report too).
    """
    ds = nmrium_data.dataset
    validation_messages = validation_messages or []

    rows = []
    for token, block, chosen_entry in nmrium_data.blocks:
        n_peaks = block["peaks"]["count"]
        rows.append(
            f"<tr><td>{_esc(token)}</td><td>{_esc(block['type'])}</td>"
            f"<td>{_esc(', '.join(block['nucleus']) if isinstance(block['nucleus'], list) else block['nucleus'])}</td>"
            f"<td>{n_peaks}</td><td class='mono'>{_esc(chosen_entry)}</td></tr>"
        )
    spectra_table = "\n".join(rows) or "<tr><td colspan='5'><em>none</em></td></tr>"

    skipped_rows = []
    for s in nmrium_data.skipped_spectra:
        skipped_rows.append(
            f"<tr><td>{_esc(s.experiment)}</td><td>{_esc(s.dimension)}D</td>"
            f"<td>{_esc(', '.join(s.nucleus))}</td></tr>"
        )
    skipped_table = "\n".join(skipped_rows)
    skipped_section = ""
    if skipped_rows:
        skipped_section = f"""
        <h2>Skipped spectra <span class="pill warn">no token mapping — not submitted</span></h2>
        <table>
          <tr><th>Experiment</th><th>Dim</th><th>Nucleus</th></tr>
          {skipped_table}
        </table>
        """

    atom_rows = []
    for rec in nmrium_data.all_atoms_info_records:
        atom_rows.append(
            f"<tr><td>{rec['atom_idx']}</td><td>{_esc(rec['atomNumber'])}</td>"
            f"<td>{_esc(rec['symbol'])}</td><td>{rec['numProtons']}</td></tr>"
        )
    atoms_table = "\n".join(atom_rows)

    if validation_ok:
        validation_html = '<p class="ok">✓ Payload passed schema + HSQC-presence validation.</p>'
    else:
        msgs = "".join(f"<li>{_esc(m)}</li>" for m in validation_messages)
        validation_html = f'<p class="err">✗ Payload FAILED validation:</p><ul>{msgs}</ul>'

    return f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>simpleNMRium preview — {_esc(ds.title or nmrium_data.file_path.name)}</title>
<style>{_CSS}</style>
</head>
<body>
  <h1>{_esc(ds.title or nmrium_data.file_path.name)}</h1>
  <p class="subtitle">{_esc(nmrium_data.file_path)}</p>

  <h2>Molecule</h2>
  <p>SMILES: <span class="mono">{_esc(nmrium_data.smiles)}</span></p>
  <p>{len(nmrium_data.all_atoms_info_records)} atoms
     ({len(nmrium_data.carbon_atoms_info_records)} carbons)</p>

  <h2>Spectra to be submitted</h2>
  <table>
    <tr><th>Token</th><th>Dim</th><th>Nucleus</th><th>Peaks/zones</th><th>chosenSpectra entry</th></tr>
    {spectra_table}
  </table>
  {skipped_section}

  <h2>Validation</h2>
  {validation_html}

  <h2>Atom numbering (RDKit order — used for allAtomsInfo/carbonAtomsInfo)</h2>
  <p class="warn">Note: NMRium's own correlation-table labels (C1, C2, ... in the
  source file) are NOT used here and are currently ignored — see project notes
  on the atom-numbering reconciliation problem.</p>
  <table>
    <tr><th>atom_idx</th><th>atomNumber</th><th>symbol</th><th>numProtons</th></tr>
    {atoms_table}
  </table>
</body>
</html>
"""
