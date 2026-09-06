"""
simplenmrium_cli.py

The runnable program (installed as the `simplenmr-nmrium` console script)
that resolves an NMRium export JSON file, converts it via json_converter.
NMRiumData, lets the person choose which detected spectra to submit and
under what experiment type via the shared spectrum-assignment dialog,
shows an HTML preview of what would be submitted in a Qt window, and —
if the person confirms — submits it to the simpleNMR server and opens
the results viewer.

Structured the same way as simpleNMRjeolTools3's jason_simpleNMR_cli.py
(orchestration/CLI only; all real conversion logic lives in
json_converter.py). Spectrum selection reuses simplenmr_builder.gui.
spectrum_assignment_dialog.SpectrumAssignmentDialog — the same dialog
Bruker/JEOL already use — pre-populated with json_converter.NMRiumData's
automatic experiment-type guesses via its suggested_assignments, exactly
mirroring jeolData.choosePeakPickedSpectaforSimpleNMR()'s dialog-then-
apply flow.
"""

from __future__ import annotations

import os

# MUST be set before any qtpy-touching import (including simplenmr_builder
# .gui.submission below, and json_converter's lazy machine_id() import) —
# qtpy auto-detects a Qt binding from whatever's importable in the
# environment, which is non-deterministic once more than one binding is
# installed. Confirmed on a real machine (2026-09): a conda env with an
# unrelated PyQt5 install present caused qtpy to resolve inconsistently,
# producing "QWidget: Must construct a QApplication before a QWidget"
# when a later qtpy-sourced widget (submission.py's QProgressDialog)
# didn't recognize an earlier qtpy-sourced QApplication as valid — the
# two calls had silently resolved to different underlying bindings.
# Forcing QT_API here makes every qtpy import in this process resolve to
# PySide6 deterministically, matching what preview_viewer.py already
# hardcodes directly (it doesn't use qtpy at all) — one binding for the
# whole process, regardless of what else happens to be pip/conda-
# installed alongside it.
os.environ.setdefault("QT_API", "pyside6")

import sys
from pathlib import Path
from typing import Optional

import fire

from .html_report import build_html_report
from .json_converter import NMRiumData

try:
    from simplenmr_builder import ContractError
    from simplenmr_builder.gui.submission import (
        SubmissionOutcome,
        check_user_registration,
        open_result_viewer_subprocess,
        submit_to_server,
    )
except ImportError as e:
    print(
        "ERROR: simplenmr_builder[gui,viewer] is not installed in this "
        "environment. Install it with:\n"
        '    pip install -e "<path-to-simpleNMRbuilder>[gui,viewer]"\n'
        f"\nUnderlying import error: {e}"
    )
    sys.exit(1)


def init_qt_app():
    from qtpy.QtWidgets import QApplication

    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    return app


def get_file_dialog():
    from qtpy.QtWidgets import QFileDialog

    file_path, _ = QFileDialog.getOpenFileName(
        None,
        "Select NMRium export JSON file",
        str(Path.home()),
        "NMRium/NOMAD files (*.json);;All Files (*)",
    )
    return Path(file_path) if file_path else None


def show_info_message(title: str, message: str, message_type=None):
    from qtpy.QtWidgets import QMessageBox

    if message_type is None:
        message_type = QMessageBox.Information
    msg_box = QMessageBox()
    msg_box.setIcon(message_type)
    msg_box.setWindowTitle(title)
    msg_box.setText(message)
    msg_box.exec_()


def validate_file(fn) -> tuple[bool, str]:
    if not fn:
        return False, "No file provided"
    fn = Path(fn)
    if not fn.exists():
        return False, f"File does not exist: {fn}"
    if fn.suffix.lower() != ".json":
        return False, f"Invalid file extension. Expected .json, got {fn.suffix}"
    return True, "File is valid"


def commandline(fn=None) -> Path:
    """Resolution order: explicit fn argument, then a file picker."""
    if fn:
        fn = Path(fn)

    is_valid, message = validate_file(fn)
    if not is_valid:
        print(f"Invalid input ({message}): opening file dialog...")
        fn = get_file_dialog()
        if not fn:
            show_info_message("No File Selected", "No file selected. Exiting.")
            sys.exit(1)
        is_valid, message = validate_file(fn)
        if not is_valid:
            show_info_message("File Validation Error", f"Selected file is invalid: {message}")
            sys.exit(1)

    return fn


def show_preview_window(html_str: str, title: str) -> bool:
    """
    Show the HTML preview with Submit/Cancel buttons, and return True if
    the person clicked Submit, False otherwise (Cancel or closing the
    window).

    Deliberately launched as a SEPARATE OS process (preview_viewer.py)
    rather than constructing PySide6 widgets in-process here — see
    preview_viewer.py's module docstring for the exact crash this avoids
    ("QWidget: Must construct a QApplication before a QWidget", from this
    process already holding a qtpy-resolved PyQt5 QApplication via
    simplenmr_builder.gui.submission's own imports). Same reasoning and
    same pattern as simplenmr_builder.gui.submission.
    open_result_viewer_subprocess() uses for the post-submission viewer.
    """
    import subprocess
    import tempfile

    with tempfile.NamedTemporaryFile(mode="w", suffix=".html", delete=False, encoding="utf-8") as f:
        f.write(html_str)
        html_path = f.name

    try:
        args = [sys.executable, "-m", "simplenmrium.preview_viewer", html_path]
        result = subprocess.run(args)
        return result.returncode == 0
    finally:
        try:
            Path(html_path).unlink()
        except OSError:
            pass


def show_spectrum_selection_dialog(nmrium_data: "NMRiumData"):
    """
    Show the shared SpectrumAssignmentDialog (simplenmr_builder.gui.
    spectrum_assignment_dialog — the same dialog Bruker/JEOL use),
    pre-populated with nmrium_data.suggested_assignments, so the person
    can confirm, override, or SKIP each detected spectrum before
    anything is submitted.

    Returns (assignments, simulated_annealing, ml_consent) if accepted,
    or None if the person cancelled. Uses qtpy (not PySide6 directly) —
    this is a plain QDialog, not QtWebEngine, so it's safe to run
    in-process alongside the rest of this file's qtpy usage (see the
    QT_API note at the top of this module); it does NOT need the
    subprocess isolation show_preview_window() requires.
    """
    from qtpy.QtWidgets import QDialog
    from simplenmr_builder.gui.spectrum_assignment_dialog import SpectrumAssignmentDialog

    entries = [a["experiment_name"] for a in nmrium_data.suggested_assignments]
    chosen_types = {a["experiment_name"]: a["experiment_type"] for a in nmrium_data.suggested_assignments}

    dialog = SpectrumAssignmentDialog(
        entries=entries,
        chosen_types=chosen_types,
        window_title="Choose which NMR experiments to submit",
    )

    if dialog.exec_() != QDialog.Accepted:
        return None

    assignments = dialog.get_assignments()
    simulated_annealing, ml_consent = dialog.get_processing_options()
    return assignments, simulated_annealing, ml_consent


def main() -> int:
    """Entry point installed as the `simplenmr-nmrium` console script."""

    local_remote = "https://test-simplenmr.pythonanywhere.com"
    simpleNMR_address = f"{local_remote}/simpleMNOVA"
    ml_address = f"{local_remote}/check_machine_learning"

    init_qt_app()

    print("SERVER ADDRESS:", local_remote)

    nmrium_fn = fire.Fire(commandline)

    if not nmrium_fn.exists():
        print("File not found:", nmrium_fn)
        show_info_message("No File Found", f"File not found: {nmrium_fn}")
        return 1

    try:
        nmrium_data = NMRiumData(nmrium_fn)
    except Exception as e:
        print(f"ERROR reading/parsing NMRium file: {e}")
        show_info_message("Read Error", f"Could not read the NMRium file:\n\n{e}")
        return 1

    dialog_result = show_spectrum_selection_dialog(nmrium_data)
    if dialog_result is None:
        print("Spectrum selection was cancelled by the person. Exiting.")
        return 0
    assignments, simulated_annealing, ml_consent = dialog_result
    nmrium_data.apply_assignments(assignments)

    validation_ok = True
    validation_messages: list[str] = []
    payload: Optional[dict] = None
    try:
        payload = nmrium_data.createJsonDict(
            simulated_annealing=simulated_annealing,
            ml_consent=ml_consent,
        )
    except ContractError as e:
        validation_ok = False
        validation_messages = [str(e)]
        print(f"\nERROR: simplenmr_builder rejected this submission before it was built: {e}\n")

    html_str = build_html_report(nmrium_data, validation_ok=validation_ok, validation_messages=validation_messages)

    # Always show the preview, even on a validation failure — the person
    # should be able to see exactly what was extracted and why it failed,
    # not just get a dialog box with no detail.
    should_submit = show_preview_window(html_str, title=f"simpleNMRium preview — {nmrium_fn.name}")

    if not validation_ok:
        show_info_message(
            "Submission Rejected",
            f"The data could not be validated for submission:\n\n{validation_messages[0]}",
        )
        return 1

    # Save the converted JSON locally regardless of whether the person
    # goes on to submit — matches Bruker/JEOL precedent of always
    # persisting the converted payload next to the source file.
    json_file_path = nmrium_fn.parent / nmrium_fn.name.replace(nmrium_fn.suffix, "_nmrium_inputdata.json")
    import json as _json

    with open(json_file_path, "w") as f:
        _json.dump(payload, f, indent=4)
    print(f"Converted payload saved to: {json_file_path}")

    if not should_submit:
        print("Person chose not to submit. Exiting.")
        return 0

    if not check_user_registration(ml_address):
        show_info_message("Registration Error", "Unable to verify registration.")
        return 1

    print("\nSubmitting to simpleNMR Server...")
    submission = submit_to_server(payload, simpleNMR_address)

    if submission.outcome in (SubmissionOutcome.SUCCESS, SubmissionOutcome.DIAGNOSTIC_HTML):
        print("Opening results viewer...")
        open_result_viewer_subprocess(submission, wait=True)
        return 0

    print(f"Server submission did not succeed: {submission.outcome.value}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
