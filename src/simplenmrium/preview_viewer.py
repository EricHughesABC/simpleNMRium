"""
preview_viewer.py

The pre-submission HTML preview window, deliberately run as a SEPARATE
OS process from simplenmrium_cli.py — NOT called in-process, for exactly
the reason documented in simplenmr_builder.gui.submission.
open_result_viewer_subprocess():

simplenmrium_cli.py already imports simplenmr_builder.gui.submission
(for check_user_registration/submit_to_server/open_result_viewer_
subprocess), which imports qtpy at module level. qtpy resolves to
whatever Qt binding is importable in the environment — confirmed
2026-09 on a real machine, it resolved to PyQt5 (Qt5) even though this
project never lists PyQt5 as a dependency, because something else in
that conda environment had it installed. This module needs PySide6
specifically for QWebEngineView (same requirement as html_viewer.py).
Constructing PySide6 widgets in a process that already has a PyQt5
QApplication crashes immediately: "QWidget: Must construct a
QApplication before a QWidget" — PySide6's QWidget does not recognize a
PyQt5-created QApplication as valid at all, even though both are
"a QApplication" conceptually. A fresh subprocess starts with nothing
Qt-related loaded, so there's nothing for PySide6 to collide with —
identical reasoning to why open_result_viewer_subprocess() exists
instead of importing html_viewer.MainWindow in-process.

Usage: python -m simplenmrium.preview_viewer <html_file_path>
Exit code: 0 if the person clicked Submit, 1 for Cancel or closing the
window (checked via subprocess.returncode by simplenmrium_cli.py — no
separate result file needed, matching html_viewer.py's own exit-code-
only convention).
"""

from __future__ import annotations

import os
import sys

from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import (
    QApplication,
    QHBoxLayout,
    QMainWindow,
    QPushButton,
    QVBoxLayout,
    QWidget,
)


class PreviewWindow(QMainWindow):
    def __init__(self, html_path: str):
        super().__init__()
        self.submitted = False

        self.setWindowTitle("simpleNMRium preview")
        self.resize(1000, 800)

        central = QWidget()
        layout = QVBoxLayout(central)

        self.view = QWebEngineView()
        with open(html_path, "r", encoding="utf-8") as f:
            html_str = f.read()
        # setHtml (not load(QUrl)) since this is generated content, not a
        # file we want QtWebEngine resolving relative links against.
        self.view.setHtml(html_str)
        layout.addWidget(self.view)

        button_row = QHBoxLayout()
        submit_button = QPushButton("Submit to simpleNMR server")
        cancel_button = QPushButton("Cancel")
        button_row.addStretch()
        button_row.addWidget(cancel_button)
        button_row.addWidget(submit_button)
        layout.addLayout(button_row)

        self.setCentralWidget(central)

        submit_button.clicked.connect(self._on_submit)
        cancel_button.clicked.connect(self._on_cancel)

    def _on_submit(self):
        self.submitted = True
        self.close()

    def _on_cancel(self):
        self.submitted = False
        self.close()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python -m simplenmrium.preview_viewer <path_to_html_file>")
        sys.exit(1)

    app = QApplication(sys.argv)
    window = PreviewWindow(sys.argv[1])
    window.show()
    app.exec()

    # os._exit(), not sys.exit(): same rationale as html_viewer.py —
    # QtWebEngine can leave background threads alive after app.exec()
    # returns, which would otherwise hang process shutdown waiting for
    # them. Nothing left to clean up once the window is closed.
    os._exit(0 if window.submitted else 1)
