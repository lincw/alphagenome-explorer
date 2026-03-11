"""
AlphaGenome Explorer — Shiny for Python app
Run with:  shiny run app.py

Project structure:
    app.py          Entry point (this file)
    config.py       Constants and configuration
    ui_layout.py    UI layout definition
    server.py       Server logic (reactive handlers)
    plot_utils.py   Plot building utilities
    score_utils.py  Scoring display utilities
"""

import matplotlib
matplotlib.use("Agg")

from shiny import App

from server import server
from ui_layout import app_ui

app = App(app_ui, server)
