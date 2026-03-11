"""
AlphaGenome Explorer — Multi-variant comparison handlers.

Registers reactive handlers for Step 8: multi-variant scoring and comparison.
Called from server.py to keep it manageable.
"""

import logging

import matplotlib.pyplot as plt
import pandas as pd
from shiny import reactive, render, ui

from multi_variant_utils import (
    build_comparison_bar,
    build_comparison_heatmap,
    parse_variant_lines,
    score_multiple_variants,
)

logger = logging.getLogger(__name__)


def register(input, state):
    """Register multi-variant comparison handlers.

    Args:
        input: Shiny input object.
        state: dict of reactive values and helper functions:
            gene_interval, multi_score_df, multi_parse_errors,
            _ensure_model, _organism, _seq_len.
    """
    gene_interval = state["gene_interval"]
    multi_score_df = state["multi_score_df"]
    multi_parse_errors = state["multi_parse_errors"]
    _ensure_model = state["_ensure_model"]
    _organism = state["_organism"]
    _seq_len = state["_seq_len"]

    @reactive.effect
    @reactive.event(input.btn_multi_score)
    def _multi_score():
        iv = gene_interval()
        if iv is None:
            ui.notification_show("Look up a gene first (Step 1).", type="error")
            return

        text = input.multi_variants().strip()
        if not text:
            ui.notification_show("Enter at least one variant.", type="error")
            return

        parsed = parse_variant_lines(text, chromosome=iv.chromosome)
        errors = [p for p in parsed if p["error"]]
        valid = [p for p in parsed if p["variant"] is not None]

        multi_parse_errors.set(errors)

        if not valid:
            ui.notification_show(
                "No valid variants found. Check the format.",
                type="error", duration=8,
            )
            return

        m = _ensure_model()
        resized = iv.resize(_seq_len())
        scorer_key = input.multi_scorer_key()

        try:
            with ui.Progress(min=0, max=1) as p:
                def progress_cb(frac, msg):
                    p.set(frac, message=msg)

                df = score_multiple_variants(
                    model=m,
                    parsed_variants=valid,
                    interval=resized,
                    scorer_key=scorer_key,
                    organism=_organism(),
                    progress_callback=progress_cb,
                )
                p.set(1.0, message="Done!")
        except Exception as e:
            logger.error("Multi-variant scoring failed: %s", e)
            ui.notification_show(
                "Multi-variant scoring failed. Check your inputs and API key.",
                type="error", duration=8,
            )
            return

        multi_score_df.set(df)
        ui.update_navs("main_tabs", selected="Multi-variant")
        ui.notification_show(
            f"Scored {len(valid)} variant(s)!"
            + (f" ({len(errors)} had errors)" if errors else ""),
            type="message", duration=5,
        )

    @render.ui
    def multi_parse_status_ui():
        errors = multi_parse_errors()
        if not errors:
            return ui.p("")
        items = [
            ui.tags.li(f"{e['label']}: {e['error']}", class_="small")
            for e in errors
        ]
        return ui.div(
            ui.p("Some variants could not be parsed:", class_="small text-warning mb-1"),
            ui.tags.ul(*items, class_="mb-2"),
            class_="bg-warning bg-opacity-10 p-2 rounded mb-2",
        )

    @render.ui
    def plot_multi_heatmap_ui():
        df = multi_score_df()
        if df is None:
            return ui.p("Compare variants to see results here.", class_="text-muted small")
        return ui.output_plot("plot_multi_heatmap", height="600px")

    @render.plot
    def plot_multi_heatmap():
        df = multi_score_df()
        if df is None:
            fig, ax = plt.subplots()
            ax.axis("off")
            return fig
        try:
            target = input.gene_symbol().strip()
        except Exception:
            target = None
        return build_comparison_heatmap(df, target_gene=target)

    @render.ui
    def plot_multi_bar_ui():
        df = multi_score_df()
        if df is None:
            return ui.p("")
        return ui.output_plot("plot_multi_bar", height="500px")

    @render.plot
    def plot_multi_bar():
        df = multi_score_df()
        if df is None:
            fig, ax = plt.subplots()
            ax.axis("off")
            return fig
        try:
            target = input.gene_symbol().strip()
        except Exception:
            target = None
        return build_comparison_bar(df, target_gene=target)

    @render.data_frame
    def table_multi_scores():
        df = multi_score_df()
        if df is None:
            return render.DataGrid(
                pd.DataFrame({"Info": ["Compare variants to see results here."]})
            )
        display_cols = [c for c in df.columns if c not in ("scored_interval",)]
        return render.DataGrid(df[display_cols].reset_index(drop=True), filters=True)
