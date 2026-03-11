"""
AlphaGenome Explorer — Batch variant scoring handlers.

Registers reactive handlers for Step 9: batch variant scoring with VCF input.
Called from server.py to keep it manageable.
"""

import logging

import pandas as pd
from shiny import reactive, render, ui

from batch_utils import (
    ALL_SCORER_KEYS,
    get_selected_scorers,
    parse_vcf_text,
    run_batch_scoring,
)
from config import SEQUENCE_LENGTHS
from download_utils import score_df_to_csv

logger = logging.getLogger(__name__)


def register(input, state):
    """Register batch scoring handlers.

    Args:
        input: Shiny input object.
        state: dict of reactive values and helper functions:
            batch_score_df, batch_errors, _ensure_model, _organism.
    """
    batch_score_df = state["batch_score_df"]
    batch_errors = state["batch_errors"]
    _ensure_model = state["_ensure_model"]
    _organism = state["_organism"]

    @reactive.effect
    @reactive.event(input.btn_batch_score)
    def _batch_score():
        text = input.batch_vcf_text().strip()
        if not text:
            ui.notification_show("Paste VCF data first.", type="error")
            return

        # Parse VCF input
        variants, parse_errors = parse_vcf_text(text)
        if not variants:
            batch_errors.set(parse_errors)
            ui.notification_show(
                "No valid variants parsed. " + (parse_errors[0] if parse_errors else ""),
                type="error", duration=8,
            )
            return

        # Get selected scorers
        selected_keys = list(input.batch_scorers())
        if not selected_keys:
            ui.notification_show("Select at least one scorer.", type="error")
            return

        organism = _organism()
        scorers, excluded = get_selected_scorers(selected_keys, organism)
        if not scorers:
            ui.notification_show(
                "No scorers available for this organism. " +
                (f"Excluded: {', '.join(excluded)}" if excluded else ""),
                type="error", duration=8,
            )
            return

        # Get sequence length
        seq_len = SEQUENCE_LENGTHS[input.batch_seq_length()]

        m = _ensure_model()

        try:
            with ui.Progress(min=0, max=1) as p:
                def progress_cb(frac, msg):
                    p.set(frac, message=msg)

                df, scoring_errors = run_batch_scoring(
                    model=m,
                    variants=variants,
                    scorers=scorers,
                    sequence_length=seq_len,
                    organism=organism,
                    progress_callback=progress_cb,
                )
                p.set(1.0, message="Done!")
        except Exception as e:
            logger.error("Batch scoring failed: %s", e)
            ui.notification_show(
                "Batch scoring failed. Check your inputs and API key.",
                type="error", duration=8,
            )
            return

        all_errors = parse_errors + scoring_errors
        if excluded:
            all_errors.append(f"Excluded scorers: {', '.join(excluded)}")
        batch_errors.set(all_errors)
        batch_score_df.set(df)
        ui.update_navs("main_tabs", selected="Batch scores")

        n_variants = df["variant_id"].nunique() if "variant_id" in df.columns else "?"
        n_scorers = len(scorers)
        ui.notification_show(
            f"Batch scoring complete: {n_variants} variants x {n_scorers} scorers "
            f"= {len(df):,} rows.",
            type="message", duration=5,
        )

    @render.ui
    def batch_errors_ui():
        errors = batch_errors()
        if not errors:
            return ui.p("")
        items = [ui.tags.li(e, class_="small") for e in errors[:10]]
        return ui.div(
            ui.p("Warnings / errors:", class_="small text-warning mb-1"),
            ui.tags.ul(*items, class_="mb-2"),
            class_="bg-warning bg-opacity-10 p-2 rounded mb-2",
        )

    @render.ui
    def batch_summary_ui():
        df = batch_score_df()
        if df is None:
            return ui.p("Paste VCF data and run batch scoring to see results.",
                        class_="text-muted small")

        n_variants = df["variant_id"].nunique() if "variant_id" in df.columns else "?"
        n_scorers = df["scorer"].nunique() if "scorer" in df.columns else "?"
        n_genes = df["gene_name"].nunique() if "gene_name" in df.columns else "N/A"
        n_tissues = df["biosample_name"].nunique() if "biosample_name" in df.columns else "?"

        return ui.div(
            ui.h6("Batch scoring summary"),
            ui.tags.ul(
                ui.tags.li(f"Variants: {n_variants}", class_="small"),
                ui.tags.li(f"Scorers: {n_scorers}", class_="small"),
                ui.tags.li(f"Genes: {n_genes}", class_="small"),
                ui.tags.li(f"Cell types / tissues: {n_tissues}", class_="small"),
                ui.tags.li(f"Total rows: {len(df):,}", class_="small"),
                class_="mb-0",
            ),
            class_="bg-light p-3 rounded mb-3",
        )

    @render.data_frame
    def table_batch_scores():
        df = batch_score_df()
        if df is None:
            return render.DataGrid(
                pd.DataFrame({"Info": ["Run batch scoring to see results here."]})
            )
        display_cols = [c for c in df.columns if c not in ("scored_interval",)]
        return render.DataGrid(df[display_cols].reset_index(drop=True), filters=True)

    # Download
    @render.download(filename="batch_variant_scores.csv")
    def dl_batch_csv():
        df = batch_score_df()
        yield score_df_to_csv(df) if df is not None else ""
