"""
AlphaGenome Explorer — Server logic.

All reactive handlers, input validation, and rendering functions.
"""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from shiny import reactive, render, ui

logger = logging.getLogger(__name__)

from alphagenome.data import gene_annotation, genome
from alphagenome.data import transcript as transcript_utils
from alphagenome.models import dna_client, variant_scorers
from alphagenome.visualization import plot_components

from config import OUTPUT_TYPES, SEQUENCE_LENGTHS
from plot_utils import (
    build_interval_components,
    build_variant_components,
    count_interval_tracks,
    count_variant_tracks,
)
from ism_utils import ISM_WIDTHS, build_ism_plot, extract_ism_matrix, run_ism
from score_utils import build_score_comparison_plot, build_score_summary, count_score_rows
from variant_utils import fetch_ref_sequence, lookup_rsid
import server_batch
import server_downloads
import server_multi
import server_sequence


def server(input, output, session):

    # ==========================================================================
    # Shared reactive values
    # ==========================================================================
    model = reactive.value(None)
    gtf = reactive.value(None)
    transcript_ext = reactive.value(None)
    gene_interval = reactive.value(None)
    interval_output = reactive.value(None)
    variant_output = reactive.value(None)
    score_df = reactive.value(None)
    ism_data = reactive.value(None)
    multi_score_df = reactive.value(None)
    multi_parse_errors = reactive.value([])
    batch_score_df = reactive.value(None)
    batch_errors = reactive.value([])
    seq_output = reactive.value(None)
    onto_choices = reactive.value({})

    # ==========================================================================
    # Lazy loaders & helpers
    # ==========================================================================

    current_api_key = reactive.value("")

    def _ensure_model():
        key = input.api_key().strip()
        if not key:
            raise RuntimeError("Enter your Google AI Studio API key first.")
        # Re-create model if key changed
        if key != current_api_key():
            m = dna_client.create(key)
            model.set(m)
            current_api_key.set(key)
        return model()

    @render.ui
    def api_key_status_ui():
        key = input.api_key().strip()
        if not key:
            return ui.p("Required to use the app.", class_="small text-warning mb-0")
        if not key.startswith("AIza") or len(key) != 39:
            return ui.p("Invalid key format. Google AI keys start with 'AIza' and are 39 characters.",
                        class_="small text-danger mb-0")
        return ui.p("Key format OK.", class_="small text-success mb-0")

    def _ensure_gtf():
        if gtf() is None:
            cache = Path(__file__).parent / "gencode_v46_annotation.feather"
            if cache.exists():
                g = pd.read_feather(cache)
            else:
                g = pd.read_feather(
                    "https://storage.googleapis.com/alphagenome/reference/gencode/"
                    "hg38/gencode.v46.annotation.gtf.gz.feather"
                )
                g.to_feather(cache)
            gtf.set(g)
            gt = gene_annotation.filter_protein_coding(g)
            gt = gene_annotation.filter_to_mane_select_transcript(gt)
            transcript_ext.set(transcript_utils.TranscriptExtractor(gt))
        return gtf()

    def _organism():
        return (
            dna_client.Organism.HOMO_SAPIENS
            if input.organism() == "human"
            else dna_client.Organism.MUS_MUSCULUS
        )

    def _seq_len():
        return SEQUENCE_LENGTHS[input.seq_length()]

    # ------------------------------------------------------------------
    # Variant input helpers (rsID lookup + auto-REF)
    # ------------------------------------------------------------------
    rsid_status = reactive.value("")
    _rsid_just_set = reactive.value(False)  # flag to prevent auto-fetch overwrite
    @reactive.effect
    @reactive.event(input.btn_lookup_rsid)
    def _on_lookup_rsid():
        rsid_str = input.var_rsid().strip()
        if not rsid_str:
            ui.notification_show("Enter an rsID first.", type="warning")
            return

        rsid_status.set("Looking up... (querying NCBI, may take a few seconds)")
        try:
            result = lookup_rsid(rsid_str)
        except Exception as e:
            rsid_status.set(f"Failed: {e}")
            ui.notification_show(str(e), type="error", duration=8)
            return

        # Check chromosome matches the gene
        iv = gene_interval()
        if iv and result["chrom"] != iv.chromosome:
            rsid_status.set(
                f"{rsid_str} is on {result['chrom']}, but gene is on {iv.chromosome}. "
                f"Check that the gene and variant match."
            )
            ui.notification_show(
                f"Chromosome mismatch: {rsid_str} is on {result['chrom']}, "
                f"gene is on {iv.chromosome}.",
                type="warning", duration=8,
            )

        # Fill in position, REF, ALT
        _rsid_just_set.set(True)
        ui.update_text("var_pos", value=str(result["pos"]))
        ui.update_text("var_ref", value=result["ref"])

        alts = result["alts"]
        if len(alts) == 1:
            ui.update_text("var_alt", value=alts[0]["alt"])
            vtype = alts[0]["variant_type"]
            rsid_status.set(
                f"{rsid_str}: {result['chrom']}:{result['pos']} "
                f"{result['ref']}>{alts[0]['alt']} ({vtype})"
            )
        else:
            # Multiple ALT alleles — use the first one, notify user
            ui.update_text("var_alt", value=alts[0]["alt"])
            alt_strs = [f"{a['alt']} ({a['variant_type']})" for a in alts]
            rsid_status.set(
                f"{rsid_str}: {result['chrom']}:{result['pos']} "
                f"REF={result['ref']}  |  ALTs: {', '.join(alt_strs)}  "
                f"(first ALT selected)"
            )
        ui.notification_show(
            f"Filled variant from {rsid_str}.", type="message", duration=3
        )

    @render.ui
    def rsid_status_ui():
        msg = rsid_status()
        if not msg:
            return ui.p("")
        return ui.p(msg, class_="small text-muted mb-1")

    # ------------------------------------------------------------------
    # Auto-fetch REF from hg38 when position is entered manually
    # ------------------------------------------------------------------
    @reactive.effect
    @reactive.event(input.var_pos)
    def _auto_fetch_ref():
        # Skip if position was just set by rsID lookup
        if _rsid_just_set():
            _rsid_just_set.set(False)
            return
        iv = gene_interval()
        if iv is None:
            return
        try:
            pos_str = str(input.var_pos()).strip().replace(",", "")
            if not pos_str:
                return
            pos = int(pos_str)
        except (ValueError, TypeError):
            return
        base = fetch_ref_sequence(iv.chromosome, pos - 1, pos)
        if base:
            ui.update_text("var_ref", value=base)

    def _parse_variant_fields():
        """Parse and validate variant input fields. Returns (variant, error_msg)."""
        iv = gene_interval()
        if iv is None:
            return None, "Look up a gene first (Step 1)."

        try:
            pos_str = str(input.var_pos()).strip()
        except Exception:
            pos_str = ""
        try:
            ref_str = str(input.var_ref()).strip()
        except Exception:
            ref_str = ""
        try:
            alt_str = str(input.var_alt()).strip()
        except Exception:
            alt_str = ""

        missing = []
        if not pos_str:
            missing.append("Position")
        if not ref_str:
            missing.append("REF (check position is valid)")
        if not alt_str:
            missing.append("ALT")
        if missing:
            return None, f"Missing: {', '.join(missing)}."

        try:
            pos = int(pos_str.replace(",", ""))
        except ValueError:
            return None, f"Position '{pos_str}' is not a valid number."

        variant = genome.Variant(
            chromosome=iv.chromosome,
            position=pos,
            reference_bases=ref_str,
            alternate_bases=alt_str,
        )
        return variant, None

    def _get_ontology_terms():
        """Safely get selected ontology terms."""
        try:
            val = input.ontology_terms()
            if val:
                return list(val)
        except Exception:
            pass
        return []

    # ==========================================================================
    # 1. Gene lookup
    # ==========================================================================

    @reactive.effect
    @reactive.event(input.btn_lookup)
    def _lookup_gene():
        g = _ensure_gtf()
        try:
            iv = gene_annotation.get_gene_interval(
                g, gene_symbol=input.gene_symbol().strip()
            )
            gene_interval.set(iv)
            ui.update_navs("main_tabs", selected="Gene context")
        except Exception as e:
            logger.error("Gene lookup failed: %s", e)
            ui.notification_show("Gene not found. Check the symbol and try again.",
                                 type="error", duration=5)

    @render.text
    def gene_info():
        iv = gene_interval()
        if iv is None:
            return "No gene loaded yet."
        return f"{iv.chromosome}:{iv.start:,}-{iv.end:,} ({iv.strand}) | width={iv.width:,} bp"

    @render.ui
    def var_gene_info_ui():
        iv = gene_interval()
        if iv is None:
            return ui.p("Look up a gene first (Step 1).", class_="text-warning small")
        return ui.div(
            ui.p(
                ui.strong(f"{input.gene_symbol().strip()} "),
                f"  {iv.chromosome}:{iv.start:,}-{iv.end:,} ({iv.strand})",
                class_="small mb-1",
            ),
            ui.p(
                f"Valid positions: {iv.start:,} to {iv.end:,}",
                class_="text-muted small mb-0",
            ),
            class_="bg-light p-2 rounded mb-2",
        )

    @render.plot
    def plot_gene_context():
        iv = gene_interval()
        if iv is None:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, "Click 'Look up gene' first",
                    ha="center", va="center", fontsize=14)
            ax.axis("off")
            return fig
        _ensure_gtf()
        transcripts = transcript_ext().extract(iv)
        fig = plot_components.plot(
            [plot_components.TranscriptAnnotation(transcripts)],
            interval=iv,
            title=f"Gene context: {input.gene_symbol().strip()}",
        )
        return fig

    # ==========================================================================
    # 3. Ontology / cell type search
    # ==========================================================================

    @reactive.effect
    @reactive.event(input.btn_search_onto)
    def _search_ontology():
        try:
            m = _ensure_model()
            keyword = input.ontology_search().strip()
            if not keyword:
                ui.notification_show("Enter a search keyword first.", type="warning")
                return
            meta = m.output_metadata(_organism()).concatenate()
            df = meta[["ontology_curie", "biosample_name", "output_type"]].drop_duplicates()
            mask = df["biosample_name"].str.contains(keyword, case=False, na=False, regex=False)
            matches = df[mask]
            if matches.empty:
                ui.notification_show(f"No matches for '{keyword}'.", type="warning")
                onto_choices.set({})
                return
            grouped = (
                matches.groupby(["ontology_curie", "biosample_name"])["output_type"]
                .apply(lambda x: ", ".join(sorted(str(v) for v in x)))
                .rename("output_types")
                .reset_index()
            )
            grouped = grouped.head(30)
            choices = {
                row.ontology_curie: f"{row.biosample_name} ({row.ontology_curie}) [{row.output_types}]"
                for row in grouped.itertuples()
            }
            onto_choices.set(choices)
        except Exception as e:
            logger.error("Ontology search failed: %s", e)
            ui.notification_show("Cell type search failed. Try a different keyword.",
                                 type="error", duration=8)

    @render.ui
    def ontology_select_ui():
        choices = onto_choices()
        if not choices:
            return ui.p("Search for a cell type above.", class_="text-muted small")
        return ui.input_select(
            "ontology_terms", "Select cell type(s)",
            choices=choices,
            multiple=True,
            selected=[list(choices.keys())[0]],
            size=min(len(choices), 6),
        )

    @render.data_frame
    def table_ontology():
        m = _ensure_model()
        meta = m.output_metadata(_organism()).concatenate()
        df = meta[["ontology_curie", "biosample_name", "output_type"]].drop_duplicates()
        df["output_type"] = df["output_type"].apply(lambda x: str(x).replace("OutputType.", ""))
        df["cell_type_tissue"] = df["biosample_name"] + " (" + df["ontology_curie"] + ")"
        df = df[["cell_type_tissue", "ontology_curie", "biosample_name", "output_type"]]
        return render.DataGrid(df.reset_index(drop=True), filters=True)

    # ==========================================================================
    # 4. Predict interval
    # ==========================================================================

    @reactive.effect
    @reactive.event(input.btn_predict)
    def _predict_interval():
        iv = gene_interval()
        if iv is None:
            ui.notification_show("Look up a gene first (Step 1).", type="error")
            return
        selected_outputs = input.output_types()
        if not selected_outputs:
            ui.notification_show("Select at least one output type (Step 2).", type="error")
            return
        onto = _get_ontology_terms()
        if not onto:
            ui.notification_show("Search and select a cell type first (Step 3).", type="error")
            return

        m = _ensure_model()
        req_outputs = {OUTPUT_TYPES[k] for k in selected_outputs}
        resized = iv.resize(_seq_len())

        # Warn about unavailable output types
        meta = m.output_metadata(_organism()).concatenate()
        available = [str(v) for v in meta[meta["ontology_curie"].isin(onto)]["output_type"].unique()]
        missing = [k for k in selected_outputs if k not in available]
        if missing:
            ui.notification_show(
                f"Warning: {', '.join(missing)} not available for the selected cell type. "
                f"Available: {', '.join(available)}",
                type="warning", duration=8,
            )

        try:
            with ui.Progress(min=0, max=1) as p:
                p.set(0.1, message="Running AlphaGenome prediction...")
                result = m.predict_interval(
                    interval=resized,
                    requested_outputs=req_outputs,
                    ontology_terms=onto,
                    organism=_organism(),
                )
                p.set(1.0, message="Done!")
        except Exception as e:
            logger.error("Interval prediction failed: %s", e)
            ui.notification_show("Prediction failed. Check your inputs and API key.",
                                 type="error", duration=8)
            return
        interval_output.set((result, iv, selected_outputs))
        ui.update_navs("main_tabs", selected="Predicted tracks")
        ui.notification_show("Interval prediction complete!", type="message", duration=3)

    @render.ui
    def plot_interval_ui():
        data = interval_output()
        n_tracks = count_interval_tracks(data)
        h = max(400, 120 * n_tracks)
        return ui.output_plot("plot_interval", height=f"{h}px")

    @render.plot
    def plot_interval():
        data = interval_output()
        if data is None:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, "Click 'Run prediction' to generate tracks",
                    ha="center", va="center", fontsize=14)
            ax.axis("off")
            return fig

        result, iv, selected_outputs = data
        _ensure_gtf()
        transcripts = transcript_ext().extract(iv)
        components = build_interval_components(result, selected_outputs, transcripts)

        if len(components) == 1:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5,
                    "No tracks returned for the selected output type + cell type.\n\n"
                    "Check the 'Available cell types' tab to find valid combinations.",
                    ha="center", va="center", fontsize=12, wrap=True)
            ax.axis("off")
            return fig

        fig = plot_components.plot(
            components, interval=iv,
            title=f"Predictions for {input.gene_symbol().strip()}",
        )
        return fig

    # ==========================================================================
    # 5. Variant effect
    # ==========================================================================

    @reactive.effect
    @reactive.event(input.btn_variant)
    def _predict_variant():
        variant, err = _parse_variant_fields()
        if err:
            ui.notification_show(err, type="error")
            return
        selected_outputs = input.output_types()
        if not selected_outputs:
            ui.notification_show("Select at least one output type (Step 2).", type="error")
            return
        onto = _get_ontology_terms()
        if not onto:
            ui.notification_show("Search and select a cell type first (Step 3).", type="error")
            return

        iv = gene_interval()
        m = _ensure_model()
        req_outputs = {OUTPUT_TYPES[k] for k in selected_outputs}
        resized = iv.resize(_seq_len())

        try:
            with ui.Progress(min=0, max=1) as p:
                p.set(0.1, message="Running variant prediction...")
                result = m.predict_variant(
                    interval=resized,
                    variant=variant,
                    requested_outputs=req_outputs,
                    ontology_terms=onto,
                    organism=_organism(),
                )
                p.set(1.0, message="Done!")
        except Exception as e:
            logger.error("Variant prediction failed: %s", e)
            ui.notification_show("Variant prediction failed. Check your inputs and API key.",
                                 type="error", duration=8)
            return
        variant_output.set((result, iv, selected_outputs, variant))
        ui.update_navs("main_tabs", selected="Variant effect")
        ui.notification_show("Variant prediction complete!", type="message", duration=3)

    @render.ui
    def plot_variant_ui():
        data = variant_output()
        n_tracks = count_variant_tracks(data)
        h = max(400, 120 * n_tracks)
        return ui.output_plot("plot_variant", height=f"{h}px")

    @render.plot
    def plot_variant():
        data = variant_output()
        if data is None:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, "Click 'Predict variant effect' to compare REF vs ALT",
                    ha="center", va="center", fontsize=14)
            ax.axis("off")
            return fig

        try:
            result, iv, selected_outputs, variant = data
            _ensure_gtf()

            # Determine display interval (zoom around variant)
            display_window = int(input.var_display_window())
            if display_window > 0:
                display_iv = genome.Interval(
                    chromosome=variant.chromosome,
                    start=variant.position,
                    end=variant.position,
                ).resize(display_window)
            else:
                display_iv = iv

            transcripts = transcript_ext().extract(display_iv)
            components = build_variant_components(result, selected_outputs, transcripts)

            fig = plot_components.plot(
                components, interval=display_iv,
                annotations=[plot_components.VariantAnnotation([variant])],
                title=(
                    f"Variant effect: {variant}  (display: {display_window // 1000} KB)"
                    if display_window > 0
                    else f"Variant effect: {variant}"
                ),
            )
            return fig

        except Exception as e:
            logger.error("Variant effect plot error: %s", e)
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, "Could not render variant effect plot.\n"
                    "Check your inputs and try again.",
                    ha="center", va="center", fontsize=11, color="red", wrap=True)
            ax.axis("off")
            return fig

    # ==========================================================================
    # 6. Variant scoring
    # ==========================================================================

    @reactive.effect
    @reactive.event(input.btn_score)
    def _score_variant():
        variant, err = _parse_variant_fields()
        if err:
            ui.notification_show(err, type="error")
            return

        iv = gene_interval()
        m = _ensure_model()
        resized = iv.resize(_seq_len())
        scorer_key = input.scorer_key()
        scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS.get(scorer_key)
        if scorer is None:
            ui.notification_show(f"Scorer '{scorer_key}' not found.", type="error")
            return

        try:
            with ui.Progress(min=0, max=1) as p:
                p.set(0.1, message="Scoring variant...")
                scores = m.score_variant(
                    interval=resized,
                    variant=variant,
                    variant_scorers=[scorer],
                    organism=_organism(),
                )
                p.set(1.0, message="Done!")
        except Exception as e:
            logger.error("Variant scoring failed: %s", e)
            ui.notification_show("Scoring failed. Check your inputs and API key.",
                                 type="error", duration=8)
            return

        df = variant_scorers.tidy_scores([scores], match_gene_strand=True)
        score_df.set(df)
        ui.update_navs("main_tabs", selected="Variant scores")
        ui.notification_show("Variant scoring complete!", type="message", duration=3)

    @render.ui
    def score_summary_ui():
        try:
            target = input.gene_symbol().strip()
        except Exception:
            target = None
        return build_score_summary(score_df(), target_gene=target)

    @render.ui
    def plot_scores_ui():
        df = score_df()
        if df is None or "raw_score" not in (df.columns if df is not None else []):
            return ui.p("")
        # Chart shows max 20 bars (top 10 + bottom 10)
        n_bars = min(count_score_rows(df), 20)
        h = max(400, 30 * n_bars + 150)
        return ui.output_plot("plot_scores_comparison", height=f"{h}px")

    @render.plot
    def plot_scores_comparison():
        return build_score_comparison_plot(score_df())

    @render.data_frame
    def table_scores():
        df = score_df()
        if df is None:
            return render.DataGrid(
                pd.DataFrame({"Info": ["Score a variant to see results here."]})
            )
        display_cols = [c for c in df.columns if c not in ("scored_interval",)]
        return render.DataGrid(df[display_cols].reset_index(drop=True), filters=True)

    # ==========================================================================
    # 7. In Silico Mutagenesis (ISM)
    # ==========================================================================

    @reactive.effect
    @reactive.event(input.btn_ism)
    def _run_ism():
        iv = gene_interval()
        if iv is None:
            ui.notification_show("Look up a gene first (Step 1).", type="error")
            return

        onto = _get_ontology_terms()
        if not onto:
            ui.notification_show("Search and select a cell type first (Step 3).", type="error")
            return

        # Determine center position
        center_str = input.ism_center().strip().replace(",", "")
        if not center_str:
            # Fall back to variant position if available
            try:
                center_str = str(input.var_pos()).strip().replace(",", "")
            except Exception:
                center_str = ""
        if not center_str:
            # Fall back to gene center
            center_str = str((iv.start + iv.end) // 2)

        try:
            center_pos = int(center_str)
        except ValueError:
            ui.notification_show(f"Invalid center position: '{center_str}'.", type="error")
            return

        ism_width = ISM_WIDTHS[input.ism_width()]
        output_type_key = input.ism_output_type()

        m = _ensure_model()

        try:
            with ui.Progress(min=0, max=1) as p:
                def progress_cb(frac, msg):
                    p.set(frac, message=msg)

                variant_scores, ism_interval, seq_interval = run_ism(
                    dna_model=m,
                    center_position=center_pos,
                    chromosome=iv.chromosome,
                    ism_width=ism_width,
                    output_type_key=output_type_key,
                    organism=_organism(),
                    progress_callback=progress_cb,
                )

                p.set(0.9, message="Building ISM matrix...")

                # Extract for the first selected ontology term
                ontology_curie = onto[0]
                ism_result = extract_ism_matrix(variant_scores, ontology_curie)

                # Get a human-readable label for the ontology
                try:
                    meta = m.output_metadata(_organism()).concatenate()
                    label_match = meta[meta["ontology_curie"] == ontology_curie]["biosample_name"]
                    onto_label = label_match.iloc[0] if len(label_match) > 0 else ontology_curie
                except Exception:
                    onto_label = ontology_curie

                p.set(1.0, message="Done!")

        except Exception as e:
            logger.error("ISM failed: %s", e)
            ui.notification_show(
                "ISM failed. Check your inputs and API key.", type="error", duration=8,
            )
            return

        ism_data.set((ism_result, ism_interval, output_type_key, onto_label))
        ui.update_navs("main_tabs", selected="ISM")
        ui.notification_show("ISM complete!", type="message", duration=3)

    @render.ui
    def plot_ism_ui():
        data = ism_data()
        if data is None:
            return ui.output_plot("plot_ism", height="300px")
        return ui.output_plot("plot_ism", height="400px")

    @render.plot
    def plot_ism():
        data = ism_data()
        if data is None:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, "Click 'Run ISM' to perform in silico mutagenesis",
                    ha="center", va="center", fontsize=14)
            ax.axis("off")
            return fig

        ism_result, ism_interval, output_type_key, onto_label = data
        return build_ism_plot(ism_result, ism_interval, output_type_key, onto_label)

    # ==========================================================================
    # 8. Multi-variant comparison (delegated to server_multi.py)
    # ==========================================================================
    server_multi.register(input, state={
        "gene_interval": gene_interval,
        "multi_score_df": multi_score_df,
        "multi_parse_errors": multi_parse_errors,
        "_ensure_model": _ensure_model,
        "_organism": _organism,
        "_seq_len": _seq_len,
    })

    # ==========================================================================
    # Downloads (delegated to server_downloads.py)
    # ==========================================================================
    # ==========================================================================
    # 9. Batch variant scoring (delegated to server_batch.py)
    # ==========================================================================
    server_batch.register(input, state={
        "batch_score_df": batch_score_df,
        "batch_errors": batch_errors,
        "_ensure_model": _ensure_model,
        "_organism": _organism,
    })

    # ==========================================================================
    # 10. Predict from sequence (delegated to server_sequence.py)
    # ==========================================================================
    server_sequence.register(input, state={
        "seq_output": seq_output,
        "_ensure_model": _ensure_model,
        "_organism": _organism,
        "_get_ontology_terms": _get_ontology_terms,
    })

    # ==========================================================================
    # Downloads (delegated to server_downloads.py)
    # ==========================================================================
    server_downloads.register(input, state={
        "interval_output": interval_output,
        "variant_output": variant_output,
        "score_df": score_df,
        "ism_data": ism_data,
        "multi_score_df": multi_score_df,
        "_ensure_gtf": _ensure_gtf,
        "transcript_ext": transcript_ext,
    })
