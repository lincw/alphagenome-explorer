"""
AlphaGenome Explorer — Download handlers.

Registers @render.download handlers for exporting data and plots.
Called from server.py to keep it manageable.
"""

from alphagenome.data import genome
from alphagenome.visualization import plot_components
from shiny import render

from download_utils import fig_to_png_bytes, interval_tracks_to_csv, score_df_to_csv
from ism_utils import build_ism_plot
from multi_variant_utils import build_comparison_bar, build_comparison_heatmap
from plot_utils import build_interval_components, build_variant_components
from score_utils import build_score_comparison_plot


def register(input, state):
    """Register all download handlers.

    Args:
        input: Shiny input object.
        state: dict of reactive values and helper functions:
            interval_output, variant_output, score_df, ism_data, multi_score_df,
            _ensure_gtf, transcript_ext.
    """
    interval_output = state["interval_output"]
    variant_output = state["variant_output"]
    score_df = state["score_df"]
    ism_data = state["ism_data"]
    multi_score_df = state["multi_score_df"]
    _ensure_gtf = state["_ensure_gtf"]
    transcript_ext = state["transcript_ext"]

    # -- Predicted tracks --
    @render.download(filename="predicted_tracks.csv")
    def dl_tracks_csv():
        data = interval_output()
        if data is None:
            yield ""
            return
        result, _, selected_outputs = data
        yield interval_tracks_to_csv(result, selected_outputs)

    @render.download(filename="predicted_tracks.png")
    def dl_tracks_png():
        data = interval_output()
        if data is None:
            return
        result, iv, selected_outputs = data
        _ensure_gtf()
        transcripts = transcript_ext().extract(iv)
        components = build_interval_components(result, selected_outputs, transcripts)
        fig = plot_components.plot(
            components, interval=iv,
            title=f"Predictions for {input.gene_symbol().strip()}",
        )
        yield fig_to_png_bytes(fig)

    # -- Variant effect --
    @render.download(filename="variant_effect.png")
    def dl_variant_png():
        data = variant_output()
        if data is None:
            return
        result, iv, selected_outputs, variant = data
        _ensure_gtf()
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
            title=f"Variant effect: {variant}",
        )
        yield fig_to_png_bytes(fig)

    # -- Variant scores --
    @render.download(filename="variant_scores.csv")
    def dl_scores_csv():
        df = score_df()
        yield score_df_to_csv(df) if df is not None else ""

    @render.download(filename="variant_scores.png")
    def dl_scores_png():
        df = score_df()
        if df is None:
            return
        fig = build_score_comparison_plot(df)
        yield fig_to_png_bytes(fig)

    # -- ISM --
    @render.download(filename="ism_seqlogo.png")
    def dl_ism_png():
        data = ism_data()
        if data is None:
            return
        ism_result, ism_interval, output_type_key, onto_label = data
        fig = build_ism_plot(ism_result, ism_interval, output_type_key, onto_label)
        yield fig_to_png_bytes(fig)

    # -- Multi-variant --
    @render.download(filename="multi_variant_scores.csv")
    def dl_multi_csv():
        df = multi_score_df()
        yield score_df_to_csv(df) if df is not None else ""

    @render.download(filename="multi_variant_heatmap.png")
    def dl_multi_heatmap_png():
        df = multi_score_df()
        if df is None:
            return
        try:
            target = input.gene_symbol().strip()
        except Exception:
            target = None
        fig = build_comparison_heatmap(df, target_gene=target)
        yield fig_to_png_bytes(fig)

    @render.download(filename="multi_variant_barplot.png")
    def dl_multi_bar_png():
        df = multi_score_df()
        if df is None:
            return
        try:
            target = input.gene_symbol().strip()
        except Exception:
            target = None
        fig = build_comparison_bar(df, target_gene=target)
        yield fig_to_png_bytes(fig)
