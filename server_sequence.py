"""
AlphaGenome Explorer — Predict-from-sequence handlers.

Registers reactive handlers for the Sequence prediction feature.
Called from server.py to keep it manageable.
"""

import logging

import matplotlib.pyplot as plt
from shiny import reactive, render, ui

from alphagenome.data import genome
from alphagenome.visualization import plot_components

from config import OUTPUT_ATTR_MAP, OUTPUT_TYPES, SEQUENCE_LENGTHS
from download_utils import fig_to_png_bytes
from sequence_utils import get_sequence_info, pad_sequence, validate_sequence

logger = logging.getLogger(__name__)


def register(input, state):
    """Register sequence prediction handlers.

    Args:
        input: Shiny input object.
        state: dict of reactive values and helper functions:
            seq_output, _ensure_model, _organism, _get_ontology_terms.
    """
    seq_output = state["seq_output"]
    _ensure_model = state["_ensure_model"]
    _organism = state["_organism"]
    _get_ontology_terms = state["_get_ontology_terms"]

    @reactive.effect
    @reactive.event(input.btn_predict_seq)
    def _predict_sequence():
        raw = input.seq_input()
        seq, err = validate_sequence(raw)
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

        # Determine target length
        target_key = input.seq_pad_length()
        if target_key == "auto":
            padded, actual_len, pad_err = pad_sequence(seq)
        else:
            target_len = SEQUENCE_LENGTHS[target_key]
            padded, actual_len, pad_err = pad_sequence(seq, target_len)

        if pad_err:
            ui.notification_show(pad_err, type="error")
            return

        m = _ensure_model()
        req_outputs = {OUTPUT_TYPES[k] for k in selected_outputs}

        # Create a synthetic interval so TrackData objects have .interval set
        # (required by the plotting library). The API accepts an optional
        # interval= parameter and propagates it to all output TrackData.
        synth_interval = genome.Interval(chromosome="seq", start=0, end=actual_len)

        try:
            with ui.Progress(min=0, max=1) as p:
                p.set(0.1, message=f"Predicting from sequence ({actual_len:,} bp)...")
                result = m.predict_sequence(
                    sequence=padded,
                    requested_outputs=req_outputs,
                    ontology_terms=onto,
                    organism=_organism(),
                    interval=synth_interval,
                )
                p.set(1.0, message="Done!")
        except Exception as e:
            logger.error("Sequence prediction failed: %s", e)
            ui.notification_show(
                "Sequence prediction failed. Check your inputs and API key.",
                type="error", duration=8,
            )
            return

        info = get_sequence_info(seq)
        seq_output.set((result, selected_outputs, info, actual_len))
        ui.update_navs("main_tabs", selected="Sequence prediction")
        ui.notification_show("Sequence prediction complete!", type="message", duration=3)

    @render.ui
    def seq_info_ui():
        data = seq_output()
        if data is None:
            return ui.p("")
        _, _, info, actual_len = data
        return ui.div(
            ui.p(
                f"Input: {info['non_n_length']:,} bp "
                f"(GC: {info['gc_content']:.1f}%), "
                f"padded to {actual_len:,} bp",
                class_="small text-muted mb-0",
            ),
            class_="bg-light p-2 rounded mb-2",
        )

    @render.ui
    def plot_seq_ui():
        data = seq_output()
        if data is None:
            return ui.output_plot("plot_seq_tracks", height="300px")
        result, selected_outputs, _, _ = data
        n_tracks = 0
        for key in selected_outputs:
            attr = OUTPUT_ATTR_MAP.get(key)
            tdata = getattr(result, attr, None) if attr else None
            if tdata is not None and tdata.values.shape[-1] > 0:
                n_tracks += tdata.values.shape[-1]
        h = max(400, 120 * max(n_tracks, 1))
        return ui.output_plot("plot_seq_tracks", height=f"{h}px")

    @render.plot
    def plot_seq_tracks():
        data = seq_output()
        if data is None:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, "Paste a DNA sequence and click 'Predict from sequence'",
                    ha="center", va="center", fontsize=14)
            ax.axis("off")
            return fig

        result, selected_outputs, _, _ = data
        components = []

        for key in selected_outputs:
            attr = OUTPUT_ATTR_MAP.get(key)
            tdata = getattr(result, attr, None) if attr else None
            if tdata is None or tdata.values.shape[-1] == 0:
                continue

            if key == "CONTACT_MAPS":
                components.append(plot_components.ContactMaps(
                    tdata=tdata, ylabel_template=f"{key}: {{name}}",
                ))
            elif key == "SPLICE_JUNCTIONS":
                components.append(plot_components.Sashimi(
                    junction_track=tdata,
                    ylabel_template=f"{key}: {{name}} ({{strand}})",
                ))
            else:
                components.append(plot_components.Tracks(
                    tdata=tdata,
                    ylabel_template=f"{key}: {{name}} ({{strand}})\n{{name}}",
                    filled=True,
                ))

        if not components:
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5,
                    "No tracks returned for the selected output type + cell type.",
                    ha="center", va="center", fontsize=12)
            ax.axis("off")
            return fig

        first_attr = OUTPUT_ATTR_MAP.get(selected_outputs[0])
        first_tdata = getattr(result, first_attr, None) if first_attr else None
        interval = getattr(first_tdata, "interval", None) if first_tdata else None

        kwargs = {"title": "Predictions from raw sequence"}
        if interval is not None:
            kwargs["interval"] = interval

        fig = plot_components.plot(components, **kwargs)
        return fig

    # Download
    @render.download(filename="sequence_prediction.png")
    def dl_seq_png():
        data = seq_output()
        if data is None:
            return
        result, selected_outputs, _, _ = data
        components = []
        for key in selected_outputs:
            attr = OUTPUT_ATTR_MAP.get(key)
            tdata = getattr(result, attr, None) if attr else None
            if tdata is None or tdata.values.shape[-1] == 0:
                continue
            if key == "CONTACT_MAPS":
                components.append(plot_components.ContactMaps(
                    tdata=tdata, ylabel_template=f"{key}: {{name}}",
                ))
            elif key == "SPLICE_JUNCTIONS":
                components.append(plot_components.Sashimi(
                    junction_track=tdata,
                    ylabel_template=f"{key}: {{name}} ({{strand}})",
                ))
            else:
                components.append(plot_components.Tracks(
                    tdata=tdata,
                    ylabel_template=f"{key}: {{name}} ({{strand}})\n{{name}}",
                    filled=True,
                ))
        if not components:
            return
        first_attr = OUTPUT_ATTR_MAP.get(selected_outputs[0])
        first_tdata = getattr(result, first_attr, None) if first_attr else None
        interval = getattr(first_tdata, "interval", None) if first_tdata else None
        kwargs = {"title": "Predictions from raw sequence"}
        if interval is not None:
            kwargs["interval"] = interval
        fig = plot_components.plot(components, **kwargs)
        yield fig_to_png_bytes(fig)
