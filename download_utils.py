"""
AlphaGenome Explorer — Export / download utilities.

Helpers for converting prediction results and figures to downloadable formats.
"""

import io

import pandas as pd

from config import OUTPUT_ATTR_MAP


def interval_tracks_to_csv(result, selected_outputs):
    """Convert interval prediction tracks to a tidy CSV string.

    Each row: output_type, track_name, strand, position, value.
    To keep file size manageable, values are sampled at 128bp resolution
    for sequences longer than 100KB.

    Returns:
        CSV string.
    """
    rows = []
    for key in selected_outputs:
        attr = OUTPUT_ATTR_MAP.get(key)
        tdata = getattr(result, attr, None) if attr else None
        if tdata is None or tdata.values.shape[-1] == 0:
            continue

        n_positions = tdata.values.shape[0]
        # Downsample for large sequences
        step = 128 if n_positions > 100_000 else 1
        interval = tdata.interval
        positions = range(interval.start, interval.end, step)

        for track_idx in range(tdata.values.shape[-1]):
            name = tdata.metadata.iloc[track_idx].get("name", f"track_{track_idx}")
            strand = tdata.metadata.iloc[track_idx].get("strand", ".")
            values = tdata.values[::step, track_idx]
            for pos, val in zip(positions, values):
                rows.append({
                    "output_type": key,
                    "track_name": name,
                    "strand": strand,
                    "position": pos,
                    "value": float(val),
                })

    df = pd.DataFrame(rows)
    return df.to_csv(index=False)


def score_df_to_csv(df):
    """Convert score DataFrame to CSV, dropping internal columns."""
    if df is None:
        return ""
    drop_cols = [c for c in ["scored_interval"] if c in df.columns]
    return df.drop(columns=drop_cols).to_csv(index=False)


def fig_to_png_bytes(fig, dpi=150):
    """Save a matplotlib Figure to PNG bytes."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    buf.seek(0)
    return buf.read()
