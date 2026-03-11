"""
AlphaGenome Explorer — Variant scoring display utilities.

Functions for building score summary UI and comparison plots.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from shiny import render, ui


def _build_gene_section(gene, rows, is_target=False):
    """Build a UI section for one gene's score rows."""
    # Sort rows within gene by absolute LFC descending
    rows.sort(key=lambda r: -abs(r.get("raw_score", 0)))

    row_elements = []
    for row in rows[:5]:  # Show top 5 cell types per gene
        cell = row.get("biosample_name", None) or row.get("track_name", "all")
        assay = row.get("Assay title", None) or row.get("output_type", "")
        strand = row.get("track_strand", "")
        lfc = row.get("raw_score", 0)
        fc = 2 ** lfc
        direction = "higher" if lfc > 0 else "lower"
        abs_pct = abs(fc - 1) * 100
        q = row.get("quantile_score", None)
        q_str = f"  |  quantile: {q:.3f}" if q is not None else ""
        color = "green" if lfc > 0 else "red" if lfc < 0 else "gray"

        row_elements.append(
            ui.p(
                f"  {cell} ({assay}, {strand}): ",
                ui.span(
                    ui.strong(f"ALT is {abs_pct:.1f}% {direction}"),
                    style=f"color: {color};",
                ),
                f" than REF  (LFC: {lfc:.4f}, FC: {fc:.4f}{q_str})",
                class_="mb-1 small ps-3",
            )
        )
    if len(rows) > 5:
        row_elements.append(
            ui.p(f"  ... and {len(rows) - 5} more cell types",
                 class_="mb-1 small text-muted ps-3")
        )

    # Gene header with top-level summary
    top_lfc = rows[0].get("raw_score", 0)
    top_dir = "higher" if top_lfc > 0 else "lower"
    top_pct = abs(2 ** top_lfc - 1) * 100
    target_badge = " (target gene)" if is_target else ""

    border = "border border-primary rounded" if is_target else ""
    bg = "bg-primary bg-opacity-10 p-2" if is_target else ""

    return ui.div(
        ui.h6(
            ui.strong(gene),
            ui.span(target_badge, class_="text-primary") if is_target else "",
            ui.span(
                f" — strongest effect: ALT is {top_pct:.1f}% {top_dir}",
                style=f"color: {'green' if top_lfc > 0 else 'red'};",
            ),
            class_="mb-1",
        ),
        *row_elements,
        class_=f"mb-3 {border} {bg}",
    )


def build_score_summary(df, target_gene=None):
    """Build a UI summary card from a scoring DataFrame, grouped by gene.

    Args:
        df: scoring DataFrame from variant_scorers.tidy_scores().
        target_gene: gene symbol from Step 1 to always show first.

    Returns:
        A shiny UI element with plain-English interpretation of scores.
    """
    if df is None:
        return ui.p("")
    if "raw_score" not in df.columns:
        return ui.p("No raw_score column found.", class_="text-muted")

    # Group rows by gene
    gene_groups = {}
    for _, row in df.iterrows():
        gene = row.get("gene_name", None) or row.get("output_type", "?")
        gene = str(gene) if gene else "?"
        if gene not in gene_groups:
            gene_groups[gene] = []
        gene_groups[gene].append(row)

    # Sort genes: largest absolute mean LFC first
    def _gene_sort_key(gene_name):
        rows = gene_groups[gene_name]
        return -max(abs(r.get("raw_score", 0)) for r in rows)
    sorted_genes = sorted(gene_groups.keys(), key=_gene_sort_key)

    # Always show target gene first if present
    target_key = None
    if target_gene:
        tg_upper = target_gene.strip().upper()
        for g in sorted_genes:
            if g.upper() == tg_upper:
                target_key = g
                break

    sections = []
    total_rows = 0

    # Target gene first (highlighted)
    if target_key:
        sorted_genes = [g for g in sorted_genes if g != target_key]
        rows = gene_groups[target_key]
        sections.append(_build_gene_section(target_key, rows, is_target=True))
        total_rows += min(len(rows), 5)
    elif target_gene:
        # Target gene not in results (region-based scorer)
        sections.append(
            ui.div(
                ui.p(
                    ui.strong(target_gene),
                    " — not found in scoring results. This scorer may not be gene-specific.",
                    class_="text-muted small",
                ),
                class_="mb-3 border border-secondary rounded p-2 bg-light",
            )
        )

    # Remaining genes
    for gene in sorted_genes:
        rows = gene_groups[gene]
        sections.append(_build_gene_section(gene, rows))
        total_rows += min(len(rows), 5)
        if total_rows >= 25:
            remaining = len(sorted_genes) - (len(sections) - (1 if target_key else 0))
            if remaining > 0:
                sections.append(
                    ui.p(f"... and {remaining} more genes", class_="text-muted small")
                )
            break

    return ui.div(
        ui.h6("Gene-level expression change summary"),
        *sections,
        class_="bg-light p-3 rounded mb-3",
    )


def build_score_comparison_plot(df):
    """Build a horizontal bar chart comparing variant effect across cell types.

    Returns:
        A matplotlib Figure.
    """
    if df is None or "raw_score" not in df.columns:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "Score a variant to see comparison", ha="center", va="center")
        ax.axis("off")
        return fig

    # Build labels from available columns
    label_parts = []
    for col in ["biosample_name", "track_name", "track_strand", "Assay title", "output_type"]:
        if col in df.columns:
            label_parts.append(col)
    keep_cols = ["gene_name", "raw_score", "quantile_score"] + label_parts
    keep_cols = [c for c in keep_cols if c in df.columns]
    plot_df = df[keep_cols].copy()

    name_col = (
        "biosample_name" if "biosample_name" in plot_df.columns
        else "track_name" if "track_name" in plot_df.columns
        else None
    )
    strand_col = "track_strand" if "track_strand" in plot_df.columns else None
    assay_col = (
        "Assay title" if "Assay title" in plot_df.columns
        else "output_type" if "output_type" in plot_df.columns
        else None
    )

    parts = []
    if name_col:
        parts.append(plot_df[name_col].fillna("").astype(str))
    if strand_col:
        parts.append(" (" + plot_df[strand_col].fillna("").astype(str) + ")")
    if assay_col:
        parts.append(" " + plot_df[assay_col].fillna("").astype(str))
    if parts:
        label = parts[0]
        for p in parts[1:]:
            label = label + p
        plot_df["label"] = label
    else:
        plot_df["label"] = "track"
    plot_df = plot_df.sort_values("raw_score", ascending=True)

    # Show only top 10 most increased + top 10 most decreased
    n_show = 10
    total = len(plot_df)
    if total > n_show * 2:
        top = plot_df.tail(n_show)
        bottom = plot_df.head(n_show)
        plot_df = pd.concat([bottom, top])
        subtitle = f"(showing top {n_show} increased + top {n_show} decreased out of {total})"
    else:
        subtitle = ""

    # Color by direction
    colors = ["#E63946" if x < 0 else "#2A9D8F" for x in plot_df["raw_score"]]

    fig, ax = plt.subplots(figsize=(10, max(4, len(plot_df) * 0.35 + 1.5)))
    y_pos = range(len(plot_df))
    ax.barh(y_pos, plot_df["raw_score"], color=colors, height=0.7, edgecolor="none")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_df["label"], fontsize=9)
    ax.axvline(0, color="black", linewidth=0.5)
    ax.set_xlabel("Log fold-change (ALT vs REF)", fontsize=10)

    gene_val = plot_df["gene_name"].iloc[0] if ("gene_name" in plot_df.columns and len(plot_df) > 0) else None
    gene = gene_val if (gene_val and str(gene_val) != "nan" and str(gene_val) != "None") else (
        df["output_type"].iloc[0] if "output_type" in df.columns else "variant"
    )
    variant_obj = df["variant_id"].iloc[0] if "variant_id" in df.columns else ""
    title = f"Variant effect on {gene} across cell types / tissues\n{variant_obj}"
    if subtitle:
        title += f"\n{subtitle}"
    ax.set_title(title, fontsize=11)

    # Add quantile annotations
    for i, (_, row) in enumerate(plot_df.iterrows()):
        q = row["quantile_score"]
        if not np.isnan(q):
            ax.text(
                row["raw_score"], i, f" q={q:.2f}",
                va="center",
                ha="left" if row["raw_score"] >= 0 else "right",
                fontsize=7, color="gray",
            )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return fig


def count_score_rows(df):
    """Count rows for dynamic plot height of score comparison chart."""
    if df is None or "raw_score" not in df.columns:
        return 0
    # Use available columns for counting unique entries
    if "biosample_name" in df.columns and "track_strand" in df.columns:
        return df[["biosample_name", "track_strand"]].drop_duplicates().shape[0]
    return len(df)
