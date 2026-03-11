"""
AlphaGenome Explorer — Multi-variant comparison utilities.

Parse variant lists, score in batch, and build comparison visualizations.
"""

import logging
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from alphagenome.data import genome
from alphagenome.models import variant_scorers

from variant_utils import fetch_ref_sequence, lookup_rsid

logger = logging.getLogger(__name__)

MAX_VARIANTS = 10


def parse_variant_lines(text, chromosome=None):
    """Parse multi-line variant input.

    Supported formats (one per line):
        chr12:112919388:G:A  disease_snp
        rs35482426  risk_variant
        112919388:G:A  label          (uses chromosome from gene lookup)
        # comment lines are ignored

    Args:
        text: Multi-line string of variants.
        chromosome: Default chromosome from gene lookup (e.g. "chr12").

    Returns:
        List of dicts with keys: variant (genome.Variant or None),
        label (str), error (str or None).
    """
    results = []
    for i, line in enumerate(text.strip().splitlines(), 1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if len(results) >= MAX_VARIANTS:
            results.append({
                "variant": None,
                "label": f"line {i}",
                "error": f"Maximum {MAX_VARIANTS} variants allowed.",
            })
            break

        # Split line into variant spec + optional label
        parts = line.split(None, 1)
        spec = parts[0]
        label = parts[1] if len(parts) > 1 else spec

        result = _parse_single_variant(spec, chromosome, label, i)
        results.append(result)

    return results


def _parse_single_variant(spec, chromosome, label, line_num):
    """Parse a single variant specification string."""
    # rsID format
    if re.match(r"^rs\d+$", spec, re.IGNORECASE):
        try:
            info = lookup_rsid(spec)
            alts = info["alts"]
            variant = genome.Variant(
                chromosome=info["chrom"],
                position=info["pos"],
                reference_bases=info["ref"],
                alternate_bases=alts[0]["alt"],
                name=spec,
            )
            return {"variant": variant, "label": label, "error": None}
        except Exception as e:
            return {"variant": None, "label": label, "error": f"rsID lookup failed: {e}"}

    # Colon-separated: chr:pos:ref:alt or pos:ref:alt
    colon_parts = spec.split(":")
    if len(colon_parts) == 4:
        chrom, pos_str, ref, alt = colon_parts
        try:
            pos = int(pos_str.replace(",", ""))
            variant = genome.Variant(
                chromosome=chrom, position=pos,
                reference_bases=ref, alternate_bases=alt, name=label,
            )
            return {"variant": variant, "label": label, "error": None}
        except (ValueError, TypeError) as e:
            return {"variant": None, "label": label, "error": str(e)}

    if len(colon_parts) == 3:
        # pos:ref:alt — use default chromosome
        pos_str, ref, alt = colon_parts
        if not chromosome:
            return {"variant": None, "label": label,
                    "error": "No chromosome — use chr:pos:ref:alt format or look up a gene first."}
        try:
            pos = int(pos_str.replace(",", ""))
            variant = genome.Variant(
                chromosome=chromosome, position=pos,
                reference_bases=ref, alternate_bases=alt, name=label,
            )
            return {"variant": variant, "label": label, "error": None}
        except (ValueError, TypeError) as e:
            return {"variant": None, "label": label, "error": str(e)}

    return {"variant": None, "label": label,
            "error": f"Unrecognized format on line {line_num}. Use chr:pos:ref:alt or rsID."}


def score_multiple_variants(model, parsed_variants, interval, scorer_key,
                            organism, progress_callback=None):
    """Score multiple variants and return a combined DataFrame.

    Args:
        model: AlphaGenome model instance.
        parsed_variants: List of dicts from parse_variant_lines (valid ones only).
        interval: Resized gene interval for scoring.
        scorer_key: Key from RECOMMENDED_VARIANT_SCORERS.
        organism: dna_client.Organism.
        progress_callback: Optional callable(fraction, message).

    Returns:
        pd.DataFrame with an added "variant_label" column.
    """
    scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS.get(scorer_key)
    if scorer is None:
        raise ValueError(f"Scorer '{scorer_key}' not found.")

    all_results = []
    labels = []
    n = len(parsed_variants)

    for i, pv in enumerate(parsed_variants):
        variant = pv["variant"]
        label = pv["label"]

        if progress_callback:
            progress_callback(
                (i + 0.5) / n,
                f"Scoring variant {i + 1}/{n}: {label}...",
            )

        try:
            scores = model.score_variant(
                interval=interval,
                variant=variant,
                variant_scorers=[scorer],
                organism=organism,
            )
            all_results.append(scores)
            labels.append(label)
        except Exception as e:
            logger.error("Failed to score variant %s: %s", label, e)
            # Skip this variant but continue with others

    if not all_results:
        raise RuntimeError("All variant scoring attempts failed.")

    df = variant_scorers.tidy_scores(all_results, match_gene_strand=True)

    # Add variant_label column based on variant_id mapping
    variant_ids = []
    for res in all_results:
        # Each result is a list (one per scorer); get variant_id from first
        vid = res[0].uns.get("variant", "?") if hasattr(res[0], "uns") else "?"
        variant_ids.append(str(vid))

    label_map = dict(zip(variant_ids, labels))
    df["variant_label"] = df["variant_id"].map(label_map).fillna(df["variant_id"])

    return df


def build_comparison_heatmap(df, target_gene=None):
    """Build a heatmap: cell types (rows) x variants (columns), colored by LFC.

    Focuses on the target gene (or gene with strongest effect).

    Returns:
        matplotlib Figure.
    """
    if df is None or df.empty:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        ax.axis("off")
        return fig

    # Pick the gene to display
    gene_col = "gene_name" if "gene_name" in df.columns else None
    if gene_col:
        if target_gene:
            tg_upper = target_gene.strip().upper()
            match = [g for g in df[gene_col].dropna().unique()
                     if str(g).upper() == tg_upper]
            focus_gene = match[0] if match else None
        else:
            focus_gene = None

        if focus_gene is None:
            # Pick gene with largest mean absolute LFC
            gene_effects = df.groupby(gene_col)["raw_score"].apply(
                lambda x: x.abs().mean()
            )
            focus_gene = gene_effects.idxmax()

        sub = df[df[gene_col] == focus_gene].copy()
    else:
        sub = df.copy()
        focus_gene = "region"

    # Build label for cell type
    name_col = "biosample_name" if "biosample_name" in sub.columns else None
    if name_col:
        sub = sub.copy()
        sub["cell_label"] = sub[name_col].fillna("").astype(str)
    else:
        sub["cell_label"] = "all"

    # Pivot: cell_label x variant_label → raw_score (mean if duplicates)
    pivot = sub.pivot_table(
        index="cell_label", columns="variant_label",
        values="raw_score", aggfunc="mean",
    )

    # Keep top 20 most variable cell types
    row_var = pivot.var(axis=1).fillna(0)
    top_cells = row_var.nlargest(min(20, len(row_var))).index
    pivot = pivot.loc[top_cells]

    if pivot.empty:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No data to compare", ha="center", va="center")
        ax.axis("off")
        return fig

    fig, ax = plt.subplots(figsize=(max(6, len(pivot.columns) * 1.5 + 2),
                                     max(4, len(pivot) * 0.35 + 2)))
    vmax = max(abs(pivot.values.min()), abs(pivot.values.max()), 0.001)
    im = ax.imshow(pivot.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")

    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=8)

    # Annotate cells
    for i in range(len(pivot.index)):
        for j in range(len(pivot.columns)):
            val = pivot.values[i, j]
            if not np.isnan(val):
                ax.text(j, i, f"{val:.4f}", ha="center", va="center",
                        fontsize=7, color="white" if abs(val) > vmax * 0.6 else "black")

    fig.colorbar(im, ax=ax, label="LFC (log fold-change)", shrink=0.8)
    ax.set_title(f"Multi-variant LFC comparison — {focus_gene}\n"
                 f"(top {len(pivot)} most variable cell types)", fontsize=11)
    fig.tight_layout()
    return fig


def build_comparison_bar(df, target_gene=None):
    """Build a grouped bar chart: LFC per variant for top cell types.

    Returns:
        matplotlib Figure.
    """
    if df is None or df.empty:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        ax.axis("off")
        return fig

    gene_col = "gene_name" if "gene_name" in df.columns else None
    if gene_col:
        if target_gene:
            tg_upper = target_gene.strip().upper()
            match = [g for g in df[gene_col].dropna().unique()
                     if str(g).upper() == tg_upper]
            focus_gene = match[0] if match else None
        else:
            focus_gene = None

        if focus_gene is None:
            gene_effects = df.groupby(gene_col)["raw_score"].apply(
                lambda x: x.abs().mean()
            )
            focus_gene = gene_effects.idxmax()
        sub = df[df[gene_col] == focus_gene].copy()
    else:
        sub = df.copy()
        focus_gene = "region"

    name_col = "biosample_name" if "biosample_name" in sub.columns else None
    if name_col:
        sub["cell_label"] = sub[name_col].fillna("").astype(str)
    else:
        sub["cell_label"] = "all"

    # Pick top 10 cell types by max absolute LFC across variants
    cell_max = sub.groupby("cell_label")["raw_score"].apply(
        lambda x: x.abs().max()
    ).nlargest(10)
    top_cells = cell_max.index.tolist()
    sub = sub[sub["cell_label"].isin(top_cells)]

    variants = sub["variant_label"].unique()
    n_variants = len(variants)
    n_cells = len(top_cells)

    if n_cells == 0 or n_variants == 0:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        ax.axis("off")
        return fig

    fig, ax = plt.subplots(figsize=(10, max(4, n_cells * 0.5 + 2)))
    bar_height = 0.7 / n_variants
    colors = plt.cm.Set2(np.linspace(0, 1, max(n_variants, 3)))

    for v_idx, variant in enumerate(variants):
        v_data = sub[sub["variant_label"] == variant]
        # Mean LFC per cell type for this variant
        means = v_data.groupby("cell_label")["raw_score"].mean()
        y_positions = [top_cells.index(c) + v_idx * bar_height - 0.35 + bar_height / 2
                       for c in means.index if c in top_cells]
        vals = [means[c] for c in means.index if c in top_cells]
        ax.barh(y_positions, vals, height=bar_height, label=variant,
                color=colors[v_idx], edgecolor="none")

    ax.set_yticks(range(n_cells))
    ax.set_yticklabels(top_cells, fontsize=9)
    ax.axvline(0, color="black", linewidth=0.5)
    ax.set_xlabel("Log fold-change (ALT vs REF)", fontsize=10)
    ax.set_title(f"Variant comparison — {focus_gene}\n"
                 f"(top {n_cells} most affected cell types)", fontsize=11)
    ax.legend(fontsize=8, loc="best")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return fig
