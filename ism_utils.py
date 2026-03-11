"""
AlphaGenome Explorer — In Silico Mutagenesis (ISM) utilities.

Helpers for running ISM analysis and building SeqLogo visualizations.
"""

import logging

import numpy as np

from alphagenome.data import genome
from alphagenome.interpretation import ism
from alphagenome.models import dna_client, variant_scorers
from alphagenome.visualization import plot_components

logger = logging.getLogger(__name__)

# ISM widths available in the UI
ISM_WIDTHS = {
    "128": 128,
    "256": 256,
    "512": 512,
}

# Output types supported for ISM scoring (CenterMaskScorer-compatible)
ISM_OUTPUT_TYPES = {
    "DNASE": dna_client.OutputType.DNASE,
    "ATAC": dna_client.OutputType.ATAC,
    "CAGE": dna_client.OutputType.CAGE,
    "CHIP_HISTONE": dna_client.OutputType.CHIP_HISTONE,
    "CHIP_TF": dna_client.OutputType.CHIP_TF,
    "RNA_SEQ": dna_client.OutputType.RNA_SEQ,
    "PROCAP": dna_client.OutputType.PROCAP,
}


def run_ism(dna_model, center_position, chromosome, ism_width, output_type_key,
            organism, progress_callback=None):
    """Run ISM analysis and return (variant_scores, ism_interval, sequence_interval).

    Args:
        dna_model: The AlphaGenome model instance.
        center_position: 1-based genomic position to center ISM on.
        chromosome: e.g. "chr12".
        ism_width: Width of ISM interval in bp (128, 256, or 512).
        output_type_key: Key from ISM_OUTPUT_TYPES (e.g. "DNASE").
        organism: dna_client.Organism enum value.
        progress_callback: Optional callable(fraction, message).

    Returns:
        (variant_scores, ism_interval, sequence_interval)
    """
    # Build intervals — use 16KB context for speed
    sequence_interval = genome.Interval(
        chromosome=chromosome,
        start=center_position,
        end=center_position,
    ).resize(dna_client.SEQUENCE_LENGTH_16KB)

    ism_interval = sequence_interval.resize(ism_width)

    # Build scorer
    output_type = ISM_OUTPUT_TYPES[output_type_key]
    scorer = variant_scorers.CenterMaskScorer(
        requested_output=output_type,
        width=501,
        aggregation_type=variant_scorers.AggregationType.DIFF_MEAN,
    )

    if progress_callback:
        progress_callback(0.1, "Running ISM (this may take a minute)...")

    variant_scores = dna_model.score_ism_variants(
        interval=sequence_interval,
        ism_interval=ism_interval,
        variant_scorers=[scorer],
    )

    if progress_callback:
        progress_callback(0.9, "Processing results...")

    return variant_scores, ism_interval, sequence_interval


def extract_ism_matrix(variant_scores, ontology_curie):
    """Extract ISM matrix for a specific cell type / ontology term.

    Args:
        variant_scores: List of (scorer_results,) from score_ism_variants.
        ontology_curie: Ontology term to extract (e.g. "EFO:0002067").

    Returns:
        numpy array of shape (ism_width, 4) — the ISM contribution scores.
    """
    def extract_ontology(adata):
        mask = adata.var["ontology_curie"] == ontology_curie
        if mask.sum() == 0:
            # Fall back to averaging all tracks
            logger.warning(
                "Ontology %s not found in ISM results, averaging all tracks.",
                ontology_curie,
            )
            return float(np.mean(adata.X))
        values = adata.X[:, mask]
        return float(np.mean(values))

    ism_result = ism.ism_matrix(
        [extract_ontology(x[0]) for x in variant_scores],
        variants=[v[0].uns["variant"] for v in variant_scores],
    )
    return ism_result


def build_ism_plot(ism_result, ism_interval, output_type_key, ontology_label):
    """Build a SeqLogo figure from ISM results.

    Args:
        ism_result: numpy array (ism_width, 4) from extract_ism_matrix.
        ism_interval: genome.Interval for the ISM region.
        output_type_key: e.g. "DNASE".
        ontology_label: Human-readable label for the cell type.

    Returns:
        matplotlib Figure.
    """
    fig = plot_components.plot(
        [
            plot_components.SeqLogo(
                scores=ism_result,
                scores_interval=ism_interval,
                ylabel=f"ISM {output_type_key}\n{ontology_label}",
            )
        ],
        interval=ism_interval,
        fig_width=35,
        title=f"In Silico Mutagenesis — {output_type_key} ({ontology_label})",
    )
    return fig


def get_available_ontologies(variant_scores):
    """Get list of ontology terms available in the ISM results.

    Returns:
        List of (ontology_curie, name) tuples.
    """
    if not variant_scores or len(variant_scores) == 0:
        return []
    adata = variant_scores[0][0]
    meta = adata.var
    if "ontology_curie" not in meta.columns:
        return []
    pairs = meta[["ontology_curie", "name"]].drop_duplicates()
    return [(row.ontology_curie, row["name"]) for row in pairs.itertuples()]
