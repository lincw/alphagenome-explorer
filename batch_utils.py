"""
AlphaGenome Explorer — Batch variant scoring utilities.

Parse VCF-style input, score variants with multiple scorers, and return results.
"""

import logging
from io import StringIO

import pandas as pd

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

logger = logging.getLogger(__name__)

MAX_BATCH_VARIANTS = 20

# All recommended scorers, keyed by lowercase name
ALL_SCORER_KEYS = list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.keys())


def parse_vcf_text(text):
    """Parse VCF-style tab-separated text into a list of variant dicts.

    Expected columns: variant_id, CHROM, POS, REF, ALT
    (header required; tab or comma separated)

    Returns:
        (variants: list[dict], errors: list[str])
        Each variant dict has keys: variant (genome.Variant), label (str).
    """
    text = text.strip()
    if not text:
        return [], ["Empty input."]

    # Auto-detect separator
    first_line = text.splitlines()[0]
    sep = "\t" if "\t" in first_line else ","

    try:
        df = pd.read_csv(StringIO(text), sep=sep)
    except Exception as e:
        return [], [f"Could not parse input: {e}"]

    # Normalize column names
    df.columns = df.columns.str.strip()

    required = ["CHROM", "POS", "REF", "ALT"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        return [], [f"Missing required columns: {', '.join(missing)}. "
                    f"Found: {', '.join(df.columns)}"]

    if len(df) > MAX_BATCH_VARIANTS:
        return [], [f"Too many variants ({len(df)}). Maximum is {MAX_BATCH_VARIANTS}."]

    variants = []
    errors = []
    for i, row in df.iterrows():
        try:
            chrom = str(row["CHROM"]).strip()
            pos = int(row["POS"])
            ref = str(row["REF"]).strip()
            alt = str(row["ALT"]).strip()
            name = str(row.get("variant_id", f"{chrom}:{pos}:{ref}>{alt}")).strip()

            variant = genome.Variant(
                chromosome=chrom,
                position=pos,
                reference_bases=ref,
                alternate_bases=alt,
                name=name,
            )
            variants.append({"variant": variant, "label": name})
        except Exception as e:
            errors.append(f"Row {i + 1}: {e}")

    return variants, errors


def get_selected_scorers(scorer_keys, organism):
    """Build list of scorer objects from selected keys, filtering unsupported ones.

    Args:
        scorer_keys: List of scorer key strings (e.g. ["RNA_SEQ", "DNASE"]).
        organism: dna_client.Organism enum value.

    Returns:
        (selected_scorers: list, excluded: list[str])
    """
    all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    selected = []
    excluded = []

    for key in scorer_keys:
        scorer = all_scorers.get(key)
        if scorer is None:
            excluded.append(f"{key} (not found)")
            continue

        # Check organism support
        try:
            supported = variant_scorers.SUPPORTED_ORGANISMS.get(
                scorer.base_variant_scorer, {}
            )
            if organism.value not in supported:
                excluded.append(f"{key} (not supported for {organism.name})")
                continue
        except Exception:
            pass

        # Check PROCAP for mouse
        try:
            if (scorer.requested_output == dna_client.OutputType.PROCAP
                    and organism == dna_client.Organism.MUS_MUSCULUS):
                excluded.append(f"{key} (PROCAP unavailable for mouse)")
                continue
        except Exception:
            pass

        selected.append(scorer)

    return selected, excluded


def run_batch_scoring(model, variants, scorers, sequence_length, organism,
                      progress_callback=None):
    """Score a batch of variants with multiple scorers.

    Each variant uses its own interval (centered on the variant position).

    Args:
        model: AlphaGenome model instance.
        variants: List of dicts with "variant" (genome.Variant) and "label".
        scorers: List of variant scorer objects.
        sequence_length: Sequence length constant (e.g. SEQUENCE_LENGTH_1MB).
        organism: dna_client.Organism.
        progress_callback: Optional callable(fraction, message).

    Returns:
        pd.DataFrame from tidy_scores with all variants and scorers.
    """
    all_results = []
    n = len(variants)
    errors = []

    for i, v in enumerate(variants):
        variant = v["variant"]
        label = v["label"]

        if progress_callback:
            progress_callback(
                (i + 0.5) / (n + 1),
                f"Scoring variant {i + 1}/{n}: {label}...",
            )

        interval = variant.reference_interval.resize(sequence_length)

        try:
            scores = model.score_variant(
                interval=interval,
                variant=variant,
                variant_scorers=scorers,
                organism=organism,
            )
            all_results.append(scores)
        except Exception as e:
            logger.error("Batch scoring failed for %s: %s", label, e)
            errors.append(f"{label}: {e}")

    if not all_results:
        raise RuntimeError(
            "All variants failed to score. "
            + (" ".join(errors[:3]) if errors else "")
        )

    if progress_callback:
        progress_callback(n / (n + 1), "Building results table...")

    df = variant_scorers.tidy_scores(all_results, match_gene_strand=True)
    return df, errors
