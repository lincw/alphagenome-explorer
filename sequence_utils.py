"""
AlphaGenome Explorer — Sequence prediction utilities.

Helpers for raw DNA sequence input: validation, padding, and length selection.
"""

import re

from alphagenome.models import dna_client

# Valid sequence lengths in ascending order
VALID_LENGTHS = [
    ("16 KB", dna_client.SEQUENCE_LENGTH_16KB),
    ("100 KB", dna_client.SEQUENCE_LENGTH_100KB),
    ("500 KB", dna_client.SEQUENCE_LENGTH_500KB),
    ("1 MB", dna_client.SEQUENCE_LENGTH_1MB),
]

# Regex for valid DNA characters (case-insensitive)
_DNA_PATTERN = re.compile(r"^[ACGTNacgtn]+$")


def validate_sequence(seq):
    """Validate a raw DNA sequence string.

    Returns:
        (cleaned_seq, error_msg) — cleaned_seq is uppercase, error_msg is None if OK.
    """
    seq = seq.strip()
    if not seq:
        return None, "Enter a DNA sequence."

    # Remove whitespace and newlines within the sequence
    seq = re.sub(r"\s+", "", seq)
    seq = seq.upper()

    if not _DNA_PATTERN.match(seq):
        bad_chars = set(seq) - set("ACGTN")
        return None, f"Invalid characters: {', '.join(sorted(bad_chars))}. Only A, C, G, T, N allowed."

    if len(seq) > dna_client.SEQUENCE_LENGTH_1MB:
        return None, f"Sequence too long ({len(seq):,} bp). Maximum is {dna_client.SEQUENCE_LENGTH_1MB:,} bp."

    return seq, None


def pad_sequence(seq, target_length=None):
    """Pad a DNA sequence with N to a valid model length.

    If target_length is None, uses the smallest valid length >= len(seq).

    Returns:
        (padded_seq, actual_length, error_msg)
    """
    if target_length is None:
        # Auto-select smallest valid length that fits
        for _, length in VALID_LENGTHS:
            if len(seq) <= length:
                target_length = length
                break
        if target_length is None:
            return None, 0, "Sequence exceeds maximum length."

    if len(seq) > target_length:
        return None, 0, (
            f"Sequence ({len(seq):,} bp) exceeds target length ({target_length:,} bp). "
            f"Choose a larger length or trim the sequence."
        )

    # Center the sequence with N padding
    padded = seq.center(target_length, "N")
    return padded, target_length, None


def get_sequence_info(seq):
    """Get basic statistics about a DNA sequence.

    Returns:
        dict with keys: length, gc_content, n_count, non_n_length.
    """
    length = len(seq)
    n_count = seq.count("N")
    non_n = length - n_count
    gc = seq.count("G") + seq.count("C")
    gc_pct = (gc / non_n * 100) if non_n > 0 else 0.0

    return {
        "length": length,
        "gc_content": gc_pct,
        "n_count": n_count,
        "non_n_length": non_n,
    }
