"""
AlphaGenome Explorer — Variant lookup utilities.

Functions for fetching reference sequence from UCSC and looking up rsIDs
from the NCBI Variation API (dbSNP).
"""

import json
import logging
import urllib.parse
import urllib.request

logger = logging.getLogger(__name__)

# Map NCBI RefSeq accessions (GRCh38) to UCSC chromosome names
_ACCESSION_TO_CHROM = {
    f"NC_0000{i:02d}.{v}": f"chr{i}"
    for i, v in [
        (1, 11), (2, 12), (3, 12), (4, 12), (5, 10), (6, 12),
        (7, 14), (8, 11), (9, 12), (10, 11), (11, 10), (12, 12),
        (13, 11), (14, 9), (15, 10), (16, 10), (17, 11), (18, 10),
        (19, 10), (20, 11), (21, 9), (22, 12), (23, 11), (24, 10),
    ]
}
# Also support unversioned lookup
for _acc, _chrom in list(_ACCESSION_TO_CHROM.items()):
    _ACCESSION_TO_CHROM[_acc.split(".")[0]] = _chrom


def fetch_ref_sequence(chrom, start_0based, end_0based):
    """Fetch reference sequence from UCSC REST API.

    Args:
        chrom: UCSC-style chromosome name (e.g. "chr3").
        start_0based: 0-based start position (inclusive).
        end_0based: 0-based end position (exclusive).

    Returns:
        Uppercase DNA string, or "" on failure.
    """
    try:
        params = urllib.parse.urlencode({
            "genome": "hg38",
            "chrom": str(chrom),
            "start": int(start_0based),
            "end": int(end_0based),
        })
        url = f"https://api.genome.ucsc.edu/getData/sequence?{params}"
        with urllib.request.urlopen(url, timeout=5) as resp:
            data = json.loads(resp.read())
        return data["dna"].upper()
    except Exception as e:
        logger.warning("Failed to fetch ref sequence %s:%s-%s: %s",
                       chrom, start_0based, end_0based, e)
        return ""


def lookup_rsid(rsid_str):
    """Look up an rsID from the NCBI Variation API (dbSNP).

    Parses the GRCh38 placement and converts SPDI coordinates to VCF-style
    (1-based position with padding base for indels).

    Args:
        rsid_str: e.g. "rs35482426" or "35482426".

    Returns:
        dict with keys {chrom, pos, ref, alts} where alts is a list of
        {chrom, pos, ref, alt, variant_type} dicts.  Returns None on failure.
    """
    rsid_str = rsid_str.strip().lower()
    if rsid_str.startswith("rs"):
        rsid_str = rsid_str[2:]
    try:
        int(rsid_str)
    except ValueError:
        raise ValueError(f"'{rsid_str}' is not a valid rsID number.")

    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{rsid_str}"
    try:
        with urllib.request.urlopen(url, timeout=20) as resp:
            data = json.loads(resp.read())
    except urllib.error.HTTPError as e:
        raise ConnectionError(f"NCBI API returned HTTP {e.code} for rs{rsid_str}.")
    except urllib.error.URLError as e:
        raise ConnectionError(f"Cannot reach NCBI API: {e.reason}")
    except Exception as e:
        raise ConnectionError(f"NCBI API request failed: {e}")

    # Some rsIDs have been merged into another
    if "merged_snapshot_data" in data and "primary_snapshot_data" not in data:
        merged_into = data["merged_snapshot_data"].get("merged_into", [])
        if merged_into:
            raise ValueError(
                f"rs{rsid_str} has been merged into rs{merged_into[0]}. "
                f"Try rs{merged_into[0]} instead."
            )
        raise ValueError(f"rs{rsid_str} has been merged and has no current data.")

    # Find GRCh38 top-level placement (NC_ accession)
    placements = data.get("primary_snapshot_data", {}).get(
        "placements_with_allele", []
    )
    grch38_placement = None
    for p in placements:
        seq_id = p.get("seq_id", "")
        if not seq_id.startswith("NC_"):
            continue
        annot = p.get("placement_annot", {})
        for asm in annot.get("seq_id_traits_by_assembly", []):
            if asm.get("assembly_name", "").startswith("GRCh38"):
                grch38_placement = p
                break
        if grch38_placement:
            break

    if grch38_placement is None:
        raise ValueError(f"No GRCh38 placement found for rs{rsid_str}.")

    seq_id = grch38_placement["seq_id"]
    chrom = _ACCESSION_TO_CHROM.get(seq_id) or _ACCESSION_TO_CHROM.get(
        seq_id.split(".")[0]
    )
    if not chrom:
        # Fallback: derive from accession number
        try:
            num = int(seq_id.split("_")[1].split(".")[0])
            if num <= 22:
                chrom = f"chr{num}"
            elif num == 23:
                chrom = "chrX"
            elif num == 24:
                chrom = "chrY"
        except Exception:
            raise ValueError(f"Cannot determine chromosome from accession {seq_id}.")

    # Parse SPDI alleles
    alleles = grch38_placement.get("alleles", [])
    ref_pos = None
    alt_variants = []

    for a in alleles:
        spdi = a.get("allele", {}).get("spdi", {})
        if not spdi:
            continue
        pos_0 = spdi.get("position", 0)
        deleted = spdi.get("deleted_sequence", "")
        inserted = spdi.get("inserted_sequence", "")

        if deleted == inserted:
            ref_pos = pos_0
        else:
            alt_variants.append({
                "pos_0": pos_0,
                "deleted": deleted,
                "inserted": inserted,
            })

    if ref_pos is None or not alt_variants:
        raise ValueError(f"Could not parse alleles for rs{rsid_str}.")

    # Convert SPDI to VCF-style coordinates
    results = []
    for alt in alt_variants:
        deleted = alt["deleted"]
        inserted = alt["inserted"]

        if len(deleted) > 0 and len(inserted) > 0 and len(deleted) == len(inserted):
            # Simple substitution (SNP or MNP) — no padding needed
            vcf_pos = ref_pos + 1  # 1-based
            vcf_ref = deleted
            vcf_alt = inserted
            vtype = "SNP" if len(deleted) == 1 else "MNP"
        else:
            # Indel — need padding base before the variant
            pad_base = fetch_ref_sequence(chrom, ref_pos - 1, ref_pos)
            if not pad_base:
                pad_base = "N"
            vcf_pos = ref_pos  # 1-based (the padding base position)
            vcf_ref = pad_base + deleted
            vcf_alt = pad_base + inserted
            if not inserted:
                vtype = "DEL"
            elif not deleted:
                vtype = "INS"
            else:
                vtype = "INDEL"

        results.append({
            "chrom": chrom,
            "pos": vcf_pos,
            "ref": vcf_ref,
            "alt": vcf_alt,
            "variant_type": vtype,
        })

    if not results:
        return None

    return {
        "chrom": results[0]["chrom"],
        "pos": results[0]["pos"],
        "ref": results[0]["ref"],
        "alts": results,
    }
