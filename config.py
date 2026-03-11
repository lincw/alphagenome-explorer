"""
AlphaGenome Explorer — Configuration and constants.
"""

from alphagenome.models import dna_client

DEFAULT_API_KEY = ""  # Users provide their own key via the UI

SEQUENCE_LENGTHS = {
    "16 KB": dna_client.SEQUENCE_LENGTH_16KB,
    "100 KB": dna_client.SEQUENCE_LENGTH_100KB,
    "500 KB": dna_client.SEQUENCE_LENGTH_500KB,
    "1 MB": dna_client.SEQUENCE_LENGTH_1MB,
}

OUTPUT_TYPES = {
    "RNA_SEQ": dna_client.OutputType.RNA_SEQ,
    "CAGE": dna_client.OutputType.CAGE,
    "DNASE": dna_client.OutputType.DNASE,
    "ATAC": dna_client.OutputType.ATAC,
    "CHIP_HISTONE": dna_client.OutputType.CHIP_HISTONE,
    "CHIP_TF": dna_client.OutputType.CHIP_TF,
    "SPLICE_SITES": dna_client.OutputType.SPLICE_SITES,
    "SPLICE_SITE_USAGE": dna_client.OutputType.SPLICE_SITE_USAGE,
    "SPLICE_JUNCTIONS": dna_client.OutputType.SPLICE_JUNCTIONS,
    "CONTACT_MAPS": dna_client.OutputType.CONTACT_MAPS,
    "PROCAP": dna_client.OutputType.PROCAP,
}

SCORER_KEYS = [
    "RNA_SEQ", "CAGE", "PROCAP", "ATAC", "DNASE",
    "CHIP_TF", "CHIP_HISTONE",
    "SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS",
    "POLYADENYLATION", "CONTACT_MAPS",
]

# Descriptions for the Welcome tab (key, full name, description, output_metric)
OUTPUT_TYPE_INFO = [
    ("1. RNA_SEQ", "RNA sequencing",
     "Predicts mRNA abundance from total RNA-seq and polyA RNA-seq across 285 biosamples "
     "(667 tracks). Reflects transcription activity across the gene body.",
     "Output: normalized read signal, log(1+x) transformed, at 1 bp resolution. "
     "For variant scoring, reported as relative log-fold-change between alt and ref."),
    ("2. CAGE", "Cap Analysis of Gene Expression",
     "Predicts transcription start site (TSS) usage at single-nucleotide resolution "
     "across 264 biosamples (546 tracks). Captures promoter-level activity.",
     "Output: normalized read signal at TSSs, log(1+x) transformed, at 1 bp resolution."),
    ("3. PROCAP", "Precision Run-On with Cap",
     "Captures active, nascent transcription initiation using PRO-cap sequencing across "
     "6 biosamples (12 tracks). Higher resolution than CAGE for short-lived transcripts.",
     "Output: normalized nascent RNA signal at TSSs, log(1+x) transformed, at 1 bp resolution."),
    ("4. DNASE", "DNase-seq",
     "Predicts chromatin accessibility via DNase I hypersensitivity across 305 biosamples "
     "(305 tracks). Open chromatin indicates active regulatory elements (promoters, enhancers).",
     "Output: normalized insertion signal (untransformed), at 1 bp resolution."),
    ("5. ATAC", "ATAC-seq",
     "Predicts chromatin accessibility using Tn5 transposase insertion across 167 biosamples "
     "(167 tracks). Identifies open regulatory regions.",
     "Output: normalized insertion signal (untransformed), at 1 bp resolution."),
    ("6. CHIP_HISTONE", "ChIP-seq (histone modifications)",
     "Predicts 24 histone modification marks (e.g. H3K27ac, H3K4me1, H3K9me3, H3K27me3, "
     "H3K36me3) across 219 biosamples (1,116 tracks). Marks indicate active promoters, "
     "enhancers, or repressed regions.",
     "Output: fold-change over input/control, log(1+x) transformed, aggregated in 128 bp bins."),
    ("7. CHIP_TF", "ChIP-seq (transcription factors)",
     "Predicts binding profiles of 43 transcription factors (e.g. CTCF, SPI1) across "
     "163 biosamples (1,617 tracks).",
     "Output: fold-change over input/control, log(1+x) transformed, aggregated in 128 bp bins."),
    ("8. SPLICE_SITES", "Splice site strength",
     "Predicts the probability of each position being a splice donor or acceptor site "
     "(4 tracks, both strands). Sequence-intrinsic — no biosample dimension. "
     "Useful for assessing how variants disrupt splicing signals.",
     "Output: probability score (0–1) at 1 bp resolution (untransformed)."),
    ("9. SPLICE_SITE_USAGE", "Splice site usage",
     "Predicts how frequently each splice site is used across 282 biosamples (734 tracks). "
     "Analogous to percent-spliced-in (PSI); captures cell-type-specific splicing.",
     "Output: fraction of transcripts using a given splice site (0–1) at 1 bp resolution."),
    ("10. SPLICE_JUNCTIONS", "Splice junctions",
     "Predicts exon-exon junction read counts (sashimi-style) across 282 biosamples "
     "(734 tracks). Shows which donor-acceptor pairs are connected by splicing.",
     "Output: normalized junction signal, log(1+x) transformed, for all donor-acceptor pairs."),
    ("11. CONTACT_MAPS", "3D contact maps",
     "Predicts chromatin 3D organization (Micro-C / Hi-C) across 12 biosamples (28 tracks). "
     "Reveals enhancer-promoter loops and topologically associating domains (TADs).",
     "Output: log-fold-change over the genomic-distance expectation (power-law decay removed), "
     "at 2,048 bp resolution."),
]

# Maps output type key -> attribute name on prediction result objects
OUTPUT_ATTR_MAP = {
    "RNA_SEQ": "rna_seq",
    "CAGE": "cage",
    "DNASE": "dnase",
    "ATAC": "atac",
    "CHIP_HISTONE": "chip_histone",
    "CHIP_TF": "chip_tf",
    "SPLICE_SITES": "splice_sites",
    "SPLICE_SITE_USAGE": "splice_site_usage",
    "SPLICE_JUNCTIONS": "splice_junctions",
    "CONTACT_MAPS": "contact_maps",
    "PROCAP": "procap",
}
