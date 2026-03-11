"""
AlphaGenome Explorer — UI layout definition.
"""

from shiny import ui

from batch_utils import ALL_SCORER_KEYS
from config import OUTPUT_TYPE_INFO, OUTPUT_TYPES, SCORER_KEYS, SEQUENCE_LENGTHS
from ism_utils import ISM_OUTPUT_TYPES, ISM_WIDTHS


def _download_bar(*buttons):
    """Wrap download/action buttons in a consistent toolbar."""
    return ui.div(*buttons, class_="mb-2 d-flex gap-2 flex-wrap")


def create_sidebar():
    """Build the sidebar with all input controls."""
    return ui.sidebar(
        # --- API Key ---
        ui.input_password("api_key", "Google AI Studio API Key", placeholder="AIza..."),
        ui.p(
            ui.tags.a(
                "Get API key",
                href="https://deepmind.google.com/science/alphagenome/",
                target="_blank",
            ),
            class_="small text-muted mb-0",
        ),
        ui.output_ui("api_key_status_ui"),

        # --- 1. Gene / Interval ---
        ui.h5("1. Gene or Region"),
        ui.input_text("gene_symbol", "Gene symbol", value="OAS1"),
        ui.input_select(
            "organism", "Organism",
            {"human": "Human (hg38)", "mouse": "Mouse (mm10)"},
        ),
        ui.input_select(
            "seq_length", "Sequence length",
            list(SEQUENCE_LENGTHS.keys()),
            selected="1 MB",
        ),
        ui.input_action_button("btn_lookup", "Look up gene", class_="btn-primary btn-sm w-100"),
        ui.output_text("gene_info"),

        # --- 2. Output types ---
        ui.h5("2. Output types"),
        ui.input_checkbox_group(
            "output_types", None,
            choices=list(OUTPUT_TYPES.keys()),
            selected=["RNA_SEQ"],
        ),

        # --- 3. Ontology ---
        ui.h5("3. Cell type / Tissue"),
        ui.input_text(
            "ontology_search", "Search keyword", placeholder="e.g. monocyte, lung, K562"
        ),
        ui.input_action_button("btn_search_onto", "Search", class_="btn-sm btn-outline-primary w-100"),
        ui.output_ui("ontology_select_ui"),

        # --- 4. Predict ---
        ui.h5("4. Predict tracks"),
        ui.p(
            "Predict signal across the gene region. "
            "Requires cell type from Step 3.",
            class_="text-muted small",
        ),
        ui.input_action_button("btn_predict", "Run prediction", class_="btn-success btn-sm w-100"),

        # --- 5. Variant ---
        ui.h5("5. Variant analysis (optional)"),
        ui.p(
            "Enter an rsID to auto-fill from dbSNP, or enter manually. "
            "Requires cell type from Step 3.",
            class_="text-muted small",
        ),
        ui.output_ui("var_gene_info_ui"),
        ui.div(
            ui.input_text("var_rsid", "rsID", placeholder="e.g. rs35482426"),
            ui.input_action_button(
                "btn_lookup_rsid", "Look up rsID",
                class_="btn-sm btn-outline-secondary w-100 mb-2",
            ),
            ui.output_ui("rsid_status_ui"),
        ),
        ui.hr(class_="my-2"),
        ui.input_text("var_pos", "Chromosomal position", placeholder="e.g. 112919388"),
        ui.input_text("var_ref", "REF (auto-filled, editable)", placeholder="auto-detected from hg38"),
        ui.input_text("var_alt", "ALT", placeholder="A"),
        ui.input_select(
            "var_display_window", "Display window (zoom)",
            choices={
                "5000": "5 KB",
                "10000": "10 KB",
                "20000": "20 KB",
                "50000": "50 KB",
                "100000": "100 KB",
                "0": "Full gene region",
            },
            selected="10000",
        ),
        ui.input_action_button("btn_variant", "Predict variant effect", class_="btn-warning btn-sm w-100"),

        # --- 6. Variant scoring ---
        ui.h5("6. Variant scoring (optional)"),
        ui.p(
            "Scores the variant across ALL cell types automatically "
            "(does not use the cell type selected in Step 3).",
            class_="text-muted small",
        ),
        ui.input_select(
            "scorer_key", "Scorer",
            choices=SCORER_KEYS,
            selected="RNA_SEQ",
        ),
        ui.input_action_button("btn_score", "Score variant", class_="btn-info btn-sm w-100"),

        # --- 7. ISM ---
        ui.h5("7. In silico mutagenesis (ISM)"),
        ui.p(
            "Systematically mutate every base in a short region to find "
            "functionally important positions. Uses 16 KB context for speed. "
            "Requires cell type from Step 3.",
            class_="text-muted small",
        ),
        ui.input_text(
            "ism_center", "Center position",
            placeholder="e.g. 112919388 (or leave blank to use variant pos)",
        ),
        ui.input_select(
            "ism_width", "ISM width (bp)",
            choices={k: f"{v} bp" for k, v in ISM_WIDTHS.items()},
            selected="256",
        ),
        ui.input_select(
            "ism_output_type", "Output type to score",
            choices=list(ISM_OUTPUT_TYPES.keys()),
            selected="DNASE",
        ),
        ui.input_action_button("btn_ism", "Run ISM", class_="btn-secondary btn-sm w-100"),

        # --- 8. Multi-variant comparison ---
        ui.h5("8. Multi-variant comparison"),
        ui.p(
            "Enter multiple variants to compare their scoring results. "
            "One variant per line: chr:pos:ref:alt or rsID. "
            "Optionally add a label after a space. Max 10 variants.",
            class_="text-muted small",
        ),
        ui.input_text_area(
            "multi_variants", "Variants (one per line)",
            placeholder=(
                "chr12:112919388:G:A disease\n"
                "chr12:112919400:C:T benign\n"
                "rs35482426 risk"
            ),
            rows=5,
        ),
        ui.input_select(
            "multi_scorer_key", "Scorer",
            choices=SCORER_KEYS,
            selected="RNA_SEQ",
        ),
        ui.input_action_button(
            "btn_multi_score", "Compare variants",
            class_="btn-danger btn-sm w-100",
        ),

        # --- 9. Batch variant scoring ---
        ui.h5("9. Batch variant scoring"),
        ui.p(
            "Score many variants with multiple scorers at once. "
            "Paste VCF-style data (tab or comma separated). "
            "Each variant uses its own interval. Max 20 variants.",
            class_="text-muted small",
        ),
        ui.input_text_area(
            "batch_vcf_text", "VCF data (with header)",
            placeholder=(
                "variant_id\tCHROM\tPOS\tREF\tALT\n"
                "snp1\tchr12\t112919388\tG\tA\n"
                "snp2\tchr3\t58394738\tA\tT"
            ),
            rows=6,
        ),
        ui.input_checkbox_group(
            "batch_scorers", "Scorers",
            choices=ALL_SCORER_KEYS,
            selected=["RNA_SEQ", "CAGE", "DNASE"],
        ),
        ui.input_select(
            "batch_seq_length", "Sequence length",
            choices=list(SEQUENCE_LENGTHS.keys()),
            selected="1 MB",
        ),
        ui.input_action_button(
            "btn_batch_score", "Run batch scoring",
            class_="btn-dark btn-sm w-100",
        ),

        # --- 10. Predict from sequence ---
        ui.h5("10. Predict from sequence"),
        ui.p(
            "Predict signals from a raw DNA sequence (no genomic coordinates). "
            "Useful for synthetic or edited sequences. "
            "Requires output types (Step 2) and cell type (Step 3).",
            class_="text-muted small",
        ),
        ui.input_text_area(
            "seq_input", "DNA sequence (A/C/G/T/N)",
            placeholder="Paste DNA sequence here (up to 1 MB)...",
            rows=5,
        ),
        ui.input_select(
            "seq_pad_length", "Pad to length",
            choices={
                "auto": "Auto (smallest valid length)",
                **{k: k for k in SEQUENCE_LENGTHS.keys()},
            },
            selected="auto",
        ),
        ui.input_action_button(
            "btn_predict_seq", "Predict from sequence",
            class_="btn-outline-success btn-sm w-100",
        ),

        width=340,
    )


def _welcome_tab():
    """Build the Welcome / Documentation tab."""
    rows = []
    for name, full_name, description, output_metric in OUTPUT_TYPE_INFO:
        rows.append(
            ui.div(
                ui.div(
                    ui.strong(name), ui.span(f" ({full_name})", class_="text-muted"),
                    class_="mb-1",
                ),
                ui.p(description, class_="small mb-1"),
                ui.p(
                    ui.span("⬡ ", class_="text-primary"),
                    output_metric,
                    class_="small mb-0 text-secondary fst-italic",
                ),
                class_="border-bottom py-2",
            )
        )

    return ui.nav_panel(
        "Welcome",
        ui.div(
            ui.h4("Welcome to AlphaGenome Explorer"),
            ui.p(
                "AlphaGenome is a DNA foundation model by Google DeepMind that predicts "
                "how the genome encodes regulatory information. Given a DNA sequence "
                "(up to 1 MB, human or mouse), it predicts epigenomic and transcriptomic "
                "signals across hundreds of cell types and tissues.",
                class_="lead",
            ),
            ui.h5("What can AlphaGenome predict?"),
            ui.p(
                "The model produces 11 output modalities covering gene expression, "
                "chromatin accessibility, histone modifications, transcription factor binding, "
                "splicing, and 3D genome organization:",
            ),
            ui.div(*rows, class_="mb-4"),
            ui.h5("Quick-start guide"),
            ui.tags.ol(
                ui.tags.li(ui.strong("API key: "),
                           "Paste your Google AI Studio API key in the sidebar. ",
                           ui.tags.a("Get a free key here.",
                                     href="https://deepmind.google.com/science/alphagenome/",
                                     target="_blank")),
                ui.tags.li(ui.strong("Step 1 — Gene lookup: "),
                           "Enter a gene symbol (e.g. OAS1) and click 'Look up gene'."),
                ui.tags.li(ui.strong("Step 2 — Output types: "),
                           "Select which modalities to predict (RNA_SEQ is selected by default)."),
                ui.tags.li(ui.strong("Step 3 — Cell type: "),
                           "Search for a cell type or tissue (e.g. 'monocyte', 'K562') and select from results."),
                ui.tags.li(ui.strong("Step 4 — Predict tracks: "),
                           "Click 'Run prediction' to generate signal tracks across the gene region."),
                ui.tags.li(ui.strong("Step 5 — Variant analysis (optional): "),
                           "Enter a variant (position, REF, ALT) to compare predicted signals between alleles."),
                ui.tags.li(ui.strong("Step 6 — Variant scoring (optional): "),
                           "Quantify gene-level expression change (log-fold-change) and compare across cell types."),
                ui.tags.li(ui.strong("Step 7 — ISM (optional): "),
                           "Run in silico mutagenesis to find functionally important bases near a position of interest."),
                ui.tags.li(ui.strong("Step 8 — Multi-variant comparison (optional): "),
                           "Enter multiple variants to compare their effects side-by-side."),
                ui.tags.li(ui.strong("Step 9 — Batch scoring (optional): "),
                           "Paste VCF data to score many variants with multiple scorers. "
                           "Produces a comprehensive results table for downstream analysis."),
                ui.tags.li(ui.strong("Step 10 — Predict from sequence (optional): "),
                           "Paste a raw DNA sequence (A/C/G/T/N) to predict signals "
                           "without genomic coordinates. Useful for synthetic or edited sequences."),
                class_="mb-4",
            ),
            ui.p(
                "Use the sidebar on the left to follow these steps. Results will appear "
                "in the corresponding tabs above. Each results tab has download buttons "
                "for exporting data (CSV) and plots (PNG).",
                class_="text-muted",
            ),
            class_="p-4",
        ),
    )


def create_main_panel():
    """Build the main tabbed panel."""
    return ui.navset_tab(
        _welcome_tab(),
        ui.nav_panel(
            "Gene context",
            ui.div(
                ui.h6("How to read this plot", class_="mb-2"),
                ui.p(
                    "This shows the transcript structure (exons and introns) for the looked-up gene "
                    "and its neighbors within the prediction window. ",
                    ui.strong("Thick blocks"), " are exons (protein-coding regions), ",
                    ui.strong("thin lines"), " are introns, and ",
                    ui.strong("arrows"), " indicate the direction of transcription. "
                    "Multiple rows may appear if MANE Select transcripts exist for multiple "
                    "genes in the window.",
                    class_="small mb-0",
                ),
                class_="bg-info bg-opacity-10 p-3 rounded mb-3",
            ),
            ui.output_plot("plot_gene_context", height="500px"),
        ),
        ui.nav_panel(
            "Predicted tracks",
            ui.div(
                ui.h6("How to read this plot", class_="mb-2"),
                ui.p(
                    "Each row is a predicted signal track for the selected cell type and output modality. "
                    "The x-axis is the genomic position; the y-axis is predicted signal intensity. ",
                    class_="small mb-1",
                ),
                ui.tags.ul(
                    ui.tags.li("RNA-seq / CAGE / PROCAP: peaks indicate active transcription or TSS usage.",
                               class_="small"),
                    ui.tags.li("DNase / ATAC: peaks mark open chromatin (regulatory elements like promoters and enhancers).",
                               class_="small"),
                    ui.tags.li("ChIP (histone/TF): peaks show histone modifications or transcription factor binding sites.",
                               class_="small"),
                    ui.tags.li("Splice junctions (sashimi): arcs connect spliced exon pairs; height = junction read count.",
                               class_="small"),
                    ui.tags.li("Contact maps: heatmap of predicted 3D chromatin interactions (enhancer-promoter loops).",
                               class_="small"),
                    class_="mb-0",
                ),
                class_="bg-info bg-opacity-10 p-3 rounded mb-3",
            ),
            _download_bar(
                ui.download_button("dl_tracks_csv", "CSV", class_="btn-outline-secondary btn-sm"),
                ui.download_button("dl_tracks_png", "PNG", class_="btn-outline-secondary btn-sm"),
            ),
            ui.output_ui("plot_interval_ui"),
        ),
        ui.nav_panel(
            "Variant effect",
            ui.div(
                ui.h6("How to read this plot", class_="mb-2"),
                ui.p(
                    "This compares the predicted signal between REF and ALT alleles for your variant. "
                    "Each output type shows two rows:",
                    class_="small mb-1",
                ),
                ui.tags.ul(
                    ui.tags.li(
                        ui.strong("Overlaid tracks: "),
                        "Gray = REF prediction, Red = ALT prediction. Where the red line "
                        "diverges from gray, the variant changes the predicted signal.",
                        class_="small",
                    ),
                    ui.tags.li(
                        ui.strong("Difference track (ALT - REF): "),
                        "Positive peaks (above zero) mean ALT has stronger signal; "
                        "negative peaks (below zero) mean REF is stronger. "
                        "Peaks near the variant position indicate local effects; "
                        "distant peaks suggest long-range regulatory changes.",
                        class_="small",
                    ),
                    class_="mb-1",
                ),
                ui.p(
                    "The vertical dashed line marks the variant position. "
                    "Use the 'Display window' dropdown to zoom in/out around the variant.",
                    class_="small text-muted mb-0",
                ),
                class_="bg-info bg-opacity-10 p-3 rounded mb-3",
            ),
            _download_bar(
                ui.download_button("dl_variant_png", "PNG", class_="btn-outline-secondary btn-sm"),
            ),
            ui.output_ui("plot_variant_ui"),
        ),
        ui.nav_panel(
            "Variant scores",
            ui.div(
                ui.h6("How to read these results", class_="mb-2"),
                ui.p(
                    "Variant scoring quantifies the overall effect of the variant on each gene, "
                    "summed over the entire gene body (not per-base). Results are shown in three parts:",
                    class_="small mb-1",
                ),
                ui.tags.ul(
                    ui.tags.li(
                        ui.strong("Gene-level summary: "),
                        "Grouped by gene, showing how the variant changes expression in each cell type. "
                        "Genes are sorted by the strongest effect first. "
                        "Green = ALT is higher, Red = REF is higher.",
                        class_="small",
                    ),
                    ui.tags.li(
                        ui.strong("Bar chart: "),
                        "Compares the variant's effect on the most-affected gene across all cell types. "
                        "Shows the top 10 most increased + top 10 most decreased cell types. "
                        "Longer bars = larger effect.",
                        class_="small",
                    ),
                    ui.tags.li(
                        ui.strong("LFC (log-fold-change): "),
                        "log2(ALT / REF) of the summed gene-body signal. "
                        "LFC = 0.01 means ~1% increase; LFC = 0.1 means ~7% increase.",
                        class_="small",
                    ),
                    ui.tags.li(
                        ui.strong("Quantile: "),
                        "How extreme this change is compared to ~348,000 common SNPs. "
                        "Quantile = 1.000 means this variant's effect is larger than almost all "
                        "common variants; quantile near 0 means it's a typical (small) effect.",
                        class_="small",
                    ),
                    class_="mb-0",
                ),
                class_="bg-info bg-opacity-10 p-3 rounded mb-3",
            ),
            _download_bar(
                ui.download_button("dl_scores_csv", "CSV", class_="btn-outline-secondary btn-sm"),
                ui.download_button("dl_scores_png", "PNG", class_="btn-outline-secondary btn-sm"),
            ),
            ui.output_ui("score_summary_ui"),
            ui.output_ui("plot_scores_ui"),
            ui.output_data_frame("table_scores"),
        ),
        ui.nav_panel(
            "ISM",
            ui.div(
                ui.h6("How to read this plot", class_="mb-2"),
                ui.p(
                    "In silico mutagenesis (ISM) systematically tests the effect of every possible "
                    "single-nucleotide change in a short region. The result is displayed as a ",
                    ui.strong("sequence logo"),
                    ":",
                    class_="small mb-1",
                ),
                ui.tags.ul(
                    ui.tags.li(
                        "Each position shows up to 4 letters (A, C, G, T). "
                        "Letter height = magnitude of that mutation's effect on the scored signal.",
                        class_="small",
                    ),
                    ui.tags.li(
                        ui.strong("Tall letters "),
                        "indicate positions where mutations cause large changes — "
                        "these are functionally important bases.",
                        class_="small",
                    ),
                    ui.tags.li(
                        ui.strong("Positive (upward) letters "),
                        "mean the mutation increases the signal; ",
                        ui.strong("negative (downward) "),
                        "means it decreases it.",
                        class_="small",
                    ),
                    ui.tags.li(
                        "Clusters of tall letters often correspond to transcription factor "
                        "binding motifs or other regulatory elements.",
                        class_="small",
                    ),
                    class_="mb-1",
                ),
                ui.p(
                    "ISM uses 16 KB context and a CenterMaskScorer (501 bp window, DIFF_MEAN). "
                    "Wider ISM regions take longer to compute.",
                    class_="small text-muted mb-0",
                ),
                class_="bg-info bg-opacity-10 p-3 rounded mb-3",
            ),
            _download_bar(
                ui.download_button("dl_ism_png", "PNG", class_="btn-outline-secondary btn-sm"),
            ),
            ui.output_ui("plot_ism_ui"),
        ),
        ui.nav_panel(
            "Multi-variant",
            ui.div(
                ui.h6("How to read these results", class_="mb-2"),
                ui.p(
                    "This tab compares the scoring results of multiple variants side-by-side.",
                    class_="small mb-1",
                ),
                ui.tags.ul(
                    ui.tags.li(
                        ui.strong("Heatmap: "),
                        "Rows are cell types, columns are variants. Color intensity shows the "
                        "log fold-change (LFC) — blue = REF higher, red = ALT higher. "
                        "Shows the 20 most variable cell types for the target gene.",
                        class_="small",
                    ),
                    ui.tags.li(
                        ui.strong("Grouped bar chart: "),
                        "Compares LFC across variants for the top 10 most affected cell types. "
                        "Different colors represent different variants.",
                        class_="small",
                    ),
                    ui.tags.li(
                        ui.strong("Data table: "),
                        "Full scoring results for all variants, filterable and sortable.",
                        class_="small",
                    ),
                    class_="mb-0",
                ),
                class_="bg-info bg-opacity-10 p-3 rounded mb-3",
            ),
            ui.output_ui("multi_parse_status_ui"),
            _download_bar(
                ui.download_button("dl_multi_csv", "CSV", class_="btn-outline-secondary btn-sm"),
                ui.download_button("dl_multi_heatmap_png", "Heatmap PNG", class_="btn-outline-secondary btn-sm"),
                ui.download_button("dl_multi_bar_png", "Bar chart PNG", class_="btn-outline-secondary btn-sm"),
            ),
            ui.output_ui("plot_multi_heatmap_ui"),
            ui.output_ui("plot_multi_bar_ui"),
            ui.output_data_frame("table_multi_scores"),
        ),
        ui.nav_panel(
            "Batch scores",
            ui.div(
                ui.h6("How to read these results", class_="mb-2"),
                ui.p(
                    "Batch scoring runs multiple scorers on multiple variants at once, "
                    "producing a comprehensive results table. Key columns:",
                    class_="small mb-1",
                ),
                ui.tags.ul(
                    ui.tags.li(
                        ui.strong("variant_id: "),
                        "Identifier for each variant (from VCF input or auto-generated).",
                        class_="small",
                    ),
                    ui.tags.li(
                        ui.strong("scorer: "),
                        "Which scorer was used (RNA_SEQ, DNASE, CAGE, etc.).",
                        class_="small",
                    ),
                    ui.tags.li(
                        ui.strong("gene_name: "),
                        "Which gene the score applies to (for gene-specific scorers). "
                        "Region-based scorers (ATAC, DNASE, CHIP) do not have gene annotations.",
                        class_="small",
                    ),
                    ui.tags.li(
                        ui.strong("raw_score (LFC): "),
                        "log2(ALT / REF). Positive = ALT higher, negative = REF higher.",
                        class_="small",
                    ),
                    ui.tags.li(
                        ui.strong("quantile_score: "),
                        "Rank vs ~348K common SNPs. Values near 1.0 = extreme effect.",
                        class_="small",
                    ),
                    class_="mb-1",
                ),
                ui.p(
                    "Use the table filters to explore results by variant, scorer, gene, "
                    "or cell type. Download the full table as CSV for further analysis.",
                    class_="small text-muted mb-0",
                ),
                class_="bg-info bg-opacity-10 p-3 rounded mb-3",
            ),
            ui.output_ui("batch_errors_ui"),
            ui.output_ui("batch_summary_ui"),
            _download_bar(
                ui.download_button("dl_batch_csv", "CSV", class_="btn-outline-secondary btn-sm"),
            ),
            ui.output_data_frame("table_batch_scores"),
        ),
        ui.nav_panel(
            "Sequence prediction",
            ui.div(
                ui.h6("How to read these results", class_="mb-2"),
                ui.p(
                    "This shows predicted signal tracks from a raw DNA sequence "
                    "(no genomic coordinates). Useful for synthetic or edited sequences.",
                    class_="small mb-1",
                ),
                ui.tags.ul(
                    ui.tags.li(
                        "The input sequence is center-padded with N bases to reach a "
                        "valid model length (16 KB, 100 KB, 500 KB, or 1 MB).",
                        class_="small",
                    ),
                    ui.tags.li(
                        "Track interpretation is the same as for interval predictions "
                        "(see 'Predicted tracks' tab for details).",
                        class_="small",
                    ),
                    ui.tags.li(
                        "Since there are no genomic coordinates, the x-axis shows "
                        "position within the padded sequence.",
                        class_="small",
                    ),
                    class_="mb-0",
                ),
                class_="bg-info bg-opacity-10 p-3 rounded mb-3",
            ),
            ui.output_ui("seq_info_ui"),
            _download_bar(
                ui.download_button("dl_seq_png", "PNG", class_="btn-outline-secondary btn-sm"),
            ),
            ui.output_ui("plot_seq_ui"),
        ),
        ui.nav_panel(
            "Available cell types",
            ui.p(
                "Browse all cell types / tissues available in the model. "
                "Use the filter to search.",
                class_="text-muted small",
            ),
            ui.output_data_frame("table_ontology"),
        ),
        id="main_tabs",
    )


app_ui = ui.page_sidebar(
    create_sidebar(),
    create_main_panel(),
    title="AlphaGenome Explorer",
)
