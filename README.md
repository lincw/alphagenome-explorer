# AlphaGenome Explorer

A [Shiny for Python](https://shiny.posit.co/py/) web app for interactively exploring [AlphaGenome](https://deepmind.google/discover/blog/alphagenome-a-unified-model-of-the-human-genome/), Google DeepMind's DNA foundation model that predicts regulatory genomic signals from sequence alone.

## What it does

Given a gene symbol and a cell type, the app lets you:

1. **View gene context** -- transcript structure (exons, introns, strand) from GENCODE v46.
2. **Predict signal tracks** -- predicted RNA-seq, CAGE, DNase, ATAC, ChIP, splicing, and 3D contact maps across the gene region for any supported cell type.
3. **Analyze variant effects** -- compare REF vs ALT allele predictions side-by-side with overlaid tracks and ALT-REF difference tracks. REF allele is auto-detected from hg38 via the UCSC REST API.
4. **Score variants** -- quantify gene-level expression change (log-fold-change) and compare the variant's impact across all cell types with a bar chart showing the most and least affected tissues.
5. **In silico mutagenesis (ISM)** -- systematically mutate every base in a short region to identify functionally important positions, visualized as a sequence logo.
6. **Multi-variant comparison** -- enter multiple variants (disease vs. benign, etc.) and compare their scoring results side-by-side with heatmaps and grouped bar charts.
7. **Batch variant scoring** -- paste VCF-style data to score many variants with multiple scorers at once, producing a comprehensive downloadable results table.
8. **Export results** -- download predicted tracks (CSV), score tables (CSV), and all plots (PNG) from each results tab.

The app supports all 11 AlphaGenome output modalities:

| # | Key | Description |
|---|-----|-------------|
| 1 | RNA_SEQ | Gene expression (total & polyA RNA-seq read coverage) |
| 2 | CAGE | TSS usage at single-nucleotide resolution |
| 3 | PROCAP | Nascent RNA TSS mapping |
| 4 | DNASE | Chromatin accessibility (DNase I hypersensitivity) |
| 5 | ATAC | Chromatin accessibility (Tn5 transposase) |
| 6 | CHIP_HISTONE | Histone modification signals (H3K27ac, H3K4me3, etc.) |
| 7 | CHIP_TF | Transcription factor binding (CTCF, etc.) |
| 8 | SPLICE_SITES | Splice donor/acceptor site probability |
| 9 | SPLICE_SITE_USAGE | Cell-type-specific splice site usage frequency |
| 10 | SPLICE_JUNCTIONS | Exon-exon junction counts (sashimi plots) |
| 11 | CONTACT_MAPS | 3D chromatin organization (Micro-C / Hi-C-like) |

## Project structure

```
colabs/
  app.py            Entry point (run with: shiny run app.py)
  config.py         Constants, API key, output type definitions
  ui_layout.py      UI layout (sidebar + tabbed main panel)
  server.py           Server logic (reactive handlers, Steps 1–7)
  server_downloads.py Download handlers (CSV + PNG export)
  server_multi.py     Multi-variant comparison handlers (Step 8)
  server_batch.py     Batch variant scoring handlers (Step 9)
  server_sequence.py  Predict-from-sequence handlers (Step 10)
  batch_utils.py      VCF parsing, multi-scorer batch scoring
  sequence_utils.py   DNA sequence validation, padding, statistics
  plot_utils.py       Plot component builders for interval & variant predictions
  score_utils.py      Scoring display (gene-grouped summary, comparison bar chart)
  variant_utils.py    rsID lookup (NCBI dbSNP) & reference sequence fetch (UCSC)
  ism_utils.py        In silico mutagenesis helpers (ISM scoring + SeqLogo)
  multi_variant_utils.py  Multi-variant parsing, batch scoring, heatmap/bar charts
  download_utils.py   CSV/PNG export helpers
  environment.yml     Conda environment specification
```

## Setup

### 1. Create the conda environment

```bash
conda env create -f environment.yml
conda activate alphagenome
```

### 2. Get an API key

The app uses the AlphaGenome API. You'll be prompted to enter your key when the app starts.

1. Go to [AlphaGenome](https://www.alphagenomedocs.com/index.html) and click "Get API key".
2. Paste it into the "API Key" field in the app sidebar.

### 3. Run the app

```bash
shiny run app.py
```

The app will open at `http://localhost:8000` by default.

**Note:** On the first run, the app will automatically download a ~300MB GENCODE v46 annotation file (`gencode_v46_annotation.feather`) from Google Cloud Storage. This file is cached locally for future use.

## Usage workflow

1. **Step 1 -- Gene lookup:** Enter a gene symbol (e.g. `OAS1`, `DLG1`, `BRCA1`) and click "Look up gene".
2. **Step 2 -- Output types:** Select which modalities to predict (RNA_SEQ is selected by default).
3. **Step 3 -- Cell type:** Search for a cell type or tissue (e.g. `monocyte`, `lung`, `K562`) and select from results.
4. **Step 4 -- Predict tracks:** Click "Run prediction" to generate signal tracks across the gene region.
5. **Step 5 -- Variant analysis (optional):** Enter an **rsID** (e.g. `rs35482426`) and click "Look up rsID" to auto-fill position, REF, and ALT from dbSNP. Or enter coordinates manually -- the REF allele is auto-detected from hg38. Click "Predict variant effect" to see overlaid REF/ALT tracks and difference tracks.
6. **Step 6 -- Variant scoring (optional):** Choose a scorer (e.g. RNA_SEQ) and click "Score variant" to get gene-level log-fold-change values and a cross-cell-type comparison.
7. **Step 7 -- ISM (optional):** Enter a center position (or leave blank to use the variant position / gene center), choose an ISM width and output type, then click "Run ISM" to see a sequence logo of base-level importance.
8. **Step 8 -- Multi-variant comparison (optional):** Enter multiple variants (one per line, `chr:pos:ref:alt` or `rsID`, optionally followed by a label) and click "Compare variants" to see a heatmap and grouped bar chart of their effects side-by-side.
9. **Step 9 -- Batch scoring (optional):** Paste VCF-style data (tab or comma separated, with header: `variant_id, CHROM, POS, REF, ALT`). Select which scorers to run. Each variant is scored independently with its own interval. Download the full results table as CSV.
10. **Step 10 -- Predict from sequence (optional):** Paste a raw DNA sequence (A/C/G/T/N, up to 1 MB) to predict signals without genomic coordinates. The sequence is center-padded with N to the nearest valid model length (16 KB / 100 KB / 500 KB / 1 MB). Useful for synthetic or edited sequences. Requires output types (Step 2) and cell type (Step 3).

## Interpreting results

### Predicted tracks
- Peaks in RNA-seq / CAGE / PROCAP indicate active transcription or TSS usage.
- Peaks in DNase / ATAC mark open chromatin (promoters, enhancers).
- ChIP peaks show histone modifications or TF binding sites.
- Sashimi arcs connect spliced exon pairs; arc height = junction read count.
- Contact maps show predicted 3D chromatin interactions.

### Variant effect
- **Overlaid tracks:** Gray = REF, Red = ALT. Divergence = variant-induced change.
- **Difference tracks (ALT - REF):** Positive = ALT stronger, Negative = REF stronger.

### Variant scores
- **LFC (log-fold-change):** log2(ALT / REF) of gene-body signal. LFC = 0.01 means ~1% increase.
- **Quantile:** How extreme the change is vs. ~348K common SNPs. Quantile = 1.000 means this variant's effect is larger than nearly all common variants.
- The summary groups results by gene, sorted by strongest effect first.
- The bar chart shows the top 10 most increased + top 10 most decreased cell types.

### ISM (sequence logo)
- Each position shows letters (A, C, G, T) whose height reflects the magnitude of that mutation's effect.
- Tall letters = functionally important bases (mutations there cause large signal changes).
- Clusters of tall letters often correspond to transcription factor binding motifs.
- Uses 16 KB context and CenterMaskScorer (501 bp window) for efficiency.

### Multi-variant comparison
- **Heatmap:** Rows = cell types, columns = variants. Color = LFC (blue = REF higher, red = ALT higher).
- **Grouped bar chart:** Compares LFC across variants for the top most-affected cell types.
- Supports up to 10 variants per comparison. Input format: `chr:pos:ref:alt label` or `rsID label`.

### Predict from sequence
- Input sequence is center-padded with N bases to a valid model length.
- The info bar shows the original sequence length, GC content, and padded length.
- Track interpretation is identical to interval predictions (peaks = active signal).
- The x-axis shows position within the padded sequence (no genomic coordinates).

### Batch scoring
- Scores many variants (up to 20) with multiple scorers simultaneously.
- Each variant uses its own genomic interval (centered on the variant position).
- Variants can be on different chromosomes — not tied to a gene lookup.
- Results include raw LFC scores and quantile scores for all variant/scorer/gene/cell-type combinations.
- Use the table filters to explore, or download the full CSV for external analysis.

### Exporting results
- Each results tab has download buttons for CSV (data) and PNG (plots).
- Track CSV is downsampled to 128 bp resolution for sequences > 100 KB.

## Dependencies

- Python 3.12
- [alphagenome](https://pypi.org/project/alphagenome/) (Google DeepMind)
- [shiny](https://shiny.posit.co/py/) (Posit)
- numpy, pandas, matplotlib, pyarrow
