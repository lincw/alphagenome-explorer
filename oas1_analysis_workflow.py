# %%
# =============================================================================
# OAS1 Locus Analysis Workflow using AlphaGenome
# =============================================================================
# Goal: Explore how genomic variants near the OAS1 locus affect gene expression.
#
# OAS1 (2'-5'-Oligoadenylate Synthetase 1) is a key innate immune gene on
# chromosome 12 (chr12:112,906,962-112,933,219, + strand, GRCh38). It encodes
# an enzyme activated by double-stranded RNA during viral infection, playing a
# critical role in antiviral defense. Variants in OAS1 have been associated with
# COVID-19 severity, susceptibility to viral infections, and differential
# interferon response.
#
# Key activities:
#   1. Visualize the genomic context of known OAS1 variants.
#   2. Predict the functional impact of a specific variant on gene expression,
#      accessibility, and histone marks.
#   3. Systematically compare the predicted effects of known disease-associated
#      variants to a set of background variants.
#
# Based on the TAL1 analysis workflow from AlphaGenome:
# https://www.alphagenomedocs.com/colabs/example_analysis_workflow.html
# =============================================================================

# %% Install AlphaGenome (uncomment if running in Colab)
# from IPython.display import clear_output
# ! pip install alphagenome
# clear_output()

# %% Imports
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend: save only, no preview window

import io
import itertools

from alphagenome import colab_utils
from alphagenome.data import gene_annotation
from alphagenome.data import genome
from alphagenome.data import transcript as transcript_utils
from alphagenome.models import dna_client
from alphagenome.models import variant_scorers
from alphagenome.visualization import plot_components
from IPython.display import clear_output
import numpy as np
import pandas as pd
import plotnine as gg

# %% Initialize the DNA model
secret = ''  # Paste your API key here (get one from https://deepmind.google.com/science/alphagenome/)
dna_model = dna_client.create(secret)

# %% Create output directory for figures
import os
output_dir = os.path.join(os.path.dirname(__file__), 'oas1_figures')
os.makedirs(output_dir, exist_ok=True)
print(f"Figures will be saved to: {output_dir}")

# %% Load gene annotations from GENCODE (cache locally to avoid re-downloading)
gtf_cache_path = os.path.join(os.path.dirname(__file__), 'gencode_v46_annotation.feather')
if os.path.exists(gtf_cache_path):
    print("Loading GTF from local cache...")
    gtf = pd.read_feather(gtf_cache_path)
else:
    print("Downloading GTF from GENCODE (first time only)...")
    gtf = pd.read_feather(
        'https://storage.googleapis.com/alphagenome/reference/gencode/'
        'hg38/gencode.v46.annotation.gtf.gz.feather'
    )
    gtf.to_feather(gtf_cache_path)
    print(f"Saved GTF cache to {gtf_cache_path}")

gtf_transcript = gene_annotation.filter_protein_coding(gtf)
gtf_transcript = gene_annotation.filter_to_mane_select_transcript(
    gtf_transcript
)
transcript_extractor = transcript_utils.TranscriptExtractor(gtf_transcript)


# =============================================================================
# Utility functions
# =============================================================================

# %% Utilities
def generate_background_variants(
    variant: genome.Variant, max_number: int = 100
) -> pd.DataFrame:
    """Generates background variants by creating random sequences of the same
    length as the alternate allele at the same position. This lets us test
    whether the specific variant sequence has a greater effect than a random
    sequence of the same length."""
    nucleotides = np.array(list('ACGT'), dtype='<U1')

    def generate_unique_strings(n, max_number, random_seed=42):
        rng = np.random.default_rng(random_seed)
        if 4**n < max_number:
            raise ValueError(
                'Cannot generate that many unique strings for the given length.'
            )
        generated_strings = set()
        while len(generated_strings) < max_number:
            indices = rng.integers(0, 4, size=n)
            new_string = ''.join(nucleotides[indices])
            if new_string != variant.alternate_bases:
                generated_strings.add(new_string)
        return list(generated_strings)

    permutations = []
    if 4 ** len(variant.alternate_bases) < max_number:
        for p in itertools.product(
            nucleotides, repeat=len(variant.alternate_bases)
        ):
            permutations.append(''.join(p))
    else:
        permutations = generate_unique_strings(
            len(variant.alternate_bases), max_number
        )

    ism_candidates = pd.DataFrame({
        'ID': ['mut_' + str(variant.position) + '_' + x for x in permutations],
        'CHROM': variant.chromosome,
        'POS': variant.position,
        'REF': variant.reference_bases,
        'ALT': permutations,
        'output': 0.0,
        'original_variant': variant.name,
    })
    return ism_candidates


def disease_and_background_variants(
    input_sequence_length: int, number_of_background_variants: int = 20
) -> pd.DataFrame:
    """Generates a dataframe combining disease-associated and shuffled
    background variants for OAS1."""
    disease_variants = oas1_disease_variants()

    variants = []
    for vcf_row in disease_variants.itertuples():
        variants.append(
            genome.Variant(
                chromosome=str(vcf_row.CHROM),
                position=int(vcf_row.POS),
                reference_bases=vcf_row.REF,
                alternate_bases=vcf_row.ALT,
                name=vcf_row.ID,
            )
        )

    background_variants = pd.concat([
        generate_background_variants(variant, number_of_background_variants)
        for variant in variants
    ])
    all_variants = pd.concat([disease_variants, background_variants])
    return inference_df(
        all_variants, input_sequence_length=input_sequence_length
    )


def vcf_row_to_variant(vcf_row: pd.Series) -> genome.Variant:
    """Parse a row of a VCF dataframe into a genome.Variant."""
    variant = genome.Variant(
        chromosome=str(vcf_row.CHROM),
        position=int(vcf_row.POS),
        reference_bases=vcf_row.REF,
        alternate_bases=vcf_row.ALT,
        name=vcf_row.ID,
    )
    return variant


def inference_df(
    qtl_df: pd.DataFrame,
    input_sequence_length: int,
) -> pd.DataFrame:
    """Returns a DataFrame with variants and intervals ready for inference."""
    df = []
    for _, row in qtl_df.iterrows():
        variant = vcf_row_to_variant(row)
        interval = genome.Interval(
            chromosome=row['CHROM'], start=row['POS'], end=row['POS']
        ).resize(input_sequence_length)
        df.append({
            'interval': interval,
            'variant': variant,
            'output': row['output'],
            'variant_id': row['ID'],
            'POS': row['POS'],
            'REF': row['REF'],
            'ALT': row['ALT'],
            'CHROM': row['CHROM'],
        })
    return pd.DataFrame(df)


def variant_region_groups(eval_df):
    """Groups variants by genomic region (promoter, splice site, coding,
    3'UTR) for comparative visualization."""
    grp = []
    for row in eval_df.itertuples():
        pos = row.POS
        if pos < 112907000:
            grp.append('upstream')
        elif pos < 112910000:
            grp.append('promoter_5utr')
        elif pos < 112920000:
            grp.append('intronic_splice')
        elif pos < 112930000:
            grp.append('coding')
        else:
            grp.append('3utr_downstream')
    grp = pd.Series(grp)
    return pd.Categorical(grp, categories=sorted(grp.unique()), ordered=True)


# =============================================================================
# Assembling variant data
# =============================================================================

# %% Define OAS1 disease-associated variants
def oas1_disease_variants() -> pd.DataFrame:
    """Returns a dataframe of known disease-associated variants affecting OAS1.

    These variants have been identified through GWAS, eQTL studies, and
    functional genomics analyses. Key references:
    - Zeberg & Paabo 2021 (PNAS): Neanderthal OAS1 haplotype protective
      against COVID-19
    - Pairo-Castineira et al. 2021 (Nature): COVID-19 GWAS
    - Zhou et al. 2021 (Nature Medicine): COVID-19 severity
    - GTEx eQTL data: Expression quantitative trait loci in whole blood

    Variant positions are on GRCh38/hg38.
    NOTE: Verify exact positions against current dbSNP before use.
    """
    variant_data = """\
ID	CHROM	POS	REF	ALT	output	Category	Source
rs10774671	chr12	112919388	G	A	1	splice_site	Zeberg_Paabo_2021
rs1131454	chr12	112919637	G	A	1	coding	COVID19_HGI_2021
rs2660	chr12	112930876	A	G	1	3utr	GTEx_eQTL
rs4767027	chr12	112907944	T	C	1	intronic_eQTL	GTEx_eQTL
rs2285933	chr12	112909893	G	A	1	promoter_regulatory	Pairo-Castineira_2021
rs2057778	chr12	112918360	C	T	1	intronic	Zhou_2021
rs1042856	chr12	112920171	C	T	1	coding_missense	COVID19_HGI_2021
rs2249243	chr12	112908498	A	G	1	intronic_eQTL	GTEx_eQTL
rs3741981	chr12	112919524	G	A	1	coding_synonymous	Zeberg_Paabo_2021
rs2854600	chr12	112906270	G	A	1	upstream_regulatory	ENCODE_regulatory
"""
    return pd.read_table(io.StringIO(variant_data), sep='\t')


print("OAS1 disease-associated variants:")
print(oas1_disease_variants().head(10))


# =============================================================================
# Visualize variant positions in genomic context
# =============================================================================

# %% Visualise variant positions
# Define the OAS1 interval (gene + flanking regions for context)
oas1_interval = genome.Interval(
    chromosome='chr12', start=112904000, end=112935000, strand='+'
)

# Gather unique variant positions
unique_positions = oas1_disease_variants()['POS'].unique()
unique_positions.sort()

# Define labels (manually curated to avoid overplotting)
labels = [
    'rs2854600\n(upstream)',
    'rs4767027\n(intronic eQTL)',
    'rs2249243',
    'rs2285933\n(promoter)',
    'rs2057778\n(intronic)',
    'rs10774671\n(splice)',
    'rs3741981',
    'rs1131454\n(coding)',
    'rs1042856\n(missense)',
    'rs2660\n(3\'UTR)',
]

# Build visualization
fig_variant_positions = plot_components.plot(
    [
        plot_components.TranscriptAnnotation(
            transcript_extractor.extract(oas1_interval)
        ),
    ],
    annotations=[
        plot_components.VariantAnnotation(
            [
                genome.Variant(
                    chromosome='chr12',
                    position=x,
                    reference_bases='N',
                    alternate_bases='N',
                )
                for x in unique_positions
            ],
            labels=labels,
            use_default_labels=False,
        )
    ],
    interval=oas1_interval,
    title='Positions of disease-associated variants near OAS1',
)
fig_variant_positions.savefig(
    os.path.join(output_dir, '01_variant_positions.png'), dpi=200, bbox_inches='tight'
)
print(f"Saved: 01_variant_positions.png")


# =============================================================================
# Exploring individual variant outcomes
# =============================================================================

# %% Define the variant of interest
# rs10774671 (G>A) is the key splice-site variant associated with COVID-19
# severity. It alters the splice acceptor site of exon 7, producing a shorter
# p42 isoform (A allele) instead of the longer, more enzymatically active p46
# isoform (G allele).
variant = vcf_row_to_variant(oas1_disease_variants().iloc[0])
print(f"Variant of interest: {variant}")

# %% Discover available ontology terms for RNA-seq
# List available RNA-seq ontology terms so we can pick a relevant cell type.
rna_seq_metadata = dna_model.output_metadata(
    dna_client.Organism.HOMO_SAPIENS
).rna_seq
print("Available RNA-seq ontology terms (first 30):")
print(rna_seq_metadata[['ontology_curie', 'biosample_name']].drop_duplicates().head(30).to_string())

# %% Define the ontology context
# Pick an immune-relevant cell type available in the model.
# Preferred order: monocyte > whole blood > B cell > CD34+ HSC
preferred_ontologies = [
    'CL:0000576',   # monocyte
    'CL:0000236',   # B cell
    'CL:0000624',   # CD4-positive T cell
    'CL:0000084',   # T cell
    'CL:0001059',   # common myeloid progenitor, CD34-positive
    'UBERON:0000178',  # blood
]
available_curies = set(rna_seq_metadata['ontology_curie'].unique())
ontology_terms = None
for curie in preferred_ontologies:
    if curie in available_curies:
        ontology_terms = [curie]
        break
if ontology_terms is None:
    # Fallback: use the first available term
    ontology_terms = [rna_seq_metadata['ontology_curie'].iloc[0]]
print(f"\nSelected ontology term: {ontology_terms[0]}")

# %% Make predictions for reference and alternate alleles
output = dna_model.predict_variant(
    interval=oas1_interval.resize(2**20),
    variant=variant,
    requested_outputs={
        dna_client.OutputType.RNA_SEQ,
        dna_client.OutputType.CHIP_HISTONE,
        dna_client.OutputType.DNASE,
    },
    ontology_terms=ontology_terms,
)

# %% Visualize variant effect on predicted tracks
transcripts = transcript_extractor.extract(oas1_interval)
fig_variant_effect = plot_components.plot(
    [
        plot_components.TranscriptAnnotation(transcripts),
        # RNA-seq tracks (ALT - REF difference)
        plot_components.Tracks(
            tdata=output.alternate.rna_seq.filter_to_positive_strand()
            - output.reference.rna_seq.filter_to_positive_strand(),
            ylabel_template='{biosample_name} ({strand})\n{name}',
            filled=True,
        ),
        # DNase tracks
        plot_components.Tracks(
            tdata=output.alternate.dnase.filter_to_positive_strand()
            - output.reference.dnase.filter_to_positive_strand(),
            ylabel_template='{biosample_name} ({strand})\n{name}',
            filled=True,
        ),
        # ChIP histone tracks
        plot_components.Tracks(
            tdata=output.alternate.chip_histone.filter_to_positive_strand()
            - output.reference.chip_histone.filter_to_positive_strand(),
            ylabel_template='{biosample_name} ({strand})\n{name}',
            filled=True,
        ),
    ],
    annotations=[plot_components.VariantAnnotation([variant])],
    interval=oas1_interval,
    title=(
        'Effect of rs10774671 (G>A) on predicted RNA Expression, DNase, and'
        f' ChIP-Histone in monocytes.\n{variant=}'
    ),
)
fig_variant_effect.savefig(
    os.path.join(output_dir, '02_rs10774671_variant_effect.png'), dpi=200, bbox_inches='tight'
)
print(f"Saved: 02_rs10774671_variant_effect.png")


# =============================================================================
# Batch variant scoring: comparing disease variants vs. background
# =============================================================================

# %% Prepare variant groups with background variants
eval_df = disease_and_background_variants(
    input_sequence_length=2**20, number_of_background_variants=3
)

# Add annotations
eval_df['ALT_len'] = eval_df['ALT'].str.len()
eval_df['variant_group'] = (
    eval_df['POS'].astype(str) + '_' + eval_df['ALT_len'].astype(str)
)
eval_df['output'] = eval_df['output'].fillna(0) != 0
eval_df['region_group'] = variant_region_groups(eval_df)

# %% Score all variants using RNA-seq predictions
scores = dna_model.score_variants(
    intervals=eval_df['interval'].to_list(),
    variants=eval_df['variant'].to_list(),
    variant_scorers=[variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ']],
    max_workers=2,
)
clear_output()

# %% Extract OAS1 expression scores
# Find the index corresponding to the OAS1 gene
gene_index = scores[0][0].obs.query('gene_name == "OAS1"').index[0]

# Find the index for our selected cell type
selected_curie = ontology_terms[0]
scorer_var = scores[0][0].var
print(f"Available ontology_curie values in scorer output:\n{scorer_var['ontology_curie'].unique()}")

cell_type_matches = scorer_var.query(f'ontology_curie == "{selected_curie}"')
if len(cell_type_matches) == 0:
    # Fallback: use the first available cell type in the scorer output
    print(f"Warning: {selected_curie} not found in scorer output. Using first available.")
    cell_type_index = scorer_var.index[0]
else:
    cell_type_index = cell_type_matches.index[0]
print(f"Using cell type index: {cell_type_index} ({scorer_var.loc[cell_type_index].to_dict()})")


def get_oas1_score(score_data):
    """Extracts the OAS1 expression score for the selected cell type."""
    return score_data[gene_index, cell_type_index].X[0, 0]


eval_df['oas1_diff'] = [get_oas1_score(x[0]) for x in scores]

print("\nScore summary:")
print(eval_df[['variant_id', 'output', 'oas1_diff', 'region_group']].head(20))


# =============================================================================
# Comparative visualization: disease variants vs. background
# =============================================================================

# %% Prepare plotting data
plot_df = eval_df.loc[eval_df.REF != eval_df.ALT].copy()
plot_df['variant'] = plot_df['variant'].astype(str)

plot_df = plot_df.loc[
    :,
    [
        'variant',
        'output',
        'oas1_diff',
        'region_group',
    ],
].drop_duplicates()

facet_title_by_group = {
    '3utr_downstream': 'chr12: 3\'UTR /\ndownstream',
    'coding': 'chr12: coding\nregion',
    'intronic_splice': 'chr12: intronic /\nsplice site',
    'promoter_5utr': 'chr12: promoter /\n5\'UTR',
    'upstream': 'chr12: upstream\nregulatory',
}

# %% Generate rain-cloud style plots for each region group
plt_dict = {}
for group in plot_df.region_group.unique():
    subplot_df = pd.concat(
        [plot_df.assign(plot_group='density'), plot_df.assign(plot_group='rain')]
    )
    subplot_df = subplot_df[subplot_df.region_group == group]
    subplot_df = subplot_df[
        ~((subplot_df.plot_group == 'density') & (subplot_df.output))
    ]

    score_range = np.ptp(subplot_df.oas1_diff)
    col_width = score_range / 200 if score_range > 0 else 0.001
    subplot_df['col_width'] = subplot_df['output'].map(
        {True: 1.5 * col_width, False: 1.25 * col_width}
    )

    plt_ = (
        gg.ggplot(subplot_df)
        + gg.aes(x='oas1_diff')
        + gg.geom_col(
            gg.aes(
                y=1,
                width='col_width',
                fill='output',
                x='oas1_diff',
                alpha='output',
            ),
            data=subplot_df[subplot_df['plot_group'] == 'rain'],
        )
        + gg.geom_density(
            gg.aes(
                x='oas1_diff',
                fill='output',
            ),
            data=subplot_df[subplot_df['plot_group'] == 'density'],
            color='white',
        )
        + gg.facet_wrap('~output + plot_group', nrow=1, scales='free_x')
        + gg.scale_alpha_manual({True: 1, False: 0.3})
        + gg.scale_fill_manual({True: '#E63946', False: 'gray'})
        + gg.labs(title=facet_title_by_group.get(group, group))
        + gg.theme_minimal()
        + gg.geom_vline(xintercept=0, linetype='dotted')
        + gg.theme(
            figure_size=(1.5, 3),
            legend_position='none',
            axis_text_x=gg.element_blank(),
            panel_grid_major_x=gg.element_blank(),
            panel_grid_minor_x=gg.element_blank(),
            strip_text=gg.element_blank(),
            axis_title_y=gg.element_blank(),
            axis_title_x=gg.element_blank(),
            plot_title=gg.element_text(size=9),
        )
        + gg.scale_y_reverse()
        + gg.coord_flip()
    )
    plt_dict[group] = plt_

clear_output()

# %% Display plots
print("Displaying variant effect plots by genomic region:")
print("Each plot shows disease variants (red) vs. background (gray)")
print("Y-axis = predicted OAS1 expression change (ALT vs REF) in monocytes\n")

for group_name, plt_ in plt_dict.items():
    fname = f"03_raincloud_{group_name}.png"
    plt_.save(os.path.join(output_dir, fname), dpi=200)
    print(f"Saved: {fname}")

print(f"\nAll figures saved to: {output_dir}")
