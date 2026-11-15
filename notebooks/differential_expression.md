# Differential gene-expression analysis

## Objective
Identify genes with different estimated expression between sensitive and resistant *Candida auris* samples.

## Quantification and import
Salmon 2.8.0 quantified transcripts. tximport 1.38.2 imported transcript estimates and summarized them to gene-level data. Length-scaled tximport objects were used for differential-expression analysis. Sample conditions are recorded in `data/metadata/sample_metadata.csv`.

Two strategies were compared:
- baseline Salmon index with exact duplicate transcript sequences collapsed;
- duplicate-aware Salmon index built with `--keepDuplicates`, with four exact duplicate pairs grouped into synthetic ambiguous groups.

The grouped analysis was selected as primary; baseline results were retained for sensitivity assessment.

## Methods
- edgeR 4.10.7 with a quasi-likelihood framework;
- limma 3.68.5 with voom.

Methods were compared using FDR, log fold change, significance overlap, direction concordance, and effect-size correlation. The reported significance threshold is FDR < 0.05.

## Grouped-analysis results
| Metric | edgeR | limma-voom |
|---|---:|---:|
| Input groups | 5,415 | 5,415 |
| Retained groups | 5,376 | 5,376 |
| Significant groups | 1,325 | 1,214 |
| Up in resistant | 707 | 672 |
| Down in resistant | 618 | 542 |

Comparison:
- significant in both: 1,210;
- edgeR only: 115;
- limma only: 4;
- neither: 4,047;
- direction concordance: 5,274/5,376 (98.1%);
- Spearman correlation of log fold changes: 0.997304.

## Baseline-versus-grouped sensitivity analysis
The two branches shared 5,372 non-ambiguous biological genes. For these genes, edgeR and limma-voom log-fold-change Spearman correlations were both 1.0, direction concordance was 100%, and no significance-classification changes were observed. The maximum absolute log-fold-change difference was approximately 0.000075.

The four ambiguous transcript groups are computational groups, not additional biological genes.

## Interpretation
Differential expression identifies associations with phenotype labels in this dataset; it does not establish causality. The small number of biological replicates limits power and generalizability. Independent datasets and experimental validation are needed for robust signatures.

## Outputs and scripts
Results are stored in `results/differential_expression/` and `results/differential_expression_grouped/`.

Relevant scripts: `02_tximport.R`, `04_tximport_length_scaled.R`, `05_edgeR_DE.R`, `07_limma_voom_DE.R`, `08_compare_DE_methods.R`, and grouped counterparts `20`–`24` in `scripts/`.
