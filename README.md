# Fungal Bulk RNA-seq Transcriptomics
## Identification of transcriptional signatures associated with antifungal resistance in *Candida auris*

[![Bioinformatics](https://img.shields.io/badge/Field-Bioinformatics-2b6cb0)]()
[![Data](https://img.shields.io/badge/Data-Bulk%20RNA--seq-2f855a)]()
[![Status](https://img.shields.io/badge/Status-Analysis%20completed%3B%20documentation%20in%20progress-d69e2e)]()

## Overview

This project investigates transcriptional differences associated with antifungal resistance in *Candida auris* using bulk RNA-seq data and a computational analysis workflow.

The analysis combines transcript quantification, differential gene-expression testing, pathway-level analysis, alternative-splicing analysis, and functional annotation. A particular focus is the integration of gene-expression and splicing evidence to characterize candidate resistance-associated transcriptional patterns.

This repository describes the bulk-transcriptomics component of a broader *C. auris* bioinformatics portfolio. It provides analytical outputs that may inform later single-cell, integrative, machine-learning, and deep-learning projects. Those later analyses are **not** part of the completed work described here.

## Research objective

The main objective is to identify gene-expression and alternative-splicing patterns associated with resistant and sensitive phenotypes in the analysed *C. auris* samples.

The project addresses five questions:

1. Which genes show differential expression between resistant and sensitive groups?
2. Which biological pathways differ in enrichment or estimated activity between the groups?
3. Which alternative-splicing events are detected, and which meet the statistical significance threshold?
4. Which candidate genes are supported by multiple analytical layers, including differential expression, pathway analysis, and functional annotation?
5. How consistent are differential-expression results when transcript-quantification duplicate handling is changed?

## Dataset

The analysis uses six bulk RNA-seq samples:

| Phenotype | SRA run accessions | Samples |
|---|---|---:|
| Sensitive | `SRR10059687`, `SRR10059688`, `SRR10059689` | 3 |
| Resistant | `SRR10059690`, `SRR10059691`, `SRR10059692` | 3 |
| **Total** | **6 runs** | **6** |

The sample metadata used in the analysis are stored in `data/metadata/sample_metadata.csv`.

The sample size is small (three samples per group). Results should be interpreted as findings from this dataset, not as independently validated biomarkers or a clinically validated resistance classifier.

## Analysis workflow

### 1. Reference preparation and transcript quantification

- Reference genome and transcript-annotation preparation.
- STAR splice-aware alignment of RNA-seq reads.
- Salmon transcript quantification.
- `tximport` import and aggregation of transcript-level estimates to gene-level data.
- Construction of gene-level expression matrices for downstream analysis.

Two transcript-mapping strategies were evaluated:

- a baseline Salmon index with duplicate transcript sequences collapsed;
- a second index built with `--keepDuplicates`, followed by aggregation of exact duplicate transcript records into four explicit ambiguous transcript groups.

This comparison assessed whether duplicate transcript handling materially changed downstream differential-expression results.

### 2. Differential gene-expression analysis

Differential expression was analysed using two complementary R/Bioconductor approaches:

- **edgeR** with a quasi-likelihood framework;
- **limma-voom**.

The methods were compared using adjusted p-values, log fold changes, significance overlap, direction concordance, and correlation of effect-size estimates.

The duplicate-aware grouped analysis was retained as the primary analysis strategy, while the baseline analysis was preserved as a sensitivity comparison.

### 3. Pathway analysis

Pathway-level analysis used:

- **fgsea** for gene-set enrichment analysis;
- **GSVA** for pathway activity estimation.

Results from the two methods were compared in terms of pathway direction, statistical significance, and agreement between methods.

### 4. Alternative-splicing analysis

Alternative splicing was analysed with **rMATS** using STAR alignment files.

The analysis considered five standard event classes:

- skipped exon (SE);
- mutually exclusive exons (MXE);
- alternative 5′ splice site (A5SS);
- alternative 3′ splice site (A3SS);
- retained intron (RI).

Both junction-count (JC) and junction-count-plus-exon-count (JCEC) result sets were considered. Significant events were integrated with gene-expression results for joint interpretation.

### 5. Functional annotation and gene-set interpretation

Functional interpretation included:

- eggNOG-mapper annotation;
- KEGG pathway analysis;
- Gene Ontology (GO) over-representation analysis;
- directional GO analysis;
- integration of pathway leading-edge genes with differential-expression results;
- prioritization and annotation audit of a core candidate-gene set.

The core-gene analysis produced two non-overlapping groups: energy/carbon-metabolism genes and ribosome/RNA-related genes.

## Key results

The following results summarize the current analysis outputs. They are dataset-specific and are not independent biological validation.

### Differential expression

In the duplicate-aware grouped analysis:

| Method | Significant gene groups (FDR < 0.05) | Up in resistant | Down in resistant |
|---|---:|---:|---:|
| edgeR | 1,325 | 707 | 618 |
| limma-voom | 1,214 | 672 | 542 |

The two methods identified 1,210 significant groups in common. Their log-fold-change estimates were highly correlated (Spearman's rho = 0.9973), and the direction of change agreed for 5,274 of 5,376 tested groups (98.1%).

The comparison between baseline and duplicate-aware analyses showed that the 5,372 non-ambiguous biological genes had unchanged significance classification and direction. This supports the robustness of the main differential-expression findings to the duplicate-transcript handling strategy used here.

### Pathway analysis

In the grouped analysis:

- fgsea tested 228 pathways and identified 13 significant pathways (FDR < 0.05);
- GSVA tested 228 pathways and identified 6 significant pathways;
- five pathways were significant in both methods;
- pathway-level effect estimates were positively correlated (Spearman's rho = 0.8554).

Agreement between methods supports prioritization of shared pathway-level patterns, while method-specific findings remain exploratory.

### Alternative splicing

The rMATS analysis detected:

| Event type | Events reported | Significant events (FDR < 0.05) |
|---|---:|---:|
| SE | 11 | 1 |
| MXE | 2 | 2 |
| A3SS | 0 | 0 |
| A5SS | 0 | 0 |
| RI | 0 | 0 |

The significant-event set therefore contains three biological events across two genes. JC and JCEC outputs represent two counting approaches for the same event set and should not be counted as six independent biological events.

Interpretation is cautious: the annotation contains one transcript per gene, which limits isoform-level interpretation, and the observed effect sizes and junction support require careful consideration. The splicing results are candidate findings rather than experimentally validated splice switches.

### Core candidate genes

The integrated GO-based prioritization produced **43 unique core genes**:

- 15 genes associated with energy/carbon metabolism;
- 28 genes associated with ribosome biogenesis and RNA-related processes.

All 43 genes were significant in both edgeR and limma-voom and had concordant log-fold-change directions in the audited grouped analysis. This is an internal consistency check of the selected gene set, not independent external validation.

The annotation audit found named functional product annotations for 25 genes; 18 had hypothetical or uncharacterized product descriptions in the initial annotation. A subsequent eggNOG/protein-reference audit recovered additional annotation evidence for the 18 genes requiring further annotation.

## Software and methods

The analysis uses Linux/WSL, Bash, Python, R, and command-line bioinformatics tools.

### Main tools

| Tool | Role | Recorded version |
|---|---|---:|
| FastQC | Read-quality assessment | 0.12.1 |
| fastp | Read preprocessing / quality filtering | 0.23.4 |
| STAR | Splice-aware alignment executable | 2.7.11b |
| STAR index | Reference index build version | 2.7.4a |
| Salmon | Transcript quantification | 2.8.0 |
| tximport | Transcript-to-gene summarization | 1.38.2 |
| edgeR | Differential expression | 4.10.7 |
| limma | Differential expression with voom | 3.68.5 |
| fgsea | Gene-set enrichment | 1.36.2 |
| GSVA | Pathway activity analysis | 2.4.9 |
| rMATS | Alternative-splicing analysis | 4.4.0 |
| SAMtools | BAM inspection and validation | 1.24 |
| BEDTools | Genomic interval operations | 2.31.1 |
| MultiQC | Aggregated quality-control reporting | 1.33 |
| R | Statistical analysis environment | 4.5.3 |
| Python | Scripting environment | 3.11.x |

The STAR executable version and the STAR index build version are reported separately because the existing reference index was generated with STAR 2.7.4a, while the currently recorded executable is STAR 2.7.11b.

### Main R packages

`tximport`, `edgeR`, `limma`, `fgsea`, `GSVA`, `ComplexHeatmap`, `circlize`, `ggplot2`, `data.table`, `GenomicRanges`, `AnnotationDbi`, and `GO.db`.

### Main Python libraries

`pandas`, `numpy`, and `pysam`, alongside Python standard-library modules.

Exact environment files and installation instructions are being finalized. The versions above describe the recorded analysis environment and should not be taken as a guarantee that every package is installed in a fresh environment.

## Repository structure

The repository is being prepared to separate source code, metadata, reference-file documentation, analysis outputs, figures, and project documentation.

```text
.
├── README.md
├── .gitignore
├── environment.yml
├── requirements.txt
├── config.yaml
├── data/
│   ├── metadata/
│   └── reference/
├── scripts/
├── results/
│   ├── differential_expression/
│   ├── differential_expression_grouped/
│   ├── pathway_enrichment/
│   ├── pathway_enrichment_grouped/
│   ├── go_enrichment/
│   ├── integrated_analysis/
│   └── rmats/
├── figures/
└── docs/
```

Large raw sequencing files, processed FASTQ files, BAM/BAI files, large reference databases, and aligner/quantifier indices are not intended to be version-controlled. They are excluded from the public repository to keep it manageable. Reference sources and acquisition instructions should be documented separately.

The final tracked-file list may differ from this illustrative structure as repository curation is completed.

## Reproducibility notes

- Run scripts from the project root unless a script explicitly documents otherwise.
- Sample-to-phenotype assignments are defined in `data/metadata/sample_metadata.csv`.
- Reference files and indexes are not all distributed with the repository; see `data/reference/README.md` when available.
- Some scripts in the current working copy still contain machine-specific absolute paths. These must be converted to portable project-relative paths before the repository is considered fully reproducible.
- Environment specifications and dependency versions are being finalized. Do not assume that a clean environment can reproduce every output until those files and execution instructions have been validated.
- The analysis uses six samples only. No independent validation cohort is included in this project.
- Pathway enrichment and GO analyses are downstream interpretations of the same expression-derived data and are not independent experimental validation.
- rMATS results are exploratory and require cautious interpretation given the annotation structure and limited sample size.

## Limitations

1. **Small sample size:** three samples per phenotype limit statistical power and generalizability.
2. **No independent validation cohort:** candidate signatures have not been validated on an independent dataset in this project.
3. **Annotation limitations:** the reference annotation contains one transcript per gene, constraining transcript-isoform interpretation.
4. **Alternative-splicing evidence:** significant rMATS events are computational candidates and have not been experimentally validated.
5. **Functional annotation uncertainty:** a proportion of candidate genes remain hypothetical or incompletely characterized.
6. **No ML/DL model training in this project:** predictive modelling is planned for later portfolio stages and is not claimed as a result here.

## Project status

The main bulk RNA-seq computational analyses and internal result audits have been performed. Repository documentation, portable paths, environment specifications, and GitHub curation are being finalized.

This repository is a research portfolio project and computational analysis, not a clinically validated diagnostic or resistance-prediction tool.

## Acknowledgements and data provenance

The analysis uses public sequencing run accessions listed in the Dataset section. The final repository should include the original study accession, source links, reference-genome accession/version, annotation provenance, and any applicable data-use or citation requirements once these have been verified and documented.

## License

A license has not yet been selected. Until a license is added, reuse and redistribution permissions should not be assumed.
