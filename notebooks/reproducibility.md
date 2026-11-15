# Reproducibility and computational environment

## Environment
The project was analysed in Ubuntu on WSL using Conda, Bash, Python, R/Bioconductor, and command-line bioinformatics tools.

Two principal Conda environments were used:
- `p2_ngs`: primary NGS, quantification, and R analysis environment;
- `p2_rmats`: separate rMATS environment.

Recorded versions:
- Python 3.11.16 in `p2_ngs`;
- Python 3.11.18 in `p2_rmats`;
- R 4.5.3;
- Nextflow 26.04.6 (installed, but no `main.nf` or `nextflow.config` was found in the project audit);
- rMATS 4.4.0.

## Command-line tool versions
| Tool | Version |
|---|---:|
| FastQC | 0.12.1 |
| fastp | 0.23.4 |
| Cutadapt | 5.2 |
| HISAT2 | 2.2.3 |
| Salmon | 2.8.0 |
| STAR executable | 2.7.11b |
| STAR index build | 2.7.4a |
| SAMtools | 1.24 |
| BEDTools | 2.31.1 |
| featureCounts | 2.0.6 |
| MultiQC | 1.33 |
| SRA Toolkit | 3.4.1 |
| BLASTN | 2.12.0+ |
| SeqKit | 2.3.0 |
| rMATS | 4.4.0 |

The STAR index was built with 2.7.4a; the later executable was 2.7.11b. These versions must remain distinguished.

## R package versions
| Package | Version |
|---|---:|
| tximport | 1.38.2 |
| edgeR | 4.10.7 |
| limma | 3.68.5 |
| fgsea | 1.36.2 |
| GSVA | 2.4.9 |
| ComplexHeatmap | 2.26.1 |
| circlize | 0.4.18 |
| ggplot2 | 4.0.3 |
| data.table | 1.18.6.1 |
| GenomicRanges | 1.62.1 |
| AnnotationDbi | 1.72.0 |
| GO.db | 3.22.0 |

DESeq2 and rtracklayer were not installed in the audited `p2_ngs` R environment. Do not list them as dependencies of the completed workflow unless the scripts are changed.

## Python libraries
Third-party imports used by the audited scripts include pandas 3.0.6, numpy 2.4.6, and pysam 0.24.1 in `p2_ngs` (pysam 0.24.0 in `p2_rmats`). Scripts also use standard-library modules such as `collections`, `csv`, `json`, `pathlib`, and `re`.

A raw `pip freeze` export contained local file URLs and should not be copied unedited into a portable `requirements.txt`.

## Reproduction prerequisites
A fresh clone does not contain all data and reference assets required to reproduce every result. A reproduction attempt will require:
1. the six source SRA runs;
2. the matching *C. auris* reference genome and annotation;
3. the required eggNOG-mapper database;
4. locally rebuilt STAR and Salmon indices;
5. Conda environments and R/Python dependencies;
6. sample metadata and transcript-to-gene mapping files;
7. scripts executed in dependency order as documented in `pipeline_overview.md`.

The exact reference accession/release, database release, download URLs, and checksums must be completed from original acquisition records, not guessed.

## Portability status
The audit found machine-specific absolute paths in several scripts, including scripts 02–10 and 20–24. These must be replaced with project-root discovery or a documented configuration mechanism before claiming that a fresh clone can run the workflow.

Before release:
- verify all referenced inputs exist in the documented layout;
- create output directories safely;
- avoid overwriting existing results without an explicit decision;
- test scripts from the project root in the intended environment;
- document dependencies on data excluded from Git.

## Data and index distribution
Large raw/trimmed FASTQ files, BAM/BAI files, eggNOG databases, and generated STAR/Salmon indices are not intended for Git version control. Document their sources and regeneration steps instead.

## Limitations
The sample size is three per phenotype; the repository is not yet a fully tested one-command workflow; some scripts depend on previous outputs; reference/database provenance and checksums need completion; machine-specific paths need remediation; and no independent validation cohort is included.
