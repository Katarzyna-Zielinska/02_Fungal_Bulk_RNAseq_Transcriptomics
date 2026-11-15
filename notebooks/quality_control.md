# Quality control

## Scope
QC covered read-level reports, reference/annotation consistency, STAR/BAM integrity, and sample-level expression diagnostics. These checks do not replace review of every per-sample metric and plot.

## Read-level QC
FastQC 0.12.1 assessed read quality; MultiQC 1.33 aggregated reports. Preprocessing used fastp 0.23.4 and Cutadapt 5.2.

One inspected raw read file from SRR26185932 had PASS for Basic Statistics, Per Base Sequence Quality, and Adapter Content, and FAIL for Per Tile Sequence Quality and Sequence Duplication Levels. This observation applies to that inspected file, not automatically to all six samples in this project.

## Reference checks
Recorded validation:
- 15 genome contigs;
- 12,365,959 total bases;
- 15 annotation contigs;
- 5,419 transcript records and 5,890 exon records in the mRNA-only GTF;
- no missing annotation contigs;
- no invalid or out-of-bounds coordinates.

## STAR and BAM checks
The STAR index was built with STAR 2.7.4a; the later executable was 2.7.11b.

All six STAR BAM headers matched the reference FASTA:
- 15 contigs per BAM;
- zero contig-length mismatches;
- zero BAM-only or FASTA-only contigs.

All six BAMs passed `samtools quickcheck`, and all had BAI indices.

| SRA run | Properly paired |
|---|---:|
| SRR10059687 | 99.84% |
| SRR10059688 | 99.76% |
| SRR10059689 | 99.77% |
| SRR10059690 | 99.81% |
| SRR10059691 | 99.73% |
| SRR10059692 | 99.82% |

These checks establish file/header integrity and basic pairing statistics. They do not prove that every alignment is biologically correct or replace full alignment QC.

## Sample-level diagnostics
Generated diagnostic outputs include library-size comparison, sample PCA, sample-correlation heatmap, edgeR BCV plot, edgeR MDS plot, edgeR quasi-likelihood dispersion plot, and limma-voom mean-variance plot.

## Limitations
- There are three samples per phenotype.
- A failed FastQC module must be interpreted in context; it does not automatically make a sample unusable.
- `samtools quickcheck` tests basic BAM integrity, not full alignment quality.
- Proper-pair percentages do not replace inspection of mapping rate, coverage, strandedness, contamination, and other relevant metrics.
- The STAR index-build and executable versions differ and are documented separately.
