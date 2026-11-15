# Methodology

## Project scope
This project investigates bulk RNA-seq transcriptional differences associated with antifungal-sensitive and antifungal-resistant phenotypes in *Candida auris*. It combines read QC and preprocessing, splice-aware alignment, transcript quantification, differential-expression analysis, pathway analysis, alternative-splicing analysis, and functional annotation.

This is a computational, exploratory analysis. It does not establish causal mechanisms, clinically validated biomarkers, or a validated resistance-prediction model. Single-cell, machine-learning, and deep-learning analyses belong to later portfolio projects and are not claimed as completed here.

## Dataset
Six public SRA runs were analysed:

| Phenotype | SRA runs | n |
|---|---|---:|
| Sensitive | SRR10059687, SRR10059688, SRR10059689 | 3 |
| Resistant | SRR10059690, SRR10059691, SRR10059692 | 3 |

Sample metadata are stored in `data/metadata/sample_metadata.csv`. Three samples per group is a small sample size; findings are dataset-specific and require independent validation.

## Reference and annotation
The project used a *C. auris* B8441 reference genome and an mRNA-focused annotation. Recorded validation found 15 genome contigs, 12,365,959 total bases, 5,419 transcript records and 5,890 exon records in the mRNA-only annotation. All annotation contigs matched the genome, with no invalid or out-of-bounds coordinates.

The mRNA-only annotation contains one transcript per gene in the records used here. This limits transcript-isoform-level interpretation. Exact source URL, assembly accession/release, annotation release, download date, and checksums must be verified from the original acquisition records before public release; they are not guessed here.

## Read processing, alignment, and quantification
Read QC used FastQC and MultiQC. Preprocessing used fastp and Cutadapt. STAR was used for splice-aware alignment, Salmon for transcript quantification, and tximport for transcript-to-gene summarization.

The STAR index was built with STAR 2.7.4a using `candida_auris_genome.fna`, `candida_auris_mRNA_only.gtf`, `sjdbOverhang=149`, `genomeSAindexNbases=10`, and four threads. The later recorded STAR executable was 2.7.11b; these versions are distinct.

Salmon 2.8.0 indices used k=21. Two strategies were compared: a baseline index with duplicate transcript sequences collapsed (5,415 references), and an index built with `--keepDuplicates` (5,419 references).

## Duplicate-transcript handling
The baseline `tx2gene.tsv` contains 5,419 transcript records and 5,419 unique gene identifiers. The duplicate-aware `tx2gene_keepDuplicates_grouped.tsv` contains 5,419 transcript records and 5,415 gene groups. Four pairs of exact duplicate sequences were grouped:

| Group | Transcript identifiers |
|---|---|
| AMBIGUOUS_GROUP_001 | B9J08_004680 / B9J08_005536 |
| AMBIGUOUS_GROUP_002 | B9J08_004683 / B9J08_005540 |
| AMBIGUOUS_GROUP_003 | B9J08_004684 / B9J08_005541 |
| AMBIGUOUS_GROUP_004 | B9J08_000329 / B9J08_005587 |

These synthetic groups represent ambiguous transcript-to-locus assignments, not additional biological genes. The duplicate-aware grouped analysis was selected as primary; baseline results were retained as a sensitivity comparison.

## Differential expression
Gene-level differential expression used edgeR with a quasi-likelihood framework and limma-voom. Methods were compared using FDR, log fold changes, significance overlap, direction concordance, and correlation of effect-size estimates. The grouped analysis started with 5,415 groups and retained 5,376 for testing.

## Pathway analysis
fgsea was used for ranked gene-set enrichment and GSVA for pathway activity estimation. Their results were compared for direction, significance, and agreement. Both analyses derive from the same expression dataset and are not independent validation.

## Alternative splicing
rMATS was run on STAR-aligned BAM files. Event classes were SE, MXE, A5SS, A3SS, and RI. Both JC and JCEC outputs were considered; these are alternative counting modes and must not be counted as independent biological events.

## Functional interpretation
Functional annotation used eggNOG-mapper and protein/reference annotation resources. Interpretation included KEGG pathway gene sets, GO over-representation, directional GO analysis, leading-edge integration, and core-gene prioritization.

The core set contained 43 unique genes: 15 associated with energy/carbon metabolism and 28 associated with ribosome biogenesis/RNA-related processes. Internal consistency checks are not external biological validation.

## Statistical interpretation
The reported significance threshold is FDR < 0.05 unless stated otherwise. Effect sizes, sample size, annotation limitations, and method agreement must be considered alongside adjusted p-values. GO and pathway analyses are downstream interpretations of the same expression evidence.
