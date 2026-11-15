# Reference data and annotations

This directory contains the reference sequences, gene/transcript mappings, annotations, and locally generated indices used by the *Candida auris* bulk RNA-seq analysis.

The reference assets support STAR alignment, Salmon transcript quantification, transcript-to-gene summarization with `tximport`, and downstream functional interpretation.

> **Important:** Large reference assets and generated indices are not intended to be committed to the public Git repository. This README documents their expected role and layout. Obtain source files from their original providers and regenerate software-specific indices locally.

## Reference assets

The current working project contains the following reference-related files. File availability in a public clone may differ because large files and generated indices are excluded from version control.

| File / directory | Description | Use |
|---|---|---|
| `star/candida_auris_genome.fna` | Reference genome FASTA used for STAR alignment | Genome reference |
| `star/candida_auris_genome.fna.fai` | FASTA index for the reference genome | Random access to FASTA |
| `star/candida_auris_mRNA_only.gtf` | mRNA-focused gene annotation used to build the STAR index | STAR splice-junction annotation |
| `star/candida_auris_annotation.gtf` | Broader gene annotation file | Annotation/reference comparison |
| `mrna_only.gff` | mRNA-focused GFF annotation | Transcript/gene feature reference |
| `candida_auris_mrna.fa` | mRNA/transcript FASTA used for the baseline Salmon reference | Transcript quantification |
| `candida_auris_transcripts.fa` | Alternative transcript FASTA retained in the working reference directory | Reference comparison / audit |
| `proteins/candida_auris_proteins.faa` | Protein FASTA | Protein-to-gene mapping and functional annotation |
| `tx2gene.tsv` | Transcript-to-gene mapping for the baseline quantification | `tximport` |
| `tx2gene_keepDuplicates_grouped.tsv` | Transcript-to-gene mapping with four explicit ambiguous duplicate groups | Duplicate-aware `tximport` analysis |
| `star/index/` | Locally generated STAR genome index | STAR alignment |
| `salmon_index_mrna/` | Baseline Salmon index | Baseline Salmon quantification |
| `salmon_index_mrna_keepDuplicates/` | Salmon index generated with duplicate transcript records retained | Duplicate-aware Salmon quantification |
| `eggnog_db/` | eggNOG-mapper database files | Functional annotation |

Some of the files above may be stored in nested subdirectories. The table describes their logical roles; it is not a promise that every asset will be distributed with the repository.

## Reference and annotation characteristics

The recorded genome reference is the *C. auris* B8441 assembly used by the project.

Reference validation performed during the analysis reported:

- 15 genome contigs;
- 12,365,959 total reference bases;
- 15 contigs represented in the mRNA-only GTF;
- 5,419 transcript records and 5,890 exon records in the mRNA-only annotation;
- no annotation contigs absent from the genome FASTA;
- no invalid or out-of-bounds annotation coordinates.

The mRNA-only annotation contains one transcript per gene in the records used for this analysis. This is an important limitation when interpreting transcript-isoform-level alternative-splicing results.

The broader annotation file and the mRNA-only annotation have different record counts. Use the specific annotation file documented for each analysis step rather than substituting one for the other.

## Transcript-to-gene mappings and duplicate handling

Two mapping files are maintained because the project compared alternative handling of duplicate transcript sequences.

### `tx2gene.tsv`

Baseline transcript-to-gene mapping used with the baseline Salmon index.

- 5,419 transcript records;
- 5,419 unique transcript identifiers;
- 5,419 unique gene identifiers;
- no duplicated transcript identifiers;
- no synthetic ambiguous groups.

### `tx2gene_keepDuplicates_grouped.tsv`

Mapping used with the Salmon index built using `--keepDuplicates`.

- 5,419 transcript records;
- 5,419 unique transcript identifiers;
- 5,415 gene groups;
- four synthetic groups, each combining a pair of exact duplicate transcript sequences.

The synthetic groups are:

| Synthetic group | Transcript identifiers |
|---|---|
| `AMBIGUOUS_GROUP_001` | `B9J08_004680`, `B9J08_005536` |
| `AMBIGUOUS_GROUP_002` | `B9J08_004683`, `B9J08_005540` |
| `AMBIGUOUS_GROUP_003` | `B9J08_004684`, `B9J08_005541` |
| `AMBIGUOUS_GROUP_004` | `B9J08_000329`, `B9J08_005587` |

These synthetic identifiers represent ambiguous transcript-to-locus assignments for the purpose of grouped quantification. They are **not** additional biological genes and should not be interpreted as such.

The baseline and duplicate-aware analyses were compared as a sensitivity analysis. For the 5,372 non-ambiguous biological genes, the differential-expression significance classification and direction were unchanged between the two approaches.

## Software-specific indices

The following indices were generated locally and should be rebuilt for the software versions and reference files used in a new environment:

### STAR

The STAR index is stored locally under `star/index/`.

Recorded index-build parameters include:

- genome FASTA: `candida_auris_genome.fna`;
- annotation: `candida_auris_mRNA_only.gtf`;
- `sjdbOverhang`: 149;
- `genomeSAindexNbases`: 10;
- build-time STAR version: 2.7.4a.

The STAR executable recorded for the subsequent analysis environment was version 2.7.11b. The index-build version and executable version are therefore intentionally documented separately.

### Salmon

Two Salmon indices were generated with Salmon 2.8.0 and k-mer length 21:

- `salmon_index_mrna/`: baseline index; 5,415 references; duplicate sequences collapsed;
- `salmon_index_mrna_keepDuplicates/`: duplicate-aware index; 5,419 references; duplicate sequences retained.

The two indices were used for separate quantification branches. Do not interchange their associated quantification outputs or transcript-to-gene mapping files.

## Rebuilding and obtaining reference assets

A fresh clone may not contain the genome FASTA, transcript/protein FASTA files, annotation files, eggNOG database, or STAR/Salmon indices.

Before reproducing the analysis:

1. Obtain the *C. auris* B8441 genome assembly and matching annotation from the original public reference provider.
2. Verify the assembly accession, annotation release, file checksums, and provider citation against the original download records.
3. Place the files in the corresponding paths described in this README.
4. Confirm that FASTA sequence identifiers match the annotation identifiers.
5. Rebuild the STAR and Salmon indices using the documented input files and parameters.
6. Keep the baseline and duplicate-aware Salmon branches separate.
7. Record software versions and checksums for downloaded reference files in the project documentation.

**Provenance fields to complete before public release:** the exact source URL, assembly accession/version, annotation release/date, download date, and checksums for each reference FASTA/GFF/GTF. These details are not asserted here because they must be checked against the original acquisition records.

## Data and licensing

Reference genome, annotation, protein, and eggNOG database files may have their own source-specific terms of use and citation requirements. Cite the original providers and publications when using these resources. Do not assume that the project's software license, if one is later selected, changes the terms attached to third-party reference data.

Large databases and generated indices should remain outside Git version control. The public repository should contain this documentation and reproducible instructions rather than multi-gigabyte reference databases or index directories.
