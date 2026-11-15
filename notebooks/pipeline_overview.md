# Pipeline overview

## Overview
The project uses R and Python scripts with command-line bioinformatics tools. Scripts are stored in `scripts/` and should be run from the project root unless otherwise documented.

The workflow contains connected branches for read QC/alignment/quantification, baseline expression analysis, duplicate-aware sensitivity analysis, pathway analysis, alternative-splicing analysis, and GO/core-gene interpretation.

This is a script-based workflow. Although Nextflow is installed in the environment, the audit found no `main.nf` or `nextflow.config`; this repository does not claim to provide an implemented Nextflow pipeline.

## Samples
- Sensitive: SRR10059687, SRR10059688, SRR10059689
- Resistant: SRR10059690, SRR10059691, SRR10059692

Metadata: `data/metadata/sample_metadata.csv`.

## Script groups

### Baseline expression
| Script | Role |
|---|---|
| `02_tximport.R` | Import baseline Salmon quantifications |
| `03_sample_qc.R` | Sample-level QC |
| `04_tximport_length_scaled.R` | Create length-scaled tximport object |
| `05_edgeR_DE.R` | Baseline edgeR analysis |
| `06_edgeR_diagnostics.R` | edgeR diagnostics |
| `07_limma_voom_DE.R` | Baseline limma-voom |
| `08_compare_DE_methods.R` | Compare DE methods |

### Annotation and baseline pathways
| Script | Role |
|---|---|
| `09_build_gene_annotation.py` | Build gene annotation table |
| `10_analyze_eggnog_annotations.py` | Process eggNOG annotations |
| `10_validate_protein_gene_mapping.py` | Validate protein-to-gene mapping |
| `11_integrate_DE_annotation.py` | Integrate DE and annotation |
| `12_functional_characterization.py` | Functional characterization |
| `13_fgsea_pathway_enrichment.R` | Baseline fgsea |
| `14_gsva_pathway_analysis.R` | Baseline GSVA |
| `15_compare_fgsea_gsva.R` | Compare pathway methods |

### Alternative splicing
| Script | Role |
|---|---|
| `16_integrate_splicing_expression.py` | Integrate rMATS and expression |
| `17_validate_se_bam.py` | BAM-based SE investigation |
| `18_plot_se_validation.R` | Plot SE validation |
| `19_splicing_event_deduplication.py` | Summarize unique biological events |

rMATS execution used the separate `p2_rmats` environment; `rmats.py` was not on the `p2_ngs` PATH.

### Duplicate-aware grouped analysis
| Script | Role |
|---|---|
| `20_tximport_keepDuplicates_grouped.R` | Import duplicate-aware Salmon quantifications |
| `21_tximport_keepDuplicates_length_scaled.R` | Length-scaled grouped tximport |
| `22_edgeR_DE_grouped.R` | Grouped edgeR |
| `23_limma_voom_DE_grouped.R` | Grouped limma-voom |
| `24_compare_DE_methods_grouped.R` | Compare grouped DE methods |
| `25_fgsea_grouped.R` | Grouped fgsea |
| `26_gsva_grouped.R` | Grouped GSVA |
| `27_compare_fgsea_gsva_grouped.R` | Compare grouped pathway methods |
| `28_integrate_shared_pathway_genes.R` | Integrate leading-edge genes |
| `29_rank_leading_edge_genes.R` | Rank leading-edge genes |

### GO and core genes
| Script | Role |
|---|---|
| `31_go_enrichment_ora.R` | GO over-representation |
| `32_go_directional_ora.R` | Directional GO |
| `33_plot_go_directional.R` | Plot directional GO |
| `34_integrate_go_gene_expression.R` | Integrate GO and expression |
| `35_identify_core_go_genes.R` | Identify core genes |
| `36_audit_core_gene_annotations.R` | Audit annotations |
| `37_merge_core_gene_eggnog_annotations.py` | Refine eggNOG annotations |

There is no script numbered 30 in the audited list.

## Important dependencies
- Salmon `quant.sf` + `tx2gene.tsv` -> baseline tximport.
- Baseline tximport + metadata -> baseline edgeR and limma-voom.
- Duplicate-aware `quant.sf` + `tx2gene_keepDuplicates_grouped.tsv` -> grouped tximport.
- Grouped tximport + metadata -> grouped edgeR and limma-voom.
- Grouped DE + annotation -> fgsea/GSVA and leading-edge integration.
- rMATS outputs + STAR BAMs + annotation/expression -> splicing integration and validation.
- Tier 1 pathway-leading-edge genes -> GO scripts 31–34 -> script 35 core-gene selection -> annotation audits.

## Distinct gene sets
Scripts 31–34 use the Tier 1 50-gene functional-characterization input. Script 35 subsequently selects a 43-gene core set from selected GO terms. These are different sets and must not be conflated.

## Portability
Several scripts (including scripts 02–10 and 20–24) contain machine-specific absolute paths. These must be converted to portable project-root discovery/configuration before claiming that a fresh clone can run the complete workflow. Check required inputs and outputs at each step.
