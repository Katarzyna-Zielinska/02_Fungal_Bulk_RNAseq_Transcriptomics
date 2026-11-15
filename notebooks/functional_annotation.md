# Functional annotation and GO interpretation

## Objective
Functional annotation adds biological context to differentially expressed genes and pathway-leading-edge genes. The workflow combines protein/gene annotation, KEGG pathway interpretation, Gene Ontology analysis, and candidate-gene prioritization.

## Resources and processing
The project used eggNOG-mapper annotations, a *C. auris* protein FASTA, gene/transcript annotation, KEGG pathway gene sets, and GO resources accessed through Bioconductor packages including GO.db and AnnotationDbi.

The exact eggNOG database release and citation must be recorded from the original acquisition records before publication.

## Leading-edge genes
The grouped fgsea leading-edge integration contained 99 gene-pathway records representing 91 unique genes. Tier 1 contained 50 genes significant in both edgeR and limma-voom with concordant direction.

GO scripts 31–34 use the Tier 1 50-gene functional-characterization input. Script 35 subsequently selects a separate 43-gene core set.

## GO over-representation
Recorded settings:
- background: 4,007 genes with valid GO annotation;
- foreground: 50 Tier 1 genes, 49 with GO annotation;
- `gene-B9J08_003082` had no GO annotation;
- GO.db 3.22.0;
- Benjamini-Hochberg adjustment applied separately by ontology;
- minimum gene-set size 5; maximum 500.

| Ontology | Tested terms | Significant terms |
|---|---:|---:|
| BP | 618 | 144 |
| MF | 111 | 13 |
| CC | 86 | 24 |

Examples:
- CC nucleolus: 27/258, fold enrichment 8.56, FDR 2.20e-18;
- BP rRNA metabolic process: 26/271, fold enrichment 7.85, FDR 1.06e-15.

GO enrichment indicates over-representation in the selected gene set. It does not by itself establish pathway activation or causality.

## Directional GO
The Tier 1 set contained 18 genes upregulated in resistant and 32 downregulated. All 18 upregulated genes and 31 downregulated genes had GO annotations. Separate enrichment results and dot plots were generated for up- and downregulated genes.

Because the foreground was selected from expression/pathway results, directional GO is an interpretation of the same evidence, not independent validation.

## Core-gene prioritization
Script 35 selected 16 representative significant GO terms and identified 43 unique genes:

| Group | Genes |
|---|---:|
| Energy/carbon metabolism | 15 |
| Ribosome biogenesis/RNA-related | 28 |
| **Total** | **43** |

The groups did not overlap; their union exactly matched the 43-gene set. All 43 genes were significant in both edgeR and limma-voom with concordant directions:
- 15/15 energy/carbon-metabolism genes were up in resistant;
- 28/28 ribosome/RNA-related genes were down in resistant.

This is an internal consistency check, not external validation.

## Annotation audit
The initial audit found named functional product annotations for 25 genes and hypothetical/uncharacterized product descriptions for 18 genes. The 18 genes requiring further annotation were uniquely mapped to the protein FASTA and eggNOG records:
- 18/18 had PFAM evidence;
- 18/18 had GO annotation;
- 17/18 had a preferred name.

Recovered evidence included families/domains associated with fructose-bisphosphatase, hexokinase, WD40 proteins, U3 snoRNA-associated proteins, nucleolar pre-rRNA processing, GTPases, RNase/PIN domains, and RNA polymerase subunits. Interpret these descriptions only at the level supported by the source annotation.

## Outputs and scripts
Main output areas: `results/pathway_enrichment_grouped/`, `results/go_enrichment/`, and `results/integrated_analysis/`.

Scripts: `31_go_enrichment_ora.R`, `32_go_directional_ora.R`, `33_plot_go_directional.R`, `34_integrate_go_gene_expression.R`, `35_identify_core_go_genes.R`, `36_audit_core_gene_annotations.R`, and `37_merge_core_gene_eggnog_annotations.py`.
