# Pathway analysis

## Objective and methods
Pathway analysis interpreted gene-level expression patterns at pathway level. Two approaches were used:
- fgsea 1.36.2 for ranked gene-set enrichment;
- GSVA 2.4.9 for pathway activity estimation.

These analyses use the same expression dataset as differential-expression testing. They are complementary interpretations, not independent experimental validation.

## Grouped-analysis results
| Metric | fgsea | GSVA |
|---|---:|---:|
| Pathways tested | 228 | 228 |
| Significant (FDR < 0.05) | 13 | 6 |

Five pathways were significant in both methods; eight were fgsea-only and one was GSVA-only. Across 228 common pathways, directional agreement was 192/228 (84.21%), and Spearman correlation of pathway effect estimates was 0.8554112.

### Pathways significant in both methods
| KEGG ID | fgsea NES | fgsea FDR | GSVA logFC | GSVA FDR | Higher in |
|---|---:|---:|---:|---:|---|
| ko03008 | -2.231802 | 4.460259e-06 | -0.6797348 | 0.04473045 | Sensitive |
| ko03020 | -2.045073 | 0.004086842 | -0.6829869 | 0.04473045 | Sensitive |
| ko04066 | 1.931130 | 0.01667009 | 0.7555071 | 0.04473045 | Resistant |
| ko00710 | 1.931308 | 0.01820716 | 0.8632716 | 0.03706494 | Resistant |
| ko05230 | 1.921070 | 0.02175015 | 0.7488413 | 0.04473045 | Resistant |

Positive direction indicates higher in resistant for the reported contrast; negative indicates higher in sensitive.

## Leading-edge prioritization
The fgsea leading-edge integration contained 99 gene-pathway records representing 91 unique genes. Tier 1 contained 50 genes significant in both differential-expression methods with concordant direction.

GO scripts 31–34 use the Tier 1 50-gene functional-characterization input. Script 35 later derives a distinct 43-gene core set from selected GO terms.

## Interpretation and limitations
- fgsea and GSVA address related but different questions.
- Pathway gene sets overlap and may be hierarchical.
- A significant pathway does not mean every member gene is activated.
- Results are exploratory with three samples per group.
- GO and pathway results derive from the same expression evidence and are not independent replication.

## Outputs and scripts
Grouped results: `results/pathway_enrichment_grouped/fgsea/`, `gsva/`, and `method_comparison/`.

Scripts: `25_fgsea_grouped.R`, `26_gsva_grouped.R`, `27_compare_fgsea_gsva_grouped.R`, `28_integrate_shared_pathway_genes.R`, and `29_rank_leading_edge_genes.R`.
