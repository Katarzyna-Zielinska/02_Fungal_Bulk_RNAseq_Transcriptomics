# Alternative-splicing analysis

## Objective
Identify candidate alternative-splicing differences between sensitive and resistant *Candida auris* samples and integrate them with gene-expression results.

## Method
rMATS 4.4.0 was run using STAR-aligned BAM files. The five event classes were:
- SE: skipped exon;
- MXE: mutually exclusive exons;
- A5SS: alternative 5-prime splice site;
- A3SS: alternative 3-prime splice site;
- RI: retained intron.

Both JC (junction counts) and JCEC (junction plus exon counts) outputs were considered. They are alternative counting modes and must not be counted as independent biological events.

## Results
| Event type | Events reported | Significant (FDR < 0.05) |
|---|---:|---:|
| SE | 11 | 1 |
| MXE | 2 | 2 |
| A3SS | 0 | 0 |
| A5SS | 0 | 0 |
| RI | 0 | 0 |

There were three unique significant biological events across two genes. The integrated table contains six records because both JC and JCEC records are represented.

### Significant SE event
The significant SE event was associated with `gene-B9J08_004599` (event ID 8). FDR was 0.0001906325 for JC and 0.0001920031 for JCEC. Reported delta PSI was approximately -0.032. The gene was also upregulated in the resistant group in both differential-expression methods.

The delta PSI is modest. BAM-derived PSI values in this project are approximate proxies, not formal replacements for rMATS PSI.

### Significant MXE events
Two MXE events were associated with `gene-B9J08_000840`. JC/JCEC FDR values were approximately 0.0400366 and 0.0065021. The gene was also upregulated in resistant in both differential-expression methods. Junction support was weak or uneven, and the two paths were not equally supported. These remain exploratory candidate events, not a robustly confirmed splice switch.

## Annotation limitation
The mRNA-only annotation contains 5,419 transcript records and 5,419 distinct gene IDs (one annotated transcript per gene). This limits transcript-isoform-level interpretation. Event calls do not prove that all isoforms are represented in the annotation.

## BAM-based investigation
For the strict paired-read audit of the selected locus in SRR10059691:
- 624 query names were inspected;
- 481 met the defined complete reciprocal-pair criteria in both loci;
- 143 had no valid complete pair under those criteria;
- no valid complete pair was assigned exclusively to only one of the two duplicate loci.

This is an interval-based audit, not a genome-wide proof that all reads are ambiguous or uniquely assignable.

## Interpretation
The findings are computational candidates. Small sample size, one-transcript-per-gene annotation, modest SE delta PSI, and uneven MXE junction support limit strong biological claims. Experimental validation and additional datasets are required.

## Outputs and scripts
Integrated outputs include `results/integrated_analysis/splicing/splicing_expression_integrated.csv`, `splicing_gene_level_summary.csv`, and `unique_splicing_events_summary.csv`.

Scripts: `16_integrate_splicing_expression.py`, `17_validate_se_bam.py`, `18_plot_se_validation.R`, and `19_splicing_event_deduplication.py`. rMATS execution used the separate `p2_rmats` environment.
