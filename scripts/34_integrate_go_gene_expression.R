#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
})

cat("\n========== STEP 34: GO-GENE EXPRESSION INTEGRATION ==========\n")

input_go <- "results/go_enrichment/directional_ORA/GO_ORA_directional_significant.csv"
annotation_file <- "results/integrated_analysis/integrated_DE_functional_annotation.csv"
tier1_file <- "results/pathway_enrichment_grouped/biological_interpretation/tier1_50_genes_functional_characterization.csv"

out_dir <- "results/go_enrichment/biological_interpretation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

go <- fread(input_go)
ann <- fread(annotation_file)
tier1 <- fread(tier1_file)

required_go <- c(
    "direction", "ONTOLOGY", "GO", "TERM",
    "foreground_count", "FoldEnrichment", "FDR", "genes"
)

required_ann <- c(
    "gene_id", "gene_name", "Preferred_name", "product",
    "COG_category", "KEGG_Pathway",
    "edgeR_logFC", "edgeR_FDR",
    "limma_logFC", "limma_FDR",
    "consensus_logFC", "consensus_direction"
)

stopifnot(all(required_go %in% names(go)))
stopifnot(all(required_ann %in% names(ann)))

# --------------------------------------------------
# 1. EXPAND GO TERM GENE LISTS
# --------------------------------------------------

go_long <- go[
    ,
    .(gene_id = unlist(strsplit(genes, ";", fixed = TRUE))),
    by = .(
        direction, ONTOLOGY, GO, TERM,
        foreground_count, FoldEnrichment, pvalue, FDR
    )
]

go_long <- unique(go_long)

# --------------------------------------------------
# 2. JOIN FUNCTIONAL ANNOTATION
# --------------------------------------------------

annotation_cols <- c(
    "gene_id", "gene_name", "Preferred_name", "product",
    "COG_category", "KEGG_Pathway",
    "edgeR_logFC", "edgeR_FDR",
    "limma_logFC", "limma_FDR",
    "consensus_logFC", "consensus_direction"
)

ann_small <- unique(
    ann[, ..annotation_cols],
    by = "gene_id"
)

integrated <- merge(
    go_long,
    ann_small,
    by = "gene_id",
    all.x = TRUE
)

# --------------------------------------------------
# 3. VALIDATE DIRECTION
# --------------------------------------------------

integrated[
    ,
    direction_check := fifelse(
        direction == "up_in_resistant" &
            consensus_direction == "up_in_resistant",
        "CONCORDANT",
        fifelse(
            direction == "down_in_resistant" &
                consensus_direction == "down_in_resistant",
            "CONCORDANT",
            "CHECK"
        )
    )
]

# --------------------------------------------------
# 4. ADD TIER 1 PRIORITY INFORMATION
# --------------------------------------------------

tier1_small <- unique(
    tier1[
        ,
        .(
            gene_id,
            rank,
            pathways,
            pathway_count,
            pathway_directions
        )
    ],
    by = "gene_id"
)

integrated <- merge(
    integrated,
    tier1_small,
    by = "gene_id",
    all.x = TRUE
)

setorder(integrated, direction, ONTOLOGY, FDR, GO, gene_id)

# --------------------------------------------------
# 5. SAVE GENE-LEVEL TABLE
# --------------------------------------------------

fwrite(
    integrated,
    file.path(out_dir, "GO_significant_terms_gene_expression_integration.csv")
)

fwrite(
    integrated[direction == "up_in_resistant"],
    file.path(out_dir, "GO_upregulated_genes_integrated.csv")
)

fwrite(
    integrated[direction == "down_in_resistant"],
    file.path(out_dir, "GO_downregulated_genes_integrated.csv")
)

# --------------------------------------------------
# 6. TERM-LEVEL SUMMARY
# --------------------------------------------------

term_summary <- integrated[
    ,
    .(
        genes_in_term = uniqueN(gene_id),
        genes_common_significant_DE = uniqueN(
            gene_id[
                !is.na(edgeR_FDR) &
                !is.na(limma_FDR) &
                edgeR_FDR < 0.05 &
                limma_FDR < 0.05
            ]
        ),
        genes_edgeR_significant = uniqueN(
            gene_id[!is.na(edgeR_FDR) & edgeR_FDR < 0.05]
        ),
        genes_limma_significant = uniqueN(
            gene_id[!is.na(limma_FDR) & limma_FDR < 0.05]
        ),
        concordant_genes = sum(direction_check == "CONCORDANT"),
        genes_with_named_products = uniqueN(
            gene_id[
                !is.na(product) &
                product != "" &
                product != "hypothetical protein"
            ]
        )
    ),
    by = .(
        direction, ONTOLOGY, GO, TERM,
        foreground_count, FoldEnrichment, FDR
    )
]

setorder(term_summary, direction, FDR)

fwrite(
    term_summary,
    file.path(out_dir, "GO_term_gene_expression_summary.csv")
)

# --------------------------------------------------
# 7. VALIDATION AND REPORT
# --------------------------------------------------

cat("\n========== INTEGRATION SUMMARY ==========\n")
cat("Significant GO term-gene records:", nrow(integrated), "\n")
cat("Unique genes:", uniqueN(integrated$gene_id), "\n")
cat("Unique GO terms:", uniqueN(integrated$GO), "\n")

cat("\nDirection consistency:\n")
print(
    integrated[
        ,
        .(
            records = .N,
            concordant = sum(direction_check == "CONCORDANT"),
            check_required = sum(direction_check == "CHECK")
        ),
        by = direction
    ]
)

cat("\nTop terms with gene-level expression support:\n")
print(
    term_summary[
        order(FDR)
    ][
        ,
        .(
            direction, ONTOLOGY, TERM,
            genes_in_term,
            genes_common_significant_DE,
            FoldEnrichment, FDR
        )
    ][1:min(20, .N)]
)

cat("\n========== OUTPUT FILES ==========\n")
cat(out_dir, "\n")

cat("\nPASS: GO-gene expression integration completed.\n")
