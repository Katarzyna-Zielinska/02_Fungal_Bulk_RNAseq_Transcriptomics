#!/usr/bin/env Rscript

# ============================================================
# Project 2: Candida auris Bulk RNA-seq Transcriptomics
# Step 29: Ranking of shared pathway leading-edge genes
# ============================================================

suppressPackageStartupMessages({
    library(data.table)
})

input_file <- paste0(
    "results/pathway_enrichment_grouped/biological_interpretation/",
    "shared_pathways_leading_edge_gene_integration.csv"
)

output_dir <- "results/pathway_enrichment_grouped/biological_interpretation"

if (!file.exists(input_file)) {
    stop("Missing input file: ", input_file)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n=== LOADING INTEGRATED GENE DATA ===\n")

x <- fread(input_file)

required <- c(
    "gene_id", "pathway", "NES", "fgsea_direction",
    "gene_name", "product", "Description", "Preferred_name",
    "COG_category", "GOs", "has_functional_annotation",
    "edgeR_logFC", "edgeR_FDR",
    "limma_logFC", "limma_FDR",
    "DE_class", "consensus_direction"
)

missing_columns <- setdiff(required, names(x))

if (length(missing_columns) > 0) {
    stop(
        "Missing columns: ",
        paste(missing_columns, collapse = ", ")
    )
}

cat("Gene-pathway records:", nrow(x), "\n")
cat("Unique genes:", uniqueN(x$gene_id), "\n")

# ------------------------------------------------------------
# 1. Evaluate agreement between pathway and gene direction
# ------------------------------------------------------------

cat("\n=== EVALUATING DIRECTIONAL AGREEMENT ===\n")

x[, pathway_direction_concordant := (
    (fgsea_direction == "Resistant" &
     consensus_direction == "up_in_resistant") |
    (fgsea_direction == "Sensitive" &
     consensus_direction == "down_in_resistant")
)]

# ------------------------------------------------------------
# 2. Create one record per gene
# ------------------------------------------------------------

cat("\n=== AGGREGATING GENES ===\n")

gene_info <- x[
    ,
    .(
        gene_name = first(na.omit(gene_name)),
        product = first(na.omit(product)),
        Description = first(na.omit(Description)),
        Preferred_name = first(na.omit(Preferred_name)),
        COG_category = first(na.omit(COG_category)),
        GOs = first(na.omit(GOs)),
        has_functional_annotation = any(
            has_functional_annotation %in% TRUE,
            na.rm = TRUE
        ),
        edgeR_logFC = first(edgeR_logFC),
        edgeR_FDR = first(edgeR_FDR),
        limma_logFC = first(limma_logFC),
        limma_FDR = first(limma_FDR),
        DE_class = first(DE_class),
        consensus_direction = first(consensus_direction),
        consensus_logFC = mean(
            c(first(edgeR_logFC), first(limma_logFC)),
            na.rm = TRUE
        ),
        min_DE_FDR = min(
            c(first(edgeR_FDR), first(limma_FDR)),
            na.rm = TRUE
        ),
        pathway_count = uniqueN(pathway),
        pathways = paste(sort(unique(pathway)), collapse = ";"),
        pathway_directions = paste(
            sort(unique(fgsea_direction)),
            collapse = ";"
        ),
        all_pathway_directions_match = all(
            pathway_direction_concordant
        ),
        pathway_direction_matches = sum(
            pathway_direction_concordant,
            na.rm = TRUE
        )
    ),
    by = gene_id
]

# ------------------------------------------------------------
# 3. Define transparent priority classes
# ------------------------------------------------------------

cat("\n=== ASSIGNING PRIORITY CLASSES ===\n")

gene_info[, priority_class := fcase(
    DE_class == "common_significant" &
        all_pathway_directions_match,
    "Tier 1 - significant in both DE methods; pathway direction concordant",

    DE_class == "common_significant" &
        !all_pathway_directions_match,
    "Tier 2 - significant in both DE methods; direction requires review",

    DE_class %in% c("edgeR_only", "limma_only"),
    "Tier 3 - significant in one DE method",

    default = "Tier 4 - not significant in either DE method"
)]

gene_info[, priority_order := fcase(
    grepl("^Tier 1", priority_class), 1L,
    grepl("^Tier 2", priority_class), 2L,
    grepl("^Tier 3", priority_class), 3L,
    default = 4L
)]

# ------------------------------------------------------------
# 4. Rank genes
# ------------------------------------------------------------

cat("\n=== RANKING GENES ===\n")

gene_info[, abs_consensus_logFC := abs(consensus_logFC)]

setorder(
    gene_info,
    priority_order,
    -pathway_count,
    min_DE_FDR,
    -abs_consensus_logFC
)

gene_info[, abs_consensus_logFC := NULL]

gene_info[, rank := seq_len(.N)]

setcolorder(
    gene_info,
    c("rank", "gene_id", "priority_class", "pathway_count",
      "pathways", "pathway_directions",
      "pathway_direction_matches", "all_pathway_directions_match",
      "gene_name", "Preferred_name", "product", "Description",
      "COG_category", "GOs", "has_functional_annotation",
      "edgeR_logFC", "edgeR_FDR",
      "limma_logFC", "limma_FDR",
      "consensus_logFC", "consensus_direction",
      "min_DE_FDR", "DE_class", "priority_order")
)

# ------------------------------------------------------------
# 5. Export complete gene ranking
# ------------------------------------------------------------

ranking_file <- file.path(
    output_dir,
    "leading_edge_gene_priority_ranking.csv"
)

fwrite(gene_info, ranking_file)

# ------------------------------------------------------------
# 6. Export high-priority genes
# ------------------------------------------------------------

high_priority <- gene_info[
    priority_order == 1
]

high_priority_file <- file.path(
    output_dir,
    "leading_edge_high_priority_genes.csv"
)

fwrite(high_priority, high_priority_file)

# ------------------------------------------------------------
# 7. Export pathway-level gene ranking
# ------------------------------------------------------------

pathway_ranking <- copy(x)

pathway_ranking[, pathway_direction_concordant := as.logical(
    pathway_direction_concordant
)]

setorder(
    pathway_ranking,
    pathway,
    DE_class,
    edgeR_FDR,
    limma_FDR
)

pathway_ranking_file <- file.path(
    output_dir,
    "leading_edge_gene_ranking_by_pathway.csv"
)

fwrite(pathway_ranking, pathway_ranking_file)

# ------------------------------------------------------------
# 8. Priority summary
# ------------------------------------------------------------

priority_summary <- gene_info[
    ,
    .(
        genes = .N,
        genes_with_annotation = sum(
            has_functional_annotation,
            na.rm = TRUE
        ),
        genes_in_multiple_pathways = sum(pathway_count > 1)
    ),
    by = priority_class
]

summary_file <- file.path(
    output_dir,
    "leading_edge_gene_priority_summary.csv"
)

fwrite(priority_summary, summary_file)

# ------------------------------------------------------------
# 9. Validation
# ------------------------------------------------------------

cat("\n=== FINAL VALIDATION ===\n")

stopifnot(
    nrow(gene_info) == uniqueN(x$gene_id),
    !anyDuplicated(gene_info$gene_id),
    all(!is.na(gene_info$rank)),
    all(gene_info$rank == seq_len(nrow(gene_info))),
    all(gene_info$priority_order %in% 1:4),
    all(high_priority$DE_class == "common_significant"),
    all(high_priority$all_pathway_directions_match)
)

cat("Unique genes ranked:", nrow(gene_info), "\n")
cat("Tier 1 genes:", nrow(high_priority), "\n")
cat("Genes in multiple pathways:", sum(gene_info$pathway_count > 1), "\n")
cat("Genes with functional annotation:", sum(
    gene_info$has_functional_annotation
), "\n")

cat("\n=== PRIORITY SUMMARY ===\n")
print(priority_summary)

cat("\n=== TOP 20 GENES ===\n")

print(
    gene_info[
        1:min(20, .N),
        .(
            rank,
            gene_id,
            Preferred_name,
            product,
            pathway_count,
            edgeR_logFC,
            edgeR_FDR,
            limma_logFC,
            limma_FDR,
            consensus_direction,
            priority_class
        )
    ]
)

cat("\n========================================\n")
cat("STEP 29 COMPLETED SUCCESSFULLY\n")
cat("========================================\n")

cat("\nGenerated files:\n")
cat(ranking_file, "\n")
cat(high_priority_file, "\n")
cat(pathway_ranking_file, "\n")
cat(summary_file, "\n")
