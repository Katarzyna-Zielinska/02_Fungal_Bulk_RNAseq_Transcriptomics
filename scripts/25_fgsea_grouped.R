#!/usr/bin/env Rscript

# ============================================================
# Project 2: Candida auris Bulk RNA-seq Transcriptomics
# Step 25: KEGG pathway enrichment using fgsea
# Analysis: Grouped transcript/gene mapping
# Contrast: Resistant vs Sensitive
# ============================================================


# ============================================================
# 1. CONFIGURATION
# ============================================================

project_dir <- normalizePath(".", mustWork = TRUE)

edger_file <- file.path(
    project_dir,
    "results/differential_expression_grouped/edgeR_all_genes_grouped.csv"
)

integrated_file <- file.path(
    project_dir,
    "results/integrated_analysis/integrated_DE_functional_annotation.csv"
)

output_dir <- file.path(
    project_dir,
    "results/pathway_enrichment_grouped/fgsea"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

FDR_THRESHOLD <- 0.05
MIN_SIZE <- 10
MAX_SIZE <- 500

EXPECTED_GENES <- 5376

# Original IDs belonging to four identical-sequence pairs.
# They are represented by synthetic groups in grouped analysis
# and must not be assigned independently to KEGG pathways.

duplicate_locus_ids <- c(
    "gene-B9J08_004680",
    "gene-B9J08_005536",
    "gene-B9J08_004683",
    "gene-B9J08_005540",
    "gene-B9J08_004684",
    "gene-B9J08_005541",
    "gene-B9J08_000329",
    "gene-B9J08_005587"
)


# ============================================================
# 2. LOAD PACKAGES
# ============================================================

cat("\n=== LOADING PACKAGES ===\n")

required_packages <- c("fgsea", "data.table")

for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Required package is not installed: ", pkg)
    }

    cat(
        pkg,
        ":",
        as.character(packageVersion(pkg)),
        "\n"
    )
}

suppressPackageStartupMessages({
    library(fgsea)
    library(data.table)
})


# ============================================================
# 3. CHECK INPUT FILES
# ============================================================

cat("\n=== CHECKING INPUT FILES ===\n")

if (!file.exists(edger_file)) {
    stop("Grouped edgeR input file not found: ", edger_file)
}

if (!file.exists(integrated_file)) {
    stop("Integrated annotation file not found: ", integrated_file)
}

cat("[FOUND]", edger_file, "\n")
cat("[FOUND]", integrated_file, "\n")


# ============================================================
# 4. LOAD DATA
# ============================================================

cat("\n=== LOADING DATA ===\n")

edger <- fread(
    edger_file,
    data.table = FALSE
)

integrated <- fread(
    integrated_file,
    data.table = FALSE
)

cat("Grouped edgeR genes:", nrow(edger), "\n")
cat("Integrated annotation rows:", nrow(integrated), "\n")


# ============================================================
# 5. VALIDATE INPUT DATA
# ============================================================

cat("\n=== VALIDATING DATA ===\n")

required_edger <- c(
    "gene_id",
    "logFC",
    "F",
    "PValue",
    "FDR"
)

required_integrated <- c(
    "gene_id",
    "KEGG_Pathway"
)

if (!all(required_edger %in% colnames(edger))) {
    stop("Required columns missing in grouped edgeR results.")
}

if (!all(required_integrated %in% colnames(integrated))) {
    stop("Required columns missing in integrated annotation.")
}

if (anyDuplicated(edger$gene_id) > 0) {
    stop("Duplicated gene IDs in grouped edgeR results.")
}

if (anyDuplicated(integrated$gene_id) > 0) {
    stop("Duplicated gene IDs in integrated annotation.")
}

if (anyNA(edger$gene_id) || anyNA(integrated$gene_id)) {
    stop("Missing gene IDs detected.")
}

if (anyNA(edger$logFC) || anyNA(edger$F)) {
    stop("Missing logFC or F statistic.")
}

if (
    any(!is.finite(edger$logFC)) ||
    any(!is.finite(edger$F))
) {
    stop("Non-finite edgeR statistics detected.")
}

if (any(edger$F < 0)) {
    stop("Negative F statistics detected.")
}

if (nrow(edger) != EXPECTED_GENES) {
    stop(
        "Expected ",
        EXPECTED_GENES,
        " grouped genes, found ",
        nrow(edger)
    )
}

if (any(edger$gene_id %in% duplicate_locus_ids)) {
    stop("Original duplicate locus IDs found in grouped edgeR results.")
}

cat("Grouped edgeR input validation passed.\n")


# ============================================================
# 6. CREATE SIGNED GENE RANKING
# ============================================================

cat("\n=== CREATING SIGNED GENE RANKING ===\n")

edger$ranking_statistic <- (
    sign(edger$logFC) * sqrt(edger$F)
)

ranking <- edger[
    ,
    c(
        "gene_id",
        "logFC",
        "F",
        "PValue",
        "FDR",
        "ranking_statistic"
    )
]

ranking <- ranking[
    order(
        ranking$ranking_statistic,
        decreasing = TRUE
    ),
]

stats <- ranking$ranking_statistic

names(stats) <- ranking$gene_id

stats <- sort(
    stats,
    decreasing = TRUE
)

if (length(stats) != EXPECTED_GENES) {
    stop("Unexpected number of genes in ranking.")
}

if (anyDuplicated(names(stats)) > 0) {
    stop("Duplicate gene IDs in ranking.")
}

if (any(!is.finite(stats))) {
    stop("Non-finite ranking statistics.")
}

tie_count <- sum(duplicated(stats))

cat("Genes in ranking:", length(stats), "\n")
cat("Positive ranking values:", sum(stats > 0), "\n")
cat("Negative ranking values:", sum(stats < 0), "\n")
cat("Zero ranking values:", sum(stats == 0), "\n")
cat("Tied ranking values:", tie_count, "\n")


# ============================================================
# 7. EXTRACT KEGG PATHWAY ASSIGNMENTS
# ============================================================

cat("\n=== EXTRACTING KEGG PATHWAYS ===\n")

# Remove the eight original locus IDs from the annotation.
# Synthetic ambiguous groups have no KEGG assignments.
# No pathway membership is inferred for those groups.

integrated_grouped <- integrated[
    !integrated$gene_id %in% duplicate_locus_ids,
    ,
    drop = FALSE
]

pathway_rows <- list()

for (i in seq_len(nrow(integrated_grouped))) {

    gene_id <- integrated_grouped$gene_id[i]
    pathway_value <- integrated_grouped$KEGG_Pathway[i]

    if (
        is.na(pathway_value) ||
        pathway_value == ""
    ) {
        next
    }

    pathway_ids <- unlist(
        strsplit(
            as.character(pathway_value),
            ",",
            fixed = TRUE
        )
    )

    pathway_ids <- trimws(pathway_ids)

    # Keep koXXXXX identifiers only.
    # Exclude mapXXXXX identifiers.

    pathway_ids <- pathway_ids[
        grepl(
            "^ko[0-9]{5}$",
            pathway_ids
        )
    ]

    pathway_ids <- unique(pathway_ids)

    if (length(pathway_ids) == 0) {
        next
    }

    pathway_rows[[length(pathway_rows) + 1]] <- data.frame(
        gene_id = gene_id,
        pathway_id = pathway_ids,
        stringsAsFactors = FALSE
    )
}

if (length(pathway_rows) == 0) {
    stop("No valid KEGG pathway assignments found.")
}

pathway_membership_all <- unique(
    do.call(rbind, pathway_rows)
)

cat(
    "All gene-pathway assignments after duplicate-locus exclusion:",
    nrow(pathway_membership_all),
    "\n"
)

cat(
    "Unique KEGG pathway IDs:",
    length(unique(pathway_membership_all$pathway_id)),
    "\n"
)

fwrite(
    pathway_membership_all,
    file.path(
        output_dir,
        "KEGG_ko_gene_pathway_membership_all_genes_grouped.csv"
    )
)


# ============================================================
# 8. FILTER MEMBERSHIP TO TESTED GENES
# ============================================================

cat("\n=== FILTERING TO GROUPED TESTED GENES ===\n")

tested_gene_ids <- names(stats)

pathway_membership <- pathway_membership_all[
    pathway_membership_all$gene_id %in% tested_gene_ids,
    ,
    drop = FALSE
]

pathway_membership <- unique(pathway_membership)

excluded_assignments <- (
    nrow(pathway_membership_all) -
    nrow(pathway_membership)
)

cat(
    "Assignments retained:",
    nrow(pathway_membership),
    "\n"
)

cat(
    "Assignments excluded because genes were not tested:",
    excluded_assignments,
    "\n"
)

if (nrow(pathway_membership) == 0) {
    stop("No KEGG assignments remain after filtering.")
}

fwrite(
    pathway_membership,
    file.path(
        output_dir,
        "KEGG_ko_gene_pathway_membership_tested_genes_grouped.csv"
    )
)


# ============================================================
# 9. BUILD PATHWAY GENE SETS
# ============================================================

cat("\n=== BUILDING PATHWAY GENE SETS ===\n")

pathways <- split(
    pathway_membership$gene_id,
    pathway_membership$pathway_id
)

pathways <- lapply(
    pathways,
    unique
)

pathway_sizes <- data.frame(
    pathway_id = names(pathways),
    gene_count = lengths(pathways),
    stringsAsFactors = FALSE
)

pathway_sizes <- pathway_sizes[
    order(
        pathway_sizes$gene_count,
        decreasing = TRUE
    ),
]

eligible <- (
    lengths(pathways) >= MIN_SIZE &
    lengths(pathways) <= MAX_SIZE
)

pathways_tested <- pathways[eligible]

cat("All KEGG pathway sets:", length(pathways), "\n")
cat("Pathways meeting size criteria:", length(pathways_tested), "\n")
cat("Minimum pathway size:", MIN_SIZE, "\n")
cat("Maximum pathway size:", MAX_SIZE, "\n")

if (length(pathways_tested) == 0) {
    stop("No pathways meet the selected size criteria.")
}

fwrite(
    pathway_sizes,
    file.path(
        output_dir,
        "KEGG_pathway_gene_set_sizes_grouped.csv"
    )
)


# ============================================================
# 10. VALIDATE PATHWAY GENE SETS
# ============================================================

cat("\n=== VALIDATING PATHWAY GENE SETS ===\n")

if (!all(pathway_membership$gene_id %in% tested_gene_ids)) {
    stop("Unranked genes remain in filtered pathway membership.")
}

if (anyDuplicated(pathway_membership) > 0) {
    stop("Duplicated gene-pathway assignments detected.")
}

for (pathway_id in names(pathways_tested)) {

    genes <- pathways_tested[[pathway_id]]

    if (anyDuplicated(genes) > 0) {
        stop("Duplicated genes in pathway: ", pathway_id)
    }

    if (!all(genes %in% names(stats))) {
        stop("Unranked genes found in pathway: ", pathway_id)
    }
}

cat("Pathway gene sets validated.\n")


# ============================================================
# 11. RUN FGSEA MULTILEVEL
# ============================================================

cat("\n=== RUNNING FGSEA MULTILEVEL ===\n")

set.seed(123)

fgsea_results <- fgseaMultilevel(
    pathways = pathways_tested,
    stats = stats,
    minSize = MIN_SIZE,
    maxSize = MAX_SIZE,
    eps = 0
)

fgsea_results <- as.data.frame(fgsea_results)

if (nrow(fgsea_results) == 0) {
    stop("fgsea returned no results.")
}

fgsea_results <- fgsea_results[
    order(
        fgsea_results$padj,
        fgsea_results$pval
    ),
]

significant_results <- fgsea_results[
    !is.na(fgsea_results$padj) &
    fgsea_results$padj < FDR_THRESHOLD,
]

positive_results <- significant_results[
    significant_results$NES > 0,
]

negative_results <- significant_results[
    significant_results$NES < 0,
]

cat("Pathways tested:", nrow(fgsea_results), "\n")
cat(
    "Pathways with adjusted P-value < 0.05:",
    nrow(significant_results),
    "\n"
)


# ============================================================
# 12. PREPARE OUTPUT TABLES
# ============================================================

cat("\n=== PREPARING OUTPUT TABLES ===\n")

fgsea_export <- fgsea_results

fgsea_export$leadingEdge <- vapply(
    fgsea_export$leadingEdge,
    function(x) {
        paste(x, collapse = ";")
    },
    character(1)
)

significant_export <- fgsea_export[
    !is.na(fgsea_export$padj) &
    fgsea_export$padj < FDR_THRESHOLD,
]

positive_export <- fgsea_export[
    !is.na(fgsea_export$padj) &
    fgsea_export$padj < FDR_THRESHOLD &
    fgsea_export$NES > 0,
]

negative_export <- fgsea_export[
    !is.na(fgsea_export$padj) &
    fgsea_export$padj < FDR_THRESHOLD &
    fgsea_export$NES < 0,
]


# ============================================================
# 13. GENERATE SUMMARY
# ============================================================

summary <- data.frame(
    metric = c(
        "genes_ranked",
        "positive_ranking_values",
        "negative_ranking_values",
        "tied_ranking_values",
        "all_gene_pathway_assignments",
        "tested_gene_pathway_assignments",
        "excluded_unranked_assignments",
        "KEGG_pathway_sets",
        "pathways_tested",
        "significant_pathways_FDR05",
        "pathways_enriched_in_resistant",
        "pathways_enriched_in_sensitive",
        "FDR_threshold",
        "minimum_pathway_size",
        "maximum_pathway_size"
    ),
    value = c(
        length(stats),
        sum(stats > 0),
        sum(stats < 0),
        tie_count,
        nrow(pathway_membership_all),
        nrow(pathway_membership),
        excluded_assignments,
        length(pathways),
        nrow(fgsea_results),
        nrow(significant_results),
        nrow(positive_results),
        nrow(negative_results),
        FDR_THRESHOLD,
        MIN_SIZE,
        MAX_SIZE
    )
)


# ============================================================
# 14. FINAL VALIDATION
# ============================================================

cat("\n=== FINAL VALIDATION ===\n")

stopifnot(length(stats) == EXPECTED_GENES)

stopifnot(
    all(pathway_membership$gene_id %in% names(stats))
)

stopifnot(
    !anyDuplicated(pathway_membership)
)

stopifnot(
    all(lengths(pathways_tested) >= MIN_SIZE)
)

stopifnot(
    all(lengths(pathways_tested) <= MAX_SIZE)
)

stopifnot(
    all(fgsea_results$pathway %in% names(pathways_tested))
)

stopifnot(
    nrow(significant_results) ==
    sum(fgsea_results$padj < FDR_THRESHOLD, na.rm = TRUE)
)

cat("All validation checks passed.\n")


# ============================================================
# 15. EXPORT RESULTS
# ============================================================

cat("\n=== EXPORTING RESULTS ===\n")

fwrite(
    ranking,
    file.path(
        output_dir,
        "edgeR_signed_gene_ranking_grouped.csv"
    )
)

fwrite(
    fgsea_export,
    file.path(
        output_dir,
        "fgsea_all_pathways_grouped.csv"
    )
)

fwrite(
    significant_export,
    file.path(
        output_dir,
        "fgsea_significant_pathways_FDR05_grouped.csv"
    )
)

fwrite(
    positive_export,
    file.path(
        output_dir,
        "fgsea_enriched_in_resistant_grouped.csv"
    )
)

fwrite(
    negative_export,
    file.path(
        output_dir,
        "fgsea_enriched_in_sensitive_grouped.csv"
    )
)

fwrite(
    summary,
    file.path(
        output_dir,
        "fgsea_summary_grouped.tsv"
    ),
    sep = "\t"
)


# ============================================================
# 16. FINAL REPORT
# ============================================================

cat("\n========================================\n")
cat("GROUPED FGSEA ANALYSIS COMPLETED\n")
cat("========================================\n")

cat("Genes ranked:", length(stats), "\n")
cat("Pathways tested:", nrow(fgsea_results), "\n")
cat("Significant pathways:", nrow(significant_results), "\n")
cat("Enriched in resistant:", nrow(positive_results), "\n")
cat("Enriched in sensitive:", nrow(negative_results), "\n")

cat("\nOutput directory:\n")
cat(output_dir, "\n")

cat("\nAll files were generated and validated.\n")
