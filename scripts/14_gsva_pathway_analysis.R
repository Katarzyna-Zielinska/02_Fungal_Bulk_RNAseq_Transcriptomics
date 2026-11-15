
#!/usr/bin/env Rscript

# ============================================================
# Project 2: Candida auris Bulk RNA-seq Transcriptomics
# Step 14: KEGG pathway activity analysis using GSVA
# Contrast: Resistant vs Sensitive
# Method: TMM-normalized logCPM + GSVA + limma
# ============================================================


# ============================================================
# 1. CONFIGURATION
# ============================================================

project_dir <- normalizePath(".", mustWork = TRUE)

dge_file <- file.path(
    project_dir,
    "results/differential_expression/edgeR_DGEList_TMM.rds"
)

integrated_file <- file.path(
    project_dir,
    "results/integrated_analysis/integrated_DE_functional_annotation.csv"
)

metadata_file <- file.path(
    project_dir,
    "data/metadata/sample_metadata.csv"
)

output_dir <- file.path(
    project_dir,
    "results/pathway_enrichment/gsva"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

FDR_THRESHOLD <- 0.05
MIN_SIZE <- 10
MAX_SIZE <- 500

set.seed(123)


# ============================================================
# 2. LOAD PACKAGES
# ============================================================

cat("\n=== LOADING PACKAGES ===\n")

required_packages <- c(
    "GSVA",
    "edgeR",
    "limma",
    "data.table"
)

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
    library(GSVA)
    library(edgeR)
    library(limma)
    library(data.table)
})


# ============================================================
# 3. CHECK INPUT FILES
# ============================================================

cat("\n=== CHECKING INPUT FILES ===\n")

input_files <- c(
    dge_file,
    integrated_file,
    metadata_file
)

for (f in input_files) {

    if (!file.exists(f)) {
        stop("Input file not found: ", f)
    }

    cat("[FOUND]", f, "\n")
}


# ============================================================
# 4. LOAD DATA
# ============================================================

cat("\n=== LOADING DATA ===\n")

y <- readRDS(dge_file)

integrated <- fread(
    integrated_file,
    data.table = FALSE
)

metadata <- read.csv(
    metadata_file,
    stringsAsFactors = FALSE
)

cat("DGEList dimensions:", paste(dim(y), collapse = " x "), "\n")
cat("Integrated annotation:", nrow(integrated), "genes\n")
cat("Metadata samples:", nrow(metadata), "\n")


# ============================================================
# 5. VALIDATE INPUT DATA
# ============================================================

cat("\n=== VALIDATING INPUT DATA ===\n")

if (!inherits(y, "DGEList")) {
    stop("Input object is not an edgeR DGEList.")
}

if (nrow(y) != 5376 || ncol(y) != 6) {
    stop(
        "Expected a DGEList with 5376 genes and 6 samples."
    )
}

required_integrated <- c(
    "gene_id",
    "KEGG_Pathway"
)

if (!all(required_integrated %in% colnames(integrated))) {
    stop("Required columns missing in integrated annotation.")
}

required_metadata <- c(
    "sample",
    "condition"
)

if (!all(required_metadata %in% colnames(metadata))) {
    stop("Required columns missing in sample metadata.")
}

if (anyDuplicated(integrated$gene_id) > 0) {
    stop("Duplicated gene IDs in integrated annotation.")
}

if (anyDuplicated(metadata$sample) > 0) {
    stop("Duplicated sample IDs in metadata.")
}

if (anyNA(integrated$gene_id)) {
    stop("Missing gene IDs in integrated annotation.")
}

metadata$condition <- factor(
    metadata$condition,
    levels = c("Sensitive", "Resistant")
)

if (anyNA(metadata$condition)) {
    stop("Unexpected or missing sample condition.")
}

rownames(metadata) <- metadata$sample

if (!setequal(colnames(y), rownames(metadata))) {
    stop("Sample IDs do not match between DGEList and metadata.")
}

metadata <- metadata[
    colnames(y),
    ,
    drop = FALSE
]

if (!identical(colnames(y), rownames(metadata))) {
    stop("Sample order mismatch.")
}

if (!identical(
    as.character(y$samples$group),
    as.character(metadata$condition)
)) {
    stop("DGEList group assignments differ from metadata.")
}

if (!all(table(metadata$condition) == c(3, 3))) {
    stop("Expected 3 Sensitive and 3 Resistant samples.")
}

cat("Gene and sample validation passed.\n")


# ============================================================
# 6. CREATE NORMALIZED LOGCPM MATRIX
# ============================================================

cat("\n=== CREATING NORMALIZED LOGCPM MATRIX ===\n")

cat("Normalization method: TMM\n")
cat("Using saved edgeR normalization factors.\n")

print(
    y$samples[, c("group", "lib.size", "norm.factors")]
)

logCPM <- edgeR::cpm(
    y,
    log = TRUE,
    prior.count = 2,
    normalized.lib.sizes = TRUE
)

if (!is.matrix(logCPM)) {
    stop("logCPM output is not a matrix.")
}

if (anyNA(logCPM) || any(!is.finite(logCPM))) {
    stop("Missing or non-finite values in logCPM matrix.")
}

if (nrow(logCPM) != 5376 || ncol(logCPM) != 6) {
    stop("Unexpected logCPM matrix dimensions.")
}

if (!identical(rownames(logCPM), rownames(y))) {
    stop("Gene ID order changed during logCPM transformation.")
}

if (!identical(colnames(logCPM), colnames(y))) {
    stop("Sample order changed during logCPM transformation.")
}

cat("Genes:", nrow(logCPM), "\n")
cat("Samples:", ncol(logCPM), "\n")
cat("logCPM matrix successfully generated.\n")

write.csv(
    data.frame(
        gene_id = rownames(logCPM),
        logCPM,
        check.names = FALSE
    ),
    file.path(output_dir, "TMM_normalized_logCPM.csv"),
    row.names = FALSE
)


# ============================================================
# 7. EXTRACT KEGG PATHWAY ASSIGNMENTS
# ============================================================

cat("\n=== EXTRACTING KEGG PATHWAYS ===\n")

pathway_rows <- list()

for (i in seq_len(nrow(integrated))) {

    gene_id <- integrated$gene_id[i]

    pathway_value <- integrated$KEGG_Pathway[i]

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
    # Exclude corresponding mapXXXXX identifiers.

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
    "All gene-pathway assignments:",
    nrow(pathway_membership_all),
    "\n"
)

cat(
    "Unique KEGG pathway IDs:",
    length(unique(pathway_membership_all$pathway_id)),
    "\n"
)

write.csv(
    pathway_membership_all,
    file.path(
        output_dir,
        "KEGG_ko_gene_pathway_membership_all_genes.csv"
    ),
    row.names = FALSE
)


# ============================================================
# 8. FILTER MEMBERSHIP TO TESTED GENES
# ============================================================

cat("\n=== FILTERING TO TESTED GENES ===\n")

tested_gene_ids <- rownames(logCPM)

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
    "Assignments excluded:",
    excluded_assignments,
    "\n"
)

if (nrow(pathway_membership) == 0) {
    stop("No KEGG assignments remain after filtering.")
}

write.csv(
    pathway_membership,
    file.path(
        output_dir,
        "KEGG_ko_gene_pathway_membership_tested_genes.csv"
    ),
    row.names = FALSE
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

cat("All pathway sets:", length(pathways), "\n")
cat("Pathways meeting size criteria:", length(pathways_tested), "\n")
cat("Minimum pathway size:", MIN_SIZE, "\n")
cat("Maximum pathway size:", MAX_SIZE, "\n")

if (length(pathways_tested) == 0) {
    stop("No pathways meet the selected size criteria.")
}

write.csv(
    pathway_sizes,
    file.path(
        output_dir,
        "KEGG_pathway_gene_set_sizes.csv"
    ),
    row.names = FALSE
)


# ============================================================
# 10. VALIDATE PATHWAY GENE SETS
# ============================================================

cat("\n=== VALIDATING PATHWAY GENE SETS ===\n")

if (!all(pathway_membership$gene_id %in% rownames(logCPM))) {
    stop("Unmatched gene IDs remain in pathway membership.")
}

for (pathway_id in names(pathways_tested)) {

    genes <- pathways_tested[[pathway_id]]

    if (anyDuplicated(genes) > 0) {
        stop("Duplicated genes in pathway: ", pathway_id)
    }

    if (!all(genes %in% rownames(logCPM))) {
        stop("Unmatched genes in pathway: ", pathway_id)
    }
}

cat("All pathway gene sets validated.\n")


# ============================================================
# 11. RUN GSVA
# ============================================================

cat("\n=== RUNNING GSVA ===\n")

cat("GSVA version:", as.character(packageVersion("GSVA")), "\n")
cat("Method: GSVA\n")
cat("Input: TMM-normalized logCPM\n")
cat("Kernel: Gaussian\n")

gsva_parameters <- GSVA::gsvaParam(
    exprData = logCPM,
    geneSets = pathways_tested,
    minSize = MIN_SIZE,
    maxSize = MAX_SIZE,
    kcdf = "Gaussian",
    maxDiff = TRUE,
    verbose = TRUE
)

gsva_scores <- GSVA::gsva(
    gsva_parameters,
    verbose = TRUE
)

gsva_scores <- as.matrix(gsva_scores)

if (nrow(gsva_scores) == 0 || ncol(gsva_scores) != 6) {
    stop("Unexpected GSVA score matrix dimensions.")
}

if (anyNA(gsva_scores) || any(!is.finite(gsva_scores))) {
    stop("Missing or non-finite GSVA scores.")
}

if (!identical(colnames(gsva_scores), colnames(logCPM))) {
    stop("Sample order mismatch in GSVA results.")
}

cat("Pathways with GSVA scores:", nrow(gsva_scores), "\n")
cat("Samples:", ncol(gsva_scores), "\n")

write.csv(
    data.frame(
        pathway_id = rownames(gsva_scores),
        gsva_scores,
        check.names = FALSE
    ),
    file.path(
        output_dir,
        "GSVA_pathway_scores.csv"
    ),
    row.names = FALSE
)

saveRDS(
    gsva_scores,
    file.path(
        output_dir,
        "GSVA_pathway_scores.rds"
    )
)

saveRDS(
    gsva_parameters,
    file.path(
        output_dir,
        "GSVA_parameters.rds"
    )
)


# ============================================================
# 12. DIFFERENTIAL PATHWAY ACTIVITY USING LIMMA
# ============================================================

cat("\n=== DIFFERENTIAL PATHWAY ACTIVITY ===\n")

design <- model.matrix(
    ~ condition,
    data = metadata
)

rownames(design) <- rownames(metadata)

if (qr(design)$rank != ncol(design)) {
    stop("Design matrix is not full rank.")
}

cat("Experimental design:\n")
print(design)

contrast <- limma::makeContrasts(
    conditionResistant,
    levels = design
)

fit <- limma::lmFit(
    gsva_scores,
    design
)

fit <- limma::contrasts.fit(
    fit,
    contrast
)

fit <- limma::eBayes(
    fit,
    robust = TRUE
)

pathway_results <- limma::topTable(
    fit,
    coef = 1,
    number = Inf,
    sort.by = "P"
)

pathway_results$pathway_id <- rownames(pathway_results)

pathway_results <- pathway_results[
    ,
    c(
        "pathway_id",
        "logFC",
        "AveExpr",
        "t",
        "P.Value",
        "adj.P.Val",
        "B"
    )
]

colnames(pathway_results) <- c(
    "pathway_id",
    "logFC",
    "AveExpr",
    "t",
    "PValue",
    "FDR",
    "B"
)

pathway_results <- pathway_results[
    order(pathway_results$FDR, pathway_results$PValue),
]

significant_results <- pathway_results[
    !is.na(pathway_results$FDR) &
    pathway_results$FDR < FDR_THRESHOLD,
]

positive_results <- significant_results[
    significant_results$logFC > 0,
]

negative_results <- significant_results[
    significant_results$logFC < 0,
]

cat("Pathways tested:", nrow(pathway_results), "\n")
cat(
    "Significant pathways (FDR < 0.05):",
    nrow(significant_results),
    "\n"
)

cat(
    "Higher GSVA scores in Resistant:",
    nrow(positive_results),
    "\n"
)

cat(
    "Higher GSVA scores in Sensitive:",
    nrow(negative_results),
    "\n"
)

saveRDS(
    fit,
    file.path(
        output_dir,
        "GSVA_limma_fit.rds"
    )
)


# ============================================================
# 13. PREPARE OUTPUT TABLES
# ============================================================

cat("\n=== PREPARING OUTPUT TABLES ===\n")

write.csv(
    pathway_results,
    file.path(
        output_dir,
        "GSVA_limma_all_pathways.csv"
    ),
    row.names = FALSE
)

write.csv(
    significant_results,
    file.path(
        output_dir,
        "GSVA_significant_pathways_FDR05.csv"
    ),
    row.names = FALSE
)

write.csv(
    positive_results,
    file.path(
        output_dir,
        "GSVA_higher_activity_in_resistant.csv"
    ),
    row.names = FALSE
)

write.csv(
    negative_results,
    file.path(
        output_dir,
        "GSVA_higher_activity_in_sensitive.csv"
    ),
    row.names = FALSE
)


# ============================================================
# 14. VISUALIZE GSVA SCORES
# ============================================================

cat("\n=== GENERATING GSVA VISUALIZATIONS ===\n")

# Select up to 30 pathways with the lowest adjusted P-values.
# If no pathway is significant, display the 30 lowest raw P-values.

if (nrow(significant_results) > 0) {

    plot_pathways <- significant_results$pathway_id

} else {

    plot_pathways <- pathway_results$pathway_id
}

plot_pathways <- head(plot_pathways, 30)

plot_matrix <- gsva_scores[
    plot_pathways,
    ,
    drop = FALSE
]

# Row-wise Z-score transformation for visualization only.

plot_matrix_z <- t(
    scale(t(plot_matrix))
)

plot_matrix_z[!is.finite(plot_matrix_z)] <- 0

pdf(
    file.path(
        output_dir,
        "GSVA_pathway_activity_heatmap.pdf"
    ),
    width = 11,
    height = max(7, 0.25 * length(plot_pathways) + 3)
)

if (nrow(plot_matrix_z) >= 2) {

    heatmap(
        plot_matrix_z,
        scale = "none",
        Rowv = NA,
        Colv = NA,
        margins = c(9, 12),
        main = "GSVA pathway activity",
        xlab = "Samples",
        ylab = "KEGG pathway"
    )

} else {

    plot(
        1,
        type = "n",
        main = "Insufficient pathways for heatmap"
    )

}

dev.off()


# ============================================================
# 15. GENERATE SAMPLE-LEVEL SUMMARY
# ============================================================

cat("\n=== GENERATING SAMPLE SUMMARY ===\n")

sample_summary <- data.frame(
    sample = colnames(gsva_scores),
    condition = metadata[colnames(gsva_scores), "condition"],
    mean_GSVA_score = colMeans(gsva_scores),
    median_GSVA_score = apply(
        gsva_scores,
        2,
        median
    ),
    stringsAsFactors = FALSE
)

write.csv(
    sample_summary,
    file.path(
        output_dir,
        "GSVA_sample_summary.csv"
    ),
    row.names = FALSE
)


# ============================================================
# 16. GENERATE ANALYSIS SUMMARY
# ============================================================

cat("\n=== GENERATING ANALYSIS SUMMARY ===\n")

summary_table <- data.frame(
    metric = c(
        "GSVA_version",
        "normalization",
        "expression_transformation",
        "GSVA_kernel",
        "genes_tested",
        "samples",
        "sensitive_samples",
        "resistant_samples",
        "all_gene_pathway_assignments",
        "tested_gene_pathway_assignments",
        "excluded_unmatched_assignments",
        "all_KEGG_pathway_sets",
        "pathways_meeting_size_criteria",
        "pathways_scored_by_GSVA",
        "pathways_tested_by_limma",
        "significant_pathways_FDR05",
        "higher_activity_in_resistant",
        "higher_activity_in_sensitive",
        "FDR_threshold",
        "minimum_pathway_size",
        "maximum_pathway_size"
    ),
    value = c(
        as.character(packageVersion("GSVA")),
        "TMM",
        "logCPM",
        "Gaussian",
        nrow(logCPM),
        ncol(logCPM),
        sum(metadata$condition == "Sensitive"),
        sum(metadata$condition == "Resistant"),
        nrow(pathway_membership_all),
        nrow(pathway_membership),
        excluded_assignments,
        length(pathways),
        length(pathways_tested),
        nrow(gsva_scores),
        nrow(pathway_results),
        nrow(significant_results),
        nrow(positive_results),
        nrow(negative_results),
        FDR_THRESHOLD,
        MIN_SIZE,
        MAX_SIZE
    ),
    stringsAsFactors = FALSE
)

write.table(
    summary_table,
    file.path(
        output_dir,
        "GSVA_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 17. FINAL VALIDATION
# ============================================================

cat("\n=== FINAL VALIDATION ===\n")

stopifnot(nrow(logCPM) == 5376)
stopifnot(ncol(logCPM) == 6)

stopifnot(
    all(pathway_membership$gene_id %in% rownames(logCPM))
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
    all(rownames(gsva_scores) %in% names(pathways_tested))
)

stopifnot(
    identical(colnames(gsva_scores), colnames(logCPM))
)

stopifnot(
    nrow(pathway_results) == nrow(gsva_scores)
)

stopifnot(
    nrow(significant_results) ==
    sum(pathway_results$FDR < FDR_THRESHOLD, na.rm = TRUE)
)

cat("All validation checks passed.\n")


# ============================================================
# 18. FINAL REPORT
# ============================================================

cat("\n========================================\n")
cat("GSVA ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("========================================\n")

cat("Genes analyzed:", nrow(logCPM), "\n")
cat("Samples analyzed:", ncol(logCPM), "\n")
cat("Pathways scored:", nrow(gsva_scores), "\n")
cat("Pathways tested by limma:", nrow(pathway_results), "\n")
cat("Significant pathways:", nrow(significant_results), "\n")
cat("Higher activity in Resistant:", nrow(positive_results), "\n")
cat("Higher activity in Sensitive:", nrow(negative_results), "\n")

cat("\n=== TOP 10 PATHWAYS ===\n")
print(head(pathway_results, 10))

cat("\nOutput directory:\n")
cat(output_dir, "\n")

cat("\nAll files were generated and validated.\n")

