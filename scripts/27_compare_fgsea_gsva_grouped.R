#!/usr/bin/env Rscript

# ============================================================
# Project 2: Candida auris Bulk RNA-seq Transcriptomics
# Step 27: Comparison of FGSEA and GSVA pathway analyses
# Analysis: Grouped transcript/gene mapping
# ============================================================


# ============================================================
# 1. CONFIGURATION
# ============================================================

project_dir <- normalizePath(".", mustWork = TRUE)

fgsea_file <- file.path(
    project_dir,
    "results/pathway_enrichment_grouped/fgsea/fgsea_all_pathways_grouped.csv"
)

gsva_file <- file.path(
    project_dir,
    "results/pathway_enrichment_grouped/gsva/GSVA_limma_all_pathways_grouped.csv"
)

output_dir <- file.path(
    project_dir,
    "results/pathway_enrichment_grouped/method_comparison"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

FDR_THRESHOLD <- 0.05


# ============================================================
# 2. LOAD PACKAGES
# ============================================================

cat("\n=== LOADING PACKAGES ===\n")

if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package data.table is not installed.")
}

suppressPackageStartupMessages({
    library(data.table)
})


# ============================================================
# 3. CHECK INPUT FILES
# ============================================================

cat("\n=== CHECKING INPUT FILES ===\n")

if (!file.exists(fgsea_file)) {
    stop("FGSEA results file not found: ", fgsea_file)
}

if (!file.exists(gsva_file)) {
    stop("GSVA results file not found: ", gsva_file)
}

cat("[FOUND] FGSEA results.\n")
cat("[FOUND] GSVA results.\n")


# ============================================================
# 4. LOAD DATA
# ============================================================

cat("\n=== LOADING RESULTS ===\n")

fgsea <- fread(
    fgsea_file,
    data.table = FALSE
)

gsva <- fread(
    gsva_file,
    data.table = FALSE
)

cat("FGSEA pathways:", nrow(fgsea), "\n")
cat("GSVA pathways:", nrow(gsva), "\n")


# ============================================================
# 5. VALIDATE INPUT DATA
# ============================================================

cat("\n=== VALIDATING INPUT DATA ===\n")

required_fgsea <- c(
    "pathway",
    "NES",
    "pval",
    "padj"
)

required_gsva <- c(
    "pathway_id",
    "logFC",
    "PValue",
    "FDR"
)

if (!all(required_fgsea %in% colnames(fgsea))) {
    stop("Required columns missing in FGSEA results.")
}

if (!all(required_gsva %in% colnames(gsva))) {
    stop("Required columns missing in GSVA results.")
}

if (anyDuplicated(fgsea$pathway) > 0) {
    stop("Duplicated pathway IDs in FGSEA results.")
}

if (anyDuplicated(gsva$pathway_id) > 0) {
    stop("Duplicated pathway IDs in GSVA results.")
}

if (
    anyNA(fgsea$NES) ||
    anyNA(fgsea$padj) ||
    anyNA(gsva$logFC) ||
    anyNA(gsva$FDR)
) {
    stop("Missing values detected in required statistics.")
}

if (
    any(!is.finite(fgsea$NES)) ||
    any(!is.finite(fgsea$padj)) ||
    any(!is.finite(gsva$logFC)) ||
    any(!is.finite(gsva$FDR))
) {
    stop("Non-finite values detected.")
}

if (any(fgsea$padj < 0 | fgsea$padj > 1)) {
    stop("Invalid FGSEA adjusted P-values.")
}

if (any(gsva$FDR < 0 | gsva$FDR > 1)) {
    stop("Invalid GSVA adjusted P-values.")
}

cat("Input validation passed.\n")


# ============================================================
# 6. MATCH PATHWAYS
# ============================================================

cat("\n=== MATCHING PATHWAYS ===\n")

common_pathways <- intersect(
    fgsea$pathway,
    gsva$pathway_id
)

fgsea_only <- setdiff(
    fgsea$pathway,
    gsva$pathway_id
)

gsva_only <- setdiff(
    gsva$pathway_id,
    fgsea$pathway
)

cat("Common pathways:", length(common_pathways), "\n")
cat("Only in FGSEA:", length(fgsea_only), "\n")
cat("Only in GSVA:", length(gsva_only), "\n")

if (length(common_pathways) == 0) {
    stop("No common pathways found.")
}

all_pathways <- merge(
    fgsea,
    gsva,
    by.x = "pathway",
    by.y = "pathway_id",
    all = TRUE,
    sort = FALSE
)

if (
    nrow(all_pathways) !=
    length(union(fgsea$pathway, gsva$pathway_id))
) {
    stop("Unexpected number of pathways after merging.")
}


# ============================================================
# 7. STANDARDIZE AND CLASSIFY RESULTS
# ============================================================

cat("\n=== CLASSIFYING PATHWAY RESULTS ===\n")

all_pathways$fgsea_significant <- (
    !is.na(all_pathways$padj) &
    all_pathways$padj < FDR_THRESHOLD
)

all_pathways$GSVA_significant <- (
    !is.na(all_pathways$FDR) &
    all_pathways$FDR < FDR_THRESHOLD
)

all_pathways$fgsea_direction <- ifelse(
    is.na(all_pathways$NES),
    NA_character_,
    ifelse(
        all_pathways$NES > 0,
        "Resistant",
        ifelse(
            all_pathways$NES < 0,
            "Sensitive",
            "Neutral"
        )
    )
)

all_pathways$GSVA_direction <- ifelse(
    is.na(all_pathways$logFC),
    NA_character_,
    ifelse(
        all_pathways$logFC > 0,
        "Resistant",
        ifelse(
            all_pathways$logFC < 0,
            "Sensitive",
            "Neutral"
        )
    )
)

all_pathways$direction_concordant <- ifelse(
    is.na(all_pathways$fgsea_direction) |
    is.na(all_pathways$GSVA_direction),
    NA,
    all_pathways$fgsea_direction ==
        all_pathways$GSVA_direction
)

all_pathways$significance_class <- ifelse(
    all_pathways$fgsea_significant &
        all_pathways$GSVA_significant,
    "Significant in both methods",
    ifelse(
        all_pathways$fgsea_significant &
            !all_pathways$GSVA_significant,
        "Significant only in FGSEA",
        ifelse(
            !all_pathways$fgsea_significant &
                all_pathways$GSVA_significant,
            "Significant only in GSVA",
            "Not significant in either method"
        )
    )
)

all_pathways$comparison_class <- ifelse(
    all_pathways$fgsea_significant &
        all_pathways$GSVA_significant &
        all_pathways$direction_concordant,
    "Significant in both - concordant",
    ifelse(
        all_pathways$fgsea_significant &
            all_pathways$GSVA_significant &
            !all_pathways$direction_concordant,
        "Significant in both - discordant",
        ifelse(
            all_pathways$fgsea_significant,
            "Significant only in FGSEA",
            ifelse(
                all_pathways$GSVA_significant,
                "Significant only in GSVA",
                ifelse(
                    all_pathways$direction_concordant,
                    "Concordant - not significant in both",
                    "Not significant - discordant"
                )
            )
        )
    )
)

all_pathways <- all_pathways[
    order(
        all_pathways$padj,
        all_pathways$FDR,
        na.last = TRUE
    ),
]


# ============================================================
# 8. CALCULATE DIRECTIONAL AGREEMENT
# ============================================================

cat("\n=== DIRECTIONAL AGREEMENT ===\n")

matched <- all_pathways[
    !is.na(all_pathways$NES) &
    !is.na(all_pathways$logFC),
]

matched_nonzero <- matched[
    matched$NES != 0 &
    matched$logFC != 0,
]

n_matched <- nrow(matched)
n_directional <- nrow(matched_nonzero)

n_concordant <- sum(
    matched_nonzero$direction_concordant
)

n_discordant <- n_directional - n_concordant

direction_agreement <- if (n_directional > 0) {
    n_concordant / n_directional
} else {
    NA_real_
}

cat("Matched pathways:", n_matched, "\n")
cat("Non-zero directional comparisons:", n_directional, "\n")
cat("Concordant directions:", n_concordant, "\n")
cat("Discordant directions:", n_discordant, "\n")

cat(
    "Directional agreement:",
    round(direction_agreement * 100, 2),
    "%\n"
)


# ============================================================
# 9. CALCULATE SPEARMAN CORRELATION
# ============================================================

cat("\n=== SPEARMAN CORRELATION ===\n")

if (n_matched >= 3) {

    correlation <- cor.test(
        matched$NES,
        matched$logFC,
        method = "spearman",
        exact = FALSE
    )

    spearman_rho <- unname(correlation$estimate)
    spearman_pvalue <- correlation$p.value

} else {

    spearman_rho <- NA_real_
    spearman_pvalue <- NA_real_
}

cat("Spearman rho:", spearman_rho, "\n")
cat("P-value:", spearman_pvalue, "\n")


# ============================================================
# 10. SUMMARIZE SIGNIFICANCE OVERLAP
# ============================================================

cat("\n=== SIGNIFICANCE OVERLAP ===\n")

n_both <- sum(
    matched$fgsea_significant &
    matched$GSVA_significant
)

n_fgsea_only <- sum(
    matched$fgsea_significant &
    !matched$GSVA_significant
)

n_gsva_only <- sum(
    !matched$fgsea_significant &
    matched$GSVA_significant
)

n_neither <- sum(
    !matched$fgsea_significant &
    !matched$GSVA_significant
)

significance_summary <- data.frame(
    category = c(
        "Significant in both methods",
        "Significant only in FGSEA",
        "Significant only in GSVA",
        "Not significant in either method"
    ),
    pathways = c(
        n_both,
        n_fgsea_only,
        n_gsva_only,
        n_neither
    )
)

print(significance_summary)

direction_summary <- data.frame(
    category = c(
        "Concordant",
        "Discordant"
    ),
    pathways = c(
        n_concordant,
        n_discordant
    )
)


# ============================================================
# 11. GENERATE SUMMARY REPORT
# ============================================================

summary_table <- data.frame(
    metric = c(
        "FGSEA_pathways",
        "GSVA_pathways",
        "common_pathways",
        "FGSEA_only_pathways",
        "GSVA_only_pathways",
        "FGSEA_significant_FDR05",
        "GSVA_significant_FDR05",
        "significant_in_both",
        "significant_only_FGSEA",
        "significant_only_GSVA",
        "not_significant_in_either",
        "matched_pathways",
        "nonzero_directional_comparisons",
        "directionally_concordant",
        "directionally_discordant",
        "directional_agreement_percent",
        "spearman_rho",
        "spearman_pvalue",
        "FDR_threshold"
    ),
    value = c(
        nrow(fgsea),
        nrow(gsva),
        length(common_pathways),
        length(fgsea_only),
        length(gsva_only),
        sum(fgsea$padj < FDR_THRESHOLD),
        sum(gsva$FDR < FDR_THRESHOLD),
        n_both,
        n_fgsea_only,
        n_gsva_only,
        n_neither,
        n_matched,
        n_directional,
        n_concordant,
        n_discordant,
        direction_agreement * 100,
        spearman_rho,
        spearman_pvalue,
        FDR_THRESHOLD
    ),
    stringsAsFactors = FALSE
)


# ============================================================
# 12. EXPORT COMPARISON TABLES
# ============================================================

cat("\n=== EXPORTING TABLES ===\n")

fwrite(
    all_pathways,
    file.path(
        output_dir,
        "fgsea_GSVA_comparison_all_pathways_grouped.csv"
    )
)

fwrite(
    matched,
    file.path(
        output_dir,
        "fgsea_GSVA_comparison_common_pathways_grouped.csv"
    )
)

fwrite(
    matched[
        matched$comparison_class ==
            "Significant in both - concordant",
    ],
    file.path(
        output_dir,
        "pathways_significant_both_concordant_grouped.csv"
    )
)

fwrite(
    matched[
        matched$comparison_class ==
            "Significant in both - discordant",
    ],
    file.path(
        output_dir,
        "pathways_significant_both_discordant_grouped.csv"
    )
)

fwrite(
    significance_summary,
    file.path(
        output_dir,
        "pathway_significance_overlap_grouped.csv"
    )
)

fwrite(
    direction_summary,
    file.path(
        output_dir,
        "pathway_directional_concordance_grouped.csv"
    )
)

fwrite(
    summary_table,
    file.path(
        output_dir,
        "fgsea_GSVA_comparison_summary_grouped.tsv"
    ),
    sep = "\t"
)


# ============================================================
# 13. GENERATE VISUALIZATIONS
# ============================================================

cat("\n=== GENERATING VISUALIZATIONS ===\n")

# 13.1 NES vs GSVA logFC scatterplot

pdf(
    file.path(
        output_dir,
        "fgsea_GSVA_direction_correlation_grouped.pdf"
    ),
    width = 8,
    height = 7
)

plot(
    matched$NES,
    matched$logFC,
    pch = 19,
    xlab = "FGSEA Normalized Enrichment Score (NES)",
    ylab = "GSVA logFC",
    main = "FGSEA vs GSVA pathway direction - grouped"
)

abline(
    h = 0,
    v = 0,
    lty = 2
)

if (!is.na(spearman_rho)) {
    legend(
        "topleft",
        legend = paste(
            "Spearman rho =",
            round(spearman_rho, 3)
        ),
        bty = "n"
    )
}

dev.off()


# 13.2 Significance overlap barplot

pdf(
    file.path(
        output_dir,
        "fgsea_GSVA_significance_overlap_grouped.pdf"
    ),
    width = 9,
    height = 6
)

barplot(
    significance_summary$pathways,
    names.arg = significance_summary$category,
    las = 2,
    ylab = "Number of pathways",
    main = "Significance overlap: FGSEA and GSVA - grouped",
    cex.names = 0.8
)

dev.off()


# 13.3 Directional concordance

pdf(
    file.path(
        output_dir,
        "fgsea_GSVA_directional_concordance_grouped.pdf"
    ),
    width = 7,
    height = 6
)

barplot(
    direction_summary$pathways,
    names.arg = direction_summary$category,
    ylab = "Number of pathways",
    main = "Directional concordance - grouped"
)

dev.off()


# ============================================================
# 14. FINAL VALIDATION
# ============================================================

cat("\n=== FINAL VALIDATION ===\n")

stopifnot(length(common_pathways) > 0)

stopifnot(
    nrow(all_pathways) ==
    length(union(fgsea$pathway, gsva$pathway_id))
)

stopifnot(
    n_both + n_fgsea_only + n_gsva_only + n_neither ==
    n_matched
)

stopifnot(
    n_concordant + n_discordant ==
    n_directional
)

stopifnot(
    nrow(matched) == length(common_pathways)
)

cat("All validation checks passed.\n")


# ============================================================
# 15. FINAL REPORT
# ============================================================

cat("\n========================================\n")
cat("GROUPED FGSEA-GSVA COMPARISON COMPLETED\n")
cat("========================================\n")

cat("Common pathways:", length(common_pathways), "\n")
cat("Significant in both:", n_both, "\n")
cat("Significant only in FGSEA:", n_fgsea_only, "\n")
cat("Significant only in GSVA:", n_gsva_only, "\n")
cat("Not significant in either:", n_neither, "\n")
cat(
    "Directional agreement:",
    round(direction_agreement * 100, 2),
    "%\n"
)
cat("Spearman rho:", spearman_rho, "\n")

cat("\nOutput directory:\n")
cat(output_dir, "\n")

cat("\nAll files were generated and validated.\n")
