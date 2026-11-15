# ============================================================
# PROJECT 2: Candida auris Bulk RNA-seq
# Comparison: grouped edgeR vs grouped limma-voom
# ============================================================

project <- "/mnt/e/Fungal_Bioinformatics_Portfolio/02_Fungal_Bulk_RNAseq_Transcriptomics"

de_dir <- file.path(project, "results/differential_expression_grouped")

edgeR_file <- file.path(de_dir, "edgeR_all_genes_grouped.csv")
limma_file <- file.path(de_dir, "limma_voom_all_genes_grouped.csv")

stopifnot(file.exists(edgeR_file), file.exists(limma_file))

edgeR <- read.csv(edgeR_file, check.names = FALSE)
limma <- read.csv(limma_file, check.names = FALSE)

edgeR <- edgeR[, c("gene_id", "logFC", "PValue", "FDR")]
limma <- limma[, c("gene_id", "logFC", "P.Value", "adj.P.Val")]

colnames(edgeR) <- c(
  "gene_id", "edgeR_logFC", "edgeR_PValue", "edgeR_FDR"
)

colnames(limma) <- c(
  "gene_id", "limma_logFC", "limma_PValue", "limma_FDR"
)

comparison <- merge(edgeR, limma, by = "gene_id")

if (nrow(comparison) != 5376) {
  stop("Unexpected number of genes in comparison: ", nrow(comparison))
}

comparison$edgeR_significant <- comparison$edgeR_FDR < 0.05
comparison$limma_significant <- comparison$limma_FDR < 0.05

comparison$edgeR_direction <- ifelse(
  comparison$edgeR_logFC > 0, "Up",
  ifelse(comparison$edgeR_logFC < 0, "Down", "Zero")
)

comparison$limma_direction <- ifelse(
  comparison$limma_logFC > 0, "Up",
  ifelse(comparison$limma_logFC < 0, "Down", "Zero")
)

comparison$direction_concordant <-
  comparison$edgeR_direction == comparison$limma_direction

comparison$significance_group <- ifelse(
  comparison$edgeR_significant & comparison$limma_significant,
  "Significant_both",
  ifelse(
    comparison$edgeR_significant,
    "edgeR_only",
    ifelse(
      comparison$limma_significant,
      "limma_only",
      "Not_significant"
    )
  )
)

write.csv(
  comparison,
  file.path(de_dir, "edgeR_limma_comparison_all_genes_grouped.csv"),
  row.names = FALSE
)

common <- comparison[
  comparison$edgeR_significant & comparison$limma_significant,
]

write.csv(
  common,
  file.path(de_dir, "common_significant_genes_edgeR_limma_grouped.csv"),
  row.names = FALSE
)

ambiguous_ids <- sprintf("AMBIGUOUS_GROUP_%03d", 1:4)

ambiguous <- comparison[
  comparison$gene_id %in% ambiguous_ids,
]

write.csv(
  ambiguous,
  file.path(de_dir, "ambiguous_groups_edgeR_limma_grouped.csv"),
  row.names = FALSE
)

rho <- cor(
  comparison$edgeR_logFC,
  comparison$limma_logFC,
  method = "spearman",
  use = "complete.obs"
)

cat("\n=== GROUPED edgeR vs limma-voom ===\n")
cat("Groups compared:", nrow(comparison), "\n")
cat("Significant in edgeR:", sum(comparison$edgeR_significant), "\n")
cat("Significant in limma-voom:", sum(comparison$limma_significant), "\n")
cat("Significant in both:", nrow(common), "\n")
cat("edgeR only:", sum(comparison$significance_group == "edgeR_only"), "\n")
cat("limma only:", sum(comparison$significance_group == "limma_only"), "\n")
cat("Not significant in either:",
    sum(comparison$significance_group == "Not_significant"), "\n")

cat("\nDirection concordance:\n")
print(table(comparison$direction_concordant))

cat("\nSpearman correlation of logFC:\n")
print(rho)

cat("\n=== FOUR AMBIGUOUS GROUPS ===\n")
print(ambiguous)

cat("\nGrouped comparison completed successfully.\n")
