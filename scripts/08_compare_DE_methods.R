project <- "/mnt/e/Fungal_Bioinformatics_Portfolio/02_Fungal_Bulk_RNAseq_Transcriptomics"

de_dir <- file.path(project, "results/differential_expression")

edgeR_file <- file.path(de_dir, "edgeR_all_genes.csv")
limma_file <- file.path(de_dir, "limma_voom_all_genes.csv")

edgeR <- read.csv(edgeR_file, check.names = FALSE)
limma <- read.csv(limma_file, check.names = FALSE)

# Select and rename relevant columns
edgeR <- edgeR[, c("gene_id", "logFC", "PValue", "FDR")]
limma <- limma[, c("gene_id", "logFC", "P.Value", "adj.P.Val")]

colnames(edgeR) <- c(
  "gene_id", "edgeR_logFC", "edgeR_PValue", "edgeR_FDR"
)

colnames(limma) <- c(
  "gene_id", "limma_logFC", "limma_PValue", "limma_FDR"
)

# Merge results by gene ID
comparison <- merge(edgeR, limma, by = "gene_id")

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
    ifelse(comparison$limma_significant, "limma_only", "Not_significant")
  )
)

# Save complete comparison
write.csv(
  comparison,
  file.path(de_dir, "edgeR_limma_comparison_all_genes.csv"),
  row.names = FALSE
)

# Genes significant in both methods
common <- comparison[
  comparison$edgeR_significant & comparison$limma_significant,
]

write.csv(
  common,
  file.path(de_dir, "common_significant_genes_edgeR_limma.csv"),
  row.names = FALSE
)

# Common genes with concordant direction
common_concordant <- common[common$direction_concordant, ]

write.csv(
  common_concordant,
  file.path(de_dir, "common_significant_concordant_genes.csv"),
  row.names = FALSE
)

# Summary
cat("\n=== edgeR vs limma-voom COMPARISON ===\n")
cat("Genes compared:", nrow(comparison), "\n")
cat("Significant in edgeR:", sum(comparison$edgeR_significant), "\n")
cat("Significant in limma-voom:", sum(comparison$limma_significant), "\n")
cat("Significant in both:", nrow(common), "\n")
cat("Significant in both, same direction:", nrow(common_concordant), "\n")
cat("edgeR only:", sum(comparison$significance_group == "edgeR_only"), "\n")
cat("limma-voom only:", sum(comparison$significance_group == "limma_only"), "\n")
cat("Not significant in either:", sum(comparison$significance_group == "Not_significant"), "\n")

cat("\nDirection concordance across all genes:\n")
print(table(comparison$direction_concordant))

cat("\nSpearman correlation of logFC:\n")
print(cor(
  comparison$edgeR_logFC,
  comparison$limma_logFC,
  method = "spearman",
  use = "complete.obs"
))

cat("\nTop 15 genes significant in both methods:\n")
common <- common[order(common$edgeR_FDR, common$limma_FDR), ]
print(head(common, 15))

cat("\n=== COMPARISON COMPLETED ===\n")
