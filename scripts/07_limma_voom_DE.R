suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

project <- "/mnt/e/Fungal_Bioinformatics_Portfolio/02_Fungal_Bulk_RNAseq_Transcriptomics"

count_file <- file.path(
  project, "results/counts/tximport_gene_level_length_scaled.rds"
)

metadata_file <- file.path(
  project, "data/metadata/sample_metadata.csv"
)

out_dir <- file.path(project, "results/differential_expression")
fig_dir <- file.path(project, "results/figures")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== LIMMA-VOOM ANALYSIS ===\n")

# 1. Load length-scaled tximport data
txi <- readRDS(count_file)
counts <- txi$counts

# 2. Load metadata
metadata <- read.csv(
  metadata_file,
  row.names = 1,
  check.names = FALSE
)

metadata$condition <- factor(
  metadata$condition,
  levels = c("Sensitive", "Resistant")
)

# 3. Verify sample correspondence
stopifnot(setequal(colnames(counts), rownames(metadata)))

metadata <- metadata[colnames(counts), , drop = FALSE]

stopifnot(identical(colnames(counts), rownames(metadata)))

cat("Samples:\n")
print(metadata)

# 4. Create design matrix
design <- model.matrix(~ condition, data = metadata)

cat("\nDesign matrix:\n")
print(design)

# 5. Create DGEList and filter low-expression genes
y <- DGEList(
  counts = counts,
  group = metadata$condition
)

keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]

cat("\nGenes retained:", nrow(y), "\n")
cat("Genes filtered:", sum(!keep), "\n")

# 6. TMM normalization
y <- calcNormFactors(y, method = "TMM")

cat("\nTMM normalization factors:\n")
print(y$samples[, c("group", "lib.size", "norm.factors")])

# 7. voom transformation
pdf(
  file.path(fig_dir, "limma_voom_mean_variance.pdf"),
  width = 8,
  height = 6
)

v <- voom(y, design, plot = TRUE)

dev.off()

# 8. Linear model
fit <- lmFit(v, design)

# 9. Empirical Bayes moderation
fit <- eBayes(fit, robust = TRUE)

# 10. Extract results: Resistant vs Sensitive
results <- topTable(
  fit,
  coef = "conditionResistant",
  number = Inf,
  sort.by = "P"
)

results$gene_id <- rownames(results)
results <- results[, c("gene_id", setdiff(colnames(results), "gene_id"))]

# 11. Save all results
write.csv(
  results,
  file.path(out_dir, "limma_voom_all_genes.csv"),
  row.names = FALSE
)

# 12. Significant genes
significant <- results[!is.na(results$adj.P.Val) &
                       results$adj.P.Val < 0.05, ]

write.csv(
  significant,
  file.path(out_dir, "limma_voom_significant_genes_FDR05.csv"),
  row.names = FALSE
)

# 13. Upregulated and downregulated genes
up <- significant[significant$logFC > 0, ]
down <- significant[significant$logFC < 0, ]

write.csv(
  up,
  file.path(out_dir, "limma_voom_upregulated_resistant.csv"),
  row.names = FALSE
)

write.csv(
  down,
  file.path(out_dir, "limma_voom_downregulated_resistant.csv"),
  row.names = FALSE
)

# 14. Save model objects
saveRDS(y, file.path(out_dir, "limma_voom_DGEList_TMM.rds"))
saveRDS(v, file.path(out_dir, "limma_voom_object.rds"))
saveRDS(fit, file.path(out_dir, "limma_voom_fit.rds"))

# 15. Summary
cat("\n=== LIMMA-VOOM SUMMARY ===\n")
cat("Genes tested:", nrow(results), "\n")
cat("Significant genes (FDR < 0.05):", nrow(significant), "\n")
cat("Upregulated in Resistant:", nrow(up), "\n")
cat("Downregulated in Resistant:", nrow(down), "\n")

cat("\nTop 10 genes:\n")
print(head(results, 10))

cat("\n=== LIMMA-VOOM COMPLETED ===\n")
