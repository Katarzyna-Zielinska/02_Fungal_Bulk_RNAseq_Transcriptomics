# ============================================================
# PROJECT 2: Candida auris Bulk RNA-seq Transcriptomics
# Script: 05_edgeR_DE.R
# Purpose: Differential gene expression with edgeR
# Contrast: Resistant vs Sensitive
# Method: TMM + quasi-likelihood
# ============================================================

suppressPackageStartupMessages({
  library(edgeR)
})

# 1. Paths
project <- "/mnt/e/Fungal_Bioinformatics_Portfolio/02_Fungal_Bulk_RNAseq_Transcriptomics"

counts_dir <- file.path(project, "results/counts")
de_dir <- file.path(project, "results/differential_expression")
qc_dir <- file.path(project, "results/qc")

dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Load length-scaled tximport object
txi <- readRDS(
  file.path(counts_dir, "tximport_gene_level_length_scaled.rds")
)

# 3. Load metadata
metadata <- read.csv(
  file.path(project, "data/metadata/sample_metadata.csv"),
  stringsAsFactors = FALSE
)

metadata$condition <- factor(
  metadata$condition,
  levels = c("Sensitive", "Resistant")
)

rownames(metadata) <- metadata$sample

# 4. Verify sample identities and order
stopifnot(
  identical(colnames(txi$counts), rownames(metadata))
)

stopifnot(
  all(table(metadata$condition) == c(3, 3))
)

# 5. Create DGEList with explicitly assigned groups
y <- DGEList(
  counts = txi$counts,
  group = metadata$condition
)

# 6. Design matrix
design <- model.matrix(~ condition, data = metadata)
rownames(design) <- rownames(metadata)

cat("\n=== EXPERIMENTAL DESIGN ===\n")
print(design)

if (qr(design)$rank != ncol(design)) {
  stop("ERROR: Design matrix is not full rank.")
}

# 7. Filter lowly expressed genes
keep <- filterByExpr(y, design = design)

cat("\n=== GENE FILTERING ===\n")
cat("Genes before filtering:", nrow(y), "\n")
cat("Genes retained:", sum(keep), "\n")
cat("Genes filtered:", sum(!keep), "\n")

y <- y[keep, , keep.lib.sizes = FALSE]

# 8. TMM normalization
y <- calcNormFactors(y, method = "TMM")

cat("\n=== TMM NORMALIZATION ===\n")
print(y$samples[, c("group", "lib.size", "norm.factors")])

# 9. Estimate dispersion
y <- estimateDisp(
  y,
  design,
  robust = TRUE
)

cat("\n=== DISPERSION ESTIMATION ===\n")
cat("Common dispersion:", y$common.dispersion, "\n")
cat("Tagwise dispersions:", length(y$tagwise.dispersion), "\n")

# 10. Quasi-likelihood model
fit <- glmQLFit(
  y,
  design,
  robust = TRUE
)

# 11. Resistant vs Sensitive contrast
contrast <- makeContrasts(
  conditionResistant,
  levels = design
)

qlf <- glmQLFTest(
  fit,
  contrast = contrast
)

# 12. Extract all results
de_results <- topTags(
  qlf,
  n = Inf,
  sort.by = "PValue"
)$table

de_results$gene_id <- rownames(de_results)

de_results <- de_results[
  ,
  c("gene_id", "logFC", "logCPM", "F", "PValue", "FDR")
]

# 13. Significant genes
significant <- subset(
  de_results,
  FDR < 0.05
)

upregulated <- subset(
  significant,
  logFC > 0
)

downregulated <- subset(
  significant,
  logFC < 0
)

# 14. Save results
write.csv(
  de_results,
  file.path(de_dir, "edgeR_all_genes.csv"),
  row.names = FALSE
)

write.csv(
  significant,
  file.path(de_dir, "edgeR_significant_genes_FDR05.csv"),
  row.names = FALSE
)

write.csv(
  upregulated,
  file.path(de_dir, "edgeR_upregulated_resistant.csv"),
  row.names = FALSE
)

write.csv(
  downregulated,
  file.path(de_dir, "edgeR_downregulated_resistant.csv"),
  row.names = FALSE
)

saveRDS(
  y,
  file.path(de_dir, "edgeR_DGEList_TMM.rds")
)

saveRDS(
  fit,
  file.path(de_dir, "edgeR_QL_fit.rds")
)

saveRDS(
  qlf,
  file.path(de_dir, "edgeR_QL_test.rds")
)

write.csv(
  y$samples,
  file.path(de_dir, "edgeR_sample_normalization.csv")
)

write.table(
  data.frame(
    gene_id = rownames(txi$counts),
    retained = rownames(txi$counts) %in% rownames(y)
  ),
  file.path(de_dir, "edgeR_gene_filtering.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 15. Summary
cat("\n=== edgeR DIFFERENTIAL EXPRESSION SUMMARY ===\n")
cat("Genes tested:", nrow(de_results), "\n")
cat("Significant genes (FDR < 0.05):", nrow(significant), "\n")
cat("Upregulated in Resistant:", nrow(upregulated), "\n")
cat("Downregulated in Resistant:", nrow(downregulated), "\n")

cat("\n=== TOP 10 GENES ===\n")
print(head(de_results, 10))

cat("\n=== edgeR ANALYSIS COMPLETED ===\n")
