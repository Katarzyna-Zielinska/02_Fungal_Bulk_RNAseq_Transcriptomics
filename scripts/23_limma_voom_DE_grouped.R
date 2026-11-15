# ============================================================
# PROJECT 2: Candida auris Bulk RNA-seq
# Grouped duplicate transcript analysis - limma-voom
# Contrast: Resistant vs Sensitive
# ============================================================

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

project <- "/mnt/e/Fungal_Bioinformatics_Portfolio/02_Fungal_Bulk_RNAseq_Transcriptomics"

input_file <- file.path(
  project,
  "results/counts_keepDuplicates_grouped",
  "tximport_gene_level_length_scaled_grouped.rds"
)

metadata_file <- file.path(
  project,
  "data/metadata/sample_metadata.csv"
)

output_dir <- file.path(
  project,
  "results/differential_expression_grouped"
)

figure_dir <- file.path(
  project,
  "results/figures_grouped"
)

txi <- readRDS(input_file)
counts <- txi$counts

metadata <- read.csv(
  metadata_file,
  row.names = 1,
  check.names = FALSE
)

metadata$condition <- factor(
  metadata$condition,
  levels = c("Sensitive", "Resistant")
)

stopifnot(setequal(colnames(counts), rownames(metadata)))

metadata <- metadata[colnames(counts), , drop = FALSE]

stopifnot(
  identical(colnames(counts), rownames(metadata)),
  nrow(counts) == 5415,
  all(is.finite(counts)),
  all(counts >= 0)
)

design <- model.matrix(~ condition, data = metadata)

if (qr(design)$rank != ncol(design)) {
  stop("Design matrix is not full rank.")
}

cat("=== GROUPED LIMMA-VOOM ANALYSIS ===\n")
cat("\nExperimental design:\n")
print(design)

y <- DGEList(
  counts = counts,
  group = metadata$condition
)

keep <- filterByExpr(y, design = design)

cat("\n=== GENE FILTERING ===\n")
cat("Before:", nrow(y), "\n")
cat("Retained:", sum(keep), "\n")
cat("Filtered:", sum(!keep), "\n")

y <- y[keep, , keep.lib.sizes = FALSE]

y <- calcNormFactors(y, method = "TMM")

cat("\n=== TMM NORMALIZATION ===\n")
print(y$samples[, c("group", "lib.size", "norm.factors")])

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

pdf(
  file.path(figure_dir, "limma_voom_mean_variance_grouped.pdf"),
  width = 8,
  height = 6
)

v <- voom(y, design, plot = TRUE)

dev.off()

fit <- lmFit(v, design)
fit <- eBayes(fit, robust = TRUE)

results <- topTable(
  fit,
  coef = "conditionResistant",
  number = Inf,
  sort.by = "P"
)

results$gene_id <- rownames(results)
results <- results[
  ,
  c("gene_id", setdiff(colnames(results), "gene_id"))
]

significant <- results[
  !is.na(results$adj.P.Val) & results$adj.P.Val < 0.05,
  ,
  drop = FALSE
]

up <- significant[significant$logFC > 0, , drop = FALSE]
down <- significant[significant$logFC < 0, , drop = FALSE]

write.csv(
  results,
  file.path(output_dir, "limma_voom_all_genes_grouped.csv"),
  row.names = FALSE
)

write.csv(
  significant,
  file.path(output_dir, "limma_voom_significant_genes_FDR05_grouped.csv"),
  row.names = FALSE
)

write.csv(
  up,
  file.path(output_dir, "limma_voom_upregulated_resistant_grouped.csv"),
  row.names = FALSE
)

write.csv(
  down,
  file.path(output_dir, "limma_voom_downregulated_resistant_grouped.csv"),
  row.names = FALSE
)

saveRDS(
  y,
  file.path(output_dir, "limma_voom_DGEList_TMM_grouped.rds")
)

saveRDS(
  v,
  file.path(output_dir, "limma_voom_object_grouped.rds")
)

saveRDS(
  fit,
  file.path(output_dir, "limma_voom_fit_grouped.rds")
)

cat("\n=== GROUPED LIMMA-VOOM SUMMARY ===\n")
cat("Genes tested:", nrow(results), "\n")
cat("Significant (FDR < 0.05):", nrow(significant), "\n")
cat("Upregulated in Resistant:", nrow(up), "\n")
cat("Downregulated in Resistant:", nrow(down), "\n")

cat("\n=== FOUR AMBIGUOUS GROUPS ===\n")
ambiguous_ids <- sprintf("AMBIGUOUS_GROUP_%03d", 1:4)

print(
  results[
    results$gene_id %in% ambiguous_ids,
    ,
    drop = FALSE
  ]
)

cat("\nGROUPED LIMMA-VOOM COMPLETED SUCCESSFULLY.\n")
