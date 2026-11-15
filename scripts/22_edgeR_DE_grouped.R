# ============================================================
# PROJECT 2: Candida auris Bulk RNA-seq
# Grouped duplicate transcript analysis - edgeR
# ============================================================

suppressPackageStartupMessages({
  library(edgeR)
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

if (dir.exists(output_dir) &&
    length(list.files(output_dir, all.files = TRUE, no.. = TRUE)) > 0) {
  stop("Output directory is not empty. Refusing to overwrite results.")
}

txi <- readRDS(input_file)

metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE
)

metadata$condition <- factor(
  metadata$condition,
  levels = c("Sensitive", "Resistant")
)

rownames(metadata) <- metadata$sample

stopifnot(
  identical(colnames(txi$counts), rownames(metadata)),
  nrow(txi$counts) == 5415,
  all(table(metadata$condition) == c(3, 3)),
  all(is.finite(txi$counts)),
  all(txi$counts >= 0)
)

y <- DGEList(
  counts = txi$counts,
  group = metadata$condition
)

design <- model.matrix(~ condition, data = metadata)
rownames(design) <- rownames(metadata)

if (qr(design)$rank != ncol(design)) {
  stop("Design matrix is not full rank.")
}

cat("\n=== EXPERIMENTAL DESIGN ===\n")
print(design)

keep <- filterByExpr(y, design = design)

cat("\n=== GENE FILTERING ===\n")
cat("Before:", nrow(y), "\n")
cat("Retained:", sum(keep), "\n")
cat("Filtered:", sum(!keep), "\n")

y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")

cat("\n=== TMM NORMALIZATION ===\n")
print(y$samples[, c("group", "lib.size", "norm.factors")])

y <- estimateDisp(
  y,
  design,
  robust = TRUE
)

fit <- glmQLFit(
  y,
  design,
  robust = TRUE
)

contrast <- makeContrasts(
  conditionResistant,
  levels = design
)

qlf <- glmQLFTest(
  fit,
  contrast = contrast
)

results <- topTags(
  qlf,
  n = Inf,
  sort.by = "PValue"
)$table

results$gene_id <- rownames(results)

results <- results[
  ,
  c("gene_id", "logFC", "logCPM", "F", "PValue", "FDR")
]

significant <- subset(results, FDR < 0.05)
up <- subset(significant, logFC > 0)
down <- subset(significant, logFC < 0)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  results,
  file.path(output_dir, "edgeR_all_genes_grouped.csv"),
  row.names = FALSE
)

write.csv(
  significant,
  file.path(output_dir, "edgeR_significant_genes_FDR05_grouped.csv"),
  row.names = FALSE
)

write.csv(
  up,
  file.path(output_dir, "edgeR_upregulated_resistant_grouped.csv"),
  row.names = FALSE
)

write.csv(
  down,
  file.path(output_dir, "edgeR_downregulated_resistant_grouped.csv"),
  row.names = FALSE
)

write.table(
  data.frame(
    gene_id = rownames(txi$counts),
    retained = rownames(txi$counts) %in% rownames(y)
  ),
  file.path(output_dir, "edgeR_gene_filtering_grouped.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(y, file.path(output_dir, "edgeR_DGEList_TMM_grouped.rds"))
saveRDS(fit, file.path(output_dir, "edgeR_QL_fit_grouped.rds"))
saveRDS(qlf, file.path(output_dir, "edgeR_QL_test_grouped.rds"))

write.csv(
  y$samples,
  file.path(output_dir, "edgeR_sample_normalization_grouped.csv")
)

cat("\n=== GROUPED edgeR SUMMARY ===\n")
cat("Genes tested:", nrow(results), "\n")
cat("Significant (FDR < 0.05):", nrow(significant), "\n")
cat("Upregulated in Resistant:", nrow(up), "\n")
cat("Downregulated in Resistant:", nrow(down), "\n")

cat("\n=== FOUR AMBIGUOUS GROUPS ===\n")
print(
  results[
    results$gene_id %in% sprintf("AMBIGUOUS_GROUP_%03d", 1:4),
    ,
    drop = FALSE
  ]
)

cat("\nGROUPED edgeR COMPLETED SUCCESSFULLY.\n")
