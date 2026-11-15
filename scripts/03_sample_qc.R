
# ============================================================
# PROJECT 2: Fungal Bulk RNA-seq Transcriptomics
# Script: 03_sample_qc.R
# Purpose: Sample-level quality control
# Organism: Candida auris
# Comparison: Sensitive vs Resistant
# ============================================================

suppressPackageStartupMessages({
  library(tximport)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
})

# 1. Project paths
project_dir <- "/mnt/e/Fungal_Bioinformatics_Portfolio/02_Fungal_Bulk_RNAseq_Transcriptomics"

results_dir <- file.path(project_dir, "results")
qc_dir <- file.path(results_dir, "qc")
figures_dir <- file.path(results_dir, "figures")

dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Load tximport object and metadata
txi <- readRDS(file.path(
  results_dir, "counts", "tximport_gene_level.rds"
))

metadata <- read.csv(
  file.path(project_dir, "data/metadata/sample_metadata.csv"),
  stringsAsFactors = FALSE
)

metadata$condition <- factor(
  metadata$condition,
  levels = c("Sensitive", "Resistant")
)

rownames(metadata) <- metadata$sample

# 3. Verify sample order
stopifnot(
  identical(colnames(txi$counts), rownames(metadata))
)

# 4. Create edgeR DGEList
y <- DGEList(
  counts = txi$counts,
  samples = metadata
)

# 5. Filter lowly expressed genes
design <- model.matrix(~ condition, data = metadata)

keep <- filterByExpr(y, design = design)

cat("\n=== GENE FILTERING ===\n")
cat("Genes before filtering:", nrow(y), "\n")
cat("Genes retained:", sum(keep), "\n")
cat("Genes filtered:", sum(!keep), "\n")

write.table(
  data.frame(
    gene_id = rownames(y),
    retained = keep
  ),
  file.path(qc_dir, "gene_filtering_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

y <- y[keep, , keep.lib.sizes = FALSE]

# 6. TMM normalization
y <- calcNormFactors(y, method = "TMM")

cat("\n=== TMM NORMALIZATION ===\n")
print(y$samples[, c("group", "lib.size", "norm.factors")])

# 7. Calculate log2 CPM
logCPM <- cpm(
  y,
  log = TRUE,
  prior.count = 2,
  normalized.lib.sizes = TRUE
)

# 8. Sample correlation matrix
sample_cor <- cor(
  logCPM,
  method = "pearson",
  use = "pairwise.complete.obs"
)

write.csv(
  sample_cor,
  file.path(qc_dir, "sample_correlation_matrix.csv")
)

# 9. Correlation heatmap
annotation <- data.frame(
  Condition = metadata$condition
)

rownames(annotation) <- metadata$sample

ha <- HeatmapAnnotation(
  df = annotation,
  col = list(
    Condition = c(
      Sensitive = "#377EB8",
      Resistant = "#E41A1C"
    )
  )
)

pdf(
  file.path(figures_dir, "sample_correlation_heatmap.pdf"),
  width = 8,
  height = 7
)

Heatmap(
  sample_cor,
  name = "Pearson\ncorrelation",
  col = colorRamp2(
    c(0.5, 0.75, 1),
    c("#2166AC", "#F7F7F7", "#B2182B")
  ),
  top_annotation = ha,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_title = "Sample correlation",
  row_names_side = "left"
)

dev.off()

# 10. PCA
pca <- prcomp(
  t(logCPM),
  center = TRUE,
  scale. = FALSE
)

percent_var <- 100 * pca$sdev^2 / sum(pca$sdev^2)

pca_data <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  condition = metadata[rownames(pca$x), "condition"]
)

p_pca <- ggplot(
  pca_data,
  aes(x = PC1, y = PC2, color = condition, label = sample)
) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
  scale_color_manual(
    values = c(
      Sensitive = "#377EB8",
      Resistant = "#E41A1C"
    )
  ) +
  labs(
    title = "PCA of Candida auris RNA-seq samples",
    x = paste0("PC1 (", round(percent_var[1], 2), "%)"),
    y = paste0("PC2 (", round(percent_var[2], 2), "%)"),
    color = "Condition"
  ) +
  theme_classic(base_size = 12)

ggsave(
  file.path(figures_dir, "PCA_samples.png"),
  p_pca,
  width = 8,
  height = 6,
  dpi = 300
)

# 11. Library size comparison
library_data <- data.frame(
  sample = colnames(y),
  condition = metadata[colnames(y), "condition"],
  library_size = y$samples$lib.size,
  normalized_library_size =
    y$samples$lib.size * y$samples$norm.factors
)

p_library <- ggplot(
  library_data,
  aes(x = sample, y = library_size, fill = condition)
) +
  geom_col() +
  scale_fill_manual(
    values = c(
      Sensitive = "#377EB8",
      Resistant = "#E41A1C"
    )
  ) +
  labs(
    title = "Library sizes after gene filtering",
    x = "Sample",
    y = "Library size",
    fill = "Condition"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(figures_dir, "library_size_comparison.png"),
  p_library,
  width = 9,
  height = 6,
  dpi = 300
)

# 12. Save normalized data and QC objects
saveRDS(
  y,
  file.path(qc_dir, "dgeList_TMM_filtered.rds")
)

saveRDS(
  logCPM,
  file.path(qc_dir, "logCPM_TMM_filtered.rds")
)

write.csv(
  pca_data,
  file.path(qc_dir, "PCA_coordinates.csv"),
  row.names = FALSE
)

write.csv(
  library_data,
  file.path(qc_dir, "library_size_summary.csv"),
  row.names = FALSE
)

cat("\n=== QC COMPLETED ===\n")
cat("Retained genes:", nrow(y), "\n")
cat("Samples:", ncol(y), "\n")
cat("QC outputs saved in results/qc and results/figures\n")
