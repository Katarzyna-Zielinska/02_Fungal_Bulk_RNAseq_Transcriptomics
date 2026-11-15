
suppressPackageStartupMessages({
  library(edgeR)
})

project <- "/mnt/e/Fungal_Bioinformatics_Portfolio/02_Fungal_Bulk_RNAseq_Transcriptomics"

de_dir <- file.path(project, "results/differential_expression")
fig_dir <- file.path(project, "results/figures")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

y <- readRDS(file.path(de_dir, "edgeR_DGEList_TMM.rds"))
fit <- readRDS(file.path(de_dir, "edgeR_QL_fit.rds"))

# 1. Biological coefficient of variation
pdf(file.path(fig_dir, "edgeR_BCV_plot.pdf"), width = 8, height = 6)
plotBCV(y, main = "edgeR: Biological coefficient of variation")
dev.off()

# 2. Mean-variance trend
pdf(file.path(fig_dir, "edgeR_QL_dispersion.pdf"), width = 8, height = 6)
plotQLDisp(fit)
dev.off()

# 3. MDS sample distances
pdf(file.path(fig_dir, "edgeR_MDS_plot.pdf"), width = 8, height = 6)
plotMDS(
  y,
  labels = colnames(y),
  col = c(rep("#377EB8", 3), rep("#E41A1C", 3)),
  main = "edgeR: MDS of RNA-seq samples"
)
legend(
  "topright",
  legend = c("Sensitive", "Resistant"),
  col = c("#377EB8", "#E41A1C"),
  pch = 16
)
dev.off()

cat("\n=== edgeR DIAGNOSTICS COMPLETED ===\n")
cat("BCV, QL dispersion and MDS plots saved.\n")
