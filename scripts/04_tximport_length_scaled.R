
# ============================================================
# PROJECT 2: Candida auris Bulk RNA-seq
# Script: 04_tximport_length_scaled.R
# Purpose: Length-scaled transcript-to-gene import
# ============================================================

suppressPackageStartupMessages({
  library(tximport)
})

project <- "/mnt/e/Fungal_Bioinformatics_Portfolio/02_Fungal_Bulk_RNAseq_Transcriptomics"

samples <- c(
  "SRR10059687", "SRR10059688", "SRR10059689",
  "SRR10059690", "SRR10059691", "SRR10059692"
)

conditions <- factor(
  c("Sensitive", "Sensitive", "Sensitive",
    "Resistant", "Resistant", "Resistant"),
  levels = c("Sensitive", "Resistant")
)

groups <- ifelse(conditions == "Sensitive", "sensitive", "resistant")

files <- file.path(
  project, "results/quantification",
  groups, samples, "quant.sf"
)

names(files) <- samples

stopifnot(all(file.exists(files)))

tx2gene <- read.delim(
  file.path(project, "data/reference/tx2gene.tsv"),
  header = FALSE,
  col.names = c("TXNAME", "GENEID"),
  stringsAsFactors = FALSE
)

txi_scaled <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM"
)

stopifnot(identical(colnames(txi_scaled$counts), samples))

cat("\n=== LENGTH-SCALED TXIMPORT ===\n")
cat("Genes:", nrow(txi_scaled$counts), "\n")
cat("Samples:", ncol(txi_scaled$counts), "\n")
cat("NA:", sum(is.na(txi_scaled$counts)), "\n")
cat("Infinite:", sum(!is.finite(txi_scaled$counts)), "\n")

saveRDS(
  txi_scaled,
  file.path(project, "results/counts/tximport_gene_level_length_scaled.rds")
)

cat("\nLength-scaled tximport completed successfully.\n")
