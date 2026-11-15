# ============================================================
# PROJECT 2: Candida auris Bulk RNA-seq
# Script: 21_tximport_keepDuplicates_length_scaled.R
# Purpose: Length-scaled tximport for grouped duplicate loci
# Method: Salmon + tximport + lengthScaledTPM
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
  project,
  "results/quantification_keepDuplicates",
  groups,
  samples,
  "quant.sf"
)

names(files) <- samples

if (!all(file.exists(files))) {
  stop(
    "Missing quant.sf files:\n",
    paste(files[!file.exists(files)], collapse = "\n")
  )
}

tx2gene <- read.delim(
  file.path(project, "data/reference/tx2gene_keepDuplicates_grouped.tsv"),
  header = FALSE,
  col.names = c("TXNAME", "GENEID"),
  stringsAsFactors = FALSE
)

if (nrow(tx2gene) != 5419 || length(unique(tx2gene$GENEID)) != 5415) {
  stop("Unexpected tx2gene mapping dimensions.")
}

txi_scaled <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM",
  ignoreTxVersion = FALSE
)

stopifnot(identical(colnames(txi_scaled$counts), samples))
stopifnot(nrow(txi_scaled$counts) == 5415)
stopifnot(sum(is.na(txi_scaled$counts)) == 0)
stopifnot(all(is.finite(txi_scaled$counts)))

output_dir <- file.path(
  project,
  "results/counts_keepDuplicates_grouped"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(
  txi_scaled,
  file.path(output_dir, "tximport_gene_level_length_scaled_grouped.rds")
)

write.csv(
  txi_scaled$counts,
  file.path(output_dir, "gene_counts_length_scaled_grouped.csv"),
  quote = FALSE
)

cat("\n=== GROUPED LENGTH-SCALED TXIMPORT ===\n")
cat("Transcripts mapped:", nrow(tx2gene), "\n")
cat("Gene groups:", nrow(txi_scaled$counts), "\n")
cat("Samples:", ncol(txi_scaled$counts), "\n")
cat("NA values:", sum(is.na(txi_scaled$counts)), "\n")
cat("Infinite values:", sum(!is.finite(txi_scaled$counts)), "\n")
cat("Method: lengthScaledTPM\n")

cat("\n=== AMBIGUOUS GROUPS ===\n")
ambiguous_ids <- grep(
  "^AMBIGUOUS_GROUP_",
  rownames(txi_scaled$counts),
  value = TRUE
)
print(txi_scaled$counts[ambiguous_ids, , drop = FALSE])

cat("\nSaved files:\n")
cat("tximport_gene_level_length_scaled_grouped.rds\n")
cat("gene_counts_length_scaled_grouped.csv\n")
cat("\nGROUPED LENGTH-SCALED TXIMPORT COMPLETED.\n")
