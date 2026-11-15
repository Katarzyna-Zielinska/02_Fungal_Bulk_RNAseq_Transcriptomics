
# ============================================================
# PROJECT 2: Fungal Bulk RNA-seq Transcriptomics
# Script: 02_tximport.R
# Purpose: Import Salmon quantifications using tximport
# Organism: Candida auris
# Comparison: Sensitive vs Resistant
# ============================================================

suppressPackageStartupMessages({
  library(tximport)
})

# 1. Project directory
project_dir <- "/mnt/e/Fungal_Bioinformatics_Portfolio/02_Fungal_Bulk_RNAseq_Transcriptomics"

# 2. Sample metadata
samples <- data.frame(
  sample = c(
    "SRR10059687",
    "SRR10059688",
    "SRR10059689",
    "SRR10059690",
    "SRR10059691",
    "SRR10059692"
  ),
  condition = factor(
    c(
      "Sensitive",
      "Sensitive",
      "Sensitive",
      "Resistant",
      "Resistant",
      "Resistant"
    ),
    levels = c("Sensitive", "Resistant")
  )
)

rownames(samples) <- samples$sample

# 3. Locate Salmon quant.sf files
quant_files <- file.path(
  project_dir,
  "results",
  "quantification",
  ifelse(samples$condition == "Sensitive", "sensitive", "resistant"),
  samples$sample,
  "quant.sf"
)

names(quant_files) <- samples$sample

if (!all(file.exists(quant_files))) {
  stop("ERROR: One or more Salmon quant.sf files are missing.")
}

# 4. Load transcript-to-gene mapping
tx2gene <- read.delim(
  file.path(project_dir, "data/reference/tx2gene.tsv"),
  header = FALSE,
  stringsAsFactors = FALSE
)

colnames(tx2gene) <- c("TXNAME", "GENEID")

# 5. Import Salmon quantifications
txi <- tximport(
  quant_files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "no",
  ignoreTxVersion = FALSE
)

# 6. Check imported data
stopifnot(identical(colnames(txi$counts), samples$sample))
stopifnot(identical(colnames(txi$abundance), samples$sample))
stopifnot(identical(colnames(txi$length), samples$sample))

cat("\n=== TXIMPORT SUMMARY ===\n")
cat("Number of genes:", nrow(txi$counts), "\n")
cat("Number of samples:", ncol(txi$counts), "\n")
cat("Sensitive samples:", sum(samples$condition == "Sensitive"), "\n")
cat("Resistant samples:", sum(samples$condition == "Resistant"), "\n")

cat("\n=== COUNT MATRIX DIMENSIONS ===\n")
print(dim(txi$counts))

# 7. Create output directories
output_dir <- file.path(project_dir, "results/counts")
metadata_dir <- file.path(project_dir, "data/metadata")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

# 8. Save imported data
saveRDS(
  txi,
  file.path(output_dir, "tximport_gene_level.rds")
)

write.csv(
  txi$counts,
  file.path(output_dir, "gene_counts.csv"),
  quote = FALSE
)

write.csv(
  txi$abundance,
  file.path(output_dir, "gene_TPM.csv"),
  quote = FALSE
)

write.csv(
  txi$length,
  file.path(output_dir, "gene_effective_lengths.csv"),
  quote = FALSE
)

write.csv(
  samples,
  file.path(metadata_dir, "sample_metadata.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("\n=== FILES SAVED ===\n")
cat("results/counts/tximport_gene_level.rds\n")
cat("results/counts/gene_counts.csv\n")
cat("results/counts/gene_TPM.csv\n")
cat("results/counts/gene_effective_lengths.csv\n")
cat("data/metadata/sample_metadata.csv\n")

cat("\nTXIMPORT COMPLETED SUCCESSFULLY.\n")
