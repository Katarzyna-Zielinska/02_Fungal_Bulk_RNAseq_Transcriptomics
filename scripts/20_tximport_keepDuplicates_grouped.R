# ============================================================
# PROJECT 2: Candida auris Bulk RNA-seq
# Script: 20_tximport_keepDuplicates_grouped.R
# Purpose: Group exact duplicate transcripts before gene-level
#          differential expression analysis
# ============================================================

suppressPackageStartupMessages({
  library(tximport)
})

project_dir <- "/mnt/e/Fungal_Bioinformatics_Portfolio/02_Fungal_Bulk_RNAseq_Transcriptomics"

samples <- data.frame(
  sample = c(
    "SRR10059687", "SRR10059688", "SRR10059689",
    "SRR10059690", "SRR10059691", "SRR10059692"
  ),
  condition = factor(
    c(
      "Sensitive", "Sensitive", "Sensitive",
      "Resistant", "Resistant", "Resistant"
    ),
    levels = c("Sensitive", "Resistant")
  )
)

rownames(samples) <- samples$sample

quant_files <- file.path(
  project_dir,
  "results",
  "quantification_keepDuplicates",
  ifelse(samples$condition == "Sensitive", "sensitive", "resistant"),
  samples$sample,
  "quant.sf"
)

names(quant_files) <- samples$sample

if (!all(file.exists(quant_files))) {
  stop("ERROR: Missing quant.sf files:\n",
       paste(quant_files[!file.exists(quant_files)], collapse = "\n"))
}

tx2gene <- read.delim(
  file.path(project_dir, "data/reference/tx2gene_keepDuplicates_grouped.tsv"),
  header = FALSE,
  col.names = c("TXNAME", "GENEID"),
  stringsAsFactors = FALSE
)

# Check that every Salmon transcript has exactly one mapping.
quant_names <- lapply(quant_files, function(f) {
  read.delim(f, nrows = 1, check.names = FALSE)$Name
})

for (sample_name in names(quant_files)) {
  q <- read.delim(quant_files[[sample_name]], check.names = FALSE)

  missing_mapping <- setdiff(q$Name, tx2gene$TXNAME)
  extra_mapping <- setdiff(tx2gene$TXNAME, q$Name)

  if (length(missing_mapping) > 0 || length(extra_mapping) > 0) {
    stop("Transcript ID mismatch in sample: ", sample_name,
         "\nMissing mappings: ", length(missing_mapping),
         "\nExtra mappings: ", length(extra_mapping))
  }

  if (anyDuplicated(q$Name)) {
    stop("Duplicated transcript names in quant.sf: ", sample_name)
  }
}

cat("\n=== INPUT CHECK ===\n")
cat("Samples:", length(quant_files), "\n")
cat("Mapped transcripts:", nrow(tx2gene), "\n")
cat("Unique gene groups:", length(unique(tx2gene$GENEID)), "\n")

txi <- tximport(
  quant_files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "no",
  ignoreTxVersion = FALSE
)

stopifnot(identical(colnames(txi$counts), samples$sample))
stopifnot(identical(colnames(txi$abundance), samples$sample))
stopifnot(identical(colnames(txi$length), samples$sample))

expected_genes <- length(unique(tx2gene$GENEID))

if (nrow(txi$counts) != expected_genes) {
  stop("Unexpected number of gene groups after tximport.")
}

output_dir <- file.path(
  project_dir,
  "results/counts_keepDuplicates_grouped"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(txi, file.path(output_dir, "tximport_gene_level_grouped.rds"))

write.csv(
  txi$counts,
  file.path(output_dir, "gene_counts_grouped.csv"),
  quote = FALSE
)

write.csv(
  txi$abundance,
  file.path(output_dir, "gene_TPM_grouped.csv"),
  quote = FALSE
)

write.csv(
  txi$length,
  file.path(output_dir, "gene_effective_lengths_grouped.csv"),
  quote = FALSE
)

write.csv(
  samples,
  file.path(output_dir, "sample_metadata.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("\n=== TXIMPORT SUMMARY ===\n")
cat("Gene groups:", nrow(txi$counts), "\n")
cat("Samples:", ncol(txi$counts), "\n")
cat("Count matrix dimensions:\n")
print(dim(txi$counts))

cat("\n=== AMBIGUOUS GROUP COUNTS ===\n")
ambiguous_ids <- grep("^AMBIGUOUS_GROUP_", rownames(txi$counts), value = TRUE)
print(txi$counts[ambiguous_ids, , drop = FALSE])

cat("\n=== OUTPUT FILES ===\n")
cat(output_dir, "\n")
cat("TXIMPORT GROUPED COMPLETED SUCCESSFULLY.\n")
