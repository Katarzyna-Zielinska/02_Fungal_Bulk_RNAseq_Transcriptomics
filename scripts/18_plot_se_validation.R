# ============================================================
# Project 2: Candida auris bulk RNA-seq
# Validation and visualization of SE event
# Gene: gene-B9J08_004599
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

input_file <- "results/rmats/se_bam_validation_gene-B9J08_004599.csv"
output_dir <- "figures"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(input_file, stringsAsFactors = FALSE)

# Basic validation
required_columns <- c(
  "sample", "group", "rMATS_PSI_JC",
  "rMATS_PSI_JCEC", "BAM_PSI_proxy"
)

missing_columns <- setdiff(required_columns, colnames(df))

if (length(missing_columns) > 0) {
  stop("Brak wymaganych kolumn: ", paste(missing_columns, collapse = ", "))
}

if (nrow(df) != 6) {
  stop("Oczekiwano 6 próbek, znaleziono: ", nrow(df))
}

df$group <- factor(
  df$group,
  levels = c("Sensitive", "Resistant")
)

if (any(is.na(df$group))) {
  stop("Nieprawidłowe przypisanie próbek do grup.")
}

# ------------------------------------------------------------
# 1. rMATS PSI (JC)
# ------------------------------------------------------------

p_rmats <- ggplot(
  df,
  aes(x = group, y = rMATS_PSI_JC)
) +
  geom_point(
    position = position_jitter(width = 0.07, height = 0, seed = 42),
    size = 3
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 95,
    size = 12
  ) +
  scale_y_continuous(
    limits = c(0.90, 1.01),
    breaks = seq(0.90, 1.00, by = 0.02)
  ) +
  labs(
    title = "SE event validation: rMATS PSI",
    subtitle = "gene-B9J08_004599 | Event ID 8 | JC counts",
    x = "Experimental group",
    y = "PSI",
    caption = "Points represent biological replicates (n = 3/group); horizontal marks represent group means."
  ) +
  theme_classic(base_size = 12)

ggsave(
  file.path(output_dir, "SE_gene-B9J08_004599_rMATS_PSI.png"),
  p_rmats, width = 7, height = 5, dpi = 300
)

ggsave(
  file.path(output_dir, "SE_gene-B9J08_004599_rMATS_PSI.pdf"),
  p_rmats, width = 7, height = 5
)

# ------------------------------------------------------------
# 2. BAM-derived PSI proxy
# ------------------------------------------------------------

p_bam <- ggplot(
  df,
  aes(x = group, y = BAM_PSI_proxy)
) +
  geom_point(
    position = position_jitter(width = 0.07, height = 0, seed = 42),
    size = 3
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 95,
    size = 12
  ) +
  scale_y_continuous(
    limits = c(0.90, 1.01),
    breaks = seq(0.90, 1.00, by = 0.02)
  ) +
  labs(
    title = "SE junction audit: BAM PSI proxy",
    subtitle = "gene-B9J08_004599 | Independent junction-count audit",
    x = "Experimental group",
    y = "Approximate PSI proxy",
    caption = "Exploratory proxy based on unique junction read names; not a formal PSI estimate."
  ) +
  theme_classic(base_size = 12)

ggsave(
  file.path(output_dir, "SE_gene-B9J08_004599_BAM_PSI_proxy.png"),
  p_bam, width = 7, height = 5, dpi = 300
)

ggsave(
  file.path(output_dir, "SE_gene-B9J08_004599_BAM_PSI_proxy.pdf"),
  p_bam, width = 7, height = 5
)

cat("Wykresy zapisane w katalogu:", output_dir, "\n")
cat("1. SE_gene-B9J08_004599_rMATS_PSI.png/pdf\n")
cat("2. SE_gene-B9J08_004599_BAM_PSI_proxy.png/pdf\n")
cat("\nLiczba próbek:", nrow(df), "\n")
cat("Średnie rMATS PSI:\n")
print(aggregate(rMATS_PSI_JC ~ group, data = df, FUN = mean))
cat("\nŚrednie BAM PSI proxy:\n")
print(aggregate(BAM_PSI_proxy ~ group, data = df, FUN = mean))
cat("\n========== SE VISUALIZATION COMPLETE ==========\n")
