#!/usr/bin/env Rscript

# ============================================================
# Project 2: Candida auris Bulk RNA-seq Transcriptomics
# Step 28: Integration of shared pathway leading-edge genes
# Grouped analysis: FGSEA + eggNOG + edgeR + limma-voom
# ============================================================

suppressPackageStartupMessages({
    library(data.table)
})

project_dir <- normalizePath(".", mustWork = TRUE)

fgsea_file <- file.path(
    project_dir,
    "results/pathway_enrichment_grouped/fgsea/fgsea_all_pathways_grouped.csv"
)

annotation_file <- file.path(
    project_dir,
    "results/integrated_analysis/integrated_DE_functional_annotation.csv"
)

edgeR_file <- file.path(
    project_dir,
    "results/differential_expression_grouped/edgeR_all_genes_grouped.csv"
)

limma_file <- file.path(
    project_dir,
    "results/differential_expression_grouped/limma_voom_all_genes_grouped.csv"
)

output_dir <- file.path(
    project_dir,
    "results/pathway_enrichment_grouped/biological_interpretation"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

FDR_THRESHOLD <- 0.05

shared_pathways <- c(
    "ko03008",
    "ko03020",
    "ko04066",
    "ko00710",
    "ko05230"
)

# ------------------------------------------------------------
# 1. Check input files
# ------------------------------------------------------------

cat("\n=== CHECKING INPUT FILES ===\n")

input_files <- c(
    fgsea_file,
    annotation_file,
    edgeR_file,
    limma_file
)

if (!all(file.exists(input_files))) {
    stop(
        "Missing input files:\n",
        paste(input_files[!file.exists(input_files)], collapse = "\n")
    )
}

cat("All four input files found.\n")

# ------------------------------------------------------------
# 2. Load data
# ------------------------------------------------------------

cat("\n=== LOADING DATA ===\n")

fgsea <- fread(fgsea_file)
annotation <- fread(annotation_file)
edgeR <- fread(edgeR_file)
limma <- fread(limma_file)

stopifnot(
    all(c("pathway", "NES", "padj", "size", "leadingEdge") %in% names(fgsea)),
    all(c(
        "gene_id", "gene_name", "product", "Description",
        "Preferred_name", "GOs", "COG_category",
        "KEGG_ko", "KEGG_Pathway"
    ) %in% names(annotation)),
    all(c("gene_id", "logFC", "logCPM", "PValue", "FDR") %in% names(edgeR)),
    all(c("gene_id", "logFC", "AveExpr", "P.Value", "adj.P.Val") %in% names(limma))
)

# ------------------------------------------------------------
# 3. Select the five shared pathways
# ------------------------------------------------------------

cat("\n=== SELECTING SHARED PATHWAYS ===\n")

fgsea_selected <- fgsea[pathway %in% shared_pathways]

if (nrow(fgsea_selected) != length(shared_pathways)) {
    stop("Expected exactly five shared pathways in FGSEA results.")
}

cat("Selected pathways:", nrow(fgsea_selected), "\n")

# ------------------------------------------------------------
# 4. Expand leading-edge genes
# ------------------------------------------------------------

cat("\n=== EXPANDING LEADING-EDGE GENES ===\n")

leading_edge <- fgsea_selected[
    ,
    .(gene_id = unlist(strsplit(leadingEdge, ";", fixed = TRUE))),
    by = .(pathway, NES, fgsea_padj = padj, pathway_size = size)
]

leading_edge <- unique(leading_edge)

leading_edge[, fgsea_direction := ifelse(
    NES > 0, "Resistant",
    ifelse(NES < 0, "Sensitive", "No direction")
)]

cat("Unique gene-pathway records:", nrow(leading_edge), "\n")
cat("Unique genes:", uniqueN(leading_edge$gene_id), "\n")

# ------------------------------------------------------------
# 5. Prepare annotation data
# ------------------------------------------------------------

cat("\n=== PREPARING EGGNOG ANNOTATIONS ===\n")

annotation_selected <- annotation[
    ,
    .(
        gene_id,
        gene_name,
        product,
        Description,
        Preferred_name,
        COG_category,
        GOs,
        KEGG_ko,
        KEGG_Pathway,
        KEGG_Module,
        has_functional_annotation
    )
]

annotation_selected <- unique(annotation_selected, by = "gene_id")

# ------------------------------------------------------------
# 6. Prepare grouped edgeR and limma-voom results
# ------------------------------------------------------------

cat("\n=== PREPARING GROUPED DIFFERENTIAL EXPRESSION ===\n")

edgeR_selected <- edgeR[
    ,
    .(
        gene_id,
        edgeR_logFC = logFC,
        edgeR_logCPM = logCPM,
        edgeR_PValue = PValue,
        edgeR_FDR = FDR
    )
]

limma_selected <- limma[
    ,
    .(
        gene_id,
        limma_logFC = logFC,
        limma_AveExpr = AveExpr,
        limma_PValue = P.Value,
        limma_FDR = adj.P.Val
    )
]

# ------------------------------------------------------------
# 7. Integrate all data
# ------------------------------------------------------------

cat("\n=== INTEGRATING DATASETS ===\n")

integrated <- merge(
    leading_edge,
    annotation_selected,
    by = "gene_id",
    all.x = TRUE
)

integrated <- merge(
    integrated,
    edgeR_selected,
    by = "gene_id",
    all.x = TRUE
)

integrated <- merge(
    integrated,
    limma_selected,
    by = "gene_id",
    all.x = TRUE
)

# ------------------------------------------------------------
# 8. Classify differential expression
# ------------------------------------------------------------

cat("\n=== CLASSIFYING DIFFERENTIAL EXPRESSION ===\n")

integrated[, edgeR_significant := !is.na(edgeR_FDR) & edgeR_FDR < FDR_THRESHOLD]
integrated[, limma_significant := !is.na(limma_FDR) & limma_FDR < FDR_THRESHOLD]

integrated[, DE_class := fcase(
    edgeR_significant & limma_significant, "common_significant",
    edgeR_significant & !limma_significant, "edgeR_only",
    !edgeR_significant & limma_significant, "limma_only",
    default = "not_significant"
)]

integrated[, direction_concordant := sign(edgeR_logFC) == sign(limma_logFC)]

integrated[, consensus_direction := fcase(
    edgeR_logFC > 0 & limma_logFC > 0, "up_in_resistant",
    edgeR_logFC < 0 & limma_logFC < 0, "down_in_resistant",
    default = "discordant_or_unavailable"
)]

integrated[, annotation_status := ifelse(
    is.na(has_functional_annotation),
    "annotation_not_found",
    ifelse(has_functional_annotation, "annotated", "no_functional_annotation")
)]

# ------------------------------------------------------------
# 9. Order and export detailed table
# ------------------------------------------------------------

setorder(integrated, pathway, fgsea_padj, edgeR_FDR)

detail_file <- file.path(
    output_dir,
    "shared_pathways_leading_edge_gene_integration.csv"
)

fwrite(integrated, detail_file)

# ------------------------------------------------------------
# 10. Create pathway-level summary
# ------------------------------------------------------------

cat("\n=== CREATING PATHWAY SUMMARY ===\n")

summary_table <- integrated[
    ,
    .(
        leading_edge_genes = uniqueN(gene_id),
        annotated_genes = uniqueN(gene_id[annotation_status == "annotated"]),
        genes_common_DE = uniqueN(gene_id[DE_class == "common_significant"]),
        genes_edgeR_only = uniqueN(gene_id[DE_class == "edgeR_only"]),
        genes_limma_only = uniqueN(gene_id[DE_class == "limma_only"]),
        genes_not_significant = uniqueN(gene_id[DE_class == "not_significant"]),
        genes_up_in_resistant = uniqueN(gene_id[consensus_direction == "up_in_resistant"]),
        genes_down_in_resistant = uniqueN(gene_id[consensus_direction == "down_in_resistant"])
    ),
    by = .(pathway, NES, fgsea_padj, fgsea_direction)
]

setorder(summary_table, fgsea_padj)

summary_file <- file.path(
    output_dir,
    "shared_pathways_leading_edge_summary.csv"
)

fwrite(summary_table, summary_file)

# ------------------------------------------------------------
# 11. Export genes significant in both DE methods
# ------------------------------------------------------------

common_DE <- integrated[
    DE_class == "common_significant"
]

common_file <- file.path(
    output_dir,
    "shared_pathways_common_DE_genes.csv"
)

fwrite(common_DE, common_file)

# ------------------------------------------------------------
# 12. Validation
# ------------------------------------------------------------

cat("\n=== FINAL VALIDATION ===\n")

stopifnot(
    all(integrated$pathway %in% shared_pathways),
    !anyDuplicated(integrated[, .(pathway, gene_id)]),
    all(integrated$DE_class %in% c(
        "common_significant",
        "edgeR_only",
        "limma_only",
        "not_significant"
    )),
    nrow(summary_table) == 5
)

cat("Pathways:", uniqueN(integrated$pathway), "\n")
cat("Gene-pathway records:", nrow(integrated), "\n")
cat("Unique genes:", uniqueN(integrated$gene_id), "\n")
cat("Annotated genes:", uniqueN(integrated[annotation_status == "annotated", gene_id]), "\n")
cat("Common DE genes:", uniqueN(common_DE$gene_id), "\n")

cat("\n=== PATHWAY SUMMARY ===\n")
print(summary_table)

cat("\n============================================\n")
cat("STEP 28 COMPLETED SUCCESSFULLY\n")
cat("============================================\n")
cat("Detailed table:\n", detail_file, "\n")
cat("Pathway summary:\n", summary_file, "\n")
cat("Common DE genes:\n", common_file, "\n")
