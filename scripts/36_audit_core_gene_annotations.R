#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
})

cat("\n========== STEP 36: CORE GENE ANNOTATION AUDIT ==========\n")

input_file <- "results/go_enrichment/biological_interpretation/core_genes/core_GO_gene_summary.csv"
out_dir <- "results/go_enrichment/biological_interpretation/annotation_audit"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
    stop("Missing input file: ", input_file)
}

genes <- fread(input_file)

required <- c(
    "gene_id", "Preferred_name", "product",
    "COG_category", "KEGG_Pathway",
    "consensus_logFC", "consensus_direction",
    "GO_term_count"
)

if (!all(required %in% names(genes))) {
    stop(
        "Missing columns: ",
        paste(setdiff(required, names(genes)), collapse = ", ")
    )
}

# --------------------------------------------------
# 1. CLASSIFY PRODUCT ANNOTATION QUALITY
# --------------------------------------------------

genes[
    ,
    product_status := fcase(
        is.na(product) | trimws(product) == "",
        "Missing product annotation",

        grepl(
            "hypothetical protein|uncharacterized protein",
            product,
            ignore.case = TRUE
        ),
        "Hypothetical or uncharacterized protein",

        default = "Named product annotation"
    )
]

genes[
    ,
    preferred_name_status := fifelse(
        is.na(Preferred_name) |
            trimws(Preferred_name) == "",
        "Missing preferred name",
        "Preferred name available"
    )
]

genes[
    ,
    COG_status := fifelse(
        is.na(COG_category) |
            trimws(COG_category) == "",
        "No COG category",
        "COG category available"
    )
]

genes[
    ,
    KEGG_status := fifelse(
        is.na(KEGG_Pathway) |
            trimws(KEGG_Pathway) == "",
        "No KEGG pathway",
        "KEGG pathway available"
    )
]

# This is an annotation evidence classification,
# not experimental validation of protein function.
genes[
    ,
    annotation_evidence := fcase(
        product_status == "Named product annotation" &
            COG_status == "COG category available" &
            KEGG_status == "KEGG pathway available",
        "Named product + COG + KEGG",

        product_status == "Named product annotation" &
            (COG_status == "COG category available" |
             KEGG_status == "KEGG pathway available"),
        "Named product + additional annotation",

        product_status == "Named product annotation",
        "Named product only",

        product_status != "Named product annotation" &
            (COG_status == "COG category available" |
             KEGG_status == "KEGG pathway available"),
        "Hypothetical/missing product + COG/KEGG",

        default = "Limited functional annotation"
    )
]

# --------------------------------------------------
# 2. SAVE COMPLETE AUDIT
# --------------------------------------------------

setorder(
    genes,
    biological_group,
    -GO_term_count,
    gene_id
)

fwrite(
    genes,
    file.path(out_dir, "core_gene_annotation_quality_audit.csv")
)

# --------------------------------------------------
# 3. SUMMARIZE ANNOTATION QUALITY
# --------------------------------------------------

product_summary <- genes[
    ,
    .(genes = .N),
    by = product_status
][order(-genes)]

evidence_summary <- genes[
    ,
    .(genes = .N),
    by = annotation_evidence
][order(-genes)]

cog_summary <- genes[
    ,
    .(genes = .N),
    by = COG_category
][order(-genes)]

group_summary <- genes[
    ,
    .(
        total_genes = .N,
        named_products = sum(
            product_status == "Named product annotation"
        ),
        hypothetical_or_missing = sum(
            product_status != "Named product annotation"
        ),
        with_COG = sum(COG_status == "COG category available"),
        with_KEGG = sum(KEGG_status == "KEGG pathway available")
    ),
    by = biological_group
]

fwrite(
    product_summary,
    file.path(out_dir, "product_annotation_summary.csv")
)

fwrite(
    evidence_summary,
    file.path(out_dir, "annotation_evidence_summary.csv")
)

fwrite(
    cog_summary,
    file.path(out_dir, "COG_category_summary.csv")
)

fwrite(
    group_summary,
    file.path(out_dir, "annotation_quality_by_biological_group.csv")
)

# --------------------------------------------------
# 4. EXPORT GENES REQUIRING FURTHER ANNOTATION
# --------------------------------------------------

needs_annotation <- genes[
    product_status != "Named product annotation"
]

fwrite(
    needs_annotation,
    file.path(out_dir, "genes_requiring_further_annotation.csv")
)

# --------------------------------------------------
# 5. VALIDATION
# --------------------------------------------------

cat("\n========== AUDIT SUMMARY ==========\n")
cat("Total core genes:", nrow(genes), "\n")
cat("Unique gene IDs:", uniqueN(genes$gene_id), "\n")

cat("\nProduct annotation:\n")
print(product_summary)

cat("\nAnnotation evidence:\n")
print(evidence_summary)

cat("\nCOG categories:\n")
print(cog_summary)

cat("\nAnnotation quality by biological group:\n")
print(group_summary)

cat("\nGenes requiring further annotation:",
    nrow(needs_annotation), "\n")

cat("\n========== OUTPUT FILES ==========\n")
cat(out_dir, "\n")

if (nrow(genes) != 43 || uniqueN(genes$gene_id) != 43) {
    stop("Validation failed: expected exactly 43 unique core genes.")
}

cat("\nPASS: Core gene annotation audit completed.\n")
