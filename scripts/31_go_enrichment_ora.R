#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(GO.db)
    library(AnnotationDbi)
})

cat("\n========== PROJECT 2: GO ENRICHMENT ORA ==========\n")

# --------------------------------------------------
# 1. INPUT FILES
# --------------------------------------------------

annotation_file <- "results/integrated_analysis/integrated_DE_functional_annotation.csv"

tier1_file <- "results/pathway_enrichment_grouped/biological_interpretation/tier1_50_genes_functional_characterization.csv"

filter_file <- "results/differential_expression_grouped/edgeR_gene_filtering_grouped.tsv"

output_dir <- "results/go_enrichment/ORA"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------
# 2. LOAD DATA
# --------------------------------------------------

ann <- fread(annotation_file)
tier1 <- fread(tier1_file)
filtering <- fread(filter_file)

required_ann <- c("gene_id", "GOs", "Preferred_name", "product")
required_tier1 <- c("gene_id", "GOs")

stopifnot(all(required_ann %in% names(ann)))
stopifnot(all(required_tier1 %in% names(tier1)))

ann[, gene_id := as.character(gene_id)]
tier1[, gene_id := as.character(gene_id)]

# --------------------------------------------------
# 3. DEFINE BACKGROUND AND FOREGROUND
# --------------------------------------------------

universe <- unique(
    filtering[
        retained == TRUE,
        as.character(gene_id)
    ]
)

universe <- universe[
    !grepl("^AMBIGUOUS_GROUP_", universe)
]

foreground <- unique(tier1$gene_id)

# --------------------------------------------------
# 4. BUILD GENE-GO MAPPING
# --------------------------------------------------

ann_go <- ann[
    !is.na(GOs) & trimws(GOs) != "",
    .(gene_id, GOs, Preferred_name, product)
]

ann_go <- unique(ann_go)

go_map <- ann_go[
    ,
    .(GO = unlist(strsplit(GOs, ",", fixed = TRUE))),
    by = .(gene_id, Preferred_name, product)
]

go_map[, GO := trimws(GO)]

go_map <- unique(
    go_map[GO != "" & !is.na(GO)]
)

# Restrict to genes in the tested universe
go_map <- go_map[gene_id %in% universe]

# Keep only GO identifiers recognized by GO.db
valid_go_ids <- AnnotationDbi::keys(
    GO.db,
    keytype = "GOID"
)

go_map <- go_map[GO %in% valid_go_ids]

# Remove duplicate gene-GO assignments
go_map <- unique(go_map, by = c("gene_id", "GO"))

annotated_universe <- intersect(
    universe,
    unique(go_map$gene_id)
)

annotated_foreground <- intersect(
    foreground,
    annotated_universe
)

cat("\n========== DATASET VALIDATION ==========\n")
cat("Tested biological genes:", length(universe), "\n")
cat("GO-annotated background genes:", length(annotated_universe), "\n")
cat("Tier 1 genes:", length(foreground), "\n")
cat("Tier 1 genes with valid GO:", length(annotated_foreground), "\n")
cat("Unique valid GO identifiers:",
    uniqueN(go_map$GO), "\n")

stopifnot(length(annotated_foreground) == 49)
stopifnot(length(annotated_universe) > length(annotated_foreground))

# --------------------------------------------------
# 5. ADD GO TERM DESCRIPTIONS AND ONTOLOGY
# --------------------------------------------------

go_info <- AnnotationDbi::select(
    GO.db,
    keys = unique(go_map$GO),
    keytype = "GOID",
    columns = c("TERM", "ONTOLOGY")
)

go_info <- unique(
    as.data.table(go_info),
    by = "GOID"
)

setnames(go_info, "GOID", "GO")

go_map <- merge(
    go_map,
    go_info,
    by = "GO",
    all = FALSE
)

# --------------------------------------------------
# 6. PREPARE HYPERGEOMETRIC TEST
# --------------------------------------------------

N <- length(annotated_universe)
n <- length(annotated_foreground)

min_size <- 5
max_size <- 500

ontology_names <- c(
    BP = "Biological_Process",
    MF = "Molecular_Function",
    CC = "Cellular_Component"
)

all_results <- list()

for (ont in names(ontology_names)) {

    cat("\n========== ONTOLOGY:", ont, "==========\n")

    ont_map <- go_map[ONTOLOGY == ont]

    term_stats <- ont_map[
        ,
        .(
            background_count = uniqueN(gene_id),
            foreground_count = uniqueN(
                gene_id[gene_id %in% annotated_foreground]
            ),
            genes = paste(
                sort(unique(gene_id[gene_id %in% annotated_foreground])),
                collapse = ";"
            ),
            gene_names = paste(
                sort(unique(
                    Preferred_name[
                        gene_id %in% annotated_foreground &
                        !is.na(Preferred_name) &
                        Preferred_name != ""
                    ]
                )),
                collapse = ";"
            ),
            products = paste(
                sort(unique(
                    product[
                        gene_id %in% annotated_foreground &
                        !is.na(product) &
                        product != ""
                    ]
                )),
                collapse = ";"
            )
        ),
        by = .(GO, TERM)
    ]

    term_stats <- term_stats[
        background_count >= min_size &
        background_count <= max_size
    ]

    term_stats <- term_stats[foreground_count > 0]

    if (nrow(term_stats) == 0) {
        cat("No terms passed the size and overlap criteria.\n")
        next
    }

    # One-sided hypergeometric test:
    # probability of observing at least k foreground genes
    # in a GO term of size K.

    term_stats[
        ,
        pvalue := phyper(
            foreground_count - 1,
            background_count,
            N - background_count,
            n,
            lower.tail = FALSE
        )
    ]

    term_stats[, FDR := p.adjust(pvalue, method = "BH")]

    term_stats[
        ,
        GeneRatio := paste0(foreground_count, "/", n)
    ]

    term_stats[
        ,
        BgRatio := paste0(background_count, "/", N)
    ]

    term_stats[
        ,
        FoldEnrichment := (
            foreground_count / n
        ) / (
            background_count / N
        )
    ]

    term_stats[, ONTOLOGY := ont]

    setorder(term_stats, FDR, pvalue)

    all_results[[ont]] <- term_stats

    cat("GO terms tested:", nrow(term_stats), "\n")
    cat("Significant terms FDR < 0.05:",
        sum(term_stats$FDR < 0.05), "\n")
}

# --------------------------------------------------
# 7. COMBINE RESULTS
# --------------------------------------------------

if (length(all_results) == 0) {
    stop("No GO terms were tested. Check annotations and thresholds.")
}

results <- rbindlist(all_results, fill = TRUE)

setcolorder(
    results,
    c(
        "ONTOLOGY",
        "GO",
        "TERM",
        "foreground_count",
        "background_count",
        "GeneRatio",
        "BgRatio",
        "FoldEnrichment",
        "pvalue",
        "FDR",
        "genes",
        "gene_names",
        "products"
    )
)

significant <- results[FDR < 0.05]

# --------------------------------------------------
# 8. SAVE RESULTS
# --------------------------------------------------

fwrite(
    results,
    file.path(output_dir, "GO_ORA_all_terms.csv")
)

fwrite(
    significant,
    file.path(output_dir, "GO_ORA_significant_FDR05.csv")
)

for (ont in names(ontology_names)) {

    subset_results <- results[ONTOLOGY == ont]

    fwrite(
        subset_results,
        file.path(
            output_dir,
            paste0("GO_ORA_", ont, "_all.csv")
        )
    )

    fwrite(
        subset_results[FDR < 0.05],
        file.path(
            output_dir,
            paste0("GO_ORA_", ont, "_significant.csv")
        )
    )
}

# --------------------------------------------------
# 9. SUMMARY
# --------------------------------------------------

summary_table <- results[
    ,
    .(
        tested_terms = .N,
        significant_terms = sum(FDR < 0.05),
        minimum_FDR = min(FDR),
        minimum_pvalue = min(pvalue)
    ),
    by = ONTOLOGY
]

fwrite(
    summary_table,
    file.path(output_dir, "GO_ORA_summary.csv")
)

cat("\n========== GO ORA SUMMARY ==========\n")
print(summary_table)

cat("\n========== TOP SIGNIFICANT GO TERMS ==========\n")

if (nrow(significant) > 0) {
    print(
        significant[
            order(FDR)
        ][
            ,
            .(
                ONTOLOGY,
                GO,
                TERM,
                foreground_count,
                background_count,
                FoldEnrichment,
                pvalue,
                FDR
            )
        ][1:min(20, .N)]
    )
} else {
    cat("No terms reached FDR < 0.05.\n")
}

cat("\n========== OUTPUT FILES ==========\n")
cat(output_dir, "\n")
cat("GO_ORA_all_terms.csv\n")
cat("GO_ORA_significant_FDR05.csv\n")
cat("GO_ORA_BP_all.csv / GO_ORA_BP_significant.csv\n")
cat("GO_ORA_MF_all.csv / GO_ORA_MF_significant.csv\n")
cat("GO_ORA_CC_all.csv / GO_ORA_CC_significant.csv\n")
cat("GO_ORA_summary.csv\n")

cat("\nPASS: GO over-representation analysis completed.\n")
