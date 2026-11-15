#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(GO.db)
    library(AnnotationDbi)
})

base_dir <- "results/go_enrichment/ORA"
out_dir <- "results/go_enrichment/directional_ORA"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ann <- fread("results/integrated_analysis/integrated_DE_functional_annotation.csv")
tier1 <- fread("results/pathway_enrichment_grouped/biological_interpretation/tier1_50_genes_functional_characterization.csv")
filtering <- fread("results/differential_expression_grouped/edgeR_gene_filtering_grouped.tsv")

universe <- unique(filtering[retained == TRUE, as.character(gene_id)])
universe <- universe[!grepl("^AMBIGUOUS_GROUP_", universe)]

ann_go <- ann[!is.na(GOs) & trimws(GOs) != "",
              .(gene_id, GOs)]
ann_go <- unique(ann_go)

go_map <- ann_go[
    ,
    .(GO = unlist(strsplit(GOs, ",", fixed = TRUE))),
    by = gene_id
]
go_map[, GO := trimws(GO)]
go_map <- unique(go_map[GO != "" & gene_id %in% universe])

valid_go <- AnnotationDbi::keys(GO.db, keytype = "GOID")
go_map <- go_map[GO %in% valid_go]

go_info <- as.data.table(AnnotationDbi::select(
    GO.db,
    keys = unique(go_map$GO),
    keytype = "GOID",
    columns = c("TERM", "ONTOLOGY")
))
setnames(go_info, "GOID", "GO")
go_info <- unique(go_info, by = "GO")

go_map <- merge(go_map, go_info, by = "GO")

N <- uniqueN(go_map$gene_id)

tier1[, gene_id := as.character(gene_id)]

directions <- list(
    up_in_resistant = tier1[consensus_direction == "up_in_resistant", gene_id],
    down_in_resistant = tier1[consensus_direction == "down_in_resistant", gene_id]
)

cat("GO background:", N, "\n")
cat("Tier 1 up in Resistant:", length(directions$up_in_resistant), "\n")
cat("Tier 1 down in Resistant:", length(directions$down_in_resistant), "\n")

all_results <- list()

for (direction in names(directions)) {

    foreground <- intersect(
        unique(directions[[direction]]),
        unique(go_map$gene_id)
    )

    n <- length(foreground)

    cat("\nDirection:", direction, "\n")
    cat("Genes with GO:", n, "\n")

    if (n < 2) next

    for (ont in c("BP", "MF", "CC")) {

        subset_map <- go_map[ONTOLOGY == ont]

        stats <- subset_map[
            ,
            .(
                background_count = uniqueN(gene_id),
                foreground_count = uniqueN(gene_id[gene_id %in% foreground]),
                genes = paste(
                    sort(unique(gene_id[gene_id %in% foreground])),
                    collapse = ";"
                )
            ),
            by = .(GO, TERM)
        ]

        stats <- stats[
            background_count >= 5 &
            background_count <= 500 &
            foreground_count > 0
        ]

        if (nrow(stats) == 0) next

        stats[, pvalue := phyper(
            foreground_count - 1,
            background_count,
            N - background_count,
            n,
            lower.tail = FALSE
        )]

        stats[, FDR := p.adjust(pvalue, method = "BH")]
        stats[, FoldEnrichment := (foreground_count / n) /
                                      (background_count / N)]
        stats[, GeneRatio := paste0(foreground_count, "/", n)]
        stats[, BgRatio := paste0(background_count, "/", N)]
        stats[, direction := direction]
        stats[, ONTOLOGY := ont]

        setorder(stats, FDR, pvalue)

        fwrite(
            stats,
            file.path(out_dir, paste0("GO_ORA_", direction, "_", ont, ".csv"))
        )

        all_results[[paste(direction, ont, sep = "_")]] <- stats

        cat(
            ont,
            "- tested:", nrow(stats),
            "- significant:", sum(stats$FDR < 0.05), "\n"
        )
    }
}

combined <- rbindlist(all_results, fill = TRUE)

fwrite(
    combined,
    file.path(out_dir, "GO_ORA_directional_all.csv")
)

fwrite(
    combined[FDR < 0.05],
    file.path(out_dir, "GO_ORA_directional_significant.csv")
)

cat("\n========== TOP DIRECTIONAL GO TERMS ==========\n")

if (nrow(combined[FDR < 0.05]) > 0) {
    print(
        combined[FDR < 0.05][order(FDR)][
            ,
            .(direction, ONTOLOGY, GO, TERM,
              foreground_count, FoldEnrichment, FDR, genes)
        ][1:min(30, .N)]
    )
} else {
    cat("No directional GO terms passed FDR < 0.05.\n")
}

cat("\nPASS: Directional GO ORA completed.\n")
