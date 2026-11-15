#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
})

input_file <- "results/go_enrichment/biological_interpretation/GO_significant_terms_gene_expression_integration.csv"
out_dir <- "results/go_enrichment/biological_interpretation/core_genes"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

x <- fread(input_file)

# Selected representative terms from directional GO analysis.
core_terms <- data.table(
    GO = c(
        "GO:0006090",
        "GO:0005996",
        "GO:0006096",
        "GO:0042866",
        "GO:0046364",
        "GO:0006006",
        "GO:0006094",
        "GO:0005730",
        "GO:0016072",
        "GO:0006364",
        "GO:0030684",
        "GO:0032040",
        "GO:0030686",
        "GO:0042254",
        "GO:0022613",
        "GO:0006396"
    ),
    biological_group = c(
        rep("Energy_and_carbon_metabolism", 7),
        rep("Ribosome_biogenesis_and_RNA", 9)
    )
)

selected <- merge(
    x,
    core_terms,
    by = "GO",
    all = FALSE
)

# One record per gene and GO term.
selected <- unique(selected, by = c("gene_id", "GO"))

# Gene-level summary.
gene_summary <- selected[
    ,
    .(
        biological_group = paste(
            sort(unique(biological_group)),
            collapse = ";"
        ),
        GO_term_count = uniqueN(GO),
        GO_terms = paste(sort(unique(TERM)), collapse = "; "),
        GO_IDs = paste(sort(unique(GO)), collapse = ";"),
        gene_name = first(gene_name),
        Preferred_name = first(Preferred_name),
        product = first(product),
        COG_category = first(COG_category),
        KEGG_Pathway = first(KEGG_Pathway),
        edgeR_logFC = first(edgeR_logFC),
        edgeR_FDR = first(edgeR_FDR),
        limma_logFC = first(limma_logFC),
        limma_FDR = first(limma_FDR),
        consensus_logFC = first(consensus_logFC),
        consensus_direction = first(consensus_direction)
    ),
    by = gene_id
]

setorder(
    gene_summary,
    biological_group,
    -GO_term_count,
    gene_id
)

# Term-level gene membership.
term_membership <- selected[
    ,
    .(
        biological_group = first(biological_group),
        gene_count = uniqueN(gene_id),
        genes = paste(sort(unique(gene_id)), collapse = ";"),
        gene_names = paste(
            sort(unique(
                Preferred_name[
                    !is.na(Preferred_name) &
                    Preferred_name != ""
                ]
            )),
            collapse = ";"
        ),
        FDR = first(FDR),
        FoldEnrichment = first(FoldEnrichment)
    ),
    by = .(GO, TERM, direction, ONTOLOGY)
]

setorder(term_membership, biological_group, FDR)

# Save outputs.
fwrite(
    selected,
    file.path(out_dir, "core_GO_term_gene_membership.csv")
)

fwrite(
    gene_summary,
    file.path(out_dir, "core_GO_gene_summary.csv")
)

fwrite(
    gene_summary[
        grepl("Energy_and_carbon", biological_group)
    ],
    file.path(out_dir, "energy_metabolism_core_genes.csv")
)

fwrite(
    gene_summary[
        grepl("Ribosome_biogenesis", biological_group)
    ],
    file.path(out_dir, "ribosome_RNA_core_genes.csv")
)

fwrite(
    term_membership,
    file.path(out_dir, "core_GO_term_summary.csv")
)

cat("\n========== CORE GO GENES ==========\n")
cat("Unique genes:", uniqueN(gene_summary$gene_id), "\n")
cat("Energy/metabolism genes:",
    sum(grepl("Energy_and_carbon", gene_summary$biological_group)),
    "\n")
cat("Ribosome/RNA genes:",
    sum(grepl("Ribosome_biogenesis", gene_summary$biological_group)),
    "\n")

cat("\n========== TOP GENES ==========\n")
print(
    gene_summary[
        ,
        .(
            gene_id,
            Preferred_name,
            product,
            GO_term_count,
            consensus_logFC,
            consensus_direction
        )
    ][1:min(30, .N)]
)

cat("\n========== OUTPUT DIRECTORY ==========\n")
cat(out_dir, "\n")

cat("\nPASS: Core GO gene identification completed.\n")
