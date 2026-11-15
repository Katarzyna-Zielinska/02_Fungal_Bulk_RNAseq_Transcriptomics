#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

input_file <- "results/go_enrichment/directional_ORA/GO_ORA_directional_all.csv"
output_dir <- "results/go_enrichment/directional_ORA/figures"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
    stop("Missing directional GO ORA results: ", input_file)
}

go <- fread(input_file)

required <- c(
    "direction", "ONTOLOGY", "GO", "TERM",
    "foreground_count", "FoldEnrichment", "FDR"
)

stopifnot(all(required %in% names(go)))

go <- go[
    !is.na(FDR) &
    FDR > 0 &
    is.finite(FoldEnrichment) &
    FoldEnrichment > 0
]

go[, neg_log10_FDR := -log10(FDR)]

ontology_labels <- c(
    BP = "Biological Process",
    MF = "Molecular Function",
    CC = "Cellular Component"
)

direction_labels <- c(
    up_in_resistant = "Upregulated in Resistant",
    down_in_resistant = "Downregulated in Resistant"
)

make_plot <- function(direction_value, output_prefix) {

    plot_data <- go[
        direction == direction_value & FDR < 0.05
    ]

    if (nrow(plot_data) == 0) {
        cat("No significant terms for:", direction_value, "\n")
        return(NULL)
    }

    # Select 10 most significant terms per ontology.
    plot_data <- plot_data[
        order(FDR),
        head(.SD, 10),
        by = ONTOLOGY
    ]

    plot_data[, ONTOLOGY_LABEL :=
        factor(
            ontology_labels[ONTOLOGY],
            levels = ontology_labels[c("BP", "MF", "CC")]
        )
    ]

    # Order terms by enrichment within the displayed dataset.
    plot_data[, TERM_ORDER := paste(
        ONTOLOGY,
        TERM,
        sep = "__"
    )]

    plot_data[, TERM_ORDER := factor(
        TERM_ORDER,
        levels = unique(TERM_ORDER[order(FoldEnrichment)])
    )]

    p <- ggplot(
        plot_data,
        aes(
            x = FoldEnrichment,
            y = TERM_ORDER
        )
    ) +
        geom_point(
            aes(
                size = foreground_count,
                colour = neg_log10_FDR
            ),
            alpha = 0.9
        ) +
        facet_wrap(
            ~ ONTOLOGY_LABEL,
            scales = "free_y",
            ncol = 1
        ) +
        scale_y_discrete(
            labels = function(x) sub("^.*__", "", x)
        ) +
        scale_size_continuous(
            name = "Number of genes",
            range = c(2.5, 8)
        ) +
        scale_colour_viridis_c(
            name = expression(-log[10]("FDR"))
        ) +
        labs(
            title = direction_labels[direction_value],
            subtitle = "GO over-representation analysis; top 10 terms per ontology",
            x = "Fold Enrichment",
            y = NULL,
            caption = paste(
                "Foreground: Tier 1 genes with GO annotations.",
                "Background: tested genes with GO annotations.",
                "Significance threshold: FDR < 0.05."
            )
        ) +
        theme_bw(base_size = 11) +
        theme(
            plot.title = element_text(
                face = "bold",
                size = 15
            ),
            plot.subtitle = element_text(size = 10),
            strip.text = element_text(face = "bold"),
            axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 9),
            legend.position = "right",
            panel.grid.minor = element_blank(),
            plot.caption = element_text(
                size = 8,
                hjust = 0
            )
        )

    ggsave(
        filename = file.path(
            output_dir,
            paste0(output_prefix, ".pdf")
        ),
        plot = p,
        width = 11,
        height = 13
    )

    ggsave(
        filename = file.path(
            output_dir,
            paste0(output_prefix, ".png")
        ),
        plot = p,
        width = 11,
        height = 13,
        dpi = 300
    )

    fwrite(
        plot_data,
        file.path(
            output_dir,
            paste0(output_prefix, "_plotted_terms.csv")
        )
    )

    cat("\nCreated:", output_prefix, "\n")
    cat("Terms plotted:", nrow(plot_data), "\n")

    return(p)
}

p_up <- make_plot(
    "up_in_resistant",
    "GO_dotplot_upregulated_resistant"
)

p_down <- make_plot(
    "down_in_resistant",
    "GO_dotplot_downregulated_resistant"
)

cat("\n========== OUTPUT VALIDATION ==========\n")

files <- list.files(
    output_dir,
    full.names = TRUE
)

print(files)

stopifnot(
    file.exists(file.path(
        output_dir,
        "GO_dotplot_upregulated_resistant.pdf"
    )),
    file.exists(file.path(
        output_dir,
        "GO_dotplot_downregulated_resistant.pdf"
    ))
)

cat("\nPASS: GO directional dot plots generated.\n")
