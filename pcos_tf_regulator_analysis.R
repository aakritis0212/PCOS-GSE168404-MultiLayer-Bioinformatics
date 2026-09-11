# ============================================================
# PCOS CORE-GENE TF / REGULATOR ENRICHMENT
# ============================================================

options(stringsAsFactors = FALSE)

cat("\n============================================================\n")
cat("        PCOS TF / REGULATOR ENRICHMENT ANALYSIS\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# DIRECTORIES
# ------------------------------------------------------------

input_dir <- "/Users/aakritisharma/newpipelin/dataset/PCOS_regulatory_analysis"

output_dir <- "/Users/aakritisharma/newpipelin/dataset/PCOS_TF_regulator_analysis"

if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

cat("Input directory:\n", input_dir, "\n\n")
cat("Output directory:\n", output_dir, "\n\n")


# ------------------------------------------------------------
# REQUIRED PACKAGES
# ------------------------------------------------------------

required <- c(
    "dplyr",
    "readr",
    "ggplot2",
    "tibble"
)

missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]

if (length(missing) > 0) {

    cat("Installing missing packages:\n")
    print(missing)

    install.packages(
        missing,
        repos = "https://cloud.r-project.org"
    )
}

library(dplyr)
library(readr)
library(ggplot2)
library(tibble)


# ------------------------------------------------------------
# FIND CORE GENE FILE
# ------------------------------------------------------------

core_file <- file.path(
    input_dir,
    "PCOS_core_genes_for_TF_enrichment.csv"
)

if (!file.exists(core_file)) {

    cat("\nERROR: Core gene file not found.\n")
    cat("Expected:\n")
    cat(core_file, "\n\n")

    cat("Files currently available:\n")
    print(list.files(input_dir))

    quit(save = "no", status = 1)
}


cat("Core gene file found:\n")
cat(basename(core_file), "\n\n")


# ------------------------------------------------------------
# LOAD CORE GENES
# ------------------------------------------------------------

core <- read_csv(
    core_file,
    show_col_types = FALSE
)

cat("Core gene rows:", nrow(core), "\n")
cat("Core gene columns:\n")
print(colnames(core))


# ------------------------------------------------------------
# DETECT GENE COLUMN
# ------------------------------------------------------------

gene_candidates <- c(
    "Gene",
    "GeneSymbol",
    "gene",
    "gene_symbol",
    "SYMBOL",
    "symbol"
)

gene_col <- gene_candidates[
    gene_candidates %in% colnames(core)
][1]

if (is.na(gene_col)) {

    cat("\nERROR: Could not identify gene column.\n")
    cat("Available columns:\n")
    print(colnames(core))

    quit(save = "no", status = 1)
}

cat("\nDetected gene column:", gene_col, "\n")


# ------------------------------------------------------------
# CLEAN GENE LIST
# ------------------------------------------------------------

genes <- core[[gene_col]]

genes <- as.character(genes)

genes <- trimws(genes)

genes <- genes[
    !is.na(genes) &
    genes != "" &
    genes != "NA"
]

genes <- unique(genes)

cat("Unique core genes:", length(genes), "\n\n")


# ------------------------------------------------------------
# SAVE CLEAN GENE LIST
# ------------------------------------------------------------

clean_gene_table <- tibble(
    GeneSymbol = genes
)

write_csv(
    clean_gene_table,
    file.path(
        output_dir,
        "PCOS_TF_enrichment_core_gene_list.csv"
    )
)


# ============================================================
# TF / REGULATOR DATABASE ANALYSIS
# ============================================================

cat("============================================================\n")
cat("        TF / REGULATOR ENRICHMENT PREPARATION\n")
cat("============================================================\n\n")

cat("The core gene set contains", length(genes), "genes.\n\n")

cat("Important:\n")
cat("This stage prepares the validated PCOS core-gene set for\n")
cat("transcription-factor / regulator enrichment.\n\n")


# ------------------------------------------------------------
# HIGH-CENTRALITY GENES
# ------------------------------------------------------------

high_centrality_file <- file.path(
    input_dir,
    "PCOS_high_centrality_genes_for_TF_enrichment.csv"
)

if (file.exists(high_centrality_file)) {

    high <- read_csv(
        high_centrality_file,
        show_col_types = FALSE
    )

    cat("High-centrality gene file found.\n")
    cat("Rows:", nrow(high), "\n")

    high_gene_col <- gene_candidates[
        gene_candidates %in% colnames(high)
    ][1]

    if (!is.na(high_gene_col)) {

        high_genes <- unique(
            trimws(as.character(high[[high_gene_col]]))
        )

        high_genes <- high_genes[
            !is.na(high_genes) &
            high_genes != ""
        ]

        write_csv(
            tibble(
                GeneSymbol = high_genes
            ),
            file.path(
                output_dir,
                "PCOS_high_centrality_gene_list.csv"
            )
        )

        cat(
            "High-centrality genes:",
            length(high_genes),
            "\n"
        )
    }
}


# ------------------------------------------------------------
# CORE GENE SUMMARY
# ------------------------------------------------------------

summary_table <- tibble(
    Metric = c(
        "Total core genes",
        "High-centrality genes"
    ),
    Value = c(
        length(genes),
        if (exists("high_genes")) length(high_genes) else NA
    )
)

write_csv(
    summary_table,
    file.path(
        output_dir,
        "PCOS_TF_ANALYSIS_INPUT_SUMMARY.csv"
    )
)


# ============================================================
# TF ENRICHMENT USING g:Profiler
# ============================================================

cat("\n============================================================\n")
cat("        TF ENRICHMENT\n")
cat("============================================================\n\n")

cat("For the actual TF enrichment, the recommended next step is\n")
cat("to query a regulator database such as TRRUST, ChEA, or\n")
cat("ENCODE-derived TF-target datasets.\n\n")

cat("This script therefore DOES NOT fabricate TF enrichment\n")
cat("statistics from the pathway data.\n\n")


# ------------------------------------------------------------
# WRITE GENE LIST AS TXT
# ------------------------------------------------------------

writeLines(
    genes,
    file.path(
        output_dir,
        "PCOS_core_genes_for_TF_analysis.txt"
    )
)


# ============================================================
# SIMPLE CORE-GENE CENTRALITY PLOT
# ============================================================

if (
    "Pathway_Count" %in% colnames(core)
) {

    plot_data <- core %>%
        mutate(
            GeneSymbol = as.character(.data[[gene_col]]),
            Pathway_Count = as.numeric(Pathway_Count)
        ) %>%
        filter(
            !is.na(Pathway_Count)
        ) %>%
        arrange(
            desc(Pathway_Count)
        ) %>%
        slice_head(n = 30)

    if (nrow(plot_data) > 0) {

        p <- ggplot(
            plot_data,
            aes(
                x = reorder(GeneSymbol, Pathway_Count),
                y = Pathway_Count
            )
        ) +
            geom_col() +
            coord_flip() +
            labs(
                title = "Top 30 PCOS Core Genes by Pathway Centrality",
                x = "Gene",
                y = "Number of Associated Pathways"
            ) +
            theme_minimal()

        ggsave(
            file.path(
                output_dir,
                "PCOS_top30_core_gene_centrality.png"
            ),
            p,
            width = 9,
            height = 8,
            dpi = 300
        )

        ggsave(
            file.path(
                output_dir,
                "PCOS_top30_core_gene_centrality.pdf"
            ),
            p,
            width = 9,
            height = 8
        )
    }
}


# ============================================================
# MASTER SUMMARY
# ============================================================

master <- tibble(
    Analysis = c(
        "Core gene set",
        "High-centrality gene set",
        "TF analysis"
    ),
    Result = c(
        paste(length(genes), "unique core genes"),
        if (exists("high_genes"))
            paste(length(high_genes), "high-centrality genes")
        else
            "Not available",
        "Gene list prepared for regulator enrichment"
    )
)

write_csv(
    master,
    file.path(
        output_dir,
        "PCOS_TF_REGULATOR_MASTER_SUMMARY.csv"
    )
)


# ============================================================
# COMPLETE
# ============================================================

cat("\n============================================================\n")
cat("        TF / REGULATOR ANALYSIS INPUT READY\n")
cat("============================================================\n\n")

cat("Core genes:", length(genes), "\n")

if (exists("high_genes")) {
    cat("High-centrality genes:", length(high_genes), "\n")
}

cat("\nOutput directory:\n")
cat(output_dir, "\n\n")

cat("Generated files:\n")
cat("1. PCOS_TF_enrichment_core_gene_list.csv\n")
cat("2. PCOS_high_centrality_gene_list.csv\n")
cat("3. PCOS_core_genes_for_TF_analysis.txt\n")
cat("4. PCOS_TF_ANALYSIS_INPUT_SUMMARY.csv\n")
cat("5. PCOS_TF_REGULATOR_MASTER_SUMMARY.csv\n")
cat("6. PCOS_top30_core_gene_centrality.png/pdf\n")

cat("\n============================================================\n")
cat("             STAGE COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")