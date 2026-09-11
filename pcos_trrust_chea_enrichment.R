############################################################
# PCOS TRRUST + ChEA TF/REGULATOR ENRICHMENT
#
# Purpose:
#   Identify transcription factors/regulators enriched among
#   PCOS GSEA-derived core genes and high-centrality genes.
#
# Databases:
#   1. TRRUST human
#   2. ChEA/TF-target data supplied as tftargets.rda
#
# IMPORTANT:
#   No regulator statistics are fabricated.
############################################################

options(stringsAsFactors = FALSE)

cat("\n============================================================\n")
cat(" PCOS TRRUST + ChEA REGULATOR ENRICHMENT\n")
cat("============================================================\n\n")

############################################################
# 1. DIRECTORIES
############################################################

base_dir <- "/Users/aakritisharma/newpipelin/dataset"

input_dir <- file.path(
    base_dir,
    "PCOS_TF_regulator_analysis"
)

db_dir <- file.path(
    base_dir,
    "PCOS_TF_regulator_enrichment",
    "databases"
)

output_dir <- file.path(
    base_dir,
    "PCOS_TF_regulator_enrichment",
    "results"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Input directory:\n", input_dir, "\n\n")
cat("Database directory:\n", db_dir, "\n\n")
cat("Output directory:\n", output_dir, "\n\n")


############################################################
# 2. PACKAGES
############################################################

required <- c(
    "dplyr",
    "ggplot2",
    "readr",
    "tidyr",
    "stringr"
)

missing <- required[
    !sapply(required, requireNamespace, quietly = TRUE)
]

if(length(missing) > 0){

    cat("Installing missing packages:\n")
    print(missing)

    install.packages(
        missing,
        repos = "https://cloud.r-project.org"
    )
}

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(readr)
    library(tidyr)
    library(stringr)
})


############################################################
# 3. LOCATE CORE GENE FILE
############################################################

core_candidates <- c(
    file.path(
        input_dir,
        "PCOS_TF_enrichment_core_gene_list.csv"
    ),
    file.path(
        input_dir,
        "PCOS_core_genes_for_TF_enrichment.csv"
    )
)

core_file <- core_candidates[
    file.exists(core_candidates)
][1]

if(is.na(core_file)){
    stop(
        "Could not find the PCOS core-gene file."
    )
}

cat("Core gene file:\n", basename(core_file), "\n\n")

core_df <- read_csv(
    core_file,
    show_col_types = FALSE
)

cat("Core-gene columns:\n")
print(names(core_df))
cat("\n")


############################################################
# 4. DETECT GENE COLUMN
############################################################

gene_candidates <- c(
    "GeneSymbol",
    "Gene",
    "gene",
    "Symbol",
    "SYMBOL",
    "GeneID"
)

gene_col <- gene_candidates[
    gene_candidates %in% names(core_df)
][1]

if(is.na(gene_col)){
    stop(
        "Could not detect gene column in core-gene file."
    )
}

core_genes <- unique(
    toupper(
        trimws(
            as.character(core_df[[gene_col]])
        )
    )
)

core_genes <- core_genes[
    !is.na(core_genes) &
    core_genes != ""
]

cat("Detected gene column:", gene_col, "\n")
cat("Unique core genes:", length(core_genes), "\n\n")


############################################################
# 5. HIGH-CENTRALITY GENES
############################################################

high_candidates <- c(
    file.path(
        input_dir,
        "PCOS_high_centrality_genes_for_TF_enrichment.csv"
    ),
    file.path(
        input_dir,
        "PCOS_high_centrality_gene_list.csv"
    )
)

high_file <- high_candidates[
    file.exists(high_candidates)
][1]

high_genes <- character(0)

if(!is.na(high_file)){

    high_df <- read_csv(
        high_file,
        show_col_types = FALSE
    )

    high_col <- gene_candidates[
        gene_candidates %in% names(high_df)
    ][1]

    if(!is.na(high_col)){

        high_genes <- unique(
            toupper(
                trimws(
                    as.character(high_df[[high_col]])
                )
            )
        )

        high_genes <- high_genes[
            !is.na(high_genes) &
            high_genes != ""
        ]
    }
}

cat("High-centrality genes:", length(high_genes), "\n\n")


############################################################
# 6. LOAD TRRUST
############################################################

cat("============================================================\n")
cat(" LOADING TRRUST\n")
cat("============================================================\n\n")

trrust_file <- file.path(
    db_dir,
    "trrust_rawdata.human.tsv"
)

if(!file.exists(trrust_file)){
    stop(
        "TRRUST file not found:\n",
        trrust_file
    )
}

trrust <- read.delim(
    trrust_file,
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
)

cat("TRRUST dimensions:\n")
print(dim(trrust))
cat("\n")

cat("First rows:\n")
print(head(trrust))
cat("\n")


############################################################
# 7. STANDARDIZE TRRUST
############################################################

if(ncol(trrust) < 2){
    stop("TRRUST file does not contain enough columns.")
}

trrust_std <- data.frame(
    TF = toupper(trimws(as.character(trrust[[1]]))),
    Target = toupper(trimws(as.character(trrust[[2]]))),
    stringsAsFactors = FALSE
)

trrust_std <- trrust_std[
    trrust_std$TF != "" &
    trrust_std$Target != "" &
    !is.na(trrust_std$TF) &
    !is.na(trrust_std$Target),
]

trrust_std <- unique(trrust_std)

cat("TRRUST TF-target interactions:", nrow(trrust_std), "\n")
cat("TRRUST regulators:", length(unique(trrust_std$TF)), "\n\n")


############################################################
# 8. LOAD tftargets.rda
############################################################

cat("============================================================\n")
cat(" LOADING ChEA / TF TARGET DATABASE\n")
cat("============================================================\n\n")

rda_file <- file.path(
    db_dir,
    "tftargets.rda"
)

if(!file.exists(rda_file)){
    stop(
        "tftargets.rda not found:\n",
        rda_file
    )
}

rda_env <- new.env()

load(
    rda_file,
    envir = rda_env
)

rda_objects <- ls(rda_env)

cat("Objects found in tftargets.rda:\n")
print(rda_objects)
cat("\n")


############################################################
# 9. FIND DATA FRAME / TABLE IN RDA
############################################################

extract_tables <- function(env){

    objects <- ls(env)

    result <- list()

    for(obj in objects){

        x <- get(obj, envir = env)

        if(
            is.data.frame(x) ||
            is.matrix(x)
        ){

            result[[obj]] <- as.data.frame(x)
        }
    }

    result
}

tables <- extract_tables(rda_env)

if(length(tables) == 0){

    stop(
        paste(
            "No data.frame or matrix object was found",
            "inside tftargets.rda."
        )
    )
}

cat("Tabular objects available:\n")

for(nm in names(tables)){

    cat(
        "\nObject:",
        nm,
        "\nDimensions:",
        paste(dim(tables[[nm]]), collapse = " x "),
        "\nColumns:\n"
    )

    print(names(tables[[nm]]))
}

cat("\n")


############################################################
# 10. AUTOMATICALLY IDENTIFY TF-TARGET TABLE
############################################################

find_tf_target_table <- function(tables){

    for(nm in names(tables)){

        x <- tables[[nm]]

        cn <- tolower(names(x))

        tf_candidates <- which(
            cn %in% c(
                "tf",
                "transcription_factor",
                "transcriptionfactor",
                "regulator",
                "source",
                "factor"
            )
        )

        target_candidates <- which(
            cn %in% c(
                "target",
                "target_gene",
                "targetgene",
                "gene",
                "genes"
            )
        )

        if(
            length(tf_candidates) >= 1 &&
            length(target_candidates) >= 1
        ){

            return(
                list(
                    name = nm,
                    data = x,
                    tf_col = tf_candidates[1],
                    target_col = target_candidates[1]
                )
            )
        }
    }

    return(NULL)
}

chea_info <- find_tf_target_table(tables)


############################################################
# 11. IF COLUMN NAMES ARE NOT STANDARD,
#     TRY FIRST TWO COLUMNS
############################################################

if(is.null(chea_info)){

    cat(
        "Could not automatically identify TF/target columns.\n"
    )

    cat(
        "Using the first suitable two-column table.\n\n"
    )

    for(nm in names(tables)){

        x <- tables[[nm]]

        if(ncol(x) >= 2){

            chea_info <- list(
                name = nm,
                data = x,
                tf_col = 1,
                target_col = 2
            )

            break
        }
    }
}

if(is.null(chea_info)){
    stop(
        "Could not identify a usable TF-target table."
    )
}

cat("Selected TF-target object:\n")
cat(chea_info$name, "\n\n")

chea_raw <- chea_info$data

chea_std <- data.frame(
    TF = toupper(
        trimws(
            as.character(
                chea_raw[[chea_info$tf_col]]
            )
        )
    ),
    Target = toupper(
        trimws(
            as.character(
                chea_raw[[chea_info$target_col]]
            )
        )
    ),
    stringsAsFactors = FALSE
)

chea_std <- chea_std[
    chea_std$TF != "" &
    chea_std$Target != "" &
    !is.na(chea_std$TF) &
    !is.na(chea_std$Target),
]

chea_std <- unique(chea_std)

cat(
    "TF-target interactions detected:",
    nrow(chea_std),
    "\n"
)

cat(
    "Unique regulators:",
    length(unique(chea_std$TF)),
    "\n\n"
)


############################################################
# 12. ENRICHMENT FUNCTION
############################################################

run_regulator_enrichment <- function(
    genes,
    interactions,
    database_name
){

    genes <- unique(
        toupper(
            trimws(
                as.character(genes)
            )
        )
    )

    genes <- genes[
        !is.na(genes) &
        genes != ""
    ]

    universe <- unique(
        c(
            interactions$Target,
            genes
        )
    )

    regulators <- unique(
        interactions$TF
    )

    results <- vector(
        "list",
        length(regulators)
    )

    for(i in seq_along(regulators)){

        tf <- regulators[i]

        targets <- unique(
            interactions$Target[
                interactions$TF == tf
            ]
        )

        N <- length(universe)

        K <- length(
            intersect(
                targets,
                universe
            )
        )

        n <- length(genes)

        k <- length(
            intersect(
                genes,
                targets
            )
        )

        if(K == 0 || k == 0){
            p <- 1
        } else {

            p <- phyper(
                k - 1,
                K,
                N - K,
                n,
                lower.tail = FALSE
            )
        }

        results[[i]] <- data.frame(
            Database = database_name,
            Regulator = tf,
            Target_Count = K,
            Core_Gene_Count = k,
            Universe_Size = N,
            Input_Gene_Count = n,
            Fold_Enrichment =
                ifelse(
                    K > 0,
                    (k / n) / (K / N),
                    NA
                ),
            P_Value = p,
            stringsAsFactors = FALSE
        )
    }

    result <- bind_rows(results)

    result$FDR <- p.adjust(
        result$P_Value,
        method = "BH"
    )

    result <- result %>%
        arrange(FDR, P_Value)

    result
}


############################################################
# 13. RUN TRRUST - CORE GENES
############################################################

cat("============================================================\n")
cat(" TRRUST ENRICHMENT - CORE GENES\n")
cat("============================================================\n\n")

trrust_core <- run_regulator_enrichment(
    core_genes,
    trrust_std,
    "TRRUST"
)

write.csv(
    trrust_core,
    file.path(
        output_dir,
        "PCOS_TRRUST_core_gene_enrichment.csv"
    ),
    row.names = FALSE
)

cat(
    "Significant TRRUST regulators FDR < 0.05:",
    sum(trrust_core$FDR < 0.05, na.rm = TRUE),
    "\n\n"
)


############################################################
# 14. RUN TRRUST - HIGH CENTRALITY
############################################################

if(length(high_genes) > 0){

    trrust_high <- run_regulator_enrichment(
        high_genes,
        trrust_std,
        "TRRUST"
    )

    write.csv(
        trrust_high,
        file.path(
            output_dir,
            "PCOS_TRRUST_high_centrality_enrichment.csv"
        ),
        row.names = FALSE
    )

} else {

    trrust_high <- data.frame()
}


############################################################
# 15. RUN ChEA/TF TARGETS - CORE
############################################################

cat("============================================================\n")
cat(" ChEA / TF TARGET ENRICHMENT - CORE GENES\n")
cat("============================================================\n\n")

chea_core <- run_regulator_enrichment(
    core_genes,
    chea_std,
    "ChEA_TF_Targets"
)

write.csv(
    chea_core,
    file.path(
        output_dir,
        "PCOS_ChEA_core_gene_enrichment.csv"
    ),
    row.names = FALSE
)

cat(
    "Significant ChEA/TF regulators FDR < 0.05:",
    sum(chea_core$FDR < 0.05, na.rm = TRUE),
    "\n\n"
)


############################################################
# 16. ChEA - HIGH CENTRALITY
############################################################

if(length(high_genes) > 0){

    chea_high <- run_regulator_enrichment(
        high_genes,
        chea_std,
        "ChEA_TF_Targets"
    )

    write.csv(
        chea_high,
        file.path(
            output_dir,
            "PCOS_ChEA_high_centrality_enrichment.csv"
        ),
        row.names = FALSE
    )

} else {

    chea_high <- data.frame()
}


############################################################
# 17. COMBINE RESULTS
############################################################

all_results <- bind_rows(
    trrust_core,
    trrust_high,
    chea_core,
    chea_high
)

write.csv(
    all_results,
    file.path(
        output_dir,
        "PCOS_ALL_REGULATOR_ENRICHMENT_RESULTS.csv"
    ),
    row.names = FALSE
)


############################################################
# 18. SIGNIFICANT REGULATORS
############################################################

significant <- all_results %>%
    filter(
        !is.na(FDR),
        FDR < 0.05
    ) %>%
    arrange(
        FDR,
        desc(Fold_Enrichment)
    )

write.csv(
    significant,
    file.path(
        output_dir,
        "PCOS_SIGNIFICANT_REGULATORS_FDR05.csv"
    ),
    row.names = FALSE
)


############################################################
# 19. TOP REGULATORS
############################################################

top_regulators <- significant %>%
    group_by(Database) %>%
    slice_min(
        order_by = FDR,
        n = 30,
        with_ties = FALSE
    ) %>%
    ungroup()

write.csv(
    top_regulators,
    file.path(
        output_dir,
        "PCOS_TOP30_REGULATORS.csv"
    ),
    row.names = FALSE
)


############################################################
# 20. REGULATOR CONSENSUS
############################################################

consensus <- significant %>%
    group_by(Regulator) %>%
    summarise(
        Database_Count = n_distinct(Database),
        Best_FDR = min(FDR, na.rm = TRUE),
        Best_P_Value = min(P_Value, na.rm = TRUE),
        Max_Fold_Enrichment =
            max(
                Fold_Enrichment,
                na.rm = TRUE
            ),
        Total_Targets =
            max(
                Target_Count,
                na.rm = TRUE
            ),
        Total_Input_Target_Hits =
            max(
                Core_Gene_Count,
                na.rm = TRUE
            ),
        .groups = "drop"
    ) %>%
    arrange(
        Best_FDR,
        desc(Database_Count)
    )

write.csv(
    consensus,
    file.path(
        output_dir,
        "PCOS_REGULATOR_CONSENSUS_RANKING.csv"
    ),
    row.names = FALSE
)


############################################################
# 21. TOP REGULATOR PLOT
############################################################

plot_data <- significant %>%
    arrange(FDR) %>%
    group_by(Database) %>%
    slice_head(n = 15) %>%
    ungroup()

if(nrow(plot_data) > 0){

    plot_data$Regulator <- factor(
        plot_data$Regulator,
        levels =
            rev(
                unique(
                    plot_data$Regulator
                )
            )
    )

    p <- ggplot(
        plot_data,
        aes(
            x = Regulator,
            y = -log10(FDR),
            size = Fold_Enrichment
        )
    ) +
        geom_point() +
        coord_flip() +
        facet_wrap(
            ~Database,
            scales = "free_y"
        ) +
        theme_bw() +
        labs(
            title =
                "PCOS Core-Gene Regulator Enrichment",
            x = "Regulator",
            y = "-log10(FDR)",
            size = "Fold enrichment"
        )

    ggsave(
        file.path(
            output_dir,
            "PCOS_top_regulators.png"
        ),
        p,
        width = 11,
        height = 8,
        dpi = 300
    )

    ggsave(
        file.path(
            output_dir,
            "PCOS_top_regulators.pdf"
        ),
        p,
        width = 11,
        height = 8
    )
}


############################################################
# 22. DATABASE COMPARISON
############################################################

database_summary <- all_results %>%
    group_by(Database) %>%
    summarise(
        Regulators_Tested =
            n(),
        Significant_FDR05 =
            sum(
                FDR < 0.05,
                na.rm = TRUE
            ),
        Best_FDR =
            min(
                FDR,
                na.rm = TRUE
            ),
        .groups = "drop"
    )

write.csv(
    database_summary,
    file.path(
        output_dir,
        "PCOS_REGULATOR_DATABASE_SUMMARY.csv"
    ),
    row.names = FALSE
)


############################################################
# 23. MASTER SUMMARY
############################################################

master_summary <- data.frame(
    Metric = c(
        "Core genes",
        "High-centrality genes",
        "TRRUST interactions",
        "TRRUST regulators",
        "ChEA/TF-target interactions",
        "ChEA/TF regulators",
        "TRRUST significant regulators",
        "ChEA significant regulators",
        "Total significant regulator records",
        "Unique significant regulators"
    ),
    Value = c(
        length(core_genes),
        length(high_genes),
        nrow(trrust_std),
        length(unique(trrust_std$TF)),
        nrow(chea_std),
        length(unique(chea_std$TF)),
        sum(
            trrust_core$FDR < 0.05,
            na.rm = TRUE
        ),
        sum(
            chea_core$FDR < 0.05,
            na.rm = TRUE
        ),
        nrow(significant),
        length(
            unique(
                significant$Regulator
            )
        )
    )
)

write.csv(
    master_summary,
    file.path(
        output_dir,
        "PCOS_REGULATOR_ENRICHMENT_MASTER_SUMMARY.csv"
    ),
    row.names = FALSE
)


############################################################
# 24. PRINT TOP RESULTS
############################################################

cat("\n============================================================\n")
cat(" TOP SIGNIFICANT REGULATORS\n")
cat("============================================================\n\n")

if(nrow(significant) > 0){

    print(
        significant %>%
            select(
                Database,
                Regulator,
                Core_Gene_Count,
                Fold_Enrichment,
                P_Value,
                FDR
            ) %>%
            head(30)
    )

} else {

    cat(
        "No regulator passed FDR < 0.05.\n"
    )
}


############################################################
# 25. FINAL STATUS
############################################################

cat("\n============================================================\n")
cat(" PCOS REGULATOR ENRICHMENT COMPLETE\n")
cat("============================================================\n\n")

cat(
    "Core genes:",
    length(core_genes),
    "\n"
)

cat(
    "High-centrality genes:",
    length(high_genes),
    "\n"
)

cat(
    "TRRUST interactions:",
    nrow(trrust_std),
    "\n"
)

cat(
    "ChEA/TF-target interactions:",
    nrow(chea_std),
    "\n"
)

cat(
    "Significant TRRUST regulators:",
    sum(
        trrust_core$FDR < 0.05,
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "Significant ChEA regulators:",
    sum(
        chea_core$FDR < 0.05,
        na.rm = TRUE
    ),
    "\n"
)

cat(
    "\nResults saved to:\n",
    output_dir,
    "\n\n"
)

cat("============================================================\n")
cat(" ANALYSIS FINISHED SUCCESSFULLY\n")
cat("============================================================\n")