# ============================================================
# PCOS FINAL REGULATORY ANALYSIS
# TRRUST + ENCODE/TF TARGET DATABASE
# Core-gene regulatory architecture
# ============================================================

options(stringsAsFactors = FALSE)
options(scipen = 999)

cat("\n============================================================\n")
cat("        PCOS FINAL REGULATORY ANALYSIS\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------

required_packages <- c(
  "dplyr",
  "ggplot2",
  "tidyr",
  "readr"
)

missing <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)
]

if (length(missing) > 0) {
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
  library(tidyr)
  library(readr)
})

# ------------------------------------------------------------
# 2. DIRECTORIES
# ------------------------------------------------------------

base_dir <- "/Users/aakritisharma/newpipelin/dataset"

input_dir <- file.path(
  base_dir,
  "PCOS_TF_regulator_analysis"
)

database_dir <- file.path(
  base_dir,
  "PCOS_TF_regulator_enrichment",
  "databases"
)

output_dir <- file.path(
  base_dir,
  "PCOS_FINAL_REGULATORY_RESULTS"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Input directory:\n", input_dir, "\n\n")
cat("Database directory:\n", database_dir, "\n\n")
cat("Output directory:\n", output_dir, "\n\n")

# ------------------------------------------------------------
# 3. INPUT FILES
# ------------------------------------------------------------

core_file <- file.path(
  input_dir,
  "PCOS_TF_enrichment_core_gene_list.csv"
)

high_centrality_file <- file.path(
  input_dir,
  "PCOS_high_centrality_gene_list.csv"
)

trrust_file <- file.path(
  database_dir,
  "trrust_rawdata.human.tsv"
)

tftargets_file <- file.path(
  database_dir,
  "tftargets.rda"
)

# ------------------------------------------------------------
# 4. CHECK FILES
# ------------------------------------------------------------

files_required <- c(
  core_file,
  trrust_file,
  tftargets_file
)

missing_files <- files_required[
  !file.exists(files_required)
]

if (length(missing_files) > 0) {

  cat("ERROR: Required files are missing:\n")
  print(missing_files)

  stop(
    "\nPlease make sure the required database/input files are present."
  )
}

# ------------------------------------------------------------
# 5. LOAD CORE GENES
# ------------------------------------------------------------

cat("============================================================\n")
cat("                    LOADING CORE GENES\n")
cat("============================================================\n\n")

core <- read.csv(
  core_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Core-gene columns:\n")
print(names(core))

# Automatically detect gene column

gene_candidates <- c(
  "GeneSymbol",
  "Gene",
  "gene",
  "Gene.symbol",
  "Symbol"
)

gene_col <- gene_candidates[
  gene_candidates %in% names(core)
][1]

if (is.na(gene_col)) {
  stop("Could not identify the gene column in core-gene file.")
}

cat("\nDetected gene column:", gene_col, "\n")

core_genes <- unique(
  toupper(
    trimws(
      as.character(core[[gene_col]])
    )
  )
)

core_genes <- core_genes[
  !is.na(core_genes) &
  core_genes != ""
]

cat("Unique core genes:", length(core_genes), "\n\n")

# ------------------------------------------------------------
# 6. LOAD HIGH CENTRALITY GENES
# ------------------------------------------------------------

high_centrality_genes <- character(0)

if (file.exists(high_centrality_file)) {

  hc <- read.csv(
    high_centrality_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  hc_candidates <- c(
    "GeneSymbol",
    "Gene",
    "gene",
    "Symbol"
  )

  hc_col <- hc_candidates[
    hc_candidates %in% names(hc)
  ][1]

  if (!is.na(hc_col)) {

    high_centrality_genes <- unique(
      toupper(
        trimws(
          as.character(hc[[hc_col]])
        )
      )
    )

    high_centrality_genes <-
      high_centrality_genes[
        !is.na(high_centrality_genes) &
        high_centrality_genes != ""
      ]
  }
}

cat("High-centrality genes:", length(high_centrality_genes), "\n\n")

# ------------------------------------------------------------
# 7. LOAD TRRUST
# ------------------------------------------------------------

cat("============================================================\n")
cat("                       LOADING TRRUST\n")
cat("============================================================\n\n")

trrust <- read.delim(
  trrust_file,
  header = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

if (ncol(trrust) < 2) {
  stop("TRRUST file does not contain the expected columns.")
}

colnames(trrust)[1:4] <- c(
  "TF",
  "Target",
  "Effect",
  "PMID"
)

trrust$TF <- toupper(trimws(trrust$TF))
trrust$Target <- toupper(trimws(trrust$Target))

trrust <- trrust[
  trrust$TF != "" &
  trrust$Target != "",
]

cat("TRRUST interactions:", nrow(trrust), "\n")
cat("TRRUST regulators:", length(unique(trrust$TF)), "\n\n")

# ------------------------------------------------------------
# 8. TRRUST ENRICHMENT
# ------------------------------------------------------------

cat("============================================================\n")
cat("                     TRRUST ENRICHMENT\n")
cat("============================================================\n\n")

# Background = genes represented as targets in TRRUST

background_trrust <- unique(trrust$Target)

N <- length(background_trrust)

K <- length(core_genes[
  core_genes %in% background_trrust
])

cat("TRRUST background genes:", N, "\n")
cat("Core genes represented in TRRUST:", K, "\n\n")

regulators <- unique(trrust$TF)

trrust_results <- lapply(
  regulators,
  function(tf) {

    targets <- unique(
      trrust$Target[
        trrust$TF == tf
      ]
    )

    M <- length(targets)

    overlap_genes <- intersect(
      core_genes,
      targets
    )

    k <- length(overlap_genes)

    if (k == 0) {
      return(NULL)
    }

    # Hypergeometric enrichment
    pval <- phyper(
      q = k - 1,
      m = M,
      n = N - M,
      k = K,
      lower.tail = FALSE
    )

    odds_ratio <- (
      (k + 0.5) *
      (N - M - K + k + 0.5)
    ) /
      (
        (M - k + 0.5) *
        (K - k + 0.5)
      )

    data.frame(
      Regulator = tf,
      Target_Count = M,
      Overlap = k,
      OddsRatio = odds_ratio,
      PValue = pval,
      OverlapGenes = paste(
        overlap_genes,
        collapse = ";"
      ),
      stringsAsFactors = FALSE
    )
  }
)

trrust_results <- bind_rows(trrust_results)

if (nrow(trrust_results) > 0) {

  trrust_results$FDR <- p.adjust(
    trrust_results$PValue,
    method = "BH"
  )

  trrust_results <- trrust_results %>%
    arrange(FDR, PValue)

}

cat(
  "TRRUST regulators tested:",
  nrow(trrust_results),
  "\n"
)

cat(
  "Significant TRRUST regulators FDR < 0.05:",
  sum(trrust_results$FDR < 0.05),
  "\n\n"
)

write.csv(
  trrust_results,
  file.path(
    output_dir,
    "PCOS_TRRUST_enrichment_all.csv"
  ),
  row.names = FALSE
)

write.csv(
  trrust_results %>%
    filter(FDR < 0.05),
  file.path(
    output_dir,
    "PCOS_TRRUST_significant_regulators.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. LOAD TFTARGETS DATABASE
# ------------------------------------------------------------

cat("============================================================\n")
cat("                 LOADING TF TARGET DATABASE\n")
cat("============================================================\n\n")

tf_env <- new.env()

load(
  tftargets_file,
  envir = tf_env
)

objects <- ls(tf_env)

cat("Objects found:\n")
print(objects)
cat("\n")

# ------------------------------------------------------------
# 10. EXTRACT ENCODE
# ------------------------------------------------------------

encode_object <- NULL

if ("ENCODE" %in% objects) {

  encode_object <- get(
    "ENCODE",
    envir = tf_env
  )

  cat(
    "ENCODE object detected:",
    class(encode_object),
    "\n"
  )

} else {

  cat(
    "ENCODE object not found. Skipping ENCODE enrichment.\n"
  )
}

# ------------------------------------------------------------
# 11. GENERIC TF DATABASE EXTRACTION
# ------------------------------------------------------------

extract_pairs <- function(
    obj,
    regulator_name = "Regulator"
) {

  pairs <- data.frame(
    Regulator = character(),
    Target = character(),
    stringsAsFactors = FALSE
  )

  # Case 1: data.frame

  if (is.data.frame(obj) || is.matrix(obj)) {

    x <- as.data.frame(obj)

    if (ncol(x) >= 2) {

      pairs <- data.frame(
        Regulator = as.character(x[[1]]),
        Target = as.character(x[[2]]),
        stringsAsFactors = FALSE
      )
    }
  }

  # Case 2: list

  else if (is.list(obj)) {

    temp <- lapply(
      names(obj),
      function(nm) {

        element <- obj[[nm]]

        # Matrix/data.frame element

        if (
          is.data.frame(element) ||
          is.matrix(element)
        ) {

          element <- as.data.frame(element)

          if (ncol(element) >= 2) {

            return(
              data.frame(
                Regulator = nm,
                Target = as.character(
                  element[[2]]
                ),
                stringsAsFactors = FALSE
              )
            )
          }
        }

        # Character vector

        if (is.character(element)) {

          return(
            data.frame(
              Regulator = nm,
              Target = element,
              stringsAsFactors = FALSE
            )
          )
        }

        NULL
      }
    )

    pairs <- bind_rows(temp)
  }

  pairs$Regulator <- toupper(
    trimws(
      as.character(
        pairs$Regulator
      )
    )
  )

  pairs$Target <- toupper(
    trimws(
      as.character(
        pairs$Target
      )
    )
  )

  pairs <- pairs[
    pairs$Regulator != "" &
    pairs$Target != "" &
    !is.na(pairs$Regulator) &
    !is.na(pairs$Target),
  ]

  unique(pairs)
}

# ------------------------------------------------------------
# 12. ENCODE ENRICHMENT
# ------------------------------------------------------------

encode_results <- data.frame()

if (!is.null(encode_object)) {

  cat("\n============================================================\n")
  cat("                    ENCODE ENRICHMENT\n")
  cat("============================================================\n\n")

  encode_pairs <- extract_pairs(
    encode_object
  )

  cat(
    "Extracted ENCODE TF-target pairs:",
    nrow(encode_pairs),
    "\n"
  )

  if (nrow(encode_pairs) > 0) {

    background <- unique(
      encode_pairs$Target
    )

    N2 <- length(background)

    K2 <- length(
      intersect(
        core_genes,
        background
      )
    )

    encode_regs <- unique(
      encode_pairs$Regulator
    )

    encode_list <- lapply(
      encode_regs,
      function(tf) {

        targets <- unique(
          encode_pairs$Target[
            encode_pairs$Regulator == tf
          ]
        )

        M <- length(targets)

        overlap_genes <- intersect(
          core_genes,
          targets
        )

        k <- length(overlap_genes)

        if (k == 0) {
          return(NULL)
        }

        pval <- phyper(
          q = k - 1,
          m = M,
          n = N2 - M,
          k = K2,
          lower.tail = FALSE
        )

        data.frame(
          Regulator = tf,
          Target_Count = M,
          Overlap = k,
          PValue = pval,
          OverlapGenes = paste(
            overlap_genes,
            collapse = ";"
          ),
          stringsAsFactors = FALSE
        )
      }
    )

    encode_results <- bind_rows(
      encode_list
    )

    if (nrow(encode_results) > 0) {

      encode_results$FDR <- p.adjust(
        encode_results$PValue,
        method = "BH"
      )

      encode_results <- encode_results %>%
        arrange(FDR, PValue)

      write.csv(
        encode_results,
        file.path(
          output_dir,
          "PCOS_ENCODE_enrichment_all.csv"
        ),
        row.names = FALSE
      )

      write.csv(
        encode_results %>%
          filter(FDR < 0.05),
        file.path(
          output_dir,
          "PCOS_ENCODE_significant_regulators.csv"
        ),
        row.names = FALSE
      )
    }
  }
}

# ------------------------------------------------------------
# 13. COMBINE REGULATORY RESULTS
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("                 COMBINING REGULATORY RESULTS\n")
cat("============================================================\n\n")

trrust_combined <- trrust_results %>%
  mutate(
    Database = "TRRUST"
  ) %>%
  select(
    Database,
    Regulator,
    Target_Count,
    Overlap,
    OddsRatio,
    PValue,
    FDR,
    OverlapGenes
  )

if (nrow(encode_results) > 0) {

  encode_combined <- encode_results %>%
    mutate(
      Database = "ENCODE"
    ) %>%
    mutate(
      OddsRatio = NA_real_
    ) %>%
    select(
      Database,
      Regulator,
      Target_Count,
      Overlap,
      OddsRatio,
      PValue,
      FDR,
      OverlapGenes
    )

} else {

  encode_combined <- data.frame()
}

master <- bind_rows(
  trrust_combined,
  encode_combined
) %>%
  arrange(FDR, PValue)

write.csv(
  master,
  file.path(
    output_dir,
    "PCOS_FINAL_REGULATOR_MASTER_TABLE.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 14. SIGNIFICANT MASTER RESULTS
# ------------------------------------------------------------

significant_master <- master %>%
  filter(
    !is.na(FDR),
    FDR < 0.05
  ) %>%
  arrange(
    FDR,
    PValue
  )

write.csv(
  significant_master,
  file.path(
    output_dir,
    "PCOS_FINAL_SIGNIFICANT_REGULATORS.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 15. HIGH-CENTRALITY REGULATOR OVERLAP
# ------------------------------------------------------------

if (length(high_centrality_genes) > 0) {

  high_centrality_regulators <- master %>%
    filter(
      Regulator %in% high_centrality_genes
    ) %>%
    arrange(FDR)

  write.csv(
    high_centrality_regulators,
    file.path(
      output_dir,
      "PCOS_HIGH_CENTRALITY_REGULATORS.csv"
    ),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 16. TOP REGULATORS
# ------------------------------------------------------------

top_regulators <- master %>%
  filter(
    !is.na(FDR)
  ) %>%
  arrange(FDR, PValue) %>%
  slice_head(n = 30)

write.csv(
  top_regulators,
  file.path(
    output_dir,
    "PCOS_TOP30_REGULATORS.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 17. REGULATOR-GENE NETWORK
# ------------------------------------------------------------

network_edges <- master %>%
  filter(
    !is.na(FDR),
    FDR < 0.05
  ) %>%
  select(
    Regulator,
    OverlapGenes
  ) %>%
  separate_rows(
    OverlapGenes,
    sep = ";"
  ) %>%
  rename(
    Gene = OverlapGenes
  ) %>%
  filter(
    Gene != "",
    !is.na(Gene)
  ) %>%
  distinct()

write.csv(
  network_edges,
  file.path(
    output_dir,
    "PCOS_significant_regulator_gene_network_edges.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 18. REGULATOR CENTRALITY
# ------------------------------------------------------------

regulator_centrality <- network_edges %>%
  count(
    Regulator,
    name = "Core_Gene_Targets"
  ) %>%
  arrange(
    desc(Core_Gene_Targets)
  )

write.csv(
  regulator_centrality,
  file.path(
    output_dir,
    "PCOS_regulator_target_centrality.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 19. PLOT: TOP REGULATORS
# ------------------------------------------------------------

if (nrow(top_regulators) > 0) {

  plot_data <- top_regulators %>%
    slice_head(n = 20) %>%
    mutate(
      Regulator = reorder(
        Regulator,
        -log10(PValue)
      )
    )

  p1 <- ggplot(
    plot_data,
    aes(
      x = Regulator,
      y = -log10(PValue)
    )
  ) +
    geom_col() +
    coord_flip() +
    labs(
      title = "Top Regulatory Enrichment Signals in PCOS Core Genes",
      x = "Regulator",
      y = "-log10(P-value)"
    ) +
    theme_minimal(
      base_size = 12
    )

  ggsave(
    file.path(
      output_dir,
      "PCOS_top20_regulators.png"
    ),
    p1,
    width = 9,
    height = 7,
    dpi = 300
  )

  ggsave(
    file.path(
      output_dir,
      "PCOS_top20_regulators.pdf"
    ),
    p1,
    width = 9,
    height = 7
  )
}

# ------------------------------------------------------------
# 20. PLOT: REGULATOR TARGET COUNTS
# ------------------------------------------------------------

if (nrow(regulator_centrality) > 0) {

  plot_data2 <- regulator_centrality %>%
    slice_head(n = 20) %>%
    mutate(
      Regulator = reorder(
        Regulator,
        Core_Gene_Targets
      )
    )

  p2 <- ggplot(
    plot_data2,
    aes(
      x = Regulator,
      y = Core_Gene_Targets
    )
  ) +
    geom_col() +
    coord_flip() +
    labs(
      title = "Regulator–Core Gene Connectivity",
      x = "Regulator",
      y = "Number of Core-Gene Targets"
    ) +
    theme_minimal(
      base_size = 12
    )

  ggsave(
    file.path(
      output_dir,
      "PCOS_regulator_target_centrality.png"
    ),
    p2,
    width = 9,
    height = 7,
    dpi = 300
  )

  ggsave(
    file.path(
      output_dir,
      "PCOS_regulator_target_centrality.pdf"
    ),
    p2,
    width = 9,
    height = 7
  )
}

# ------------------------------------------------------------
# 21. FINAL SUMMARY
# ------------------------------------------------------------

summary_table <- data.frame(
  Metric = c(
    "Core genes",
    "High-centrality genes",
    "TRRUST interactions",
    "TRRUST regulators tested",
    "TRRUST significant regulators FDR < 0.05",
    "ENCODE results available",
    "Total regulator results",
    "Significant regulators FDR < 0.05"
  ),

  Value = c(
    length(core_genes),
    length(high_centrality_genes),
    nrow(trrust),
    nrow(trrust_results),
    sum(trrust_results$FDR < 0.05),
    nrow(encode_results) > 0,
    nrow(master),
    sum(master$FDR < 0.05, na.rm = TRUE)
  )
)

write.csv(
  summary_table,
  file.path(
    output_dir,
    "PCOS_FINAL_REGULATORY_ANALYSIS_SUMMARY.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 22. PRINT FINAL RESULTS
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("          FINAL REGULATORY ANALYSIS COMPLETE\n")
cat("============================================================\n\n")

cat("Core genes:", length(core_genes), "\n")
cat(
  "High-centrality genes:",
  length(high_centrality_genes),
  "\n"
)

cat(
  "TRRUST interactions:",
  nrow(trrust),
  "\n"
)

cat(
  "TRRUST significant regulators:",
  sum(trrust_results$FDR < 0.05),
  "\n"
)

cat(
  "Total regulator results:",
  nrow(master),
  "\n"
)

cat(
  "Significant regulators FDR < 0.05:",
  sum(master$FDR < 0.05, na.rm = TRUE),
  "\n\n"
)

cat("TOP REGULATORS:\n\n")

if (nrow(significant_master) > 0) {

  print(
    significant_master %>%
      select(
        Database,
        Regulator,
        Overlap,
        PValue,
        FDR
      ) %>%
      slice_head(n = 20)
  )

} else {

  cat(
    "No regulators reached FDR < 0.05.\n"
  )
}

cat("\n============================================================\n")
cat("OUTPUT DIRECTORY:\n")
cat(output_dir, "\n")
cat("============================================================\n\n")

cat("Major files generated:\n\n")

cat("1. PCOS_TRRUST_enrichment_all.csv\n")
cat("2. PCOS_TRRUST_significant_regulators.csv\n")
cat("3. PCOS_ENCODE_enrichment_all.csv [if available]\n")
cat("4. PCOS_ENCODE_significant_regulators.csv [if available]\n")
cat("5. PCOS_FINAL_REGULATOR_MASTER_TABLE.csv\n")
cat("6. PCOS_FINAL_SIGNIFICANT_REGULATORS.csv\n")
cat("7. PCOS_HIGH_CENTRALITY_REGULATORS.csv\n")
cat("8. PCOS_TOP30_REGULATORS.csv\n")
cat("9. PCOS_significant_regulator_gene_network_edges.csv\n")
cat("10. PCOS_regulator_target_centrality.csv\n")
cat("11. PCOS_top20_regulators.png/pdf\n")
cat("12. PCOS_regulator_target_centrality.png/pdf\n")
cat("13. PCOS_FINAL_REGULATORY_ANALYSIS_SUMMARY.csv\n\n")

cat("============================================================\n")
cat("                  ANALYSIS FINISHED\n")
cat("============================================================\n")
