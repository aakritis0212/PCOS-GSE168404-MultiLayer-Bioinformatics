############################################################
# PCOS REGULATORY LAYER ANALYSIS
# GSE168404
#
# Input:
#   PCOS_post_GSEA_results/
#     PCOS_core_genes_3_or_more_pathways.csv
#     PCOS_leading_edge_gene_pathway_table.csv
#     PCOS_significant_GSEA_pathways.csv
#
# Output:
#   PCOS_regulatory_analysis/
#
# Purpose:
#   Identify pathway-central core genes and candidate
#   transcriptional regulators associated with the
#   PCOS transcriptomic phenotype.
############################################################

cat("\n============================================================\n")
cat("       PCOS CORE-GENE REGULATORY ANALYSIS\n")
cat("============================================================\n\n")

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------

required <- c(
  "dplyr",
  "readr",
  "stringr",
  "ggplot2"
)

missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]

if (length(missing) > 0) {
  cat("Missing packages:", paste(missing, collapse = ", "), "\n")
  cat("Install with:\n")
  cat(
    "install.packages(c(",
    paste0("'", missing, "'", collapse = ", "),
    "))\n"
  )
  quit(status = 1)
}

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)

# ------------------------------------------------------------
# 2. DIRECTORIES
# ------------------------------------------------------------

base_dir <- getwd()

gsea_dir <- file.path(
  base_dir,
  "PCOS_post_GSEA_results"
)

out_dir <- file.path(
  base_dir,
  "PCOS_regulatory_analysis"
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Input directory:\n", gsea_dir, "\n\n")
cat("Output directory:\n", out_dir, "\n\n")

# ------------------------------------------------------------
# 3. INPUT FILES
# ------------------------------------------------------------

core_file <- file.path(
  gsea_dir,
  "PCOS_core_genes_3_or_more_pathways.csv"
)

leading_file <- file.path(
  gsea_dir,
  "PCOS_leading_edge_gene_pathway_table.csv"
)

pathway_file <- file.path(
  gsea_dir,
  "PCOS_significant_GSEA_pathways.csv"
)

files <- c(
  core_file,
  leading_file,
  pathway_file
)

for (f in files) {
  if (!file.exists(f)) {
    stop(
      paste0(
        "\nERROR: Required file not found:\n",
        f
      )
    )
  }
}

# ------------------------------------------------------------
# 4. LOAD DATA
# ------------------------------------------------------------

cat("Loading core genes...\n")

core <- read_csv(
  core_file,
  show_col_types = FALSE
)

cat("Core-gene rows:", nrow(core), "\n")
cat("Core-gene columns:\n")
print(names(core))

cat("\nLoading leading-edge associations...\n")

leading <- read_csv(
  leading_file,
  show_col_types = FALSE
)

cat("Leading-edge rows:", nrow(leading), "\n")
cat("Leading-edge columns:\n")
print(names(leading))

cat("\nLoading significant pathways...\n")

pathways <- read_csv(
  pathway_file,
  show_col_types = FALSE
)

cat("Pathway rows:", nrow(pathways), "\n")
cat("Pathway columns:\n")
print(names(pathways))

# ------------------------------------------------------------
# 5. DETECT GENE COLUMN
# ------------------------------------------------------------

find_gene_column <- function(df) {

  candidates <- c(
    "GeneSymbol",
    "Gene",
    "gene",
    "Symbol",
    "symbol",
    "GeneID"
  )

  hit <- candidates[candidates %in% names(df)]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[1]
}

core_gene_col <- find_gene_column(core)

leading_gene_col <- find_gene_column(leading)

if (is.na(core_gene_col)) {
  stop("Could not identify gene column in core-gene file.")
}

if (is.na(leading_gene_col)) {
  stop("Could not identify gene column in leading-edge file.")
}

cat("\nCore gene column:", core_gene_col, "\n")
cat("Leading-edge gene column:", leading_gene_col, "\n")

# ------------------------------------------------------------
# 6. STANDARDIZE CORE GENE TABLE
# ------------------------------------------------------------

core_std <- core %>%
  mutate(
    GeneSymbol = as.character(.data[[core_gene_col]])
  ) %>%
  filter(
    !is.na(GeneSymbol),
    GeneSymbol != "",
    GeneSymbol != "NA"
  ) %>%
  distinct(GeneSymbol, .keep_all = TRUE)

cat(
  "\nUnique core genes:",
  nrow(core_std),
  "\n"
)

# ------------------------------------------------------------
# 7. STANDARDIZE LEADING-EDGE DATA
# ------------------------------------------------------------

leading_std <- leading %>%
  mutate(
    GeneSymbol = as.character(.data[[leading_gene_col]])
  ) %>%
  filter(
    !is.na(GeneSymbol),
    GeneSymbol != "",
    GeneSymbol != "NA"
  )

# ------------------------------------------------------------
# 8. IDENTIFY PATHWAY COLUMN
# ------------------------------------------------------------

find_pathway_column <- function(df) {

  candidates <- c(
    "Term",
    "pathway",
    "Pathway",
    "Description"
  )

  hit <- candidates[candidates %in% names(df)]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[1]
}

leading_pathway_col <- find_pathway_column(leading)

if (is.na(leading_pathway_col)) {
  stop("Could not identify pathway column.")
}

cat(
  "Leading-edge pathway column:",
  leading_pathway_col,
  "\n"
)

leading_std <- leading_std %>%
  mutate(
    Pathway = as.character(
      .data[[leading_pathway_col]]
    )
  ) %>%
  filter(
    !is.na(Pathway),
    Pathway != ""
  )

# ------------------------------------------------------------
# 9. CORE-GENE PATHWAY INTERSECTION
# ------------------------------------------------------------

core_leading <- leading_std %>%
  semi_join(
    core_std %>% select(GeneSymbol),
    by = "GeneSymbol"
  )

cat(
  "\nCore-gene pathway associations:",
  nrow(core_leading),
  "\n"
)

# ------------------------------------------------------------
# 10. PATHWAY FREQUENCY PER CORE GENE
# ------------------------------------------------------------

gene_pathway_frequency <- core_leading %>%
  group_by(GeneSymbol) %>%
  summarise(
    Pathway_Count = n_distinct(Pathway),
    Pathways = paste(
      unique(Pathway),
      collapse = "; "
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(Pathway_Count))

write_csv(
  gene_pathway_frequency,
  file.path(
    out_dir,
    "PCOS_core_gene_pathway_centrality.csv"
  )
)

cat(
  "Saved: PCOS_core_gene_pathway_centrality.csv\n"
)

# ------------------------------------------------------------
# 11. TOP PATHWAY-CENTRAL GENES
# ------------------------------------------------------------

top_genes <- gene_pathway_frequency %>%
  slice_max(
    Pathway_Count,
    n = 30,
    with_ties = FALSE
  )

write_csv(
  top_genes,
  file.path(
    out_dir,
    "PCOS_top30_pathway_central_genes.csv"
  )
)

# ------------------------------------------------------------
# 12. PATHWAY CENTRALITY
# ------------------------------------------------------------

pathway_gene_frequency <- core_leading %>%
  group_by(Pathway) %>%
  summarise(
    Core_Gene_Count = n_distinct(GeneSymbol),
    Core_Genes = paste(
      unique(GeneSymbol),
      collapse = ";"
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(Core_Gene_Count))

write_csv(
  pathway_gene_frequency,
  file.path(
    out_dir,
    "PCOS_pathway_core_gene_centrality.csv"
  )
)

# ------------------------------------------------------------
# 13. CORE GENE × PATHWAY MATRIX
# ------------------------------------------------------------

matrix_df <- core_leading %>%
  distinct(GeneSymbol, Pathway) %>%
  mutate(Value = 1) %>%
  tidyr::pivot_wider(
    names_from = Pathway,
    values_from = Value,
    values_fill = 0
  )

write_csv(
  matrix_df,
  file.path(
    out_dir,
    "PCOS_core_gene_pathway_matrix.csv"
  )
)

# ------------------------------------------------------------
# 14. IDENTIFY MULTI-PATHWAY GENES
# ------------------------------------------------------------

multi_pathway <- gene_pathway_frequency %>%
  filter(Pathway_Count >= 5)

write_csv(
  multi_pathway,
  file.path(
    out_dir,
    "PCOS_multi_pathway_core_genes.csv"
  )
)

cat(
  "\nGenes present in >=5 pathways:",
  nrow(multi_pathway),
  "\n"
)

# ------------------------------------------------------------
# 15. FUNCTIONAL THEME INFORMATION
# ------------------------------------------------------------

theme_file <- file.path(
  gsea_dir,
  "PCOS_pathway_functional_theme_assignment.csv"
)

if (file.exists(theme_file)) {

  themes <- read_csv(
    theme_file,
    show_col_types = FALSE
  )

  cat(
    "\nFunctional theme file loaded:",
    nrow(themes),
    "rows\n"
  )

  print(names(themes))

  theme_pathway_col <- find_pathway_column(themes)

  if (!is.na(theme_pathway_col)) {

    core_theme <- core_leading %>%
      left_join(
        themes,
        by = setNames(
          theme_pathway_col,
          "Pathway"
        )
      )

    write_csv(
      core_theme,
      file.path(
        out_dir,
        "PCOS_core_gene_functional_theme_links.csv"
      )
    )

    cat(
      "Saved: PCOS_core_gene_functional_theme_links.csv\n"
    )
  }
}

# ------------------------------------------------------------
# 16. TOP CORE-GENE PLOT
# ------------------------------------------------------------

plot_df <- top_genes %>%
  arrange(Pathway_Count)

if (nrow(plot_df) > 0) {

  p1 <- ggplot(
    plot_df,
    aes(
      x = Pathway_Count,
      y = reorder(
        GeneSymbol,
        Pathway_Count
      )
    )
  ) +
    geom_col() +
    labs(
      title = "Top PCOS Pathway-Central Core Genes",
      x = "Number of Significant Pathways",
      y = "Gene"
    ) +
    theme_minimal()

  ggsave(
    file.path(
      out_dir,
      "PCOS_top30_pathway_central_genes.png"
    ),
    p1,
    width = 9,
    height = 8,
    dpi = 300
  )

  ggsave(
    file.path(
      out_dir,
      "PCOS_top30_pathway_central_genes.pdf"
    ),
    p1,
    width = 9,
    height = 8
  )
}

# ------------------------------------------------------------
# 17. TOP PATHWAYS
# ------------------------------------------------------------

top_pathways <- pathway_gene_frequency %>%
  slice_max(
    Core_Gene_Count,
    n = 30,
    with_ties = FALSE
  ) %>%
  arrange(Core_Gene_Count)

write_csv(
  top_pathways,
  file.path(
    out_dir,
    "PCOS_top30_core_gene_rich_pathways.csv"
  )
)

p2 <- ggplot(
  top_pathways,
  aes(
    x = Core_Gene_Count,
    y = reorder(
      Pathway,
      Core_Gene_Count
    )
  )
) +
  geom_col() +
  labs(
    title = "PCOS Core-Gene Enrichment Across Significant Pathways",
    x = "Number of Core Genes",
    y = "Pathway"
  ) +
  theme_minimal()

ggsave(
  file.path(
    out_dir,
    "PCOS_top30_core_gene_rich_pathways.png"
  ),
  p2,
  width = 11,
  height = 9,
  dpi = 300
)

ggsave(
  file.path(
    out_dir,
    "PCOS_top30_core_gene_rich_pathways.pdf"
  ),
  p2,
  width = 11,
  height = 9
)

# ------------------------------------------------------------
# 18. REGULATORY CANDIDATE INPUT
# ------------------------------------------------------------

write_csv(
  core_std %>%
    select(GeneSymbol),
  file.path(
    out_dir,
    "PCOS_core_genes_for_TF_enrichment.csv"
  )
)

write_csv(
  multi_pathway %>%
    select(GeneSymbol),
  file.path(
    out_dir,
    "PCOS_high_centrality_genes_for_TF_enrichment.csv"
  )
)

# ------------------------------------------------------------
# 19. MASTER SUMMARY
# ------------------------------------------------------------

summary_df <- data.frame(
  Metric = c(
    "Total core genes",
    "Core gene-pathway associations",
    "Unique pathways represented",
    "Genes in >=5 pathways",
    "Top gene pathway count",
    "Top pathway core-gene count"
  ),
  Value = c(
    nrow(core_std),
    nrow(core_leading),
    n_distinct(core_leading$Pathway),
    nrow(multi_pathway),
    max(
      gene_pathway_frequency$Pathway_Count,
      na.rm = TRUE
    ),
    max(
      pathway_gene_frequency$Core_Gene_Count,
      na.rm = TRUE
    )
  )
)

write_csv(
  summary_df,
  file.path(
    out_dir,
    "PCOS_REGULATORY_LAYER_SUMMARY.csv"
  )
)

# ------------------------------------------------------------
# 20. FINAL REPORT
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("        REGULATORY LAYER ANALYSIS COMPLETE\n")
cat("============================================================\n\n")

cat(
  "Core genes:",
  nrow(core_std),
  "\n"
)

cat(
  "Core gene-pathway associations:",
  nrow(core_leading),
  "\n"
)

cat(
  "Unique pathways represented:",
  n_distinct(core_leading$Pathway),
  "\n"
)

cat(
  "Genes in >=5 pathways:",
  nrow(multi_pathway),
  "\n"
)

cat(
  "\nTop pathway-central genes:\n"
)

print(
  head(
    top_genes,
    20
  )
)

cat(
  "\nOutput directory:\n",
  out_dir,
  "\n"
)

cat("\n============================================================\n")
cat("NEXT STAGE: TF / REGULATOR ENRICHMENT\n")
cat("============================================================\n")
