# ============================================================
# PCOS GSEA CORE BIOLOGY + LEADING-EDGE ANALYSIS
# GSE168404
# ============================================================

cat("\n============================================================\n")
cat(" PCOS GSEA CORE BIOLOGY ANALYSIS\n")
cat("============================================================\n\n")

# -----------------------------
# 1. REQUIRED PACKAGES
# -----------------------------

required <- c(
  "dplyr",
  "readr",
  "stringr",
  "tidyr",
  "ggplot2"
)

for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)

# -----------------------------
# 2. INPUT FILES
# -----------------------------

gsea_file <- "PCOS_GSEA_GO_BP_results.csv"
leading_file <- "PCOS_GSEA_leading_edge_genes.csv"

if (!file.exists(gsea_file)) {
  stop(
    paste0(
      "\nERROR: Cannot find ", gsea_file,
      "\nMake sure it is in the current dataset folder.\n"
    )
  )
}

if (!file.exists(leading_file)) {
  stop(
    paste0(
      "\nERROR: Cannot find ", leading_file,
      "\nMake sure it is in the current dataset folder.\n"
    )
  )
}

# -----------------------------
# 3. LOAD GSEA RESULTS
# -----------------------------

gsea <- read_csv(gsea_file, show_col_types = FALSE)

cat("GSEA rows loaded:", nrow(gsea), "\n")
cat("GSEA columns:\n")
print(names(gsea))

# -----------------------------
# 4. STANDARDISE COLUMN NAMES
# -----------------------------

# Handle common naming differences

if ("FDR q-val" %in% names(gsea)) {
  gsea <- gsea %>%
    rename(FDR = `FDR q-val`)
}

if ("padj" %in% names(gsea) && !"FDR" %in% names(gsea)) {
  gsea <- gsea %>%
    rename(FDR = padj)
}

if (!"Term" %in% names(gsea)) {
  stop("ERROR: GSEA results do not contain a 'Term' column.")
}

if (!"NES" %in% names(gsea)) {
  stop("ERROR: GSEA results do not contain an 'NES' column.")
}

# -----------------------------
# 5. SIGNIFICANT PATHWAYS
# -----------------------------

sig <- gsea %>%
  filter(!is.na(FDR), FDR < 0.05) %>%
  mutate(
    Direction = ifelse(NES > 0, "Positive", "Negative"),
    Abs_NES = abs(NES)
  ) %>%
  arrange(desc(Abs_NES))

cat("\n============================================================\n")
cat(" SIGNIFICANT PATHWAYS\n")
cat("============================================================\n\n")

cat("Total pathways tested:", nrow(gsea), "\n")
cat("FDR < 0.05:", nrow(sig), "\n")
cat("Positive pathways:", sum(sig$NES > 0), "\n")
cat("Negative pathways:", sum(sig$NES < 0), "\n\n")

# -----------------------------
# 6. SAVE ALL SIGNIFICANT PATHWAYS
# -----------------------------

write_csv(
  sig,
  "PCOS_GSEA_significant_pathways.csv"
)

# -----------------------------
# 7. TOP POSITIVE / NEGATIVE
# -----------------------------

top_positive <- sig %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  slice_head(n = 20)

top_negative <- sig %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  slice_head(n = 20)

cat("\nTOP POSITIVE PATHWAYS\n")
print(
  top_positive %>%
    select(Term, NES, FDR) %>%
    as.data.frame()
)

cat("\nTOP NEGATIVE PATHWAYS\n")
print(
  top_negative %>%
    select(Term, NES, FDR) %>%
    as.data.frame()
)

write_csv(
  top_positive,
  "PCOS_GSEA_top20_positive.csv"
)

write_csv(
  top_negative,
  "PCOS_GSEA_top20_negative.csv"
)

# -----------------------------
# 8. LOAD LEADING EDGE GENES
# -----------------------------

leading <- read_csv(
  leading_file,
  show_col_types = FALSE
)

cat("\n============================================================\n")
cat(" LEADING-EDGE DATA\n")
cat("============================================================\n\n")

cat("Rows:", nrow(leading), "\n")
cat("Columns:\n")
print(names(leading))

# -----------------------------
# 9. FLEXIBLE LEADING-EDGE PARSER
# -----------------------------

# Detect pathway column

term_col <- names(leading)[
  names(leading) %in% c(
    "Term",
    "term",
    "Pathway",
    "pathway",
    "Description"
  )
][1]

if (is.na(term_col)) {
  stop(
    "\nERROR: Could not identify pathway column in leading-edge file."
  )
}

# Detect gene column

gene_col <- names(leading)[
  names(leading) %in% c(
    "Gene",
    "gene",
    "GeneSymbol",
    "gene_symbol",
    "GeneID",
    "gene_id"
  )
][1]

if (is.na(gene_col)) {

  # Try common multi-gene column

  gene_col <- names(leading)[
    str_detect(
      tolower(names(leading)),
      "lead|gene"
    )
  ][1]
}

if (is.na(gene_col)) {
  stop(
    "\nERROR: Could not identify gene column in leading-edge file."
  )
}

cat("\nPathway column:", term_col, "\n")
cat("Gene column:", gene_col, "\n")

# -----------------------------
# 10. CREATE PATHWAY-GENE TABLE
# -----------------------------

leading_clean <- leading %>%
  transmute(
    Term = .data[[term_col]],
    Gene = .data[[gene_col]]
  ) %>%
  filter(
    !is.na(Term),
    !is.na(Gene),
    Gene != ""
  )

# If genes are semicolon-separated,
# split them automatically

leading_long <- leading_clean %>%
  separate_rows(
    Gene,
    sep = "[;,|]"
  ) %>%
  mutate(
    Term = str_trim(Term),
    Gene = str_trim(Gene)
  ) %>%
  filter(
    Gene != "",
    !is.na(Gene)
  ) %>%
  distinct()

cat("\nUnique pathway-gene relationships:",
    nrow(leading_long), "\n")

# -----------------------------
# 11. COUNT GENE PARTICIPATION
# -----------------------------

gene_frequency <- leading_long %>%
  count(
    Gene,
    name = "Pathway_Count"
  ) %>%
  arrange(desc(Pathway_Count))

cat("\n============================================================\n")
cat(" CORE LEADING-EDGE GENES\n")
cat("============================================================\n\n")

print(
  head(gene_frequency, 30)
)

write_csv(
  gene_frequency,
  "PCOS_core_leading_edge_genes.csv"
)

# -----------------------------
# 12. TOP CORE GENES
# -----------------------------

top_core_genes <- gene_frequency %>%
  filter(Pathway_Count >= 3)

write_csv(
  top_core_genes,
  "PCOS_core_genes_shared_pathways.csv"
)

cat(
  "\nGenes appearing in >=3 significant pathways:",
  nrow(top_core_genes),
  "\n"
)

# -----------------------------
# 13. PATHWAY REDUNDANCY
# -----------------------------

pathway_gene_counts <- leading_long %>%
  count(
    Term,
    name = "Leading_Edge_Genes"
  ) %>%
  arrange(desc(Leading_Edge_Genes))

write_csv(
  pathway_gene_counts,
  "PCOS_pathway_leading_edge_counts.csv"
)

# -----------------------------
# 14. METABOLIC PATHWAY SCREEN
# -----------------------------

metabolic_keywords <- paste(
  c(
    "mitochond",
    "oxidative phosphorylation",
    "electron transport",
    "respiration",
    "ATP",
    "pyruvate",
    "carbohydrate",
    "nucleotide",
    "purine",
    "metabolic",
    "lipid",
    "fatty acid",
    "cholesterol",
    "sterol",
    "energy"
  ),
  collapse = "|"
)

metabolic <- sig %>%
  filter(
    str_detect(
      tolower(Term),
      metabolic_keywords
    )
  ) %>%
  arrange(NES)

cat("\n============================================================\n")
cat(" METABOLIC / MITOCHONDRIAL PATHWAYS\n")
cat("============================================================\n\n")

cat("Matching pathways:", nrow(metabolic), "\n\n")

print(
  metabolic %>%
    select(Term, NES, FDR) %>%
    head(50) %>%
    as.data.frame()
)

write_csv(
  metabolic,
  "PCOS_metabolic_mitochondrial_pathways.csv"
)

# -----------------------------
# 15. PATHWAY THEMES
# -----------------------------

sig <- sig %>%
  mutate(
    Theme = case_when(

      str_detect(
        tolower(Term),
        "mitochond|oxidative phosphorylation|electron transport|respiration|atp"
      )
      ~ "Mitochondrial / Energy Metabolism",

      str_detect(
        tolower(Term),
        "purine|nucleotide|nucleoside|ribose phosphate"
      )
      ~ "Nucleotide Metabolism",

      str_detect(
        tolower(Term),
        "pyruvate|carbohydrate|glucose|monocarboxylic|glycol"
      )
      ~ "Carbohydrate / Central Metabolism",

      str_detect(
        tolower(Term),
        "lipid|fatty acid|cholesterol|sterol"
      )
      ~ "Lipid / Sterol Metabolism",

      str_detect(
        tolower(Term),
        "translation|ribosome|ribonucleoprotein"
      )
      ~ "Translation / Ribosome",

      str_detect(
        tolower(Term),
        "protein folding|protein targeting|protein localization"
      )
      ~ "Protein Processing",

      TRUE
      ~ "Other"
    )
  )

# -----------------------------
# 16. THEME SUMMARY
# -----------------------------

theme_summary <- sig %>%
  group_by(
    Theme,
    Direction
  ) %>%
  summarise(
    Pathways = n(),
    Mean_NES = mean(NES),
    Strongest_NES = ifelse(
      Direction == "Positive",
      max(NES),
      min(NES)
    ),
    .groups = "drop"
  ) %>%
  arrange(
    Theme,
    desc(Pathways)
  )

write_csv(
  theme_summary,
  "PCOS_GSEA_theme_summary.csv"
)

cat("\n============================================================\n")
cat(" BIOLOGICAL THEME SUMMARY\n")
cat("============================================================\n\n")

print(theme_summary)

# -----------------------------
# 17. THEME PLOT
# -----------------------------

theme_plot <- sig %>%
  count(Theme, Direction) %>%
  ggplot(
    aes(
      x = reorder(Theme, n),
      y = n,
      fill = Direction
    )
  ) +
  geom_col() +
  coord_flip() +
  labs(
    title = "PCOS GSEA Significant Pathway Themes",
    x = "Biological Theme",
    y = "Number of Significant Pathways"
  ) +
  theme_minimal()

ggsave(
  "PCOS_GSEA_theme_summary.png",
  theme_plot,
  width = 10,
  height = 6,
  dpi = 300
)

# -----------------------------
# 18. TOP PATHWAYS PLOT
# -----------------------------

top_plot_data <- sig %>%
  arrange(desc(Abs_NES)) %>%
  slice_head(n = 25) %>%
  mutate(
    Term = factor(
      Term,
      levels = rev(Term)
    )
  )

p <- ggplot(
  top_plot_data,
  aes(
    x = NES,
    y = Term
  )
) +
  geom_col() +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = "Top GSEA Pathways in PCOS",
    x = "Normalized Enrichment Score (NES)",
    y = "Biological Process"
  ) +
  theme_minimal()

ggsave(
  "PCOS_GSEA_top25_pathways.png",
  p,
  width = 11,
  height = 9,
  dpi = 300
)

# -----------------------------
# 19. FINAL SUMMARY
# -----------------------------

cat("\n============================================================\n")
cat(" ANALYSIS COMPLETE\n")
cat("============================================================\n\n")

cat("Generated files:\n\n")

cat("1. PCOS_GSEA_significant_pathways.csv\n")
cat("2. PCOS_GSEA_top20_positive.csv\n")
cat("3. PCOS_GSEA_top20_negative.csv\n")
cat("4. PCOS_core_leading_edge_genes.csv\n")
cat("5. PCOS_core_genes_shared_pathways.csv\n")
cat("6. PCOS_pathway_leading_edge_counts.csv\n")
cat("7. PCOS_metabolic_mitochondrial_pathways.csv\n")
cat("8. PCOS_GSEA_theme_summary.csv\n")
cat("9. PCOS_GSEA_theme_summary.png\n")
cat("10. PCOS_GSEA_top25_pathways.png\n\n")

cat("NEXT STEP:\n")
cat("Use the core leading-edge genes for regulatory-network analysis.\n\n")
