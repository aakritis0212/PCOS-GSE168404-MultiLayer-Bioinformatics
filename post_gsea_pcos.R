############################################################
# PCOS POST-GSEA FUNCTIONAL ANALYSIS
# GSE168404 - Granulosa Cell Transcriptomics
#
# Input:
#   PCOS_GSEA_GO_BP_results.csv
#
# Outputs:
#   - Pathway summaries
#   - Functional theme summaries
#   - Leading-edge gene analysis
#   - Core genes
#   - Pathway-gene network
#   - Publication-quality PNG/PDF plots
#   - CSV tables
############################################################

cat("\n")
cat("============================================================\n")
cat("       PCOS POST-GSEA FUNCTIONAL ANALYSIS\n")
cat("============================================================\n\n")


############################################################
# 1. PACKAGE MANAGEMENT
############################################################

required_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "stringr",
  "forcats",
  "readr",
  "igraph"
)

missing_packages <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_packages) > 0) {

  cat("Missing packages:\n")
  print(missing_packages)

  cat("\nInstalling missing packages...\n")

  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org"
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(readr)
  library(igraph)
})


############################################################
# 2. CREATE OUTPUT DIRECTORY
############################################################

output_dir <- "PCOS_post_GSEA_results"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

cat("Output directory:", output_dir, "\n\n")


############################################################
# 3. LOAD GSEA RESULTS
############################################################

gsea_file <- "PCOS_GSEA_GO_BP_results.csv"

if (!file.exists(gsea_file)) {

  stop(
    paste0(
      "\nERROR: Cannot find ", gsea_file,
      "\nMake sure this script is inside:\n",
      getwd(), "\n"
    )
  )
}

cat("Loading GSEA results...\n")

gsea <- read_csv(
  gsea_file,
  show_col_types = FALSE
)

cat("GSEA rows:", nrow(gsea), "\n")
cat("GSEA columns:\n")
print(colnames(gsea))


############################################################
# 4. VERIFY REQUIRED COLUMNS
############################################################

required_columns <- c(
  "Term",
  "NES",
  "pval",
  "padj",
  "size",
  "LeadingEdgeGenes"
)

missing_columns <- setdiff(
  required_columns,
  colnames(gsea)
)

if (length(missing_columns) > 0) {

  stop(
    paste(
      "\nERROR: Required columns are missing:",
      paste(missing_columns, collapse = ", ")
    )
  )
}

cat("\nRequired GSEA columns successfully detected.\n")


############################################################
# 5. CLEAN GSEA DATA
############################################################

gsea <- gsea %>%

  mutate(
    Term = as.character(Term),
    NES = as.numeric(NES),
    pval = as.numeric(pval),
    padj = as.numeric(padj),
    size = as.numeric(size),
    LeadingEdgeGenes = as.character(LeadingEdgeGenes)
  ) %>%

  filter(
    !is.na(Term),
    !is.na(NES),
    !is.na(padj)
  )


############################################################
# 6. DEFINE SIGNIFICANT PATHWAYS
############################################################

sig <- gsea %>%
  filter(padj < 0.05) %>%
  mutate(
    Direction = case_when(
      NES > 0 ~ "Positive",
      NES < 0 ~ "Negative",
      TRUE ~ "Neutral"
    ),
    AbsNES = abs(NES)
  )

cat("\n")
cat("============================================================\n")
cat("                    GSEA SUMMARY\n")
cat("============================================================\n")

cat("Total valid pathways:", nrow(gsea), "\n")
cat("Significant pathways:", nrow(sig), "\n")
cat(
  "Positive pathways:",
  sum(sig$Direction == "Positive"),
  "\n"
)
cat(
  "Negative pathways:",
  sum(sig$Direction == "Negative"),
  "\n"
)


############################################################
# 7. SAVE COMPLETE SIGNIFICANT PATHWAY TABLE
############################################################

write_csv(
  sig,
  file.path(
    output_dir,
    "PCOS_significant_GSEA_pathways.csv"
  )
)


############################################################
# 8. TOP POSITIVE PATHWAYS
############################################################

positive_pathways <- sig %>%
  filter(NES > 0) %>%
  arrange(desc(NES))

write_csv(
  positive_pathways,
  file.path(
    output_dir,
    "PCOS_positive_GSEA_pathways.csv"
  )
)


############################################################
# 9. TOP NEGATIVE PATHWAYS
############################################################

negative_pathways <- sig %>%
  filter(NES < 0) %>%
  arrange(NES)

write_csv(
  negative_pathways,
  file.path(
    output_dir,
    "PCOS_negative_GSEA_pathways.csv"
  )
)


############################################################
# 10. TOP 20 PATHWAYS BY ABSOLUTE NES
############################################################

top20 <- sig %>%
  arrange(desc(AbsNES)) %>%
  slice_head(n = 20)

write_csv(
  top20,
  file.path(
    output_dir,
    "PCOS_top20_GSEA_pathways.csv"
  )
)


############################################################
# 11. GSEA TOP PATHWAY PLOT
############################################################

plot_data <- sig %>%
  arrange(desc(AbsNES)) %>%
  slice_head(n = 20) %>%
  mutate(
    Term = fct_reorder(Term, NES)
  )

p1 <- ggplot(
  plot_data,
  aes(
    x = NES,
    y = Term
  )
) +
  geom_col(
    aes(fill = NES > 0)
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  scale_fill_manual(
    values = c(
      "TRUE" = "#D95F02",
      "FALSE" = "#1B9E77"
    ),
    labels = c(
      "FALSE" = "Negative NES",
      "TRUE" = "Positive NES"
    ),
    name = "Direction"
  ) +
  labs(
    title = "Top Significant PCOS Biological Processes",
    subtitle = "GSEA pathways ranked by absolute Normalized Enrichment Score",
    x = "Normalized Enrichment Score (NES)",
    y = "GO Biological Process"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 8)
  )

ggsave(
  file.path(
    output_dir,
    "PCOS_top20_GSEA_NES.png"
  ),
  p1,
  width = 11,
  height = 8,
  dpi = 300
)

ggsave(
  file.path(
    output_dir,
    "PCOS_top20_GSEA_NES.pdf"
  ),
  p1,
  width = 11,
  height = 8
)


############################################################
# 12. PATHWAY BUBBLE PLOT
############################################################

bubble_data <- sig %>%
  arrange(desc(AbsNES)) %>%
  slice_head(n = 25) %>%
  mutate(
    Term = fct_reorder(Term, NES)
  )

p2 <- ggplot(
  bubble_data,
  aes(
    x = NES,
    y = Term,
    size = size,
    color = -log10(padj)
  )
) +
  geom_point(alpha = 0.85) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  scale_color_viridis_c() +
  labs(
    title = "PCOS GSEA Pathway Bubble Plot",
    x = "Normalized Enrichment Score (NES)",
    y = "Biological Process",
    size = "Gene-set size",
    color = "-log10(FDR)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 8)
  )

ggsave(
  file.path(
    output_dir,
    "PCOS_GSEA_bubble_plot.png"
  ),
  p2,
  width = 11,
  height = 9,
  dpi = 300
)

ggsave(
  file.path(
    output_dir,
    "PCOS_GSEA_bubble_plot.pdf"
  ),
  p2,
  width = 11,
  height = 9
)


############################################################
# 13. FUNCTIONAL THEME SUMMARISATION
############################################################

cat("\n")
cat("============================================================\n")
cat("             FUNCTIONAL THEME SUMMARISATION\n")
cat("============================================================\n")


# Keyword-based biological theme classification.
# A pathway can belong to more than one biological theme.

assign_theme <- function(term) {

  term_low <- tolower(term)

  themes <- c()

  if (
    str_detect(
      term_low,
      "mitochond|oxidative phosphorylation|electron transport|cellular respiration|atp|energy"
    )
  ) {
    themes <- c(themes, "Mitochondrial / Energy Metabolism")
  }

  if (
    str_detect(
      term_low,
      "nucleotide|purine|pyrimidine|nucleoside|ribose phosphate"
    )
  ) {
    themes <- c(themes, "Nucleotide Metabolism")
  }

  if (
    str_detect(
      term_low,
      "carbohydrate|glucose|glycolysis|pyruvate|monocarboxylic|lipid|fatty acid"
    )
  ) {
    themes <- c(themes, "Metabolic Reprogramming")
  }

  if (
    str_detect(
      term_low,
      "ribosome|translation|ribosomal|ribonucleoprotein"
    )
  ) {
    themes <- c(themes, "Protein Synthesis / Translation")
  }

  if (
    str_detect(
      term_low,
      "immune|inflamm|cytokine|leukocyte|lymphocyte|interferon"
    )
  ) {
    themes <- c(themes, "Immune / Inflammatory Signaling")
  }

  if (
    str_detect(
      term_low,
      "cell cycle|mitotic|chromosome|dna replication|replication"
    )
  ) {
    themes <- c(themes, "Cell Cycle / Genome Maintenance")
  }

  if (
    str_detect(
      term_low,
      "axon|neuron|synapse|neuronal|cell morphogenesis"
    )
  ) {
    themes <- c(themes, "Neuronal / Morphogenesis Processes")
  }

  if (
    str_detect(
      term_low,
      "transcription|rna processing|mrna|gene expression|chromatin"
    )
  ) {
    themes <- c(themes, "Gene Expression / RNA Regulation")
  }

  if (
    str_detect(
      term_low,
      "protein folding|proteolysis|ubiquitin|protein modification"
    )
  ) {
    themes <- c(themes, "Protein Homeostasis")
  }

  if (
    str_detect(
      term_low,
      "membrane|transport|vesicle|endocyt"
    )
  ) {
    themes <- c(themes, "Membrane / Cellular Transport")
  }

  if (
    str_detect(
      term_low,
      "development|differentiation|morphogenesis"
    )
  ) {
    themes <- c(themes, "Development / Differentiation")
  }

  if (length(themes) == 0) {
    themes <- "Other Biological Processes"
  }

  unique(themes)
}


############################################################
# 14. APPLY FUNCTIONAL THEMES
############################################################

theme_table <- sig %>%
  rowwise() %>%
  mutate(
    Theme = list(assign_theme(Term))
  ) %>%
  ungroup() %>%
  unnest(Theme)


############################################################
# 15. THEME SUMMARY
############################################################

theme_summary <- theme_table %>%
  group_by(Theme) %>%
  summarise(
    Significant_Pathways = n(),
    Mean_NES = mean(NES, na.rm = TRUE),
    Median_NES = median(NES, na.rm = TRUE),
    Strongest_NES = max(abs(NES), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Significant_Pathways))


write_csv(
  theme_summary,
  file.path(
    output_dir,
    "PCOS_functional_theme_summary.csv"
  )
)

write_csv(
  theme_table,
  file.path(
    output_dir,
    "PCOS_pathway_functional_theme_assignment.csv"
  )
)

cat("\nFunctional themes identified:\n")
print(theme_summary)


############################################################
# 16. FUNCTIONAL THEME PLOT
############################################################

theme_plot_data <- theme_summary %>%
  slice_head(n = 12) %>%
  mutate(
    Theme = fct_reorder(
      Theme,
      Significant_Pathways
    )
  )

p3 <- ggplot(
  theme_plot_data,
  aes(
    x = Significant_Pathways,
    y = Theme
  )
) +
  geom_col() +
  labs(
    title = "Major Functional Themes in PCOS",
    subtitle = "Number of significantly enriched GO Biological Processes",
    x = "Number of Significant Pathways",
    y = "Functional Theme"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(
    output_dir,
    "PCOS_functional_theme_summary.png"
  ),
  p3,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(
    output_dir,
    "PCOS_functional_theme_summary.pdf"
  ),
  p3,
  width = 10,
  height = 7
)


############################################################
# 17. LEADING-EDGE ANALYSIS
############################################################

cat("\n")
cat("============================================================\n")
cat("                LEADING-EDGE ANALYSIS\n")
cat("============================================================\n")


# IMPORTANT:
# Explicitly use LeadingEdgeGenes.
# NEVER substitute Term.

leading_edge_col <- "LeadingEdgeGenes"

if (!leading_edge_col %in% colnames(sig)) {

  stop(
    "\nERROR: LeadingEdgeGenes column is missing.\n",
    "The script will NOT use Term as a fallback.\n"
  )
}

cat(
  "Confirmed leading-edge column:",
  leading_edge_col,
  "\n"
)


############################################################
# 18. EXTRACT LEADING-EDGE GENES
############################################################

leading_edge_long <- sig %>%

  select(
    Term,
    NES,
    padj,
    LeadingEdgeGenes
  ) %>%

  filter(
    !is.na(LeadingEdgeGenes),
    LeadingEdgeGenes != "",
    LeadingEdgeGenes != "NA"
  ) %>%

  mutate(
    LeadingEdgeGenes = str_replace_all(
      LeadingEdgeGenes,
      "\\[|\\]|'|\"",
      ""
    )
  ) %>%

  separate_rows(
    LeadingEdgeGenes,
    sep = "[;,]"
  ) %>%

  mutate(
    Gene = str_trim(LeadingEdgeGenes)
  ) %>%

  filter(
    Gene != "",
    Gene != "NA",
    !is.na(Gene)
  ) %>%

  select(
    Term,
    NES,
    padj,
    Gene
  )


############################################################
# 19. REMOVE DUPLICATE GENE-PATHWAY PAIRS
############################################################

leading_edge_long <- leading_edge_long %>%
  distinct(
    Term,
    Gene,
    .keep_all = TRUE
  )


cat(
  "Leading-edge pathway-gene associations:",
  nrow(leading_edge_long),
  "\n"
)

cat(
  "Unique leading-edge genes:",
  n_distinct(leading_edge_long$Gene),
  "\n"
)


############################################################
# 20. SAVE LEADING-EDGE TABLE
############################################################

write_csv(
  leading_edge_long,
  file.path(
    output_dir,
    "PCOS_leading_edge_gene_pathway_table.csv"
  )
)


############################################################
# 21. CORE GENE FREQUENCY
############################################################

core_genes <- leading_edge_long %>%

  group_by(Gene) %>%

  summarise(
    Pathway_Count = n_distinct(Term),
    Mean_NES = mean(NES, na.rm = TRUE),
    Strongest_NES = max(abs(NES), na.rm = TRUE),
    .groups = "drop"
  ) %>%

  arrange(desc(Pathway_Count), desc(Strongest_NES))


write_csv(
  core_genes,
  file.path(
    output_dir,
    "PCOS_leading_edge_gene_frequency.csv"
  )
)


############################################################
# 22. CORE GENES IN >= 3 PATHWAYS
############################################################

core_genes_3 <- core_genes %>%
  filter(Pathway_Count >= 3)


cat(
  "Core genes occurring in >=3 pathways:",
  nrow(core_genes_3),
  "\n"
)

write_csv(
  core_genes_3,
  file.path(
    output_dir,
    "PCOS_core_genes_3_or_more_pathways.csv"
  )
)


############################################################
# 23. TOP 30 CORE GENES
############################################################

top_core_genes <- core_genes %>%
  slice_head(n = 30)


write_csv(
  top_core_genes,
  file.path(
    output_dir,
    "PCOS_top30_core_genes.csv"
  )
)


############################################################
# 24. CORE GENE BARPLOT
############################################################

if (nrow(top_core_genes) > 0) {

  core_plot_data <- top_core_genes %>%
    mutate(
      Gene = fct_reorder(
        Gene,
        Pathway_Count
      )
    )

  p4 <- ggplot(
    core_plot_data,
    aes(
      x = Pathway_Count,
      y = Gene
    )
  ) +
    geom_col() +
    labs(
      title = "Top Leading-Edge Core Genes in PCOS",
      subtitle = "Genes recurring across enriched biological processes",
      x = "Number of Enriched Pathways",
      y = "Gene"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold")
    )

  ggsave(
    file.path(
      output_dir,
      "PCOS_top_core_genes.png"
    ),
    p4,
    width = 9,
    height = 8,
    dpi = 300
  )

  ggsave(
    file.path(
      output_dir,
      "PCOS_top_core_genes.pdf"
    ),
    p4,
    width = 9,
    height = 8
  )
}


############################################################
# 25. PATHWAY-GENE NETWORK
############################################################

cat("\n")
cat("============================================================\n")
cat("             PATHWAY-GENE NETWORK ANALYSIS\n")
cat("============================================================\n")


if (nrow(core_genes_3) > 0) {

  network_genes <- core_genes_3$Gene

  network_edges <- leading_edge_long %>%
    filter(
      Gene %in% network_genes
    ) %>%
    select(
      Term,
      Gene
    ) %>%
    distinct()


  # Keep only pathways connected to >=1 core gene

  network_pathways <- unique(
    network_edges$Term
  )

  network_edges <- network_edges %>%
    filter(
      Term %in% network_pathways
    )


  ##########################################################
  # Create unique node names
  ##########################################################

  pathway_nodes <- data.frame(
    name = unique(network_edges$Term),
    type = "Pathway",
    stringsAsFactors = FALSE
  )

  gene_nodes <- data.frame(
    name = unique(network_edges$Gene),
    type = "Gene",
    stringsAsFactors = FALSE
  )

  nodes <- bind_rows(
    pathway_nodes,
    gene_nodes
  ) %>%
    distinct(name, .keep_all = TRUE)


  ##########################################################
  # Edges
  ##########################################################

  edges <- network_edges %>%
    rename(
      from = Term,
      to = Gene
    )


  ##########################################################
  # SAFETY CHECK FOR DUPLICATED VERTICES
  ##########################################################

  nodes <- nodes %>%
    filter(
      !is.na(name),
      name != ""
    ) %>%
    distinct(
      name,
      .keep_all = TRUE
    )


  edges <- edges %>%
    filter(
      from %in% nodes$name,
      to %in% nodes$name
    ) %>%
    distinct(
      from,
      to,
      .keep_all = TRUE
    )


  ##########################################################
  # Build graph
  ##########################################################

  g <- graph_from_data_frame(
    d = edges,
    vertices = nodes,
    directed = FALSE
  )


  cat(
    "Network nodes:",
    vcount(g),
    "\n"
  )

  cat(
    "Network edges:",
    ecount(g),
    "\n"
  )


  ##########################################################
  # SAVE NETWORK EDGE TABLE
  ##########################################################

  write_csv(
    edges,
    file.path(
      output_dir,
      "PCOS_core_gene_pathway_network_edges.csv"
    )
  )


  write_csv(
    nodes,
    file.path(
      output_dir,
      "PCOS_core_gene_pathway_network_nodes.csv"
    )
  )


  ##########################################################
  # SAVE GRAPH OBJECT
  ##########################################################

  saveRDS(
    g,
    file.path(
      output_dir,
      "PCOS_core_gene_pathway_network.rds"
    )
  )


  ##########################################################
  # NETWORK PDF
  ##########################################################

  pdf(
    file.path(
      output_dir,
      "PCOS_core_gene_pathway_network.pdf"
    ),
    width = 12,
    height = 10
  )

  set.seed(123)

  plot(
    g,
    vertex.size = 5,
    vertex.label.cex = 0.55,
    vertex.label.color = "black",
    edge.arrow.size = 0.2,
    layout = layout_with_fr(g)
  )

  dev.off()


} else {

  cat(
    "\nNo genes occurred in >=3 pathways.\n"
  )

  cat(
    "Network was therefore not generated.\n"
  )
}


############################################################
# 26. TOP PATHWAY × CORE GENE MATRIX
############################################################

if (nrow(core_genes_3) > 0) {

  top_network_genes <- core_genes_3 %>%
    slice_head(n = 30) %>%
    pull(Gene)


  pathway_gene_matrix <- leading_edge_long %>%

    filter(
      Gene %in% top_network_genes
    ) %>%

    mutate(
      Present = 1
    ) %>%

    distinct(
      Term,
      Gene,
      .keep_all = TRUE
    ) %>%

    select(
      Term,
      Gene,
      Present
    ) %>%

    pivot_wider(
      names_from = Gene,
      values_from = Present,
      values_fill = 0
    )


  write_csv(
    pathway_gene_matrix,
    file.path(
      output_dir,
      "PCOS_pathway_core_gene_matrix.csv"
    )
  )
}


############################################################
# 27. DIRECTIONAL FUNCTIONAL SUMMARY
############################################################

direction_summary <- sig %>%

  group_by(Direction) %>%

  summarise(
    Pathways = n(),
    Mean_NES = mean(NES, na.rm = TRUE),
    Median_NES = median(NES, na.rm = TRUE),
    Mean_FDR = mean(padj, na.rm = TRUE),
    .groups = "drop"
  )


write_csv(
  direction_summary,
  file.path(
    output_dir,
    "PCOS_GSEA_direction_summary.csv"
  )
)


############################################################
# 28. MASTER SUMMARY TABLE
############################################################

master_summary <- data.frame(

  Metric = c(
    "Total pathways tested",
    "Significant pathways FDR < 0.05",
    "Positive pathways",
    "Negative pathways",
    "Unique leading-edge genes",
    "Core genes >=3 pathways"
  ),

  Value = c(
    nrow(gsea),
    nrow(sig),
    sum(sig$Direction == "Positive"),
    sum(sig$Direction == "Negative"),
    n_distinct(leading_edge_long$Gene),
    nrow(core_genes_3)
  )
)


write_csv(
  master_summary,
  file.path(
    output_dir,
    "PCOS_POST_GSEA_MASTER_SUMMARY.csv"
  )
)


############################################################
# 29. FINAL REPORT
############################################################

cat("\n")
cat("============================================================\n")
cat("                  POST-GSEA COMPLETE\n")
cat("============================================================\n\n")

cat(
  "Total pathways tested: ",
  nrow(gsea),
  "\n",
  sep = ""
)

cat(
  "Significant pathways: ",
  nrow(sig),
  "\n",
  sep = ""
)

cat(
  "Positive pathways: ",
  sum(sig$Direction == "Positive"),
  "\n",
  sep = ""
)

cat(
  "Negative pathways: ",
  sum(sig$Direction == "Negative"),
  "\n",
  sep = ""
)

cat(
  "Unique leading-edge genes: ",
  n_distinct(leading_edge_long$Gene),
  "\n",
  sep = ""
)

cat(
  "Core genes >=3 pathways: ",
  nrow(core_genes_3),
  "\n",
  sep = ""
)

cat("\n")
cat("Output directory:\n")
cat(
  normalizePath(output_dir),
  "\n\n"
)

cat("Important files generated:\n\n")

cat("1. PCOS_significant_GSEA_pathways.csv\n")
cat("2. PCOS_positive_GSEA_pathways.csv\n")
cat("3. PCOS_negative_GSEA_pathways.csv\n")
cat("4. PCOS_top20_GSEA_pathways.csv\n")
cat("5. PCOS_functional_theme_summary.csv\n")
cat("6. PCOS_pathway_functional_theme_assignment.csv\n")
cat("7. PCOS_leading_edge_gene_pathway_table.csv\n")
cat("8. PCOS_leading_edge_gene_frequency.csv\n")
cat("9. PCOS_core_genes_3_or_more_pathways.csv\n")
cat("10. PCOS_top30_core_genes.csv\n")
cat("11. PCOS_core_gene_pathway_network_edges.csv\n")
cat("12. PCOS_core_gene_pathway_network_nodes.csv\n")
cat("13. PCOS_pathway_core_gene_matrix.csv\n")
cat("14. PCOS_POST_GSEA_MASTER_SUMMARY.csv\n")

cat("\nPlots:\n\n")

cat("15. PCOS_top20_GSEA_NES.png/pdf\n")
cat("16. PCOS_GSEA_bubble_plot.png/pdf\n")
cat("17. PCOS_functional_theme_summary.png/pdf\n")
cat("18. PCOS_top_core_genes.png/pdf\n")
cat("19. PCOS_core_gene_pathway_network.pdf\n")

cat("\n============================================================\n")
cat("             ANALYSIS FINISHED SUCCESSFULLY\n")
cat("============================================================\n")