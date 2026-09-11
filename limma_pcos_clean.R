# ============================================================
# GSE168404 PCOS mRNA DIFFERENTIAL EXPRESSION ANALYSIS
# Clean limma pipeline
# ============================================================

cat("\n")
cat("============================================================\n")
cat(" GSE168404 PCOS mRNA DIFFERENTIAL EXPRESSION ANALYSIS\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Load package
# ------------------------------------------------------------

if (!requireNamespace("limma", quietly = TRUE)) {
    stop(
        "limma is not installed.\n",
        "Install it with BiocManager::install('limma')"
    )
}

library(limma)

# ------------------------------------------------------------
# 2. Input file
# ------------------------------------------------------------

input_file <- "GSE168404_mRNA_10samples_clean.csv"

if (!file.exists(input_file)) {
    stop(
        "Input file not found: ",
        input_file,
        "\nMake sure the R script is running in the dataset directory."
    )
}

# ------------------------------------------------------------
# 3. Read data
# ------------------------------------------------------------

dat <- read.csv(
    input_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

cat("Input dimensions:", nrow(dat), ncol(dat), "\n\n")

cat("Input columns:\n")
print(colnames(dat))

# ------------------------------------------------------------
# 4. Identify GeneID and GeneSymbol
# ------------------------------------------------------------

# Based on your actual file:
# id = GeneID
# Unnamed: 1 = GeneSymbol

if (!"id" %in% colnames(dat)) {
    stop("Column 'id' was not found.")
}

if (!"Unnamed: 1" %in% colnames(dat)) {
    stop("Column 'Unnamed: 1' was not found.")
}

gene_id <- as.character(dat$id)
gene_symbol <- as.character(dat$`Unnamed: 1`)

# Replace missing symbols with GeneID
gene_symbol[
    is.na(gene_symbol) |
    gene_symbol == "" |
    gene_symbol == "NA"
] <- gene_id[
    is.na(gene_symbol) |
    gene_symbol == "" |
    gene_symbol == "NA"
]

# ------------------------------------------------------------
# 5. Define samples
# ------------------------------------------------------------

control_samples <- c(
    "C1", "C2", "C3", "C4", "C5"
)

pcos_samples <- c(
    "P1", "P2", "P3", "P4", "P5"
)

sample_names <- c(
    control_samples,
    pcos_samples
)

# Check samples
missing_samples <- setdiff(sample_names, colnames(dat))

if (length(missing_samples) > 0) {
    stop(
        "Missing sample columns: ",
        paste(missing_samples, collapse = ", ")
    )
}

# ------------------------------------------------------------
# 6. Build expression matrix
# ------------------------------------------------------------

expr <- dat[, sample_names, drop = FALSE]

# Convert every expression column to numeric
expr <- as.data.frame(
    lapply(expr, function(x) as.numeric(as.character(x))),
    check.names = FALSE
)

expr <- as.matrix(expr)

# Check numeric values
if (any(is.na(expr))) {
    cat("\nWARNING: NA values detected after numeric conversion.\n")
    cat("Number of NA values:", sum(is.na(expr)), "\n")

    expr[is.na(expr)] <- 0
}

# ------------------------------------------------------------
# 7. Remove duplicate GeneIDs
# ------------------------------------------------------------

duplicate_ids <- duplicated(gene_id)

cat("\nDuplicate GeneIDs:", sum(duplicate_ids), "\n")

if (any(duplicate_ids)) {

    # Keep first occurrence
    keep_unique <- !duplicate_ids

    expr <- expr[keep_unique, , drop = FALSE]
    gene_id <- gene_id[keep_unique]
    gene_symbol <- gene_symbol[keep_unique]
}

# Make GeneIDs unique for row names
gene_id_unique <- make.unique(gene_id)

rownames(expr) <- gene_id_unique

# ------------------------------------------------------------
# 8. Log2 transformation
# ------------------------------------------------------------

expr_log2 <- log2(expr + 1)

cat("\nLog2 transformation completed.\n")

# ------------------------------------------------------------
# 9. Expression filtering
# ------------------------------------------------------------

# Keep genes with log2(expression + 1) > 1
# in at least 3 samples

keep <- rowSums(expr_log2 > 1) >= 3

expr_filtered <- expr_log2[keep, , drop = FALSE]

gene_id_filtered <- gene_id[keep]
gene_symbol_filtered <- gene_symbol[keep]

cat("\nGenes before filtering:", nrow(expr_log2), "\n")
cat("Genes after filtering:", nrow(expr_filtered), "\n")

# ------------------------------------------------------------
# 10. Sample information
# ------------------------------------------------------------

group <- factor(
    c(
        rep("Control", length(control_samples)),
        rep("PCOS", length(pcos_samples))
    ),
    levels = c("Control", "PCOS")
)

sample_info <- data.frame(
    Sample = sample_names,
    Group = group
)

print(sample_info)

# ------------------------------------------------------------
# 11. PCA
# ------------------------------------------------------------

pca <- prcomp(
    t(expr_filtered),
    scale. = TRUE
)

pca_percent <- 100 * (pca$sdev^2 / sum(pca$sdev^2))

pca_df <- data.frame(
    Sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    Group = group
)

pdf(
    "PCOS_PCA.pdf",
    width = 7,
    height = 6
)

plot(
    pca_df$PC1,
    pca_df$PC2,
    pch = 19,
    cex = 1.5,
    xlab = paste0(
        "PC1 (",
        round(pca_percent[1], 1),
        "%)"
    ),
    ylab = paste0(
        "PC2 (",
        round(pca_percent[2], 1),
        "%)"
    ),
    main = "GSE168404 PCOS vs Control PCA"
)

text(
    pca_df$PC1,
    pca_df$PC2,
    labels = pca_df$Sample,
    pos = 3,
    cex = 0.8
)

dev.off()

# ------------------------------------------------------------
# 12. Sample correlation
# ------------------------------------------------------------

cor_matrix <- cor(
    expr_filtered,
    method = "pearson"
)

pdf(
    "PCOS_sample_correlation.pdf",
    width = 7,
    height = 6
)

heatmap(
    cor_matrix,
    main = "Sample Pearson Correlation",
    margins = c(8, 8)
)

dev.off()

# ------------------------------------------------------------
# 13. Limma design matrix
# ------------------------------------------------------------

design <- model.matrix(
    ~ group
)

colnames(design) <- c(
    "Intercept",
    "PCOS_vs_Control"
)

cat("\nDesign matrix:\n")
print(design)

# ------------------------------------------------------------
# 14. Fit limma model
# ------------------------------------------------------------

fit <- lmFit(
    expr_filtered,
    design
)

fit <- eBayes(
    fit
)

# ------------------------------------------------------------
# 15. Extract differential expression
# ------------------------------------------------------------

results_limma <- topTable(
    fit,
    coef = "PCOS_vs_Control",
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
)

cat("\nLimma result columns:\n")
print(colnames(results_limma))

# ------------------------------------------------------------
# 16. Add GeneID
# ------------------------------------------------------------

# topTable retains the expression matrix row names.
# Use those directly as GeneID.

results_limma$GeneID <- rownames(results_limma)

# ------------------------------------------------------------
# 17. Add GeneSymbol safely
# ------------------------------------------------------------

symbol_match <- match(
    results_limma$GeneID,
    make.unique(gene_id_filtered)
)

results_limma$GeneSymbol <- gene_symbol_filtered[
    symbol_match
]

# If annotation is missing, use GeneID
missing_symbol <- is.na(results_limma$GeneSymbol) |
                  results_limma$GeneSymbol == ""

results_limma$GeneSymbol[missing_symbol] <-
    results_limma$GeneID[missing_symbol]

# ------------------------------------------------------------
# 18. Add significance classification
# ------------------------------------------------------------

results_limma$Significance <- "Not Significant"

results_limma$Significance[
    results_limma$adj.P.Val < 0.05 &
    results_limma$logFC > 1
] <- "Upregulated"

results_limma$Significance[
    results_limma$adj.P.Val < 0.05 &
    results_limma$logFC < -1
] <- "Downregulated"

# ------------------------------------------------------------
# 19. Reorder columns
# ------------------------------------------------------------

results_final <- results_limma[
    ,
    c(
        "GeneID",
        "GeneSymbol",
        "logFC",
        "AveExpr",
        "t",
        "P.Value",
        "adj.P.Val",
        "B",
        "Significance"
    ),
    drop = FALSE
]

# ------------------------------------------------------------
# 20. Save complete results
# ------------------------------------------------------------

write.csv(
    results_final,
    "PCOS_limma_results.csv",
    row.names = FALSE
)

# ------------------------------------------------------------
# 21. Significant genes
# ------------------------------------------------------------

significant <- results_final[
    !is.na(results_final$adj.P.Val) &
    results_final$adj.P.Val < 0.05,
    ,
    drop = FALSE
]

upregulated <- significant[
    significant$logFC > 1,
    ,
    drop = FALSE
]

downregulated <- significant[
    significant$logFC < -1,
    ,
    drop = FALSE
]

write.csv(
    significant,
    "PCOS_significant_genes.csv",
    row.names = FALSE
)

write.csv(
    upregulated,
    "PCOS_upregulated_genes.csv",
    row.names = FALSE
)

write.csv(
    downregulated,
    "PCOS_downregulated_genes.csv",
    row.names = FALSE
)

# ------------------------------------------------------------
# 22. Volcano plot
# ------------------------------------------------------------

volcano_x <- results_final$logFC

volcano_y <- -log10(
    pmax(
        results_final$adj.P.Val,
        1e-300
    )
)

plot(
    volcano_x,
    volcano_y,
    pch = 19,
    cex = 0.5,
    col = "grey70",
    xlab = "Log2 Fold Change",
    ylab = "-Log10 Adjusted P-value",
    main = "GSE168404 PCOS Differential Expression"
)

# Highlight upregulated
up_idx <- which(
    results_final$adj.P.Val < 0.05 &
    results_final$logFC > 1
)

points(
    volcano_x[up_idx],
    volcano_y[up_idx],
    pch = 19,
    cex = 0.7,
    col = "red"
)

# Highlight downregulated
down_idx <- which(
    results_final$adj.P.Val < 0.05 &
    results_final$logFC < -1
)

points(
    volcano_x[down_idx],
    volcano_y[down_idx],
    pch = 19,
    cex = 0.7,
    col = "blue"
)

abline(
    v = c(-1, 1),
    lty = 2
)

abline(
    h = -log10(0.05),
    lty = 2
)

pdf(
    "PCOS_limma_volcano.pdf",
    width = 7,
    height = 6
)

plot(
    volcano_x,
    volcano_y,
    pch = 19,
    cex = 0.5,
    col = "grey70",
    xlab = "Log2 Fold Change",
    ylab = "-Log10 Adjusted P-value",
    main = "GSE168404 PCOS Differential Expression"
)

points(
    volcano_x[up_idx],
    volcano_y[up_idx],
    pch = 19,
    cex = 0.7,
    col = "red"
)

points(
    volcano_x[down_idx],
    volcano_y[down_idx],
    pch = 19,
    cex = 0.7,
    col = "blue"
)

abline(
    v = c(-1, 1),
    lty = 2
)

abline(
    h = -log10(0.05),
    lty = 2
)

dev.off()

# ------------------------------------------------------------
# 23. Top genes
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat(" DIFFERENTIAL EXPRESSION SUMMARY\n")
cat("============================================================\n\n")

cat(
    "Total genes tested:",
    nrow(results_final),
    "\n"
)

cat(
    "FDR < 0.05:",
    nrow(significant),
    "\n"
)

cat(
    "Upregulated (FDR < 0.05, logFC > 1):",
    nrow(upregulated),
    "\n"
)

cat(
    "Downregulated (FDR < 0.05, logFC < -1):",
    nrow(downregulated),
    "\n\n"
)

cat("Top 20 genes by adjusted P-value:\n\n")

print(
    head(
        results_final[
            order(results_final$adj.P.Val),
            c(
                "GeneSymbol",
                "logFC",
                "P.Value",
                "adj.P.Val",
                "Significance"
            )
        ],
        20
    )
)

# ------------------------------------------------------------
# 24. Output files
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat(" ANALYSIS COMPLETE\n")
cat("============================================================\n\n")

cat("Generated files:\n")
cat("  PCOS_limma_results.csv\n")
cat("  PCOS_significant_genes.csv\n")
cat("  PCOS_upregulated_genes.csv\n")
cat("  PCOS_downregulated_genes.csv\n")
cat("  PCOS_PCA.pdf\n")
cat("  PCOS_sample_correlation.pdf\n")
cat("  PCOS_limma_volcano.pdf\n\n")

cat("Controls: C1 C2 C3 C4 C5\n")
cat("PCOS:     P1 P2 P3 P4 P5\n\n")

cat("============================================================\n")
