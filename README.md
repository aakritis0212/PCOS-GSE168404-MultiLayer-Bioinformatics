# Multi-Layer Bioinformatics Analysis of PCOS Transcriptomic Data

### GSE168404 | Differential Expression → GSEA → Core Genes → Regulatory Network

[![R](https://img.shields.io/badge/R-4.5.3-blue)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-compatible-green)](https://www.bioconductor.org/)
[![Streamlit](https://img.shields.io/badge/Interactive-Streamlit-red)](https://streamlit.io/)

---

## Overview

This repository contains a reproducible multi-layer bioinformatics analysis of transcriptomic data from **NCBI GEO accession GSE168404**, comparing control and PCOS samples.

The project integrates differential expression analysis, pathway-level enrichment, leading-edge analysis, core-gene prioritization and transcriptional regulator enrichment to identify biological processes and regulatory factors associated with PCOS.

---

## Study Design

| Group   | Samples    |
| ------- | ---------- |
| Control | C1–C5      |
| PCOS    | P1–P5      |
| Total   | 10 samples |
| Dataset | GSE168404  |

---

## Analytical Workflow

```text
GSE168404
    │
    ▼
Expression Matrix QC
    │
    ▼
Expression Filtering
    │
    ▼
limma Differential Expression
    │
    ▼
Ranked Gene List
    │
    ▼
Pre-ranked GSEA
    │
    ▼
Functional Theme Summarisation
    │
    ▼
Leading-edge Analysis
    │
    ▼
Core Gene Identification
    │
    ▼
Pathway Centrality
    │
    ▼
TRRUST Regulator Enrichment
    │
    ▼
Regulatory Network
```

---

## Main Results

### Differential Expression

The analysis retained 14,749 genes after expression filtering.

Differential expression was performed using the **limma** framework with a PCOS-versus-control design.

---

### GSEA

A pre-ranked gene-set enrichment analysis was performed using Gene Ontology Biological Process gene sets.

Key enriched biological themes included:

* mitochondrial gene expression
* mitochondrial translation
* oxidative phosphorylation
* cellular respiration
* nucleotide metabolism
* purine metabolism
* carbohydrate metabolism
* pyruvate metabolism
* developmental processes
* neuronal/morphogenetic processes

---

### Leading-Edge Analysis

The GSEA analysis identified:

* **446 significant pathways**
* **2,956 unique leading-edge genes**
* **1,774 core genes occurring in ≥3 pathways**

These recurrent genes were subsequently used for pathway-centrality and regulatory analysis.

---

### Regulatory Analysis

TRRUST enrichment identified statistically significant regulatory associations among the PCOS core-gene set.

The strongest FDR-significant regulators were:

| Regulator | Target overlap | Odds ratio |     FDR |
| --------- | -------------: | ---------: | ------: |
| TP53      |             51 |       2.44 | 0.00078 |
| HIF1A     |             29 |       2.83 | 0.00723 |
| HDAC1     |             24 |       2.67 |  0.0423 |

These regulators represent candidate upstream regulatory nodes associated with the transcriptional architecture identified in this dataset.

---

## Reproducibility

All major computational stages are implemented as R scripts.

Run the scripts in order:

```bash
Rscript scripts/01_QC.R
Rscript scripts/02_limma_DEG.R
Rscript scripts/03_GSEA.R
Rscript scripts/04_post_GSEA_functional_analysis.R
Rscript scripts/05_core_gene_regulatory_analysis.R
Rscript scripts/06_TRRUST_regulator_enrichment.R
Rscript scripts/07_final_regulatory_analysis.R
```

---

## Software

The analysis uses R/Bioconductor and associated packages including:

* limma
* fgsea
* dplyr
* ggplot2
* igraph
* ggraph
* forcats

---

## Data Sources

Transcriptomic data:

**NCBI Gene Expression Omnibus — GSE168404**

Regulatory analysis:

**TRRUST human transcriptional regulatory network**

Additional TF-target resources were used where applicable.

Raw and third-party database files are not redistributed in this repository. Source information and processing instructions are provided separately.

---

## Important Interpretation Note

This project represents a computational analysis and hypothesis-generation framework.

Enrichment and network associations should not be interpreted as experimental proof of regulatory causality. Candidate regulators require independent validation using additional datasets and/or experimental approaches.

---

## Interactive Dashboard

An interactive Streamlit dashboard is provided for exploring:

* differential-expression results
* GSEA pathways
* functional themes
* core genes
* regulator enrichment
* regulatory networks

**Live dashboard:** Coming soon

---

## Author

**Aakriti Sharma**

Molecular Medicine | Genomics | Bioinformatics | NGS Data Analysis

---

## License

This project is intended for research and educational use. Please see the repository license for details.

