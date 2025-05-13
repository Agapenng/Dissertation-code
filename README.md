# Alzheimer's Disease Bioinformatics Analysis

This R script performs a comprehensive bioinformatics analysis of two Alzheimer's disease gene expression datasets ([GSE5281](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5281) and [GSE48350](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48350)) from GEO. The workflow includes data preprocessing, differential expression analysis, visualization, machine learning classification, pathway enrichment, and network analysis.

## Features

- **Data Preprocessing:** Loads and normalizes GEO microarray datasets.
- **Differential Expression Analysis:** Identifies differentially expressed genes (DEGs) using the `limma` package.
- **Visualization:** Generates volcano plots, heatmaps, UMAP, PCA, and t-SNE plots for exploratory data analysis.
- **Machine Learning:** Trains and evaluates multiple classifiers (Logistic Regression, Random Forest, SVM, XGBoost, Lasso, Neural Network) for disease classification.
- **Model Evaluation:** Plots ROC curves, computes AUC, and summarizes confusion matrices for all models.
- **Feature Importance:** Assesses feature importance and correlation among selected genes.
- **Pathway and GO Enrichment:** Performs KEGG and Gene Ontology enrichment analysis on DEGs.
- **Network Analysis:** Visualizes protein-protein interaction networks using STRING-DB and `ggraph`.
- **Model Explainability:** Uses SHAP and DALEX for model interpretation.

## Requirements

- R (version 4.0 or higher recommended)
- R packages: `BiocManager`, `GEOquery`, `limma`, `clusterProfiler`, `org.Hs.eg.db`, `KEGGREST`, `preprocessCore`, `Rtsne`, `glmnet`, `pheatmap`, `dplyr`, `randomForest`, `pROC`, `e1071`, `tensorflow`, `keras`, `caret`, `nnet`, `VennDiagram`, `ggfortify`, `igraph`, `STRINGdb`, `httr`, `tidyverse`, `jsonlite`, `networkD3`, `ggraph`, `enrichplot`, `AnnotationDbi`, `hgu133plus2.db`, `ROCR`, `RColorBrewer`, `corrplot`, `shapper`, `DALEX`, `gridExtra`, `xgboost`

Install dependencies using:
```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "clusterProfiler", "org.Hs.eg.db", "KEGGREST", "igraph", "STRINGdb", "httr", "enrichplot", "AnnotationDbi", "hgu133plus2.db"))
install.packages(c("preprocessCore", "Rtsne", "glmnet", "pheatmap", "dplyr", "randomForest", "pROC", "e1071", "tensorflow", "keras", "caret", "nnet", "VennDiagram", "ggfortify", "tidyverse", "jsonlite", "networkD3", "ggraph", "ROCR", "RColorBrewer", "corrplot", "shapper", "DALEX", "gridExtra", "xgboost"))
