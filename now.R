##Install & Load Required Libraries

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "clusterProfiler", "org.Hs.eg.db", "KEGGREST"))
install.packages(c("preprocessCore", "Rtsne", "glmnet", "pheatmap", "dplyr", "randomForest", "pROC", 
                   "e1071", "tensorflow", "keras", "caret", "nnet", "VennDiagram"))


# Install and load the necessary package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("VennDiagram")
install.packages("ggfortify")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "igraph", "STRINGdb", "httr"))

install.packages(c("tidyverse", "jsonlite", "networkD3", "ggraph"))
install.packages("clusterProfiler")
install.packages("org.Hs.eg.db")
install.packages("enrichplot")

#Load Necessary libraries
library(GEOquery)
library(limma)
library(umap)
library(caret)
library(glmnet)
library(randomForest)
library(e1071)
library(xgboost)
library(ggplot2)
library(ggfortify)
library(reshape2)
library(pheatmap)
library(Rtsne)
library(dplyr)
library(VennDiagram)
library(AnnotationDbi)
library(hgu133plus2.db)
library(caret)
library(nnet)
library(pROC)
library(ROCR)
library(igraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(AnnotationDbi)
library(cnetplot)
library(emapplot)
library(STRINGdb)
library(httr)
library(jsonlite)
library(tidyverse)
library(networkD3)
library(ggraph)
library(httr)
library(jsonlite)
library(RColorBrewer)
library(corrplot)

#Preprocessing GEO dataset
gset_5281 <- getGEO("GSE5281", GSEMatrix = TRUE, AnnotGPL = TRUE)

# load series and platform data from GEO

if (length(gset_5281) > 1) idx_5281 <- grep("GPL570", attr(gset_5281, "names")) else idx_5281 <- 1
gset_5281 <- gset_5281[[idx_5281]]

# make proper column names to match toptable
fvarLabels(gset_5281) <- make.names(fvarLabels(gset_5281))

# group membership for all samples
gsms_5281 <- paste0("00000000000000000000000000000000000000000000000000",
                    "00000000000000000000000011111111111111111111111111",
                    "11111111111111111111111111111111111111111111111111",
                    "11111111111")
sml_5281 <- strsplit(gsms_5281, split="")[[1]]

# log2 transformation
ex_5281 <- exprs(gset_5281)
qx_5281<- as.numeric(quantile(ex_5281, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC_5281 <- (qx_5281[5] > 100) ||
  (qx_5281[6]-qx_5281[1] > 50 && qx_5281[2] > 0)
if (LogC_5281) { ex_5281[which(ex_5281<= 0)] <- NaN
exprs(gset_5281) <- log2(ex_5281) }


# assign samples to groups and set up design matrix
gs_5281 <- factor(sml_5281)  # Only keep selected samples
groups_5281 <- make.names(c("Normal","Alzheimers Disease"))
levels(gs_5281) <- groups_5281
gset_5281$group <- gs_5281
table(gset_5281$group)
design_5281 <- model.matrix(~0 + gs_5281)  
colnames(design_5281) <- levels(gs_5281)

#gset_5281<- gset_5281[complete.cases(exprs(gset_5281)), ] # skip missing values
# Recreate design matrix
fit_5281 <- lmFit(gset_5281, design_5281)  # fit linear model



# set up contrasts of interest and recalculate model coefficients
cts_5281 <- paste(groups_5281[1], groups_5281[2], sep="-")
cont.matrix_5281 <- makeContrasts(contrasts=cts_5281, levels=design_5281)
fit2_5281 <- contrasts.fit(fit_5281, cont.matrix_5281)

# compute statistics and table of top significant genes
fit2_5281 <- eBayes(fit2_5281)
tT_5281 <- topTable(fit2_5281, adjust="fdr", sort.by="B", number=Inf)


filtered_data_5281 <- tT_5281[tT_5281$adj.P.Val < 0.05 & (tT_5281$logFC)> 1,]
tT_5281 <- subset(tT_5281, select=c("ID","adj.P.Val","P.Value","logFC","GI","Gene.symbol","Gene.title", "Chromosome.location","Nucleotide.Title" ))


# Assume 'tT_5281' is your data frame
tT_5281 <- as.data.frame(tT_5281)

tT_5281$Significant <- ifelse(
  tT_5281$adj.P.Val < 0.05 & tT_5281$logFC > 1, "Up",
  ifelse(tT_5281$adj.P.Val < 0.05 & tT_5281$logFC < -1, "Down", "Not")
)


# Create the Significant column with specified categories
tT_5281 <- tT_5281 %>%
  mutate(Significant = case_when(
    adj.P.Val < 0.05 & logFC> 1  ~ "Up",
    adj.P.Val < 0.05 & logFC < -1 ~ "Down",
    TRUE ~ "Not"
  ))

# Create the volcano plot with color
ggplot(tT_5281, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
  geom_point() +
  scale_color_manual(values = c("Down" = "blue", "Not" = "grey", "Up" = "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Genes GSE5281 for Normal and Alzheimer's Disease ", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value")

#Pheatmap for Gse5281
# Categorizing all genes in the toptable (tT)
tT_heat_5281 <- tT_5281 %>%
  mutate(Significant = case_when(
    adj.P.Val < 0.05 & logFC > 1  ~ "Up",
    adj.P.Val < 0.05 & logFC < -1 ~ "Down",
    TRUE ~ "Not"
  ))
# Filter only significant genes
significant_genes_5281 <- tT_heat_5281 %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

# Select the top 10 significant genes based on adjusted p-value
top_10_genes_5281 <- significant_genes_5281 %>%
  arrange(adj.P.Val) %>%
  head(10)

# Extract expression data for top 10 genes
top_10_expression_5281 <- exprs(gset_5281)[rownames(gset_5281) %in% rownames(top_10_genes_5281), ]
top_gene_names_5281 <- rownames(top_10_genes_5281)

print(significant_genes_5281)
# Set the seed for reproducibility
set.seed(123)

# Specify the number of samples you want to keep
sample_size <- 20

# Randomly select sample indices
selected_samples5281<- sample(ncol(top_10_expression_5281), sample_size)

# Subset the expression matrix to include only the selected samples
reduced_top_10_expression5281 <- top_10_expression_5281[, selected_samples5281]


# Create the heatmap with adjusted font sizes and rotation for better readability
pheatmap(reduced_top_10_expression5281,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8, # Adjust row font size
         fontsize_col = 8, # Adjust column font size
         angle_col = 90)  # Rotate column names for better readability


# Map Probe IDs to Gene Symbols
probe_ids_5281 <- rownames(top_10_genes_5281)
gene_symbols_5281 <- mapIds(hgu133plus2.db,
                            keys = probe_ids_5281,
                            column = "SYMBOL",
                            keytype = "PROBEID",
                            multiVals = "first")

# Update the row names with gene symbols
rownames(reduced_top_10_expression5281) <- gene_symbols_5281[rownames(reduced_top_10_expression5281)]

# Create the heatmap with adjusted font sizes and rotation for better readability

pheatmap(reduced_top_10_expression5281,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 10,  # Adjust row font size
         fontsize_col = 10,  # Adjust column font size
         angle_col = 90,
         main = "Heatmap of Top 10 Differentially Expressed Genes")     # Rotate column names for better readability


# Convert expression data to a matrix (genes as rows, samples as columns)
expression_data_5281 <- t(exprs(gset_5281))  # Transpose for UMAP (samples as rows)

# Perform UMAP dimensionality reduction
set.seed(123)  # For reproducibility
umap_result_5281 <- umap(expression_data_5281, n_neighbors = 15, min_dist = 0.1, metric = "euclidean")

# Create a dataframe for plotting
umap_df_5281 <- data.frame(
  UMAP1 = umap_result_5281$layout[,1],  # UMAP x-axis
  UMAP2 = umap_result_5281$layout[,2],  # UMAP y-axis
  Group = gset_5281$group  # Group labels (ND vs AD)
)

# Generate UMAP scatter plot
ggplot(umap_df_5281, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("Normal" = "green", "Alzheimers.Disease" = "red")) +
  theme_minimal() +
  labs(title = "GSE5281: UMAP(nbrs=15)", x = "UMAP1", y = "UMAP2") +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )



#Preprocessing GEO dataset 48350
gset <- getGEO("GSE48350", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
# Define sample labels correctly
gsms <- paste0("00000000000000000000000000000000000000000000000000",
               "00000000000000000000000000000000000000000000000000",
               "00000000000000000000000000000000000000000000000000",
               "00000000000000000000000111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "111")
sml <- strsplit(gsms, split = "")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("Normal","Alzheimers Disease"))
levels(gs) <- groups
gset$group <- gs
table(gset$group)
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model


#set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
colnames(tT)  # Print available columns


# Assume 'tT_48350' is your data frame
tT <- as.data.frame(tT)

tT$Significant <- ifelse(
  tT$adj.P.Val < 0.05 & tT$logFC > 1, "Up",
  ifelse(tT$adj.P.Val < 0.05 & tT$logFC < -1, "Down", "Not")
)


# Create the Significant column with specified categories
tT <- tT %>%
  mutate(Significant = case_when(
    adj.P.Val < 0.05 & logFC> 1  ~ "Up",
    adj.P.Val < 0.05 & logFC < -1 ~ "Down",
    TRUE ~ "Not"
  ))

# Create the volcano plot with color
ggplot(tT, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
  geom_point() +
  scale_color_manual(values = c("Down" = "blue", "Not" = "grey", "Up" = "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Genes GSE48350 for Normal and Alzheimer's Disease ", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value")


# Categorizing all genes in the toptable (tT)
tT_heat_48350 <- tT %>%
  mutate(Significant = case_when(
    adj.P.Val < 0.05 & logFC > 1  ~ "Up",
    adj.P.Val < 0.05 & logFC < -1 ~ "Down",
    TRUE ~ "Not"
  ))
# Filter only significant genes
significant_genes_48350 <- tT_heat_48350 %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

# Select the top 10 significant genes based on adjusted p-value
top_10_genes_48350 <- significant_genes_48350 %>%
  arrange(adj.P.Val) %>%
  head(10)

# Extract expression data for top 10 genes
top_10_expression_48350 <- exprs(gset)[rownames(gset) %in% rownames(top_10_genes_48350), ]
top_gene_names_48350 <- rownames(top_10_genes_48350)
# Set the seed for reproducibility
set.seed(123)

# Specify the number of samples you want to keep
sample_size <- 20

# Randomly select sample indices
selected_samples48350<- sample(ncol(top_10_expression_48350), sample_size)

# Subset the expression matrix to include only the selected samples
reduced_top_10_expression48350 <- top_10_expression_48350[, selected_samples48350]


# Create the heatmap with adjusted font sizes and rotation for better readability
pheatmap(reduced_top_10_expression48350,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8, # Adjust row font size
         fontsize_col = 8, # Adjust column font size
         angle_col = 90)  # Rotate column names for better readability


# Map Probe IDs to Gene Symbols
probe_ids_48350 <- rownames(top_10_genes_48350)
gene_symbols_48350 <- mapIds(hgu133plus2.db,
                            keys = probe_ids_48350,
                            column = "SYMBOL",
                            keytype = "PROBEID",
                            multiVals = "first")

# Update the row names with gene symbols
rownames(reduced_top_10_expression48350) <- gene_symbols_48350[rownames(reduced_top_10_expression48350)]

# Create the heatmap with adjusted font sizes and rotation for better readability

pheatmap(reduced_top_10_expression48350,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 10,  # Adjust row font size
         fontsize_col = 10,  # Adjust column font size
         angle_col = 90,
         main = "Heatmap of Top 10 Differentially Expressed Genes")     # Rotate column names for better readability


# Convert expression data to a matrix (genes as rows, samples as columns)
expression_data_48350 <- t(exprs(gset))  # Transpose for UMAP (samples as rows)

# Perform UMAP dimensionality reduction
set.seed(123)  # For reproducibility
umap_result_48350 <- umap(expression_data_48350, n_neighbors = 15, min_dist = 0.1, metric = "euclidean")

# Create a dataframe for plotting
umap_df_48350 <- data.frame(
  UMAP1 = umap_result_48350$layout[,1],  # UMAP x-axis
  UMAP2 = umap_result_48350$layout[,2],  # UMAP y-axis
  Group = gset$group  # Group labels (Normal vs AD)
)

# Generate UMAP scatter plot
ggplot(umap_df_48350, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("Normal" = "green", "Alzheimers.Disease" = "red")) +
  theme_minimal() +
  labs(title = "GSE48350: UMAP(nbrs=15)", x = "UMAP1", y = "UMAP2") +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


# Merge DEGs from GSE5281 and GSE48350
# Merge DEGs from GSE5281 and GSE48350
DEGs_5281 <- rownames(significant_genes_5281)
DEGs_48350 <- rownames(significant_genes_48350)

# Union of DEGs without duplicates
common_DEGs <- unique(c(DEGs_5281, DEGs_48350))

# Venn Diagram
venn.plot <- draw.pairwise.venn(
  area1 = length(DEGs_5281),
  area2 = length(DEGs_48350),
  cross.area = length(intersect(DEGs_5281, DEGs_48350)),
  category = c("GSE5281", "GSE48350"),
  fill = c("red", "blue"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5
)
grid.draw(venn.plot)

print(common_DEGs)
# Find common DEGs
common_DEGs <- intersect(DEGs_5281, DEGs_48350)
print(common_DEGs)

# Define the common differentially expressed genes (DEGs) as probe IDs
common_DEGs <- c("223913_s_at", "212833_at", "206552_s_at")

# Retrieve gene symbols from the annotation database
gene_symbols <- mapIds(hgu133plus2.db,
                       keys = common_DEGs,
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")

# Print the mapped gene symbols
print(gene_symbols)


# PCA Analysis
gset_merged <- cbind(exprs(gset_5281)[common_DEGs, ], exprs(gset)[common_DEGs, ])
gset_merged <- gset_merged[!duplicated(rownames(gset_merged)), ]  # Remove duplicate genes
gset_merged <- gset_merged[, !duplicated(t(gset_merged))]  # Remove duplicate samples
group_labels <- c(gset_5281$group, gset$group)

# PCA Analysis
pca_merged <- prcomp(t(gset_merged), scale. = TRUE)

# Plot PCA
autoplot(pca_merged, data = data.frame(group = group_labels), colour = 'group') +
  ggtitle("PCA of Merged GSE5281 and GSE48350")

# t-SNE Analysis
set.seed(123)
tsne_merged <- Rtsne(t(gset_merged), perplexity = 30, theta = 0.0)

tsne_df_merged <- data.frame(
  X = tsne_merged$Y[,1],
  Y = tsne_merged$Y[,2],
  Group = group_labels
)

# Plot t-SNE
ggplot(tsne_df_merged, aes(x = X, y = Y, color = Group)) +
  geom_point(size = 3) +
  ggtitle("t-SNE of Merged GSE5281 and GSE48350")


# Recursive Feature Elimination (RFE)
# Recursive Feature Elimination (RFE)
set.seed(123)
control <- rfeControl(functions = rfFuncs, method = "cv", number = 10)

gset_merged_df <- as.data.frame(t(gset_merged))

# Run RFE
rfe_results_merged <- rfe(gset_merged_df, group_labels, sizes = c(1:10), rfeControl = control)

# Plot RFE performance
plot(rfe_results_merged, type = c("g", "o"), main = "RFE Performance - Merged GSE5281 & GSE48350")

#MACHINE LEARNING
# Prepare Data for ML
merged_expr <- cbind(exprs(gset_5281)[common_DEGs, ], exprs(gset)[common_DEGs, ])
group_labels <- c(gset_5281$group, gset$group)
trainIndex <- createDataPartition(group_labels, p = 0.7, list = FALSE)
trainData <- t(merged_expr[, trainIndex])
testData <- t(merged_expr[, -trainIndex])
trainLabels <- group_labels[trainIndex]
testLabels <- group_labels[-trainIndex]

# Feature Selection (RFE)
control <- rfeControl(functions = rfFuncs, method = "cv", number = 10)
rfe_results <- rfe(trainData, trainLabels, sizes = c(5, 10, 20), rfeControl = control)
selected_features <- predictors(rfe_results)
trainData <- trainData[, selected_features]
testData <- testData[, selected_features]


# Ensure labels are factors first
trainLabels <- as.factor(trainLabels)
testLabels <- as.factor(testLabels)

# Convert factor levels to numeric (0 and 1)
trainLabels <- as.numeric(trainLabels) - 1
testLabels <- as.numeric(testLabels) - 1

# Now train the logistic regression model


# Train Models
## Train Logistic Regression Model

log_model <- glm(trainLabels ~ ., data = as.data.frame(trainData), family = binomial)

testData <- as.data.frame(testData)
dim(trainData)  # Check the dimensions of trainData
dim(testData)   # Check the dimensions of testData
length(trainLabels)  # Ensure the number of labels matches trainData rows
length(testLabels)   # Ensure the number of labels matches testData rows

testLabels <- testLabels[1:nrow(testData)]  # Ensure testLabels matches testData rows
length(testLabels) == nrow(testData)  # Should return TRUE

unique(testLabels)  # Should return only 0 and 1
valid_indices <- !is.na(testLabels)  # Identify non-NA rows
testLabels <- testLabels[valid_indices]
testData <- testData[valid_indices, ]  # Remove NA rows from testData



log_preds <- predict(log_model, as.data.frame(testData), type = "response")
length(log_preds)    # Ensure predictions match testLabels
sum(is.na(log_preds))  # Should return 0 if no NA values
# Step 3: Check Length Consistency
if (length(log_preds) != length(testLabels)) {
  testLabels <- testLabels[1:length(log_preds)]
}

# Step 4: Compute ROC
roc_log <- roc(testLabels, log_preds)
auc_log <- auc(roc_log)
print(paste("AUC Logistic Regression:", auc_log))

 plot(roc_log, col = "blue", main = "ROC Curve - Logistic Regression")

#RANDOM FOREST
trainLabels <- as.factor(trainLabels)
rf_model <- randomForest(x = trainData, y = trainLabels)
rf_preds <- predict(rf_model, testData, type = "prob")[, 2]
roc_rf <- roc(testLabels, rf_preds)
auc_rf <- auc(roc_rf)
plot(roc_rf, col = "red", main = "ROC Curve - Random Forest")
print(auc_rf)

#SVM
svm_model <- svm(trainLabels ~ ., data = trainData, probability = TRUE)
svm_preds <- attr(predict(svm_model, testData, probability = TRUE), "probabilities")[, 2]
roc_svm <- roc(testLabels, svm_preds)
auc_svm <- auc(roc_svm)
plot(roc_svm, col = "blue", main = "ROC Curve - SVM")
print(auc_svm)

#svm_model <- svm(trainLabels ~ ., data = as.data.frame(trainData), probability = TRUE)
xgb_train <- xgb.DMatrix(data = as.matrix(trainData), label = as.numeric(trainLabels)-1)
xgb_model <- xgboost(data = xgb_train, nrounds = 100, objective = "binary:logistic", eval_metric = "auc")
xgb_preds <- predict(xgb_model, xgb.DMatrix(data = as.matrix(testData)))
auc_xgb <- auc(testLabels, xgb_preds)


## Train Lasso Regression
lasso_model <- cv.glmnet(as.matrix(trainData), as.numeric(trainLabels) - 1, family = "binomial", alpha = 1)
pred_lasso <- predict(lasso_model, newx = as.matrix(testData), s = "lambda.min", type = "response")
# Ensure pred_lasso is a numeric vector
pred_lasso <- as.numeric(pred_lasso)  # Convert matrix to numeric vector

# Compute ROC and AUC
roc_lasso <- roc(testLabels, pred_lasso)
auc_lasso <- auc(roc_lasso)

# Print AUC value
print(paste("Lasso Regression AUC:", auc_lasso))


auc_lasso <- auc(roc_lasso)
plot(roc_lasso, col = "orange", main = "ROC Curve - Lasso")

## Train Neural Network Model
library(nnet)
library(caret)

# Reduce number of features using Recursive Feature Elimination (RFE)
set.seed(123)
control <- rfeControl(functions = rfFuncs, method = "cv", number = 5)
rfe_results <- rfe(trainData, trainLabels, sizes = c(10, 20, 50), rfeControl = control)

# Select top features
selected_features <- predictors(rfe_results)
trainData_reduced <- trainData[, selected_features]
testData_reduced <- testData[, selected_features]

# Scale the data for better convergence
trainData_scaled <- as.data.frame(scale(trainData_reduced))
testData_scaled <- as.data.frame(scale(testData_reduced))

# Convert labels to factors for classification
trainLabels <- as.factor(trainLabels)
testLabels <- as.factor(testLabels)

# Train Neural Network with fewer neurons
nn_model <- nnet(trainLabels ~ ., data = trainData_scaled, size = 3, decay = 0.1, maxit = 200)

# Make Predictions
nn_preds <- predict(nn_model, testData_scaled, type = "raw")

# Convert predictions to class labels
nn_class_preds <- ifelse(nn_preds > 0.5, 1, 0)

# Compute ROC and AUC
roc_nn <- roc(testLabels, as.numeric(nn_preds))
auc_nn <- auc(roc_nn)

# Print AUC
print(paste("Neural Network AUC:", auc_nn))
plot(roc_nn, col = "green", main = "ROC Curve - Neural Network")


# Load required libraries
library(pROC)
library(ggplot2)

# Define a common sequence of FPR points (0 to 1 with small intervals)
common_fpr <- seq(0, 1, length.out = 100)

# Function to interpolate sensitivity at common FPRs
extract_tpr <- function(roc_obj, fpr_seq) {
  coords(roc_obj, x = fpr_seq, input = "specificity", ret = "sensitivity", transpose = TRUE, as.list = FALSE)
}

# Compute ROC for each model
roc_log <- roc(testLabels, log_preds)
roc_rf <- roc(testLabels, rf_preds)
roc_svm <- roc(testLabels, svm_preds)
roc_xgb <- roc(testLabels, xgb_preds)
roc_lasso <- roc(testLabels, pred_lasso)
roc_nn <- roc(testLabels, nn_preds)

# Compute AUC for each model
auc_log <- auc(roc_log)
auc_rf <- auc(roc_rf)
auc_svm <- auc(roc_svm)
auc_xgb <- auc(roc_xgb)
auc_lasso <- auc(roc_lasso)
auc_nn <- auc(roc_nn)

# Print AUC Scores
print(paste("AUC Logistic Regression:", auc_log))
print(paste("AUC Random Forest:", auc_rf))
print(paste("AUC SVM:", auc_svm))
print(paste("AUC XGBoost:", auc_xgb))
print(paste("AUC Lasso:", auc_lasso))
print(paste("AUC Neural Network:", auc_nn))

# Extract TPR values at common FPRs
tpr_log <- extract_tpr(roc_log, common_fpr)
tpr_rf <- extract_tpr(roc_rf, common_fpr)
tpr_svm <- extract_tpr(roc_svm, common_fpr)
tpr_xgb <- extract_tpr(roc_xgb, common_fpr)
tpr_lasso <- extract_tpr(roc_lasso, common_fpr)
tpr_nn <- extract_tpr(roc_nn, common_fpr)

# Create a data frame with aligned values
roc_data <- data.frame(
  FPR = rep(common_fpr, 6),
  TPR = c(tpr_log, tpr_rf, tpr_svm, tpr_xgb, tpr_lasso, tpr_nn),
  Model = rep(c("Logistic Regression", "Random Forest", "SVM", "XGBoost", "Lasso", "Neural Network"), each = length(common_fpr))
)

# Plot ROC Curves
ggplot(roc_data, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(size = 1.2) +
  geom_abline(linetype = "dashed", color = "gray") +  # Reference diagonal
  theme_minimal() +
  labs(
    title = "ROC Curves for Multiple Models",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  )


## Load Required Libraries
library(caret)
library(pROC)
library(dplyr)
library(igraph)
library(ggplot2)
library(reshape2)
library(ROCR)
log_preds <- log_preds[1:length(testLabels)]
rf_preds <- rf_preds[1:length(testLabels)]
svm_preds <- svm_preds[1:length(testLabels)]
xgb_preds <- xgb_preds[1:length(testLabels)]
pred_lasso <- pred_lasso[1:length(testLabels)]
nn_preds <- nn_preds[1:length(testLabels)]




## Collecting Confusion Matrices from Models
conf_matrices <- list(
  "Logistic Regression" = confusionMatrix(as.factor(ifelse(log_preds > 0.5, 1, 0)), as.factor(testLabels)),
  "Random Forest" = confusionMatrix(as.factor(ifelse(rf_preds > 0.5, 1, 0)), as.factor(testLabels)),
  "SVM" = confusionMatrix(as.factor(ifelse(svm_preds > 0.5, 1, 0)), as.factor(testLabels)),
  "XGBoost" = confusionMatrix(as.factor(ifelse(xgb_preds > 0.5, 1, 0)), as.factor(testLabels)),
  "Lasso Regression" = confusionMatrix(as.factor(ifelse(pred_lasso > 0.5, 1, 0)), as.factor(testLabels)),
  "Neural Network" = confusionMatrix(as.factor(ifelse(nn_preds > 0.5, 1, 0)), as.factor(testLabels))
)

## Extracting Confusion Matrices into a Table
conf_matrix_df <- do.call(rbind, lapply(names(conf_matrices), function(model) {
  cm <- conf_matrices[[model]]
  data.frame(
    Model = model,
    Accuracy = cm$overall["Accuracy"],
    Sensitivity = cm$byClass["Sensitivity"],
    Specificity = cm$byClass["Specificity"],
    Precision = cm$byClass["Precision"],
    F1_Score = cm$byClass["F1"]
  )
}))

print(conf_matrix_df)


# Train the Random Forest Model (if not already trained)
rf_model <- randomForest(x = trainData, y = trainLabels, importance = TRUE)

# Extract Feature Importance
importance_df <- as.data.frame(importance(rf_model))

# Add Probe IDs (Feature Names)
importance_df$ProbeID <- rownames(importance_df)

# Map Probe IDs to Gene Symbols
gene_symbols <- mapIds(hgu133plus2.db,
                       keys = importance_df$ProbeID,
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")

# Replace missing values with "Unknown"
gene_symbols[is.na(gene_symbols)] <- "Unknown"

# Create a new column with Probe IDs
importance_df$Feature <- paste(importance_df$ProbeID, "(", gene_symbols, ")", sep="")

# Sort by importance (MeanDecreaseGini)
importance_df <- importance_df[order(-importance_df$MeanDecreaseGini), ]

# Plot Feature Importance with Probe IDs
ggplot(importance_df, aes(x = reorder(Feature, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Feature Importance - Random Forest",
       x = "Probe ID (Gene Symbol)",
       y = "Mean Decrease in Gini")


print(importance_df)

# Predictions
log_preds <- predict(log_model, as.data.frame(testData), type = "response")
rf_preds <- predict(rf_model, testData, type = "prob")[,2]
svm_preds <- attr(predict(svm_model, as.data.frame(testData), probability = TRUE), "probabilities")[,2]
xgb_preds <- predict(xgb_model, xgb.DMatrix(data = as.matrix(testData)))
lasso_preds <- predict(lasso_model, newx = as.matrix(testData), s = "lambda.min", type = "response")
nn_preds <- predict(nn_model, scale(as.data.frame(testData)), type = "raw")

# Compute Mean Gini Importance
mean_gini_importance <- mean(importance(rf_model)[,1])
print(paste("Mean Gini Importance:", mean_gini_importance))



# Compute Spearman Correlation Matrix
correlation_matrix <- cor(trainData, method = "spearman")

# Extract Probe IDs from row names
probe_ids <- rownames(correlation_matrix)

# Map Probe IDs to Gene Symbols
gene_symbols <- mapIds(hgu133plus2.db,
                       keys = probe_ids,
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")

# Replace missing gene symbols with "Unknown"
gene_symbols[is.na(gene_symbols)] <- "Unknown"

# Rename row and column names with Probe IDs and Gene Symbols
rownames(correlation_matrix) <- paste(probe_ids, "(", gene_symbols, ")", sep="")
colnames(correlation_matrix) <- paste(probe_ids, "(", gene_symbols, ")", sep="")

# Plot Correlation Matrix using Corrplot (Color Method)
corrplot(correlation_matrix, method = "color", type = "lower",
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         col = colorRampPalette(c("blue", "white", "red"))(200))

# Plot Correlation Matrix using Corrplot (Circle Method)
corrplot(correlation_matrix, method = "circle")

# Plot Correlation Matrix using Pheatmap
pheatmap(correlation_matrix, method = "color",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         show_rownames = TRUE, show_colnames = TRUE,
         fontsize_row = 6, fontsize_col = 6,
         main = "Correlation Matrix of Features (Probe IDs)")




# Predictions
log_preds <- predict(log_model, as.data.frame(testData), type = "response")
rf_preds <- predict(rf_model, testData, type = "prob")[,2]
svm_preds <- attr(predict(svm_model, as.data.frame(testData), probability = TRUE), "probabilities")[,2]
xgb_preds <- predict(xgb_model, xgb.DMatrix(data = as.matrix(testData)))
lasso_preds <- predict(lasso_model, newx = as.matrix(testData), s = "lambda.min", type = "response")
nn_preds <- predict(nn_model, scale(as.data.frame(testData)), type = "raw")

# Compute Mean Gini Importance
mean_gini_importance <- mean(importance(rf_model)[,1])
print(paste("Mean Gini Importance:", mean_gini_importance))



## Prediction Distribution for each model
par(mfrow=c(1,3))
hist(log_preds, main="Logistic regression Prediction Probabilities", col="black", xlab="Predicted Probabilities")
hist(rf_preds, main="Random Forest Prediction Probabilities", col="blue", xlab="Predicted Probabilities")
hist(xgb_preds, main="XGBoost Prediction Probabilities", col="red", xlab="Predicted Probabilities")
hist(svm_preds, main="SVM Prediction Probabilities", col="green", xlab="Predicted Probabilities")
hist(pred_lasso, main="Lasso Prediction Probabilities", col="orange", xlab="Predicted Probabilities")
hist(nn_preds, main="Neural Network Prediction Probabilities", col="purple", xlab="Predicted Probabilities")
par(mfrow=c(1,3))

install.packages("shapper")
install.packages("DALEX")
install.packages("ggplot2")
library(shapper)
library(DALEX)
library(ggplot2)
library(gridExtra)

explainer_log <- DALEX::explain(
  model = log_model, 
  data = trainData, 
  y = trainLabels, 
  label = "Logistic Regression"
)


explainer_rf <- DALEX::explain(
  model = rf_model, 
  data = trainData, 
  y = trainLabels, 
  label = "Random Forest Model"
)

explainer_svm <- DALEX::explain(
  model = svm_model, 
  data = trainData, 
  y = trainLabels, 
  label = "SVM Model"
)

explainer_xgb <- DALEX::explain(
  model = xgb_model, 
  data = trainData, 
  y = trainLabels, 
  label = "XGBoost Model"
)

explainer_lasso <- DALEX::explain(
  model = lasso_model, 
  data = trainData, 
  y = trainLabels, 
  label = "Lasso Regression"
)

explainer_nn <- DALEX::explain(
  model = nn_model, 
  data = trainData, 
  y = trainLabels, 
  label = "Neural Network"
)

shap_log <- shap(explainer_log, new_observation = testData[1:10, ])
shap_rf <- shap(explainer_rf, new_observation = testData[1:10, ])
shap_svm <- shap(explainer_svm, new_observation = testData[1:10, ])
shap_xgb <- shap(explainer_xgb, new_observation = testData[1:10, ])
shap_lasso <- shap(explainer_lasso, new_observation = testData[1:10, ])
shap_nn <- shap(explainer_nn, new_observation = testData[1:10, ])


plot_rf <- plot(shap_rf)
plot_log <- plot(shap_log)
plot_svm <- plot(shap_svm)
plot_xgb <- plot(shap_xgb)
plot_lasso <- plot(shap_lasso)
plot_nn <- plot(shap_nn)

# Display all models together
grid.arrange(plot_rf, plot_log, plot_svm, plot_xgb, plot_lasso, plot_nn, ncol = 2)


shap_long_rf <- shap_rf$shap_score_long
shap_long_xgb <- shap_xgb$shap_score_long

ggplot(shap_long_rf, aes(value, contribution)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~ variable, scales = "free") +
  theme_minimal() +
  labs(title = "SHAP Dependence Plot - Random Forest", x = "Feature Value", y = "SHAP Contribution")

ggplot(shap_long_xgb, aes(value, contribution)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~ variable, scales = "free") +
  theme_minimal() +
  labs(title = "SHAP Dependence Plot - XGBoost", x = "Feature Value", y = "SHAP Contribution")

plot(shap_rf, row_id = 1)  # Random Forest Prediction Explanation
plot(shap_xgb, row_id = 1)  # XGBoost Prediction Explanation

grid.arrange(
  plot(shap_rf), 
  plot(shap_xgb), 
  plot(shap_svm), 
  plot(shap_log), 
  plot(shap_lasso), 
  plot(shap_nn), 
  ncol = 2
)




# Get gene symbols for DEGs (ensure that the IDs are converted to Entrez IDs)
gene_symbols_5281 <- significant_genes_5281$Gene.symbol
gene_symbols_48350 <- significant_genes_48350$Gene.symbol

# Convert gene symbols to Entrez IDs using org.Hs.eg.db
gene_ids_5281 <- bitr(gene_symbols_5281, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_ids_48350 <- bitr(gene_symbols_48350, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Merge the results
all_gene_ids <- unique(c(gene_ids_5281$ENTREZID, gene_ids_48350$ENTREZID))

# KEGG Enrichment Analysis
kegg_results <- enrichKEGG(gene = all_gene_ids,
                           organism = "hsa", # Human species
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.2)
# Barplot for KEGG enrichment results
barplot(kegg_results, showCategory = 10, main = "KEGG Pathway Enrichment")

# Dotplot for KEGG pathways
dotplot(kegg_results, showCategory = 10, title = "KEGG Pathway Enrichment")

# KEGG Network Plot (Connections between genes and pathways)
cnetplot(kegg_results, showCategory = 5, title = "KEGG Network Plot")

# Convert Entrez IDs to Gene Symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = all_gene_ids,  # Your list of Entrez IDs
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")


# Ensure no missing values
gene_symbols[is.na(gene_symbols)] <- "Unknown"

# Rename the names in KEGG results
kegg_results@result$geneID <- sapply(strsplit(kegg_results@result$geneID, "/"), 
                                     function(ids) paste(gene_symbols[ids], collapse="/"))

# Generate cnetplot with gene symbols instead of Entrez IDs
cnetplot(kegg_results, 
         showCategory = 5,   # Show top 5 enriched pathways
         node_label = "all",  # Show gene symbols on nodes
         foldChange = NULL,   # Optional: Color based on logFC values
         title = "KEGG Pathway Network with Gene Symbols")



# GO Enrichment Analysis (for Biological Process, Molecular Function, and Cellular Component)
go_results_BP <- enrichGO(gene = all_gene_ids,
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID",
                          ont = "BP", # Biological Process
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)



# Visualize the GO enrichment results
barplot(go_results_BP, showCategory = 10)
dotplot(go_results_BP, showCategory = 10)


# Extract gene symbols from DEGs
gene_symbols_5281 <- significant_genes_5281$Gene.symbol
gene_symbols_48350 <- significant_genes_48350$Gene.symbol

# Convert gene symbols to Entrez IDs
gene_ids_5281 <- bitr(gene_symbols_5281, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_ids_48350 <- bitr(gene_symbols_48350, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Merge gene lists and remove duplicates
all_gene_ids <- unique(c(gene_ids_5281$ENTREZID, gene_ids_48350$ENTREZID))

# Biological Process (BP)
go_results_BP <- enrichGO(gene = all_gene_ids,
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID",
                          ont = "BP", # Biological Process
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)

# Molecular Function (MF)
go_results_MF <- enrichGO(gene = all_gene_ids,
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID",
                          ont = "MF", # Molecular Function
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)

# Cellular Component (CC)
go_results_CC <- enrichGO(gene = all_gene_ids,
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID",
                          ont = "CC", # Cellular Component
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)


# Barplots for GO categories
barplot(go_results_BP, showCategory = 10, main = "GO Biological Process")
barplot(go_results_MF, showCategory = 10, main = "GO Molecular Function")
barplot(go_results_CC, showCategory = 10, main = "GO Cellular Component")

# Dotplots for better visualization
dotplot(go_results_BP, showCategory = 10, title = "GO Biological Process")
dotplot(go_results_MF, showCategory = 10, title = "GO Molecular Function")
dotplot(go_results_CC, showCategory = 10, title = "GO Cellular Component")

# Remove NAs (replace missing symbols with "Unknown")
gene_symbols[is.na(gene_symbols)] <- "Unknown"

# Replace Entrez IDs in GO enrichment results with Gene Symbols
go_results_BP@result$geneID <- sapply(strsplit(go_results_BP@result$geneID, "/"), 
                                      function(ids) paste(gene_symbols[ids], collapse="/"))

# Network plot of enriched GO terms
cnetplot(go_results_BP, showCategory = 5, title = "GO Biological Process Network")

go_results_MF@result$geneID <- sapply(strsplit(go_results_MF@result$geneID, "/"), 
                                      function(ids) paste(gene_symbols[ids], collapse="/"))


cnetplot(go_results_MF, showCategory = 5, title = "GO Molecular Function Network")

go_results_CC@result$geneID <- sapply(strsplit(go_results_CC@result$geneID, "/"), 
                                      function(ids) paste(gene_symbols[ids], collapse="/"))

cnetplot(go_results_CC, showCategory = 5, title = "GO Cellular Component Network")

# Pathway-Gene Graph using enrichplot

# Compute pairwise similarity between enriched GO terms
go_results_BP <- pairwise_termsim(go_results_BP)

# Generate the pathway graph using `emapplot()`
emapplot(go_results_BP, showCategory = 10) + ggtitle("GO BP Pathway Graph")

go_results_MF <- pairwise_termsim(go_results_MF)

# Generate the pathway graph using `emapplot()`
emapplot(go_results_MF, showCategory = 10) + ggtitle("GO MF Pathway Graph")

go_results_CC <- pairwise_termsim(go_results_CC)

# Generate the pathway graph using `emapplot()`
emapplot(go_results_CC, showCategory = 10) + ggtitle("GO CC Pathway Graph")



# Define target genes
target_genes <- c("SLC25A46", "TAC1", "MIR7-3HG")

# Function to fetch STRING-DB interactions
fetch_string_interactions <- function(gene) {
  url <- paste0("https://string-db.org/api/json/interaction_partners?identifiers=", gene, "&limit=10")
  
  response <- GET(url)
  
  if (http_status(response)$category == "Success") {
    response_content <- content(response, as = "text", encoding = "UTF-8")
    
    # Print raw JSON response for debugging
    print(paste("Raw API Response for:", gene))
    print(response_content)
    
    # Convert response to JSON
    content_data <- fromJSON(response_content, flatten = TRUE)
    
    # Check if response contains valid data
    if (!is.null(content_data) && length(content_data) > 0) {
      if ("preferredName_B" %in% names(content_data)) {
        interaction_df <- data.frame(
          Gene = gene,
          Interaction_Partner = content_data$preferredName_B,
          Confidence_Score = content_data$score
        )
        return(interaction_df)
      } else {
        print(paste("No valid interactions found for", gene))
      }
    }
  } else {
    print(paste("API request failed for", gene))
  }
  return(NULL)
}

# Fetch STRING-DB interactions
all_string_interactions <- lapply(target_genes, fetch_string_interactions)

# Combine results into a dataframe
string_interaction_df <- do.call(rbind, all_string_interactions)

# Print extracted STRING-DB interactions
print("STRING-DB Drug-Protein Interactions:")
print(string_interaction_df)



# Define nodes (genes + their interactors)
nodes <- unique(c(string_interaction_df$Gene, string_interaction_df$Interaction_Partner))

# Create a dataframe of nodes with type (Target Gene vs Interacting Protein)
node_df <- data.frame(
  Name = nodes,
  Type = ifelse(nodes %in% target_genes, "Target Gene", "Interacting Protein")
)

# Create edges (interactions)
edges <- data.frame(
  from = string_interaction_df$Gene,
  to = string_interaction_df$Interaction_Partner,
  weight = string_interaction_df$Confidence_Score
)

# Create a graph object
graph <- graph_from_data_frame(d = edges, vertices = node_df, directed = FALSE)

# Assign colors based on node type
node_colors <- ifelse(V(graph)$name %in% target_genes, "red", "blue")

# Generate a colorful layout for better visibility
set.seed(123)  # Ensures reproducibility
plot_layout <- create_layout(graph, layout = "fr")  # Fruchterman-Reingold layout

# Define a color gradient for edges
edge_colors <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)

# Plot the graph using ggraph
ggraph(plot_layout) +
  geom_edge_link(aes(edge_alpha = weight, color = weight), show.legend = TRUE) +
  geom_node_point(aes(color = Type), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, size = 5) +  # Fixed column name
  scale_color_manual(values = c("Target Gene" = "red", "Interacting Protein" = "blue")) +
  scale_edge_color_gradientn(colors = edge_colors) +
  theme_void() +
  ggtitle("Protein Interaction Network")

