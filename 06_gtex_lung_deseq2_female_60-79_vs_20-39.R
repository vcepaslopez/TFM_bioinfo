####################################################################################################
####### ANALYSIS OF GTEx LUNG RNAseq DATA
####### From Count Matrix to Differentially Expressed Genes Using DESeq2 
####################################################################################################

####################################################################################################
##### 0. PRELIMINARY INFORMATION 
####################################################################################################

# This script was developed in RStudio using R version 4.4.1 and Bioconductor version 3.19.
# It analyzes GTEx lung RNA-seq data from female subjects aged 20–39 vs. 60–79 using DESeq2.

# Set the working directory (uncomment and modify as needed).
#setwd("/Users/vanesacepaslopez/Documents/TFM/GTEx_lung_RNAseq/")

####################################################################################################
###### 1. INSTALLING REQUIRED PACKAGES  
####################################################################################################

# Installing the necessary R packages if not already installed

#### DESeq2
# https://www.bioconductor.org/packages/release/bioc/html/DESeq2.html
# Statistical package for differential gene expression analysis based on negative binomial
# distribution, designed for count data from RNA-seq experiments.
if (!("DESeq2" %in% installed.packages())) { 
  BiocManager::install("DESeq2");
}

#### org.Hs.eg.db
# https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html
# Human genome annotation package: maps gene identifiers (e.g., ENSEMBL, Entrez, SYMBOL)
# to metadata like gene names, functions, and pathways.
if (!("org.Hs.eg.db" %in% installed.packages())) { 
  BiocManager::install("org.Hs.eg.db");
}

#### AnnotationDbi
# https://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html
# Provides a framework for interfacing with Bioconductor annotation packages (e.g., org.Hs.eg.db),
# enables querying and mapping of gene identifiers.
if (!("AnnotationDbi" %in% installed.packages())) { 
  BiocManager::install("AnnotationDbi");
}

#### clusterProfiler
# https://www.bioconductor.org/packages/release/bioc/html/clusterProfiler.html
# Performs statistical analysis and visualization of functional profiles (e.g., GO, KEGG)
# of gene clusters or lists from differential expression results.
if (!("clusterProfiler" %in% installed.packages())) { 
  BiocManager::install("clusterProfiler");
}

####  EnhancedVolcano
# https://www.bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html
# Simplifies the creation of Volcano plots for visualizing differential expression results,
# highlighting significance and fold change. 
if (!("EnhancedVolcano" %in% installed.packages())) { 
  BiocManager::install("EnhancedVolcano");
}

#### ggplot2
# https://ggplot2.tidyverse.org/
# Provides a flexible grammar of graphics for data visualization in R.
# Used for creating publication-quality figures, including PCA plots,
# volcano plots, and enrichment visualizations.

if (!("ggplot2" %in% installed.packages())) {
  install.packages("ggplot2")
}

####################################################################################################
##### 2. LOADING REQUIRED PACKAGES 
####################################################################################################

# Loading the necessary packages
library("DESeq2")
library("org.Hs.eg.db")
library("AnnotationDbi")
library("clusterProfiler")
library("EnhancedVolcano")
library("ggplot2")

####################################################################################################
#### 3. READING AND PREPROCESSING DATA 
####################################################################################################

# We analyze GTEx lung RNA-seq data comparing female subjects aged 20–39 vs. 60–79.

# Load the count matrix (gene expression counts with genes as rows and samples as columns).
count_matrix <- read.delim("gene_reads_lung_female_20-39_60-79.gct", header = TRUE)

# Load lung subject metadata with phenotypic annotations. 
# Each row represents a donor (subject) and includes information about sex, age and dthhrdy.
subjects_metadata <- read.delim("GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS_lung_women_20-39_60-79.txt", 
                                header = TRUE, stringsAsFactors = FALSE)

# Check dimensions of the count matrix:
dim(count_matrix)
# [1] 59033   97
# 59033 genes x 97 columns (2 gene information columns + 95 samples)
# 70 donors 60-79 and 25 donors 20-39

# View data if needed:
# View(count_matrix)
# View(subjects_metadata)

# All samples are from females (SEX == 2), and only 20–39 and 60–79 age groups are included.
# DTHHRDY encodes cause of death using Hardy scale:
# 0: Ventilator Case; 1: Violent and fast death; 2: Fast death of natural causes; 3: Intermediate death; 4: Slow death

# Set AGE_GROUP as a factor, ordering levels to define the comparison direction: 20–39 (baseline) vs. 60–79
# so overexpressed genes (positive log2FoldChange) will have a higher expression value in older subjects (60-79).
subjects_metadata$AGE_GROUP <- factor(subjects_metadata$AGE_GROUP, levels = c("20-39", "60-79"))

####################################################################################################
#### 4. DIFFERENTIAL EXPRESSION ANALYSIS
####################################################################################################

# Differential expression analysis consists on three steps:
# 1. Filtering lowly expressed genes.
# 2. Normalization of data.
# 3. Running statistical models and tests for differential expression.


### 4.1. Create DESeqDataSet object
## ------------------------------------------------------------------------
# First, we create the DESeqDataSet object, specifying the experimental design.
# In this analysis we want to identify genes that are differentially expressed
# between samples of young women (ages 20-39) and old women (ages 60-79).
# With the function DESeqDataSetFromMatrix, we provide the count matrix (countData), 
# samples information (colData), and the experimental design (design).
# ?DESeq

# Extract sample IDs from count matrix (removing first two columns: Name and Description)
sample_ids <- colnames(count_matrix)[-(1:2)]

# Extract subject IDs from sample names (e.g., GTEX.1122O)
extract_subject_id <- function(x) paste(unlist(strsplit(x, "\\."))[1:2], collapse = ".")
sample_subject_ids <- sapply(sample_ids, extract_subject_id)

# Make sure metadata has SUBJID_dot column (subject IDs with dots)
# Then match and reorder rows based on the sample_subject_ids
colData <- subjects_metadata[match(sample_subject_ids, subjects_metadata$SUBJID_dot), ]

rownames(colData) <- sample_ids

# Set gene names as rownames and isolate count data, excluding the first two columns (Gene ID and Description)
rownames(count_matrix) <- count_matrix$Name
count_data <- count_matrix[, -(1:2)]

# Create DESeq2 dataset with AGE_GROUP as design variable
dds <- DESeqDataSetFromMatrix(countData = as.matrix(count_data),
                              colData = colData,
                              design = ~ AGE_GROUP)

# View experimental matrix to verify structure
# In the column AGE_GROUP60-79, a value of "1" indicates that the subject belongs to the 60-79 age group,
# while a value of "0" indicates that the subject belongs to the 20-39 age group.
model.matrix(~ AGE_GROUP, data = colData)

### 4.2. Filtering lowly expressed genes
## ------------------------------------------------------------------------

# We filter out genes with low expression across all samples
# to improve the statistical power of the differential expression test.
# Specifically, we keep genes that have at least 25 reads in total across all samples.
# That is, we retain rows in the count matrix where the sum of counts across all columns is 25 or higher.

keep <- rowSums(counts(dds)) >= 25
dds <- dds[keep,]

# Show how many genes are kept after filtering
table(keep)
dim(dds)
# Of the initial 59033 genes, 40036 genes are retained and the remaining 18997 are discarded.

### 4.3. Running DESeq and assessing dispersion
## ------------------------------------------------------------------------

# Now we need to identify which genes show significant changes in expression levels 
# between the two experimental conditions: older subjects (ages 60-79) and younger subjects (ages 20-39).

# Perform DESeq2 analysis using Wald test to identify differentially expressed genes
dds <- DESeq(dds)
res <- results(dds)
res
summary(res)

# We have fitted the data to a model and want to check that the data points are not too dispersed from it.
# The estimated curve should follow the general trend of the majority of the data
# and should not leave many points far from it or outside of it (outliers).
# To assess this, we are going to plot the dispersion. 

# Plot raw and normalized dispersion estimates to assess model fit.
DESeq2::plotSparsity(dds) # raw dispersion
DESeq2::plotDispEsts(dds) # normalized dispersion

# The red adjustment line represents the expected trend based on the data.
# In the estimation of dispersion, it is assumed that genes with similar expression levels
# should also have similar dispersion. Based on this assumption, the adjustment is performed. 
# During this process, values that are relatively far from the trend are pulled closer to the curve 
# using a method called shrinkage.

# Principal Component Analysis (PCA)

# Transform counts using VST ((variance stabilizing transformation) for downstream visualization (e.g., PCA)
# ?vst
vst_dds <- vst(dds)

# We perform a principal component analysis (PCA) on the data using the DESeq object that has already been transformed.
# Given a dataset with multiple variables, PCA reduces its dimensionality by identifying combinations of variables
# (principal components) that explain the highest possible proportion of variance. 

# Visualize separation between age groups via PCA

pca_data_age <- plotPCA(vst_dds, intgroup = "AGE_GROUP", ntop = nrow(dds),returnData = TRUE)

# Percentage of variance explained
percentVar <- round(100 * attr(pca_data_age, "percentVar"))

# Convert age group to factor (important for discrete colors)
pca_data_age$AGE_GROUP <- factor(
  pca_data_age$AGE_GROUP,
  levels = c("20-39", "60-79")
)

# Create PCA plot with ggplot2
ggplot(pca_data_age, aes(PC1, PC2, color = AGE_GROUP)) +
  geom_point(size = 3, alpha = 0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(
    color = "Age Group",
    title = "PCA by Age Group") +
  theme_classic() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 20, hjust = 0.5),   # plot title
    axis.title.x = element_text(size = 18),              # x-axis title
    axis.title.y = element_text(size = 18),              # y-axis title
    axis.text.x = element_text(size = 16),               # x-axis labels
    axis.text.y = element_text(size = 16),               # y-axis labels
    legend.title = element_text(size = 16),              # legend title size
    legend.text = element_text(size = 16)                # legend text
  )
# When we plot the PCA, we observe that the first principal component (X-axis) separates
# the samples by age group (20-39 and 60-79).


# Visualize separation between Hardy scale via PCA

# Compute PCA manually (same logic as plotPCA)
pca_data_dthhrdy <- plotPCA(vst_dds, intgroup = "DTHHRDY", ntop = nrow(dds), returnData = TRUE)

# Percentage of variance explained
percentVar <- round(100 * attr(pca_data_dthhrdy, "percentVar"))

# Convert Hardy scale to factor (important for discrete colors)
pca_data_dthhrdy$DTHHRDY <- factor(
  pca_data_dthhrdy$DTHHRDY,
  levels = c(0, 1, 2, 3, 4)
)

# Create PCA plot with ggplot2
ggplot(pca_data_dthhrdy, aes(PC1, PC2, color = DTHHRDY)) +
  geom_point(size = 3, alpha = 0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(
    color = "Hardy Scale",
    title = "PCA by Hardy Scale") +
  theme_classic() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 20, hjust = 0.5),   # plot title
    axis.title.x = element_text(size = 18),              # x-axis title
    axis.title.y = element_text(size = 18),              # y-axis title
    axis.text.x = element_text(size = 16),               # x-axis labels
    axis.text.y = element_text(size = 16),               # y-axis labels
    legend.title = element_text(size = 16),              # legend title size
    legend.text = element_text(size = 16)                # legend text
  )

# When we plot the PCA, we observe that the first principal component (X-axis) separates
# the samples by Hardy scale (0, ventilator and 1-4, the rest types of death).

### 4.4. Extracting differential expression results
## ------------------------------------------------------------------------

# Next, we obtain the matrix of differentially expressed genes.
# The first step is to decide which conditions we want to compare
# To do this, we visualize the internal names that DESeq has assigned to the possible comparisons,
# based on the factor levels we previously defined.

# List available result contrasts
resultsNames(dds)
# [1] "Intercept"                "AGE_GROUP_60.79_vs_20.39"

# Extract differential expression results for the comparison: 60–79 vs. 20–39,
# applying p-value adjustment using the FDR method.
res_trt <- results(dds, name = "AGE_GROUP_60.79_vs_20.39", pAdjustMethod = "fdr")

# We will calculate the number of up- and down-regulated genes.
# Define thresholds for statistical significance.
logfc_threshold <- 0
padj_threshold <- 0.1 

# log2FC = log2(A/B) = 0 means that the expression in A is the same as in B.
# We consider a gene to be differentially expressed if the adjusted p-value is less than 0.1.

# 4.4.1 IDENTIFYING UPREGULATED GENES (higher in 60–79 group)
# ------------------------------------------------------------------------

# We keep genes with a log2FC > 0 and an adjusted p-value < 0.1

table(res_trt$log2FoldChange > logfc_threshold & res_trt$padj < padj_threshold)
# We observe that 1593 genes are overexpressed in the 60-79 age group vs. 20-39 age group (logfc = 0).
# We observe that 1028 genes are overexpressed in the 60-79 age group vs. 20-39 age group (logfc = 0.5).
# We observe that 278 genes are overexpressed in the 60-79 age group vs. 20-39 age group (logfc = 1).

# Filter results for up-regulated genes (log2FC > 0 and padj < 0.1)
upregulated_genes <- res_trt %>% as.data.frame() %>% filter(log2FoldChange > logfc_threshold & padj < padj_threshold)

# View the top rows
head(upregulated_genes)

# Save the table to a file
write.table(upregulated_genes, "upregulated_genes_60-79_vs_20-39_logfc0.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# 4.4.2 IDENTIFYING DOWNREGULATED GENES (higher in 20–39 group)
# ------------------------------------------------------------------------

# We keep genes with a log2FC < 0 and an adjusted p-value < 0.1

table(res_trt$log2FoldChange < -logfc_threshold & res_trt$padj < padj_threshold)
# We observe that 923 genes are underexpressed in the 60-79 age group vs. 20-39 age group (logfc = 0).
# We observe that 444 genes are underexpressed in the 60-79 age group vs. 20-39 age group (logfc = 0.5).
# We observe that 103 genes are underexpressed in the 60-79 age group vs. 20-39 age group (logfc = 1).

# Filter results for down-regulated genes (log2FC > 0 and padj < 0.1)
downregulated_genes <- res_trt %>% as.data.frame() %>% filter(log2FoldChange < -logfc_threshold & padj < padj_threshold)

# View the top rows
head(downregulated_genes)

# Save the table to a file
write.table(downregulated_genes, "downregulated_genes_60-79_vs_20-39_logfc0.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# 4.4.3 IDENTIFYING HIGHLY UPREGULATED GENES (example)
# ------------------------------------------------------------------------

# Count highly upregulated genes (log2FC > 2, i.e., ≥ 4-fold change)

# Now, we calculate as an example, how many of the overexpressed genes are expressed at least 4 times more in 60-79 subjects (FC = 4).
# The logarithm base 2 of 4 is 2, so we create a variable to set the threshold:

logfc_threshold_fc4 <- 2

# We keep genes with a log2FC > 2 and an adjusted p-value < 0.1
table(res_trt$log2FoldChange > logfc_threshold_fc4 & res_trt$padj < padj_threshold)
# We observe that 17 genes are expressed at least 4 times more in the 60-69 age group vs. 20-39 age group.

# We keep genes with a log2FC < 2 and an adjusted p-value < 0.1
table(res_trt$log2FoldChange < -logfc_threshold_fc4 & res_trt$padj < padj_threshold)
# We observe that 10 genes are expressed at least 4 times less in the 60-69 age group vs. 20-39 age group.

####################################################################################################
#### 5. GENE ANNOTATION AND FUNCTIONAL ENRICHMENT
####################################################################################################

### 5.1. Gene Annotation
## ------------------------------------------------------------------------

# We perform annotation using the org.Hs.eg.db package from Bioconductor, which provides human gene annotations.
# We check the available column names for annotation:
# columns(org.Hs.eg.db)

# If we visualize the row names of the matrix, 
# we see that they correspond to the Ensembl gene IDs ("ENSEMBL") with the version number. 
# rownames(res_trt)
head(res_trt)

# We strip the version from the Ensembl IDs and add the Ensembl gene ID as a new column.
res_trt$ensembl <- gsub("\\.\\d+$", "", rownames(res_trt))
# head(res_trt)

# Using the mapIds function, we add different annotations to the data.
# We use the ENSEMBL identifier, which we already have, to retrieve the corresponding SYMBOL identifier from the 
# database. With the argument "multivals" set to "first", we specify that if multiple matches are found, 
# only the first match will be returned.
res_trt$symbol <- mapIds(org.Hs.eg.db, keys=res_trt$ensembl, column="SYMBOL", keytype="ENSEMBL", multiVals="first")

# Add gene name
res_trt$genename <- mapIds(org.Hs.eg.db, keys=res_trt$ensembl, column="GENENAME", keytype="ENSEMBL", multiVals="first")

# Add Entrez ID
res_trt$entrezid <- mapIds(org.Hs.eg.db, keys=res_trt$ensembl, column="ENTREZID", keytype="ENSEMBL", multiVals="first")

# We also add identifiers of all the pathways from the KEGG database in which
# each gene is involved, using "list"as the value for the multiVals argument.
res_trt$KEGG<-mapIds(org.Hs.eg.db, keys=res_trt$ensembl, column="PATH", keytype="ENSEMBL", multiVals="list")

# Sort by adjusted p-value
res_trt <- res_trt[order(res_trt$padj), ]
# head(res_trt)
# tail(res_trt)

# We convert KEGG column from list to character to write table
res_df <- as.data.frame(res_trt)
res_df[] <- lapply(res_df, as.character)

write.table(res_df, file = "rnaseq_women_60-79_vs_20-39.txt", row.names = FALSE, sep = "\t", quote = FALSE)

### 5.2. Functional Enrichment Analysis: Identifying Biological Processes, Molecular Functions, and Cellular Components
## ------------------------------------------------------------------------

# We perform a biological enrichment analysis to identify the enriched pathways or biological functions
# in the list of genes differentially expressed between the 60-79 and 20-39 age groups. 

# We set two thresholds to define the upregulated and downregulated genes:
# 1. The fold change must be greater than 2, which corresponds to a log2 fold change threshold of 1.
# 2. The adjusted p-value (padj) must be smaller than 0.1.

logfc_threshold <- 0.5
padj_threshold <- 0.1

# We create a list of Entrez IDs for upregulated and downregulated genes based on the defined thresholds:
deg <- list(up = res_trt$entrezid[which(res_trt$log2FoldChange > logfc_threshold & res_trt$padj < padj_threshold)], 
            down = res_trt$entrezid[which(res_trt$log2FoldChange < -logfc_threshold & res_trt$padj < padj_threshold)])

# We remove any "NA" entries from the list of Entrez IDs.
deg <- lapply(deg, na.omit)

# Now, we perform the Gene Ontology (GO) enrichment analysis for three categories: Biological Process (BP),
# Molecular Function (MF), and Cellular Component (CC) using the clusterProfiler package.

# Gene Ontology (GO) Enrichment: Biological Process (BP)
fea_GO_BP <- compareCluster(deg, 
                            ont = "BP",
                            fun = "enrichGO", 
                            OrgDb = org.Hs.eg.db)

# Gene Ontology (GO) Enrichment: Molecular Function (MF)
fea_GO_MF <- compareCluster(deg, 
                            ont = "MF",
                            fun = "enrichGO", 
                            OrgDb = org.Hs.eg.db)

# Gene Ontology (GO) Enrichment: Cellular Component (CC)
fea_GO_CC <- compareCluster(deg, 
                            ont = "CC",
                            fun = "enrichGO", 
                            OrgDb = org.Hs.eg.db)

# Visualize the top 10 GO terms enriched for each category (BP, MF, CC)
dotplot(fea_GO_BP, showCategory = 10, title="BP: Biological Process")
dotplot(fea_GO_MF, showCategory = 10, title="MF: Molecular Function")
dotplot(fea_GO_CC, showCategory = 10, title="CC: Cellular component")

# Save the plots into a single PDF file

pdf("GO_functional_enrichment_women_60-79_vs_20-39_logfc0.5_padj0.1.pdf", 
    width = 6, # width of the PDf in inches 
    height = 10 # height of the PDF in inches
)
dotplot(fea_GO_BP, showCategory = 10, title="BP: Biological Process")
dotplot(fea_GO_MF, showCategory = 10, title="MF: Molecular Function")
dotplot(fea_GO_CC, showCategory = 10, title="CC: Cellular component")

dev.off()

# Save the detailed results of the functional enrichment analysis for each GO category to text files
write.table(x = data.frame(fea_GO_BP), 
            file = "deseq2_results_GO_BP_enrichment_women_60-79_vs_20-39_logfc0.5_padj0.1.txt", 
            row.names = FALSE, sep = "\t", quote = FALSE)

write.table(x = data.frame(fea_GO_MF), 
            file = "deseq2_results_GO_MF_enrichment_women_60-79_vs_20-39_logfc0.5_padj0.1.txt", 
            row.names = FALSE, sep = "\t", quote = FALSE)

write.table(x = data.frame(fea_GO_CC), 
            file = "deseq2_results_GO_CC_enrichment_women_60-79_vs_20-39_logfc0.5_padj0.1.txt", 
            row.names = FALSE, sep = "\t", quote = FALSE)

# KEGG Pathway Enrichment
fea_KEGG <- compareCluster(geneCluster = deg, 
                           fun = "enrichKEGG", 
                           organism = "hsa")

# Visualize the top 10 enriched KEGG Pathways 
dotplot(fea_KEGG, showCategory = 10, title="KEGG Pathway Enrichment")

# Save the plots into a single PDF file

pdf("KEGG_functional_enrichment_women_60-79_vs_20-39_logfc0.5_padj0.1.pdf", 
    width = 6, # width of the PDf in inches 
    height = 10 # height of the PDF in inches
)
dotplot(fea_KEGG, showCategory = 10, title="KEGG Pathway Enrichment")

dev.off()

# Save the detailed results of the KEGG enrichment analysis to a text file
write.table(x = data.frame(fea_KEGG), 
            file = "deseq2_results_KEGG_enrichment_women_60-79_vs_20-39_logfc0.5_padj0.1.txt", 
            row.names = FALSE, sep = "\t", quote = FALSE)

####################################################################################################
# 6. VISUALIZATION: VOLCANO PLOT
####################################################################################################

# We will generate a volcano plot to visualize the differential expression analysis results.
# This plot shows the genes in terms of their log2 fold change (logFC) and adjusted p-value (padj),
# helping to identify significantly upregulated or downregulated genes.

# The EnhancedVolcano package is used to create the plot, and we annotate the genes by their symbol.
# ?EnhancedVolcano

pdf("volcano_plot_women_60-79_vs_20-39_logfc0.5_padj1e-5.pdf",
    width = 7, # width of the PDf in inches
    height = 6 # height of the PDF in inches
)
EnhancedVolcano(res_trt, # table of results
                lab = res_trt$symbol, # name of the gene to print
                x = "log2FoldChange", # logFC
                y = "pvalue", # p-value non adjusted
                pCutoff = 1e-5,  # threshold of non adjusted p-value
                FCcutoff = 0.5, # threshold of logFC,
                title = "Volcano plot", # title of the plot
                subtitle = "", # subtitle
                caption = paste0("total = ", nrow(res_trt), " genes"), # text for the caption
                labSize = 4.5, # size of the name
                pointSize = 2, # size of the point
                drawConnectors = TRUE # if we need to draw a line connecting name and point
)
dev.off()

#######################################################################################################
