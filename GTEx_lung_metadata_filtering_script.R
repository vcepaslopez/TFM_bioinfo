####################################################################################################
####### CREATION OF GTEx LUNG METATADA FILE 
####### From the original GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt metadata file
####################################################################################################

####################################################################################################
##### 0. PRELIMINARY INFORMATION 
####################################################################################################

# This script was developed in RStudio using R version 4.4.1 and Bioconductor version 3.19.

# Set the working directory (uncomment and modify as needed).
# setwd("/Users/vanesacepaslopez/Documents/TFM/GTEx/")

# The original file 'GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt' contains metadata 
# for all GTEx subjects. However, not all subjects have corresponding lung RNA-seq data.
# The goal of this script is to create a new metadata file containing only those subjects 
# for whom lung RNA-seq data is available. 

####################################################################################################
##### 1. READING DATA 
####################################################################################################

# Load RNA-seq count matrix from the .gct file.
# This file contains gene expression counts with genes as rows and samples as columns.
# Skip the first two lines (metadata headers), and treat the third line as the column header.
gct <- read.delim("gene_reads_lung.gct", skip = 2, header = TRUE)

# Load subject metadata containing phenotypic annotations (e.g., sex, age).
# Load metadata file with phenotypic information for GTEx subjects.
# Each row represents a donor (subject) and includes information about sex, age and dthhrdy.
# We disable stringsAsFactors to keep character columns as strings (not factors).
subjects <- read.delim("GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt", header = TRUE, stringsAsFactors = FALSE)

####################################################################################################
##### 2. PROCESSING IDENTIFIERS FOR MATCHING 
####################################################################################################

# The metadata column SUBJID contains Subject IDs in the format: GTEX-111CU.
# In the count matrix, sample column names look like: GTEX.111CU.0326.SM.5GZXO.
# To match these formats, create a new column in metadata replacing '-' with '.' in SUBJID.
subjects$SUBJID_dot <- gsub("-", ".", subjects$SUBJID)

# Extract sample IDs from the count matrix.
# These are the column names from the 3rd column onward (first two columns are gene ID info).
gct_sample_ids <- colnames(gct)[-(1:2)]

# Define a function to extract the subject ID from the sample ID by keeping only
# the first two components (e.g., "GTEX.111CU" from "GTEX.111CU.0326.SM.5GZXO").
extract_subject_id <- function(x) {
  parts <- unlist(strsplit(x, "\\."))
  paste(parts[1], parts[2], sep = ".")
}

# Apply the function to all sample IDs from the count matrix
# to get extract corresponding subject IDs for each sample.
gct_subject_ids <- sapply(gct_sample_ids, extract_subject_id)

####################################################################################################
##### 3. FILTERING METADATA TO INCLUDE ONLY MATCHED SUBJECTS
####################################################################################################

# Subset the metadata to include only those subjects that are present in the lung RNAseq data.
# This ensures we are working only with metadata for subjects for whom we have gene expression data.
matched_metadata <- subjects[subjects$SUBJID_dot %in% gct_subject_ids, ]

####################################################################################################
##### 4. SAVING FILTERED METADATA
####################################################################################################

# Save the filtered metadata (with matched subjects only) to a tab-delimited text file.
write.table(matched_metadata, "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS_lung.txt", sep = "\t", row.names = FALSE, quote = FALSE)
