####################################################################################################
####### GTEx LUNG RNAseq PREPROCESSING AND FILTERING PIPELINE
####### Filtering for female subjects, selected age groups (20-39 and 60-79) and Hardy scale (0 = Ventilator)
####################################################################################################

####################################################################################################
##### 0. PRELIMINARY INFORMATION 
####################################################################################################
# This script was developed in RStudio using R version 4.4.1 and Bioconductor version 3.19.

# Set the working directory (uncomment and modify as needed).
# setwd("/Users/vanesacepaslopez/Documents/TFM/GTEx/")

# The file 'GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS_lung.txt' contains metadata 
# for both men and women. We want to use only women for RNAseq analysis.
# Subjects are grouped by age in 20-29, 30-39, 40-49, 50-59, 60-69 and 70-79.
# Since the number of young subjects is low, we will create a new group 20-39.
# We will also create a new group 60-79.
# We want to use only 20-39 and 60-79 age groups.

# Objective:
# - Filter GTEx metadata and RNA-seq count matrix to include only female subjects.
# - Further restrict to specific age groups: 20–39 and 60–79.
# - Save processed metadata and count matrix for downstream RNA-seq analysis.

####################################################################################################
##### 1. LOADING REQUIRED PACKAGES 
####################################################################################################

library(tidyverse)

####################################################################################################
##### 2. READING DATA 
####################################################################################################

# Load RNA-seq count matrix from the .gct file.
# This file contains gene expression counts with genes as rows and samples as columns.
# Skip the first two lines (metadata headers), and treat the third line as the column header.
count_matrix <- read.delim("gene_reads_lung.gct", skip = 2, header = TRUE)

# Load subject metadata containing phenotypic annotations (e.g., sex, age).
# Load metadata file with phenotypic information for GTEx subjects.
# Each row represents a donor (subject) and includes information about sex, age and dthhrdy.
# We disable stringsAsFactors to keep character columns as strings (not factors).
subjects <- read.delim("GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS_lung.txt", header = TRUE, stringsAsFactors = FALSE)

####################################################################################################
##### 3. FILTERING METADATA BY SEX (ONLY FEMALES)
####################################################################################################

# Filter to include only female subjects (SEX == 2).
subjects_female <- subjects %>% filter(SEX == 2)

####################################################################################################
##### 4. SAVING FILTERED METADATA (ONLY WOMEN)
####################################################################################################

# Save the filtered metadata (with female subjects only) to a tab-delimited text file.
write.table(subjects_female, "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS_lung_women.txt", sep = "\t", row.names = FALSE, quote = FALSE)

####################################################################################################
##### 5. RECODING AGE GROUPS AND FILTERING METADATA BY AGE
####################################################################################################

# Recode AGE to group 20–29 and 30–39 into "20-39", and 60–69 and 70–79 into "60-79".
subjects_female <- subjects_female %>%
  mutate(AGE_GROUP = case_when(
    AGE %in% c("20-29", "30-39") ~ "20-39",
    AGE == "40-49"              ~ "40-49",
    AGE == "50-59"              ~ "50-59",
    AGE %in% c("60-69", "70-79") ~ "60-79",
    TRUE                        ~ NA_character_  # catch any others
  ))

# Filter metadata to include only subjects in the 20–39 and 60–79 age groups.
subjects_female_age <- subjects_female %>% filter(AGE_GROUP %in% c("20-39", "60-79"))

####################################################################################################
##### 6. FILTERING METADATA BY HARDY SCALE (Cause of death)
####################################################################################################

# Filter metadata to include only subjects with a 0 value in the Hardy scale = ventilator case.
subjects_female_age_death <- subjects_female_age %>% filter(DTHHRDY %in% c("0"))


####################################################################################################
##### 7. SAVING FILTERED METADATA (FEMALES 20-39 AND 60-79)
####################################################################################################

# Save the filtered metadata (with 20-39 and 60-79 years old subjects and Hardy sclae = 0 only) to a tab-delimited text file.
write.table(subjects_female_age, "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS_lung_women_20-39_60-79_Hardy0.txt", sep = "\t", row.names = FALSE, quote = FALSE)

####################################################################################################
##### 8. EXTRACT SUBJECT IDs FROM SAMPLE NAMES IN COUNT MATRIX
####################################################################################################

# Sample names in the count matrix are formatted as: GTEX.111CU.0326.SM.5GZXO
# We extract the subject portion: GTEX.111CU

# Extract sample IDs from the count matrix.
# These are the column names from the 3rd column onward (first two columns are gene ID info).
count_matrix_sample_ids <- colnames(count_matrix)[-(1:2)]

# Define a function to extract the subject ID from the sample ID by keeping only
# the first two components (e.g., "GTEX.111CU" from "GTEX.111CU.0326.SM.5GZXO").
extract_subject_id <- function(x) {
  parts <- unlist(strsplit(x, "\\."))
  paste(parts[1], parts[2], sep = ".")
}

# Apply the function to all sample IDs from the count matrix
# to get extract corresponding subject IDs for each sample.
count_matrix_subject_ids <- sapply(count_matrix_sample_ids, extract_subject_id)

####################################################################################################
##### 9. FILTERING COUNT MATRIX TO MATCH SELECTED SUBJECTS
####################################################################################################

# Get vector of subject IDs to keep from the filtered metadata
subjects_to_keep <- subjects_female_age_death$SUBJID_dot

# Keep only sample columns from the count matrix that match the filtered subjects
keep_sample_cols <- count_matrix_sample_ids[count_matrix_subject_ids %in% subjects_to_keep]

# Subset the count matrix to include only selected columns
filtered_counts <- count_matrix[, c("Name", "Description", keep_sample_cols)]

####################################################################################################
##### 10. SAVE FILTERED COUNT MATRIX
####################################################################################################

# Save the filtered RNA-seq count matrix to a new file
write.table(filtered_counts, "gene_reads_lung_female_20-39_60-79_Hardy0.gct", sep = "\t", row.names = FALSE, quote = FALSE)
