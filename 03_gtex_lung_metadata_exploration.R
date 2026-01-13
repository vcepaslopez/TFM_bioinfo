####################################################################################################
####### GTEx LUNG DONOR METADATA EXPLORATION
####### Distribution of Sex, Age, Hardy Scale, and Ventilator Use
####################################################################################################

####################################################################################################
##### 0. PRELIMINARY INFORMATION
####################################################################################################

# This script was developed in RStudio using R version 4.4.1.
# It explores GTEx lung donor metadata to understand the distribution of sex, age, cause of death (Hardy scale),
# and ventilator usage, including visualizations and basic statistical comparisons.

# Set the working directory (uncomment and modify as needed).
setwd("/Users/vanesacepaslopez/Documents/TFM/GTEx/")

####################################################################################################
###### 1. INSTALLING REQUIRED PACKAGES  
####################################################################################################

# Installing the necessary R packages if not already installed
#### ggplot2
# https://ggplot2.tidyverse.org/
# Provides a flexible grammar of graphics for data visualization in R.
# Used for creating publication-quality figures, including PCA plots,
# volcano plots, and enrichment visualizations.

if (!("ggplot2" %in% installed.packages())) {
  install.packages("ggplot2")
}

#### dplyr
# https://dplyr.tidyverse.org/
# Provides a grammar of data manipulation in R.
# Used for filtering, selecting, mutating, summarizing, and grouping data frames,
# which simplifies data preprocessing and preparation for analysis and visualization.

if (!("dplyr" %in% installed.packages())) {
  install.packages("dplyr")
}

#### tidyr
# https://tidyr.tidyverse.org/
# Provides tools for tidying data in R.
# Used for reshaping data frames, including pivoting, unnesting, and filling missing combinations,
# making datasets easier to work with in downstream analyses and plotting.

if (!("tidyr" %in% installed.packages())) {
  install.packages("tidyr")
}

####################################################################################################
##### 2. LOADING REQUIRED PACKAGES 
####################################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)

####################################################################################################
##### 3. READING AND PREPROCESSING SUBJECT METADATA
####################################################################################################

# Load subject metadata containing phenotypic annotations (e.g., sex, age).
# Load metadata file with phenotypic information for GTEx subjects.
# Each row represents a donor (subject) and includes information about sex, age and dthhrdy.
# We disable stringsAsFactors to keep character columns as strings (not factors).

subjects <- read.delim("GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt", header = TRUE, stringsAsFactors = FALSE)

# The column SUBJID is Subject ID, GTEx Public Donor ID and looks like GTEX-111CU.
# In the .gct file I have one column per sample and the name looks like GTEX.111CU.0326.SM.5GZXO.
# The first part of the name describes the SubjectID.
# So I will now create a new column in the Subject metadata file to include the SubjectID with a . instead of a -
  
subjects$SubjectID_dot <- gsub("-", ".", subjects$SUBJID)

####################################################################################################
##### 4. READING EXPRESSION DATA AND MATCHING TO SUBJECTS
####################################################################################################

# Now I load bulk RNA-seq count matrix (.gct format)
# The first two columns are gene information (Name, Description), columns 3+ are samples

gct <- read.delim("gene_reads_lung.gct", skip = 2, header = TRUE)

# I create a variable with the sample IDs which are from column 3 onwards:
gct_sample_ids <- colnames(gct)[-(1:2)]

# Function to extract subject ID from sample name
extract_subject_id <- function(x) {
  parts <- unlist(strsplit(x, "\\."))
  paste(parts[1], parts[2], sep = ".")
}

# Apply function to get subject IDs corresponding to each sample
gct_subject_ids <- sapply(gct_sample_ids, extract_subject_id)

# Identify samples and subjects present in both expression data and metadata
matched_samples <- gct_sample_ids[gct_subject_ids %in% subjects$SubjectID_dot]

matched_subjects <- unique(gct_subject_ids[gct_subject_ids %in% subjects$SubjectID_dot])

# Subset metadata to matched subjects
matched_metadata <- subjects[subjects$SubjectID_dot %in% matched_subjects, ]

# Save matched metadata for later use
write.table(matched_metadata, "matched_metadata_sex_age.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# I load the cleaned metadata, i.e. the file I just created
meta <- read.delim("matched_metadata_sex_age.txt")

# Remove any rows with missing data
meta <- meta[complete.cases(meta), ]

####################################################################################################
##### 5. EXPLORATORY ANALYSIS: SEX, AGE, AND HARDY SCALE
####################################################################################################

# 5.1 SEX

# SEX: 1 = Male, 2 = Female (standard GTEx coding)
table(meta$SEX)
#Result
#1   2 
#414 189

# Percentages of each sex
prop.table(table(meta$SEX))
#1         2 
#0.6865672 0.3134328

# Convert SEX to factor with descriptive labels for plotting
meta$SEX <- factor(meta$SEX, levels = c(1, 2), labels = c("Male", "Female"))

# Create the bar plot: distribution of sex in GTEx lung donors
ggplot(meta, aes(x = SEX, fill = SEX)) +
  geom_bar(width = 0.6) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    vjust = -0.4,
    size = 5
  ) +
  scale_fill_manual(
    values = c(
      Male = "#00BFC4",
      Female = "#F8766D"
    )
  ) +
  labs(
    title = "Sex Distribution of GTEx Lung Donors",
    x = "Sex",
    y = "Number of donors",
    fill = "Sex"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 16),
    plot.title = element_text(size = 20, hjust = 0.5),  # plot title
    axis.title.x = element_text(size = 18),              # x-axis title
    axis.title.y = element_text(size = 18),              # y-axis title
    axis.text.x = element_text(size = 16),               # x-axis labels
    axis.text.y = element_text(size = 16),               # y-axis labels
    #legend.text = element_text(size = 16)               # legend text (if used)
  )

# 5.2 AGE

table(meta$AGE)
#Result
#20-29 30-39 40-49 50-59 60-69 70-79 
#35    44    93   211   198    22 

# Percentages per age group
prop.table(table(meta$AGE))  # Percentages
#20-29      30-39      40-49      50-59      60-69      70-79 
#0.05804312 0.07296849 0.15422886 0.34991708 0.32835821 0.03648425  

# Make sure AGE_GROUP is a factor with all original bins in the desired order
meta$AGE <- factor(meta$AGE, levels = c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79"))

# Define a palette of blue tones (light to dark) for age groups
blue_palette <- c("lightblue", "skyblue", "deepskyblue", "dodgerblue", "royalblue", "navy")

# Create the bar plot: age distribution in GTEx lung donors
ggplot(meta, aes(x = AGE, fill = AGE)) +
  geom_bar(width = 0.6) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    vjust = -0.4,
    size = 5
  ) +
  scale_fill_manual(values = blue_palette) + 
  labs(
    title = "Age Distribution of GTEx Lung Donors",
    x = "Age group",
    y = "Number of donors",
    fill = "Age group"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20, hjust = 0.5),  # plot title
    axis.title.x = element_text(size = 18),              # x-axis title
    axis.title.y = element_text(size = 18),              # y-axis title
    axis.text.x = element_text(size = 16),               # x-axis labels
    axis.text.y = element_text(size = 16),               # y-axis labels
    #legend.text = element_text(size = 16)               # legend text (if used)
  )

# 5.3 SEX BY AGE GROUP

# Contingency table: count of males/females in each age group
table_sex_age <- table(meta$SEX, meta$AGE)
rownames(table_sex_age) <- c("Male", "Female")

# Grouped bar plot: age distribution by sex
ggplot(meta, aes(x = AGE, fill = SEX)) +
  geom_bar(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    position = position_dodge(width = 0.7),
    vjust = -0.4,
    size = 5
  ) +
  scale_fill_manual(
    values = c(
      Male = "#00BFC4",
      Female = "#F8766D"
    )
  ) +
  labs(
    title = "Age Distribution of GTEx Lung Donors by Sex",
    x = "Age Group",
    y = "Count",
    fill = "Sex"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 20, hjust = 0.5),  # plot title
    axis.title.x = element_text(size = 18),              # x-axis title
    axis.title.y = element_text(size = 18),              # y-axis title
    axis.text.x = element_text(size = 16),               # x-axis labels
    axis.text.y = element_text(size = 16),               # y-axis labels
    legend.title = element_text(size = 16),              # legend title size
    legend.text = element_text(size = 14)               # legend text (if used)
  )

#5.4 HARDY SCALE (Cause/Time of Death)

# DTHHRDY:  
# "0" = "Ventilator",
# "1" = "Fast death\n(violent injury)",
# "2" = "Fast death\n(natural causes)",
# "3" = "Intermediate death",
# "4" = "Slow death"

table(meta$DTHHRDY)
#Result
#0   1   2   3   4 
#317  26 165  32  63

# Percentages per Hardy scale
prop.table(table(meta$DTHHRDY))
#0          1          2          3          4 
#0.52570481 0.04311774 0.27363184 0.05306799 0.10447761 

# Convert DTHHRDY to factor
meta$DTHHRDY <- factor(meta$DTHHRDY, levels = c(0, 1, 2, 3, 4))
#labels = c("0 – Ventilator case", "1 – Violent & fast death", "2 – Fast natural death", "3 – Intermediate death", "4 – Slow chronic death")
#)

# Define green palette for Hardy scale
green_palette <- c("lightgreen", "palegreen3", "mediumseagreen", "forestgreen", "darkgreen")

# Create bar plot: Hardy scale distribution
ggplot(meta, aes(x = DTHHRDY, fill = DTHHRDY)) +
  geom_bar(width = 0.6) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    vjust = -0.4,
    size = 5
  ) +
  scale_fill_manual(values = green_palette) +
  labs(
    title = "Hardy Scale Distribution of GTEx Lung Donors",
    x = "Hardy scale",
    y = "Number of donors",
    fill = "Hardy scale"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20, hjust = 0.5),  # plot title
    axis.title.x = element_text(size = 18),              # x-axis title
    axis.title.y = element_text(size = 18),              # y-axis title
    axis.text.x = element_text(size = 16),               # x-axis labels
    axis.text.y = element_text(size = 16),               # y-axis labels
    #legend.text = element_text(size = 16)               # legend text (if used)
  )

# Create plot: Hardy scale by sex
ggplot(meta, aes(x = DTHHRDY, fill = SEX)) +
  geom_bar(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    position = position_dodge(width = 0.7),
    vjust = -0.4,
    size = 5
  ) +
  scale_fill_manual(
    values = c(
      Male = "#00BFC4",
      Female = "#F8766D"
    )
  ) +
  labs(
    title = "Hardy Scale Distribution of GTEx Lung Donors by Sex",
    x = "Hardy Scale",
    y = "Count",
    fill = "Sex"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 20, hjust = 0.5),  # plot title
    axis.title.x = element_text(size = 18),              # x-axis title
    axis.title.y = element_text(size = 18),              # y-axis title
    axis.text.x = element_text(size = 16),               # x-axis labels
    axis.text.y = element_text(size = 16),               # y-axis labels
    legend.title = element_text(size = 16),              # legend title size
    legend.text = element_text(size = 14)               # legend text (if used)
  )

####################################################################################################
##### 6. FOCUS ON FEMALE DONORS
####################################################################################################

# Subset metadata to only female donors
meta_female <- meta[meta$SEX == "Female", ]

# Bar plot: age distribution of female donors
ggplot(meta_female, aes(x = AGE, fill = AGE)) +
  geom_bar(width = 0.6) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    vjust = -0.4,
    size = 5
  ) +
  scale_fill_manual(values = blue_palette) + 
  labs(
    title = "Age Distribution of GTEx Female Lung Donors",
    x = "Age group",
    y = "Number of donors",
    fill = "Age group"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20, hjust = 0.5),  # plot title
    axis.title.x = element_text(size = 18),              # x-axis title
    axis.title.y = element_text(size = 18),              # y-axis title
    axis.text.x = element_text(size = 16),               # x-axis labels
    axis.text.y = element_text(size = 16),               # y-axis labels
    #legend.text = element_text(size = 16)               # legend text (if used)
  )

# Bar plot: Hardy scale distribution of female donors
ggplot(meta_female, aes(x = DTHHRDY, fill = DTHHRDY)) +
  geom_bar(width = 0.6) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    vjust = -0.4,
    size = 5
  ) +
  scale_fill_manual(values = green_palette) +
  labs(
    title = "Hardy Scale Distribution of GTEx Female Lung Donors",
    x = "Hardy scale",
    y = "Number of donors",
    fill = "Hardy scale"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20, hjust = 0.5),  # plot title
    axis.title.x = element_text(size = 18),              # x-axis title
    axis.title.y = element_text(size = 18),              # y-axis title
    axis.text.x = element_text(size = 16),               # x-axis labels
    axis.text.y = element_text(size = 16),               # y-axis labels
    #legend.text = element_text(size = 16)               # legend text (if used)
  )   

# I create a bar plot: Hardy scale distribution by age group for female donors
ggplot(meta_female, aes(x = AGE, fill = DTHHRDY)) +
  geom_bar(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    position = position_dodge(width = 0.7),
    vjust = -0.4,
    size = 4
  ) +
  scale_fill_manual(values = green_palette) +
  labs(
    title = "Hardy Scale Distribution of Female GTEx Lung Donors by Age",
    x = "Age Group",
    y = "Count",
    fill = "Hardy Scale"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 20, hjust = 0.5),  # plot title
    axis.title.x = element_text(size = 18),              # x-axis title
    axis.title.y = element_text(size = 18),              # y-axis title
    axis.text.x = element_text(size = 16),               # x-axis labels
    axis.text.y = element_text(size = 16),               # y-axis labels
    legend.title = element_text(size = 16),              # legend title size
    legend.text = element_text(size = 14)               # legend text (if used)
  )

# I create agin the plot differently so this time a 0 is plotted when there are no samples in one group.                        
# Make sure AGE and DTHHRDY are factors with all levels
meta_female$AGE <- factor(meta_female$AGE, 
                          levels = c("20-29","30-39","40-49","50-59","60-69","70-79"))
meta_female$DTHHRDY <- factor(meta_female$DTHHRDY, levels = 0:4)

# Compute counts for all AGE x DTHHRDY combinations (fill missing with 0)
meta_female_complete <- meta_female %>%
  group_by(AGE, DTHHRDY) %>%
  summarise(n = n(), .groups = "drop") %>%
  complete(AGE, DTHHRDY, fill = list(n = 0))

# Plot age distribution of female donors by Hardy scale using precomputed counts geom_col()
ggplot(meta_female_complete, aes(x = AGE, y = n, fill = DTHHRDY)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = n),
            position = position_dodge(width = 0.7),
            vjust = -0.4,
            size = 5) +
  scale_fill_manual(values = green_palette) +
  labs(
    title = "Age Distribution of Female GTEx Lung Donors by Hardy Scale",
    x = "Age Group",
    y = "Count",
    fill = "Hardy Scale"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

####################################################################################################
##### 7. VENTILATOR USE IN FEMALE DONORS
####################################################################################################

# Create new binary variable: Ventilator use
meta_female$Ventilator <- ifelse(meta_female$DTHHRDY == 0, "Ventilator", "No ventilator")

# Bar plot: overall ventilator usage
barplot(table(meta_female$Ventilator), 
        main = "Ventilator Distribution",
        col = c("forestgreen", "palegreen3"),
        ylab = "Count")

# Ventilator use by age group
ventilator_age_female <- table(meta_female$Ventilator, meta_female$AGE)

# Bar plot:  ventilator use by age
barplot(ventilator_age_female,
        beside = TRUE,
        col = c("forestgreen", "palegreen3"),
        main = "Ventilator Distribution",
        xlab = "Age",
        ylab = "Count",
        legend = rownames(ventilator_age_female),
        las = 2)
