**Master’s Thesis Title**: “Characterising the Transcriptional Landscape of Lung Ageing in Women”

**Student**: Vanesa Cepas López

**Tutors**: Daniel Rico, Patricia Altea Manzano


**Abstract**:

This project investigates age-associated transcriptional changes in the human female lung using bulk RNA-seq data from the GTEx project. By comparing lung samples from young (20–39 years) and older (60–79 years) female donors, we identify widespread age-related gene expression changes affecting extracellular matrix organisation, immune and inflammatory signalling, and cellular proliferation. These changes suggest age-related remodelling of the lung microenvironment with potential relevance to metastatic colonisation. The study also demonstrates that mechanical ventilation prior to death strongly influences lung transcriptomic profiles, highlighting the importance of clinical metadata when analysing post-mortem transcriptomic datasets.


**Repository contents**:

**Curated reference tables**:

01_lung_scRNAseq_public_datasets_metadata.csv
Curated table of publicly available human lung single-cell RNA-seq datasets, integrating multiple studies and summarising key information including study name, experimental condition, age range, sex, sample type, number of donors, and data access details.
https://github.com/vcepaslopez/TFM_bioinfo/blob/main/01_lung_scRNAseq_public_datasets_metadata.csv 

**GTEx lung data preprocessing and metadata exploration**

02_gtex_lung_metadata_filter.R
Script to generate a filtered GTEx metadata file containing only donors with available lung bulk RNA-seq data.
https://github.com/vcepaslopez/TFM_bioinfo/blob/main/02_gtex_lung_metadata_filter.R

03_gtex_lung_metadata_exploration.R
Script to explore and visualise GTEx lung donor metadata, including distributions of sex, age, Hardy scale, and ventilator use.
https://github.com/vcepaslopez/TFM_bioinfo/blob/main/03_gtex_lung_metadata_exploration.R

**GTEx lung RNA-seq data filtering**

04_gtex_lung_filter_counts_female_age.R
Script to filter GTEx lung RNA-seq count matrices and metadata for female donors aged 20–39 and 60–79, and to save processed files for downstream analyses.
https://github.com/vcepaslopez/TFM_bioinfo/blob/main/04_gtex_lung_filter_counts_female_age.R

05_gtex_lung_filter_counts_female_age_hardy0.R
Script to filter GTEx lung RNA-seq count matrices and metadata for female donors aged 20–39 and 60–79 with Hardy scale 0 (mechanical ventilation prior to death), and to save processed files for downstream analyses.
https://github.com/vcepaslopez/TFM_bioinfo/blob/main/05_gtex_lung_filter_counts_female_age_hardy0.R

**Differential expression analyses**

06_gtex_lung_deseq2_female_60-79_vs_20-39.R
Script to perform differential gene expression analysis using DESeq2 on GTEx lung RNA-seq data from female donors aged 60–79 vs 20–39.
https://github.com/vcepaslopez/TFM_bioinfo/blob/main/06_gtex_lung_deseq2_female_60-79_vs_20-39.R

07_gtex_lung_deseq2_female_60-79_vs_20-39_hardy0.R
Script to perform differential gene expression analysis using DESeq2 on GTEx lung RNA-seq data from female donors with Hardy scale 0, comparing ages 60–79 vs 20–39.
https://github.com/vcepaslopez/TFM_bioinfo/blob/main/07_gtex_lung_deseq2_female_60-79_vs_20-39_hardy0.R

