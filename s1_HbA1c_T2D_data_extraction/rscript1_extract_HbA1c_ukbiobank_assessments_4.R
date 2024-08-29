############################################################################
### 1.4. HbA1c extraction UK Biobank biomarker assessments
############################################################################
### In this script we:
### 1. Extract some basic covariate information,
### 2. Extract HbA1c measurements from UK Biobank assessments and
### 3. Apply correction to HbA1c measurements (so that they can be used with
### primary care data).
### 4. Restructure extracted HbA1c biomarker data to be in long format
############################################################################
### Authors: Saskia Hagenaars, Alexandra C Gillett
############################################################################
# User add library path where R packages are stored if required:
#.libPaths(new = "your_rlibrary_path")
# Data extraction uses ukbkings package
# https://github.com/kenhanscombe/ukbkings
# ukbkings requires data to be formatted using ukbproject
# https://github.com/kenhanscombe/ukbproject
############################################################################
### Add R packages:
library(ukbkings)
library(tidyverse)
library(data.table)
### Set up directories required:
### For ukbkings, set project directory where UKB data is stored:
out_path <- "/path_where_you_to_store_extracted_data/hba1c_data/"
project_dir <- "/path_to_your_ukb_data_dir/"
############################################################################
### 1. Extract some basic covariate info from UKB data using ukbkings
############################################################################
f <- bio_field(project_dir)
############################################################################
### sex ("31-0.0"), "assessment_centre"="54-0.0", "age_at_assessment"="21003-0.0"
############################################################################
f %>%
select(field, name) %>%
filter(str_detect(field, c("54-0.0"))) %>%
slice(1) %>%
bio_field_add(paste(out_path, "covariates.txt", sep=""))
f %>%
select(field, name) %>%
filter(str_detect(field, c("31-0.0"))) %>%
slice(1) %>%
bio_field_add(paste(out_path, "covariates.txt", sep=""))
f %>%
select(field, name) %>%
filter(str_detect(field, c("21003-0.0"))) %>%
slice(1) %>%
bio_field_add(paste(out_path, "covariates.txt", sep=""))
### bio_phen will extract the required covariate data in covariates.txt and store it as a rds file:
bio_phen(
 project_dir,
 field = paste(out_path, "covariates.txt", sep=""),
 out = paste(out_path, "covariates", sep="")
)
############################################################################
### 2. Extract HbA1c biomarker data from UK Biobank assessments
### {Prior extractions were from primary care}
############################################################################
### a) Extract the biomarker fields
f %>%
select(field, name) %>%
filter(str_detect(name, "hba1c")) %>%
slice(1:4) %>%
bio_field_add(paste(out_path, "HbA1c_field_subset.txt", sep=""))
### b) Extract HbA1c biomarker data from UKB and store
bio_phen(
 project_dir,
 field = paste(out_path, "HbA1c_field_subset.txt", sep=""),
 out = paste(out_path, "HbA1c_phenotype_subset", sep="")
)
### c) Replace biomarker assay date with date of assessment (This is the date the sample was taken)
### Extract assessment date field and data
f %>%
select(field, name) %>%
filter(str_detect(field, c("53-0.0"))) %>%
slice(1) %>%
bio_field_add(paste(out_path, "assessment_date_field.txt", sep=""))
f %>%
select(field, name) %>%
filter(str_detect(field, c("53-1.0"))) %>%
slice(1) %>%
bio_field_add(paste(out_path, "assessment_date_field.txt", sep=""))
bio_phen(
 project_dir,
 field = paste(out_path, "assessment_date_field.txt", sep=""),
 out = paste(out_path, "assessment_date_subset", sep="")
)
#system("head assessment_date_subset.rds")
### Read in assessment date data:
date_asses <- readRDS(paste(out_path, "assessment_date_subset.rds", sep=""))
### Read in HbA1c biomarker data:
HbA1c_biomarker <- readRDS(paste(out_path, "HbA1c_phenotype_subset.rds", sep=""))
### Replace appropriate dates:
HbA1c_biomarker <- HbA1c_biomarker %>%
left_join(date_asses) %>%
select(1,2,3,6,7)%>%
rename(hba1c_bio_0 = 2,
  hba1c_bio_1 = 3,
  hba1c_date_0 = 4,
  hba1c_date_1 = 5)
############################################################################
### 3. Apply correction to HbA1c measurements (so that they can be used with
### primary care data).
############################################################################
## add correction to biomarker data as suggested by Katie Young (Exeter)
## correction: calibrated HbA1c = -0.9696 * original Hba1c + 3.3595
HbA1c_biomarker <- HbA1c_biomarker %>%
mutate(hba1c_bio_0_corrected = (0.9696 * hba1c_bio_0) + 3.3595,
  hba1c_bio_1_corrected = (0.9696 * hba1c_bio_1) + 3.3595)
############################################################################
### 4. Restructure extracted HbA1c biomarker data to be in long format and save
############################################################################
## split data to add the repeat visit as extra rows instead of columns
### Create a dataset for repeat visit. Remove NA values
HbA1c_biomarker_repeat_visit <- HbA1c_biomarker[!is.na(HbA1c_biomarker$hba1c_bio_1), ]
### Reduce columns to EID and repeat visit data and rename
HbA1c_biomarker_repeat_visit <- HbA1c_biomarker_repeat_visit %>%
select(eid, hba1c_bio_1, hba1c_bio_1_corrected, hba1c_date_1) %>%
rename(hba1c_bio_raw = hba1c_bio_1, hba1c = hba1c_bio_1_corrected, event_dt = hba1c_date_1)
### In original dataset remove NA for initial visit biomarkers, reduce to EID and initial visit columns and rename
### Then bind rows with repeat data to create long format
HbA1c_biomarker_separate <- HbA1c_biomarker %>%
filter(!is.na(hba1c_bio_0)) %>%
select(eid, hba1c_bio_0, hba1c_bio_0_corrected, hba1c_date_0) %>%
rename(hba1c_bio_raw = hba1c_bio_0, hba1c = hba1c_bio_0_corrected, event_dt = hba1c_date_0) %>%
bind_rows(HbA1c_biomarker_repeat_visit) %>%
mutate(hba1c_type="biomarker")

### Remove biomarker values > 195
HbA1c_biomarker_separate <- HbA1c_biomarker_separate[HbA1c_biomarker_separate$hba1c <= 195, ]

fwrite(HbA1c_biomarker_separate, file=paste(out_path, "/data/HbA1c_biomarker_rowwise.txt", sep=""),col.names=TRUE, row.names=FALSE, quote=FALSE)
