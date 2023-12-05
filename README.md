# MDD_T2D_HbA1c_LMM
This repository contains code, typically R-scripts, to accompany research titled 'The impact of major depressive disorder on glycaemic control in type 2 diabetes: a longitudinal cohort study using UK Biobank primary care records'.

MDD = major depressive disorder. T2D = type 2 diabetes. HbA1c = glycated haemaglobin. LMM = linear mixed effects model. LDA = longitudinal data analysis. GP = general practioner.

## Overview of content

### s1_HbA1c_T2D_data_extraction 

R code to extract: HbA1c measurements from UKBiobank assessments and GP records, prescribed T2D medications from prescription data, and, diabetes diagnosis codes from GP records. Also contains an Inputs folder with files providing diabetes (Diabetes_read2_codes, Diabetes_read3_codes), HbA1c (mapping_read2_read3.txt, read3_code_term) and drug codes used in data extraction (drug_names.txt, drug_read2.txt and read2drugs.txt).

### s2_HbA1c_T2D_data_management

Post-extraction data management.

### s3_T2D_prs

Code used to perform T2D polygenic (risk) score analysis to validate the T2D phenotype.

### s4_covariate_extraction_LDA

R code to extract additional covariates used in longitudinal data analysis.

### s5_MI_Amelia

R code for multiple imputation using the Amelia package and for post-imputation data management.

### s6_LDA_run_models

R code for running the linear mixed effects models for the longitudinal data analysis work. This is for: (a) unadjusted model, (b) adjusted model and (c) within-individual variation analysis.

### s7_LDA_summary

R code for summarising and presenting the output from the linear mixed effects models.

## Description of content

Content is described in the order in which code should be run.

### 1. s1_HbA1c_T2D_data_extraction

All extraction is done using the ukbkings package in R. Please adapt scripts to your preferred UK Biobank data extraction/ management approach.

1.0. step1_mkdir_0: Script with code to: a) make directory to store inputs and outputs from data extraction, b) make files where fields to extracted are stored and c) make directory to store analysis output. 

1.1. rscript1_extract_diabetes_gpdata_1.R: Extract diabetes diagnosis codes from GP data. Inputs: Diabetes_read2_codes, Diabetes_read3_codes.

1.2. rscript1_extract_diabetes_meds_gpdata_2.R: Medication extraction from prescription records using read2 codes. Inputs: read2drugs.txt.

1.3. rscript1_extract_HbA1c_gpdata_3.R: (a) Extracts HbA1c measurements from GP records. Inputs: mapping_read2_read3.txt, read3_code_term. (b) Additional medication extraction code- extraction using drug name, and then drug code if drug name is unavailable. Inputs: drug_names.txt, drug_read2.txt.

1.4. rscript1_extract_HbA1c_ukbiobank_assessments_4.R: Extract HbA1c measurements from UK Biobank biomarker assessment data. Script also includes extraction of some basic covariate information (sex, assessment centre, age at assessment(s), assessment dates)

1.5. Before running data management scripts, MDD phenotypes needs to be extracted and defined. Please follow instructions found here: https://github.com/chiarafabbri/depression_phenotypes

### 2. s2_HbA1c_T2D_data_management

2.0. step2_basic_cov_file_0.sh: Script to create files for ukbkings package to store requested basic covariate data used in this step.

2.1. rscript2_HbA1c_T2D_data_management_1.R: Initial HbA1c and T2D phenotype data management script. Takes raw extracted HbA1c datasets, merges these to create a main long-formatted HbA1c dataset. Creates the T2D phenotype. Exclusion criteria for T2D case-control data applied. Phenotype dataset for use in PRSice (for phenotype validation) created.

### 3. s3_T2D_prs

3.1. PRSice_T2D_for_HbA1c_LMM.sh: Example batch script for PRSice with T2D case-control phenotype in UK Biobank data.

### 4. s4_covariate_extraction_LDA

4.0. step4_basic_cov_file_0.sh: Script to create files for ukbkings package to store requested basic covariate data used in this step.

4.1. rscript4_ukb_covariate_extraction_1.R: Variables used in imputation and/ or longitudinal data analysis are extracted and/ or created. These include BMI at T2D diagnosis, MDD exposures, etc. The dataset to use in imputation is also created. This is a long script!

### 5. s5_MI_Amelia

5.1. rscript5_amelia_imputation_1.R: Script to perform mulitple imputation using Amelia package.

5.2. rscript5_postMI_data_management_2.R: Script to perform post-imputation data management in preparation for longitudinal data analysis.

### 6. s6_LDA_run_models

#### LDA_unadjusted

Folder containing the script (R) to run the unadjusted linear mixed effects models. Analysis does not use imputed data. Code can be submitted as a batch job (user to write own job submission script). Models are run with both maximum likelihood (ML) estimation (for fixed effect model p-values) and restricted maximum likelihood (REML) estimation (for model parameter estimates).

#### LDA_adjusted 

Folder containing the script (R) to run the adjusted linear mixed effects model. This is the selected unadjusted model plus additional covariates. Analysis uses imputed data. Code can be submitted as batch jobs (user to write own job submission script). Model-code split into 5 scripts- each script runs the same models but in different imputation datasets in subgroups of 10 MI datasets per script. Models are run with both maximum likelihood (ML) estimation (for fixed effect model p-values) and restricted maximum likelihood (REML) estimation (for model parameter estimates).

#### LDA_variance

Folder containing the script (R) to run the adjusted linear mixed effects model with MDD allowed to impact within-individual variation. Analysis uses imputed data. Three models considered: 1) pre-T2D MDD (binary) impacting residual HbA1c variation (REMLvar_pret2d_diff scripts), 2) time-varying post-T2D MDD (binary) impacting residual HbA1c variation (REMLvar_postt2dTV_diff scripts), and 3) retrospective post-T2D MDD (binary) impacting residual HbA1c variation (REMLvar_postt2d_retro_diff scripts). Code can be submitted as batch jobs (user to write own job submission script). Model-code split into 5 scripts- each script runs the same models but in different imputation datasets in subgroups of 10 MI datasets per script. Models are run using restricted maximum likelihood (REML) estimation only since variance parameters are those of interest here. 


### 7. s7_LDA_summary

7.1. rscript7_unadjusted_summary_1.R: Script selecting and summarising the unadjusted model. Contains code to generate plots used in manuscript (HbA1c over time, using model, for example pre- and post-T2D MDD individuals).

7.2. rscript7_adjusted_summary_2.R: Script to summarise pooled output from the adjusted model run in 50 MI dataset.

7.3. rscript7_MDD_variance_3.R: Script to summarise pooled variance parameter estimates, and perform hypothesis testing, for the three models considered: 1) pre-T2D MDD (binary) impacting residual HbA1c variation, 2) time-varying post-T2D MDD (binary) impacting residual HbA1c variation, and 3) retrospective post-T2D MDD (binary) impacting residual HbA1c variation. 
