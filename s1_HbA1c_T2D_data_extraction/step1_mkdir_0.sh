############################################################################
### Step 1.0
############################################################################
### Set up directory where you want to store extracted data
############################################################################
cd path_where_you_to_store_extracted_data
mkdir hba1c_dat
cd hba1c_dat
### Recommended to move content from github Input folder to subdirectory inputs
mkdir inputs
### Make subdirectory to store data
mkdir data
############################################################################
### Create files within hba1c_dat folder for the ukbkings package to save
### UK Biobank field identifiers in- these are then used for data extraction
############################################################################
### Create a file for some basic covariate fields:
touch covariates.txt
### Create a file for HbA1c biomarker fields:
touch HbA1c_field_subset.txt
### Create a file to store assessment date field info:
touch assessment_date_field.txt
############################################################################
### Additionally create directory space for Analysis output
############################################################################
cd ..
mkdir HbA1c_LMM
cd HbA1c_LMM
mkdir data
mkdir desc_tabs
mkdir PRSice_out
mkdir model_output
mkdir plots
mkdir scripts
