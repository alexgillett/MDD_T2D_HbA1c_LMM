############################################################################
### Step 4.0: Create a files in the directory where extracted HbA1c/ T2D data
### is being stored, so that UK Biobank assessment data fields can be
### saved for use in data extraction via ukbkings
############################################################################
cd path_to_where_you_store_your_analysis_output
cd hba1c_LMM
cd data
### Create file to store field identifiers from UK Biobank for required covariates
touch covariate_fields.txt
### Create file to store field identifiers from UK Biobank for initial blood pressure measurements
touch ukb_BPfield
### Create file to store field identifiers from UK Biobank for initial height measurements
touch ukb_hgt_field.txt
