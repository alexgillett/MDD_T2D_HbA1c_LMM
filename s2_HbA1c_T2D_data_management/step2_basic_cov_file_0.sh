############################################################################
### Step 2.0: Create a files in the directory where extracted HbA1c/ T2D data
### is being stored, so that UK Biobank assessment data fields can be
### saved for use in data extraction via ukbkings
############################################################################
cd path_where_you_to_store_extracted_data
cd hba1c_dat
### Create file to store month and year of birth field identifiers from UK Biobank
touch month_year_birth.txt
### Create file to store field info for self-reported phenotypes
touch self_report_diabetes_phenotypes.txt




