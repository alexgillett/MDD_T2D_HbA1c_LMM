############################################################################
### 1.3. HbA1c extraction from GP records, plus prescription extractions
############################################################################
### In this script we:
### 1. Extract HbA1c data from the GP records, and,
### 2. Extract T2D medication data from prescription records with read3.
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
library(ukbkings)
library(tidyverse)
library(data.table)
### Set up directories required:
### For ukbkings, set project directory where UKB data is stored:
project_dir <- "/path_to_your_ukb_data_dir/"
inputs_dir <- "/path_to_provided_input_data/"
hba1c_data_output <- "/path_to_where_to_save_extracted_data/"
############################################################################
### 1. Extract HbA1c data from GP records
############################################################################
gp_clinical <- bio_record(project_dir, "gp_clinical")

readv2=fread(paste(inputs_dir, "mapping_read2_read3.txt", sep=""),header=T)
readv3=fread(paste(inputs_dir, "read3_code_term", sep=""),header=T)
### Previously readv2 rows 21:22 were blank.
### Just to be sure remove rows before starting, and save.
readv2 <-readv2[1:20, ]
fwrite(readv2, paste(inputs_dir, "mapping_read2_read3.txt", sep=""))

# separately extract read 2 and then map these to read3 code before extracting read3 and binding rows
# extract read 2 codes and add mapped read3 codes and terms
gp_HbA1c_read2 <- gp_clinical %>%
filter(read_2 %in% readv2$READV2_CODE) %>%
  left_join(readv2, by=c("read_2"="READV2_CODE")) %>%
  mutate(raw_HbA1c_results = ifelse(value1 == "OPR003", value2, value1)) %>%
  select(eid, event_dt, raw_HbA1c_results, read_2) %>%
  collect()
### ifelse(testexpression, value if true, value if false)

gp_HbA1c_read3 <- gp_clinical %>%
filter(read_3 %in% readv3$READV3_CODE) %>%
  left_join(readv3, by=c("read_3"="READV3_CODE")) %>%
  mutate(raw_HbA1c_results = ifelse(value1 == "OPR003", value2, value1)) %>%
  select(eid, event_dt, raw_HbA1c_results, read_3) %>%
  collect()

### save output:
fwrite(gp_HbA1c_read2, paste(hba1c_data_output, "gp_HbA1c_read2_extracted.txt", sep=""), col.names=T,row.names=F,quote=F,sep="\t")

fwrite(gp_HbA1c_read3, paste(hba1c_data_output, "gp_HbA1c_read3_extracted.txt"), sep="", col.names=T,row.names=F,quote=F,sep="\t")

############################################################################
### 2. Extract medication data (using drug names, and read codes if drug
### name is missing
############################################################################
gp_scripts <- bio_record(project_dir, "gp_scripts")

# drug names/codes
med_names = fread(paste(inputs_dir, "drug_names.txt", sep=""), header=T)
med_codes = fread(paste(inputs_dir, "drug_read2.txt", sep=""),header=T)

# extract medication based on drug name (this does not include Welsh data, as they have no drugname)
drug_data <- data.frame()
for (i in med_names$drug_name){
  print(paste("started ",i, "at", Sys.time()))
  temp <- grep(paste0("^",i,"."), gp_scripts[,drug_name],ignore.case=T) #extract rows with complete drugname
  temp2 <-  gp_scripts[temp,]
  drug_data <- rbind(drug_data,temp2)
  print(dim(drug_data))
  }
#Note: grep looks for pattern.
# pattern here is paste0("^", i, ".") =
# ^.drugname
#paste0("^", "vipidia", ".")
#case is ignored

drug_data <- drug_data %>%
left_join(med_names, by="drug_name") 

fwrite(drug_data, paste(hba1c_data_output, "gp_HbA1c_drug_data_names_extracted.txt", sep=""), col.names=T,row.names=F,quote=F,sep="\t")

rm(drug_data)
# extract medication data with missing drugnames (dataprovider = 4)
drug_read2 <- gp_scripts %>%
filter(read_2 %in% med_codes$read_code & data_provider==4)  %>%
left_join(med_codes, by=c("read_2"="read_code")) %>%
collect()

fwrite(drug_read2, paste(hba1c_data_output, "gp_HbA1c_drug_data_codes_extracted.txt", sep=""),col.names=T,row.names=F,quote=F,sep="\t")

