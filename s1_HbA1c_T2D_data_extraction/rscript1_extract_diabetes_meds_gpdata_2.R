############################################################################
### 1.2. Medication extraction for read2 codes from prescription records
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
project_dir <- "/path_to_your_ukb_data_dir/"
inputs_dir <- "/path_to_provided_input_data/"
hba1c_data_output <- "/path_to_where_to_save_extracted_data/"
############################################################################
### 1. Extract medication data
############################################################################
### Extracted using read2 codes

gp_scripts <- bio_record(project_dir, "gp_scripts")

readv2=fread(paste(inputs_dir, "read2drugs.txt", sep=""), header=T)

# separately extract read 2
gp_DM_meds_read2 <- gp_scripts %>%
  filter(read_2 %in% readv2$read_code) %>%
  left_join(readv2, by=c("read_2"="read_code")) %>%
  collect()

print(table(gp_DM_meds_read2$read_2))

fwrite(gp_DM_meds_read2, paste(hba1c_data_output, "gp_DM_meds_read2_extracted.txt", sep="",col.names=T,row.names=F,quote=F,sep="\t")
