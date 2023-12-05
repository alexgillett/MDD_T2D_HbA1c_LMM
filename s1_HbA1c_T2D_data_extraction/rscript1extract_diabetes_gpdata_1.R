############################################################################
### 1.1. Code to extract diabetes cases in primary care
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
### Extract individuals with a diabetes diagnosis code in the
### GP records
############################################################################
gp_clinical <- bio_record(project_dir, "gp_clinical")

readv2=fread(paste(inputs_dir, "Diabetes_read2_codes", sep=""), header=T)
readv3=fread(paste(inputs_dir, "Diabetes_read3_codes", sep=""), header=T)

# separately extract read 2 and then map these to read3 code before extracting read3 and binding rows
# extract read 2 codes and add mapped read3 codes and terms
gp_DM_read2 <- gp_clinical %>%
  filter(read_2 %in% readv2$read_2) %>% 
  left_join(readv2, by=c("read_2"="read_2")) %>%
  select(!c(value2,value3)) %>%
  collect()

print(table(gp_DM_read2$read_2))

gp_DM_read3 <- gp_clinical %>%
  filter(read_3 %in% readv3$read_3) %>%
  left_join(readv3, by=c("read_3"="read_3")) %>%
  select(!c(value2,value3)) %>%
  collect()

print(table(gp_DM_read3$read_3))

fwrite(gp_DM_read2, paste(hba1c_data_output, "gp_DM_read2_extracted.txt", sep=""),col.names=T,row.names=F,quote=F,sep="\t")

fwrite(gp_DM_read3, paste(hba1c_data_output, "gp_DM_read3_extracted.txt", sep="") ,col.names=T,row.names=F,quote=F,sep="\t")
