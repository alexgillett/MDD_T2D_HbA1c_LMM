############################################################################
### 4.1. Covariate extraction and data reformating for LDA/ MI
### LDA = longitudinal data analysis, MI = multiple imputation
############################################################################
### Authors: Alexandra C Gillett
############################################################################
### In script we do a lot... Strap in!
### 1. Extract required UK Biobank assessment data
### 2. Extract BMI from primary care records
### 3. Reformat BMI from UK Biobank assessments (extracted above)
### 4. Merge UK Biobank BMI with primary care BMI dataset
### 5. Extract weight from GP data
### 6. Extract height from GP data
### 7. Extract height from UK Biobank initial assessment
### 8. Create additional BMI measurements with weight (GP) and height (UKB)
### 9. Merge BMI datasets (remove duplicates)
### 10. Create basic HbA1c dataset for LDA
### 11. Apply pre-imputation exclusion critera
### 12. Create HbA1c observation count variables
### 13. Create BMI at T2D diagnosis (baseline) covariate
### 14. Create MDD status at baseline and medications at baseline
### 15. Add additional UK Biobank assessment variables to HbA1c dataset
### 16. Update coding for the medicine variables
### 17. Catch all section...
### 18. Create mean SBP and DBP variables using UK Biobank assessments
### 19. Create BMI measurment count columns (GP BMI only)
### 20. Extract BP measurements from GP records
### 21. Create BP at T2D diagnosis (baseline) covariates
### 22. Create additional BP variables (i.e. count of observations)
### 23. Create retrospective MDD grouping variable
### 24. Standardise continuous imputation and LDA variables
### 25. MDD exposures (duration variables created)
### 26. Maximum follow-up calculated
### 27. Alcohol never created
### 28. Add initial BMI at UK Biobank assessment into LDA dataset
### 29. Sum of HbA1c, BP and BMI measurements before T2D diagnosis
### 30. Standardise pret2dmdd_time (robust scaler standardising)
### 31. Subgroup variable for pre-T2D MDD (exploratory)
### 32. Splines for time since post-T2D MDD diagnosis
### 33. Plot of time between pre-T2D MDD diagnosis and T2D diagnosis
### 34. Percentage diagnosed outside of 2006 and 2010
############################################################################
# User add library path where R packages are stored if required:
#.libPaths(new = "your_rlibrary_path")
# Data extraction uses ukbkings package
# https://github.com/kenhanscombe/ukbkings
# ukbkings requires data to be formatted using ukbproject
# https://github.com/kenhanscombe/ukbproject
############################################################################
### Add R packages:
library(tidyverse)
library(data.table)
library(lubridate)
library(readxl)
library(knitr)
library(xtable)
library(ggplot2)
library(readr)
library(zoo)
library(ukbkings)
### Assign paths to relevant directories
hba1c_dat_path <- "/path_where_you_to_store_extracted_data/hba1c_data/"
project_dir <- "/path_to_your_ukb_data_dir/"
depression_dir <- "/path_to_extracted_depression_data/"
analysis_dir <- "/path_to_analysis_dir/"

############################################################################
### 1. Extract covariates from UK Biobank assessment data
############################################################################
f <- bio_field(project_dir)
### Covariates may want are:
### covariates data set with sex ("31-0.0"), "assessment_centre"="54-0.0", genetic batch == "22000",
### genetic ethnic background == "22006", self report ethnic background == "21000",
### bmi ("21001"), dates (for BMI/ instances) = 53
### Physical activity (measured at initial assessment). Overall MET mins summed per week == 22040,
### walking = 22037, mod exercise = 22038, vig exer = 22039
### Qualifications: 6138, age completed FT education = 845
### Comparative body size at age 10 = 1687
### smoking status "current" = 20116, ever smoked: 20160
### Alcohol intake frequency: 1558, alc. drinking status 20117
### tdi at UKB recruitment: 189
### Average total household income before tax 738
### Cooked vegetable intake 1289, salad/raw veg intake 1299, fresh fruit intake 1309, processed meat intake 1349
f %>%
    select(field, name) %>%
    filter(str_detect(name, "townsend"))
f %>%
    select(field, name) %>%
    filter(str_detect(field, "1289"))
covs_f <- c("31-0.0", "54-0.0", "22000-0.0", "22006-0.0", "21000-0.0", "21001-0.0", "21001-1.0", "21001-2.0",
"21001-3.0", "53-0.0", "53-1.0", "53-2.0", "53-3.0", "189-0.0", "22040-0.0", "22037-0.0", "22038-0.0", "22039-0.0", "6138-0.0", "6138-0.1", "6138-0.2",
"6138-0.3", "6138-0.4", "6138-0.5", "845-0.0", "845-1.0", "845-2.0", "1687-0.0", "1687-1.0", "1687-2.0", "20116-0.0", "20160-0.0",
"1558-0.0", "20117-0.0", "739-0.0", "1289-0.0", "1299-0.0", "1309-0.0", "1349-0.0")

### Add fields to file:
for(i in 1:length(covs_f)){
f %>%
select(field, name) %>%
filter(str_detect(field, covs_f[i])) %>%
slice(1) %>%
bio_field_add(paste(analysis_dir, "data/covariate_fields.txt", sep = ""))
}
### Extract fields from UK Biobank phenotype data:
bio_phen(
 project_dir,
 field = paste(analysis_dir, "data/covariate_fields.txt", sep=""),
 out = paste(analysis_dir, "data/covariate_lmm", sep="")
)
### Read in data frame of extracted covariate data
df <- readRDS(paste(analysis_dir, "data/covariate_lmm.rds", sep=""))
### Rename from fields to covariate names
df <- bio_rename(df, f)

### Additionally extract blood pressure measurements at initial visit
f %>%
    select(field, name) %>%
    filter(str_detect(name, "systolic"))
## 93-0.0 93-0.1 4080-0.0 4080-0.1 == systolic at initial
## 94-0.0 94-0.1 4079-0.0 4079-0.1 == diastolic at initial
BPfields <- c("93-0.0", "93-0.1", "4080-0.0", "4080-0.1", "94-0.0", "94-0.1", "4079-0.0", "4079-0.1")
for(i in 1:length(BPfields)){
f %>%
select(field, name) %>%
filter(str_detect(field, BPfields[i])) %>%
slice(1) %>%
bio_field_add(paste(analysis_dir, "data/ukb_BPfield.txt", sep=""))
}
bio_phen(
 project_dir,
 field = paste(analysis_dir, "data/ukb_BPfield.txt", sep=""),
 out = paste(analysis_dir, "data/ukb_BP", sep="")
)

############################################################################
### 2. Extract BMI from primary care records
############################################################################
t2d_wrk <- readRDS(paste(hba1c_dat_path, "data/hba1c_t2d_mdd_medication_clean.rds", sep=""))

sum(is.na(t2d_wrk$dob)) # 0
sum(is.na(t2d_wrk$first_t2d_date)) # 0
### Create age at T2D diagnosis
t2d_wrk <- t2d_wrk %>% mutate(age_t2d_diag = time_length(interval(as.Date(dob), as.Date(first_t2d_date)), "years"))
### Restrict to over 18
t2d_wrk <- t2d_wrk[t2d_wrk$age_t2d_diag > 18, ]
length(unique(t2d_wrk$eid)) # 17544

### Double check number of rows per individual...
#t2d_dt <- data.table(t2d_wrk, key="eid")
#out <- t2d_dt[, .(rowCount = .N), by=eid]
#any(out$rowCount < 2) # FALSE
#rm(t2d_dt)

### List of EIDs included in dataset
eids_t2d <- unique(t2d_wrk$eid)

### Start BMI extraction
gp_clinical <- bio_record(project_dir, "gp_clinical")
head(gp_clinical)
### Main code for BMI == 22K.. (read v2 and v3):
### Search v3
gp_bmi_read3_22k <- gp_clinical %>% 
    filter(str_detect(read_3, "22K..")) %>%
    collect()
dim(gp_bmi_read3_22k) # 1129676       8
gp_bmi_read3_22k <- gp_bmi_read3_22k[gp_bmi_read3_22k$eid %in% eids_t2d, ]
dim(gp_bmi_read3_22k)
# [1] 229645      8

### Search v2
gp_bmi_read2_22k <- gp_clinical %>% 
    filter(str_detect(read_2, "22K..")) %>%
    collect()
dim(gp_bmi_read2_22k) #  473277      8
gp_bmi_read2_22k <- gp_bmi_read2_22k[gp_bmi_read2_22k$eid %in% eids_t2d, ]
dim(gp_bmi_read2_22k)
# [1] 113049      8

### 2 additional read 3 codes mentioned online: Xa7wG and X76CO
#gp_bmi_read3_Xa7wG <- gp_clinical %>% 
#    filter(str_detect(read_3, "Xa7wG")) %>%
#    collect()
#dim(gp_bmi_read3_Xa7wG) # [1] 72  8 ### Values appear blank
#sum(gp_bmi_read3_Xa7wG$value1 == "") # 72
#sum(gp_bmi_read3_Xa7wG$value2 == "") # 72
#sum(gp_bmi_read3_Xa7wG$value3 == "") # 72
#rm(gp_bmi_read3_Xa7wG)
### Just double check for read 2...
#gp_bmi_read2_Xa7wG <- gp_clinical %>% 
#    filter(str_detect(read_2, "Xa7wG")) %>%
#    collect()
#dim(gp_bmi_read2_Xa7wG) # [1] 0  8 ### nothing there
#rm(gp_bmi_read2_Xa7wG)
#
#gp_bmi_read3_X76CO <- gp_clinical %>% 
#    filter(str_detect(read_3, "X76CO")) %>%
#    collect()
#dim(gp_bmi_read3_X76CO)
## [1] 0 8
###
### What are 'viable'/ reasonable values of BMI. What may be data input errors...
### The following article dropped BMI outside of range 5-200 kg/m2
### Bhaskaran K, Forbes HJ, Douglas I, et al Representativeness and optimal use of body mass index (BMI) in the UK Clinical Practice Research Datalink (CPRD) BMJ Open 2013;3:e003389. doi: 10.1136/bmjopen-2013-003389
gp_bmi_read2_22k
gp_bmi_read3_22k

### BMI data cleaning
### Start with read 2...
range(gp_bmi_read2_22k$value1) # ""       "98.000"
table(gp_bmi_read2_22k$read_2)

table(gp_bmi_read2_22k$value1)[1:10] ### has "^" included- set to ""
table(gp_bmi_read2_22k$value2) ### Value 2 can be ignored in this round... But check in future extractions.
table(gp_bmi_read2_22k$value3) ### "K/M2"    "kg/m2"  "Unknown" set to ""

gp_bmi_read2_22k$value1[gp_bmi_read2_22k$value1 == "^"] <- ""
gp_bmi_read2_22k$value3[gp_bmi_read2_22k$value3 == "K/M2"] <- ""
gp_bmi_read2_22k$value3[gp_bmi_read2_22k$value3 == "kg/m2"] <- ""
gp_bmi_read2_22k$value3[gp_bmi_read2_22k$value3 == "Unknown"] <- ""

gp_bmi_read2_22k_1 <- gp_bmi_read2_22k %>% 
    mutate(miss1 = as.numeric(value1 == ""), miss2 = as.numeric(value2 == ""),miss3 = as.numeric(value3 == ""))

### ignore value 2- so interested in miss1 and miss3
gp_bmi_read2_22k_1 <- gp_bmi_read2_22k_1 %>%
    mutate(sum_miss = miss1+miss3)

table(gp_bmi_read2_22k_1$sum_miss)
### both missing -> remove
#check <- gp_bmi_read2_22k_1[gp_bmi_read2_22k_1$sum_miss == 2, ]
#check
#rm(check)
gp_bmi_read2_22k_1 <- gp_bmi_read2_22k_1[gp_bmi_read2_22k_1$sum_miss != 2, ]
### 1 == keep as this means just one column has a value...
#check <- gp_bmi_read2_22k_1[gp_bmi_read2_22k_1$sum_miss == 0, ]
#check
# rm(check)
### 3 individuals with 2 values- hard to know which to use so rm
gp_bmi_read2_22k_1 <- gp_bmi_read2_22k_1[gp_bmi_read2_22k_1$sum_miss != 0, ]

### Make BMI numeric and combine BMI values from value1 and value3
gp_bmi_read2_22k_1 %>% mutate(bmi = as.numeric(value1)) -> gp_bmi_read2_22k_1
gp_bmi_read2_22k_1 %>% mutate(bmi3 = as.numeric(value3)) -> gp_bmi_read2_22k_1
sum(is.na(gp_bmi_read2_22k_1$bmi)) # 18621
sum(!is.na(gp_bmi_read2_22k_1$bmi3)) #18621
gp_bmi_read2_22k_1 %>% 
    mutate(bmi = ifelse(is.na(bmi), bmi3, bmi)) -> gp_bmi_read2_22k_1
sum(is.na(gp_bmi_read2_22k_1$bmi)) # 0

gp_bmi_read2_22k_1 %>% select(eid, data_provider, event_dt, bmi) -> gp_bmi_read2_22k_1
dim(gp_bmi_read2_22k_1) # 95970     4
length(unique(gp_bmi_read2_22k_1$eid)) # 5291
gp_bmi_read2_22k_1 <- gp_bmi_read2_22k_1[gp_bmi_read2_22k_1$bmi >= 5, ]
gp_bmi_read2_22k_1 <- gp_bmi_read2_22k_1[gp_bmi_read2_22k_1$bmi <=200, ]
dim(gp_bmi_read2_22k_1) # 95525     4
length(unique(gp_bmi_read2_22k_1$eid)) # 5291
any(is.na(gp_bmi_read2_22k_1$event_dt)) # F
any(gp_bmi_read2_22k_1$event_dt == "") # T
sum(gp_bmi_read2_22k_1$event_dt == "") # 4
gp_bmi_read2_22k_1 <- gp_bmi_read2_22k_1[gp_bmi_read2_22k_1$event_dt != "", ]
dim(gp_bmi_read2_22k_1) # 95521     4
length(unique(gp_bmi_read2_22k_1$eid)) # 5291

### move onto read3...

table(gp_bmi_read3_22k$value1)[1:100] ### all "number" or ""
table(gp_bmi_read3_22k$value2) ### Value 2 all ""... But check in future extractions.
table(gp_bmi_read3_22k$value3) ### Value 3 all ""... But check in future extractions.

gp_bmi_read3_22k_1 <- gp_bmi_read3_22k[gp_bmi_read3_22k$value1 != "", ]
gp_bmi_read3_22k_1 %>% mutate(bmi = as.numeric(value1)) -> gp_bmi_read3_22k_1
sum(is.na(gp_bmi_read3_22k_1$bmi)) # 0
range(gp_bmi_read3_22k_1$bmi) # 0.000 135.087

gp_bmi_read3_22k_1 %>% select(eid, data_provider, event_dt, bmi) -> gp_bmi_read3_22k_1

dim(gp_bmi_read3_22k_1) # 221878      4
length(unique(gp_bmi_read3_22k_1$eid)) # 12263

gp_bmi_read3_22k_1 <- gp_bmi_read3_22k_1[gp_bmi_read3_22k_1$bmi >= 5, ]
gp_bmi_read3_22k_1 <- gp_bmi_read3_22k_1[gp_bmi_read3_22k_1$bmi <=200, ]
dim(gp_bmi_read3_22k_1) # 218376      4
length(unique(gp_bmi_read3_22k_1$eid)) # 12261

any(is.na(gp_bmi_read3_22k_1$event_dt)) # F
any(gp_bmi_read3_22k_1$event_dt == "") # T
sum(gp_bmi_read3_22k_1$event_dt == "") # 91
gp_bmi_read3_22k_1 <- gp_bmi_read3_22k_1[gp_bmi_read3_22k_1$event_dt != "", ]
dim(gp_bmi_read3_22k_1) #  218285      4
length(unique(gp_bmi_read3_22k_1$eid)) # 12261

### Combine read2 and read3 BMI observations...
gp_bmi <- gp_bmi_read3_22k_1 %>% bind_rows(gp_bmi_read2_22k_1)
dim(gp_bmi) # 313806      4
length(unique(gp_bmi$eid)) # 17469
gp_bmi %>% mutate(event_dt = as.Date(event_dt, "%d/%m/%Y")) %>%
    distinct() -> gp_bmi1

dim(gp_bmi1) # 302422      4
length(unique(gp_bmi1$eid)) # 17469

gp_bmi2 <- gp_bmi1 %>% select(-data_provider)
rm(gp_bmi1)
saveRDS(gp_bmi2, file= paste(hba1c_dat_path, "data/GP_t2d_bmi.rds", sep=""))

############################################################################
### 3. Reformat BMI from UK Biobank assessments (extracted above)
############################################################################
df_bmi0 <- df %>% select(eid, date_of_attending_assessment_centre_f53_0_0, "body_mass_index_(bmi)_f21001_0_0")
df_bmi1 <- df %>% select(eid, date_of_attending_assessment_centre_f53_1_0, "body_mass_index_(bmi)_f21001_1_0")
df_bmi2 <- df %>% select(eid, date_of_attending_assessment_centre_f53_2_0, "body_mass_index_(bmi)_f21001_2_0")
df_bmi3 <- df %>% select(eid, date_of_attending_assessment_centre_f53_3_0, "body_mass_index_(bmi)_f21001_3_0")

df_bmi0 <- df_bmi0 %>% rename(event_dt = date_of_attending_assessment_centre_f53_0_0, bmi = "body_mass_index_(bmi)_f21001_0_0")
df_bmi1 <- df_bmi1 %>% rename(event_dt = date_of_attending_assessment_centre_f53_1_0, bmi = "body_mass_index_(bmi)_f21001_1_0")
df_bmi2 <- df_bmi2 %>% rename(event_dt = date_of_attending_assessment_centre_f53_2_0, bmi = "body_mass_index_(bmi)_f21001_2_0")
df_bmi3 <- df_bmi3 %>% rename(event_dt = date_of_attending_assessment_centre_f53_3_0, bmi = "body_mass_index_(bmi)_f21001_3_0")

df_bmi <- df_bmi0 %>%  
    bind_rows(df_bmi1) %>%  bind_rows(df_bmi2) %>%  bind_rows(df_bmi3)

df_bmi <- df_bmi[!is.na(df_bmi$bmi), ]

any(df_bmi$event_dt == "") # F
any(is.na(df_bmi$event_dt)) # F
### All BMI
df_bmi <- df_bmi%>% distinct()
### BMI restricted to T2D individuals
df_bmi1 <- df_bmi[df_bmi$eid %in% eids_t2d,]

saveRDS(df_bmi1, file=paste(hba1c_dat_path, "data/UKBassessments_t2d_bmi.rds", sep=""))
saveRDS(df_bmi, file=paste(hba1c_dat_path, "data/UKBassessments_bmi.rds", sep=""))

df_bmi1 <- df_bmi1 %>% mutate(event_dt = as.Date(event_dt, "%Y-%m-%d"))

############################################################################
### 4. Merge UK Biobank BMI with primary care BMI dataset
############################################################################
### merge gp and UKB assessment BMI info...
gp_bmi2 <- gp_bmi2 %>% mutate(type="gp_bmi")
df_bmi1$type <- "ukb_bmi"

bmi_df <- gp_bmi2 %>% bind_rows(df_bmi1)
bmi_df <- bmi_df %>% distinct()
### Save merged BMI dataset
saveRDS(bmi_df, file=paste(hba1c_dat_path, "data/t2d_bmi.rds", sep=""))
rm(gp_bmi2, df_bmi1, df_bmi)
### Make merged dataset a data.table with keys eids and event_dt
bmi_dt <- data.table(bmi_df, key= c("eid", "event_dt"))

############################################################################
### 5. Extract weight from GP data
############################################################################
gp_clinical <- bio_record(project_dir, "gp_clinical")
head(gp_clinical)
t2d_wrk <- readRDS(paste(hba1c_dat_path, "data/hba1c_t2d_mdd_medication_clean_over18.rds", sep=""))
eids_t2d <- unique(t2d_wrk$eid)

### Search read3 column for main code
### Main code for wgt observation == 22A.. (read v2 and v3):
### Search v3
gp_wgt_read3_22a <- gp_clinical %>% 
    filter(str_detect(read_3, "22A..")) %>%
    collect()
dim(gp_wgt_read3_22a) # 1160956       8

gp_wgt_read3_22a <- gp_wgt_read3_22a[gp_wgt_read3_22a$eid %in% eids_t2d, ]
dim(gp_wgt_read3_22a) ### [1] 237540      8

### Search read2 column for main code
gp_wgt_read2_22a <- gp_clinical %>%
    filter(str_detect(read_2, "22A..")) %>%
    collect()
dim(gp_wgt_read2_22a) # 492581       8

### Restrict to T2D EIDs.
gp_wgt_read2_22a <- gp_wgt_read2_22a[gp_wgt_read2_22a$eid %in% eids_t2d, ]
dim(gp_wgt_read2_22a)

### Additional weigt codes:
### Check for codes 1622, X76CG, XM01G, XE1h4, Xa7wI
### Check read2
gp_wgt_read2_1622 <- gp_clinical %>%
    filter(str_detect(read_2, "X76CG")) %>%
    collect()
dim(gp_wgt_read2_1622) # [1] 701   8
gp_wgt_read2_X76CG <- gp_clinical %>% 
    filter(str_detect(read_2, "X76CG")) %>%
    collect()
dim(gp_wgt_read2_X76CG) # 0 8
gp_wgt_read2_XM01G <- gp_clinical %>% 
    filter(str_detect(read_2, "XM01G")) %>%
    collect()
dim(gp_wgt_read2_XM01G) # 0 8 
gp_wgt_read2_XE1h4 <- gp_clinical %>% 
    filter(str_detect(read_2, "XE1h4")) %>%
    collect()
dim(gp_wgt_read2_XE1h4) # 0 8
gp_wgt_read2_Xa7wI <- gp_clinical %>% 
    filter(str_detect(read_2, "Xa7wI")) %>%
    collect()
dim(gp_wgt_read2_Xa7wI) # 0 8
rm(gp_wgt_read2_X76CG, gp_wgt_read2_Xa7wI, gp_wgt_read2_XE1h4, gp_wgt_read2_XM01G)
### Check read3
gp_wgt_read3_1662 <- gp_clinical %>%
    filter(str_detect(read_3, "1662")) %>%
    collect()
dim(gp_wgt_read3_1662) # [1] 48  8
gp_wgt_read3_X76CG <- gp_clinical %>% 
    filter(str_detect(read_3, "X76CG")) %>%
    collect()
dim(gp_wgt_read3_X76CG) # [1] 1 8
gp_wgt_read3_XM01G <- gp_clinical %>% 
    filter(str_detect(read_3, "XM01G")) %>%
    collect()
dim(gp_wgt_read3_XM01G) # [1] 1874    8
gp_wgt_read3_XE1h4 <- gp_clinical %>% 
    filter(str_detect(read_3, "XE1h4")) %>%
    collect()
dim(gp_wgt_read3_XE1h4) # [1] 290   8
gp_wgt_read3_Xa7wI <- gp_clinical %>% 
    filter(str_detect(read_3, "Xa7wI")) %>%
    collect()
dim(gp_wgt_read3_Xa7wI) # [1] 20  8
### Bind extracted data from additonal weight codes together
gp_wgt_other <- rbind(gp_wgt_read2_1622, gp_wgt_read3_1622, gp_wgt_read3_X76CG,
gp_wgt_read3_XM01G, gp_wgt_read3_XE1h4, gp_wgt_read3_Xa7wI)

### Little glance at them if you want...
#gp_wgt_read2_22a
#gp_wgt_read3_22a
#gp_wgt_other

### Initial data cleaning for weight data:
### Start with read 2...
range(gp_wgt_read2_22a$value1) # ""       "998"
table(gp_wgt_read2_22a$read_2)

table(gp_wgt_read2_22a$value1)[1:10] ### has "^" included- set to ""
table(gp_wgt_read2_22a$value2) ### Value 2 can be ignored in this round... But check in future extractions.
table(gp_wgt_read2_22a$value3) ### "K/M2"    "kg/m2"  "Unknown" set to "" Think ignore for now as some look like BMI measurements rather than weight.
### Implement the above
gp_wgt_read2_22a$value1[gp_wgt_read2_22a$value1 == "^"] <- ""
gp_wgt_read2_22a_1 <- gp_wgt_read2_22a %>% 
    mutate(miss1 = as.numeric(value1 == ""))

### Ignore value 2 and 3- so interested in miss1
table(gp_wgt_read2_22a_1$miss1)
gp_wgt_read2_22a_1 <- gp_wgt_read2_22a_1[gp_wgt_read2_22a_1$miss1 != 1, ]

### Make weight numeric
gp_wgt_read2_22a_1 %>% mutate(wgt = as.numeric(value1)) -> gp_wgt_read2_22a_1

gp_wgt_read2_22a_1 %>% select(eid, data_provider, event_dt, wgt) -> gp_wgt_read2_22a_1

### What boundaries to set to determine 'valid' weight?
### https://www.cell.com/ajhg/fulltext/S0002-9297(22)00049-0#supplementaryMaterial
### only includes those between 30 and 200kg

gp_wgt_read2_22a_1 <- gp_wgt_read2_22a_1[gp_wgt_read2_22a_1$wgt >= 30, ]
gp_wgt_read2_22a_1 <- gp_wgt_read2_22a_1[gp_wgt_read2_22a_1$wgt <=200, ]

### Check for valid event_dt:
any(is.na(gp_wgt_read2_22a_1$event_dt)) # F
any(gp_wgt_read2_22a_1$event_dt == "") # T
sum(gp_wgt_read2_22a_1$event_dt == "") # 10
gp_wgt_read2_22a_1 <- gp_wgt_read2_22a_1[gp_wgt_read2_22a_1$event_dt != "", ]

### move onto read3...
### Explore possible content of columns
table(gp_wgt_read3_22a$value1)[1:100] ### all "number" or ""
table(gp_wgt_read3_22a$value2) ### Value 2 all ""... But check in future extractions.
table(gp_wgt_read3_22a$value3) ### Value 3 all ""... But check in future extractions.

### Restrict to those without missing in value1, make weight numberic
gp_wgt_read3_22a_1 <- gp_wgt_read3_22a[gp_wgt_read3_22a$value1 != "", ]
gp_wgt_read3_22a_1 %>% mutate(wgt = as.numeric(value1)) -> gp_wgt_read3_22a_1
sum(is.na(gp_wgt_read3_22a_1$wgt)) # 0
### Restrict to 'valid' weights
gp_wgt_read3_22a_1 <- gp_wgt_read3_22a_1[gp_wgt_read3_22a_1$wgt >= 30, ]
gp_wgt_read3_22a_1 <- gp_wgt_read3_22a_1[gp_wgt_read3_22a_1$wgt <=200, ]
### Select required columns
gp_wgt_read3_22a_1 %>% select(eid, data_provider, event_dt, wgt) -> gp_wgt_read3_22a_1

### Check dates:
any(is.na(gp_wgt_read3_22a_1$event_dt)) # F
any(gp_wgt_read3_22a_1$event_dt == "") # T
sum(gp_wgt_read3_22a_1$event_dt == "") # 122
gp_wgt_read3_22a_1 <- gp_wgt_read3_22a_1[gp_wgt_read3_22a_1$event_dt != "", ]

### Combine read2 and read3 weight observations...
gp_wgt <- gp_wgt_read3_22a_1 %>% bind_rows(gp_wgt_read2_22a_1)
### Make sure date is in a 'date' format using as.Date
### Remove repeat rows
gp_wgt %>% mutate(event_dt = as.Date(event_dt, "%d/%m/%Y")) %>%
    distinct() -> gp_wgt1

### Remove any rows with incorrect dates
gp_wgt1 <- gp_wgt1[gp_wgt1$event_dt != "1902-02-02", ]
gp_wgt1 <- gp_wgt1[gp_wgt1$event_dt != "2037-07-07", ]
gp_wgt1 <- gp_wgt1[gp_wgt1$event_dt != "1903-03-03", ]

gp_wgt2 <- gp_wgt1 %>% select(-data_provider)
rm(gp_wgt1)

### Move onto weight data extracted with additional codes
### other
### Check the columns
table(gp_wgt_other$value1)[1:100] ### all "number" or ""
table(gp_wgt_other$value2) ### Value 2 all ""... But check in future extractions.
table(gp_wgt_other$value3) ### Value 3. Some numbers but ignore in favour of value1

### Restrict to rows with non-missing values for the value1 column
gp_wgt_other <- gp_wgt_other[gp_wgt_other$value1 != "", ]
### Make weight numeric
gp_wgt_other %>% mutate(wgt = as.numeric(value1)) -> gp_wgt_other
### Check for missing data
sum(is.na(gp_wgt_other$wgt)) # 0
### Remove 'invalid' weight measurements
gp_wgt_other <- gp_wgt_other[gp_wgt_other$wgt >= 30, ]
gp_wgt_other <- gp_wgt_other[gp_wgt_other$wgt <=200, ]
range(gp_wgt_other$wgt) ### 34.5 200

gp_wgt_other %>% select(eid, data_provider, event_dt, wgt) -> gp_wgt_other

### Check for missing event dates
any(is.na(gp_wgt_other$event_dt)) # F
any(gp_wgt_other$event_dt == "") # F

### Make the date into a date format R recognises using as.Date
### Remove repeat rows
gp_wgt_other %>% mutate(event_dt = as.Date(event_dt, "%d/%m/%Y")) %>%
    distinct() -> gp_wgt_other

### Check date range, and remove any with invalid dates
range(gp_wgt_other$event_dt) ### "1990-07-02" "2017-04-05"

gp_wgt_other <- gp_wgt_other %>% select(-data_provider)

### Merge with primary weight dataset, removing repeat rows
gp_wgt2 <- gp_wgt2 %>% bind_rows(gp_wgt_other)
gp_wgt2 <- gp_wgt2 %>% distinct()
gp_wgt3 <- gp_wgt2[gp_wgt2$eid %in% eids_t2d, ]
### Save data
saveRDS(gp_wgt3, file=paste(hba1c_dat_path, "data/GP_t2d_wgt.rds", sep=""))
rm(gp_wgt2)

############################################################################
### 6. Extract height from GP data
############################################################################
### Extract height: 229
gp_hgt_read3_229 <- gp_clinical %>% 
    filter(str_detect(read_3, "229")) %>%
    collect()
dim(gp_hgt_read3_229) ### [1] 737799      8
gp_hgt_read2_229 <- gp_clinical %>% 
    filter(str_detect(read_2, "229")) %>%
    collect()
dim(gp_hgt_read2_229) ### [1] 261435      8
### Restrict to T2D EIDs
gp_hgt_read2_229 <- gp_hgt_read2_229[gp_hgt_read2_229$eid %in% eids_t2d, ]
dim(gp_hgt_read2_229) ### [1] 47819     8
gp_hgt_read3_229 <- gp_hgt_read3_229[gp_hgt_read3_229$eid %in% eids_t2d, ]
dim(gp_hgt_read3_229) ### [1] 124225      8

### Explore contents of columns
table(gp_hgt_read2_229$value1)[1:100] ### has '^'
table(gp_hgt_read2_229$value2) ### Can ignore
table(gp_hgt_read2_229$value3) ### can ignore

table(gp_hgt_read3_229$value1)[1:100]
table(gp_hgt_read3_229$value2) ### all ""
table(gp_hgt_read3_229$value3) ### all ""

### For both read2 and read3 extractions focus on value1, remove missing rows
gp_hgt_read2_229 <- gp_hgt_read2_229[gp_hgt_read2_229$value1 != "^",  ]
gp_hgt_read2_229 <- gp_hgt_read2_229[gp_hgt_read2_229$value1 != "",  ]
gp_hgt_read3_229 <- gp_hgt_read3_229[gp_hgt_read3_229$value1 != "",  ]

### Bind rows of read2 and read3 extracted height data
### Restrict to columns of interest
gp_hgt <- gp_hgt_read2_229 %>% bind_rows(gp_hgt_read3_229) %>% select(-value2, -value3, -data_provider)
### Make hgt column numeric
gp_hgt <- gp_hgt %>% mutate(hgt = as.numeric(value1)) %>% select(-value1)
range(gp_hgt$hgt)
### Some hgt measurements could be in Ms, some in CM (possibly imperial measurements too?)
### Make a 'CM' column
gp_hgt <- gp_hgt %>% mutate(hgt2 = hgt*100)
### Make a rule prioritising original hgt column within a range of CM values.
### If hgt column is outside of this range check hgt2
### If hgt2 outside of range = NA
gp_hgt <- gp_hgt %>% mutate(hgt3 = ifelse(hgt <= 210 & hgt >= 125, hgt, ifelse(hgt2 <= 210 & hgt2 >= 125, hgt2, NA)))
gp_hgt <- gp_hgt[!is.na(gp_hgt$hgt3), ]
dim(gp_hgt) # [1] 166816      7
length(unique(gp_hgt$eid)) # 17463
### Tidy to just hgt3 and rename this as hgt
gp_hgt <- gp_hgt %>% select(-read_2, -read_3, -hgt, -hgt2) %>% rename(hgt = hgt3)
### Focus on dates
### Reformat
### Remove missing and incorrect dates
gp_hgt$event_dt <- as.Date(gp_hgt$event_dt, "%d/%m/%Y")
gp_hgt <- gp_hgt[!is.na(gp_hgt$event_dt), ]
range(gp_hgt$event_dt)
gp_hgt <- gp_hgt[gp_hgt$event_dt != "1902-02-02", ]
gp_hgt <- gp_hgt[gp_hgt$event_dt != "1903-03-03", ]
gp_hgt <- gp_hgt[gp_hgt$event_dt != "2037-07-07", ]
### Make dataset a data.table otdered by eid and event_dt and save
gp_hgt <- data.table(gp_hgt, key=c("eid", "event_dt"))
saveRDS(gp_hgt, file=paste(hba1c_dat_path, "data/GP_t2d_hgt.rds", sep=""))
### It is hard to know if height is correct with confidence so I will extract
### UK Biobank assessment height measuresments

############################################################################
### 7. Extract height from UK Biobank initial assessment
############################################################################
f %>%
    select(field, name) %>%
    filter(str_detect(field, "50-0.0"))

# standing_height_f50_0_0
### Save field to field file
f %>%
select(field, name) %>%
filter(str_detect(field, "50-0.0")) %>%
slice(1) %>%
bio_field_add(paste(analysis_dir, "data/ukb_hgt_field.txt", sep=""))
### Extract height data
bio_phen(
 project_dir,
 field = paste(analysis_dir, "data/ukb_hgt_field.txt", sep=""),
 out = paste(analysis_dir, "data/ukb_hgt", sep="")
)
### Read in extracted UKB height data
hgt_ukb <- readRDS(paste(analysis_dir, "data/ukb_hgt.rds", sep=""))
### Rename columns
hgt_ukb <- bio_rename(hgt_ukb, f)
colnames(hgt_ukb)
### Make height into Ms
hgt_ukb <- hgt_ukb %>% mutate(hgt = standing_height_f50_0_0/100) %>% select(-standing_height_f50_0_0)
sum(is.na(hgt_ukb$hgt)) ### 2583

############################################################################
### 8. Create additional BMI measurements with weight (GP) and height (UKB)
############################################################################
### Add height column into weight data set so that an individuals height is repeated for each row with their EID
bmi_ds2 <- gp_wgt3 %>% left_join(hgt_ukb)
### Calculate BMI
bmi_ds2 <- bmi_ds2 %>% mutate(bmi = wgt/(hgt^2))
### Remove height and weight and save BMI data
bmi_ds2 <- bmi_ds2 %>% select(-wgt, -hgt)
bmi_ds2 <- data.table(bmi_ds2, key=c("eid", "event_dt"))

saveRDS(bmi_ds2, file=paste(hba1c_dat_path, "data/GP_t2d_BMIfrom_wgt.rds", sep=""))

############################################################################
### 9. Merge BMI datasets (remove duplicates)
############################################################################
bmi_dt <- readRDS(paste(hba1c_dat_path, "data/t2d_bmi.rds", sep=""))
bmi_dt2 <- bmi_dt %>% bind_rows(bmi_ds2) %>% distinct()
dim(bmi_dt2) # [1] 650707      4
length(unique(bmi_dt2$eid)) ### 17542
bmi_dt2$type[is.na(bmi_dt2$type)] <- "gp_bmi_wgt"

saveRDS(bmi_dt2, file=paste(hba1c_dat_path, "data/t2d_bmi_bmiwgt.rds", sep=""))

############################################################################
### 10. Create basic HbA1c dataset for LDA
############################################################################
### How many T2D individuals have a BMI measurement 'close' to their baseline?
### Need to restructure data for LDA to find this...
### Read in HbA1c dataset:
df.s <- readRDS(paste(hba1c_dat_path, "data/hba1c_t2d_mdd_medication_clean_over18.rds", sep=""))
### Make data.table
dt.s <- data.table(df.s, key=c("eid", "event_dt"))
length(unique(dt.s$eid)) ### 17544

### Create time variable for analysis
### Just want to check that everyone has a unique T2D diagnosis date:
### For each EID extract their T2D diagnosis date
summ6 <- dt.s[, (unique(as.Date(first_t2d_date), na.rm=T)), by=eid]
sum(is.na(summ6$V1)) ### 0
### Create a data.table counting the number of rows for each EID in dt.s
count6 <- dt.s[, .(rowCount = .N), by=eid]

### Decided to rename T2D diagnosis date variable
df.s$t2d_diag_date <- df.s$first_t2d_date
dt.s <- data.table(df.s, key=c("eid", "event_dt"))
df.s <- data.frame(dt.s)

### Create a time_base variable measuring the time between an observation and T2D diagnosis
### Negative dates -> observation occurred prior to T2D diagnosis
df.s$time_base <- time_length(interval(as.Date(df.s$t2d_diag_date), as.Date(df.s$event_dt)),"years")
range(df.s$time_base) ### -19.03005  56.16667
dt.s <- data.table(df.s, key=c("eid", "event_dt"))
df.s <- data.frame(dt.s)
### How many people have a baseline measurement?
sum(as.numeric(df.s$time_base == 0)) ## 8575

### There are some duplicates still in data atm...
### I.e. multiple HbA1c measurements on same date
### Restrict to fewer key columns
dt.s1 <- dt.s %>% select(eid, event_dt, hba1c, dob, type, depression, metformin:time_base)
dim(dt.s1) # 355216     32
length(unique(dt.s1$eid)) # 17544
### remove duplicates part 1:
dt.s1 <- dt.s1%>% distinct()
dim(dt.s1) # 303132     32
length(unique(dt.s1$eid)) # 17544
dt.s <- dt.s1
df.s <- data.frame(dt.s)
### Check for remaining duplicates
sum(as.numeric(df.s$time_base == 0)) ### 7727 
### Check for duplicated baselines...
### restrict to baseline measurements
test <- df.s[df.s$time_base == 0,]
test.dt <- data.table(test)
### Count how many rows each ID has (should be 1)
dup <- test.dt[, .(rowCount = .N), by=eid]
table(dup$rowCount)
dup[dup$rowCount == 2, ] 
dim(dup[dup$rowCount == 2, ] ) ### 107 2
### There are still date-duplicates (same date, different HbA1c values)
### order by eid and then time_base
dt.s <- data.table(df.s, key=c("eid", "time_base")) 
df.s <- data.frame(dt.s)
### NB: There are probably other duplicates too at t != 0 too...
### Assume if there are multiple on same day, then can take the mean...
### Prioritise hba1c from gp (type == hba1c) rather than UKB (type == biomarker)
eids_hba1c <- unique(dt.s$eid)
### Loop checks the following
### For each EID restricts to their HbA1c data and counts the number of times
### a time_base measurement is observed.
### If dim of dataset equals length of unqiue time_base -> no duplicates and
### we output the original hba1c data for this EID. Else
### Extract unique dates for EID and for each date
### Restrict to hba1c data for that date. If there is one observation that is
### added to the output
### If there are more we restrict to GP measurements (if there are any) and take the mean OR
### if there are no GP HbA1c measurements we take the mean of the UKB measurements

hba1c_dat <- dt.s
test <- NULL
for(i in 1:length(eids_hba1c)){
    x <- eids_hba1c[i]
    dati <- hba1c_dat[hba1c_dat$eid == x, ]
    dup_check_i <- dati[, length(unique(time_base))]
    if(dim(dati)[1] == dup_check_i){
        out <- data.table(dati)
    }else{
        out <- NULL
        datesi <- unique(dati$event_dt)
        for(j in 1:length(datesi)){
            X <- datesi[j]
            dat_datei <- dati[dati$event_dt == X, ]
            if(dim(dat_datei)[1] == 1){
                out_dat <- dat_datei
            }else{
                if(any(dat_datei$type == "hba1c")){
                    dat_datei <- dat_datei[dat_datei$type == "hba1c", ]
                    out_dat <- dat_datei[1,]
                    out_dat$hba1c <- mean(dat_datei[, hba1c])
                }else{
                    out_dat <- dat_datei[1,]
                    out_dat$hba1c <- mean(dat_datei[, hba1c])
                }      
            }
            out <- rbind(out, out_dat)
            rm(out_dat)
        }
        rm(datesi)
        rm(X)
    }
    out$event_dt <- as.Date(out$event_dt)
    out$first_t2d_date <- as.Date(out$first_t2d_date)
    out$dob <- as.Date(out$dob)
    out$dep_first_date <- as.Date(out$dep_first_date)
    out$t2d_diag_date <- as.Date(out$t2d_diag_date)
    test <- rbind(test, out)
    rm(out)
    }
length(unique(test$eid)) ## 17544
dim(test) ### 297141     32

dt.s <- data.table(test, key=c("eid", "event_dt"))
df.s <- data.frame(dt.s)
sum(as.numeric(df.s$time_base == 0)) ### has dropped 7620

### save interim duplicate removed file
saveRDS(dt.s,  file= paste(analysis_dir, "data/t2d_hba1c_o18_dups_removed.rds", sep=""))

### For individuals with no baseline obs add in a t = 0 row with NA
### populate this additional row with other info

### Initial idea is to identify the row closest to t=0...
### Which IDs have a HbA1c observation at T2D diagnosis
ids_with_t0 <- dt.s$eid[dt.s$time_base == 0]
### hba1c dataset for those WITHOUT baseline (not in the ID list)
base0addin <- dt.s[!(eid %in% ids_with_t0), ]
length(unique(base0addin$eid)) ### 9924
### Find the minimum absolute time_base value (i.e. values closest to 0)
abs_min <- base0addin[ , min(abs(time_base)), by=eid]
dim(abs_min) ### 9924    2

### Count the rows for each EID
count.tab <- base0addin[, .(rowCount = .N), by=eid]
### Create a column with the abs_min time_base value:
base0addin$abs_min_t <- unlist(mapply(rep, x= abs_min$V1, times= count.tab$rowCount))
### Try to restrict to rows where time_base == abs_min_t
test <- base0addin[abs(base0addin$time_base) == base0addin$abs_min_t, ]
### Check row counts in this new dataset by EID
count.tab <- test[, .(rowCount = .N), by=eid]
count.tab[rowCount == 2, ] ### 2 IDs where there is a +/- issue
### That is, there are 2 IDs where there is a +abs_min_value and a -abs_min_value
### Identify these EIDs:
id_list <- count.tab[rowCount == 2, eid]
### Remove the osbervation occuring AFTER baseline (i.e. priorise the negative)
rm_bin <- rep(0, dim(test)[1])
rm_bin[(test$eid %in% id_list) & test$time_base > 0] <- 1
test <- test[rm_bin == 0, ]
count.tab <- test[, .(rowCount = .N), by=eid]
table(count.tab$rowCount)
#   1 
#9924
### all 1 so no double ups now :-)
### Rename this dataset...
base0addin <- test
rm(test, id_list, rm_bin, count.tab)
### Update values in this base
base0addin <- base0addin %>% select(-abs_min_t)
base0addin$type <- "missing_baseline"
base0addin$hba1c <- NA
base0addin$event_dt <- base0addin$first_t2d_date
base0addin$time_base <- 0
base0addin$metformin <- NA
base0addin$sulfonylurea <- NA
base0addin$TZD <- NA
base0addin$acarbose <- NA
base0addin$DPP4i <- NA
base0addin$GLP1Ra <- NA
base0addin$insulin <- NA
base0addin$meglitinide <- NA
base0addin$SGLT2i <- NA
base0addin$sum_meds <- NA
base0addin$medication_coded <- NA
### Update medication codes for these individuals
### Recall, the following are the functions required for medication info to
### populate dataset:
med_f_internal <- function(x, med_x){
  out <- rep(0, 9)
  x <- as.Date(x)
  xb <- x - days(100) 
  med_in_j <- as.numeric(med_x$event_dt <= x)
  med_in_j[med_x$event_dt < xb] <- 0
  if(sum(med_in_j) > 0){
    meds_i <- med_x[med_in_j == 1, ]
    out[1] <- as.numeric((as.numeric(any(meds_i$drug == "metformin")) + as.numeric(any(meds_i$drug == "metformin_DPP4i")) + as.numeric(any(meds_i$drug == "metformin_TZD"))) > 0)
    out[2] <- as.numeric((any(meds_i$drug == "sulfonylurea")))
    out[3] <- as.numeric((as.numeric(any(meds_i$drug == "TZD")) + as.numeric(any(meds_i$drug == "metformin_TZD"))) > 0)
    out[4] <- as.numeric((any(meds_i$drug == "meglitinide")))
    out[5] <- as.numeric((any(meds_i$drug == "acarbose")))
    out[6] <- as.numeric((as.numeric(any(meds_i$drug == "DPP4i")) + as.numeric(any(meds_i$drug == "metformin_DPP4i"))) > 0)
    out[7] <- as.numeric((as.numeric(any(meds_i$drug == "GLP1Ra")) + as.numeric(any(meds_i$drug == "GLP1Ra_insulin"))) > 0)
    out[8] <- as.numeric((any(meds_i$drug == "SGLT2i")))
    out[9] <- as.numeric((as.numeric(any(meds_i$drug == "insulin")) + as.numeric(any(meds_i$drug == "GLP1Ra_insulin"))) > 0)
  }
  out <- matrix(out)
}
med_f <- function(x, hba1c_data = hba1c_v3dt, med_data = med_dt){
  hba1c_x <- hba1c_data[hba1c_data$eid == x, ]
  med_x <- med_data[med_data$eid == x, ]
  if(dim(med_x)[1] > 0){
    tmp_out <- t(mapply(med_f_internal,x = hba1c_x$event_dt, MoreArgs = list(med_x=med_x)))
  }else{
    tmp_out <- matrix(NA, nrow=dim(hba1c_x)[1], ncol=9)
  }
  (tmp_out)
}
### Read in the medication dataset:
drug_data_hba1csubset_v2 <- fread(paste(hba1c_dat_path, "data/gp_diabetes_medication_data.txt", sep=""))
### Make sure date formatting is ok
drug_data_hba1csubset_v2$event_dt <- as.Date(drug_data_hba1csubset_v2$event_dt)
drug_data_hba1csubset_v2$prev_drug_date <- as.Date(drug_data_hba1csubset_v2$prev_drug_date)
drug_data_hba1csubset_v2$med_start <- as.Date(drug_data_hba1csubset_v2$med_start)
drug_data_hba1csubset_v2$med_end <- as.Date(drug_data_hba1csubset_v2$med_end)
drug_data_hba1csubset_v2 <- drug_data_hba1csubset_v2[(drug_data_hba1csubset_v2$eid %in% unique(base0addin$eid)), ]
### Make a data.table version to use in functions
med_dt <- data.table(drug_data_hba1csubset_v2, key=c("eid", "event_dt"))
rm(drug_data_hba1csubset_v2)
### Apply function:
meds_out <- mapply(x=unique(base0addin$eid), med_f, MoreArgs = list(hba1c_data = base0addin, med_data = med_dt),SIMPLIFY = TRUE)
meds_out <- t(meds_out)
dim(meds_out)
### Add medication info into dataset:
base0addin$metformin <- meds_out[,1]
base0addin$sulfonylurea <- meds_out[,2]
base0addin$TZD <- meds_out[,3]
base0addin$meglitinide <- meds_out[,4]
base0addin$acarbose <- meds_out[,5]
base0addin$DPP4i <- meds_out[,6]
base0addin$GLP1Ra <- meds_out[,7]
base0addin$SGLT2i <- meds_out[,8]
base0addin$insulin <- meds_out[,9]
rm(med_dt)
### Sum the meds an individual is on at baseline
base0addin[, sum_meds:= (metformin + sulfonylurea + TZD + meglitinide + acarbose + DPP4i + GLP1Ra + SGLT2i + insulin)]
base0addin[, medication_coded:= sum_meds]
### Create the four level medication variable
base0addin$medication_coded[base0addin$insulin == 1] <- 3
base0addin$medication_coded[base0addin$sum_meds >= 3] <- 3
table(base0addin$sum_meds) ### 10 people in this subgroup on 2 meds at baseline -> these individuals need to be removed later.

### join with orig dataset:
dt.s <- dt.s %>% bind_rows(base0addin)
### Then make data.table again with key EID and event_dt to re-order the data table
dt.s <- data.table(dt.s, key=c("eid", "event_dt"))
df.s <- data.frame(dt.s)

### For the depression variable- is it prevalent or time-varying?
check.dep <- dt.s[, unique(depression), by=eid]
dim(check.dep)[1] ### 17544 -> not time varying- is ever-never by end of study period...
### Create a time-varying version...
dt.s$event_dt <- as.Date(dt.s$event_dt, "%Y-%m-%d")
dt.s$first_t2d_date <- as.Date(dt.s$first_t2d_date, "%Y-%m-%d")
dt.s$dob <- as.Date(dt.s$dob, "%Y-%m-%d")
dt.s$dep_first_date <- as.Date(dt.s$dep_first_date, "%Y-%m-%d")
dt.s <- dt.s %>% mutate(dep_tv = ifelse(is.na(dep_first_date), NA, as.numeric(dep_first_date <= event_dt)))
check.dep <- dt.s[, unique(dep_tv), by=eid]
dim(check.dep)[1] ### 18121 >> 17544 - time varying now
rm(check.dep)

### save intermediate file
saveRDS(dt.s, file=paste(analysis_dir, "data/t2d_hba1c_traj_analysis_prelim.rds", sep=""))

############################################################################
### 11. Apply pre-imputation exclusion critera
############################################################################
### Check some exlusion criteria
### Prescriptions of insulin within 1 year
### Prescriptions at baseline...
### Age at T2D (< 18 removed)- done above
### at least one measurement >= 39, one of which must occur within a 6 month window of diagnosis
############################################################################
length(unique(dt.s$eid)) #  17544
### Prescriptions at baseline check
### Since adding in a baseline row for those missing HbA1c at baseline we identify 12 more individuals with multi-medications at baseline
dt0 <- dt.s[time_base ==0, ]
rm_bin <- as.numeric(dt0$sum_meds > 1)
sum(rm_bin, na.rm=T) ### 12
rm_bin[is.na(rm_bin)] <- 0
rm_eids <- dt0$eid[rm_bin == 1]

dt.s <- dt.s[!(dt.s$eid %in% rm_eids), ]
df.s <- data.frame(dt.s)
length(unique(dt.s$eid)) ### 17532

### at least one measurement >= 39, one of which must occur within a 6 month window of diagnosis
max_obs <- dt.s[, max(hba1c, na.rm=T), by=eid]
range(max_obs$V1) ### 28 to 187
rm_eids <- max_obs$eid[max_obs$V1 < 39]
length(rm_eids) ### 53
dt.s <- dt.s[!(dt.s$eid %in% rm_eids), ]
df.s <- data.frame(dt.s)
length(unique(dt.s$eid)) ### 17479
# one of which in 6 month window of diagnosis...
wrk <- dt.s
wrk[, min_int := t2d_diag_date %m-% months(6)]
wrk[, max_int := t2d_diag_date %m+% months(6)]
bin1 <- as.numeric(wrk$event_dt <= wrk$max_int)
bin1[wrk$event_dt < wrk$min_int] <- 0
wrk <- wrk[bin1 == 1, ]
length(unique(wrk$eid)) ### 17479
### need to remove those with no measurement within a 6 month window of diagnosis...
wrk <- wrk[!is.na(hba1c), ]
length(unique(wrk$eid)) ### 12942 have measurements within 6 month window of diagnosis...
rm_eid <- unique(dt.s$eid[!(dt.s$eid %in% wrk$eid)])
length(rm_eid) ### 4537

### Remove those with no measurement > 38 within 6 month window...
max_obs <- wrk[, max(hba1c, na.rm=T), by=eid]
range(max_obs$V1) # 15 184
sum(max_obs$V1 < 39) ### 255
rm_eids2 <- max_obs$eid[max_obs$V1 < 39]

rm_eid_combi <- unique(c(rm_eid, rm_eids2))
length(rm_eid_combi) ### 4792

dt.s <- dt.s[!(eid %in% rm_eid_combi), ]
length(unique(dt.s$eid)) # 12687
df.s <- data.frame(dt.s)

### Prescriptions of insulin within 1 year
#dt.s <- dt.s[time_base >= 0 & time_base <= 1, ]
#med_counts <- dt.s[, max(insulin, na.rm=T), by=eid]
#table(med_counts$V1)
#any(med_counts$insulin == 1) # FALSE - all insulin 1 year now removed

length(unique(df.s$eid)) ### 12687
max(df.s$time_base) ### 26.05191
### remove measurements after 10 years:
df.s <- df.s[df.s$time_base <= 10, ]
dim(df.s) ### 176980     35
length(unique(df.s$eid)) ### 12687
dt.s <- data.table(df.s, key=c("eid", "event_dt"))
count.tab <- dt.s[ , .(rowCount = .N), by=eid]
table(count.tab$rowCount)
saveRDS(dt.s, file=paste(analysis_dir, "data/MDDhba1c_imputeprep.txt", sep=""))

############################################################################
### 12. Create HbA1c observation count variables
############################################################################
### Create variable with the number of observations occuring before baseline
### Number occuring after baseline (up to 10 years)
### Number of visits in year prior to diagnosis...
###########################################################################
summ.pre <- dt.s[, sum(as.numeric(time_base < 0)), by=eid]
table(summ.pre$V1)
length(summ.pre$eid) ### 12687
summ.post <- dt.s[, sum(as.numeric(time_base > 0)), by=eid]
table(summ.post$V1) ### there are 228 individuals with 0 measurements after diag, and 532 with only 1. Keep for imputation though
length(summ.post$eid) ### 12687
count.tab <- dt.s[ , .(rowCount = .N), by=eid]
df.s <- data.frame(dt.s)
obs_pre <- unlist(mapply(rep, x= summ.pre$V1, times= count.tab$rowCount))
obs_post <- unlist(mapply(rep, x= summ.post$V1, times= count.tab$rowCount))
### total number of observations overall, including if baseline is measured.
total_obs <- unlist(mapply(rep, x= count.tab$rowCount, times= count.tab$rowCount))
dt.s$obs_pre <- obs_pre
dt.s$obs_post <- obs_post
dt.s$total_obs <- total_obs
df.s <- data.frame(dt.s)
summ1 <- dt.s[ , sum(as.numeric(time_base < 0 & time_base >= -1)), by=eid]
table(summ1$V1)
obs_pre1yr <- unlist(mapply(rep, x= summ1$V1, times= count.tab$rowCount))
dt.s$obs_pre1yr <- obs_pre1yr
df.s <- data.frame(dt.s)
saveRDS(dt.s, file=paste(analysis_dir, "data/MDDhba1c_imputeprep.txt", sep=""))

############################################################################
### 13. Create BMI at T2D diagnosis (baseline) covariate
############################################################################
### Read in required data:
dt.s <- readRDS(file=paste(analysis_dir, "data/MDDhba1c_imputeprep.txt", sep=""))
bmi_dt <- readRDS(paste(hba1c_dat_path, "data/t2d_bmi_bmiwgt.rds", sep=""))
# is.Date(dt.s$event_dt) ### TRUE

diag_date <- dt.s[, unique(t2d_diag_date), by=eid]
### restrict BMI dataset to those in the current HbA1c trajectory dataset
bmi_dt <- bmi_dt[bmi_dt$eid %in% unique(dt.s$eid), ]
length(unique(bmi_dt$eid)) # 12685 (2 less than the hba1c dataset)

colnames(diag_date)[2] <- "t2d_diag_date"
### Double check BMI dates (as perhaps didn't with weight and height version?)
range(bmi_dt$event_dt) # "1902-02-02" "2022-02-03"
bmi_dt <- bmi_dt[bmi_dt$event_dt != "1902-02-02", ]

### Add T2D diagnosis date info to BMI dataset
bmi_dt <- bmi_dt %>% left_join(diag_date)
### Order BMI dataset by EID and event_dt
bmi_dt <- data.table(bmi_dt, key=c("eid", "event_dt"))
### Create a time_diff variable
bmi_dt[, time_diff:= time_length(interval(as.Date(t2d_diag_date), as.Date(event_dt)), "years")]
#closest_bmi <- bmi_dt[, min(abs(time_diff)), by=eid]
#range(closest_bmi$V1) # 0.00000 13.37534
### Ok event range for baseline BMI is [-5, 2/12]...
### Restrict to individuals with BMI measurement(s) available between
### 5 years before and 2 months after a T2D diagnosis
### Know this is a WIDE window.
bmi_dt2 <- bmi_dt[time_diff <= 2/12, ]
length(unique(bmi_dt2$eid)) ### 12376
range(bmi_dt2$time_diff) #  -35.8821918   0.1666667
bmi_dt2 <- bmi_dt2[time_diff >= -5, ]
length(unique(bmi_dt2$eid)) ### 12039
### Now need to select 1 obs per individual...
### The closest BMI observation to baseline
### Prioritise obsevation prior to T2D diagnosis if there is a absmin tie
closest_bmi <- bmi_dt2[, min(abs(time_diff)), by=eid]
colnames(closest_bmi)[2] <- "absmin"
bmi_dt2 <- bmi_dt2 %>% left_join(closest_bmi)
bmi_dt2 <- bmi_dt2 %>% mutate(keep = as.numeric(abs(time_diff) == absmin))
### Restrict to BMI data that equals the closest to BMI diagnosis date
bmi_dt2 <- bmi_dt2[keep==1, ]
bmi_dt2 <- data.table(bmi_dt2, key=c("eid", "event_dt"))
count_tab <- bmi_dt2[, .(rowCount = .N), by=eid]
### Can see there are ties...
table(count_tab$rowCount)
bmi_dt2 <- bmi_dt2 %>% select(-keep)
bmi_dt2 <- bmi_dt2 %>% mutate(negtime_bin = time_diff <=0, )
ids_bmi <- unique(bmi_dt2$eid)
### Write a function of what rows to keep for each EID in the BMI dataset
### If there is only 1 BMI observation - keep that. Else,
### If there are no BMI observations prior to T2D diagnosis then:
### Priorise UKB BMI records, then GP BMI records then GP weight derived BMI records
### Keep the SMALLEST BMI measurement available (conservative)
### If there are tied closest BMI measurements occuring prior to T2D restrict to these and again priorise data source and smallest observation as above:
to_keep_f <- function(x, bmi_dt=bmi_dt2){
    dati <- bmi_dt[bmi_dt$eid == x, ]
    if(dim(dati)[1] == 1){
        keep <- 1
    }else{
        time_bin <- as.numeric(dati$time_diff <= 0)
        ukb_bin <- as.numeric(dati$type == "ukb_bmi")
        gp_bin <- as.numeric(dati$type == "gp_bmi")
        gpwgt_bin <- as.numeric(dati$type == "gp_bmi_wgt")
        if(sum(time_bin) == 0){           
            if(any(ukb_bin == 1)){
                if(sum(ukb_bin) == 1){
                    keep = ukb_bin
                }else{
                    keep = ukb_bin
                    keep[dati$bmi != min(dati$bmi[keep==1])] <- 0
                }
            }else{
                if(any(gp_bin == 1)){
                    if(sum(gp_bin) == 1){
                        keep = gp_bin
                    }else{
                        keep = gp_bin
                        keep[dati$bmi != min(dati$bmi[keep==1])] <- 0
                    }
                }else{
                    if(sum(gpwgt_bin) == 1){
                        keep = gpwgt_bin
                    }else{
                        keep = gpwgt_bin
                        keep[dati$bmi != min(dati$bmi[keep==1])] <- 0
                    }
                }
            }
        }else{
            if(sum(time_bin) == 1){
                keep <- time_bin
            }else{
                ukb_bin2 <- as.numeric((time_bin + ukb_bin) > 1)
                gp_bin2 <- as.numeric((time_bin + gp_bin) > 1)
                gpwgt_bin2 <- as.numeric((time_bin + gpwgt_bin) > 1)
                if(any(ukb_bin2==1)){
                    if(sum(ukb_bin2)==1){
                        keep = ukb_bin2
                    }else{
                        keep = ukb_bin2
                        keep[dati$bmi != min(dati$bmi[keep==1])] <- 0
                    }
                }else{
                    if(any(gp_bin2 == 1)){
                        if(sum(gp_bin2) == 1){
                            keep = gp_bin2
                        }else{
                            keep = gp_bin2
                            keep[dati$bmi != min(dati$bmi[keep==1])] <- 0
                        }
                    }else{
                        if(sum(gpwgt_bin2) == 1){
                            keep = gpwgt_bin2
                        }else{
                            keep = gpwgt_bin2
                            keep[dati$bmi != min(dati$bmi[keep==1])] <- 0
                        }
                    }
                }
            }
        }
    }
    keep
}
### Apply horribly long function!
test <- mapply(to_keep_f,x = ids_bmi, MoreArgs = list(bmi_dt=bmi_dt2))
test1 <- do.call(c, test)
length(unique(bmi_dt2$eid))
sum(test1)

bmi_dt2$keep <- test1

length(unique(bmi_dt2$eid))
length(unique(bmi_dt2$eid[bmi_dt2$keep==1]))
### Restrict to those rows we are keeping
base_bmi <- bmi_dt2[bmi_dt2$keep == 1, ]
### Remove extra columns
base_bmi <- base_bmi %>% select(-keep, -absmin, -time_diff)
### save interim BMI baseline dataset
saveRDS(base_bmi, file=paste(analysis_dir, "data/baseline_bmi.txt", sep=""))

### join
dt.s <- readRDS(file=paste(analysis_dir, "data/MDDhba1c_imputeprep.txt", sep=""))
base_bmi <- readRDS(file=paste(analysis_dir, "data/baseline_bmi.txt", sep=""))
dt.s <- data.table(dt.s, key=c("eid", "event_dt"))
base_bmi <- data.table(base_bmi, key="eid")
base_bmi <- base_bmi %>% select(eid, bmi)

dt.s <- dt.s %>% left_join(base_bmi)

############################################################################
### 14. Create MDD status at baseline and medications at baseline
############################################################################
### Now have joined BMI, create depression_base, meds_base for later use...
############################################################################
### Update dep_tv slightly...
dt.s$dep_tv[dt.s$depression == 0] <- 0
### Restrict to baseline (t==0) 
dt0 <- dt.s[time_base == 0, ]
length(unique(dt.s$eid)) # [1] 12687
length(unique(dt0$eid))  # [1] 12687
dt0 <- dt0 %>% select(eid, metformin:medication_coded, dep_tv)
table(dt0$dep_tv)
#    0     1 
# 11474  1213 
table(dt0$medication_coded)
#    0    1 
# 8941  449 
colnames(dt0) <- c("eid", "metformin_base", "sulfonylurea_base", "TZD_base", "meglitinide_base", "acarbose_base", 
"DPP4i_base", "GLP1Ra_base", "SGLT2i_base", "insulin_base", "sum_meds_base", "medication_coded_base", "dep_base")
dt.s <- dt.s %>% left_join(dt0)

############################################################################
### 15. Add additional UK Biobank assessment variables to HbA1c dataset
############################################################################
### Join with additional UKB variables... At least for imputation...
############################################################################
### Read in covariate data
df <- readRDS(paste(analysis_dir, "data/covariate.rds", sep=""))
### Rename the columns
df <- bio_rename(df, f)
### Make a data.table version
dt_covs <- data.table(df, key="eid")
### Select covariates want to include
dt_covs <- dt_covs %>% select(eid, uk_biobank_assessment_centre_f54_0_0, townsend_deprivation_index_at_recruitment_f189_0_0, age_completed_full_time_education_f845_0_0,
age_completed_full_time_education_f845_1_0, age_completed_full_time_education_f845_2_0, cooked_vegetable_intake_f1289_0_0,
"salad_/_raw_vegetable_intake_f1299_0_0", fresh_fruit_intake_f1309_0_0, processed_meat_intake_f1349_0_0, alcohol_intake_frequency__f1558_0_0,
comparative_body_size_at_age_10_f1687_0_0, comparative_body_size_at_age_10_f1687_1_0, comparative_body_size_at_age_10_f1687_2_0,
qualifications_f6138_0_0, smoking_status_f20116_0_0, alcohol_drinker_status_f20117_0_0, ever_smoked_f20160_0_0, liking_for_vegetables_f20739_0_0,
ethnic_background_f21000_0_0, genotype_measurement_batch_f22000_0_0, genetic_ethnic_grouping_f22006_0_0, met_minutes_per_week_for_walking_f22037_0_0,
met_minutes_per_week_for_moderate_activity_f22038_0_0, met_minutes_per_week_for_vigorous_activity_f22039_0_0, summed_met_minutes_per_week_for_all_activity_f22040_0_0) %>%
rename(assessment_centre=uk_biobank_assessment_centre_f54_0_0, tdi=townsend_deprivation_index_at_recruitment_f189_0_0, 
age_ft_edu1=age_completed_full_time_education_f845_0_0, age_ft_edu2=age_completed_full_time_education_f845_1_0, 
age_ft_edu3=age_completed_full_time_education_f845_2_0, cooked_veg_intake=cooked_vegetable_intake_f1289_0_0, 
salad_raw_veg_intake="salad_/_raw_vegetable_intake_f1299_0_0", fresh_fruit_intake=fresh_fruit_intake_f1309_0_0,
processed_meat_intake=processed_meat_intake_f1349_0_0, alc_intake_freq=alcohol_intake_frequency__f1558_0_0,
comp_body_size1=comparative_body_size_at_age_10_f1687_0_0, comp_body_size2=comparative_body_size_at_age_10_f1687_1_0,
comp_body_size3=comparative_body_size_at_age_10_f1687_2_0, qualifications=qualifications_f6138_0_0, smk_status=smoking_status_f20116_0_0,
alc_drinker_status=alcohol_drinker_status_f20117_0_0, ever_smk=ever_smoked_f20160_0_0, like_veg=liking_for_vegetables_f20739_0_0,
sr_ethnicity=ethnic_background_f21000_0_0, batch=genotype_measurement_batch_f22000_0_0, genetic_ethnicity=genetic_ethnic_grouping_f22006_0_0,
met_walk=met_minutes_per_week_for_walking_f22037_0_0, met_mod=met_minutes_per_week_for_moderate_activity_f22038_0_0,
met_vig=met_minutes_per_week_for_vigorous_activity_f22039_0_0, met_all=summed_met_minutes_per_week_for_all_activity_f22040_0_0)
### age full time education: combine variables across time points
# prioritise initial assessment
dt_covs <- dt_covs %>%
    mutate(age_dt_edu = ifelse(!is.na(age_ft_edu1), age_ft_edu1, ifelse(!is.na(age_ft_edu2), age_ft_edu2, age_ft_edu3)))
dt_covs <- dt_covs %>% select(-age_ft_edu1, -age_ft_edu2, -age_ft_edu3)
### Similarly for comp_body_size_age10
dt_covs <- dt_covs %>%
    mutate(comp_body_size_age10 = ifelse(!is.na(comp_body_size1), comp_body_size1, ifelse(!is.na(comp_body_size2), comp_body_size2, comp_body_size3)))
dt_covs <- dt_covs %>% select(-comp_body_size1, -comp_body_size2, -comp_body_size3)
### Cooked veg intake, salad_raw_veg_intake, fresh_fruit_intake
### -10 == < 1, make 0.5, -3 == prefer not to answer == NA, -1 == don't know == NA
### Change coding for these variables
dt_covs$cooked_veg_intake[dt_covs$cooked_veg_intake == -10] <- 0.5
dt_covs$cooked_veg_intake[dt_covs$cooked_veg_intake == -3] <- NA
dt_covs$cooked_veg_intake[dt_covs$cooked_veg_intake == -1] <- NA

dt_covs$salad_raw_veg_intake[dt_covs$salad_raw_veg_intake == -10] <- 0.5
dt_covs$salad_raw_veg_intake[dt_covs$salad_raw_veg_intake == -3] <- NA
dt_covs$salad_raw_veg_intake[dt_covs$salad_raw_veg_intake == -1] <- NA

dt_covs$fresh_fruit_intake[dt_covs$fresh_fruit_intake == -10] <- 0.5
dt_covs$fresh_fruit_intake[dt_covs$fresh_fruit_intake == -3] <- NA
dt_covs$fresh_fruit_intake[dt_covs$fresh_fruit_intake == -1] <- NA
### processed_meat_intake. Create a factor variable.
dt_covs$processed_meat_intake_f <- NA
dt_covs$processed_meat_intake_f[dt_covs$processed_meat_intake ==0] <- "Never"
dt_covs$processed_meat_intake_f[dt_covs$processed_meat_intake ==1] <- "<1_p/week"
dt_covs$processed_meat_intake_f[dt_covs$processed_meat_intake ==2] <- "1_p/week"
dt_covs$processed_meat_intake_f[dt_covs$processed_meat_intake ==3] <- "2to4_p/week"
dt_covs$processed_meat_intake_f[dt_covs$processed_meat_intake ==4] <- "5to6_p/week"
dt_covs$processed_meat_intake_f[dt_covs$processed_meat_intake ==5] <- ">=1_daily"
dt_covs <- dt_covs %>% select(-processed_meat_intake)
### alc_intake_freq. Create a factor variable.
range(dt_covs$alc_intake_freq, na.rm=T)
dt_covs$alc_intake_freq_f <- NA
dt_covs$alc_intake_freq_f[dt_covs$alc_intake_freq == 6] <- "Never"
dt_covs$alc_intake_freq_f[dt_covs$alc_intake_freq == 5] <- "Special_occasions_only"
dt_covs$alc_intake_freq_f[dt_covs$alc_intake_freq == 4] <- "1to3_p/month"
dt_covs$alc_intake_freq_f[dt_covs$alc_intake_freq == 3] <- "1to2_p/week"
dt_covs$alc_intake_freq_f[dt_covs$alc_intake_freq == 2] <- "3to4_p/week"
dt_covs$alc_intake_freq_f[dt_covs$alc_intake_freq == 1] <- "daily_almost_daily"
table(dt_covs$alc_intake_freq_f)
dt_covs <- dt_covs %>% select(-alc_intake_freq)
### qualifications. Create a factor variable.
dt_covs$qualifications_f <- NA
dt_covs$qualifications_f[dt_covs$qualifications == 6] <- "Other_prof_qual"
dt_covs$qualifications_f[dt_covs$qualifications == 5] <- "NVQ/HND/HNC/equiv"
dt_covs$qualifications_f[dt_covs$qualifications == 4] <- "CSEs/equiv"
dt_covs$qualifications_f[dt_covs$qualifications == 3] <- "Olevels/GCSEs/equiv"
dt_covs$qualifications_f[dt_covs$qualifications == 2] <- "Alevel/As_levels/equiv"
dt_covs$qualifications_f[dt_covs$qualifications == 1] <- "College_uni"
dt_covs$qualifications_f[dt_covs$qualifications == -7] <- "None_of_the_above"
dt_covs <- dt_covs %>% select(-qualifications)
### smk_status. Create a factor variable.
dt_covs$smk_status_f <- NA
dt_covs$smk_status_f[dt_covs$smk_status == 0] = "Never"
dt_covs$smk_status_f[dt_covs$smk_status == 1] = "Previous"
dt_covs$smk_status_f[dt_covs$smk_status == 2] = "Current"
dt_covs <- dt_covs %>% select(-smk_status)
### alc_drinker_status. Create a factor variable.
dt_covs$alc_drinker_status_f <- NA
dt_covs$alc_drinker_status_f[dt_covs$alc_drinker_status == 0] = "Never"
dt_covs$alc_drinker_status_f[dt_covs$alc_drinker_status == 1] = "Previous"
dt_covs$alc_drinker_status_f[dt_covs$alc_drinker_status == 2] = "Current"
dt_covs <- dt_covs %>% select(-alc_drinker_status)
### ever_smk: yes/ no. Leave as 0, 1
dt_covs <- dt_covs %>% select(-genetic_ethnicity, -like_veg)

### merge with dt.s
dt.s <- dt.s %>% left_join(dt_covs)
### Add genetic super population info
# Use ukbkings funcitons to identify super population membership
genetic_ancestry <- bio_gen_ancestry(project_dir)
dt.s <- dt.s %>% left_join(genetic_ancestry)

### Re-code self report ethicity...
dt.s$sr_ethnicity_f <- NA
dt.s$sr_ethnicity_f[dt.s$sr_ethnicity == "6"] <- "Other"
dt.s$sr_ethnicity_f[dt.s$sr_ethnicity == "1"] <- "White"
dt.s$sr_ethnicity_f[dt.s$sr_ethnicity == "1001"] <- "White"
dt.s$sr_ethnicity_f[dt.s$sr_ethnicity == "1002"] <- "White"
dt.s$sr_ethnicity_f[dt.s$sr_ethnicity == "1003"] <- "White"
dt.s$sr_ethnicity_f[dt.s$sr_ethnicity == "2"] <- "Mixed_race"
dt.s$sr_ethnicity_f[dt.s$sr_ethnicity == "2001"] <- "Mixed_race"
dt.s$sr_ethnicity_f[dt.s$sr_ethnicity == "2002"] <- "Mixed_race"
dt.s$sr_ethnicity_f[dt.s$sr_ethnicity == "2003"] <- "Mixed_race"
dt.s$sr_ethnicity_f[dt.s$sr_ethnicity == "2004"] <- "Mixed_race"
dt.s$sr_ethnicity_f[dt.s$sr_ethnicity == "3"] <- "Asian"
dt.s$sr_ethnicity_f[dt.s$sr_ethnicity == "3001"] <- "Asian"
dt.s$sr_ethnicity_f[dt.s$sr_ethnicity == "3002"] <- "Asian"
dt.s$sr_ethnicity_f[dt.s$sr_ethnicity == "3003"] <- "Asian"
dt.s$sr_ethnicity_f[dt.s$sr_ethnicity == "3004"] <- "Asian"
dt.s$sr_ethnicity_f[dt.s$sr_ethnicity == "4"] <- "Black"
dt.s$sr_ethnicity_f[dt.s$sr_ethnicity == "4001"] <- "Black"
dt.s$sr_ethnicity_f[dt.s$sr_ethnicity == "4002"] <- "Black"
dt.s$sr_ethnicity_f[dt.s$sr_ethnicity == "4003"] <- "Black"
dt.s$sr_ethnicity_f[dt.s$sr_ethnicity == "5"] <- "Chinese"
### year of diagnosis...
dt.s$yo_diag <- year(dt.s$first_t2d_date)

############################################################################
### 16. Update coding for the medicine variables
############################################################################
### NA means not present in prescription database. 
### Now assume this means no medications.
dt.s$metformin[is.na(dt.s$metformin)] <- 0
dt.s$sulfonylurea[is.na(dt.s$sulfonylurea)] <- 0
dt.s$TZD[is.na(dt.s$TZD)] <- 0
dt.s$meglitinide[is.na(dt.s$meglitinide)] <- 0
dt.s$acarbose[is.na(dt.s$acarbose)] <- 0
dt.s$DPP4i[is.na(dt.s$DPP4i)] <- 0
dt.s$GLP1Ra[is.na(dt.s$GLP1Ra)] <- 0
dt.s$SGLT2i[is.na(dt.s$SGLT2i)] <- 0
dt.s$insulin[is.na(dt.s$insulin)] <- 0
dt.s$sum_meds[is.na(dt.s$sum_meds)] <- 0
dt.s$medication_coded[is.na(dt.s$medication_coded)] <- 0

dt.s$metformin_base[is.na(dt.s$metformin_base)] <- 0
dt.s$sulfonylurea_base[is.na(dt.s$sulfonylurea_base)] <- 0
dt.s$TZD_base[is.na(dt.s$TZD_base)] <- 0
dt.s$meglitinide_base[is.na(dt.s$meglitinide_base)] <- 0
dt.s$acarbose_base[is.na(dt.s$acarbose_base)] <- 0
dt.s$DPP4i_base[is.na(dt.s$DPP4i_base)] <- 0
dt.s$GLP1Ra_base[is.na(dt.s$GLP1Ra_base)] <- 0
dt.s$SGLT2i_base[is.na(dt.s$SGLT2i_base)] <- 0
dt.s$insulin_base[is.na(dt.s$insulin_base)] <- 0
dt.s$sum_meds_base[is.na(dt.s$sum_meds_base)] <- 0
dt.s$medication_coded_base[is.na(dt.s$medication_coded_base)] <- 0

############################################################################
### 17. Catch all section...
### MDD time-varying column (check, update),
### Remove those with only 1 HbA1c obs
### Create comparative body size at age 10 factor
### Create baseline HbA1c variable
############################################################################
### Using MDD time-varying for inputation. 
### Check have removed those with missing MDD dates
dep_dt <- dt.s[depression == 1, ]
any(is.na(dep_dt$dep_tv)) # FALSE- but this should be true due to depression cases with missing dates?
dt.s$dep_tv2 <- ifelse(dt.s$depression == 0, 0, ifelse(!is.na(dt.s$dep_first_date), dt.s$dep_tv, NA))
dt.s <- dt.s %>% select(-dep_tv) %>% rename(dep_tv = dep_tv2)
dt0 <- dt.s[time_base == 0, ]
table(dt0$depression)
#     0     1 
# 11232  1451 
### Note. Have already removed those with depression but no depression date

### Create a variable to do with whether the individual will stay in the LDA.
### Everyone Need 2 observations. 
dt.s <- dt.s[total_obs != 1, ] ### Even for imputation remove those with only 1 observation ever

dt.s <- data.table(dt.s, key=c("eid", "event_dt"))
dt.s$lda_keep <- as.numeric(dt.s$obs_post > 0)

saveRDS(dt.s, paste(hba1c_dat_path, "data/t2d_imputeprep_ds.rds", sep=""))
length(unique(dt.s$eid[dt.s$lda_keep == 1])) # 12459
dt0 <- dt.s[time_base == 0, ]
table(dt0$depression[dt0$lda_keep == 1])
###     0     1 
###  11034  1425 

### Create dep_post variable
dt.s <-  dt.s %>%
    mutate(dep_post = ifelse(depression == 0, 0, ifelse(dep_base == 0, 1, 0)))

dt0 <- dt.s[time_base == 0, ]
table(dt0$dep_base[dt0$lda_keep == 1])
#    0     1 
# 11471  1212 
table(dt0$dep_post[dt0$lda_keep == 1])
#    0     1 
# 12444   239 

### Check on source of info for dep == 0
#dep0 <- dt.s[depression == 0, ]
#dep0 <- dep0[type != "missing_baseline", ]
#check <- dep0[, any(type == "hba1c"), by="eid"]
#table(check$V1)


dep_post <- dt.s[dt.s$dep_post == 1, ]
dep_post <- dep_post[lda_keep == 1, ]

dep_post0 <- dep_post[time_base == 0, ]
table(dep_post0$pop)
# AFR AMR EUR SAS 
#  6   3 187  21 
sum(is.na(dep_post0$pop)) ## 22
table(dep_post0$sr_ethnicity_f)
sum(is.na(dep_post0$sr_ethnicity_f))

sum(is.na(dep_post0$tdi)) ### 0
sum(is.na(dep_post0$qualifications_f)) ### 5
sum(is.na(dep_post0$ever_smk)) ### 2

dep_post0$dep_time_diff <- time_length(interval(as.Date(dep_post0$first_t2d_date), as.Date(dep_post0$dep_first_date)),"years")

quantile(dep_post0$dep_time_diff, probs = seq(0, 1, 0.1))
dim(dep_post0)
sum(as.numeric(dep_post0$dep_time_diff <= 10))

check <- dt.s[, any(dep_tv[lda_keep==1] > 0), by="eid"]
table(check$V1)
check2 <- dt.s[, unique(depression), by="eid"]
dim(check2)[1] - length(unique(dt.s$eid)) ### 0 so the depression variable is correct...

### Create time difference variable for each individual between T2D and MDD diagnosis
dt.s <- dt.s %>% mutate(dep_time_diff = time_length(interval(as.Date(first_t2d_date), as.Date(dep_first_date)),"years"))
range(dt.s$dep_time_diff, na.rm=T)
range(dt.s$dep_time_diff[dt.s$dep_base == 1], na.rm=T)
range(dt.s$dep_time_diff[dt.s$dep_post == 1], na.rm=T)

### Create a binary variable indicating if an individual receives their initial diagnosis during the 10 years of study follow-up
dt.s <- dt.s %>% mutate(dep_post10years = ifelse(depression == 0, 0, ifelse(dep_time_diff > 0 & dep_time_diff <= 10, 1, 0)))
dt0 <- dt.s[time_base == 0, ]
dt0 <- dt0[lda_keep == 1, ]
table(dt0$dep_post10years, dt0$dep_post)
### Note: dep_post10years indicates whether individuals time of T2D diagnosis and time of depression diag is within a 10 year interval

### Create a missing baseline value indicator value
dt.s[, missing_baseline := as.numeric(any(is.na(hba1c))), by=eid]
### Create a baseline HbA1c column for each individual
dt.s[, hba1c_base := hba1c[time_base == 0], by=eid]

### Create a factor for comparative body size at age 10
dt.s$comp_body_size_age10_f <- NA
dt.s$comp_body_size_age10_f[dt.s$comp_body_size_age10 == 1] <- "Thinner"
dt.s$comp_body_size_age10_f[dt.s$comp_body_size_age10 == 2] <- "Plumper"
dt.s$comp_body_size_age10_f[dt.s$comp_body_size_age10 == 3] <- "Average"

### Save dataset
saveRDS(dt.s, paste(hba1c_dat_path, "data/t2d_imputeprep_ds.rds", sep=""))

############################################################################
### 18. Create mean SBP and DBP variables using UK Biobank assessments
############################################################################
### Add in BP measurements...
### Blood pressure check
### Read in LDA dataset created so far
dt.s <- readRDS(paste(hba1c_dat_path, "data/t2d_imputeprep_ds.rds", sep=""))
### Read in UK Biobank BP dataset
BP_df <- readRDS(paste(analysis_dir, "data/ukb_BP.rds", sep=""))
### rename dataset
BP_df <- bio_rename(BP_df, f)
### Split into SBP and DBP
sbp <- BP_df %>%
  select(systolic_blood_pressure_manual_reading_f93_0_0, systolic_blood_pressure_manual_reading_f93_0_1, 
  systolic_blood_pressure_automated_reading_f4080_0_0, systolic_blood_pressure_automated_reading_f4080_0_1)
dbp <- BP_df %>%
  select(diastolic_blood_pressure_manual_reading_f94_0_0, diastolic_blood_pressure_manual_reading_f94_0_1, 
  diastolic_blood_pressure_automated_reading_f4079_0_0, diastolic_blood_pressure_automated_reading_f4079_0_1)
### Find the mean SBP across measurements
msbp <- rowMeans(sbp, na.rm=T)
### Find the mean SBP across measurements
mdbp <- rowMeans(dbp, na.rm=T)
### Add these the main BP dataset
BP_df$msbp <- msbp
BP_df$mdbp <- mdbp
### Remove separate dataset
rm(sbp, dbp)
head(BP_df)
### Reduce to just EID and the mean BP measurements
BP_df1 <- BP_df %>% select(eid, msbp, mdbp)

dt.s <- data.table(dt.s, key=c("eid", "event_dt"))
### Add the means into the LDA dataset and save
dt.s <- dt.s %>% left_join(BP_df1)
head(dt.s)
saveRDS(dt.s, paste(hba1c_dat_path, "data/t2d_imputeprep_ds.rds", sep=""))

############################################################################
### 19. Create BMI measurment count columns (GP BMI only)
############################################################################
### Count of BMI measurements
gp_bmi2 <- readRDS(file=paste(hba1c_dat_path, "data/GP_t2d_bmi.rds", sep=""))
gp_bmi2 <- data.table(gp_bmi2, key=c("eid", "event_dt"))
row_count <- gp_bmi2[ ,.(rowCount = .N), by=eid]
table(row_count$rowCount)

dt0work <- dt0 %>% select(eid, t2d_diag_date)

gp_bmi2 <- gp_bmi2 %>% left_join(dt0work)
gp_bmi2 <- gp_bmi2 %>%
    mutate(bin_dates = as.numeric(event_dt <= t2d_diag_date))

tmp_bmi <- gp_bmi2[gp_bmi2$bin_dates == 1, ]
row_count <- tmp_bmi[ ,.(rowCount = .N), by=eid]
table(row_count$rowCount)
dim(row_count)
colnames(row_count)[2] <- "nobs_preBMI"
dt.s <- dt.s %>% left_join(row_count)
sum(is.na(dt.s$nobs_preBMI))
dt.s$nobs_preBMI[is.na(dt.s$nobs_preBMI)] <- 0

saveRDS(dt.s, paste(hba1c_dat_path, "data/t2d_imputeprep_ds.rds", sep=""))

############################################################################
### 20. Extract BP measurements from GP records
############################################################################
### BPs...
eids_t2d <- unique(dt.s$eid)
gp_clinical <- bio_record(project_dir, "gp_clinical")
head(gp_clinical)
### Main code for BP == 246..00 (read v2 and v3):
### Search v2
gp_bp_read2_246 <- gp_clinical %>%
    filter(str_detect(read_2, "246..")) %>%
    collect()
dim(gp_bp_read2_246) # 2001170       8
### Restrict to T2D cases
gp_bp_read2_246 <- gp_bp_read2_246[gp_bp_read2_246$eid %in% eids_t2d, ]
dim(gp_bp_read2_246) #  219825      8
### Search v3
gp_bp_read3_246 <- gp_clinical %>% 
    filter(str_detect(read_3, "246..")) %>%
    collect()
dim(gp_bp_read3_246) # 6071729       8
### Restrict to T2D cases
gp_bp_read3_246 <- gp_bp_read3_246[gp_bp_read3_246$eid %in% eids_t2d, ]
dim(gp_bp_read3_246) # 676523      8
head(gp_bp_read3_246)

### Look at read 2
gp_bp_read2_246$num1 <- as.numeric(gp_bp_read2_246$value1)
gp_bp_read2_246$num2 <- as.numeric(gp_bp_read2_246$value2)
gp_bp_read2_246$num3 <- as.numeric(gp_bp_read2_246$value3)
all_na <- as.numeric((as.numeric(is.na(gp_bp_read2_246$num1)) + as.numeric(is.na(gp_bp_read2_246$num2)) + as.numeric(is.na(gp_bp_read2_246$num3))) == 3)
gp_bp_read2_246 <- gp_bp_read2_246[all_na == 0, ]
### can ignore num3
all_oob <- as.numeric((as.numeric(gp_bp_read2_246$num1 <=0) + as.numeric(gp_bp_read2_246$num2 <= 0)) == 2)
gp_bp_read2_246 <- gp_bp_read2_246[all_oob == 0, ]

both_nums_avail <- as.numeric((as.numeric(!is.na(gp_bp_read2_246$num1)) + as.numeric(!is.na(gp_bp_read2_246$num2))) == 2)
table(both_nums_avail) ### says all 1 -> all measurements have 2 entries...
gp_bp_read2_246 <- gp_bp_read2_246 %>% select(eid, event_dt, read_2, read_3, num1, num2) %>% 
    rename(value1=num1, value2=num2)
range(gp_bp_read2_246$value1)
range(gp_bp_read2_246$value2)
gp_bp_read2_246 <- gp_bp_read2_246[gp_bp_read2_246$value1 < 290]
gp_bp_read2_246 <- gp_bp_read2_246[gp_bp_read2_246$value1 > 35]
gp_bp_read2_246 <- gp_bp_read2_246[gp_bp_read2_246$value2 < 300]
gp_bp_read2_246 <- gp_bp_read2_246[gp_bp_read2_246$value2 > 35]
range(gp_bp_read2_246$value1) ### 38 260
range(gp_bp_read2_246$value2) ### 37 230
gp_bp_read2_246$event_dt <- as.Date(gp_bp_read2_246$event_dt, "%d/%m/%Y")
gp_bp_read2_246 <- data.table(gp_bp_read2_246, key=c("eid", "event_dt"))
sbp_read2 <- apply(data.frame(cbind(gp_bp_read2_246$value1, gp_bp_read2_246$value2)),1, FUN=max)
dbp_read2 <- apply(data.frame(cbind(gp_bp_read2_246$value1, gp_bp_read2_246$value2)),1, FUN=min)
gp_bp_read2_246$sbp <- sbp_read2
gp_bp_read2_246$dbp <- dbp_read2
sbp_read2 <- gp_bp_read2_246 %>% select(eid, event_dt, sbp)
dbp_read2 <- gp_bp_read2_246 %>% select(eid, event_dt, dbp)

### Look at read 3
gp_bp_read3_246$event_dt <- as.Date(gp_bp_read3_246$event_dt, "%d/%m/%Y")
table(gp_bp_read3_246$value3) ### all blank
table(gp_bp_read3_246$value2) ### all blank
### only value1 is of interest
gp_bp_read3_246$value1 <- as.numeric(gp_bp_read3_246$value1)
gp_bp_read3_246 <- gp_bp_read3_246[!is.na(gp_bp_read3_246$value1)]
gp_bp_read3_246 <- gp_bp_read3_246 %>% 
    select(eid, event_dt, value1)
range(gp_bp_read3_246$value1)
gp_bp_read3_246 <- gp_bp_read3_246[gp_bp_read3_246$value1 > 35, ]
gp_bp_read3_246 <- gp_bp_read3_246[gp_bp_read3_246$value1 < 300, ]

gp_bp_read3_246 <- data.table(gp_bp_read3_246, key=c("eid", "event_dt"))

### Use the awesome dcast function
### SBP assumed to the maximum value
sbp_read3 <- dcast(gp_bp_read3_246, eid + event_dt ~ .,
                fun.aggregate = max, 
                value.var = 'value1')
### DBP assumed to the maximum value
dbp_read3 <- dcast(gp_bp_read3_246, eid + event_dt ~ .,
                fun.aggregate = min, 
                value.var = 'value1')
colnames(dbp_read3)[3] <- "dbp"
colnames(sbp_read3)[3] <- "sbp"

### Join read2 and read3 for SBP and DBP
sbp_gp <- sbp_read2 %>% bind_rows(sbp_read3)
dbp_gp <- dbp_read2 %>% bind_rows(dbp_read3)
### Remove repeat rows
sbp_gp <- sbp_gp %>% distinct()
dbp_gp <- dbp_gp %>% distinct()

### Remove those with no date info
dbp_gp <- dbp_gp[!is.na(dbp_gp$event_dt), ]
sbp_gp <- sbp_gp[!is.na(sbp_gp$event_dt), ]
### Remove invalid dates
dbp_gp <- dbp_gp[dbp_gp$event_dt != "1902-02-02", ]
dbp_gp <- dbp_gp[dbp_gp$event_dt != "2037-07-07", ]
dbp_gp <- dbp_gp[dbp_gp$event_dt != "1958-12-15", ]
dbp_gp <- dbp_gp[dbp_gp$event_dt >= "1960-01-01", ]

sbp_gp <- sbp_gp[sbp_gp$event_dt >= "1960-01-01", ]
sbp_gp <- sbp_gp[sbp_gp$event_dt != "2037-07-07", ]

### Save BP GP datasets
saveRDS(sbp_gp, file=paste(hba1c_dat_path, "data/GP_t2d_SBP.rds", sep=""))
saveRDS(dbp_gp, file=paste(hba1c_dat_path, "data/GP_t2d_DBP.rds", sep=""))

############################################################################
### 21. Create BP at T2D diagnosis (baseline) covariates
############################################################################
### Restrict HbA1c LDA dataset created so far to baseline dataset
dt0 <- dt.s[time_base == 0, ]
### Mean ukb values...
ukb_sbp <- dt0 %>% select(eid, msbp)
ukb_dbp <- dt0 %>% select(eid, mdbp)
### Read in UK Biobank assessment covariate dataset
df <- readRDS(paste(analysis_dir, "data/covariate.rds", sep=""))
### Rename columns
df <- bio_rename(df, f)
### Pull out assessment centre dates
dates_df <- df %>% select(eid, date_of_attending_assessment_centre_f53_0_0)
rm(df)
### Add dates to the 2 BP datasets
ukb_sbp <- ukb_sbp %>% left_join(dates_df)
ukb_dbp <- ukb_dbp %>% left_join(dates_df)

ukb_sbp <- ukb_sbp %>% rename(sbp = msbp, event_dt = date_of_attending_assessment_centre_f53_0_0) %>%
    select(eid, event_dt, sbp)
ukb_dbp <- ukb_dbp %>% rename(dbp = mdbp, event_dt = date_of_attending_assessment_centre_f53_0_0) %>%
    select(eid, event_dt, dbp)
### Treating mean like a single observation from UK Biobank assessment data
ukb_sbp$event_dt <- as.Date(ukb_sbp$event_dt)
ukb_dbp$event_dt <- as.Date(ukb_dbp$event_dt)

### Add UKB BP datasets to GP datasets
### Remove repeat rows
sbp_all <- sbp_gp %>% bind_rows(ukb_sbp) %>% distinct()
dbp_all <- dbp_gp %>% bind_rows(ukb_dbp) %>% distinct()

t2d_diag <- dt0 %>% select(eid, t2d_diag_date)
### Work out baseline BP measurements
### Add T2D diagnosis date to BP datasets
sbp_all <- sbp_all %>% left_join(t2d_diag)
dbp_all <- dbp_all %>% left_join(t2d_diag)
### Create time difference column
dbp_all <- dbp_all %>%
    mutate(time_diff = time_length(interval(as.Date(t2d_diag_date), as.Date(event_dt)), "years"))
sbp_all <- sbp_all %>%
    mutate(time_diff = time_length(interval(as.Date(t2d_diag_date), as.Date(event_dt)), "years"))

### Restrict to time window [-5 years, +2 months]
sbp_tmp <- sbp_all[time_diff >= -5, ]
sbp_tmp <- sbp_tmp[time_diff <= 2/12, ]
length(unique(sbp_tmp$eid)) ### 12140

dbp_tmp <- dbp_all[time_diff >= -5, ]
dbp_tmp <- dbp_tmp[time_diff <= 2/12, ]
length(unique(dbp_tmp$eid)) ### 12140
### Create data.table ordering by EID and event date
dbp_tmp <- data.table(dbp_tmp, key=c("eid", "event_dt"))
sbp_tmp <- data.table(sbp_tmp, key=c("eid", "event_dt"))
### Find the closest observation time to T2D diagnosis within this window
absmin_sbp <- sbp_tmp[, min(abs(time_diff)), by="eid"]
absmin_dbp <- dbp_tmp[, min(abs(time_diff)), by="eid"]
### Find the observation(s) matching this closest time (abs)
sbp_tmp <- sbp_tmp %>% left_join(absmin_sbp) %>% rename(absmin = V1)
dbp_tmp <- dbp_tmp %>% left_join(absmin_dbp) %>% rename(absmin = V1)
sbp_tmp <- sbp_tmp[sbp_tmp$absmin == abs(sbp_tmp$time_diff), ]
sbp_tmp <- sbp_tmp %>% distinct()
dbp_tmp <- dbp_tmp[dbp_tmp$absmin == abs(dbp_tmp$time_diff), ]
dbp_tmp <- dbp_tmp %>% distinct()
### Count the rows to check for ties in data
sbp_count <- sbp_tmp[, .(rowCount = .N), by=eid]
dbp_count <- dbp_tmp[, .(rowCount = .N), by=eid]
sbp_tmp <- sbp_tmp %>% left_join(sbp_count)
dbp_tmp <- dbp_tmp %>% left_join(dbp_count) 
sbp_tmp2 <- sbp_tmp %>% select(eid, event_dt, sbp)
dbp_tmp2 <- dbp_tmp %>% select(eid, event_dt, dbp)
dbp_tmp2 <- data.table(dbp_tmp2, key=c("eid", "event_dt"))
sbp_tmp2 <- data.table(sbp_tmp2, key=c("eid", "event_dt"))
### Mean of ties used
sbp_tmp2 <- dcast(sbp_tmp2, eid ~ ., fun.aggregate = mean,
                value.var = 'sbp')
dbp_tmp2 <- dcast(dbp_tmp2, eid ~ ., fun.aggregate = mean, 
                value.var = 'dbp')

colnames(sbp_tmp2)[2] <- "sbp_base"
colnames(dbp_tmp2)[2] <- "dbp_base"
### Save BP baseline data
saveRDS(sbp_tmp2, file=paste(analysis_dir, "data/SBP_baseline.rds", sep=""))
saveRDS(dbp_tmp2, file=paste(analysis_dir, "data/DBP_baseline.rds", sep=""))

### Add baseline BP information to HbA1c LDA dataset
dt.s <- data.table(dt.s, key=c("eid", "event_dt"))
dt.s <- dt.s %>% left_join(sbp_tmp2)
dt.s <- dt.s %>% left_join(dbp_tmp2)

############################################################################
### 22. Create additional BP variables
# count of measurements pre T2D, median SBP, median DBP, max SBP, max DBP,
# mean SBP, mean DBP
### Some are used to explore missing data
############################################################################
sbp_all <- data.table(sbp_all, key=c("eid", "event_dt"))
dbp_all <- data.table(dbp_all, key=c("eid", "event_dt"))
sbp_pre <- sbp_all[time_diff <= 0, ]
sbprow_pre <- sbp_pre[, .(rowCount = .N), by=eid]
colnames(sbprow_pre)[2] <- "nobs_pre_SBP"
dbp_pre <- dbp_all[time_diff <= 0, ]
dbprow_pre <- dbp_pre[, .(rowCount = .N), by=eid]
colnames(dbprow_pre)[2] <- "nobs_pre_DBP"
### Add to main dataset
dt.s <- data.table(dt.s, key=c("eid", "event_dt"))
dt.s <- dt.s %>% left_join(sbprow_pre)
dt.s <- dt.s %>% left_join(dbprow_pre)

sbp_median <- dcast(sbp_all, eid ~ ., fun.aggregate = median, 
                value.var = 'sbp')
colnames(sbp_median)[2] <- "sbp_median"

dbp_median <- dcast(dbp_all, eid ~ ., fun.aggregate = median, 
                value.var = 'dbp')
colnames(dbp_median)[2] <- "dbp_median"

sbp_max <- dcast(sbp_all, eid ~ ., fun.aggregate = max, 
                value.var = 'sbp')
colnames(sbp_max)[2] <- "sbp_max"

dbp_max <- dcast(dbp_all, eid ~ ., fun.aggregate = max, 
                value.var = 'dbp')
colnames(dbp_max)[2] <- "dbp_max"

sbp_mean <- dcast(sbp_all, eid ~ ., fun.aggregate = mean, 
                value.var = 'sbp')
colnames(sbp_mean)[2] <- "sbp_mean"

dbp_mean <- dcast(dbp_all, eid ~ ., fun.aggregate = mean, 
                value.var = 'dbp')
colnames(dbp_mean)[2] <- "dbp_mean"

dt.s <- data.table(dt.s, key=c("eid", "event_dt"))
dt.s <- dt.s %>% left_join(sbp_median)
dt.s <- dt.s %>% left_join(dbp_median)
dt.s <- dt.s %>% left_join(sbp_max)
dt.s <- dt.s %>% left_join(dbp_max)
dt.s <- dt.s %>% left_join(sbp_mean)
dt.s <- dt.s %>% left_join(dbp_mean)

sum(is.na(dt.s$nobs_pre_SBP))
sum(is.na(dt.s$nobs_pre_DBP))
dt.s$nobs_pre_SBP[is.na(dt.s$nobs_pre_SBP)] <- 0
dt.s$nobs_pre_DBP[is.na(dt.s$nobs_pre_DBP)] <- 0

############################################################################
### 23. Create retrospective MDD grouping variable
############################################################################
dt.s <- dt.s %>%
    mutate(dep_group = ifelse(dep_post == 1, "t2d_dep", ifelse(dep_base == 1, "dep_t2d", "control")))
### Save dataset
saveRDS(dt.s, paste(hba1c_dat_path, "data/t2d_imputeprep_ds.rds", sep=""))

############################################################################
### 24. Standardise continuous imputation and LDA variables
############################################################################
### For continuous variables to be used in LDA (and imputation model) create baseline standardised versions...
dt.s <- readRDS(paste(hba1c_dat_path, "data/t2d_imputeprep_ds.rds", sep=""))
### Make data.table with ordering by EID and event date
dt.s <- data.table(dt.s, key=c("eid", "event_dt"))
### Create baseline dataset
dt0 <- dt.s[time_base ==0, ]
### Create a standardised parameter dataset to save for later use
std_params <- rbind(c(mean(dt0$hba1c, na.rm=T), sd(dt0$hba1c, na.rm=T)),
c(mean(dt0$yo_diag, na.rm=T), sd(dt0$yo_diag, na.rm=T)),
c(mean(dt0$age_t2d_diag, na.rm=T), sd(dt0$age_t2d_diag, na.rm=T)), 
c(mean(dt0$bmi, na.rm=T), sd(dt0$bmi, na.rm=T)),
c(mean(dt0$tdi, na.rm=T), sd(dt0$tdi, na.rm=T)),
c(mean(dt0$sbp_base, na.rm=T), sd(dt0$sbp_base, na.rm=T)),
c(mean(dt0$dbp_base, na.rm=T), sd(dt0$dbp_base, na.rm=T)),
c(mean(dt0$sbp_max, na.rm=T), sd(dt0$sbp_max, na.rm=T)), 
c(mean(dt0$mdbp, na.rm=T), sd(dt0$mdbp, na.rm=T)))
rownames(std_params) <- c("hba1c", "yo_diag", "age_t2d_diag", "bmi", "tdi", "sbp_base", "dbp_base", "sbp_max", "mdbp")
colnames(std_params) <- c("Mean_baseline", "SD_baseline")
### save standardising parameters
saveRDS(std_params, paste(hba1c_dat_path, "data/cont_vars_mean_SD_baseline.rds", sep=""))

### Then use these to standardise columns accordingly
dt.s <- dt.s %>%
    mutate(hba1c_std = (hba1c-std_params[1,1])/std_params[1,2],
    yo_diag_std = (yo_diag - std_params[2,1])/std_params[2,2],
    age_t2d_diag_std = (age_t2d_diag - std_params[3,1])/std_params[3,2],
    bmi_std = (bmi - std_params[4,1])/std_params[4,2],
    tdi_std = (tdi - std_params[5,1])/std_params[5,2],
    sbp_base_std = (sbp_base - std_params[6,1])/std_params[6,2],
    dbp_base_std = (dbp_base - std_params[7,1])/std_params[7,2],
    sbp_max_std = (sbp_max - std_params[8,1])/std_params[8,2],
    mdbp_std = (mdbp - std_params[9,1])/std_params[9,2])

dt.s <- dt.s %>% mutate(hba1c_base_std = (hba1c_base-std_params[1,1])/std_params[1,2])

dt.s <- data.table(dt.s, key=c("eid", "event_dt"))
saveRDS(dt.s, paste(hba1c_dat_path, "data/t2d_imputeprep_ds.rds", sep=""))

############################################################################
### 25. MDD exposures (duration variables created)
############################################################################
### dep_change_base: time-varying depression minus depression at baseline
### depression duration variables for pre and post T2D MDD individuals
dt.s <- readRDS(paste(hba1c_dat_path, "data/t2d_imputeprep_ds.rds", sep=""))
dt.s <- data.table(dt.s, key=c("eid", "time_base"))
dt.s <- dt.s %>% mutate(dep_change_base = dep_tv - dep_base)

saveRDS(dt.s, paste(hba1c_dat_path, "data/t2d_imputeprep_ds.rds", sep=""))

### Add in some extra depression duration variables for pre and post T2D MDD individuals...
### pret2dmdd_time: time between MDD diagnosis and T2D diagnosis for pret2dmdd individuals- 0 otherwise. Time invariant.
### postt2dmdd_time: time between MDD diagnosis and T2D diagnosis for postt2dmdd individuals- 0 otherwise. Time varying.
### Read in main LDA dataset
dt.s <- readRDS(paste(hba1c_dat_path, "data/t2d_imputeprep_ds.rds", sep=""))
dt.s <- data.table(dt.s, key=c("eid", "time_base"))

### dep_time_og is each HbA1c event relative to depression
dt.s$dep_time_og <- time_length(interval(as.Date(dt.s$dep_first_date), as.Date(dt.s$event_dt)), "years")
### Will be NA for those without depression- set to 0
dt.s$dep_time_og[is.na(dt.s$dep_time_og)] <- 0
### Want a post-T2D MDD time variable which is this og depression time multipled by dep_change_base
dt.s$postt2dmdd_time <- dt.s$dep_time_og*dt.s$dep_change_base
### Ensure data.table ordering is maintained
dt.s <- data.table(dt.s, key=c("eid", "time_base"))
range(dt.s$postt2dmdd_time) ### [1] 0.000000 9.592896

### Want a pre-T2D MDD time variable
### This will be time between MDD and T2D diagnosis
dt.s$pret2dmdd_time <- time_length(interval(as.Date(dt.s$dep_first_date), as.Date(dt.s$t2d_diag_date)), "years")
dt.s$pret2dmdd_time[is.na(dt.s$pret2dmdd_time)] <- 0
dt.s$pret2dmdd_time <- dt.s$pret2dmdd_time*dt.s$dep_base
range(dt.s$pret2dmdd_time)

dt.s <- data.table(dt.s, key=c("eid", "time_base"))
saveRDS(dt.s, paste(hba1c_dat_path, "data/t2d_imputeprep_ds.rds", sep=""))

############################################################################
### 26. Maximum follow-up calculated
############################################################################
### Add max followup
dt.s <- data.table(dt.s, key=c("eid", "event_dt"))
max_fu <- dt.s[, max(time_base), by="eid"]
colnames(max_fu)[2] <- "followup"
dt.s <- dt.s %>% left_join(max_fu)
dt.s <- data.table(dt.s, key=c("eid", "event_dt"))
saveRDS(dt.s, paste(hba1c_dat_path, "data/t2d_imputeprep_ds.rds", sep=""))

############################################################################
### 27. Alcohol never created
############################################################################
### alcohol never
dt.s <- dt.s[, alcohol_never := as.numeric(alc_drinker_status_f == "Never")]
saveRDS(dt.s, paste(hba1c_dat_path, "data/t2d_imputeprep_ds.rds", sep=""))

############################################################################
### 28. Add initial BMI at UK Biobank assessment into LDA dataset 
### (auxiliary MI variable)
############################################################################
### Read in UK Biobank assessment covariates
df <- readRDS(paste(analysis_dir, "data/covariate.rds", sep=""))
### Rename columns
df <- bio_rename(df, f)
### Ensure is a data frame
df <- data.frame(df)
### Select first instance of BMI (i.e. at initial UK Biobank assessment
df_bmi0 <- df %>% dplyr::select(eid, "body_mass_index_.bmi._f21001_0_0")
dt.s <- data.table(dt.s, key=c("eid", "event_dt"))
colnames(df_bmi0)[2] <- "bmi_ukb"
### Add into main dataset
dt.s <- dt.s %>% left_join(df_bmi0)
## Standardised
dt.s$bmi_ukb_std <- (dt.s$bmi_ukb - mean(dt.s$bmi_ukb, na.rm=T))/sd(dt.s$bmi_ukb, na.rm=T)

dt.s <- data.table(dt.s, key=c("eid", "event_dt"))
saveRDS(dt.s, paste(hba1c_dat_path, "data/t2d_imputeprep_ds.rds", sep=""))

############################################################################
### 29. Sum of HbA1c, BP and BMI measurements before T2D diagnosis
############################################################################
### Read in HbA1c LDA dataset so far:
dt.s <- readRDS(paste(hba1c_dat_path, "data/t2d_imputeprep_ds.rds", sep=""))
## Make data.table ordered by EID and time_base (so event_dt)
dt.s <- data.table(dt.s, key=c("eid", "time_base"))
dt.s <- dt.s %>% mutate(previsits_tot = nobs_pre_DBP + nobs_preBMI + obs_pre)
dt0 <- dt.s[time_base ==0, ]
meanpre_tot <- mean(dt0$previsits_tot) ### 20.48758
sdpre_tot <- sd(dt0$previsits_tot) ### 20.03423
dt.s <- dt.s %>% mutate(previsits_tot_std = (previsits_tot - meanpre_tot)/sdpre_tot)

############################################################################
### 30. Standardise pret2dmdd_time (robust scaler standardising)
############################################################################
### standardise pret2dmdd_time
dt.s <- data.table(dt.s, key=c("eid", "time_base"))
dt0pre <- dt0[dep_base == 1, ]
premean <- mean(dt0pre$pret2dmdd_time) ### 12.90697
presd <- sd(dt0pre$pret2dmdd_time) ### 10.05384
dt.s <- dt.s %>% mutate(pret2dmdd_time_std = ifelse(dep_base == 1, (pret2dmdd_time - premean)/presd, 0))

### percentile for pret2dmdd_time
dt.s <- dt.s %>%
        mutate(percentile_pretime = ifelse(dep_base == 0, 0 , pnorm(-1*pret2dmdd_time, mean=-1*premean, sd=presd)))
### Robust-scaler for pret2dmdd_time
Q1 <- quantile(dt0pre$pret2dmdd_time, probs =0.25) ### 5.026712
Q3 <- quantile(dt0pre$pret2dmdd_time, probs =0.75) ### 17.61943
dt.s <- dt.s %>%
        mutate(robust_pretime = ifelse(dep_base == 0, 0 , (pret2dmdd_time - Q1)/(Q3-Q1)))
dt.s <- dt.s %>%
        mutate(uncentered_robust_pretime = ifelse(dep_base == 0, 0 , (pret2dmdd_time)/(Q3-Q1)))
dt.s <- data.table(dt.s, key=c("eid", "time_base"))
range(dt.s$uncentered_robust_pretime[dt.s$dep_base == 1])
dt.s <- dt.s %>%
        mutate(uncentered_robust_pretime2 = ifelse(dep_base == 0, 0 , (pret2dmdd_time + (2/52))/(Q3-Q1)))
dt.s <- data.table(dt.s, key=c("eid", "time_base"))

############################################################################
### 31. Subgroup variable for pre-T2D MDD (exploratory)
############################################################################
### Subgroup for pre_time...
### Not used in main analysis
dt.s <- dt.s %>%
        mutate(dep_base_sg = ifelse(pret2dmdd_time > 15.865753, "gt75pc", "25pc_75pc")) %>%
        mutate(dep_base_sg = ifelse(pret2dmdd_time <= 4.314208, "10pc_25pc", dep_base_sg)) %>%
        mutate(dep_base_sg = ifelse(pret2dmdd_time <= 1.786301, "le10pc", dep_base_sg)) %>%
        mutate(dep_base_sg = ifelse(dep_base == 0, "depbase0", dep_base_sg))
dt.s <- data.table(dt.s, key=c("eid", "time_base"))
dt.s <- dt.s %>%
        mutate(dep_base_sg2 = ifelse(pret2dmdd_time > 15.865753, 1, 2)) %>%
        mutate(dep_base_sg2 = ifelse(pret2dmdd_time <= 4.314208, 3, dep_base_sg2)) %>%
        mutate(dep_base_sg2 = ifelse(pret2dmdd_time <= 1.786301, 4, dep_base_sg2)) %>%
        mutate(dep_base_sg2 = ifelse(dep_base == 0, 0, dep_base_sg2))
dt.s <- dt.s %>% mutate(dep_base_sg_f = as.factor(dep_base_sg2))
dt.s <- data.table(dt.s, key=c("eid", "time_base"))

############################################################################
### 32. Splines for time since post-T2D MDD diagnosis
############################################################################
### Splines for postt2dmdd_time, use knot = 4.
dtchange <- dt.s[dep_change_base == 1, ]
library(rms)
postspline <- rcs(dtchange$postt2dmdd_time, 4) 
attributes(postspline)$parms ### [1] 0.2745205 1.9666532 4.0706849 7.8142466
postspline <- rcs(dt.s$postt2dmdd_time, parms=c(0.2745205, 1.9666532, 4.0706849, 7.8142466))
dt.s$postspline1 <- as.vector(postspline[,1])
dt.s$postspline2 <- as.vector(postspline[,2])
dt.s$postspline3 <- as.vector(postspline[,3])
dt.s <- data.table(dt.s, key=c("eid", "time_base"))

# table(dt.s$postt2dmdd_time - dt.s$postspline1) ### all zero as expected

saveRDS(dt.s, paste(hba1c_dat_path, "data/t2d_imputeprep_ds.rds", sep=""))

############################################################################
### 33. Plot of time between pre-T2D MDD diagnosis and T2D diagnosis
### Various standardisations of pret2dmdd_time
############################################################################
datw <- dt.s[dt.s$dep_base == 1, ]
pdf(paste(analysis_dir, "plots/pret2d_time_hba1c_loess.pdf", sep=""))
   ggplot(data = datw, aes(x = pret2dmdd_time, y = hba1c))+ geom_point() + geom_smooth(method = "loess")
   ggplot(data = datw, aes(x = pret2dmdd_time_std, y = hba1c))+ geom_point() + geom_smooth(method = "loess") 
   ggplot(data = datw, aes(x = percentile_pretime, y = hba1c))+ geom_point() + geom_smooth(method = "loess") 
   ggplot(data = datw, aes(x = robust_pretime, y = hba1c))+ geom_point() + geom_smooth(method = "loess") 
   ggplot(data = datw, aes(x = uncentered_robust_pretime2, y = hba1c))+ geom_point() + geom_smooth(method = "loess")    
dev.off()

############################################################################
### 34. Percentage diagnosed outside of 2006 and 2010
############################################################################
dt.s <- readRDS(paste(hba1c_dat_path, "data/t2d_imputeprep_ds.rds", sep=""))
dt.s <- data.table(dt.s, key=c("eid", "time_base"))
dt.s1 <- dt.s[lda_keep == 1, ]
dt.s2 <- dt.s1[, unique(t2d_diag_date), by="eid"]
dt.s2a <- dt.s2[dt.s2$V1 < "2006-01-01", ]
dt.s2b <- dt.s2[dt.s2$V1 > "2010-12-31", ]
(dim(dt.s2a)[1] + dim(dt.s2b)[1])/dim(dt.s2)[1] ### 0.5869274
