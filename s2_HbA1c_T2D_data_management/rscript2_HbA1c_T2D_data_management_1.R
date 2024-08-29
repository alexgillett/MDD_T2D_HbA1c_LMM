############################################################################
### 2.1. HbA1c data management
############################################################################
### Authors: Saskia Hagenaars, Alexandra C Gillett
############################################################################
### In script we:
### 1. Merge GP extracted HbA1c datasets and UK Biobank assessment HbA1c dataset.
### 2. Merge main HbA1c dataset with MDD case/ control data.
### 3. Extract and add in date of birth information to main HbA1c dataset.
### 4. Prepare diabetic medication data.
### 4b. OPTIONAL: Merge diabetic medication data and main HbA1c dataset.
### 5. In main HbA1c dataset, create medication at HbA1c observation columns
### 6. Create the T2D case-control phenotype
### 6a. HbA1c < 48 mmol/mol
### 6b. Self-reported phenotype from UK Biobank assessment data
### 6c. HES diagnostic code for T2D
### 6d. GP T2D diagnosis
### 6e. GP Medication for T2D
### 6f. Define T2D phenotype using 6a-e
### 6g. First date of recorded diabetes
### 7. T2D case exclusions
### 8. Add in MDD diagnosis date
### 9. Additional exclusion criteria
### 10. Create T2D related datasets for PRS analysis
############################################################################
# User add library path where R packages are stored if required:
#.libPaths(new = "your_rlibrary_path")
# Data extraction uses ukbkings package
# https://github.com/kenhanscombe/ukbkings
# ukbkings requires data to be formatted using ukbproject
# https://github.com/kenhanscombe/ukbproject
############################################################################
### Add R packages:
library(data.table)
library(ukbkings)
library(tidyverse)
library(lubridate)
library(readxl)
library(knitr)
library(xtable)
library(ggplot2)
library(readr)
library(zoo)
library(DescTools)
### Assign paths to relevant directories
out_path <- "/path_to_where_you_store_extracted_data/hba1c_data/"
project_dir <- "/path_to_your_ukb_data_dir/"
depression_dir <- "/path_to_extracted_depression_data/"
analysis_dir <- "/path_to_analysis_dir/"
############################################################################
### 1. Merge GP extracted HbA1c datasets and UK Biobank assessment HbA1c dataset
############################################################################
### a) Merge GP extracted HbA1c datasets
### Read in hba1c data extracted from primary care data
### This is separate read2 and read3 data atm
hba1c_read2<-fread(paste(out_path, "gp_HbA1c_read2_extracted.txt", sep=""))
hba1c_read3<-fread(paste(out_path, "gp_HbA1c_read3_extracted.txt", sep=""))
### Merge read2 and read3, replace missing with NA
### NAs for read2- will get an error message about NA by coersion...
hba1c_read2$raw_HbA1c_results[hba1c_read2$raw_HbA1c_results == "^"] <- ""
hba1c_read2_NA <- hba1c_read2 %>%
select(eid,event_dt, raw_HbA1c_results) %>%
mutate(raw_HbA1c_results = na_if(raw_HbA1c_results, "")) %>%
mutate(raw_HbA1c_results=as.numeric(raw_HbA1c_results))
hba1c_read3_NA <- hba1c_read3 %>% 
select(eid,event_dt, raw_HbA1c_results)
### Merge
hba1c_tab <-  hba1c_read2_NA %>%
bind_rows(hba1c_read3_NA) %>% rename(raw_HbA1c_result=raw_HbA1c_results)
### Remove duplicates and incorrect dates
### remove any duplicates that have been added due to the way we extract the data
### remove values with incorrect dates (see UKB primary care documentation for more detail)
dates_rm <- c("01/01/1901", "02/02/1902", "03/03/1903", "07/07/2037", "")
hba1c_tab <- hba1c_tab[!(hba1c_tab$event_dt %in% dates_rm), ]
hba1c_tab$event_dt <- as.Date(hba1c_tab$event_dt, "%d/%m/%Y")
range(hba1c_tab$event_dt) ### "1950-07-28" "2017-09-27"
hba1c_tab <- hba1c_tab[hba1c_tab$event_dt != "1950-07-28", ]
range(hba1c_tab$event_dt) ### "1989-02-20" "2017-09-27"
hba1c_tab <- hba1c_tab%>%
    distinct()


### Transform HbA1c values to IFCC values instead of DCCT values:
### remove values between 15 and 20 as we can't be sure whether these values are very low IFCC or very high DCCT values
hba1c_tab <- hba1c_tab[hba1c_tab$raw_HbA1c_result >= 3.5, ]
hba1c_tab <- hba1c_tab[hba1c_tab$raw_HbA1c_result <= 195, ]
keep <- rep(0, dim(hba1c_tab)[1])
keep[hba1c_tab$raw_HbA1c_result < 15]<-1
keep[hba1c_tab$raw_HbA1c_result >= 20]<-1
hba1c_tab <- hba1c_tab[keep == 1, ]
dim(hba1c_tab) ### For reference, values we have are:  [1] 581446      3
length(unique(hba1c_tab$eid)) ### For reference, values we have are:  106506
rm(keep)
### Use appropriate transformation (see SuppMethods of paper for details)
hba1c_tab <- hba1c_tab %>%
mutate(hba1c_v2 = ifelse((raw_HbA1c_result >= 3.5 & raw_HbA1c_result <20), round((10.93*raw_HbA1c_result)-23.50), raw_HbA1c_result),
	type="hba1c")
### Save GP merged and tranformed data:
fwrite(hba1c_tab, paste(out_path, "data/hba1c_read2_read3_original_transformed_values.txt", sep=""),
  col.names=TRUE, row.names=FALSE, quote=FALSE)

### remove original read2, read3 HbA1c files:
rm(hba1c_read2, hba1c_read3)

### b) Merge UK Biobank biomarker data to main HbA1c dataset
### Read in UK Biobank assessment biomarker data for HbA1c
HbA1c_biomarker_separate <- fread(paste(out_path, "data/HbA1c_biomarker_rowwise.txt", sep=""), header=T)

### Create list of ids to filter down biomarker data to people who have GP HbA1c
### Make sure columns names are suitable for merging
### Make sure dates are in correct format in both datasets
hba1c_ids <- hba1c_tab %>%
distinct(eid, .keep_all=FALSE)

HbA1c_biomarker_separate <- HbA1c_biomarker_separate %>%
filter(eid %in% hba1c_ids$eid) %>%
rename(type = hba1c_type)%>%
mutate(event_dt = as.Date(event_dt, format="%Y-%m-%d"))

hba1c_tab <- hba1c_tab %>%
  mutate(event_dt = as.Date(event_dt, format="%Y-%m-%d"))

### Merge
hba1c_tab <- hba1c_tab %>%
rename(hba1c = hba1c_v2)%>%
bind_rows(HbA1c_biomarker_separate)

############################################################################
### 2. Merge main HbA1c dataset with depression case/ control data
############################################################################
### Add depression cases/controls to HbA1c data
### Depression EIDs:
dep_eid <- fread(paste(depression_dir, "dep_AD_eids.txt", sep=""), header=T)
### Binary yes/no depression variable for those in HbA1c main dataset:
depression <- as.numeric(hba1c_tab$eid %in% dep_eid$eid)
### Add binary depression variable to main HbA1c dataset
hba1c_tab_dep <- hba1c_tab
hba1c_tab_dep$depression <- depression
rm(depression)
############################################################################
### 3. Add in approximate date of birth
############################################################################
### add date of birth
### get date of birth (month/year only)
### Store required data fields to txt file (created in step0_basic_cov_file.sh)
f %>%
select(field, name) %>% 
filter(str_detect(field, "^34-0.0|^52-0.0")) %>%
bio_field_add(paste(out_path, "month_year_birth.txt", sep=""))
### Extract fields and store data
bio_phen(project_dir,
	field=paste(out_path, "month_year_birth.txt", sep=""),
	out=paste(out_path, "month_year_birth", sep=""))
### Read in data
dob=readRDS(paste(out_path, "month_year_birth.rds", sep=""))
### Rename columns, create approximate date of birth with birth day being 1st of month
### This is a common approach in UK Biobank data
dob <- dob %>%
rename("year_of_birth"=2,
	"month_of_birth"=3) %>%
mutate(dob=paste0("1-",month_of_birth,"-",year_of_birth),
	dob=as.Date(dob, "%d-%m-%Y"))
### Add dob to main HbA1c dataset
hba1c_tab_dep <- hba1c_tab_dep %>%
group_by(eid) %>%
left_join(dob[,c(1,4)]) %>%
ungroup() 
### Store intermediate output:
fwrite(hba1c_tab_dep, paste(out_path, "data/HbA1c_data_clean_dep.txt", sep=""),col.names=TRUE,row.names=FALSE,quote=FALSE, sep="\t")

############################################################################
### 4. Prepare diabetic medication data
############################################################################
### load extracted medication data
drug_data <- fread(paste(out_path, "gp_HbA1c_drug_data_names_extracted.txt", sep=""))
drug_read2 <- fread(paste(out_path, "gp_HbA1c_drug_data_codes_extracted.txt", sep=""))
med_names <- fread(paste(out_path, "inputs/drug_names.txt", sep=""),header=T)

### Because the drug name in the gp_scripts file has the actual name of the drug it doesn't match with drug_name in mapping file
### split drug name to extract only first word of the name (and check if that works)
names <- as.character(drug_data$drug_name)
a <- unique(names)

### Start with bnf/dmd extracted data:
### check drug_data...
drug_data %>%
select(-drug_substance,-drug_class1,-drug_class2)%>%
separate(drug_name, c("drug_name","rest"), sep="([\\s/-2,])", extra="merge") -> test #separators are: any white space, /, -,2, ","
test %>% 
separate(drug_name, c("drug_name","rest2"), sep="-", extra="merge") -> test2 # for some reason i couldn't sep by just the '-' (for glibenclamide-p42) in previous line
###  NB. rest2 column contains no additional drugs... Focus on 'rest' column
### Following loop explores 'rest' column create above to see if additional drug
### names are contained in this info
### If there is this is put into a column called 'drug_name2' in 'test2' dataset
test2$drug_name2 <- NA
names_w <- as.character(med_names$drug_name)
for(i in 1:length(names_w)){
  check_i <- as.numeric(sapply(as.character(tolower(test2$rest)), grepl, pattern = names_w[i]))
  if(sum(check_i) > 0){
    test2$drug_name2[check_i == 1] <- names_w[i] 
  }
}
### Create a test2 dataset which contains non-NA drug_name2 rows
test3 <- test2[!is.na(test2$drug_name2), ]
test3 %>% 
  mutate(drug_name=drug_name2) %>%
  select(-drug_name2) -> test3

colnames(test2)
colnames(test3)
### In test2 remove drug_name2, covert drug_name to lower case, and bind rows with test3 (which has additonal drug names
test2 %>% select(-drug_name2) %>%
  mutate(drug_name=tolower(drug_name)) %>%
  bind_rows(test3) -> test2
rm(test3)

### Add in prescription column to say these are prescription data
### Rename issue_date to event_dt, in-keeping with date naming in HbA1c dataset
### Format event_dt as date.
test2 %>%
	mutate(drug_name=tolower(drug_name),
    type="prescription") %>%
    rename(event_dt=issue_date) %>%
    mutate(event_dt=as.Date(event_dt, "%d/%m/%Y")) -> test2

### Add in additional medication information specific to each drug contained in med_names
test2 %>%
left_join(med_names, by="drug_name") %>%
select(-read_2,-bnf_code,-dmd_code,-rest2,-rest) -> drug_data
#head(drug_data)
### Remove rows with missing dates
drug_data <- drug_data[!is.na(drug_data$event_dt), ]
### Remove rows with incorrect dates
drug_data <- drug_data[drug_data$event_dt != "1902-02-02", ]

### Read2 extracted drug data
drug_read2 %>%
rename(event_dt=issue_date) %>%
mutate(event_dt=as.Date(event_dt, "%d/%m/%Y"),
	type="prescription") %>%
select(eid, data_provider, event_dt, drug_name, quantity,type, drug_substance,drug_class1, drug_class2) -> drug_read2

### Combine drug name (bnf/dmd) and read2 files and make subset based on people who have hba1c values
drug_data %>%
bind_rows(drug_read2) %>%
distinct() %>%
filter(eid %in% hba1c_ids$eid) %>%
mutate(event_dt=as.Date(event_dt, "%Y-%m-%d")) -> drug_data_hba1csubset

drug_data_hba1csubset$drug_class2[which(drug_data_hba1csubset$drug_class2=="")] <- NA

drug_data_hba1csubset %>%
mutate(drug=ifelse(!is.na(drug_class2), paste0(drug_class1,"_",drug_class2), drug_class1))->drug_data_hba1csubset

table(drug_data_hba1csubset$drug)
sum(is.na(drug_data_hba1csubset$drug)) # 0

### Create prescription episodes for each drug (max 6 months between prescriptions)

drug_data_hba1csubset %>%
  group_by(eid) %>%
  arrange(event_dt) %>%
  mutate(prev_drug = lag(drug)) %>%
  ungroup() %>%
  group_by(eid, drug) %>%
  mutate(prev_drug_date = dplyr::lag(event_dt, n = 1, default = NA),
    diff_weeks_drug = as.numeric(difftime(event_dt, prev_drug_date, units = "weeks")),
    prescription_episode = ifelse(diff_weeks_drug > 26, seq_along(diff_weeks_drug), 1),
    prescription_episode = ifelse(is.na(prescription_episode), 999, prescription_episode),
    prescription_episode = ifelse( prescription_episode == 1,
  NA, prescription_episode)) %>%
  fill(prescription_episode) %>%
  mutate(prescription_episode = ifelse (prescription_episode == 999, 1, prescription_episode)) %>%
  ungroup() %>%
  group_by(eid,drug,prescription_episode) %>%
  arrange(event_dt) %>%
  mutate(med_start = min(event_dt, na.rm=T)) %>%
  mutate(med_end = max(event_dt, na.rm=T)) %>%
  ungroup() -> drug_data_hba1csubset_v2


### Save medication file
fwrite(drug_data_hba1csubset_v2, paste(out_path, "data/gp_diabetes_medication_data.txt", sep=""),
  col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

############################################################################
### 4b. Optional: Merge diabetic medication data and main HbA1c dataset
############################################################################
### READ IN hba1c_tab_dep
### READ IN drug_data_hba1csubset_v2
hba1c_tab_dep <- fread(paste(out_path, "data/HbA1c_data_clean_dep.txt", sep=""))
drug_data_hba1csubset_v2 <- fread(paste(out_path, "data/gp_diabetes_medication_data.txt", sep=""))
### Make sure date columns are formatted a dates using as.Date:
drug_data_hba1csubset_v2$event_dt <- as.Date(drug_data_hba1csubset_v2$event_dt)
drug_data_hba1csubset_v2$prev_drug_date <- as.Date(drug_data_hba1csubset_v2$prev_drug_date)
drug_data_hba1csubset_v2$med_start <- as.Date(drug_data_hba1csubset_v2$med_start)
drug_data_hba1csubset_v2$med_end <- as.Date(drug_data_hba1csubset_v2$med_end)
hba1c_tab_dep$event_dt <- as.Date(hba1c_tab_dep$event_dt)
hba1c_tab_dep$dob <- as.Date(hba1c_tab_dep$dob)

### EIDs with HbA1c data available:
hba1c_ids <- unique(hba1c_tab_dep$eid)


### Restrict drug dataset to HbA1c EIDs:
drug_data_hba1csubset_v2 <- drug_data_hba1csubset_v2[(drug_data_hba1csubset_v2$eid %in% hba1c_ids), ]
dim(drug_data_hba1csubset_v2)


### Bind rows for HbA1c data and prescription data
hba1c_tab_dep %>%
bind_rows(drug_data_hba1csubset_v2) -> hba1c_v2


### save combined medication hba1c file
fwrite(hba1c_v2, paste(out_path, "data/gp_hba1c_medication_data.txt", sep=""),
  col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
rm(drug_data_hba1csubset)
rm(dep_eid, dob)

############################################################################
### 5. In main HbA1c dataset, create medication at HbA1c observation columns
############################################################################
### create medication combinations based on what medication people were taking at the time of HbA1c measurement
### Create hba1c interval (hba1c - 3 months) within a new dataset hba1c_v3
hba1c_tab_dep %>%
  mutate(hba1c_3m_bound = event_dt %m-% months(3)) -> hba1c_v3

### Create (blank for now) columns to fill with 1 = yes and 0 = no for drug prescription in the 3 months prior to HbA1c observation
hba1c_v3$metformin <- NA
hba1c_v3$sulfonylurea <- NA
hba1c_v3$TZD <- NA
hba1c_v3$meglitinide <- NA
hba1c_v3$acarbose <- NA
hba1c_v3$DPP4i <- NA
hba1c_v3$GLP1Ra <- NA
hba1c_v3$SGLT2i <- NA
hba1c_v3$insulin <- NA

### Create a data.table version of working dataset
hba1c_v3dt <- data.table(hba1c_v3, key=c("eid", "event_dt"))
hba1c_eids <- unique(hba1c_v3dt$eid)
### Create a data.table version of the drug dataset
med_dt <- data.table(drug_data_hba1csubset_v2, key=c("eid", "event_dt"))

### Functions to determine if individual received a prescription for a drug
### in the 3 months prior to HbA1c observation.
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

### Apply the above functions(s)
### Essentially, for each unique EID (i) we extract hba1c and medication data for that individual. Then for each HbA1c event_dt for individual i, we search for the drug prescriptions within the 3 months prior to that date.
### The output contains a column for each medication (metformin = col1,
### sulfonylurea = col2, TZD = col3, meglitinide = col4, acarbose = col5,
### DPP4i = col6, GLP1Ra = col7, SGLT2i = col8, insulin = col9.
meds_out <- mapply(x=hba1c_eids, med_f, MoreArgs = list(hba1c_data = hba1c_v3dt, med_data = med_dt),SIMPLIFY = TRUE)
meds_out <- do.call(rbind, meds_out)
### Add this drug data into the data.table for hba1c (which is ordered by eid and event_dt to match output in meds_out
hba1c_v3dt$metformin <- meds_out[,1]
hba1c_v3dt$sulfonylurea <- meds_out[,2]
hba1c_v3dt$TZD <- meds_out[,3]
hba1c_v3dt$meglitinide <- meds_out[,4]
hba1c_v3dt$acarbose <- meds_out[,5]
hba1c_v3dt$DPP4i <- meds_out[,6]
hba1c_v3dt$GLP1Ra <- meds_out[,7]
hba1c_v3dt$SGLT2i <- meds_out[,8]
hba1c_v3dt$insulin <- meds_out[,9]
### save data:
saveRDS(hba1c_v3dt, file=paste(out_path, "data/hba1c_meds_date3mth.rds", sep=""))
### Create column which sums the number of medications being prescribed at each HbA1c event_dt
hba1c_v3dt[, sum_meds:= (metformin + sulfonylurea + TZD + meglitinide + acarbose + DPP4i + GLP1Ra + SGLT2i + insulin)]
### Create medication_coded variable which is a categorical variable with 4 levels:
### 0 = no prescribed medications, 1 = 1 prescribed medication,
### 2 = 2 prescribed medications, and 3 = >2 medications or insulin.
hba1c_v3dt[, medication_coded:= sum_meds]
hba1c_v3dt$medication_coded[hba1c_v3dt$insulin == 1] <- 3
hba1c_v3dt$medication_coded[hba1c_v3dt$sum_meds >= 3] <- 3

### Save dataset:
saveRDS(hba1c_v3dt, file=paste(out_path, "data/HbA1c_full_data_medication_coded.rds", sep=""))

############################################################################
### 6. Create the T2D case-control phenotype
############################################################################
### 6a. HbA1c > 48
############################################################################
hba1c_v9 <- readRDS(file=paste(out_path, "data/HbA1c_full_data_medication_coded.rds", sep=""))

hba1c_v9[, hba1c_over48:= as.numeric(any(hba1c > 48)), by=eid]
hba1c_tab_v2 <- hba1c_v9
############################################################################
### 6b. Self-reported phenotype from UK Biobank assessment data
############################################################################
### Extract relevant field identifiers
f %>%
select(field, name) %>%
filter(str_detect(field, "20002")) -> sr_illness
f %>%
select(field, name) %>%
filter(str_detect(field, "21003")) -> age_assessment
### Add fields to created txt file
f %>%
select(field, name) %>%
filter(str_detect(name, "diabetes")) %>%
slice(1:16) %>%
bind_rows(sr_illness) %>%
bind_rows(age_assessment) %>%
bio_field_add(paste(out_path, "self_report_diabetes_phenotypes.txt", sep=""))
### Extract those fields from UK Biobank data
bio_phen(
 project_dir,
 field = paste(out_path, "self_report_diabetes_phenotypes.txt", sep=""),
 out = paste(out_path, "self_report_diabetes_subset", sep="")
)
### Read in the extracted data:
sr_diabetes <- readRDS(paste(out_path, "self_report_diabetes_subset.rds", sep=""))
dim(sr_diabetes) ### FYI: [1] 502411    158
### What are T2D related codes:
# a) t2d self-report interview (variable 20002, t2d (1223) or generic diabetes(1220))
t2d_codes <- c("1223","1220")
### check for 20002- columns...
t2d_wrk <- sr_diabetes[, c(18:153)]
t2d_sr_1223 <- apply(X=t2d_wrk, FUN = function(X){as.numeric(any(X == "1223", na.rm=T))}, MARGIN=1)
t2d_sr_1220 <- apply(X=t2d_wrk, FUN = function(X){as.numeric(any(X == "1220", na.rm=T))}, MARGIN=1)
t2d_sr_Combined <- as.numeric(t2d_sr_1223+t2d_sr_1220 > 0)
sr_diabetes$t2d_sr_Combined<-t2d_sr_Combined
print(table(sr_diabetes$t2d_sr_Combined))

# b) diabetes self-report interview touchscreen (variable 2443)
diab_touch_code <- "1"
### Identify 2443 columns
t2d_wrk <- sr_diabetes[, 2:5]
diab_touch_Combined <- apply(X=t2d_wrk, FUN = function(X){as.numeric(any(X == "1", na.rm=T))}, MARGIN=1)
sr_diabetes$diab_touch_Combined<-diab_touch_Combined
print(table(sr_diabetes$diab_touch_Combined))

### combine touchscreen with nurse interview diagnosis

sr_diabetes$t2d_sr_diagnosis <- as.numeric((sr_diabetes$t2d_sr_Combined + sr_diabetes$diab_touch_Combined) > 0)

print(summary(as.factor(sr_diabetes$t2d_sr_diagnosis)))

### exclusion criteria
# other diabetes diagnosis (type 1, gestational, insipidus)
diabetes_other_codes <- c("1222", "1221", "1521")
diab_other <- apply(sr_diabetes[,grep("20002-", colnames(sr_diabetes))], 1, function(row) diabetes_other_codes %in% row)
diab_other_Combined <- numeric()
for(indiv in 1:dim(sr_diabetes)[1]){
  diab_other_Combined[indiv]<-ifelse(sum(diab_other[(((indiv-1)*length(diabetes_other_codes))+1):(indiv*length(diabetes_other_codes))]) == 0, 0, 1)
}
sr_diabetes$diab_other_Combined<-diab_other_Combined
sr_diabetes$diab_other_Combined[c(diab_other_Combined == 1)] <- 1
print(table(sr_diabetes$diab_other_Combined))


# add touchscreen gestational diabetes to this variable
sr_diabetes$diab_other_Combined[which(sr_diabetes$"4041-0.0" == 1 | sr_diabetes$"4041-1.0" == 1 | sr_diabetes$"4041-2.0" == 1 | sr_diabetes$"4041-3.0" == 1)] <- 1
print(table(sr_diabetes$diab_other_Combined))


# insulin within first year of diagnosis
sr_diabetes$insulin_first_year <- NA
sr_diabetes$insulin_first_year[which(sr_diabetes$"2986-0.0" == 1 | sr_diabetes$"2986-1.0" == 1 | sr_diabetes$"2986-2.0" == 1 | sr_diabetes$"2986-3.0" == 1)] <- 1
print(table(sr_diabetes$insulin_first_year))


# individuals with diagnosis under age 35 or missing age at diagnosis (for self-report)
sr_diabetes$under35 <- NA
sr_diabetes$under35[which(sr_diabetes$"2976-0.0" <= 35 | sr_diabetes$"2976-1.0" <= 35 | sr_diabetes$"2976-2.0" <= 35 | sr_diabetes$"2976-3.0" <= 35)] <- 1
print(table(sr_diabetes$under35))

# excluded individuals diagnosed with diabetes within the year prior to the baseline study visit
# unable to determine whether they were using insulin within the first year
sr_diabetes$diff_age_1 <- NA
sr_diabetes$diff_age_1[which((sr_diabetes$"21003-0.0" - sr_diabetes$"2976-0.0" == 1) | 
  (sr_diabetes$"21003-1.0" - sr_diabetes$"2976-1.0" == 1) | 
  (sr_diabetes$"21003-2.0" - sr_diabetes$"2976-2.0" == 1) | 
  (sr_diabetes$"21003-3.0" - sr_diabetes$"2976-3.0" == 1))] <- 1
print(table(sr_diabetes$diff_age_1))

# combine exclusion criteria
sr_diabetes$exclusion_sr <- NA
sr_diabetes$exclusion_sr[which(sr_diabetes$diab_other_Combined == 1 | sr_diabetes$insulin_first_year == 1 | sr_diabetes$under35 == 1 | sr_diabetes$diff_age_1 == 1)] <- 1
print(table(sr_diabetes$exclusion_sr))


## remove exclusion from cases
sr_diabetes$t2d_sr_diagnosis_clean <- sr_diabetes$t2d_sr_diagnosis
sr_diabetes$t2d_sr_diagnosis_clean[which(sr_diabetes$exclusion_sr==1)]<- NA

### save outpt
saveRDS(sr_diabetes, paste(out_path, "data/self_report_diabetes_subset_diagnosis_pheno.rds", sep=""))

sr_diabetes <- readRDS(paste(out_path, "data/self_report_diabetes_subset_diagnosis_pheno.rds", sep=""))
#self_report_diab <- sr_diabetes[, c("eid", "t2d_sr_diagnosis_clean")]
sr_diabetes %>% 
select(eid, t2d_sr_diagnosis_clean) -> self_report_diab

table(self_report_diab$t2d_sr_diagnosis_clean)

sum(is.na(self_report_diab$t2d_sr_diagnosis_clean))

############################################################################
### 6c. HES diagnostic code
############################################################################
### HES DATA
hesin_diag <- bio_record(project_dir, "hesin_diag")

icd9_t2d_hes_code <- c("25000","25010","25020","25090")
icd10_t2d_hes_code <- c("E110","E111","E112","E113","E114","E115","E116","E117","E118","E119")

hes_diab_all <- hesin_diag %>%
filter(diag_icd9 %in% icd9_t2d_hes_code | diag_icd10 %in% icd10_t2d_hes_code) %>%
mutate(hes_diag_dm = 1) %>% collect()

### save
saveRDS(hes_diab_all, paste(out_path, "data/HES_diabetes_subset_diagnosis_pheno.rds", sep=""))

hes_diab_all <- readRDS(paste(out_path, "data/HES_diabetes_subset_diagnosis_pheno.rds", sep=""))

hes_diab_uniqueID <- hes_diab_all %>%
distinct(eid, .keep_all=TRUE) %>%
select(eid, hes_diag_dm) 

############################################################################
### 6d. GP diagnosis
############################################################################
### define people with type 2 diabetes using gp data
### have EIDs...
dm_read2=fread(paste(out_path, "gp_DM_read2_extracted.txt", sep=""))
dm_read3=fread(paste(out_path, "gp_DM_read3_extracted.txt", sep=""))

# in excel I have made a new variable indicating the codes that are likely type 2 codes,
# several (such as 'diabetes mellitus') are unclear if type 1 or type 2 -> check with other data but for now use this file

read2=fread(paste(out_path, "inputs/table_diabetes_diagnosis_read2.txt", sep=""))
read3=fread(paste(out_path, "inputs/table_diabetes_diagnosis_read3.txt", sep=""))

dm_read2 %>%
left_join(read2[,c(1,4)]) -> dm_read2

dm_read3 %>%
left_join(read3[,c(1,4)]) %>%
bind_rows(dm_read2) %>%
distinct()-> dm_read_data

dm_read_data$diabetes_read_type[dm_read_data$read_2_term=="Type 1 diabetes mellitus with ketoacidosis"]<-0
dm_read_data %>%
filter(!is.na(read_3_term) | !is.na(read_2_term)) %>%
mutate(read_3_term = na_if(read_3_term, ""),
    read_2_term = na_if(read_2_term, ""),
    read_term = coalesce(read_3_term,read_2_term)) %>%
select(eid, read_term, diabetes_read_type) -> diabetes_only

### create gp diagnosis based on which code (for which type of diabetes) is most common wihtin an individual
### Make data.table ordered by eid
### Create some summary statistics to explore
diabetes_only_dt <- data.table(diabetes_only, key="eid")
diabetes_likely_pheno <- diabetes_only_dt[, .(mode_type=min(as.numeric(Mode(diabetes_read_type)))), by="eid"]
dat.min <- diabetes_only_dt[, .(min_type=min(diabetes_read_type)), by="eid"]
dat.max <- diabetes_only_dt[, .(max_type=max(diabetes_read_type)), by="eid"]
dat.n <- diabetes_only_dt[, .(n=length(diabetes_read_type)), by="eid"]
dat2count <- diabetes_only_dt[, .(count2=sum(diabetes_read_type == 2)), by="eid"]

diabetes_likely_pheno$min_type <- dat.min$min_type
diabetes_likely_pheno$max_type <- dat.max$max_type
diabetes_likely_pheno$n <- dat.n$n
diabetes_likely_pheno$count2 <- dat2count$count2
diabetes_likely_pheno$prop2 <- diabetes_likely_pheno$count2/diabetes_likely_pheno$n

### Likely diabetes phenotype created:
diabetes_likely_pheno$mode_type[is.na(diabetes_likely_pheno$mode_type)] <- diabetes_likely_pheno$min_type[is.na(diabetes_likely_pheno$mode_type)]
rm(dat.min, dat.max, dat.n, dat2count)

### Then binary variable indicating if T2D is most likely GP diabetes phenotype
diabetes_likely_pheno %>%
mutate(gp_t2d = ifelse(mode_type==0, NA, 
        ifelse(mode_type == 2, 1, 
            ifelse(mode_type==1,0, NA)))) -> diabetes_likely_pheno

### save output
fwrite(dm_read_data, paste(out_path, "data/GP_t2d_read_data.txt", sep=""), col.names=TRUE, row.names=FALSE, quote=FALSE)
fwrite(diabetes_likely_pheno, paste(out_path, "data/GP_t2d_phenotype.txt", sep=""),col.names=TRUE, row.names=FALSE, quote=FALSE)

diabetes_likely_pheno <- fread(paste(out_path, "data/GP_t2d_phenotype.txt", sep=""))

############################################################################
### 6e. GP Medication for T2D
############################################################################
hba1c_tab_v2[, medication_gp:= ifelse(is.na(medication_coded), NA, as.numeric(max(medication_coded, na.rm=T)> 0)), by=eid]

############################################################################
### 6f. Define T2D phenotype using 6a-e
############################################################################
### combine all t2d phenotype files into one to create the actual phenotype

hba1c_tab_v2 %>%
left_join(diabetes_likely_pheno[,c("eid","gp_t2d")]) %>%
left_join(hes_diab_uniqueID) %>%
left_join(self_report_diab) -> hba1c_tab_v3

any(is.na(hba1c_tab_v3$t2d_sr_diagnosis_clean)) # TRUE
any(is.na(hba1c_tab_v3$gp_t2d)) # TRUE

hba1c_tab_v3$t2d_sr_diagnosis_clean[is.na(hba1c_tab_v3$t2d_sr_diagnosis_clean)] <- -999
hba1c_tab_v3$gp_t2d[is.na(hba1c_tab_v3$gp_t2d)] <- -999
hba1c_tab_v3$medication_gp[is.na(hba1c_tab_v3$medication_gp)] <- -999
hba1c_tab_v3$hba1c_over48[is.na(hba1c_tab_v3$hba1c_over48)] <- -999
hba1c_tab_v3$hes_diag_dm[is.na(hba1c_tab_v3$hes_diag_dm)] <- -999
t2d_sr_diagnosis_clean <- hba1c_tab_v3$t2d_sr_diagnosis_clean
gp_t2d <- hba1c_tab_v3$gp_t2d
medication_gp <- hba1c_tab_v3$medication_gp
hba1c_over48 <- hba1c_tab_v3$hba1c_over48
hes_diag_dm <- hba1c_tab_v3$hes_diag_dm
t2d_sr_diagnosis_clean[t2d_sr_diagnosis_clean == -999] <- 0
gp_t2d[gp_t2d == -999] <- 0
medication_gp[medication_gp == -999] <- 0
medication_gp <- as.numeric(medication_gp > 0)
hba1c_over48[hba1c_over48 == -999] <- 0
hes_diag_dm[hes_diag_dm == -999] <- 0

temp <- t2d_sr_diagnosis_clean + gp_t2d + medication_gp + hba1c_over48 + hes_diag_dm

temp2 <- as.numeric(temp > 1)

hba1c_tab_v3$t2d_diagnosis_all_tmp <- temp2
hba1c_tab_v3[, t2d_diagnosis_all_v2:= as.numeric(any(t2d_diagnosis_all_tmp > 0)), by=eid]


############################################################################
### 6g. First date of recorded diabetes
############################################################################
### get list of ids for t2d cases
hba1c_tab_v3 %>%
filter(t2d_diagnosis_all_v2 ==1) %>%
distinct(eid, .keep_all=TRUE) %>%
select(eid, t2d_diagnosis_all_v2) -> ids_t2d_cases
### HES dates:
hesin_diag <- bio_record(project_dir, "hesin_diag")
hesin_main <- bio_record(project_dir, "hesin")
icd9_t2d_hes_code <- c("25000","25010","25020","25090")
icd10_t2d_hes_code <- c("E110","E111","E112","E113","E114","E115","E116","E117","E118","E119")

hes_t2d_cases <- hesin_diag %>%
filter(eid %in% ids_t2d_cases$eid & (diag_icd9 %in% icd9_t2d_hes_code | diag_icd10 %in% icd10_t2d_hes_code)) %>%
mutate(eid_index = paste0(eid,"_",ins_index)) %>%
select(eid, ins_index, eid_index, diag_icd9, diag_icd10) %>%
rename(eid_1 = eid,
ins_index_1 = ins_index) %>% collect()

## to get date of episode in hospital you need the main hesin table.
## the combination of eid and ins_index give a unique identifier for each diagnostic code
## combine them and then extract the date of episode start for each person

hes_t2d_datex_temp <- hesin_main %>%
filter(eid %in% ids_t2d_cases$eid) %>%
mutate(eid_index = paste0(eid,"_",ins_index)) %>%
select(eid, ins_index, eid_index, epistart, epiend) %>% collect()

hes_t2d_dates <- hes_t2d_datex_temp %>%
right_join(hes_t2d_cases, by="eid_index") %>%
mutate(epistart = as.Date(epistart, "%d/%m/%Y"),
    epiend = as.Date(epiend, "%d/%m/%Y"))

rm(hes_t2d_datex_temp)

hes_t2d_dates %>%
arrange(eid, epistart) %>%
group_by(eid) %>%
mutate(first_t2d_date = first(epistart)) %>%
ungroup() %>%
distinct(eid, .keep_all=TRUE) %>%
select(eid, first_t2d_date) %>%
mutate(type_first_date = "hes") -> t2d_hes_first_code

### get t2d first date from GP data
t2d_gp_first_code <-fread(paste(out_path, "data/GP_t2d_read_data.txt", sep=""))
t2d_gp_first_code %>%
mutate(event_dt=as.Date(event_dt, "%d/%m/%Y")) -> t2d_gp_first_code
t2d_gp_first_code <- t2d_gp_first_code[!is.na(t2d_gp_first_code$event_dt), ]
t2d_gp_first_code <- t2d_gp_first_code[t2d_gp_first_code$event_dt > "1902-02-02", ]

t2d_gp_first_code %>%
filter(eid %in% ids_t2d_cases$eid) %>%
arrange(eid, event_dt) %>%
group_by(eid) %>%
mutate(first_t2d_date = first(event_dt)) %>%
ungroup()%>%
distinct(eid, .keep_all=TRUE) %>%
select(eid, first_t2d_date) %>%
mutate(type_first_date = "gp") -> t2d_gp_first_code

t2d_gp_first_code %>%
bind_rows(t2d_hes_first_code) -> t2d_first_code

### get t2d first date based on medication

t2d_meds_first_code<-fread(paste(out_path, "data/gp_diabetes_medication_data.txt", sep=""),header=T)
t2d_meds_first_code <- t2d_meds_first_code %>%
mutate(event_dt=as.Date(event_dt, "%Y-%m-%d")) %>%
filter(eid %in% ids_t2d_cases$eid & event_dt>"1902-02-02") %>%
arrange(eid, event_dt) %>%
group_by(eid) %>%
mutate(first_t2d_date = first(event_dt)) %>%
ungroup()%>%
distinct(eid, .keep_all=TRUE) %>%
select(eid, first_t2d_date) %>%
mutate(type_first_date = "meds") 

t2d_first_code %>%
bind_rows(t2d_meds_first_code) %>%
arrange(eid, first_t2d_date)-> t2d_first_code


### get t2d first date based on hba1c > 48

hba1c_tab_v3_tmp <- hba1c_tab_v3[hba1c_tab_v3$eid %in% ids_t2d_cases$eid, ]
hba1c_tab_v3_tmp <- hba1c_tab_v3_tmp[hba1c_tab_v3_tmp$hba1c > 48, ]
hba1c_tab_v3_tmp <- hba1c_tab_v3_tmp[, min(event_dt), by = eid]
colnames(hba1c_tab_v3_tmp)[2] <- "first_t2d_date"
t2d_hba1c_first_code <- hba1c_tab_v3_tmp %>% mutate(type_first_date = "hba1c_48") 
t2d_hba1c_first_code$first_t2d_date <- as.Date(t2d_hba1c_first_code$first_t2d_date, "%Y-%m-%d")
rm(hba1c_tab_v3_tmp)

t2d_first_code %>%
bind_rows(t2d_hba1c_first_code) %>%
arrange(eid, first_t2d_date)  %>%
full_join(ids_t2d_cases, by="eid") -> t2d_first_code_v2


t2d_first_code_v2 <- data.table(t2d_first_code_v2, key=c("eid", "first_t2d_date"))
t2d_first_code_v3 <- t2d_first_code_v2[, min(first_t2d_date), by=eid]
t2d_first_code_v4 <- t2d_first_code_v2[, type_first_date[first_t2d_date == min(first_t2d_date)], by=eid]

colnames(t2d_first_code_v3)[2] <- "first_t2d_date"
t2d_first_code_v3$type_first_date <- NA
for(i in 1:dim(t2d_first_code_v3)[1]){
  eid_i <- t2d_first_code_v3$eid[i]
  t2d_first_code_v4i <- t2d_first_code_v4[t2d_first_code_v4$eid == eid_i, ]
  if(dim(t2d_first_code_v4i)[1] == 1){
    t2d_first_code_v3$type_first_date[i] <- t2d_first_code_v4i$V1
  }else{
    a <- t2d_first_code_v4i$V1[1]
    for(j in 2:dim(t2d_first_code_v4i)[1]){
      a <- paste(a, t2d_first_code_v4i$V1[j], sep="/")
    }
    t2d_first_code_v3$type_first_date[i] <- a
  }
}
table(t2d_first_code_v3$type_first_date)
t2d_first_code_v3$type_first_date[t2d_first_code_v3$type_first_date == "NA/NA/NA"] <- NA
t2d_first_code_v3$type_first_date[t2d_first_code_v3$type_first_date == "NA/NA/NA/NA"] <- NA

saveRDS(t2d_first_code_v3, file=paste(out_path, "data/t2d_first_occurrence.rds", sep=""))

rm(t2d_first_code_v2, t2d_first_code_v4, t2d_first_code_v4i)

hba1c_tab_v3 %>%
left_join(t2d_first_code_v3) -> hba1c_tab_v3

############################################################################
### 7. T2D case exclusions
############################################################################
### REMOVE THESE PEOPLE
# 1. people who have medication as type first date and >1 medication as these might be in middle of treatment and not original diagnosis

rm_dat <- (as.numeric(hba1c_tab_v3$type_first_date == "meds") + as.numeric(hba1c_tab_v3$type_first_date == "gp/meds") + 
as.numeric(hba1c_tab_v3$type_first_date == "gp/meds/hba1c_48") + as.numeric(hba1c_tab_v3$type_first_date == "meds/hba1c_48"))
rm_dat[is.na(rm_dat)] <- 0
hba1c_tmp <- hba1c_tab_v3[rm_dat == 1, ]

hba1c_tmp[, difftime := difftime(first_t2d_date, event_dt, units = "days")]
hba1c_tmp[ , min_event_dt := min(event_dt), by=eid]
keep <- as.numeric(hba1c_tmp$event_dt == hba1c_tmp$min_event_dt)
hba1c_tmp <- hba1c_tmp[keep==1, ]
hba1c_tmp[, rm_data:= as.numeric(difftime <= 0), by=eid]

ids_remove <- unique(hba1c_tmp$eid[hba1c_tmp$rm_data == 1])

length(ids_remove)

# insulin within first year of diagnosis
hba1c_tab_v3$event_dt <- as.Date(hba1c_tab_v3$event_dt, "%Y-%m-%d")
insulin_dt <- hba1c_tab_v3[insulin == 1, ]
insulin_dt[, min_ins_date:= min(event_dt), by=eid]
insulin_dt <- insulin_dt[insulin_dt$t2d_diagnosis_all_v2 == 1, ]
insulin_dt[, followup:= time_length(interval(as.Date(first_t2d_date), as.Date(min_ins_date)),"years")]
range(insulin_dt$followup, na.rm=T)
insulin_dt <- insulin_dt[insulin_dt$followup <= 1, ]
ids_insulin_rm <- unique(insulin_dt$eid)

## combine ids to remove total of 1165 people
ids_remove <- c(ids_remove, ids_insulin_rm)
ids_remove <- unique(ids_remove)

# remove these people from main dataframe and save data
# add sex

hba1c_tab_v4 <- hba1c_tab_v3[!(hba1c_tab_v3$eid %in% ids_remove), ]

saveRDS(hba1c_tab_v4, paste(out_path, "HbA1c_full_data_medication_coded_t2d_first_occ.rds", sep=""))

covariates <- readRDS(paste(out_path, "covariates.rds", sep="")) %>%
rename("sex"="31-0.0",
    "assessment_centre"="54-0.0",
    "age_at_assessment"="21003-0.0") %>%
select(eid, sex)

hba1c_tab_v4 %>%
left_join(covariates) -> hba1c_tab_v4

hba1c_tab_v4 <- hba1c_tab_v4[!is.na(hba1c_tab_v4$sex), ]
 
saveRDS(hba1c_tab_v4, paste(out_path, "HbA1c_full_data_medication_coded_t2d_first_occ.rds", sep=""))

############################################################################
### 8. Add in MDD diagnosis date
############################################################################
gp_dep<- fread(paste(depression_dir, "dep_ADs_data.txt", sep=""),header=T)

# get list of ids
hba1c_tab_v4 %>%
distinct(eid, .keep_all=TRUE) %>%
filter(depression==1) %>%
select(eid, depression) -> ids_depression_cases

gp_dep %>%
filter(eid %in% ids_depression_cases$eid & type=="depression") %>%
mutate(event_dt = as.Date(event_dt, "%Y-%m-%d")) %>%
arrange(eid, event_dt) %>%
group_by(eid) %>%
mutate(first_date_code = first(event_dt)) %>%
ungroup() %>%
distinct(eid, .keep_all=TRUE) %>%
select(eid, first_date_code) -> dep_first_code

fwrite(dep_first_code, paste(out_path, "data/depression_first_occurrence.txt", sep=""),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")

hba1c_tab_v4 %>%
left_join(dep_first_code) -> hba1c_tab_v4

hba1c_tab_v4 %>% rename(dep_first_date = first_date_code) -> hba1c_tab_v4

############################################################################
### 9. Additional exclusion criteria
############################################################################
## remove people without t2d first date
## remove people with <2 hba1c

hba1c_tab_v5 <- hba1c_tab_v4[hba1c_tab_v4$t2d_diagnosis_all_v2 == 1, ]
hba1c_tab_v5[ , n_hba1c:=length(event_dt), by=eid]
hba1c_tab_v5 <- hba1c_tab_v5[hba1c_tab_v5$n_hba1c > 1, ]
hba1c_tab_v5 <- hba1c_tab_v5[!is.na(hba1c_tab_v5$first_t2d_date), ]


## save final file
saveRDS(hba1c_tab_v5, file=paste(out_path, "data/hba1c_t2d_mdd_medication_clean.rds", sep=""))

hba1c_tab_v5[, age_t2d:= time_length(interval(as.Date(dob), as.Date(first_t2d_date)),"years")]
range(hba1c_tab_v5$age_t2d)

hba1c_tab_v5 <- hba1c_tab_v5[hba1c_tab_v5$age_t2d > 18, ]
length(unique(hba1c_tab_v5$eid)) # 17544 # removed 168 with age at onset <= 18

saveRDS(hba1c_tab_v5, file=paste(out_path, "data/hba1c_t2d_mdd_medication_clean_over18.rds", sep=""))

############################################################################
### 10. Create phenotype and European .fam file for PRS validation analysis
############################################################################
### create phenotype file for PRS validation
t2d_data = readRDS(paste(out_path, "HbA1c_full_data_medication_coded_t2d_first_occ.rds", sep=""))

### Create a dataset wit EID and case-control status:
dt1 <- t2d_data[, unique(t2d_diagnosis_all_v2), by=eid]

### Create a file for the case-control phenotype where missing data for UKB participants who are not cases/ controls in this analysis marked with NA.
### Order of EIDs in this file should match the order of the stored UKB genetic data you hold.
### We did this by using a pre-existing covariate file for out application.
### This file contains: FID, IID, PC1, PC2, PC3, PC4, PC5, PC6, genetic batch (batch) and assessment centre for all UKB participants, and it's rows (EIDs) are ordered according to the genetic data we hold.
### Read this file in:
covar=fread("path_to_the_above_described_covariates/ukb_prj_number_covariates.txt")

t2d_data %>% select(eid, t2d_diagnosis_all_v2) %>%
distinct(eid, .keep_all=TRUE)  -> t2d_pheno_tmp 

### In T2D case-control data remove eid to FID and create IID (== FID)
### Select FID, IID and t2d phenotype columns only
t2d_pheno_tmp %>%
rename(FID = eid) %>%
mutate(IID = FID) %>%
select(FID, IID, t2d_diagnosis_all_v2) -> t2d_pheno_tmp
### covar is according the genetic data so left join our T2D data to this
covar %>% select(FID, IID) %>%
left_join(t2d_pheno_tmp) -> t2d_pheno_PRSice


### The above includes is all ancestries and all degree of genetic relatedness.

### Using QC laid out for individuals here:
### https://opain.github.io/UKB-GenoPrep/
### We read in a list of the rows in our UKB .fam file that are individuals of
### European ancestry.
### GenoPrep is written and maintained by Oliver Pain.

EUR_qc_rows <- read.table(file="/path_to_ukbiobank_QC_dir/UKB.postQC.EUR.keep", header=F)

### Read in your applications fam file
fam_file <- read.table(file="/path_to_your_ukbiobank_genotyped_data/applicationfamfile.fam", header=F)
### Extract European rows of fam file
### This is the European fam file
EUR_qc_fam <- fam_file[EUR_qc_rows$V1, ]

t2d_pheno_PRSice2 <- t2d_pheno_PRSice 
### Set the phenotypes of non-european IDs all equal to NA so they are excluded from the analysis:
t2d_pheno_PRSice2$t2d_diagnosis_all_v2[!(t2d_pheno_PRSice2$IID %in% EUR_qc_fam$V1)] <- NA
table(t2d_pheno_PRSice2$t2d_diagnosis_all_v2)

### European IDs to use within the greed relatedness algorithm to identify related individuals for removal
keep_ukb <- t2d_pheno_PRSice2$IID[!is.na(t2d_pheno_PRSice2$t2d_diagnosis_all_v2)]
### Greedy relatedness software can be downloaded from here:
### https://gitlab.com/choishingwan/GreedyRelated
rel_rm <- bio_gen_related_remove(project_dir=project_dir,
greedy_related="/path_to_software/GreedyRelated-master/greedy_related",
keep = keep_ukb,
seed = 1234)
### Set phenotypes of related IDs to NA so they are excluded.
t2d_pheno_PRSice2$t2d_diagnosis_all_v2[t2d_pheno_PRSice2$IID %in% rel_rm$eid] <- NA

### Write pheno file to output to use in PRSice
write.table(t2d_pheno_PRSice2, paste(analysis_dir, "data/EURt2d_pheno_file_for_PRSice.txt", sep=""),col.names=TRUE, row.names=FALSE, quote=FALSE)

### Write the European post-QC fam file to output
write.table(EUR_qc_fam, file=paste(analysis_dir, "data/EUR_postQC.fam", sep=""),col.names=FALSE, row.names=FALSE, quote=FALSE)

### Note. PRSice needs other files included:
### The summary statistics from the base/ discovery GWAS.
### cov file. Ours has columns: FID, IID, PC1 to PC6, genetic batch, assessment centre.
### extract file. This is a file of SNPs to use. Ours follows the SNP QC outlined in the GenoPrep pipeline (https://opain.github.io/UKB-GenoPrep/)
### target files (your ukbiobank genetic data).
### See PRSice documentation for details.


