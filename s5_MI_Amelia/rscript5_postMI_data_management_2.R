############################################################################
### 5.2. Data management script to run AFTER Amelia (MI)
### Preps dataset(s) for LDA
############################################################################
### Authors: Alexandra C Gillett
############################################################################
### In script we:
### 1. Load in data from MI (Amelia)
### 2. Add columns to imputed and original datasets
### 3. Add in some variables to the imputated datasets from the original
### 4. Remove measurements occuring prior to T2D diagnosis
### 5. Add in additional depression variables
### 6. Variable prep for linear mixed effects models
# includes removing those with missing self-reported ethnicity, and
# saving data
############################################################################
# User add library path where R packages are stored if required:
#.libPaths(new = "your_rlibrary_path")
############################################################################
### Add R packages:
library(data.table)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(nlme)
library(MASS)
library(Amelia)
library(AICcmodavg)
library(lmeInfo)
library(msm)
library(mitml)
library(rms)
### Assign paths to relevant directories
hba1c_dat_path <- "/path_to_where_you_store_extracted_data/hba1c_data/"
analysis_dir <- "/path_to_analysis_dir/"
############################################################################
### 1. Load in data from MI (Amelia)
############################################################################
load(paste(analysis_dir, "data/hba1c_amelia50.RData", sep=""))
############################################################################
### 2. Add columns to imputed and original datasets:
### imputed need: hba1c_std, hba1c_base and hba1c_base_std
### original needs: hba1c_base_std
############################################################################
### NOTE: It is likely that I would be able to do this with output from Amelia still in it's list form, but I'm not sure...
### So am going to split the data into different datasets, and work from there.
rm(df.imp, dt.imp, dt.imp2, df.imp2, bound.mat2, dat.base, data.comp.dens.t0,data.comp.dens.t0std, hba1c_observed, mean_hba1c_imp2, mean_hba1c_imp, s)
### Read in the standardisation parameters:
std_params <- readRDS(paste(hba1c_dat_path, "data/cont_vars_mean_SD_baseline.rds", sep=""))
### Create hba1c_std for each imputation, then create hba1c_base_std variable (and hba1c_base)
for(i in 1:50){
    hba1c.imp$imputations[[i]] <- data.table(hba1c.imp$imputations[[i]], key=c("eid", "time_base"))
    dti <- hba1c.imp$imputations[[i]]
    dti$hba1c_std <- (dti$hba1c - std_params[1,1])/std_params[1,2]
    dti0 <- data.table(dti[time_base == 0, ], key="eid")
    dti0 <- dti0[, c("eid", "hba1c", "hba1c_std")]
    dti0 <- dti0 %>% rename(hba1c_base = hba1c, hba1c_base_std = hba1c_std)
    dti0 <- data.table(dti0, key="eid")
    dti <- dti %>% left_join(dti0)
    dti <- data.table(dti, key=c("eid", "time_base"))
    hba1c.imp$imputations[[i]] <- dti
}
### dt.s has hba1c_std but no baseline value
### it does have hba1c_base though
dt.s$hba1c_base_std <- (dt.s$hba1c_base - std_params[1,1])/std_params[1,2]
dt.s <- data.table(dt.s, key=c("eid", "time_base"))
############################################################################
### 3. Add in some variables to the imputated datasets from the original
############################################################################
table(dt.s$time_base - hba1c.imp$imputations[[1]]$time_base) ### all 0 as expected -> order is same...
for(i in 1:50){
    hba1c.imp$imputations[[i]] <- data.table(hba1c.imp$imputations[[i]], key=c("eid", "time_base"))
    dti <- hba1c.imp$imputations[[i]]
    dti$pret2dmdd_time <- dt.s$pret2dmdd_time
    dti$postt2dmdd_time <- dt.s$postt2dmdd_time
    dti$pret2dmdd_sg10to25pc <- dt.s$pret2dmdd_sg10to25pc
    dti <- data.table(dti, key=c("eid", "time_base"))
    hba1c.imp$imputations[[i]] <- dti
}
############################################################################
### 4. Remove measurements occuring prior to T2D diagnosis
### This is to create LDA datasets
############################################################################
dt.s <- dt.s[time_base > 0, ]
dim(dt.s)
length(unique(dt.s$eid))
### Count the number of rows per ID:
count.tab <- dt.s[, .(rowCount = .N), by=eid]
table(count.tab$rowCount) ### There are individuals with only 1 post T2D row - remove
eids_rm <- count.tab$eid[count.tab$rowCount == 1]
dt.s <- dt.s[!(eid %in% eids_rm), ]
dt.s <- data.table(dt.s, key=c("eid", "time_base"))

for(i in 1:50){
    dti <- hba1c.imp$imputations[[i]]
    dti <- data.table(dti, key=c("eid", "time_base"))
    dti <- dti[time_base > 0, ]
    dti <- dti[!(eid %in% eids_rm), ]
    hba1c.imp$imputations[[i]] <- dti
    rm(dti)
}
rm(dti0, i)
############################################################################
### 5. Add in additional depression variables
### dep_time_diff: time between T2D and MDD diagnoses
### dep_post10years: individual newly diagnosed with MDD during the 10 years
### of follow-up
### dep_group: pre, post, no MDD (called control in code)
############################################################################
### Couple of bonus columns that should be in dt.s for LDA...
dt.s <- dt.s %>% mutate(dep_time_diff = time_length(interval(as.Date(first_t2d_date), as.Date(dep_first_date)),"years"))
dt.s <- dt.s %>% mutate(dep_post10years = ifelse(depression == 0, 0, ifelse(dep_time_diff > 0 & dep_time_diff <= 10, 1, 0)))
dt.s <- data.table(dt.s, key=c("eid", "time_base"))
dt.s <- dt.s %>% 
    mutate(dep_group = ifelse(dep_post10years == 1, "t2d_dep", ifelse(dep_base == 1, "dep_t2d", "control")))
depgrouptab <- dt.s[, unique(dep_group), by=eid]
colnames(depgrouptab)[2] <- "dep_group"
for(i in 1:50){
    dti <- hba1c.imp$imputations[[i]]
    dti <- data.table(dti, key=c("eid", "time_base"))
    dti <- dti %>% left_join(depgrouptab)
    dti <- data.table(dti, key=c("eid", "time_base"))
    hba1c.imp$imputations[[i]] <- dti
    rm(dti)
}
rm(depgrouptab)
############################################################################
### 6. Variable prep for linear mixed effects models
### Includes removing those missing self-reported ethnicity
############################################################################
### Add in splines, and standardised time_base variable...
#table(dt.s$eid - hba1c.imp$imputations[[1]]$eid)
#table(dt.s$time_base - hba1c.imp$imputations[[1]]$time_base)
### order is therefore the same across dt.s and imputation datasets...
### using 4 knots due to interactions tested and sample size
spline <- rcs(dt.s$time_base, 4)
dt.s$time_spline1 <- as.vector(spline[,1])
dt.s$time_spline2 <- as.vector(spline[,2])
dt.s$time_spline3 <- as.vector(spline[,3])
time_mean <- mean(dt.s$time_base) ### 4.059199
time_sd <- sd(dt.s$time_base) ### 2.727724
dt.s$time_base_std <- (dt.s$time_base -time_mean)/time_sd
dt.s <- data.table(dt.s, key=c("eid", "time_base"))
for(i in 1:50){
    dti <- hba1c.imp$imputations[[i]]
    dti <- data.table(dti, key=c("eid", "time_base"))
    dti$time_spline1 <- dt.s$time_spline1
    dti$time_spline2 <- dt.s$time_spline2
    dti$time_spline3 <- dt.s$time_spline3
    dti$time_base_std <- (dti$time_base -time_mean)/time_sd
    dti <- data.table(dti, key=c("eid", "time_base"))
    hba1c.imp$imputations[[i]] <- dti
    rm(dti)
}
### The main aim of the imputation algorithm was to impute HbA1c at diagnosis.
### The aim was not to impute ethnicity or genetic ancestry, etc- which I do not think
### could be imputed from the info provided tbh...
### Introduce NAs back in for these at least...
eids_na_ga <- unique(dt.s$eid[is.na(dt.s$pop)])
length(eids_na_ga)
eids_na_srethn <- unique(dt.s$eid[is.na(dt.s$sr_ethnicity_f)])
length(eids_na_srethn)
for(i in 1:50){
    dti <- hba1c.imp$imputations[[i]]
    dti <- data.table(dti, key=c("eid", "time_base"))
    dti$pop_f[dti$eid %in% eids_na_ga] <- NA
    dti$sr_ethnicity_f[dti$eid %in% eids_na_srethn] <- NA
    dti <- data.table(dti, key=c("eid", "time_base"))
    hba1c.imp$imputations[[i]] <- dti
    rm(dti)
}

### Make factors and other preps...
dt.s <- data.table(dt.s, key=c("eid", "time_base"))
dt.s$med_f1 <- as.numeric(dt.s$medication_coded == 1)
dt.s$med_f2 <- as.numeric(dt.s$medication_coded == 2)
dt.s$med_f3 <- as.numeric(dt.s$medication_coded == 3)

dt.s$med_f1_ts1 <- as.numeric(dt.s$medication_coded == 1)*dt.s$time_spline1
dt.s$med_f1_ts2 <- as.numeric(dt.s$medication_coded == 1)*dt.s$time_spline2
dt.s$med_f1_ts3 <- as.numeric(dt.s$medication_coded == 1)*dt.s$time_spline3

dt.s$med_f2_ts1 <- as.numeric(dt.s$medication_coded == 2)*dt.s$time_spline1
dt.s$med_f2_ts2 <- as.numeric(dt.s$medication_coded == 2)*dt.s$time_spline2
dt.s$med_f2_ts3 <- as.numeric(dt.s$medication_coded == 2)*dt.s$time_spline3

dt.s$med_f3_ts1 <- as.numeric(dt.s$medication_coded == 3)*dt.s$time_spline1
dt.s$med_f3_ts2 <- as.numeric(dt.s$medication_coded == 3)*dt.s$time_spline2
dt.s$med_f3_ts3 <- as.numeric(dt.s$medication_coded == 3)*dt.s$time_spline3

dt.s$med_base_ts1 <- as.numeric(dt.s$medication_coded_base == 1)*dt.s$time_spline1
dt.s$med_base_ts2 <- as.numeric(dt.s$medication_coded_base == 1)*dt.s$time_spline2
dt.s$med_base_ts3 <- as.numeric(dt.s$medication_coded_base == 1)*dt.s$time_spline3

dt.s$dep_base_ts1 <- dt.s$dep_base*dt.s$time_spline1
dt.s$dep_base_ts2 <- dt.s$dep_base*dt.s$time_spline2
dt.s$dep_base_ts3 <- dt.s$dep_base*dt.s$time_spline3

dt.s$dep_change_ts1 <- dt.s$dep_change_base*dt.s$time_spline1
dt.s$dep_change_ts2 <- dt.s$dep_change_base*dt.s$time_spline2
dt.s$dep_change_ts3 <- dt.s$dep_change_base*dt.s$time_spline3

dt.s$dep_t2d <- as.numeric(dt.s$dep_group == "dep_t2d")
dt.s$t2d_dep <- as.numeric(dt.s$dep_group == "t2d_dep")  

dt.s$dep_t2d_ts1 <- dt.s$dep_t2d*dt.s$time_spline1
dt.s$dep_t2d_ts2 <- dt.s$dep_t2d*dt.s$time_spline2
dt.s$dep_t2d_ts3 <- dt.s$dep_t2d*dt.s$time_spline3

dt.s$t2d_dep_ts1 <- dt.s$t2d_dep*dt.s$time_spline1
dt.s$t2d_dep_ts2 <- dt.s$t2d_dep*dt.s$time_spline2
dt.s$t2d_dep_ts3 <- dt.s$t2d_dep*dt.s$time_spline3

### uncentered_robust_pretime2 interacting with time splines...
dt.s$robust_pre_ts1 <- dt.s$uncentered_robust_pretime2*dt.s$time_spline1
dt.s$robust_pre_ts2 <- dt.s$uncentered_robust_pretime2*dt.s$time_spline2
dt.s$robust_pre_ts3 <- dt.s$uncentered_robust_pretime2*dt.s$time_spline3

### Subgroups interacting with time_splines...
dt.s$pret2dmdd_sgo75pc_ts1 <- dt.s$pret2dmdd_sgo75pc*dt.s$time_spline1
dt.s$pret2dmdd_sgo75pc_ts2 <- dt.s$pret2dmdd_sgo75pc*dt.s$time_spline2
dt.s$pret2dmdd_sgo75pc_ts3 <- dt.s$pret2dmdd_sgo75pc*dt.s$time_spline3

dt.s$pret2dmdd_sg25to75pc_ts1 <- dt.s$pret2dmdd_sg25to75pc*dt.s$time_spline1
dt.s$pret2dmdd_sg25to75pc_ts2 <- dt.s$pret2dmdd_sg25to75pc*dt.s$time_spline2
dt.s$pret2dmdd_sg25to75pc_ts3 <- dt.s$pret2dmdd_sg25to75pc*dt.s$time_spline3

dt.s$pret2dmdd_sg10to25pc_ts1 <- dt.s$pret2dmdd_sg10to25pc*dt.s$time_spline1
dt.s$pret2dmdd_sg10to25pc_ts2 <- dt.s$pret2dmdd_sg10to25pc*dt.s$time_spline2
dt.s$pret2dmdd_sg10to25pc_ts3 <- dt.s$pret2dmdd_sg10to25pc*dt.s$time_spline3

dt.s$pret2dmdd_sgl10pc_ts1 <- dt.s$pret2dmdd_sgl10pc*dt.s$time_spline1
dt.s$pret2dmdd_sgl10pc_ts2 <- dt.s$pret2dmdd_sgl10pc*dt.s$time_spline2
dt.s$pret2dmdd_sgl10pc_ts3 <- dt.s$pret2dmdd_sgl10pc*dt.s$time_spline3

dt.s <- data.table(dt.s, key=c("eid", "time_base"))
for(i in 1:50){
    dti <- hba1c.imp$imputations[[i]]
    dti <- data.table(dti, key=c("eid", "time_base"))
    dti$eid_f <- as.factor(dti$eid)
    dti$dep_group_f <- as.factor(dti$dep_group)
    dti$med_f <- as.factor(dti$med_f)
    ### medication and time splines interactions...
    dti$med_f1_ts1 <- dt.s$med_f1_ts1
    dti$med_f1_ts2 <- dt.s$med_f1_ts2
    dti$med_f1_ts3 <- dt.s$med_f1_ts3
    
    dti$med_f2_ts1 <- dt.s$med_f2_ts1
    dti$med_f2_ts2 <- dt.s$med_f2_ts2
    dti$med_f2_ts3 <- dt.s$med_f2_ts3
    
    dti$med_f3_ts1 <- dt.s$med_f3_ts1
    dti$med_f3_ts2 <- dt.s$med_f3_ts2
    dti$med_f3_ts3 <- dt.s$med_f3_ts3

    dti$med_base_ts1 <- dt.s$med_base_ts1
    dti$med_base_ts2 <- dt.s$med_base_ts2
    dti$med_base_ts3 <- dt.s$med_base_ts3
    ### depression with splines...
    dti$dep_base_ts1 <- dt.s$dep_base_ts1
    dti$dep_base_ts2 <- dt.s$dep_base_ts2
    dti$dep_base_ts3 <- dt.s$dep_base_ts3
    
    dti$dep_change_ts1 <- dt.s$dep_change_ts1
    dti$dep_change_ts2 <- dt.s$dep_change_ts2
    dti$dep_change_ts3 <- dt.s$dep_change_ts3
    ### broad depression groups with splines...
    dti$dep_t2d <- as.numeric(dti$dep_group == "dep_t2d")
    dti$t2d_dep <- as.numeric(dti$dep_group == "t2d_dep") 

    dti$dep_t2d_ts1 <- dti$dep_t2d*dti$time_spline1
    dti$dep_t2d_ts2 <- dti$dep_t2d*dti$time_spline2
    dti$dep_t2d_ts3 <- dti$dep_t2d*dti$time_spline3

    dti$t2d_dep_ts1 <- dti$t2d_dep*dti$time_spline1
    dti$t2d_dep_ts2 <- dti$t2d_dep*dti$time_spline2
    dti$t2d_dep_ts3 <- dti$t2d_dep*dti$time_spline3

    ### pret2dmdd time variable with splines...
    dti$robust_pre_ts1 <- dt.s$robust_pre_ts1
    dti$robust_pre_ts2 <- dt.s$robust_pre_ts2
    dti$robust_pre_ts3 <- dt.s$robust_pre_ts3

    ### Subgroups interacting with time_splines...
    dti$pret2dmdd_sgo75pc_ts1 <- dt.s$pret2dmdd_sgo75pc_ts1
    dti$pret2dmdd_sgo75pc_ts2 <- dt.s$pret2dmdd_sgo75pc_ts2
    dti$pret2dmdd_sgo75pc_ts3 <- dt.s$pret2dmdd_sgo75pc_ts3
    
    dti$pret2dmdd_sg25to75pc_ts1 <- dt.s$pret2dmdd_sg25to75pc_ts1
    dti$pret2dmdd_sg25to75pc_ts2 <- dt.s$pret2dmdd_sg25to75pc_ts2
    dti$pret2dmdd_sg25to75pc_ts3 <- dt.s$pret2dmdd_sg25to75pc_ts3
    
    dti$pret2dmdd_sg10to25pc_ts1 <- dt.s$pret2dmdd_sg10to25pc_ts1
    dti$pret2dmdd_sg10to25pc_ts2 <- dt.s$pret2dmdd_sg10to25pc_ts2
    dti$pret2dmdd_sg10to25pc_ts3 <- dt.s$pret2dmdd_sg10to25pc_ts3
    
    dti$pret2dmdd_sgl10pc_ts1 <- dt.s$pret2dmdd_sgl10pc_ts1
    dti$pret2dmdd_sgl10pc_ts2 <- dt.s$pret2dmdd_sgl10pc_ts2
    dti$pret2dmdd_sgl10pc_ts3 <- dt.s$pret2dmdd_sgl10pc_ts3

    dti <- data.table(dti, key=c("eid", "time_base"))
    hba1c.imp$imputations[[i]] <- dti
    rm(dti)
}
### baseline hba1c standardised x time_splines...
for(i in 1:50){
    dti <- hba1c.imp$imputations[[i]]
    dti <- data.table(dti, key=c("eid", "time_base"))
    dti$hba1c_ts1 <- dti$hba1c_base_std*dti$time_spline1
    dti$hba1c_ts2 <- dti$hba1c_base_std*dti$time_spline2
    dti$hba1c_ts3 <- dti$hba1c_base_std*dti$time_spline3
    dti <- data.table(dti, key=c("eid", "time_base"))
    hba1c.imp$imputations[[i]] <- dti
    rm(dti)
}
hba1c.imp$imputations[[1]]
dt.s$hba1c_ts1 <- dt.s$hba1c_base_std*dt.s$time_spline1
dt.s$hba1c_ts2 <- dt.s$hba1c_base_std*dt.s$time_spline2
dt.s$hba1c_ts3 <- dt.s$hba1c_base_std*dt.s$time_spline3
### Might need depression tv and baseline as factors...
dt.s$dep_change_base_f <- as.factor(dt.s$dep_change_base)
dt.s$dep_base_f <- as.factor(dt.s$dep_base)
for(i in 1:50){
    dti <- hba1c.imp$imputations[[i]]
    dti <- data.table(dti, key=c("eid", "time_base"))
    dti$dep_change_base_f <- as.factor(dti$dep_change_base)
    dti$dep_base_f <- as.factor(dti$dep_base)
    dti <- data.table(dti, key=c("eid", "time_base"))
    hba1c.imp$imputations[[i]] <- dti
    rm(dti)   
}
### BMI, SBP and DBP time-spline interactions...
for(i in 1:50){
    dti <- hba1c.imp$imputations[[i]]
    dti <- data.table(dti, key=c("eid", "time_base"))
    ### bmi
    dti$bmi_ts1 <- dti$bmi_std*dti$time_spline1
    dti$bmi_ts2 <- dti$bmi_std*dti$time_spline2
    dti$bmi_ts3 <- dti$bmi_std*dti$time_spline3
    ### sbp
    dti$sbp_ts1 <- dti$sbp_base_std*dti$time_spline1
    dti$sbp_ts2 <- dti$sbp_base_std*dti$time_spline2
    dti$sbp_ts3 <- dti$sbp_base_std*dti$time_spline3
    ### dbp
    dti$dbp_ts1 <- dti$dbp_base_std*dti$time_spline1
    dti$dbp_ts2 <- dti$dbp_base_std*dti$time_spline2
    dti$dbp_ts3 <- dti$dbp_base_std*dti$time_spline3
    ### tdi
    dti$tdi_ts1 <- dti$tdi_std*dti$time_spline1
    dti$tdi_ts2 <- dti$tdi_std*dti$time_spline2
    dti$tdi_ts3 <- dti$tdi_std*dti$time_spline3

    dti <- data.table(dti, key=c("eid", "time_base"))
    hba1c.imp$imputations[[i]] <- dti
    rm(dti)  
}
### bmi
dt.s$bmi_ts1 <- dt.s$bmi_std*dt.s$time_spline1
dt.s$bmi_ts2 <- dt.s$bmi_std*dt.s$time_spline2
dt.s$bmi_ts3 <- dt.s$bmi_std*dt.s$time_spline3
### sbp
dt.s$sbp_ts1 <- dt.s$sbp_base_std*dt.s$time_spline1
dt.s$sbp_ts2 <- dt.s$sbp_base_std*dt.s$time_spline2
dt.s$sbp_ts3 <- dt.s$sbp_base_std*dt.s$time_spline3
### dbp
dt.s$dbp_ts1 <- dt.s$dbp_base_std*dt.s$time_spline1
dt.s$dbp_ts2 <- dt.s$dbp_base_std*dt.s$time_spline2 
dt.s$dbp_ts3 <- dt.s$dbp_base_std*dt.s$time_spline3
### tdi
dt.s$tdi_ts1 <- dt.s$tdi_std*dt.s$time_spline1
dt.s$tdi_ts2 <- dt.s$tdi_std*dt.s$time_spline2 
dt.s$tdi_ts3 <- dt.s$tdi_std*dt.s$time_spline3
dt.s <- data.table(dt.s, key=c("eid", "time_base"))

### Save all created in this R space
save.image(paste(analysis_dir, "data/LDAhba1c_amelia.RData", sep=""))
### Save the original unimputed dataset separately
saveRDS(dt.s, file=paste(analysis_dir, "data/LDA_dataset.rds", sep=""))
