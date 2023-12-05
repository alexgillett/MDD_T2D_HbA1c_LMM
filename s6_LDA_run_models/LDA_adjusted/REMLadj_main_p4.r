# User add library path where R packages are stored if required:
#.libPaths(new = "your_rlibrary_path")
### Add R packages
library(data.table)
library(lubridate)
library(ggplot2)
library(nlme)
library(lme4)
library(MASS)
library(rms)
library(Amelia)
library(mitml)
### Assign paths to relevant directories
analysis_dir <- "/path_to_analysis_dir/"

### Load data
load(paste(analysis_dir, "data/LDAhba1c_amelia.RData", sep=""))
### Remove un-needed datasets
rm(count.tab, eids_rm, eids_na_ga, eids_na_srethn, i, spline, dt.s)
### Rename main imputation dataset
hba1c.imp_srethn <- hba1c.imp
### remove original imputation dataset
rm(hba1c.imp, p1)
### Remove individuals with missing self-reported ethnicity data from each
### MI dataset
for(i in 1:50){
    dti <- hba1c.imp_srethn$imputations[[i]]
    dti <- data.table(dti, key=c("eid", "time_base"))
    dti <- dti[!is.na(dti$sr_ethnicity_f), ]
    dti <- data.table(dti, key=c("eid", "time_base"))
    hba1c.imp_srethn$imputations[[i]] <- dti
    rm(dti)
}
### Add in a MI dataset identifier into each dataset
for(i in 1:50){
    dti <- hba1c.imp_srethn$imputations[[i]]
    dti <- data.table(dti, key=c("eid", "time_base"))
    dti$mi_id <- i
    dti <- data.table(dti, key=c("eid", "time_base"))
    hba1c.imp_srethn$imputations[[i]] <- dti
    rm(dti)
}

### Extract the 10 imputed datasets required (fourth set of 10 here)
for(i in 1:10){
  # assign function within loop
  assign(paste0("df_", i), hba1c.imp_srethn$imputations[[i+30]])
}
### Bind them for use with mitml package
df_mi <- rbind(df_1, df_2, df_3, df_4, df_5, df_6, df_7, df_8, df_9, df_10)
rm(df_1, df_2, df_3, df_4, df_5, df_6, df_7, df_8, df_9, df_10)
mi_list4 <- long2mitml.list(df_mi, split = "mi_id", exclude = 0)
rm(df_mi)

### Run, using REML, the (main/ only) adjusted model
fit_adj_mainREML4 <- with(mi_list4,
lme(hba1c_std ~ sex + ac_f + sr_ethnicity_f + age_t2d_diag_std + hba1c_base_std +
ever_smk + alcohol_never + qualifications_f + previsits_tot_std + tdi_std + bmi_std + sbp_base_std + dbp_base_std +
med_base_f + med_f +
time_spline1 + time_spline2 + time_spline3 +
hba1c_ts1 + hba1c_ts2 + hba1c_ts3 +
med_base_ts1 + med_base_ts2 + med_base_ts3 +
med_f1_ts1 + med_f1_ts2 + med_f1_ts3 + med_f2_ts1 + med_f2_ts2 + med_f2_ts3 + med_f3_ts1 + med_f3_ts2 + med_f3_ts3 +
bmi_ts1 + bmi_ts2 + bmi_ts3 +
dep_base + uncentered_robust_pretime2 + dep_change_base + postt2dmdd_time, random = ~ time_base_std|eid_f,
correlation = corCAR1(form=~time_base|eid_f), method="REML"),
include.data=TRUE)

rm(hba1c.imp_srethn, mi_list4)

### save output
save.image(file=paste(analysis_dir, "model_output/adjusted_mainREML4.RData", sep=""))

