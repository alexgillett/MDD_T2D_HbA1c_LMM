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
### Extract a single dataset (start by extracting 10 and then just use df_1,
### in this case (all will be the same for the unadjusted analysis)
for(i in 1:10){
  # assign function within loop
  assign(paste0("df_", i), hba1c.imp_srethn$imputations[[i]])
}
df_mi <- df_1
### Remove additional MI datasets created in loop above
rm(df_1, df_2, df_3, df_4, df_5, df_6, df_7, df_8, df_9, df_10)

### Run 'Main' model:
# This is the time splines, pre-T2D MDD (dep_base, = binary),
# A centred version of pre-T2D MDD duration up to baseline/ index
# (uncentered_robust_pretime2, semi-continuous variable),
# time-varying post-T2D MDD (dep_change_base),
# and time since post-T2D MDD (postspline1)
fit_unadj_mainREML <- lme(hba1c_std ~ time_spline1 + time_spline2 + time_spline3 +
dep_base + uncentered_robust_pretime2 + dep_change_base + postspline1,
random = ~ time_base_std|eid_f,
correlation = corCAR1(form=~time_base|eid_f), method="REML", data=df_mi)

### Only need to run the 'final' fixed effects model for REML (to get parameter estimates)

rm(hba1c.imp_srethn)
save.image(file=paste(analysis_dir, "model_output/unadj_REML.RData", sep=""))
