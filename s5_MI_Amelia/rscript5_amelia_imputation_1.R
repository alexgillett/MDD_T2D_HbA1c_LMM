############################################################################
### 5.1. AMELIA MULTIPLE IMPUTATION SCRIPT
############################################################################
### Authors: Alexandra C Gillett
############################################################################
### In script we:
### 1. Create dataset with just the variables required for imputation
### 2. Indentify bounds for HbA1c at T2D diagnosis
### 3. Test run of imputation models
### 4. Plots for test imputation models
### 5. Imputation model to generate data to use in LDA analysis 2 and 3
### 6. Plots for imputation mdoel
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
### Assign paths to relevant directories
hba1c_dat_path <- "/path_where_you_to_store_extracted_data/hba1c_data/"
analysis_dir <- "/path_to_analysis_dir/"
############################################################################
### 1. Create dataset with just the variables required for imputation
############################################################################
### Read in dataset to use in imputation script
dt.s <- readRDS(paste(hba1c_dat_path, "data/t2d_imputeprep_ds.rds", sep=""))
############################################################################
### What columns to include in imputation...
### eid, hba1c, type_first_date (medication_gp, gp_t2d, hes_diag_dm, hba1c_over48, etc), medication_coded, med at baseline,
### sex, age_t2d_diag, yo_diag, time_base, assessment_centre, tdi, ever_smk, alcohol_never, sr_ethnicity_f, 
### dep_change_base, qualifications_f, pop_f, 
### obs_pre, obs_post, bmi, sbp_base, dbp_base, sbp_max, mdbp, nobs_pre_DBP, followup, previsits_tot_std,
### uncentered_robust_pretime2, 
### postspline1 (postt2dmdd_time), postspline2, postspline3
### Need to check the above for correlation though as variances are non-invertible when there are 
### highly co-linear variables.
### Cannot have dep_base_sg_f and dep_base in the model. Note, dep_base  == one of the factor levels.
### Just have dep_base_sg_f as this includes dep_base
############################################################################
### Ensure factor variables are factors
dt.s$ac_f <- as.factor(dt.s$assessment_centre)
dt.s$med_f <- as.factor(dt.s$medication_coded)
dt.s$pop_f <- as.factor(dt.s$pop)
dt.s$sr_ethnicity_f <- as.factor(dt.s$sr_ethnicity_f)
dt.s$qualifications_f <- as.factor(dt.s$qualifications)
dt.s$med_base_f <- as.factor(dt.s$medication_coded_base)
dt.s <- data.table(dt.s, key=c("eid", "event_dt"))
dt.s$obs_pre <- as.numeric(dt.s$obs_pre)
### Tidy T2D 'type' variables
dt.s$hba1c_over48[is.na(dt.s$hba1c_over48)] <- 0
dt.s$medication_gp[is.na(dt.s$medication_gp)] <- 0
dt.s$medication_gp[dt.s$medication_gp == -999] <- 0
dt.s$gp_t2d[is.na(dt.s$gp_t2d)] <- 0
dt.s$gp_t2d[dt.s$gp_t2d == -999] <- 0
dt.s$hes_diag_dm[is.na(dt.s$hes_diag_dm)] <- 0
dt.s$hes_diag_dm[dt.s$hes_diag_dm == -999] <- 0
dt.s$t2d_sr_diagnosis_clean[is.na(dt.s$t2d_sr_diagnosis_clean)] <- 0
dt.s$t2d_sr_diagnosis_clean[dt.s$t2d_sr_diagnosis_clean == -999] <- 0
### Make pre-T2D MDD duration subgroup binary variables
# This is so that we can exclude one factor level to avoid
# multicollinearity with the pre-T2D MDD binary variable
dt.s$pret2dmdd_sgo75pc <- as.numeric(dt.s$dep_base_sg2 == 1)
dt.s$pret2dmdd_sg25to75pc <- as.numeric(dt.s$dep_base_sg2 == 2)
dt.s$pret2dmdd_sg10to25pc <- as.numeric(dt.s$dep_base_sg2 == 3)
dt.s$pret2dmdd_sgl10pc <- as.numeric(dt.s$dep_base_sg2 == 4)
### Create imputation dataset...
dt.imp <- dt.s[, .(eid, time_base, hba1c_std, med_base_f, med_f, sex, ac_f, age_t2d_diag_std, yo_diag_std, tdi_std,
   dep_change_base, pop_f, sr_ethnicity_f, bmi_std, bmi_ukb_std, sbp_base_std, dbp_base_std, sbp_max_std, mdbp_std, ever_smk, alcohol_never, qualifications_f,
   obs_pre, obs_post, nobs_pre_DBP, previsits_tot_std, followup, hba1c_over48, medication_gp, gp_t2d, hes_diag_dm, t2d_sr_diagnosis_clean,
   uncentered_robust_pretime2, postspline1, postspline2, postspline3, dep_base, pret2dmdd_sgo75pc, pret2dmdd_sg25to75pc, pret2dmdd_sgl10pc)]
### Make it into a data.frame (required by Amelia)
### Make an equivalent version with unstandardised HbA1c to use in log-Normal
### imputation model
df.imp <- data.frame(dt.imp)
dt.imp2 <- dt.s[, .(eid, time_base, hba1c, med_base_f, med_f, sex, ac_f, age_t2d_diag_std, yo_diag_std, tdi_std,
   dep_change_base, pop_f, sr_ethnicity_f, bmi_std, bmi_ukb_std, sbp_base_std, dbp_base_std, sbp_max_std, mdbp_std, ever_smk, alcohol_never, qualifications_f,
   obs_pre, obs_post, nobs_pre_DBP, previsits_tot_std, followup, hba1c_over48, medication_gp, gp_t2d, hes_diag_dm, t2d_sr_diagnosis_clean,
   uncentered_robust_pretime2, postspline1, postspline2, postspline3, dep_base, pret2dmdd_sgo75pc, pret2dmdd_sg25to75pc, pret2dmdd_sgl10pc)]
df.imp2 <- data.frame(dt.imp2)

############################################################################
### 2. Indentify bounds for HbA1c at T2D diagnosis
############################################################################
# Using a bound on HbA1c to accomodate that we are only imputing t=0 values
# which do not follow a standard normal distribution
############################################################################
dat.base <- dt.s[dt.s$time_base ==0, ]
hba1c_observed <- dat.base$hba1c[!is.na(dat.base$hba1c)]
range(hba1c_observed) ### 30 184
quantile(hba1c_observed, probs = c(0, 0.005, 0.01, 0.02))

quantile(hba1c_observed, probs = c(0, 0.005, 0.01, 0.012, 0.013, 0.0131, 0.015))

unique(dt.s$hba1c_std[dt.s$hba1c == 43])
### Corresponds to a standardised HbA1c of -0.9865925
unique(dt.s$hba1c_std[dt.s$hba1c == 184])
### 44.5: approx 1.3%
### Corresponds to a standardised HbA1c of 7.121568
bound.mat0 <- t(matrix(c(3, -0.9003352, 7.121568)))
bound.mat <- t(matrix(c(3, -0.9865925, 7.121568)))
bound.mat2 <- t(matrix(c(3, log(44.5), log(184))))
bound.mat3 <- t(matrix(c(3, log(43), log(184))))
############################################################################
### 3. Test run of imputation models
############################################################################
### a) Run a test with HbA1c ~ Normal (hba1c.impnormal)
### b) Run a test with HbA1c ~ doubly truncated Normal 
# with lower trunc == 43 (hba1c.imp)
### c) Run a test with HbA1c ~ doubly truncated Normal
# with lower trunc == 44.5 (hba1c.imp44_5)
### d) Run a test with HbA1c ~ doubly truncated Log-Normal
# with lower trunc == log(43) (hba1c.implog_43)
### e) Run a test with HbA1c ~ doubly truncated Log-Normal
# with lower trunc == log(44.5) (hba1c.implog)
# Note- because assessment_centre has so many categories/ levels Amelia will
# give you a warning message to tell you this. You can ignore this.
s <- 665
set.seed(s) ### setting a seed in case need to reproduce.
hba1c.impnormal <- amelia(df.imp, m=10, ts="time_base", cs="eid", splinetime=3, lags="hba1c_std", leads="hba1c_std", ords=c("med_f"),
    noms=c("med_base_f", "ac_f", "sex", "dep_change_base", "pop_f", "sr_ethnicity_f",
    "ever_smk", "alcohol_never", "qualifications_f", "hba1c_over48", "medication_gp", "gp_t2d", "hes_diag_dm",
    "t2d_sr_diagnosis_clean", "dep_base", "pret2dmdd_sgo75pc", "pret2dmdd_sg25to75pc", "pret2dmdd_sgl10pc"),
    p2s = 1)
s <- 666
set.seed(s)
hba1c.imp <- amelia(df.imp, m=10, ts="time_base", cs="eid", splinetime=3, lags="hba1c_std", leads="hba1c_std", ords=c("med_f"),
    noms=c("med_base_f", "ac_f", "sex", "dep_change_base", "pop_f", "sr_ethnicity_f",
    "ever_smk", "alcohol_never", "qualifications_f", "hba1c_over48", "medication_gp", "gp_t2d", "hes_diag_dm", 
    "t2d_sr_diagnosis_clean", "dep_base", "pret2dmdd_sgo75pc", "pret2dmdd_sg25to75pc", "pret2dmdd_sgl10pc"), 
    p2s = 1, bounds=bound.mat)
s <- 666
set.seed(s)
hba1c.imp44_5 <- amelia(df.imp, m=10, ts="time_base", cs="eid", splinetime=3, lags="hba1c_std", leads="hba1c_std", ords=c("med_f"),
    noms=c("med_base_f", "ac_f", "sex", "dep_change_base", "pop_f", "sr_ethnicity_f",
    "ever_smk", "alcohol_never", "qualifications_f", "hba1c_over48", "medication_gp", "gp_t2d", "hes_diag_dm",
    "t2d_sr_diagnosis_clean", "dep_base", "pret2dmdd_sgo75pc", "pret2dmdd_sg25to75pc", "pret2dmdd_sgl10pc"),
    p2s = 1, bounds=bound.mat0)

set.seed(s+1)
hba1c.implog <- amelia(df.imp2, m=10, ts="time_base", cs="eid", splinetime=3, lags="hba1c", leads="hba1c", ords=c("med_f"), 
    noms=c("med_base_f", "ac_f", "sex", "dep_change_base", "pop_f", "sr_ethnicity_f",
    "ever_smk", "alcohol_never", "qualifications_f", "hba1c_over48", "medication_gp", "gp_t2d", "hes_diag_dm", 
    "t2d_sr_diagnosis_clean", "dep_base", "pret2dmdd_sgo75pc", "pret2dmdd_sg25to75pc", "pret2dmdd_sgl10pc"), 
    p2s = 1, logs=c("hba1c"), bounds=bound.mat2)

set.seed(s+1)
hba1c.implog_43 <- amelia(df.imp2, m=10, ts="time_base", cs="eid", splinetime=3, lags="hba1c", leads="hba1c", ords=c("med_f"),
    noms=c("med_base_f", "ac_f", "sex", "dep_change_base", "pop_f", "sr_ethnicity_f",
    "ever_smk", "alcohol_never", "qualifications_f", "hba1c_over48", "medication_gp", "gp_t2d", "hes_diag_dm",
    "t2d_sr_diagnosis_clean", "dep_base", "pret2dmdd_sgo75pc", "pret2dmdd_sg25to75pc", "pret2dmdd_sgl10pc"),
    p2s = 1, logs=c("hba1c"), bounds=bound.mat3)
############################################################################
### 4. Plots for test imputation models
############################################################################
### Density plots for missing baseline data.
############################################################################
### Functions to generate output for plots of HbA1c at baseline/ index
data.comp.dens.t0 <- function(amelia.dat, orig.dat, m=10){
    bin.na <- as.numeric(is.na(orig.dat$hba1c))
    n_rows <- sum(bin.na)
    mat.work <- matrix(0, nrow=n_rows, ncol=m)
    for(i in 1:m){
        imp_i <- amelia.dat$imputations[[i]]
        imp_i <- imp_i[bin.na == 1, ]
        mat.work[,i] <- imp_i$hba1c
        rm(imp_i)
    }
    imp.mean <- rowMeans(mat.work)
    
    dat.base <- orig.dat[orig.dat$time_base ==0, ]
    hba1c_observed <- dat.base$hba1c[!is.na(dat.base$hba1c)]
    out.df <- data.frame(c(hba1c_observed, imp.mean), c(rep("Yes", length(hba1c_observed)), rep("No", length(imp.mean))))
    colnames(out.df) <- c("HbA1c", "Observed")
    out.df
}
data.comp.dens.t0std <- function(amelia.dat, orig.dat, m=10){
    bin.na <- as.numeric(is.na(orig.dat$hba1c_std))
    n_rows <- sum(bin.na)
    mat.work <- matrix(0, nrow=n_rows, ncol=m)
    for(i in 1:m){
        imp_i <- amelia.dat$imputations[[i]]
        imp_i <- imp_i[bin.na == 1, ]
        mat.work[,i] <- imp_i$hba1c_std
        rm(imp_i)
    }
    imp.mean <- rowMeans(mat.work)
    
    dat.base <- orig.dat[orig.dat$time_base ==0, ]
    hba1c_observed <- dat.base$hba1c_std[!is.na(dat.base$hba1c_std)]
    out.df <- data.frame(c(hba1c_observed, imp.mean), c(rep("Yes", length(hba1c_observed)), rep("No", length(imp.mean))))
    colnames(out.df) <- c("HbA1c", "Observed")
    out.df
}
### Density plot for mean imputed versus observed only considering baseline
mean_hba1c_imp0 <- data.comp.dens.t0std(amelia.dat=hba1c.imp44_5, orig.dat=dt.s)
p0 <- ggplot(data=mean_hba1c_imp0, aes(x=HbA1c, color=Observed, fill=Observed)) + geom_density(alpha=0.6)
pdf(paste(analysis_dir, "plots/ameliatest_dtnormal_min44pt5.pdf", sep=""))
    p0
dev.off()

mean_hba1c_imp <- data.comp.dens.t0std(amelia.dat=hba1c.imp, orig.dat=dt.s)
p1 <- ggplot(data=mean_hba1c_imp, aes(x=HbA1c, color=Observed, fill=Observed)) + geom_density(alpha=0.6)
pdf(paste(analysis_dir, "plots/ameliatest_dtnormal.pdf", sep=""))
    p1
dev.off()

mean_hba1c_imp2 <- data.comp.dens.t0(amelia.dat=hba1c.implog, orig.dat=dt.s)
p2 <- ggplot(data=mean_hba1c_imp2, aes(x=HbA1c, color=Observed, fill=Observed)) + geom_density(alpha=0.6)
pdf(paste(analysis_dir, "plots/ameliatest_dtlognormal.pdf", sep=""))
    p2
dev.off()

mean_hba1c_imp3 <- data.comp.dens.t0std(amelia.dat=hba1c.impnormal, orig.dat=dt.s)
p3 <- ggplot(data=mean_hba1c_imp3, aes(x=HbA1c, color=Observed, fill=Observed)) + geom_density(alpha=0.6)
pdf(paste(analysis_dir, "plots/ameliatest_normal_distn.pdf", sep=""))
    p3
dev.off()

mean_hba1c_imp4 <- data.comp.dens.t0(amelia.dat=hba1c.implog_43, orig.dat=dt.s)
p4 <- ggplot(data=mean_hba1c_imp4, aes(x=HbA1c, color=Observed, fill=Observed)) + geom_density(alpha=0.6)
pdf(paste(analysis_dir, "plots/ameliatest_dtlognormal_min43.pdf", sep=""))
    p4
dev.off()

### Normal distibution  = no good. Imputed values lower on average and more spread.
### Double truncated normal, with lower trunc = 43. Good, central point slightly higher HbA1c for imputed compared to observed. Spread good.
### Double truncated normal, with lower trunc = 44.5. Not good,
# central tendency too high.
### Double truncated log-normal with lower trunc = 43. Pretty good,
# central tendency slightly too low.
### Double truncated log-normal with lower trunc = 44.5. Very good. Selected

### Looking at plots, using a log gives rise to an average imputed HbA1c_diag value more alined with the observed values...
### Move forwards with log...
###
rm(hba1c.imp, hba1c.implog, hba1c.impnormal, hba1c.impnormal44_5, hba1c.implog_43)
rm(bound.mat, p1, p2, bound.mat0, bound.mat3, p3, p4, p0)
rm(mean_hba1c_imp4, mean_hba1c_imp3, mean_hba1c_imp2, mean_hba1c_imp, mean_hba1c_imp0)

############################################################################
### 5. Imputation model to generate data to use in LDA analysis 2 and 3
############################################################################
### Set number of imputations at 50
s <- 668
set.seed(s) ### setting a seed in case need to reproduce.
hba1c.imp <- amelia(df.imp2, m=50, ts="time_base", cs="eid", splinetime=3, lags="hba1c", leads="hba1c", ords=c("med_f"), 
    noms=c("med_base_f", "ac_f", "sex", "dep_change_base", "pop_f", "sr_ethnicity_f",
    "ever_smk", "alcohol_never", "qualifications_f", "hba1c_over48", "medication_gp", "gp_t2d", "hes_diag_dm", 
    "t2d_sr_diagnosis_clean", "dep_base", "pret2dmdd_sgo75pc", "pret2dmdd_sg25to75pc", "pret2dmdd_sgl10pc"), 
    p2s = 1, , logs=c("hba1c"), bounds=bound.mat2)

save.image(paste(analysis_dir, "data/hba1c_amelia50.RData", sep=""))
load(paste(analysis_dir, "data/hba1c_amelia50.RData", sep=""))
### Note: could also write output to file using write.amelia().

############################################################################
### Plots...
############################################################################
### (1) mean imputed versus observed only considering baseline
mean_hba1c_imp <- data.comp.dens.t0(amelia.dat=hba1c.imp, orig.dat=dt.s, m=50)
p1 <- ggplot(data=mean_hba1c_imp, aes(x=HbA1c, color=Observed, fill=Observed)) + geom_density(alpha=0.6)
pdf(paste(analysis_dir, "plots/amelia_comp_density_baselineHbA1c.pdf", sep=""))
    p1
dev.off()
### (2) mean imputed versus observed-  default from Amelia considering all time points
pdf(paste(analysis_dir, "plots/amelia_comp_density_default.pdf", sep=""))
compare.density(hba1c.imp, var = "hba1c")
compare.density(hba1c.imp, var = "bmi_std")
compare.density(hba1c.imp, var = "bmi_ukb_std")
compare.density(hba1c.imp, var = "tdi_std")
compare.density(hba1c.imp, var = "sbp_base_std")
compare.density(hba1c.imp, var = "dbp_base_std")
compare.density(hba1c.imp, var = "sbp_max_std")
compare.density(hba1c.imp, var = "mdbp_std")
dev.off()

### Default over-imputation available in Amelia not considered. Firstly, it takes a long time to run. Secondly, it is not designed to accurately describe the imputation quality when only imputing part of a variable (baseline HbA1c here) which is known to follow a different distribution (more skewed) compared to HbA1c at other time points (since for some individuals it is the value at this baseline HbA1c measure that determines T2D status and so selection into the LDA part of the study).
