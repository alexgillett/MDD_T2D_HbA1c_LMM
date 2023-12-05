############################################################################
### 7.3. Does MDD impact residual within-person HbA1c variability?
############################################################################
### Authors: Alexandra C Gillett
############################################################################
### In script we do the following:
### 1. Load in all relevant model output data
### 2. Test if pre-T2D MDD individuals have different residual variation to
### all others.
### 3. Test if post-T2D MDD individuals have different residual variation to
### all others after their MDD diagnosis.
### 4. Test if post-T2D MDD individuals have different residual variation to
### all others looking across all follow-up.
############################################################################
# User add library path where R packages are stored if required:
#.libPaths(new = "your_rlibrary_path")
############################################################################
### Add R packages:
.libPaths(new = "/scratch/prj/ukbiobank/usr/alex_gillett/software/Rpackages/")
# Load libraries
library(data.table)
library(lubridate)
library(ggplot2)
library(nlme)
library(lme4) 
library(MASS)
library(rms) 
library(Amelia)
library(AICcmodavg)
library(lmeInfo)
library(msm)
library(mitml)
### Assign paths to relevant directories
hba1c_dat_path <- "/path_where_you_to_store_extracted_data/hba1c_data/"
analysis_dir <- "/path_to_analysis_dir/"

############################################################################
### 1. Load in all relevant model output data
############################################################################
### REML, original data without variance being a function of MDD
############################################################################
load(file=paste(analysis_dir,"model_output/adjusted_mainREML1.RData", sep=""))
load(file=paste(analysis_dir,"model_output/adjusted_mainREML2.RData", sep=""))
load(file=paste(analysis_dir,"model_output/adjusted_mainREML3.RData", sep=""))
load(file=paste(analysis_dir,"model_output/adjusted_mainREML4.RData", sep=""))
load(file=paste(analysis_dir,"model_output/adjusted_mainREML5.RData", sep=""))
### Combine and remove individual model output datasets
fit_adj_mainREML <- c(fit_adj_mainREML1, fit_adj_mainREML2, fit_adj_mainREML3, fit_adj_mainREML4, fit_adj_mainREML5)
rm(fit_adj_mainREML1, fit_adj_mainREML2, fit_adj_mainREML3, fit_adj_mainREML4, fit_adj_mainREML5)
############################################################################
### REML, data with pre-T2D MDD individuals allowed to have differing residual HbA1c variation
############################################################################
load(file=paste(analysis_dir, "model_output/var_pret2d_diffREML1.RData", sep=""))
load(file=paste(analysis_dir, "model_output/var_pret2d_diffREML2.RData", sep=""))
load(file=paste(analysis_dir, "model_output/var_pret2d_diffREML3.RData", sep=""))
load(file=paste(analysis_dir, "model_output/var_pret2d_diffREML4.RData", sep=""))
load(file=paste(analysis_dir, "model_output/var_pret2d_diffREML5.RData", sep=""))
### Combine and remove individual model output datasets
fit_pret2d_diff_REML <- c(fit_pret2d_diff_REML1, fit_pret2d_diff_REML2, fit_pret2d_diff_REML3, fit_pret2d_diff_REML4, fit_pret2d_diff_REML5)
rm(fit_pret2d_diff_REML1, fit_pret2d_diff_REML2, fit_pret2d_diff_REML3, fit_pret2d_diff_REML4, fit_pret2d_diff_REML5)
############################################################################
### REML, data with post-T2D MDD (time-varying) individuals allowed to have differing residual HbA1c variation
############################################################################
load(file=paste(analysis_dir, "model_output/var_postt2dTV_diffREML1.RData", sep=""))
load(file=paste(analysis_dir, "model_output/var_postt2dTV_diffREML2.RData", sep=""))
load(file=paste(analysis_dir, "model_output/var_postt2dTV_diffREML3.RData", sep=""))
load(file=paste(analysis_dir, "model_output/var_postt2dTV_diffREML4.RData", sep=""))
load(file=paste(analysis_dir, "model_output/var_postt2dTV_diffREML5.RData", sep=""))
### Combine and remove individual model output datasets
fit_postt2dTV_diff_REML <- c(fit_postt2dTV_diff_REML1, fit_postt2dTV_diff_REML2, fit_postt2dTV_diff_REML3, fit_postt2dTV_diff_REML4, fit_postt2dTV_diff_REML5)
rm(fit_postt2dTV_diff_REML1, fit_postt2dTV_diff_REML2, fit_postt2dTV_diff_REML3, fit_postt2dTV_diff_REML4, fit_postt2dTV_diff_REML5)
############################################################################
### REML, data with (retro) post-T2D MDD individuals allowed to have differing residual HbA1c variation
############################################################################
load(file=paste(analysis_dir, "model_output/var_postt2dretro_diffREML1.RData", sep=""))
load(file=paste(analysis_dir, "model_output/var_postt2dretro_diffREML2.RData", sep=""))
load(file=paste(analysis_dir, "model_output/var_postt2dretro_diffREML3.RData", sep=""))
load(file=paste(analysis_dir, "model_output/var_postt2dretro_diffREML4.RData", sep=""))
load(file=paste(analysis_dir, "model_output/var_postt2dretro_diffREML5.RData", sep=""))
### Combine and remove individual model output datasets
fit_retro_postt2d_diff_REML <- c(fit_postt2dretro_diff_REML1, fit_postt2dretro_diff_REML2, fit_postt2dretro_diff_REML3, fit_postt2dretro_diff_REML4, fit_postt2dretro_diff_REML5)
rm(fit_postt2dretro_diff_REML1, fit_postt2dretro_diff_REML2, fit_postt2dretro_diff_REML3, fit_postt2dretro_diff_REML4, fit_postt2dretro_diff_REML5)

############################################################################
### 2. Test if pre-T2D MDD individuals have different residual variation to
### all others.
############################################################################
### Compare pre-T2D MDD variance model to original
############################################################################
### Note. Non-convergence for models from MI dataset = 3, 7, 9, 15, 20, 25, 28, 36, 41, 45, 46, 50.
### 12 models with non-convergence. 24% of models
exclude_m <- c(3, 7, 9, 15, 20, 25, 28, 36, 41, 45, 46, 50)
out.p1 <- NULL
for(i in 1:50){
    pi <- anova(fit_pret2d_diff_REML[[i]], fit_adj_mainREML[[i]])$p[2]
    out.p1 <- c(out.p1, pi)
    rm(pi)
}
out.p1 <- out.p1[-exclude_m]
median(out.p1) 
### Parameter estimates
summ_reml_pret2d <- testEstimates(fit_pret2d_diff_REML, extra.pars=T)

### Function to get mean and 95% CI...
m <- 50
x.reint <- NULL
se.reint <- NULL
x.reslope <- NULL
se.reslope <- NULL
x.phi <- NULL
se.phi <- NULL
x.dept2d <- NULL
se.dept2d <- NULL
#x.t2ddep <- NULL
#se.t2ddep <- NULL
x.sigma <- NULL
se.sigma <- NULL
for(i in 1:m){
    if(!any(exclude_m == i)){
    var <- fit_pret2d_diff_REML[[i]]$apVar
    par<-attr(var, "Pars")
    x.reint <- c(x.reint, exp(par[1]))
    x.reslope <- c(x.reslope, exp(par[2]))
    x.phi <- c(x.phi, exp(par[4]))
    x.dept2d <- c(x.dept2d, exp(par[5]))
 #   x.t2ddep <- c(x.t2ddep, exp(par[6]))
    x.sigma <- c(x.sigma, exp(par[6]))
    se.reint <- c(se.reint, deltamethod(~ exp(x1), par, var))
    se.reslope <- c(se.reslope, deltamethod(~ exp(x2), par, var))
    se.phi <- c(se.phi, deltamethod(~ exp(x4), par, var))
    se.dept2d <- c(se.dept2d, deltamethod(~ exp(x5), par, var))
 #   se.t2ddep <- c(se.t2ddep, deltamethod(~ exp(x6), par, var))
    se.sigma <- c(se.sigma, deltamethod(~ exp(x6), par, var))
    rm(par, var)
    }
}
point_est_rubin_f2 <- function(x, se){
    ### x here is a vector of point estimates, se = SE
    pev <- mean(x)
    m <- length(x)
    varpe <- se^2
    Ubar <- mean(varpe)
    B <- sum((x - pev)^2)/(m-1)
    T <- Ubar + ((1 + (1/m))*B)
    seT <- sqrt(T)
    lower <- pev - 1.96*seT
    upper <- pev + 1.96*seT
    out <- c(pev, seT, lower, upper)
    out
}
point_est_rubin_f2(x= x.reint, se=se.reint) # 0.491534445 0.004372927 0.482963508 0.500105382
point_est_rubin_f2(x= x.reslope, se=se.reslope) # 0.230364960 0.003412321 0.223676810 0.237053109
point_est_rubin_f2(x= x.phi, se=se.phi) # 0.046592117 0.001425617 0.043797908 0.049386326
point_est_rubin_f2(x= x.dept2d, se=se.dept2d) # 0.983232207 0.007660051 0.968218507 0.998245907
#point_est_rubin_f2(x= x.t2ddep, se=se.t2ddep) # 
point_est_rubin_f2(x= x.sigma, se=se.sigma) # 0.542833276 0.001498478 0.539896259 0.545770293

############################################################################
### 3. Test if post-T2D MDD individuals have different residual variation to
### all others after their MDD diagnosis.
############################################################################
### Compare post-T2D MDD (time-varying) variance model to original
############################################################################
### Note. Non-convergence for models from MI dataset = 2, 5, 6, 7, 9, 10, 14, 17, 18, 25, 28, 29, 31, 34, 35, 36, 37, 38, 40, 42, 45, 49, 50
### 24 models with non-convergence. 48% of models
exclude_m2 <- c(2, 5, 6, 7, 9, 10, 14, 17, 18, 25, 28, 29, 31, 32, 34, 35, 36, 37, 38, 40, 42, 45, 49, 50)
out.p2 <- NULL
for(i in 1:50){
    pi <- anova(fit_postt2dTV_diff_REML[[i]], fit_adj_mainREML[[i]])$p[2]
    out.p2 <- c(out.p2, pi)
    rm(pi)
}
out.p2 <- out.p2[-exclude_m2]
median(out.p2) 

### Function to get mean and 95% CI...
m <- 50
x.reint <- NULL
se.reint <- NULL
x.reslope <- NULL
se.reslope <- NULL
x.phi <- NULL
se.phi <- NULL
#x.dept2d <- NULL
#se.dept2d <- NULL
x.t2ddep <- NULL
se.t2ddep <- NULL
x.sigma <- NULL
se.sigma <- NULL 
### non-convergence for models in datasets = 2, 5, 6, 7, 9, 10, 14, 17, 18, 25, 28, 29, 31,32, 34, 35, 36
### 37, 38, 40, 42, 45, 49, 50
for(i in 1:m){
    if(!any(exclude_m2 == i)){
    var <- fit_postt2dTV_diff_REML[[i]]$apVar
    par<-attr(var, "Pars")
    x.reint <- c(x.reint, exp(par[1]))
    x.reslope <- c(x.reslope, exp(par[2]))
    x.phi <- c(x.phi, exp(par[4]))
    #x.dept2d <- c(x.dept2d, exp(par[5]))
    x.t2ddep <- c(x.t2ddep, exp(par[5]))
    x.sigma <- c(x.sigma, exp(par[6]))
    se.reint <- c(se.reint, deltamethod(~ exp(x1), par, var))
    se.reslope <- c(se.reslope, deltamethod(~ exp(x2), par, var))
    se.phi <- c(se.phi, deltamethod(~ exp(x4), par, var))
    #se.dept2d <- c(se.dept2d, deltamethod(~ exp(x5), par, var))
    se.t2ddep <- c(se.t2ddep, deltamethod(~ exp(x5), par, var))
    se.sigma <- c(se.sigma, deltamethod(~ exp(x6), par, var))
    rm(par, var)
    }
}

point_est_rubin_f2(x= x.reint, se=se.reint) # 0.491399818 0.004444909 0.482687797 0.500111839
point_est_rubin_f2(x= x.reslope, se=se.reslope) # 0.230421860 0.003415363 0.223727748 0.237115973
point_est_rubin_f2(x= x.phi, se=se.phi) # 0.04657307 0.00142255 0.04378487 0.04936127
#point_est_rubin_f2(x= x.dept2d, se=se.dept2d) # 
point_est_rubin_f2(x= x.t2ddep, se=se.t2ddep) # 0.9732649 0.0113030 0.9511110 0.9954188
point_est_rubin_f2(x= x.sigma, se=se.sigma) # 0.542205992 0.001432907 0.539397495 0.545014490

############################################################################
### 4. Test if post-T2D MDD individuals have different residual variation to
### all others looking across all follow-up.
############################################################################
### Compare post-T2D MDD (retrospective) variance model to original
############################################################################
out.p3 <- NULL
for(i in 1:50){
    pi <- anova(fit_retro_postt2d_diff_REML[[i]], fit_adj_mainREML[[i]])$p[2]
    out.p3 <- c(out.p3, pi)
    rm(pi)
}
median(out.p3) 

### Function to get mean and 95% CI...
m <- 50
x.reint <- NULL
se.reint <- NULL
x.reslope <- NULL
se.reslope <- NULL
x.phi <- NULL
se.phi <- NULL
#x.dept2d <- NULL
#se.dept2d <- NULL
x.t2ddep <- NULL
se.t2ddep <- NULL
x.sigma <- NULL
se.sigma <- NULL 

for(i in 1:m){
    #if(!any(exclude_m2 == i)){
    var <- fit_retro_postt2d_diff_REML[[i]]$apVar
    par<-attr(var, "Pars")
    x.reint <- c(x.reint, exp(par[1]))
    x.reslope <- c(x.reslope, exp(par[2]))
    x.phi <- c(x.phi, exp(par[4]))
    #x.dept2d <- c(x.dept2d, exp(par[5]))
    x.t2ddep <- c(x.t2ddep, exp(par[5]))
    x.sigma <- c(x.sigma, exp(par[6]))
    se.reint <- c(se.reint, deltamethod(~ exp(x1), par, var))
    se.reslope <- c(se.reslope, deltamethod(~ exp(x2), par, var))
    se.phi <- c(se.phi, deltamethod(~ exp(x4), par, var))
    #se.dept2d <- c(se.dept2d, deltamethod(~ exp(x5), par, var))
    se.t2ddep <- c(se.t2ddep, deltamethod(~ exp(x5), par, var))
    se.sigma <- c(se.sigma, deltamethod(~ exp(x6), par, var))
    rm(par, var)
   # }
}

point_est_rubin_f2(x= x.reint, se=se.reint) # 0.491399818 0.004444909 0.482687797 0.500111839
point_est_rubin_f2(x= x.reslope, se=se.reslope) # 0.230421860 0.003415363 0.223727748 0.237115973
point_est_rubin_f2(x= x.phi, se=se.phi) # 0.04657307 0.00142255 0.04378487 0.04936127
#point_est_rubin_f2(x= x.dept2d, se=se.dept2d) # 
point_est_rubin_f2(x= x.t2ddep, se=se.t2ddep) # 1.16444759 0.01506678 1.13491670 1.19397848
point_est_rubin_f2(x= x.sigma, se=se.sigma) # 0.542205992 0.001432907 0.539397495 0.545014490
