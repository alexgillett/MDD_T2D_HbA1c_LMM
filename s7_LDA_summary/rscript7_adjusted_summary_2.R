############################################################################
### 7.1. Summary of the adjusted LDA model
############################################################################
### Authors: Alexandra C Gillett
############################################################################
### In script we:
### 1. Load in model output
### 2. P-values for fixed effects
### 3. Pooled summary of model output- tables
############################################################################
# User add library path where R packages are stored if required:
#.libPaths(new = "your_rlibrary_path")
############################################################################
### Add R packages:
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
### 1. Load in model output
############################################################################
### ML
load(file= paste(analysis_dir, "model_output/adjusted_mainML1.RData", sep=""))
load(file= paste(analysis_dir, "model_output/adjusted_mainML2.RData", sep=""))
load(file= paste(analysis_dir, "model_output/adjusted_mainML3.RData", sep=""))
load(file= paste(analysis_dir, "model_output/adjusted_mainML4.RData", sep=""))
load(file= paste(analysis_dir, "model_output/adjusted_mainML5.RData", sep=""))
### Combine:
fit_adj_mainML <- c(fit_adj_mainML1, fit_adj_mainML2, fit_adj_mainML3, fit_adj_mainML4, fit_adj_mainML5)
### Remove individual ML output datasets
rm(fit_adj_mainML1, fit_adj_mainML2, fit_adj_mainML3, fit_adj_mainML4, fit_adj_mainML5)
### REML
load(file= paste(analysis_dir, "model_output/adjusted_mainREML1.RData", sep=""))
load(file= paste(analysis_dir, "model_output/adjusted_mainREML2.RData", sep=""))
load(file= paste(analysis_dir, "model_output/adjusted_mainREML3.RData", sep=""))
load(file= paste(analysis_dir, "model_output/adjusted_mainREML4.RData", sep=""))
load(file= paste(analysis_dir, "model_output/adjusted_mainREML5.RData", sep=""))
### Combine:
fit_adj_mainREML <- c(fit_adj_mainREML1, fit_adj_mainREML2, fit_adj_mainREML3, fit_adj_mainREML4, fit_adj_mainREML5)
### Remove individual REML output datasets
rm(fit_adj_mainREML1, fit_adj_mainREML2, fit_adj_mainREML3, fit_adj_mainREML4, fit_adj_mainREML5)
############################################################################
### 2. P-values for fixed effects
############################################################################
### Test for pre-T2D MDD variables:
############################################################################
pret2dmdd_c <- c("dep_base", "uncentered_robust_pretime2")
testConstraints(fit_adj_mainML, constraints = pret2dmdd_c, method="D1")
testConstraints(fit_adj_mainML, constraints = pret2dmdd_c, method="D2")[[3]][4] ### 0.1884773
### Not significant after adjusting for additional covariates
### Just dep_base == 0
testConstraints(fit_adj_mainML, constraints = "dep_base", method="D2")[[3]][4] ### 0.06734393
############################################################################
### Test for post-T2D MDD variables:
############################################################################
postt2dmdd_c <- c("dep_change_base", "postt2dmdd_time")
testConstraints(fit_adj_mainML, constraints = postt2dmdd_c, method="D2")
testConstraints(fit_adj_mainML, constraints = postt2dmdd_c, method="D2")[[3]][4] ### 0.01523378
### Still marginally significant after adjusting for additional covariates
testConstraints(fit_adj_mainML, constraints = "dep_change_base", method="D2")[[3]][4] ### 0.2554125
testConstraints(fit_adj_mainML, constraints = "postt2dmdd_time", method="D2")[[3]][4] ### 0.003832691
### driven by postt2dmdd_time

############################################################################
### 3. Pooled summary of model output- tables
############################################################################
### Tables to summarise results mitml::testEstimates
############################################################################
### Fixed effects summary from ML (want p-values from this):
basic.summ.outML <- testEstimates(fit_adj_mainML, extra.pars=T)
out.tabML <- basic.summ.outML$estimates
out.tabML
### Summary from REML:
summ.reml <- testEstimates(fit_adj_mainREML, extra.pars=T)
### Fixed effects from REML
out.tab.reml <- summ.reml$estimates
### Random effects from REML
out.re.reml <- summ.reml$extra.pars
CI_REML <- confint(testEstimates(fit_adj_mainREML, extra.pars=T))
### Write files to output
write.csv(CI_REML, file=paste(analysis_dir, "model_output/adj_REML_CIfe_pooled.csv", sep=""))
write.csv(out.tab.reml, file=paste(analysis_dir, "model_output/adj_REMLfe_pooled.csv", sep=""))
write.csv(out.re.reml, file=paste(analysis_dir, "model_output/adj_REMLre_pooled.csv", sep=""))

### Create file with Estimates (REML), 95% CI (REML) and pvalues (ML) for fixed effects
CI_REML2 <- paste(paste("(", as.character(round(CI_REML[,1], 2)), sep=""), 
paste(as.character(round(CI_REML[,2], 2)), ")", sep=""), sep = ", ")
fe.out.tab <- cbind(out.tab.reml[,1], CI_REML2, out.tabML[, 5])
colnames(fe.out.tab) <- c("Estimate", "95CI", "pvalue")
write.csv(fe.out.tab, file=paste(analysis_dir, "model_output/adj_fe_pooled_tabSM.csv", sep=""))
rm(fe.out.tab)

### Unstandardised updated params for exposures...
Q3 <- 17.61943
Q1 <- 5.026712
CIexp <- CI_REML[(dim(CI_REML)[1]-3):dim(CI_REML)[1], ]
coef_exp <- out.tab.reml[(dim(CI_REML)[1]-3):dim(CI_REML)[1],1]
p <- out.tab.reml[(dim(CI_REML)[1]-3):dim(CI_REML)[1],5]
sd_hba1c <- std_params[1,2]
summ_exp <- cbind(coef_exp,CIexp)*sd_hba1c
summ_exp[2,] <- summ_exp[2,]/(Q3-Q1)
fe_summ_exp <- data.frame(summ_exp[,1], 
paste(paste("(", as.character(round(summ_exp[,2], 2)), sep=""), 
paste(as.character(round(summ_exp[,3], 2)), ")", sep=""), sep=", "), 
p)
colnames(fe_summ_exp) <- c("Estimate", "95CI", "p_value")
write.csv(fe_summ_exp, file=paste(analysis_dir, "model_output/adj_unstdHbA1c_exposures_fe_summ.csv", sep=""))

### Do not have intervals for RE or variance-covariance params...
### Function to get mean and 95% CI...
m <- 50
x.reint <- NULL
se.reint <- NULL
x.reslope <- NULL
se.reslope <- NULL
x.phi <- NULL
se.phi <- NULL
x.sigma <- NULL
se.sigma <- NULL
for(i in 1:m){
    var <- fit_adj_mainREML[[i]]$apVar
    par<-attr(var, "Pars")
    x.reint <- c(x.reint, exp(par[1]))
    x.reslope <- c(x.reslope, exp(par[2]))
    x.phi <- c(x.phi, exp(par[4]))
    x.sigma <- c(x.sigma, exp(par[5]))
    se.reint <- c(se.reint, deltamethod(~ exp(x1), par, var))
    se.reslope <- c(se.reslope, deltamethod(~ exp(x2), par, var))
    se.phi <- c(se.phi, deltamethod(~ exp(x4), par, var))
    se.sigma <- c(se.sigma, deltamethod(~ exp(x5), par, var))
    rm(par, var)
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
### Manually add the following to pooled model summary table:
point_est_rubin_f2(x= x.reint, se=se.reint) # 0.491570653 0.004370099 0.483005259 0.500136047
point_est_rubin_f2(x= x.reslope, se=se.reslope) # 0.230406263 0.003393407 0.223755186 0.237057341
point_est_rubin_f2(x= x.phi, se=se.phi) # 0.04651344 0.00141974 0.04373075 0.04929613
point_est_rubin_f2(x= x.sigma, se=se.sigma) # 0.541960458 0.001430698 0.539156291 0.544764626
#attributes(intervals(fit_adj_mainREML[[1]]))
corr.re <- matrix(0, nrow=50, ncol=3)
for(i in 1:m){
    corr.re[i,] <- unlist(intervals(fit_adj_mainREML[[i]])$reStruct[[1]][3,])
}
colMeans(corr.re) # 0.4135670 0.4412037 0.4680273
