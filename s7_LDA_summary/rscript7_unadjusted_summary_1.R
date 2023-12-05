############################################################################
### 7.1. Summary of the unadjusted LDA models
############################################################################
### Authors: Alexandra C Gillett
############################################################################
### In script we do/ obtain the following:
### 1. P-values for fixed effects
### 2. Summarising selected model output
### 2.1. Summarising model output into tables
### 2.2. Summarising model output using example individuals- plots
### 2.3. Summarising model output using example individuals- tables
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
### 1. P-values for fixed effects
############################################################################
### Load in ML models
load(file= paste(analysis_dir, "model_output/unadj_ML.RData", sep=""))
############################################################################
### Test for pre-T2D duration interaction
############################################################################
anova(fit_unadj_pret2d_int_allML, fit_unadj_mainML)
#                            Model df      AIC      BIC    logLik   Test  L.Ratio
#fit_unadj_pret2d_int_allML     1 19 247790.8 247979.4 -123876.4                
#fit_unadj_mainML               2 13 247786.0 247915.0 -123880.0 1 vs 2 7.149783
#                           p-value
#fit_unadj_pret2d_int_allML        
#fit_unadj_mainML            0.3072 No evidence for interaction
############################################################################
### Test any pre-T2D MDD effect
############################################################################
anova(fit_unadj_mainML, fit_unadj_nopreML)
#                   Model df      AIC      BIC    logLik   Test  L.Ratio p-value
#fit_unadj_mainML      1 13 247786.0 247915.0 -123880.0                        
#fit_unadj_nopreML     2 11 247792.6 247901.7 -123885.3 1 vs 2 10.61108   0.005
anova(fit_unadj_mainML, fit_unadj_nopreML)[2,9]
### p = 0.004964008
### Evidence for pre-T2D MDD effect
summary(fit_unadj_mainML)$tTable
############################################################################
### Test for pre-T2D MDD duration (until baseline) effect
############################################################################
anova(fit_unadj_mainML, fit_unadj_nopredurML)
# anova(fit_unadj_mainML, fit_unadj_nopredurML)
#                    Model df    AIC      BIC  logLik   Test  L.Ratio p-value
#fit_unadj_mainML         1 13 247786 247915.0 -123880                        
#fit_unadj_nopredurML     2 12 247794 247913.1 -123885 1 vs 2 10.04175  0.0015
############################################################################
### Test post-T2D MDD duration/ time effect:
############################################################################
anova(fit_unadj_nopostdurML, fit_unadj_mainML)
#                      Model df      AIC      BIC    logLik   Test  L.Ratio
#fit_unadj_nopostdurML     1 12 247793.5 247912.6 -123884.7                
#fit_unadj_mainML          2 13 247786.0 247915.0 -123880.0 1 vs 2 9.519033
#                      p-value
#fit_unadj_nopostdurML        
#fit_unadj_mainML        0.002
############################################################################
### Test for any post-T2D MDD effect
############################################################################
anova(fit_unadj_mainML, fit_unadj_nopostML)
#                    Model df    AIC      BIC  logLik   Test  L.Ratio p-value
#fit_unadj_mainML       1 13 247786 247915.0 -123880                        
#fit_unadj_nopostML     2 11 247792 247901.2 -123885 1 vs 2 10.05671  0.0065
anova(fit_unadj_mainML, fit_unadj_nopostML)[2,9]
# p = 0.006549583
### Evidence for post-T2D MDD effect on HbA1c- this comes through the semi-continuous variable for time since depression.
anova(fit_unadj_nopostbaseML, fit_unadj_mainML)
############################################################################
### Test for post-T2D MDD splines (Supplemetary materials)
############################################################################
anova(fit_unadj_mainML, fit_unadj_postsplineML)
#                        Model df      AIC      BIC    logLik   Test  L.Ratio
#fit_unadj_mainML           1 13 247786.0 247915.0 -123880.0                
#fit_unadj_postsplineML     2 15 247786.8 247935.7 -123878.4 1 vs 2 3.110491
 #                      p-value
#fit_unadj_mainML              
#fit_unadj_postsplineML  0.2111
anova(fit_unadj_mainML, fit_unadj_postsplineML)[2,9]
### p = 0.2111375
### No evidence for non-linear post-T2D splines...
############################################################################
### Test for pre-T2D MDD interacting with time-splines...
############################################################################
### interaction with binary pre-T2D MDD variable only
anova(fit_unadj_mainML, fit_unadj_pret2d_int_baseML)[2,9] ###  0.1583923
### interaction with pre-T2D MDD duration variable only
anova(fit_unadj_mainML, fit_unadj_pret2d_int_durML)[2,9] ### 0.2265567
### interaction with any pre-T2D MDD variables
anova(fit_unadj_mainML, fit_unadj_pret2d_int_allML)[2,9] ### 0.3072173
### No evidence that pre-T2D interacts with time-splines. 
# anova(fit_unadj_pret2d_int_linML, fit_unadj_mainML)[2,9] ### 0.1679341
### No evidence for linear interaction either...

############################################################################
### 2. Summarising selected model output
############################################################################
### Primary unadjusted model ML version = fit_unadj_mainML, REML version == fit_unadj_mainREML
############################################################################
### 2.1. Summarising model output into tables
############################################################################
### Load in REML model
load(paste(analysis_dir, "model_output/unadj_REML.RData", sep=""))
### Extract Fixed effects summary from ML model
feML <- summary(fit_unadj_mainML)$tTable
### Extract Fixed effects summary from REML model
feREML <- summary(fit_unadj_mainREML)$tTable
### Get intervals for fixed effects from REML model
int_fe <- intervals(fit_unadj_mainREML)$fixed
### Get intervals for random effects from REML model
int_re <- intervals(fit_unadj_mainREML)$reStruct
### Get intervals for CAR1 parameter from REML model
int_phi <- intervals(fit_unadj_mainREML)$corStruct
### Get intervals for residual variance parameter from REML model
int_sigma <- intervals(fit_unadj_mainREML)$sigma

### Create a data frame for fixed effects summary
### Contains: estimate (REML), lower 95% CI (REML), upper 95% CI (REML), p-value (ML)
fe_summ <- data.frame(int_fe[,2],
paste(paste("(", as.character(round(int_fe[,1], 2)), sep=""), 
paste(as.character(round(int_fe[,3], 2)), ")", sep=""), sep=", "), 
feML[,5])
colnames(fe_summ) <- c("Estimate", "95CI", "p_value")

### Create a data frame for summary info for:
### random effects and variance-covariance parameters
### Contains: estimate (REML), lower 95% CI (REML), upper 95% CI (REML)
out_add <- rbind(int_re$eid_f, int_sigma, int_phi)
add_params_summ <- data.frame(out_add[,2],
paste(paste("(", as.character(round(out_add[,1], 2)), sep=""), 
paste(as.character(round(out_add[,3], 2)), ")", sep=""), sep=", ")
)
colnames(add_params_summ) <- c("Estimate", "95CI")
rownames(add_params_summ) <- c("sd((Intercept))", "sd(time_base_std) ", "cor((Intercept),time_base_std)", 
"sigma", "phi")

### Write summaries to output
write.csv(fe_summ, file= paste(analysis_dir,  "model_output/unadj_fixed_effects_summary.csv", sep=""))
write.csv(add_params_summ, file= paste(analysis_dir, "model_output/unadj_re_sigma_phi_summ.csv", sep=""))

### The above fixed effects parameters are:
# 1) for standardised HbA1c, and,
# 2) for a centred version of the pre-T2D MDD duration variable.
### Create FE summaries for unstandardised and uncentred variables:
### Unstandardised updated params for exposures...
Q3 <- 17.61943
Q1 <- 5.026712
summ_exp <- int_fe[5:8,]*sd_hba1c
summ_exp[2,] <- summ_exp[2,]/(Q3-Q1)
fe_summ_exp <- data.frame(summ_exp[,2], 
paste(paste("(", as.character(round(summ_exp[,1], 2)), sep=""), 
paste(as.character(round(summ_exp[,3], 2)), ")", sep=""), sep=", "), 
feML[5:8,5])
colnames(fe_summ_exp) <- c("Estimate", "95CI", "p_value")
### Write to output
write.csv(fe_summ_exp, file= paste(analysis_dir, "model_output/unadj_unstdHbA1c_exposures_fe_summ.csv", sep=""))

############################################################################
### 2.2. Summarising model output using example individuals- plots
############################################################################
### Want to create some 'example' individuals to plot and discuss model results
### for.
### For pre- and post-T2D MDD duration variables find the 10th, 50th (median)
### and 90th percentile and use these as examples, alongside no MDD.
############################################################################
### Median and 10th and 90th decile for pre-T2D MDD time:
############################################################################
datw <- fit_unadj_mainREML$data
datw <- datw[datw$dep_base == 1, ]
datwdt <- data.table(datw, key="eid")
length(unique(datw$eid)) ### 1119
upredur <- datwdt[, unique(pret2dmdd_time), by="eid"]

median(upredur$V1) ### 10.69041 (approx 10 years 8 months) Use 10.5 years
quantile(upredur$V1, probs=c(0.10)) ### 2.232877 (approx 2 year 3 month) 
quantile(upredur$V1, probs=c(0.90)) ### 27.45973 (approx 27 years 6 months) 
m_pretime <- median(upredur$V1)/(17.61943 - 5.026712) ### 0.8489359
lower_pretime <- quantile(upredur$V1, probs=c(0.10))/(17.61943 - 5.026712) ### 0.1773149
upper_pretime <- quantile(upredur$V1, probs=c(0.90))/(17.61943 - 5.026712) ### 2.180604 
############################################################################
### Median post-T2D MDD when onset:
############################################################################
datw <- fit_unadj_mainREML$data
datw <- datw[datw$dep_change_base == 1, ]
datwdt <- data.table(datw, key="eid")
datpost <- datwdt[, min(time_base), by="eid"]
m_posttime <- median(datpost$V1) ### 3.767123 (approx 3 years 9 months) 
lower_posttime <- quantile(datpost$V1, probs=c(0.10)) ### 0.9622996 (approx 11.5 months) 
upper_posttime <- quantile(datpost$V1, probs=c(0.90)) ### 7.126027 (approx 7 years, 2 months) 
rm(datw)
############################################################################
### For example individuals want to find their expected HbA1c trajectories
### from the model. We need to give the model their 'variables'.
### That is the time splines, and exposure variables for our example
### individuals
############################################################################
### Find time splines for some time values
############################################################################
datw <- fit_unadj_mainREML$data
get_time <- quantile(datw$time_base, probs = seq(from=0, to=1, by=0.002))
time_base <- NULL
time_spline1 <- NULL
time_spline2 <- NULL
time_spline3 <- NULL
for(i in 1:length(get_time)){
    dfi <- datw[datw$time_base == get_time[i], ]
    time_base <- c(time_base, dfi$time_base[1])
    time_spline1 <- c(time_spline1, dfi$time_spline1[1])
    time_spline2 <- c(time_spline2, dfi$time_spline2[1])
    time_spline3 <- c(time_spline3, dfi$time_spline3[1])
}
time_spline1 <- time_spline1[!is.na(time_base)]
time_spline2 <- time_spline2[!is.na(time_base)]
time_spline3 <- time_spline3[!is.na(time_base)]

time_base <- time_base[!is.na(time_base)]
n_time <- length(time_base) ##482

### 1 week, 1 month, 6 months, 1 year, 1.5 years, 2-10 years:
time_base_tab <- c(0.019178082, 0.079452055, 0.501369863, 0.997267760, 1.501369863, 1.997260274, 2.991780822,
3.994520548, 4.983606557, 6.000000000, 6.991803279, 8.010958904, 8.961643836, 10)

############################################################################
### Create 'new' dataset...
############################################################################
### no MDD
############################################################################
timedf <- data.frame(time_base, time_spline1, time_spline2, time_spline3)
dep_base <- rep(0, n_time)
dep_change_base <- rep(0, n_time)
uncentered_robust_pretime2 <- rep(0, n_time)
postspline1 <- rep(0, n_time)
Example <- rep("No MDD", n_time)
dep_group <- rep("No MDD", n_time)
noMDD <- data.frame(timedf, dep_base, dep_change_base, uncentered_robust_pretime2, postspline1, Example, dep_group)
############################################################################
### pre-T2D MDD median
############################################################################
dep_base <- rep(1, n_time)
dep_change_base <- rep(0, n_time)
uncentered_robust_pretime2 <- rep(m_pretime, n_time)
postspline1 <- rep(0, n_time)
Example <- rep("Pre-T2D MDD: 10.7 years before T2D", n_time)
dep_group <- rep("Pre-T2D MDD", n_time)
preMDDmed <- data.frame(timedf, dep_base, dep_change_base, uncentered_robust_pretime2, postspline1, Example, dep_group)
############################################################################
### pre-T2D MDD lower
############################################################################
dep_base <- rep(1, n_time)
dep_change_base <- rep(0, n_time)
uncentered_robust_pretime2 <- rep(lower_pretime, n_time)
postspline1 <- rep(0, n_time)
Example <- rep("Pre-T2D MDD: 2.2 years before T2D", n_time)
dep_group <- rep("Pre-T2D MDD", n_time)
preMDDlower <- data.frame(timedf, dep_base, dep_change_base, uncentered_robust_pretime2, postspline1, Example, dep_group)
############################################################################
### pre-T2D MDD upper
############################################################################
dep_base <- rep(1, n_time)
dep_change_base <- rep(0, n_time)
uncentered_robust_pretime2 <- rep(upper_pretime, n_time)
postspline1 <- rep(0, n_time)
Example <- rep("Pre-T2D MDD: 27.5 years before T2D", n_time)
dep_group <- rep("Pre-T2D MDD", n_time)
preMDDupper <- data.frame(timedf, dep_base, dep_change_base, uncentered_robust_pretime2, postspline1, Example, dep_group)
############################################################################
### post-T2D MDD median
############################################################################
dep_base <- rep(0, n_time)
dep_change_base <- as.numeric(time_base >= m_posttime)
uncentered_robust_pretime2 <- rep(0, n_time)
postspline1 <- time_base - m_posttime
postspline1[postspline1 < 0] <- 0
Example <- rep("Post-T2D MDD: 3.8 years after T2D", n_time)
dep_group <- rep("Post-T2D MDD", n_time)
postMDDmed <- data.frame(timedf, dep_base, dep_change_base, uncentered_robust_pretime2, postspline1, Example, dep_group)
############################################################################
### post-T2D MDD lower
############################################################################
dep_base <- rep(0, n_time)
dep_change_base <- as.numeric(time_base >= lower_posttime)
uncentered_robust_pretime2 <- rep(0, n_time)
postspline1 <- time_base - lower_posttime
postspline1[postspline1 < 0] <- 0
Example <- rep("Post-T2D MDD: 1 year after T2D", n_time)
dep_group <- rep("Post-T2D MDD", n_time)
postMDDlower <- data.frame(timedf, dep_base, dep_change_base, uncentered_robust_pretime2, postspline1, Example, dep_group)
############################################################################
### post-T2D MDD upper
############################################################################
dep_base <- rep(0, n_time)
dep_change_base <- as.numeric(time_base >= upper_posttime)
uncentered_robust_pretime2 <- rep(0, n_time)
postspline1 <- time_base - upper_posttime
postspline1[postspline1 < 0] <- 0
Example <- rep("Post-T2D MDD: 7.1 years after T2D", n_time)
dep_group <- rep("Post-T2D MDD", n_time)
postMDDupper <- data.frame(timedf, dep_base, dep_change_base, uncentered_robust_pretime2, postspline1, Example, dep_group)
############################################################################
### Overall dataframe:
############################################################################
dfnew <- rbind(noMDD, preMDDlower, preMDDmed, preMDDupper, postMDDlower, postMDDmed, postMDDupper)
rm(noMDD, preMDDlower, preMDDmed, preMDDupper, postMDDlower, postMDDmed, postMDDupper)
############################################################################
### Calculated expected HbA1c (standardised) using model for these examples
############################################################################
prednew <- predictSE(fit_unadj_mainREML, newdata=dfnew, level=0, se.fit=T)
dfnew$fit <- prednew$fit
dfnew$se <- prednew$se.fit
### Find 95% CI
dfnew$lowerCI <- dfnew$fit - (1.96*dfnew$se)
dfnew$upperCI <- dfnew$fit + (1.96*dfnew$se)
### Standardisation parameters for HbA1c
mean_hba1c <-  std_params[1,1]
sd_hba1c <- std_params[1,2]
### Use these to create unstandardised outcomes and CIs
dfnew$fit_unstd <- (sd_hba1c*dfnew$fit) + mean_hba1c
#dfnew$lCIunstd <- (sd_hba1c*dfnew$lowerCI) + mean_hba1c
#dfnew$uCIunstd <- (sd_hba1c*dfnew$upperCI) + mean_hba1c
dfnew$lCIunstd <- dfnew$fit_unstd - 1.96*sd_hba1c*dfnew$se
dfnew$uCIunstd <- dfnew$fit_unstd + 1.96*sd_hba1c*dfnew$se
############################################################################
### Use expected outcome and 95% CIs for example individuals to create
### plots for paper
############################################################################
### All together (too busy)
pdf(paste(analysis_dir, "plots/unadj_example_plots_all.pdf", sep=""))
ggplot(data = dfnew, aes(x = time_base, y = fit, group = Example, color=Example)) +
geom_ribbon(aes(ymin = lowerCI, ymax = upperCI, color=Example, group=Example, fill = Example), alpha = 0.25, linetype="dashed") +
  geom_line() + facet_wrap(~Example) +
  labs(y = "HbA1c (standardised)", x="T2D disease duration (years)")

ggplot(data = dfnew, aes(x = time_base, y = fit_unstd, group = Example, color=Example)) + 
geom_ribbon(aes(ymin = lCIunstd, ymax = uCIunstd, color=Example, group=Example, fill = Example), alpha = 0.25, linetype="dashed") +
  geom_line() + facet_wrap(~Example)  +
  labs(y = "HbA1c (mmol/mol)", x="T2D disease duration (years)") 
dev.off()

### Reduce to a pre-T2D MDD + no MDD dataset and
### a post-T2D MDD + no MDD dataset
dfpre <- dfnew[dfnew$dep_group != "Post-T2D MDD", ]
dfpost <- dfnew[dfnew$dep_group != "Pre-T2D MDD", ]
### Plots for these datasets separately:
pdf(paste(analysis_dir, "plots/unadj_preT2D_noMDD_plots.pdf", sep=""))

ggplot(data = dfpre, aes(x = time_base, y = fit_unstd, color=Example)) + 
geom_ribbon(aes(ymin = lCIunstd, ymax = uCIunstd, color=Example, fill=Example), alpha = 0.25, linetype="dashed") +
  geom_line() + facet_wrap(~Example) + theme(legend.position = "none")  +
  labs(y = "HbA1c (mmol/mol)", x="T2D disease duration (years)") 
dev.off()

pdf(paste(analysis_dir, "plots/unadj_postT2D_noMDD_plots.pdf", sep=""))

ggplot(data = dfpost, aes(x = time_base, y = fit_unstd, color=Example)) + 
geom_ribbon(aes(ymin = lCIunstd, ymax = uCIunstd, color=Example, fill=Example), alpha = 0.25, linetype="dashed") +
  geom_line() + facet_wrap(~Example)+ theme(legend.position = "none")  +
  labs(y = "HbA1c (mmol/mol)", x="T2D disease duration (years)") 
dev.off()
############################################################################
### Additional variations tried for plots
#nomdd <- dfpre[dfpre$dep_group == "No MDD", ]
#dfprepost <- dfnew[dfnew$dep_group != "No MDD", ]
#
#pdf(paste(analysis_dir, "plots/unadj_noMDD_plotsv2.pdf", sep=""))
#ggplot(data = nomdd, aes(x = time_base, y = fit_unstd, color=Example)) +
#geom_ribbon(aes(ymin = lCIunstd, ymax = uCIunstd, color=Example, fill=Example), alpha = 0.25, linetype="dashed") +
#  geom_line() + facet_wrap(~Example)+ theme(legend.position = "none")  +
#  labs(y = "HbA1c (mmol/mol)", x="T2D disease duration (years)")
#dev.off()

#tmp <- strsplit(dfprepost$Example, split= ":")
#Example2 <- NULL
#for(i in 1:dim(dfprepost)[1]){
#  outi <- tmp[[i]]
#  Example2 <- c(Example2, outi[2])
#}
#dfprepost$Example2 <- Example2
#dfprepost$Example3 <- factor(Example2, levels=c(" 2.2 years before T2D", " 10.7 years before T2D", " 27.5 years before T2D",
#  " 1 year after T2D", " 3.8 years after T2D", " 7.1 years after T2D"))
#dfprepost$dep_group2 <- factor(dfprepost$dep_group, levels=c("Pre-T2D MDD", "Post-T2D MDD"))

#pdf(paste(analysis_dir, "plots/unadj_dfpre_postMDD_plotsv2.pdf", sep=""))
#ggplot(data = dfprepost, aes(x = time_base, y = fit_unstd, color=Example)) +
#geom_ribbon(aes(ymin = lCIunstd, ymax = uCIunstd, color=Example, fill=Example), alpha = 0.25, linetype="dashed") +
 # geom_line() + facet_wrap(dep_group2 ~ Example3)+ theme(legend.position = "none")  +
 # labs(y = "HbA1c (mmol/mol)", x="T2D disease duration (years)")
#dev.off()

#tmp <- strsplit(dfnew$Example, split= ":")
#Example2 <- NULL
#for(i in 1:dim(dfnew)[1]){
#  outi <- tmp[[i]]
#  Example2 <- c(Example2, outi[2])
#}
#Example2[is.na(Example2)] <- " "
#dfnew$Example2 <- Example2
#dfnew$Example3 <- factor(Example2, levels=c(" 2.2 years before T2D", " 10.7 years before T2D", " 27.5 years before T2D",
#  " 1 year after T2D", " 3.8 years after T2D", " 7.1 years after T2D", " "))
#dfnew$dep_group2 <- factor(dfnew$dep_group, levels=c("Pre-T2D MDD", "Post-T2D MDD", "No MDD"))

#pdf(paste(analysis_dir, "plots/unadj_all_plotsv2.pdf", sep=""))
#ggplot(data = dfnew, aes(x = time_base, y = fit_unstd, color=Example)) +
#geom_ribbon(aes(ymin = lCIunstd, ymax = uCIunstd, color=Example, fill=Example), alpha = 0.25, linetype="dashed") +
 # geom_line() + facet_wrap(dep_group2 ~ Example3)+ theme(legend.position = "none")  +
#  labs(y = "HbA1c (mmol/mol)", x="T2D disease duration (years)")
#dev.off()

############################################################################
### 2.3. Summarising model output using example individuals- tables
############################################################################
### Tables: pre-T2D MDD examples
############################################################################
### HbA1c at time points for pre-T2D MDD individuals and no MDD individual
############################################################################
dfpre <- dfnew[dfnew$dep_group != "Post-T2D MDD", ]
dfpost <- dfnew[dfnew$dep_group != "Pre-T2D MDD", ]
library(tidyverse)
dfpre_tab <- dfpre[round(dfpre$time_base, 5) %in% round(time_base_tab, 5), ]
dfpre_tab <- dfpre_tab %>% select(time_base, Example, fit_unstd, lCIunstd, uCIunstd)
dfpre_tab$time_base <- round(dfpre_tab$time_base, 1)
dfpre_tab$HbA1c <- paste(as.character(round(dfpre_tab$fit_unstd, 2)), 
paste(paste("(", as.character(round(dfpre_tab$lCIunstd, 2)), sep=""), paste(as.character(round(dfpre_tab$uCIunstd, 2)), ")", sep=""),sep=", "),sep=" ")
dfpre_tab2 <- dfpre_tab %>% select(time_base, Example, HbA1c)
dfpre_tab2 <- reshape(dfpre_tab2, idvar = "time_base", timevar = "Example", direction = "wide")
### Write to output
write.csv(dfpre_tab2, file= paste(analysis_dir, "desc_tabs/unadjusted_pre_t2d_tab.csv", sep=""))
############################################################################
### Difference in HbA1c (from model) between pre-T2D MDD examples and no MDD
############################################################################
dfpre_tab$temp1 <- rep(dfpre_tab$fit_unstd[dfpre_tab$Example == "No MDD"], 4)
dfpre_tab$diff <- dfpre_tab$fit_unstd - dfpre_tab$temp1

dfpre_tab$duration <- c(rep(0, length(time_base_tab)), rep(lower_pretime, length(time_base_tab)), 
rep(m_pretime , length(time_base_tab)), rep(upper_pretime, length(time_base_tab)))

### Using parameters extracted from model:
sebeta1 <- 0.02702733*sd_hba1c
sebeta_dur <- 0.02098304*sd_hba1c
corr_pre <- -0.771
cov_pre <- corr_pre*sebeta1*sebeta_dur
dfpre_tab$se_diff <- sqrt((sebeta1^2) + ((dfpre_tab$duration^2)*(sebeta_dur^2)) + (2*dfpre_tab$duration*cov_pre))
dfpre_tab$lower_diff <- dfpre_tab$diff - 1.96*dfpre_tab$se_diff
dfpre_tab$upper_diff <- dfpre_tab$diff + 1.96*dfpre_tab$se_diff

0.0530728/-0.0665073
############################################################################
### Tables: post-T2D MDD examples
############################################################################
### HbA1c at time points for post-T2D MDD individuals and no MDD individual
############################################################################
dfpost_tab <- dfpost[round(dfpost$time_base, 5) %in% round(time_base_tab, 5), ]
dfpost_tab <- dfpost_tab %>% select(time_base, Example, fit_unstd, lCIunstd, uCIunstd)
dfpost_tab$time_base <- round(dfpost_tab$time_base, 1)
dfpost_tab$HbA1c <- paste(as.character(round(dfpost_tab$fit_unstd, 2)), 
paste(paste("(", as.character(round(dfpost_tab$lCIunstd, 2)), sep=""), paste(as.character(round(dfpost_tab$uCIunstd, 2)), ")", sep=""),sep=", "),sep=" ")
dfpost_tab2 <- dfpost_tab %>% select(time_base, Example, HbA1c)
dfpost_tab2 <- reshape(dfpost_tab2, idvar = "time_base", timevar = "Example", direction = "wide")
### Write to output
write.csv(dfpost_tab2, file= file= paste(analysis_dir, "desc_tabs/unadjusted_post_t2d_mdd_tab.csv", sep=""))

############################################################################
### Difference in HbA1c (from model) between post-T2D MDD examples and no MDD
############################################################################
dfpost_tab <- dfpost[round(dfpost$time_base, 5) %in% round(time_base_tab, 5), ]
dfpost_tab <- dfpost_tab %>% select(time_base, Example, fit_unstd, lCIunstd, uCIunstd)
dfpost_tab$time_base <- round(dfpost_tab$time_base, 1)
dfpost_tab$temp1 <- rep(dfpost_tab$fit_unstd[dfpost_tab$Example == "No MDD"], 4)
dfpost_tab$diff <- dfpost_tab$fit_unstd - dfpost_tab$temp1

time_base_tab1 <- time_base_tab - lower_posttime
time_base_tab1[time_base_tab1 < 0] <- 0

time_base_tab2 <- time_base_tab - m_posttime
time_base_tab2[time_base_tab2 < 0] <- 0

time_base_tab3 <- time_base_tab - upper_posttime
time_base_tab3[time_base_tab3 < 0] <- 0

dfpost_tab$duration <- c(rep(0, length(time_base_tab)), time_base_tab1, 
time_base_tab2, time_base_tab3)

sebeta1post <- 0.03325739*sd_hba1c
sebeta_postdur <- 0.00988531*sd_hba1c
cov_post <- -1.597178e-04*sd_hba1c*sd_hba1c

dfpost_tab$se_diff <- sqrt((sebeta1post^2) + ((dfpost_tab$duration^2)*(sebeta_postdur^2)) + (2*dfpost_tab$duration*cov_post))
dfpost_tab$lower_diff <- dfpost_tab$diff - 1.96*dfpost_tab$se_diff
dfpost_tab$upper_diff <- dfpost_tab$diff + 1.96*dfpost_tab$se_diff

diff_out <- dfpost_tab[dfpost_tab$Example != "No MDD", ]
diff_out <- diff_out %>% select(time_base, Example, diff, lower_diff, upper_diff)
diff_out$Difference <- paste(as.character(round(diff_out$diff, 2)), 
paste(paste("(", as.character(round(diff_out$lower_diff, 2)), sep=""), paste(as.character(round(diff_out$upper_diff, 2)), ")", sep=""),sep=", "),sep=" ")
diff_out2 <- diff_out %>% select(time_base, Example, Difference)
diff_out2 <- reshape(diff_out2, idvar = "time_base", timevar = "Example", direction = "wide")
### Write difference to output for post-T2D MDD group
write.csv(diff_out2, file= file= file= paste(analysis_dir, "desc_tabs/unadjusted_post_t2d_mdd_DIFFERENCEtab.csv", sep=""))


### Calculate the time at which we expect the 95% CI for the difference between
### mean HbA1c for post-T2D MDD individuals and no MDD to exclude 0
quad <- function(a, b, c)
{
  a <- as.complex(a)
  answer <- c((-b + sqrt(b^2 - 4 * a * c)) / (2 * a),
              (-b - sqrt(b^2 - 4 * a * c)) / (2 * a))
  if(all(Im(answer) == 0)) answer <- Re(answer)
  if(answer[1] == answer[2]) return(answer[1])
  answer
}

A <- ((0.0304995*sd_hba1c)^2) - ((1.96*sebeta_postdur)^2)
B <- (2*((0.0304995*sd_hba1c)*(-0.0285279*sd_hba1c))) - (2*(1.96^2)*cov_post)
C <- ((-0.0285279*sd_hba1c)^2) - ((1.96*sebeta1post)^2)
quad(a=A, b=B, c=C)
### 2.993200 -2.068524
### Approx 3 years after diagnosis we expect the 95% CI for difference between post-T2D MDD and no MDD to exclude 0.
