rm(list = ls())
library(tidyverse)
library(knitr)
library(mice)
library(miceadds)
library(mitml)
library(metafor)
library(here)
set.seed(20261001)

# Input
raw <- readRDS(here("data", "SeeData.rds"))

## Normal Approximation of the Log-Hazard Ratio
df <- mutate(raw, 
             diff = (log(ULHR) + log(LLHR)) / 2,
             se = (log(ULHR) - log(LLHR)) / (1.96*2) ,
             weight = 1/se**2) %>%
  rename(estimate = diff, std.error = se)
##  Scale
df <- mutate(df,
             eGFR=eGFR*0.1#, Stage3P=Stage3P*0.01, MaleP=MaleP*0.01, DMP=DMP*0.01, CHDP=CHDP*0.01
             )

## Age: Auxiliary variable; Used only in imputation model
indata <- df %>%
  dplyr::select(Strata, Stage3P, Age, eGFR, estimate, std.error, weight)
#  dplyr::select(Strata, FUYear, Stage3P, Age, MaleP, DMP, CHDP, eGFR, estimate, std.error, weight)
meth1 <- make.method(indata) #PMM is specified by default.
pred1 <- make.predictorMatrix(indata)
## Specify variables NOT to be used as explanatory variables in the imputation model
pred0var = c("std.error", "weight")
#pred0var = c("std.error")
macov <- "~ eGFR + Strata + Stage3P"
f_cov <- as.formula(macov)

pred1[, pred0var] <- 0 
mi_n <- 100
meth1[ meth1=="pmm" ] <- "weighted.norm"
miceseed <- 12321

# Inverse-variance weight
weights <- indata$weight

# CCA
res_cca <- rma(yi = estimate, sei = std.error, data = indata, method = "REML", 
           mods = f_cov,
           test = "knha", control=list(maxiter=10000, stepadj=0.5))
# Sensitivity analysis weight
weights_sa <- 1/(indata$std.error**2 + res_cca$tau2)

# Imputation
imp1 <- mice(data = indata, m = mi_n, predictorMatrix = pred1, method = meth1, imputationWeights=weights, seed = miceseed, printFlag = FALSE)
imp1
imp.list <- complete(imp1, "all")
# Analysis 
fit.imp1 <- lapply(imp.list, function(d) {
  rma(yi = estimate, sei = std.error, data = d, method = "REML",
      mods = f_cov,
      test = "knha", control=list(maxiter=10000, stepadj=0.5))
})
# Pooling
## Original function
metareg_pool <- function(anal.list, df, names){
  coef_list <- lapply(anal.list, function(d) {
    as.matrix(d$beta)
  })
  vcov_list <- lapply(anal.list, function(d) {
    as.matrix(d$vb)
  })
  rubin_pool <- testEstimates(qhat = coef_list, uhat = vcov_list, df.com = df)
  confint <- confint(rubin_pool)
  rubin_tab <- as.data.frame(rubin_pool$estimates) %>%
    cbind(as.data.frame(confint))
  rownames(rubin_tab) <- names
  rubin_tab
}

mr.df <- fit.imp1[[1]]$ddf
param_names <- rownames(fit.imp1[[1]]$beta)
result <- metareg_pool(fit.imp1, df = mr.df, names = param_names)
kable(result, digits = 4)

# Simple average for tau2
tau2_vals <- sapply(fit.imp1, `[[`, "tau2")
# tau2_vals
mean_tau2 <- mean(tau2_vals, na.rm = TRUE) 
sqrt(mean_tau2)

sessionInfo()
