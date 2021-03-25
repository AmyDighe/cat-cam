# evaluating predictive capacity and comparing model fit

# read in the fits

# simulated - compare fit with and without correct model assumptions
fit1_sim1 <- readRDS("fits/fit_4_1.rds")
fit4_sim4 <- readRDS("fits/fit_4.rds")
fit4_sim1 <- readRDS("fits/fit_4_sim1.rds")
fit1_sim4 <- readRDS("fits/fit_1_sim4.rds")

DIC_fit1_sim4 <- DICfunc(fit1_sim4)
DIC_fit4_sim1 <- DICfunc(fit4_sim1)

# real data
fit1_real <- readRDS("fits/fit_4_1real.rds")
fit2_real <- readRDS("fits/fit_4_2real.rds")
fit3_real <- readRDS("fits/fit_4_3real.rds")
fit4_real <- readRDS("fits/fit_4real.rds")
fit_1real_10 <- readRDS("fits/fit_4_1real_10.rds")
fit_2real_10 <- readRDS("fits/fit_4_2real_10.rds")
fit_3real_10 <- readRDS("fits/fit_4_3real_10.rds")
fit_4real_10 <- readRDS("fits/fit_4real_10.rds")

DIC_fit1_real <- DICfunc(fit = fit1_real)
DIC_fit2_real <- DICfunc(fit = fit2_real)
DIC_fit3_real <- DICfunc(fit = fit3_real)
DIC_fit3bb_real <- DICfunc(fit = fit_4_3bbreal)
DIC_fit4_real <- DICfunc(fit = fit4_real)
DIC_fit1_real_10 <- DICfunc(fit = fit_1real_10)
DIC_fit2_real_10 <- DICfunc(fit = fit_2real_10)
DIC_fit3_real_10 <- DICfunc(fit = fit_3real_10)
DIC_fit4_real_10 <- DICfunc(fit = fit_4real_10)
DIC_fit5_real_10 <- DICfunc(fit = fit_5real_10)
DIC_fit5_real_10_studspec <- DICfunc(fit = fit_5real_10_studspec)

# # extract the log likelihood from the stanfit object

WAICfunc <- function(fit){
  
ll_1s <- extract_log_lik(fit) 
ll_1s <-  ll_1s[, colSums(ll_1s != 0) > 0] # because of ragged array padding 

# WAIC and LOO calculations
WAIC1 <- waic(ll_1s)
LOO1 <- loo(ll_1s, moment_match = TRUE)
LOO_07 <- loo(ll_1s, k_threshold = 0.7, moment_match = TRUE)
return(list(waic = WAIC1,
            loo = LOO1,
            loo_07 = LOO_07))
}

WAIC1 <- WAICfunc(fit = fit1_real)
WAIC2 <-WAICfunc(fit = fit2_real)
WAIC3 <-WAICfunc(fit = fit3_real)
WAIC4 <-WAICfunc(fit = fit4_real)

WAIC1 <-WAICfunc(fit = fit1_sim1)
WAIC1 <-WAICfunc(fit = fit1_sim4)
WAIC1 <-WAICfunc(fit = fit4_sim1)
WAIC1 <-WAICfunc(fit = fit4_sim4)

compare(WAIC1$waic, WAIC2$waic)
compare(WAIC1$waic, WAIC3$waic)
compare(WAIC1$waic, WAIC4$waic)
compare(WAIC1$waic, WAIC2$waic)
compare(WAIC1$waic, WAIC2$waic)

##loo_moment_match(x = fit_4_1, loo = LOO1) #not working - need to pass user defined functions

## LOO with datasets with > 2 age classes

WAIC1 <- WAICfunc(fit = fit_1real_10)
WAIC2 <-WAICfunc(fit = fit_2real_10)
WAIC3 <-WAICfunc(fit = fit_3real_10)
WAIC4 <-WAICfunc(fit = fit_4real_10)
