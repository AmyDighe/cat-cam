library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

library(ggmcmc)
library(emdbook)
library(dplyr)
source("utils.R")
