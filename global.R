# stan
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# dependencies
library(ggmcmc)
library(emdbook)
library(dplyr)
library(gridExtra)
library(binom)
library(bayesplot)
library(loo)
source("utils.R")

# simulated data
foi <- c(0.5, 0.25, 0.05, 0.75, 1, 0.1, 1.5, 2, 0.2, 0.6)
datak0 <- readRDS("data/sim_datak0") # contain data frames for plotting
datak001 <- readRDS("data/sim_datak001")
datak01 <- readRDS("data/sim_datak01")
simk0 <- readRDS("data/sim_matk0") #contain STUDY X AGE matrices for fit
simk001 <- readRDS("data/sim_matk001")
simk01 <- readRDS("data/sim_matk01")
default_sens <- 0.999
default_spec <- 1

# real data
data_sero <- readRDS("data/data_sero.rds") # contain data frames for plotting
SEROPOS <- readRDS("data/SEROPOS.rds") #contain STUDY X AGE matrices for fit
AGE_L <- readRDS("data/AGE_L.rds")
AGE_U <- readRDS("data/AGE_U.rds")
N_CAMELS <- readRDS("data/N_CAMELS.rds")
STUDY_TEST_TYPE <- readRDS("data/STUDY_TEST_TYPE.rds")
