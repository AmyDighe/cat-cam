#########################
## FITTING TO THE DATA ##
#########################

n_datasets <- 10
n_ages <- 5

LOW_AGE <- matrix(c(0, 0.51, 3, 6, 10,
                    0, 1, 3, 4, 10,
                    0.5, 1.5, 2.5, 3.5, 4.5,
                    0, 0.51, 3, 6, 10,
                    0, 1, 3, 4, 10,
                    0, 2, 10, 10, 10,
                    0, 0.51, 3, 6, 10,
                    0, 1, 3, 4, 10,
                    0.5, 1.5, 2.5, 3.5, 4.5,
                    0, 2, 10, 10, 10), 
                  ncol = n_ages, nrow = n_datasets, byrow = T)
UPP_AGE <- matrix(c(0.5, 2.99, 5.99, 8, 20,
                    0.99, 2.99, 3.99, 5, 20,
                    1.49, 2.49, 3.49, 4.49, 5.49,
                    0.5, 2.99, 5.99, 8, 20,
                    0.99, 2.99, 3.99, 5, 20,
                    1.99, 10, 20, 20, 20,
                    0.5, 2.99, 5.99, 8, 20,
                    0.99, 2.99, 3.99, 5, 20,
                    1.49, 2.49, 3.49, 4.49, 5.49,
                    1.99, 10, 20, 20, 20),
                  ncol = n_ages, nrow = n_datasets, byrow = T)

# define number of camels per age class manually
# (could just use the ones in the real data?!)
N_CAMELS <- matrix(c(2000, 2000, 2000, 2000, 0,
                     2000, 2000, 2000, 2000, 0,
                     2000, 2000, 2000, 2000, 2000,
                     2000, 2000, 2000, 2000, 0,
                     2000, 2000, 2000, 2000, 0,
                     2000, 2000, 0, 0, 0,
                     2000, 2000, 2000, 2000, 0,
                     2000, 2000, 2000, 2000, 0,
                     2000, 2000, 2000, 2000, 2000,
                     2000, 2000, 0, 0, 0), ncol = n_ages, nrow = n_datasets, byrow = T)

##############
## BINOMIAL ##
##############

# run full model reduced to 1

fit_4_1 <- stan(
  file = here::here("stan-models/model4_reduced1b.stan"),
  data = list(
    S = nrow(simk0[[1]][[1]]$simulated),
    A =  ncol(simk0[[1]][[1]]$simulated),
    N = N_CAMELS,
    pos = simk0[[1]][[1]]$simulated,
    age1 = LOW_AGE,
    age2 = UPP_AGE,
    sigma_m = 2.1,
    sigma_r = 0,
    sens = 0.999,
    spec = 1,
    mabs = -1
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)


fit_4_2 <- stan(
  file = here::here("stan-models/model4_reduced2b.stan"),
  data = list(
    S = nrow(simk0[[2]][[1]]$simulated),
    A =  ncol(simk0[[2]][[1]]$simulated),
    N = N_CAMELS,
    pos = simk0[[2]][[1]]$simulated,
    age1 = LOW_AGE,
    age2 = UPP_AGE,
    sigma_m = 2.1,
    sens = 0.999,
    spec = 1,
    mabs = -1
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99) 
)


fit_4_3 <- stan(
  file = here::here("stan-models/model4_reduced3b.stan"),
  data = list(
    S = nrow(simk0[[3]][[1]]$simulated),
    A =  ncol(simk0[[3]][[1]]$simulated),
    N = N_CAMELS,
    pos = simk0[[3]][[1]]$simulated,
    age1 = LOW_AGE,
    age2 = UPP_AGE,
    sigma_r = 0,
    sens = sens,
    spec = spec,
    mabs = 1
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE,
  ##control = list(max_treedepth = 15)
  ##control = list(adapt_delta = 0.99)
)

fit_4 <- stan(
  file = here::here("stan-models/model4b.stan"),
  data = list(
    S = nrow(simk0[[4]][[1]]$simulated),
    A =  ncol(simk0[[4]][[1]]$simulated),
    N = N_CAMELS,
    pos = simk0[[4]][[1]]$simulated,
    age1 = LOW_AGE,
    age2 = UPP_AGE,
    sens = sens,
    spec = spec, 
    mabs = 1
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)

saveRDS(fit_4_1, "fits/fit_4_1.rds")
saveRDS(fit_4_2, "fits/fit_4_2.rds")
saveRDS(fit_4_3, "fits/fit_4_3.rds")
saveRDS(fit_4, "fits/fit_4.rds")


##################
## BETABINOMIAL ##
##################

# run full model reduced to 1

fit_4_1bb <- stan(
  file = here::here("stan-models/model4_reduced1bb.stan"),
  data = list(
    S = nrow(simk001[[1]][[1]]$simulated),
    A =  ncol(simk001[[1]][[1]]$simulated),
    N = N_CAMELS,
    pos = simk001[[1]][[1]]$simulated,
    age1 = LOW_AGE,
    age2 = UPP_AGE,
    sigma_m = 2.1,
    sigma_r = 0,
    sens = 0.999,
    spec = 1,
    mabs = -1
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)

diagnos <- ggmcmc(ggs(fit_4_1bb), here::here("diagnostics/1bb.pdf"))

fit_4_2bb <- stan(
  file = here::here("stan-models/model4_reduced2bb.stan"),
  data = list(
    S = nrow(simk001[[2]][[1]]$simulated),
    A =  ncol(simk001[[2]][[1]]$simulated),
    N = N_CAMELS,
    pos = simk001[[2]][[1]]$simulated,
    age1 = LOW_AGE,
    age2 = UPP_AGE,
    sigma_m = 2.1,
    sens = 0.999,
    spec = 1,
    mabs = -1
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99) 
)

diagnos <- ggmcmc(ggs(fit_4_2bb), here::here("diagnostics/2bb.pdf"))

fit_4_3bb <- stan(
  file = here::here("stan-models/model4_reduced3bb.stan"),
  data = list(
    S = nrow(simk001[[3]][[1]]$simulated),
    A =  ncol(simk001[[3]][[1]]$simulated),
    N = N_CAMELS,
    pos = simk001[[3]][[1]]$simulated,
    age1 = LOW_AGE,
    age2 = UPP_AGE,
    sigma_r = 0,
    sens = sens,
    spec = spec,
    mabs = 1
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE,
  ##control = list(max_treedepth = 15)
  ##control = list(adapt_delta = 0.99)
)

diagnos <- ggmcmc(ggs(fit_4_3bb), here::here("diagnostics/3bb.pdf"))

fit_4bb <- stan(
  file = here::here("stan-models/model4bb.stan"),
  data = list(
    S = nrow(simk001[[4]][[1]]$simulated),
    A =  ncol(simk001[[4]][[1]]$simulated),
    N = N_CAMELS,
    pos = simk001[[4]][[1]]$simulated,
    age1 = LOW_AGE,
    age2 = UPP_AGE,
    sens = sens,
    spec = spec, 
    mabs = 1
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)
diagnos <- ggmcmc(ggs(fit_4bb), here::here("diagnostics/4bb.pdf"))

saveRDS(fit_4_1bb, "fits/fit_4_1bb.rds")
saveRDS(fit_4_2bb, "fits/fit_4_2bb.rds")
saveRDS(fit_4_3bb, "fits/fit_4_3bb.rds")
saveRDS(fit_4bb, "fits/fit_4bb.rds")

fit_4bb <- readRDS("fits/fit_4bb.rds")




## with k = 0.1

# run full model reduced to 1

fit_4_1bb1 <- stan(
  file = here::here("stan-models/model4_reduced1bb.stan"),
  data = list(
    S = nrow(simk01[[1]][[1]]$simulated),
    A =  ncol(simk01[[1]][[1]]$simulated),
    N = N_CAMELS,
    pos = simk01[[1]][[1]]$simulated,
    age1 = LOW_AGE,
    age2 = UPP_AGE,
    sigma_m = 2.1,
    sigma_r = 0,
    sens = 0.999,
    spec = 1,
    mabs = -1
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)

diagnos <- ggmcmc(ggs(fit_4_1bb1), here::here("diagnostics/1bb1.pdf"))

fit_4_2bb1 <- stan(
  file = here::here("stan-models/model4_reduced2bb.stan"),
  data = list(
    S = nrow(simk01[[2]][[1]]$simulated),
    A =  ncol(simk01[[2]][[1]]$simulated),
    N = N_CAMELS,
    pos = simk01[[2]][[1]]$simulated,
    age1 = LOW_AGE,
    age2 = UPP_AGE,
    sigma_m = 2.1,
    sens = 0.999,
    spec = 1,
    mabs = -1,
    k = 0.1
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99) 
)

diagnos <- ggmcmc(ggs(fit_4_2bb1), here::here("diagnostics/2bb1.pdf"))

fit_4_3bb1 <- stan(
  file = here::here("stan-models/model4_reduced3bb.stan"),
  data = list(
    S = nrow(simk01[[3]][[1]]$simulated),
    A =  ncol(simk01[[3]][[1]]$simulated),
    N = N_CAMELS,
    pos = simk01[[3]][[1]]$simulated,
    age1 = LOW_AGE,
    age2 = UPP_AGE,
    sigma_r = 0,
    sens = sens,
    spec = spec,
    mabs = 1
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE,
  ##control = list(max_treedepth = 15)
  ##control = list(adapt_delta = 0.99)
)

diagnos <- ggmcmc(ggs(fit_4_3bb1), here::here("diagnostics/3bb1.pdf"))

fit_4bb1 <- stan(
  file = here::here("stan-models/model4bb.stan"),
  data = list(
    S = nrow(simk01[[4]][[1]]$simulated),
    A =  ncol(simk01[[4]][[1]]$simulated),
    N = N_CAMELS,
    pos = simk01[[4]][[1]]$simulated,
    age1 = LOW_AGE,
    age2 = UPP_AGE,
    sens = sens,
    spec = spec, 
    mabs = 1
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)
diagnos <- ggmcmc(ggs(fit_4bb1), here::here("diagnostics/4bb1.pdf"))

saveRDS(fit_4_1bb1, "fits/fit_4_1bb1.rds")
saveRDS(fit_4_2bb1, "fits/fit_4_2bb1.rds")
saveRDS(fit_4_3bb1, "fits/fit_4_3bb1.rds")
saveRDS(fit_4bb1, "fits/fit_4bb1.rds")

fit_4bb <- readRDS("fits/fit_4bb.rds")