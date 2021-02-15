# fit the models to the real data
init_fun <- function(...) list(foi = c(rep(2.911526, nrow(SEROPOS))))


fit_4_1real <- stan(
  file = here::here("stan-models/model4_reduced1b.stan"),
  data = list(
    S = nrow(SEROPOS),
    A =  ncol(SEROPOS),
    N = N_CAMELS,
    pos = SEROPOS,
    age1 = AGE_L,
    age2 = AGE_U,
    M = M_ZERO_ZERO,
    sigma_m = 12,
    sigma_r = 0
  ),
  chains = 1,
  iter = 4000,
  verbose = TRUE,
  init = init_fun
  ##control = list(adapt_delta = 0.99)
)

rstan::summary(fit_4_1real)
diagnos <- ggmcmc(ggs(fit_4_1real), here::here("diagnostics/real4_1b.pdf"))


fit_4_2real <- stan(
  file = here::here("stan-models/model4_reduced2b.stan"),
  data = list(
    S = nrow(SEROPOS),
    A =  ncol(SEROPOS),
    N = N_CAMELS,
    pos = SEROPOS,
    age1 = AGE_L,
    age2 = AGE_U,
    M = M_ZERO_ZERO,
    sigma_m = 12
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99) 
)

fit2_summary <- rstan::summary(fit_4_2real)
 diagnos <- ggmcmc(ggs(fit_4_2real), here::here("diagnostics/real4_2b.pdf"))

fit_4_3real <- stan(
  file = here::here("stan-models/model4_reduced3b.stan"),
  data = list(
    S = nrow(SEROPOS),
    A =  ncol(SEROPOS),
    N = N_CAMELS,
    pos = SEROPOS,
    age1 = AGE_L,
    age2 = AGE_U,
    M = M_ZERO,
    sigma_r = 0.00001
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE,
  ##control = list(max_treedepth = 15)
  ##control = list(adapt_delta = 0.99)
)

rstan::summary(fit_4_3real)

fit_4real <- stan(
  file = here::here("stan-models/model4b.stan"),
  data = list(
    S = nrow(SEROPOS),
    A =  ncol(SEROPOS),
    N = N_CAMELS,
    pos = SEROPOS,
    age1 = AGE_L,
    age2 = AGE_U,
    M = M_ZERO
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)

fit4_summary <- rstan::summary(fit_4real)
diagnos <- ggmcmc(ggs(fit_4real), here::here("diagnostics/real4b.pdf"))


# run full model reduced to 1

# test on sim

init_fun <- function(...) list(foi = c(rep(2.911526, n_datasets)))

fit_4_1 <- stan(
  file = here::here("stan-models/model4_reduced1b.stan"),
  data = list(
    S = nrow(simk0[[1]][[1]]$simulated),
    A =  ncol(simk0[[1]][[1]]$simulated),
    N = N_camels,
    pos = simk0[[1]][[1]]$simulated,
    age1 = age_lower,
    age2 = age_upper,
    M = simk0[[1]][[1]]$M_initial,
    sigma_m = 12, #set it greater than max foi to avoid lack of support
    sigma_r = 0
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE,
  init = init_fun
  ##control = list(adapt_delta = 0.99)
)

rstan::summary(fit_4_1)
get_inits(fit_4_1)


fit_4_2 <- stan(
  file = here::here("stan-models/model4_reduced2b.stan"),
  data = list(
    S = nrow(simk0[[2]][[1]]$simulated),
    A =  ncol(simk0[[2]][[1]]$simulated),
    N = N_camels,
    pos = simk0[[2]][[1]]$simulated,
    age1 = age_lower,
    age2 = age_upper,
    M = simk0[[2]][[1]]$M_initial,
    sigma_m = 12
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
    N = N_camels,
    pos = simk0[[3]][[1]]$simulated,
    age1 = age_lower,
    age2 = age_upper,
    M = simk0[[3]][[1]]$M_initial,
    sigma_r = 0
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE,
  ##control = list(max_treedepth = 15)
  ##control = list(adapt_delta = 0.99)
)

rstan::summary(fit_model4av)

fit_model4av <- stan(
  file = here::here("stan-models/model4b.stan"),
  data = list(
    S = nrow(simk0[[4]][[1]]$simulated),
    A =  ncol(simk0[[4]][[1]]$simulated),
    N = N_camels,
    pos = simk0[[4]][[1]]$simulated,
    age1 = age_lower,
    age2 = age_upper,
    M = simk0[[4]][[1]]$M_initial
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)

rstan::summary(fit_model4av)

saveRDS(fit_4_2real, file = "fits/fit_4_2real.rds")
saveRDS(fit_4real, file = "fits/fit_4real.rds")
