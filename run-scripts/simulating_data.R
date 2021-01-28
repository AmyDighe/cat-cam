# define age classes manually 
# (could just use the ones in the real data?!)
n_ages <- 5
n_datasets <- 3
age_lower <- matrix(c(0, 0.51, 3, 6, 10,
                      0, 1, 3, 4, 10,
                      0.5, 1.5, 2.5, 3.5, 4.5), 
                    ncol = n_ages, nrow = n_datasets, byrow = T)
age_upper <- matrix(c(0.5, 2.99, 5.99, 8, 20,
                      0.99, 2.99, 3.99, 5, 20,
                      1.49, 2.49, 3.49, 4.49, 5.49),
                    ncol = n_ages, nrow = n_datasets, byrow = T)

# define number of camels per age class manually
# (could just use the ones in the real data?!)
N_camels <- matrix(c(2000, 2000, 2000, 2000, 0,
                     2000, 2000, 2000, 2000, 0,
                     2000, 2000, 2000, 2000, 2000), ncol = n_ages, nrow = n_datasets, byrow = T)
od <- 0.1

sim4 <- sim_data_betabinom(n_datasets, n_ages,
         gamma = c(0.5, 0.25, 0.05), sigma = 0.2, omega = 2, mabs = 1,
         N_camels, age_upper, age_lower, od = od)

sim3 <- sim_data_betabinom(n_datasets, n_ages,
         gamma = c(0.5, 0.25, 0.05), sigma = 0, omega = 2, mabs = 1,
         N_camels, age_upper, age_lower, od = od)

sim2 <- sim_data_betabinom(n_datasets, n_ages,
         gamma = c(0.5, 0.25, 0.05), sigma = 0.2, omega = 2, mabs = 0,
         N_camels, age_upper, age_lower, od = od)

sim1 <- sim_data_betabinom(n_datasets, n_ages,
                 gamma = c(0.5, 0.25, 0.05), sigma = 0, omega = 2, mabs = 0,
                 N_camels, age_upper, age_lower, od = od)

# checks
# 1. if mabs = 1 pmAb and M_intial should be non-zero, if mabs = 0, pmab and M_initial should all be zero

sim2$M_initial
sim2$pmAbs

# 2. these two should be similar for those with >0 in N_camels if od is close to zero

sim2$simulated/N_camels
sim2$pmAbs + sim4$pprev

# run the full model 


fit_model4av <- stan(
  file = here::here("stan-models/model4av_hierarchical.stan"),
  data = list(
    S = nrow(sim4$simulated),
    A =  ncol(sim4$simulated),
    N = N_camels,
    pos = sim4$simulated,
    age1 = age_lower,
    age2 = age_upper,
    M = sim4$M_initial
  ),
  chains = 2,
  iter = 4000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)


# run the full model reduced to 3,2,1 on simulated data

# 4 --> 3

fit_4_3 <- stan(
  file = here::here("stan-models/model4_reduced_3.stan"),
  data = list(
    S = nrow(sim3$simulated),
    A =  ncol(sim3$simulated),
    N = N_camels,
    pos = sim3$simulated,
    age1 = age_lower,
    age2 = age_upper,
    M = sim3$M_initial,
    sigma_r = 0
  ),
  chains = 2,
  iter = 10000,
  verbose = TRUE,
  control = list(max_treedepth = 15)
  ##control = list(adapt_delta = 0.99)
)

diagnos <- ggmcmc(ggs(fit_4_3), here::here("diagnostics/4_3b.pdf"))

# run full model recuded to 2

fit_4_2 <- stan(
  file = here::here("stan-models/model4_reduced_2.stan"),
  data = list(
    S = nrow(sim2$simulated),
    A =  ncol(sim2$simulated),
    N = N_camels,
    pos = sim2$simulated,
    age1 = age_lower,
    age2 = age_upper,
    M = sim2$M_initial,
    sigma_m = 2
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)

diagnos <- ggmcmc(ggs(fit_4_2), here::here("diagnostics/4_2b.pdf"))

# run full model reduced to 1

fit_4_1 <- stan(
  file = here::here("stan-models/model4_reduced_1.stan"),
  data = list(
    S = nrow(sim1$simulated),
    A =  ncol(sim1$simulated),
    N = N_camels,
    pos = sim1$simulated,
    age1 = age_lower,
    age2 = age_upper,
    M = sim1$M_initial,
    sigma_m = 2,
    sigma_r = 0
  ),
  chains = 2,
  iter = 4000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)

diagnos <- ggmcmc(ggs(fit_4_1), here::here("diagnostics/4_1b.pdf"))


# plot the data

# merge age with N_camels, N_pos etc

  # 1 convert to long format
l_pos <- reshape::melt(sim1$simulated)
l_N_camels <- reshape::melt(N_camels)
l_age_upper <- reshape::melt(age_upper)
l_age_lower <- reshape::melt(age_lower)
l_pprev <- reshape::melt(sim1$pprev)
l_pmAbs <- reshape::melt(sim1$pmAbs)
l_ptot <- l_pprev + l_pmAbs

  # 2 stick together by X1 and X2



