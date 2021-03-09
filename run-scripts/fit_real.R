# fit the models to the real data
#init_fun <- function(...) list(foi = c(rep(2.911526, nrow(SEROPOS))))

#################
## SENS & SPEC ##
#################

TEST_TYPES <- unique(STUDY_TEST_TYPE$TEST_TYPE) #read in from global

sens_ppNT <- default_sens
sens_MN <- default_sens
sens_PM <- default_sens
sens_ELISA <- default_sens
sens_ELISA_MN <- sens_ELISA * sens_MN

spec_ppNT <- default_spec 
spec_MN <- default_spec 
spec_PM <- default_spec 
spec_ELISA <- default_spec 
spec_ELISA_MN <- spec_ELISA + spec_MN - (spec_ELISA * spec_MN)

test_sens_spec <- data.frame(TEST_TYPE = TEST_TYPES, SENS = c(sens_ppNT, sens_MN, sens_PM, sens_ELISA, sens_ELISA_MN),
                             SPEC = c(spec_ppNT, spec_MN, spec_PM, spec_ELISA, spec_ELISA_MN))
STUDY_TEST_TYPE$id <- 1:nrow(STUDY_TEST_TYPE)
STUDY_TEST_TYPE <- merge(STUDY_TEST_TYPE, test_sens_spec)
STUDY_TEST_TYPE <- STUDY_TEST_TYPE[order(STUDY_TEST_TYPE$id),]

sens <- STUDY_TEST_TYPE$SENS
spec <- STUDY_TEST_TYPE$SPEC

#############
# BINOMIAL ##
#############

fit_4_1real <- stan(
  file = here::here("stan-models/model4_reduced1b.stan"),
  data = list(
    S = nrow(SEROPOS),
    A =  ncol(SEROPOS),
    N = N_CAMELS,
    pos = SEROPOS,
    age1 = AGE_L,
    age2 = AGE_U,
    sigma_m = 12,
    sigma_r = 0,
    sens = sens,
    spec = spec,
    mabs = -1
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE#,
  ##init = init_fun
  ##control = list(adapt_delta = 0.99)
)

fit1_summary <- rstan::summary(fit_4_1real)
diagnos <- ggmcmc(ggs(fit_4_1real), here::here("diagnostics/real4_1btest.pdf"))


fit_4_2real <- stan(
  file = here::here("stan-models/model4_reduced2b.stan"),
  data = list(
    S = nrow(SEROPOS),
    A =  ncol(SEROPOS),
    N = N_CAMELS,
    pos = SEROPOS,
    age1 = AGE_L,
    age2 = AGE_U,
    sigma_m = 12,
    sens = sens,
    spec = spec,
    mabs = -1
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

 fit3_summary <- rstan::summary(fit_4_3real)
diagnos <- ggmcmc(ggs(fit_4_3real), here::here("diagnostics/4_3breal.pdf"))

 fit_4real <- stan(
  file = here::here("stan-models/model4b.stan"),
  data = list(
    S = nrow(SEROPOS),
    A =  ncol(SEROPOS),
    N = N_CAMELS,
    pos = SEROPOS,
    age1 = AGE_L,
    age2 = AGE_U,
    sens = sens,
    spec = spec,
    mabs = 1
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)

fit4_summary <- rstan::summary(fit_4real)
diagnos <- ggmcmc(ggs(fit_4real), here::here("diagnostics/real4b.pdf"))

saveRDS(fit_4_1real, "fits/fit_4_1real.rds")
saveRDS(fit_4_2real, "fits/fit_4_2real.rds")
saveRDS(fit_4_3real, "fits/fit_4_3real.rds")
saveRDS(fit_4real, "fits/fit_4real.rds")


##################
## BETABINOMIAL ##
##################

fit_4_1bbreal <- stan(
  file = here::here("stan-models/model4_reduced1bb.stan"),
  data = list(
    S = nrow(SEROPOS),
    A =  ncol(SEROPOS),
    N = N_CAMELS,
    pos = SEROPOS,
    age1 = AGE_L,
    age2 = AGE_U,
    sigma_m = 12,
    sigma_r = 0,
    sens = sens,
    spec = spec,
    mabs = -1
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE#,
  ##init = init_fun
  ##control = list(adapt_delta = 0.99)
)

fit1_summary <- rstan::summary(fit_4_1bbreal)
diagnos <- ggmcmc(ggs(fit_4_1bbreal), here::here("diagnostics/real4_1bb.pdf"))


fit_4_2bbreal <- stan(
  file = here::here("stan-models/model4_reduced2bb.stan"),
  data = list(
    S = nrow(SEROPOS),
    A =  ncol(SEROPOS),
    N = N_CAMELS,
    pos = SEROPOS,
    age1 = AGE_L,
    age2 = AGE_U,
    sigma_m = 12,
    sens = sens,
    spec = spec,
    mabs = -1
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99) 
)

fit2_summary <- rstan::summary(fit_4_2bbreal)
diagnos <- ggmcmc(ggs(fit_4_2bbreal), here::here("diagnostics/real4_2bb.pdf"))

fit_4_3bbreal <- stan(
  file = here::here("stan-models/model4_reduced3bb.stan"),
  data = list(
    S = nrow(SEROPOS),
    A =  ncol(SEROPOS),
    N = N_CAMELS,
    pos = SEROPOS,
    age1 = AGE_L,
    age2 = AGE_U,
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

fit3_summary <- rstan::summary(fit_4_3bbreal)
diagnos <- ggmcmc(ggs(fit_4_3bbreal), here::here("diagnostics/real4_3bb.pdf"))

fit_4bbreal <- stan(
  file = here::here("stan-models/model4bb.stan"),
  data = list(
    S = nrow(SEROPOS),
    A =  ncol(SEROPOS),
    N = N_CAMELS,
    pos = SEROPOS,
    age1 = AGE_L,
    age2 = AGE_U,
    sens = sens,
    spec = spec,
    mabs = 1
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)

fit4_summary <- rstan::summary(fit_4bbreal)
diagnos <- ggmcmc(ggs(fit_4bbreal), here::here("diagnostics/real4bb.pdf"))

saveRDS(fit_4_1bbreal, "fits/fit_4_1bbreal.rds")
saveRDS(fit_4_2bbreal, "fits/fit_4_2bbreal.rds")
saveRDS(fit_4_3bbreal, "fits/fit_4_3bbreal.rds")
saveRDS(fit_4bbreal, "fits/fit_4bbreal.rds")

fit_4_3bbreal <- readRDS("fits/fit_4_3bbreal.rds")
