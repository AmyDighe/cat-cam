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

# full data

test_sens_spec <- data.frame(TEST_TYPE = TEST_TYPES, SENS = c(sens_ppNT, sens_MN, sens_PM, sens_ELISA, sens_ELISA_MN),
                             SPEC = c(spec_ppNT, spec_MN, spec_PM, spec_ELISA, spec_ELISA_MN))
STUDY_TEST_TYPE$id <- 1:nrow(STUDY_TEST_TYPE)
STUDY_TEST_TYPE <- merge(STUDY_TEST_TYPE, test_sens_spec)
STUDY_TEST_TYPE <- STUDY_TEST_TYPE[order(STUDY_TEST_TYPE$id),]

sens <- STUDY_TEST_TYPE$SENS
spec <- STUDY_TEST_TYPE$SPEC

## reduced data

# test_sens_spec <- data.frame(TEST_TYPE = TEST_TYPES, SENS = c(sens_ELISA, sens_PM,sens_ppNT, sens_ELISA_MN),
#                              SPEC = c(spec_ELISA, spec_PM, spec_ppNT, spec_ELISA_MN))
# STUDY_TEST_TYPE$id <- 1:nrow(STUDY_TEST_TYPE)
# STUDY_TEST_TYPE <- merge(STUDY_TEST_TYPE, test_sens_spec)
# STUDY_TEST_TYPE <- STUDY_TEST_TYPE[order(STUDY_TEST_TYPE$id),]
# 
# sens <- STUDY_TEST_TYPE$SENS
# spec <- STUDY_TEST_TYPE$SPEC

#n_datasets <- 25
#############
# BINOMIAL ##
#############

fit_4_1real_10 <- stan(
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
  iter = 8000,
  warmup = 2000,
  verbose = TRUE#,
  ##init = init_fun
  ##control = list(adapt_delta = 0.99)
)

fit1_summary <- rstan::summary(fit_4_1real)
diagnos <- ggmcmc(ggs(fit_4_1real_10), plot = "traceplot", here::here("diagnostics/real4_1btrace10.pdf"))

saveRDS(fit_4_1real_10, "fits/fit_4_1real_10.rds")


fit_4_2real_10 <- stan(
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
  iter = 8000,
  warmup = 2000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99) 
)
saveRDS(fit_4_2real_10, "fits/fit_4_2real_10.rds")
diagnos <- ggmcmc(ggs(fit_4_2real_10), plot = "traceplot", here::here("diagnostics/real4_2btrace10.pdf"))
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
  iter = 10000,
  warmup = 2000,
  verbose = TRUE,
  ##control = list(max_treedepth = 15)
  ##control = list(adapt_delta = 0.99)
)
saveRDS(fit_4_3real, "fits/fit_4_3real.rds")
diagnos <- ggmcmc(ggs(fit_4_3real_10), plot = "traceplot", here::here("diagnostics/real4_3btrace10.pdf"))
 fit3_summary <- rstan::summary(fit_4_3real)
diagnos <- ggmcmc(ggs(fit_4_3real), here::here("diagnostics/4_3breal.pdf"))

 fit_4real_10 <- stan(
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
  iter = 8000,
  warmup = 2000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)
 saveRDS(fit_4real_9, "fits/fit_4real_9.rds")
 diagnos <- ggmcmc(ggs(fit_4real_9), plot = "traceplot", here::here("diagnostics/real4btrace9.pdf"))
 
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
  iter = 10000,
  warmup = 2000,
  verbose = TRUE#,
  ##init = init_fun
  ##control = list(adapt_delta = 0.99)
)

fit1_summary <- rstan::summary(fit_4_1bbreal)
diagnos <- ggs_traceplot(ggs(fit_4_1bbreal), family = "foi")

fit1 <- as.data.frame(fit_4_1bbreal)
fit1 <- fit1[, 1:(n_datasets + 1)]
diagnos <- ggmcmc(ggs(fit_4_1bbreal), plot = "traceplot", here::here("diagnostics/real4_1bbtrace.pdf"))

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
  iter = 10000,
  warmup = 2000,
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
  iter = 10000,
  warmup = 2000,
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
  iter = 10000,
  warmup = 2000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)

fit4_summary <- rstan::summary(fit_4bbreal)
diagnos <- ggmcmc(ggs_traceplot(fit_4bbreal, family = "foi"), here::here("diagnostics/real4bb.pdf"))
ggs_traceplot(fit4, family = "foi")

saveRDS(fit_4_1bbreal, "fits/fit_4_1bbreal.rds")
saveRDS(fit_4_2bbreal, "fits/fit_4_2bbreal.rds")
saveRDS(fit_4_3bbreal, "fits/fit_4_3bbreal.rds")
saveRDS(fit_4bbreal, "fits/fit_4bbreal.rds")

fit_4_1bbreal <- readRDS("fits/fit_4_1bbreal.rds")
fit_4_2bbreal <- readRDS("fits/fit_4_2bbreal.rds")
fit_4_3bbreal <- readRDS("fits/fit_4_3bbreal.rds")
fit_4bbreal <- readRDS("fits/fit_4bbreal.rds")

## FITTING SENSITIVITY - MODEL 5

fit_5real_9 <- stan(
  file = here::here("stan-models/model5b.stan"),
  data = list(
    S = nrow(SEROPOS),
    A =  ncol(SEROPOS),
    N = N_CAMELS,
    pos = SEROPOS,
    age1 = AGE_L,
    age2 = AGE_U,
    sigma_r = 0,
    spec = spec,
    mabs = 1
  ),
  chains = 4,
  iter = 8000,
  warmup = 2000,
  verbose = TRUE,
  ##control = list(max_treedepth = 15)
  ##control = list(adapt_delta = 0.99)
)

saveRDS(fit_5real_9, "fits/fit_5real_9.rds")
diagnos <- ggmcmc(ggs(fit_5real_9), plot = "traceplot", here::here("diagnostics/real5btrace9.pdf"))

# study specific sensitivity


fit_5real_10_studspec <- stan(
  file = here::here("stan-models/model5b_studspec.stan"),
  data = list(
    S = nrow(SEROPOS),
    A =  ncol(SEROPOS),
    N = N_CAMELS,
    pos = SEROPOS,
    age1 = AGE_L,
    age2 = AGE_U,
    sigma_r = 0,
    spec = spec,
    mabs = 1
  ),
  chains = 4,
  iter = 8000,
  warmup = 2000,
  verbose = TRUE,
  ##control = list(max_treedepth = 15)
  ##control = list(adapt_delta = 0.99)
)

saveRDS(fit_5real_10_studspec, "fits/fit_5real_10_studspec.rds")
diagnos <- ggmcmc(ggs(fit_5real_10_studspec), plot = "traceplot", here::here("diagnostics/real5btrace10_studspec.pdf"))

