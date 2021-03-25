
n_datasets <- length(unique(data_sero$STUDY_COUNTRY))
data <- data_sero

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

names(sens) <- unique(data_sero$STUDY_COUNTRY)
names(spec) <- names(sens)

fit_4_1real <- readRDS("fits/fit_4_1real.rds")
fit_4_2real <- readRDS("fits/fit_4_2real.rds")
fit_4_3real <- readRDS("fits/fit_4_3real.rds")
fit_4real <- readRDS("fits/fit_4real.rds")

############
# model 1b #
############

# chains
chains <- ggmcmc(ggs(fit_4_1real), plot = "traceplot", here::here("figs/real_1_trace.pdf"))


# fits
p1a_real <- plot_real_gq(fit = fit_4_1real, data = data, mabs = 0, sr = 0)


# parameter estimates
posterior <- as.array(fit_4_1real)
fit_4_1_df <- as.data.frame(fit_4_1real)
labs <- names(fit_4_1_df)[1:n_datasets]
color_scheme_set("red")

p1b_real <- mcmc_intervals(posterior, pars = names(fit_4_1_df)[1:n_datasets])+ 
  scale_y_discrete(breaks=names(fit_4_1_df)[1:n_datasets],
                   labels=labs,
                   limits = rev)


ggsave("figs/fit1b_real.png", plot = p1a_real)
ggsave("figs/param_est1b_real.png", plot = p1b_real)

############
# model 2b #
############

# fits
p2a_real <- plot_fit_real(fits = fit_4_2real, data, mabs = 0, sr = 1, sens, spec)


# parameter estimates
posterior <- as.array(fit_4_2real)
fit_4_2_df <- as.data.frame(fit_4_2real)
labs <- names(fit_4_2_df)[1:(n_datasets+1)]
color_scheme_set("red")

p2b_real <- mcmc_intervals(posterior, pars = names(fit_4_2_df)[1:(n_datasets+1)])+ 
  scale_y_discrete(breaks=names(fit_4_2_df)[1:(n_datasets+1)],
                   labels=labs,
                   limits = rev)


ggsave("figs/fit2b_real.png", plot = p2a_real)
ggsave("figs/param_est2b_real.png", plot = p2b_real)

############
# model 3b #
############


# fits
p3a_real <- plot_fit_real(fits = fit_4_3real, data, mabs = 1, sr = 0, sens, spec)



# parameter estimates
posterior <- as.array(fit_4_3real)
fit_4_3_df <- as.data.frame(fit_4_3real)
labs <- names(fit_4_3_df)[1:(n_datasets+1)]
color_scheme_set("red")

p3b_real <- mcmc_intervals(posterior, pars = names(fit_4_3_df)[1:(n_datasets+1)])+ 
  scale_y_discrete(breaks=names(fit_4_3_df)[1:(n_datasets+1)],
                   labels=labs,
                   limits = rev)


ggsave("figs/fit3b_real.png", plot = p3a_real)
ggsave("figs/param_est3b_real.png", plot = p3b_real)

############
# model 4b #
############

# fits
p4a_real <- plot_real_gq(fit = fit_4real, data, mabs = 1, sr = 1)

# parameter estimates
posterior <- as.array(fit_4real)
fit_4_df <- as.data.frame(fit_4real)
labs <- names(fit_4_df)[1:(n_datasets+2)]
color_scheme_set("red")

p4b_real <- mcmc_intervals(posterior, pars = names(fit_4_df)[1:(n_datasets+2)])+ 
  scale_y_discrete(breaks=names(fit_4_df)[1:(n_datasets+2)],
                   labels=labs,
                   limits = rev)

ggsave("figs/fit4b_real.png", plot = p4a_real)
ggsave("figs/param_est4b_real.png", plot = p4b_real)



# with overdispersion

fit_4_1bbreal <- readRDS("fits/fit_4_1bbreal.rds")
fit_4_2bbreal <- readRDS("fits/fit_4_2bbreal.rds")
fit_4_3bbreal <- readRDS("fits/fit_4_3bbreal.rds")
fit_4bbreal <- readRDS("fits/fit_4bbreal.rds")

############
# model 1b #
############

# chains
chains <- ggmcmc(ggs(fit_4_1bbreal), plot = "traceplot", here::here("figs/real_1bb_trace.pdf"))

# fits
p1a_real <- plot_real_gq(fit = fit_4_1bbreal, data = data, mabs = 0, sr = 0)


# parameter estimates
posterior <- as.array(fit_4_1bbreal)
fit_4_1_df <- as.data.frame(fit_4_1bbreal)
labs <- names(fit_4_1_df)[1:(n_datasets+1)]
color_scheme_set("red")

p1b_real <- mcmc_intervals(posterior, pars = names(fit_4_1_df)[1:(n_datasets+1)])+ 
  scale_y_discrete(breaks=names(fit_4_1_df)[1:(n_datasets+1)],
                   labels=labs,
                   limits = rev)


ggsave("figs/fit1bb_real.png", plot = p1a_real)
ggsave("figs/param_est1bb_real.png", plot = p1b_real)

############
# model 2b #
############

# chains
chains <- ggmcmc(ggs(fit_4_2bbreal), plot = "traceplot", here::here("figs/real_2bb_trace.pdf"))

# fits
p2a_real <- plot_real_gq(fit = fit_4_2bbreal, data, mabs = 0, sr = 1)


# parameter estimates
posterior <- as.array(fit_4_2real)
fit_4_2_df <- as.data.frame(fit_4_2real)
labs <- names(fit_4_2_df)[1:(n_datasets+1)]
color_scheme_set("red")

p2b_real <- mcmc_intervals(posterior, pars = names(fit_4_2_df)[1:(n_datasets+1)])+ 
  scale_y_discrete(breaks=names(fit_4_2_df)[1:(n_datasets+1)],
                   labels=labs,
                   limits = rev)


ggsave("figs/fit2b_real.png", plot = p2a_real)
ggsave("figs/param_est2b_real.png", plot = p2b_real)

############
# model 3b #
############

# chains
chains <- ggmcmc(ggs(fit_4_3bbreal), plot = "traceplot", here::here("figs/real_3bb_trace.pdf"))

# fits
p3a_real <- plot_fit_real(fits = fit_4_3real, data, mabs = 1, sr = 0, sens, spec)



# parameter estimates
posterior <- as.array(fit_4_3real)
fit_4_3_df <- as.data.frame(fit_4_3real)
labs <- names(fit_4_3_df)[1:(n_datasets+1)]
color_scheme_set("red")

p3b_real <- mcmc_intervals(posterior, pars = names(fit_4_3_df)[1:(n_datasets+1)])+ 
  scale_y_discrete(breaks=names(fit_4_3_df)[1:(n_datasets+1)],
                   labels=labs,
                   limits = rev)


ggsave("figs/fit3b_real.png", plot = p3a_real)
ggsave("figs/param_est3b_real.png", plot = p3b_real)

############
# model 4b #
############

# chains
chains <- ggmcmc(ggs(fit_4bbreal), plot = "traceplot", here::here("figs/real_4bb_trace.pdf"))

# fits
p4a_real <- plot_real_gq(fit = fit_4bbreal, data, mabs = 1, sr = 1)

# parameter estimates
posterior <- as.array(fit_4bbreal)
fit_4_df <- as.data.frame(fit_4bbreal)
labs <- names(fit_4_df)[1:(n_datasets+2)]
color_scheme_set("red")

p4b_real <- mcmc_intervals(posterior, pars = names(fit_4_df)[1:(n_datasets+2)])+ 
  scale_y_discrete(breaks=names(fit_4_df)[1:(n_datasets+2)],
                   labels=labs,
                   limits = rev)

ggsave("figs/fit4b_real.png", plot = p4a_real)
ggsave("figs/param_est4b_real.png", plot = p4b_real)


################################################################################
# only studies with >2 age classes                                             #
################################################################################



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

test_sens_spec <- data.frame(TEST_TYPE = TEST_TYPES, SENS = c(sens_ELISA, sens_PM,sens_ppNT, sens_ELISA_MN),
                             SPEC = c(spec_ELISA, spec_PM, spec_ppNT, spec_ELISA_MN))
STUDY_TEST_TYPE$id <- 1:nrow(STUDY_TEST_TYPE)
STUDY_TEST_TYPE <- merge(STUDY_TEST_TYPE, test_sens_spec)
STUDY_TEST_TYPE <- STUDY_TEST_TYPE[order(STUDY_TEST_TYPE$id),]

sens <- STUDY_TEST_TYPE$SENS
spec <- STUDY_TEST_TYPE$SPEC

fit_4_1real_10 <- readRDS("fits/fit_4_1real_10.rds")
fit_4_2real_10 <- readRDS("fits/fit_4_2real_10.rds")
fit_4_3real_10 <- readRDS("fits/fit_4_3real_10.rds")
fit_4real_10 <- readRDS("fits/fit_4real_10.rds")

############
# model 1b #
############

# chains
chains <- ggmcmc(ggs(fit_4_1real_10), plot = "traceplot", here::here("figs/real_10_1_trace.pdf"))


# fits
p1a_real_10 <- plot_real_gq(fit = fit_4_1real_10, data = data, mabs = 0, sr = 0)


# parameter estimates
posterior <- as.array(fit_4_1real_10)
fit_4_1_df <- as.data.frame(fit_4_1real_10)
labs <- names(fit_4_1_df)[1:n_datasets]
color_scheme_set("red")

p1b_real_10 <- mcmc_intervals(posterior, pars = names(fit_4_1_df)[1:n_datasets])+ 
  scale_y_discrete(breaks=names(fit_4_1_df)[1:n_datasets],
                   labels=labs,
                   limits = rev)


ggsave("figs/fit1b_real_10.png", plot = p1a_real_10)
ggsave("figs/param_est1b_real_10.png", plot = p1b_real_10)

############
# model 2b #
############

# chains
chains <- ggmcmc(ggs(fit_4_2real_10), plot = "traceplot", here::here("figs/real_10_2_trace.pdf"))

# fits
p2a_real_10 <- plot_real_gq(fits = fit_4_2real_10, data, mabs = 0, sr = 1, sens, spec)


# parameter estimates
posterior <- as.array(fit_4_2real_10)
fit_4_2_df <- as.data.frame(fit_4_2real_10)
labs <- names(fit_4_2_df)[1:(n_datasets+1)]
color_scheme_set("red")

p2b_real_10 <- mcmc_intervals(posterior, pars = names(fit_4_2_df)[1:(n_datasets+1)])+ 
  scale_y_discrete(breaks=names(fit_4_2_df)[1:(n_datasets+1)],
                   labels=labs,
                   limits = rev)


ggsave("figs/fit2b_real_10.png", plot = p2a_real_10)
ggsave("figs/param_est2b_real_10.png", plot = p2b_real_10)

############
# model 3b #
############

# chains
chains <- ggmcmc(ggs(fit_4_3real_10), plot = "traceplot", here::here("figs/real_10_3_trace.pdf"))

# fits
p3a_real_10 <- plot_real_gq(fits = fit_4_3real_10, data, mabs = 1, sr = 0, sens, spec)



# parameter estimates
posterior <- as.array(fit_4_3real_10)
fit_4_3_df <- as.data.frame(fit_4_3real_10)
labs <- names(fit_4_3_df)[1:(n_datasets+1)]
color_scheme_set("red")

p3b_real_10 <- mcmc_intervals(posterior, pars = names(fit_4_3_df)[1:(n_datasets+1)])+ 
  scale_y_discrete(breaks=names(fit_4_3_df)[1:(n_datasets+1)],
                   labels=labs,
                   limits = rev)


ggsave("figs/fit3b_real_10.png", plot = p3a_real_10)
ggsave("figs/param_est3b_real_10.png", plot = p3b_real_10)

############
# model 4b #
############

# chains
chains <- ggmcmc(ggs(fit_4real_10), plot = "traceplot", here::here("figs/real_10_4_trace.pdf"))

# fits
p4a_real_10 <- plot_real_gq(fit = fit_4real_10, data, mabs = 1, sr = 1)

# parameter estimates
posterior <- as.array(fit_4real_10)
fit_4_df <- as.data.frame(fit_4real_10)
labs <- names(fit_4_df)[1:(n_datasets+2)]
color_scheme_set("red")

p4b_real_10 <- mcmc_intervals(posterior, pars = names(fit_4_df)[1:(n_datasets+2)])+ 
  scale_y_discrete(breaks=names(fit_4_df)[1:(n_datasets+2)],
                   labels=labs,
                   limits = rev)

ggsave("figs/fit4b_real_10.png", plot = p4a_real_10)
ggsave("figs/param_est4b_real_10.png", plot = p4b_real_10)
