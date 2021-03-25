
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

# test_sens_spec <- data.frame(TEST_TYPE = TEST_TYPES, SENS = c(sens_ppNT, sens_MN, sens_PM, sens_ELISA, sens_ELISA_MN),
#                              SPEC = c(spec_ppNT, spec_MN, spec_PM, spec_ELISA, spec_ELISA_MN))
# STUDY_TEST_TYPE$id <- 1:nrow(STUDY_TEST_TYPE)
# STUDY_TEST_TYPE <- merge(STUDY_TEST_TYPE, test_sens_spec)
# STUDY_TEST_TYPE <- STUDY_TEST_TYPE[order(STUDY_TEST_TYPE$id),]
# 
# sens <- STUDY_TEST_TYPE$SENS
# spec <- STUDY_TEST_TYPE$SPEC

## reduced data

test_sens_spec <- data.frame(TEST_TYPE = TEST_TYPES, SENS = c(sens_ELISA, sens_PM,sens_ppNT, sens_ELISA_MN),
                             SPEC = c(spec_ELISA, spec_PM, spec_ppNT, spec_ELISA_MN))
STUDY_TEST_TYPE$id <- 1:nrow(STUDY_TEST_TYPE)
STUDY_TEST_TYPE <- merge(STUDY_TEST_TYPE, test_sens_spec)
STUDY_TEST_TYPE <- STUDY_TEST_TYPE[order(STUDY_TEST_TYPE$id),]

spec <- STUDY_TEST_TYPE$SPEC


############
# model 5b #
############
fit_5real_10 <- readRDS("fits/fit_5real_10.rds")

n_datasets <- 10

# fits
p5a_real <- plot_real5(fit = fit_5real_10, data, mabs = 1, sr = 0)



# parameter estimates
posterior <- as.array(fit_5real_10)
fit_5real_10_df <- as.data.frame(fit_5real_10)
labs <- names(fit_5real_10_df)[1:(n_datasets+2)]
color_scheme_set("red")

p5b_real <- mcmc_intervals(posterior, pars = names(fit_5real_10_df)[1:(n_datasets+2)])+ 
  scale_y_discrete(breaks=names(fit_5real_10_df)[1:(n_datasets+2)],
                   labels=labs,
                   limits = rev)


ggsave("figs/fit5b_real.png", plot = p5a_real)
ggsave("figs/param_est5b_real.png", plot = p5b_real)
