
n_datasets <- length(unique(data_sero$STUDY_COUNTRY))
data <- data_sero

############
# model 1b #
############

# fits
p1a_real <- plot_fit_real(fits = fit_4_1real, data, mabs = 0, sr = 0, sens, spec)


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
p4a_real <- plot_fit_real(fits = fit_4real, data, mabs = 1, sr = 1, sens, spec)

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
