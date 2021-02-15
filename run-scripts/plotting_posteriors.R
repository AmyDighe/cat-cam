# plotting posteriors

fit2_df <- as.data.frame(fit_4_2real)
fit4_df <- as.data.frame(fit_4real)

posterior <- as.array(fit_4real)
dim(posterior)


labs <- c(unique(STUDY_COUNTRY), "sigma_r", "sigma_m")
color_scheme_set("red")
mcmc_intervals(posterior, pars = names(fit4_df)[1:27])+ 
  scale_y_discrete(breaks=names(fit4_df)[1:27],
                   labels=labs,
                   limits = rev)


# for model 2

posterior <- as.array(fit_4_2real)
dim(posterior)


labs <- c(unique(STUDY_COUNTRY), "sigma_r")
color_scheme_set("red")
mcmc_intervals(posterior, pars = names(fit2_df)[1:26])+ 
  scale_y_discrete(breaks=names(fit2_df)[1:26],
                   labels=labs,
                   limits = rev)
