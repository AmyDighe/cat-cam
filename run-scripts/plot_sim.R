
n_datasets <- 10

fit_4_1 <- readRDS("fits/fit_4_1.rds")
fit_4_2 <- readRDS("fits/fit_4_2.rds")
fit_4_3 <- readRDS("fits/fit_4_3.rds")
fit_4 <- readRDS("fits/fit_4.rds")

fit_4_1bb <- readRDS("fits/fit_4_1bb.rds")
fit_4_2bb <- readRDS("fits/fit_4_2bb.rds")
fit_4_3bb <- readRDS("fits/fit_4_3bb.rds")
fit_4bb <- readRDS("fits/fit_4bb.rds")

############
# model 1b #
############

# fits
data <- datak0[[1]]%>% dplyr::filter(rep == "pos1")
data$seroprevalence <- data$pos/data$SERO_N
data <- data %>%
  dplyr::filter(SERO_N > 0)

p1a <- plot_fit_sim(fits = fit_4_1, data, mabs = 0, sr = 0, sens, spec)


# parameter estimates
posterior <- as.array(fit_4_1)
fit_4_1_df <- as.data.frame(fit_4_1)
true <- data.frame(params = names(fit_4_1_df)[1:n_datasets], values = foi)
labs <- names(fit_4_1_df)[1:n_datasets]
color_scheme_set("red")

p1b <- mcmc_intervals(posterior, pars = names(fit_4_1_df)[1:n_datasets])+ 
  scale_y_discrete(breaks=names(fit_4_1_df)[1:n_datasets],
                   labels=labs,
                   limits = rev)+
  geom_point(data = true, aes(x = values, y = params))+
  xlim(0,2.5)


ggsave("figs/fit1b.png", plot = p1a)
ggsave("figs/param_est1b.png", plot = p1b)

############
# model 2b #
############

# fits
data <- datak0[[2]]%>% dplyr::filter(rep == "pos1")
data$seroprevalence <- data$pos/data$SERO_N
data <- data %>%
  dplyr::filter(SERO_N > 0)

p2a <- plot_fit_sim(fits = fit_4_2, data, mabs = 0, sr = 1, sens, spec)


# parameter estimates
posterior <- as.array(fit_4_2)
fit_4_2_df <- as.data.frame(fit_4_2)
true <- data.frame(params = names(fit_4_2)[1:(n_datasets + 1)], value = c(foi, 0.2))
labs <- names(fit_4_2_df)[1:(n_datasets+1)]
color_scheme_set("red")

p2b <- mcmc_intervals(posterior, pars = names(fit_4_2_df)[1:(n_datasets+1)])+ 
  scale_y_discrete(breaks=names(fit_4_2_df)[1:(n_datasets+1)],
                   labels=labs,
                   limits = rev)+
  geom_point(data = true, aes(x = value, y = params))+
  xlim(0,2.5)


ggsave("figs/fit2b.png", plot = p2a)
ggsave("figs/param_est2b.png", plot = p2b)

############
# model 3b #
############


# fits
data <- datak0[[3]]%>% dplyr::filter(rep == "pos1")
data$seroprevalence <- data$pos/data$SERO_N
data <- data %>%
  dplyr::filter(SERO_N > 0)

p3a <- plot_fit_sim(fits = fit_4_3, data, mabs = 1, sr = 0, sens, spec)



# parameter estimates
posterior <- as.array(fit_4_3)
fit_4_3_df <- as.data.frame(fit_4_3)
true <- data.frame(params = names(fit_4_3_df)[1:(n_datasets+1)], value = c(foi, 2.1))
labs <- names(fit_4_3_df)[1:(n_datasets+1)]
color_scheme_set("red")

p3b <- mcmc_intervals(posterior, pars = names(fit_4_3_df)[1:(n_datasets+1)])+ 
  scale_y_discrete(breaks=names(fit_4_3_df)[1:(n_datasets+1)],
                   labels=labs,
                   limits = rev)+
  geom_point(data = true, aes(x = value, y = params))+
  xlim(0,2.5)


ggsave("figs/fit3b.png", plot = p3a)
ggsave("figs/param_est3b.png", plot = p3b)

############
# model 4b #
############

# fits
data <- datak0[[4]]%>% dplyr::filter(rep == "pos1")
data$seroprevalence <- data$pos/data$SERO_N
data <- data %>%
  dplyr::filter(SERO_N > 0)

p4a <- plot_fit_sim(fits = fit_4, data, mabs = 1, sr = 1, sens, spec)

# parameter estimates
posterior <- as.array(fit_4)
fit_4_df <- as.data.frame(fit_4)
true <- data.frame(params = names(fit_4_df)[1:(n_datasets+2)], value = c(foi, 0.2, 2.1))
labs <- names(fit_4_df)[1:(n_datasets+2)]
color_scheme_set("red")

p4b <- mcmc_intervals(posterior, pars = names(fit_4_df)[1:(n_datasets+2)])+ 
  scale_y_discrete(breaks=names(fit_4_df)[1:(n_datasets+2)],
                   labels=labs,
                   limits = rev)+
  geom_point(data = true, aes(x = value, y = params))+
  xlim(0,2.5)

ggsave("figs/fit4b.png", plot = p4a)
ggsave("figs/param_est4b.png", plot = p4b)


##################
## BETABINOMIAL ##
##################

#############
# model 1bb #
#############

# fits
data <- datak001[[1]]%>% dplyr::filter(rep == "pos1")
data$seroprevalence <- data$pos/data$SERO_N
data <- data %>%
  dplyr::filter(SERO_N > 0)

p1ab <- plot_fit_sim(fits = fit_4_1bb, data, mabs = 0, sr = 0, sens, spec)


# parameter estimates
posterior <- as.array(fit_4_1bb)
fit_4_1_df <- as.data.frame(fit_4_1bb)
true <- data.frame(params = names(fit_4_1_df)[1:(n_datasets + 1)], values = c(foi, 0.01))
labs <- names(fit_4_1_df)[1:(n_datasets + 1)]
color_scheme_set("red")

p1bb <- mcmc_intervals(posterior, pars = names(fit_4_1_df)[1:(n_datasets + 1)])+ 
  scale_y_discrete(breaks=names(fit_4_1_df)[1:(n_datasets + 1)],
                   labels=labs,
                   limits = rev)+
  geom_point(data = true, aes(x = values, y = params))+
  xlim(0,2.5)


ggsave("figs/fit1bb.png", plot = p1ab)
ggsave("figs/param_est1bb.png", plot = p1bb)

#############
# model 2bb #
#############

# fits
data <- datak001[[2]]%>% dplyr::filter(rep == "pos1")
data$seroprevalence <- data$pos/data$SERO_N
data <- data %>%
  dplyr::filter(SERO_N > 0)

p2ab <- plot_fit_sim(fits = fit_4_2bb, data, mabs = 0, sr = 1, sens, spec)


# parameter estimates
posterior <- as.array(fit_4_2bb)
fit_4_2_df <- as.data.frame(fit_4_2bb)
true <- data.frame(params = names(fit_4_2bb)[1:(n_datasets + 2)], value = c(foi, 0.2, 0.01))
labs <- names(fit_4_2_df)[1:(n_datasets+2)]
color_scheme_set("red")

p2bb <- mcmc_intervals(posterior, pars = names(fit_4_2_df)[1:(n_datasets+2)])+ 
  scale_y_discrete(breaks=names(fit_4_2_df)[1:(n_datasets+2)],
                   labels=labs,
                   limits = rev)+
  geom_point(data = true, aes(x = value, y = params))+
  xlim(0,2.5)


ggsave("figs/fit2bb.png", plot = p2ab)
ggsave("figs/param_est2bb.png", plot = p2bb)

#############
# model 3bb #
#############


# fits
data <- datak001[[3]]%>% dplyr::filter(rep == "pos1")
data$seroprevalence <- data$pos/data$SERO_N
data <- data %>%
  dplyr::filter(SERO_N > 0)

p3ab <- plot_fit_sim(fits = fit_4_3bb, data, mabs = 1, sr = 0, sens, spec)



# parameter estimates
posterior <- as.array(fit_4_3bb)
fit_4_3_df <- as.data.frame(fit_4_3bb)
true <- data.frame(params = names(fit_4_3_df)[1:(n_datasets+2)], value = c(foi, 2.1, 0.01))
labs <- names(fit_4_3_df)[1:(n_datasets+2)]
color_scheme_set("red")

p3bb <- mcmc_intervals(posterior, pars = names(fit_4_3_df)[1:(n_datasets+2)])+ 
  scale_y_discrete(breaks=names(fit_4_3_df)[1:(n_datasets+2)],
                   labels=labs,
                   limits = rev)+
  geom_point(data = true, aes(x = value, y = params))+
  xlim(0,2.5)


ggsave("figs/fit3bb.png", plot = p3ab)
ggsave("figs/param_est3bb.png", plot = p3bb)

#############
# model 4bb #
#############

# fits
data <- datak001[[4]]%>% dplyr::filter(rep == "pos1")
data$seroprevalence <- data$pos/data$SERO_N
data <- data %>%
  dplyr::filter(SERO_N > 0)

p4ab <- plot_fit_sim(fits = fit_4bb, data, mabs = 1, sr = 1, sens, spec)

# parameter estimates
posterior <- as.array(fit_4bb)
fit_4_df <- as.data.frame(fit_4bb)
true <- data.frame(params = names(fit_4_df)[1:(n_datasets+3)], value = c(foi, 0.2, 2.1, 0.01))
labs <- names(fit_4_df)[1:(n_datasets+3)]
color_scheme_set("red")

p4bb <- mcmc_intervals(posterior, pars = names(fit_4_df)[1:(n_datasets+3)])+ 
  scale_y_discrete(breaks=names(fit_4_df)[1:(n_datasets+3)],
                   labels=labs,
                   limits = rev)+
  geom_point(data = true, aes(x = value, y = params))+
  xlim(0,2.5)

ggsave("figs/fit4bb.png", plot = p4ab)
ggsave("figs/param_est4bb.png", plot = p4bb)
