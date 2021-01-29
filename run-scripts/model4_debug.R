# simulate data for full model

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
foi <- c(0.5, 0.25, 0.05)
sigma_r <- 0.2
sim4 <- sim_data_binom(n_datasets, n_ages,
                           gamma = foi, sigma = sigma_r, omega = 2, mabs = 1,
                           N_camels, age_upper, age_lower)

# fit the model to the simulated data

##init_fun <- function() list(sigma_m= 1.4)

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
  chains = 4,
  iter = 4000,
  verbose = TRUE,
  ##init = init_fun,
  seed = 100
  ##control = list(adapt_delta = 0.99)
)

diagnos <- ggmcmc(ggs(fit_model4av), here::here("diagnostics/4_binom.pdf"))
saveRDS(fit_model4av, file = "fit_mod4_notworking.rds")
fits <- as.data.frame(fit_model4av)

fits[which.max(fits$lp__),]

sumfit <- rstan::summary(fit_model4av)

sumfit$c_summary

inits <- get_inits(fit_model4av)

fitsarray <- as.array(fit_model4av)

chain4 <- as.data.frame(fitsarray[,4,])
chain4$iteration <- seq(1:length(chain4$sigma_r))

chain1 <- as.data.frame(fitsarray[,1,])
chain1$iteration <- seq(1:length(chain1$sigma_r))

p <- ggplot2::ggplot(chain1, aes(x = iteration, y = lp__))+
  ggplot2::geom_line()
