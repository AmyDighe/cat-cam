# run model 2 with seroreversion - hierarchical

gamma <- c(0.5, 0.25, 0.05)
omega <- 2 # maternal antibodies wane over 6 months
n_datasets <- 3
n_ages <- 5

ages_data <- matrix(c(0.25, 0.5, 2, 4, Inf,
                      0.5, 1, 3, 5, Inf,
                      1, 2, 3, 4, 5), ncol = n_ages, nrow = n_datasets, byrow = T)
N_camels <- matrix(c(1000, 1000, 1000, 1000, 0,
                     1000, 1000, 1000, 1000, 0,
                     1000, 1000, 1000, 1000, 1000), ncol = n_ages, nrow = n_datasets, byrow = T)

# need an estimate for the proportion of calves with mAbs - should be
# a function of number of adults with Abs but we don't know that yet...
# lets estimate it using the simple model as a rough ballpark figure

pprev_base <- function(gamma, age) 1 - exp(-gamma*age)

pos_data_base <- matrix(ncol = n_ages, nrow = n_datasets, byrow = T)
pred_prev_base <- matrix(ncol = n_ages, nrow = n_datasets, byrow = T)

for(s in 1:n_datasets){
  
  for(a in 1:n_ages){
    pred_prev_base[s,a] <- pprev_base(gamma = gamma[s], age = ages_data[s,a])
    pos_data_base[s,a] <- rbinom(n = 1, size = N_camels[s,a], prob = pred_prev_base[s,a])
  }
}


M_initial <- vector(length = n_datasets)

for(i in 1:n_datasets){
  M_initial[i] <- sum(pos_data_base[i,3:5])/sum(N_camels[i,3:5])  
}

# now calculate the predicted prevalence given M_initial and mAbs model

pprev <- function(gamma, age, omega, M_initial){
  1 - exp(-gamma*age) - M_initial*((gamma/(gamma-omega))*(exp(-omega*age) - exp(-gamma*age)))
}
# we want seropositives to include animals with mAbs as well as Abs
# pprev will only estimate those with Abs

pmAbs <- function(M_initial, omega, age) M_initial*exp(-omega*age)

pos_data <- matrix(ncol = n_ages, nrow = n_datasets, byrow = T)
pred_prev <- matrix(ncol = n_ages, nrow = n_datasets, byrow = T)
pred_mAb <- matrix(ncol = n_ages, nrow = n_datasets, byrow = T)


for(s in 1:n_datasets){
  
  for(a in 1:n_ages){
    pred_prev[s,a] <- pprev(gamma = gamma[s], age = ages_data[s,a], omega = omega, M_initial = M_initial[s])
    pred_mAb[s,a]<- pmAbs(M_initial = M_initial[s], omega = omega, age = ages_data[s,a])
    pos_data[s,a] <- rbinom(n = 1, size = N_camels[s,a], prob = pred_prev[s,a] + pred_mAb[s,a])
  }
}

# look at the data

pos_data_df <- as.data.frame(pos_data)



# fit the model to the simulated data

fit_model3 <- stan(
  file = here::here("stan-models/model3_hierarchical.stan"),
  data = list(
    S = nrow(ages_data),
    A =  ncol(ages_data),
    N = N_camels,
    pos = pos_data,
    age = ages_data,
    M = M_initial#,
    #sigma_m = 2
  ),
  chains = 2,
  iter = 4000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)


diagnos <- ggmcmc(ggs(fit_model3), here::here("diagnostics/test_hier_model3.pdf"))
