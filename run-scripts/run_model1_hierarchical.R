# run model on multiple test data-sets - hierarchical

gamma <- c(0.5, 0.25, 0.05)
n_datasets <- 3
n_ages <- 5

ages_data <- matrix(c(0.25, 0.5, 2, 4, 5,
                    0.5, 1, 3, 5, 6,
                    1, 2, 3, 4, 5), ncol = n_ages, nrow = n_datasets, byrow = T)
N_camels <- matrix(c(100, 100, 100, 100, 100,
                   100, 100, 100, 100, 100,
                   100, 100, 100, 100, 100), ncol = n_ages, nrow = n_datasets, byrow = T)

pprev <- function(gamma, age) 1 - exp(-gamma*age)

pos_data <- matrix(ncol = n_ages, nrow = n_datasets, byrow = T)
pred_prev <- matrix(ncol = n_ages, nrow = n_datasets, byrow = T)

for(s in 1:n_datasets){
  
  for(a in 1:n_ages){
    pred_prev[s,a] <- pprev(gamma = gamma[s], age = ages_data[s,a])
    pos_data[s,a] <- rbinom(n = 1, size = N_camels[s,a], prob = pred_prev[s,a])
  }
}

# fit the model to the simulated data

  fit_hier <- stan(
    file = here::here("stan-models/model1_hierarchical.stan"),
    data = list(
      S = nrow(ages_data),
      A =  ncol(ages_data),
      N = N_camels,
      pos = pos_data,
      age = ages_data
    ),
    chains = 2,
    iter = 5000,
    verbose = TRUE
    ##control = list(adapt_delta = 0.99)
  )

  
  diagnos <- ggmcmc(ggs(fit_hier), here::here("diagnostics/test_hier.pdf"))
  