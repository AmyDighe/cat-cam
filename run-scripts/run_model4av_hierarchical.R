# run model 4 with prev averages across age class rather than point prev for mean age

gamma <- c(0.5, 0.25, 0.05)
omega <- 2 # maternal antibodies wane over 6 months
sigma <- 0.2
n_datasets <- 3
n_ages <- 5

age_lower <- matrix(c(0, 0.51, 3, 6, 10,
                      0, 1, 3, 4, 10,
                      0.5, 1.5, 2.5, 3.5, 4.5), ncol = n_ages, nrow = n_datasets, byrow = T)
age_upper <- matrix(c(0.5, 2.99, 5.99, 8, 20,
                      0.99, 2.99, 3.99, 5, 20,
                      1.49, 2.49, 3.49, 4.49, 5.49), ncol = n_ages, nrow = n_datasets, byrow = T)
N_camels <- matrix(c(1000, 1000, 1000, 1000, 0,
                     1000, 1000, 1000, 1000, 0,
                     1000, 1000, 1000, 1000, 1000), ncol = n_ages, nrow = n_datasets, byrow = T)

# need an estimate for the proportion of calves with mAbs - should be
# a function of number of adults with Abs but we don't know that yet...
# lets estimate it using the simple model as a rough ballpark figure

pprev_base <- function(gamma, sigma, age1, age2){
  
(1/(age2 - age1))*(gamma/(gamma + sigma))*((age2 - age1) + (1/(sigma + gamma))*(exp(-(gamma + sigma)*age2) - exp(-(gamma + sigma)*age1)))
  
  }

pos_data_base <- matrix(ncol = n_ages, nrow = n_datasets, byrow = T)
pred_prev_base <- matrix(ncol = n_ages, nrow = n_datasets, byrow = T)

for(s in 1:n_datasets){
  
  for(a in 1:n_ages){
    pred_prev_base[s,a] <- pprev_base(gamma = gamma[s], age1 = age_lower[s,a], age2 = age_upper[s,a], sigma = sigma)
    pos_data_base[s,a] <- rbinom(n = 1, size = N_camels[s,a], prob = pred_prev_base[s,a])
  }
}


 M_initial <- vector(length = n_datasets)

for(i in 1:n_datasets){
  M_initial[i] <- sum(pos_data_base[i,3:5])/sum(N_camels[i,3:5])  
}

# now calculate the predicted prevalence given M_initial and mAbs model

pprev <- function(gamma, age2, age1, sigma, omega, M_initial){
  
  (1/(age2 - age1))*
    ((gamma/(gamma + sigma))*((age2 - age1) + 
                               (1/(sigma + gamma))*(exp(-(gamma + sigma)*age2) - 
                                                      exp(-(gamma + sigma)*age1)))-
     (M_initial*gamma/(gamma + sigma - omega))*
       ((1/(gamma + sigma))*(exp(-(sigma + gamma)*age2) - exp(-(sigma + gamma)*age1)) - 
       (1/omega)*(exp(-omega*age2)- exp(-omega*age1))))
}

# we want seropositives to include animals with mAbs as well as Abs
# pprev will only estimate those with Abs

pmAbs <- function(M_initial, omega, age2, age1){
 1/(age2 - age1)*((-M_initial/omega)*(exp(-omega*age2) - exp(-omega*age1))) 
} 

pos_data <- matrix(ncol = n_ages, nrow = n_datasets, byrow = T)
pred_prev <- matrix(ncol = n_ages, nrow = n_datasets, byrow = T)
pred_mAb <- matrix(ncol = n_ages, nrow = n_datasets, byrow = T)


for(s in 1:n_datasets){
  
  for(a in 1:n_ages){
    pred_prev[s,a] <- pprev(gamma = gamma[s], age2 = age_upper[s,a], age1 = age_lower[s,a], sigma = sigma, omega = omega, M_initial = M_initial[s])
    pred_mAb[s,a]<- pmAbs(M_initial = M_initial[s], omega = omega, age2 = age_upper[s,a], age1 = age_lower[s,a])
    pos_data[s,a] <- rbinom(n = 1, size = N_camels[s,a], prob = pred_prev[s,a] + pred_mAb[s,a])
  }
}

# look at the data

pos_data_df <- as.data.frame(pos_data)



# fit the model to the simulated data

fit_model4av <- stan(
  file = here::here("stan-models/model4av_hierarchical.stan"),
  data = list(
    S = nrow(pos_data),
    A =  ncol(pos_data),
    N = N_camels,
    pos = pos_data,
    age1 = age_lower,
    age2 = age_upper,
    M = M_initial
  ),
  chains = 2,
  iter = 4000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)


diagnos <- ggmcmc(ggs(fit_model4av), here::here("diagnostics/test_hier_model4av.pdf"))
