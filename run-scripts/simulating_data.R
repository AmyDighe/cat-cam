source("utils.R")


# define age classes manually 
# (could just use the ones in the real data?!)
n_ages <- 5
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

n_datasets <- 3
N_camels <- matrix(c(1000, 1000, 1000, 1000, 0,
                     1000, 1000, 1000, 1000, 0,
                     1000, 1000, 1000, 1000, 1000), ncol = n_ages, nrow = n_datasets, byrow = T)

sim4 <- sim_data(n_datasets, n_ages,
         gamma = c(0.5, 0.25, 0.05), sigma = 0.2, omega = 2, mabs = 1,
         N_camels, age_upper, age_lower)

sim3 <- sim_data(n_datasets, n_ages,
         gamma = c(0.5, 0.25, 0.05), sigma = 0, omega = 2, mabs = 1,
         N_camels, age_upper, age_lower)

sim2 <- sim_data(n_datasets, n_ages,
         gamma = c(0.5, 0.25, 0.05), sigma = 0.2, omega = 2, mabs = 0,
         N_camels, age_upper, age_lower)

sim1 <- sim_data(n_datasets, n_ages,
                 gamma = c(0.5, 0.25, 0.05), sigma = 0, omega = 2, mabs = 0,
                 N_camels, age_upper, age_lower)

# checks
# 1. if mabs = 1 pmAb and M_intial should be non-zero, if mabs = 0, pmab and M_initial should all be zero

sim1$M_initial
sim1$pmAbs

# 2. these two should be similar for those with >0 in N_camels

sim1$simulated/N_camels
sim1$pmAbs + sim1$pprev
