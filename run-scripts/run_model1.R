# run model 1 on test data

gamma <- 0.5

ages <- c(0.5,2,4,6,8) #ages in years

pprev <- function(gamma, age) 1 - exp(-gamma*age)

N_ages <- length(ages) 
predprev <- vector(length = N_ages)
N_camels <- rep(100, N_ages)
seropositives <- vector(length = N_ages)

for(i in 1:N_ages){
  
  predprev[i] <- pprev(gamma = 0.5, age = ages[i])
  
  seropositives[i] <- rbinom(n = 1, size = N_camels[i], prob = predprev[i])
}


sim_data <- data.frame(ages, seropositives, N_camels)


fit <- stan(
  file = here::here("stan-models/model1.stan"),
  data = list(
    A = nrow(sim_data),
    age = sim_data$ages,
    N = sim_data$N_camels,
    pos = sim_data$seropositives
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)

diagnos <- ggmcmc(ggs(fit), here::here("diagnostics/test.pdf"))
