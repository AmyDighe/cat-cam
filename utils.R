# UTILS

# taking age class prevalence as the point prevalence at age = midpoint of class

# simple model (NO seroreversion or mAbs)
pprev1 <- function(foi, age) 1 - exp(-foi*age)

# with seroreversion
pprev2 <- function(foi, sigma_r, age){
  (foi/(foi+sigma_r))*
    (1 - exp(-(foi + sigma_r)*age))
}

# with mAbs
pprev3 <- function(foi, sigma_m, M, age){
  1 - exp(-foi*age) - 
    M*((foi/(foi-sigma_m))*
         (exp(-sigma_m*age) - exp(-foi*age)))
}

# full model (with seroreversion AND mAbs)
pprev4 <- function(foi, sigma_r, sigma_m, M, age){
  
  pprev <- (
    (foi/(foi+sigma_r))*
      (1 - exp(-(foi + sigma_r)*age))) - 
    M*((foi/(foi + sigma_r - sigma_m))*
         (exp(-sigma_m*age) - exp( - (foi + sigma_r)*age)))

}

# taking age class prevalence as the average prevalence across the age class

# model 1
pprev1_int <- function(foi, age1, age2){
  
 1/(age2 - age1) * (
   age2 - age1 + (1/foi) * (
     exp(-foi*age2) - exp(-foi*age1))
   )

  }

# model 2

pprev2_int <- function(foi, sigma_r, age1, age2){
  
1 / (age2 - age1) * (
  foi / (foi + sigma_r) * (
      age2 - age1 +  1 / (foi + sigma_r) * (
        exp( - (foi + sigma_r) * age2) - exp( - (foi + sigma_r) * age1)))
     )

}

# model 3

pprev3_int <- function(foi, sigma_m, M, age1, age2){
  
  1/(age2 - age1) * (
  age2 - age1 + 1 / foi * (
    exp(-foi * age2) - exp(-foi * age1)) - (
  M * foi / (foi - sigma_m) * (
  1 / foi * (
    exp(-foi * age2) - exp( -foi * age1)) - 
    1 / sigma_m * (
      exp(-sigma_m*age2) - exp(-sigma_m*age1))
      )
    )
  )
}

# full model 
pprev4_int <- function(foi, sigma_r, sigma_m, M, age1, age2){
  
  1 / (age2 - age1) * (
    foi/(foi + sigma_r) * (
      (age2 - age1) + 1 / (sigma_r + foi) * (
        exp(-(foi + sigma_r)*age2) - exp(-(foi + sigma_r) * age1)
          )
      )-
       M * foi / (foi + sigma_r - sigma_m) * (
         1 / (foi + sigma_r) * (
           exp(-(sigma_r + foi) * age2) - exp(-(sigma_r + foi) * age1)
           ) - 
          1 / sigma_m * (exp(-sigma_m * age2)- exp(-sigma_m * age1))
         )
    )
}

pmAbs <- function(M, sigma_m, age2, age1){
  
  1/(age2 - age1) * (
    - M / sigma_m * (exp(-sigma_m*age2) - exp(-sigma_m*age1))
    )
  
} 

total_pprev4 <- function(foi, sigma_r, sigma_m, M, age1, age2, sens, spec){
  
 total_true <- pprev4_int(foi, sigma_r, sigma_m, M, age1, age2) + pmAbs(M, sigma_m, age2, age1)
 total_obs <- sens*total_true + (1 - spec)*total_true
}

# test reduction

test_reduction <- function(foi, sigma_r, sigma_m, age){

  M_1 <- pprev1(foi, age = 4)
  M_2 <- pprev2(foi, sigma_r, age = 4)
  
  pp4 <- pprev4(foi, sigma_r, sigma_m, M_2, age)
  pp3 <- pprev4(foi, sigma_r = 0, sigma_m, M_1, age)
  pp3_ <- pprev3(foi, sigma_m, M_1, age)
  pp2 <- pprev4(foi, sigma_r, sigma_m, M = 0, age)
  pp2_ <- pprev2(foi, sigma_r, age)
  pp1 <- pprev4(foi, sigma_r = 0, sigma_m, M = 0, age)
  pp1_ <- pprev1(foi, age)

  if(pp3==pp3_) three <- ("4-->3 success") else three <- "4-->3 FAIL"
  if(pp2==pp2_) two <- ("4-->2 success") else two <- "4-->3=2 FAIL"
  if(pp1==pp1_) one <- ("4-->1 success") else one <- "4-->1 FAIL"

    return(c(three, two , one))

}

test_reduction_int <- function(foi, sigma_r, sigma_m, age1, age2){
  
  M_1 <- pprev1_int(foi, age1 = 3.5, age2 = 4.5)
  M_2 <- pprev2_int(foi, sigma_r, age1 = 3.5, age2 = 4.5)
  
  pp4int <- pprev4_int(foi, sigma_r, sigma_m, M_2, age1, age2)
  pp3int <- pprev4_int(foi, sigma_r = 0, sigma_m, M_1, age1, age2)
  pp3int_ <- pprev3_int(foi, sigma_m, M_1, age1, age2)
  pp2int <- pprev4_int(foi, sigma_r, sigma_m, M = 0, age1, age2)
  pp2int_ <- pprev2_int(foi, sigma_r, age1, age2)
  pp1int <- pprev4_int(foi, sigma_r = 0, sigma_m, M = 0, age1, age2)
  pp1int_ <- pprev1_int(foi, age1, age2)
  
  if(pp3int==pp3int_) three <- ("4-->3 success") else three <- "4-->3 FAIL"
  if(pp2int==pp2int_) two <- ("4-->2 success") else two <- "4-->3=2 FAIL"
  if(pp1int==pp1int_) one <- ("4-->1 success") else one <- "4-->1 FAIL"
  
  return(c(three, two , one))
  
}

test_reduction(foi = 0.5, sigma_r = 0.2, sigma_m = 2, age = 2) # all working

test_reduction_int(foi = 0.5, sigma_r = 0.2, sigma_m = 2, age1 = 1.5, age2 = 2.5) # all working 


# simulate data-sets
# assign prevalence of Abs and mAbs

sim_data_betabinom <- function (n_datasets, n_ages, gamma, sigma, omega, mabs,
                      N_camels, age_upper, age_lower, od, sens, spec){
  
  M_initial <- vector(length = n_datasets)
  pred_prev <- matrix(data = NA, nrow = n_datasets, ncol = n_ages, byrow = T)
  pred_mAb <- matrix(data = NA, nrow = n_datasets, ncol = n_ages, byrow = T)
  pos_data <- matrix(data = NA, nrow = n_datasets, ncol = n_ages, byrow = T)
  obs_pprev <- matrix(data = NA, nrow = n_datasets, ncol = n_ages, byrow = T)

for(s in 1:n_datasets){
  
  if(mabs == 1){
    
    M_initial[s] <-total_pprev4(foi = gamma[s],
                               age2 = 10,
                               age1 = 4,
                               sigma_r = sigma,
                               sigma_m = omega,
                               M = 0,
                               sens = sens,
                               spec = spec) 
    
  } else {
    
    M_initial[s] <- 0
  }
  
  for(a in 1:n_ages){
    pred_prev[s,a] <- pprev4_int(foi = gamma[s],
                                 age2 = age_upper[s,a], 
                                 age1 = age_lower[s,a],
                                 sigma_r = sigma, 
                                 sigma_m = omega, 
                                 M = M_initial[s])
    pred_mAb[s,a]<- pmAbs(M = M_initial[s], 
                          sigma_m = omega, 
                          age2 = age_upper[s,a], 
                          age1 = age_lower[s,a])
    obs_pprev[s,a] <- sens * (pred_prev[s,a] + pred_mAb[s,a]) + 
                            (1 - spec) * (pred_prev[s,a] + pred_mAb[s,a])
    pos_data[s,a] <- rbetabinom(n = 1, 
                            size = N_camels[s,a], 
                            prob = obs_pprev[s,a],
                            theta = (1/od)-1)
    }
}
  
  return(list(M_initial = M_initial,
              pmAbs = pred_mAb,
              pAbs = pred_prev,
              obs_seroprev = obs_pprev,
              simulated = pos_data))
}



sim_data_binom <- function (n_datasets, n_ages, gamma, sigma, omega, mabs,
                                N_camels, age_upper, age_lower, sens, spec){
  
  M_initial <- vector(length = n_datasets)
  pred_prev <- matrix(data = NA, nrow = n_datasets, ncol = n_ages, byrow = T)
  pred_mAb <- matrix(data = NA, nrow = n_datasets, ncol = n_ages, byrow = T)
  pos_data <- matrix(data = NA, nrow = n_datasets, ncol = n_ages, byrow = T)
  obs_pprev <- matrix(data = NA, nrow = n_datasets, ncol = n_ages, byrow = T)
  
  for(s in 1:n_datasets){
    
    if(mabs == 1){
      
      M_initial[s] <- total_pprev4(foi = gamma[s],
                                 age2 = 10,
                                 age1 = 4,
                                 sigma_r = sigma,
                                 sigma_m = omega,
                                 M = 0,
                                 sens = sens,
                                 spec = spec)
    } else {
      
      M_initial[s] <- 0
    }
    
    for(a in 1:n_ages){
      pred_prev[s,a] <- pprev4_int(foi = gamma[s],
                                   age2 = age_upper[s,a], 
                                   age1 = age_lower[s,a],
                                   sigma_r = sigma, 
                                   sigma_m = omega, 
                                   M = M_initial[s])
      pred_mAb[s,a]<- pmAbs(M = M_initial[s], 
                            sigma_m = omega, 
                            age2 = age_upper[s,a], 
                            age1 = age_lower[s,a])
      obs_pprev[s,a] <- sens * (pred_prev[s,a] + pred_mAb[s,a]) + 
                              (1 - spec) * (pred_prev[s,a] + pred_mAb[s,a])
      pos_data[s,a] <- rbinom(n = 1, 
                                  size = N_camels[s,a], 
                                  prob = obs_pprev[s,a])
    }
  }
  
  return(list(M_initial = M_initial,
              pmAbs = pred_mAb,
              pAbs = pred_prev,
              obs_seroprev = obs_pprev,
              simulated = pos_data))
}


# get CIs of data points

ci_lower <- function(x, n){
  
  lower <- binom::binom.confint(x = x, n = n, method = "exact")$lower
  return(lower)
}

ci_upper <- function(x, n){
  
  upper <- binom::binom.confint(x = x, n = n, method = "exact")$upper
  return(upper)
}


##############
## PLOTTING ##
##############

plot_fit_sim <- function(fits, data, mabs, sr, sens, spec){
  
  n_datasets <- length(unique(data$STUDY_COUNTRY))
  fit <- as.data.frame((rstan::summary(fits))$summary)
  foi_df <- fit[1:n_datasets, ]
  foi_df$STUDY_COUNTRY <- factor(as.character(1:n_datasets), levels = as.character(1:n_datasets))
  
  summary_fit <- rstan::summary(fits)$summary
  
  if(sr == 1){
    sigma_r <- summary_fit["sigma_r", "mean"]
  } else{
    sigma_r <- 0
  }
  
  data_fit <- merge(foi_df, data) 
  
  if(mabs == 1){
    data_fit$m_zero <- apply(X = data_fit["mean"], MARGIN = 1, 
                             FUN = function(x){ total_pprev4(foi = x[1],
                                                             sigma_r = sigma_r,
                                                             sigma_m = 12, M = 0, 
                                                             age1 = 4, 
                                                             age2 = 10, 
                                                             sens, spec)})
    sigma_m <- fit["sigma_m", "mean"]
    
  } else {
    data_fit$m_zero = 0
    sigma_m <- 2.1
  }
  
  
  data_fit$pprev <- apply(X = data_fit[, c("LOW_AGE", "UPP_AGE", "mean", "m_zero")], MARGIN = 1,
                          FUN = function(x){total_pprev4(age1 = x[1], age2 = x[2], 
                                                         foi = x[3], sigma_r =  sigma_r,
                                                         M = x[4], sigma_m = sigma_m, sens = sens, spec = spec)})
  
  data_fit$pprev_cilow <- apply(X = data_fit[, c("LOW_AGE", "UPP_AGE", "2.5%", "m_zero")], MARGIN = 1,
                                FUN = function(x){total_pprev4(age1 = x[1], age2 = x[2], 
                                                               foi = x[3], sigma_r =  sigma_r,
                                                               M = x[4], sigma_m = sigma_m, sens = sens, spec = spec)})
  
  data_fit$pprev_ciupp <- apply(X = data_fit[, c("LOW_AGE", "UPP_AGE", "97.5%", "m_zero")], MARGIN = 1,
                                FUN = function(x){total_pprev4(age1 = x[1], age2 = x[2], 
                                                               foi = x[3], sigma_r =  sigma_r,
                                                               M = x[4], sigma_m = sigma_m, sens = sens, spec = spec)})
  
  data_fit <- data_fit%>%
    mutate(
      item = seq(1:(dim(data_fit)[1])))
  data_fit <- gather(data_fit, "BOUND", "AGE", c("LOW_AGE", "UPP_AGE"))
  
  # get continuous predicted prev
  ages <- seq(0.01,20, by = 0.1)
  STUDY_COUNTRY <- rep(unique(data_fit$STUDY_COUNTRY), each = length(ages))
  foi_mean <- rep(unique(data_fit$mean), each = length(ages))
  foi_low <- rep(unique(data_fit$"2.5%"), each = length(ages))
  foi_high <- rep(unique(data_fit$"97.5%"), each = length(ages))
  cont_pred <- data.frame(ages, STUDY_COUNTRY, foi_mean, foi_low, foi_high)
  
  if(mabs == 1){
    cont_pred$m_zero <- apply(X = cont_pred["foi_mean"], MARGIN = 1, 
                              FUN = function(x){ total_pprev4(age1 = 4, 
                                                            age2 = 10, foi = x[1],
                                                            sigma_r = sigma_r,
                                                            M = 0, sigma_m = 12, sens, spec)})
  } else{
    cont_pred$m_zero = 0
  }
  
  cont_pred$pprev_mean <- apply(X = cont_pred[,c("ages", "foi_mean", "m_zero")], MARGIN = 1, 
                                FUN = function(x){ total_pprev4(age1 = x[1] -0.01, 
                                                              age2 = x[1] + 0.01, foi = x[2],
                                                              sigma_r = sigma_r,
                                                              M = x[3], sigma_m = sigma_m, sens, spec)})
  cont_pred$pprev_low <- apply(X = cont_pred[,c("ages", "foi_low", "m_zero")], MARGIN = 1, 
                               FUN = function(x){ total_pprev4(age1 = x[1] -0.01, 
                                                             age2 = x[1] + 0.01, foi = x[2],
                                                             sigma_r = sigma_r,
                                                             M = x[3], sigma_m = sigma_m, sens, spec)})
  cont_pred$pprev_high <- apply(X = cont_pred[,c("ages", "foi_high", "m_zero")], MARGIN = 1, 
                                FUN = function(x){ total_pprev4(age1 = x[1] -0.01, 
                                                              age2 = x[1] + 0.01, foi = x[2],
                                                              sigma_r = sigma_r,
                                                              M = x[3], sigma_m = sigma_m,
                                                              sens, spec)})
  
  # PLOT
  p <- ggplot(data = data_fit, aes(y = seroprevalence, x = AGE_MID))+
    #geom_errorbar(aes(ymin = ci_low, ymax = ci_upp))+
    geom_point()+  
    geom_line(size = 2, alpha = 0.2, aes(x = AGE, group = item))+
    geom_errorbar(aes(y = pprev, ymin = pprev_cilow, ymax = pprev_ciupp), col = "firebrick")+
    geom_point(aes(y = pprev), col = "firebrick")+  
    geom_line(data = cont_pred, aes(x = ages, y = pprev_mean), col = "firebrick")+
    geom_ribbon(data = cont_pred, aes(x = ages, y = pprev_mean, ymin = pprev_low, ymax = pprev_high), 
                fill = "firebrick", alpha = 0.1)+
    facet_wrap(~STUDY_COUNTRY)+
    theme_minimal()
  
  print(p)
}


plot_fit_real <- function(fits, data, mabs, sr, sens, spec){
  
  n_datasets <- length(unique(data$STUDY_COUNTRY))
  fit <- as.data.frame((rstan::summary(fits))$summary)
  foi_df <- fit[1:n_datasets, ]
  foi_df$STUDY_COUNTRY <- unique(data$STUDY_COUNTRY)
  
  summary_fit <- rstan::summary(fits)$summary
  
  if(sr == 1){
    sigma_r <- summary_fit["sigma_r", "mean"]
  } else{
    sigma_r <- 0
  }
  
  data_fit <- merge(foi_df, data) 
  
  if(mabs == 1){
    data_fit$m_zero <- apply(X = data_fit["mean"], MARGIN = 1, 
                             FUN = function(x){ total_pprev4(foi = x[1],
                                                             sigma_r = sigma_r,
                                                             sigma_m = 12, M = 0, 
                                                             age1 = 4, 
                                                             age2 = 10, 
                                                             sens, spec)})
    sigma_m <- fit["sigma_m", "mean"]
    
  } else {
    data_fit$m_zero = 0
    sigma_m <- 12
  }
  
  
  data_fit$pprev <- apply(X = data_fit[, c("LOW_AGE", "UPP_AGE", "mean", "m_zero")], MARGIN = 1,
                          FUN = function(x){total_pprev4(age1 = x[1], age2 = x[2], 
                                                         foi = x[3], sigma_r =  sigma_r,
                                                         M = x[4], sigma_m = sigma_m, sens = sens, spec = spec)})
  
  data_fit$pprev_cilow <- apply(X = data_fit[, c("LOW_AGE", "UPP_AGE", "2.5%", "m_zero")], MARGIN = 1,
                                FUN = function(x){total_pprev4(age1 = x[1], age2 = x[2], 
                                                               foi = x[3], sigma_r =  sigma_r,
                                                               M = x[4], sigma_m = sigma_m, sens = sens, spec = spec)})
  
  data_fit$pprev_ciupp <- apply(X = data_fit[, c("LOW_AGE", "UPP_AGE", "97.5%", "m_zero")], MARGIN = 1,
                                FUN = function(x){total_pprev4(age1 = x[1], age2 = x[2], 
                                                               foi = x[3], sigma_r =  sigma_r,
                                                               M = x[4], sigma_m = sigma_m, sens = sens, spec = spec)})
  
  data_fit <- data_fit%>%
    mutate(
      item = seq(1:(dim(data_fit)[1])))
  data_fit <- gather(data_fit, "BOUND", "AGE", c("LOW_AGE", "UPP_AGE"))
  
  # get continuous predicted prev
  ages <- seq(0.01,20, by = 0.1)
  STUDY_COUNTRY <- rep(unique(data_fit$STUDY_COUNTRY), each = length(ages))
  foi_mean <- rep(unique(data_fit$mean), each = length(ages))
  foi_low <- rep(unique(data_fit$"2.5%"), each = length(ages))
  foi_high <- rep(unique(data_fit$"97.5%"), each = length(ages))
  cont_pred <- data.frame(ages, STUDY_COUNTRY, foi_mean, foi_low, foi_high)
  
  if(mabs == 1){
    cont_pred$m_zero <- apply(X = cont_pred["foi_mean"], MARGIN = 1, 
                              FUN = function(x){ total_pprev4(age1 = 4, 
                                                              age2 = 10, foi = x[1],
                                                              sigma_r = sigma_r,
                                                              M = 0, sigma_m = 12, sens, spec)})
  } else{
    cont_pred$m_zero = 0
  }
  
  cont_pred$pprev_mean <- apply(X = cont_pred[,c("ages", "foi_mean", "m_zero")], MARGIN = 1, 
                                FUN = function(x){ total_pprev4(age1 = x[1] -0.01, 
                                                                age2 = x[1] + 0.01, foi = x[2],
                                                                sigma_r = sigma_r,
                                                                M = x[3], sigma_m = sigma_m, sens, spec)})
  cont_pred$pprev_low <- apply(X = cont_pred[,c("ages", "foi_low", "m_zero")], MARGIN = 1, 
                               FUN = function(x){ total_pprev4(age1 = x[1] -0.01, 
                                                               age2 = x[1] + 0.01, foi = x[2],
                                                               sigma_r = sigma_r,
                                                               M = x[3], sigma_m = sigma_m, sens, spec)})
  cont_pred$pprev_high <- apply(X = cont_pred[,c("ages", "foi_high", "m_zero")], MARGIN = 1, 
                                FUN = function(x){ total_pprev4(age1 = x[1] -0.01, 
                                                                age2 = x[1] + 0.01, foi = x[2],
                                                                sigma_r = sigma_r,
                                                                M = x[3], sigma_m = sigma_m,
                                                                sens, spec)})
  
  # PLOT
  p <- ggplot(data = data_fit, aes(y = seroprevalence, x = AGE_MID))+
    geom_errorbar(aes(ymin = ci_low, ymax = ci_upp))+
    geom_point()+  
    geom_line(size = 2, alpha = 0.2, aes(x = AGE, group = item))+
    geom_errorbar(aes(y = pprev, ymin = pprev_cilow, ymax = pprev_ciupp), col = "firebrick")+
    geom_point(aes(y = pprev), col = "firebrick")+  
    geom_line(data = cont_pred, aes(x = ages, y = pprev_mean), col = "firebrick")+
    geom_ribbon(data = cont_pred, aes(x = ages, y = pprev_mean, ymin = pprev_low, ymax = pprev_high), 
                fill = "firebrick", alpha = 0.1)+
    facet_wrap(~STUDY_COUNTRY)+
    theme_minimal()
  
  print(p)
}