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

sim_data <- function (n_datasets, n_ages, gamma, sigma, omega, mabs,
                      N_camels, age_upper, age_lower, overdisp){
  
  M_initial <- vector(length = n_datasets)

for(s in 1:n_datasets){
  
  if(mabs == 1){
    
    M_initial[s] <- pprev4_int(foi = gamma[s],
                               age2 = 4.5,
                               age1 = 3.5,
                               sigma_r = sigma,
                               sigma_m = omega,
                               M = 0)
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
    pos_data[s,a] <- rbetabinom(n = 1, 
                            size = N_camels[s,a], 
                            prob = pred_prev[s,a] + pred_mAb[s,a],
                            theta = (1/od)-1)
    }
}
  
  return(list(M_initial = M_initial,
              pmAbs = pred_mAb,
              pprev = pred_prev,
              simulated = pos_data))
}


