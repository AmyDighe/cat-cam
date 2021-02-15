####################
# plotting results #
####################

## MODEL 2 - EVENTUALLY MAKE IT A FUNCTION 

# get predicted prevalence and CrIs using posterior parameter estimates

fit2 <- as.data.frame(fit2_summary$summary)
foi2 <- fit2[1:n_datasets,] #fit2_summary saved from fit_binomial
foi2$STUDY_COUNTRY <- study_country

data_fit2 <- merge(foi2, data_sero) 

data_fit2$pprev <- apply(X = data_fit2[, c("LOW_AGE", "UPP_AGE", "mean")], MARGIN = 1,
                         FUN = function(x){pprev4_int(age1 = x[1], age2 = x[2], 
                                                      foi = x[3], sigma_r =  fit2["sigma_r",]$mean,
                                                      M = 0, sigma_m = 12)})

data_fit2$pprev_cilow <- apply(X = data_fit2[, c("LOW_AGE", "UPP_AGE", "2.5%")], MARGIN = 1,
                         FUN = function(x){pprev4_int(age1 = x[1], age2 = x[2], 
                                                      foi = x[3], sigma_r =  fit2["sigma_r",]$mean,
                                                      M = 0, sigma_m = 12)})

data_fit2$pprev_ciupp <- apply(X = data_fit2[, c("LOW_AGE", "UPP_AGE", "97.5%")], MARGIN = 1,
                               FUN = function(x){pprev4_int(age1 = x[1], age2 = x[2], 
                                                            foi = x[3], sigma_r =  fit2["sigma_r",]$mean,
                                                            M = 0, sigma_m = 12)})
# plot

data_fit2 <- data_fit2%>%
  mutate(
    item = seq(1:(dim(data_fit2)[1])))
data_fit2 <- gather(data_fit2, "BOUND", "AGE", 14:15)

# get continuous predicted prev
ages <- seq(0.01,20, by = 0.1)
STUDY_COUNTRY <- rep(unique(data_fit2$STUDY_COUNTRY), each = length(ages))
foi_mean <- rep(unique(data_fit2$mean), each = length(ages))
foi_low <- rep(unique(data_fit2$"2.5%"), each = length(ages))
foi_high <- rep(unique(data_fit2$"97.5%"), each = length(ages))
cont_pred <- data.frame(ages, STUDY_COUNTRY, foi_mean, foi_low, foi_high)
cont_pred$pprev_mean <- apply(X = cont_pred[,c("ages", "foi_mean")], MARGIN = 1, 
                              FUN = function(x){ pprev4_int(age1 = x[1] -0.01, 
                                                            age2 = x[1] + 0.01, foi = x[2],
                                                            sigma_r = fit2["sigma_r",]$mean,
                                                            M = 0, sigma_m = 12)})
cont_pred$pprev_low <- apply(X = cont_pred[,c("ages", "foi_low")], MARGIN = 1, 
                         FUN = function(x){ pprev4_int(age1 = x[1] -0.01, 
                                                       age2 = x[1] + 0.01, foi = x[2],
                                                       sigma_r = fit2["sigma_r",]$mean,
                                                       M = 0, sigma_m = 12)})
cont_pred$pprev_high <- apply(X = cont_pred[,c("ages", "foi_high")], MARGIN = 1, 
                          FUN = function(x){ pprev4_int(age1 = x[1] -0.01, 
                                                        age2 = x[1] + 0.01, foi = x[2],
                                                        sigma_r = fit2["sigma_r",]$mean,
                                                        M = 0, sigma_m = 12)})

# PLOT
ggplot(data = data_fit2, aes(y = seroprevalence, x = AGE_MID))+
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



# MODEL 4 

# get predicted prevalence and CrIs using posterior parameter estimates

fit4 <- as.data.frame(fit4_summary$summary)
foi4 <- fit4[1:n_datasets,] #fit4_summary saved from fit_binomial
foi4$STUDY_COUNTRY <- study_country

data_fit4 <- merge(foi4, data_sero) 

data_fit4$pprev <- apply(X = data_fit4[, c("LOW_AGE", "UPP_AGE", "M_ZERO", "mean")], MARGIN = 1,
                         FUN = function(x){total_pprev4(age1 = x[1], age2 = x[2], 
                                                      foi = x[4], sigma_r =  fit4["sigma_r",]$mean,
                                                      M = x[3]/100, sigma_m = fit4["sigma_m",]$mean)})

data_fit4$pprev_cilow <- apply(X = data_fit4[, c("LOW_AGE", "UPP_AGE","M_ZERO", "2.5%")], MARGIN = 1,
                               FUN = function(x){total_pprev4(age1 = x[1], age2 = x[2], 
                                                            foi = x[4], sigma_r =  fit4["sigma_r",]$mean,
                                                            M = x[3]/100, sigma_m = fit4["sigma_m",]$mean)})

data_fit4$pprev_ciupp <- apply(X = data_fit4[, c("LOW_AGE", "UPP_AGE","M_ZERO", "97.5%")], MARGIN = 1,
                               FUN = function(x){total_pprev4(age1 = x[1], age2 = x[2], 
                                                            foi = x[4], sigma_r =  fit4["sigma_r",]$mean,
                                                            M = x[3]/100, sigma_m = fit4["sigma_m",]$mean)})
# plot

data_fit4 <- data_fit4%>%
  mutate(
    item = seq(1:(dim(data_fit4)[1])))
data_fit4 <- gather(data_fit4, "BOUND", "AGE", 14:15)

# get continuous predicted prev
ages <- seq(0.01,20, by = 0.1)
STUDY_COUNTRY <- rep(unique(data_fit4$STUDY_COUNTRY), each = length(ages))
M_ZERO <- rep(unique(data_fit4$M_ZERO), each = length(ages))
foi_mean <- rep(unique(data_fit4$mean), each = length(ages))
foi_low <- rep(unique(data_fit4$"2.5%"), each = length(ages))
foi_high <- rep(unique(data_fit4$"97.5%"), each = length(ages))
cont_pred <- data.frame(ages, STUDY_COUNTRY, M_ZERO, foi_mean, foi_low, foi_high)
cont_pred$pprev_mean <- apply(X = cont_pred[,c("ages", "M_ZERO", "foi_mean")], MARGIN = 1, 
                              FUN = function(x){ total_pprev4(age1 = x[1] -0.01, 
                                                            age2 = x[1] + 0.01, foi = x[3],
                                                            sigma_r = fit4["sigma_r",]$mean,
                                                            M = x[2]/100, sigma_m = fit4["sigma_m",]$mean)})
cont_pred$pprev_low <- apply(X = cont_pred[,c("ages", "M_ZERO", "foi_low")], MARGIN = 1, 
                             FUN = function(x){ total_pprev4(age1 = x[1] -0.01, 
                                                           age2 = x[1] + 0.01, foi = x[3],
                                                           sigma_r = fit4["sigma_r",]$mean,
                                                           M = x[2]/100, sigma_m = fit4["sigma_m",]$mean)})
cont_pred$pprev_high <- apply(X = cont_pred[,c("ages", "M_ZERO", "foi_high")], MARGIN = 1, 
                              FUN = function(x){ total_pprev4(age1 = x[1] -0.01, 
                                                            age2 = x[1] + 0.01, foi = x[3],
                                                            sigma_r = fit4["sigma_r",]$mean,
                                                            M = x[2]/100, sigma_m = fit4["sigma_m",]$mean)})

# PLOT
ggplot(data = data_fit4, aes(y = seroprevalence, x = AGE_MID))+
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
