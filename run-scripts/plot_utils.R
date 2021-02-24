

data_sero <- readRDS("data/data_sero.rds")
#data <- data_sero

data <- datak0[[1]]%>% dplyr::filter(rep == "pos1")
data$seroprevalence <- data$pos/data$SERO_N

data <- data %>%
  dplyr::filter(SERO_N > 0)

fits <- fit_4_1
mabs <- 0
sr <- 0
sens <- 0.999
spec <- 1

plot_fit_sim <- function(fits, data, mabs, sr, sens, spec){
  
  n_datasets <- length(unique(data$STUDY_COUNTRY))
  fit <- as.data.frame((rstan::summary(fits))$summary)
  foi_df <- fit[1:n_datasets, ]
  foi_df$STUDY_COUNTRY <- factor(as.character(1:n_datasets), levels = as.character(1:n_datasets))

  
if(sr == 1){
  sigma_r <- fits["sigma_r",]$mean
} else{
  sigma_r <- 0
}

data_fit <- merge(foi_df, data) 

if(mabs == 1){
  data_fit$m_zero <- apply(X = data_fit["foi_mean"], MARGIN = 1, 
                            FUN = function(x){ total_pprev4(foi = x[1],
                                                          sigma_r = sigma_r,
                                                          sigma_m = 12, M = 0, 
                                                          age1 = 4, 
                                                          age2 = 10, 
                                                          sens, spec)})
  sigma_m <- fit["sigma_m"]$mean
  
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
                           FUN = function(x){ pprev4_int(age1 = 4, 
                                                         age2 = 10, foi = x[1],
                                                         sigma_r = sigma_r,
                                                         M = 0, sigma_m = 12)})
} else{
  cont_pred$m_zero = 0
}

cont_pred$pprev_mean <- apply(X = cont_pred[,c("ages", "foi_mean", "m_zero")], MARGIN = 1, 
                              FUN = function(x){ pprev4_int(age1 = x[1] -0.01, 
                                                            age2 = x[1] + 0.01, foi = x[2],
                                                            sigma_r = sigma_r,
                                                            M = x[3], sigma_m = sigma_m)})
cont_pred$pprev_low <- apply(X = cont_pred[,c("ages", "foi_low", "m_zero")], MARGIN = 1, 
                             FUN = function(x){ pprev4_int(age1 = x[1] -0.01, 
                                                           age2 = x[1] + 0.01, foi = x[2],
                                                           sigma_r = sigma_r,
                                                           M = x[3], sigma_m = sigma_m)})
cont_pred$pprev_high <- apply(X = cont_pred[,c("ages", "foi_high", "m_zero")], MARGIN = 1, 
                              FUN = function(x){ pprev4_int(age1 = x[1] -0.01, 
                                                            age2 = x[1] + 0.01, foi = x[2],
                                                            sigma_r = sigma_r,
                                                            M = x[3], sigma_m = sigma_m)})

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

plot_fit(fits, data, mabs, sr, sens, spec)
  