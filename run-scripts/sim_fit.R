# define age classes manually 
# (could just use the ones in the real data?!)
n_ages <- 5
n_datasets <- 10
foi <- c(0.5, 0.25, 0.05, 0.75, 1, 0.1, 1.5)
age_lower <- matrix(c(0, 0.51, 3, 6, 10,
                      0, 1, 3, 4, 10,
                      0.5, 1.5, 2.5, 3.5, 4.5,
                      0, 0.51, 3, 6, 10,
                      0, 1, 3, 4, 10,
                      0, 2,
                      0, 0.51, 3, 6, 10,
                      0, 1, 3, 4, 10,
                      0.5, 1.5, 2.5, 3.5, 4.5,
                      0, 2), 
                    ncol = n_ages, nrow = n_datasets, byrow = T)
age_upper <- matrix(c(0.5, 2.99, 5.99, 8, 20,
                      0.99, 2.99, 3.99, 5, 20,
                      1.49, 2.49, 3.49, 4.49, 5.49,
                      0.5, 2.99, 5.99, 8, 20,
                      0.99, 2.99, 3.99, 5, 20,
                      1.99, 20,
                      0.5, 2.99, 5.99, 8, 20,
                      0.99, 2.99, 3.99, 5, 20,
                      1.49, 2.49, 3.49, 4.49, 5.49,
                      1.99, 20),
                    ncol = n_ages, nrow = n_datasets, byrow = T)

# define number of camels per age class manually
# (could just use the ones in the real data?!)
N_camels <- matrix(c(2000, 2000, 2000, 2000, 0,
                     2000, 2000, 2000, 2000, 0,
                     2000, 2000, 2000, 2000, 2000), ncol = n_ages, nrow = n_datasets, byrow = T)

# define pars necessary to reduce model 4 --> 1, 2, 3
sigma <- c(0, 0.2, 0, 0.2)
mabs <- c(0, 0, 1, 1)

simk0 <- list()
simk001 <- list()
simk01 <- list()
posk0 <- rep(list(list()), 4)
posk001 <- rep(list(list()), 4)
posk01 <- rep(list(list()), 4)
pprevk0 <- rep(list(list()), 4)
pprevk001 <- rep(list(list()), 4)
pprevk01 <- rep(list(list()), 4)
pmabk0 <- rep(list(list()), 4)
pmabk001 <- rep(list(list()), 4)
pmabk01 <- rep(list(list()), 4)
datak0 <- list()
p0 <- list()
datak001 <- list()
p001 <- list()
datak01 <- list()
p01 <- list()
# simulate 3 sets of datasets using each model 1:4, 
# for 3 different values of overdispersion k
sens <- 0.9999
spec <- 1

for(m in 1:4){
  
  simk0[[m]] <- replicate(3, sim_data_binom(n_datasets, n_ages,
                                             gamma = foi, sigma = sigma[m], omega = 2, mabs = mabs[m],
                                             N_camels, age_upper, age_lower, sens, spec), 
                           simplify = FALSE)
  
  simk001[[m]] <- replicate(3, sim_data_betabinom(n_datasets, n_ages,
                                  gamma = foi, sigma = sigma[m], omega = 2, mabs = mabs[m],
                                  N_camels, age_upper, age_lower, od = 0.01, sens, spec), 
            simplify = FALSE)
  
  simk01[[m]]<- replicate(3, sim_data_betabinom(n_datasets, n_ages,
                                                gamma = foi, sigma = sigma[m], omega = 2, mabs = mabs[m],
                                                N_camels, age_upper, age_lower, od = 0.1, sens, spec), 
                          simplify = FALSE)
}

# pull datasets together 
# merge age with N_camels, N_pos etc

# 1 convert to long format
l_N_camels <- reshape2::melt(N_camels, value.name = "N")
l_age_upper <- reshape2::melt(age_upper, value.name = "age_upper")
l_age_lower <- reshape2::melt(age_lower, value.name = "age_lower")

for(m in 1:4){
 for(k in 1:3){ 
  posk0[[m]][[k]] <- reshape2::melt(simk0[[m]][[k]]$simulated, value.name = paste("pos", k, sep = ""))
  posk001[[m]][[k]] <- reshape2::melt(simk001[[m]][[k]]$simulated, value.name = paste("pos", k, sep = ""))
  posk01[[m]][[k]] <- reshape2::melt(simk01[[m]][[k]]$simulated, value.name = paste("pos", k, sep = ""))
  
  pprevk0[[m]] <- reshape2::melt(simk0[[m]][[1]]$pprev, value.name = "pprev")
  pprevk001[[m]] <- reshape2::melt(simk001[[m]][[1]]$pprev, value.name = "pprev")
  pprevk01[[m]] <- reshape2::melt(simk01[[m]][[1]]$pprev, value.name = "pprev")
  
  pmabk0[[m]] <- reshape2::melt(simk0[[m]][[1]]$pmAbs, value.name = "pmAbs")
  pmabk001[[m]] <- reshape2::melt(simk001[[m]][[1]]$pmAbs, value.name = "pmAbs")
  pmabk01[[m]] <- reshape2::melt(simk01[[m]][[1]]$pmAbs, value.name = "pmAbs")
 }
  
  datak0[[m]] <- Reduce(merge, list(l_N_camels, l_age_lower,
                                    l_age_upper,
                                    posk0[[m]][[1]],
                                    posk0[[m]][[2]],
                                    posk0[[m]][[3]],
                                    pmabk0[[m]],
                                    pprevk0[[m]])
                         )
  datak0[[m]] <- datak0[[m]]%>%
    rename(study = Var1, age_class = Var2)%>%
    mutate(pprevtot = pprev + pmAbs,
           study = factor(study, levels = c("1", "2", "3")),
           av_age = age_lower + ((age_upper - age_lower)/2))
  
  datak0[[m]] <- gather(datak0[[m]], "rep", "pos", 6:8) 

  datak0[[m]] <- datak0[[m]]%>%
    mutate(
      item = seq(1:45))
  datak0[[m]] <- gather(datak0[[m]], "bound", "age", 4:5)

  p0[[m]] <- ggplot(data = datak0[[m]], aes(x = age, y = pprevtot, group = item, colour = study))+
  geom_line(size = 2, alpha = 0.5)+
  geom_line(aes(x = av_age, y = pos/N, group = interaction(study, rep), colour = study))+
  geom_point(aes(x = av_age, y = pos/N, group = interaction(study, rep), colour = study))+
  xlim(0,10)+
  ylab("seroprevalence")+
  theme_minimal()+ theme(legend.position="none")   
  
  # for k = 0.001
  
  datak001[[m]] <- Reduce(merge, list(l_N_camels, l_age_lower,
                                    l_age_upper,
                                    posk001[[m]][[1]],
                                    posk001[[m]][[2]],
                                    posk001[[m]][[3]],
                                    pmabk001[[m]],
                                    pprevk001[[m]])
  )
  datak001[[m]] <- datak001[[m]]%>%
    rename(study = Var1, age_class = Var2)%>%
    mutate(pprevtot = pprev + pmAbs,
           study = factor(study, levels = c("1", "2", "3")),
           av_age = age_lower + ((age_upper - age_lower)/2))
  
  datak001[[m]] <- gather(datak001[[m]], "rep", "pos", 6:8) 
  
  datak001[[m]] <- datak001[[m]]%>%
    mutate(
      item = seq(1:45))
  datak001[[m]] <- gather(datak001[[m]], "bound", "age", 4:5)
  
  p001[[m]] <- ggplot(data = datak001[[m]], aes(x = age, y = pprevtot, group = item, colour = study))+
    geom_line(size = 2, alpha = 0.5)+
    geom_line(aes(x = av_age, y = pos/N, group = interaction(study, rep), colour = study))+
    geom_point(aes(x = av_age, y = pos/N, group = interaction(study, rep), colour = study))+
    xlim(0,10)+
    ylab("seroprevalence")+
    theme_minimal()+ theme(legend.position="none") 
  
  # for k = 0.01
  
  datak01[[m]] <- Reduce(merge, list(l_N_camels, l_age_lower,
                                      l_age_upper,
                                      posk01[[m]][[1]],
                                      posk01[[m]][[2]],
                                      posk01[[m]][[3]],
                                      pmabk01[[m]],
                                      pprevk01[[m]])
  )
  datak01[[m]] <- datak01[[m]]%>%
    rename(study = Var1, age_class = Var2)%>%
    mutate(pprevtot = pprev + pmAbs,
           study = factor(study, levels = c("1", "2", "3")),
           av_age = age_lower + ((age_upper - age_lower)/2))
  
  datak01[[m]] <- gather(datak01[[m]], "rep", "pos", 6:8) 
  
  datak01[[m]] <- datak01[[m]]%>%
    mutate(
      item = seq(1:45))
  datak01[[m]] <- gather(datak01[[m]], "bound", "age", 4:5)
  
  p01[[m]] <- ggplot(data = datak01[[m]], aes(x = age, y = pprevtot, group = item, colour = study))+
    geom_line(size = 2, alpha = 0.5)+
    geom_line(aes(x = av_age, y = pos/N, group = interaction(study, rep), colour = study))+
    geom_point(aes(x = av_age, y = pos/N, group = interaction(study, rep), colour = study))+
    xlim(0,10)+
    ylab("seroprevalence")+
    theme_minimal()+ theme(legend.position="none") 
}

# plot in a panel

plegend <- ggplot(data = datak01[[1]], aes(x = age, y = pprevtot, group = item, colour = study))+
  geom_line(size = 2, alpha = 0.5)+
  geom_line(aes(x = av_age, y = pos/N, group = interaction(study, rep), colour = study))+
  geom_point(aes(x = av_age, y = pos/N, group = interaction(study, rep), colour = study))

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(plegend)

panel <- gridExtra::grid.arrange(arrangeGrob(p0[[1]], p001[[1]],p01[[1]],
                                 p0[[2]], p001[[2]],p01[[2]],
                                 p0[[3]], p001[[3]],p01[[3]],
                                 p0[[4]], p001[[4]],p01[[4]], 
                                 ncol = 3), mylegend, nrow=2,heights=c(10, 1))

#########################
## FITTING TO THE DATA ##
#########################


# run full model reduced to 1

fit_4_1 <- stan(
  file = here::here("stan-models/model4_reduced1bb.stan"),
  data = list(
    S = nrow(simk01[[1]][[1]]$simulated),
    A =  ncol(simk01[[1]][[1]]$simulated),
    N = N_camels,
    pos = simk01[[1]][[1]]$simulated,
    age1 = age_lower,
    age2 = age_upper,
    M = simk01[[1]][[1]]$M_initial,
    sigma_m = 2,
    sigma_r = 0
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99)
)

fits1 <- as.data.frame(fit_4_1)
fits1[which.max(fits1$lp__),]

sumfit <- rstan::summary(fit_4_1)
sumfit$c_summary

fits1_long <- reshape2::melt(fits1[,1:4])

true1 <- data.frame(variable = as.factor(c("foi[1]", "foi[2]", "foi[3]", "k")), 
                    true = c(0.5, 0.25, 0.05, 0.1))

ggplot()+
  geom_histogram(data = fits1_long, aes(x = value), binwidth = 0.01)+
  facet_wrap(~variable, scales = "free")+
  geom_vline(data = true1, aes(xintercept = true1$true))
diagnos <- ggmcmc(ggs(fit_4_1), here::here("diagnostics/4_1b.pdf"))

ggplot(data = datak01[[1]]%>%dplyr::filter(rep == "pos1"), 
       aes(x = age, y = pprevtot, group = item, colour = study))+
  geom_line(size = 2, alpha = 0.5)+
  geom_line(aes(x = av_age, y = pos/N, group = interaction(study, rep), colour = study))+
  geom_point(aes(x = av_age, y = pos/N, group = interaction(study, rep), colour = study))+
  xlim(0,10)+
  ylab("seroprevalence")+
  theme_minimal()+ theme(legend.position="none")



fit_4_2 <- stan(
  file = here::here("stan-models/model4_reduced2bb.stan"),
  data = list(
    S = nrow(simk01[[2]][[1]]$simulated),
    A =  ncol(simk01[[2]][[1]]$simulated),
    N = N_camels,
    pos = simk01[[2]][[1]]$simulated,
    age1 = age_lower,
    age2 = age_upper,
    M = simk01[[2]][[1]]$M_initial,
    sigma_m = 2
  ),
  chains = 4,
  iter = 4000,
  verbose = TRUE
  ##control = list(adapt_delta = 0.99) 
)

diagnos <- ggmcmc(ggs(fit_4_2), here::here("diagnostics/4_2b.pdf"))

fits2 <- as.data.frame(fit_4_2)
fits2[which.max(fits2$lp__),]

sumfit <- rstan::summary(fit_4_2)
sumfit$c_summary

fits2_long <- reshape2::melt(fits2[,1:5])

true2 <- data.frame(variable = as.factor(c("foi[1]", "foi[2]", "foi[3]", "sigma_r", "k")), 
                    true = c(0.5, 0.25, 0.05, 0.2, 0.1))

ggplot()+
  geom_histogram(data = fits2_long, aes(x = value), binwidth = 0.05)+
  facet_wrap(~variable, scales = "free")+
  geom_vline(data = true2, aes(xintercept = true2$true))


ggplot(data = datak01[[2]]%>%dplyr::filter(rep == "pos1"), 
       aes(x = age, y = pprevtot, group = item, colour = study))+
  geom_line(size = 2, alpha = 0.5)+
  geom_line(aes(x = av_age, y = pos/N, group = interaction(study, rep), colour = study))+
  geom_point(aes(x = av_age, y = pos/N, group = interaction(study, rep), colour = study))+
  xlim(0,10)+
  ylab("seroprevalence")+
  theme_minimal()+ theme(legend.position="none")


