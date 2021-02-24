# define age classes manually 
# (could just use the ones in the real data?!)
n_ages <- 5
n_datasets <- 10
foi <- c(0.5, 0.25, 0.05, 0.75, 1, 0.1, 1.5, 2, 0.2, 0.6)
LOW_AGE <- matrix(c(0, 0.51, 3, 6, 10,
                    0, 1, 3, 4, 10,
                    0.5, 1.5, 2.5, 3.5, 4.5,
                    0, 0.51, 3, 6, 10,
                    0, 1, 3, 4, 10,
                    0, 2, 10, 10, 10,
                    0, 0.51, 3, 6, 10,
                    0, 1, 3, 4, 10,
                    0.5, 1.5, 2.5, 3.5, 4.5,
                    0, 2, 10, 10, 10), 
                  ncol = n_ages, nrow = n_datasets, byrow = T)
UPP_AGE <- matrix(c(0.5, 2.99, 5.99, 8, 20,
                    0.99, 2.99, 3.99, 5, 20,
                    1.49, 2.49, 3.49, 4.49, 5.49,
                    0.5, 2.99, 5.99, 8, 20,
                    0.99, 2.99, 3.99, 5, 20,
                    1.99, 10, 20, 20, 20,
                    0.5, 2.99, 5.99, 8, 20,
                    0.99, 2.99, 3.99, 5, 20,
                    1.49, 2.49, 3.49, 4.49, 5.49,
                    1.99, 10, 20, 20, 20),
                  ncol = n_ages, nrow = n_datasets, byrow = T)

# define number of camels per age class manually
# (could just use the ones in the real data?!)
N_CAMELS <- matrix(c(2000, 2000, 2000, 2000, 0,
                     2000, 2000, 2000, 2000, 0,
                     2000, 2000, 2000, 2000, 2000,
                     2000, 2000, 2000, 2000, 0,
                     2000, 2000, 2000, 2000, 0,
                     2000, 2000, 0, 0, 0,
                     2000, 2000, 2000, 2000, 0,
                     2000, 2000, 2000, 2000, 0,
                     2000, 2000, 2000, 2000, 2000,
                     2000, 2000, 0, 0, 0), ncol = n_ages, nrow = n_datasets, byrow = T)

# define pars necessary to reduce model 4 --> 1, 2, 3
sigma <- c(0, 0.2, 0, 0.2)
omega <- 2.1
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
datak0L <- list()
datak01L <- list()
datak001L <- list()
# simulate 3 sets of datasets using each model 1:4, 
# for 3 different values of overdispersion k
sens <- 0.99
spec <- 1

for(m in 1:4){
  
  simk0[[m]] <- replicate(3, sim_data_binom(n_datasets, n_ages,
                                            gamma = foi, sigma = sigma[m], omega = omega, mabs = mabs[m],
                                            N_CAMELS, UPP_AGE, LOW_AGE, sens, spec), 
                          simplify = FALSE)
  
  simk001[[m]] <- replicate(3, sim_data_betabinom(n_datasets, n_ages,
                                                  gamma = foi, sigma = sigma[m], omega = omega, mabs = mabs[m],
                                                  N_CAMELS, UPP_AGE, LOW_AGE, od = 0.01, sens, spec), 
                            simplify = FALSE)
  
  simk01[[m]]<- replicate(3, sim_data_betabinom(n_datasets, n_ages,
                                                gamma = foi, sigma = sigma[m], omega = omega, mabs = mabs[m],
                                                N_CAMELS, UPP_AGE, LOW_AGE, od = 0.1, sens, spec), 
                          simplify = FALSE)
}

# pull datasets together 
# merge age with N_CAMELS, N_pos etc

# 1 convert to long format
l_N_CAMELS <- reshape2::melt(N_CAMELS, value.name = "SERO_N")
l_UPP_AGE <- reshape2::melt(UPP_AGE, value.name = "UPP_AGE")
l_LOW_AGE <- reshape2::melt(LOW_AGE, value.name = "LOW_AGE")

for(m in 1:4){
  for(k in 1:3){ 
    posk0[[m]][[k]] <- reshape2::melt(simk0[[m]][[k]]$simulated, value.name = paste("pos", k, sep = ""))
    posk001[[m]][[k]] <- reshape2::melt(simk001[[m]][[k]]$simulated, value.name = paste("pos", k, sep = ""))
    posk01[[m]][[k]] <- reshape2::melt(simk01[[m]][[k]]$simulated, value.name = paste("pos", k, sep = ""))
    
    pprevk0[[m]] <- reshape2::melt(simk0[[m]][[1]]$obs_seroprev, value.name = "pprevtot")
    pprevk001[[m]] <- reshape2::melt(simk001[[m]][[1]]$obs_seroprev, value.name = "pprevtot")
    pprevk01[[m]] <- reshape2::melt(simk01[[m]][[1]]$obs_seroprev, value.name = "pprevtot")
    
    pmabk0[[m]] <- reshape2::melt(simk0[[m]][[1]]$pmAbs, value.name = "pmAbs")
    pmabk001[[m]] <- reshape2::melt(simk001[[m]][[1]]$pmAbs, value.name = "pmAbs")
    pmabk01[[m]] <- reshape2::melt(simk01[[m]][[1]]$pmAbs, value.name = "pmAbs")
  }
  
  datak0[[m]] <- Reduce(merge, list(l_N_CAMELS, l_LOW_AGE,
                                    l_UPP_AGE,
                                    posk0[[m]][[1]],
                                    posk0[[m]][[2]],
                                    posk0[[m]][[3]],
                                    pmabk0[[m]],
                                    pprevk0[[m]])
  )
  datak0[[m]] <- datak0[[m]]%>%
    rename(STUDY_COUNTRY = Var1, age_class = Var2)%>%
    mutate(STUDY_COUNTRY = factor(STUDY_COUNTRY, levels = as.character(seq(1, n_datasets, by = 1))),
           AGE_MID = LOW_AGE + ((UPP_AGE - LOW_AGE)/2))
  
  datak0[[m]] <- gather(datak0[[m]], "rep", "pos", pos1, pos2, pos3) 
  datak0[[m]] <- datak0[[m]]%>%
    mutate(
      item = seq(1:(dim(datak0[[m]])[1])))
  
  
  datak0L[[m]] <- gather(datak0[[m]], "bound", "age", LOW_AGE, UPP_AGE)
  
  p0[[m]] <- ggplot(data = datak0L[[m]], aes(x = age, y = pprevtot, group = item, colour = STUDY_COUNTRY))+
    geom_line(size = 2, alpha = 0.5)+
    geom_line(aes(x = AGE_MID, y = pos/SERO_N, group = interaction(STUDY_COUNTRY, rep), colour = STUDY_COUNTRY))+
    geom_point(aes(x = AGE_MID, y = pos/SERO_N, group = interaction(STUDY_COUNTRY, rep), colour = STUDY_COUNTRY))+
    xlim(0,10)+
    ylab("seroprevalence")+
    theme_minimal()+ theme(legend.position="none")   
  
  # for k = 0.001
  
  datak001[[m]] <- Reduce(merge, list(l_N_CAMELS, l_LOW_AGE,
                                      l_UPP_AGE,
                                      posk001[[m]][[1]],
                                      posk001[[m]][[2]],
                                      posk001[[m]][[3]],
                                      pmabk001[[m]],
                                      pprevk001[[m]])
  )
  datak001[[m]] <- datak001[[m]]%>%
    rename(STUDY_COUNTRY = Var1, age_class = Var2)%>%
    mutate(STUDY_COUNTRY = factor(STUDY_COUNTRY, levels = as.character(seq(1, n_datasets, by = 1))),
           AGE_MID = LOW_AGE + ((UPP_AGE - LOW_AGE)/2))
  
  datak001[[m]] <- gather(datak001[[m]], "rep", "pos", pos1, pos2, pos3) 
  datak001[[m]] <- datak001[[m]]%>%
    mutate(
      item = seq(1:(dim(datak001[[m]])[1])))
  
  
  datak001L[[m]] <- gather(datak001[[m]], "bound", "age", LOW_AGE, UPP_AGE)
  
  p001[[m]] <- ggplot(data = datak001L[[m]], aes(x = age, y = pprevtot, group = item, colour = STUDY_COUNTRY))+
    geom_line(size = 2, alpha = 0.5)+
    geom_line(aes(x = AGE_MID, y = pos/SERO_N, group = interaction(STUDY_COUNTRY, rep), colour = STUDY_COUNTRY))+
    geom_point(aes(x = AGE_MID, y = pos/SERO_N, group = interaction(STUDY_COUNTRY, rep), colour = STUDY_COUNTRY))+
    xlim(0,10)+
    ylab("seroprevalence")+
    theme_minimal()+ theme(legend.position="none") 
  
  # for k = 0.01
  
  datak01[[m]] <- Reduce(merge, list(l_N_CAMELS, l_LOW_AGE,
                                     l_UPP_AGE,
                                     posk01[[m]][[1]],
                                     posk01[[m]][[2]],
                                     posk01[[m]][[3]],
                                     pmabk01[[m]],
                                     pprevk01[[m]])
  )
  datak01[[m]] <- datak01[[m]]%>%
    rename(STUDY_COUNTRY = Var1, age_class = Var2)%>%
    mutate(STUDY_COUNTRY = factor(STUDY_COUNTRY, levels = as.character(seq(1, n_datasets, by = 1))),
           AGE_MID = LOW_AGE + ((UPP_AGE - LOW_AGE)/2))
  
  datak01[[m]] <- gather(datak01[[m]], "rep", "pos", pos1, pos2, pos3) 
  datak01[[m]] <- datak01[[m]]%>%
    mutate(
      item = seq(1:(dim(datak01[[m]])[1])))
  
  datak01L[[m]] <- gather(datak01[[m]], "bound", "age", LOW_AGE, UPP_AGE)
  
  p01[[m]] <- ggplot(data = datak01L[[m]], aes(x = age, y = pprevtot, group = item, colour = STUDY_COUNTRY))+
    geom_line(size = 2, alpha = 0.5)+
    geom_line(aes(x = AGE_MID, y = pos/SERO_N, group = interaction(STUDY_COUNTRY, rep), colour = STUDY_COUNTRY))+
    geom_point(aes(x = AGE_MID, y = pos/SERO_N, group = interaction(STUDY_COUNTRY, rep), colour = STUDY_COUNTRY))+
    xlim(0,10)+
    ylab("seroprevalence")+
    theme_minimal()+ theme(legend.position="none") 
}

# save simulated data
saveRDS(datak0, file = "data/sim_datak0")
saveRDS(datak001, file = "data/sim_datak001")
saveRDS(datak01, file = "data/sim_datak01")
saveRDS(simk0, file = "data/sim_matk0")
saveRDS(simk001, file = "data/sim_matk001")
saveRDS(simk01, file = "data/sim_matk01")
# plot in a panel

plegend <- ggplot(data = datak01L[[1]], aes(x = age, y = pprevtot, group = item, colour = STUDY_COUNTRY))+
  geom_line(size = 2, alpha = 0.5)+
  geom_line(aes(x = AGE_MID, y = pos/SERO_N, group = interaction(STUDY_COUNTRY, rep), colour = STUDY_COUNTRY))+
  geom_point(aes(x = AGE_MID, y = pos/SERO_N, group = interaction(STUDY_COUNTRY, rep), colour = STUDY_COUNTRY))

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
