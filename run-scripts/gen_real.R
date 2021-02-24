
# reading in the real data

data_raw <- read.csv("data/all_rna_and_serol_new.csv")

# select columns of interest

data_min <- data_raw%>%
  dplyr::select(STUDY, COUNTRY, LOW_AGE, UPP_AGE, RNA_POS,
                RNA_N, SERO_POS, SERO_N, TEST_TYPE,
                SAMPLE)

# filter only those studies which have serology data
# remove studies which did targetted sampling linked to human cases

data_sero <- data_min %>%
  dplyr::filter(!is.na(SERO_POS))%>%
  dplyr::filter(SAMPLE != "targetted_epilink")%>%
  #dplyr::filter(STUDY != "van_doremalen")%>%
  mutate(STUDY_COUNTRY = paste(COUNTRY, STUDY, sep ="_"))#%>%
  #dplyr::filter(STUDY_COUNTRY != "kandeil_Tunisia")

# restructure the data into matrices for stan

n_datasets <- length(unique(data_sero$STUDY_COUNTRY))
n_ages <-  max(dplyr::count(data_sero, STUDY, COUNTRY)$n)

AGE_L <- matrix(data = 21, nrow = n_datasets, ncol = n_ages)
AGE_U <- matrix(data = 22, nrow = n_datasets, ncol = n_ages)
N_CAMELS <- matrix(data = 0, nrow = n_datasets, ncol = n_ages)
SEROPOS <- matrix(data = 0, nrow = n_datasets, ncol = n_ages)
row_names <- unique(data_sero$STUDY_COUNTRY)
rownames(AGE_L) <- row_names
rownames(AGE_U) <- row_names
rownames(N_CAMELS) <- row_names
rownames(SEROPOS) <- row_names
col_names <- 1:n_ages
colnames(AGE_L) <- col_names
colnames(AGE_U) <- col_names
colnames(N_CAMELS) <- col_names
colnames(SEROPOS) <- col_names

# fill these matrices with the jagged arrays 
study_country <- unique(data_sero$STUDY_COUNTRY)

for(s in 1:n_datasets){
  # lower age
c_s_agel <- (data_sero %>% filter(STUDY_COUNTRY == study_country[s]))$LOW_AGE
names(c_s_agel) <- seq(1: length(c_s_agel))

AGE_L[study_country[s], names(c_s_agel)] <- c_s_agel

  # upper age
c_s_ageu <- (data_sero %>% filter(STUDY_COUNTRY == study_country[s]))$UPP_AGE
names(c_s_ageu) <- seq(1: length(c_s_ageu))

AGE_U[study_country[s], names(c_s_ageu)] <- c_s_ageu

  # N_camels
c_s_ncam <- (data_sero %>% filter(STUDY_COUNTRY == study_country[s]))$SERO_N
names(c_s_ncam) <- seq(1: length(c_s_ncam))

N_CAMELS[study_country[s], names(c_s_ncam)] <- c_s_ncam

  # sero_pos
c_s_spos <- (data_sero %>% filter(STUDY_COUNTRY == study_country[s]))$SERO_POS
names(c_s_spos) <- seq(1: length(c_s_spos))

SEROPOS[study_country[s], names(c_s_spos)] <- c_s_spos

}

# M_ZERO <- unique(data_sero$M_ZERO)/100
# 
# M_ZERO_ZERO <- rep(0, length(M_ZERO))

  # add seroprevalnce and ci

data_sero$AGE_MID <- data_sero$LOW_AGE + (data_sero$UPP_AGE - data_sero$LOW_AGE)/2
  
data_sero$seroprevalence <- data_sero$SERO_POS / data_sero$SERO_N

data_sero$ci_low <- apply(data_sero[,c("SERO_POS", "SERO_N")], 1, 
                          function(x) ci_lower(x[1], x[2]))

data_sero$ci_upp <- apply(data_sero[,c("SERO_POS", "SERO_N")], 1, 
                          function(x) ci_upper(x[1], x[2]))

saveRDS(data_sero, file = "data/data_sero.rds")
saveRDS(SEROPOS, file = "data/SEROPOS.rds")
saveRDS(AGE_L, file = "data/AGE_L.rds")
saveRDS(AGE_U, file = "data/AGE_U.rds")
saveRDS(N_CAMELS, file = "data/N_CAMELS.rds")