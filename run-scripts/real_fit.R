
# reading in the real data

data_raw <- read.csv("data/all_rna_and_serol_new.csv")

# select columns of interest

data_min <- data_raw%>%
  dplyr::select(STUDY, COUNTRY, LOW_AGE, UPP_AGE, RNA_POS,
                RNA_N, SERO_POS, SERO_N, M_ZERO, TEST_TYPE,
                SAMPLE)

# filter only those studies which have serology data
# remove studies which did targetted sampling linked to human cases

data_sero <- data_min %>%
  dplyr::filter(!is.na(SERO_POS))%>%
  dplyr::filter(SAMPLE != "targetted_epilink")

# restructure the data into matrices for stan

n_datasets <- length(unique(data_sero$STUDY))
n_ages <-  max(dplyr::count(data_sero, STUDY, COUNTRY)$n)

