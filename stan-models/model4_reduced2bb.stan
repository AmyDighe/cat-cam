#include likelihoods.stan 
data{
  int S; //number of studies
  int A; //number of age classes
  int N[S,A]; //number of camels tested per age class per study
  int pos[S,A]; //number of seropositive camels per age class per study
  matrix[S,A] age1; //lower bound per age class per study
  matrix[S,A] age2; //upper bound per age class per study
  real sigma_m;
  real sens[S];
  real spec[S];
  real mabs;
}

parameters{
  vector<lower = 0, upper = 10>[S] foi; // force of infection parameter per study
  real<lower = 0, upper = 1> sigma_r; // rate of waning maternal Abs
  real <lower = 0.0001, upper = 0.999> k; // overdispersion
}

model{
          k ~ beta(1,100); //prior for overdispersion
  for(s in 1:S){
    for(a in 1:A){
      if(!is_inf(age1[s,a])){
        target+= model4av_bb_lpmf(pos[s,a]| N[s,a], foi[s], age1[s,a], age2[s,a], sigma_r, sigma_m, k, sens[s], spec[s], mabs);
      }
    }
  }
}

generated quantities{
  matrix <lower = 0, upper = 100> [S,A] seroprevalence;
  matrix[S,A] log_lik;

  for(s in 1:S){
    for(a in 1:A){
      seroprevalence[s,a] = seroprev(foi[s], sigma_r, sigma_m, mabs, age1[s,a], age2[s,a], sens[s], spec[s]);
      log_lik[s,a] = model4av_bb_lpmf(pos[s,a]| N[s,a], foi[s], age1[s,a], age2[s,a], sigma_r, sigma_m, k, sens[s], spec[s], mabs);
    }
  }
}

