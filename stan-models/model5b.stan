#include likelihoods.stan 
data{
  int S; //number of studies
  int A; //number of age classes
  int N[S,A]; //number of camels tested per age class per study
  int pos[S,A]; //number of seropositive camels per age class per study
  matrix[S,A] age1; //lower bound per age class per study
  matrix[S,A] age2; //upper bound per age class per study
  real sigma_r; // rate of sero-reversion
  real spec[S]; // specificity of Ab test
  real mabs;
}

parameters{
  vector<lower = 0.00001, upper = 10>[S] foi; // force of infection parameter per study
  real<lower = 0.00001, upper = 10> sigma_m; // rate of waning maternal Abs
  real<lower = 0, upper = 0.999> sens; // test sensitivity
}

model{
  for(s in 1:S){
    for(a in 1:A){
        target+= model4av_b_lpmf(pos[s,a]| N[s,a], foi[s], age1[s,a], age2[s,a], sigma_r, sigma_m, sens, spec[s], mabs);
    }
  }
}

generated quantities{
  matrix <lower = 0, upper = 100> [S,A] seroprevalence;
  matrix[S,A] log_lik;

  for(s in 1:S){
    for(a in 1:A){
      seroprevalence[s,a] = seroprev(foi[s], sigma_r, sigma_m, mabs, age1[s,a], age2[s,a], sens, spec[s]);
      log_lik[s,a] = model4av_b_lpmf(pos[s,a]| N[s,a], foi[s], age1[s,a], age2[s,a], sigma_r, sigma_m, sens, spec[s], mabs);
    }
  }
}
