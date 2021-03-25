#include likelihoods.stan 
data{
  int S; //number of studies
  int A; //number of age classes
  int N[S,A]; //number of camels tested per age class per study
  int pos[S,A]; //number of seropositive camels per age class per study
  matrix[S,A] age1; //lower bound per age class per study
  matrix[S,A] age2; //upper bound per age class per study
  real sigma_m;
  real sigma_r;
  real sens[S];
  real spec[S];
  int mabs;
}

parameters{
  vector <lower = 0.00001, upper = 10>[S] foi; // force of infection parameter per study
}


model{
  for(s in 1:S){
    for(a in 1:A){
        target+= model4av_b_lpmf(pos[s,a]| N[s,a], foi[s], age1[s,a], age2[s,a], sigma_r, sigma_m, sens[s], spec[s], mabs);
    }
  }
}

generated quantities{
  matrix <lower = 0, upper = 100> [S,A] seroprevalence;
  vector[S*A] log_lik;
  int x;
  
  x = 1;

  for(s in 1:S){
    for(a in 1:A){
      seroprevalence[s,a] = seroprev(foi[s], sigma_r, sigma_m, mabs, age1[s,a], age2[s,a], sens[s], spec[s]);
      log_lik[x] = model4av_b_lpmf(pos[s,a]| N[s,a], foi[s], age1[s,a], age2[s,a], sigma_r, sigma_m, sens[s], spec[s], mabs);
    x = x + 1;
    }
  }
}
