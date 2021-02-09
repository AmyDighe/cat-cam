#include likelihoods.stan 
data{
  int S; //number of studies
  int A; //number of age classes
  int N[S,A]; //number of camels tested per age class per study
  int pos[S,A]; //number of seropositive camels per age class per study
  matrix[S,A] age1; //lower bound per age class per study
  matrix[S,A] age2; //upper bound per age class per study
  real M[S]; //initial proportion of baby camels with maternal Abs 
  real sigma_m;
}

parameters{
  vector<lower = 0, upper = 10>[S] foi; // force of infection parameter per study
  real<lower = 0, upper = 10> sigma_r; // rate of waning maternal Abs
}

model{
  for(s in 1:S){
    for(a in 1:A){
      if(!is_inf(age1[s,a])){
        target+= model4av_b_lpmf(pos[s,a]| N[s,a], foi[s], age1[s,a], age2[s,a], sigma_r, sigma_m, M[s]);
      }
    }
  }
}
