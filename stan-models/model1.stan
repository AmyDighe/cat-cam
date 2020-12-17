#include likelihoods.stan 
data{
  int A; //number of age classes
  real age[A]; //average age
  int N[A]; //number of camels tested in age class
  int pos[A]; //number of seropositive camels in age class
}

parameters{
  real <lower = 0, upper = 10> foi; // force of infection parameter
}

model{
  for(a in 1:A){
    target+= model1_lpmf(pos[a]| N[a], foi, age[a]);
  }
}
