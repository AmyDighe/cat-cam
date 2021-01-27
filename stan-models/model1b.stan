#include likelihoods.stan 
data{
  int S; //number of studies
  int A; //number of ages
  int N[S,A]; //number of camels tested per age class per study
  int pos[S,A]; //number of seropositive camels per age class per study
  matrix[S,A] age; //average age per age class per study
}

parameters{
  vector<lower = 0, upper = 10>[S] foi; // force of infection parameter per study
  real <lower = 0, upper = 1> k; // overdispersion
}

model{
  
      k ~ beta(1,3); //prior for overdispersion
      
  for(s in 1:S){
    for(a in 1:A){
      if(!is_inf(age[s,a])){
        target+= model1b_lpmf(pos[s,a]| N[s,a], foi[s], age[s,a], k);
      }
    }
  }
}
