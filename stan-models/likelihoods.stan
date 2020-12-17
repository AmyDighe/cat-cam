functions{
  
real model1_lpmf(int seropos,
                  int N,
                  real foi, 
                  real age
                  ){
 
real pred_prev;
real loglik;

pred_prev =  1 - exp(-foi*age);

loglik = binomial_lpmf(seropos|N, pred_prev);

return loglik;
  }
}
