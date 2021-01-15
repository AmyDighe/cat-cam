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
  
  real model2_lpmf(int seropos,
                    int N,
                    real foi,
                    real age,
                    real sigma_r
                    ){
                      
    real pred_prev;
    real loglik;
    
    pred_prev = (foi/(foi+sigma_r))*(1 - exp(-(foi + sigma_r)*age));
    
    loglik = binomial_lpmf(seropos|N, pred_prev);
    
    return loglik;
    }
    
  real model3_lpmf(int seropos,
                    int N,
                    real foi,
                    real age,
                    real sigma_m,
                    real M){
    real pred_prev;
    real pred_mab;
    real pred_prev_tot;
    real loglik;
    
    pred_prev = 1 - exp(-foi*age) - M*((foi/(foi-sigma_m))*(exp(-sigma_m*age) - exp(-foi*age)));
    pred_mab = M*exp(-sigma_m*age);
    pred_prev_tot = pred_prev + pred_mab;
    print(pred_prev_tot);
    loglik = binomial_lpmf(seropos|N, pred_prev_tot);
    
    return loglik;
    }
    
      real model4_lpmf(int seropos,
                    int N,
                    real foi,
                    real age,
                    real sigma_r,
                    real sigma_m,
                    real M
                    ){
                      
    real pred_prev;
    real pred_mab;
    real pred_prev_tot;
    real loglik;
    
    pred_prev = (foi/(foi+sigma_r))*(1 - exp(-(foi + sigma_r)*age)) - M*((foi/(foi + sigma_r - sigma_m))*(exp(-sigma_m*age) - exp( - (foi + sigma_r)*age)));
    pred_mab = M*exp(-sigma_m*age);
    pred_prev_tot = pred_prev + pred_mab;
    
    loglik = binomial_lpmf(seropos|N, pred_prev_tot);
    
    return loglik;
    }
}
