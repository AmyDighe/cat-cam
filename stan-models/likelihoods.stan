functions{
  
    
real model1b_lpmf(int seropos,
                  int N,
                  real foi, 
                  real age,
                  real k
                  ){
 
real pred_prev;
real loglik;
real alpha;
real beta;

pred_prev =  1 - exp(-foi*age);
alpha = ((1/k) - 1)* pred_prev;
beta = ((1/k) - 1)* (1- pred_prev);

loglik = beta_binomial_lpmf(seropos|N, alpha, beta);

return loglik;
  }

  
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
    
   real prev4_int(real foi, 
                  real sigma_r, 
                  real sigma_m,
                  real M, 
                  real age1, 
                  real age2){
                    
    real pp;
  
    pp = 1 / (age2 - age1) * (
           foi/(foi + sigma_r) * (
              (age2 - age1) + 1 / (sigma_r + foi) * (
                exp(-(foi + sigma_r)*age2) - exp(-(foi + sigma_r) * age1)
              )
            ) -
        M * foi / (foi + sigma_r - sigma_m) * (
               1 / (foi + sigma_r) * (
                exp(-(sigma_r + foi) * age2) - exp(-(sigma_r + foi) * age1)
              ) -
               1 / sigma_m * (exp(-sigma_m * age2)- exp(-sigma_m * age1))
           )
        );
    
    return(pp);
}

    real seroprev(real foi, 
                  real sigma_r, 
                  real sigma_m,
                  real mabs, 
                  real age1, 
                  real age2,
                  real sens,
                  real spec){
                    
                    real pred_prev;
                    real pred_mab;
                    real pred_prev_tot;
                    real obs_pred_prev;
                    real m_initial_true;
                    real m_initial;
                    
              if(mabs > 0) {
       m_initial_true = prev4_int(foi, sigma_r, sigma_m, 0, 4.0000, 10.0000);
       m_initial = sens * m_initial_true + (1 - spec) * m_initial_true;
    } else {
       m_initial = 0;
    }
          pred_prev = prev4_int(foi, sigma_r, sigma_m, m_initial, age1, age2);
          pred_mab = 1/(age2 - age1)*((-m_initial/sigma_m)*(exp(-sigma_m*age2) - exp(-sigma_m*age1)));
          pred_prev_tot = pred_prev + pred_mab;  
          obs_pred_prev = sens * pred_prev_tot + (1 - spec) * pred_prev_tot;
    
    return(obs_pred_prev);
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
    
    pred_prev = ((foi/(foi+sigma_r))*(1 - exp(-(foi + sigma_r)*age))) - M*((foi/(foi + sigma_r - sigma_m))*(exp(-sigma_m*age) - exp( - (foi + sigma_r)*age)));
    pred_mab = M*exp(-sigma_m*age);
    pred_prev_tot = pred_prev + pred_mab;
    
    loglik = binomial_lpmf(seropos|N, pred_prev_tot);
    
    return loglik;
    }
    
          real model4av_bb_lpmf(int seropos,
                    int N,
                    real foi,
                    real age1,
                    real age2,
                    real sigma_r,
                    real sigma_m,
                    real k,
                    real sens, 
                    real spec,
                    real mabs
                    ){
                      
            real obs_pred_prev;
            real alpha;
            real beta;    
            real loglik;
    
    obs_pred_prev = seroprev(foi, sigma_r, sigma_m, mabs, age1, age2, sens, spec);
    alpha = ((1/k) - 1)* obs_pred_prev;
    beta = ((1/k) - 1)* (1 - obs_pred_prev);
    
    loglik = beta_binomial_lpmf(seropos|N, alpha, beta);
    
    return loglik;
    }
    
        real model4av_b_lpmf(int seropos,
                             int N,
                             real foi,
                             real age1,
                             real age2,
                             real sigma_r,
                             real sigma_m,
                             real sens,
                             real spec,
                             real mabs
                             ){
                      
                real obs_pred_prev;
                real loglik;
                      
    obs_pred_prev = seroprev(foi, sigma_r, sigma_m, mabs, age1, age2, sens, spec);
    
    loglik = binomial_lpmf(seropos|N, obs_pred_prev);
    
    return loglik;
    }
}

          real model4av_bbtrans_lpmf(int seropos,
                    int N,
                    real foi,
                    real age1,
                    real age2,
                    real sigma_r,
                    real sigma_m,
                    real gamma2,
                    real sens, 
                    real spec,
                    real mabs
                    ){
                      
            real obs_pred_prev;
            real k;
            real alpha;
            real beta;    
            real loglik;
    
    obs_pred_prev = seroprev(foi, sigma_r, sigma_m, mabs, age1, age2, sens, spec);
    k = -log(gamma2);
    alpha = ((1/k) - 1)* obs_pred_prev;
    beta = ((1/k) - 1)* (1 - obs_pred_prev);
    
    loglik = beta_binomial_lpmf(seropos|N, alpha, beta);
    
    return loglik;
    }
