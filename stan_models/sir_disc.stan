
// This Stan program defines a dual model for ILI data and for 
// hospitalizations. The ILI data is modeled according to an SIR
// compartmental model where ILI is treated as the observed 
// infected. The infected portion is the expected value of a beta
// distribution. The mean and scale parameters are independent 
// across seasons.
// A discrepancy reverse random walk is included to the mean to 
// capture any patter tha is shared by all or most flu seasons.
// Flu hospitalizations are then modeled as a linear function of 
// the ILI model, the squared ILI model, and one lagged 
// hospitalization. The parameters in the linear model between 
// the two seasons, 2022 and 2023, are independent.

functions {
  
  
  vector sir(real t, 
             vector y, 
             real beta,
             real gamma,
             int N) {
    
    real S = y[1];
    real I = y[2];
    real R = y[3];
    
    vector[3] dydt;
    
    dydt[1] = -beta * I * S / N;
    dydt[2] =  beta * I * S / N - gamma * I;
    dydt[3] =  gamma * I;
    
    return dydt;
  }
}
data {
  int M;                          //total number of data points
  int<lower=1> n_seasons;
  int<lower=1> n_weeks;
  int<lower=1> cur_yr_n_weeks;
  array[M] int<lower=0> all_seasons;                    // all_seasons
  array[n_seasons] int<lower=1> seasons;
  array[n_seasons] int<lower=1> seg_ind_start;
  array[n_seasons] int<lower=1> seg_ind_length;
  array[n_seasons] int<lower=1> seg_ind_max;
  array[M] int<lower=1> weeks;
  real S0;
  real t0;
  array[M] real<lower=0> ts;                     //array of weeks
  int N;
  array[M] real<lower=0> ili;
  int<lower=0> HM;                                // number of hospitalization points
  int<lower=0> n_seasons_hosp;
  array[HM] real<lower=0> hosp;
  array[HM] int<lower=0> hosp_seasons;
  array[HM] int<lower=0> hosp_season_weeks;
  int<lower=0> ps;
  int<lower=0> hps;
  real<lower=0> count_rate;
  
  //prior hyperparameters
  real<lower=0> sigma_sigma_gamma_w; // .01
  real<lower=0> sigma_sigma_gamma; // .05
  real<lower=0> sigma_kappa;
  real rho_mu; // .68
  real<lower=0> rho_sigma; //.08
  real beta_mu; //.8
  real<lower=0> beta_sigma; //.3
  real I0_mu; //.005
  real I0_sigma; //.03
  real<lower=0> sigma_alpha0; // 5
  real<lower=0> sigma_alpha1; // 5
  real<lower=0> sigma_alpha3; // .4
  real<lower=0> sigma_sigma_epsilon; // 4.
  
  real<lower=0> sigma_mu_alpha_w; // 2
  real<lower=0> sigma_mu_beta_w; // 2
  real<lower=0> sigma_mu_alpha; // 2
  real<lower=0> sigma_mu_beta; // .02
  real<lower=0> sigma_disc_sir; //calculated from discrepancy between MLE observation of final week across all seasons
  real<lower=0> sigma_gamma_W; // 1
  
  
}

parameters {
  array[n_seasons] real<lower=0,upper=S0> rho;
  array[n_seasons] real<lower=0> beta;
  array[n_seasons] real<lower=0,upper=1 - S0> I0;
  array[n_weeks - 1] real mu;
  array[n_seasons] real<lower=0> kappas;
  real<lower=0> sigma_gamma_w;
  real<lower=0> sigma_gamma;
  
  
}
transformed parameters{
  vector[3] y0;
  array[M] vector[3] y;
  array[M] real IS;
  array[n_seasons] real<lower=0> R0;
  array[n_seasons] real<lower=0> gamma;
  array[n_seasons,2] real theta;
  //array[HM] real ili_ps;
  array[1] real mu1; 
  array[n_weeks] real muu; 
  
  mu1[1] = -sum(mu);
  muu = append_array(mu1, mu);
  
  for (i in 1:n_seasons) {
    R0[i] = 1 - S0 - I0[i];
    gamma[i] = rho[i]*beta[i];}
  for (i in 1:n_seasons) {
    {
      
      y0[1] = S0;
      y0[2] = I0[i];
      y0[3] = R0[i];
      
      
      y[seg_ind_start[i]:seg_ind_max[i],] =
        ode_rk45(sir, y0, t0,
                 segment(ts, seg_ind_start[i], seg_ind_length[i]),
                 beta[i], gamma[i],
                 N);
      
      IS[seg_ind_start[i]:seg_ind_max[i]] = 
        y[seg_ind_start[i]:seg_ind_max[i], 2];
      
    }
    
    
  }
  
  
  
  //ili_ps = segment(ili, seg_ind_start[ps - 1], seg_ind_length[ps - 1] +
  //                                             cur_yr_n_weeks);
}

model {
  rho ~ normal(rho_mu, rho_sigma);
  beta ~ normal(beta_mu, beta_sigma);
  I0 ~ normal(I0_mu, I0_sigma);
  
  sigma_gamma_w ~ normal(sigma_disc_sir, sigma_gamma_W);
  sigma_gamma ~ normal(0, sigma_sigma_gamma);

  mu[n_weeks - 1] ~ normal(0, sqrt(sigma_gamma_w));
  for (i in 1:(n_weeks - 2)) {
	mu[i] ~ normal(mu[i + 1], sqrt(sigma_gamma));
  } 
    
   
  for (i in 1:n_seasons) kappas[i] ~ normal(0, sigma_kappa);
 

    
  for (i in 1:M) {
        
    ili[i] ~ beta_proportion(
      (IS[i]+.00001)*exp(muu[weeks[i]])/
        (1 + (exp(muu[weeks[i]]) - 1)*(IS[i]+.00001)), 
      kappas[all_seasons[i]]);
      
  }
  
    
  
  }
generated quantities {
    
  array[cur_yr_n_weeks + 5, n_seasons] real pred_ili;
  array[cur_yr_n_weeks + 5, n_seasons] real pred_ili_asg;
  array[cur_yr_n_weeks + 5] real discrepancy; //as opposed to n_weeks
    
  array[n_weeks] real tp;
  array[n_weeks] vector[3] ypm;
  array[n_weeks] real IP;
  vector[3] yp;
  yp[1] = S0;
  yp[2] = I0[ps];
  yp[3] = 1 - S0 - I0[ps];
  
  for (i in 1:n_weeks) tp[i] = i;
  
  ypm = ode_rk45(sir, yp, t0, tp, beta[ps], gamma[ps], N);
  IP = ypm[,2];
  
  for(j in 1:n_seasons) {
    for (i in 1:(cur_yr_n_weeks + 5)) {
      
      pred_ili[i,j] = beta_proportion_rng(
        (IP[i]+.00001)*exp(muu[i])/
          (1 + (exp(muu[i]) - 1)*(IP[i]+.00001)), 
        kappas[j]);
        
      pred_ili_asg[i,j] = beta_proportion_rng((IP[i]+.00001), 
                                            kappas[j]);
  
  
    }
  }
                     
 }




