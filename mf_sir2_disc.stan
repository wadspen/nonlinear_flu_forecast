//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
functions {
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {

      real S = y[1];
      real I = y[2];
      real R = y[3];
      real N = x_i[1];
      
      real beta = theta[1];
      real gamma = theta[2];
      
      real<lower=0> dS_dt = -beta * I * S / N;
      real<lower=0> dI_dt =  beta * I * S / N - gamma * I;
      real<lower=0> dR_dt =  gamma * I;
      
      return {dS_dt, dI_dt, dR_dt};
  }
  

}
data {
  int M;                          //total number of data points
  int<lower=1> n_seasons;
  int<lower=1> n_weeks;
  int<lower=1> cur_yr_n_weeks;
  int<lower=1> seasons[n_seasons];
  int<lower=1> seg_ind_start[n_seasons];
  int<lower=1> seg_ind_length[n_seasons];
  int<lower=1> seg_ind_max[n_seasons];
  int<lower=1> weeks[M];
  real S0;
  real t0;
  real ts[M];                     //array of weeks
  int N;
  real<lower=0> ili[M];
  int<lower=0> hosp[cur_yr_n_weeks];
  int ps;
  //int pos;
}
transformed data {
  real x_r[0];
  int x_i[1] = { N };
}
parameters {
  real<lower=0,upper=S0> rho[n_seasons];
  real<lower=0> beta[n_seasons];
  real<lower=0.0001,upper=.1> I0[n_seasons];
  real<lower=0> phi_inv;
  real mu[n_weeks];
  real<lower=0> sigmuT2;
  real<lower=0> sigmu2;
  
  //real alpha[n_weeks];
  //real<lower=0> sigalphaT2;
  //real<lower=0> sigalpha2;
  
  real alpha0;
  real alpha1;
  real<lower=0> xi;
  
}
transformed parameters{
  real y0[3];
  real y[M, 3];
  vector[M] IS;
  real phi = 1. / phi_inv;
  real<lower=0> R0[n_seasons];
  real<lower=0> gamma[n_seasons];
  real theta[n_seasons,2];
  real ili_ps[cur_yr_n_weeks];
  
  for (i in 1:n_seasons) {
    R0[i] = 1 - S0 - I0[i];
    gamma[i] = rho[i]*beta[i];
    //pos = 1;
    {
      
      theta[i,1] = beta[i];
      theta[i,2] = gamma[i];
      y0 = {S0, I0[i], R0[i]};
     // print(y0);
     // print(segment(ts, seg_ind_start[i], seg_ind_length[i]));
      
        //print(dims(y[seg_ind_start[i]:seg_ind_max[i],]))
        //print("seg start ", seg_ind_start[i]);
        //print("seg length ", seg_ind_length[i]);
        //print(segment(ts, seg_ind_start[i], seg_ind_length[i]));
        y[seg_ind_start[i]:seg_ind_max[i],] =
          integrate_ode_rk45(sir, y0, t0, 
                              segment(ts, seg_ind_start[i], seg_ind_length[i]),
                              theta[i,], x_r, x_i);
            //pos = pos + seg_ind_length[i];
            
        IS[seg_ind_start[i]:seg_ind_max[i]] = col(
                            to_matrix(segment(y, 
                            seg_ind_start[i], seg_ind_length[i])), 2);
                            
    }
    // print(IS);
    
  }
  // for (i in 1:M) {
  //     if (IS[i] <= 0) IS[i] = .000001;
  //   }
  
  // for (i in 1:M) {
  //   print(y[i,]);
  // }
    
  ili_ps = segment(ili, seg_ind_start[ps], seg_ind_length[ps]);
}

model {
  //priors
  rho ~ normal(0.68, 0.08);
  beta ~ normal(.8, .3);
  I0 ~ normal(.005, .03);
  phi_inv ~ exponential(5);
  sigmuT2 ~ gamma(2,2);
  sigmu2 ~ gamma(2,.02);
  
  //sigalphaT2 ~ gamma(1000,1);
  //sigalpha2 ~ gamma(1000,1);
  //sigalphaT2 ~ gamma(2,2);
  //sigalpha2 ~ gamma(2,.2);
  xi ~ exponential(.01);
  
  mu[n_weeks] ~ normal(0, sqrt(1/sigmuT2));
  
  
  for (i in 1:(n_weeks - 1)) {
    mu[i] ~ normal(mu[i + 1], sqrt(1/sigmu2));
  }
  
  //alpha[1] ~ normal(0, 1/sigalphaT2);
  //for (i in 2:n_weeks) {
  //  alpha[i] ~ normal(alpha[i - 1], 1/sigalpha2);
  //}
  
  alpha0 ~ normal(0,20000);
  alpha1 ~ normal(0,20000);
  
  
  
  //sampling distribution
  
  
  //ili ~ beta_proportion(IS, phi);
  
  for (i in 1:M) {
    
    
     ili[i] ~ beta_proportion(
                IS[i]*exp(mu[weeks[i]])/
                (1 + (exp(mu[weeks[i]]) - 1)*IS[i]), 
                              phi);
                              
    ili[i] ~ beta_proportion(inv_logit(logit(IS[i]) + mu[weeks[i]]), 
                              phi);
  }
  //print(seg)
  for (i in 1:cur_yr_n_weeks) {
  //  hosp[i] ~ neg_binomial_2(exp(alpha[i]*segment(ili, seg_ind_start[ps], 
  //                                        seg_ind_length[ps])[i]),
  //                                        xi);
    //print("i = ", i);
    //print(100*segment(ili, seg_ind_start[ps], seg_ind_length[ps])[i]);
    //hosp[i] ~ normal(alpha[i]*(ili_ps[i]*100), xi);
    // hosp[i] ~ neg_binomial_2_log(alpha0 + alpha1*(100*ili[i]), xi);
    ////hosp[i] ~ normal(alpha0 + alpha1*(100*ili[i]), xi);
    hosp[i] ~ normal(alpha0 + alpha1*logit(ili_ps[i]), xi);
  }
  
  //hosp ~ normal(alpha0 + alpha1*ili_ps, xi);
}
generated quantities {
  real r0 = beta[ps] / gamma[ps];
  real recovery_time = 1 / gamma[ps];
  //real yp[43];
  real pred_ili[n_weeks];
  real pred_ili_asg[n_weeks];
  real pred_hosp[n_weeks];
  real discrepancy[n_weeks];
  real tp[n_weeks];
  //real thetap[2];
  real ypm[n_weeks,3];
  vector[n_weeks] Ip;
  real yp[3];
  //real pred_ili_wo_disc[27];
  yp = {S0, I0[ps], 1 - S0 - I0[ps]};

  // tp = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27};
  for (i in 1:n_weeks) tp[i] = i;

  ypm = integrate_ode_rk45(sir, yp, t0, tp, theta[ps,], x_r, x_i);
  Ip = col(to_matrix(ypm) ,2);
  
  for (i in 1:n_weeks) {
    //pred_ili_wo_disc[i] = beta_proportion_rng(Ip[i], phi);
    pred_ili[i] = beta_proportion_rng(
                        Ip[i]*exp(mu[i])/
                        (1 + (exp(mu[i]) - 1)*Ip[i]), 
                        phi);
                        
    pred_ili_asg[i] = beta_proportion_rng(Ip[i], 
                        phi);
                        
    // pred_hosp[i] = normal_rng(alpha0 + alpha1*(100*pred_ili[i]), xi);
    // pred_hosp[i] = neg_binomial_2_log_rng(alpha0 + alpha1*(100*pred_ili[i]), xi);
    // pred_hosp[i] = normal_rng(alpha[i]*(100*pred_ili[i]), xi);
    pred_hosp[i] = normal_rng(alpha0 + alpha1*logit(pred_ili[i]), xi);
    discrepancy[i] = mu[i];
  }
  
  //pred_hosp = normal_rng(alpha0 + alpha1*(100*pred_ili), xi);


}




