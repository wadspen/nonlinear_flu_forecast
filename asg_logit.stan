
// This Stan program defines a model where a vector hosp is modeled as linear 
// model function of logit(ili). ili is modeled hierarchically as beta 
// distributed with mean a function ASG(week).


functions {

   real asg(row_vector theta, real beta, real x) {
    
    real eta;
    real mu; 
    real sig1;
    real sig2;
    real ASG;
    eta = exp(theta[1]);
    mu = theta[2];
    sig1 = exp(theta[3]);
    sig2 = exp(theta[4]);
    
    ASG = ((beta + (eta)*exp(-((x-mu)^2)/(2*sig1^2)))*(x < mu) +
       (beta + (eta)*exp(-((x-mu)^2)/(2*sig2^2)))*(mu <= x));
    return ASG;
  }
}

data {
int<lower=0> M;                                 // number of points
int<lower=0> cur_yr_n_weeks;                    // last observed week
int<lower=0> n_seasons;                         // number of seasons
int<lower=0> all_seasons[M];                    // all_seasons
real<lower=0> ili[M];                           // all ili data
real<lower=0> hosp[cur_yr_n_weeks];
int<lower=0> weeks[M];                          // week
int<lower=0> n_params;                          // number of parameters ASG
vector[n_params] m0;                            // prior means
matrix[n_params,n_params] C0;                   // prior sds
matrix[n_params,n_params] d_mat;                // diagonal identity matrix
real beta[n_seasons];
real nu;
real<lower=0> c;
real<lower=0> d;
int<lower=0> n_weeks;                                 // weeks in a season


int<lower=0> ps;                                // season to predict


//hyperparameters
real<lower=0> sigma_sigma_gamma_w; // .01
real<lower=0> sigma_sigma_gamma; // .05
real<lower=0> sigma_kappa; // 1000
real<lower=0> df; // 4
real<lower=0> sigma_alpha0; // 5
real<lower=0> sigma_alpha1; // 5
real<lower=0> sigma_xi; // 3
real<lower=0> sigma_eta; // .4
real<lower=0> sigma_sigma_epsilon; // 4

}
parameters {


vector[n_params] theta;
matrix[n_seasons,n_params] theta_s;
real<lower=0> sigma_gamma_w;
real<lower=0> sigma_gamma;
real<lower=0> kappas[n_seasons];


row_vector[n_params] zeta;
corr_matrix[n_params] omega;
real gamma[n_weeks - 1];

real alpha0;
real alpha1;
real<lower=0> sigma_epsilon;

}

transformed parameters {
matrix[n_params,n_params] Sigma;
real ili_ps[cur_yr_n_weeks];
real gamma1[1]; //test1
real gam[n_weeks]; //test1

Sigma = quad_form_diag(omega,zeta);
ili_ps = ili[(M - cur_yr_n_weeks + 1):M];
gamma1[1] = -sum(gamma);
gam = append_array(gamma1, gamma);
}

model {

theta ~ multi_normal(m0,C0);
sigma_gamma_w ~ normal(0, sigma_sigma_gamma_w);
sigma_gamma ~ normal(0, sigma_sigma_gamma);

zeta ~ student_t(df, c, d);
omega ~ lkj_corr(nu);

alpha0 ~ normal(0, sigma_alpha0);
alpha1 ~ normal(0, sigma_alpha1);
// sigma_epsilon ~ exponential(.81);
sigma_epsilon ~ normal(0, sigma_sigma_epsilon);


                    
for (i in 1:n_seasons) {
                          theta_s[i,] ~ multi_normal(theta, Sigma);
                          // kappas[i] ~ exponential(.001);
                          kappas[i] ~ normal(0, sigma_kappa);
}

gamma[n_weeks - 1] ~ normal(0,sigma_gamma_w);

for (i in 1:(n_weeks-2)) gamma[i] ~ normal(gamma[i+1],sigma_gamma);

for (i in 1:M) ili[i] ~ beta_proportion(
                                    inv_logit(asg(theta_s[all_seasons[i],], 
                                    beta[all_seasons[i]],
                                    weeks[i]) +
                                    gam[weeks[i]]),
                                    kappas[all_seasons[i]]);


for (i in 1:cur_yr_n_weeks) hosp[i] ~ normal(alpha0 + alpha1*logit(ili_ps[i]), 
                                             sigma_epsilon);
                                    
// for (i in 1:cur_yr_n_weeks) hosp[i] ~ normal(alpha0 + alpha1*(log(ili_ps[i]) - 
//                                               log(1 - ili_ps[i])), 
//                                              sigma_epsilon);

                          
}



generated quantities {
    real pred_ili[n_weeks, n_seasons];
    real pred_ili_asg[n_weeks, n_seasons];
    real discrepancy[n_weeks];
    real pred_hosp[n_weeks];
    
    
    for (i in 1:n_seasons) {
      for (j in 1:n_weeks) {
        pred_ili[j,i] = beta_proportion_rng(
                              inv_logit(asg(theta_s[i,], 
                              beta[i],
                              j) +
                              gam[j]),
                              kappas[i]);
        
                                  
        pred_ili_asg[j,i] = beta_proportion_rng(
                              inv_logit(asg(theta_s[i,], 
                              beta[i],
                              j)),
                              kappas[i]);
        
      }
    }
    
    for (i in 1:n_weeks) {
      
      pred_hosp[i] = normal_rng(alpha0 + alpha1*logit(pred_ili[i,ps]), 
                                sigma_epsilon);
      
      // pred_hosp[i] = normal_rng(alpha0 + alpha1*(log(pred_ili[i,ps]) - 
      //                                           log(1 - pred_ili[i,ps])), sigma_epsilon);
      
      discrepancy[i] = gam[i];
    }
    
}

