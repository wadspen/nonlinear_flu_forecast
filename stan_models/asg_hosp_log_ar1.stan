
// In this STAN model, ILI data is modeled according to the assymmetric
// Gaussian function where there is a hierarchy of parameters over flu
// seasons. The ASG function is taken as the expected value of a beta
// distribution and the scale parameter of the beta distribution is 
// considered independent over seasons.
// A discrepancy reverse random walk is included to the mean to capture
// any pattern that is shared by all or most flu seasons. 
// Flu hospitalizations are then modeled as a linear function of the ILI
// model, the squared ILI model, and one lagged hospitalization. The 
// parameters in the linear model between the two seasons, 2022 and 2023,
// are independent.

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
int<lower=0> HM;
//int<lower=0> n_seasons_hosp;
// array[HM] real<lower=0> ili_ps;                           // all ili data
array[HM] real hosp;
//array[HM] int<lower=0> hosp_seasons;
// real<lower=0> count_rate;

//hyperparameters
real<lower=0> sigma_alpha0; // 5
real<lower=0> sigma_alpha1; // 5
real<lower=0> sigma_alpha3; // .4
real<lower=0> sigma_sigma_epsilon; // 4


//asg data and hyperparameters
int<lower=0> M;                                 // number of points
int<lower=0> cur_yr_n_weeks;                    // last observed week
int<lower=0> n_seasons;                         // number of seasons
array[M] int<lower=0> all_seasons;                    // all_seasons
array[M] real<lower=0, upper=1> ili;                           // all ili data
// int<lower=0> HM;                                // number of hospitalization points
int<lower=0> n_seasons_hosp;
// array[HM] real hosp;
array[HM] int<lower=0> hosp_seasons;
array[HM] int<lower=0> hosp_season_weeks;
array[M] int<lower=0> weeks;                          // week
int<lower=0> n_params;                          // number of parameters ASG
vector[n_params] m0;                            // prior means
matrix[n_params,n_params] C0;                   // prior sds
matrix[n_params,n_params] d_mat;                // diagonal identity matrix
array[n_seasons] real beta;
real nu;
real<lower=0> c;
real<lower=0> d;
int<lower=0> n_weeks;                                 // weeks in a season


int<lower=0> ps;                                // season to predict
int<lower=0> hps;
// real<lower=0> count_rate;

//hyperparameters
real<lower=0> sigma_sigma_gamma_w; // .01
real<lower=0> sigma_sigma_gamma; // .05
real<lower=0> sigma_kappa; // 1000
real<lower=0> df; // 4
// real<lower=0> sigma_alpha0; // 5
// real<lower=0> sigma_alpha1; // 5
// real<lower=0> sigma_alpha3; // .4
real<lower=0> sigma_xi; // 3
real<lower=0> sigma_eta; // .4
// real<lower=0> sigma_sigma_epsilon; // 4
real<lower=0> sigma_disc;
real<lower=0> sigma_gamma_W;

}
parameters {

real alpha0;
real alpha1;
real<lower=-1, upper=1> alpha3;
real<lower=0> sigma_epsilon;


//asg parameters
vector[n_params] theta;
matrix[n_seasons,n_params] theta_s;
array[n_seasons] real<lower=0> kappas;


vector<lower=0>[n_params] zeta;

}


model {

alpha0 ~ normal(0, sigma_alpha0);
alpha1 ~ normal(0, sigma_alpha1);
//alpha2 ~ normal(0, sigma_alpha1);
sigma_epsilon ~ normal(0, sigma_sigma_epsilon);
alpha3 ~ normal(0,sigma_alpha3);

///////////////////////////////////////////////////
///////////////ili model part//////////////////////
//////////////////////////////////////////////////

theta ~ multi_normal(m0,C0);
zeta ~ normal(0, c);


                    
for (i in 1:n_seasons) {
                 
			  theta_s[i,] ~ multi_normal(theta, diag_matrix(square(zeta)));
                          kappas[i] ~ normal(0, sigma_kappa);
}



for (i in 1:M) ili[i] ~ beta_proportion(
                                    inv_logit(asg(theta_s[all_seasons[i],],
                                    beta[all_seasons[i]],
                                    weeks[i])),
                                    kappas[all_seasons[i]]);
                                    
////////////////////////////////////////////////
///////////////////////////////////////////////
///////////////////////////////////////////////


                                  

hosp[1] ~ normal(alpha0 +
                 alpha1*(ili[(M - HM) + 1]),// +
                 //alpha2*(ili_ps[1])^2,
			 sigma_epsilon);

for (i in 2:HM) hosp[i] ~ normal(alpha0 + 
                                 alpha1*(ili[(M - HM) + i]) +
                                 //alpha2*(ili_ps[i])^2 +
				 alpha3*hosp[i-1], 
                                             sigma_epsilon);         

                          
}

generated quantities {
    array[cur_yr_n_weeks + 5] real<lower=0,upper=1> pred_ili;
    array[cur_yr_n_weeks + 5] real pred_hosp;
            
    
    // for (i in 1:n_seasons) {
      for (j in 1:(cur_yr_n_weeks + 5)) {
        pred_ili[j] = beta_proportion_rng(
                              inv_logit(asg(theta_s[n_seasons,],
                              beta[n_seasons],
                              j)),
                              kappas[n_seasons]);
                                                         
               
        
      }
    // }
    
    pred_hosp[1] = normal_rng(alpha0 + alpha1 * (pred_ili[1]), 
                              sigma_epsilon);
    for (i in 2:(cur_yr_n_weeks + 5)) {
      pred_hosp[i] = normal_rng(alpha0 + alpha1 * (pred_ili[i]) +
                                alpha3 * pred_hosp[i - 1], sigma_epsilon);
    }
    
    
}




