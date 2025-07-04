
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
int<lower=0> M;                                 // number of points
int<lower=0> cur_yr_n_weeks;                    // last observed week
int<lower=0> n_seasons;                         // number of seasons
array[M] int<lower=0> all_seasons;                    // all_seasons
array[M] real<lower=0, upper=1> ili;                           // all ili data
int<lower=0> HM;                                // number of hospitalization points
int<lower=0> n_seasons_hosp;
array[HM] real hosp;
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
real<lower=0> count_rate;

//hyperparameters
real<lower=0> sigma_sigma_gamma_w; // .01
real<lower=0> sigma_sigma_gamma; // .05
real<lower=0> sigma_kappa; // 1000
real<lower=0> df; // 4
real<lower=0> sigma_alpha0; // 5
real<lower=0> sigma_alpha1; // 5
real<lower=0> sigma_alpha3; // .4
real<lower=0> sigma_xi; // 3
real<lower=0> sigma_eta; // .4
real<lower=0> sigma_sigma_epsilon; // 4
real<lower=0> sigma_disc;
real<lower=0> sigma_gamma_W;

}
parameters {


vector[n_params] theta;
matrix[n_seasons,n_params] theta_s;
real<lower=0> sigma_gamma_w;
real<lower=0> sigma_gamma;
array[n_seasons] real<lower=0> sigma_upsilon;
real<lower=0> sigma_sigma_upsilon;
array[n_seasons] real<lower=0> kappas;

vector<lower=0>[n_params] zeta;
array[n_weeks - 1] real gamma;

array[n_weeks - 1, n_seasons] real upsilon;

}

transformed parameters {
array[1] real gamma1;
array[n_weeks] real gam;

array[1, n_seasons] real upsilon1;
array[n_weeks, n_seasons] real ups;

gamma1[1] = -sum(gamma);
gam = append_array(gamma1, gamma);

for (i in 1:n_seasons) upsilon1[1, i] = 0;//-exp(beta[i])/(1 + exp(beta[i]));
ups = append_array(upsilon, upsilon1);

}

model {

theta ~ multi_normal(m0,C0);
sigma_gamma ~ normal(0, sigma_sigma_gamma);
sigma_gamma_w ~ normal(sigma_disc, sigma_gamma_W);
zeta ~ normal(0, c);
sigma_upsilon ~ normal(0, sigma_sigma_upsilon);
sigma_sigma_upsilon ~ normal(0, sigma_gamma_W);

                    
// for (i in 1:n_seasons) {
//                  
// 			  theta_s[i,] ~ multi_normal(theta, diag_matrix(square(zeta)));
//                           kappas[i] ~ normal(0, sigma_kappa);
// }

gamma[n_weeks - 1] ~ normal(0, sqrt(sigma_gamma_w));
for (i in 1:(n_weeks-2)) gamma[i] ~ normal(gamma[i+1], sqrt(sigma_gamma));

// gamma[1:(n_weeks-1)] ~ normal(gamma[2:(n_weeks)], sqrt(sigma_gamma));

// for (i in 1:n_seasons) {
//   upsilon[n_weeks - 1, i] ~ normal(upsilon1[1, i], sigma_upsilon[i]);
// }

upsilon[n_weeks - 1, 1:n_seasons] ~ normal(upsilon1[1, 1:n_seasons],
                                           sqrt(sigma_upsilon[1:n_seasons]));
                                           
// upsilon[1, 1:n_seasons] ~ normal(0, sqrt(sigma_upsilon[1:n_seasons]));

// for (i in 1:(n_weeks-2)) {
//   for (j in 1:n_seasons) {
//     upsilon[i, j] ~ normal(upsilon[i + 1, j], sigma_upsilon[j]);
//   }
// }

// for (i in 1:(n_weeks-2)) {
  for (j in 1:n_seasons) {
    upsilon[1:(n_weeks-2), j] ~ normal(upsilon[2:(n_weeks-1), j],
                                      sqrt(sigma_upsilon[j]));
    
    // upsilon[2:(n_weeks), j] ~ normal(to_vector(upsilon[1:(n_weeks-1), j])*alpha_up[j], 
    //                                   sqrt(sigma_upsilon[j]));
                                      
    theta_s[j,] ~ multi_normal(theta, diag_matrix(square(zeta)));
                          kappas[j] ~ normal(0, sigma_kappa);
  }
// }

// for (i in 1:(M - cur_yr_n_weeks)) ili[i] ~ beta_proportion(
//                                     inv_logit(asg(theta_s[all_seasons[i],],
//                                     beta[all_seasons[i]],
//                                     weeks[i]) +
//                                     gam[weeks[i]]),
//                                     kappas[all_seasons[i]]); 
        


for (i in 1:M) ili[i] ~ beta_proportion(
                                    inv_logit(asg(theta_s[all_seasons[i],],
                                    beta[all_seasons[i]],
                                    weeks[i]) +
                                    gam[weeks[i]] + ups[weeks[i],
                                                        all_seasons[i]]),
                                    kappas[all_seasons[i]]); 
        
}



generated quantities {
    array[cur_yr_n_weeks + 5, n_seasons] real<lower=0,upper=1> pred_ili;
    array[cur_yr_n_weeks + 5, n_seasons] real pred_ili_asg;
    // array[cur_yr_n_weeks + 5] real discrepancy; //as opposed to n_weeks
    // array[cur_yr_n_weeks + 5] real discrepancy2;    
    
    for (i in 1:n_seasons) {
      for (j in 1:(cur_yr_n_weeks + 5)) {
        pred_ili[j,i] = beta_proportion_rng(
                              inv_logit(asg(theta_s[i,],
                              beta[i],
                              j) +
                              gam[j] + ups[j, all_seasons[i]]),
                              kappas[i]);
                              
        
        
  
        pred_ili_asg[j,i] = beta_proportion_rng(
                             inv_logit(asg(theta_s[i,],
                             beta[i],
                             j)),
                             kappas[i]);
                              
        
        
        
      }
    }
    
    // discrepancy[1] = gam[1];
    // discrepancy2[1] = ups[1];

    // for (i in 2:(cur_yr_n_weeks + 5)) {
    //   discrepancy[i] = gam[i];
      // discrepancy2[i] = ups[i];
    // }
    
}

