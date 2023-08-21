
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
real<lower=0> sigma_tau; // .01
real<lower=0> sigma_taumu; // .05
real<lower=0> sigma_kappa; // 1000
real<lower=0> df; // 4
real<lower=0> sigma_alpha0; // 5
real<lower=0> sigma_alpha1; // 5
real<lower=0> sigma_xi; // 3
real<lower=0> sigma_eta; // .4

}
parameters {


vector[n_params] theta;
matrix[n_seasons,n_params] theta_s;
real<lower=0> tau;
real<lower=0> taumu;
real<lower=0> kappas[n_seasons];
real<lower=0> kappas2[n_seasons];


row_vector[n_params] zeta;
corr_matrix[n_params] omega;
real gamma[n_weeks - 1];

real alpha0;
real alpha1;
real<lower=0> sigma_epsilon;
real<lower=0> sigma_epsilon2;

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
tau ~ normal(0, .01);
taumu ~ normal(0, .05);

zeta ~ student_t(df, c, d);
omega ~ lkj_corr(nu);

alpha0 ~ normal(0, sigma_alpha0);
alpha1 ~ normal(0, sigma_alpha1);
sigma_epsilon ~ exponential(.81);
sigma_epsilon2 ~ normal(0, 4);

for (i in 1:n_seasons) {
                          
                          kappas[i] ~ normal(0, sigma_kappa);
                          kappas2[i] ~ exponential(.001);
}
}




