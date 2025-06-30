
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


data {
int<lower=0> HM;
int<lower=0> n_seasons_hosp;
array[HM] real<lower=0> ili_ps;                           // all ili data
array[HM] real hosp;
array[HM] int<lower=0> hosp_seasons;
real<lower=0> count_rate;

//hyperparameters
real<lower=0> sigma_alpha0; // 5
real<lower=0> sigma_alpha1; // 5
real<lower=0> sigma_alpha3; // .4
real<lower=0> sigma_sigma_epsilon; // 4

}
parameters {

vector[n_seasons_hosp] alpha0;
vector[n_seasons_hosp] alpha1;
vector[n_seasons_hosp] alpha2;
vector<lower=-1,upper=1>[n_seasons_hosp] alpha3;
vector<lower=0>[n_seasons_hosp] sigma_epsilon;
//real alpha0;
//real alpha1;
//real alpha2;
//real<lower=-1, upper=1> alpha3;
//real<lower=0> sigma_epsilon;

}


model {

alpha0 ~ normal(0, sigma_alpha0);
alpha1 ~ normal(0, sigma_alpha1);
alpha2 ~ normal(0, sigma_alpha1);
sigma_epsilon ~ normal(0, sigma_sigma_epsilon);
alpha3 ~ normal(0,sigma_alpha3);


                                  

hosp[1] ~ normal(alpha0[hosp_seasons[1]] +
                 alpha1[hosp_seasons[1]]*(count_rate*ili_ps[1]) +
                 alpha2[hosp_seasons[1]]*(count_rate*ili_ps[1])^2,
			count_rate*sigma_epsilon[hosp_seasons[1]]);

for (i in 2:HM) hosp[i] ~ normal(alpha0[hosp_seasons[i]] + 
                                 alpha1[hosp_seasons[i]]*(count_rate*ili_ps[i]) +
                                 alpha2[hosp_seasons[i]]*(count_rate*ili_ps[i])^2 +
				 alpha3[hosp_seasons[i]]*hosp[i-1], 
                                             count_rate*sigma_epsilon[hosp_seasons[i]]);         

                          
}




