// Optimized Stan model: focus on reducing redundant computations and using non-centered parameterization

functions {
  real asg(row_vector theta, real x) {
    real beta = theta[1];
    real eta = exp(theta[2]);
    real mu = theta[3];
    real sig1 = exp(theta[4]);
    real sig2 = exp(theta[5]);

    real dx = x - mu;
    real term = dx * dx;

    return (beta + eta * exp(-term / (2 * sig1^2))) * (x < mu) +
           (beta + eta * exp(-term / (2 * sig2^2))) * (x >= mu);
  }
}

data {
  int<lower=0> M;
  int<lower=0> cur_yr_n_weeks;
  int<lower=0> n_seasons;
  array[M] int<lower=0> all_seasons;
  array[M] real<lower=0, upper=1> ili;
  int<lower=0> HM;
  int<lower=0> n_seasons_hosp;
  array[HM] real hosp;
  array[HM] int<lower=0> hosp_seasons;
  array[HM] int<lower=0> hosp_season_weeks;
  array[M] int<lower=0> weeks;
  int<lower=0> n_params;
  vector[n_params] m0;
  matrix[n_params,n_params] C0;
  matrix[n_params,n_params] d_mat;
  array[n_seasons] real beta;
  real nu;
  real<lower=0> c;
  real<lower=0> d;
  int<lower=0> n_weeks;
  int<lower=0> ps;
  int<lower=0> hps;
  real<lower=0> count_rate;
  real<lower=0> sigma_sigma_gamma_w;
  real<lower=0> sigma_sigma_gamma;
  real<lower=0> sigma_kappa;
  real<lower=0> df;
  real<lower=0> sigma_alpha0;
  real<lower=0> sigma_alpha1;
  real<lower=0> sigma_alpha3;
  real<lower=0> sigma_xi;
  real<lower=0> sigma_eta;
  real<lower=0> sigma_sigma_epsilon;
  real<lower=0> sigma_disc;
  real<lower=0> sigma_gamma_W;
}

parameters {
  vector[n_params] theta;
  matrix[n_seasons, n_params] z_theta_s;
  real<lower=0> sigma_gamma_w;
  real<lower=0> sigma_gamma;
  array[n_seasons] real<lower=0> kappas;
  real<lower=0> sigma_upsilon;
  vector<lower=0>[n_params] zeta;
  array[n_weeks - 1] real gamma;
  array[n_weeks - 1] real upsilon;
}

transformed parameters {
  array[1] real gamma1;
  array[n_weeks] real gam;
  array[1] real upsilon1;
  array[n_weeks] real ups;
  matrix[n_seasons, n_params] theta_s;
  matrix[n_seasons, n_weeks] asg_vals;

  gamma1[1] = -sum(gamma);
  gam = append_array(gamma1, gamma);

  upsilon1[1] = 0;
  ups = append_array(upsilon, upsilon1);

  for (i in 1:n_seasons) {
    theta_s[i] = theta + zeta .* z_theta_s[i];
  }

  for (i in 1:n_seasons)
    for (j in 1:n_weeks)
      asg_vals[i, j] = asg(theta_s[i], j);
}

model {
  theta ~ multi_normal(m0, C0);
  sigma_gamma ~ normal(0, sigma_sigma_gamma);
  sigma_gamma_w ~ normal(sigma_disc, sigma_gamma_W);
  zeta ~ normal(0, c);
  sigma_upsilon ~ normal(0, sigma_gamma_W);

  to_vector(z_theta_s) ~ normal(0, 1);

  for (i in 1:n_seasons)
    kappas[i] ~ normal(0, sigma_kappa);

  gamma ~ normal(0, sqrt(sigma_gamma_w));
  for (i in 1:(n_weeks-2))
    gamma[i] ~ normal(gamma[i+1], sqrt(sigma_gamma));

  upsilon ~ normal(0, sigma_upsilon);
  for (i in 1:(n_weeks-2))
    upsilon[i] ~ normal(upsilon[i+1], sigma_upsilon);

  for (i in 1:(M - cur_yr_n_weeks))
    ili[i] ~ beta_proportion(inv_logit(asg_vals[all_seasons[i], weeks[i]] + gam[weeks[i]]),
                             kappas[all_seasons[i]]);

  for (i in (M - cur_yr_n_weeks + 1):M)
    ili[i] ~ beta_proportion(inv_logit(asg_vals[all_seasons[i], weeks[i]] + gam[weeks[i]] + ups[weeks[i]]),
                             kappas[all_seasons[i]]);
}

generated quantities {
  array[cur_yr_n_weeks + 5, n_seasons] real<lower=0,upper=1> pred_ili;

  for (i in 1:n_seasons)
    for (j in 1:(cur_yr_n_weeks + 5))
      pred_ili[j, i] = beta_proportion_rng(
        inv_logit(asg_vals[i, j] + gam[j] + ups[j]), kappas[i]);
}
