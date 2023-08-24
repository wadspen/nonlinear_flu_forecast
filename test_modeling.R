setwd("~/flu_research/nonlinear_flu_forecast/")
source('./get_data_functions.R')
source('./mle_functions.R')
library(tidyr)
ILINet <- get_ili_data()
both_flu <- comb_data()

args <- commandArgs()
# mod <- cmdstan_model(stan_file = './mf_asg_ili_link_mod.stan')
mod2 <- cmdstan_model(stan_file = "./asg_logit.stan")
mod3 <- cmdstan_model(stan_file = "./sir_logit.stan")
mod3 <- cmdstan_model(stan_file = "./mf_sir2_disc.stan")
#state <- "California"
#state <- args[6]

state <- "Texas"
print(state)
i <- 18

dat <- get_stan_data(ILINet, both_flu, m_week = i, s_region = state)
stan_dat <- dat[[1]]
info_list <- dat[[2]]

init <- list(theta = stan_dat$m0, theta_s = info_list$theta_s) #, alpha0 = info_list$alpha[1], alpha1 = info_list$alpha[2])

print(stan_dat$beta)
print(stan_dat$m0)
print(info_list$theta_s)


  
  samps <- mod3$sample(data = stan_dat,
                      chains = 1,
                      # seed = 0,
                      # adapt_delta = .999,
                      init = list(init),
                      iter_warmup = 500,
                      iter_sampling = 500)
  
  
  
  draws <- samps$draws(format = 'df')
  num <- ",13]"
  lnum <- nchar(num)
  postpr <- draws[,which(str_detect(colnames(draws),pattern = 'pred_hosp'))]
  
  # postpr <- postpr[,which(str_detect(colnames(postpr),pattern = '8]'))]
  # postpr <- postpr[,which(str_detect(str_sub(colnames(postpr),
  #                                            start = -lnum), pattern = num))]
  postpr <- postpr[,1:53]
  # bayesplot::mcmc_trace(draws, pars = 'theta.2.')
  bayesplot::mcmc_hist(draws, pars = 'sigma_epsilon')
  # bayesplot::mcmc_trace(draws, pars = 'alpha1')
  # bayesplot::mcmc_trace(draws, pars = 'epsilons.13.')
  # bayesplot::mcmc_trace(draws, pars = 'theta_s[13,1]')
  # bayesplot::mcmc_trace(draws, pars = 'epsilon')
  postpr <- postpr[is.na(postpr[,1]) == FALSE,1:53]
  pi <- as.data.frame(predictive_interval(as.matrix(postpr), prob=.95))
  colnames(pi) <- c('lower','upper')
  pi$season_week <- 1:53
  
  state <- "Texas"
  scale <- unique(both_flu$pop_2022[both_flu$region == state])
  both_flu %>%
    filter(season == 2022, region == state) %>%
    ggplot() +
    geom_point(aes(x=(season_week), y = value)) +
    # geom_point(aes(x=(season_week), y = unweighted_ili/100)) +
    # geom_line(data=pi[1:20,],aes(x=season_week, y=exp(upper)*scale[2]),colour='orange') +
    # geom_line(data=pi[1:20,],aes(x=season_week, y=exp(lower)*scale[2]),colour='orange') +
    geom_line(data=pi[,],aes(x=season_week, y=exp(upper)),colour='orange') +
    geom_line(data=pi[,],aes(x=season_week, y=exp(lower)),colour='orange') +
    geom_vline(xintercept = 25) +
    theme_bw()
  
  

  samps$draws()
  samps$sampler_diagnostics()
  file <- paste(paste("test_stan_forc_link", state, i,
                      sep = "_"), ".csv", sep = "")

  file_path <- paste("./", file, sep = "")
  # saveRDS(samps, file_path)

  draws <- samps$draws(format = 'df')
  write.csv(draws, file_path)
