library(cmdstanr)
library(lubridate)
library(ggplot2)
library(posterior)
library(bayesplot)
library(rstanarm)
library(stringr)
library(dplyr)
library(cdcfluview)
get_hosp_data <- function(lag = 1) {
  raw_flu <- read.csv(
    'https://raw.githubusercontent.com/cdcepi/Flusight-forecast-data/master/data-truth/truth-Incident%20Hospitalizations.csv')
  
  raw_flu <- raw_flu %>% 
    mutate(year = year(date),week = week(date) - lag, 
           week2 = week(date) - (lag - 1))
  
  return(raw_flu)
}

get_ili_data <- function() {
  ILINet_state <- ilinet(region = 'state')
  ILINet_us <- ilinet(region = 'national') %>% mutate(region = 'US')
  ILINet <- rbind(ILINet_state,ILINet_us)
  
  # ILINet <- ILINet %>% 
  #   mutate(year=year(week_start), week2=week(week_start),location_name=region)
  
  ILINet <- ILINet %>% 
    mutate(year=year(week_start), week2=week(week_start),location_name=region) #%>% 
  # filter(month(week_start) %in% c(1:6,8:12))
  
  # years53 <- ILINet %>% 
  #   group_by(year) %>% 
  #   summarise(
  #     m = max(week2)
  #   ) %>% 
  #   filter(m > 52) %>% 
  #   select(year) %>% 
  #   as.vector()
  # years53 <- years53$year
  
  # ILINet <- ILINet %>% 
  #   filter(week2 %in% c(38:53,1:25)) %>% 
  #   mutate(season=ifelse(week2 %in% 38:53,year,year-1)) %>% 
  #   mutate(season_week = ifelse(week2 %in% 38:53, week2-37,
  #                               ifelse(season %in% years53, 
  #                                      week2 + 16, week2 + 15))) %>% 
  #   filter(season_week < max(season_week))
  
  ILINet <- ILINet %>% 
    mutate(season = ifelse(month(week_start) <= 7, year - 1, year),
           place = 1) %>% 
    group_by(season, region) %>% 
    mutate(season_week = cumsum(place),
           unweighted_ili = ifelse(unweighted_ili == 0, .00001,
	   			   unweighted_ili)) %>%
    filter(is.na(unweighted_ili) == FALSE)
  
  return(ILINet)
}

comb_data <- function(lag = 1, s_region = NA, s_season = 2022) {
  raw_flu <- get_hosp_data(lag = lag)
  ILINet <- get_ili_data()
  
  if (is.na(s_region)) {
    both_flu <- raw_flu %>% 
      left_join(ILINet, by = c('year','week','location_name'))
  }
  
  else {
    both_flu <- raw_flu %>% 
      left_join(ILINet, by = c('year','week','location_name')) %>% 
      filter(region == s_region)
  }
 
  pops <- read.csv("./state_populations.csv")
  both_flu <- both_flu %>%
	 left_join(pops, by = "region") %>%
	 mutate(unweighted_ili = ifelse(unweighted_ili == 0, .00001,
					unweighted_ili)) %>%
 	 filter(is.na(unweighted_ili) == FALSE, is.na(value) == FALSE)

  return(both_flu)
  
}

make_stan_data <- function(dat1, dat2, ili_seasons,
                           nu=1, df = 4, c = 2, d = 4, S0 = .9, 
                           sigma_sigma_gamma_w = .01, sigma_sigma_gamma = .05,
                           sigma_kappa = 10000, sigma_alpha0 = 5,
                           sigma_alpha1 = 5, sigma_xi = 3, sigma_eta = .4,
                           sigma_sigma_epsilon = 4,
                           traj = 'asg') {
  
  start <- c(-4.2533962, 0.3416948, 22.5023919, 2.3580754, 
	     2.3519650, 100.3432742)
  pars <- asg_max_lik(dat2$unweighted_ili/100, dat2$season_week, theta = start)
  state <- unique(dat1$region)
  seasons <- unique(dat1$season)
  season <- unique(dat2$season)
  #alpha <- coef(lm(dat2$value/dat2$pop_2022 ~ dat2$unweighted_ili))
  alpha <- coef(lm(log(dat2$value + 1) ~ dat2$unweighted_ili^2))
  if (traj == 'asg') {
    m0 <- c(.3, 23, 3.69, 4.7)
    C0 <- c(.2, 5, .9, .9)
    n_params <- length(m0)
    all_season_region_mles <- read.csv("./all_season_region_mles_asg.csv") %>% 
      select(-X) %>% 
      filter(region == state, season %in% ili_seasons)
    all_season_mles <- read.csv("./all_season_mles_asg.csv") %>%
      filter(region == state)
    as_beta <- as.vector(all_season_mles$beta1)
    # print(ili_seasons)
    # beta <- all_season_region_mles$beta[all_season_region_mles$region == state
    #                                     & all_season_region_mles$region == state]
    beta <- all_season_region_mles %>%
      select(beta) %>% 
      unlist() %>% 
      as.vector()
    #beta[length(seasons)] <- mean(beta[1:(length(seasons)-1)])
    beta[length(seasons)] <- pars[1]
    theta_s <- all_season_region_mles %>% 
      select(eta, mu, sigma1, sigma2, kappa) %>% 
      as.matrix()
    theta_s[length(seasons),] <- apply(theta_s[1:(length(seasons) - 1),], MARGIN = 2, FUN = mean)
  }
  
  else if (traj == 'uln') {
    m0 <- c(-4, -2.6, 2.4, log(.4))
    C0 <- c(1, .5, 4, 4)
    n_params <- length(m0)
    all_season_region_mles <- read.csv("./all_season_region_mles_uln.csv") %>% 
      select(-X) %>% 
      filter(region == state, season %in% ili_seasons)
    beta <- all_season_region_mles %>%
      select(beta) %>% 
      unlist() %>% 
      as.vector()
    theta_s <- all_season_region_mles %>% 
      select(eta, mu, sigma, kappa) %>% 
      as.matrix()
  }
  
  
  
  week_inds <- dat1 %>% 
    group_by(season) %>% 
    summarise(
      max_week = max(season_week)
    )
  is_last <- function(x) !duplicated(x, fromLast = TRUE)
  week_inds$index <- which(is_last(dat1$season))
  dat_inds <- week_inds %>% 
    mutate(for_ind = index - max_week + 1)
  
  
  
  data_list <- list(
    #for all models
    M = nrow(dat1),
    n_seasons = length(unique(dat1$season)),
    n_weeks = max(ILINet$season_week),
    cur_yr_n_weeks = max(dat2$season_week),
    all_seasons = as.numeric(factor(dat1$season)),
    weeks = dat1$season_week,
    ps = max(as.numeric(factor(unique(dat1$season)))),
    ili = dat1$unweighted_ili/100,
    #hosp = log((dat2$value + 1)/dat2$pop_2022),
    #hosp = dat2$value,
    hosp = log(dat2$value),
    sigma_sigma_gamma_w = sigma_sigma_gamma_w,
    sigma_sigma_gamma = sigma_sigma_gamma,
    sigma_kappa = sigma_kappa,
    sigma_alpha0 = sigma_alpha0,
    sigma_alpha1 = sigma_alpha1,
    sigma_xi = sigma_xi,
    sigma_eta = sigma_eta,
    sigma_sigma_epsilon = sigma_sigma_epsilon,
    
    #asg or ulnorm
    beta = as.vector(beta),
    #beta = as_beta,
    n_params = n_params,
    m0 = m0,
    C0 = diag(C0),
    d_mat = diag(rep(1,n_params)),
    nu = nu,
    c = c,
    d = d,
    df = df,
    
    
    #sir only
    seasons = as.numeric(factor(unique(dat1$season))),
    seg_ind_start = dat_inds$for_ind,
    seg_ind_length = dat_inds$max_week,
    seg_ind_max = dat_inds$index,
    S0 = S0, t0 = 0,
    ts = dat1$season_week, #figure out how to get rid of this variable since same as weeks
    N = 1 #might go away after changing integrate_ode for ode)
    
    
    
    
    
    
    
    
    
  )
  
  info_list <- list(
    
    
    as_beta = as_beta, 
    theta_s = theta_s,
    alpha = alpha,
    all_years = seasons,
    year = season,
    region = state,
    forecast_date = ceiling_date(date(dat2$date[nrow(dat2)]), unit = 'week',
                                 week_start = 1)
    
  )
  
  return(list(data_list, info_list))
}



get_stan_data <- function(ILINet, both_flu, s_region = 'US', s_season = 2022,
                          ili_seasons = c(2003:2022), 
                          m_week = 'all', traj = 'asg') {
  m_season <- max(ili_seasons)
  if (m_week != 'all') {
    both_flu_samp <- both_flu %>% 
      filter(season == s_season, region == s_region) %>% 
      filter(season < m_season | season_week <= m_week)
    
    ILINet_samp <- ILINet %>% 
      filter(season %in% ili_seasons, region == s_region) %>% 
      filter(season < m_season | 
               (season == m_season & season_week <= m_week)) %>% 
      group_by(season) %>% 
      mutate(n_weeks = max(season_week))
  } else {
    ILINet_samp <- ILINet %>% 
      filter(season %in% ili_seasons, region == s_region) %>% 
      group_by(season) %>% 
      mutate(n_weeks = max(season_week))
    
    both_flu_samp <- both_flu %>% 
      filter(season == s_season, region == s_region)
  }
  
  stan_data <- make_stan_data(dat1 = ILINet_samp, dat2 = both_flu_samp, 
                              traj = traj, ili_seasons)
  
  return(stan_data)
  
}



get_quantile_forecasts <- function(pred_samp, forc_date, 
                                   point = 'mean',
                                   quantiles = c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99),
                                   targets = paste(1:4,'wk ahead inc flu hosp'),non_negative = TRUE,
                                   location = 'US') {
  
  forecasts <- apply(pred_samp, MARGIN=2,FUN = quantile, probs = quantiles)
  if (non_negative == TRUE) {forecasts[forecasts <= 0] <- 0}
  forecasts <- data.frame(forecasts)
  colnames(forecasts) <- targets
  forecasts$quantile <- rownames(forecasts)
  forecasts <- forecasts %>%
    mutate(quantile = as.character(as.numeric(strsplit(quantile,'%'))/100))
  
  forecasts <- forecasts %>%
    tidyr::pivot_longer(1:4,names_to = 'target') %>%
    mutate(type = 'quantile')
  
  
  mean_value <- apply(pred_samp,MARGIN=2,FUN=point)
  mean_value[mean_value < 0] <- 0
  
  targets <- paste(1:4,'wk ahead inc flu hosp')
  points <- data.frame(quantile = NA, target = targets,
                       value = mean_value, type = 'point')
  
  forecasts <- rbind(forecasts,points)
  
  forecasts <- forecasts %>%
    mutate(location = location, forecast_date = forc_date) %>%
    mutate(target_end_date=forecast_date +
             7*as.numeric(substr(target,1,1))-2) %>%
    select(forecast_date,target,target_end_date,
           location,type,quantile,value) %>%
    arrange(target,type,quantile)
  
  return(forecasts)
}
