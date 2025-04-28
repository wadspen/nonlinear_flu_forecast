library(cmdstanr)
library(lubridate)
library(ggplot2)
library(posterior)
library(bayesplot)
library(rstanarm)
library(stringr)
library(dplyr)
library(cdcfluview)
library(tidyr)

check_error <- function() {print("error is here")}
print("spencer is the coolest")
get_hosp_data <- function(lag = 1) {
  # raw_flu <- read.csv(
  #   'https://raw.githubusercontent.com/cdcepi/Flusight-forecast-data/master/data-truth/truth-Incident%20Hospitalizations.csv')
  
  last_date <- floor_date(today() - 12, unit = "week", week_start = 6)
  
  #raw_url <- paste("https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/main/auxiliary-data/target-data-archive/target-hospital-admissions_",
  #last_date, ".csv", sep = "")
  
  raw_url <- "https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/main/target-data/target-hospital-admissions.csv"
  
  # url <- "https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/main/auxiliary-data/target-data-archive/target-hospital-admissions_2023-09-23.csv"
  
  raw_flu <- read.csv(raw_url)
  raw_flu <- raw_flu %>% 
    mutate(year = year(date),week = week(date) - lag, 
           week2 = week(date) - (lag - 1)) %>% 
    mutate(week = week(date(date) - 1))
  
  return(raw_flu)
}

get_ili_data <- function(lag = 2) {
  ILINet_state <- ilinet(region = 'state')
  ILINet_us <- ilinet(region = 'national') %>% mutate(region = 'US')
  ILINet <- rbind(ILINet_state,ILINet_us)
  
  # ILINet <- ILINet %>% 
  #   mutate(year=year(week_start), week2=week(week_start),location_name=region)
  #this allows to control the lag. let's default to lag = 2
  ILINet <- ILINet %>%
	 mutate(week_start = date(week_start) + 7*(lag - 1))

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
    # group_by(region, year) %>%
    # mutate(new_week = ifelse(week == max(week) &
    #                            year < year(today()), 1, week + 1)) %>%
    # mutate(new_year = ifelse(week == max(week) &
    #                            year < year(today()),
    #                          year + 1, year)) %>%
    # ungroup() %>%
    # mutate(year = new_year, week = new_week) %>% 
    mutate(season = ifelse(month(week_start) <= 7, year - 1, year),
           place = 1) %>% 
    group_by(season, region) %>% 
    mutate(season_week = ifelse(season == 2010 & region != "US", 
                                cumsum(place) + 8, cumsum(place)),
           unweighted_ili = ifelse(unweighted_ili == 0, .00001,
                                   unweighted_ili)) %>% 
    mutate(unweighted_ili = ifelse(region == "District of Columbia" &
                                     season == 2021 & is.na(unweighted_ili),
                                   2.94468,
                                   unweighted_ili)) %>% 
    filter(is.na(unweighted_ili) == FALSE, season < 2024)
  # View(ILINet[ILINet$region == "Minnesota",c("region", "season", "season_week")])
  #max_season <- max(ILINet$season)
  #ILINet <- ILINet %>% group_by(region, season) %>%
    #mutate(new_week = ifelse(season_week == max(season_week) &
    #                           season < year(today()),
			       #season < max_season, 
    #                         1, season_week + 1)) %>%
    #mutate(new_season = ifelse(season_week == max(season_week) &
    #                             season < year(today()),
			         #season < max_season,
    #                           season + 1, season)) %>%
    #ungroup() %>%
    #mutate(season = new_season, season_week = new_week)
  
  
  pops <- read.csv("https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/main/auxiliary-data/locations.csv")
  # colnames(pops)[1] <- "region"
  pops <- pops %>% 
    select(location_name, population, count_rate1, 
           count_rate2,
           count_rate3,
           count_rate4, count_rate5)
  colnames(pops)[1] <- "region"
  
  ILINet <- ILINet %>%
    left_join(pops, by = "region")
  
  return(ILINet)
}

comb_data <- function(lag = 2, s_region = NA, s_season = 2022) {
  raw_flu <- get_hosp_data(lag = lag) %>% 
    mutate(season = ifelse(month(date) <= 7 | month(date) == 8 & mday(date) < 7, year(date) - 1, year(date)),
           place = 1) %>% 
    arrange(desc(season), date) %>% 
    group_by(season, location) %>% 
    mutate(season_week = cumsum(place)) %>%
    select(-place) 
  
  ILINet <- get_ili_data(lag = lag) %>% 
    select(region, unweighted_ili, week_start, location_name, season, place,
           season_week, population, count_rate1, 
           #count_rate2, 
           count_rate3,
           count_rate4, count_rate5) %>%
    select(-place)
  
  if (is.na(s_region)) {
    both_flu <- raw_flu %>% 
      left_join(ILINet, by = c('season','season_week','location_name'))
  }
  
  else {
    both_flu <- raw_flu %>% 
      left_join(ILINet, by = c('season','season_week','location_name')) %>% 
      filter(region == s_region)
  }
  
  # pops <- read.csv("./state_populations.csv")
  # pops <- read.csv("https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/main/auxiliary-data/locations.csv")
  # # colnames(pops)[1] <- "region"
  # pops <- pops %>% 
  #   select(location_name, population, count_rate1, count_rate2, count_rate3,
  #          count_rate4, count_rate5)
  # colnames(pops)[1] <- "region"
  both_flu <- both_flu %>%
    # left_join(pops, by = "region") %>%
    mutate(unweighted_ili = ifelse(unweighted_ili == 0, .00001,
                                   unweighted_ili)) %>% 
    filter(is.na(unweighted_ili) == FALSE, is.na(value) == FALSE) %>%
    arrange(season, season_week)
  
  # print(max(both_flu$date))
  
  return(both_flu)
  
}

make_stan_data <- function(dat1, dat2, ili_seasons,
                           nu=1, df = 4, c = 4, d = 4, S0 = .9, 
                           sigma_sigma_gamma_w = 2, sigma_sigma_gamma = .02,
                           sigma_kappa = 10000, sigma_alpha0 = 1000,
                           sigma_alpha1 = 1000, sigma_alpha3 = .5, 
			   sigma_xi = 3, sigma_eta = .4,
                           sigma_sigma_epsilon = 20,
                           rho_mu = .68, rho_sigma = .08,
                           beta_mu = .8, beta_sigma = .3,
                           I0_mu = .005, I0_sigma = .03,
                           sigma_mu_alpha_w = 2, sigma_mu_beta_w = 2,
                           sigma_mu_alpha = 2, sigma_mu_beta = .02,
			   sigma_gamma_W = 1, 
                           traj = 'asg', dist = "norm", hosp_log = FALSE) {
  
  start <- c(-4.2533962, 0.3416948, 22.5023919, 2.3580754, 
	     2.3519650, 100.3432742)
  
  state <- unique(dat1$region)
  seasons <- unique(dat1$season)
  season <- unique(dat2$season)
  #alpha <- coef(lm(dat2$value/dat2$pop_2022 ~ dat2$unweighted_ili))
  #alpha <- coef(lm(log(dat2$value + 1) ~ dat2$unweighted_ili^2))
  if (traj == 'asg') {
    
    formle <- dat1 %>%
	    ungroup() %>%
	    filter(season == max(season))
    print(dim(formle))
    pars <- asg_max_lik(formle$unweighted_ili/100, formle$season_week, 
                        theta = start)
    m0 <- c(.3, 23, 3.69, 4.7)
    C0 <- c(.2, 5, 2, 2)
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
    theta_s <- rbind(theta_s, pars[-1])
    theta_s[length(seasons),] <- apply(as.matrix(theta_s[1:(length(seasons) - 1),]), MARGIN = 2, FUN = mean)
  }
  
  else if (traj == 'uln') {
    dat2f <- dat2 %>% 
      mutate(unweighted_ili = ifelse(unweighted_ili <= .00001,
                                     .01, unweighted_ili))
    pars <- uln_max_lik(dat2f$unweighted_ili/100, dat2f$season_week, 
                        theta = start)
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
    
    theta_s[length(seasons),] <- apply(as.matrix(theta_s[1:(length(seasons) - 1),]), MARGIN = 2, FUN = mean)
  }
  
  #use the following two lines to calculate the estimated sd of the discrepancy at week W
  #asg is calculated at the MLE values for each season, and the sd is the sample sd of the
  #difference between the asg and the logit of ILI  
  dat1 <- dat1 %>%
	 left_join(all_season_region_mles, by = c("region", "season")) %>%
	 group_by(region, season) %>%
	 mutate(asg = asg(season_week, beta, eta, mu, sigma1, sigma2)) %>%
	 mutate(disc = boot::logit(unweighted_ili/100) - asg)
  print("I'm getting very frustrated...")
 
  dat1 <- dat1 %>%
	 ungroup() %>%
	 group_by(region, season) %>%
	 mutate(sir = get_sir_mle(season_week = season_week,
				  unweighted_ili = unweighted_ili/100)) %>%
	 mutate(sir_disc = boot::logit(sir) - boot::logit(unweighted_ili/100))
	  
  print("what about here?")
  disc_dat <- dat1 %>%
	  group_by(region, season) %>%
	  filter(season_week == max(season_week)) %>%
	  group_by(region) %>%
	  summarise(
		    asg_disc_sd = var(disc, na.rm = TRUE),
		    sir_disc_sd = var(sir_disc, na.rm = TRUE)
          )
  week_inds <- dat1 %>% 
    group_by(season) %>% 
    summarise(
      #max_week = max(season_week),
      max_week = n()
    )
  is_last <- function(x) !duplicated(x, fromLast = TRUE)
  week_inds$index <- which(is_last(dat1$season))
  dat_inds <- week_inds %>% 
    mutate(for_ind = index - max_week + 1) %>%
    mutate(for_ind = ifelse(for_ind < 1, 1, for_ind))
  print(dat_inds$for_ind)
  
  
  ts <- c()
  for (i in 1:length(dat_inds$max_week)) {
    temp <- 1:dat_inds$max_week[i]
    ts <- c(ts, temp)
  }
  # dat_inds$max_week[nrow(dat_inds)] <- last(dat2$season_week)
  
  #population <- read.csv("https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/main/auxiliary-data/locations.csv")
  #count_rate = population$count_rate2[population$location_name == "Minnesota"] 
  count_rate <- unique(dat1$count_rate2)

  if (hosp_log == TRUE) {hosp <- log(dat2$value + 1)}
  else {hosp <- dat2$value} 
  print("does it get here?") 
  data_list <- list(
    #for all models
    M = nrow(dat1),
    HM = nrow(dat2),
    n_seasons = length(unique(dat1$season)),
    n_weeks = max(ILINet$season_week),
    cur_yr_n_weeks = last(dat2$season_week),
    all_seasons = as.numeric(factor(dat1$season)),
    weeks = dat1$season_week,
    ps = max(as.numeric(factor(unique(dat1$season)))),
    hps = max(as.numeric(factor(unique(dat2$season)))),
    count_rate = count_rate,
    ili = dat1$unweighted_ili/100,
    # ili = boot::logit(dat1$unweighted_ili/100),
    # ili = dat2$unweighted_ili,
    #hosp = log(dat2$value + 1),
    hosp = hosp,
    n_seasons_hosp = length(unique(dat2$season)),
    #hosp = dat2$value,
    hosp_seasons = as.numeric(factor(dat2$season)),
    hosp_season_weeks = dat2$season_week,
    sigma_sigma_gamma_w = sigma_sigma_gamma_w,
    sigma_sigma_gamma = sigma_sigma_gamma,
    sigma_kappa = sigma_kappa,
    sigma_alpha0 = sigma_alpha0,
    sigma_alpha1 = sigma_alpha1,
    sigma_alpha3 = sigma_alpha3,
    sigma_xi = sigma_xi,
    sigma_eta = sigma_eta,
    sigma_sigma_epsilon = sigma_sigma_epsilon,
    sigma_disc = disc_dat$asg_disc_sd,
    sigma_disc_sir = disc_dat$sir_disc_sd,
    sigma_gamma_W = sigma_gamma_W,
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
    #ts = c(dat1$season_week, last(dat1$season_week) + 1), #figure out how to get rid of this variable since same as weeks
    ts = dat1$season_week, #the above line is for when we need to deal with the lag (on submission days since ILI is not updated yet)
    N = 1, #might go away after changing integrate_ode for ode)
    
      #prior hyperparameters
    rho_mu = rho_mu,
    rho_sigma = rho_sigma,
    beta_mu = beta_mu,
    beta_sigma = beta_sigma,
    I0_mu = I0_mu,
    I0_sigma = I0_sigma,
    sigma_mu_alpha_w = sigma_mu_alpha_w,
    sigma_mu_beta_w = sigma_mu_beta_w,
    sigma_mu_alpha = sigma_mu_alpha,
    sigma_mu_beta = sigma_mu_beta
    
    
    
    
    
    
    
    
    
  )
  
  info_list <- list(
    
    
    # as_beta = as_beta, 
    theta_s = theta_s,
    alpha = alpha,
    all_years = seasons,
    year = season,
    region = state,
    # forecast_date = ceiling_date(date(dat2$date[nrow(dat2)]), unit = 'week',
    #                              week_start = 1)
    
    forecast_date = date(dat2$date[nrow(dat2)]),
    last_count = dat2$value[nrow(dat2)],
    location = unique(dat2$location)
    
  )
  
  return(list(data_list, info_list))
}



get_stan_data <- function(ILINet, both_flu, s_region = 'US', s_seasons = 
                            2022:2023,
                          ili_seasons = c(2003:2019, 2021:2023),
                          nu=1, df = 4, c = 4, d = 4, S0 = .9, 
                          sigma_sigma_gamma_w = 2, sigma_sigma_gamma = .02,
                          sigma_kappa = 10000, sigma_alpha0 = 1000,
                          sigma_alpha1 = 1000, sigma_alpha3 = .5, 
			  sigma_xi = 3, sigma_eta = .4,
                          sigma_sigma_epsilon = 20,
			  sigma_gamma_W = 1,
                          rho_mu = .68, rho_sigma = .08,
                          beta_mu = .8, beta_sigma = .3,
                          I0_mu = .005, I0_sigma = .03,
                          m_week = 'all', traj = 'asg', dist = "norm",
			  hosp_log = FALSE) {
  m_season <- max(ili_seasons)
  ILINet <- ILINet %>%
	  filter(season != 2020)

  if (m_week != 'all') {
    both_flu_samp <- both_flu %>% 
      filter(season %in% s_seasons, region == s_region) %>% 
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
      filter(season %in% s_seasons, region == s_region) 
  }
  #print(both_flu_samp[,c("season_week", "value")])
  stan_data <- make_stan_data(dat1 = ILINet_samp, dat2 = both_flu_samp, 
                              traj = traj, ili_seasons,
                              sigma_sigma_gamma_w = sigma_sigma_gamma_w,
                              sigma_sigma_gamma = sigma_sigma_gamma, dist = dist,
  			      hosp_log = hosp_log)
  
  return(stan_data)
  
}


get_mix_params <- function(x, max.it = 30000, comp = 4) {
	fit <- mixtools::normalmixEM(x, k = comp, maxit = max.it, maxrestarts = 100)
	mix_dist <- data.frame(family = "Norm",
		       	       param1 = fit$mu,
			       param2 = fit$sigma,
			       weight = fit$lambda)

	return(mix_dist)	
}
get_quantile_forecasts <- function(pred_samp, forc_date, 
                                   point = 'mean',
                                   quantiles = c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99),
                                   targets = paste(1:5,'wk ahead inc flu hosp'),
                                   horizons = -1:3,
                                   non_negative = TRUE,
                                   location = "US",
                                   location_name = "US",
                                   m_week, last_count,
				   reference_date, dist = "norm", hosp_log = FALSE) {
  print("get quantile forecasts function")
  print(paste("m_week =", m_week))
  print(dim(pred_samp))
  pred_samp <- pred_samp[,(m_week):(m_week + 4)]
  if (hosp_log == TRUE) {pred_samp <- exp(pred_samp) - 1}
  forecasts <- apply(pred_samp, MARGIN=2,FUN = quantile, probs = quantiles)
  #if (non_negative == TRUE) {forecasts[forecasts <= 0] <- 0}
  forecasts <- data.frame(forecasts)
  colnames(forecasts) <- horizons
  forecasts$output_type_id <- rownames(forecasts)
  forecasts <- forecasts %>%
    mutate(output_type_id = 
             as.character(as.numeric(strsplit(output_type_id,'%'))/100))
  
  forecasts <- forecasts %>%
    tidyr::pivot_longer(1:5,names_to = 'horizon') %>%
    mutate(output_type = 'quantile')
  
  
  # mean_value <- apply(pred_samp,MARGIN=2,FUN=point)
  # mean_value[mean_value < 0] <- 0
  
  # targets <- paste(1:5,'wk ahead inc flu hosp')
  # points <- data.frame(quantile = NA, target = targets,
  #                      value = mean_value, type = 'point')
  
  # forecasts <- rbind(forecasts,points)
  
  # reference_date <- ceiling_date(today(), unit = "week", week_start = 6)
  forecasts <- forecasts %>%
    mutate(reference_date = reference_date,
           location = location, 
           target = "wk inc flu hosp") %>%
    mutate(target_end_date = reference_date + 7*horizons) %>%
    select(reference_date, horizon, target, target_end_date, 
           location, output_type, output_type_id, value) %>%
    arrange(horizon, output_type, output_type_id)
  
  loc <- location
  pops <- read.csv("https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/main/auxiliary-data/locations.csv")
  pop <- pops %>% 
    filter(location == loc)
  
  
  colnames(pred_samp) <- horizons
  
  long_preds <- pred_samp %>%
    tidyr::pivot_longer(1:5,names_to = 'horizon') %>%
    mutate(output_type = 'quantile')
  
  # long_preds <- cbind(long_preds, pop)
  
  pop_100k <- pop$population/100000
  yt_27 <- last_count/pop_100k
  
  long_preds_pmf <- long_preds %>% 
    mutate(stable_cut = 1/100000, cut = ifelse(as.numeric(horizon) == 3, 
                             5/100000, (as.numeric(horizon) + 3)/100000), 
           yt_h7 = value/pop_100k, yt_27 = yt_27) %>% 
    mutate(count_change = value - last_count, rate_change = yt_h7 - yt_27) %>% 
    mutate(rate_bin = ifelse(count_change < 10 | abs(rate_change) < stable_cut,
                             "stable", 
                             ifelse(abs(rate_change) < cut & rate_change < 0,
                                    "decrease",
                                    ifelse(abs(rate_change) < cut & rate_change > 0,
                                           "increase",
                                           ifelse(abs(rate_change) > cut & rate_change < 0,
                                                  "large_decrease",
                                                  ifelse(abs(rate_change) > cut & rate_change > 0,
                                                         "large_increase", NA))))))
  
  bin_preds <- long_preds_pmf %>% 
    group_by(horizon) %>% 
    reframe(
      output_type_id = names(table(rate_bin)),
      value = table(rate_bin)/n()
    ) %>% 
    mutate(reference_date = reference_date, target = "wk flu hosp rate change",
           target_end_date = reference_date + as.numeric(horizon)*7, 
           location = location,
           output_type = "pmf") %>% 
    select(reference_date, horizon, target, target_end_date, location, 
           output_type, output_type_id, value)
  
  
  pmfs <- expand.grid(horizon = as.character(-1:3), 
                      target = "wk flu hosp rate change", 
                      location = unique(forecasts$location), output_type = "pmf", 
                      output_type_id = c("large_increase", "increase", "stable", 
                                         "decrease", "large_decrease"), 
                      reference_date = reference_date) %>% 
    mutate(target_end_date = reference_date + 7*horizons)
  
  bin_preds <- pmfs %>% 
    left_join(bin_preds, by = c("horizon", "target", "location", 
                                "output_type", "output_type_id",
                                "reference_date", "target_end_date")) %>% 
    mutate(value = ifelse(is.na(value), 0, value)) %>% 
    group_by(horizon, location) %>% 
    mutate(value = value/sum(value)) %>% 
    select(reference_date, horizon, target, target_end_date, location, 
           output_type, output_type_id, value)
   #begin trouble mixture distribution thing
   # mix_dists <- long_preds %>%
   #	    group_by(horizon) %>%
   #	    reframe(
   #	    	get_mix_params(value)
   #	    )

   # dists <- expand.grid(horizon = as.character(-1:3),
   #			 target = "wk inc flu hosp",
   #			 location = unique(forecasts$location),
   #			 output_type = "mix_dist",
   #			 output_type_id = NA,
   #			 reference_date = reference_date) %>%
   # 		mutate(target_end_date = reference_date + 7*horizons,
   #		       value = NA)

   # mix_dists <- mix_dists %>%
   #	    left_join(dists, by = "horizon") %>%
   #	    select(reference_date, horizon, target, target_end_date, location,
   #		    output_type, output_type_id, value, family, param1,
   #		    param2, weight)
   #end of commented mixture distribution
	    
  # long_pred <- long_preds %>% 
    # mutate(output_type_id = ifelse(between(value, 0, count_rate1), "large decrease",
    #                                ifelse(between(value, count_rate1, count_rate2), "decrease",
    #                                       ifelse(between(value, count_rate2, count_rate3), "stable",
    #                                              ifelse(between(value, count_rate3, count_rate4)))))))
  
  forecasts <- rbind(forecasts, bin_preds) %>% 
    arrange(output_type, location, horizon, output_type_id) %>%
    mutate(family = NA, param1 = NA, param2 = NA, weight = NA)
  
  #skip the mixture distribution for now
  #forecasts <- rbind(forecasts, mix_dists)
  return(forecasts)
  #this is the file I want
}
