setwd("../")
source("./simulation/simulate_hospitalizations.R")
source("./get_data_functions_new_hosp.R") #this is a gerry-rig file. It comments out some functions that NOVA was getting stuck on. I don't think other machines will have the same problem since my machine didn't seem to have it. It gets stuck in the mle2 function in mle_functions.R and is related to the parallelzation. The issue does not occur when not paralelising. Comment out and use the line below instead on a different machine.
#source("./get_data_function_new.R")
source("./mle_functions.R")
library(parallel)
library(doParallel)
library(doMC)
library(evalcast)
n.cores <- detectCores()
my.cluster <- makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
registerDoMC(cores = n.cores)

both_flu <- comb_data(lag = 2) %>%
    filter(season >= 2010) %>% select(-X)
ILINet <- get_ili_data() %>%
    filter(season >= 2010)

args <- commandArgs()
model <- args[6]
mod_string <- paste("./", model, ".stan", sep = "")
mod <- cmdstan_model(stan_file = mod_string)
state <- "US"
both_flu_state <- both_flu %>% 
  filter(region == state)
ILINet_state <- ILINet %>%
	filter(region == state)

hosp_log <- ifelse(str_detect(model, "log"), TRUE, FALSE)

reps <- 1000
set.seed(96)
seeds <- sample(99999999, reps)

m <- 1
repeat {
   set.seed(seeds[m])
ILINet_sim <- ILINet %>% 
  filter(season != 2020, region == state) %>% 
  group_by(season, region) %>% 
  mutate(ili = count_rate2*unweighted_ili/100) %>% 
  mutate(hosp_sim = simulate_hospitalizations(ili, both_flu = both_flu_state, 
					      is.log = hosp_log)) %>%
  ungroup()
  
 

season_levels <- unique(ILINet_sim$season)
season_levels <- season_levels[which(season_levels != 2023)]

#for (j in season_levels) {
#	for (k in seq(14, 38, 5)) {
foreach(j = season_levels,
	.packages = c("tidyr", "dplyr", "evalcast")
	,.errorhandling = "remove"
	#) %:%
	) %:% 
	#k <- 14
  foreach(k = seq(14, 38, 6)) %dopar% {
  
  ILINet_arrange <- ILINet_state %>%
	  filter(season  != 2023) %>%
	  mutate(season = ifelse(season == j, 2023,
				 season)) %>%
	  arrange(season) %>%
	  filter(season < 2023 | (season == 2023 & season_week <= k))

  sims_arrange <- ILINet_sim %>%
    filter(season != 2023) %>% 
    mutate(season = ifelse(season == j, 2023, season)) %>%
    #filter(season < 2023 | season_week <= k) %>% 
    arrange(season)
    #print(unique(sims_arrange$season))
  
  both_flu_for_sim <- both_flu_state %>% 
    ungroup() %>% 
    select(region, unweighted_ili, season, season_week, value) %>% 
    filter(season != 2023)
  
  sims_to_forecast <- sims_arrange %>% 
    filter(season == 2023) %>% 
    mutate(value = hosp_sim) %>% 
    select(region, unweighted_ili, season, season_week, value)
  
  sims_both_flu <- rbind(both_flu_for_sim, sims_to_forecast) %>% 
    arrange(season) %>% 
    mutate(sim_value = value)
  
 
    #hosp_log <- ifelse(str_detect(model, "log"), TRUE, FALSE)
    ILINet_sim_hold <- sims_arrange %>% 
      filter(season < 2023 | season_week <= k)

       
    sims_both_flu_hold <- sims_both_flu %>% 
      filter(season < 2023 | season_week <= k)
    print("we getting here?") 
    dat <- get_stan_data(ILINet_arrange, sims_both_flu_hold,
                         s_region = state,
                         ili_seasons = unique(ILINet_arrange$season),
    			 hosp_log = hosp_log,
    			 m_week = k)
    print("please make it!")
    stan_dat <- dat[[1]]
    info_list <- dat[[2]]
    
    init <- list(theta = stan_dat$m0, theta_s = info_list$theta_s)
    print("does it get here?")
    samps <- mod$sample(data = stan_dat,
			chains = 1,
			adapt_delta = .9999,
			iter_warmup = 8000,
			iter_sampling = 5000)

    draws <- samps$draws(format = "df")
    
    save_name <- paste("./simulation/hosp_sim/", model,
		       "/", j, "/", "week", k, "/",
		       "rep", m, ".csv", sep = "")
	
    write.csv(draws, save_name, row.names = FALSE) 
}
    save_name_sim <- paste("./simulation/hosp_sim/", model,
			   "/sim_hosp/sims_rep", m, "_seed",
			   seeds[m], ".csv", sep = "")
    
    sims <- ILINet_sim %>%
	    select(region, season, season_week, unweighted_ili, 
		   ili, count_rate2, hosp_sim)

    write.csv(sims, save_name_sim, row.names = FALSE)
    
    m <- m + 1
    if (m >= reps) {break}   
  
}




