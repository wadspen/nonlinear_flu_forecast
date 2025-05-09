setwd("../flu_forecast_23")
#source("./simulation/simulate_hospitalizations.R")
source("./get_data_functions_new.R")
source("./mle_functions.R")
library(parallel)
library(doParallel)
library(doMC)
library(evalcast)
#n.cores <- detectCores()
n.cores <- 60
my.cluster <- makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
registerDoMC(cores = n.cores)

both_flu <- comb_data() %>%
	filter(season >= 2010)
ILINet <- get_ili_data() %>%
	filter(season >= 2010, !(season %in% c(2020, 2023)))
args <- commandArgs()
model <- args[6]#; model <- "sir"
print(model)
mod_string <- paste("../stan_models/", model, ".stan", sep = "")
mod <- cmdstan_model(stan_file = mod_string)
setwd("./flu_forecast_23")
select_regions = "US"
ILINet_state <- ILINet %>%
   filter(region == select_regions)

season_levels <- unique(ILINet$season)
season_levels <- season_levels[season_levels != 2023]

season_levels <- 2022
foreach(j = season_levels,
	.packages = c("tidyr", "dplyr", "evalcast")
	#,.errorhandling = "remove"
	) %:%
   foreach(k = 9:38) %dopar% { #c(14,20,26,32,38)) %dopar% {
#for (i in 1:(length(season_levels))) {
  ILINet_arrange <- ILINet_state %>%
	  filter(season != 2023) %>%
	  mutate(season = ifelse(season == j, 2023,
				 season)) %>%
  	  arrange(season)
  
    #ILINet_arrange_state <- ILINet_arrange %>% 
    #  filter(region == select_regions)
    
    
    #for (k in seq(14, 34, 4)) {

      ILINet_filter <- ILINet_arrange %>%
	      filter(season < 2023 | (season == 2023 & season_week <= k))
      #print(head(ILINet_filter))
      week <- max(ILINet_filter$week_start[ILINet_filter$season == 2023])
      dat <- get_stan_data(ILINet_filter, 
                           both_flu, s_region = select_regions,
                           ili_seasons = unique(ILINet_filter$season),
      			   m_week = k)
      
      stan_dat <- dat[[1]]; print(stan_dat$seg_ind_start[stan_dat$ps - 1]); 
      
      info_list <- dat[[2]]
      forecast_date <- info_list$forecast_date
      init <- list(theta = stan_dat$m0, theta_s = info_list$theta_s)
      
      samps <- mod$sample(data = stan_dat,
                          chains = 1,
                          adapt_delta = .9999,
                          init = list(init),
                          iter_warmup = 10000, 
                          iter_sampling = 50000)
 
      draws <- samps$draws(format = "df")
      pred_ili_draws <- draws[,str_detect(colnames(draws), "pred_ili.")]
      pred_draws <- pred_ili_draws[,(ncol(pred_ili_draws) - 5):ncol(pred_ili_draws)]
      save_name <- paste("./simulation/ili_fits/", model, 
			 "/", j, "/", "week", 
	    		 k,  
			 ".csv", sep = "")

      write.csv(pred_draws, save_name, row.names = FALSE)
       
    
    }
  
#}


















