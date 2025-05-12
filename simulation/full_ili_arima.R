setwd("../")
#source("./simulation/simulate_hospitalizations.R")
source("./flu_forecast_23/get_data_functions_new.R")
source("./flu_forecast_23/mle_functions.R")
library(parallel)
library(doParallel)
library(doMC)
library(evalcast)
library(forecast)
#n.cores <- detectCores()
#n.cores <- 1
n.cores <- 1
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

select_regions = "US"
ILINet_state <- ILINet %>%
   filter(region == select_regions)
print(unique(ILINet$region))
season_levels <- unique(ILINet$season)
season_levels <- season_levels[season_levels != 2023]

#season_levels <- 2022
N <- 50000
foreach(j = season_levels,
	.packages = c("tidyr", "dplyr", "evalcast")
	#,.errorhandling = "remove"
	) %:%
   foreach(k = 9:38) %dopar% { #c(14,20,26,32,38)) %dopar% {
#for (i in 1:(length(season_levels))) {
 #j <- 2013; k <- 17
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
      
      
      arima <- auto.arima(ILINet_filter$unweighted_ili)
      arimameans <- forecast(arima)$mean[1:4]
      low <- as.vector(forecast(arima)$lower[1:4,2])
      arimasds <- (arimameans - low)/qnorm(.975)


      arimaforcs <- apply(cbind(arimameans, arimasds),
			  MARGIN = 1, FUN = function(x) {
				  rnorm(N, 
					mean = x[1],
					sd = x[2])})

      
	#print(baseforcs)
      save_name <- paste("./simulation/ili_fits/", "arima", 
			 "/", j, "/", "week", 
	    		 k,  
			 ".csv", sep = "")

      write.csv(arimaforcs, save_name, row.names = FALSE)
       
    
    }
  
#}


















