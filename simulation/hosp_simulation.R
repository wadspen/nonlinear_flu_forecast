setwd("../")
source("./simulation/simulate_hospitalizations.R")
source("./get_data_functions_new.R") #this is a gerry-rig file. It comments out some functions that NOVA was getting stuck on. I don't think other machines will have the same problem since my machine didn't seem to have it. It gets stuck in the mle2 function in mle_functions.R and is related to the parallelzation. The issue does not occur when not paralelising. Comment out and use the line below instead on a different machine.
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


state <- "US"
both_flu_state <- both_flu %>% 
  filter(region == state)
ILINet_state <- ILINet %>%
	filter(region == state)

hosp_log <- TRUE
model <- "hosp_log_ar1"

reps <- 2000
set.seed(96)
seeds <- sample(99999999, reps)



rep <- 1

print("simulation starts now")
foreach(rep = 1000:reps,
	.packages = c("tidyr", "dplyr", "evalcast")
	,.errorhandling = "remove"
	) %dopar% {
   	set.seed(seeds[rep])
	ILINet_sim <- ILINet %>% 
  	filter(season != 2020, region == state) %>%
        #filter(season == 2010) %>%	
  	group_by(season, region) %>% 
  	mutate(ili = unweighted_ili) %>% 
  	mutate(hosp_sim = simulate_hospitalizations(ili, both_flu = both_flu_state, 
						      is.log = hosp_log, rep = rep)) %>%
  	ungroup()
  
  

       
    
    save_name_sim <- paste("./simulation/hosp_sim/", model,
			   "/sim_hosp/sims_rep", rep, "_seed",
			   seeds[rep], ".csv", sep = "")
    
    print(rep) 
    sims <- ILINet_sim %>%
	    select(region, season, season_week, unweighted_ili, 
		   ili, count_rate2, hosp_sim)

    write.csv(sims, save_name_sim, row.names = FALSE)
    
    rep <- rep + 1
    if (rep >= reps) {break}   
  
}




