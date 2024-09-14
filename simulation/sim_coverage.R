library(ggplot2)
library(dplyr)
library(tidyverse)
library(stringr)
library(readr)
library(evalcast)
library(scoringRules)
library(parallel)
library(doParallel)
library(doMC)
n.cores <- detectCores()
my.cluster <- makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
registerDoMC(cores = n.cores)


args <- commandArgs()
mod <- args[6]

coverage <- function(probs, values, true_value, alpha = .05) {
  probs <- round(probs, 3)
  upper <- round(1 - alpha/2, 3)
  lower <- round(alpha/2, 3)
  upper_ind <- which(probs == upper)
  lower_ind <- which(probs == lower)
  cover <- between(true_value, values[lower_ind], values[upper_ind])
  return(cover)
}

probs <- round(c(0, .0025, .005, .01, .025, seq(.05, .95, by = .05), .975, .99, .995, .9975, 1), 3)

hosp_dir <- "./hosp_sim/hosp_log_ar1/" 
sim_data_dir <- paste0(hosp_dir, "sim_hosp/")
sim_forecast_dir <- paste0(hosp_dir, "new/forecasts/", mod, "/")

all_results <- data.frame()

reps <- 1100
for (rep in 1:reps) {

#all_results <- foreach(rep = 1:150,
#		       .packages = c("dplyr", "tidyr", "stringr", "readr", "evalcast",
#				     "scoringRules")
#		       ,.errorhandling = "remove"
#		       ,.combine = rbind
#		       ) %:%
	#for (wk in seq(14, 38, 6)) {
	results <- foreach(seas = c(2010:2019, 2021:2022),
				.packages = c("dplyr", "tidyr", "stringr", "readr", "evalcast",
					      "scoringRules")
					      ,.errorhandling = "remove"
					      ,.combine = rbind
					      ) %:%
				foreach(wk = seq(14, 38, 6), .combine = rbind) %dopar% {

		     sim_file <- list.files(sim_data_dir, pattern = paste0("rep", rep, "_"))
		     sim_file_path <- paste0(sim_data_dir, sim_file)
	     	     sim_data <- read.csv(sim_file_path)
		     sim_data <- sim_data %>%
			     mutate(hosp_sim = log(hosp_sim + 1))
	#for (seas in c(2010:2019, 2021)) {
	#	for (wk in seq(14, 38, 6)) {
		
			forecast_file <- paste0(sim_forecast_dir, seas, 
						"/week", wk, "/rep", rep, ".csv")
			#if (file.exists(forecast_file) == FALSE) {next}
			sim_forecast <- read.csv(forecast_file) #%>% #select(-X)
		 	
			
                        quantiles <- apply(sim_forecast, MARGIN = 2, quantile, probs = probs) %>% data.frame()
                        quantiles$probs <- probs
                        wis_cov <- quantiles %>% 
			  pivot_longer(1:4, names_to = "week_ahead", values_to = "pred_value") %>% 
			  mutate(week_ahead = parse_number(week_ahead)) %>% 
			  arrange(week_ahead) %>% 
			  group_by(week_ahead) %>% 
			  mutate(season_week = week_ahead + wk + 1, season = seas) %>% 
			  left_join(sim_data %>% 
				      select(season, season_week, hosp_sim), by = c("season", "season_week")) %>% 
			  group_by(season, week_ahead) %>% 
			  reframe(
			    cover50 = coverage(probs, pred_value, unique(hosp_sim), alpha = .5),
			    cover60 = coverage(probs, pred_value, unique(hosp_sim), alpha = .4),
			    cover70 = coverage(probs, pred_value, unique(hosp_sim), alpha = .3),
			    cover80 = coverage(probs, pred_value, unique(hosp_sim), alpha = .2),
			    cover90 = coverage(probs, pred_value, unique(hosp_sim), alpha = .1), 
			    cover95 = coverage(probs, pred_value, unique(hosp_sim), alpha = .05),
			    cover98 = coverage(probs, pred_value, unique(hosp_sim), alpha = .02),
			    cover99 = coverage(probs, pred_value, unique(hosp_sim), alpha = .01),
			    cover99.5 = coverage(probs, pred_value, unique(hosp_sim), alpha = .005),
			    cover100 = coverage(probs, pred_value, unique(hosp_sim), alpha = 0)
			  )
		
                        results <- wis_cov %>% 
			  mutate(last_week = wk, rep = rep, model = mod) %>% 
			  select(model, season, last_week, week_ahead, rep, cover50, cover60, 
				 cover70, cover80, cover90, cover95,
				 cover95, cover98, cover99, cover99.5, cover100) %>% data.frame()
		  	#all_results <- rbind(all_results, results)
			results
                     
		}
        print(dim(results))
	print(class(results))	
	write.csv(results, "test.csv")
	all_results <- rbind(all_results, results)
	print(rep)
}
			

write.csv(all_results, paste0("./hosp_sim/cover_", mod, ".csv"), row.names = FALSE)













