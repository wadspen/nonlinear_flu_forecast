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

probs <- round(c(.01, .025, seq(.05, .95, by = .05), .975, .99), 3)

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
		 	cont_scores <- data.frame()
			for (i in 1:4) {
			  sim_value <- sim_data %>% 
			    filter(season == seas, season_week == wk + i) %>%
			    select(hosp_sim) %>% as.vector() %>% as.numeric()
		          kdens <- density(sim_forecast[,i], kernel = "gaussian", from = -12, to = 12)
			  logs2 <- -log(approx(kdens$x, kdens$y, sim_value)$y)
			  crps <- crps_sample(sim_value, as.vector(sim_forecast[,i]))
			  logs <- logs_sample(sim_value, as.vector(sim_forecast[,i]))
			  pit <- ecdf(as.vector(sim_forecast[,i]))(sim_value)
			  scores <- data.frame(season = seas, week_ahead = i - 1, 
					       crps = crps, logs = logs, logs2 = logs2, pit = pit)
			  cont_scores <- rbind(cont_scores, scores)
			}
			
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
			    wis = weighted_interval_score(probs,
							  pred_value, 
							  hosp_sim),
		            lmae = abs(median(pred_value - unique(hosp_sim))),		  
			    mae = abs(median(exp(pred_value) - unique(exp(hosp_sim)))), 
			    cover50 = coverage(probs, pred_value, unique(hosp_sim), alpha = .5), 
			    cover95 = coverage(probs, pred_value, unique(hosp_sim), alpha = .05)
			  )
		
                        results <- wis_cov %>% 
			  left_join(cont_scores, by = c("season", "week_ahead")) %>% 
			  mutate(last_week = wk, rep = rep, model = mod) %>% 
			  select(model, season, last_week, rep, week_ahead, wis, crps, 
				 logs, logs2, lmae, mae, pit, cover50, cover95) %>% data.frame()
		  	#all_results <- rbind(all_results, results)
			results
                     
		}
        print(dim(results))
	print(class(results))	
	write.csv(results, "test.csv")
	all_results <- rbind(all_results, results)
	print(rep)
}
			

write.csv(all_results, paste0("./hosp_sim/scores_", mod, ".csv"), row.names = FALSE)













