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
source("../get_data_functions_new.R")
n.cores <- detectCores()
my.cluster <- makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
registerDoMC(cores = n.cores)


args <- commandArgs()
mod <- args[6]
ILINet <- get_ili_data(lag = 1) %>%
	filter(region == "US", season %in% c(2010:2019, 2021))

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

all_results <- data.frame()

	results <- foreach(seas = c(2010:2019, 2021),
				.packages = c("dplyr", "tidyr", "stringr", "readr", "evalcast",
					      "scoringRules")
					      ,.errorhandling = "remove"
					      ,.combine = rbind
					      ) %:%
				foreach(wk = 9:38, .combine = rbind) %dopar% {
		    
			true_ili <- ILINet %>%
				filter(season == seas, season_week %in% c(wk + 1:4)) %>%
				ungroup() %>%
				select(unweighted_ili) %>% 	
				unlist() %>% as.vector()
		        
		        true_ili <- as.numeric(true_ili)/100

			#saveRDS(true_ili, "true_ili.rds")	
			forecast_file <- paste0("./ili_fits/", mod, "/", seas, 
						"/week", wk, ".csv")
			
			forecast <- read.csv(forecast_file) #%>% #select(-X)
	
			cont_scores <- data.frame()
			for (i in 2:5) {
		          true_value <- as.numeric(true_ili[i - 1])
		          kdens <- density(forecast[,i], kernel = "gaussian", from = -12, to = 12)
			  logs2 <- -log(approx(kdens$x, kdens$y, true_value)$y)
			  crps <- crps_sample(true_value, as.vector(forecast[,i]))
			  logs <- logs_sample(true_value, as.vector(forecast[,i]))
			  pit <- ecdf(as.vector(forecast[,i]))(true_value)
			  scores <- data.frame(season = seas, week_ahead = i - 1, 
					       crps = crps, logs = logs, logs2 = logs2, pit = pit)
			  cont_scores <- rbind(cont_scores, scores)
			}
			
			
		        colnames(forecast) <- paste0("week", 0:5)	
                        quantiles <- apply(forecast[,2:5], MARGIN = 2, quantile, probs = probs) %>% data.frame()
                        quantiles$probs <- probs
                        wis_cov <- quantiles %>% 
			  pivot_longer(1:4, names_to = "week_ahead", values_to = "pred_value") %>% 
			  mutate(week_ahead = parse_number(week_ahead)) %>% 
			  arrange(week_ahead) %>% 
			  group_by(week_ahead) %>% 
			  mutate(season_week = week_ahead + wk, season = seas) %>%# write.csv(wis_cov, "wis.csv") #%>% 
			  left_join(ILINet %>%
				   mutate(uw_ili = as.numeric(unweighted_ili)/100) %>% 
				      select(season, season_week, uw_ili, unweighted_ili), 
			      by = c("season", "season_week")) %>% # saveRDS(wis_cov, "wctest.rds")# %>% 
			  group_by(season, week_ahead) %>% 
			  reframe(
			    wis = weighted_interval_score(probs,
							  pred_value, 
							  uw_ili),
		            lmae = abs(median(pred_value - unique(uw_ili))),		  
			    mae = abs(median(exp(pred_value) - unique(exp(uw_ili)))), 
			    cover50 = coverage(probs, pred_value, unique(uw_ili), alpha = .5),
			    cover60 = coverage(probs, pred_value, unique(uw_ili), alpha = .4),
			    cover70 = coverage(probs, pred_value, unique(uw_ili), alpha = .3),
			    cover80 = coverage(probs, pred_value, unique(uw_ili), alpha = .2),
			    cover90 = coverage(probs, pred_value, unique(uw_ili), alpha = .1), 
			    cover95 = coverage(probs, pred_value, unique(uw_ili), alpha = .05),
			    cover98 = coverage(probs, pred_value, unique(uw_ili), alpha = .02),
			    cover99 = coverage(probs, pred_value, unique(uw_ili), alpha = .01),
			    cover99.5 = coverage(probs, pred_value, unique(uw_ili), alpha = .005),
			    cover100 = coverage(probs, pred_value, unique(uw_ili), alpha = 0)
			  )
		
                        results <- wis_cov %>% 
			  left_join(cont_scores, by = c("season", "week_ahead")) %>% 
			  mutate(last_week = wk, model = mod) %>% 
			  select(model, season, last_week, week_ahead, wis, crps, 
				 logs, logs2, lmae, mae, pit, cover50, cover60, cover70, cover80, cover90, cover95,
				 cover98, cover99, cover99.5, cover100) %>% data.frame()
		  	#all_results <- rbind(all_results, results)
			results
                     
		}
        #print(dim(results))
	print(class(results))	
	write.csv(results, "test.csv")
	all_results <- rbind(all_results, results)
	print(rep)

			

write.csv(all_results, paste0("./hosp_sim/ili_scores_", mod, ".csv"), row.names = FALSE)













