library(dplyr)
library(tidyr)
library(parallel)
library(doParallel)
library(doMC)
library(evalcast)
n.cores <- detectCores()
#n.cores <- 1
my.cluster <- makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
registerDoMC(cores = n.cores)

args <- commandArgs()
model <- args[6]
print(model)



reps <- 500

ili_dir <- paste0("ili_fits/", model, "/")
work_dir <- "hosp_sim/hosp_log_ar1/"
sim_dir <- paste0(work_dir, "sim_hosp/")
forecast_dir <- paste0(work_dir, "forecasts/")
seas <- 2022
#for (rep in 143:reps) {
foreach(rep = 1:reps,
	.packages = c("dplyr", "tidyr")
	) %:% 
	foreach(wk = seq(14, 38, 6)) %dopar% {
        #foreach (seas = c(2010:2019, 2021:2022)) %dopar% {
	file <- list.files(sim_dir, pattern = paste0("rep", rep, "_"))
	sim_rep_file <- paste0(sim_dir, file)
	sim_hosp <- read.csv(sim_rep_file)
	#count_rate <- unique(sim_hosp$count_rate2)
	count_rate <- 1
	#for (seas in c(2010:2019, 2021:2022)) {
		sim_season_dir <- paste0(work_dir, "new/", seas, "/")
		#for (wk in seq(14, 38, 6)) {
			sim_season_week_dir <- paste0(sim_season_dir, "week", wk, "/")
		
				sim_season_week_rep_file <- paste0(sim_season_week_dir, "trep",
								   rep, ".csv")
				if (file.exists(sim_season_week_rep_file) == FALSE) {next}
				sim_forecast <- read.csv(sim_season_week_rep_file)

				ili_fit_file <- paste0(ili_dir, seas, "/week", wk, ".csv")
				ili_fit <- read.csv(ili_fit_file)
				
				data <- sim_hosp %>%
				    filter(season == seas, season_week <= wk)
			  		
			    	final_hosp <- data %>%
				    filter(season_week == wk) %>%
			    	    select(hosp_sim)
				final_hosp <- log(final_hosp$hosp_sim + 1)

			        alpha0 <- sim_forecast$alpha0
				alpha1 <- sim_forecast$alpha1
				#alpha2 <- sim_forecast$alpha2
				alpha3 <- sim_forecast$alpha3
				sigma_epsilon <- sim_forecast$sigma_epsilon
					
				#forecasts <- matrix(NA, nrow = nrow(ili_fit), ncol = 4)
				forecasts <- data.frame()
				samp_size <- nrow(ili_fit)
				#ili_samp <- as.vector(apply(ili_fit, MARGIN = 2, mean))
				#ili_samp <- unlist(ili_samp)*100
				#saveRDS(ili_samp, "ili.rds")
				#saveRDS(ili_samp, "test.rds")
				for (i in 1:samp_size) {
					alpha0_samp <- sample(alpha0, 1)
					alpha1_samp <- sample(alpha1, 1)
					#alpha2_samp <- sample(alpha2, 1)
					alpha3_samp <- sample(alpha3, 1)
					sigma_epsilon_samp <- sample(sigma_epsilon, 1)
					
					ili_samp <- as.vector(ili_fit[sample(samp_size, 1), ])
					ili_samp <- unlist(ili_samp)
					ili_samp <- ili_samp*100
					#saveRDS(print(ili_samp), "test.rds")
					hosp_forecast <- c()
					#ili_samp <- as.vector(ili_fit[sample(samp_size, 1), ])
					#ili_samp <- unlist(ili_samp)
					
					for (wh in 1:4) {
						if (wh == 1) {
							hosp_lag <- final_hosp
						}
						else {
							#mean_lag <- alpha0_samp +
							#     alpha1_samp*(ili_samp[wh]*count_rate) +
							#     alpha2_samp*(ili_samp[wh]*count_rate)^2 +
							#     alpha3_samp*hosp_forecast[wh - 1]
					     		
							#hosp_lag <- rnorm(1, mean_lag, sigma_epsilon_samp)
						     	hosp_lag <- hosp_forecast[wh - 1]
						}

						mean <- alpha0_samp + 
							alpha1_samp*(ili_samp[wh + 1]*count_rate) +
							#alpha2_samp*(ili_samp[wh + 1]*count_rate)^2 +
							alpha3_samp*hosp_lag
							
						
						#hosp_forecast[wh] <- mean
						hosp_forecast[wh] <- rnorm(1, mean, count_rate*sigma_epsilon_samp)
					
					}
					forecasts <- rbind(forecasts, hosp_forecast)	
				}
				colnames(forecasts) <- paste0("ahead", 0:3)
				save_file_dir <- paste0(work_dir, "new/forecasts/", model, "/", seas, 
							"/week", wk, "/rep", rep, ".csv")
				write.csv(forecasts, save_file_dir, row.names = FALSE)
					
		#}			
	#}
}

