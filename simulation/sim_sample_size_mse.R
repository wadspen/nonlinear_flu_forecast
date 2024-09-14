library(dplyr)
library(readr)
library(stringr)
library(tidyr)


args <- commandArgs()
model <- args[6]

sim_dir <- "./hosp_sim/hosp_log_sq_ar1/sim_hosp/"
for_dir <- paste0("./hosp_sim/hosp_log_sq_ar1/new/forecasts/", model, "/")
all_vals <- data.frame()
sum_assess <- data.frame()
for (rep in 1:150) {
	file <- list.files(sim_dir, pattern = paste0("rep", rep, "_"))
	sim <- read.csv(paste0(sim_dir, file))
	for (seas in c(2010:2019, 2021)) {
		for (wk in seq(14, 38, 6)) {
			for (wh in 1:4) {
				sim_val <- sim %>%
			    	    filter(season == seas, season_week == wk + wh) %>%
			    	    select(season, season_week, hosp_sim)
				print(sim_val)
			    all_vals <- rbind(all_vals, sim_val)

			    forc_file <- paste0(for_dir, seas, "/week", wk, "/rep", rep, ".csv")
			    forc <- read.csv(forc_file)
			    #print(dim(forc))
			    forc_long <- forc %>% 
				    pivot_longer(everything(), names_to = "week_ahead", 
						 values_to = "hosp_forc") %>%
				    mutate(week_ahead = parse_number(week_ahead) + 1,
				           season = seas, forecast_week = wk) %>%
				    mutate(season_week = forecast_week + week_ahead) %>%
				    select(season, forecast_week, season_week, week_ahead, hosp_forc) %>%
				    left_join(sim, by = c("season", "season_week")) %>%
				    select(-ili, -unweighted_ili, -region, -count_rate2)

			    sum_forc_long <- forc_long %>%
				    mutate(abs_diff = abs((exp(hosp_forc) - 1) - (exp(hosp_sim) - 1)), 
						     sq_diff = ((exp(hosp_forc) - 1) - 
								(exp(hosp_sim) - 1))^2,
					   labs_diff = abs(hosp_forc - hosp_sim),
					   lsq_diff = (hosp_forc - hosp_sim)^2) %>%
			    	    ungroup() %>%
			            group_by(season, forecast_week, season_week) %>%
				    summarise(
					      mae = mean(abs_diff),
					      mse = mean(sq_diff),
					      male = mean(labs_diff),
					      msle = mean(lsq_diff)
				    )
			   sum_assess <- rbind(sum_assess, sum_forc_long) 
			}
		}
	}
}

write.csv(sum_assess, paste0("./hosp_sim/sample_size_mse_", model, ".csv"), row.names = FALSE)
