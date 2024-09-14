
library(dplyr)

sim_dir <- "./hosp_sim/hosp_log_sq_ar1/sim_hosp/"
all_vals <- data.frame()
for (rep in 1:150) {
	file <- list.files(sim_dir, pattern = paste0("rep", rep, "_"))
	sim <- read.csv(paste0(sim_dir, file))
	for (seas in c(2010:2019, 2021)) {
		for (wk in seq(14, 38, 6)) {
			for (wh in 1:4) {
				sim_val <- sim %>%
			    	    filter(season == seas, season_week == wk + wh) %>%
			    	    select(season, season_week, hosp_sim)

			    all_vals <- rbind(all_vals, sim_val)
			}
		}
	}
}




write.csv(all_vals, "./hosp_sim/all_forward_samp150.csv")

#with sample size of 1000, the max margin of error is 2.422 percent the value of the mean grouped by season and week
summarise_vals <- all_vals %>% 
	group_by(season, season_week) %>%
	summarise(mean = mean(exp(hosp_sim) - 1), sd = sd(exp(hosp_sim) - 1)) %>%
	mutate(ME = 1.96*(sd/sqrt(1000))) %>%
	mutate(perc = 100*ME/mean)

print(max(summarise_vals$perc))





