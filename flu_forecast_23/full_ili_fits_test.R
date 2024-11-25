setwd("./")
#source("./simulation/simulate_hospitalizations.R")
source("./get_data_functions_new.R")
source("./mle_functions.R")
library(parallel)
library(doParallel)
library(doMC)
library(evalcast)
n.cores <- detectCores()
n.cores <- 1
my.cluster <- makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
registerDoMC(cores = n.cores)

both_flu <- comb_data(lag = 1) %>%
	filter(season >= 2010)
ILINet <- get_ili_data(lag = 1) %>%
	filter(season >= 2010, !(season %in% c(2020)))
args <- commandArgs()
model <- args[6]
tri_st <- as.numeric(args[7])*10
tri_st <- 12
mod_string <- paste("../stan_models/", model, ".stan", sep = "")
mod <- cmdstan_model(stan_file = mod_string)

tri_wks <- tri_st:(tri_st + 9)
print(tri_wks)
select_regions <- unique(both_flu$region)
#select_regions <- "US"
season_levels <- unique(ILINet$season)
season_levels <- season_levels[season_levels != 2023]
save_dir <- paste0("./ili_fits/", model)
print(save_dir)
if (dir.exists(save_dir) == FALSE) {
	dir.create(paste0("./ili_fits/", model))
}
j <- "Alabama"
k <- 14
foreach(j = select_regions,
	.packages = c("tidyr", "dplyr", "evalcast")
	,.errorhandling = "remove"
	) %:%
   foreach(k = tri_wks) %dopar% { #c(14, 20, 26)) %dopar% { #35 is the last week
#for (i in 1:(length(season_levels))) {
#  j <- "US"
#  k <- 26
      ILINet_arrange <- ILINet %>%
	  filter(region == j)
  

      ILINet_filter <- ILINet_arrange %>%
	      filter(season < 2023 | (season == 2023 & season_week <= k))
      #print(head(ILINet_filter))
      print("does it get here?")
      dat <- get_stan_data(ILINet, 
                           both_flu, s_region = j,
                           ili_seasons = unique(ILINet$season),
      			   m = k)
      
      print("how about here?")
      stan_dat <- dat[[1]]
      info_list <- dat[[2]]
      forecast_date <- info_list$forecast_date
      init <- list(theta = stan_dat$m0, theta_s = info_list$theta_s)
      reference_date <- info_list$forecast_date + 7
      forecast_dir <- paste0(save_dir, "/reference_date_",
			     reference_date, "_week", k)
      if(dir.exists(forecast_dir) == FALSE) {
	     dir.create(forecast_dir)
      }
      num_chains <- ifelse(j == "US" & k %in% c(14, 20, 26), 4, 1)
      if (j == "US" & k %in% c(14, 20, 26)) {
		inits <- list(init, init, init, init)
      } else {inits <- list(init)} 
      samps <- mod$sample(data = stan_dat,
                          chains = num_chains,
			  parallel_chains = num_chains,
                          adapt_delta = .95,
                          init = inits,
                          iter_warmup = 10000, 
                          iter_sampling = 50000)
 
      draws <- samps$draws(format = "df")
      pred_ili_draws <- draws[,str_detect(colnames(draws), "pred_ili.")]
      pred_draws <- pred_ili_draws[,(ncol(pred_ili_draws) - 5):ncol(pred_ili_draws)]
      
      save_name <- paste("./ili_fits/", model, 
      			 "/reference_date_", reference_date, 
			 "_week", k, "/",
      			 j,  ".csv", sep = "")

      write.csv(pred_draws, save_name, row.names = FALSE)
      
      if (j == "US" & k %in% c(14, 20, 26)) {
	  summary <- samps$summary()
          diag_list <- lapply(samps$diagnostic_summary(), FUN = mean)
          diag <- diag_list %>% data.frame()
          diag$max_rhat <- max(summary$rhat, na.rm = TRUE)
          diag$min_ess <- min(summary$ess_bulk, na.rm = TRUE)
          diag$min_esst <- min(summary$ess_tail, na.rm = TRUE)
          write.csv(diag, paste0("./ili_fits/", model, "/", j, "_", k,  
		    "_diagnostics.csv"), row.names = FALSE)	
   	  
          write.csv(draws, paste0("./ili_fits/", model, "/", j, "_", k,
			"posterior.csv"), row.names = FALSE)	  
      }

       
    
 
  
}


















