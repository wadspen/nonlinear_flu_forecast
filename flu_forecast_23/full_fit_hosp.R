setwd("./")
library(tidyr)
library(readr)
source("./get_data_functions_new.R")
source("./mle_functions.R")
library(parallel)
library(doParallel)
library(doMC)
library(evalcast)
#n.cores <- detectCores()
n.cores <- 64
my.cluster <- makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
registerDoMC(cores = n.cores)

both_flu <- comb_data(lag = 1) %>%
    filter(season >= 2010) #%>% select(-X)
both_flu2 <- comb_data(lag = 2) %>%
    filter(season >= 2010) #%>% select(-X)

ILINet <- get_ili_data() %>%
    filter(season >= 2010)

args <- commandArgs()
model <- args[6]
ili_model <- args[7]
mod_string <- paste("../stan_models/", model, ".stan", sep = "")
mod <- cmdstan_model(stan_file = mod_string)

#dir.create(paste0("./hosp_fits/", model))
select_regions <- unique(ILINet$region)
#select_regions <- c("Alabama", "Alaska")#, "Idaho")
#select_regions <- c("US")
hosp_log <- ifelse(str_detect(model, "log"), TRUE, FALSE)
hosp_lst <- ifelse(str_detect(model, "lst"), TRUE, FALSE)
sqr <- ifelse(str_detect(model, "sq"), TRUE, FALSE)
if (hosp_log == TRUE) {both_flu$value <- log(both_flu$value + 1)}
if (hosp_log == TRUE) {both_flu2$value <- log(both_flu2$value + 1)}

#select_regions <- "US"
all_forecasts <- foreach(j = select_regions,
	.packages = c("tidyr", "dplyr", "evalcast", "truncnorm", "truncdist")
	,.errorhandling = "remove"
	,.combine = rbind) %:% 
	#k <- 14
  foreach(k = 
	  c(14, 20, 26)
	  #9:38
  	, .combine = rbind) %dopar% {
# j <- "Georgia"; k <- 21 
       
    both_flu_hold <- both_flu %>%
	   filter(region == j, season == 2022 | (season == 2023 & season_week <= k)) 
      
    both_flu2_hold <- both_flu2 %>%
	   filter(region == j, season == 2022 | (season == 2023 & season_week <= k))

    stan_dat <- list(HM = nrow(both_flu_hold),
		     n_seasons_hosp = length(unique(both_flu_hold$season)),
		     ili_ps = both_flu_hold$unweighted_ili/100,
		     hosp = both_flu_hold$value,
		     count_rate = unique(both_flu_hold$count_rate2),
		     hosp_seasons = as.numeric(as.factor(both_flu_hold$season)),
		     sigma_alpha0 = 5,
		     sigma_alpha1 = 5,
		     sigma_alpha3 = .4,
		     sigma_sigma_epsilon = 4,
		     sigma_nu = 15
		 )
    location <- unique(both_flu_hold$location)
    # if (j != "US" & !(k %in% c(14, 20, 26))) {
        samps <- mod$sample(data = stan_dat,
				chains = 1,
				adapt_delta = .9999,
				iter_warmup = 10000,
				iter_sampling = 50000)

    	    draws <- samps$draws(format = "df"); #print("gets here")
    # }
    # if (j == "US" & k %in% c(14, 20, 26)) {
    #         samps <- mod$sample(data = stan_dat,
    #     			chains = 4,
    #     			adapt_delta = .9999,
    #     			iter_warmup = 10000,
    #     			iter_sampling = 50000)
    #         draws <- samps$draws(format = "df")
    # 
    #         summary <- samps$summary()
    #         diag_list <- lapply(samps$diagnostic_summary(), FUN = mean)
    #         diag <- diag_list %>% data.frame()
    #         diag$max_rhat <- max(summary$rhat, na.rm = TRUE)
    #         diag$min_ess <- min(summary$ess_bulk, na.rm = TRUE)
    #         diag$min_esst <- min(summary$ess_tail, na.rm = TRUE)
    #         write.csv(diag, paste0("./hosp_fits/", ili_model, "_", model, "_", j, "_", k,
    #     			   "_diagnostics.csv"), row.names = FALSE)
    # 
    #         write.csv(draws, paste0("./hosp_fits/", ili_model, "_", model, "_", j, "_", k,
    #     			    "_posterior.csv"), row.names = FALSE)
    # } 
    #write.csv(draws, "test_draws.csv")
    reference_date <- date(max(both_flu_hold$date)) + 7
    #dir.create(paste0("./hosp_fits/", model, 
    #	       "/reference_date_", reference_date, "_week", k))
    save_name <- paste("./hosp_fits/", model,
		       "/reference_date_", reference_date, "_week", k, "/",
		       j,  ".csv", sep = "")
	
    #write.csv(draws, save_name, row.names = FALSE)

    
    final_hosp <- last(both_flu_hold$value)
    current_hosp <- last(both_flu2_hold$value)
   ili_forc <- read.csv(paste0("./ili_fits/", ili_model, "/reference_date_", date(reference_date), 
			"_week", k, "/", j, ".csv"))
   
   ili_forc <- ili_forc %>%
     filter(if_all(everything(), ~ !is.na(.)))
   
   forecasts <- data.frame()
   count_rate <- unique(both_flu_hold$count_rate2)
   count_rate_sig <- count_rate
   #print("does it get here") 
   if (hosp_log == TRUE) {
	   alpha0 <- draws$alpha0
	   alpha1 <- draws$alpha1
	   if (sqr == TRUE) {alpha2 <- draws$alpha2}
	   alpha3 <- draws$alpha3
	   sigma_epsilon <- draws$sigma_epsilon
	   count_rate <- 1
	   count_rate_sig <- 1
   } else if (hosp_log == FALSE) {
	   alpha0 <- draws$`alpha0[2]`
	   alpha1 <- draws$`alpha1[2]`
	   if (sqr == TRUE) {alpha2 <- draws$`alpha2[2]`}
	   alpha3 <- draws$`alpha3[2]`
	   sigma_epsilon <- draws$`sigma_epsilon[2]`	   
   }
   
   if (hosp_lst == TRUE) {nu <- draws$`nu[2]`}
  
   #test <- data.frame(alpha0, alpha1, alpha2, alpha3, sigma_epsilon)

   #if (draws$alpha0.2. == NULL) {write.csv(dudish)}
   #if (hosp_log == FALSE) {write.csv(test, "test2_draws.csv")}
   #write.csv(draws, "dude.csv")
   samp_size <- nrow(ili_forc)
   #print("how about here")
   for (i in 1:samp_size) {
	alpha0_samp <- sample(alpha0, 1)
   	alpha1_samp <- sample(alpha1, 1)
	if (sqr == TRUE) {alpha2_samp <- sample(alpha2, 1)}
	alpha3_samp <- sample(alpha3, 1)
	sigma_epsilon_samp <- sample(sigma_epsilon, 1)
         
	ili_samp <- as.vector(ili_forc[sample(samp_size, 1), ])
	ili_samp <- unlist(ili_samp)
	
	#print(c(alpha0_samp, alpha1_samp, alpha2_samp, alpha3_samp, sigma_epsilon_samp))
	#print(dim(ili_samp))
	hosp_forecast <- c()
	for (wh in 1:5) {
		if (wh == 1) {
			hosp_lag <- final_hosp
		} 
		else if (wh == 2) {
			hosp_lag <- current_hosp
		}
		else {
			hosp_lag <- hosp_forecast[wh - 1]
		}
		
		mean <- alpha0_samp + 
			alpha1_samp*(ili_samp[wh + 1]*count_rate) +
			#alpha2_samp*(ili_samp[wh + 1]*count_rate)^2 +
			alpha3_samp*hosp_lag
		if (sqr == TRUE) mean <- mean + alpha2_samp*(ili_samp[wh + 1]*count_rate)^2
		if (hosp_lst == TRUE) {
			nu_samp <- sample(nu, 1)
			#hsamp <- rt(1, nu_samp)*(count_rate*sigma_epsilon_samp) + mean 
			#while (hsamp < 0) {
			#	hsamp <- rt(1, nu_samp)*(count_rate*sigma_epsilon_samp) + mean
			#print(hsamp)
			#}	
			a_std <- (0 - mean)/(count_rate*sigma_epsilon_samp)
			hosp_forecast[wh] <- rtrunc(1, spec = "t", a = a_std, b = Inf, df = nu_samp)*(count_rate*sigma_epsilon_samp) + mean
			#hosp_forecast[wh] <- hsamp
			#hosp_forecast[wh] <- rt(1, nu_samp)*(count_rate*sigma_epsilon_samp) + mean
		} else if (hosp_log == TRUE) {
			hosp_forecast[wh] <- rnorm(1, mean, count_rate*sigma_epsilon_samp)

		  }	
		else {
			#hsamp <- rnorm(1, mean, count_rate*sigma_epsilon_samp)
			#while (hsamp < 0) {
			#	hsamp <- rnorm(1, mean, count_rate*sigma_epsilon_samp)
			#print(hsamp)
			#}	
			#hosp_forecast[wh] <- hsamp
			hosp_forecast[wh] <- rtruncnorm(1, a = 0, b = Inf, mean = mean, sd = count_rate*sigma_epsilon_samp)
		}
		#if (hosp_log == TRUE) {hosp_forecast <- exp(hosp_forecast) - 1}

	}
	forecasts <- rbind(forecasts, hosp_forecast); #print(paste("iteration", i))
   }
   forecasts <- forecasts %>%
     filter(if_all(everything(), ~ !is.na(.)))
   
   colnames(forecasts) <- paste0("ahead", -1:3)
   #print(head(forecasts))
   probs <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
   quantile_forecasts <- apply(forecasts, MARGIN = 2, quantile, probs = probs) %>%
	   as.data.frame() 

   quantile_forecasts$output_type_id <- rownames(quantile_forecasts)
     #print(head(quantile_forecasts))
   long_quant_forc <- quantile_forecasts %>%
	   pivot_longer(1:5, names_to = "horizon", values_to = "value") %>%
	   mutate(reference_date = date(reference_date), horizon = parse_number(horizon),
	          output_type_id = parse_number(output_type_id)/100, 
		  target = "wk inc flu hosp", 
		  target_end_date = date(reference_date) + horizon*7,
	   output_type = "quantile", location = as.character(location)) %>%
	   rowwise() %>%
	   mutate(value = ifelse(hosp_log == TRUE, exp(value) - 1, value)) %>%
	   select(reference_date, horizon, target, target_end_date, location, 
	          output_type, output_type_id, value) %>%
	   mutate(value = ifelse(value < 0, 0, value)) %>% as.data.frame()
   write.csv(long_quant_forc, "test_forecasts.csv")
   long_quant_forc
	   #if (hosp_log == TRUE) {long_quant_forc$value <- exp(lonq_quant_forc_value) - 1}
	   #print(paste("the value is", hosp_log))
   	   #write.csv(long_quant_forc, "test_forecasts.csv")
	   #print(long_quant_forc); write.csv(long_quant_forc, "test_forecasts.csv")  
   #test <- cbind(dude = alpha0, dude1 = alpha1, alpha2, alpha3, sigma_epsilon)
   #test
}


#saveRDS(all_forecasts, "test.rds")
all_forecasts <- all_forecasts %>%
	#bind_rows() %>%
	group_split(reference_date)

#saveRDS(all_forecasts, "test.rds")
forecast_dir <- paste0("./fin_hosp_forecasts/", ili_model, "_", model)
if(dir.exists(forecast_dir) == FALSE) {
	     dir.create(forecast_dir)
      }

for (i in 1:length(all_forecasts)) {
	     ref <- unique(all_forecasts[[i]]$reference_date)
	     forc_file <- paste0("./fin_hosp_forecasts/", ili_model, "_", model, "/",
				 ref, "-", ili_model, "_", model,
				 ".csv")

	     write.csv(all_forecasts[[i]], forc_file, row.names = FALSE)
}



