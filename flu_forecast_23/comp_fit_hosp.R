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
# ili_model <- args[7]
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

select_regions <- "US"
all_forecasts <- foreach(j = select_regions,
                         .packages = c("tidyr", "dplyr", "evalcast")
                         ,.errorhandling = "remove"
                         ,.combine = rbind) %:% 
  #k <- 14
  foreach(k = 
            c(14, 20, 26)
            #9:38
          , .combine = rbind) %dopar% {
            # j <- "Georgia"; k <- 21 
            
            
            dat <- get_stan_data(ILINet, 
                                 both_flu, s_region = j,
                                 ili_seasons = unique(ILINet$season),
                                 m = k)
            
            
            stan_dat <- dat[[1]]
            info_list <- dat[[2]]
            location <- j
            if (j != "US" & !(k %in% c(14, 20, 26))) {
            samps <- mod$sample(data = stan_dat,
                                chains = 1,
                                # adapt_delta = .99,
                                iter_warmup = 1000,
                                iter_sampling = 1000)
            } 
            #samps <- mod$variational(data = stan_dat)
            if (j == "US" & k %in% c(14, 20, 26)) {
            samps <- mod$sample(data = stan_dat,
        			chains = 4,
				parallel_chains = 4,
        			adapt_delta = .999,
        			iter_warmup = 1000,
        			iter_sampling = 1000)
            draws <- samps$draws(format = "df")
	    write.csv(draws, "test_draws.csv", row.names = FALSE)
            summary <- samps$summary()
	    write.csv(samps$diagnostic_summary(), "test_diagnostics.csv", 
		      row.names = FALSE)
            diag_list <- lapply(samps$diagnostic_summary(), FUN = mean)
            diag <- diag_list %>% data.frame()
            diag$max_rhat <- max(summary$rhat, na.rm = TRUE)
            diag$min_ess <- min(summary$ess_bulk, na.rm = TRUE)
            diag$min_esst <- min(summary$ess_tail, na.rm = TRUE)
            write.csv(diag, paste0("./hosp_fits/", model, "_", j, "_", k,
        			   "_full",  "_diagnostics.csv"), row.names = FALSE)

            write.csv(draws, paste0("./hosp_fits/", model, "_", j, "_", k,
        			    "_full", "_posterior.csv"), row.names = FALSE)
    } 
            draws <- samps$draws(format = "df")
            
            # reference_date <- date(max(both_flu_hold$date)) + 7
            reference_date <- info_list$forecast_date + 7
            #dir.create(paste0("./hosp_fits/", model, 
            #	       "/reference_date_", reference_date, "_week", k))
            save_name <- paste("./hosp_fits/", model, "_full",
                               "/reference_date_", reference_date, "_week", k, "/",
                               j,  ".csv", sep = "")
            
            #write.csv(draws, save_name, row.names = FALSE)
            
            
            forecasts <- draws %>% 
              select(contains("pred_hosp"))
            
            forecasts <- forecasts[,k:(k + 4)]
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
            # write.csv(long_quant_forc, "test_forecasts.csv")
            long_quant_forc
           
          }


#saveRDS(all_forecasts, "test.rds")
all_forecasts <- all_forecasts %>%
  #bind_rows() %>%
  group_split(reference_date)

#saveRDS(all_forecasts, "test.rds")
forecast_dir <- paste0("./fin_hosp_forecasts/", model, "_full")
if(dir.exists(forecast_dir) == FALSE) {
  dir.create(forecast_dir)
}

for (i in 1:length(all_forecasts)) {
  ref <- unique(all_forecasts[[i]]$reference_date)
  forc_file <- paste0("./fin_hosp_forecasts/", model, "_full/",
                      ref, "-", model,
                      "_full.csv")
  
  write.csv(all_forecasts[[i]], forc_file, row.names = FALSE)
}



