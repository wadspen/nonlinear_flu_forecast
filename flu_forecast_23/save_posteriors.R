setwd("./")
#source("./simulation/simulate_hospitalizations.R")
source("./get_data_functions_new.R")
source("./mle_functions.R")
library(parallel)
library(doParallel)
library(doMC)
library(evalcast)
#n.cores <- detectCores()
n.cores <- 10
my.cluster <- makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()
registerDoMC(cores = n.cores)

# both_flu <- comb_data(lag = 1) %>%
#   filter(season >= 2010)
# ILINet <- get_ili_data(lag = 1) %>%
#   filter(season >= 2010, !(season %in% c(2020)))

both_flu <- readRDS("./both_flu.rds")
ILINet <- readRDS("./ILINet.rds")
args <- commandArgs()
model <- args[6]
k <- as.numeric(args[7])
mod_string <- paste("../stan_models/", model, ".stan", sep = "")
mod <- cmdstan_model(stan_file = mod_string)




j <- "US"
#k <- 26
print(model)
print(j)
print(k)
    ILINet_arrange <- ILINet %>%
      filter(region == j)
    
    
    ILINet_filter <- ILINet_arrange %>%
      filter(season < 2023 | (season == 2023 & season_week <= k)) %>%
      filter(season >= 2010 & season != 2020)

    dat <- get_stan_data(ILINet, 
                         both_flu, s_region = j,
                         #ili_seasons = unique(ILINet_filter$season),
                         ili_seasons = c(2010:2019, 2021:2023),
			 m = k)
    
    stan_dat <- dat[[1]]
    info_list <- dat[[2]]
    forecast_date <- info_list$forecast_date
    init <- list(theta = stan_dat$m0, theta_s = info_list$theta_s)
    # if (j == "US") init$theta_s <- init$theta_s[-24,]
    reference_date <- info_list$forecast_date + 7
    if (model == "asg_disc2_nm") {stan_dat$m0[1] <- mean(info_list$theta_s[,1]); stan_dat$C0[1,1] <- .005}
    if (model == "asg_nm") {stan_dat$C0[1,1] <- 0.3} 
    inits <- list(init, init, init, init)
    samps <- mod$sample(data = stan_dat,
                        chains = 4,
                        parallel_chains = 4,
                        adapt_delta = .9999,
                        init = inits,
                        iter_warmup = 10000, 
                        iter_sampling = 50000)
    
    sum <- samps$summary(variables = NULL,
                         posterior::default_summary_measures()[1:4],
                         quantiles = ~ quantile2(., probs = c(0.025, 0.975)),
                         posterior::default_convergence_measures())
    
 saveRDS(sum, paste0(model, "2_", k, "_summary.rds"))
















