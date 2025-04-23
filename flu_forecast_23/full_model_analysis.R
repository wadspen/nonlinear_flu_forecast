setwd("./flu_forecast_23/")
library(tidyr)
library(readr)
source("./get_data_functions_new.R")
source("./mle_functions.R")
library(evalcast)


both_flu <- comb_data(lag = 1) %>%
  filter(season >= 2010) #%>% select(-X)
both_flu2 <- comb_data(lag = 2) %>%
  filter(season >= 2010) #%>% select(-X)

ILINet <- get_ili_data() %>%
  filter(season >= 2010)


mod <- cmdstan_model(stan_file = "../stan_models/asg_hosp_log_sq_ar1.stan")
mod2 <- cmdstan_model(stan_file = "../stan_models/hosp_log_ar1.stan")

#dir.create(paste0("./hosp_fits/", model))
select_regions <- unique(ILINet$region)
#select_regions <- c("Alabama", "Alaska")#, "Idaho")
#select_regions <- c("US")
hosp_log <- ifelse(str_detect(model, "log"), TRUE, FALSE)
hosp_lst <- ifelse(str_detect(model, "lst"), TRUE, FALSE)
sqr <- ifelse(str_detect(model, "sq"), TRUE, FALSE)
if (hosp_log == TRUE) {both_flu$value <- log(both_flu$value + 1)}
if (hosp_log == TRUE) {both_flu2$value <- log(both_flu2$value + 1)}


# all_forecasts <- foreach(j = select_regions,
#                          .packages = c("tidyr", "dplyr", "evalcast")
#                          ,.errorhandling = "remove"
#                          ,.combine = rbind) %:% 
#   #k <- 14
#   foreach(k = 
#             #c(14, 20, 26), .combine = rbind) %dopar% {
#             9:38
#           , .combine = rbind) %dopar% {
            # j <- "Georgia"; k <- 21 
j <- "US"
all_qforcs_res <- data.frame()
for (k in 13:38) {
            both_flu_hold <- both_flu %>%
              filter(region == j, season == 2022 | (season == 2023 & season_week <= k)) 
            
            both_flu2_hold <- both_flu2 %>%
              filter(region == j, season == 2022 | (season == 2023 & season_week <= k))
            
            stan_dat <- list(HM = nrow(both_flu_hold),
                             n_seasons_hosp = length(unique(both_flu_hold$season)),
                             ili_ps = both_flu_hold$unweighted_ili/100,
                             hosp = both_flu_hold$value,
                             count_rate = unique(both_flu_hold$count_rate1),
                             hosp_seasons = as.numeric(as.factor(both_flu_hold$season)),
                             sigma_alpha0 = 5,
                             sigma_alpha1 = 5,
                             sigma_alpha3 = .4,
                             sigma_sigma_epsilon = 4,
                             sigma_nu = 15
            )
            
            dat <- get_stan_data(ILINet, 
                                 both_flu, s_region = j,
                                 ili_seasons = unique(ILINet$season),
                                 m = k)
            
            stan_dat <- dat[[1]]; stan_dat$count_rate <- 12
            info_list <- dat[[2]]
            location <- unique(both_flu_hold$location)
            #if (j != "US" & !(k %in% c(14, 20, 26))) {
            init <- list(theta = stan_dat$m0, theta_s = info_list$theta_s)
            inits <- list(init)
            # samps <- mod$variational(stan_dat)
            samps <- mod$sample(data = stan_dat,
                                chains = 1,
                                # init = inits,
                                # adapt_delta = .9999,
                                iter_warmup = 5000,
                                iter_sampling = 6000)
            
            draws <- samps$draws(format = "df") %>% 
              select(contains("pred_hosp"))
            
            means <- apply(draws, MARGIN = 2, FUN = quantile, probs = .5)
            forc <- data.frame(value = means, week = 1:length(means))
            
            both_flu %>% 
              filter(season == 2023, region == "US") %>% 
              ggplot() +
              geom_point(aes(x = season_week, y = value)) +
              geom_line(data = forc, aes(x = week, y = value))
            
            
            both_flu_reg <- both_flu %>% 
              filter(region == j) %>% 
              mutate(true_value = value) %>% 
              ungroup() %>% 
              select(date, true_value) %>% 
              mutate(date = date(date))
            probs <- c(.01, .025, seq(.05, .95, by = .05), .975, .99)
            qforcs <- apply(draws[,k + 1:4], MARGIN = 2, quantile, probs = probs)
            colnames(qforcs) <- as.character(info_list$forecast_date + 1:4*7)
            qforcsl <- qforcs %>% 
              as.data.frame() %>% 
              tibble::rownames_to_column("prob") %>% 
              mutate(prob = as.numeric(str_replace(prob, "%", ""))/100) %>% 
              pivot_longer(2:5, names_to = "date") %>% 
              mutate(date = date(date)) %>% 
              arrange(date, prob) %>% 
              left_join(both_flu_reg, by = "date") %>% 
              group_by(date) %>% 
              summarise(wis = weighted_interval_score(prob, value, unique(true_value)),
                        cover = between(unique(true_value), value[prob == .025], value[prob == .975])) %>% 
              mutate(week = k)
            
            
            all_qforcs_res <- rbind(all_qforcs_res, qforcsl)
            
}

all_qforcs_res %>% 
  # group_by(week) %>% 
  # filter(date == min(date)) %>% 
  group_by(date) %>%
  summarise(mwis = mean(wis)) %>%
  ggplot() +
  geom_line(aes(x = date, y = mwis))


