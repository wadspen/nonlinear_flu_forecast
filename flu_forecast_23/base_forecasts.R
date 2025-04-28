setwd("./flu_forecast_23/")
library(tidyr)
library(readr)
source("./get_data_functions_new.R")
source("./mle_functions.R")
library(evalcast)
library(forecast)


both_flu <- comb_data(lag = 1) %>%
  filter(season >= 2010) #%>% select(-X)
both_flu2 <- comb_data(lag = 2) %>%
  filter(season >= 2010) #%>% select(-X)

ILINet <- get_ili_data() %>%
  filter(season >= 2010)


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
N <- 5000
all_qforcs_res <- data.frame()
for (j in select_regions) { #[47:length(select_regions)]) {

  for (k in 9:38) {
    both_flu_hold <- both_flu %>%
      filter(region == j, season == 2022 | (season == 2023 & season_week <= k)) 
    
    forecast_date <- date(max(both_flu_hold$date))
  
    
    probs <- c(.01, .025, seq(.05, .95, by = .05), .975, .99)
    
    both_flu_reg <- both_flu %>% 
      filter(region == j) %>% 
      mutate(true_value = value) %>% 
      ungroup() %>% 
      select(date, true_value) %>% 
      mutate(date = date(date))
    
    basemod <- naive(both_flu_hold$value)
    basemean <- unique(as.vector(basemod$mean))
    basesd <- basemod$model$sigma2
    baseforcs <- sapply(1:4, FUN = function(x) {rnorm(N, basemean, 
                                                      sqrt(basesd*x))}) 
    
    baseqforcs <- apply(baseforcs, MARGIN = 2, quantile, probs = probs) %>% 
      as.data.frame()
    
    
    colnames(baseqforcs) <- as.character(forecast_date + 1:4*7)
    
    baseqforcsl <- baseqforcs %>% 
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
      mutate(week = k, model = "base", region = j)
    
    
    
    
    arima <- auto.arima(both_flu_hold$value)
    arimameans <- forecast(arima)$mean[1:4]
    low <- as.vector(forecast(arima)$lower[1:4,2])
    arimasds <- (arimameans - low)/qnorm(.975)
    
    
  
    arimaforcs <- apply(cbind(arimameans, arimasds), 
                      MARGIN = 1, FUN = function(x) {
                        rnorm(N, 
                              mean = x[1], 
                              sd = x[2])}) 
    arimaqforcs <- apply(arimaforcs, MARGIN = 2, quantile, probs = probs) %>% 
      as.data.frame()
    
    colnames(arimaqforcs) <- as.character(forecast_date + 1:4*7)
    
    arimaqforcsl <- arimaqforcs %>% 
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
      mutate(week = k, model = "arima", region = j)
    
    
    
    
    all_qforcs_res <- rbind(all_qforcs_res, baseqforcsl, arimaqforcsl)
    
  }
  print(j)
}

all_qforcs_res %>% 
  filter(region == "Mississippi") %>% 
  # group_by(week) %>% 
  # filter(date == min(date)) %>% 
  group_by(date, model) %>%
  summarise(mwis = mean(wis)) %>%
  ggplot() +
  geom_line(aes(x = date, y = mwis, colour = model)) +
  geom_hline(yintercept = 0)

write.csv(all_qforcs_res, "res.csv", row.names = FALSE)
