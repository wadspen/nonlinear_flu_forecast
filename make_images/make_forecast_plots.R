library(gridExtra)
setwd("~/flu_research/nonlinear_flu_forecast/flu_forecast_23/fin_hosp_forecasts/")
source("~/flu_research/nonlinear_flu_forecast/flu_forecast_23/get_data_functions_new.R")
files <- list.files("./")

state <- "US"


ILINet <- get_ili_data()
both_flu <- comb_data(lag = 1)

plot_dates <- both_flu %>% ungroup() %>% filter(season_week %in% seq(13, 37, by = 6),
                                  season == 2023) %>%
  select(date) %>% unique() %>% unlist()

locations <- both_flu %>% ungroup() %>%
  select(location, location_name) %>% unique()


#everything below can be changed for each image
#uncomment everything above when first running the progam
model <- "asg_disc2_nm_hosp_ar1"
models <- c("asg_nm_hosp_sq_ar1", "asg_disc2_nm_hosp_sq_ar1", 
            "sir_hosp_sq_ar1", "sir_disc2_hosp_sq_ar1")
# models <- "FluSight-baseline"
labels <- c("ASG", "ASGD", "SIR", "SIRD")
labels <- factor(labels, levels = c("ASGD", "ASG", "SIRD", "SIR"))
# models <- "asg_disc2_nm_hosp_sq_ar1"
# lables <- "ASGD"
# labels <- "SIRD"
# models <- c("sir_hosp_log_ar1", "sir_disc_hosp_log_ar1")
# models <- "FluSight-baseline"
# labels <- c("SIR", "SIRD")
# model <- "FluSight-baseline"
forecast_file1 <- paste(models, "/", plot_dates[1], "-",
                       model, ".csv", sep = "")

read_model <- function() {
  read.csv() %>% 
    mutate(model = model)
}

# forecast_file1 <- 
#   do.call(rbind,
#           lapply(paste0(models, "/", plot_dates[1], "-",
#                         models, ".csv"), read.csv))

forecast1 <- data.frame()
for (mod in 1:length(models)) {
  file_name <- paste0(models[mod], "/", plot_dates[1], "-",
                 models[mod], ".csv") 
  file <- read.csv(file_name) %>% 
    mutate(model = labels[mod])
  forecast1 <- rbind(forecast1, file)
}

forecast2 <- data.frame()
for (mod in 1:length(models)) {
  file_name <- paste0(models[mod], "/", plot_dates[2], "-",
                      models[mod], ".csv") 
  file <- read.csv(file_name) %>% 
    mutate(model = labels[mod])
  forecast2 <- rbind(forecast2, file)
}

forecast3 <- data.frame()
for (mod in 1:length(models)) {
  file_name <- paste0(models[mod], "/", plot_dates[3], "-",
                      models[mod], ".csv") 
  file <- read.csv(file_name) %>% 
    mutate(model = labels[mod])
  forecast3 <- rbind(forecast3, file)
}

forecast4 <- data.frame()
for (mod in 1:length(models)) {
  file_name <- paste0(models[mod], "/", plot_dates[4], "-",
                      models[mod], ".csv") 
  file <- read.csv(file_name) %>% 
    mutate(model = labels[mod])
  forecast4 <- rbind(forecast4, file)
}

forecast5 <- data.frame()
for (mod in 1:length(models)) {
  file_name <- paste0(models[mod], "/", plot_dates[5], "-",
                      models[mod], ".csv") 
  file <- read.csv(file_name) %>% 
    mutate(model = labels[mod])
  forecast5 <- rbind(forecast5, file)
}


# state <- c("Alabama", "District of Columbia", "Iowa", 
#            "Montana", "Ohio", "Texas")
# state <- "Texas"
# forecast1 <- read.csv(forecast_file1) %>%
forecast1 <- forecast1 %>% 
  mutate(location = as.character(location),
         output_type_id = as.character(output_type_id)) %>%
  mutate(location = ifelse(nchar(location) < 2, paste0("0",location),
                           location)) %>%
  # filter(model == "SIRD") %>% 
  left_join(locations, by = "location") %>% 
  filter(output_type == "quantile", horizon != -1) %>%
  mutate(location = as.character(location)) %>% # %>% select(reference_date, horizon, target, target_end_date, location, output_type, output_type_id, value)
  mutate(value = ifelse(value < 0, 0, value)) %>% 
  filter( location_name %in% state, output_type_id %in% c("0.025", "0.25",
                                                         "0.75", "0.975")) %>%
  mutate(output_type_id = ifelse(output_type_id == "0.025", "lower_95",
                                        ifelse(output_type_id == "0.25", "lower_50",
                                               ifelse(output_type_id == "0.75", "upper_50", "upper_95")))) %>% 
  unique() %>% 
  pivot_wider(names_from = output_type_id)

forecast2 <- forecast2 %>%
  mutate(location = as.character(location), 
         output_type_id = as.character(output_type_id)) %>% 
  mutate(location = ifelse(nchar(location) < 2, paste0("0",location),
                           location)) %>% 
  # filter(model == "SIRD") %>% 
  left_join(locations, by = "location") %>% 
  filter(output_type == "quantile", horizon != -1) %>%
  mutate(location = as.character(location)) %>% # %>% select(reference_date, horizon, target, target_end_date, location, output_type, output_type_id, value)
  mutate(value = ifelse(value < 0, 0, value)) %>% 
  filter( location_name %in% state, output_type_id %in% c("0.025", "0.25",
                                                       "0.75", "0.975")) %>%
  mutate(output_type_id = ifelse(output_type_id == "0.025", "lower_95",
                                 ifelse(output_type_id == "0.25", "lower_50",
                                        ifelse(output_type_id == "0.75", "upper_50", "upper_95")))) %>% 
  unique() %>% 
  pivot_wider(names_from = output_type_id)


forecast3 <- forecast3 %>%
  mutate(location = as.character(location), 
         output_type_id = as.character(output_type_id)) %>%
  mutate(location = ifelse(nchar(location) < 2, paste0("0",location),
                           location)) %>% 
  # filter(model == "SIRD") %>% 
  left_join(locations, by = "location") %>% 
  filter(output_type == "quantile", horizon != -1) %>%
  mutate(location = as.character(location)) %>% # %>% select(reference_date, horizon, target, target_end_date, location, output_type, output_type_id, value)
  mutate(value = ifelse(value < 0, 0, value)) %>% 
  filter( location_name %in% state, output_type_id %in% c("0.025", "0.25",
                                                       "0.75", "0.975")) %>%
  mutate(output_type_id = ifelse(output_type_id == "0.025", "lower_95",
                                 ifelse(output_type_id == "0.25", "lower_50",
                                        ifelse(output_type_id == "0.75", "upper_50", "upper_95")))) %>% 
  unique() %>% 
  pivot_wider(names_from = output_type_id)


forecast4 <- forecast4 %>%
  mutate(location = as.character(location), 
         output_type_id = as.character(output_type_id)) %>%
  mutate(location = ifelse(nchar(location) < 2, paste0("0",location),
                           location)) %>% 
  # filter(model == "SIRD") %>% 
  left_join(locations, by = "location") %>% 
  filter(output_type == "quantile", horizon != -1) %>%
  mutate(location = as.character(location)) %>% # %>% select(reference_date, horizon, target, target_end_date, location, output_type, output_type_id, value)
  mutate(value = ifelse(value < 0, 0, value)) %>% 
  filter( location_name %in% state, output_type_id %in% c("0.025", "0.25",
                                                       "0.75", "0.975")) %>%
  mutate(output_type_id = ifelse(output_type_id == "0.025", "lower_95",
                                 ifelse(output_type_id == "0.25", "lower_50",
                                        ifelse(output_type_id == "0.75", "upper_50", "upper_95")))) %>% 
  unique() %>% 
  pivot_wider(names_from = output_type_id)

forecast5 <- forecast5 %>%
  mutate(location = as.character(location), 
         output_type_id = as.character(output_type_id)) %>%
  mutate(location = ifelse(nchar(location) < 2, paste0("0",location),
                           location)) %>% 
  # filter(model == "SIRD") %>% 
  left_join(locations, by = "location") %>% 
  filter(output_type == "quantile", horizon != -1) %>%
  mutate(location = as.character(location)) %>% # %>% select(reference_date, horizon, target, target_end_date, location, output_type, output_type_id, value)
  mutate(value = ifelse(value < 0, 0, value)) %>% 
  filter( location_name %in% state, output_type_id %in% c("0.025", "0.25",
                                                       "0.75", "0.975")) %>%
  mutate(output_type_id = ifelse(output_type_id == "0.025", "lower_95",
                                 ifelse(output_type_id == "0.25", "lower_50",
                                        ifelse(output_type_id == "0.75", "upper_50", "upper_95")))) %>% 
  unique() %>% 
  pivot_wider(names_from = output_type_id)

                                                     

both_flu %>% 
  # filter(location_name != "US") %>% 
  filter(season == 2023, location_name %in% state, date > "2023-10-07",
         date < "2024-06-01") %>% 
  # filter(season_week < 14) %>% 
  ggplot() +
  geom_ribbon(data = forecast1, aes(x = date(target_end_date), ymin = lower_95,
                  ymax = upper_95, fill = "lightpink")) +
  geom_ribbon(data = forecast1, aes(x = date(target_end_date), ymin = lower_50,
                  ymax = upper_50), fill = "pink") +
  geom_ribbon(data = forecast2, aes(x = date(target_end_date), ymin = lower_95,
                                    ymax = upper_95, fill = "lightpink")) +
  geom_ribbon(data = forecast2, aes(x = date(target_end_date), ymin = lower_50,
                                    ymax = upper_50), fill = "pink") +
  geom_ribbon(data = forecast3, aes(x = date(target_end_date), ymin = lower_95,
                                    ymax = upper_95, fill = "lightpink")) +
  geom_ribbon(data = forecast3, aes(x = date(target_end_date), ymin = lower_50,
                                    ymax = upper_50), fill = "pink") +
  geom_ribbon(data = forecast4, aes(x = date(target_end_date), ymin = lower_95,
                                    ymax = upper_95, fill = "lightpink")) +
  geom_ribbon(data = forecast4, aes(x = date(target_end_date), ymin = lower_50,
                                    ymax = upper_50), fill = "pink") +
  geom_ribbon(data = forecast5, aes(x = date(target_end_date), ymin = lower_95,
                                    ymax = upper_95, fill = "lightpink")) +
  geom_ribbon(data = forecast5, aes(x = date(target_end_date), ymin = lower_50,
                                    ymax = upper_50), fill = "pink") +
  
  geom_point(aes(x = date(date), y = value)) +
  # geom_vline(xintercept = date(plot_dates[1])) +
  # facet_wrap(~location_name, scales = "free") +
  facet_wrap(~model, scales = "free") +
  ylab("Hospitalizatios") +
  xlab("Date") +
  # ylim(c(0,25000)) +
  # coord_cartesian(ylim = c(0, 40000)) +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=18),
        legend.position = "none",
        strip.text = element_text(size = 14))
  
  

  
  