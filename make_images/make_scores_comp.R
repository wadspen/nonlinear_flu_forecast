setwd("~/forecast-hub/")
source("../flu_research/nonlinear_flu_forecast/flu_forecast_23/get_data_functions_new.R")
# setwd("./model-output/")
# source("../target-data/get_target_data.R")
library(stringr)
library(epidatr)
library(evalcast)
library(dplyr)
library(tidyr)

both_flu <- comb_data()
date_weeks <- both_flu %>%
  select(date, season_week) %>% unique
ILINet <- get_ili_data()


setwd("./FluSight-forecast-hub/")

#season_week 10 is "2023-10-14"

target_data <- read.csv("./target-data/target-hospital-admissions.csv") %>%
  mutate(true_value = value) %>% 
  dplyr::select(-value)
setwd("./model-output/")
folders <- list.files("./")
models <- folders[folders != "old_models" & folders != "target-data"]
sub_dates <- substr(list.files("./FluSight-baseline"), 1, 10)
sub_dates <- sub_dates[sub_dates < "2024-07-31"]
sub_date <- sub_dates[1]
base <- "FluSight-baseline"



non_ens <- c("CEPH-Rtrend_fluH",
             "CFA_Pyrenew-Pyrenew_HE_Flu",
             "CFA_Pyrenew-Pyrenew_H_Flu",
             "FluSight-baseline",
             "GH-model",
             "GT-FluFNP",
             "JHUAPL-DMD",
             "JHUAPL-Morris",
             "LUcompUncertLab-chimera",
             "LosAlamos_NAU-CModel_Flu",
             "MDPredict-SIRS",
             "MOBS-GLEAM_FLUH",
             "NEU_ISI-FluBcast",
             "NIH-Flu_ARIMA",
             "NU_UCSD-GLEAM_AI_FLUH",
             "PSI-PROF",
             "PSI-PROF_beta",
             "SigSci-BECAM",
             "SigSci-CREG",
             "UGA_CEID-Walk",
             "UGA_flucast-Copycat",
             "UGA_flucast-INFLAenza",
             "UGA_flucast-OKeeffe",
             "UGA_flucast-Scenariocast",
             "UGuelph-CompositeCurve",
             "UM-DeepOutbreak",
             "UMass-AR2",
             "UNC_IDD-InfluPaint",
             "UVAFluX-CESGCN",
             "VTSanghani-PRIME",
             "cfa-flumech",
             "cfarenewal-cfaepimlight"
            )

models <- non_ens




baseline_forecasts <- data.frame()
for (sub_date in sub_dates) {
  baseline_file <- paste(base, "/", sub_date, "-",
                         base, ".csv", sep = "")
  
  baseline_forecast <- read.csv(baseline_file)
  baseline_forecasts <- rbind(baseline_forecasts, baseline_forecast)
}

baseline_wis <- baseline_forecasts %>%
  filter(horizon %in% 0:3) %>%
  mutate(target_end_date = as.character(date(target_end_date))) %>%
  left_join(target_data,
            by = c("location", "target_end_date" = "date")) %>%
  mutate(output_type_id = as.numeric(output_type_id),
         value = as.numeric(value),
         true_value = as.numeric(true_value)) %>%
  filter(is.na(location_name) == FALSE) %>%
  group_by(reference_date, location_name, horizon) %>%
  summarise(
    bwis = weighted_interval_score(output_type_id,
                                   value, true_value),
    blwis = weighted_interval_score(output_type_id,
                                    log(value + 1),
                                    log(true_value + 1))
  ) %>%
  mutate(model = base)

team_scores <- data.frame()
for (sub_date in sub_dates) {
  for (i in 1:length(models)) {
    forecast_file <- paste(models[i], "/", sub_date, "-",
                           models[i], ".csv", sep = "")
    
    if (file.exists(forecast_file) == FALSE) next
    forecast <- read.csv(forecast_file) %>%
      mutate(location = as.character(location),
             output_type_id = as.character(output_type_id)) %>%
      mutate(location = ifelse(nchar(location) < 2, paste0("0",location),
                               location)) %>%
      filter(output_type == "quantile", horizon != -1) %>%
      mutate(location = as.character(location)) %>% # %>% select(reference_date, horizon, target, target_end_date, location, output_type, output_type_id, value)
      mutate(value = ifelse(value < 0, 0, value)) #%>%
    # mutate(target_end_date = as.character(date(target_end_date) - 7))
    
    for_wis <- forecast %>%
      filter(horizon %in% 0:3) %>%
      left_join(target_data,
                by = c("location", "target_end_date" = "date")) %>%
      mutate(output_type_id = as.numeric(output_type_id),
             value = as.numeric(value),
             true_value = as.numeric(true_value))
    
    # if (mean(is.na(for_wis$location_name)) == 0) next
    if (nrow(for_wis) == 0) next
    wis <- for_wis %>%
      filter(is.na(location_name) == FALSE) %>%
      group_by(reference_date, location_name, horizon) %>%
      reframe(
        wis = weighted_interval_score(output_type_id,
                                      value, true_value),
        lwis = weighted_interval_score(output_type_id,
                                       log(value + 1), log(true_value + 1)),
        mae = abs(median(value) - unique(true_value)),
        width = sort(value, partial = 17)[17] - sort(value, partial = 7)[7]
      ) %>%
      mutate(model = models[i]) %>%
      left_join(unique(ILINet[c("region", "count_rate2")]),
                by = c("location_name" = "region"))
    
    wis <- wis %>%
      left_join(baseline_wis[c("reference_date", "location_name",
                               "horizon", "bwis", "blwis")],
                by = c("reference_date", "location_name", "horizon")) %>%
      mutate(relative_wis = wis/bwis,
             relative_lwis = lwis/blwis)
    
    team_scores <- rbind(team_scores, wis)
    
  }
  print(sub_date)
}




setwd("../../../flu_research/nonlinear_flu_forecast/flu_forecast_23/fin_hosp_forecasts/")

folders <- list.files("./")
folders <- folders[!str_detect(folders, "full") & 
                   !str_detect(folders, "_disc_")]
folders <- folders[!str_detect(folders, "disc2") | 
                     str_detect(folders, "disc2_nm")]

models <- folders






ili_scores <- data.frame()
for (sub_date in sub_dates) {
  for (i in 1:length(models)) {
    forecast_file <- paste(models[i], "/", sub_date, "-",
                           models[i], ".csv", sep = "")
    
    if (file.exists(forecast_file) == FALSE) next
    forecast <- read.csv(forecast_file) %>%
      mutate(location = as.character(location),
             output_type_id = as.character(output_type_id)) %>%
      mutate(location = ifelse(nchar(location) < 2, paste0("0",location),
                               location)) %>%
      filter(output_type == "quantile", horizon != -1) %>%
      mutate(location = as.character(location)) %>% # %>% select(reference_date, horizon, target, target_end_date, location, output_type, output_type_id, value)
      mutate(value = ifelse(value < 0, 0, value)) #%>%
    # mutate(target_end_date = as.character(date(target_end_date) - 7))
    
    for_wis <- forecast %>%
      filter(horizon %in% 0:3) %>%
      left_join(target_data,
                by = c("location", "target_end_date" = "date")) %>%
      mutate(output_type_id = as.numeric(output_type_id),
             value = as.numeric(value),
             true_value = as.numeric(true_value))
    
    # if (mean(is.na(for_wis$location_name)) == 0) next
    if (nrow(for_wis) == 0) next
    wis <- for_wis %>%
      filter(is.na(location_name) == FALSE) %>%
      group_by(reference_date, location_name, horizon) %>%
      reframe(
        wis = weighted_interval_score(output_type_id,
                                      value, true_value),
        lwis = weighted_interval_score(output_type_id,
                                       log(value + 1), log(true_value + 1)),
        mae = abs(median(value) - unique(true_value)),
        width = sort(value, partial = 17)[17] - sort(value, partial = 7)[7]
      ) %>%
      mutate(model = models[i]) %>%
      left_join(unique(ILINet[c("region", "count_rate2")]),
                by = c("location_name" = "region"))
    
    wis <- wis %>%
      left_join(baseline_wis[c("reference_date", "location_name",
                               "horizon", "bwis", "blwis")],
                by = c("reference_date", "location_name", "horizon")) %>%
      mutate(relative_wis = wis/bwis,
             relative_lwis = lwis/blwis)
    
    ili_scores <- rbind(ili_scores, wis)
    
  }
  print(sub_date)
}

saveRDS(team_scores, "../hub_scores_no_ens.rds")
saveRDS(ili_scores, "../nl_scores.rds")
rbind(team_scores, ili_scores) %>% 
  filter(str_detect(model, "asg")) %>%
  filter(location_name != "Florida") %>% 
  group_by(model) %>% 
  summarise(m = mean(lwis, na.rm = TRUE),
            n = n()) %>% 
  arrange(m)


team_coverage <- data.frame()

for (sub_date in sub_dates) {
  for (i in 1:length(models)) {
    forecast_file <- paste(models[i], "/", sub_date, "-",
                           models[i], ".csv", sep = "")
    
    if (file.exists(forecast_file) == FALSE) next
    forecast <- read.csv(forecast_file) %>%
      mutate(value = ifelse(value < 0, 0, value)) %>%
      filter(output_type == "quantile") %>%
      mutate(location = as.character(location)) %>%
      filter(output_type_id %in% c(.025, .25, .75, .975)) %>%
      pivot_wider(names_from = output_type_id, values_from = value) %>%
      rename(c(lower95 = "0.025", lower50 = "0.25", upper50 = "0.75",
               upper95 = "0.975"))
    
    for_cover <- forecast %>%
      left_join(target_data,
                by = c("location", "target_end_date" = "date")) %>%
      mutate(true_value = as.numeric(true_value))
    
    cover <- for_cover %>%
      mutate(cover50 = between(true_value, lower50, upper50),
             cover95 = between(true_value, lower95, upper95)) %>%
      mutate(model = models[i])
    
    team_coverage <- rbind(team_coverage, cover)
    
    
  }
}

