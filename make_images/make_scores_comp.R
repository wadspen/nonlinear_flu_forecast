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



team_coverage <- data.frame()

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
      mutate(value = ifelse(value < 0, 0, value)) %>%
      filter(output_type == "quantile") %>%
      mutate(location = as.character(location)) %>%
      # filter(output_type_id %in% c(.025, .25, .75, .975)) %>%
      pivot_wider(names_from = output_type_id, values_from = value) %>%
      rename(c(lower98 = "0.01",
               lower95 = "0.025", 
               lower90 = "0.05",
               lower80 = "0.1",
               lower70 = "0.15",
               lower60 = "0.2",
               lower50 = "0.25",
               lower40 = "0.3",
               lower30 = "0.35",
               lower20 = "0.4",
               lower10 = "0.45",
               upper10 = "0.55",
               upper20 = "0.6",
               upper30 = "0.65",
               upper40 = "0.7",
               upper50 = "0.75",
               upper60 = "0.8",
               upper70 = "0.85",
               upper80 = "0.9",
               upper90 = "0.95",
               upper95 = "0.975",
               upper98 = "0.99"))
    
    for_cover <- forecast %>%
      left_join(target_data,
                by = c("location", "target_end_date" = "date")) %>%
      mutate(true_value = as.numeric(true_value))
    
    cover <- for_cover %>%
      mutate(cover10 = between(true_value, lower10, upper10),
             cover20 = between(true_value, lower20, upper20),
             cover30 = between(true_value, lower30, upper30),
             cover40 = between(true_value, lower40, upper40),
             cover50 = between(true_value, lower50, upper50),
             cover60 = between(true_value, lower60, upper60),
             cover70 = between(true_value, lower70, upper70),
             cover80 = between(true_value, lower80, upper80),
             cover90 = between(true_value, lower90, upper90),
             cover95 = between(true_value, lower95, upper95),
             cover98 = between(true_value, lower98, upper98)) %>%
      mutate(model = models[i])
    
    team_coverage <- rbind(team_coverage, cover)
    
    
  }
}


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


team_coverage %>% 
  dplyr::select(reference_date, horizon, location, model, contains("cover")) %>% 
  pivot_longer(contains("cover"), names_to = "coverage", 
               values_to = "covered") %>% 
  filter(horizon >= 0 & !is.na(covered)) %>% 
  mutate(coverage = .01*as.numeric(str_extract(coverage, "\\d+"))) %>% 
  group_by(model, coverage) %>% 
  summarise(pcover = mean(covered)) %>% 
  filter(model == "FluSight-baseline") %>% 
  ggplot() +
  geom_line(aes(x = coverage, y = pcover, colour = model)) +
  stat_function(fun = function(x) x + 0, color = "blue", size = 1.2)
  

team_scores <- readRDS("~/flu_research/nonlinear_flu_forecast/flu_forecast_23/hub_scores_no_ens.rds")

setwd("../../../flu_research/nonlinear_flu_forecast/flu_forecast_23/fin_hosp_forecasts/")

folders <- list.files("./")
folders <- folders[!str_detect(folders, "full") &
                   !str_detect(folders, "_disc_")]
# folders <- folders[!str_detect(folders, "disc2") |
#                      str_detect(folders, "disc2_nm")]

folders <- folders[(str_detect(folders, "asg_nm")
                   | str_detect(folders, "asg_disc2_nm")
                   | str_detect(folders, "sir_hosp")
                   | str_detect(folders, "sir_disc2"))
                   & !str_detect(folders, "logi")
                   & !str_detect(folders, "tr")]

# folders <- folders[str_detect(folders, "asg_disc2_ho") & !str_detect(folders, "tr")]
# folders <- folders[str_detect(folders, "asg_disc2_nm") & !str_detect(folders, "tr")]
models <- folders
# models <- "asg_disc2_nm_hosp_lst_ar1"
# most_coverage <- ili_coverage
ili_coverage <- data.frame()

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
      mutate(value = ifelse(value < 0, 0, value)) %>%
      filter(output_type == "quantile") %>%
      mutate(location = as.character(location)) %>%
      # filter(output_type_id %in% c(.025, .25, .75, .975)) %>%
      pivot_wider(names_from = output_type_id, values_from = value) %>%
      rename(c(lower98 = "0.01",
               lower95 = "0.025", 
               lower90 = "0.05",
               lower80 = "0.1",
               lower70 = "0.15",
               lower60 = "0.2",
               lower50 = "0.25",
               lower40 = "0.3",
               lower30 = "0.35",
               lower20 = "0.4",
               lower10 = "0.45",
               upper10 = "0.55",
               upper20 = "0.6",
               upper30 = "0.65",
               upper40 = "0.7",
               upper50 = "0.75",
               upper60 = "0.8",
               upper70 = "0.85",
               upper80 = "0.9",
               upper90 = "0.95",
               upper95 = "0.975",
               upper98 = "0.99"))
    
    for_cover <- forecast %>%
      left_join(target_data,
                by = c("location", "target_end_date" = "date")) %>%
      mutate(true_value = as.numeric(true_value))
    
    cover <- for_cover %>%
      mutate(cover10 = between(true_value, lower10, upper10),
             cover20 = between(true_value, lower20, upper20),
             cover30 = between(true_value, lower30, upper30),
             cover40 = between(true_value, lower40, upper40),
             cover50 = between(true_value, lower50, upper50),
             cover60 = between(true_value, lower60, upper60),
             cover70 = between(true_value, lower70, upper70),
             cover80 = between(true_value, lower80, upper80),
             cover90 = between(true_value, lower90, upper90),
             cover95 = between(true_value, lower95, upper95),
             cover98 = between(true_value, lower98, upper98)) %>%
      mutate(model = models[i])
    
    ili_coverage <- rbind(ili_coverage, cover)
    
    
  }
}

# saveRDS(most_scores, "~/flu_research/nonlinear_flu_forecast/flu_forecast_23/nl_scores.rds")
# saveRDS(most_coverage, "~/flu_research/nonlinear_flu_forecast/flu_forecast_23/nl_cover.rds")
# most_scores <- ili_scores
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
    if (sum(is.na(forecast$value)) > 0) {next}
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

# saveRDS(ili_scores, "~/flu_research/nonlinear_flu_forecast/flu_forecast_23/nl_scores.rds")
rbind(ili_scores, team_scores) %>%
# ili_scores %>% 
  filter(location_name == "US") %>% 
  group_by(model) %>% 
  filter(horizon != -1) %>% 
  unique() %>% 
  filter(!str_detect(model, "tr_")) %>% #, location_name != "Florida") %>% 
  filter(!is.na(wis)) %>% 
  summarise(m = mean(lwis, na.rm = TRUE),
            n = n()) %>%
  # filter(str_detect(model, "asg_nm")) %>%
  arrange(m) %>% 
  as.data.frame()


rbind(ili_scores, team_scores) %>% 
  ggplot() +
  geom_boxplot(aes(x = model, y = wis)) +
  scale_y_log10()

asgbase <- scores %>% 
  filter(model == "asg_disc2_nm_hosp_sq_ar1") %>% 
  mutate(all = paste(reference_date, location_name, horizon))

base <- ili_scores %>% 
  filter(model == "asg_disc2_hosp_sq_ar1") %>% 
  mutate(all = paste(reference_date, location_name, horizon))


scores %>% 
  mutate(all = paste(reference_date, location_name, horizon)) %>% 
  filter(all %in% base$all) %>% 
  group_by(model) %>% 
  summarise(m = mean(wis, na.rm = TRUE)) %>% 
  arrange(m)


ili_scores %>% 
  filter(str_detect(model, "asg_disc2_nm")) %>% 
  ggplot() +
  geom_boxplot(aes(x = lwis)) + 
  facet_wrap(~model, scales = "free")

ili_scores %>% 
  filter(!str_detect(model, "tr_")) %>%
  group_by(model) %>% 
  mutate(max = wis == max(wis, na.rm = TRUE)) %>% 
  filter(max != TRUE) %>% 
  group_by(model) %>% 
  summarise(m = mean(wis)) %>% 
  arrange(m) %>% View()

ili_scores <- ili_scores %>% 
  mutate(all = paste(location_name, reference_date, horizon, sep = "_"))

asg_scores <- ili_scores %>% 
  filter(model == "asg_disc2_nm_hosp_sq_ar1")

asg_lst_scores <- ili_scores %>% 
  filter(model == "asg_disc2_nm_hosp_lst_ar1")


sscores <- team_scores %>% 
  filter(model == "SigSci-CREG") %>% 
  mutate(all = paste(reference_date, location_name, horizon, sep = "_"))
# team_scores <- readRDS("./hub_scores_no_ens.rds")
# saveRDS(team_scores, "../hub_scores_no_ens.rds")
# saveRDS(ili_scores, "../nl_scores.rds")

rbind(team_scores, ili_scores) %>% # %>% filter(model == "asg_disc2_nm_hosp_sq_ar1")) %>% 
  # filter(str_detect(model, "asg")) %>%
  # filter(location_name != "Florida") %>% 
  # mutate(all = paste(reference_date, location_name, horizon, sep = "_")) %>%
  # filter(all %in% unique(sscores$all)) %>%
  group_by(model) %>% 
  summarise(m = mean(lwis, na.rm = TRUE),
            n = n()) %>% 
  arrange(m)

ili_scores %>% 
  filter(str_detect(model, "asg") #& str_detect(model, "nm")
         & !str_detect(model, "logi") & str_detect(model, "lst")) %>% 
  group_by(model) %>% 
  summarise(n = n()) %>% 
  arrange(model)


ili_scores %>% 
  # filter((str_detect(model, "sir_hosp") | str_detect(model, "sir_disc2"))
  #        & !str_detect(model, "logi")) %>% 
  filter(str_detect(model, "sir")) %>% 
  group_by(model) %>% 
  summarise(n = n()) %>% 
  arrange(model)

ili_coverage <- data.frame()

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
    
    ili_coverage <- rbind(ili_coverage, cover)
    
    
  }
}




team_scores <- readRDS("./flu_forecast_23/both_flu.rds")

team_scores %>% 
  filter(location_name == "US") %>%
  filter(!(model %in% c("FluSight-ensemble", "FluSight-baseline"))) %>% 
  mutate(forecast_date = date(reference_date) + 7*as.numeric(horizon)) %>%
  ungroup() %>% 
  left_join(both_flu %>% 
              filter(season == 2023) %>% 
              select(season_week, date) %>% 
              mutate(date = date(date)) %>% 
              unique(),
            by = join_by(forecast_date == date), multiple = "all") %>% 
  group_by(season_week, disc, traj, distq) %>% 
  summarise(mrwis = mean(lwis)) %>% 
  rowwise() %>% 
  mutate(minwis = mean(mrwis)) %>% 
  ggplot() +
  geom_hline(yintercept = 22, size = 8, colour = "red") +
  geom_tile(aes(y = season_week,
                x = disc,
                fill = mrwis)) +
  scale_fill_gradient2(
    #low = "#5ab4ac", high = "#d8b365", mid = "white",
    low = "#5ab4ac", high = "black", mid = "white",
    # low = "white", high = "black",
    # midpoint = 1, limits = c(.6, 1.5),
    na.value = NA) +
  facet_grid(traj~distq,
             labeller = labeller(distq = distq_labels)) +
  coord_flip() +
  labs(fill = "LWIS") +
  ylab("Week") +
  xlab("") +
  
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(size=16),
        axis.text.y = element_text(size = 14, angle = 90, hjust = .45),
        axis.title=element_text(size=20),
        strip.text = element_text(
          size = 13),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14),
        legend.position = "none")

