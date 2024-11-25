setwd("~/flu_research/Prelim Content/")
source("../../flu_research/nonlinear_flu_forecast/forecasts_23/get_data_functions.R")
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


setwd("~/flu_research/Prelim Content/models")

#season_week 10 is "2023-10-14"

target_data <- read.csv("./target-data/target-hospital-admissions.csv") %>%
  mutate(true_value = value) %>%
  select(-value, -X)
# setwd("~/forecast-hub/FluSight-forecast-hub/model-output/")
folders <- list.files("./")
models <- folders[folders != "old_models" & folders != "target-data"]
sub_dates <- substr(list.files("./FluSight-baseline"), 1, 10)
sub_date <- sub_dates[1]
base <- "FluSight-baseline"

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
                                   value, true_value)
  ) %>%
  mutate(model = base)

team_scores <- data.frame()
for (sub_date in sub_dates) {
  for (i in 1:length(models)) {
    forecast_file <- paste(models[i], "/", sub_date, "-",
                           models[i], ".csv", sep = "")
    
    if (file.exists(forecast_file) == FALSE) next
    forecast <- read.csv(forecast_file) %>%
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
    
    for_wis <- for_wis %>% 
      filter(is.na(location_name) == FALSE)
    if (nrow(for_wis) == 0) next
    wis <- for_wis %>%
      filter(is.na(location_name) == FALSE) %>%
      group_by(reference_date, location_name, horizon) %>%
      reframe(
        wis = weighted_interval_score(output_type_id,
                                      value, true_value),
        mae = abs(median(value) - unique(true_value)),
        width = sort(value, partial = 17)[17] - sort(value, partial = 7)[7]
      ) %>%
      mutate(model = models[i]) %>%
      left_join(unique(ILINet[c("region", "count_rate1")]),
                by = c("location_name" = "region"))
    
    wis <- wis %>%
      left_join(baseline_wis[c("reference_date", "location_name",
                               "horizon", "bwis")],
                by = c("reference_date", "location_name", "horizon")) %>%
      mutate(relative_wis = wis/bwis)
    
    team_scores <- rbind(team_scores, wis)
    
  }
}


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

team_scores %>%
  filter(location_name == "Michigan", horizon != -1,
         model != "FluSight-baseline") %>%
  mutate(horizon = paste(horizon + 1, " week(s) ahead", sep = "")) %>%
  ggplot(aes(x = reference_date, y = relative_wis, group = model,
             colour = model)) +
  geom_line() +
  facet_wrap(~horizon) +
  
  xlab("Date") +
  ylab("Relative WIS") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        axis.text=element_text(size=12),
        axis.title=element_text(size=15)) #+
scale_colour_discrete(labels=c("ASG + disc", "ASG + disc2", "ASG",
                               "CDC Baseline", "SIR + disc", "SIR + disc2",
                               "SIR"))


all_scores <- team_scores %>% 
  filter(horizon != -1, reference_date < "2024-02-31") %>%
  # filter(horizon == 3) %>% 
  group_by(model) %>%
  summarise(
    med = median(relative_wis),
    mean = mean(relative_wis),
    m = mean(relative_wis < 1),
    n = n(),
    mae = mean(mae),
    s = sum(relative_wis < 1)
  ) %>% arrange(med)


team_scores <- team_scores %>% 
  filter(!(model %in% c("ISU_NiemiLab-SIR", "ISU_NiemiLab-ENS",
                         "ISU_NiemiLab-NLH")))


weekly_ranks <- team_scores %>% 
  group_by(reference_date, model, location_name) %>% 
  summarise(
    med = mean(relative_wis),
    n = n()
  ) %>% 
  arrange(reference_date, location_name, med) %>% 
  group_by(reference_date, location_name) %>% 
  mutate(rank = rank(med)) %>% 
  mutate(max_rank = max(rank)) %>% 
  mutate(rank_perc = (max_rank - rank)/max_rank)

asg_rank <- weekly_ranks %>% 
  filter(model == "asg_disc_sq_ar1_lag") %>% 
  group_by(location_name) %>% 
  summarise(
    med_perc = mean(rank_perc)
  ) %>% arrange(desc(med_perc))

location_names <- asg_rank$location_name

weekly_ranks %>% 
  filter(model %in% c("asg_sq_ar1_lag", "asg_disc_sq_ar1_lag")) %>% 
  group_by(location_name, model) %>% 
  summarise(
    med_perc = mean(rank_perc)
  ) %>% arrange(desc(med_perc)) %>% 
  mutate(location_name = factor(location_name, 
                                levels = as.character(location_names))) %>%
  # mutate(location_name = reorder.factor(location_name, location_name)) %>% 
  # mutate(location_name = relevel(location_name, location_names)) %>% 
  ggplot() +
  geom_point(aes(x = location_name, y = med_perc,
                 group = model, colour = model), size = 2) +
  ylab("RWIS Percentile") +
  xlab("State") +
  ylim(c(0,1)) +
  labs(colour = "Model") +
  scale_colour_discrete(labels=c("ASG + disc", "ASG")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + 
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=15),
        legend.position = c(.8, .8))

rankings %>% 
  group_by(location_name) %>% 
  summarise(
    mean_perc = mean(rank_perc)
  )

state_scores <- team_scores %>% 
  filter(horizon != -1) %>%
  group_by(model, location_name) %>%
  summarise(
    med = median(relative_wis),
    mean = mean(relative_wis),
    m = mean(relative_wis < 1),
    n = n(),
    mae = mean(mae),
    s = sum(relative_wis < 1)
  ) %>% arrange(med) 

state_scores %>% 
  filter(med <= 1) %>%
  group_by(model) %>% 
  summarise(n = n())



team_scores %>%
  filter(horizon != -1) %>%
  group_by(model, location_name) %>%
  summarise(
    m = mean(relative_wis)
  ) %>%
  filter(m < 1) %>%
  group_by(model) %>%
  summarise(
    n = n()
  )


team_scores %>%
  filter(horizon != -1) %>%
  group_by(model) %>%
  summarise(
    m = max(relative_wis)
  )


team_scores %>%
  filter(horizon != -1) %>%
  ggplot(aes(x = relative_wis)) +
  geom_histogram() +
  facet_wrap(~model, scale = "free")



######show this one
team_scores %>%
  filter(model != "CADPH-FluCAT_Ensemble", horizon != -1) %>%
  # filter(model %in% c("CMU-TimeSeries", "asg_disc_rate_sq_rw_ncor")) %>%
  group_by(model, reference_date, location_name, horizon) %>%
  summarise(
    mean_wis = median(relative_wis, na.rm = TRUE)
  ) %>%
  arrange(location_name, mean_wis, reference_date) %>%
  group_by(reference_date, location_name, horizon) %>%
  mutate(rank = rank(mean_wis)) %>%
  group_by(model) %>%
  summarise(
    m_rank = median(rank)
  ) %>% arrange(m_rank)

team_scores %>%
  filter(model != "CADPH-FluCAT_Ensemble", horizon != -1) %>%
  # filter(model %in% c("CMU-TimeSeries", "asg_disc_rate_sq_rw_ncor")) %>%
  group_by(model, reference_date, location_name, horizon) %>%
  summarise(
    mean_wis = median(relative_wis, na.rm = TRUE)
  ) %>%
  arrange(location_name, mean_wis, reference_date) %>%
  group_by(reference_date, location_name, horizon) %>%
  mutate(rank = rank(mean_wis)) %>%
  filter(model %in% c("sir_disc_rate_sq_rw", "asg_disc_rate_sq_rw3_ncor")) %>%
  ggplot(aes(x = rank, fill = model)) +
  geom_density(alpha = .4)








team_coverage %>%
  group_by(model) %>%
  # filter(!(location_name %in% c("Alaska", "Nevada", "Florida", "New Hampshire",
  #                               "New Hampshire", "Utah", "Connecticut"))) %>%
  # filter(location_name == 'US') %>%
  summarise(
    mc50 = mean(cover50, na.rm = TRUE),
    mc95 = mean(cover95, na.rm = TRUE)
  ) %>% arrange(mc95)
  data.frame













##########################################
##############plots n'stuff###############
##########################################


for (sub_date in sub_dates) {
  for (i in 1:length(models)) {
    forecast_file <- paste(models[i], "/", sub_date, "-",
                           models[i], ".csv", sep = "")
    
    if (file.exists(forecast_file) == FALSE) next
    forecast <- read.csv(forecast_file) %>%
      filter(output_type == "quantile") %>%
      mutate(location = as.character(location))# %>% select(reference_date, horizon, target, target_end_date, location, output_type, output_type_id, value)
    
    
    date_forecast <- forecast %>%
      # mutate(value = (exp(value) - 1)) %>%
      filter(output_type_id %in% c("0.01", "0.025", "0.25",
                                   "0.75", "0.975", "0.99")) %>%
      mutate(output_type_id = ifelse(output_type_id == "0.01", "lower_99",
                                     ifelse(output_type_id == "0.025", "lower_95",
                                            ifelse(output_type_id == "0.25", "lower_50",
                                                   ifelse(output_type_id == "0.75", "upper_50",
                                                          ifelse(output_type_id == "0.975", "upper_95", "upper_99")))))) %>%
      unique() %>%
      pivot_wider(names_from = output_type_id) %>%
      left_join(week_dates, by = c("target_end_date" = "date"))
    
    
    date_forecast %>%
      filter(location == "US") %>%
      ggplot() +
      geom_ribbon(aes(x = season_week, ymin = lower_99,
                      ymax = upper_99, fill = "lightpink")) +
      geom_ribbon(aes(x = season_week, ymin = lower_95,
                      ymax = upper_95), fill = "pink") +
      geom_ribbon(aes(x = season_week, ymin = lower_50,
                      ymax = upper_50), fill = "red") +
      facet_wrap(~reference_date) +
      geom_point(both_flu %>%
                   filter(location == "US", season == 2023),
                 mapping = aes(x = season_week, y = value))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    minor in library science, geneaology, and german
    really did try. really struggled. practiced words with grandpa and did the best then
    Hinckleys were big on kids going to school. Bud wanted to be a doctor or engineer,
    but couldn't concentrate after being a prisoner of war. Maybe meant a lot to him
when Martin became an engineer. Bud and Ruth so wrapped up in emotional struggles
they couldn't see much passed it. Vision of education and career not really passed on
    to the kids because of it.
    Mom almost had a hard time with president Hinkley pushing education because she
    felt a woman's responsibility was in the home. 'Better to see a quote from
    President Hinckley saying marrying in the temple is the most important thing but
    a good education would facilitate it. Mom always included that quote if teaching
    about getting an education. It was to her a bridge.
    
    
    
    
    
    