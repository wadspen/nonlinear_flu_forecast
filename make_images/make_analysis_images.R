setwd("~/flu_research/nonlinear_flu_forecast/")
source("../../flu_research/nonlinear_flu_forecast/forecasts_23/get_data_functions.R")
# setwd("./model-output/")
# source("../target-data/get_target_data.R")
library(stringr)
library(epidatr)
library(evalcast)
library(dplyr)
library(tidyr)
library(scoringutils)
both_flu <- comb_data(lag = 2)
date_weeks <- both_flu %>%
  select(date, season_week) %>% unique
ILINet <- get_ili_data()


setwd("~/flu_research/Prelim Content/models")

#season_week 10 is "2023-10-14"

target_data <- read.csv("../../../forecast-hub/FluSight-forecast-hub/target-data/target-hospital-admissions.csv") %>%
  mutate(true_value = value) %>%
  select(-value, -X)
# setwd("~/forecast-hub/FluSight-forecast-hub/model-output/")
folders <- list.files("./")
models <- folders[-which(folders == "FluSight-ensemble")]

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
                                   value, true_value),
    blwis = weighted_interval_score(output_type_id,
                                    log(value + 1), log(true_value + 1))
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
}


team_scores <- team_scores %>%
  mutate(traj = ifelse(str_detect(model, "asg"), "ASG",
                       ifelse(str_detect(model, "sir"), "SIR", " ")),
         disc = ifelse(str_detect(model, "disc"), "D", " "),
         dist = ifelse(str_detect(model, "lst"), "LST",
                       ifelse(str_detect(model, "log"), "LNORM",
                              ifelse(str_detect(model, "FluSight"), " ", 
                                     "NORM"))),
         sq = ifelse(str_detect(model, "sq"), "SQ", " ")) %>% 
  mutate(mod = ifelse(str_detect(model, "FluSight"), model,
                      gsub("- ", "", paste(traj, disc, dist, sq, sep = "-"))))

team_scores <- readRDS("../team_scores.rds")
library(scales)
c_trans <- function(a, b, breaks = b$breaks, format = b$format) {
  a <- as.trans(a)
  b <- as.trans(b)
  
  name <- paste(a$name, b$name, sep = "-")
  
  trans <- function(x) a$trans(b$trans(x))
  inv <- function(x) b$inverse(a$inverse(x))
  
  trans_new(name, trans, inverse = inv, breaks = breaks, format=format)
  
}



distq_labels <- 
  c(expression(LNORM), expression(LNORM^2), expression(LST), 
    expression(LST^2), expression(NORM), expression(NORM^2))



team_scores %>% 
  # filter(dist == "LST") %>%
  mutate(distq = ifelse(sq == "SQ", paste0(dist,2),
                        dist),
         disc = ifelse(disc == "D", "Disc","No Disc")) %>% 
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


team_scores %>% 
  # filter(dist == "LST") %>%
  mutate(distq = ifelse(sq == "SQ", paste0(dist,2),
                        dist),
         disc = ifelse(disc == "D", "Disc","No Disc")) %>% 
  # filter(location_name == "US") %>%
  filter(!(model %in% c("FluSight-ensemble", "FluSight-baseline"))) %>% 
  mutate(forecast_date = date(reference_date) + 7*as.numeric(horizon)) %>%
  ungroup() %>% 
  left_join(both_flu %>% 
              filter(season == 2023) %>% 
              select(season_week, date) %>% 
              mutate(date = date(date)) %>% 
              unique(),
            by = join_by(forecast_date == date), multiple = "all") %>% 
  group_by(season_week, location_name, disc, traj) %>% 
  summarise(mrwis = mean(lwis)) %>% 
  rowwise() %>% 
  mutate(minwis = mean(mrwis)) %>% 
  ggplot() +
  geom_tile(aes(y = season_week,
                x = reorder(location_name, minwis),
                fill = mrwis)) +
  scale_fill_gradient2(
    #low = "#5ab4ac", high = "#d8b365", mid = "white",
    low = "#5ab4ac", high = "black", mid = "white",
    # low = "white", high = "black",
    # midpoint = 1, limits = c(.6, 1.5),
    na.value = NA) +
  facet_grid(disc~traj) +
  coord_flip() +
  labs(fill = "LWIS") +
  ylab("Week") +
  xlab("Location") +
  
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(size=16),
        axis.title=element_text(size=20),
        strip.text = element_text(
          size = 13),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14),
        axis.text.y = element_text(size = 7),
        legend.position = "none")





team_scores %>% 
  # filter(dist == "LST") %>%
  mutate(distq = ifelse(sq == "SQ", paste0(dist,2),
                        dist),
         disc = ifelse(disc == "D", "Disc","No Disc")) %>% 
  # filter(location_name == "US") %>%
  filter(!(model %in% c("FluSight-ensemble", "FluSight-baseline"))) %>% 
  mutate(forecast_date = date(reference_date) + 7*as.numeric(horizon)) %>%
  ungroup() %>% 
  left_join(both_flu %>% 
              filter(season == 2023) %>% 
              select(season_week, date) %>% 
              mutate(date = date(date)) %>% 
              unique(),
            by = join_by(forecast_date == date), multiple = "all") %>% 
  group_by(season_week, location_name, dist) %>% 
  summarise(mrwis = mean(lwis)) %>% 
  rowwise() %>% 
  mutate(minwis = mean(mrwis)) %>% 
  ggplot() +
  geom_tile(aes(y = season_week,
                x = reorder(location_name, minwis),
                fill = mrwis)) +
  scale_fill_gradient2(
    #low = "#5ab4ac", high = "#d8b365", mid = "white",
    low = "#5ab4ac", high = "black", mid = "white",
    # low = "white", high = "black",
    # midpoint = 1, limits = c(.6, 1.5),
    na.value = NA) +
  facet_grid(.~dist) +
  coord_flip() +
  labs(fill = "LWIS") +
  ylab("Week") +
  xlab("Location") +
  
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(size=16),
        axis.title=element_text(size=20),
        strip.text = element_text(
          size = 13),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14),
        axis.text.y = element_text(size = 7),
        legend.position = "none")





team_scores_sum <- team_scores %>% 
  # filter(dist == "LST") %>%
  mutate(distq = ifelse(sq == "SQ", paste0(dist,2),
                        dist),
         disc = ifelse(disc == "D", "Disc","No Disc")) %>% 
  # filter(location_name == "US") %>%
  filter(!(model %in% c("FluSight-ensemble", "FluSight-baseline"))) %>% 
  mutate(forecast_date = date(reference_date) + 7*as.numeric(horizon)) %>%
  ungroup() %>% 
  left_join(both_flu %>% 
              filter(season == 2023) %>% 
              select(season_week, date) %>% 
              mutate(date = date(date)) %>% 
              unique(),
            by = join_by(forecast_date == date), multiple = "all") %>% 
  group_by(dist, traj, disc, sq) %>% 
  summarise(mrwis = mean(lwis)) %>% 
  rowwise() %>% 
  mutate(minwis = mrwis) %>% 
  arrange(traj, sq, disc, dist) %>% 
  select(traj, sq, disc, dist, mrwis)















team_scores %>% 
  filter(traj %in% c("ASG", "SIR")) %>% 
  mutate(traj = paste0(traj, disc)) %>% 
  mutate(distq = ifelse(sq == "SQ", paste0(dist,2),
                        dist),
         disc = ifelse(disc == "D", "Yes","No")) %>%
  
  # filter(location_name %in% "US") %>%
  filter(location_name %in% c("Alabama", "District of Columbia",
                              "Texas", "Ohio", "Iowa", "Montana")) %>%
  group_by(sq, traj, disc, dist, location_name) %>%
  ggplot() +
  geom_boxplot(aes(x = traj, y = lwis)) +
  facet_grid(location_name ~ distq, scales = "free_y") +
  # labs(colour = "Discrepancy") +
  xlab("ILI model") +
  ylab("LWIS") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(size=12, angle = 30),
        axis.text.y = element_text(size = 12, angle = 30),
        axis.title=element_text(size=20),
        strip.text = element_text(
          size = 11),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14))


r_team_scores <- team_scores %>% 
  filter(horizon != -1) %>% 
  mutate(model = mod, interval_score = lwis) %>% 
  select(location_name, horizon, model, interval_score) %>% 
  data.table::setDT() %>% 
  summarise_scores(by = c("model"),
                   relative_skill = TRUE,
                   baseline = "FluSight-baseline",
                   relative_skill_metric = "interval_score")

r_team_scores_loc <- team_scores %>% 
  filter(horizon != -1) %>% 
  mutate(model = mod, interval_score = lwis) %>% 
  select(location_name, horizon, model, interval_score) %>% 
  data.table::setDT() %>% 
  summarise_scores(by = c("model", "location_name"),
                   relative_skill = TRUE,
                   baseline = "FluSight-baseline",
                   relative_skill_metric = "interval_score")
  

r_team_scores_loc %>% 
  group_by(model) %>% 
  summarise(beat_base = sum(scaled_rel_skill < 1)) %>% 
  arrange(beat_base) %>% tail()


mods <- c("SIR-D-NORM-SQ", "SIR-D-LNORM-SQ", "SIR-D-LNORM", "ASG-D-NORM",
          "SIR-LNORM-SQ", "SIR-LNORM")

mod_scores <- r_team_scores_loc %>%
  as.data.frame() %>% 
  mutate(relative_lwis = scaled_rel_skill) %>% 
  filter(horizon != -1) %>%
  group_by(model, location_name) %>%
  filter(model %in% mods[1]) %>% #, "FluSight-ensemble")) %>% 
  arrange(relative_lwis)

mod_scores$location_name <- factor(mod_scores$location_name, 
                                   levels = mod_scores$location_name)


mod_scores %>% 
  ggplot() +
  geom_vline(xintercept = 1, colour = "darkgrey", size = 1.5) +
  geom_point(aes(y = location_name, x = relative_lwis), size = 2.4) +
  ylab("Location") + 
  xlab("RWIS") +
  theme_bw() + 
  theme(axis.text.x=element_text(size=15),
        axis.text.y = element_text(size = 7),
        axis.title=element_text(size=23),
        legend.position = "none")


sub_dates
locations <- unique(team_scores$location_name)
r_scores_week <- data.frame()
for (i in 15:length(sub_dates)) {
  for (j in 1:length(locations)) {
    r_team <- team_scores %>%
      mutate(forecast_date = date(reference_date) + 7*as.numeric(horizon)) %>% 
      filter(horizon != -1, forecast_date == sub_dates[i],
             location_name == locations[j]) %>% 
      mutate(model = mod, interval_score = lwis) %>% 
      select(location_name, horizon, model, interval_score) %>% 
      data.table::setDT() %>% 
      summarise_scores(by = c("model"),
                       relative_skill = TRUE,
                       baseline = "FluSight-baseline",
                       relative_skill_metric = "interval_score") %>% 
      as.data.frame()
    
    r_team$reference_date <- sub_dates[i]
    r_team$location_name <- locations[j]
    
    r_scores_week <- rbind(r_scores_week, r_team)
    print(c(i,j))
  
  }
}





r_scores_date <- readRDS("./r_scores_date.rds")
r_scores_date <- r_scores_date %>% 
  mutate(traj = ifelse(str_detect(model, "ASG"), "ASG",
                       ifelse(str_detect(model, "SIR"), "SIR", "")),
         disc = ifelse(str_detect(model, "-D-"), "Disc", "No Disc"),
         sq = ifelse(str_detect(model, "-SQ"), "SQ", ""),
         dist = ifelse(str_detect(model, "-NORM"), "NORM",
                       ifelse(str_detect(model, "-LNORM"), "LNORM",
                              ifelse(str_detect(model, "-LST"), "LST", ""))),
         relative_wis = scaled_rel_skill)





r_scores_date %>% 
  filter(!(model %in% c("FluSight-ensemble", "FluSight-baseline"))) %>% 
  mutate(forecast_date = date(reference_date)) %>% 
  ungroup() %>% 
  left_join(both_flu %>% 
              filter(season == 2023) %>% 
              select(season_week, date) %>% 
              mutate(date = date(date)) %>% 
              unique(),
            by = join_by(forecast_date == date), multiple = "all") %>% 
  group_by(location_name, season_week, dist) %>% 
  summarise(mrwis = mean(relative_wis)) %>% 
  rowwise() %>% 
  mutate(mrwis = min(2, mrwis)) %>% 
  group_by(location_name) %>% 
  mutate(minwis = mean(mrwis)) %>% 
  ggplot() +
  geom_tile(aes(y = season_week, x = reorder(location_name, minwis),
                fill = mrwis)) +
  scale_fill_gradient2(low = "#5ab4ac", high = "#d8b365", mid = "white",
                       midpoint = 1,
                       na.value = NA) +
  # scale_y_reverse()
  # scale_y_continuous(transform = "reverse") +
  facet_wrap(~dist) +
  coord_flip() + 
  labs(fill = "RLWIS") +
  ylab("Week") +
  xlab("Location") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(size=15),
        axis.text.y = element_text(size = 7),
        axis.title=element_text(size=23),
        strip.text = element_text(
          size = 12),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14))


r_scores_date %>% 
  filter(!(model %in% c("FluSight-ensemble", "FluSight-baseline"))) %>% 
  mutate(forecast_date = date(reference_date)) %>% 
  ungroup() %>% 
  left_join(both_flu %>% 
              filter(season == 2023) %>% 
              select(season_week, date) %>% 
              mutate(date = date(date)) %>% 
              unique(),
            by = join_by(forecast_date == date), multiple = "all") %>% 
  group_by(location_name, forecast_date, season_week, traj, disc, dist, sq) %>% 
  summarise(mrwis = mean(relative_wis)) %>% 
  mutate(mrwis = min(mrwis, 2)) %>% 
  # filter(between(season_week, 19, 22)) %>%
  ggplot() +
  geom_tile(aes(y = season_week, x = reorder(location_name, mrwis),
                fill = mrwis)) +
  scale_fill_gradient2(low = "#5ab4ac", high = "#d8b365", mid = "white",
                       midpoint = 1, 
                       na.value = NA) +
  
  labs(fill = "RLWIS") +
  ylab("Week") +
  xlab("Location") +
  # scale_y_reverse()
  coord_flip() +
  # scale_y_continuous(transform = "reverse") +
  facet_grid(disc~traj) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x=element_text(size=15),
        axis.text.y = element_text(size = 5),
        axis.title=element_text(size=23),
        strip.text = element_text(
          size = 12),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14))




r_scores_sum <- r_scores_date %>% 
  # filter(dist == "LST") %>%
  mutate(distq = ifelse(sq == "SQ", paste0(dist,2),
                        dist),
         disc = ifelse(disc == "D", "Disc","No Disc")) %>% 
  # filter(location_name == "US") %>%
  filter(!(model %in% c("FluSight-ensemble", "FluSight-baseline"))) %>% 
  # mutate(forecast_date = date(reference_date) + 7*as.numeric(horizon)) %>%
  mutate(reference_date = date(reference_date)) %>% 
  ungroup() %>% 
  left_join(both_flu %>% 
              filter(season == 2023) %>% 
              select(season_week, date) %>% 
              mutate(date = date(date)) %>% 
              unique(),
            by = join_by(reference_date == date), multiple = "all") %>% 
  group_by(dist, traj, disc, sq) %>% 
  summarise(mrwis = mean(relative_wis)) %>% 
  rowwise() %>% 
  mutate(minwis = mrwis) %>% 
  arrange(traj, sq, disc, dist) %>% 
  select(traj, sq, disc, dist, mrwis)


