source("../flu_forecast_23/get_data_functions_new.R")
library(datasets)
library(ggtext)
state.name <- c(state.name, "Puerto Rico", "District of Columbia", "US")
state.abb <- c(state.abb, "PR", "DC", "US")
states <- data.frame(loc_abb = state.abb, location_name = state.name)


team_scores <- readRDS("../flu_forecast_23/hub_scores_no_ens.rds") %>% 
  left_join(states, by = "location_name") %>% 
  mutate(model = ifelse(model == "FluSight-baseline", 
                        # "<b><i>FluSight-baseline</i></b>", 
                        "<span style='font-size:12pt'><b><i>FluSight-baseline</i></b></span>",
                        model)) %>% 
  unique()
ili_scores <- readRDS("../flu_forecast_23/nl_scores.rds") %>% 
  left_join(states, by = "location_name") %>% 
  unique()
  
team_scores <- rbind(team_scores, ili_scores %>% 
                       filter(model == "sir_disc2_hosp_sq_ar1")) %>% 
  mutate(model = ifelse(model == "sir_disc2_hosp_sq_ar1", 
                        "<span style='font-size:12pt'><b><i>SIRD_NORM2</i></b></span>",
                        model))






team_scores <- rbind(team_scores, 
                     ili_scores %>% 
                       filter(model == "asg_disc2_nm_hosp_sq_ar1")) %>% 
  unique()

best_nl <- ili_scores %>% 
  filter(model == "asg_disc2_nm_hosp_sq_ar1") %>% 
  mutate(nl_wis = wis, nl_lwis = lwis, nl_mae = mae) %>% 
  dplyr::select(reference_date, loc_abb, horizon, nl_wis, nl_lwis, nl_mae)

team_scores %>% 
  left_join(best_nl, by = c("reference_date", "loc_abb", "horizon")) %>% 
  filter(model %in% c( "FluSight-baseline")) %>% 
  group_by(model, reference_date, horizon) %>% 
  summarise(mwis = mean(wis, na.rm = TRUE),
            nl_mwis = mean(nl_wis, na.rm = TRUE)) %>% 
  # filter(model == "FluSight-baseline") %>%
  ggplot() +
  geom_point(aes(x = nl_mwis, y = mwis, colour = model)) +
  stat_function(fun = function(x) x + 0, color = "blue3") +
  scale_y_log10() +
  scale_x_log10() +
  theme(legend.position = "none")


team_scores %>% 
  left_join(best_nl, by = c("reference_date", "loc_abb", "horizon")) %>% 
  group_by(model) %>% 
  summarise(m = mean(lwis < nl_lwis, na.rm = TRUE)) %>% 
  arrange(desc(m))

team_diffs <- team_scores %>% 
  filter(model != "asg_disc2_nm_hosp_sq_ar1") %>% 
  left_join(best_nl, by = c("reference_date", "loc_abb", "horizon")) %>% 
  group_by(model, loc_abb) %>% 
  # filter(model != "GH-model") %>% 
  summarise(mlwis = mean(lwis, na.rm = TRUE),
            nl_mlwis = mean(nl_lwis, na.rm = TRUE),
            mwis = mean(wis, na.rm = TRUE),
            nl_mwis = mean(nl_wis, na.rm = TRUE)) %>% 
  arrange(mwis) %>% 
  mutate(sldiff = nl_mlwis/mlwis,
         sdiff = nl_mwis/mwis)  

scores_ord <- team_scores %>% 
  filter(model != "asg_disc2_nm_hosp_sq_ar1") %>%
  left_join(best_nl, by = c("reference_date", "loc_abb", "horizon")) %>% 
  group_by(model) %>% 
  # filter(model != "GH-model") %>% 
  summarise(mlwis = mean(lwis, na.rm = TRUE),
            nl_mlwis = mean(nl_lwis, na.rm = TRUE),
            mwis = mean(wis, na.rm = TRUE),
            nl_mwis = mean(nl_wis, na.rm = TRUE)) %>% 
  arrange(mwis) %>% 
  mutate(sldiff = nl_mlwis/mlwis,
         sdiff = nl_mwis/mwis) %>% 
  arrange(sldiff)

scores_ord$model

state_ord <- team_scores %>% 
  left_join(best_nl, by = c("reference_date", "loc_abb", "horizon")) %>% 
  group_by(loc_abb) %>% 
  # filter(model != "GH-model") %>% 
  summarise(mlwis = mean(lwis, na.rm = TRUE),
            nl_mlwis = mean(nl_lwis, na.rm = TRUE),
            mwis = mean(wis, na.rm = TRUE),
            nl_mwis = mean(nl_wis, na.rm = TRUE)) %>% 
  arrange(mwis) %>% 
  mutate(sldiff = nl_mlwis/mlwis,
         sdiff = nl_mwis/mwis) %>% 
  arrange(sldiff)

state_ldiff <- team_diffs %>% 
  filter(!(model %in% c("UNC_IDD-InfluPaint", "NU_UCSD-GLEAM_AI_FLUH",
                        "SigSci-CREG", "PSI-PROF_beta", "NIH-Flu_ARIMA",
                        "GH-model", "UGA_flucast-OKeeffe"))) %>% 
  filter(model != "asg_disc2_nm_hosp_sq_ar1") %>%
  mutate(model = factor(model, levels = scores_ord$model),
         loc_abb = factor(loc_abb, levels = state_ord$loc_abb)) %>% 
  mutate(diff_lim = ifelse(sdiff < -1, 
                           #-max(team_diffs$sdiff), 
                           -1, 
                           #-max(team_diffs$sdiff), 
                           sdiff)) %>% 
  ggplot() +
  geom_tile(aes(y = model,
                x = loc_abb,
                fill = sldiff)) +
  scale_fill_gradient2(
    low = "blue3", high = "red1", mid = "white",
    midpoint = 1) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  xlab("Location") +
  ylab("") +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey50", color = NA),
        # legend.position = "none",
        axis.text.x = element_text(size = 11, angle = 90),
        legend.title = element_blank(),
        # legend.position = "none",
        # axis.text.y = element_text(size = 10, angle = 20),
        axis.title = element_text(size = 16),
        axis.text.y = element_markdown())



state_diff <- team_diffs %>% 
  filter(!(model %in% c("UNC_IDD-InfluPaint", "NU_UCSD-GLEAM_AI_FLUH",
                        "SigSci-CREG", "PSI-PROF_beta", "NIH-Flu_ARIMA",
                        "GH-model", "UGA_flucast-OKeeffe"))) %>% 
  filter(model != "asg_disc2_nm_hosp_sq_ar1") %>%
  mutate(model = factor(model, levels = scores_ord$model),
         loc_abb = factor(loc_abb, levels = state_ord$loc_abb)) %>% 
  mutate(diff_lim = ifelse(sdiff < -1, 
                           #-max(team_diffs$sdiff), 
                           -1, 
                           #-max(team_diffs$sdiff), 
                           sdiff)) %>% 
  ggplot() +
  geom_tile(aes(y = model,
                x = loc_abb,
                fill = sdiff)) +
  scale_fill_gradient2(
    low = "blue3", high = "red1", mid = "white",
    midpoint = 1) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  xlab("Location") +
  ylab("") +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey50", color = NA),
        # legend.position = "none",
        axis.text.x = element_text(size = 11, angle = 90),
        legend.title = element_blank(),
        # legend.position = "none",
        # axis.text.y = element_text(size = 10, angle = 20),
        axis.title = element_text(size = 16),
        axis.text.y = element_markdown())














team_ddiffs <- team_scores %>% 
  filter(!(model %in% c("UNC_IDD-InfluPaint", "NU_UCSD-GLEAM_AI_FLUH",
                        "SigSci-CREG", "PSI-PROF_beta", "NIH-Flu_ARIMA",
                        "GH-model", "UGA_flucast-OKeeffe"))) %>% 
  filter(model != "asg_disc2_nm_hosp_sq_ar1") %>% 
  left_join(best_nl, by = c("reference_date", "loc_abb", "horizon")) %>% 
  group_by(model, reference_date) %>% 
  # filter(model != "GH-model") %>% 
  summarise(mlwis = mean(lwis, na.rm = TRUE),
            nl_mlwis = mean(nl_lwis, na.rm = TRUE),
            mwis = mean(wis, na.rm = TRUE),
            nl_mwis = mean(nl_wis, na.rm = TRUE)) %>% 
  arrange(mwis) %>% 
  mutate(sldiff = nl_mlwis/mlwis,
         sdiff = nl_mwis/mwis) %>% 
  filter(model != "asg_disc2_nm_hosp_sq_ar1") %>%
  mutate(model = factor(model, levels = scores_ord$model)) %>% 
  mutate(diff_lim = ifelse(sdiff < -1, 
                           #-max(team_diffs$sdiff), 
                           -1, 
                           #-max(team_diffs$sdiff), 
                           sdiff))

library(scales)
library(ggtext)
# labels <- levels(team_ddiffs$model)
# labels[labels == "FluSight-baseline"] <- "<b>FluSight-baseline</b>"
date_ldiff <- team_ddiffs %>%
  filter(!(model %in% c("UNC_IDD-InfluPaint", "NU_UCSD-GLEAM_AI_FLUH",
                        "SigSci-CREG", "PSI-PROF_beta", "NIH-Flu_ARIMA",
                        "GH-model", "UGA_flucast-OKeeffe"))) %>% 
  mutate(model = factor(model, levels = scores_ord$model)) %>%  
  ggplot() +
  geom_tile(aes(y = model,
                x = date(reference_date),
                fill = sldiff)) +
  scale_fill_gradient2(
    low = "blue3", high = "red1", mid = "white",
    midpoint = 1) +
  scale_x_date(expand = c(0,0),
               breaks = pretty_breaks(n = 6)) +
  scale_y_discrete(expand = c(0,0)) +
  ylab("") +
  xlab("Date") +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey50", color = NA),
        # legend.position = "none",
        axis.text.x = element_text(size = 11),
        legend.title = element_blank(),
        # legend.position = "none",
        # axis.text.y = element_text(size = 10, angle = 20),
        axis.title = element_text(size = 16),
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank()
        axis.text.y = element_markdown())


% UNC\_IDD-InfluPaint         & 107 & 0.35 & 16 & 0.24 & 0.043 & 78 & 0.58 & 0.95 \\ 
% NU\_UCSD-GLEAM\_AI\_FLUH    & 98 & 0.39 & 22 & 0.24 & 0.014 & 74 & 0.68 & 0.83 \\ 
% SigSci-CREG                 & 59 & 0.49 & 12 & 0.31 & 0.018 & 71 & 0.62 & 0.74 \\ 
% PSI-PROF\_beta              & 77 & 0.32 & 18 & 0.20 & 0.004 & 67 & 0.87 & 0.93 \\ 
% NIH-Flu\_ARIMA              & 230 & 0.34 & 28 & 0.23 & 0.001 & 20 & 0.57 & 0.88 \\ 
% GH-model                    & 161 & 1.52 & 28 & 1.47 & 0.373 & 17 & 0.18 & 0.27 \\ 
% UGA\_flucast-OKeeffe        & 27 & 0.59 & 8 & 0.43 & 0.003 & 10 & 0.60 & 0.58 \\ 

date_diff <- team_ddiffs %>%
  filter(!(model %in% c("UNC_IDD-InfluPaint", "NU_UCSD-GLEAM_AI_FLUH",
                        "SigSci-CREG", "PSI-PROF_beta", "NIH-Flu_ARIMA",
                        "GH-model", "UGA_flucast-OKeeffe"))) %>% 
  mutate(model = factor(model, levels = scores_ord$model)) %>%  
  ggplot() +
  geom_tile(aes(y = model,
                x = date(reference_date),
                fill = sdiff)) +
  scale_fill_gradient2(
    low = "blue3", high = "red1", mid = "white",
    midpoint = 1) +
  scale_x_date(expand = c(0,0),
               breaks = pretty_breaks(n = 6)) +
  scale_y_discrete(expand = c(0,0)) +
  ylab("") +
  xlab("Date") +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey50", color = NA),
        # legend.position = "none",
        axis.text.x = element_text(size = 11),
        legend.title = element_blank(),
        # legend.position = "none",
        # axis.text.y = element_text(size = 10, angle = 20),
        axis.title = element_text(size = 16),
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank()
        axis.text.y = element_markdown())




library(gridExtra)
library(grid)

grid.arrange(state_ldiff, date_ldiff, 
             nrow = 2, 
             left = textGrob("LWIS",
                             rot = 90, vjust = 1,
                             gp = gpar(fontsize = 16)))

grid.arrange(state_diff, date_diff, 
             nrow = 2, 
             left = textGrob("WIS",
                             rot = 90, vjust = 1,
                             gp = gpar(fontsize = 16)))

state_ldiff / state_diff
date_ldiff / date_diff
"lightred"

state_ldiff / date_ldiff
state_diff / date_diff

##################################
#############Coverage#############
##################################
# saveRDS(ili_coverage, "./nl_cover.rds")
team_coverage <- readRDS("../flu_forecast_23/hub_cover_no_ens.rds")
ili_coverage <- readRDS("../flu_forecast_23/nl_cover.rds")
# team_coverage <- rbind(team_coverage, ili_coverage %>% 
#                          filter(model == "sir_disc2_hosp_sq_ar1")) %>% 
#   mutate(model = ifelse(model == "sir_disc2_hosp_sq_ar1", "SIRD_NORM2",
#                         model))
pcover <- team_coverage %>% 
  filter(!(model %in% c("UNC_IDD-InfluPaint", "NU_UCSD-GLEAM_AI_FLUH",
                        "SigSci-CREG", "PSI-PROF_beta", "NIH-Flu_ARIMA",
                        "GH-model", "UGA_flucast-OKeeffe"))) %>%
  filter(!str_detect(model, "FluSight-baseline")) %>% 
  dplyr::select(reference_date, horizon, location, model, contains("cover")) %>% 
  pivot_longer(contains("cover"), names_to = "coverage", 
               values_to = "covered") %>% 
  filter(horizon >= 0 & !is.na(covered)) %>% 
  mutate(coverage = .01*as.numeric(str_extract(coverage, "\\d+"))) %>% 
  group_by(model, coverage) %>% 
  summarise(pcover = mean(covered)) %>% 
  # filter(model == "FluSight-baseline") %>% 
  ggplot(aes(x = coverage, y = pcover, line_type = model)) +
  geom_line(colour = "grey", size = .9) +
  geom_abline(intercept = 0, slope = 1, size = 1.2) +
  geom_line(data = ili_coverage %>%
              filter(model == "asg_disc2_nm_hosp_sq_ar1") %>%
              mutate(model = "ASGD_NORM2") %>% 
              dplyr::select(reference_date, horizon, location, model, contains("cover")) %>% 
              pivot_longer(contains("cover"), names_to = "coverage", 
               values_to = "covered") %>% 
              filter(horizon >= 0 & !is.na(covered)) %>% 
              mutate(coverage = .01*as.numeric(str_extract(coverage, "\\d+"))) %>% 
              group_by(model, coverage) %>% 
              summarise(pcover = mean(covered)),
            aes(x = coverage, y = pcover, linetype = model, colour = model), 
            size = 1.2) +
  geom_line(data = ili_coverage %>%
              filter(model == "sir_disc2_hosp_sq_ar1") %>%
              mutate(model = "SIRD_NORM2") %>% 
              dplyr::select(reference_date, horizon, location, model, contains("cover")) %>% 
              pivot_longer(contains("cover"), names_to = "coverage", 
                           values_to = "covered") %>% 
              filter(horizon >= 0 & !is.na(covered)) %>% 
              mutate(coverage = .01*as.numeric(str_extract(coverage, "\\d+"))) %>% 
              group_by(model, coverage) %>% 
              summarise(pcover = mean(covered)),
            aes(x = coverage, y = pcover, linetype = model, colour = model), 
            size = 1.2) +
  geom_line(data = team_coverage %>%
              filter(str_detect(model, "FluSight-baseline")) %>%
              mutate(model = "FluSight-baseline") %>% 
              dplyr::select(reference_date, horizon, location, model, contains("cover")) %>% 
              pivot_longer(contains("cover"), names_to = "coverage", 
                           values_to = "covered") %>% 
              filter(horizon >= 0 & !is.na(covered)) %>% 
              mutate(coverage = .01*as.numeric(str_extract(coverage, "\\d+"))) %>% 
              group_by(model, coverage) %>% 
              summarise(pcover = mean(covered)),
            aes(x = coverage, y = pcover, linetype = model, colour = model), 
            size = 1.2) +
  scale_color_manual(values = c("ASGD_NORM2" = "blue3", 
                                "FluSight-baseline" = "red1",
                                "SIRD_NORM2" = "forestgreen")) +
  scale_linetype_manual(values = c("ASGD_NORM2" = "dashed", 
                                   "FluSight-baseline" = "twodash",
                                   "SIRD_NORM2" = "longdash")) +
  
  xlab("Theoretical Coverage") +
  ylab("Empirical Coverage") +
  theme_bw() +
  theme(legend.position = "none",
        # legend.position = c(.2, .8),
        legend.title = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 12))


plwis <- team_scores %>% 
  filter(str_detect(model, "asg") | str_detect(model, "FluSight-baseline")
         | str_detect(model, "SIRD")) %>% 
  mutate(model = ifelse(str_detect(model, "asg"), "ASGD_NORM2", 
                                   ifelse(str_detect(model, "SIRD"), 
                                   "SIRD_NORM2", "FluSight-baseline"))) %>% 
  mutate(model = factor(model, levels = c("ASGD_NORM2", "SIRD_NORM2",
                        "FluSight-baseline"))) %>% 
  group_by(model, reference_date) %>% 
  summarise(mlwis = mean(lwis, na.rm = TRUE),
            mwis = mean(wis, na.rm = TRUE)) %>% 
  ggplot() +
  geom_line(data = team_scores %>% 
              
              filter(!(model %in% c("UNC_IDD-InfluPaint", "NU_UCSD-GLEAM_AI_FLUH",
                                    "SigSci-CREG", "PSI-PROF_beta", "NIH-Flu_ARIMA",
                                    "GH-model", "UGA_flucast-OKeeffe"))) %>%
              filter(!str_detect(model, "asg") & !str_detect(model, "FluSight-baseline")) %>% 
              group_by(model, reference_date) %>% 
              summarise(mlwis = mean(lwis, na.rm = TRUE),
                        mwis = mean(wis, na.rm = TRUE)),
            aes(x = date(reference_date), y = mlwis, group = model),
            colour = "grey", size = .9) + 
  geom_line(aes(x = date(reference_date), y = mlwis, 
                colour = model, linetype = model), size = 1.2) +
  scale_color_manual(values = c("ASGD_NORM2" = "blue3", 
                                "FluSight-baseline" = "red1",
                                "SIRD_NORM2" = "forestgreen")) +
  scale_linetype_manual(values = c("ASGD_NORM2" = "dashed", 
                                   "FluSight-baseline" = "twodash",
                                   "SIRD_NORM2" = "longdash")) +
  xlab("Date") +
  ylab("LWIS") +
  theme_bw() +
  theme(legend.position = c(.75, .83),
        legend.title = element_blank(),
        axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12))


pwis <- team_scores %>% 
  filter(str_detect(model, "asg") | str_detect(model, "FluSight-baseline")
         | str_detect(model, "SIRD")) %>% 
  mutate(model = ifelse(str_detect(model, "asg"), "ASGD_NORM2", 
                        ifelse(str_detect(model, "SIRD"), 
                               "SIRD_NORM2", "FluSight-baseline"))) %>% 
  mutate(model = factor(model, levels = c("ASGD_NORM2", "SIRD_NORM2",
                        "FluSight-baseline"))) %>%
  group_by(model, reference_date) %>% 
  summarise(mlwis = mean(lwis, na.rm = TRUE),
            mwis = mean(wis, na.rm = TRUE)) %>% 
  ggplot() +
  geom_line(data = team_scores %>% 
              
              filter(!(model %in% c("UNC_IDD-InfluPaint", "NU_UCSD-GLEAM_AI_FLUH",
                                    "SigSci-CREG", "PSI-PROF_beta", "NIH-Flu_ARIMA",
                                    "GH-model", "UGA_flucast-OKeeffe"))) %>% 
              filter(!str_detect(model, "asg") & !str_detect(model, "FluSight-baseline")) %>% 
              group_by(model, reference_date) %>% 
              summarise(mlwis = mean(lwis, na.rm = TRUE),
                        mwis = mean(wis, na.rm = TRUE)),
            aes(x = date(reference_date), y = mwis, group = model),
            colour = "grey", size = .9) + 
  geom_line(aes(x = date(reference_date), y = mwis, 
                colour = model, linetype = model), size = 1.2) +
  scale_color_manual(values = c("ASGD_NORM2" = "blue3", 
                                "FluSight-baseline" = "red1",
                                "SIRD_NORM2" = "forestgreen")) +
  scale_linetype_manual(values = c("ASGD_NORM2" = "dashed", 
                                   "FluSight-baseline" = "twodash",
                                   "SIRD_NORM2" = "longdash")) +
  xlab("Date") +
  ylab("WIS") +
  theme_bw() +
  theme(legend.position = c(.75, .83),
        legend.title = element_blank(),
        axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12))





library(patchwork)

plbox <- rbind(team_scores, ili_scores) %>% 
  filter(model %in% c("asg_disc2_nm_hosp_sq_ar1") |
           str_detect(model, "FluSight") |
           str_detect(model, "SIRD")) %>% 
  mutate(model = ifelse(str_detect(model, "asg"), "ASGD_NORM2",
                        ifelse(str_detect(model, "SIRD"), "SIRD_NORM2",
                        "FluSight-baseline"))) %>% 
  mutate(model = factor(model, levels = c("ASGD_NORM2", "SIRD_NORM2",
                                          "FluSight-baseline"))) %>% 
  select(model, horizon, lwis, wis) %>% 
  pivot_longer(cols = c("wis", "lwis"), 
               names_to = "score", values_to = "wis") %>% 
  mutate(score = ifelse(score == "wis", "WIS", "LWIS")) %>% 
  select(model, score, wis, horizon) %>% unique() %>% 
  filter(score == "LWIS") %>%
  ggplot() +
  geom_boxplot(aes(x = factor(horizon + 1), y = wis, colour = model)) +
  scale_colour_manual(values = c("ASGD_NORM2" = "blue3", "FluSight-baseline" = "red1",
                                 "SIRD_NORM2" = "forestgreen")) +
  # scale_y_log10() +
  xlab("Weeks Ahead") +
  ylab("") +
  # facet_wrap(~score, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        # legend.position = c(.15, .75),
        legend.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))



pbox <- rbind(team_scores, ili_scores) %>% 
  filter(model %in% c("asg_disc2_nm_hosp_sq_ar1") |
           str_detect(model, "FluSight") |
           str_detect(model, "SIRD")) %>% 
  mutate(model = ifelse(str_detect(model, "asg"), "ASGD_NORM2",
                        ifelse(str_detect(model, "SIRD"), "SIRD_NORM2",
                               "FluSight-baseline"))) %>% 
  mutate(model = factor(model, levels = c("ASGD_NORM2", "SIRD_NORM2",
                                          "FluSight-baseline"))) %>%  
  select(model, horizon, lwis, wis) %>% 
  pivot_longer(cols = c("wis", "lwis"), 
               names_to = "score", values_to = "wis") %>% 
  mutate(score = ifelse(score == "wis", "WIS", "LWIS")) %>% 
  select(model, score, wis, horizon) %>% unique() %>% 
  filter(score == "WIS") %>%
  ggplot() +
  geom_boxplot(aes(x = factor(horizon + 1), y = wis, colour = model)) +
  scale_colour_manual(values = c("ASGD_NORM2" = "blue3", 
                                 "FluSight-baseline" = "red1",
                                 "SIRD_NORM2" = "forestgreen")) +
  scale_y_log10() +
  xlab("Weeks Ahead") +
  ylab("") +
  # facet_wrap(~score, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        # legend.position = "none", #c(.15, .75),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))



# (p1 / p2)

grid.arrange(plwis, plbox, pcover, nrow = 2)
plot_grid(plwis, plbox, pcover, nrow = 2, labels = c("a)", "b)", "c)"))
grid.arrange(pwis, pbox, nrow = 1)

library(cowplot)
plot_grid(pwis, pbox, labels = c("a)", "b)"))

library(xtable)
cover_sum <- rbind(team_coverage, ili_coverage %>% unique() %>% 
        filter(model %in% c("asg_disc2_nm_hosp_sq_ar1",
                            "sir_disc2_hosp_sq_ar1"))) %>% 
  mutate(model = ifelse(str_detect(model, "asg"), "ASGD_NORM2", model)) %>% 
  mutate(model = ifelse(str_detect(model, "FluSight"), 
                        "FluSight-baseline", 
                        ifelse(str_detect(model, "sir_d"), "SIRD_NORM2",
                               model))) %>% 
  # filter(!str_detect(model, "FluSight-baseline")) %>% 
  dplyr::select(reference_date, horizon, location, model, contains("cover")) %>% 
  pivot_longer(contains("cover"), names_to = "coverage", 
               values_to = "covered") %>% 
  filter(horizon >= 0 & !is.na(covered)) %>% 
  mutate(coverage = .01*as.numeric(str_extract(coverage, "\\d+"))) %>% 
  group_by(model, coverage) %>% 
  summarise(pcover = mean(covered)) %>%
  group_by(model) %>% 
  summarise(mcov = round(mean((coverage - pcover)^2), 3)) %>% 
  arrange(mcov)

wis_sum <- team_scores %>% unique() %>% 
  filter(horizon != -1) %>% 
  mutate(model = ifelse(str_detect(model, "asg"), "ASGD_NORM2", model)) %>% 
  mutate(model = ifelse(str_detect(model, "FluSight"), 
                        "FluSight-baseline", model)) %>% 
  mutate(model = ifelse(str_detect(model, "SIR"), "SIRD_NORM2", model)) %>% 
  group_by(model) %>% 
  summarise(mdwis = median(wis, na.rm = TRUE),
            mdlwis = median(lwis, na.rm = TRUE),
            mwis = mean(wis, na.rm = TRUE),
            mlwis = mean(lwis, na.rm = TRUE),
            n = n()) %>% 
  mutate(n = n/max(n)*100) %>% 
  arrange(mwis)


rel_scores <- team_scores %>% 
  # filter(model != "asg_disc2_nm_hosp_sq_ar1") %>%
  mutate(model = ifelse(str_detect(model, "asg"), "ASGD_NORM2", model)) %>% 
  mutate(model = ifelse(str_detect(model, "FluSight"), 
                        "FluSight-baseline", model)) %>% 
  mutate(model = ifelse(str_detect(model, "SIRD"), 
                        "SIRD_NORM2", model)) %>% 
  left_join(best_nl, by = c("reference_date", "loc_abb", "horizon")) %>% 
  group_by(model) %>% 
  # filter(model != "GH-model") %>% 
  summarise(mlwis = mean(lwis, na.rm = TRUE),
            nl_mlwis = mean(nl_lwis, na.rm = TRUE),
            mwis = mean(wis, na.rm = TRUE),
            nl_mwis = mean(nl_wis, na.rm = TRUE)) %>% 
  arrange(mwis) %>% 
  mutate(sldiff = round(nl_mlwis/mlwis, 2),
         sdiff = round(nl_mwis/mwis, 2)) %>% 
  select(model, sdiff, sldiff) %>% 
  arrange(sldiff)

wis_sum %>% 
  left_join(cover_sum, by = "model") %>% 
  left_join(rel_scores, by = "model") %>% 
  select(model, sldiff, sdiff, mlwis, mwis, mcov, n) %>% # mdwis, mdlwis, mcov, n) %>% 
  arrange(desc(n), ) %>% 
  as.data.frame() %>% 
  xtable(digits = c(0, 2, 2, 2, 2, 0, 3, 0)) %>% 
  print(include.rownames = FALSE)

#Know about fMRI data analysis. Cross correlation. MIT stuff. 

