source("../flu_forecast_23/get_data_functions_new.R")
library(datasets)
state.name <- c(state.name, "Puerto Rico", "District of Columbia", "US")
state.abb <- c(state.abb, "PR", "DC", "US")
states <- data.frame(loc_abb = state.abb, location_name = state.name)


team_scores <- readRDS("../flu_forecast_23/hub_scores_no_ens.rds") %>% 
  left_join(states, by = "location_name") %>% 
  mutate(model = ifelse(model == "FluSight-baseline", 
                        # "<b><i>FluSight-baseline</i></b>", 
                        "<span style='font-size:12pt'><b><i>FluSight-baseline</i></b></span>",
                        model))
ili_scores <- readRDS("../flu_forecast_23/nl_scores.rds") %>% 
  left_join(states, by = "location_name")
  


rbind(team_scores, ili_scores) %>% 
  # filter(str_detect(model, "asg")) %>%
  # filter(location_name != "Florida") %>% 
  group_by(model) %>% 
  summarise(m = mean(lwis, na.rm = TRUE),
            n = n()) %>% 
  arrange(m)


rbind(team_scores, ili_scores) %>% 
  filter(model %in% c("asg_disc2_nm_hosp_sq_ar1", "FluSight-baseline")) %>% 
  ggplot() +
  geom_boxplot(aes(x = factor(horizon), y = lwis, colour = model)) #+
  scale_y_log10() 

team_scores <- rbind(team_scores, 
                     ili_scores %>% 
                       filter(model == "asg_disc2_nm_hosp_sq_ar1")) 

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
  stat_function(fun = function(x) x + 0, color = "blue") +
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
    low = "royalblue", high = "red", mid = "white",
    midpoint = 1) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  xlab("") +
  ylab("LWIS") +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey50", color = NA),
        # legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        # legend.position = "none",
        # axis.text.y = element_text(size = 10, angle = 20),
        axis.title = element_text(size = 16),
        axis.text.y = element_markdown())



state_diff <- team_diffs %>% 
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
    low = "royalblue", high = "red", mid = "white",
    midpoint = 1) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  xlab("Location") +
  ylab("WIS") +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey50", color = NA),
        # legend.position = "none",
        axis.text.x = element_text(size = 7, angle = 90),
        legend.title = element_blank(),
        # legend.position = "none",
        # axis.text.y = element_text(size = 10, angle = 20),
        axis.title = element_text(size = 16),
        axis.text.y = element_markdown())














team_ddiffs <- team_scores %>% 
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
  mutate(model = factor(model, levels = scores_ord$model)) %>%  
  ggplot() +
  geom_tile(aes(y = model,
                x = date(reference_date),
                fill = sldiff)) +
  scale_fill_gradient2(
    low = "royalblue", high = "red", mid = "white",
    midpoint = 1) +
  scale_x_date(expand = c(0,0),
               breaks = pretty_breaks(n = 6)) +
  scale_y_discrete(expand = c(0,0)) +
  ylab("LWIS") +
  xlab("") +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey50", color = NA),
        legend.title = element_blank(),
        # legend.position = "none",
        # axis.text.x = element_text(size = 11),
        # legend.position = "none",
        # axis.text.y = element_text(size = 10, angle = 20),
        axis.title = element_text(size = 16),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_markdown())
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank())



date_diff <- team_ddiffs %>%
  mutate(model = factor(model, levels = scores_ord$model)) %>%  
  ggplot() +
  geom_tile(aes(y = model,
                x = date(reference_date),
                fill = sdiff)) +
  scale_fill_gradient2(
    low = "royalblue", high = "red", mid = "white",
    midpoint = 1) +
  scale_x_date(expand = c(0,0),
               breaks = pretty_breaks(n = 6)) +
  scale_y_discrete(expand = c(0,0)) +
  ylab("WIS") +
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

grid.arrange(state_ldiff, date_ldiff, state_diff, date_diff, 
             ncol = 2, 
             # left = textGrob("Model", 
             #                 rot = 90, vjust = 1, 
             #                 gp = gpar(fontsize = 16)),
             widths = c(1.5, 1))

state_ldiff / state_diff
date_ldiff / date_diff

##################################
#############Coverage#############
##################################
# saveRDS(ili_coverage, "./nl_cover.rds")
team_coverage <- readRDS("./hub_cover_no_ens.rds")
# ili_coverage <- readRDS("./nl_cover.rds")
p3 <- team_coverage %>% 
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
              mutate(model = "ASG") %>% 
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
  scale_color_manual(values = c("ASG" = "blue", "FluSight-baseline" = "red")) +
  scale_linetype_manual(values = c("ASG" = "dashed", 
                                   "FluSight-baseline" = "twodash")) +
  
  xlab("Theoretical Coverage") +
  ylab("Empirical Coverage") +
  theme_bw() +
  theme(legend.position = c(.2, .8),
        legend.title = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 12))


team_scores %>% 
  filter(str_detect(model, "asg") | str_detect(model, "FluSight-baseline")) %>% 
  mutate(model = ifelse(str_detect(model, "asg"), "ASG", 
                                   "FluSight-baseline")) %>% 
  group_by(model, reference_date) %>% 
  summarise(mlwis = mean(lwis, na.rm = TRUE),
            mwis = mean(wis, na.rm = TRUE)) %>% 
  ggplot() +
  geom_line(aes(x = date(reference_date), y = mlwis, 
                colour = model, linetype = model), size = .9) +
  scale_colour_manual(values = c("ASG" = "blue", 
                                 "FluSight-baseline" = "red")) +
  xlab("") +
  ylab("LWIS") +
  theme_bw() +
  theme(legend.position = c(.2, .89),
        legend.title = element_blank(),
        axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12))

team_scores %>% 
  filter(str_detect(model, "asg") | str_detect(model, "FluSight-baseline")) %>% 
  mutate(model = ifelse(str_detect(model, "asg"), "ASG", 
                        "FluSight-baseline")) %>% 
  group_by(model, reference_date) %>% 
  summarise(mlwis = mean(lwis, na.rm = TRUE),
            mwis = mean(wis, na.rm = TRUE)) %>% 
  ggplot() +
  geom_line(aes(x = date(reference_date), y = mlwis, 
                colour = model, linetype = model), size = .9) +
  scale_colour_manual(values = c("ASG" = "blue", 
                                 "FluSight-baseline" = "red")) +
  xlab("Date") +
  ylab("WIS") +
  theme_bw() +
  theme(legend.position = c(.2, .89),
        legend.title = element_blank(),
        axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 12),
        legend.text = element_text(size = 12))




library(patchwork)

p1 <- rbind(team_scores, ili_scores) %>% 
  filter(model %in% c("asg_disc2_nm_hosp_sq_ar1") |
           str_detect(model, "FluSight")) %>% 
  mutate(model = ifelse(str_detect(model, "asg"), "ASG",
                                   "FluSight-baseline")) %>% 
  select(model, horizon, lwis, wis) %>% 
  pivot_longer(cols = c("wis", "lwis"), 
               names_to = "score", values_to = "wis") %>% 
  mutate(score = ifelse(score == "wis", "WIS", "LWIS")) %>% 
  select(model, score, wis, horizon) %>% unique() %>% 
  filter(score == "LWIS") %>%
  ggplot() +
  geom_boxplot(aes(x = factor(horizon + 1), y = wis, colour = model)) +
  scale_colour_manual(values = c("ASG" = "blue", "FluSight-baseline" = "red")) +
  # scale_y_log10() +
  xlab("") +
  ylab("LWIS") +
  # facet_wrap(~score, scales = "free_y") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(.15, .75),
        legend.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))



p2 <- rbind(team_scores, ili_scores) %>% 
  filter(model %in% c("asg_disc2_nm_hosp_sq_ar1") |
           str_detect(model, "FluSight")) %>% 
  mutate(model = ifelse(str_detect(model, "asg"), "ASG",
                        "FluSight-baseline")) %>% 
  select(model, horizon, lwis, wis) %>% 
  pivot_longer(cols = c("wis", "lwis"), 
               names_to = "score", values_to = "wis") %>% 
  mutate(score = ifelse(score == "wis", "WIS", "LWIS")) %>% 
  select(model, score, wis, horizon) %>% unique() %>% 
  filter(score == "WIS") %>%
  ggplot() +
  geom_boxplot(aes(x = factor(horizon + 1), y = wis, colour = model)) +
  scale_colour_manual(values = c("ASG" = "blue", "FluSight-baseline" = "red")) +
  scale_y_log10() +
  xlab("Horizon") +
  ylab("WIS") +
  # facet_wrap(~score, scales = "free_y") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "none", #c(.15, .75),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))



(p1 / p2)



rbind(team_coverage, ili_coverage) %>% 
  # filter(!str_detect(model, "FluSight-baseline")) %>% 
  dplyr::select(reference_date, horizon, location, model, contains("cover")) %>% 
  pivot_longer(contains("cover"), names_to = "coverage", 
               values_to = "covered") %>% 
  filter(horizon >= 0 & !is.na(covered)) %>% 
  mutate(coverage = .01*as.numeric(str_extract(coverage, "\\d+"))) %>% 
  group_by(model, coverage) %>% 
  summarise(pcover = mean(covered)) %>%
  group_by(model) %>% 
  summarise(m = mean((coverage - pcover)^2)) %>% 
  arrange(m)



