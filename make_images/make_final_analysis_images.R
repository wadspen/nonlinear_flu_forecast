source("./flu_forecast_23/get_data_functions_new.R")
team_scores <- readRDS("./flu_forecast_23/hub_scores_no_ens.rds")
ili_scores <- readRDS("./flu_forecast_23/nl_scores.rds")

rbind(team_scores, ili_scores) %>% 
  filter(str_detect(model, "asg")) %>%
  filter(location_name != "Florida") %>% 
  group_by(model) %>% 
  summarise(m = mean(lwis, na.rm = TRUE),
            n = n()) %>% 
  arrange(m)


rbind(team_scores, ili_scores) %>% 
  filter(model %in% c("asg_disc2_nm_hosp_sq_ar1", "FluSight-baseline")) %>% 
  ggplot() +
  geom_boxplot(aes(x = factor(horizon), y = wis, colour = model)) +
  scale_y_log10() 

team_scores <- rbind(team_scores, ili_scores %>% filter(model == "asg_disc2_nm_hosp_sq_ar1"))

best_nl <- ili_scores %>% 
  filter(model == "asg_disc2_nm_hosp_sq_ar1") %>% 
  mutate(nl_wis = wis, nl_lwis = lwis, nl_mae = mae) %>% 
  dplyr::select(reference_date, location_name, horizon, nl_wis, nl_lwis, nl_mae)

team_scores %>% 
  left_join(best_nl, by = c("reference_date", "location_name", "horizon")) %>% 
  group_by(model, reference_date) %>% 
  summarise(mwis = mean(wis, na.rm = TRUE),
            nl_mwis = mean(nl_wis, na.rm = TRUE)) %>% 
  # filter(model == "FluSight-baseline") %>%
  ggplot() +
  geom_line(aes(x = nl_mwis, y = mwis, colour = model)) +
  stat_function(fun = function(x) x + 0, color = "blue") +
  scale_y_log10() +
  scale_x_log10()


team_scores %>% 
  left_join(best_nl, by = c("reference_date", "location_name", "horizon")) %>% 
  group_by(model, reference_date) %>% 
  summarise(mwis = mean(wis, na.rm = TRUE),
            nl_mwis = mean(nl_wis, na.rm = TRUE)) %>% 
  filter(model %in% c("asg_disc2_nm_hosp_sq_ar1", "FluSight-baseline")) %>%
  ggplot() +
  geom_line(aes(x = date(reference_date), y = mwis, colour = model))

m_scores <- team_scores %>% 
  # filter(location_name == "US") %>%
  group_by(model, location_name) %>% 
  summarise(mrwis = median(wis, na.rm = TRUE)) 


m_scores %>% 
  filter(model != "GH-model") %>% 
  filter(mrwis < 200) %>% 
  ggplot() +
  geom_tile(aes(y = model,
                x = location_name,
                fill = mrwis)) +
  scale_fill_gradient2(
    #low = "#5ab4ac", high = "#d8b365", mid = "white",
    # low = "#5ab4ac", high = "black", mid = "white",
    low = "navyblue", high = "red", mid = "white",
    midpoint = 1, #mean(m_scores$mrwis, na.rm = TRUE),
    # limits = c(min(m_scores$mrwis), max(m_scores$mrwis)),
    na.value = NA) #+
  # coord_flip() +
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
