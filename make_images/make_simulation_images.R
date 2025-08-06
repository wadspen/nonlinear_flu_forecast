library(ggplot2)
library(dplyr)
library(stringr)
library(ggforce)
library(tidyr)



score_loc <- "../simulation/hosp_sim/scores_"
models <- c("arima_full", "base_full", "asg_nm", 
            "asg_disc2_nm", "sir", "sir_disc2")


paths <- paste0(score_loc, models, ".csv")
names <- c("ARIMA", "BASE", "ASG", "ASGD", "SIR", "SIRD")

scores <- data.frame()
for (p in 1:length(paths)) {
  score <- read.csv(paths[p])
  score$model <- names[p]
  
  scores <- rbind(scores, score)
}


scores <- scores %>% 
  mutate(model = factor(model, levels = c("ASGD", "SIRD", "ASG", "SIR",
                                          "ARIMA", "BASE"))) %>% 
  mutate(week_ahead = factor(week_ahead + 1))# %>% 
  # mutate(week_ahead = ifelse(week_ahead == 1, paste(week_ahead, 
  #                                                    "week", sep = " "),
  #                             paste(week_ahead, "weeks", sep = " "))) 


sim_cov <- scores %>% 
  select(model, contains("cover")) %>% 
  pivot_longer(contains("cover"), , names_to = "coverage", 
               values_to = "covered") %>% 
  mutate(coverage = .01*as.numeric(str_extract(coverage, "\\d+"))) %>% 
  group_by(model, coverage) %>% 
  summarise(pcover = mean(covered)) %>% 
  # filter(model == "FluSight-baseline") %>% 
  ggplot(aes(x = coverage, y = pcover, linetype = model, colour = model)) +
  geom_line(size = 1.4) +
  geom_abline(intercept = 0, slope = 1, size = 1.2) +
  scale_colour_manual(values = c("ASGD" = "#E69F00",  # orange
                                 "BASE" = "#E66100",
                                 "ASG" = "#56B4E9",  # sky blue
                                 "ARIMA" = "#009E73",  # bluish green
                                 "SIR" = "#D6C840",  # yellow
                                 "SIRD" = "#0072B2"   # blue
  )) +
  scale_linetype_manual(values = c("ASGD" = "solid",  # orange
                                   "BASE" = "12345678",
                                   "ASG" = "longdash",  # sky blue
                                   "ARIMA" = "dotdash",  # bluish green
                                   "SIR" = "twodash",  # yellow
                                   "SIRD" = "F1"   # blue
  )) +
  labs(colour = "Model", linetype = "Model") +
  xlab("Theoretical Coverage") +
  ylab("Empirical Coverage") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, hjust = 1),
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank(),
        # legend.position = c(.1, .75),
        # legend.title = "none",
        legend.position = "none",
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 17),
        legend.title = element_text(size = 13),
        strip.text = element_text(size = 11))


# scores %>% 
#   ggplot() +
#   geom_boxplot(aes(x = model, y = wis)) +
#   ggh4x::facet_grid2(last_week~week_ahead, 
#                      scales = "free_y", independent = "y") +
#   ylab("WIS") +
#   xlab("Model") +
#   theme_bw() +
#   theme(axis.text.x = element_text(size = 12, angle = 35, hjust = 1),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.title = element_text(size = 17),
#         strip.text = element_text(size = 11))


sim_wis <- scores %>% 
  group_by(model, last_week) %>% 
  summarise(mwis = mean(wis)) %>% 
  ggplot() +
  geom_line(aes(x = last_week, y = mwis, colour = model, linetype = model),
            size = 1.4) +
  scale_colour_manual(values = c("ASGD" = "#E69F00",  # orange
                                 "BASE" = "#E66100",
                                 "ASG" = "#56B4E9",  # sky blue
                                 "ARIMA" = "#009E73",  # bluish green
                                 "SIR" = "#D6C840",  # yellow
                                 "SIRD" = "#0072B2"   # blue
  )) +
  scale_linetype_manual(values = c("ASGD" = "solid",  # orange
                                   "BASE" = "12345678",
                                   "ASG" = "longdash",  # sky blue
                                   "ARIMA" = "dotdash",  # bluish green
                                   "SIR" = "twodash",  # yellow
                                   "SIRD" = "F1"   # blue
  )) +
  labs(colour = "Model", linetype = "Model") +
  ylab("WIS") +
  xlab("Week") +
  guides(color = guide_legend(direction = "horizontal", nrow = 2)) +
  theme_bw() +
  theme(legend.position = c(.69, .85),
        # legend.title = element_blank(),
        axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12))
  


sim_box <- scores %>% 
  ggplot() +
  geom_boxplot(aes(x = week_ahead, y = wis, colour = model)) +
  # facet_wrap(~week_ahead, scales = "free_y") +
  # coord_cartesian(ylim = c(0, .4)) +
  scale_colour_manual(values = c("ASGD" = "#E69F00",  # orange
                                 "BASE" = "#E66100",
                                 "ASG" = "#56B4E9",  # sky blue
                                 "ARIMA" = "#009E73",  # bluish green
                                 "SIR" = "#D6C840",  # yellow
                                 "SIRD" = "#0072B2"   # blue
  )) +
  labs(colour = "Model") +
  ylab("") +
  xlab("Weeks Ahead") +
  theme_bw() +
  theme(legend.position = c(.14, .72),
        # legend.title = element_blank(),
        axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12))

library(cowplot)
# grid.arrange(plwis, plbox, pcover, nrow = 2)
plot_grid(sim_wis, sim_box, sim_cov, nrow = 2, labels = c("a)", "b)", "c)"))
  

library(xtable)
scores %>% 
  group_by(model) %>% 
  summarise(WIS = mean(wis),
            # MAE = mean(lmae),
            c50 = mean(cover50),
            c95 = mean(cover95)) %>% 
  arrange(WIS) %>% 
  as.data.frame() %>% 
  xtable(digits = c()) %>% 
  print(include.rownames = FALSE)


install.packages("devtools")  # if needed
remove.packages("ggforce") 
devtools::install_github("thomasp85/ggforce")




























