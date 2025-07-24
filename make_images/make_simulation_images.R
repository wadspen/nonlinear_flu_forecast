library(ggplot2)
library(dplyr)
library(stringr)
library(ggforce)



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


scores %>% 
  mutate(model = factor(model, levels = c("ASGD", "ASG", "ARIMA", 
                                   "BASE", "SIRD", "SIR"))) %>% 
  mutate(week_ahead = week_ahead + 1) %>% 
  mutate(week_ahead = ifelse(week_ahead == 1, paste(week_ahead, 
                                                     "week", sep = " "),
                              paste(week_ahead, "weeks", sep = " "))) %>% 
  ggplot() +
  geom_boxplot(aes(x = model, y = wis)) +
  ggh4x::facet_grid2(last_week~week_ahead, 
                     scales = "free_y", independent = "y") +
  ylab("WIS") +
  xlab("Model") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, angle = 35, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 17),
        strip.text = element_text(size = 11))

library(xtable)
scores %>% 
  group_by(model) %>% 
  summarise(WIS = mean(wis),
            MAE = mean(lmae),
            c50 = mean(cover50),
            c95 = mean(cover95)) %>% 
  arrange(WIS) %>% 
  as.data.frame() %>% 
  xtable(digits = c()) %>% 
  print(include.rownames = FALSE)


install.packages("devtools")  # if needed
remove.packages("ggforce") 
devtools::install_github("thomasp85/ggforce")




























