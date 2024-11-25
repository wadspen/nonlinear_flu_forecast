setwd("~/flu_research/Prelim Content/")
source("~/flu_research/nonlinear_flu_forecast/forecasts_23/get_data_functions.R")
source("~/flu_research/nonlinear_flu_forecast/forecasts_23/mle_functions.R")
library(ggplot2)
library(dplyr)
library(extraDistr)
library(tidyr)
library(ggpubr)




################################################
###############alphas, sigma####################
################################################

all_posts <- data.frame()
for (week in c(14, 20, 26)) {
  for (mod in c("", "_log", "_lst")) {
    for (sq in c("", "_sq")) {
      file_name <- paste0("hosp_fits/asg_disc_hosp", mod, sq, "_ar1_US_",
                          week, "_posterior.csv")
      
      post <- read.csv(file_name) 
      if (mod %in% c("", "_lst")) {
        post <- post %>% 
          select(contains("alpha") & contains(".2."), 
                 contains("sigma_epsilon") & contains(".2.")) %>% 
          select(-contains("pha2"))
        colnames(post) <- c("alpha0", "alpha1", "alpha3", "sigma_epsilon")
      } else if (mod == "_log") {
        post <- post %>% 
          select(contains("alpha"), contains("sigma")) %>% 
          select(-contains("pha2"))
      }
      
      post$model <- mod
      post$sq <- sq
      post$season_week <- week
      
      all_posts <- rbind(all_posts, post)
      
    }
  }
  
}


all_posts <- all_posts %>% 
  mutate(modwk = paste0(model, season_week))


posts_95 <- all_posts %>% 
  group_by(modwk, model, sq, season_week) %>% 
  reframe(alpha0_low95 = quantile(alpha0, probs = .025),
          alpha0_upp95 = quantile(alpha0, probs = .975),
          alpha1_low95 = quantile(alpha1, probs = .025),
          alpha1_upp95 = quantile(alpha1, probs = .975),
          phi_low95 = quantile(alpha3, probs = .025),
          phi_upp95 = quantile(alpha3, probs = .975),
          sigma_low95 = quantile(sigma_epsilon, probs = 0),
          sigma_upp95 = quantile(sigma_epsilon, probs = .95))

posts_95 <- posts_95 %>% 
  mutate(model = ifelse(model == "", "NORM", 
                        ifelse(model == "_log", "LNORM", "LST")),
         sq = ifelse(sq == "", "w/out square ILI", "w/ square ILI"))

library(truncnorm)
posts_95 <- posts_95 %>%
  bind_rows(data.frame(modwk = "prior",
                       model = "Prior",
                       expand.grid(season_week = c(14, 20, 26),
                            sq = c("w/out square ILI", "w/ square ILI")),
                       alpha0_low95 = qnorm(.025, 0, 5),
                       alpha0_upp95 = qnorm(.975, 0, 5),
                       alpha1_low95 = qnorm(.025, 0, 5),
                       alpha1_upp95 = qnorm(.975, 0, 5),
                       phi_low95 = qtruncnorm(.025, a = -1, b = 1, 0, .4),
                       phi_upp95 = qtruncnorm(.975, a = -1, b = 1, 0, .4),
                       sigma_low95 = qhnorm(0, 4),
                       sigma_upp95 = qhnorm(.95, 4)))

posts_95 <- posts_95 %>% 
  mutate(model = factor(model, levels = c("Prior", "LNORM", "LST", "NORM")))


posts_95 %>% 
  ggplot() +
  geom_segment(aes(y = model, x = alpha0_low95, 
                   xend = alpha0_upp95), 
               size = 2) +
  facet_grid(sq ~ season_week) +
  xlab(expression(alpha[0])) +
  ylab("") +
  theme_bw() + 
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=23),
        strip.text = element_text(
          size = 12),
        legend.position = "none")


posts_95 %>% 
  ggplot() +
  geom_segment(aes(y = model, x = alpha1_low95, 
                   xend = alpha1_upp95), 
               size = 2) +
  facet_grid(sq ~ season_week) +
  xlab(expression(alpha[1])) +
  ylab("") +
  theme_bw() + 
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=23),
        strip.text = element_text(
          size = 12),
        legend.position = "none")


posts_95 %>% 
  ggplot() +
  geom_segment(aes(y = model, x = phi_low95, 
                   xend = phi_upp95), 
               size = 2) +
  facet_grid(sq ~ season_week) +
  xlab(expression(phi)) +
  ylab("") +
  theme_bw() + 
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=23),
        strip.text = element_text(
          size = 12),
        legend.position = "none")


posts_95 %>% 
  # filter(model != "Prior") %>% 
  ggplot() +
  geom_segment(aes(y = model, x = sigma_low95, 
                   xend = sigma_upp95), 
               size = 2) +
  facet_grid(sq ~ season_week) +
  xlab(expression(sigma[epsilon])) +
  ylab("") +
  coord_cartesian(xlim=c(0, .5)) +
  theme_bw() + 
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=23),
        strip.text = element_text(
          size = 12),
        legend.position = "none")




















################################################
###############alphas, sigma####################
################################################

all_posts2 <- data.frame()
for (week in c(14, 20, 26)) {
  for (mod in c("", "_log", "_lst")) {
    for (sq in c("_sq")) {
      file_name <- paste0("hosp_fits/asg_disc_hosp", mod, sq, "_ar1_US_",
                          week, "_posterior.csv")
      
      post <- read.csv(file_name) 
      if (mod %in% c("", "_lst")) {
        post <- post %>% 
          select(contains("alpha") & contains(".2."), 
                 contains("sigma_epsilon") & contains(".2.")) %>% 
          select(contains("pha2"))
        colnames(post) <- c("alpha2")
      } else if (mod == "_log") {
        post <- post %>%  
          select(contains("pha2"))
        colnames(post) <- c("alpha2")
      }
      
      post$model <- mod
      post$sq <- sq
      post$season_week <- week
      
      all_posts2 <- rbind(all_posts2, post)
      
    }
  }
  
}


posts_95 <- all_posts2 %>% 
  group_by(model, sq, season_week) %>% 
  reframe(alpha2_low95 = quantile(alpha2, probs = .025),
          alpha2_upp95 = quantile(alpha2, probs = .975))

posts_95 <- posts_95 %>% 
  mutate(model = ifelse(model == "", "NORM", 
                        ifelse(model == "_log", "LNORM", "LST")),
         sq = ifelse(sq == "", "w/out square ILI", "w/ square ILI"))


posts_95 <- posts_95 %>%
  bind_rows(data.frame(modwk = "prior",
                       model = "Prior",
                       expand.grid(season_week = c(14, 20, 26),
                                   sq = c("w/out square ILI", "w/ square ILI")),
                       alpha2_low95 = qnorm(.025, 0, 5),
                       alpha2_upp95 = qnorm(.975, 0, 5)))

posts_95 <- posts_95 %>% 
  mutate(model = factor(model, levels = c("Prior", "LNORM", "LST", "NORM")))


# library(truncnorm)
# posts_95 <- posts_95 %>% 
#   bind_rows(data.frame(modwk = "prior",
#                        model = "prior",
#                        season_week = 8,
#                        alpha0_low95 = qnorm(.025, 0, 1000),
#                        alpha0_upp95 = qnorm(.975, 0, 1000)))


posts_95 %>% 
  # filter(model != "LNORM") %>% 
  ggplot() +
  geom_segment(aes(y = model, x = alpha2_low95, 
                   xend = alpha2_upp95), 
               size = 2) +
  facet_grid(season_week~.) +
  xlab(expression(alpha[2])) +
  coord_cartesian(xlim=c(0, 1)) +
  ylab("") +
  theme_bw() + 
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=23),
        strip.text = element_text(
          size = 12),
        legend.position = "none")


























all_posts3 <- data.frame()
for (week in c(14, 20, 26)) {
  for (mod in c("_lst")) {
    for (sq in c("", "_sq")) {
      file_name <- paste0("hosp_fits/asg_disc_hosp", mod, sq, "_ar1_US_",
                          week, "_posterior.csv")
      
      post <- read.csv(file_name) 
      if (mod %in% c("", "_lst")) {
        post <- post %>% 
          select(contains("nu") & contains(".2."))
        colnames(post) <- c("nu")
      } else if (mod == "_log") {
        post <- post %>%  
          select(contains("pha2"))
        colnames(post) <- c("alpha2")
      }
      
      post$model <- mod
      post$sq <- sq
      post$season_week <- week
      
      all_posts3 <- rbind(all_posts3, post)
      
    }
  }
  
}


posts_95 <- all_posts3 %>% 
  group_by(model, sq, season_week) %>% 
  reframe(nu_low95 = quantile(nu, probs = .025),
          nu_upp95 = quantile(nu, probs = .975))

posts_95 <- posts_95 %>% 
  mutate(model = ifelse(model == "", "NORM", 
                        ifelse(model == "_log", "LNORM", "LST")),
         sq = ifelse(sq == "", "w/out square ILI", "w/ square ILI"))

posts_95 <- posts_95 %>%
  bind_rows(data.frame(modwk = "prior",
                       model = "Prior",
                       expand.grid(season_week = c(14, 20, 26),
                                   sq = c("w/out square ILI", "w/ square ILI")),
                       nu_low95 = qhnorm(0, 15),
                       nu_upp95 = qhnorm(.95, 15)))

posts_95 <- posts_95 %>% 
  mutate(model = factor(model, levels = c("Prior", "LNORM", "LST", "NORM")))


posts_95 %>% 
  # filter(model != "LNORM") %>% 
  ggplot() +
  geom_segment(aes(y = model, x = nu_low95, 
                   xend = nu_upp95), 
               size = 2) +
  facet_grid(sq~season_week) +
  xlab(expression(omega)) +
  # coord_cartesian(xlim=c(0, 1)) +
  ylab("") +
  theme_bw() + 
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=23),
        strip.text = element_text(
          size = 12),
        legend.position = "none")

