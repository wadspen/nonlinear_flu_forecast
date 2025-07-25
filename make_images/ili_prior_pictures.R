library(dplyr)
library(ggplot2)
library(stringr)
library(truncnorm)

weeks <- c(14, 20, 26)
mods <- c("asg_nm", "asg_disc2_nm", "sir", "sir_disc2")


asg14 <- readRDS("../flu_forecast_23/asg_nm_summary.rds") %>% 
  mutate(model = "ASG")
asgd14 <- readRDS("../flu_forecast_23/asg_disc2_nm_summary.rds") %>% 
  mutate(model = "ASGD")

sir14 <- readRDS("../flu_forecast_23/sir_summary.rds") %>% 
  mutate(model = "SIR")

sird14 <- readRDS("../flu_forecast_23/sir_summary.rds") %>% 
  mutate(model = "SIRD")

asgd14 %>% filter(rhat > 1.1)

asg_model_post <- rbind(asg14, asgd14) %>% 
  select(model, variable, q2.5, q97.5)

sir_model_post <- rbind(sir14, sird14) %>% 
  select(model, variable, q2.5, q97.5) %>% 
  mutate(variable = as.character(variable)) 
  

variables <- c("I0", "beta", "rho")
sm <- c(.005, .8, .68)
sc <- c(.03, .3, .08)
st <- c(.1, Inf, Inf)

sir_samps <- apply(cbind(st, sm, sc), MARGIN = 1, FUN = function(x) {
                                                          rtruncnorm(2000, 0, 
                                                                     x[1], 
                                                                     x[2], 
                                                                     sqrt(x[3]))
                                                                    })

sir_bounds <- apply(sir_samps, MARGIN = 2, 
                    FUN = quantile, probs = c(.025, .975)) %>% 
  t

colnames(sir_bounds) <- c("q2.5", "q97.5")
sir_bounds <- as.data.frame(sir_bounds)
sir_priors <- cbind(model = "prior", variable = variables, sir_bounds)


one <- rbind(sir_model_post); one$week <- "14"
two <- rbind(sir_model_post); two$week <- "20"
three <- rbind(sir_model_post, sir_priors); three$week <- "26"

sir_plot_data <- rbind(one, two, three) %>%
  filter(variable %in% c("I0[26]", "beta[26]", "rho[26]")) %>% 
  mutate(variable = str_replace(variable, "\\[26\\]", "")) %>% 
  mutate(
    week = ifelse(model == "prior", "prior", week),
    week = factor(week, levels = c("prior", "14", "20", "26")),
    model = factor(model, levels = c("prior", "SIR", "SIRD")),
    
    # Convert factor to numeric for offsetting
    week_num = as.numeric(week),
    
    # Offset by model
    y_pos = case_when(
      model == "SIR" ~ week_num + 0.04,
      model == "SIRD" ~ week_num - 0.04,
      TRUE ~ week_num
    )
  ) %>%
  filter(str_detect(variable, "I0") | str_detect(variable, "beta") |
           str_detect(variable, "rho")) %>% 
  mutate(variable = factor(variable, levels = c("I0", "beta", "rho")))





sir_plot_data %>% 
  ggplot() +
  geom_segment(aes(x = q2.5, xend = q97.5, y = y_pos, yend = y_pos, colour = model),
               size = 1.2) +
  scale_y_continuous(
    breaks = sort(unique(plot_data$week_num)),
    labels = levels(plot_data$week)
  ) +
  facet_wrap(~variable, scales = "free_x", 
             labeller = labeller(variable = label_parsed),
             nrow = 1) +
  labs(y = "Week", x = "", colour = "Model") +
  scale_colour_manual(values = c("prior" = "#E69F00",  # orange
                                 "ASG" = "#56B4E9",  # sky blue
                                 "ASGD" = "#009E73",  # bluish green
                                 "SIR" = "#F0E442",  # yellow
                                 "SIRD" = "#0072B2"   # blue
  )) + 
  theme_bw() +
  theme(legend.position = c(.95,.65),
        strip.text = element_text(size = 12))

























#######################################
#############ASG#######################
#######################################






variables <- paste0("theta[", 1:5, "]")
m0 <- c(-4.5, .3, 23, 3.69, 4.7)
c0 <- c(.3, .2, 5, 2, 2)

theta_samps <- apply(cbind(m0, c0), MARGIN = 1, FUN = function(x) {
                                                    rnorm(2000, x[1], sqrt(x[2]))
                                                        })


theta_bounds <- apply(theta_samps, MARGIN = 2, FUN = quantile, probs = c(.025, .975)) %>% 
  t

colnames(theta_bounds) <- c("q2.5", "q97.5")
theta_bounds <- as.data.frame(theta_bounds)
theta_priors <- cbind(model = "prior", variable = variables, theta_bounds)

zeta_bounds <- rtruncnorm(2000, 0, Inf, 0, 4) %>% 
  quantile(., probs = c(.025, .975)) %>% 
  t %>% 
  data.frame()

colnames(zeta_bounds) <- c("q2.5", "q97.5")
variables <- paste0("zeta[", 1:5, "]")
zeta_priors <- cbind(model = "prior", variable = variables, zeta_bounds)

priors <- rbind(theta_priors, zeta_priors)


one <- rbind(asg_model_post); one$week <- "14"
two <- rbind(asg_model_post); two$week <- "20"
three <- rbind(asg_model_post, priors); three$week <- "26"

asg_plot_data <- rbind(one, two, three) %>%
  mutate(
    week = ifelse(model == "prior", "prior", week),
    week = factor(week, levels = c("prior", "14", "20", "26")),
    model = factor(model, levels = c("prior", "ASG", "ASGD")),
    
    # Convert factor to numeric for offsetting
    week_num = as.numeric(week),
    
    # Offset by model
    y_pos = case_when(
      model == "ASG" ~ week_num + 0.04,
      model == "ASGD" ~ week_num - 0.04,
      TRUE ~ week_num
    )
  ) %>%
  filter((str_detect(variable, "theta") & !str_detect(variable, "theta_s")) |
         str_detect(variable, "zeta")) %>% 
  mutate(variable_label = recode(variable,
                                 "theta[1]" = "lambda",
                                 "theta[2]" = "log(eta)",
                                 "theta[3]" = "mu",
                                 "theta[4]" = "log(sigma[1]^2)",
                                 "theta[5]" = "log(sigma[2]^2)",
                                 "zeta[1]" = "zeta[1]",
                                 "zeta[2]" = "zeta[2]",
                                 "zeta[3]" = "zeta[3]",
                                 "zeta[4]" = "zeta[4]",
                                 "zeta[5]" = "zeta[5]"
  )) %>% 
  mutate(variable_label = factor(variable_label, 
                                 levels = c("lambda",
                                 "log(eta)",
                                 "mu",
                                 "log(sigma[1]^2)",
                                 "log(sigma[2]^2)",
                                 "zeta[1]",
                                 "zeta[2]",
                                 "zeta[3]",
                                 "zeta[4]",
                                 "zeta[5]")))




asg_plot_data %>% 
  ggplot() +
  geom_segment(aes(x = q2.5, xend = q97.5, y = y_pos, yend = y_pos, colour = model),
               size = 1.2) +
  scale_y_continuous(
    breaks = sort(unique(plot_data$week_num)),
    labels = levels(plot_data$week)
  ) +
  facet_wrap(~variable_label, scales = "free_x", 
             labeller = labeller(variable_label = label_parsed),
             nrow = 2) +
  labs(y = "Week", x = "", colour = "Model") +
  scale_colour_manual(values = c("prior" = "#E69F00",  # orange
                                 "ASG" = "#56B4E9",  # sky blue
                                 "ASGD" = "#009E73",  # bluish green
                                 "SIR" = "#F0E442",  # yellow
                                 "SIRD" = "#0072B2"   # blue
                                 )) + 
  theme_bw() +
  theme(legend.position = c(.95,.85),
        strip.text = element_text(size = 12))









#####################################
#############cOMB####################
#####################################



comb_model_post <- rbind(asg_model_post, sir_model_post) %>% 
  filter(variable %in% c("sigma_gamma", "sigma_upsilon",
                         "kappas[26]")) %>% 
  mutate(variable = str_replace(variable, "s\\[26\\]", ""))

variables <- c("kappa", "sigma_gamma", "sigma_upsilon")
sds <- c(10000, .02, .1) %>% as.matrix()
samps <- apply(sds, MARGIN = 1, FUN = function(x) {
  rtruncnorm(2000, 0, 
             Inf, 
             0, 
             x)})


bounds <- apply(samps, MARGIN = 2, 
                    FUN = quantile, probs = c(.025, .975)) %>% 
  t

colnames(bounds) <- c("q2.5", "q97.5")
bounds <- as.data.frame(bounds)
priors <- cbind(model = "prior", variable = variables, bounds, week = "prior")

one <- comb_model_post %>%
  mutate(week = "14")
two <- comb_model_post %>%
  mutate(week = "20")
three <- comb_model_post %>%
  mutate(week = "26")





comb_plot_data <- rbind(one, two, three, priors) %>%
  mutate(
    week = ifelse(model == "prior", "prior", week),
    week = factor(week, levels = c("prior", "14", "20", "26")),
    model = factor(model, levels = c("prior", "ASG", "ASGD", "SIR", "SIRD")),
    
    # Convert factor to numeric for offsetting
    week_num = as.numeric(week),
    
    # Offset by model
    y_pos = case_when(
      model == "ASG" ~ week_num + 0.12,
      model == "ASGD" ~ week_num + 0.04,
      model == "SIR" ~ week_num - 0.12,
      model == "SIRD" ~ week_num - 0.04,
      TRUE ~ week_num
    )
  ) %>%
  filter(str_detect(variable, "kappa") | str_detect(variable, "sigma_gamma") |
           str_detect(variable, "sigma_upsilon")) %>% 
  mutate(variable_label = recode(variable,
                                 "kappa" = "kappa",
                                 "sigma_gamma" = "sigma[gamma]^2",
                                 "sigma_upsilon" = "sigma[upsilon]"
  )) %>% 
  mutate(variable_label = factor(variable_label, levels = c("kappa",
                                                            "sigma[gamma]^2",
                                                            "sigma[upsilon]")))




comb_plot_data %>% 
  ggplot() +
  geom_segment(aes(x = q2.5, xend = q97.5, y = y_pos, yend = y_pos, colour = model),
               size = 1.2) +
  scale_y_continuous(
    breaks = sort(unique(plot_data$week_num)),
    labels = levels(plot_data$week)
  ) +
  facet_wrap(~variable_label, scales = "free_x", 
             labeller = labeller(variable_label = label_parsed),
             nrow = 2) +
  labs(y = "Week", x = "", colour = "Model") +
  scale_colour_manual(values = c("prior" = "#E69F00",  # orange
                                 "ASG" = "#56B4E9",  # sky blue
                                 "ASGD" = "#009E73",  # bluish green
                                 "SIR" = "#F0E442",  # yellow
                                 "SIRD" = "#0072B2"   # blue
  )) + 
  theme_bw() +
  theme(legend.position = c(.95,.85),
        strip.text = element_text(size = 12),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank())






