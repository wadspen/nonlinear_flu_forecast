
setwd("~/flu_research/nonlinear_flu_forecast/")
source('~/flu_research/nonlinear_flu_forecast/functions_data/get_data_functions_new.R')
source('~/flu_research/nonlinear_flu_forecast/functions_data/mle_functions.R')



both_flu <- comb_data(lag = 2) %>% 
  mutate(uw_odd = (unweighted_ili/100)/(1 - (unweighted_ili/100)))
ILINet <- get_ili_data()

state <- "US"
sample(unique(ILINet$region), 6)
dat <- get_stan_data(ILINet, both_flu, s_region = "Texas")
stan_dat <- dat[[1]]

mean_ili <- ILINet %>%
  filter(season %in% c(2010:2022)) %>%
  group_by(region, season_week) %>%
  summarise(
    mean_ili = mean(unweighted_ili, na.rm = TRUE)
  ) %>%
  filter(region %in% c("District of Columbia", "Alabama", "Texas",
                       "Ohio", "Iowa", "Montana")) %>%
  mutate(region = factor(region, levels=c("District of Columbia", "Alabama",
                                         "Texas", "Ohio", "Iowa", "Montana")))

ILINet %>% filter(region == state, season >= 2010) %>% 
  ggplot() +
  geom_line(aes(x = season_week, y = unweighted_ili),
            show.legend = FALSE, size = 1) +
  facet_wrap(~season) +
  ylab("ILI %") +
  xlab("Week") +
  # geom_vline(aes(xintercept = 22)) +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))

both_flu %>% filter(region == "US", season >= 2022) %>% 
  ggplot() +
  geom_line(aes(x = season_week, y = value),
            show.legend = FALSE, size = 1.4) +
  facet_wrap(~season) +
  ylab("Hospitalizations") +
  xlab("Week") +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))

ILINet %>%
  filter(region %in% c("Alabama", "District of Columbia", "Texas",
                       "Ohio", "Iowa", "Montana"),
         season %in% c(2010:2023)) %>%
  mutate(region = factor(region, levels=c("Alabama", "District of Columbia",
                                          "Iowa", "Montana",
                                          "Ohio", "Texas"))) %>%
  ggplot() +
  geom_line(aes(x = season_week, y = unweighted_ili, group = season),
            colour = "darkgrey", size = .7) +
  geom_line(data = mean_ili, aes(x = season_week, y = mean_ili), size = 1.4) +
  facet_wrap(~region, scales = "free") +
  ylab("ILI %") +
  xlab("Week") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))

both_flu %>%
  filter(region %in% c("District of Columbia", "Alabama", "Texas",
                       "Ohio", "Iowa", "Montana"), 
         season %in% 2022:2023) %>%
  ggplot() +
  geom_line(aes(x = season_week, y = value, group = season,
                colour = factor(season)), 
            size = 1) +
  facet_wrap(~region, scales = "free") +
  scale_color_manual(values = c("darkgrey", "black")) +
  labs(colour = "Season") +
  ylab("Hospitalizations") +
  xlab("Week") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))


both_flu %>%
  filter(region %in% c("District of Columbia", "Alabama", "Texas",
                       "Ohio", "Iowa", "Montana"), season == 2022) %>%
  ggplot() +
  geom_line(aes(x = season_week, y = sqrt(value)), size = 1.3) +
  facet_wrap(~region, scales = "free") +
  ylab(expression(sqrt(Hospitalizations))) +
  xlab("Week") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))


ILINet %>%
  filter(region %in% c("District of Columbia", "Alabama", "Texas",
                       "Ohio", "Iowa", "Montana"),
         season == 2022) %>%
  mutate(region = factor(region, levels=c("District of Columbia", "Alabama",
                                          "Texas", "Ohio",
                                          "Iowa", "Montana"))) %>%
  ggplot() +
  geom_line(aes(x = season_week, y = unweighted_ili, group = season),
            size = 1.3) +
  facet_wrap(~region, scales = "free") +
  ylab("ILI %") +
  xlab("Week") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))


both_flu %>%
  filter(region %in% c("Alabama", "District of Columbia", "Texas",
                       "Ohio", "Iowa", "Montana"), season %in% 2022:2023) %>%
  mutate(region = factor(region, levels=c("Alabama", "District of Columbia",
                                          "Texas", "Ohio",
                                          "Iowa", "Montana"))) %>%
  ggplot(aes(x = unweighted_ili, y = value)) +
  geom_point() +
  # geom_path() +
  # geom_smooth(method = "lm", se = FALSE) +
  # geom_smooth(se = FALSE) +
  facet_wrap(~region, scales = "free") +
  # ylab(expression(sqrt(Hospitalizations))) +
  ylab("Hospitalizations") +
  xlab("ILI %") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))


both_flu %>%
  filter(region %in% c("District of Columbia", "Alabama", "Texas",
                       "Ohio", "Iowa", "Montana"), season == 2022) %>%
  mutate(region = factor(region, levels=c("District of Columbia", "Alabama",
                                          "Texas", "Ohio",
                                          "Iowa", "Montana")),
         week = factor(season_week)) %>%
  ggplot(aes(x = unweighted_ili, y = value,
             label = week)) +
  
  geom_path(colour = "orange") +
  geom_point() +
  geom_text(vjust = 2, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~region, scales = "free") +
  # ylab(expression(sqrt(Hospitalizations))) +
  ylab("Hospitalizations") +
  xlab("ILI %") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))

both_flu <- both_flu %>% 
  mutate(uwili2 = count_rate2*unweighted_ili^2, 
         cruili = unweighted_ili*count_rate2)

sample_flu <- both_flu %>% filter(region == "Iowa", season == 2022:2023)
y <- sample_flu$value; x = sample_flu$unweighted_ili*sample_flu$count_rate1
x2 <- x^2
mod <- lm(y ~ x)
ressq <- resid(mod)^2

modres <- lm(ressq ~ x)
N <- nobs(modres)
# gmodres <- glance(modres)

get_fit_flu <- function(x, y, square = FALSE) {
  x2 <- x^2
  x3 <- x^3
  y <- log(y + 1)
  if (square == FALSE) {
    mod <- lm(y ~ x + lag(y))
    preds <- predict(mod)
    preds <- c(y[1], preds)
  }
  else if (square == 3) {
    mod <- lm(y ~ x + x2 + x3)
    preds <- predict(mod)
  }
  else {
    mod <- lm(y ~ x + x2 + lag(y))
    preds <- predict(mod)
    preds <- c(y[1], preds)
  }
  print(anova(mod))
  return(exp(preds) - 1)
}



get_lag_fits <- function(x, y, square = FALSE, log = FALSE) {
  x2 <- x^2
  x3 <- x^3
  if (log == TRUE) {y <- log(y + 1)}
  if (square == FALSE) {
    mod <- lm(y ~ lag(y) + x)
    ar1 <- coef(mod)[2]
  }
  # else if (square == 3) {
  #   mod <- lm(y ~ x + x2 + x3)
  #   preds <- predict(mod)
  # }
  else {
    mod <- lm(y ~ lag(y) + x + x2)
    mod <- lm(y ~ lag(y) + x)
    ar1 <- coef(mod)[2]
  }
  # ar1 <- .64
  lag_val <- y - ar1*lag(y)
  print(ar1)
  return(lag_val)
}

both_flu1 <- comb_data(lag = 1)
both_flu2 <- comb_data(lag = 2)

both_flu1 %>% 
  filter(region == "US", season == 2023) %>% 
  ggplot(aes(x = season_week, y = value)) +
  geom_point()

both_flu1 %>% 
  filter(region == "US", season == 2023) %>% 
  ggplot(aes(x = season_week, y = unweighted_ili)) +
  geom_point()

both_flu2 %>% 
  filter(region == "US", season == 2023) %>% 
  ggplot(aes(x = season_week, y = value)) +
  geom_point()

both_flu2 %>% 
  filter(region == "US", season == 2023) %>% 
  ggplot(aes(x = season_week, y = unweighted_ili)) +
  geom_point()

both_flu1 %>% 
  filter(region == "US", season == 2023) %>% 
  ggplot(aes(x  = unweighted_ili, y = value)) +
  geom_point()

both_flu2 %>% 
  filter(region == "US", season == 2022) %>% 
  # filter(season_week != max(season_week)) %>% 
  ggplot(aes(x  = unweighted_ili, y = value)) +
  geom_point()
  
both_test <- both_flu1 %>% 
  filter(region == "US", season == 2023) %>% 
  mutate(lag_val = lag(value), uili2 = unweighted_ili^2) %>% 
  mutate(ar1_val = value- .64*lag_val)
mod <- lm(value ~ lag_val + unweighted_ili, data = both_test)

plot(ar1_val ~ unweighted_ili^2, data = both_test)
plot(value ~ unweighted_ili, data = both_test)


library(latex2exp)
both_flu1 %>% 
  filter(region %in% c("US"),
         # c("Idaho", "Utah", "Oregon", "Minnesota", "California",
         #             "North Carolina", "Georgia"),
         # c("District of Columbia", "Alabama", "Texas",
         #             "Ohio", "Iowa", "Montana"),
         season %in% c(2022:2023)) %>% 
  group_by(region) %>% 
  mutate(lag_val_log = 
           get_lag_fits(unweighted_ili, value, square = FALSE,
                                 log = TRUE)) %>% 
  mutate(lag_val = 
           get_lag_fits(unweighted_ili, value, square = FALSE,
                        log = FALSE)) %>% 
  pivot_longer(cols = contains("lag_val"),
               values_to = "lag_vals",
               names_to = "lag_type") %>% 
  mutate(lag_type = ifelse(lag_type == "lag_val",
                           "Hospitalizations", 
                           "log Hospitalizations")) %>% 
  ggplot() +
  geom_point(aes(y = lag_vals, x = unweighted_ili, 
                 colour = factor(season)), size = 2) +
  # facet_wrap(~season) +
  scale_color_manual(values = c("darkgrey", "black")) +
  xlab("ILI %") +
  ylab(TeX("$H_{s,w} - \\phi H_{s, w-1}")) +
  labs(colour = "Season") +
  facet_wrap(~lag_type, scales = "free_y") +
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=19),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=14),
        strip.text = element_text(size = 14))


both_flu1 %>% 
  filter(region %in% c("US"),
           # c("Idaho", "Utah", "Oregon", "Minnesota", "California",
           #             "North Carolina", "Georgia"),
           # c("District of Columbia", "Alabama", "Texas",
           #             "Ohio", "Iowa", "Montana"),
         season %in% 2023) %>% 
  group_by(region, season) %>% 
  mutate(lm_preds = get_fit_flu(unweighted_ili, value, square = TRUE)) %>% 
  ggplot() +
  # geom_point(aes(y = value - .57*lag(value), x = unweighted_ili)) +
  geom_point(aes(x = season_week, y = value)) +
  geom_line(aes(x = season_week, y = lm_preds), colour = "red") +
  # facet_grid(region~season, scale = "free") + 
  # geom_point(aes(x = lm_preds, y = value)) +
  ylab("Hospitalizations") +
  xlab("Week") +
  facet_wrap(~region, scale = "free") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))


both_flu_errors <- both_flu1 %>% 
  # filter(region %in% c("Idaho", "Utah", "Oregon", "Minnesota", "California",
  #                      "North Carolina", "Georgia", "US")) %>%
  group_by(region, season) %>% 
  mutate(lm_preds = get_fit_flu(unweighted_ili*count_rate1, 
                                value, square =TRUE)) %>% 
  mutate(errors = lm_preds - value) #%>% 
  

# ggplot(aes(x = season_week, y = errors)) +
  # geom_line() +
  # facet_wrap(~region, scale = "free")

state <- "US"
y <- log(both_flu_errors$value[both_flu_errors$region == state & 
                             both_flu_errors$season >= 2022] + 1)

x <- both_flu_errors$unweighted_ili[both_flu_errors$region == state & 
                             both_flu_errors$season >= 2022]

s <- both_flu_errors$season[both_flu_errors$region == state & 
                                      both_flu_errors$season >= 2022]

x2 <- x^2
x3 <- x^3
lx2 <- lag(x)^2
mod <- lm(y ~ x + x2 + lag(y))
dude <- resid(mod)
mod2 <- lm(dude ~ x[-c(1)])
# mod <- lm(y[-length(y)] ~ x[-length(y)] + x2[-length(y)] + dude + lag(y[-length(y)]))
preds <- predict(mod)
plot(preds ~ y[-c(1)])
abline(a = 0, b = 1)

predssq <- preds^2
mod2 <- lm(y[-1] ~ preds + lag(y[-1]))
preds2 <- predict(mod2)
plot(preds2 ~ y[-c(1:2)])

error <- both_flu_errors$errors[both_flu_errors$region == "US" & 
                                      both_flu_errors$season == 2022]

weeks <- both_flu_errors$season_week[both_flu_errors$region == "US" & 
                                  both_flu_errors$season == 2022]


yb <- y[1:(length(y) - 1)]
xb <- x[1:(length(y) - 1)]
errorb <- error[1:(length(y) - 1)]
yf <- y[2:length(y)]
weeksf <- weeks[2:length(y)]

xb2 <- xb^2
xb3 <- xb^3
ucar <- yb - (xb + xb2 + xb3)
mod <- lm(yf ~ xb + xb2 + xb3)
mod2 <- lm(yf ~ ucar + xb + xb2 + xb3)

plot(resid(mod) ~ weeksf, type = "l")

x2 <- x^2
x3 <- x^3
mod <- lm(y ~ x + x2 + lag(y))

library(forecast)
AR <- auto.arima(error[1:16])
ts.plot(error, xlim = c(0,21))
AR_forecast <- predict(AR, n.ahead = 10)$pred
AR_forecast_se <- predict(AR, n.ahead = 10)$se
points(AR_forecast, type = "l", col = 2)
points(AR_forecast - 2*AR_forecast_se, type = "l", col = 2, lty = 2)
points(AR_forecast + 2*AR_forecast_se, type = "l", col = 2, lty = 2)


both_flu %>% 
  filter(region %in% c("District of Columbia", "Alabama", "Texas",
                       "Ohio", "Iowa", "Montana", "Idaho"), season %in% 2023) %>% 
  group_by(region) %>% 
  mutate(lm_preds = get_fit_flu(unweighted_ili*count_rate1, value)) %>% 
  filter(season == 2023) %>% 
  ggplot() +
  # facet_grid(region~season, scale = "free") + 
  geom_point(aes(x = season_week, y = value)) +
  geom_line(aes(x = season_week, y = lm_preds), colour = "red") +
  ylab("Hospitalizations") +
  xlab("Week") +
  facet_wrap(~region, scale = "free") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))


both_flu %>% 
  filter(region %in% c("US"), season %in% 2022:2023) %>% 
  group_by(region) %>% 
  mutate(lm_preds = get_fit_flu(unweighted_ili*count_rate1, value)) %>% 
  filter(season_week < 17) %>%
  ggplot() +
  # facet_grid(region~season, scale = "free") + 
  geom_line(aes(x = season_week, y = lm_preds - value)) +
  # geom_line(aes(x = season_week, y = lm_preds), colour = "red") +
  ylab("Hospitalizations") +
  xlab("Week") +
  facet_wrap(~season, scale = "free") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))




#########################################
######For 75th anniversary poster########
#########################################

mean_ili_us <- ILINet %>%
  filter(season %in% c(2010:2022)) %>%
  group_by(region, season_week) %>%
  summarise(
    mean_ili = mean(unweighted_ili, na.rm = TRUE)
  ) %>%
  filter(region == "US")

ILINet %>%
  filter(region == "Oregon",
         season %in% 2010:2022) %>%
  # mutate(region = factor(region, levels=c("District of Columbia", "Alabama",
  #                                         "Texas", "Ohio",
  #                                         "Iowa", "Montana"))) %>%
  ggplot() +
  geom_line(aes(x = season_week, y = unweighted_ili, group = season,
            colour = factor(season)), size = .4) +
  # geom_line(data = mean_ili_us, aes(x = season_week, y = mean_ili), size = 1.8) +
  # facet_wrap(~region, scales = "free") +
  ylab("ILI %") +
  xlab("Week") +
  theme_bw() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=23), 
        legend.position = "none")


both_flu %>%
  filter(region == "Iowa", season == 2022) %>%
  # mutate(region = factor(region, levels=c("District of Columbia", "Alabama",
  #                                         "Texas", "Ohio",
  #                                         "Iowa", "Montana"))) %>%
  ggplot(aes(x = unweighted_ili, 
             y = value, label = season_week)) +
  geom_point(size = 2.3) +
  geom_text(colour = "orange", vjust = 2) +
  # geom_smooth(se = FALSE, method = "lm") +
  # facet_wrap(~region, scales = "free") +
  # ylab(expression(sqrt(Hospitalizations))) +
  ylab("log(Hospitalizations)") +
  xlab("logit(ILI)") +
  theme_bw() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=23))


both_flu %>%
  filter(region == "US", season == 2022) %>%
  # mutate(region = factor(region, levels=c("District of Columbia", "Alabama",
  #                                         "Texas", "Ohio",
  #                                         "Iowa", "Montana"))) %>%
  ggplot(aes(x = season_week, 
             y = value)) +
  geom_line(size = 1.7) +
  # geom_smooth(se = FALSE, method = "lm") +
  # facet_wrap(~region, scales = "free") +
  # ylab(expression(sqrt(Hospitalizations))) +
  ylab("Hospitalizations") +
  xlab("Week") +
  theme_bw() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=23))



ILINet_us <- ILINet %>%
  filter(region == "US",
         season %in% c(2010:2023)) %>%
  group_by(season) %>% 
  mutate(asg_pred = asg_fit_calc(unweighted_ili/100, season_week,
                                 start = start <- c(-4.2533962, 0.3416948, 22.5023919, 2.3580754, 
                                                    2.3519650, 100.3432742)))
ILINet_select <- ILINet %>%
  filter(region %in% c("District of Columbia", "Alabama",
                     "Texas", "Ohio",
                     "Iowa", "Montana"),
         season %in% c(2010:2022)) %>%
  group_by(region, season) %>% 
  mutate(asg_pred = asg_fit_calc(unweighted_ili/100, season_week,
                                 start = start <- c(-4.2533962, 0.3416948, 
                                                    22.5023919, 2.3580754, 
                                                    2.3519650, 100.3432742)))


#need this one to show best asg fits
ILINet_us %>% 
  filter(season %in% c(2010:2019, 2021:2022)) %>%
  ggplot() + 
  geom_line(aes(x = season_week, y = unweighted_ili),
            colour = "darkgrey", size = 1.2) +
  geom_line(aes(x = season_week, y = asg_pred*100), size = 1.2) +
  facet_wrap(~season, scales = "free_y") +
  ylab("ILI %") +
  xlab("Week") +
  theme_bw() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=23), 
        legend.position = "none")


ILINet_select %>% 
  filter(season == 2018) %>% 
  ggplot() + 
  geom_point(aes(x = season_week, y = unweighted_ili),
             size = 1.9) +
  geom_line(aes(x = season_week, y = asg_pred*100), 
            colour = "blue",
            size = 1.5) +
  ylab("ILI %") +
  xlab("Week") +
  facet_wrap(~region) +
  theme_bw() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=23), 
        legend.position = "none")


ILINet_us <- ILINet_us %>% 
  mutate(disc = unweighted_ili - asg_pred*100)

ILINet_us %>% 
  ggplot() +
  geom_line

ILINet_select <- ILINet_select %>% 
  mutate(disc = unweighted_ili - asg_pred*100)

mean_disc_us <- ILINet_us %>%
  filter(season %in% c(2010:2019,2021:2022)) %>%
  group_by(region, season_week) %>%
  summarise(
    mean_disc = mean(disc, na.rm = TRUE)
  ) 

# mean_disc_select <- ILINet_select %>%
#   filter(season %in% c(2010:2023)) %>%
#   group_by(region, season_week) %>%
#   summarise(
#     mean_disc = mean(disc, na.rm = TRUE)
#   ) 



ILINet_us %>%
  filter(region == "US",
         season %in% c(2010:2019,2021:2022)) %>%
  # mutate(region = factor(region, levels=c("District of Columbia", "Alabama",
  #                                         "Texas", "Ohio",
  #                                         "Iowa", "Montana"))) %>%
  ggplot() +
  geom_line(aes(x = season_week, y = disc, group = season
                # ,colour = factor(season)
                ), 
            colour = "darkgrey",
            size = .8) +
  geom_line(data = mean_disc_us, aes(x = season_week, y = mean_disc), size = 2.4) +
  # facet_wrap(~region, scales = "free") +
  ylab("Discrepancy") +
  xlab("Week") +
  theme_bw() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=23), 
        legend.position = "none")


ILINet_select %>%
  filter(region == "Ohio",
         season %in% c(2010:2022)) %>%
  # mutate(region = factor(region, levels=c("District of Columbia", "Alabama",
  #                                         "Texas", "Ohio",
  #                                         "Iowa", "Montana"))) %>%
  ggplot() +
  geom_line(aes(x = season_week, y = disc, group = season,
                colour = factor(season)), size = .4) +
  geom_line(data = mean_disc_select[mean_disc_select$region == "Ohio",], aes(x = season_week, y = mean_disc), size = 1.8) +
  # facet_wrap(~region, scales = "free") +
  ylab("Discrepancy %") +
  xlab("Week") +
  theme_bw() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=23), 
        legend.position = "none")



ILINet_us %>% 
  ggplot() +
  geom_line(aes(x = season_week, y = disc, group = season))




# 
# fit <- stl(tsus, s.window="periodic",
#            robust=TRUE)
# fit %>% seasadj() %>% 
#   autoplot() + ylab("New orders index") +
#   ggtitle("Naive forecasts of seasonally adjusted data")
# 
# 
# 
# fcast <- stlf(tsus, method='naive')

x <- runif(23, 2, 100)
vars <- 2 + 3006*x + rnorm(length(x), 0, 1000)
for (i in 1:length(y)) {y[i] <- -2 + 7*x[i] + rnorm(1, 0, sqrt(vars[i]))}
y <- -2 + 7*x + rnorm(length(x), 0, sqrt(vars))
plot(y~x)
