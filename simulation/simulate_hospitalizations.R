#source('./get_data_functions.R')




get_post_lm_default <- function(value, scale_uwili, region, poly = 1, 
                                is.log = FALSE) {
  if (is.log == TRUE) {value = log(value + 1)}
  y <- value[-1]
  ylag <- value[-length(value)]
  scale_uwili2 <- scale_uwili^2
  x <- scale_uwili[-1]
  x2 <- scale_uwili2[-length(scale_uwili2)]
  
  if (poly == 2) {
    fit <- lm(y ~ x + x2 + ylag)
    
  }
  else {
    fit <- lm(y ~ x + ylag)
  }
  #print("we get to this place")
  X <- model.matrix(fit)
  beta_hat <- as.vector(coef(fit))
  n_k <- nrow(X) - ncol(X)
  
  #print("i'm guessing we get here?")
  s2 <- as.vector((1/n_k)*t(y - X%*%beta_hat)%*%(y - X%*%beta_hat))
  #print(X)
  V_beta <- solve(t(X)%*%X)
  #print("but not here")
  # return(data.frame(region, beta0 = beta_hat[1], beta1 = beta_hat[2], 
  #                   beta2 = beta_hat[3],
  #                   n_k = n_k, s2 = s2, Vbeta0 = V_beta[1,1], Vbeta1 = V_beta[2,2], 
  #                   Vbetad = V_beta[1,2]))
  
  return(list(region = region, beta_hat = beta_hat, V_beta = V_beta, 
              n_k = n_k, s2 = s2))
}


simulate_hospitalizations <- function(x, state = "US", both_flu, 
                                      is.log = FALSE, rep) {
  
  for_posts <- both_flu %>%
    group_by(season) %>%
    filter(location_name == state,
           season == 2022 | (season == 2021 & season_week == max(season_week))) %>%
    mutate(x = unweighted_ili)
  
  if (is.log == TRUE) {
    posterior <- get_post_lm_default(for_posts$value,
                                for_posts$x, unique(for_posts$region),
                                is.log = is.log)
    print(posterior) 
    sig_draw <- LaplacesDemon::rinvchisq(1, posterior$n_k, posterior$s2)
    beta_draw <- mixtools::rmvnorm(1, mu = posterior$beta_hat, 
                                   sig_draw*posterior$V_beta)
    
    #sig_draw <- (posterior$n_k*posterior$s2)/(2 + posterior$n_k)
    #beta_draw <- posterior$beta_hat
    write.csv(cbind(betas = beta_draw, sig = sqrt(sig_draw)), paste0("test_param_rep", rep, ".csv"), row.names =FALSE)
    y <- c()
    y[1] <- beta_draw[1] + beta_draw[2]*x[1] + #beta_draw[3]*x[1]^2 +
      rnorm(1, 0, sqrt(sig_draw))
    for (i in 2:length(x)) {
      y[i] <- beta_draw[1] + beta_draw[2]*x[i] + #beta_draw[3]*x[i]^2 +
        beta_draw[3]*y[i - 1] +
        rnorm(1, 0, sqrt(sig_draw))
    }
    
    return(exp(y) - 1)
  }
  else {
    posterior <- get_post_lm_default(for_posts$value,
                                     for_posts$x, unique(for_posts$region),
                                     is.log = is.log)
     
    sig_draw <- LaplacesDemon::rinvchisq(1, posterior$n_k, posterior$s2)
    
    beta_draw <- mixtools::rmvnorm(1, mu = posterior$beta_hat, 
                                   sig_draw*posterior$V_beta) 
    y <- c()
    y[1] <- beta_draw[1] + beta_draw[2]*x[1] + beta_draw[3]*x[1]^2 +
      rnorm(1, 0, sqrt(sig_draw))
    for (i in 2:length(x)) {
      y[i] <- beta_draw[1] + beta_draw[2]*x[i] + beta_draw[3]*x[i]^2 +
        beta_draw[4]*y[i - 1] +
        rnorm(1, 0, sqrt(sig_draw))
    }
    
    return(y)
  }
  
  
  
  
}


get_post_lm <- function(value, scale_uwili, region, is.log = FALSE,
                        m0 = rep(0,2), C0 = diag(rep(1,2)), v0 = 2, s02 = 2) {
  if (is.log == TRUE) {value = log(value + 1)}
  p <- length(m0)
  n <- length(value)
  I0 <- diag(rep(1, p))
  In <- diag(rep(1, n))
  
  fit <- lm(value ~ scale_uwili)
  X <- model.matrix(fit)
  # beta_hat <- as.vector(coef(fit))
  # n_k <- nrow(X) - ncol(X)
  # s2 <- as.vector((1/n_k)*t(value - X%*%beta_hat)%*%(value - X%*%beta_hat))
  # V_beta <- solve(t(X)%*%X)
  
  mn <- m0 + C0%*%t(X)%*%solve(X%*%C0%*%t(X) + 
                                In)%*%(value - X%*%m0)
  
  Cn <- C0 - C0%*%t(X)%*%solve(X%*%C0%*%t(X) + 
                                In)%*%X%*%C0
  
  vn <- v0 + n
  vnsn2 = v0*s02 + 
    t(value - X%*%m0)%*%solve(X%*%C0%*%t(X) + In)%*%(value - X%*%m0)
  
  sn2 <- vnsn2/vn
  
  
  return(data.frame(region, beta0 = mn[1], beta1 = mn[2], 
                    n_k = vn, s2 = sn2, Vbeta0 = Cn[1,1], Vbeta1 = Cn[2,2], 
                    Vbetad = Cn[1,2]))
}



post_simulate <- function(df, is.log = FALSE) {
  post <- get_post_lm_default
  N <- nrow(df)
  x <- (df$unweighted_ili/100)*df$count_rate2
  ests <- unique(df[, c("beta0", "beta1", "n_k", "s2", "Vbeta0", "Vbeta1",
                        "Vbetad")])
  
  Vbeta <- matrix(c(ests$Vbeta0, ests$Vbetad, ests$Vbetad, ests$Vbeta1), 
                  nrow = 2,
                  ncol = 2)
  
  
  sigma2 <- LaplacesDemon::rinvchisq(1, ests$n_k, ests$s2)
  
  Sigma <- sigma2*Vbeta
  mu <- c(unique(ests$beta0), unique(ests$beta1))
  
  
  beta <- mixtools::rmvnorm(1, mu = mu, sigma = Sigma)
  
  
  
  y <- beta[1] + beta[2]*x + rnorm(N, sd = sqrt(sigma2))
  # mixtools::rmvnorm(1, rep(0,N), 500*matrix(.8, nrow = N, ncol = N) +
  #                     200*diag(rep(.2, N)))#rnorm(N, sd = sqrt(sigma2))
  
  if (is.log == FALSE) {
    y[y < 0] <- 0
  } else {
    y <- exp(y)
  }
  
  
  return(data.frame(region = unique(df$region), season = unique(df$season),
                    season_week = df$season_week, sim_value = as.vector(y)))
  
}



simulate_post_pred_lm <- function(df, is.log = FALSE) {
  N <- nrow(df)
  x <- (df$unweighted_ili/100)*df$count_rate2
  ests <- unique(df[, c("beta0", "beta1", "n_k", "s2", "Vbeta0", "Vbeta1",
                        "Vbetad")])
  
  Vbeta <- matrix(c(ests$Vbeta0, ests$Vbetad, ests$Vbetad, ests$Vbeta1), 
                  nrow = 2,
                  ncol = 2)
  
  
  sigma2 <- LaplacesDemon::rinvchisq(1, ests$n_k, ests$s2)
  
  Sigma <- sigma2*Vbeta
  mu <- c(unique(ests$beta0), unique(ests$beta1))
 
  
  beta <- mixtools::rmvnorm(1, mu = mu, sigma = Sigma)
  
  
  
  y <- beta[1] + beta[2]*x + rnorm(N, sd = sqrt(sigma2))
    # mixtools::rmvnorm(1, rep(0,N), 500*matrix(.8, nrow = N, ncol = N) +
    #                     200*diag(rep(.2, N)))#rnorm(N, sd = sqrt(sigma2))
  
  if (is.log == FALSE) {
    y[y < 0] <- 0
  } else {
    y <- exp(y)
  }
  
  
  return(data.frame(region = unique(df$region), season = unique(df$season),
                    season_week = df$season_week, sim_value = as.vector(y)))
  
}




get_posterior <- function(region, season, value, unweighted_ili) {
  
  value <- log(value + 1)
  region <- unique(region)
  season <- unique(season)
  lili <- boot::logit(unweighted_ili/100)
  dat <- data.frame(value, lili)
  mod <- rstanarm::stan_glm(value ~ lili, data = dat)
  b0 <- mod$coefficients[1]
  b1 <- mod$coefficients[2]
  b0sd <- summary(mod)[1,3]
  b1sd <- summary(mod)[2,3]
  sigma <- summary(mod)[3,1]
  sigmasd <- summary(mod)[3,3]
  
  dist <- data.frame(region = region, season = season, b0 = b0, b1 = b1, 
                     b0sd = b0sd, b1sd = b1sd, sigma = sigma, 
                     sigmasd = sigmasd)
  
  return(dist)
  
}


simulate_post_pred <- function(df) {
  N <- nrow(df)
  x <- boot::logit(df$unweighted_ili/100)
  ests <- unique(df[, c("b0", "b1", "b0sd", "b1sd", "sigma", "sigmasd")])
  b0 <- ests$b0
  b1 <- ests$b1
  b0sd <- ests$b0sd
  b1sd <- ests$b1sd
  sigma <- ests$sigma
  sigmasd <- ests$sigmasd
  
  y <- rnorm(1, b0, b0sd) + rnorm(1, b1, b1sd)*x + 
    rnorm(N, 0, rnorm(N, sigma, sigmasd))
  
  
  return(data.frame(region = unique(df$region), 
                    season_week = df$season_week, sim_value = exp(y)))
}

get_simulated_data <- function(ILINet, s_season = 2022, posts,
                               select_regions = c("Alabama", 
                                                  "District of Columbia", 
                                                  "Iowa","Montana",
                                                  "Ohio", "Texas"),
                               is.log = FALSE) {
  
  posts <- posts %>% 
    filter(region %in% select_regions)
  
  
  flu_data_pp <- ILINet %>% 
    filter(region %in% select_regions) %>% 
    left_join(posts, by = c("region")) %>% 
    select(region, unweighted_ili, season, season_week, count_rate1, 
           count_rate2, count_rate3, count_rate4, count_rate5, beta0,
           beta1, n_k, s2, Vbeta0, Vbeta1, Vbetad)
  
  regions <- unique(flu_data_pp$region)
  seasons <- unique(flu_data_pp$season)
  
  sim_dat <- data.frame()
  for (i in 1:length(regions)) {
    for (j in 1:length(seasons)) {
      sim_part <- flu_data_pp %>%
        filter(region == regions[i], season == seasons[j]) %>%
        simulate_post_pred_lm(is.log = is.log)
      
      sim_dat <- rbind(sim_dat, sim_part)
    }
  }
  
  
  
  # data <- cbind(flu_data_pp, sim_value = sim_dat$sim_value) %>% 
  #   filter(season != s_season)
  
  data <- flu_data_pp %>% 
    left_join(sim_dat, by = c("region", "season", "season_week")) 
  
  return(data)
  
  
}


