library(extraDistr)
library(ggplot2)
post <- readRDS("./flu_research/Prelim Content/US2024-02-03.rds")


hist(post$`kappas[20]`)
hist(post$sigma_gamma)
hist(post$sigma_gamma_w) #N+(sigmahat, 1)
hist(post$`alpha2[2]`)
hist(post$`theta[1]`)


#sigma_gamma_w, note the mean value sigma_gamma_hat is not included though should be
x <- seq(0.0001, 4, length.out = 1001)
y <- dhnorm(x, sigma=1)


post %>% 
  ggplot() +
  geom_histogram(aes(x = sigma_gamma_w, y = after_stat(density)),
                 fill = "lightgrey", colour = "black") +
  geom_line(data = data.frame(x,y), aes(x = x, y = y), col = "red", size = 1) +
  ylab("") +
  xlab(~ sigma[gamma[W]]^2) +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))


#sigma_gamma
x <- seq(0.0001, .015, length.out = 1001)
y <- dhnorm(x, sigma=.02)

post %>% 
  ggplot() +
  geom_histogram(aes(x = `sigma_gamma`, y = after_stat(density)),
                 fill = "lightgrey", colour = "black") +
  geom_line(data = data.frame(x,y), aes(x = x, y = y), col = "red", size = 1) +
  ylab("") +
  xlab(~ sigma[gamma] ^ 2) +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))

#kappa for 2023
x <- seq(0.1, 31000, length.out = 1001)
y <- dhnorm(x, sigma=10000)

post %>% 
  ggplot() +
  geom_histogram(aes(x = `kappas[20]`, y = after_stat(density)),
                 fill = "lightgrey", colour = "black") +
  geom_line(data = data.frame(x,y), aes(x = x, y = y), col = "red", size = 1) +
  ylab("") +
  xlab(~ kappa[23]) +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))



#beta0 2023

x <- seq(-2700, 2700, length.out = 1001)
y <- dnorm(x, 0, 1000)

post %>% 
  ggplot() +
  geom_histogram(aes(x = `alpha0[2]`, y = after_stat(density)),
                 fill = "lightgrey", colour = "black") +
  geom_line(data = data.frame(x,y), aes(x = x, y = y), col = "red", size = 1) +
  ylab("") +
  xlab(~ beta[0[23]]) +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))


#beta1 2023

x <- seq(-2000, 2000, length.out = 1001)
y <- dnorm(x, 0, 1000)

post %>% 
  ggplot() +
  geom_histogram(aes(x = `alpha1[2]`, y = after_stat(density)),
                 fill = "lightgrey", colour = "black") +
  geom_line(data = data.frame(x,y), aes(x = x, y = y), col = "red", size = 1) +
  ylab("") +
  xlab(~ beta[1[23]]) +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))


#beta2 2023

x <- seq(-500, 500, length.out = 1001)
y <- dnorm(x, 0, 1000)

post %>% 
  ggplot() +
  geom_histogram(aes(x = `alpha2[2]`, y = after_stat(density)),
                 fill = "lightgrey", colour = "black") +
  geom_line(data = data.frame(x,y), aes(x = x, y = y), col = "red", size = 1) +
  ylab("") +
  xlab(~ beta[1[23]]) +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))




#phi

x <- seq(-.9, .9, length.out = 1001)
y <- distr::d(distr::Truncate(distr::Norm(0,.5), -1, 1))(x)

post %>% 
  ggplot() +
  geom_histogram(aes(x = `alpha3[2]`, y = after_stat(density)),
                 fill = "lightgrey", colour = "black") +
  geom_line(data = data.frame(x,y), aes(x = x, y = y), col = "red", size = 1) +
  ylab("") +
  xlab(~ phi) +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))




#sigma epsilon 23

x <- seq(.0001, 13, length.out = 1001)
y <- distr::d(distr::Truncate(distr::Norm(0,20), 0, Inf))(x)

post %>% 
  ggplot() +
  geom_histogram(aes(x = `sigma_epsilon[2]`, y = after_stat(density)),
                 fill = "lightgrey", colour = "black") +
  geom_line(data = data.frame(x,y), aes(x = x, y = y), col = "red", size = 1) +
  ylab("") +
  xlab(~ sigma[epsilon[23]]) +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))


#eta

x <- seq(.0001, 1, length.out = 1001)
y <- dnorm(x, .3, .2)

post %>% 
  ggplot() +
  geom_histogram(aes(x = `theta[1]`, y = after_stat(density)),
                 fill = "lightgrey", colour = "black") +
  geom_line(data = data.frame(x,y), aes(x = x, y = y), col = "red", size = 1) +
  ylab("") +
  xlab(~ paste("log(", eta, ")")) +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))


#mu

x <- seq(10, 35, length.out = 1001)
y <- dnorm(x, 23, 5)

post %>% 
  ggplot() +
  geom_histogram(aes(x = `theta[2]`, y = after_stat(density)),
                 fill = "lightgrey", colour = "black") +
  geom_line(data = data.frame(x,y), aes(x = x, y = y), col = "red", size = 1) +
  ylab("") +
  xlab(~ mu) +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))




#sig1

x <- seq(0, 6, length.out = 1001)
y <- dnorm(x, 3.69, 2)

post %>% 
  ggplot() +
  geom_histogram(aes(x = `theta[3]`, y = after_stat(density)),
                 fill = "lightgrey", colour = "black") +
  geom_line(data = data.frame(x,y), aes(x = x, y = y), col = "red", size = 1) +
  ylab("") +
  xlab(~ paste("log(", sigma[1] ^2, ")")) +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))



#sig2

x <- seq(0, 7, length.out = 1001)
y <- dnorm(x, 4.7, 2)

post %>% 
  ggplot() +
  geom_histogram(aes(x = `theta[3]`, y = after_stat(density)),
                 fill = "lightgrey", colour = "black") +
  geom_line(data = data.frame(x,y), aes(x = x, y = y), col = "red", size = 1) +
  ylab("") +
  xlab(~ paste("log(", sigma[2] ^ 2, ")")) +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))



#zeta

x <- seq(0, 12, length.out = 1001)
y <- dhnorm(x, 4)

post %>% 
  ggplot() +
  geom_histogram(aes(x = `zeta[2]`, y = after_stat(density)),
                 fill = "lightgrey", colour = "black") +
  geom_line(data = data.frame(x,y), aes(x = x, y = y), col = "red", size = 1) +
  ylab("") +
  xlab(~ zeta[2]) +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))




























