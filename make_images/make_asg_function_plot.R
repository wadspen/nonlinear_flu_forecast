
source("~/flu_research/flu_models/mle_functions.R")
library(ggplot2)
library(deSolve)
library(dplyr)


asg <- function(x,lambda1,eta,mu,sigma1,sigma2) {
  # eta <- lambda1 + exp(eta)
  # eta <- exp(eta)
  asg <- (lambda1 + (eta - lambda1)*exp(-(x - mu)^2/(2*exp(sigma1)^2)))*(x < mu) +
    (lambda1 + (eta - lambda1)*exp(-(x - mu)^2/(2*exp(sigma2)^2)))*(x >= mu)

  return(asg)
}

lambda1 = .2
eta = .42
mu = 24
sigma1 = 1.5 #this will be the log of the standard deviance
sigma2 = 2.3 #also the log

x <- seq(3, 53, length.out = 1001)
y <- asg(x, lambda1 = lambda1, eta = eta, mu = mu,
         sigma1 = sigma1, sigma2 = sigma2)

inf_left <- asg(mu - exp(sigma1), lambda1 = lambda1, eta = eta, mu = mu,
                sigma1 = sigma1, sigma2 = sigma2)

inf_right <- asg(mu + exp(sigma2), lambda1 = lambda1, eta = eta, mu = mu,
                 sigma1 = sigma1, sigma2 = sigma2)

plot(y ~ x, type = "l")

data.frame(x, y) %>%
  ggplot() +
  geom_line(aes(x = x, y = y), size = 1.4) +
  xlim(c(0, 53)) +
  geom_segment(aes(x = mu - exp(sigma1),
                   xend = mu + exp(sigma2),
                   y = inf_left, yend = inf_right), linetype = 2,
               size = 1.2) +

  geom_segment(aes(x = mu, xend = mu, y = lambda1 + .01, yend = eta),
               linetype = 2,
               size = 1.2) +

  geom_segment(aes(x = 7, xend = mu, y = eta, yend = eta),
               linetype = 2,
               size = 1.2) +

  annotate("text", x = mu + exp(sigma2)/2, y = inf_right,
           label = expression(sigma[2]), vjust = 1.6,
           size = 10) +
  annotate("text", x = mu - exp(sigma1)/2, y = inf_left,
           label = expression(sigma[1]), vjust = 1.6,
           size = 10) +

  annotate("text", x = 3, y = eta,
           label = expression(eta + lambda), hjust = .5,
           size = 10) +

  annotate("text", x = 3, y = lambda1,
           label = expression(lambda), hjust = 2,
           size = 10) +

  annotate("text", x = mu, y = .202,
           label = expression(mu),
           size = 10) +
  ylab(expression(f[theta](w))) +
  xlab("w") +
  theme_bw() +
  theme(axis.text=element_text(size=15),
        axis.title = element_text(size=21))

# Code for plotting the SIR compartment trajectories. Here I got some help
# from chatGPT


# Define the SIR model function
sir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I
    dI <- beta * S * I - gamma * I
    dR <- gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

# Set parameters
parameters <- c(beta = 0.85,  # Infection rate
                gamma = 0.5) # Recovery rate

# Set initial values
initial_state <- c(S = 0.9,  # Susceptible
                   I = 0.01,  # Infected
                   R = 0.09)     # Recovered

# Set time vector
times <- seq(0, 200, by = 1)

# Solve the differential equations
out <- ode(y = initial_state, times = times, func = sir_model, parms = parameters)

# Plot the results
# plot(out, main = "SIR Model", xlab = "Time", ylab = "Proportion", type = "l", col = c("blue", "red", "green"), lwd = 2)
# legend("topright", legend = c("Susceptible", "Infected", "Recovered"), col = c("blue", "red", "green"), lty = 1, lwd = 2)



# Convert the output from ode to a data frame
df <- as.data.frame(out)
# Melt the data frame to long format for ggplot
df <- reshape2::melt(df, id.vars = "time")

df %>% 
  ggplot(aes(x = time, y = value)) +
  geom_line(size = 1.4) +
  xlim(c(0, 37)) +
  ylab("Proportion") +
  xlab("t") +
  facet_wrap(~variable, scale = "free_y",) +
  theme_bw() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        strip.text = element_text(size = 17))







