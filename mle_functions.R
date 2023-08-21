

asg <- function(x,beta1,eta,mu,sigma1,sigma2) {
  # eta <- beta1 + exp(eta)
  eta <- exp(eta)
  asg <- (beta1 + (eta)*exp(-(x - mu)^2/(2*exp(sigma1)^2)))*(x < mu) +
    (beta1 + (eta)*exp(-(x - mu)^2/(2*exp(sigma2)^2)))*(x >= mu)
  
  return(asg)
}







lik_asg <- function(y,x,beta1,eta,mu,sigma1,sigma2,kappa=4) {
  asg <- asg(x,beta1,eta,mu,sigma1,sigma2)
  mean <- boot::inv.logit(asg)
  bdens <- dbetaprop(y=y,mu=mean,kappa=kappa)
  lik <- -sum(log(bdens))
  
  return(lik)
}






asg_max_lik <- function(y,x,theta,kappa,M=500,m=0) {
  
  beta1 <- theta[1]
  eta <- theta[2]
  mu <- theta[3]
  sigma1 <- theta[4]
  sigma2 <- theta[5]
  kappa <- theta[6]
  
  repeat{
    
    beta1 <- optim(par=beta1, lik_asg, x=x, y=y, eta = eta,
                   mu = mu, sigma1 = sigma1, sigma2 = sigma2,
                   kappa = kappa)$par
    
  
    
    eta <- optim(par=eta, lik_asg, x=x, y=y, beta1 = beta1,
                 mu = mu, sigma1 = sigma1, sigma2 = sigma2,
                 kappa = kappa)$par
    
    mu <- optim(par=mu, lik_asg, x=x, y=y, beta1 = beta1, 
                eta = eta, sigma1 = sigma1, sigma2 = sigma2,
                kappa = kappa)$par
    
    sigma1 <- optim(par=sigma1, lik_asg, x=x, y=y, beta1 = beta1, 
                    eta = eta, mu = mu, sigma2 = sigma2,
                    kappa = kappa)$par
    
    sigma2 <- optim(par=sigma2, lik_asg, x=x, y=y, beta1 = beta1, 
                    eta = eta, mu = mu, sigma1 = sigma1,
                    kappa = kappa)$par
    
    kappa <- optim(par=kappa, lik_asg, x=x, y=y, beta1 = beta1, 
                   eta = eta, mu = mu, sigma1 = sigma1,
                   sigma2 = sigma2)$par
    
    if (m > M) {break}
    m <- m + 1
    
    
    
  }
  print(c(beta1,eta,mu,sigma1,sigma2,kappa))
  return(c(beta1,eta,mu,sigma1,sigma2,kappa))
  
}








#####################################################
##############try least squares instead##############
#####################################################

ls_asg <- function(y,x,beta1,eta,mu,sigma1,sigma2) {
  asg <- asg(x,beta1,eta,mu,sigma1,sigma2)
  mean <- boot::inv.logit(asg)
  sse <- sum((mean - y)^2)
  
  return(sse)
}

asg_min_ls <- function(y,x,theta,M=500,m=0) {
  
  beta1 <- theta[1]
  eta <- theta[2]
  mu <- theta[3]
  sigma1 <- theta[4]
  sigma2 <- theta[5]
  
  repeat{
    
    beta1 <- optim(par=beta1, ls_asg, x=x, y=y, eta = eta,
                   mu = mu, sigma1 = sigma1, sigma2 = sigma2)$par
    
    
    eta <- optim(par=eta, ls_asg, x=x, y=y, beta1 = beta1, 
                 mu = mu, sigma1 = sigma1, sigma2 = sigma2)$par
    
    mu <- optim(par=mu, ls_asg, x=x, y=y, beta1 = beta1, 
                eta = eta, sigma1 = sigma1, sigma2 = sigma2)$par
    
    sigma1 <- optim(par=sigma1, ls_asg, x=x, y=y, beta1 = beta1,  
                    eta = eta, mu = mu, sigma2 = sigma2)$par
    
    sigma2 <- optim(par=sigma2, ls_asg, x=x, y=y, beta1 = beta1,  
                    eta = eta, mu = mu, sigma1 = sigma1)$par
    
    
    if (m > M) {break}
    m <- m + 1
    
    
    
  }
  print(c(beta1,eta,mu,sigma1,sigma2,kappa))
  return(c(beta1 = beta1,eta = eta,
                    mu = mu,sigma1 = sigma1,sigma2 = sigma2,
                    kappa = kappa))
  
}


###########################################################
#############end of least squares trash####################
###########################################################


asg_calc <- function(x,f_theta) {
  
  y_pred <- boot::inv.logit(asg(x,f_theta[1],f_theta[2],f_theta[3],f_theta[4],
                          f_theta[5],f_theta[6]))
  
  return(y_pred)
}

asg_fit_calc <- function(y,x,start = c(-2, -1, 125, 30,30, 1000) ,M=500) {
  beta1 <- start[1]
  eta <- start[2]
  mu <- start[3]
  sigma1 <- start[4]
  sigma2 <- start[5]
  kappa <- start[6]
  
  
  f_theta <- asg_max_lik(y,x,theta = start,M=M)
  f_theta <- unlist(as.vector(f_theta))
  y_pred <- asg_calc(x,f_theta)
  
  return(y_pred)
  
}





######################################################
#########unnormalized likelihood below here###########
######################################################


dbetaprop <- function(y,mu,kappa) {
  alpha <- mu*kappa
  beta <- (1 - mu)*kappa
  num <- (y^(alpha - 1))*((1 - y)^(beta - 1))
  dens <- num/beta(alpha,beta)
  # dens <- dbeta(x,alpha,beta)
  return(dens)
}

uln <- function(x,beta,eta,mu,sigma) {
  uln <- beta + (exp(eta))*exp(-((log(x) - mu)^2)/(2*exp(sigma)^2))
  return(uln)
}



lik_uln <- function(y,x,beta,eta,mu,sigma,kappa) {
  uln <- uln(x,beta,eta,mu,sigma)
  mean <- boot::inv.logit(uln)
  # mean <- uln
  bdens <- dbetaprop(y=y,mu=mean,kappa=kappa)
  lik <- -sum(log(bdens))
  return(lik)
}


uln_max_lik <- function(y,x, theta = c(-3.7, 0.18, 1.2, 1, 1000), 
                        M=500, m=0) {
  
  beta <- theta[1]
  eta <- theta[2]
  mu <- theta[3]
  sigma <- theta[4]
  kappa <- theta[5]
  
  repeat{
    
    beta <- optim(par=beta, lik_uln, x=x, y=y, eta = eta,
                   mu = mu, sigma = sigma, kappa = kappa)$par
    
    
    eta <- optim(par=eta, lik_uln, x=x, y=y, beta = beta,
                 mu = mu, sigma = sigma, kappa = kappa)$par
    
    mu <- optim(par=mu, lik_uln, x=x, y=y, beta = beta,  
                eta = eta, sigma = sigma, kappa = kappa)$par
    
    sigma <- optim(par=sigma, lik_uln, x=x, y=y, beta = beta, 
                    eta = eta, mu = mu, kappa = kappa)$par
    
    
    kappa <- optim(par=kappa, lik_uln, x=x, y=y, beta = beta, 
                   eta = eta, mu = mu, sigma = sigma)$par
    
    if (m > M) {break}
    m <- m + 1
    
    
    
  }
  print(c(beta,eta,mu,sigma,kappa))
  return(c(beta = beta,eta = eta,mu = mu,sigma = sigma,kappa = kappa))
  
}








uln_calc <- function(x,f_theta) {
  
  y_pred <- boot::inv.logit(uln(x,f_theta[1],f_theta[2],f_theta[3],f_theta[4]))
  
  return(y_pred)
}

uln_fit_calc <- function(y,x,start,M=500) {
  beta <- start[1]
  eta <- start[2]
  mu <- start[3]
  sigma <- start[4]
  kappa <- start[5]
  
  
  f_theta <- uln_max_lik(y,x,beta,eta,mu,sigma,kappa,M=M)
  y_pred <- uln_calc(x,f_theta)
  
  return(y_pred)
  
}











