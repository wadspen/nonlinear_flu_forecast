library(deSolve)
library(bbmle)



####################################################
####################SIR MLE Stuff###################
####################################################

sir <- function(t,x,parms){
	S <- x[1]
	I <- x[2]
	R <- x[3]
	with(as.list(parms),
	     {
	       dS <- -beta*S*I
	       dI <- beta*S*I - gamma*I
	       dR <- gamma*I
	       res <- c(dS,dI,dR)
	       list(res)
	     })
}

sirLL <- function(lbeta, lgamma, S = .9, logI0, logkappa, season_week,
		  unweighted_ili) {
	parms <- c(beta = exp(lbeta), gamma = exp(lgamma))
	I <- exp(logI0)
	R <- 1 - (S + exp(logI0))
	x0 <- c(S = S, I = I, R = R)
	out <- ode(y = x0, season_week, sir, parms)
	bdens <- dbetaprop(y = unweighted_ili, mu = out[,3], kappa = exp(logkappa))
	-sum(log(bdens))
}

get_sir_mle <- function(season_week, unweighted_ili, lbeta = log(1.1),
			lgamma = log(.2), logI0 = log(.05),
			logkappa = log(4), S = .9, mit = 1e5) {
	#print("here baby!")
	fit <- mle2(sirLL,
		    start = list(lbeta = lbeta,
				 lgamma = lgamma,
				 logI0 = logI0,
				 logkappa = logkappa),
		    data = list(season_week = season_week,
				unweighted_ili = unweighted_ili,
				S = S),
		    method = "Nelder-Mead",
		    control = list(maxit = mit, trace = 0),
		    trace = FALSE)
	
        #print("finally here!")
	params <- coef(fit)
	beta <- exp(params[1])
	gamma <- exp(params[2])
	I0 <- exp(params[3])
	kappa <- exp(params[4])

	S0 <- S
	R0 <- 1 - (S0 + I0)
	
	N <- 1
	parms <- c(N = N, beta = as.numeric(beta), gamma = as.numeric(gamma))
	x0 <- c(S0, I0, R0)
	stateMatrix <- ode(y = x0, season_week, sir, parms)
	colnames(stateMatrix) <- c("season_week", "S", "I", "R")
        print("totally here!")
	return(stateMatrix[,"I"])
}


###############################################################
#######################ASG MLE Stuff###########################
###############################################################






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
  #print(c(beta1,eta,mu,sigma1,sigma2,kappa))
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
  #print(c(beta1,eta,mu,sigma1,sigma2,kappa))
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
  #print(c(beta,eta,mu,sigma,kappa))
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











