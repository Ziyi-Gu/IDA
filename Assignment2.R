# Q2.(b)
load("dataex2.Rdata")
x <- dataex2$X
r <- dataex2$R
likelihood <- function(miu,x,r){
  like <- sum(r*log(dnorm(x,miu,1.5)) + (1-r)*log(pnorm(x,miu,1.5)))
  return(-like)
}
mle <- optim(par = 1, fn = likelihood, x = x, r = r)$par
mle


# Q4
load('dataex4.Rdata')
ind <- which(is.na(dataex4$Y) == FALSE)
# split the data into observed and missing 
X1 <- dataex4$X[ind]
X2 <- dataex4$X[-ind]
Y <- dataex4$Y[ind]

# E-step: Q(beta|beta(t))
Q <- function(beta,X1,X2,Y,beta_old){
  q <- sum(Y*(beta[1]+X1*beta[2]) - log(1+exp(beta[1]+X1*beta[2]))) + 
    sum((beta[1]+X2*beta[2])*exp(beta_old[1]+X2*beta_old[2])/(1+exp(beta_old[1]+X2*beta_old[2])) - 
    log(1+exp(beta[1]+X2*beta[2])))
  return(-q)
}

# M-step
M <- function(beta0, eps){
  diff <- 1
  beta <- beta0
  while(diff> eps){
    beta_old <- beta
    beta <- optim(par = c(1,-2), fn = Q, X1 = X1, X2 = X2, Y = Y, beta_old = beta_old)$par
    diff <-sum(abs(beta - beta_old))
  }
  
  return(beta)
}
M(c(10,10),0.0005)


# Q5
load("dataex5.Rdata")
EM <- function(y,theta0,eps){
  n <- length(y)
  theta <- theta0
  p <- theta[1]
  lamda <- theta[2]
  mu <- theta[3]
  diff <- 1
  while(diff > eps){
    theta_old <- theta
    # E-step
    pt1 <- p*lamda/y^(lamda+1)
    pt2 <- (1-p)*mu/y^(mu+1)
    pt <- pt1/(pt1+pt2)
    # M-step
    p <- mean(pt)
    lamda <- sum(pt)/sum(pt*log(y))
    mu <- sum(1-pt)/sum((1-pt)*log(y))
    theta <- c(p,lamda,mu)
    diff <- sum(abs(theta-theta_old))
  }
  return(theta)
}
y <- dataex5
res <- EM(y=y, c(0.3,0.3,0.4),0.0001)
p <- res[1]
lamda <- res[2]
mu <- res[3]
p; lamda; mu
hist(y, breaks = "Freedman-Diaconis", main = "Dataex5", freq= F,xlim = c(1,6))
curve(p*lamda/x^(lamda+1) + (1-p)*mu/x^(mu+1), add = TRUE, col="red")

