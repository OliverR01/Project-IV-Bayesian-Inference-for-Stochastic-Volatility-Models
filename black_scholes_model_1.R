# Black-Scholes Model Code File 1

# This file simulates synthetic data and uses Metropolis-Hastings and Gibbs Schemes
# to perform inference on parameters mu and sigma under the Black-Scholes model

################################################################################
## Libraries:

library(MASS)
library(latex2exp)
library(scales)
library(coda)

################################################################################
## Inputs Required:

# No prior inputs required

################################################################################
## Functions Required:

# black_scholes_model_1.R
# - data_gen ---- Simulate Synthetic Data
# - logprior ---- Calculate Log-Prior Density
# - loglikeA ---- Calculate Log-Likelihood (Analytic Scheme)
# - loglikeEM --- Calculate Log-Likelihood (Approximate Scheme)
# - logpostA ---- Calculate Log-Posterior Density (Analytic Scheme)
# - logpostEM --- Calculate Log-Posterior Density (Analytic Scheme)
# - MHA --------- Metropolis-Hastings Scheme (Analytic)
# - MHtheta ----- Metropolis-Hastings Scheme (Approximate)
# - MHx --------- Metropolis-Hastings Scheme (Imputations)
# - gibbs ------- Gibbs Sampler for Parameter and Augmented Data updates
# - initSVf ----- Initial Value Function

################################################################################
## New Functions:

data_gen <- function(T, deltat, S0, mu, sigma){
  n <- floor(T/deltat)
  Svec <- rep(NA, n+1)
  Svec[1] <- S0
  for (i in 1:n){
    Svec[i+1] <- rlnorm(1,
                        log(Svec[i]) + (mu-0.5*(sigma^2))*deltat,
                        sigma*sqrt(deltat))
  }
  return(Svec)
}

logprior <- function(theta,b){
  lprior <- dnorm(theta[1],0,b[1],log=TRUE) + dnorm(theta[2],0,b[2],log=TRUE)
  return(lprior)
}

loglikeA <- function(y,theta,deltat){
  n <- length(y)-1
  theta1 <- theta[1]; theta2 <- theta[2]
  Xvec <- log(y)
  mu <- exp(theta1); sigma <- exp(theta2)
  m <- (mu-sigma^2/2)*deltat
  s <- sigma*sqrt(deltat)
  diff <- Xvec[2:(n+1)] - Xvec[1:n]
  return(sum(dnorm(diff,m,s,log=TRUE)))
}

loglikeEM <- function(y,theta,deltat){
  n <- length(y)-1
  mu <- exp(theta[1]); sigma <- exp(theta[2])
  y1 <- y[1:n]
  y2 <- y[2:(n+1)]
  m <- y1 + mu*y1*deltat
  s <- sigma*y1*sqrt(deltat)
  llike <- dnorm(y2, m, s, log = TRUE)
  return(sum(llike))
}

logpostA <- function(y,theta,b,deltat){
  lpost <- logprior(theta,b) + loglikeA(y,theta,deltat)
  return(lpost)
}

logpostEM <- function(y,theta,b,deltat){
  lpost <- logprior(theta,b) + loglikeEM(y,theta,deltat)
  return(lpost)
}

MHA <- function(N, y, deltat, init, V, b){
  theta.mat <- matrix(NA, nrow = N, ncol = 2)
  # Initialise
  theta.mat[1,] <- init
  
  curr <- init
  for (j in 2:N){
    # Propose
    can <- mvrnorm(1,curr,V)
    # Acceptance Probability
    laprob <- logpostA(y,can,b,deltat) - logpostA(y,curr,b,deltat)
    # Acceptance Stage
    if (log(runif(1)) < laprob){
      curr <- can
    }
    theta.mat[j,] <- curr
  }
  return(theta.mat)
}

MHtheta <- function(N, y, deltat, init, V, b){
  theta.mat <- matrix(NA, nrow = N, ncol = 2)
  # Initialise
  theta.mat[1,] <- init
  
  curr <- init
  for (j in 2:N){
    # Propose
    can <- mvrnorm(1,curr,V)
    # Acceptance Probability
    laprob <- logpostEM(y,can,b,deltat) - logpostEM(y,curr,b,deltat)
    # Acceptance Stage
    if (log(runif(1)) < laprob){
      curr <- can
    }
    theta.mat[j,] <- curr
  }
  return(theta.mat)
}

MHx <- function(N, initS, Pvec, m){
  nn <- length(initS)
  
  Smat <- matrix(NA, nrow = N, ncol = nn)
  Smat[1,] <- initS
  
  deltatau <- 1/m
  Tvec <- c(exp(Pvec[1]), exp(Pvec[2]))
  
  for (j in 1:nn){
    if (j%%m==1){
      Smat[,j] <- initS[j]
    }
  }
  
  for (k in 2:N){
    curr <- Smat[(k-1),]
    can <- Smat[k,]
    for (i in 2:(nn-1)){
      if (i%%m!=1){
        currval <- curr[i]
        canval <- rnorm(1, 0.5*(can[i-1]+curr[i+1]), can[i-1]*Tvec[2]*sqrt(deltatau))
        
        num1 <- dnorm(canval, can[i-1] + can[i-1]*Tvec[1]*deltatau, can[i-1]*Tvec[2]*sqrt(deltatau), log = TRUE)
        den1 <- dnorm(currval, can[i-1] + can[i-1]*Tvec[1]*deltatau, can[i-1]*Tvec[2]*sqrt(deltatau), log = TRUE)
        num2 <- dnorm(curr[i+1], canval + canval*Tvec[1]*deltatau, canval*Tvec[2]*sqrt(deltatau), log = TRUE)
        den2 <- dnorm(curr[i+1], currval + currval*Tvec[1]*deltatau, currval*Tvec[2]*sqrt(deltatau), log = TRUE)
        numq <- dnorm(currval, 0.5*(can[i-1]+curr[i+1]), can[i-1]*Tvec[2]*sqrt(deltatau), log = TRUE)
        denq <- dnorm(canval, 0.5*(can[i-1]+curr[i+1]), can[i-1]*Tvec[2]*sqrt(deltatau), log = TRUE)
        
        laprob <- num1+num2+numq-den1-den2-denq
        if(log(runif(1))<laprob){
          currval <- canval
        }
        can[i] <- currval
      }
    }
    Smat[k,] <- can
  }
  return(Smat)
}

gibbs <- function(N,initP,initS,m,V,b){
  nn <- length(initS)
  
  Pmat <- matrix(NA, nrow = N, ncol = 2)
  Smat <- matrix(NA, nrow = N, ncol = nn)
  
  Pmat[1,] <- initP 
  Smat[1,] <- initS
  
  deltatau <- 1/m
  for (s in 2:N){
    data <- Smat[s-1,]
    
    # Update Theta using MHtheta
    Pmat[s,] <- MHtheta(2,data,deltatau,Pmat[s-1,],V,b)[2,]
    
    # Update S_{t} using MHx
    Smat[s,] <- MHx(2,data,Pmat[s,],m)[2,]
  }
  return(list(P = Pmat,
              S = Smat))
}

initSVf <- function(Sobs, Vobs, m){  
  T <- length(Sobs)-1
  nn <- (m*T)+1
  initSV <- matrix(NA, nrow = nn, ncol = 2)
  for (i in 1:T){
    initSV[(1+(i-1)*m):(i*m),1] <- seq(Sobs[i],Sobs[(i+1)],length.out = m+1)[1:m]
  }
  initSV[nn,1] <- Sobs[(T+1)]
  if (length(Vobs)==1){
    if (is.na(Vobs)){
      return(initSV[,1])    
    }
    initSV[,2] <- Vobs
  }
  if (length(Vobs)==2){
    initSV[,2] <- seq(Vobs[1],Vobs[2],length=nn)
  }
  if (length(Vobs) == length(Sobs)){
    for (i in 1:T){
      initSV[(1+(i-1)*m):(i*m),2] <- seq(Vobs[i],Vobs[(i+1)],length.out = m+1)[1:m]
    }
    initSV[nn,2] <- Vobs[(T+1)]
  }
  return(initSV)
}

################################################################################
## Output:

# Data Generation

S0 <- 1
mu <- 0.05
sigma <- 0.1
deltat <- 0.01
set.seed(1)
x <- data_gen(100,deltat,S0,mu,sigma)
par(mfrow=c(1,1))
plot(ts(x,start=0,deltat=deltat), col = 1)

Sobs <- x[seq(1, 10001, by = 100)]
points(seq(0,100,by=1), Sobs, lwd = 2, pch = 4, cex = 0.5, col = 4)

# Inference: Analytic

b <- c(5,5)
initP <- c(0,0)

ll <- 2.42/sqrt(length(Sobs))
Sigma <- diag(rep(ll/sqrt(2), 2))
system.time(ALout1 <- MHA(100000, Sobs, 1, initP, Sigma, b))

postvar <- var(ALout1)
lambda <- 0.8
Sigmatune <- lambda*(2.42^2)*postvar/2
system.time(ALout1.2 <- MHA(1000000, Sobs, 1, initP, Sigmatune, b))
# CPU time approx: 125s
length(unique(ALout1.2[,1]))/1000000
# Acceptance rate target: approx 35%

par(mfrow = c(2,3))
plot(ts(exp(ALout1.2[-(1:200000),1])), main = 'Trace Plot', ylab = TeX(r'($\mu$)'))
acf(exp(ALout1.2[-(1:200000),1]), main = 'ACF Plot')
mudensA <- density(exp(ALout1.2[-(1:200000),1]))
plot(mudensA, main = 'Posterior Density Plot', xlab = TeX(r'($\mu$)'))
abline(v=quantile(exp(ALout1.2[-(1:200000),1]), c(0.025, 0.975)),col=c(3,3))
abline(v=mudensA$x[which.max(mudensA$y)],col = 2)
points(mudensA$x[which.max(mudensA$y)], mudensA$y[which.max(mudensA$y)], col = 2, lwd = 2, pch = 13, cex = 2)
abline(v=mu, col = 4, lty = 2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds", "True value"),
       col=c(2, 3, 4), 
       lty=c(1,1,2), 
       bty = "n",
       cex=0.65)

plot(ts(exp(ALout1.2[-(1:200000),2])), main = 'Trace Plot', ylab = TeX(r'($\sigma$)'))
acf(exp(ALout1.2[-(1:200000),2]), main = 'ACF Plot')
sigmadensA <- density(exp(ALout1.2[-(1:200000),2]))
plot(sigmadensA, main = 'Posterior Density Plot', xlab = TeX(r'($\sigma$)'))
abline(v=quantile(exp(ALout1.2[-(1:200000),2]), c(0.025, 0.975)),col=c(3,3))
abline(v=sigmadensA$x[which.max(sigmadensA$y)],col = 2)
points(sigmadensA$x[which.max(sigmadensA$y)], sigmadensA$y[which.max(sigmadensA$y)], col = 2, lwd = 2, pch = 13, cex = 2)
abline(v=sigma,col=4,lty=2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds", "True value"),
       col=c(2, 3, 4), 
       lty=c(1,1,2), 
       bty = "n",
       cex=0.65)

effectiveSize(exp(ALout1.2[-(1:200000),]))

# Inference: Approximate

ll <- 2.42/sqrt(length(Sobs))
Sigma <- diag(rep(ll/sqrt(2), 2))
system.time(EMout1 <- MHtheta(100000, Sobs, 1, initP, Sigma, b))

postvar <- var(EMout1)
lambda <- 0.7
Sigmatune <- lambda*(2.42^2)*postvar/2
system.time(EMout1.2 <- MHtheta(1000000, Sobs, 1, initP, Sigmatune, b))
# CPU time approx: 120s
length(unique(EMout1.2[,1]))/1000000
# Acceptance rate target: approx 35% 

par(mfrow = c(2,3))
plot(ts(exp(EMout1.2[-(1:200000),1])), main = 'Trace Plot', ylab = TeX(r'($\mu$)'))
acf(exp(EMout1.2[-(1:200000),1]), main = 'ACF Plot')
mudensEM <- density(exp(EMout1.2[-(1:200000),1]))
plot(mudensEM, main = 'Posterior Density Plot', xlab = TeX(r'($\mu$)'))
abline(v=quantile(exp(EMout1.2[-(1:200000),1]), c(0.025, 0.975)),col=c(3,3))
abline(v=mudensEM$x[which.max(mudensEM$y)],col = 2)
points(mudensEM$x[which.max(mudensEM$y)], mudensEM$y[which.max(mudensEM$y)], col = 2, lwd = 2, pch = 13, cex = 2)
abline(v=mu, col = 4, lty = 2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds", "True value"),
       col=c(2, 3, 4), 
       lty=c(1,1,2), 
       bty = "n",
       cex=0.65)

plot(ts(exp(EMout1.2[-(1:200000),2])), main = 'Trace Plot', ylab = TeX(r'($\sigma$)'))
acf(exp(EMout1.2[-(1:200000),2]), main = 'ACF Plot')
sigmadensEM <- density(exp(EMout1.2[-(1:200000),2]))
plot(sigmadensEM, main = 'Posterior Density Plot', xlab = TeX(r'($\sigma$)'))
abline(v=quantile(exp(EMout1.2[-(1:200000),2]), c(0.025, 0.975)),col=c(3,3))
abline(v=sigmadensEM$x[which.max(sigmadensEM$y)],col = 2)
points(sigmadensEM$x[which.max(sigmadensEM$y)], sigmadensEM$y[which.max(sigmadensEM$y)], col = 2, lwd = 2, pch = 13, cex = 2)
abline(v=sigma,col=4,lty=2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds", "True value"),
       col=c(2, 3, 4), 
       lty=c(1,1,2), 
       bty = "n",
       cex=0.65)

effectiveSize(exp(EMout1.2[-(1:200000),]))

# Inference: Approximate with Imputations

m <- 10
initS <- initSVf(Sobs, NA, m)
initP <- c(-1,-1)
plot(ts(x,start=0,deltat=deltat))
lines(ts(Sobs, start = 0, deltat=1), col = 2, lwd=2)

ll <- 2.42/sqrt(length(initS))
Sigma <- diag(rep(ll/sqrt(2), 2))
system.time(IMPout1 <- gibbs(10000, initP, initS, m, Sigma, b))

postvar <- var(IMPout1$P)
lambda <- 0.31
Sigmatune <- lambda*(2.42^2)*postvar/2
system.time(IMPout1.2 <- gibbs(1000, initP, initS, m, Sigmatune, b))
# CPU time: 1500s
length(unique(IMPout1.2$P))/100000
# Acceptance rate target: approx 35%

par(mfrow = c(2,3))
plot(ts(exp(IMPout1.2$P[-(1:20000),1])), main = 'Trace Plot', ylab = TeX(r'($\mu$)'))
acf(exp(IMPout1.2$P[-(1:20000),1]), main = 'ACF Plot')
mudensIMP <- density(exp(IMPout1.2$P[-(1:20000),1]))
plot(mudensIMP, main = 'Posterior Density Plot', xlab = TeX(r'($\mu$)'))
abline(v=quantile(exp(IMPout1.2$P[-(1:20000),1]), c(0.025, 0.975)),col=c(3,3))
abline(v=mudensIMP$x[which.max(mudensIMP$y)],col = 2)
points(mudensIMP$x[which.max(mudensIMP$y)], mudensIMP$y[which.max(mudensIMP$y)], col = 2, lwd = 2, pch = 13, cex = 2)
abline(v=mu, col = 4, lty = 2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds", "True value"),
       col=c(2, 3, 4), 
       lty=c(1,1,2), 
       bty = "n",
       cex=0.65)

plot(ts(exp(IMPout1.2$P[-(1:20000),2])), main = 'Trace Plot', ylab = TeX(r'($\sigma$)'))
acf(exp(IMPout1.2$P[-(1:20000),2]), main = 'ACF Plot')
sigmadensIMP <- density(exp(IMPout1.2$P[-(1:20000),2]))
plot(sigmadensIMP, main = 'Posterior Density Plot', xlab = TeX(r'($\sigma$)'))
abline(v=quantile(exp(IMPout1.2$P[-(1:20000),2]), c(0.025, 0.975)),col=c(3,3))
abline(v=sigmadensIMP$x[which.max(sigmadensIMP$y)],col = 2)
points(sigmadensIMP$x[which.max(sigmadensIMP$y)], sigmadensIMP$y[which.max(sigmadensIMP$y)], col = 2, lwd = 2, pch = 13, cex = 2)
abline(v=sigma,col=4,lty=2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds", "True value"),
       col=c(2, 3, 4), 
       lty=c(1,1,2), 
       bty = "n",
       cex=0.65)

par(mfrow=c(1,1))
plot(ts(x, start = 0, frequency = 100), main = 'Every 100th iteration of the Augmented Dataset', ylab = TeX(r'($S_{t}$)'), ylim = c(0, 55), xlim = c(0,100))
abline(v=seq(0,100,by=1), col = alpha(1,0.2))
for (i in 1:800){
  lines(ts(IMPout1.2$S[(100*i)+20000,], start = 0, frequency = m), col = i)
}

plot(ts(x, start = 0, frequency = 100), xlim = c(20, 25), ylim = c(1.5,2.6), ylab = TeX(r'($S_{t}$)'), main = TeX(r'($\[20,25\]$)'))
abline(v=seq(20,25,by=1), col = alpha(1, 0.2))
for (i in 1:800){
  lines(ts(IMPout1.2$S[(100*i)+20000,], start = 0, frequency = m), col = i)
}
lines(ts(x, start = 0, frequency = 100), lwd = 3, col = 1)

plot(ts(x, start = 0, frequency = 100), xlim = c(80, 90), ylim = c(15,32), ylab = TeX(r'($S_{t}$)'), main = TeX(r'($\[80,90\]$)'))
abline(v=seq(80,90,by=1), col = alpha(1, 0.2))
for (i in 1:800){
  lines(ts(IMPout1.2$S[(100*i)+20000,], start = 0, frequency = m), col = i)
}
lines(ts(x, start = 0, frequency = 100), lwd = 3, col = 1)

plot(ts(x, start = 0, deltat = deltat), ylab = TeX(r'($S_{t}$)'), main = 'MAP and 95% Credible interval of Augmented Dataset' )
abline(v = seq(0, 100, by=1), col = alpha(1,0.2))
SdensIMP <- rep(NA, ncol(IMPout1.2$S))
for (i in 1:ncol(IMPout1.2$S)){
  dens <- density(IMPout1.2$S[-(1:20000),i])
  SdensIMP[i] <- dens$x[which.max(dens$y)]
}
lines(ts(SdensIMP, start = 0, frequency = m), col = 2, lwd = 2)
store <- matrix(NA, nrow = 3, ncol = ncol(IMPout1.2$S))
for (i in 1:ncol(IMPout1.2$S)){
  store[c(1,3),i] <- quantile(IMPout1.2$S[-(1:20000),i], c(0.025, 0.5, 0.975))[c(1,3)]
  store[2,i] <- mean(IMPout1.2$S[-(1:20000),i])
}
lines(ts(store[1,], start = 0, frequency = m), col = 3, lty = 1)
lines(ts(store[3,], start = 0, frequency = m), col = 3, lty = 1)
legend("topleft", 
       c("MAP estimate", "95% Credible bounds", "True value"),
       col=c(2, 3, 1), 
       lty=c(1,1,1),
       cex=1.5)

plot(ts(x, start = 0, deltat = deltat), xlim = c(20, 25), ylim = c(1.5,2.6), main = TeX(r'($\[20,25\]$)'), ylab = TeX(r'($S_{t}$)'))
abline(v=seq(20,25,by=1), col = alpha(1, 0.2))
lines(ts(SdensIMP, start = 0, frequency = m), col = 2, lty=1)
lines(ts(store[1,], start = 0, frequency = m), col = 3, lty = 1)
lines(ts(store[3,], start = 0, frequency = m), col = 3, lty = 1)
legend("topleft", 
       c("MAP estimate", "95% Credible bounds", "True value"),
       col=c(2, 3, 1), 
       lty=c(1,1,1),
       cex=1.5)

plot(ts(x, start = 0, deltat = deltat), xlim = c(80, 90), ylim = c(15,32), main = TeX(r'($\[80,90\]$)'), ylab = TeX(r'($S_{t}$)'))
abline(v=seq(80,90,by=1), col = alpha(1, 0.2))
lines(ts(SdensIMP, start = 0, frequency = m), col = 2, lty=1)
lines(ts(store[1,], start = 0, frequency = m), col = 3, lty = 1)
lines(ts(store[3,], start = 0, frequency = m), col = 3, lty = 1)
legend("topleft", 
       c("MAP estimate", "95% Credible bounds", "True value"),
       col=c(2, 3, 1), 
       lty=c(1,1,1),
       cex=1.5)

effectiveSize(IMPout1.2$P[-(1:20000),])

# Summary

plot(density(exp(ALout1.2[-(1:200000),1])), main = expression(mu), xlim = c(0,0.08), ylim = c(0, 50))
lines(density(exp(EMout1.2[-(1:200000),1])), col = 2)
lines(density(exp(IMPout1.2$P[-(1:20000),1])), col = 3)
abline(v = mudensA$x[which.max(mudensA$y)], lty = 2, col = 1)
abline(v = mudensEM$x[which.max(mudensEM$y)], lty = 2, col = 2)
abline(v = mudensIMP$x[which.max(mudensIMP$y)], lty = 2, col = 3)
abline(v = mu, lty = 2, col = 4)
points(mudensA$x[which.max(mudensA$y)], mudensA$y[which.max(mudensA$y)], col = 1, lwd = 2, pch = 13, cex = 2)
points(mudensEM$x[which.max(mudensEM$y)], mudensEM$y[which.max(mudensEM$y)], col = 2, lwd = 2, pch = 13, cex = 2)
points(mudensIMP$x[which.max(mudensIMP$y)], mudensIMP$y[which.max(mudensIMP$y)], col = 3, lwd = 2, pch = 13, cex = 2)
legend("topleft", 
       c("Analytic Likelihood", "Approximate Likelihood", "Approximate Likelihood with Imputations", "True value"),
       col=c(1, 2, 3, 4), 
       lty=c(1,1,1, 2),
       cex=1.5)

plot(density(exp(ALout1.2[-(1:200000),2])), main = expression(sigma), xlim = c(0.05,0.15), ylim = c(0, 65))
lines(density(exp(EMout1.2[-(1:200000),2])), col = 2)
lines(density(exp(IMPout1.2$P[-(1:20000),2])), col = 3)
abline(v = sigmadensA$x[which.max(sigmadensA$y)], lty = 2, col = 1)
abline(v = sigmadensEM$x[which.max(sigmadensEM$y)], lty = 2, col = 2)
abline(v = sigmadensIMP$x[which.max(sigmadensIMP$y)], lty = 2, col = 3)
abline(v = sigma, lty = 2, col = 4)
points(sigmadensA$x[which.max(sigmadensA$y)], sigmadensA$y[which.max(sigmadensA$y)], col = 1, lwd = 2, pch = 13, cex = 2)
points(sigmadensEM$x[which.max(sigmadensEM$y)], sigmadensEM$y[which.max(sigmadensEM$y)], col = 2, lwd = 2, pch = 13, cex = 2)
points(sigmadensIMP$x[which.max(sigmadensIMP$y)], sigmadensIMP$y[which.max(sigmadensIMP$y)], col = 3, lwd = 2, pch = 13, cex = 2)
legend("topleft", 
       c("Analytic Likelihood", "Approximate Likelihood", "Approximate Likelihood with Imputations", "True value"),
       col=c(1, 2, 3, 4), 
       lty=c(1,1,1, 2),
       cex=1.5)

# Analytic
mean(exp(ALout1.2[-(1:200000),1]))
median(exp(ALout1.2[-(1:200000),1]))
mudensA$x[which.max(mudensA$y)]
quantile(exp(ALout1.2[-(1:200000),1]), c(0.025, 0.975))

mean(exp(ALout1.2[-(1:200000),2]))
median(exp(ALout1.2[-(1:200000),2]))
sigmadensA$x[which.max(sigmadensA$y)]
quantile(exp(ALout1.2[-(1:200000),2]), c(0.025, 0.975))

# Approximate
mean(exp(EMout1.2[-(1:200000),1]))
median(exp(EMout1.2[-(1:200000),1]))
mudensEM$x[which.max(mudensEM$y)]
quantile(exp(EMout1.2[-(1:200000),1]), c(0.025, 0.975))

mean(exp(EMout1.2[-(1:200000),2]))
median(exp(EMout1.2[-(1:200000),2]))
sigmadensEM$x[which.max(sigmadensEM$y)]
quantile(exp(EMout1.2[-(1:200000),2]), c(0.025, 0.975))

# Approx. with Imputations
mean(exp(IMPout1.2$P[-(1:20000),1]))
median(exp(IMPout1.2$P[-(1:20000),1]))
mudensIMP$x[which.max(mudensIMP$y)]
quantile(exp(IMPout1.2$P[-(1:20000),1]), c(0.025, 0.975))

mean(exp(IMPout1.2$P[-(1:20000),2]))
median(exp(IMPout1.2$P[-(1:20000),2]))
sigmadensIMP$x[which.max(sigmadensIMP$y)]
quantile(exp(IMPout1.2$P[-(1:20000),2]), c(0.025, 0.975))