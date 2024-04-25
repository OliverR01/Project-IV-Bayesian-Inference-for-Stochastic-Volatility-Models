# Heston Model Code File 1 

# This file simulates synthetic data and uses Metropolis-Hastings to perform
# inference on parameters mu, kappa, theta, sigma and rho under the Heston model.
# Augmented Data not yet introduced.

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
# - initSVf ---------- Initial Value Function

# heston_model_1.R
# - deltaSV ---------- Proposes next step in (S,V) process
# - data.generator --- Simulates Synthetic Data
# - lprior ----------- Calculates log-prior density
# - lgpdf ------------ Calculates log of the multivariate normal density
# - llikehood -------- Calculates log-likelihood 
# - lpost ------------ Calculates log-Posterior
# - MH.P ------------- Metropolis-Hastings Scheme for Parameters


################################################################################
## New Functions:

deltaSV <- function(S, V, deltat, Tvec){
  # Takes the initial state of the SV process, and, given a time-step deltat and 
  # parameters in Tvec, gives a simulation of the next step of the process by Euler-Maruyama
  
  mu <- Tvec[1]; kappa <- Tvec[2]; theta <- Tvec[3]; sigma <- Tvec[4]; rho <- Tvec[5]
  
  m <- c(S + mu*S*deltat, V + kappa*(theta - V)*deltat)
  Sig <- V*deltat*matrix(c(S*S, rho*sigma*S, rho*sigma*S, sigma*sigma), nrow = 2, ncol = 2)
  
  return(mvrnorm(1,m,Sig))
}

data.generator <- function(T, deltat, S0, V0, Tvec){
  # Uses deltaSV to simulate a skeleton path upto time T of the SV process given fixed 
  # parameters and an initial state
  
  n <- floor(T/deltat)
  SVmat <- matrix(NA, nrow = n+1, ncol = 2)
  SVmat[1,] <- c(S0, V0)
  
  for (i in 1:n){
    SVmat[(i+1),] <- deltaSV(SVmat[i,1], SVmat[i,2], deltat, Tvec)
    if(SVmat[(i+1),2]<=0){
      # practical fix to keep V_{t} > 0
      SVmat[(i+1),2] <- 0.0001
    }
  }
  return(SVmat)
}

lprior <- function(Pvec, a, b){
  # let a be a vector of means, b a vector of sdevs (only for the first four parameters)
  # Provides log prior for parameters in Pvec= c(log(mu),log(kappa),log(theta),log(sigma),rho)
  # first four parameters have normal prior
  lprior <- rep(NA, length(Pvec))
  for (i in 1:4){
    lprior[i] <-dnorm(Pvec[i],a[i],b[i], log = TRUE)
  }
  # rho has uniform prior
  lprior[5] <- dunif(Pvec[5], -1, 1, log = TRUE)
  return(sum(lprior))
}

lgpdf <- function(x, m, V){
  d=length(m)
  return(-0.5*log(det(V)) -0.5*t(x-m)%*%solve(V)%*%(x-m)-0.5*d*log(2*pi))
}


llikehood <- function(Pvec, data, deltat){
  # data must be a matrix, 2 columns, rows from t=0 to t=N
  # Pvec must be a vector of 5 elements
  mu <- exp(Pvec[1]); kappa <- exp(Pvec[2]); theta <- exp(Pvec[3]); sigma <- exp(Pvec[4]); rho <- Pvec[5]
  n <- nrow(data)-1
  
  llike <- rep(NA, n)
  for (i in 1:n){
    s <- data[i,1]; v <- data[i,2]
    M <- c(s + mu*s*deltat, v + kappa*(theta-v)*deltat)
    Sig <- deltat*v*matrix(c(s*s, rho*sigma*s, rho*sigma*s, sigma*sigma), ncol = 2, nrow = 2)
    llike[i] <- lgpdf(data[(i+1),], M, Sig)
  }
  return(sum(llike))
}

lpost <- function(Pvec, data, deltat, a, b){
  return(lprior(Pvec, a, b) + llikehood(Pvec, data, deltat))
}

MH.P <- function(N, initP, data, deltat, Sigma, a, b){
  Pmat <- matrix(NA, ncol = length(initP), nrow = N)
  Pmat[1,] <- initP
  
  count <- 0
  # Initialise:
  curr <- Pmat[1,]
  for (i in 2:N){
    # Propose:
    can <- mvrnorm(1, curr, Sigma)
    if(can[5] <= 1 & can[5] >= -1){
      # Calculate Acceptance probability:
      laprob <- lpost(can, data, deltat, a, b) - lpost(curr, data, deltat, a, b)
      if(log(runif(1)) < laprob){
        # Accept move:
        curr <- can
        count <- count+1
      }
    }
    Pmat[i,] <- curr
  }
  #print(count/(N-1))
  return(Pmat)
}

################################################################################
## Output:

# Data Generation

# We set the following parameters
T <- 10 # Final time step
n <- 100 # Generated Discretisation
m <- 10 # Desired Discretisation

mu <- 0.4 # Drift constant
kappa <- 0.92 # Mean-reversion Rate
theta <- 0.35 # Long-run Variance
sigma <- 0.18 # Vol. of Vol
rho <- 0 # Correlation
Tvec <- c(mu, kappa, theta, sigma, rho)

S0 <- 1 # Initial Stock Price
V0 <- 2 # Initial Volatility

set.seed(45)
Y <- data.generator(T,1/n,S0,V0,Tvec)
Yobs <- Y[seq(1,(n*T)+1, by = n),]
initSV <- initSVf(Yobs[,1],Yobs[,2],m)
Sobs <- Yobs[,1]
Vobs <- Yobs[,2]
nn <- (m*T)+1

par(mfrow=c(1,2))
plot(ts(Y[,1], start = 0, frequency = n), ylim = range(Y[,1]), xlab = 'Time', main = TeX(r'(The Stock Price Process $S_{t}$)'), ylab = '')
abline(v=seq(0,T,by=1), col = alpha(1, 0.2))
axis(side=1, at=seq(0,10,by=1))
points(seq(0,T,by=1),Yobs[,1], lwd = 3, col = 2, pch = 1, cex = 1.5)
#lines(seq(0,T,by=1/m),initSV[,1],col=2)

plot(ts(Y[,2], start = 0, frequency = n), ylim = range(Y[,2]), col = 1, xlab = 'Time', main = TeX(r'(The Volatility Process $V_{t}$)'), ylab = '')
abline(v=seq(0,T,by=1), col = alpha(2, 0.2))
axis(side=1, at=seq(0,10,by=1))
points(seq(0,T,by=1),Yobs[,2], lwd = 3, col = 2, pch = 1, cex = 1.5)
#abline(h=theta, lty = 2, col=3)
#lines(seq(0,T,by=1/m),initSV[,2],col=2)

Pvec <- c(log(Tvec[1]),log(Tvec[2]),log(Tvec[3]),log(Tvec[4]),Tvec[5])

# Initial Value and Priors
initP <- rep(0,5)
a <- rep(-1.2,4)
b <- rep(0.5,4)

# Prior distribution plots
par(mfrow = c(1,1))
plot(exp(seq(a[1]-2,a[1]+2,length=1000)), dnorm(seq(a[1]-2,a[1]+2,length=1000),a[1],b[1]), 'l', main = TeX(r'(Prior for $\mu$)'), ylab = '', xlab = TeX(r'($\mu$)'))
plot(exp(seq(a[2]-2,a[2]+2,length=1000)), dnorm(seq(a[2]-2,a[2]+2,length=1000),a[2],b[2]), 'l', main = TeX(r'(Prior for $\kappa$)'), ylab = '', xlab = TeX(r'($\kappa)'))
plot(exp(seq(a[3]-2,a[3]+2,length=1000)), dnorm(seq(a[3]-2,a[3]+2,length=1000),a[3],b[3]), 'l', main = TeX(r'(Prior for $\theta$)'), ylab = '', xlab = TeX(r'($\theta$)'))
plot(exp(seq(a[4]-2,a[4]+2,length=1000)), dnorm(seq(a[4]-2,a[4]+2,length=1000),a[4],b[4]), 'l', main = TeX(r'(Prior for $\sigma$)'), ylab = '', xlab = TeX(r'($\sigma$)'))

# MCMC scheme for Pvec
ll <- 2.40/sqrt(11)
Sigma <- diag(rep(ll/sqrt(5), 5))
system.time(output1 <- MH.P(100000,initP,Yobs,1,Sigma,a,b))

postvar <- var(output1)
lambda <- 0.8
Sigmatune <- lambda*(2.40^2)*postvar/5
system.time(output2 <- MH.P(1000000, initP, Yobs, 1, Sigmatune, a, b))
length(unique(output2[,1]))/1000000
# Target AP: 28.39%
# CPU approx: 1500s  

out <- output2
par(mfrow = c(5,3))
plot(ts(exp(out[-(1:200000),1])), main = 'Trace Plot', ylab = TeX(r'($\mu$)'))
acf(exp(out[-(1:200000),1]), main = 'ACF Plot')
mudens <- density(exp(out[-(1:200000),1]))
plot(mudens, main = 'Posterior Density Plot', xlab = TeX(r'($\mu$)'))
abline(v=quantile(exp(out[-(1:200000),1]), c(0.025, 0.975)),col=c(3,3))
abline(v=mudens$x[which.max(mudens$y)],col = 2)
points(mudens$x[which.max(mudens$y)], mudens$y[which.max(mudens$y)], col = 2, lwd = 2, pch = 13, cex = 2)
abline(v=mu, col = 4, lty = 2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds", "True value"),
       col=c(2, 3, 4), 
       lty=c(1,1,2), 
       bty = "n",
       cex=0.65)

plot(ts(exp(out[-(1:200000),2])), main = 'Trace Plot', ylab = TeX(r'($\kappa$)'))
acf(exp(out[-(1:200000),2]), main = 'ACF Plot')
kappadens <- density(exp(out[-(1:200000),2]))
plot(kappadens, main = 'Posterior Density Plot', xlab = TeX(r'($\kappa$)'))
abline(v=quantile(exp(out[-(1:200000),2]), c(0.025, 0.975)),col=c(3,3))
abline(v=kappadens$x[which.max(kappadens$y)],col = 2)
points(kappadens$x[which.max(kappadens$y)], kappadens$y[which.max(kappadens$y)], col = 2, lwd = 2, pch = 13, cex = 2)
abline(v=kappa, col = 4, lty = 2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds", "True value"),
       col=c(2, 3, 4), 
       lty=c(1,1,2), 
       bty = "n",
       cex=0.65)

plot(ts(exp(out[-(1:200000),3])), main = 'Trace Plot', ylab = TeX(r'($\theta$)'))
acf(exp(out[-(1:200000),3]), main = 'ACF Plot')
thetadens <- density(exp(out[-(1:200000),3]))
plot(thetadens, main = 'Posterior Density Plot', xlab = TeX(r'($\theta$)'))
abline(v=quantile(exp(out[-(1:200000),3]), c(0.025, 0.975)),col=c(3,3))
abline(v=thetadens$x[which.max(thetadens$y)],col = 2)
points(thetadens$x[which.max(thetadens$y)], thetadens$y[which.max(thetadens$y)], col = 2, lwd = 2, pch = 13, cex = 2)
abline(v=theta, col = 4, lty = 2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds", "True value"),
       col=c(2, 3, 4), 
       lty=c(1,1,2), 
       bty = "n",
       cex=0.65)

plot(ts(exp(out[-(1:200000),4])), main = 'Trace Plot', ylab = TeX(r'($\sigma$)'))
acf(exp(out[-(1:200000),4]), main = 'ACF Plot')
sigmadens <- density(exp(out[-(1:200000),4]))
plot(sigmadens, main = 'Posterior Density Plot', xlab = TeX(r'($\sigma$)'))
abline(v=quantile(exp(out[-(1:200000),4]), c(0.025, 0.975)),col=c(3,3))
abline(v=sigmadens$x[which.max(sigmadens$y)],col = 2)
points(sigmadens$x[which.max(sigmadens$y)], sigmadens$y[which.max(sigmadens$y)], col = 2, lwd = 2, pch = 13, cex = 2)
abline(v=sigma, col = 4, lty = 2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds", "True value"),
       col=c(2, 3, 4), 
       lty=c(1,1,2), 
       bty = "n",
       cex=0.65)

plot(ts(out[-(1:200000),5]), main = 'Trace Plot', ylab = TeX(r'($\rho$)'))
acf(out[-(1:200000),5], main = 'ACF Plot')
rhodens <- density(out[-(1:200000),5])
plot(rhodens, main = 'Posterior Density Plot', xlab = TeX(r'($\rho$)'))
abline(v=quantile(out[-(1:200000),5], c(0.025, 0.975)),col=c(3,3))
abline(v=rhodens$x[which.max(rhodens$y)],col = 2)
points(rhodens$x[which.max(rhodens$y)], rhodens$y[which.max(rhodens$y)], col = 2, lwd = 2, pch = 13, cex = 2)
abline(v=rho, col = 4, lty = 2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds", "True value"),
       col=c(2, 3, 4), 
       lty=c(1,1,2), 
       bty = "n",
       cex=0.65)

effectiveSize(out[-(1:200000),])

# mu summary statistics:
mean(exp(out[-(1:200000),1]))
median(exp(out[-(1:200000),1]))
mudens$x[which.max(mudens$y)]
quantile(exp(out[-(1:200000),1]), c(0.025, 0.975))

# kappa summary statistics:
mean(exp(out[-(1:200000),2]))
median(exp(out[-(1:200000),2]))
kappadens$x[which.max(kappadens$y)]
quantile(exp(out[-(1:200000),2]), c(0.025, 0.975))

# theta summary statistics:
mean(exp(out[-(1:200000),3]))
median(exp(out[-(1:200000),3]))
thetadens$x[which.max(thetadens$y)]
quantile(exp(out[-(1:200000),3]), c(0.025, 0.975))

# sigma summary statistics:
mean(exp(out[-(1:200000),4]))
median(exp(out[-(1:200000),4]))
sigmadens$x[which.max(sigmadens$y)]
quantile(exp(out[-(1:200000),4]), c(0.025, 0.975))

# rho summary statistics:
mean(out[-(1:200000),5])
median(out[-(1:200000),5])
rhodens$x[which.max(rhodens$y)]
quantile(out[-(1:200000),5], c(0.025, 0.975))