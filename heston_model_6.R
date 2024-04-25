# Heston Model Code File 6

# This file uses the modified diffusion bridge scheme on the fully observed dataset.

################################################################################
## Libraries:

library(MASS)
library(latex2exp)
library(scales)
library(coda)

################################################################################
## Inputs Required:

# heston_model_1.R
# - Yobs
# - m
# - Y
# - n

################################################################################
## Functions Required:

# black_scholes_model_1.R
# - initSVf ---------- Initial Value Function

# heston_model_1.R
# - lprior ----------- Calculates log-prior density
# - lgpdf ------------ Calculates log of the multivariate normal density
# - llikehood -------- Calculates log-likelihood 
# - lpost ------------ Calculates log-Posterior
# - MH.P ------------- Metropolis-Hastings Scheme for Parameters

# heston_model_2.R
# - EM.mu ------------ Calculates one-step mean for likelihood function
# - EM.Sig ----------- Calculates one-step variance for likelihood function

# heston_model_6.R
# - mean.fun --------- Restructured EM.mu
# - var.fun ---------- Restructured EM.Sig
# - varf3 ------------ Restructured EM.Sig
# - block.update1 ---- Performs the block update under the MDB scheme
# - gibbs.mdb.full --- Gibbs Sampler for fully observed data, modified diffusion bridge

################################################################################
## New Functions:

mean.fun <- function(s,v,Tvec,dt) c(s*(1 + Tvec[1]*dt), v + Tvec[2]*(Tvec[3]-v)*dt)
var.fun <- function(s,v,Tvec,dt) v*dt*matrix(c(s*s,Tvec[4]*Tvec[5]*s,Tvec[4]*Tvec[5]*s,Tvec[4]*Tvec[4]),nrow=2,ncol=2)
varf3 <- function(vec) var.fun(vec[1],vec[2],Tvec,1)

block.update1 <- function(curr, Tvec){
  m <- ncol(curr)-1
  deltatau <- 1/m
  tau <- seq(0,1,length = m+1)
  
  start <- curr[,1]; end <- curr[,(m+1)]
  
  store <- matrix(NA, nrow = 2, ncol = m+1)
  store[,1] <- start; store[,(m+1)] <- end
  
  can <- store
  num1 <- rep(NA, m); den1 <- rep(NA, m)
  numq <- rep(NA, m-1); denq <- rep(NA, m-1)
  
  meanf1 <- function(vec) EM.mu(vec, m, Tvec)
  varf1 <- function(vec) EM.Sig(vec, m, Tvec)
  
  for (k in 1:m){
    if (k!=m){
      beta.k <- varf3(can[,k])
      mu.k <- ( (end - can[,k])/(1 - tau[k]) )
      Phi.k <- ( (1 - tau[k+1])/(1 - tau[k]) )*beta.k
      
      beta.k1 <- varf3(curr[,k])
      mu.k1 <- ( (end - curr[,k])/(1 - tau[k]) )
      Phi.k1 <- ( (1 - tau[k+1])/(1 - tau[k]) )*beta.k1
      
      canvec <- mvrnorm(1, can[,k] + mu.k*deltatau, Phi.k*deltatau)
      
      can[,(k+1)] <- canvec
      
      numq[k] <- lgpdf(curr[,(k+1)], curr[,k] + mu.k1*deltatau, Phi.k1*deltatau)
      denq[k] <- lgpdf(can[,(k+1)], can[,k] + mu.k*deltatau, Phi.k*deltatau)
    }
    
    num1[k] <- lgpdf(can[,(k+1)], meanf1(can[,k]), varf1(can[,k]))
    den1[k] <- lgpdf(curr[,(k+1)], meanf1(curr[,k]), varf1(curr[,k]))
    
  }
  
  laprob <- sum(num1)+sum(numq)-sum(den1)-sum(denq)
  
  final <- curr
  if(log(runif(1)) < laprob){
    final <- can
  }
  return(final)
}


gibbs.mdb.full <- function(N, initP, initS, initV, Sigma, a, b, m){
  nn <- length(initS)
  n <- (nn-1)/m
  
  Pmat <- matrix(NA, nrow = N, ncol = length(Pvec))
  Smat <- matrix(NA, nrow = N, ncol = nn)
  Vmat <- matrix(NA, nrow = N, ncol = nn)
  
  Pmat[1,] <- initP
  Smat[1,] <- initS
  Vmat[1,] <- initV
  
  store <- matrix(NA, nrow = 2, ncol = nn)
  deltatau <- 1/m                
  for (k in 2:N){
    data <- store
    data[1,] <- Smat[(k-1),]
    data[2,] <- Vmat[(k-1),]
    
    # Update Pmat
    Pmat[k,] <- MH.P(2, Pmat[(k-1),], t(data), deltatau, Sigma, a, b)[2,]
    Tvec <- c(exp(Pmat[k,1]), exp(Pmat[k,2]), exp(Pmat[k,3]), exp(Pmat[k,4]), Pmat[k,5])
    
    # Update Smat and Vmat, block i
    for (i in 1:n){
      block <- data[,(1+(i-1)*m):(1+i*m)]
      SVblock <- block.update1(block, Tvec)
      Smat[k,(1+(i-1)*m):(1+i*m)] <- SVblock[1,]
      Vmat[k,(1+(i-1)*m):(1+i*m)] <- SVblock[2,]
    }
  }
  return(list(S = Smat,
              V = Vmat,
              P = Pmat))
}

################################################################################
## Output:
Sobs <- Yobs[,1]
Vobs <- Yobs[,2]

init <- initSVf(Sobs, Vobs, m)
initS <- init[,1]
initV <- init[,2]
initP <- rep(0,5)

a <- rep(-1.2,4)
b <- rep(0.5,4)

ll <- 2.40/sqrt(101)
Sigma <- diag(rep(ll/sqrt(5), 5))
output.mdb.full1 <- gibbs.mdb.full(2000, initP, initS, initV, Sigma, a, b, m)

postvar <- var(output.mdb.full1$P)
lambda <- 0.35
Sigmatune <- lambda*(2.40^2)*postvar/5
system.time(output.mdb.full2 <- gibbs.mdb.full(20000, initP, initS, initV, Sigmatune, a, b, m))
length(unique(output.mdb.full2$P[,1]))/20000
# Target AP: 28.39%
# CPU: 900s

par(mfrow=c(1,1))
plot(ts(Y[,1],start=0,frequency=n),ylim=c(min(output.mdb.full2$S),max(output.mdb.full2$S)), ylab = 'S', main = 'MDB Augmented Time Series, S')
for (i in 0:1600){
  lines(ts(output.mdb.full2$S[(4000+10*i),],start=0,frequency=m), col = i)
}
lines(ts(Y[,1],start=0,frequency=n), lwd=2)

plot(ts(Y[,2],start=0,frequency=n),ylim=c(min(output.mdb.full2$V),max(output.mdb.full2$V)), ylab = 'V', main = 'MDB Augmented Time Series, V')
for (i in 0:1600){
  lines(ts(output.mdb.full2$V[(400+10*i),],start=0,frequency=m), col = i)
}
lines(ts(Y[,2],start=0,frequency=n), lwd=2)

out <- output.mdb.full2
par(mfrow = c(5,3))
plot(ts(exp(out$P[-(1:4000),1])), main = 'Trace Plot', ylab = TeX(r'($\mu$)'))
acf(exp(out$P[-(1:4000),1]), main = 'ACF Plot')
mudens <- density(exp(out$P[-(1:4000),1]))
plot(mudens, main = 'Posterior Density Plot', xlab = TeX(r'($\mu$)'))
abline(v=quantile(exp(out$P[-(1:4000),1]), c(0.025, 0.975)),col=c(3,3))
abline(v=mudens$x[which.max(mudens$y)],col = 2)
points(mudens$x[which.max(mudens$y)], mudens$y[which.max(mudens$y)], col = 2, lwd = 2, pch = 13, cex = 2)
abline(v=mu, col = 4, lty = 2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds", "True value"),
       col=c(2, 3, 4), 
       lty=c(1,1,2), 
       bty = "n",
       cex=0.65)

plot(ts(exp(out$P[-(1:4000),2])), main = 'Trace Plot', ylab = TeX(r'($\kappa$)'))
acf(exp(out$P[-(1:4000),2]), main = 'ACF Plot')
kappadens <- density(exp(out$P[-(1:4000),2]))
plot(kappadens, main = 'Posterior Density Plot', xlab = TeX(r'($\kappa$)'))
abline(v=quantile(exp(out$P[-(1:4000),2]), c(0.025, 0.975)),col=c(3,3))
abline(v=kappadens$x[which.max(kappadens$y)],col = 2)
points(kappadens$x[which.max(kappadens$y)], kappadens$y[which.max(kappadens$y)], col = 2, lwd = 2, pch = 13, cex = 2)
abline(v=kappa, col = 4, lty = 2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds", "True value"),
       col=c(2, 3, 4), 
       lty=c(1,1,2), 
       bty = "n",
       cex=0.65)

plot(ts(exp(out$P[-(1:4000),3])), main = 'Trace Plot', ylab = TeX(r'($\theta$)'))
acf(exp(out$P[-(1:4000),3]), main = 'ACF Plot')
thetadens <- density(exp(out$P[-(1:4000),3]))
plot(thetadens, main = 'Posterior Density Plot', xlab = TeX(r'($\theta$)'))
abline(v=quantile(exp(out$P[-(1:4000),3]), c(0.025, 0.975)),col=c(3,3))
abline(v=thetadens$x[which.max(thetadens$y)],col = 2)
points(thetadens$x[which.max(thetadens$y)], thetadens$y[which.max(thetadens$y)], col = 2, lwd = 2, pch = 13, cex = 2)
abline(v=theta, col = 4, lty = 2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds", "True value"),
       col=c(2, 3, 4), 
       lty=c(1,1,2), 
       bty = "n",
       cex=0.65)

plot(ts(exp(out$P[-(1:4000),4])), main = 'Trace Plot', ylab = TeX(r'($\sigma$)'))
acf(exp(out$P[-(1:4000),4]), main = 'ACF Plot')
sigmadens <- density(exp(out$P[-(1:4000),4]))
plot(sigmadens, main = 'Posterior Density Plot', xlab = TeX(r'($\sigma$)'))
abline(v=quantile(exp(out$P[-(1:4000),4]), c(0.025, 0.975)),col=c(3,3))
abline(v=sigmadens$x[which.max(sigmadens$y)],col = 2)
points(sigmadens$x[which.max(sigmadens$y)], sigmadens$y[which.max(sigmadens$y)], col = 2, lwd = 2, pch = 13, cex = 2)
abline(v=sigma, col = 4, lty = 2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds", "True value"),
       col=c(2, 3, 4), 
       lty=c(1,1,2), 
       bty = "n",
       cex=0.65)

plot(ts(out$P[-(1:4000),5]), main = 'Trace Plot', ylab = TeX(r'($\rho$)'))
acf(out$P[-(1:4000),5], main = 'ACF Plot')
rhodens <- density(out$P[-(1:4000),5])
plot(rhodens, main = 'Posterior Density Plot', xlab = TeX(r'($\rho$)'))
abline(v=quantile(out$P[-(1:4000),5], c(0.025, 0.975)),col=c(3,3))
abline(v=rhodens$x[which.max(rhodens$y)],col = 2)
points(rhodens$x[which.max(rhodens$y)], rhodens$y[which.max(rhodens$y)], col = 2, lwd = 2, pch = 13, cex = 2)
abline(v=rho, col = 4, lty = 2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds", "True value"),
       col=c(2, 3, 4), 
       lty=c(1,1,2), 
       bty = "n",
       cex=0.65)

effectiveSize(out$P[-(1:4000),])

# mu summary statistics:
mean(exp(out$P[-(1:4000),1]))
median(exp(out$P[-(1:4000),1]))
mudens$x[which.max(mudens$y)]
quantile(exp(out$P[-(1:4000),1]), c(0.025, 0.975))

# kappa summary statistics:
mean(exp(out$P[-(1:4000),2]))
median(exp(out$P[-(1:4000),2]))
kappadens$x[which.max(kappadens$y)]
quantile(exp(out$P[-(1:4000),2]), c(0.025, 0.975))

# theta summary statistics:
mean(exp(out$P[-(1:4000),3]))
median(exp(out$P[-(1:4000),3]))
thetadens$x[which.max(thetadens$y)]
quantile(exp(out$P[-(1:4000),3]), c(0.025, 0.975))

# sigma summary statistics:
mean(exp(out$P[-(1:4000),4]))
median(exp(out$P[-(1:4000),4]))
sigmadens$x[which.max(sigmadens$y)]
quantile(exp(out$P[-(1:4000),4]), c(0.025, 0.975))

# rho summary statistics:
mean(out$P[-(1:4000),5])
median(out$P[-(1:4000),5])
rhodens$x[which.max(rhodens$y)]
quantile(out$P[-(1:4000),5], c(0.025, 0.975))