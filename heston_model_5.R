# Heston Model Code File 5

# This file uses a single-site update scheme on the unobserved observed dataset,
# where volatility is completely unobserved. 

################################################################################
## Libraries:

library(MASS)
library(latex2exp)
library(scales)
library(coda)

################################################################################
## Inputs Required:

# heston_model_1.R
# - Sobs
# - m 
# - Y
# - n
# - kappadens
# - sigmadens

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

# heston_model_5.R
# - MH.imp.none ------ Metropolis-Hastings Scheme for unobserved volatility data
# - gibbs.ssu.none --- Gibbs Sampler for unobserved volatility data, single-site updates

################################################################################
## New Functions:

MH.imp.none <- function(N, initS, initV, Pvec, m, sigw=0.75, alpha=3, beta=2){
  nn <- length(initS)
  
  Smat <- matrix(NA, nrow = N, ncol = nn)
  Vmat <- matrix(NA, nrow = N, ncol = nn)
  
  Smat[1,] <- initS
  Vmat[1,] <- initV
  V0 <- rgamma(1, shape = alpha, rate = beta)
  Vmat[1,1] <- V0
  
  deltatau <- 1/m
  Tvec <- c(exp(Pvec[1]), exp(Pvec[2]), exp(Pvec[3]), exp(Pvec[4]), Pvec[5])
  
  store <- matrix(NA, nrow = 2, ncol = nn)
  
  meanf1 <- function(vec) EM.mu(vec, m, Tvec)
  varf1 <- function(vec) EM.Sig(vec, m, Tvec)
  
  for (k in 2:N){
    curr <- store
    curr[1,] <- Smat[(k-1),]
    curr[2,] <- Vmat[(k-1),]
    
    can <- store
    for (i in 1:nn){
      currvec <- curr[,i]
      if(i==1){
        canvec <- c(currvec[1],exp(rnorm(1,log(currvec[1]),sigw)))
        if(canvec[2]<0){
          canvec[2] <- 0.0001
        }
        num1 <- dgamma(canvec[2], shape = alpha, rate = beta, log = TRUE)
        den1 <- dgamma(currvec[2], shape = alpha, rate = beta, log = TRUE)
        num2 <- lgpdf(curr[,(i+1)],meanf1(canvec),varf1(canvec))
        den2 <- lgpdf(curr[,(i+1)],meanf1(currvec),varf1(currvec))
        numq <- canvec[2]
        denq <- currvec[2]
        
      }else if(i==nn){
        canvec <- c(currvec[1],rnorm(1, can[2,(i-1)],Tvec[4]*sqrt(deltatau*can[2,(i-1)])))
        if(canvec[2]<0){
          canvec[2] <- 0.0001
        }
        num1 <- lgpdf(canvec, meanf1(can[,(i-1)]),varf1(can[,(i-1)]))
        den1 <- lgpdf(currvec, meanf1(can[,(i-1)]),varf1(can[,(i-1)]))
        num2 <- 0
        den2 <- 0
        numq <- dnorm(currvec[2],can[2,(i-1)],Tvec[4]*sqrt(deltatau*can[2,(i-1)]),log=TRUE)
        denq <- dnorm(canvec[2],can[2,(i-1)],Tvec[4]*sqrt(deltatau*can[2,(i-1)]),log=TRUE)
        
      }else if(i%%m==1 | m==1){
        canvec <- c(currvec[1],rnorm(1,0.5*(can[2,(i-1)]+curr[2,(i+1)]),Tvec[4]*sqrt(deltatau*can[2,(i-1)])))
        if(canvec[2]<0){
          canvec[2] <- 0.0001
        }
        num1 <- lgpdf(canvec,  meanf1(can[,(i-1)]),varf1(can[,(i-1)]))
        den1 <- lgpdf(currvec, meanf1(can[,(i-1)]),varf1(can[,(i-1)]))
        num2 <- lgpdf(curr[,(i+1)],meanf1(canvec), varf1(canvec) )
        den2 <- lgpdf(curr[,(i+1)],meanf1(currvec),varf1(currvec))
        numq <- dnorm(currvec[2],0.5*(can[2,(i-1)]+curr[2,(i+1)]),Tvec[4]*sqrt(deltatau*can[2,(i-1)]),log=TRUE)
        denq <- dnorm(canvec[2],0.5*(can[2,(i-1)]+curr[2,(i+1)]),Tvec[4]*sqrt(deltatau*can[2,(i-1)]),log=TRUE)
        
      }else{
        canvec <- mvrnorm(1, 0.5*(can[,(i-1)]+curr[,(i+1)]), varf1(can[,(i-1)]))
        if(canvec[2]<0){
          canvec[2] <- 0.0001
        }
        num1 <- lgpdf(canvec,  meanf1(can[,(i-1)]),varf1(can[,(i-1)]))
        den1 <- lgpdf(currvec, meanf1(can[,(i-1)]),varf1(can[,(i-1)]))
        num2 <- lgpdf(curr[,(i+1)],meanf1(canvec), varf1(canvec) )
        den2 <- lgpdf(curr[,(i+1)],meanf1(currvec),varf1(currvec))
        numq <- lgpdf(currvec,0.5*(can[,(i-1)]+curr[,(i+1)]),varf1(can[,(i-1)]))
        denq <- lgpdf(canvec, 0.5*(can[,(i-1)]+curr[,(i+1)]),varf1(can[,(i-1)]))
        
      }
      laprob <- num1+num2+numq-den1-den2-denq
      if(log(runif(1))<laprob){
        currvec <- canvec
      }
      can[,i] <- currvec
    }
    Smat[k,] <- can[1,]
    Vmat[k,] <- can[2,]
  }
  return(list(S = Smat,
              V = Vmat))
}

gibbs.ssu.none <- function(N, initP, initS, Sigma, a, b, m, sigw, alpha, beta){
  nn <- length(initS)
  
  Pmat <- matrix(NA, nrow = N, ncol = length(initP))
  Smat <- matrix(NA, nrow = N, ncol = nn)
  Vmat <- matrix(NA, nrow = N, ncol = nn)
  
  Pmat[1,] <- initP
  Smat[1,] <- initS
  V0 <- rgamma(1, shape = alpha, rate = beta)
  Vmat[1,] <- rep(V0,nn)
  
  store <- matrix(NA, nrow = 2, ncol = nn)
  deltatau <- 1/m                
  for (k in 2:N){
    data <- store
    data[1,] <- Smat[(k-1),]
    data[2,] <- Vmat[(k-1),]
    
    # Update Pmat
    Pmat[k,] <- MH.P(2, Pmat[(k-1),], t(data), deltatau, Sigma, a, b)[2,]
    
    # Update Smat and Vmat
    SVmat <- MH.imp.none(2, data[1,], data[2,], Pmat[k,], m, sigw, alpha, beta)
    Smat[k,] <- SVmat$S[2,]
    Vmat[k,] <- SVmat$V[2,]
  }
  return(list(S = Smat,
              V = Vmat,
              P = Pmat))
}

################################################################################
## Output:

# MCMC

initS <- initSVf(Sobs, NA, m)
initP <- rep(0,5)

a <- rep(-1.2,4)
b <- rep(0.5,4)
a[2] <- log(kappadens$x[which.max(kappadens$y)])
b[2] <- 0.05
a[4] <- log(sigmadens$x[which.max(sigmadens$y)])
b[4] <- 0.05

# Investigating V_0 prior
plot(seq(0,10, length=100), dgamma(seq(0,10,length=100), shape = 1, rate = 2), 'l')

ll <- 2.4/sqrt(101)
Sigma <- diag(rep(ll/sqrt(5), 5))
output.ssu.none1 <- gibbs.ssu.none(2000, initP, initS, Sigma, a, b, m, sigw=0.5, alpha=1, beta=2)

postvar <- var(output.ssu.none1$P)
lambda <- 0.4
Sigmatune <- lambda*(2.40^2)*postvar/5
system.time(output.ssu.none2 <- gibbs.ssu.none(20000, initP, initS, Sigmatune, a, b, m, sigw=0.5, alpha=1, beta=2))
length(unique(output.ssu.none2$P[,1]))/20000
# Target AP: 28.39%
# CPU approx: 1200s

par(mfrow=c(1,1))
plot(ts(Y[,1],start=0,frequency=n),ylim=c(min(output.ssu.none2$S),max(output.ssu.none2$S)), ylab = 'S', main = 'SSU Augmented Time Series, S')
for (i in 0:1600){
  lines(ts(output.ssu.none2$S[(4000+10*i),],start=0,frequency=m), col = i)
}
lines(ts(Y[,1],start=0,frequency=n), lwd=2)

plot(ts(Y[,2],start=0,frequency=n),ylim=c(min(output.ssu.none2$V),max(output.ssu.none2$V)), ylab = 'V', main = 'SSU Augmented Time Series, V')
for (i in 0:1600){
  lines(ts(output.ssu.none2$V[(400+10*i),],start=0,frequency=m), col = i)
}
lines(ts(Y[,2],start=0,frequency=n), lwd=2)

out <- output.ssu.none2
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