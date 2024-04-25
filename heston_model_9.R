# Heston Model Code File 9

# This file uses the modified diffusion bridge scheme on the unobserved observed dataset,
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

# heston_model_7.R
# - pMDB.mu ---------- Calculates partial MDB mu
# - pMDB.psi --------- Calculates partial MDB psi
# - CFS.mu ----------- Calculates Conditional Forward Simulation mu
# - CFS.Sig ---------- Calculates Conditional Forward Simulation Sigma
# - block.update2 ---- Performs the block update under the partial MDB Scheme, F = (1,0)

# heston_model_9.R
# - block.update3 --- Performs the block update where no volatility is observed (block 1 only)
# - gibbs.mdb.none -- Gibbs Sampler for unobserved volatility data, modified diffusion bridge

################################################################################
## New Functions:

block.update3 <- function(block, Tvec, sigw, alpha, beta){
  m <- ncol(block)-1
  deltatau <- 1/m
  
  current <- block
  candidate <- matrix(NA, nrow = 2, ncol = m+1)
  sT <- block[1,(m+1)]
  candidate[1,(m+1)] <- sT
  
  # First: V0
  cur <- current[2,1]
  can <- exp(rnorm(1, log(current), sigw))
  numq <- dgamma(cur, shape = alpha, rate = beta, log=TRUE)*can
  denq <- dgamma(can, shape = alpha, rate = beta, log=TRUE)*cur
  can <- c(block[1,1], can)
  candidate[,1] <- can
  num1 <- lgpdf(current[,2], EM.mu(can, m, Tvec), EM.Sig(can, m, Tvec))
  den1 <- lgpdf(current[,2], EM.mu(current[,1], m, Tvec), EM.Sig(current[,1], m, Tvec))
  
  laprob <- num1 - den1 + numq - denq
  
  # Second: Vt, t>0
  for (j in 1:m){
    precur <- current[,j]
    cur <- current[,(j+1)]
    precan <- candidate[,j]
    s <- precan[1]; v <- precan[2]
    if (j!=m){
      # Update S,V using Partial MDB
      M <- precan + pMDB.mu(precan, j-1, m, sT, Tvec)*deltatau
      S <- pMDB.psi(precan, j-1, m, sT, Tvec)*deltatau
      can <- mvrnorm(1, M, S)
      if(can[2]<0.0001){
        can[2] <- 0.0001
      }
      numq <- lgpdf(cur, precur + pMDB.mu(precur, j-1, m, sT, Tvec)*deltatau, pMDB.psi(precur, j-1, m, sT, Tvec)*deltatau)
      denq <- lgpdf(can, M, S)
    }else{
      # Update using Conditional Forward Step
      M <- CFS.mu(precan, m, sT, Tvec)
      S <- CFS.Sig(precan, m, sT, Tvec)
      can <- rnorm(1, M, sqrt(S))
      can <- c(sT, can)
      if(can[2]<0.0001){
        can[2] <- 0.0001
      }
      numq <- dnorm(cur[2], CFS.mu(precur, m, sT, Tvec), sqrt(CFS.Sig(precur, m, sT, Tvec)), log = TRUE)
      denq <- dnorm(can[2], M, sqrt(S), log = TRUE)
    }
    num1 <- lgpdf(can, EM.mu(precan, m, Tvec), EM.Sig(precan, m, Tvec))
    den1 <- lgpdf(cur, EM.mu(precur, m, Tvec), EM.Sig(precur, m, Tvec))
    laprob <- laprob + num1 - den1 + numq - denq
    candidate[,(j+1)] <- can
  }
  if(log(runif(1))<laprob){
    current <- candidate
  }
  return(current)
}

gibbs.mdb.none <- function(N, initP, initS, Sigma, a, b, m, sigw, alpha, beta){
  nn <- length(initS)
  n <- (nn-1)/m
  
  Pmat <- matrix(NA, nrow = N, ncol = length(initP))
  Smat <- matrix(NA, nrow = N, ncol = nn)
  Vmat <- matrix(NA, nrow = N, ncol = nn)
  
  Pmat[1,] <- initP
  Smat[1,] <- initS
  V0 <- rgamma(1, shape = alpha, rate = beta)
  Vmat[1,] <- seq(V0, exp(initP[3]), length.out = nn)
  
  store <- matrix(NA, nrow = 2, ncol = nn)
  deltatau <- 1/m                
  for (k in 2:N){
    data <- store
    data[1,] <- Smat[(k-1),]
    data[2,] <- Vmat[(k-1),]
    
    # Update Pmat
    Pmat[k,] <- MH.P(2, Pmat[(k-1),], t(data), deltatau, Sigma, a, b)[2,]
    Tvec <- c(exp(Pmat[k,1]), exp(Pmat[k,2]), exp(Pmat[k,3]), exp(Pmat[k,4]), Pmat[k,5])

    
    # Update Smat and Vmat
    
    # Block 1, (start unknown, end unknown)
    block <- data[,1:(m+1)]
    SVblock <- block.update3(block, Tvec, sigw, alpha, beta)
    Smat[k,1:(1+m)] <- SVblock[1,]
    Vmat[k,1:(1+m)] <- SVblock[2,]
    data[2,(1+m)] <- SVblock[2,(m+1)]
    
    # Block i, (start known, end unknown)
    for (i in 2:n){
      block <- data[,(1+(i-1)*m):(1+i*m)]
      SVblock <- block.update2(block, Tvec)
      Smat[k,(1+(i-1)*m):(1+i*m)] <- SVblock[1,]
      Vmat[k,(1+(i-1)*m):(1+i*m)] <- SVblock[2,]
      data[2,(1+i*m)] <- SVblock[2,(m+1)]
    }
  }
  return(list(S = Smat,
              V = Vmat,
              P = Pmat))
}

################################################################################
## Output:

initP <- c(rep(-2, 4),0)
initS <- initSVf(Sobs, NA, m)

a <- rep(-1.2,4)
b <- rep(0.5,4)
a[2] <- log(kappadens$x[which.max(kappadens$y)])
b[2] <- 0.05
a[4] <- log(sigmadens$x[which.max(sigmadens$y)])
b[4] <- 0.05

ll <- 2.40/sqrt(101)
Sigma <- diag(rep(ll/sqrt(5), 5))
output.mdb.none1 <- gibbs.mdb.none(2000, initP, initS, Sigma, a, b, m, 0.5, 1, 2)

postvar <- var(output.mdb.none1$P)
lambda <- 0.20
Sigmatune <- lambda*(2.40^2)*postvar/5
system.time(output.mdb.none2 <- gibbs.mdb.none(20000, initP, initS, Sigmatune, a, b, m, 0.5, 1, 2))
length(unique(output.mdb.none2$P[,1]))/20000
# Target AP: 28.39%
# CPU approx: 900s

par(mfrow=c(1,1))
plot(ts(Y[,1],start=0,frequency=n),ylim=c(min(output.mdb.none2$S),max(output.mdb.none2$S)), ylab = 'S', main = 'MDB Augmented Time Series, S')
for (i in 0:1600){
  lines(ts(output.mdb.none3$S[(4000+10*i),],start=0,frequency=m), col = i)
}
lines(ts(Y[,1],start=0,frequency=n), lwd=2)

plot(ts(Y[,2],start=0,frequency=n),ylim=c(min(output.mdb.none2$V),max(output.mdb.none2$V)), ylab = 'V', main = 'MDB Augmented Time Series, V')
for (i in 0:1600){
  lines(ts(output.mdb.none3$V[(400+10*i),],start=0,frequency=m), col = i)
}
lines(ts(Y[,2],start=0,frequency=n), lwd=2)

out <- output.mdb.none3
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