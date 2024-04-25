# Heston Model Code File 7

# This file uses the modified diffusion bridge scheme on the partially observed dataset,
# where volatility is observed only at the start and end of the interval. 

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
# - nn
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

# heston_model_6.R
# - block.update1 ---- Performs the block update under the MDB scheme, F = I

# heston_model_7.R
# - pMDB.mu ---------- Calculates partial MDB mu
# - pMDB.psi --------- Calculates partial MDB psi
# - CFS.mu ----------- Calculates Conditional Forward Simulation mu
# - CFS.Sig ---------- Calculates Conditional Forward Simulation Sigma
# - block.update2 ---- Performs the block update under the partial MDB Scheme, F = (1,0)
# - gibbs.mdb.part1 -- Gibbs Sampler for partially observed data (start and end), modified diffusion bridge

################################################################################
## New Functions:

pMDB.mu <- function(vec, j, m, sT, Tvec){
  mu <- Tvec[1]; kappa <- Tvec[2]; theta <- Tvec[3]; sigma <- Tvec[4]; rho <- Tvec[5]
  s <- vec[1]; v <- vec[2]
  
  deltatau <- 1/m
  delta <- deltatau*(m-j)
  
  vec1 <- (sT-s)/delta
  vec2 <- kappa*(theta-v) - mu*rho*sigma + rho*sigma*(sT-s)/(delta*s)
  return(c(vec1, vec2))
}

pMDB.psi <- function(vec, j, m, sT, Tvec){
  mu <- Tvec[1]; kappa <- Tvec[2]; theta <- Tvec[3]; sigma <- Tvec[4]; rho <- Tvec[5]
  s <- vec[1]; v <- vec[2]
  
  ratio <- (m - j - 1)/(m - j)
  Sig11 <- s*s*v*ratio
  Sig12 <- rho*sigma*s*v*ratio
  Sig22 <- sigma*sigma*v*(1 - rho*rho/(m-j))
  mat1 <- matrix(c(Sig11, Sig12,
                   Sig12, Sig22), nrow = 2, ncol = 2)
  return(mat1)
}

CFS.mu <- function(vec, m, sT, Tvec){
  mu <- Tvec[1]; kappa <- Tvec[2]; theta <- Tvec[3]; sigma <- Tvec[4]; rho <- Tvec[5]
  s <- vec[1]; v <- vec[2]
  
  deltatau <- 1/m
  val1 <- v + kappa*(theta-v)*deltatau + rho*sigma*(sT-s)/s - mu*rho*sigma*deltatau
  return(val1)
}

CFS.Sig <- function(vec, m, sT, Tvec){
  mu <- Tvec[1]; kappa <- Tvec[2]; theta <- Tvec[3]; sigma <- Tvec[4]; rho <- Tvec[5]
  s <- vec[1]; v <- vec[2]
  
  deltatau <- 1/m
  val1 <- 1-rho^2
  val2 <- sigma*sigma*v*val1*deltatau
  return(val2)
}

block.update2 <- function(block, Tvec){
  m <- ncol(block)-1
  deltatau <- 1/m
  
  current <- block
  candidate <- matrix(NA, nrow = 2, ncol = m+1)
  candidate[,1] <- block[,1]
  sT <- block[1,(m+1)]
  candidate[1,(m+1)] <- sT
  
  laprob <- 0
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

gibbs.mdb.part1 <- function(N, initP, initS, initV, Sigma, a, b, m){
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
    
    # Update Smat and Vmat
    
    # Block i, (end of block unknown)
    for (i in 1:(n-1)){
      block <- data[,(1+(i-1)*m):(1+i*m)]
      SVblock <- block.update2(block, Tvec)
      Smat[k,(1+(i-1)*m):(1+i*m)] <- SVblock[1,]
      Vmat[k,(1+(i-1)*m):(1+i*m)] <- SVblock[2,]
      data[2,(1+i*m)] <- SVblock[2,(m+1)]
    }
    
    # Block n, (end of block known)
    block <- data[,(1+(n-1)*m):(1+n*m)]
    SVblock <- block.update1(block, Tvec)
    Smat[k,(1+(n-1)*m):(1+n*m)] <- SVblock[1,]
    Vmat[k,(1+(n-1)*m):(1+n*m)] <- SVblock[2,]
  }
  
  return(list(S = Smat,
              V = Vmat,
              P = Pmat))
}

################################################################################
## Output:

Vobs <- c(Y[1,2], Y[nrow(Y),2])
initP <- c(rep(-0.7, 4),0)
initS <- initSVf(Sobs, NA, m)
initV <- seq(Vobs[1], Vobs[2], length.out=nn)

a <- rep(-1.2,4)
b <- rep(0.5,4)
a[2] <- log(kappadens$x[which.max(kappadens$y)])
b[2] <- 0.02
a[4] <- log(sigmadens$x[which.max(sigmadens$y)])
b[4] <- 0.02

ll <- 2.40/sqrt(101)
Sigma <- diag(rep(ll/sqrt(5), 5))
output.mdb.part1.1 <- gibbs.mdb.part1(2000, initP, initS, initV, Sigma, a, b, m)

postvar <- var(output.mdb.part1.1$P)
lambda <- 0.21
Sigmatune <- lambda*(2.40^2)*postvar/5
system.time(output.mdb.part1.2 <- gibbs.mdb.part1(20000, initP, initS, initV, Sigmatune, a, b, m))
length(unique(output.mdb.part1.2$P[,1]))/20000
# Target AP: 28.39%
# CPU approx: 900s

par(mfrow=c(1,1))
plot(ts(Y[,1],start=0,frequency=n),ylim=c(min(output.mdb.part1.2$S),max(output.mdb.part1.2$S)), ylab = 'S', main = 'MDB Augmented Time Series, S')
for (i in 0:1600){
  lines(ts(output.mdb.part1.2$S[(4000+10*i),],start=0,frequency=m), col = i)
}
lines(ts(Y[,1],start=0,frequency=n), lwd=2)

plot(ts(Y[,2],start=0,frequency=n),ylim=c(min(output.mdb.part1.2$V),max(output.mdb.part1.2$V)), ylab = 'V', main = 'MDB Augmented Time Series, V')
for (i in 0:1600){
  lines(ts(output.mdb.part1.2$V[(400+10*i),],start=0,frequency=m), col = i)
}
lines(ts(Y[,2],start=0,frequency=n), lwd=2)

out <- output.mdb.part1.2
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