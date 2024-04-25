# Real Data Application Code File 1 

# This file performs inference on a real dataset

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

# No new functions

################################################################################
## Output:

# Mercedes-Benz
# This file, including the URL and exact dates recorded are found in the report
MBG.DE <- read.csv("MBG.DE.csv", header = TRUE)
plot(ts(MBG.DE$Open), ylim = c(55, 80), main = "Mercedes-Benz Group AG", ylab = "", xlab = 'Time')
S <- MBG.DE$Open
plot(seq(0, 255, by=1), S)

################################################################################
# Under BSM
# Analytic
b <- c(5,5)
initP <- c(0,0)
ll <- 2.42/sqrt(length(S))
Sigma <- diag(rep(ll/sqrt(2), 2))
MBG.out1.1 <- MHA(100000, S, deltat=1, initP, Sigma, b)
postvar <- var(MBG.out1.1)
lambda <- 0.37
Sigmatune <- lambda*(2.42^2)*postvar/2
# Target AP: 35%
# CPU approx: 160s

# Run 1 
system.time(MBG.out1.2 <- MHA(1000000, S, deltat=1, initP, Sigmatune, b))
length(unique(MBG.out1.2[,1]))/1000000
# Run 2
system.time(MBG.out1.3 <- MHA(1000000, S, deltat=1, initP, Sigmatune, b))
length(unique(MBG.out1.3[,1]))/1000000
# Run 3
system.time(MBG.out1.4 <- MHA(1000000, S, deltat=1, initP, Sigmatune, b))
length(unique(MBG.out1.4[,1]))/1000000

################################################################################
# Results 
# Try MBG.out1.2, MBG.out1.3, MBG.out1.4
out <- MBG.out1.2

par(mfrow = c(2,3))
plot(ts(exp(out[-(1:200000),1])), main = 'Trace Plot', ylab = TeX(r'($\mu$)'))
acf(exp(out[-(1:200000),1]), main = 'ACF Plot')
mudens <- density(exp(out[-(1:200000),1]))
plot(mudens, main = 'Posterior Density Plot', xlab = TeX(r'($\mu$)'), xlim = c(0,0.0005))
abline(v=quantile(exp(out[-(1:200000),1]), c(0.025, 0.975)),col=c(3,3))
abline(v=mudens$x[which.max(mudens$y)],col = 2)
points(mudens$x[which.max(mudens$y)], mudens$y[which.max(mudens$y)], col = 2, lwd = 2, pch = 13, cex = 2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds"),
       col=c(2, 3), 
       lty=c(1,1), 
       bty = "n",
       cex=0.65)
plot(ts(exp(out[-(1:200000),2])), main = 'Trace Plot', ylab = TeX(r'($\sigma$)'))
acf(exp(out[-(1:200000),2]), main = 'ACF Plot')
sigmadens <- density(exp(out[-(1:200000),2]))
plot(sigmadens, main = 'Posterior Density Plot', xlab = TeX(r'($\sigma$)'))
abline(v=quantile(exp(out[-(1:200000),2]), c(0.025, 0.975)),col=c(3,3))
abline(v=sigmadens$x[which.max(sigmadens$y)],col = 2)
points(sigmadens$x[which.max(sigmadens$y)], sigmadens$y[which.max(sigmadens$y)], col = 2, lwd = 2, pch = 13, cex = 2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds"),
       col=c(2, 3), 
       lty=c(1,1), 
       bty = "n",
       cex=0.65)

effectiveSize(exp(out[-(1:200000),]))

# mu summary statistics:
mean(exp(out[-(1:200000),1]))
median(exp(out[-(1:200000),1]))
mudens$x[which.max(mudens$y)]
quantile(exp(out[-(1:200000),1]), c(0.025, 0.975))

# sigma summary statistics:
mean(exp(out[-(1:200000),2]))
median(exp(out[-(1:200000),2]))
sigmadens$x[which.max(sigmadens$y)]
quantile(exp(out[-(1:200000),2]), c(0.025, 0.975))

set.seed(1)
# MAP estimates from the three runs
MBG.sigma <- (0.01317105 + 0.01321263 + 0.01310707)/3

system.time(mc_MBG <- mc_BS(N = 10000,
                            t = 0,
                            S0 = S[256],
                            r = 0.05,
                            sigma = MBG.sigma,
                            K = 65,
                            T = 30,
                            euro_call))

mc_MBG$price
mc_MBG$sd

d1 <- (log(S[256]/65) + (0.05 + 0.5*MBG.sigma*MBG.sigma)*(30))/(MBG.sigma*sqrt(30))
d2 <- (log(S[256]/65) + (0.05 - 0.5*MBG.sigma*MBG.sigma)*(30))/(MBG.sigma*sqrt(30))
act_MBG <- S[256]*pnorm(d1) - 65*exp(-0.05*30)*pnorm(d2)
act_MBG

################################################################################
# Under Heston
a <- c(-11, -1.2, -4.3, -1.2)
b <- c(2, 0.5, 1, 0.5)
initP <- c(rnorm(1,a[1],b[1]),rnorm(1,a[2],b[2]),rnorm(1,a[3],b[3]),rnorm(1,a[4],b[4]),0)
alpha <- 2
beta <- 180

# No volatility is observed, so cannot use MH scheme as is
# Of the schemes used, only the unobserved volatility schemes can be used
# Start with SSU, then use with priors to give MDB (due to convergence problem)

par(mfrow = c(1,1))
plot(exp(seq(a[1]-2,a[1]+2,length=1000)), dnorm(seq(a[1]-2,a[1]+2,length=1000),a[1],b[1]), 'l', main = TeX(r'(Prior for $\mu$)'), ylab = '', xlab = TeX(r'($\mu$)'))
plot(exp(seq(a[2]-2,a[2]+2,length=1000)), dnorm(seq(a[2]-2,a[2]+2,length=1000),a[2],b[2]), 'l', main = TeX(r'(Prior for $\kappa$)'), ylab = '', xlab = TeX(r'($\kappa$)'))
plot(exp(seq(a[3]-2,a[3]+2,length=1000)), dnorm(seq(a[3]-2,a[3]+2,length=1000),a[3],b[3]), 'l', main = TeX(r'(Prior for $\theta$)'), ylab = '', xlab = TeX(r'($\theta$)'))
plot(exp(seq(a[4]-2,a[4]+2,length=1000)), dnorm(seq(a[4]-2,a[4]+2,length=1000),a[4],b[4]), 'l', main = TeX(r'(Prior for $\sigma$)'), ylab = '', xlab = TeX(r'($\sigma$)'))
plot(seq(0,0.1,length.out=100), dgamma(seq(0,0.1,length.out=100), shape = alpha, rate = beta), 'l', main = TeX(r'(Prior for $V_0$)'), ylab = '', xlab = TeX(r'($V_0$)'))

# SSU
m <- 1
initS <- initSVf(S, NA, m)
ll <- 2.4/sqrt(length(initS))
Sigma <- diag(rep(ll/sqrt(5), 5))
system.time(MBG.out2.1 <- gibbs.ssu.none(1000, initP, initS, Sigma, a, b, m, sigw=0.5, alpha=2, beta=180))
postvar <- var(MBG.out6.1$P)
lambda <- 0.005
Sigmatune <- lambda*(2.40^2)*postvar/5
# Target AP: 28.39%
# CPU approx: 11000s

# Run 1
system.time(MBG.out2.2 <- gibbs.ssu.none(100000, initP, initS, Sigmatune, a, b, m, sigw=0.5, alpha=2, beta=180))
length(unique(MBG.out2.2$P[,1]))/100000
# Run 2
system.time(MBG.out2.3 <- gibbs.ssu.none(100000, initP, initS, Sigmatune, a, b, m, sigw=0.5, alpha=2, beta=180))
length(unique(MBG.out2.3$P[,1]))/100000
# Run 3
system.time(MBG.out2.4 <- gibbs.ssu.none(100000, initP, initS, Sigmatune, a, b, m, sigw=0.5, alpha=2, beta=180))
length(unique(MBG.out2.4$P[,1]))/100000
# Run 4
system.time(MBG.out2.5 <- gibbs.ssu.none(100000, initP, initS, Sigmatune, a, b, m, sigw=0.5, alpha=2, beta=180))
length(unique(MBG.out2.5$P[,1]))/100000

# Try MBG.out2.2, MBG.out2.3, MBG.out2.4, MBG.out2.5
out <- MBG.out2.2
length(unique(out$P[,1]))/100000

par(mfrow=c(1,1))
plot(ts(S, start = 0), ylim = c(55, 80), main = 'Proposed Stock Process', ylab = TeX(r'($S_{t}$)'))
for (i in 0:800){
  lines(ts(out$S[(20000+100*i),], start = 0, frequency = m), col = i)
}
lines(ts(S, start = 0), lwd = 2)

plot(ts(out$V[200,], start = 0, frequency = m), ylim = c(0, 0.02), main = 'Proposed Volatility Process', ylab = TeX(r'($V_{t}$)'))
for (i in 1:800){
  lines(ts(out$V[(20000+100*i),], start = 0, frequency = m), col = i)
}

par(mfrow = c(5,3), mar = c(5,4,4,2)+0.1) # default
par(mfrow = c(5,3), mar = c(2,4,2,2), omi = c(0.1,0.1,0.3,0.1))

plot(ts(exp(out$P[-(1:20000),1])), main = '', ylab = TeX(r'($\mu$)'))
acf(exp(out$P[-(1:20000),1]), main = '')
mudens <- density(exp(out$P[-(1:20000),1]))
plot(mudens, main = '', xlab = TeX(r'($\mu$)'))
mtext(TeX(r'(Mercedes-Benz Group SSU Results)'), side = 3, line = , outer = TRUE, cex = 1.5)
abline(v=quantile(exp(out$P[-(1:20000),1]), c(0.025, 0.975)),col=c(3,3))
abline(v=mudens$x[which.max(mudens$y)],col = 2)
points(mudens$x[which.max(mudens$y)], mudens$y[which.max(mudens$y)], col = 2, lwd = 2, pch = 13, cex = 2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds"),
       col=c(2, 3), 
       lty=c(1,1), 
       bty = "n",
       cex=0.65)
plot(ts(exp(out$P[-(1:20000),2])), main = '', ylab = TeX(r'($\kappa$)'))
acf(exp(out$P[-(1:20000),2]), main = '')
kappadens <- density(exp(out$P[-(1:20000),2]))
plot(kappadens, main = '', xlab = TeX(r'($\kappa$)'))
abline(v=quantile(exp(out$P[-(1:20000),2]), c(0.025, 0.975)),col=c(3,3))
abline(v=kappadens$x[which.max(kappadens$y)],col = 2)
points(kappadens$x[which.max(kappadens$y)], kappadens$y[which.max(kappadens$y)], col = 2, lwd = 2, pch = 13, cex = 2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds"),
       col=c(2, 3), 
       lty=c(1,1), 
       bty = "n",
       cex=0.65)
plot(ts(exp(out$P[-(1:20000),3])), main = '', ylab = TeX(r'($\theta$)'))
acf(exp(out$P[-(1:20000),3]), main = '')
thetadens <- density(exp(out$P[-(1:20000),3]))
plot(thetadens, main = '', xlab = TeX(r'($\theta$)'))
abline(v=quantile(exp(out$P[-(1:20000),3]), c(0.025, 0.975)),col=c(3,3))
abline(v=thetadens$x[which.max(thetadens$y)],col = 2)
points(thetadens$x[which.max(thetadens$y)], thetadens$y[which.max(thetadens$y)], col = 2, lwd = 2, pch = 13, cex = 2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds"),
       col=c(2, 3), 
       lty=c(1,1), 
       bty = "n",
       cex=0.65)
plot(ts(exp(out$P[-(1:20000),4])), main = '', ylab = TeX(r'($\sigma$)'))
acf(exp(out$P[-(1:20000),4]), main = '')
sigmadens <- density(exp(out$P[-(1:20000),4]))
plot(sigmadens, main = '', xlab = TeX(r'($\sigma$)'))
abline(v=quantile(exp(out$P[-(1:20000),4]), c(0.025, 0.975)),col=c(3,3))
abline(v=sigmadens$x[which.max(sigmadens$y)],col = 2)
points(sigmadens$x[which.max(sigmadens$y)], sigmadens$y[which.max(sigmadens$y)], col = 2, lwd = 2, pch = 13, cex = 2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds"),
       col=c(2, 3), 
       lty=c(1,1), 
       bty = "n",
       cex=0.65)
plot(ts(out$P[-(1:20000),5]), main = '', ylab = TeX(r'($\rho$)'))
acf(out$P[-(1:20000),5], main = '')
rhodens <- density(out$P[-(1:20000),5])
plot(rhodens, main = '', xlab = TeX(r'($\rho$)'))
abline(v=quantile(out$P[-(1:20000),5], c(0.025, 0.975)),col=c(3,3))
abline(v=rhodens$x[which.max(rhodens$y)],col = 2)
points(rhodens$x[which.max(rhodens$y)], rhodens$y[which.max(rhodens$y)], col = 2, lwd = 2, pch = 13, cex = 2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds"),
       col=c(2, 3), 
       lty=c(1,1), 
       bty = "n",
       cex=0.65)

effectiveSize(exp(out$P[-(1:20000),]))

# mu summary statistics:
mean(exp(out$P[-(1:20000),1]))
median(exp(out$P[-(1:20000),1]))
mudens$x[which.max(mudens$y)]
quantile(exp(out$P[-(1:20000),1]), c(0.025, 0.975))

# kappa summary statistics:
mean(exp(out$P[-(1:20000),2]))
median(exp(out$P[-(1:20000),2]))
kappadens$x[which.max(kappadens$y)]
quantile(exp(out$P[-(1:20000),2]), c(0.025, 0.975))

# theta summary statistics:
mean(exp(out$P[-(1:20000),3]))
median(exp(out$P[-(1:20000),3]))
thetadens$x[which.max(thetadens$y)]
quantile(exp(out$P[-(1:20000),3]), c(0.025, 0.975))

# sigma summary statistics:
mean(exp(out$P[-(1:20000),4]))
median(exp(out$P[-(1:20000),4]))
sigmadens$x[which.max(sigmadens$y)]
quantile(exp(out$P[-(1:20000),4]), c(0.025, 0.975))

# rho summary statistics:
mean(out$P[-(1:20000),5])
median(out$P[-(1:20000),5])
rhodens$x[which.max(rhodens$y)]
quantile(out$P[-(1:20000),5], c(0.025, 0.975))

######################
# MDB
m <- 1
initS <- initSVf(S, NA, m)
ll <- 2.4/sqrt(length(initS))
Sigma <- diag(rep(ll/sqrt(5), 5))
system.time(MBG.out3.1 <- gibbs.mdb.none(1000, initP, initS, Sigma, a, b, m, sigw=0.5, alpha=2, beta=180))
postvar <- var(MBG.out3.1$P)
lambda <- 0.5
Sigmatune <- lambda*(2.40^2)*postvar/5
# Target AP: 28.39%
# CPU approx: 7000s

# Run1
system.time(MBG.out3.2 <- gibbs.mdb.none(10000, initP, initS, Sigmatune, a, b, m, sigw=0.5, alpha=2, beta=180))
length(unique(MBG.out4.2$P[,1]))/10000
# Run 2
system.time(MBG.out3.3 <- gibbs.mdb.none(100000, initP, initS, Sigmatune, a, b, m, sigw=0.5, alpha=2, beta=180))
length(unique(MBG.out3.3$P[,1]))/100000
# Run 3
system.time(MBG.out3.4 <- gibbs.mdb.none(100000, initP, initS, Sigmatune, a, b, m, sigw=0.5, alpha=2, beta=180))
length(unique(MBG.out3.4$P[,1]))/100000
# Run 4
system.time(MBG.out3.5 <- gibbs.mdb.none(100000, initP, initS, Sigmatune, a, b, m, sigw=0.5, alpha=2, beta=180))
length(unique(MBG.out3.5$P[,1]))/100000

# Try MBG.out3.2, MBG.out3.3, MBG.out3.4, MBG.out3.5
out <- MBG.out3.2
length(unique(out$P[,1]))/100000
par(mfrow=c(1,1))
plot(ts(S, start = 0), ylim = c(0, 120), main = 'Proposed Stock Process', ylab = TeX(r'($S_{t}$)'))
for (i in 0:800){
  lines(ts(out$S[(20000+100*i),], start = 0, frequency = m), col = i)
}
lines(ts(S, start = 0), lwd = 2)

plot(ts(out$V[20000,], start = 0, frequency = m), ylim = c(0, 1), main = 'Proposed Volatility Process', ylab = TeX(r'($V_{t}$)'))
for (i in 0:800){
  lines(ts(out$V[(20000+100*i),], start = 0, frequency = m), col = i)
}

par(mfrow = c(5,3), mar = c(5,4,4,2)+0.1) # default
par(mfrow = c(5,3), mar = c(2,4,2,2), omi = c(0.1,0.1,0.3,0.1))
plot(ts(exp(out$P[-(1:20000),1])), main = '', ylab = TeX(r'($\mu$)'))
acf(exp(out$P[-(1:20000),1]), main = '')
mudens <- density(exp(out$P[-(1:20000),1]))
plot(mudens, main = '', xlab = TeX(r'($\mu$)'))
mtext(TeX(r'(Mercedes-Benz Group MDB Results)'), side = 3, line = , outer = TRUE, cex = 1.5)
abline(v=quantile(exp(out$P[-(1:20000),1]), c(0.025, 0.975)),col=c(3,3))
abline(v=mudens$x[which.max(mudens$y)],col = 2)
points(mudens$x[which.max(mudens$y)], mudens$y[which.max(mudens$y)], col = 2, lwd = 2, pch = 13, cex = 2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds"),
       col=c(2, 3), 
       lty=c(1,1), 
       bty = "n",
       cex=0.65)
plot(ts(exp(out$P[-(1:20000),2])), main = '', ylab = TeX(r'($\kappa$)'))
acf(exp(out$P[-(1:20000),2]), main = '')
kappadens <- density(exp(out$P[-(1:20000),2]))
plot(kappadens, main = '', xlab = TeX(r'($\kappa$)'))
abline(v=quantile(exp(out$P[-(1:20000),2]), c(0.025, 0.975)),col=c(3,3))
abline(v=kappadens$x[which.max(kappadens$y)],col = 2)
points(kappadens$x[which.max(kappadens$y)], kappadens$y[which.max(kappadens$y)], col = 2, lwd = 2, pch = 13, cex = 2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds"),
       col=c(2, 3), 
       lty=c(1,1), 
       bty = "n",
       cex=0.65)
plot(ts(exp(out$P[-(1:20000),3])), main = '', ylab = TeX(r'($\theta$)'))
acf(exp(out$P[-(1:20000),3]), main = '')
thetadens <- density(exp(out$P[-(1:20000),3]))
plot(thetadens, main = '', xlab = TeX(r'($\theta$)'))
abline(v=quantile(exp(out$P[-(1:20000),3]), c(0.025, 0.975)),col=c(3,3))
abline(v=thetadens$x[which.max(thetadens$y)],col = 2)
points(thetadens$x[which.max(thetadens$y)], thetadens$y[which.max(thetadens$y)], col = 2, lwd = 2, pch = 13, cex = 2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds"),
       col=c(2, 3), 
       lty=c(1,1), 
       bty = "n",
       cex=0.65)
plot(ts(exp(out$P[-(1:20000),4])), main = '', ylab = TeX(r'($\sigma$)'))
acf(exp(out$P[-(1:20000),4]), main = '')
sigmadens <- density(exp(out$P[-(1:20000),4]))
plot(sigmadens, main = '', xlab = TeX(r'($\sigma$)'))
abline(v=quantile(exp(out$P[-(1:20000),4]), c(0.025, 0.975)),col=c(3,3))
abline(v=sigmadens$x[which.max(sigmadens$y)],col = 2)
points(sigmadens$x[which.max(sigmadens$y)], sigmadens$y[which.max(sigmadens$y)], col = 2, lwd = 2, pch = 13, cex = 2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds"),
       col=c(2, 3), 
       lty=c(1,1), 
       bty = "n",
       cex=0.65)
plot(ts(out$P[-(1:20000),5]), main = '', ylab = TeX(r'($\rho$)'))
acf(out$P[-(1:20000),5], main = '')
rhodens <- density(out$P[-(1:20000),5])
plot(rhodens, main = '', xlab = TeX(r'($\rho$)'))
abline(v=quantile(out$P[-(1:20000),5], c(0.025, 0.975)),col=c(3,3))
abline(v=rhodens$x[which.max(rhodens$y)],col = 2)
points(rhodens$x[which.max(rhodens$y)], rhodens$y[which.max(rhodens$y)], col = 2, lwd = 2, pch = 13, cex = 2)
legend("topright", 
       c("MAP estimate", "95% Credible bounds"),
       col=c(2, 3), 
       lty=c(1,1), 
       bty = "n",
       cex=0.65)

effectiveSize(exp(out$P[-(1:20000),]))

# mu summary statistics:
mean(exp(out$P[-(1:20000),1]))
median(exp(out$P[-(1:20000),1]))
mudens$x[which.max(mudens$y)]
quantile(exp(out$P[-(1:20000),1]), c(0.025, 0.975))

# kappa summary statistics:
mean(exp(out$P[-(1:20000),2]))
median(exp(out$P[-(1:20000),2]))
kappadens$x[which.max(kappadens$y)]
quantile(exp(out$P[-(1:20000),2]), c(0.025, 0.975))

# theta summary statistics:
mean(exp(out$P[-(1:20000),3]))
median(exp(out$P[-(1:20000),3]))
thetadens$x[which.max(thetadens$y)]
quantile(exp(out$P[-(1:20000),3]), c(0.025, 0.975))

# sigma summary statistics:
mean(exp(out$P[-(1:20000),4]))
median(exp(out$P[-(1:20000),4]))
sigmadens$x[which.max(sigmadens$y)]
quantile(exp(out$P[-(1:20000),4]), c(0.025, 0.975))

# rho summary statistics:
mean(out$P[-(1:20000),5])
median(out$P[-(1:20000),5])
rhodens$x[which.max(rhodens$y)]
quantile(out$P[-(1:20000),5], c(0.025, 0.975))