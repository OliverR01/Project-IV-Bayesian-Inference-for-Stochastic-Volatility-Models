# Monte Carlo Option Pricing Code File 1 

# This file performs Monte Carlo Option Pricing as seen in Chapter 5. 

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
# - data_gen --------- Simulate Synthetic Data

# heston_model_1.R
# - deltaSV ---------- Proposes next step in (S,V) process
# - data.generator --- Simulates Synthetic Data

# monte_carlo_op.R
# - euro_call -------- Calculates the payoff of a European Call option
# - asian_call ------- Calculates the payoff of an Asian Call option
# - mc_BS ------------ Calculates risk-neutral price under Black-Scholes model
# - mc_H ------------- Calculates risk-neutral price under Heston model

################################################################################
## New Functions:

euro_call <- function(Svec, K, Timevec){
  T <- Timevec[1]
  ST <- Svec[T+1]
  payoff <- max(ST-K,0)
  return(payoff)
}

asian_call <- function(Svec, K, Timevec){
  T <- Timevec[1]; T0 <- Timevec[2]
  avrg_period <- Svec[(T0+1):(T+1)]
  A <- sum(avrg_period)/(T-T0)
  payoff <- max(A-K, 0)
  return(payoff)
}

mc_BS <- function(N, t, S0, r, sigma, K, Timevec, phi, deltat = 1){
  T <- Timevec[1]; T0 <- Timevec[2]
  val <- rep(NA, N)
  for (i in 1:N){
    sim <- data_gen(T-t, deltat, S0, r, sigma)
    val[i] <- phi(sim, K, Timevec)
  }
  phibar <- mean(val)
  price <- exp(-r*(T-t))*phibar
  sd <- sqrt(sum((val[i]-phibar)^2)/(N-1))
  return(list(price = price,
              sd = sd))
}

mc_H <- function(N, t, S0, V0, r, Tvec, K, Timevec, phi, deltat = 1){
  kappa <- Tvec[2]; theta <- Tvec[3]; sigma <- Tvec[4]; rho <- Tvec[5]; lambda <- Tvec[6]
  T <- Timevec[1]; T0 <- Timevec[2]
  RNTvec <- c(r, kappa+lambda, kappa*theta/(kappa+lambda), sigma, rho)
  val <- rep(NA, N)
  for (i in 1:N){
    sim <- data.generator(T-t, deltat, S0, V0, RNTvec)
    val[i] <- phi(sim[,1], K, Timevec)
  }
  phibar <- mean(val)
  price <- exp(-r*(T-t))*phibar
  sd <- sqrt(sum((val[i]-phibar)^2)/(N-1))
  return(list(price = price,
              sd = sd))
}

################################################################################
## Output:

# Example 1:
N = 10000
t = 0
S0 = 10
r = 0.05
sigma = 0.02
K = 50
T = 60
deltat = 1

set.seed(1)
system.time(mc1 <- mc_BS(N = 10000,
                         t = 0,
                         S0 = 10,
                         r = 0.05,
                         sigma = 0.02,
                         K = 50,
                         T = 60,
                         euro_call))

mc1$price
mc1$sd

d1 <- (log(S0/K) + (r + 0.5*sigma*sigma)*(T-t))/(sigma*sqrt(T-t))
d2 <- (log(S0/K) + (r - 0.5*sigma*sigma)*(T-t))/(sigma*sqrt(T-t))
act1 <- S0*pnorm(d1) - K*exp(-r*(T-t))*pnorm(d2)

act1
100*abs(act1-mc1$price)/act1

set.seed(1)
system.time(mc2 <- mc_BS(10000, 0, 10, 0.05, 0.02, 50, c(60, 10), asian_call))

mc2$price
mc2$sd

################################################################################
# Example 2:
N = 10000
t = 0
S0 = 10
V0 = 1
r = 0.05
Tvec = c(0.1, 0.75, 0.02, 0.2, 0.5, 0) # lambda included as the 6th term 
K = 50
Timevec = 60
deltat = 1

set.seed(1)
system.time(mc3 <- mc_H(N, t, S0, V0, r, Tvec, K, Timevec, euro_call, deltat = 1))

mc3$price
mc3$sd

Timevec = c(60,10)
set.seed(1)
system.time(mc4 <- mc_H(N, t, S0, V0, r, Tvec, K, Timevec, asian_call, deltat = 1))

mc4$price
mc4$sd

########
t = 0
S0 = 10
V0 = 1
r = 0.05
Tvec = c(0.1, 0.75, 0.02, 0.2, 0.5, 0)
K = 50
deltat = 1
sigma = 0.02
T = 60
RNTvec = c(0.05, 0.75, 0.02, 0.5, 0)

set.seed(1)
par(mfrow = c(1,2))
plot(ts(data_gen(T-t, deltat, S0, r, sigma), start = 0), main = TeX(r'(Black-Scholes $S_t$ Processes)'), ylab = TeX(r'($S_t$)'))
for (i in 1:100){
  lines(ts(data_gen(T-t, deltat, S0, r, sigma), start = 0), col = i)
}
abline(h=K, col = 'red', lty = 2, lwd = 2)
plot(ts(data.generator(T-t, deltat, S0, V0, RNTvec)[,1], start = 0), main = TeX(r'(Heston $S_t$ Processes)'), ylab = TeX(r'($S_t$)'), ylim = c(0,800))
for (i in 1:150){
  data <- data.generator(T-t, deltat, S0, V0, RNTvec)
  x <- 0
  for (j in 1:nrow(data)){
    if(data[j,1]<0){
      x <- x+1
    }
  }
  if(x==0){
    lines(ts(data[,1], start = 0), col = i)
  }
}
abline(h=K, col = 'red', lty = 2, lwd = 2)
