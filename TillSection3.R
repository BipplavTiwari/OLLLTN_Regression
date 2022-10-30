library(pracma)
library(stats)
library(copula)
library(conflicted)
# Functions for pdf, quantile, sampling, moment, skewness and kurtosis
G <- function(mu, sigma, y){
  q <- (logit(y)-logit(mu))/sigma
  eita <- pnorm(q, 0.0, 1.0)
  return (eita)
}

Pdf <- function(mu, sigma, v, y){
  p1 <- v/(y*(1-y)*(2*pi*(sigma)^2)^0.5)
  p2 <- exp(-((logit(y)-logit(mu))^2/(2*sigma^2)))
  eita <- G(mu, sigma, y)
  p3 <- (eita*(1-eita))^(v-1)
  p4 <- (eita^v + (1-eita)^v)^(-2)
  return (p1*p2*p3*p4)
}
Quantile <-function(u, mu, sigma, v){
  u1 <- u^(1/v)
  u2 <- (1-u)^(1/v)
  u <- u1/(u1+u2)
  x <- logit(mu) + sqrt(2)*sigma*erfinv(2*u-1)
  q_ollltn <- 1/(1+exp(-x))
  return (q_ollltn)
}

moment <-function(k, mu, sigma ,v ) {
  integrand <- function(u) {
    return (Quantile(u,mu,sigma,v)^k)
  }
  return (integrate(integrand, 0, 1)$value)
}

skewness <-function(mu, sigma, v) {
  num<- Quantile(1/4,mu,sigma,v) + Quantile(3/4,mu,sigma,v) - 2*Quantile(1/2,mu,sigma,v)
  den <- Quantile(3/4, mu, sigma, v) - Quantile(1/4, mu, sigma, v)
  return (num/den)
}

kurtosis <- function(mu, sigma, v) {
  num <- Quantile(7/8, mu, sigma , v) - Quantile(5/8, mu, sigma , v) +Quantile(3/8, mu, sigma , v) -Quantile(1/8, mu, sigma , v)
  den <- Quantile(6/8, mu, sigma , v) - quantile(2/8, mu, sigma , v)
  return (num/den)
}
sample <-function(mu, sigma, v) {
  u <- runif(1)
  return (Quantile(u,mu,sigma,v))
}


print(sample(0.5,0.5,0.25))
# num_samples <- 2
# mu <- 0.5
# sigma <- 0.5
# v <- 0.25
# den.vec <- rep(NA, num_samples-1)
# y.vec <- rep(NA, num_samples-1)
# quantile.vec <-rep(NA, num_samples-1)
# skewness.mat <-matrix(NA, num_samples-1, num_samples-1)
# skewness.vec <-rep(NA, num_samples-1)
# kurtosis.vec <-rep(NA, num_samples-1)
# kurtosis.mat <-matrix(NA, num_samples-1, num_samples-1)
# sigma.vec <- seq(0.1,1, length.out = num_samples-1)
# nu.vec <- seq(0.1,1, length.out = num_samples-1)
# strt <- 0
# ed <- 0.5
# diff <- 1/num_samples
# for(i in 2:num_samples){ 
#   y.vec[i-1] <- strt + (i-1)*diff
# } 

# for(i in 2:num_samples){
#   den.vec[i-1] <- Pdf(mu, sigma, v, y.vec[i-1])
# }
# for(i in 2:num_samples) {
#   quantile.vec[i-1] <-Quantile(runif(1), mu, sigma, v)
# }
# sample.vec  <- rep(NA, num_samples-1)
# for(i in 2:num_samples) {
#   sample.vec[i-1] <-sample(mu, sigma, v)
# }
# for(i in 2:num_samples) {
#   for(j in 2:num_samples)
#     skewness.mat[i-1][j-1] <-skewness(mu, sigma.vec[i-1], nu.vec[j-1])
#     kurtosis.mat[i-1][j-1] <-kurtosis(mu, sigma.vec[i-1], nu.vec[j-1])
  
# }
# z <- outer(sigma.vec, nu.vec, function(a,b) skewness(0.1,a,b))
# # print(y.vec)
# # print(den.vec)
# plot(y.vec,den.vec, type='l')

