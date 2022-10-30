library(pracma)

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

num_samples <- 10000
mu <- 0.15
sigma <- 0.25
v <- 0.20
den.vec <- rep(NA, num_samples-1)
y.vec <- rep(NA, num_samples-1)
strt <- 0
ed <- 0.5
diff <- 1/num_samples
for(i in 2:num_samples){ 
    y.vec[i-1] <- strt + (i-1)*diff
} 

for(i in 2:num_samples){
    den.vec[i-1] <- Pdf(mu, sigma, v, y.vec[i-1])
}

# print(y.vec)
# print(den.vec)
plot(y.vec,den.vec, type='l')