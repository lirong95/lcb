library(EnvStats)
library(MASS)
############################################################################
################### Example 1 in Section 3.1 ###############################
############################################################################
n <- 500
set.seed(123456)
# data generating process
x <- runif(n, min=-6, max=6)
euclid.dist <- dist(x)
euclid.dist <- as.matrix(euclid.dist)
eta <- sin(1.5*x-1) + 0.05*exp(x) + 0.2*x
mu <- pnorm(eta)
y <- rbinom(n, 1, mu)

# plot y with jittering
y_j <- y + runif(n, min=-0.04, max=0.04)
par(mfrow=c(1,4))
plot(x[which(y==0)], y_j[which(y==0)], pch=20, ylab="Y", xlab="x", col="black",
     xlim=c(-6,6),  ylim=c(0,1.3), main="Orignal data", cex = 0.5, cex.lab=1.5)
points(x[which(y==1)], y_j[which(y==1)], pch=1, col="red", cex = 0.5)
lines(x[order(x)], mu[order(x)], pch=20, col="green", lwd=2, cex = 0.5)
legend("topright", c("Y=0", "Y=1"), col=c("black", "red"), pch=c(20,1), cex=1.2, 
       text.width = 1, y.intersp = 0.8, x.intersp = 0.2)

# fit with the assumed model
fit <- glm(y ~ x, family = binomial(link = "probit"))
beta_hat <- fit$coefficients
eta_hat <- -(beta_hat[1] + beta_hat[2] * x)

# surrogate residual
ss <- NULL
for(i in 1:n){
  if(y[i]==0){
    zmean <- fit$coefficients[1] + fit$coefficients[2]*x[i]
    s <- rnormTrunc(1, mean = zmean, sd = 1, min = -Inf, max = 0)
  }
  else{
    zmean <- fit$coefficients[1] + fit$coefficients[2]*x[i]
    s <- rnormTrunc(1, mean = zmean, sd = 1, min = 0, max = Inf)
  }
  ss <- c(ss,s)
}
r <- ss - (fit$coefficients[1] + fit$coefficients[2]*x)
plot(x[which(y==0)], r[which(y==0)], pch=20, ylab="Surrogate residuals", xlab="x", col="black",
     xlim=c(-6,6), ylim=c(-4,4), main="Surrogate residual", cex = 0.5, cex.lab=1.5)
points(x[which(y==1)], r[which(y==1)], pch=1, col="red", cex = 0.5)
lines(x, eta_hat, col = "blue", lwd=2, cex = 0.5)
legend("topright", c("Y=0", "Y=1"), col=c("black", "red"), pch=c(20,1), cex=1.2, 
       text.width = 1, y.intersp = 0.8, x.intersp = 0.2)

# Pearson residual with jittering
pi_hat <- pnorm(fit$coefficients[1]+fit$coefficients[2]*x)
pr <- (y-pi_hat)/sqrt(pi_hat*(1-pi_hat))
pr_j <- NULL
pr_j[which(y==1)] <- pr[which(y==1)] + runif(length(which(y==1)), min=0, max=0.6)
pr_j[which(y==0)] <- pr[which(y==0)] + runif(length(which(y==0)), min=-0.6, max=0)
plot(x[which(y==0)], pr_j[which(y==0)], pch=20, ylab="Pearson residuals", xlab="x", col="black",
     xlim=c(-6,6), ylim=c(-4,4), main="Pearson residual", cex = 0.5, cex.lab=1.5)
points(x[which(y==1)], pr_j[which(y==1)], pch=1, col="red", cex = 0.5)
legend("topright", c("Y=0", "Y=1"), col=c("black", "red"), pch=c(20,1), cex=1.2, 
       text.width = 1, y.intersp = 0.8, x.intersp = 0.2)


# SBS residual with jittering
library(PResiduals)
rsbs <- presid(fit)
rsbs_j <- NULL
rsbs_j[which(y==1)] <- rsbs[which(y==1)] + runif(length(which(y==1)), min=0, max=0.6)
rsbs_j[which(y==0)] <- rsbs[which(y==0)] + runif(length(which(y==0)), min=-0.6, max=0)
plot(x[which(y==0)], rsbs_j[which(y==0)], pch=20, ylab="SBS residuals", xlab="x", col="black",
     xlim=c(-6,6), ylim=c(-4,4), main="SBS residual", cex = 0.5, cex.lab=1.5)
points(x[which(y==1)], rsbs_j[which(y==1)], pch=1, col="red", cex = 0.5)
legend("topright", c("Y=0", "Y=1"), col=c("black", "red"), pch=c(20,1), cex=1.2, 
       text.width = 1, y.intersp = 0.8, x.intersp = 0.2)

# Parametric bootstrapped repsonse with jittering
eta_para <- fit$coefficients[1] + fit$coefficients[2] * x
mu_para <- pnorm(eta_para)
y_para <- rbinom(n, 1, mu_para)
y_para_j2 <- y_para + runif(n, min=-0.04, max=0.04)
plot(x[which(y==0)], y_para_j2[which(y==0)], pch=17, ylab="Bootstrapped Y*", xlab="x",
     xlim=c(-6,6), ylim=c(0,1.3), col="black", main="Bootstrapped data \n Parametric bootstrap", cex = 0.7, cex.lab=1.5)
points(x[which(y==1)], y_para_j2[which(y==1)], pch=2, col="red", cex=0.7)
lines(x[order(x)], mu_para[order(x)], pch=20, col="orange", lwd=2, cex = 0.5)
legend("topright", c("Y=0", "Y=1"), col=c("black", "red"), pch=c(17,2), cex=1.2, 
       text.width = 1, y.intersp = 0.8, x.intersp = 0.2)


# Local residual bootstrapped surrogate residuals and response
delta <- 4
r_b <- NULL
r_mean <- r 
for(i in 1:n){
  ind <- order(euclid.dist[i,])[1:delta]
  r_b[i] <- sample(r_mean[ind], size = 1, replace = T)
}
s_b <- r_b + fit$coefficients[1] + fit$coefficients[2]*x
y_b <- NULL
for(i in 1:n){
  if(s_b[i]<=0){
    yy = 0}else{
      yy = 1
    }
  y_b <- c(y_b, yy)
}
plot(x[which(y==0)], r_b[which(y==0)], pch=15, ylab="Bootstrapped residuals", xlab="x", xlim=c(-6,6), ylim=c(-4,4),
     col="black", main="Bootstrapped residuals \n LRB with Surrogate residual", cex=0.7, cex.lab=1.5)
points(x[which(y==1)], r_b[which(y==1)], pch=0, col="red", cex=0.7)
lines(x, eta_hat, col = "blue", lwd=2, cex = 0.5)
legend("topright", c("Y=0", "Y=1"), col=c("black", "red"), pch=c(15,0), cex=1.2, 
       text.width = 1, y.intersp = 0.8, x.intersp = 0.2)
y_b_j <- y_b + runif(n, min=-0.04, max=0.04)
plot(x[which(y==0)], y_b_j[which(y==0)], pch=17, ylab="Bootstrapped Y*", xlab="x",
     xlim=c(-6,6), ylim=c(0,1.3), col="black", main="Bootstrapped data \n LRB with Surrogate residual", cex = 0.7, cex.lab=1.5)
points(x[which(y==1)], y_b_j[which(y==1)], pch=2, col="red", cex=0.7)
legend("topright", c("Y=0", "Y=1"), col=c("black", "red"), pch=c(17,2), cex=1.2, 
       text.width = 1, y.intersp = 0.8, x.intersp = 0.2)

# Local residual bootstrapped Pearson residuals and response
r_mean <- pr_j
delta <- 4
r_b <- NULL
for(i in 1:n){
  ind <- order(euclid.dist[i,])[1:delta]
  r_b[i] <- sample(r_mean[ind], size = 1, replace = T)
}
y_pear <- pi_hat + sqrt(pi_hat*(1-pi_hat))*r_b
y_p <- round(y_pear)
plot(x[which(y==0)], r_b[which(y==0)], pch=15, ylab="Bootstrapped residuals", xlab="x", xlim=c(-6,6), ylim=c(-4,4),
     col="black", main="Bootstrapped residuals \n LRB with Pearson residual", cex=0.7, cex.lab=1.5)
points(x[which(y==1)], r_b[which(y==1)], pch=0, col="red", cex=0.7)
legend("topright", c("Y=0", "Y=1"), col=c("black", "red"), pch=c(15,0), cex=1.2, 
       text.width = 1, y.intersp = 0.8, x.intersp = 0.2)
y_p_j <- y_p + runif(n, min=-0.04, max=0.04)
plot(x[which(y==0)], y_p_j[which(y==0)], pch=17, ylab="Bootstrapped Y*", xlab="x",
     xlim=c(-6,6), ylim=c(0,1.3), col="black", main="Bootstrapped data \n LRB with Pearson residual", cex = 0.7, cex.lab=1.5)
points(x[which(y==1)], y_p_j[which(y==1)], pch=2, col="red", cex=0.7)
legend("topright", c("Y=0", "Y=1"), col=c("black", "red"), pch=c(17,2), cex=1.2, 
       text.width = 1, y.intersp = 0.8, x.intersp = 0.2)

# Local residual bootstrapped SBS residuals and response
s_mean <- rsbs_j
delta <- 4
r_b <- NULL
for(i in 1:n){
  ind <- order(euclid.dist[i,])[1:delta]
  #ind <- which(y==y[i])
  r_b[i] <- sample(s_mean[ind], size = 1, replace = T)
}
y_sbs <- NULL
for(i in 1:n){
  if(r_b[i]<=0) yy = 0
  else yy = 1
  y_sbs <- c(y_sbs, yy)
}
plot(x[which(y==0)], r_b[which(y==0)], pch=15, ylab="Bootstrapped residuals", xlab="x", xlim=c(-6,6), ylim=c(-4,4),
     col="black", main="Bootstrapped residuals \n LRB with SBS residual", cex=0.7, cex.lab=1.5)
points(x[which(y==1)], r_b[which(y==1)], pch=0, col="red", cex=0.7)
legend("topright", c("Y=0", "Y=1"), col=c("black", "red"), pch=c(15,0), cex=1.2, 
       text.width = 1, y.intersp = 0.8, x.intersp = 0.2)
y_sbs_j <- y_sbs + runif(n, min=-0.04, max=0.04)
plot(x[which(y==0)], y_sbs_j[which(y==0)], pch=17, ylab="Bootstrapped Y*", xlab="x",
     xlim=c(-6,6), ylim=c(0,1.3), col="black", main="Bootstrapped data \n LRB with SBS residual", cex = 0.7, cex.lab=1.5)
points(x[which(y==1)], y_sbs_j[which(y==1)], pch=2, col="red", cex=0.7)
legend("topright", c("Y=0", "Y=1"), col=c("black", "red"), pch=c(17,2), cex=1.2, 
       text.width = 1, y.intersp = 0.8, x.intersp = 0.2)

############################################################################
######### Example in Section S3.3 in supplementary materials ###############
############################################################################
source("functions.R")
# data generating process
n <- 2000
set.seed(1111)
times <- 100
B = 500
x <- runif(n, min=-6, max=6)
beta1 = 2
alpha = 12
beta2 = -2
eta <- alpha + beta1 * x + beta2 * x^2
mu <- pnorm(eta)
y <- rbinom(n, 1, mu)

# fit with the assumed model
euclid.dist <- dist(x)
euclid.dist <- as.matrix(euclid.dist)
fit.mis <- glm(y ~ x, family = binomial(link = "probit"))
# setting intial for all samples 
l_itera <- round(n^(1/3))
res_whole <- surr.boot(y, x, B, fit.mis, l_itera, euclid.dist)
sd.hat <- apply(res_whole$beta.surr, 2, sd)

# set candidate size and calculate the mse under each candidate
times <- 5
candidate.size <- c(2,4,6,8,10,12,14,16)
beta_cand <- array(NA, dim=c(times, 2, length(candidate.size)))
for(m in 1:times){
  ind.sample <- sample(1:n, size=0.9*n, replace = F)
  y.sample <- y[ind.sample]
  x.sample <- x[ind.sample]
  euclid.dist.sample <- dist(x.sample)
  euclid.dist.sample <- as.matrix(euclid.dist.sample)
  fit.mis.sample <- glm(y.sample ~ x.sample, family = binomial(link = "probit"))
  beta_hat.sample <- fit.mis.sample$coefficients
  ss <- NULL
  for(i in 1:length(ind.sample)){
    if(y.sample[i]==0){
      zmean <- beta_hat.sample[1] + beta_hat.sample[2] * x.sample[i]
      s <- rnormTrunc(1, mean = zmean, sd = 1, min = -Inf, max = 0)
    }
    else{
      zmean <- beta_hat.sample[1] + beta_hat.sample[2] * x.sample[i]
      s <- rnormTrunc(1, mean = zmean, sd = 1, min = 0, max = Inf)
    }
    ss <- c(ss,s)
  }
  r <- ss - (beta_hat.sample[1] + beta_hat.sample[2] * x.sample)
  r_mean <- r
  j = 1
  for(l_cand in candidate.size){
    print(c(m, l_cand))
    beta.surr <- NULL
    for(b in 1:B){
      r_b <- apply(euclid.dist.sample, 1, function(c) sample(r_mean[order(c)[1:l_cand]], size = 1))
      s_b <- r_b + beta_hat.sample[1] + beta_hat.sample[2] * x.sample
      y_surr <- ifelse(s_b>0, 1, 0)
      fit_surr <- glm(y_surr ~ x.sample, family = binomial(link = "probit"))
      beta.surr <- rbind(beta.surr, fit_surr$coefficients)
    }
    beta_cand[m,,j] <- apply(beta.surr, 2, sd)
    j = j + 1
  }
}
beta_mse <- matrix(NA, nrow=times, ncol=length(candidate.size))
for(j in 1:length(candidate.size)){
  for(i in 1:times){
    beta_mse[i,j] <- crossprod(beta_cand[i,,j] - sd.hat)
  }
}
l_itera_update <- candidate.size[which.min(apply(beta_mse, 2, mean))]


############################################################################
########################### Example in Section 5 ###########################
############################################################################
source("functions.R")
n <- 2000
set.seed(1111)
B = 500
x <- runif(n, min=-6, max=6)
euclid.dist <- dist(x)
euclid.dist <- as.matrix(euclid.dist)
beta1 = 2
alpha = 12
beta2 = -2
eta <- alpha + beta1 * x + beta2 * x^2
mu <- pnorm(eta)
y <- rbinom(n, 1, mu)
fit.mis <- glm(y ~ x, family = binomial(link = "probit"))
size.cand <- c(2,4,6,10)
delta = selec_neigsize(y, x, B, fit.mis, euclid.dist, size.cand, subsample=5, m.prop=0.9)
surr.res <- surr.boot(y, x, B, fit.mis, delta, euclid.dist)
