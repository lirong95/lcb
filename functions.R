library(EnvStats)
library(MASS)
#### local residual bootstrap with surrogate residual for probit model #####
#x: predictors
#y: reponse
#B: the times of bootstrap
#fit.mis: the fitted model (e.g. output from glm function)
#delta: neighborhood size
#euclid.dist: Euclidean distance of predictors
surr.boot <- function(y, x, B, fit.mis, delta, euclid.dist){
  beta_hat <- fit.mis$coefficients
  nn <- length(y)
  ##get surrogate residual
  ss <- NULL
  for(i in 1:nn){
    if(y[i]==0){
      zmean <- beta_hat[1] + beta_hat[2] * x[i]
      s <- rnormTrunc(1, mean = zmean, sd = 1, min = -Inf, max = 0)
    }
    else{
      zmean <- beta_hat[1] + beta_hat[2] * x[i]
      s <- rnormTrunc(1, mean = zmean, sd = 1, min = 0, max = Inf)
    }
    ss <- c(ss,s)
  }
  r <- ss - (beta_hat[1] + beta_hat[2] * x)
  r_mean <- r
  beta.surr <- y.surr <- NULL
  for(b in 1:B){
    r_b <- apply(euclid.dist, 1, function(c) sample(r_mean[order(c)[1:delta]], size = 1))
    ##generate bootstrap y from resampling surrogate residual
    s_b <- r_b + beta_hat[1] + beta_hat[2] * x
    y_surr <- ifelse(s_b>0, 1, 0)
    fit_surr <- glm(y_surr ~ x, family = binomial(link = "probit"))
    beta.surr <- rbind(beta.surr, fit_surr$coefficients)
    y.surr <- cbind(y.surr, y_surr)
  }
  return(list(beta.surr=beta.surr, y.surr=y.surr))
}
#beta.surr: B*p matrix of bootstrapped coefficient estimate
#y.surr: n*B matrix of bootstrapped response

#### local residual bootstrap with surrogate residual for logistic model #####
surr.log.boot <- function(y, x, B, fit.mis, delta, euclid.dist){
  beta_hat <- fit.mis$coefficients
  ss <- NULL
  for(i in 1:n){
    if(y[i]==0){
      zmean <- beta_hat[1] + beta_hat[2]*x[i]
      usamp <- runif(1, min = 0, max = exp(-zmean)/(1+exp(-zmean)))
      s <- log(usamp/(1-usamp)) + zmean
    }
    else{
      zmean <- beta_hat[1] + beta_hat[2]*x[i]
      usamp <- runif(1, min = exp(-zmean)/(1+exp(-zmean)), max = 1)
      s <- log(usamp/(1-usamp)) + zmean
    }
    ss <- c(ss,s)
  }
  r <- ss - (beta_hat[1] + beta_hat[2]*x)
  r_mean <- r
  beta.surr <- y.surr <- NULL
  for(b in 1:B){
    r_b <- apply(euclid.dist, 1, function(c) sample(r_mean[order(c)[1:delta]], size = 1))
    s_b <- r_b + (beta_hat[1] + beta_hat[2] * x) 
    y_surr <- ifelse(s_b>0, 1, 0)
    fit_surr <- glm(y_surr ~ x, family = binomial(link = "logit"))
    beta.surr <- rbind(beta.surr, fit_surr$coefficients)
    y.surr <- cbind(y.surr, y_surr)
  }
  return(list(beta.surr=beta.surr, y.surr=y.surr))
}


#### local residual bootstrap with pearson residual for poisson model #####
pearson.poi.boot <- function(y, x, B, fit.mis, delta, euclid.dist){
  beta_hat <- fit.mis$coefficients
  pi_hat <- exp(beta_hat[1] + beta_hat[2]*x)
  pi_hat <- as.vector(pi_hat)
  pr <- (y-pi_hat)/sqrt(pi_hat)
  r_mean <- pr
  beta.pearson <- y.pearson <- NULL
  for(b in 1:B){
    r_b <- apply(euclid.dist, 1, function(c) sample(r_mean[order(c)[1:delta]], size=1, replace=T))
    y_pear <- pi_hat + sqrt(pi_hat)*r_b
    y_p <- round(y_pear)
    y_p[which(y_p < 0)] <- 0 
    fit_pearson <- glm(y_p ~ x, family = poisson(link = "log"))
    beta.pearson <- rbind(beta.pearson, fit_pearson$coefficients)
    y.pearson <- cbind(y.pearson, y_p)
  }
  return(list(beta.pearson=beta.pearson, y.pearson=y.pearson))
}

#### local residual bootstrap with pearson residual for gamma model #####
pearson.gam.boot <- function(y, x, B, fit.mis, delta, euclid.dist){
  beta_hat <- fit.mis$coefficients
  eta_hat <- -(beta_hat[1] + beta_hat[2] * x)
  pi_hat <- -1/(eta_hat)
  pi_hat <- as.vector(pi_hat)
  pr <- (y-pi_hat)/pi_hat
  r_mean <- pr
  beta.pearson <- y.pearson <- NULL
  for(b in 1:B){
    r_b <- apply(euclid.dist, 1, function(c) sample(r_mean[order(c)[1:delta]], size=1, replace=T))
    y_pear <- pi_hat + pi_hat*r_b
    y_pear[which(y_pear <= 0)] <- 1e-4
    fit_pearson <- glm(y_pear ~ x, family = Gamma(link = "inverse"))
    beta.pearson <- rbind(beta.pearson, fit_pearson$coefficients)
    y.pearson <- cbind(y.pearson, y_pear)
  }
  return(list(beta.pearson=beta.pearson, y.pearson=y.pearson))
}


#### selection of neighborhood size 
selec_neigsize <- function(y, x, B, fit.mis, euclid.dist, size.cand, subsample, m.prop){
  #x: predictors
  #y: reponse
  #B: the times of bootstrap
  #fit.mis: the fitted model (e.g. output from glm function)
  #euclid.dist: Euclidean distance of predictors
  #subsample: number of subsamples
  #size.cand: candidate of neighborhood size
  #m.prop: the observations proportion sampled from the entire data set
  beta.cand <- array(NA, dim=c(length(size.cand), 2, subsample))
  for(m in 1:subsample){
    print(paste("This the subsample", m))
    ind.sample <- sample(1:length(y), size=m.prop*length(y), replace = F)
    y.sample <- y[ind.sample]
    x.sample <- x[ind.sample]
    euclid.dist.sample <- dist(x.sample)
    euclid.dist.sample <- as.matrix(euclid.dist.sample)
    fit.mis.sample <- glm(y.sample ~ x.sample, family = binomial(link = "probit"))
    j <- 1
    for(l_cand in size.cand){
      res.surr.sample <- surr.boot(y.sample, x.sample, B, fit.mis.sample, delta=l_cand, euclid.dist.sample)
      beta.cand[j,,m] <- apply(res.surr.sample$beta.surr, 2, sd)
      j <- j +1
    }
  }
  diff <- 1; l_iter <- round(length(y)^(1/3))
  if(diff != 0){
    res_whole <- surr.boot(y, x, B, fit.mis, l_iter, euclid.dist)
    sd.hat <- apply(res_whole$beta.surr, 2, sd)
    beta.diff <- NULL
    for(m in 1:subsample){
      beta.diff<- rbind(beta.diff, apply(beta.cand[,,m], 1, function(c) crossprod(c-sd.hat)))
    }
    l_new <- round((1/m.prop)^(1/3)*size.cand[which(apply(beta.diff, 2, mean)==min(apply(beta.diff, 2, mean)))])
    diff <- abs(l_new - l_iter)
    l_iter <- l_new
  }
  return(neig_size=l_iter)
}

