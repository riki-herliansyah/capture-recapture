require(TMB)
require(TMBhelper)
require(statmod)
require(boot)
require(microbenchmark)
path <- "your working directory"
setwd(path)
#make sure you put all files including .cpp (C++ files) within the same folder

#Mh Type models
compile(paste0('Mh_type', '.cpp'))
dyn.load(paste0('Mh_type'))
source("mtbh.R")

#preparing data 
x <- read.table("golftees.txt", quote="\"", comment.char="")
x <- as.matrix(x) 
#fitting Mh model using Laplace approximations and GHQ
#model consists of Mh, Mth, Mbh and Mtbh
mtbh(x, method = 'LA2', model = 'Mh', niter = 5)
mtbh(x, method = 'LA4', model = 'Mh', niter = 5)
mtbh(x, method = 'GHQ', model = 'Mh', niter = 5)

#cI for golf tees data using non-paramyeric bootstrap
B=1000
pN1 = data.frame(alpha=rep(0,B), log_sigma=rep(0,B), N=rep(0,B)) 
pN2 = data.frame(alpha=rep(0,B), log_sigma=rep(0,B), N=rep(0,B)) 
pN3 = data.frame(alpha=rep(0,B), log_sigma=rep(0,B), N=rep(0,B)) 
for (p in 1:B) {
  sampel = sample(1:n,replace=TRUE)
  xb = x[sampel, ]
  pN1[p,] =  mtbh(xb, method = 'LA2', model = 'Mh', niter = 5)$par[1:3]
  pN2[p,] =  mtbh(xb, method = 'LA4', model = 'Mh', niter = 5)$par[1:3]
  pN3[p,] =  mtbh(xb, method = 'GHQ', model = 'Mh', npoints = 50)$par[1:3]
}
sapply(pN1, quantile, 0.025)
sapply(pN1, quantile, 0.975)
sapply(pN2, quantile, 0.025)
sapply(pN2, quantile, 0.975)
sapply(pN3, quantile, 0.025)
sapply(pN3, quantile, 0.975)

#the first simulation study in Section 4.1 
#we use a nonparametric Bootstrap to compute the confidence interval
source("sim_closed.R")
nsim = 100
lower.la2 = data.frame(log_sigma=rep(0,nsim), N=rep(0,nsim)) 
upper.la2 = data.frame(log_sigma=rep(0,nsim), N=rep(0,nsim)) 
param.la2 = data.frame(log_sigma=rep(0,nsim), N=rep(0,nsim)) 
lower.la4 = data.frame(log_sigma=rep(0,nsim), N=rep(0,nsim)) 
upper.la4 = data.frame(log_sigma=rep(0,nsim), N=rep(0,nsim)) 
param.la4 = data.frame(log_sigma=rep(0,nsim), N=rep(0,nsim)) 
lower.ghq = data.frame(log_sigma=rep(0,nsim), N=rep(0,nsim)) 
upper.ghq = data.frame(log_sigma=rep(0,nsim), N=rep(0,nsim)) 
param.ghq = data.frame(log_sigma=rep(0,nsim), N=rep(0,nsim))

set.seed(2345)
parameter = list (alpha=rep(-1.0, 6), N=100, sigma=0.75, beta=0.75)
for (k in 1:nsim) {
  x <- sim_closed(parameter = parameter, T=6, behaviour=FALSE)
  n <- dim(x)[1]
  B=1000
  pN1 = data.frame(log_sigma=rep(0,B), N=rep(0,B)) 
  pN2 = data.frame(log_sigma=rep(0,B), N=rep(0,B)) 
  pN3 = data.frame(log_sigma=rep(0,B), N=rep(0,B)) 
  for (p in 1:B) {
    sampel = sample(1:n, replace=TRUE)
    xb = x[sampel, ]
    pN1[p,] =  mtbh(xb, method = 'LA2', model = 'Mh', niter = 5)$par[2:3]
    pN2[p,] =  mtbh(xb, method = 'LA4', model = 'Mh', niter = 5)$par[2:3]
    pN3[p,] =  mtbh(xb, method = 'GHQ', model = 'Mh', npoints = 50)$par[2:3]
  }
  lower.la2[k,] = sapply(pN1, quantile, 0.025)
  upper.la2[k,] = sapply(pN1, quantile, 0.975)
  lower.la4[k,] = sapply(pN2, quantile, 0.025)
  upper.la4[k,] = sapply(pN2, quantile, 0.975)
  lower.ghq[k,] = sapply(pN3, quantile, 0.025)
  upper.ghq[k,] = sapply(pN3, quantile, 0.975)
  param.la2[k,] = mtbh(x, method = 'LA2', model = 'Mh', niter = 5)$par[2:3]
  param.la4[k,] = mtbh(x, method = 'LA4', model = 'Mh', niter = 5)$par[2:3]
  param.ghq[k,] = mtbh(x, method = 'GHQ', model = 'Mh', npoints = 50)$par[2:3]
}

mean((param.la2[,2]- N)/N)
mean((param.la4[,2]- N)/N)
mean((param.ghq[,2]- N)/N)
mean((exp(param.la2[,1])- sigma)/sigma)
mean((exp(param.la4[,1])- sigma)/sigma)
mean((exp(param.ghq[,1])- sigma)/sigma)
cp = 0
for (i in 1:100) {
  if (N >= lower.la4[i,2] && N <= upper.la4[i,2]){
    cp = cp + 1
  } else {
    cp = cp + 0
  }
} 
cp
cp = 0
for (i in 1:100) {
  if (sigma >= exp(lower.la4[i,1]) && sigma <= exp(upper.la4[i,1])){
    cp = cp + 1
  } else {
    cp = cp + 0
  }
} 
cp

#the second simulation study in Section 4.1 
#we use a nonparametric Bootstrap to compute the confidence interval
#we set heterogeneity variability as sigma=1, 1.5 and 2.0
set.seed(2345)
parameter = list (alpha=rep(-1.5, 8), N=250, sigma=1.0, beta=NULL)
for (k in 1:nsim) {
  x <- sim_closed(parameter = parameter, T=8, behaviour=FALSE)
  n <- dim(x)[1]
  B=1000
  pN1 = data.frame(log_sigma=rep(0,B), N=rep(0,B)) 
  pN2 = data.frame(log_sigma=rep(0,B), N=rep(0,B)) 
  pN3 = data.frame(log_sigma=rep(0,B), N=rep(0,B)) 
  for (p in 1:B) {
    sampel = sample(1:n, replace=TRUE)
    xb = x[sampel, ]
    pN1[p,] =  mtbh(xb, method = 'LA2', model = 'Mh', niter = 5)$par[2:3]
    pN2[p,] =  mtbh(xb, method = 'LA4', model = 'Mh', niter = 5)$par[2:3]
    pN3[p,] =  mtbh(xb, method = 'GHQ', model = 'Mh', npoints = 50)$par[2:3]
  }
  lower.la2[k,] = sapply(pN1, quantile, 0.025)
  upper.la2[k,] = sapply(pN1, quantile, 0.975)
  lower.la4[k,] = sapply(pN2, quantile, 0.025)
  upper.la4[k,] = sapply(pN2, quantile, 0.975)
  lower.ghq[k,] = sapply(pN3, quantile, 0.025)
  upper.ghq[k,] = sapply(pN3, quantile, 0.975)
  param.la2[k,] = mtbh(x, method = 'LA2', model = 'Mh', niter = 5)$par[2:3]
  param.la4[k,] = mtbh(x, method = 'LA4', model = 'Mh', niter = 5)$par[2:3]
  param.ghq[k,] = mtbh(x, method = 'GHQ', model = 'Mh', npoints = 50)$par[2:3]
}

mean((param.la2[,2]- N)/N)
mean((param.la4[,2]- N)/N)
mean((param.ghq[,2]- N)/N)
mean((exp(param.la2[,1])- sigma)/sigma)
mean((exp(param.la4[,1])- sigma)/sigma)
mean((exp(param.ghq[,1])- sigma)/sigma)
sapply(upper.la2-lower.la2, mean)
sapply(upper.la4-lower.la4, mean)
sapply(upper.ghq-lower.ghq, mean)

cp = 0
for (i in 1:iter) {
  if (N >= lower.la4[i,2] && N <= upper.la4[i,2]){
    cp = cp + 1
  } else {
    cp = cp + 0
  }
} 
cp
cp = 0
for (i in 1:iter) {
  if (sigma >= exp(lower.la4[i,1]) && sigma <= exp(upper.la4[i,1])){
    cp = cp + 1
  } else {
    cp = cp + 0
  }
} 
cp

#computational comparisons to execute the likelihood of Mh model from the last simulation
microbenchmark('LA2' = {
  mtbh(x, method = 'LA2', model = 'Mh', niter = 5)
},'LA4' = {
  mtbh(x, method = 'LA4', model = 'Mh', niter = 5)
},'GHQ' = {
  mtbh(x, method = 'GHQ', model = 'Mh', npoints = 50)
}, times = 100)

#CJS with continuous covariates
compile(paste0('CJSc_HMM', '.cpp'))
dyn.load(paste0('CJSc_HMM'))
compile(paste0('CJSc_LA', '.cpp'))
dyn.load(paste0('CJSc_LA'))
compile(paste0('CJSc_MI', '.cpp'))
dyn.load(paste0('CJSc_MI'))

source("CJSc.R")
xvole <- read.table("volecap.dat", quote="\"", comment.char="")
yvole <- read.table("voleweight.dat", quote="\"", comment.char="")


#since we model the capture probability as a function of covariates (covp=TRUE), 
#we can ignore timep
#timecov indicates the number of parameters assigned to covariate means, 
#e.g. 0 indicating mu_1 since C++ starts indexing from 0
#number of mus is T-1 which 3
#fitting the dataset via 20-HMM
CJSc(xvole, yvole, method = "HMM", covp=TRUE, timecov = seq(0,2,1), m=20) 

#fitting the dataset via the Laplace's method
CJSc(xvole, yvole, method = "LA2", covp=TRUE, timecov = seq(0,2,1)) 

#fitting the dataset using model 3 via multiple imputation
#first we need to create a function to estimate the imputation model for covariates
ll <- function(param, y, fi, timecov){
  k = length(param);  y <- as.matrix(y)
  n = dim(y)[1]; t = dim(y)[2]
  mut = param[-k]; sigma = exp(param[k])
  mu = rep(0,t-1)
  for (j in 1:t-1) {
    mu[j] = mut[timecov[j]]
  }
  nll = 0.0
  for (i in 1:n) {
    for (j in (fi[i]+1):t) {
      if (y[i,j] > 0 && y[i,j-1] > 0){
        nll = nll - dnorm(y[i,j], y[i,j-1] + mu[j-1], sigma, TRUE)
      } else {
        nll = nll - 0
      }
    }
  }
  return(nll)
}
#we consider 20 repeatition 
par.mi = data.frame(beta0=rep(0,20), beta1=rep(0,20), gamma0=rep(0,20), gamma1=rep(0,20))
n <- dim(yvole)[1]; T<-dim(yvole)[2]
timecov <- c(1,2,3) #since R starts from 1 for indexing
param <- rep(0, max(timecov)+1)
fi = rep(0,n); la = rep(0,n)
for (i in 1:n) {
  for (j in 1:t) {
    fi[i] = min(which(xvole[i,]==1))
    la[i] = max(which(xvole[i,]==1))
  }
}
#find the estimates of imputation model parameters
par = optim(param, ll, y=yvole, fi=fi, timecov=timecov)$par
mu = par[-4]; sm = exp(par[4])
#imputing covariates using the estmated parameters obtained above
for (q in 1:20) {
  wnew = as.matrix(yvole)
  for (i in 1:n) {
    for (j in (fi[i]+1):T) {
      if (wnew[i,j] == 0 && j < la[i]) {
        k = min(which(wnew[i,j:T] > 0))-1
        m = (k*(wnew[i,j-1] + mu[j-1] + wnew[i,j+k]) - sum(mu[(j):(j+k-1)]))/(k+1) 
        s = sqrt(k*sm^2/(k+1))
        wnew[i,j] = rnorm(1, mean = m, sd = s) 
      } else if (wnew[i,j] == 0 && j > la[i]) {
        wnew[i,j] = wnew[i,j-1] + mu[j-1] + sm*rnorm(1)
      } else {
        wnew[i,j] = wnew[i,j]
      }
    }
  }
  par.mi[q,] <- CJSc(xvole, wnew, method = "MI", covp=TRUE, timecov = seq(0,2,1))$par
}
sapply(par.mi, mean)

#use a nonparametric Bootstrap to obtain the confidence intervals
B=1000
param.mi = data.frame(beta0=rep(0,B), beta1=rep(0,B), gamma0=rep(0,B), gamma1=rep(0,B))
for (r in 1:B) {
  sampel = sample(1:n,replace=TRUE)
  xb = xvole[sampel, ]; yb = yvole[sampel, ]
  for (i in 1:n) {
    for (j in 1:t) {
      fi[i] = min(which(xb[i,]==1))
      la[i] = max(which(xb[i,]==1))
    }
  }
  par = optim(param, ll, y=yb, fi=fi, timecov=timecov)$par
  mu = par[-4]; sm = exp(par[4])
  wnew = yb
  for (q in 1:20) {
    for (i in 1:n) {
      for (j in (fi[i]+1):T) {
        if (wnew[i,j] == 0 && j < la[i]) {
          k = min(which(wnew[i,j:T] > 0))-1
          m = (k*(wnew[i,j-1] + mu[j-1] + wnew[i,j+k]) - sum(mu[(j):(j+k-1)]))/(k+1) 
          s = sqrt(k*sm^2/(k+1))
          wnew[i,j] = rnorm(1, mean = m, sd = s) 
        } else if (wnew[i,j] == 0 && j > la[i]) {
          wnew[i,j] = wnew[i,j-1] + mu[j-1] + sm*rnorm(1)
        } else {
          wnew[i,j] = wnew[i,j]
        }
      }
    }
    par.mi[q,] <- CJSc(xb, wnew, method = "MI", covp=TRUE, timecov = seq(0,2,1))$par 
  }
  param.mi[r,] <- sapply(par.mi, mean)
}
sapply(param.mi, quantile, 0.025)
sapply(param.mi, quantile, 0.975)

#simulation studies
#note here we exclude the multiple imputation since the only way to compute the confidence interval is via a nonparametric Bootstrap
#the computational time is extremely intensvie
#here we repeat several scenarios as explained in the draft
#set p=0.25, 0.5 and 0.75 for model 1 and 2 (covp = FALSE)
# timep = rep(0,9), because p is constant over time, if p is time-varying, then
#we set timep=timecov=seq(0,8,1)
source("sim_open.R")
nsim = 100
par_l1 = data.frame(beta0 = rep(0,nsim), beta1=rep(0,nsim), se_b0=rep(0,nsim), se_b1=rep(0,nsim))
par_l2 = data.frame(beta0 = rep(0,nsim), beta1=rep(0,nsim), se_b0=rep(0,nsim), se_b1=rep(0,nsim))
set.seed(3456)
parameters <- list(beta = c(-3, 0.2), mu = NULL, p=0.25, sigma=1.2, gamma=NULL)
for (i in 1:100) {
  sim <- sim_open(parameters = parameters, N=500, T=10, initial.cov = list(mu=15, sd=2))
  x <- sim$x 
  w <- sim$w
  hmm <- CJSc(x, w, method = "HMM", covp=FALSE, temporal=FALSE, timep = rep(0,9), m=20)
  par_l1[i,1:2] <- hmm$par[1:2]
  par_l1[i,3:4] <- hmm$se[1:2]
  
  la2 <- CJSc(x, w, method = "LA2", covp=FALSE, timecov = seq(0,8,1), timep = rep(0,9))
  par_l2[i,1:2] <- la2$par[1:2]
  par_l2[i,3:4] <- la2$se[1:2]
}

mean((par_l1[,1]-(-3.0))/(-3.0))
mean((par_l2[,1]-(-3.0))/(-3.0))

mean((par_l1[,2]-(0.2))/(0.2))
mean((par_l2[,2]-(0.2))/(0.2))

#for example we want to compute the covergae probabilites of the HMM
cp = 0
lower = par_l2[,2] - qnorm(0.975)*par_l2[,4]
upper = par_l2[,2] + qnorm(0.975)*par_l2[,4]
for (i in 1:nsim) {
  if (0.2 >= lower[i] && 0.2 <= upper[i]){
    cp = cp + 1
  } else {
    cp = cp + 0
  }
} 
cp
mean(upper-lower)
#computational comparisons
sim <- sim_open(parameters = parameters, N=500, T=10, initial.cov = list(mu=15, sd=2))
x <- sim$x 
w <- sim$w
microbenchmark('HMM' = {
    CJSc(x, w, method = "HMM", covp=FALSE, timecov = seq(0,8,1), timep = rep(0,9), m=20)
  }, 
  'LA2' = {
    CJSc(x, w, method = "LA2", covp=FALSE, timecov = seq(0,8,1), timep = rep(0,9))
}, 
times = 100)
