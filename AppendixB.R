require(TMB)
require(TMBhelper)
require(statmod)
require(boot)
require(microbenchmark)
path <- "C:/Users/s1898267/OneDrive - University of Edinburgh/PhD/First Year/Annual Review/Appendix codes"
setwd(path)

#Mh Type models
compile(paste0('Mh_type', '.cpp'))
dyn.load(paste0('Mh_type'))
source("mtbh.R")

#preparing data 
x <- read.table("golftees.txt", quote="\"", comment.char="")
x <- as.matrix(x) 
#fitting Mh model using Laplace approximations and GHQ
#model consists of Mh, Mth, Mbh and Mtbh
mtbh(x, method = 'LA2', model = 'Mh', niter = 20)
mtbh(x, method = 'LA4', model = 'Mh', niter = 5)
mtbh(x, method = 'GHQ', model = 'Mbh', niter = 5)

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
iter = 100
lower.la2 = data.frame(log_sigma=rep(0,iter), N=rep(0,iter)) 
upper.la2 = data.frame(log_sigma=rep(0,iter), N=rep(0,iter)) 
param.la2 = data.frame(log_sigma=rep(0,iter), N=rep(0,iter)) 
lower.la4 = data.frame(log_sigma=rep(0,iter), N=rep(0,iter)) 
upper.la4 = data.frame(log_sigma=rep(0,iter), N=rep(0,iter)) 
param.la4 = data.frame(log_sigma=rep(0,iter), N=rep(0,iter)) 
lower.ghq = data.frame(log_sigma=rep(0,iter), N=rep(0,iter)) 
upper.ghq = data.frame(log_sigma=rep(0,iter), N=rep(0,iter)) 
param.ghq = data.frame(log_sigma=rep(0,iter), N=rep(0,iter))

set.seed(2345)
parameter = list (alpha=rep(-1.0, 8), N=250, sigma=0.75, beta=0.75)
for (k in 1:iter) {
  x <- sim_closed(parameter = parameter, T=8, behaviour=TRUE)
  n <- dim(x)[1]
  B=1000
  pN1 = data.frame(log_sigma=rep(0,B), N=rep(0,B)) 
  pN2 = data.frame(log_sigma=rep(0,B), N=rep(0,B)) 
  pN3 = data.frame(log_sigma=rep(0,B), N=rep(0,B)) 
  for (p in 1:B) {
    sampel = sample(1:n, replace=TRUE)
    xb = x[sampel, ]
    pN1[p,] =  mtbh(xb, method = 'LA2', model = 'Mbh', niter = 5)$par[2:3]
    pN2[p,] =  mtbh(xb, method = 'LA4', model = 'Mbh', niter = 5)$par[2:3]
    pN3[p,] =  mtbh(xb, method = 'GHQ', model = 'Mbh', npoints = 50)$par[2:3]
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
for (k in 1:iter) {
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
  mtbh(x, method = 'LA2', model = 'Mh', npoints = 50)
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
xv = as.matrix(xvole)
yv = as.matrix(yvole)
t = dim(yv)[2]; n = dim(yv)[1]
fi = rep(0,n); l = rep(0,n)
#assign 0 if the capture history is missing
#assuming no zero value in covariate
for (i in 1:n) {
  for (j in 1:t) {
    if(is.na(xv[i,j]) == TRUE){
      xv[i,j] = 0.0
    } else {
      xv[i,j] = xv[i,j]
    }
  }
}

for (i in 1:n) {
  for (j in 1:t) {
    fi[i] = min(which(xv[i,]==1))
    l[i] = max(which(xv[i,]==1))
  }
}
#assign 0 if the covariate is missing
for (i in 1:n) {
  for (j in 1:t) {
    if(is.na(yv[i,j]) == TRUE){
      yv[i,j] = 0.0
    } else {
      yv[i,j] = yv[i,j]
    }
  }
}
#counting the total number of missing covariates need to be imputed for Laplace approximation
nk = 0.0
for (i in 1:n) {
  for (j in fi[i]:t) {
    if (fi[i] == t){
      yv[i,j] = yv[i,j]
    } else {
      if(yv[i,j] == 0){
        nk = nk + 1
      } 
    }
  }
}

#fitting the dataset 3 via the HMM
#find the lower an dthe upper limit of covariate values
mx = max(c(yv)[c(yv)>0])
mn = min(c(yv)[c(yv)>0])
#for this dataset, m=20 is adequate
m=20
max.w  <-  round(12/10*mx) # essential range upper limit
min.w  <-  round(8/10*mn)  # essential range lower limit
K      <-m+1
#discretizing into m spaces
Bj     <- seq(min.w, max.w,length=K)
#finding the mid points
Bs  <-(Bj[-1]+Bj[-K])*0.5
#prepapring inputs
#since we model the capture probability as a function of covariates, we can ignore time_a
#time_mu indicates the number of parameters assigned to covariate means, e.g. 0 indicating mu_1
#C++ starts indexing from 0

data = list(x=xv, y=yv, Bj=Bj, Bs=Bs, fi=fi, la=l,timep=c(0,1,2), timecov=c(0,1,2), model=3)
parameters = list(beta=c(0,0), gamma=c(0,0), mu= c(0,0,0), log_sigma=0)
hm <- MakeADFun(data, parameters, DLL='CJSc_HMM', silent = TRUE)
#we consider to use the fit_tmb for optimization from TMBhelper
fit_tmb(hm)
se<-sdreport(hm)

#fitting the dataset via the Laplace's method
data = list(x=xv, y=yv, e=fi, l=l, time_a=c(0,1,2), time_mu=c(0,1,2), model=3)
parameters = list(beta=c(0,0), gamma=c(0,0), mu = rep(0,3), log_sigma=0, u=rep(0,nk))
la <- MakeADFun(data, parameters, DLL='CJSc_LA', random='u', silent = TRUE)
fit_tmb(la)
sdreport(la)


#fitting the dataset using model 3 via multiple imputation
#first we need to create a function to estimate the imputation model for covariates
ll <- function(param, w, e, tm){
  k = length(param); n = dim(w)[1]
  t = dim(w)[2]
  mut = param[-k]; sigma = exp(param[k])
  mu = rep(0,t-1)
  for (j in 1:t-1) {
    mu[j] = mut[tm[j]]
  }
  nll = 0.0
  for (i in 1:n) {
    for (j in (e[i]+1):t) {
      if (w[i,j] > 0 && w[i,j-1] > 0){
        nll = nll - dnorm(w[i,j], w[i,j-1] + mu[j-1], sigma, TRUE)
      } else {
        nll = nll - 0
      }
    }
  }
  return(nll)
}
#we consider 20 repeatition 
pl4 = data.frame(beta0=rep(0,20), beta1=rep(0,20), gamma0=rep(0,20), gamm1=rep(0,20))
param =  rep(0,4)
#find the estimates of imputation model parameters
par = optim(param, ll, w=yv, e=fi, tm=c(1,2,3))$par
mut = par[-4]; sm = exp(par[4]); tm=c(1,2,3)
mu = rep(0,t-1)
for (j in 1:t-1) {
  mu[j] = mut[tm[j]]
}
#imputing covariates using the estmated parameters obtained above
for (q in 1:m) {
    wnew = yv
    for (i in 1:n) {
      for (j in (fi[i]+1):t) {
        if (wnew[i,j] == 0 && j < l[i]) {
          k = min(which(wnew[i,j:t] > 0))-1
          m = (k*(wnew[i,j-1] + mu[j-1] + wnew[i,j+k]) - sum(mu[(j):(j+k-1)]))/(k+1) 
          s = sqrt(k*sm^2/(k+1))
          wnew[i,j] = rnorm(1, mean = m, sd = s) 
        } else if (wnew[i,j] == 0 && j > l[i]) {
          wnew[i,j] = wnew[i,j-1] + mu[j-1] + sm*rnorm(1)
        } else {
          wnew[i,j] = wnew[i,j]
        }
      }
    }
    data = list(x=xv, y=wnew, l=l, e=fi, time_a=tm, model=2)
    parameters = list(gamma=c(0,0), beta=c(0,0))
    l4 <- MakeADFun(data, parameters, DLL='CJSc_MI', silent = TRUE)
    fit <- optim(l4$par, l4$fn, l4$gr, l4$he, method = 'BFGS', control=list(maxit=100))
    pl4[q,] <- fit$par 
}
sapply(pl4, mean)

#use a nonparametric Bootstrap to obtain the confidence intervals
B=1000
param.mi = data.frame(beta0=rep(0,B), beta1=rep(0,B), gamma0=rep(0,B), gamm1=rep(0,B))
for (r in 1:B) {
  sampel = sample(1:n,replace=TRUE)
  xb = xv[sampel, ]; yb = yv[sampel, ]
  for (i in 1:n) {
    for (j in 1:t) {
      e[i] = min(which(xb[i,]==1))
      l[i] = max(which(xb[i,]==1))
    }
  }
  par = optim(param, ll, w=yb, e=e, tm=c(1,2,3))$par
  mut = par[-4]; sm = exp(par[4])
  mu = rep(0,t-1)
  for (j in 1:t-1) {
    mu[j] = mut[tm[j]]
  }
  wnew = yb
  for (q in 1:20) {
    for (i in 1:n) {
      for (j in (e[i]+1):t) {
        if (wnew[i,j] == 0 && j < l[i]) {
          k = min(which(wnew[i,j:t] > 0))-1
          m = (k*(wnew[i,j-1] + mu[j-1] + wnew[i,j+k]) - sum(mu[(j):(j+k-1)]))/(k+1) 
          s = sqrt(k*sm^2/(k+1))
          wnew[i,j] = rnorm(1, mean = m, sd = s) 
        } else if (wnew[i,j] == 0 && j > l[i]) {
          wnew[i,j] = wnew[i,j-1] + mu[j-1] + sm*rnorm(1)
        } else {
          wnew[i,j] = wnew[i,j]
        }
      }
    }
    data = list(x=xb, y=wnew, l=l, e=e)
    parameters = list(gamma=c(0,0), beta=c(0,0))
    l4 <- MakeADFun(data, parameters, DLL='CJSc_MI', silent = TRUE)
    fit <- optim(l4$par, l4$fn, l4$gr, l4$he, method = 'BFGS', control=list(maxit=1000))
    pl4[q,] <- fit$par 
  }
  param.mi[r,] <- sapply(pl4, mean)
}
sapply(param.mi, quantile, 0.025)
sapply(param.mi, quantile, 0.975)

#simulation studies
#note here we exclude the multiple imputation since the only way to compute the confidence interval is via a nonparametric Bootstrap
#the computational time is huge and unreasonable
nsim = 100
par_l1 = data.frame(beta0 = rep(0,nsim), beta1=rep(0,nsim), se_b0=rep(0,nsim), se_b1=rep(0,nsim))
par_l2 = data.frame(beta0 = rep(0,nsim), beta1=rep(0,nsim), se_b0=rep(0,nsim), se_b1=rep(0,nsim))

set.seed(2345)
source("sim_open.R")
parameters <- list(beta = c(-3, 0.2), mu = rnorm(9), p=0.25, sigma=1.2, gamma=NULL)
sim <- sim_open(parameters = parameters, 500, 10, initial.cov = list(mu=15,sd=2))
CJSc(sim$x, sim$w, covp=FALSE, timecov = seq(0,8,1), timep = rep(0,9))
for (r in 1:nsim) {
  N=500
  T=10
  x<-matrix(NA,N,T)
  w<-matrix(NA,N,T)
  alive<-matrix(NA,N,T)
  mu <- rnorm(9) #the time-invariant constant of the covariate process 
  beta<-c(-3.0, 0.2)   
  p <- 0.25
  sigma<-1.2        
  
  ## simulate survival, covariate and observation process
  for (k in 1:N){
    fi<-sample(1:(T-1), size=1)
    alive[k,fi]<-1
    x[k,fi]<-1
    w[k,fi]<-rnorm(1,15,2) # initial body mass distribution
    for (i in (fi+1):T){
      w[k,i]<- w[k,i-1]+mu[i-1]+sigma*rnorm(1) 
      ifelse(alive[k,i-1]==1,alive[k,i]<-rbinom(1,size=1,prob=inv.logit(beta[1]+beta[2]*w[k,i-1])),alive[k,i]<-0)
      if (alive[k,i]==1) {x[k,i]<-rbinom(1,size=1,prob=p)} else {ifelse(alive[k,i-1]==0,x[k,i]<-0,x[k,i]<-0)}
    }
  }
  
  ## set missing covariate values equal to NA 
  for (i in 1:N){
    for (j in 1:T){
      if (!is.na(x[i,j]) & x[i,j]==0) w[i,j]<-NA
      if (is.na(x[i,j]))              w[i,j]<-NA
    }}
  
  t = dim(w)[2]; n = dim(w)[1]
  fi = rep(0,n); l = rep(0,n)
  
  for (i in 1:n) {
    for (j in 1:t) {
      fi[i] = min(which(x[i,]==1))
      l[i] = max(which(x[i,]==1))
    }
  }
  #set missing value back to 0
  for (i in 1:n) {
    for (j in 1:t) {
      if(is.na(w[i,j]) == TRUE){
        w[i,j] = 0.0
      } else {
        w[i,j] = w[i,j]
      }
    }
  }
  
  nk = 0.0
  for (i in 1:n) {
    for (j in fi[i]:t) {
      if (fi[i] == T){
        w[i,j] = w[i,j]
      } else {
        if(w[i,j] == 0){
          nk = nk + 1
        } 
      }
    }
  }
  
  time_mu = seq(0,8,1)
  time_a = rep(0,9)
  
  #HMM
  mx = max(c(w)[c(w)>0])
  mn = min(c(w)[c(w)>0])
  L=30
  max.w  <-  round(12/10*mx) # essential range upper limit
  min.w  <-  round(8/10*mn)  # essential range lower limit
  K      <-L+1
  Bj     <- seq(min.w, max.w,length=K)
  Bs  <-(Bj[-1]+Bj[-K])*0.5
  data = list(x=x, y=w, Bj=Bj, Bs=Bs, fi=fi, la=l, timep=time_a, timecov=time_mu, model=1)
  parameters = list(beta=c(0,0), alpha=c(0), mu=rep(0,9), log_sigma=0)
  l1 <- MakeADFun(data, parameters, DLL='CJSc_HMM', silent = TRUE)
  try(fit_tmb(l1), silent = TRUE)
  se <- try(sdreport(l1), silent = TRUE)
  par_l1[r,1:2] <- se$par.fixed[1:2]
  par_l1[r,3:4] <- try(sqrt(diag(se$cov.fixed)), silent = T)[1:2]
  
  #Laplace
  data = list(x=x, y=w, fi=fi, la=l, timep=time_a, timecov=time_mu, model=1)
  parameters = list(beta=c(0,0), alpha=c(0), mu=rep(0,9), log_sigma=0, u=rep(0,nk))
  l2 <- MakeADFun(data, parameters, DLL='CJSc_LA', random='u', silent = TRUE)
  try(fit_tmb(l2),silent = TRUE)
  se <- try(sdreport(l2), silent = TRUE)
  par_l2[r,1:2] <- se$par.fixed[1:2]
  par_l2[r,3:4] <- try(sqrt(diag(se$cov.fixed)), silent = T)[1:2]
  
}

mean((par_l1[,1]-(-3.0))/(-3.0))
mean((par_l2[,1]-(-3.0))/(-3.0))

mean((par_l1[,2]-(0.2))/(0.2))
mean((par_l2[,2]-(0.2))/(0.2))

#for example we want to compute the covergae probabilites of the HMM
cp = 0
lower = par_l1[,1] - qnorm(0.975)*par_l1[,3]
upper = par_l1[,3] + qnorm(0.975)*par_l1[,3]
for (i in 1:nsim) {
  if (-3 >= lower[i] && -3 <= upper[i]){
    cp = cp + 1
  } else {
    cp = cp + 0
  }
} 
cp
mean(upper-lower)
#computational comparisons
microbenchmark('HMM' = {
  data = list(x=x, y=w, Bj=Bj, Bs=Bs, e=fi, la=l, time_a=time_a, time_mu=time_mu, model=1)
  parameters = list(beta=c(0,0), alpha=c(0), mu=rep(0,9), log_sigma=0)
  l1 <- MakeADFun(data, parameters, DLL='CJSc_HMM', silent = TRUE)
  try(fit_tmb(l1), silent = TRUE)
  se <- try(sdreport(l1), silent = TRUE)
},
'LA2' = {
  data = list(x=x, y=w, e=fi, l=l, time_a=time_a, time_mu=time_mu, model=1)
  parameters = list(beta=c(0,0), alpha=c(0), mu=rep(0,9), log_sigma=0, u=rep(0,nk))
  l2 <- MakeADFun(data, parameters, DLL='CJSc_LA', random='u', silent = TRUE)
  try(fit_tmb(l2),silent = TRUE)
  se <- try(sdreport(l2), silent = TRUE)
}, 
times = 100)
