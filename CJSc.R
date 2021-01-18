CJSc <- function(x, y, method = "LA2", covp=TRUE, temporal=TRUE, m=20, timep=NULL, timecov=NULL){
  x = as.matrix(x); y = as.matrix(y)
  T = dim(y)[2]; n = dim(y)[1]
  fi = rep(0,n); la = rep(0,n)
  for (i in 1:n) {
    fi[i] = min(which(x[i,]==1))
    la[i] = max(which(x[i,]==1))
  }
  #assign 0 if the covariate is missing
  #assuming no zero value in covariate
  for (i in 1:n) {
    for (j in 1:T) {
      if(is.na(y[i,j]) == TRUE){
        y[i,j] = 0.0
      } else {
        y[i,j] = y[i,j]
      }
    }
  }
  #counting the total number of missing covariates need to be imputed for Laplace approximation
  nk = 0.0
  for (i in 1:n) {
    for (j in fi[i]:T) {
      if (fi[i] == T){
        y[i,j] = y[i,j]
      } else {
        if(y[i,j] == 0){
          nk = nk + 1
        } 
      }
    }
  }

  #fitting the dataset 3 via the HMM
  #find the lower an dthe upper limit of covariate values
  mx = max(c(y)[c(y)>0])
  mn = min(c(y)[c(y)>0])
  #for this dataset, m=20 is adequate

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
  if (covp == TRUE){
    ifelse (temporal == TRUE, model<-3, model<-4)
    parmu <- max(timecov)+1
    timep <- 0
    data = list(x=x, y=y, fi=fi, la=la, timep=timep, timecov=timecov, model=model, Bj=Bj, Bs=Bs)
    if (method == "HMM"){
      parameters = list(beta=c(0,0), gamma=c(0,0), mu= rep(0,parmu), log_sigma=0)
      fittmb <- MakeADFun(data=data, parameters=parameters, DLL="CJSc_HMM", silent = TRUE)
      
    } else if (method == "LA2"){
      parameters = list(beta=c(0,0), gamma=c(0,0), mu = rep(0,parmu), log_sigma=0, u=rep(0,nk))
      fittmb <- MakeADFun(data, parameters, DLL="CJSc_LA", random='u', silent = TRUE)
    } else if (method == "MI"){
      parameters = list(gamma=c(0,0), beta=c(0,0))
      fittmb <- MakeADFun(data=data, parameters=parameters, DLL='CJSc_MI', silent = TRUE)
    }
  } else {
    ifelse(temporal == TRUE, model<-1, model<-2)
    ifelse(temporal == TRUE, timecov<-timecov, timecov<-0)
    parmu <- max(timecov)+1
    parp <- max(timep)+1
    data = list(x=x, y=y, fi=fi, la=la, timep=timep, timecov=timecov, model=model, Bj=Bj, Bs=Bs)
    if (method == "HMM"){
      parameters = list(beta=c(0,0), gamma=c(0,0), alpha=rep(0,parp), mu=rep(0,parmu), log_sigma=0)
      fittmb <- MakeADFun(data=data, parameters=parameters, DLL="CJSc_HMM", silent = TRUE)
      
    } else if (method == "LA2"){
      parameters = list(beta=c(0,0), gamma=c(0,0), alpha=rep(0,parp), mu=rep(0,parmu), log_sigma=0, u=rep(0,nk))
      fittmb <- MakeADFun(data=data, parameters=parameters, DLL="CJSc_LA", random='u', silent = TRUE)
    } else if (method == "MI"){
      parameters = list(gamma=c(0,0), alpha=rep(0,parp))
      fittmb <- MakeADFun(data=data, parameters=parameters, DLL='CJSc_MI', silent = TRUE)
    }
  }
  if (method == "binomial"){
    parp <- max(timep)+1
    data = lis(x=x, y=y, timep=timep)
    parameters = list(beta=c(0,0), alpha = rep(0,parp))
    fittmb <- MakeADFun(data=data, parameters=parameters, DLL="CJSc_Binomial", silent = TRUE)
  }
  #we consider to use the fit_tmb for optimization from TMBhelper
  try(fit_tmb(fittmb), silent = TRUE)
  se <- try(sdreport(fittmb), silent = TRUE)
  param <- se$par.fixed
  separam <- try(sqrt(diag(se$cov.fixed)), silent = TRUE)
  k <- length(param)
  AIC <- 2*k + 2*fittmb$fn(param)
  AICc <- AIC + ((2*k^2+2*k)/(n-k-1))
  LogL <- data.frame(LogL = -fittmb$fn(param), AIC, AICc)
  
  output <- list()
  output$par <- param
  output$se <- separam
  output$LogL <- LogL
  return(output)
}
