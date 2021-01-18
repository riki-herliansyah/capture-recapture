sim_open <- function(parameters, N, T, initial.cov=list(mu=15, sd=2)){
  x<-matrix(NA,N,T)
  w<-matrix(NA,N,T)
  alive<-matrix(NA,N,T)
  mu <- parameters$mu #the constant of the covariate process 
  beta<- parameters$beta #regression parameters of the survival probability
  gam <- parameters$gamma #regression parameters of the capture probability, assuming p is specified as a function of covariates
  p <- parameters$p
  ifelse(length(p)==1, p<-rep(p,T), p<-p) 
  sigma<- parameters$sigma        
  mu.in <- initial.cov$mu; sd.in <- initial.cov$sd
  
  ## simulate survival, covariate and observation process
  for (k in 1:N){
    fi<-sample(1:(T-1), size=1)
    alive[k,fi]<-1
    x[k,fi]<-1
    w[k,fi]<-rnorm(1, mu.in, sd.in) # initial body mass distribution
    for (i in (fi+1):T){
      w[k,i]<- w[k,i-1]+mu[i-1]+sigma*rnorm(1) 
      ifelse(alive[k,i-1]==1,alive[k,i]<-rbinom(1,size=1,prob=inv.logit(beta[1]+beta[2]*w[k,i-1])),alive[k,i]<-0)
      if (alive[k,i]==1) {
        ifelse(is.null(gam)==TRUE, x[k,i]<-rbinom(1,size=1,prob=p[i]), 
               x[k,i]<-rbinom(1,size=1, prob=inv.logit(gam[1]+gam[2]*w[k,i])))
      } 
      else {
        ifelse(alive[k,i-1]==0,x[k,i]<-0,x[k,i]<-0)
      }
    }
  }
  
  ## set missing covariate values equal to NA 
  for (i in 1:N){
    for (j in 1:T){
      if (!is.na(x[i,j]) & x[i,j]==0) w[i,j]<-NA
      if (is.na(x[i,j]))              w[i,j]<-NA
    }
  }
  for (i in 1:N) {
    for (j in 1:T) {
      if(is.na(x[i,j]) == TRUE){
        x[i,j] = 0.0
      } else {
        x[i,j] = x[i,j]
      }
    }
  }
  output <- list()
  output$x <- x
  output$w <- w
  return(output)
}

