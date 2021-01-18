#create a function for general Mh-type models
mtbh <-function(x, method = 'LA2', model = 'Mh', npoints=50, maxit=1000, niter=5){
  x <- as.matrix(x); n <- dim(x)[1]; T <- dim(x)[2]
  ys = rowSums(x);  ys = ys[ys!=0]
  nk = rep(0,T)
  for(k in 1:T){
    nk[k] = length(ys[ys==k])
  }
  y = matrix(0, nrow=n, ncol=T)
  for (i in 1:n) {
    fi = min(which(x[i,]==1))
    if (fi == T) {
      y[i,T] = 0
    } else {
      for (j in fi:(T-1)) {
        y[i,(j+1)] = 1
      }
    }
  }
  
  w <- gauss.quad(npoints,"hermite")$weights
  node <- gauss.quad(npoints,"hermite")$nodes
  data <- list()
  parameters <- list()
  if(model == 'Mh' | model == 'Mbh'){
    if (model == "Mh") {model = 1}
    else {model = 3}
    parameters <- list(alphat=rep(0,1), log_sigma=0, N=n, b=0)
    if(method == 'LA2'){data <- list(x=x, y=y, nk=nk, T=T, w=w, node=node, method=1, model=model, niter=niter)} 
    else if (method == 'LA4'){data <- list(x=x, y=y, nk=nk, T=T, w=w, node=node, method=2, model=model, niter=niter)}
    else if (method == 'GHQ'){data <- list(x=x, y=y, nk=nk, T=T, w=w, node=node, method=3, model=model, niter=niter)}
    else { print('Select the correct method')}
  }
  else if (model == 'Mth' | model == "Mtbh"){
    if (model == "Mth") {model = 2}
    else {model = 4}
    parameters <- list(alphat=rep(0,T), log_sigma=0, N=n, b=0)
    if(method == 'LA2'){data <- list(x=x, y=y, nk=nk, T=T, w=w, node=node, method=1, model=model, niter=niter)} 
    else if (method == 'LA4'){data <- list(x=x, y=y, nk=nk, T=T, w=w, node=node, method=2, model=model, niter=niter)}
    else if (method == 'GHQ'){data <- list(x=x, y=y, nk=nk, T=T, w=w, node=node, method=3, model=model, niter=niter)}
    else { print('Select the correct method')}
  } else {print('Select the correct model')}
  
  fittmb <- MakeADFun(data=data, parameters=parameters, DLL='Mh_type', silent = TRUE)
  fit <- try(optim(fittmb$par, fittmb$fn, fittmb$gr, fittmb$he, method = 'BFGS', 
                   control = list(reltol=1e-8, maxit=maxit), hessian = FALSE), silent = TRUE)
  k <- length(fit$par)
  AIC <- 2*k + 2*fittmb$fn(fit$par)
  AICc <- AIC + ((2*k^2+2*k)/(n-k-1))
  LogL <- data.frame(LogL = -fittmb$fn(fit$par), AIC, AICc)
  
  output <- list()
  output$par <- fit$par
  output$LogL <- LogL
  output$message <- fit$convergence
  return(output)
}