#create a function for general Mh-type models
mtbh <-function(x, method = 'LA2', model = 'Mh', npoints=50, maxit=1e+03, niter=5, 
                starting.value = NULL){
  x <- as.matrix(x); n <- dim(x)[1]; T <- dim(x)[2]
  ys = rowSums(x);  ys = ys[ys!=0]
  nk = rep(0,T)
  for(k in 1:T){
    nk[k] = length(ys[ys==k])
  }
  y = matrix(0, nrow=n, ncol=T)
  for (i in 1:n) {
    for (j in 1:(T-1)) {
      if (x[i,j] == 1){
        y[i,c((j+1):T)] = c(rep(1,T-j))
        y[i,j] = 0
        break
      } else if(x[i,T] == 1){
        y[i,T] = 0
      } else {
        y[i,j] = 0
      }
    }
  }
  alpha <- NULL; 
  w <- gauss.quad(npoints,"hermite")$weights
  node <- gauss.quad(npoints,"hermite")$nodes
  if(model == 'Mh'){
    y <- matrix(0, nrow = n, ncol = T)
    if(method == 'LA2'){
      data <- list(x=x, y=y, nk=nk, T=T, w=w, node=node, method=1, model=1, niter=niter)
      parameters <- list(alphat=rep(0,1), log_sigma=0, N=n)
      la.tmb <- MakeADFun(data, parameters, DLL='Mh_type', silent = TRUE)
    } 
    else if (method == 'LA4'){
      data <- list(x=x, y=y, nk=nk, T=T, w=w, node=node, method=2, model=1, niter=niter)
      parameters <- list(alphat=rep(0,1), log_sigma=0, N=n)
      la.tmb <- MakeADFun(data, parameters, DLL='Mh_type', silent = TRUE)
    }
    else if (method == 'GHQ'){
      data <- list(x=x, y=y, nk=nk, T=T, w=w, node=node, method=3, model=1, niter=niter)
      parameters <- list(alphat=rep(0,1), log_sigma=0, N=n)
      la.tmb <- MakeADFun(data, parameters, DLL='Mh_type', silent = TRUE)
    }
    else { print('Select the correct method')}
  }
  else if (model == 'Mth'){
    y <- matrix(0, nrow = n, ncol = T)
    if(method == 'LA2'){
      data <- list(x=x, y=y, nk=nk, T=T, w=w, node=node, method=1, model=2, niter=niter)
      parameters <- list(alphat=rep(0,T), log_sigma=0, N=n)
      la.tmb <- MakeADFun(data, parameters, DLL='Mh_type', silent = TRUE)
    } 
    else if (method == 'LA4'){
      data <- list(x=x, y=y, nk=nk, T=T, w=w, node=node, method=2, model=2, niter=niter)
      parameters <- list(alphat=rep(0,T), log_sigma=0, N=n)
      la.tmb <- MakeADFun(data, parameters, DLL='Mh_type', silent = TRUE)
    }
    else if (method == 'GHQ'){
      data <- list(x=x, y=y, nk=nk, T=T, w=w, node=node, method=3, model=2, niter=niter)
      parameters <- list(alphat=rep(0,T), log_sigma=0, N=n)
      la.tmb <- MakeADFun(data, parameters, DLL='Mh_type', silent = TRUE)
    }
    else { print('Select the correct method')}
  }
  else if (model == 'Mbh'){
    if(method == 'LA2'){
      data <- list(x=x, y=y, nk=nk, T=T, w=w, node=node, method=1, model=3, niter=niter)
      parameters <- list(alphat=rep(0,1), log_sigma=0, N=n, b=0)
      la.tmb <- MakeADFun(data, parameters, DLL='Mh_type', silent = TRUE)
    } 
    else if (method == 'LA4'){
      data <- list(x=x, y=y, nk=nk, T=T, w=w, node=node, method=2, model=3, niter=niter)
      parameters <- list(alphat=rep(0,1), log_sigma=0, N=n, b=0)
      la.tmb <- MakeADFun(data, parameters, DLL='Mh_type', silent = TRUE)
    }
    else if (method == 'GHQ'){
      data <- list(x=x, y=y, nk=nk, T=T, w=w, node=node, method=3, model=3, niter=niter)
      parameters <- list(alphat=rep(0,1), log_sigma=0, N=n, b=0)
      la.tmb <- MakeADFun(data, parameters, DLL='Mh_type', silent = TRUE)
    }
    else { print('Select the correct method')}
  }
  else if (model == 'Mtbh'){
    if(method == 'LA2'){
      data <- list(x=x, y=y, nk=nk, T=T, w=w, node=node, method=1, model=4, niter=niter)
      parameters <- list(alphat=rep(0,T), log_sigma=0, N=n, b=0)
      la.tmb <- MakeADFun(data, parameters, DLL='Mh_type', silent = TRUE)
    } 
    else if (method == 'LA4'){
      data <- list(x=x, y=y, nk=nk, T=T, w=w, node=node, method=2, model=4, niter=niter)
      parameters <- list(alphat=rep(0,T), log_sigma=0, N=n, b=0)
      la.tmb <- MakeADFun(data, parameters, DLL='Mh_type', silent = TRUE)
    }
    else if (method == 'GHQ'){
      data <- list(x=x, y=y, nk=nk, T=T, w=w, node=node, method=3, model=4, niter=niter)
      parameters <- list(alphat=rep(0,T), log_sigma=0, N=n, b=0)
      la.tmb <- MakeADFun(data, parameters, DLL='Mh_type', silent = TRUE)
    }
    else { print('Select the correct method')}
  }
  fit <- try(optim(la.tmb$par, la.tmb$fn, la.tmb$gr, la.tmb$he, method = 'BFGS', 
               control = list(reltol=1e-8, maxit=maxit), hessian = FALSE), silent = TRUE)
  k <- length(fit$par)
  AIC <- 2*k + 2*la.tmb$fn(fit$par)
  AICc <- AIC + ((2*k^2+2*k)/(n-k-1))
  LogL <- data.frame(LogL = -la.tmb$fn(fit$par), AIC, AICc)
  
  output <- list()
  output$par <- fit$par
  output$LogL <- LogL
  output$message <- fit$convergence
  return(output)
}