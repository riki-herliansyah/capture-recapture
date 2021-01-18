sim_closed <- function(parameter, T, behaviour= FALSE){
  alpha = parameter$alpha; N = parameter$N; sigma = parameter$sigma
  p = matrix(0, ncol = T, nrow = N)
  y = matrix(0, ncol = T+1, nrow = N)
  z = sigma*rnorm(N)
  for (i in 1:N){
    for (j in 1:T) {
        p[i,j] = inv.logit(alpha[j] + z[i])
        y[i,j] = rbinom(1, size=1, prob=p[i,j])
    }
    if (behaviour == TRUE){
      beta = parameter$beta
      if (sum(y[i,c(1:T-1)]) > 0){
        fi = min(which(y[i,c(1:T-1)]==1))
        for (j in (fi+1):T) {
          p[i,j] = inv.logit(alpha[j] + beta + z[i])
          y[i,j] = rbinom(1, size=1, prob=p[i,j])
        }
      }
    }
    y[i, T+1] = sum(y[i,])
  }
  x <- y[y[,(T+1)]!=0,-(T+1)]
  return(x)
}
  