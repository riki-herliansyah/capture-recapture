#include <TMB.hpp>
#include<math.h>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(x);
  DATA_MATRIX(y);
  DATA_IVECTOR(la);
  DATA_IVECTOR(fi);
  DATA_IVECTOR(timep);
  DATA_INTEGER(model);

  PARAMETER_VECTOR(beta);

  int n = x.rows();
  int T = x.cols();

  matrix<Type> eta(n,T-1);
  matrix<Type> phi(n,T-1);
   for (int i = 0; i < n; i++){
        for (int j = 0; j < T-1; j++){
            eta(i,j) = beta[0] + beta[1]*y(i,j);
            phi(i,j) = invlogit(eta(i,j));
        }
  }
  Type nll = 0.0;
  vector<Type> chi(T);
  matrix<Type> mu(n,T-1);
  matrix<Type> p(n,T-1);

  if (model == 2 || model == 1){
  PARAMETER_VECTOR(alpha);
  vector<Type> alphat(T-1);
  for (int j = 0; j < T-1; j++){
    alphat[j] = alpha(timep[j]);
  }
  for (int i = 0; i < n; i++){
        for (int j = 0; j < T-1; j++){
            mu(i,j) = alphat[j];
            p(i,j) = invlogit(mu(i,j));
        }
  }
  }

  else if (model == 4 || model == 3){
  PARAMETER_VECTOR(gamma);
  for (int i = 0; i < n; i++){
        for (int j = 0; j < T-1; j++){
            mu(i,j) = gamma[0] + gamma[1]*y(i,j+1);
            p(i,j) = invlogit(mu(i,j));
        }
  }
  }
    for (int i = 0; i < n; i++){
        for (int j = 1; j < T; j++){
            chi[T-1] = 1;
            chi[T-j-1] = (1-phi(i,T-j-1)) + phi(i,T-j-1)*(1-p(i,T-j-1))*chi[T-j];
        }
        if(fi[i] == la[i]){
            nll -= log(chi(la[i]-1));
        } else {
            Type logp = 0.0;
            for (int j = fi[i]; j < la[i]; j++){
                logp += x(i, j)*log(p(i,j-1)) + (1-x(i,j))*log(1-p(i,j-1));
            }
            Type logph = 0.0;
            for (int j = (fi[i]-1); j < (la[i]-1); j++){
                logph += log(phi(i,j));
            }
            nll -= logph + logp + log(chi(la[i]-1));
        }
    }
    return nll;
}
