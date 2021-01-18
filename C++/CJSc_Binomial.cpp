#include <TMB.hpp>
#include<math.h>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(x);
  DATA_MATRIX(y);
  DATA_IVECTOR(timep);

  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(alpha);

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
  vector<Type> alphat(T-1);
  for (int j = 0; j < T-1; j++){
    alphat[j] = alpha(timep[j]);
  }
  vector<Type> p = invlogit(alphat);
  for (int i = 0; i < n; i++){
       for (int j = 0; j < T-1; j++){
            if (x(i,j+1) == 1) {
                if (x(i,j) == 1){
                    nll -= log(phi(i,j)) + log(p[j]);
                } else if (x(i,j) == 0) {
                    nll -= 0.0;
                }
            } else if (x(i,j+1) == 0) {
               if (x(i,j) == 1){
                    nll -= log(1 - (phi(i,j)*p[j]));
                } else if (x(i,j) == 0) {
                    nll -= 0.0;
                }
            }
       }
  }
    return nll;
}
