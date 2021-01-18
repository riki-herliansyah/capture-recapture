#include <TMB.hpp>
#include<math.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(x);
  DATA_MATRIX(y);
  DATA_VECTOR(Bj);
  DATA_VECTOR(Bs);
  DATA_IVECTOR(fi);
  DATA_IVECTOR(la);
  DATA_IVECTOR(timep);
  DATA_IVECTOR(timecov);
  DATA_INTEGER(model);

  PARAMETER_VECTOR(beta);
  PARAMETER(log_sigma);

  int n = x.rows();
  int T = x.cols();
  int m = Bs.size();
  int ind = 0.0;
  Type sigma = exp(log_sigma);

  matrix<Type> one(m+2,1);
  one.fill(1.0);
  Type nll = 0.0;

  if (model == 1){
  PARAMETER_VECTOR(alpha);
  PARAMETER_VECTOR(mu);
  vector<Type> alphat(T-1);
  for (int j = 0; j < T-1; j++){
    alphat[j] = alpha(timep[j]);
  }
  vector<Type> mut(T-1);
  for (int j = 0; j < T-1; j++){
    mut[j] = mu(timecov[j]);
  }
  vector<Type> p = invlogit(alphat);
  for (int i = 0; i < n; i++){
        matrix<Type> delta(1,m+2);
        delta.fill(0.0);
        for (int k = 1; k < m+1; k++) {
            if (y(i,fi[i]-1) < Bj[k]) {
                delta(0,k-1) = 1;
                break;
            }
        }
        matrix<Type> GQ(m+2,m+2);
        GQ.fill(0.0);
        for (int k = 0; k < m+2; k++) {
            GQ(k,k) = 1.0;
        }
        for (int j = fi[i]; j < T; j++){
            matrix<Type> Q(m+2,m+2);
            Q.fill(0.0);
            for (int k = 0; k < m; k++) {
                if (x(i,j) == 0) {
                    Q(k,k) = 1 - p[j-1];
                } else {
                    Q(k,k) = p[j-1];
                }
            }
            matrix<Type> gamma(m+2,m+2);
            gamma.fill(0.0);
            gamma(m,m+1) = 1.0;
            gamma(m+1,m+1) = 1.0;
            if(x(i,j) == 0){
                Q(m,m) = 1.0;
                Q(m+1,m+1) = 1.0;
                if (y(i,j-1) == 0){
                    for (int k = 0; k < m; k++){
                        for (int l = 0; l < m; l++){
                        Type psi = pnorm(Bj[k+1], Bs[l] + mut[j-1], sigma) - pnorm(Bj[k], Bs[l] + mut[j-1], sigma);
                        Type phi = invlogit(beta[0] + beta[1]*Bs[k]);
                        gamma(k,l) = psi*phi;
                        gamma(k,m) = 1 - phi;
                        }
                    }
                    if (j < la[i]-1){
                        gamma.col(m).fill(0.0);
                        gamma.col(m+1).fill(0.0);
                    }
                }
                else if (y(i,j-1) != 0){
                    for (int k = 1; k < m; k++) {
                        if (y(i,j-1) < Bj[k]) {
                            ind = k-1;
                            break;
                        }
                    }
                    for (int k = 0; k < m; k++){
                        Type phi = invlogit(beta[0] + beta[1]*y(i,j-1));
                        Type psi = pnorm(Bj[k+1], y(i,j-1) + mut[j-1], sigma) - pnorm(Bj[k], y(i,j-1) + mut[j-1], sigma);
                        gamma(ind,k) = phi*psi;
                        gamma(ind,m) = 1 - phi;
                    }
                    if (j < la[i]-1){
                        gamma.col(m).fill(0.0);
                        gamma.col(m+1).fill(0.0);
                    }
                }
            }
            else if(x(i,j) == 1) {
                Q(m,m) = 0.0;
                Q(m+1,m+1) = 0.0;
                gamma.fill(0.0);
                for (int k = 1; k < m; k++) {
                    if (y(i,j) < Bj[k]) {
                        ind = k-1;
                        break;
                    }
                }
                if (y(i,j-1) == 0){
                    for (int k = 0; k < m; k++){
                        Type phi = invlogit(beta[0] + beta[1]*Bs[k]);
                        Type psi = dnorm(y(i,j), Bs[k] + mut[j-1], sigma, false);
                        gamma(k,ind) = phi*psi;
                    }
                }
                else if (y(i,j-1) != 0) {
                    for (int k = 0; k < m; k++){
                        Type phi = invlogit(beta[0] + beta[1]*y(i,j-1));
                        Type psi = dnorm(y(i,j), y(i,j-1) + mut[j-1], sigma, false);
                        gamma(k,ind) = phi*psi;
                    }
                }
            }

            matrix<Type> foo = gamma*Q;
            GQ *= foo;
        }
        vector<Type> nl = GQ*one;
        nll -= log((delta*nl).sum());
  }
  }

  if (model == 2){
  PARAMETER_VECTOR(alpha);
  vector<Type> alphat(T-1);
  for (int j = 0; j < T-1; j++){
    alphat[j] = alpha(timep[j]);
  }
  vector<Type> p = invlogit(alphat);
  for (int i = 0; i < n; i++){
        matrix<Type> delta(1,m+2);
        delta.fill(0.0);
        for (int k = 1; k < m+1; k++) {
            if (y(i,fi[i]-1) < Bj[k]) {
                delta(0,k-1) = 1;
                break;
            }
        }
        matrix<Type> GQ(m+2,m+2);
        GQ.fill(0.0);
        for (int k = 0; k < m+2; k++) {
            GQ(k,k) = 1.0;
        }
        for (int j = fi[i]; j < T; j++){
            matrix<Type> Q(m+2,m+2);
            Q.fill(0.0);
            for (int k = 0; k < m; k++) {
                if (x(i,j) == 0) {
                    Q(k,k) = 1 - p[j-1];
                } else {
                    Q(k,k) = p[j-1];
                }
            }
            matrix<Type> gamma(m+2,m+2);
            gamma.fill(0.0);
            gamma(m,m+1) = 1.0;
            gamma(m+1,m+1) = 1.0;
            if(x(i,j) == 0){
                Q(m,m) = 1.0;
                Q(m+1,m+1) = 1.0;
                if (y(i,j-1) == 0){
                    for (int k = 0; k < m; k++){
                        for (int l = 0; l < m; l++){
                        Type psi = pnorm(Bj[k+1], Bs[l], sigma) - pnorm(Bj[k], Bs[l], sigma);
                        Type phi = invlogit(beta[0] + beta[1]*Bs[k]);
                        gamma(k,l) = psi*phi;
                        gamma(k,m) = 1 - phi;
                        }
                    }
                    if (j < la[i]-1){
                        gamma.col(m).fill(0.0);
                        gamma.col(m+1).fill(0.0);
                    }
                }
                else if (y(i,j-1) != 0){
                    for (int k = 1; k < m; k++) {
                        if (y(i,j-1) < Bj[k]) {
                            ind = k-1;
                            break;
                        }
                    }
                    for (int k = 0; k < m; k++){
                        Type phi = invlogit(beta[0] + beta[1]*y(i,j-1));
                        Type psi = pnorm(Bj[k+1], y(i,j-1), sigma) - pnorm(Bj[k], y(i,j-1), sigma);
                        gamma(ind,k) = phi*psi;
                        gamma(ind,m) = 1 - phi;
                    }
                    if (j < la[i]-1){
                        gamma.col(m).fill(0.0);
                        gamma.col(m+1).fill(0.0);
                    }
                }
            }
            else if(x(i,j) == 1) {
                Q(m,m) = 0.0;
                Q(m+1,m+1) = 0.0;
                gamma.fill(0.0);
                for (int k = 1; k < m; k++) {
                    if (y(i,j) < Bj[k]) {
                        ind = k-1;
                        break;
                    }
                }
                if (y(i,j-1) == 0){
                    for (int k = 0; k < m; k++){
                        Type phi = invlogit(beta[0] + beta[1]*Bs[k]);
                        Type psi = dnorm(y(i,j), Bs[k], sigma, false);
                        gamma(k,ind) = phi*psi;
                    }
                }
                else if (y(i,j-1) != 0) {
                    for (int k = 0; k < m; k++){
                        Type phi = invlogit(beta[0] + beta[1]*y(i,j-1));
                        Type psi = dnorm(y(i,j), y(i,j-1), sigma, false);
                        gamma(k,ind) = phi*psi;
                    }
                }
            }

            matrix<Type> foo = gamma*Q;
            GQ *= foo;
        }
        vector<Type> nl = GQ*one;
        nll -= log((delta*nl).sum());
  }
  }

  else if (model == 3) {
  PARAMETER_VECTOR(gamma);
  PARAMETER_VECTOR(mu);
  vector<Type> mut(T-1);
  for (int j = 0; j < T-1; j++){
    mut[j] = mu(timecov[j]);
  }
  for (int i = 0; i < n; i++){
        matrix<Type> delta(1,m+2);
        delta.fill(0.0);
        for (int k = 1; k < m+1; k++) {
            if (y(i,fi[i]-1) < Bj[k]) {
                delta(0,k-1) = 1;
                break;
            }
        }
        matrix<Type> GQ(m+2,m+2);
        GQ.fill(0.0);
        for (int k = 0; k < m+2; k++) {
            GQ(k,k) = 1.0;
        }
        for (int j = fi[i]; j < T; j++){
            matrix<Type> Q(m+2,m+2);
            Q.fill(0.0);
            for (int k = 0; k < m; k++) {
                if (x(i,j) == 0) {
                    Q(k,k) = 1 - invlogit(gamma[0] + gamma[1]*Bs[k]);;
                } else {
                    Q(k,k) = invlogit(gamma[0] + gamma[1]*y(i,j));;
                }
            }
            matrix<Type> gamma(m+2,m+2);
            gamma.fill(0.0);
            gamma(m,m+1) = 1.0;
            gamma(m+1,m+1) = 1.0;
            if(x(i,j) == 0){
                Q(m,m) = 1.0;
                Q(m+1,m+1) = 1.0;
                if (y(i,j-1) == 0){
                    for (int k = 0; k < m; k++){
                        for (int l = 0; l < m; l++){
                        Type psi = pnorm(Bj[k+1], Bs[l] + mut[j-1], sigma) - pnorm(Bj[k], Bs[l] + mut[j-1], sigma);
                        Type phi = invlogit(beta[0] + beta[1]*Bs[k]);
                        gamma(k,l) = psi*phi;
                        gamma(k,m) = 1 - phi;
                        }
                    }
                    if (j < la[i]-1){
                        gamma.col(m).fill(0.0);
                        gamma.col(m+1).fill(0.0);
                    }
                }
                else if (y(i,j-1) != 0){
                    for (int k = 1; k < m; k++) {
                        if (y(i,j-1) < Bj[k]) {
                            ind = k-1;
                            break;
                        }
                    }
                    for (int k = 0; k < m; k++){
                        Type phi = invlogit(beta[0] + beta[1]*y(i,j-1));
                        Type psi = pnorm(Bj[k+1], y(i,j-1) + mut[j-1], sigma) - pnorm(Bj[k], y(i,j-1) + mut[j-1], sigma);
                        gamma(ind,k) = phi*psi;
                        gamma(ind,m) = 1 - phi;
                    }
                    if (j < la[i]-1){
                        gamma.col(m).fill(0.0);
                        gamma.col(m+1).fill(0.0);
                    }
                }
            }
            else if(x(i,j) == 1) {
                Q(m,m) = 0.0;
                Q(m+1,m+1) = 0.0;
                gamma.fill(0.0);
                for (int k = 1; k < m; k++) {
                    if (y(i,j) < Bj[k]) {
                        ind = k-1;
                        break;
                    }
                }
                if (y(i,j-1) == 0){
                    for (int k = 0; k < m; k++){
                        Type phi = invlogit(beta[0] + beta[1]*Bs[k]);
                        Type psi = dnorm(y(i,j), Bs[k] + mut[j-1], sigma, false);
                        gamma(k,ind) = phi*psi;
                    }
                }
                else if (y(i,j-1) != 0) {
                    for (int k = 0; k < m; k++){
                        Type phi = invlogit(beta[0] + beta[1]*y(i,j-1));
                        Type psi = dnorm(y(i,j), y(i,j-1) + mut[j-1], sigma, false);
                        gamma(k,ind) = phi*psi;
                    }
                }
            }

            matrix<Type> foo = gamma*Q;
            GQ *= foo;
        }
        vector<Type> nl = GQ*one;
        nll -= log((delta*nl).sum());
  }
  }

  else if (model == 4) {
  PARAMETER_VECTOR(gamma);
  for (int i = 0; i < n; i++){
        matrix<Type> delta(1,m+2);
        delta.fill(0.0);
        for (int k = 1; k < m+1; k++) {
            if (y(i,fi[i]-1) < Bj[k]) {
                delta(0,k-1) = 1;
                break;
            }
        }
        matrix<Type> GQ(m+2,m+2);
        GQ.fill(0.0);
        for (int k = 0; k < m; k++) {
            GQ(k,k) = 1.0;
        }
        for (int j = fi[i]; j < T; j++){
            matrix<Type> Q(m+2,m+2);
            Q.fill(0.0);
            for (int k = 0; k < m; k++) {
                if (x(i,j) == 0) {
                    Q(k,k) = 1 - invlogit(gamma[0] + gamma[1]*Bs[k]);
                } else {
                    Q(k,k) = invlogit(gamma[0] + gamma[1]*y(i,j));
                }
            }
            matrix<Type> gamma(m+2,m+2);
            gamma.fill(0.0);
            gamma(m,m+1) = 1.0;
            gamma(m+1,m+1) = 1.0;
            if(x(i,j) == 0){
                Q(m,m) = 1.0;
                Q(m+1,m+1) = 1.0;
                if (y(i,j-1) == 0){
                    for (int k = 0; k < m; k++){
                        for (int l = 0; l < m; l++){
                        Type psi = pnorm(Bj[k+1], Bs[l], sigma) - pnorm(Bj[k], Bs[l], sigma);
                        Type phi = invlogit(beta[0] + beta[1]*Bs[k]);
                        gamma(k,l) = psi*phi;
                        gamma(k,m) = 1 - phi;
                        }
                    }
                    if (j < la[i]-1){
                        gamma.col(m).fill(0.0);
                        gamma.col(m+1).fill(0.0);
                    }
                }
                else if (y(i,j-1) != 0){
                    for (int k = 1; k < m; k++) {
                        if (y(i,j-1) < Bj[k]) {
                            ind = k-1;
                            break;
                        }
                    }
                    for (int k = 0; k < m; k++){
                        Type phi = invlogit(beta[0] + beta[1]*y(i,j-1));
                        Type psi = pnorm(Bj[k+1], y(i,j-1), sigma) - pnorm(Bj[k], y(i,j-1), sigma);
                        gamma(ind,k) = phi*psi;
                        gamma(ind,m) = 1 - phi;
                    }
                    if (j < la[i]-1){
                        gamma.col(m).fill(0.0);
                        gamma.col(m+1).fill(0.0);
                    }
                }
            }
            else if(x(i,j) == 1) {
                Q(m,m) = 0.0;
                Q(m+1,m+1) = 0.0;
                gamma.fill(0.0);
                for (int k = 1; k < m; k++) {
                    if (y(i,j) < Bj[k]) {
                        ind = k-1;
                        break;
                    }
                }
                if (y(i,j-1) == 0){
                    for (int k = 0; k < m; k++){
                        Type phi = invlogit(beta[0] + beta[1]*Bs[k]);
                        Type psi = dnorm(y(i,j), Bs[k], sigma, false);
                        gamma(k,ind) = phi*psi;
                    }
                }
                else if (y(i,j-1) != 0) {
                    for (int k = 0; k < m; k++){
                        Type phi = invlogit(beta[0] + beta[1]*y(i,j-1));
                        Type psi = dnorm(y(i,j), y(i,j-1), sigma, false);
                        gamma(k,ind) = phi*psi;
                    }
                }
            }

            matrix<Type> foo = gamma*Q;
            GQ *= foo;
        }
        vector<Type> nl = GQ*one;
        nll -= log((delta*nl).sum());
  }
  }
    return nll;
}
