#include <TMB.hpp>
#include <math.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(x);
  DATA_MATRIX(y);
  DATA_VECTOR(nk);
  DATA_INTEGER(T);

  DATA_VECTOR(w);
  DATA_VECTOR(node);
  DATA_INTEGER(niter); // Newton step for inner optimisation

  DATA_INTEGER(method); //1 = LA 2nd order, 2 = LA 4th Order, 3 = GH quadrature
  DATA_INTEGER(model); // 1 = Mh, 2 = Mth, 3 = Mbh, 4 = Mtbh

  PARAMETER_VECTOR(alphat);
  PARAMETER(log_sigma);
  PARAMETER(N);

  int n = x.rows();
  int q = w.size();

  Type sigma = exp(log_sigma);
  Type nll = 0.0;

  if (model == 1){
        Type alpha = alphat[0];
        if ( method == 1){
        Type z0 = 0;
        //solving the inner problem of the second integral
        for (int i=0; i<niter; i++){
                    Type mua = alpha + z0;
                    Type pa = invlogit(mua);
                    Type gr0 = T*pa + z0/pow(sigma,2);
                    Type he0 = T*pa*(1-pa) + pow(1/sigma,2);
                    z0 = z0 - gr0/he0;
        }
        Type mu0 = alpha + z0;
        Type p0 = invlogit(mu0);
        //compute the negative likelihood for the integral
        Type L1 = 0.0;
        L1 -= T*log(1-p0) - 0.5*pow(z0/sigma,2) - 0.5*log(pow(sigma,2)) - 0.5*log(T*p0*(1-p0) + pow(1/sigma,2));

        ADREPORT(z0);
        REPORT(z0);

        vector<Type> z(T);
        z.setZero();
        //solving the inner problem of the first integral
        for (int k = 1; k < T+1; k++){
            for (int i=0; i<niter; i++){
                Type mub = alpha + z[k-1];
                Type pb = invlogit(mub);
                Type gr = -k*(1-pb) + (T-k)*pb + z[k-1]/pow(sigma,2);
                Type he = T*pb*(1-pb) + pow(1/sigma,2);
                z[k-1] = z[k-1] - gr/he;
            }
        }

        vector<Type> mu = alpha + z;
        vector<Type> p = invlogit(mu);
        //computing the negative likelihood of the first integral
        Type L2 = 0.0;
        for (int k = 1; k < T+1; k++){
            L2 -= nk[k-1]*(k*log(p[k-1]) + (T-k)*log(1-p[k-1]) - 0.5*pow(z[k-1]/sigma,2) - 0.5*log(pow(sigma,2)) - 0.5*log(T*p[k-1]*(1-p[k-1])+pow(1/sigma,2)));
        }
        ADREPORT(z);
        REPORT(z);

        //finally, combining the first and the second integral
        nll += (N-n)*L1 + L2 - lgamma(N+1) + lgamma(N-n+1);
        }

        else if ( method == 2){
        //basically it's similar to the second-order Laplace above with additional higher order derivatives
        Type z0 = 0;
        for (int i=0; i<niter; i++){
            Type mua = alpha + z0;
            Type pa = invlogit(mua);
            Type gr0 = T*pa + z0/pow(sigma,2);
            Type he0 = T*pa*(1-pa) + pow(1/sigma,2);
            z0 = z0 - gr0/he0;
        }

        Type mu0 = alpha + z0;
        Type p0 = invlogit(mu0);
        Type he0 = T*p0*(1-p0) + pow(1/sigma,2);
        Type h30 = T*(pow(1-p0,2)*p0-pow(p0,2)*(1-p0));
        Type h40 = T*(pow(1-p0,3)*p0-4*pow(p0,2)*pow(1-p0,2)+pow(p0,3)*(1-p0));
        Type r30 = pow(h30,2)/pow(he0,3);
        Type r40 = h40/pow(he0,2);
        Type r50 = 1 + ((5*r30-3*r40)/(24));
        Type L1 = 0.0;
        L1 -= T*log(1-p0) - 0.5*pow(z0/sigma,2) - 0.5*log(pow(sigma,2)) - 0.5*log(he0) + log(r50);

        ADREPORT(z0);
        REPORT(z0);

        vector<Type> z(T);
        z.setZero();
        for (int k = 1; k < T+1; k++){
            for (int i=0; i<niter; i++){
                Type mub = alpha + z[k-1];
                Type pb = invlogit(mub);
                Type gr = -k*(1-pb) + (T-k)*pb + z[k-1]/pow(sigma,2);
                Type he = T*pb*(1-pb) + pow(1/sigma,2);
                z[k-1] = z[k-1] - gr/he;
            }
        }

        vector<Type> mu = alpha + z;
        vector<Type> p = invlogit(mu);

        Type L2 = 0.0;
        for (int k = 1; k < T+1; k++){
            Type he = T*p[k-1]*(1-p[k-1])+pow(1/sigma,2);
            Type h3 = T*(pow(1-p[k-1],2)*p[k-1]-pow(p[k-1],2)*(1-p[k-1]));
            Type h4 = T*(pow(1-p[k-1],3)*p[k-1]-4*pow(p[k-1],2)*pow(1-p[k-1],2)+pow(p[k-1],3)*(1-p[k-1]));
            Type r3 = pow(h3,2)/pow(he,3);
            Type r4 = h4/pow(he,2);
            Type r5 = 1 + ((5*r3-3*r4)/(24));
            L2 -= nk[k-1]*(k*log(p[k-1]) + (T-k)*log(1-p[k-1]) - 0.5*pow(z[k-1]/sigma,2) - 0.5*log(pow(sigma,2)) - 0.5*log(he) + log(r5));
        }
        ADREPORT(z);
        REPORT(z);

        nll += (N-n)*L1 + L2 - lgamma(N+1) + lgamma(N-n+1);
        }

        else if (method == 3){
        //computing the approximation of the first integral
        Type w1 = 0.0;
        for (int i = 0; i < q; i++){
                w1 += w[i]/sqrt(M_PI)*(1/pow((1+exp(sqrt(2)*node[i]*sigma+alpha)),T));
        }
        Type L1 = -log(w1);

        //computing the approximation of the second integral
        Type L2 = 0;
        for (int k = 1; k < T+1; k++){
                Type w2 = 0.0;
                for (int i = 0; i < q; i++){
                    w2 += w[i]/sqrt(M_PI)*(1/pow((1+exp(sqrt(2)*node[i]*sigma+alpha)),T-k))*(1/pow((1+exp(-sqrt(2)*node[i]*sigma-alpha)),k));
                }
        L2 -= nk[k-1]*log(w2);
        }
        //combining the first and the second integral
        nll += (N-n)*L1 + L2 - lgamma(N+1) + lgamma(N-n+1);
        }
  }

  //the rest of coding is similar coded for different models.
  else if (model == 2){
        vector<Type> alpha = alphat;
        if(method == 1){
        Type z0 = 0.0;
        vector<Type> mua(T);
        vector<Type> pa(T);
        vector<Type> qa(T);
        for (int i=0; i<niter; i++){
            for (int k = 0; k < T; k++){
                mua[k] = alpha[k] + z0;
                pa[k] = invlogit(mua[k]);
                qa[k] = (1- pa[k])*pa[k];
            }
        Type gr0 = pa.sum() + z0/pow(sigma,2);
        Type he0 = qa.sum() + pow(1/sigma,2);
        z0 = z0 - gr0/he0;
        }

        vector<Type> mu0 = alpha + z0;
        vector<Type> p0 = invlogit(mu0);
        vector<Type> h0(T);
        vector<Type> q0(T);
        Type L1 = 0.0;
        for (int k = 0; k < T; k++){
            h0[k] = (1-p0[k])*p0[k];
            q0[k] = log(1-p0[k]);
        }
        L1 -= q0.sum() - 0.5*pow(z0/sigma,2) - 0.5*log(pow(sigma,2)) - 0.5*log(h0.sum() + pow(1/sigma,2));

        ADREPORT(z0);
        REPORT(z0);

        vector<Type> z(n);
        z.setZero();
        vector<Type> mb(T);
        vector<Type> pb(T);
        vector<Type> gb(T);
        vector<Type> hb(T);
        for (int i = 0; i < n; i++){
            for (int j=0; j<niter; j++){
                for (int k = 0; k < T; k++){
                    mb[k] = alpha[k] + z[i];
                    pb[k] = invlogit(mb[k]);
                    gb[k] = -x(i,k)*(1-pb[k]) + (1-x(i,k))*pb[k];
                    hb[k] = (1-pb[k])*pb[k];
                }
            Type gr = gb.sum() + z[i]/pow(sigma,2);
            Type he = hb.sum() + pow(1/sigma,2);
            z[i] = z[i] - gr/he;
            }
        }

        matrix<Type> mu(n,T);
        matrix<Type> p(n,T);
        for (int i = 0; i < n; i++){
            for (int k = 0; k < T; k++){
                mu(i,k)= alpha[k] + z[i];
                p(i,k) = invlogit(mu(i,k));
            }
        }

        vector<Type> h1(T);
        vector<Type> q1(T);
        Type L2 = 0.0;
        for (int i = 0; i < n; i++){
            for (int k = 0; k < T; k++){
                h1[k] = p(i,k)*(1-p(i,k));
                q1[k] = x(i,k)*log(p(i,k)) + (1-x(i,k))*log(1-p(i,k));
            }
        L2 -= q1.sum() - 0.5*pow(z[i]/sigma,2) - 0.5*log(pow(sigma,2)) - 0.5*log(h1.sum()+pow(1/sigma,2));
        }

        ADREPORT(z);
        REPORT(z);

        nll += (N-n)*L1 + L2 - lgamma(N+1) + lgamma(N-n+1);
        }

        else if (method == 2) {
        Type z0 = 0.0;
        vector<Type> mua(T);
        vector<Type> pa(T);
        vector<Type> qa(T);
        for (int i=0; i<niter; i++){
            for (int k = 0; k < T; k++){
                mua[k] = alpha[k] + z0;
                pa[k] = invlogit(mua[k]);
                qa[k] = (1- pa[k])*pa[k];
            }
        Type gr0 = pa.sum() + z0/pow(sigma,2);
        Type he0 = qa.sum() + pow(1/sigma,2);
        z0 = z0 - gr0/he0;
        }
        vector<Type> mu0 = alpha + z0;
        vector<Type> p0 = invlogit(mu0);
        vector<Type> h0(T);
        vector<Type> q0(T);
        vector<Type> h30(T);
        vector<Type> h40(T);

        Type L1 = 0.0;
        for (int k = 0; k < T; k++){
            h0[k] = (1-p0[k])*p0[k];
            q0[k] = log(1-p0[k]);
            h30[k] = pow(1-p0[k],2)*p0[k]-pow(p0[k],2)*(1-p0[k]);
            h40[k] = pow(1-p0[k],3)*p0[k]-4*pow(p0[k],2)*pow(1-p0[k],2)+pow(p0[k],3)*(1-p0[k]);
        }
        Type h20 = h0.sum() + pow(1/sigma,2);
        Type r30 = pow(h30.sum(),2)/pow(h20,3);
        Type r40 = h40.sum()/pow(h20,2);
        Type r50 = 1 + ((5*r30-3*r40)/(24));

        L1 -= q0.sum() - 0.5*pow(z0/sigma,2) - 0.5*log(pow(sigma,2)) - 0.5*log(h20) + log(r50);

        ADREPORT(z0);
        REPORT(z0);

        vector<Type> z(n);
        z.setZero();
        vector<Type> mb(T);
        vector<Type> pb(T);
        vector<Type> gb(T);
        vector<Type> hb(T);
        for (int i = 0; i < n; i++){
            for (int j=0; j<niter; j++){
                for (int k = 0; k < T; k++){
                    mb[k] = alpha[k] + z[i];
                    pb[k] = invlogit(mb[k]);
                    gb[k] = -x(i,k)*(1-pb[k]) + (1-x(i,k))*pb[k];
                    hb[k] = (1-pb[k])*pb[k];
                }
            Type gr = gb.sum() + z[i]/pow(sigma,2);
            Type he = hb.sum() + pow(1/sigma,2);
            z[i] = z[i] - gr/he;
            }
        }
        matrix<Type> mu(n,T);
        matrix<Type> p(n,T);
        for (int i = 0; i < n; i++){
            for (int k = 0; k < T; k++){
                mu(i,k)= alpha[k] + z[i];
                p(i,k) = invlogit(mu(i,k));
            }
        }

        vector<Type> h1(T);
        vector<Type> q1(T);
        vector<Type> h3(T);
        vector<Type> h4(T);

        Type L2 = 0.0;
        for (int i = 0; i < n; i++){
            for (int k = 0; k < T; k++){
                q1[k] = x(i,k)*log(p(i,k)) + (1-x(i,k))*log(1-p(i,k));
                h1[k] = p(i,k)*(1-p(i,k));
                h3[k] = p(i,k)*pow(1-p(i,k),2)-pow(p(i,k),2)*(1-p(i,k));
                h4[k] = p(i,k)*pow(1-p(i,k),3)-4*pow(p(i,k),2)*pow(1-p(i,k),2)+pow(p(i,k),3)*(1-p(i,k));
            }
        Type h2 = h1.sum() + pow(1/sigma,2);
        Type r3 = pow(h3.sum(),2)/pow(h2,3);
        Type r4 = h4.sum()/pow(h2,2);
        Type r5 = 1 + ((5*r3-3*r4)/(24));
        L2 -= q1.sum() - 0.5*pow(z[i]/sigma,2) - 0.5*log(pow(sigma,2)) - 0.5*log(h2) + log(r5);
        }
        ADREPORT(z);
        REPORT(z);

        nll += (N-n)*L1 + L2 - lgamma(N+1) + lgamma(N-n+1);
        }

        else if (method == 3){
        Type w1 = 0.0;
        for (int i = 0; i < q; i++){
                Type q1 = 1.0;
                for (int k = 0; k < T; k++){
                    q1 *= 1/(1+exp(sqrt(2)*node[i]*sigma+alpha[k]));
                }
                w1 += w[i]/sqrt(M_PI)*q1;
            }
        Type L1 = -log(w1);

        Type L2 = 0;
        for (int i = 0; i < n; i++){
            Type w2 = 0.0;
            for (int j = 0; j < q; j++){
                Type q2 = 1.0;
                for (int k = 0; k < T; k++){
                    q2 *= (1/pow((1+exp(sqrt(2)*node[j]*sigma+alpha[k])),1-x(i,k)))*(1/pow((1+exp(-sqrt(2)*node[j]*sigma-alpha[k])),x(i,k)));
                }
                w2 += w[j]/sqrt(M_PI)*q2;
            }
        L2 -= log(w2);
        }

        nll += (N-n)*L1 + L2 - lgamma(N+1) + lgamma(N-n+1);
        }

  }
  else if (model == 3){
        PARAMETER(b);
        Type alpha = alphat[0];
        if(method == 1){
        Type z0 = 0;
        for (int i=0; i<niter; i++){
                    Type mua = alpha + z0;
                    Type pa = invlogit(mua);
                    Type gr0 = T*pa + z0/pow(sigma,2);
                    Type he0 = T*(1-pa)*pa + pow(1/sigma,2);
                    z0 = z0 - gr0/he0;
        }
        Type mu0 = alpha + z0;
        Type p0 = invlogit(mu0);
        Type h20 = T*(1-p0)*p0 + pow(1/sigma,2);
        Type L1 = 0.0;
        L1 -= T*log(1-p0) - 0.5*pow(z0/sigma,2) - 0.5*log(pow(sigma,2)) - 0.5*log(h20);

        ADREPORT(z0);
        REPORT(z0);

        vector<Type> z(n);
        z.setZero();
        vector<Type> mb(T);
        vector<Type> pb(T);
        vector<Type> gb(T);
        vector<Type> hb(T);
        for (int i = 0; i < n; i++){
            for (int j=0; j<niter; j++){
                for (int k = 0; k < T; k++){
                    mb[k] = alpha + z[i] + b*y(i,k);
                    pb[k] = invlogit(mb[k]);
                    gb[k] = -x(i,k)*(1-pb[k]) + (1-x(i,k))*pb[k];
                    hb[k] = (1-pb[k])*pb[k];
                }
            Type gr = gb.sum() + z[i]/pow(sigma,2);
            Type he = hb.sum() + pow(1/sigma,2);
            z[i] = z[i] - gr/he;
            }
        }

        matrix<Type> mu(n,T);
        matrix<Type> p(n,T);
        for (int i = 0; i < n; i++){
            for (int k = 0; k < T; k++){
                mu(i,k)= alpha + z[i] + b*y(i,k);
                p(i,k) = invlogit(mu(i,k));
            }
        }

        vector<Type> h1(T);
        vector<Type> q1(T);
        Type L2 = 0.0;
        for (int i = 0; i < n; i++){
            for (int k = 0; k < T; k++){
                h1[k] = p(i,k)*(1-p(i,k));
                q1[k] = x(i,k)*log(p(i,k)) + (1-x(i,k))*log(1-p(i,k));
            }
        L2 -= q1.sum() - 0.5*pow(z[i]/sigma,2) - 0.5*log(pow(sigma,2)) - 0.5*log(h1.sum()+pow(1/sigma,2));
        }

        ADREPORT(z);
        REPORT(z);

        nll += (N-n)*L1 + L2 - lgamma(N+1) + lgamma(N-n+1);
        }

        else if (method == 2) {
        Type z0 = 0;
        for (int i=0; i<niter; i++){
            Type mua = alpha + z0;
            Type pa = invlogit(mua);
            Type gr0 = T*pa + z0/pow(sigma,2);
            Type he0 = T*pa*(1-pa) + pow(1/sigma,2);
            z0 = z0 - gr0/he0;
        }

        Type mu0 = alpha + z0;
        Type p0 = invlogit(mu0);
        Type he0 = T*p0*(1-p0) + pow(1/sigma,2);
        Type h30 = T*(pow(1-p0,2)*p0-pow(p0,2)*(1-p0));
        Type h40 = T*(pow(1-p0,3)*p0-4*pow(p0,2)*pow(1-p0,2)+pow(p0,3)*(1-p0));
        Type r30 = pow(h30,2)/pow(he0,3);
        Type r40 = h40/pow(he0,2);
        Type r50 = 1 + ((5*r30-3*r40)/(24));
        Type L1 = 0.0;
        L1 -= T*log(1-p0) - 0.5*pow(z0/sigma,2) - 0.5*log(pow(sigma,2)) - 0.5*log(he0) + log(r50);

        ADREPORT(z0);
        REPORT(z0);

        vector<Type> z(n);
        z.setZero();
        vector<Type> mb(T);
        vector<Type> pb(T);
        vector<Type> gb(T);
        vector<Type> hb(T);
        for (int i = 0; i < n; i++){
            for (int j=0; j<niter; j++){
                for (int k = 0; k < T; k++){
                    mb[k] = alpha + z[i] + b*y(i,k);
                    pb[k] = invlogit(mb[k]);
                    gb[k] = -x(i,k)*(1-pb[k]) + (1-x(i,k))*pb[k];
                    hb[k] = (1-pb[k])*pb[k];
                }
            Type gr = gb.sum() + z[i]/pow(sigma,2);
            Type he = hb.sum() + pow(1/sigma,2);
            z[i] = z[i] - gr/he;
            }
        }
        matrix<Type> mu(n,T);
        matrix<Type> p(n,T);
        for (int i = 0; i < n; i++){
            for (int k = 0; k < T; k++){
                mu(i,k)= alpha + z[i] + b*y(i,k);
                p(i,k) = invlogit(mu(i,k));
            }
        }

        vector<Type> h1(T);
        vector<Type> q1(T);
        vector<Type> h3(T);
        vector<Type> h4(T);

        Type L2 = 0.0;
        for (int i = 0; i < n; i++){
            for (int k = 0; k < T; k++){
                q1[k] = x(i,k)*log(p(i,k)) + (1-x(i,k))*log(1-p(i,k));
                h1[k] = p(i,k)*(1-p(i,k));
                h3[k] = p(i,k)*pow(1-p(i,k),2)-pow(p(i,k),2)*(1-p(i,k));
                h4[k] = p(i,k)*pow(1-p(i,k),3)-4*pow(p(i,k),2)*pow(1-p(i,k),2)+pow(p(i,k),3)*(1-p(i,k));
            }
        Type h2 = h1.sum() + pow(1/sigma,2);
        Type r3 = pow(h3.sum(),2)/pow(h2,3);
        Type r4 = h4.sum()/pow(h2,2);
        Type r5 = 1 + ((5*r3-3*r4)/(24));
        L2 -= q1.sum() - 0.5*pow(z[i]/sigma,2) - 0.5*log(pow(sigma,2)) - 0.5*log(h2) + log(r5);
        }
        ADREPORT(z);
        REPORT(z);

        nll += (N-n)*L1 + L2 - lgamma(N+1) + lgamma(N-n+1);
        }

        else if (method == 3){
        Type w1 = 0.0;
        for (int i = 0; i < q; i++){
                w1 += w[i]/sqrt(M_PI)*(1/pow((1+exp(sqrt(2)*node[i]*sigma+alpha)),T));
        }
        Type L1 = -log(w1);

        Type L2 = 0;
        for (int i = 0; i < n; i++){
            Type w2 = 0.0;
            for (int j = 0; j < q; j++){
                Type q2 = 1.0;
                for (int k = 0; k < T; k++){
                    q2 *= (1/pow((1+exp(sqrt(2)*node[j]*sigma+alpha+b*y(i,k))),1-x(i,k)))*(1/pow((1+exp(-sqrt(2)*node[j]*sigma-alpha-b*y(i,k))),x(i,k)));
                }
                w2 += w[j]/sqrt(M_PI)*q2;
            }
        L2 -= log(w2);
        }

        nll += (N-n)*L1 + L2 - lgamma(N+1) + lgamma(N-n+1);
        }
  }

  else if (model == 4){
        PARAMETER(b);
        vector<Type> alpha = alphat;
        if(method == 1){
         Type z0 = 0.0;
        vector<Type> mua(T);
        vector<Type> pa(T);
        vector<Type> qa(T);
        for (int i=0; i<niter; i++){
            for (int k = 0; k < T; k++){
                mua[k] = alpha[k] + z0;
                pa[k] = invlogit(mua[k]);
                qa[k] = (1- pa[k])*pa[k];
            }
        Type gr0 = pa.sum() + z0/pow(sigma,2);
        Type he0 = qa.sum() + pow(1/sigma,2);
        z0 = z0 - gr0/he0;
        }

        vector<Type> mu0 = alpha + z0;
        vector<Type> p0 = invlogit(mu0);
        vector<Type> h0(T);
        vector<Type> q0(T);
        Type L1 = 0.0;
        for (int k = 0; k < T; k++){
            h0[k] = (1-p0[k])*p0[k];
            q0[k] = log(1-p0[k]);
        }
        L1 -= q0.sum() - 0.5*pow(z0/sigma,2) - 0.5*log(pow(sigma,2)) - 0.5*log(h0.sum() + pow(1/sigma,2));

        ADREPORT(z0);
        REPORT(z0);

        vector<Type> z(n);
        z.setZero();
        vector<Type> mb(T);
        vector<Type> pb(T);
        vector<Type> gb(T);
        vector<Type> hb(T);
        for (int i = 0; i < n; i++){
            for (int j=0; j<niter; j++){
                for (int k = 0; k < T; k++){
                    mb[k] = alpha[k] + z[i] + b*y(i,k);
                    pb[k] = invlogit(mb[k]);
                    gb[k] = -x(i,k)*(1-pb[k]) + (1-x(i,k))*pb[k];
                    hb[k] = (1-pb[k])*pb[k];
                }
            Type gr = gb.sum() + z[i]/pow(sigma,2);
            Type he = hb.sum() + pow(1/sigma,2);
            z[i] = z[i] - gr/he;
            }
        }

        matrix<Type> mu(n,T);
        matrix<Type> p(n,T);
        for (int i = 0; i < n; i++){
            for (int k = 0; k < T; k++){
                mu(i,k)= alpha[k] + z[i] + b*y(i,k);
                p(i,k) = invlogit(mu(i,k));
            }
        }

        vector<Type> h1(T);
        vector<Type> q1(T);
        Type L2 = 0.0;
        for (int i = 0; i < n; i++){
            for (int k = 0; k < T; k++){
                h1[k] = p(i,k)*(1-p(i,k));
                q1[k] = x(i,k)*log(p(i,k)) + (1-x(i,k))*log(1-p(i,k));
            }
        L2 -= q1.sum() - 0.5*pow(z[i]/sigma,2) - 0.5*log(pow(sigma,2)) - 0.5*log(h1.sum()+pow(1/sigma,2));
        }

        ADREPORT(z);
        REPORT(z);

        nll += (N-n)*L1 + L2 - lgamma(N+1) + lgamma(N-n+1);
        }

        else if (method == 2) {
        Type z0 = 0.0;
        vector<Type> mua(T);
        vector<Type> pa(T);
        vector<Type> qa(T);
        for (int i=0; i<niter; i++){
            for (int k = 0; k < T; k++){
                mua[k] = alpha[k] + z0;
                pa[k] = invlogit(mua[k]);
                qa[k] = (1- pa[k])*pa[k];
            }
        Type gr0 = pa.sum() + z0/pow(sigma,2);
        Type he0 = qa.sum() + pow(1/sigma,2);
        z0 = z0 - gr0/he0;
        }
        vector<Type> mu0 = alpha + z0;
        vector<Type> p0 = invlogit(mu0);
        vector<Type> h0(T);
        vector<Type> q0(T);
        vector<Type> h30(T);
        vector<Type> h40(T);

        Type L1 = 0.0;
        for (int k = 0; k < T; k++){
            h0[k] = (1-p0[k])*p0[k];
            q0[k] = log(1-p0[k]);
            h30[k] = pow(1-p0[k],2)*p0[k]-pow(p0[k],2)*(1-p0[k]);
            h40[k] = pow(1-p0[k],3)*p0[k]-4*pow(p0[k],2)*pow(1-p0[k],2)+pow(p0[k],3)*(1-p0[k]);
        }
        Type h20 = h0.sum() + pow(1/sigma,2);
        Type r30 = pow(h30.sum(),2)/pow(h20,3);
        Type r40 = h40.sum()/pow(h20,2);
        Type r50 = 1 + ((5*r30-3*r40)/(24));

        L1 -= q0.sum() - 0.5*pow(z0/sigma,2) - 0.5*log(pow(sigma,2)) - 0.5*log(h20) + log(r50);

        ADREPORT(z0);
        REPORT(z0);

        vector<Type> z(n);
        z.setZero();
        vector<Type> mb(T);
        vector<Type> pb(T);
        vector<Type> gb(T);
        vector<Type> hb(T);
        for (int i = 0; i < n; i++){
            for (int j=0; j<niter; j++){
                for (int k = 0; k < T; k++){
                    mb[k] = alpha[k] + z[i] + b*y(i,k);
                    pb[k] = invlogit(mb[k]);
                    gb[k] = -x(i,k)*(1-pb[k]) + (1-x(i,k))*pb[k];
                    hb[k] = (1-pb[k])*pb[k];
                }
            Type gr = gb.sum() + z[i]/pow(sigma,2);
            Type he = hb.sum() + pow(1/sigma,2);
            z[i] = z[i] - gr/he;
            }
        }
        matrix<Type> mu(n,T);
        matrix<Type> p(n,T);
        for (int i = 0; i < n; i++){
            for (int k = 0; k < T; k++){
                mu(i,k)= alpha[k] + z[i] + b*y(i,k);
                p(i,k) = invlogit(mu(i,k));
            }
        }

        vector<Type> h1(T);
        vector<Type> q1(T);
        vector<Type> h3(T);
        vector<Type> h4(T);

        Type L2 = 0.0;
        for (int i = 0; i < n; i++){
            for (int k = 0; k < T; k++){
                q1[k] = x(i,k)*log(p(i,k)) + (1-x(i,k))*log(1-p(i,k));
                h1[k] = p(i,k)*(1-p(i,k));
                h3[k] = p(i,k)*pow(1-p(i,k),2)-pow(p(i,k),2)*(1-p(i,k));
                h4[k] = p(i,k)*pow(1-p(i,k),3)-4*pow(p(i,k),2)*pow(1-p(i,k),2)+pow(p(i,k),3)*(1-p(i,k));
            }
        Type h2 = h1.sum() + pow(1/sigma,2);
        Type r3 = pow(h3.sum(),2)/pow(h2,3);
        Type r4 = h4.sum()/pow(h2,2);
        Type r5 = 1 + ((5*r3-3*r4)/(24));
        L2 -= q1.sum() - 0.5*pow(z[i]/sigma,2) - 0.5*log(pow(sigma,2)) - 0.5*log(h2) + log(r5);
        }
        ADREPORT(z);
        REPORT(z);

        nll += (N-n)*L1 + L2 - lgamma(N+1) + lgamma(N-n+1);
        }

        else if (method == 3){
        Type w1 = 0.0;
        for (int i = 0; i < q; i++){
                Type q1 = 1.0;
                for (int k = 0; k < T; k++){
                    q1 *= 1/(1+exp(sqrt(2)*node[i]*sigma+alpha[k]));
                }
                w1 += w[i]/sqrt(M_PI)*q1;
            }
        Type L1 = -log(w1);

        Type L2 = 0;
        for (int i = 0; i < n; i++){
            Type w2 = 0.0;
            for (int j = 0; j < q; j++){
                Type q2 = 1.0;
                for (int k = 0; k < T; k++){
                    q2 *= (1/pow((1+exp(sqrt(2)*node[j]*sigma+alpha[k]+b*y(i,k))),1-x(i,k)))*(1/pow((1+exp(-sqrt(2)*node[j]*sigma-alpha[k]-b*y(i,k))),x(i,k)));
                }
                w2 += w[j]/sqrt(M_PI)*q2;
            }
        L2 -= log(w2);
        }

        nll += (N-n)*L1 + L2 - lgamma(N+1) + lgamma(N-n+1);
        }
  }

  ADREPORT(sigma);
  REPORT(sigma);

  return (nll);
}
