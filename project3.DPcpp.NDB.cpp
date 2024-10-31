#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <Rmath.h>
#include <iostream>
#include <algorithm>
#include <math.h> 
#include <random>
//#include "utilFunctions.h"
  
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;
using namespace stats;

void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);  
}

// double sigmoid(double x) {
//   return( ( 1 / (1 + exp(- x) ) ) );
// }


// func 01.
uvec histC(uvec x) {
  int xmax = max(x);
  return(hist( x, xmax ));
}

// func 02.
vec removeelement(vec x, int n) {
  uvec n1(1);
  n1(0) = n;
  mat y(x.size(),1);
  y.col(0) = x;
  y.shed_rows(n1);
  return(y.col(0));
}

// func 03.
vec appendelement(vec x, double a) {
  vec a1(1);
  a1(0) = a;
  mat y(x.size(),1);
  y.col(0) = x;
  y.insert_rows(x.size(),a1);
  return(y.col(0));
}

// func 04.
double dlogsknorm(double s, vec x, rowvec beta, double sig2, double xi) {
  double mu = dot(x,beta);
  double sig = sqrt(sig2);
  double z = (log(s)-mu)/sig;
  return(2/(s*sig)*normpdf(z)*normcdf(xi*z));
}

double dlogsknorm_log(double s, vec x, rowvec beta, double sig2, double xi) {
  double mu = dot(x,beta);
  double sig = sqrt(sig2);
  double z = (log(s)-mu)/sig;
  return(log(2)-log(s*sig) + log_normpdf(z) + log(normcdf(xi*z)));
}

// func 05.
int rmultinomF(vec const& p) {
  vec csp = cumsum(p/sum(p));
  double rnd = runif(1)[0];
  int res = 0;
  int psize = p.size();
  for(int i = 0; i < psize; i++) {
    if(rnd>csp(i)) res = res+1;
  }
  return(res+1);
}

int rmultinomFlog(vec const& logp) {
  double logptemp = logp(0);
  for(int i=1; i<logp.size(); i++) {
    if(logp(i)>logptemp) {
      logptemp = logp(i)+log(1+exp(logptemp-logp(i)));
    } else {
      logptemp = logptemp+log(1+exp(logp(i)-logptemp));
    }
  }
  vec p = exp(logp - logptemp);
  vec csp = cumsum(p);
  double rnd = runif(1)[0];
  int res = 0;
  int psize = p.size();
  for(int i = 0; i < psize; i++) {
    if(rnd>csp(i)) res = res+1;
  }
  return(res+1);
}

// func 07.
double rinvgamma(double u0, double v0) {
  double x = R::rgamma(u0, 1/v0);
  return(1/x);
}




//////////////////////////////////////////////////////////////////////////////// + Kappaparam and its friends 

// [[Rcpp::export]]
List clusterDP(vec S, vec X, vec Z, uvec cl_membership,
               vec piparam, mat kappaparam, vec lambda2param, 
               mat betaparam, vec sig2param, vec xiparam, double alphaparam,
               vec f0X, vec f0S, double u0, double v0, vec m0, mat SIG_b0, double nu0,
               double g0, double h0, vec kappa_v0, mat SIG_k0, mat SIG_inv_k0, double c0, double d0, 
               double gamma0, double psi0, double varinf) {
  
  int n = S.size();
  uvec cluster_si;
  uvec nj;
  int j;
  int J; // total number of clusters
  J = max(cl_membership);
  vec x_i(3);
  
  x_i(0) = 1; // for convenience
  
  double X_bar;
  X_bar = mean(X);
  int newclust;  

  // vec beta_old_j;
  // double sig2_old_j;
  // double xi_old_j;
  // vec beta_p;
  // double sig2_p;
  // double xi_p;
  // double numerator;
  // double denominator;
  // double ratio;
  // vec compare(2);
  // double U;

  mat KK1j(1,2); // (1,Z)
  KK1j(0,0) = 1;
  mat KK2j(2,1); // (X*, X*z)
  

  for(int i = 0; i < n; i++) {
    //Rprintf("i=%d,",i);
    
    cluster_si = cl_membership;
    nj = histC(cluster_si); // counts of each cluster
    
    // if observation i is the only obs in its cluster, we remove that obs from 
    // cluster and take out cluster from list of available clusters
    if( nj(cluster_si(i)-1)==1 ){
      j = cluster_si(i);
      
      for(int k = 0; k < n; k++) {
        if(cluster_si(k)>j) {
          cluster_si(k) = cluster_si(k)-1;
        }
      }
      
      // for Z
      piparam = removeelement(piparam, j-1);
      
      // for X
      lambda2param = removeelement(lambda2param, j-1);
      
      uvec j1(1);
      j1(0) = j-1;
      kappaparam.shed_rows(j1);
      
      // for S
      betaparam.shed_rows(j1);
      sig2param = removeelement(sig2param, j-1);
      xiparam = removeelement(xiparam, j-1);
      
      J = J-1;
    }
    
    //(a) put the value 0 in the cluster membership vector for observation i
    cluster_si(i) = 0; //% replace the membership value (the first obv i) with "0"!
    nj = histC( cluster_si(find(cluster_si>0)) ); // counts from each cluster without obv i
    vec probs(nj.size()+1, fill::zeros); // vector of length nj.size()+

    
    //(b)
    // xi = (1, X, Z)
    x_i(1) = X(i);
    x_i(2) = Z(i);

    for(int j=0; j<nj.size(); j++) {
      probs(j) = nj(j)/(n-1+alphaparam)*dlogsknorm(S(i), x_i, betaparam.row(j), sig2param(j), xiparam(j) )*
        //dbinom R::dbinom
        (x_i(2)*piparam(j)+(1-x_i(2))*(1-piparam(j)))*
        normpdf(x_i(1), kappaparam(j,0) + kappaparam(j,1)*x_i(2), sqrt(lambda2param(j)));
      //Rprintf("n.logSN.beta=%.8f,", dlogsknorm(S(i), x_i, betaparam.row(j), sig2param(j), xiparam(j) ) );
      //Rprintf("prob val.solid=%.8f,", probs(j) );
    }
    
    //(c)
    j=nj.size();
    probs(j) = alphaparam/(n-1+alphaparam)*f0S(i)*f0X(i);
    //Rprintf("prob val.liquid=%.8f,", probs(j) );
    
    //(d) assigned a new cluster to observation i
    newclust = rmultinomF(probs); //sample new cluster

    cl_membership = cluster_si;
    cl_membership(i) = newclust;
    nj = histC(cl_membership); // cluster counts including observation i

    //(e)
    if( nj(cl_membership(i)-1)==1 ){
      Rprintf("new-cl");

      piparam = appendelement(piparam, R::rbeta(g0, h0));
      //Rprintf("piparam.size=%d,", piparam.size() );
      
      lambda2param = appendelement(lambda2param, rinvgamma(c0/2, d0/2));
      //Rprintf("lamb2.size=%d,", lambda2param.size() );
      
      KK1j(0,1) = x_i(2);
      KK2j(0,0) = x_i(1);
      KK2j(1,0) = x_i(1)*x_i(2);
      mat SIG_j = SIG_inv_k0+1/lambda2param(nj.size()-1)*(KK1j.t() * KK1j);
      vec kappatemp = SIG_j.i() * (SIG_inv_k0 * kappa_v0 + 1/lambda2param(nj.size()-1)*KK2j);
      kappaparam.insert_rows(kappaparam.n_rows, mvnrnd(kappatemp, SIG_j.i()).t());
      
      betaparam.insert_rows(betaparam.n_rows, mvnrnd(m0, SIG_b0*varinf).t());
      sig2param = appendelement(sig2param, rinvgamma(u0/2, v0/2));
      xiparam = appendelement(xiparam, R::rt(nu0));
      
      J = J+1;
    

      /////////////////////////// we won't need MH /////////////////////////////
      // beta_old_j = mvnrnd(m0, SIG_b0*varinf);
      // sig2_old_j = rinvgamma(u0, v0);
      // xi_old_j = R::rt(nu0); 
      //   
      // beta_p = mvnrnd(m0, SIG_b0*varinf);
      // sig2_p = rinvgamma(u0, v0);
      // xi_p = R::rt(nu0); 
      // 
      // //Rprintf("sig2_old_j=%.4f,sig2_p=%.4f,",sig2_old_j,sig2_p);
      // //Rprintf("x_i(1)=%.4f,",x_i(1));
      // numerator = dlogsknorm_log(S(i), x_i, beta_p.t(), sig2_p, xi_p);
      // denominator = dlogsknorm_log(S(i), x_i, beta_old_j.t(), sig2_old_j, xi_old_j);
      // 
      // // compute the ratio
      // compare(0) = exp(numerator-denominator);
      // compare(1) = 1;
      // ratio = min(compare);
      // //Rprintf("ratio=%.4f,",ratio);
      //   
      // U = randu();
      // if(U < ratio) {
      //   betaparam.insert_rows(betaparam.n_rows, beta_p.t());
      //   sig2param = appendelement(sig2param,sig2_p);
      //   xiparam = appendelement(xiparam,xi_p);
      // } else {
      //   betaparam.insert_rows(betaparam.n_rows, beta_old_j.t());
      //   sig2param = appendelement(sig2param,sig2_old_j);
      //   xiparam = appendelement(xiparam,xi_old_j);
      // }
      //////////////////////////////////////////////////////////////////////////
    }
  }
  
  return List::create(_["cl_membership"]=cl_membership,
                      _["piparam"]=piparam,
                      _["kappaparam"]=kappaparam,
                      _["lambdaparam"]=lambda2param,
                      _["betaparam"]=betaparam,
                      _["sig2param"]=sig2param,
                      _["xiparam"]=xiparam);
}

