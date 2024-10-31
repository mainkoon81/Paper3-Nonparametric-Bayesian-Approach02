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


// [[Rcpp::export]]
// func 01.
uvec histC(uvec x) {
  int xmax = max(x);
  return(hist( x, xmax ));
}

// [[Rcpp::export]]
// func 02.
vec removeelement(vec x, int n) {
  uvec n1(1);
  n1(0) = n;
  mat y(x.size(),1);
  y.col(0) = x;
  y.shed_rows(n1);
  return(y.col(0));
}

// [[Rcpp::export]]
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



// [[Rcpp::export]]
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

// func 06.
double rinvgamma(double u0, double v0) {
  double x = R::rgamma(u0, 1/v0);
  return(1/x);
}

// func 07.
// double Mean(const std::vector<double>& v) {
//   return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
// }






////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List clusterDPBasic(vec S, vec X, vec Z, uvec cl_membership,
               vec piparam, vec lambda2param, 
               mat betaparam, vec sig2param, vec xiparam, double alphaparam,
               vec f0X, vec f0S, double u0, double v0, vec m0, mat SIG_b0, double nu0,
               double g0, double h0, double c0, double d0, 
               double gamma0, double psi0, double varinf) {

  int n = S.size();
  uvec cluster_si;
  uvec nj;
  int j;
  int J; // total number of clusters
  J = max(cl_membership);
  vec x_i(3);
  x_i(0) = 1;
  
  double X_bar;
  X_bar = mean(X);
  int newclust;


  for(int i = 0; i < n; i++) {
    //Rprintf("index=%d", i);
    
    cluster_si = cl_membership;
    nj = histC(cluster_si); // counts of each cluster
  
    
    if( nj(cluster_si(i)-1)==1 ){
      j = cluster_si(i);
      
      for(int k = 0; k < n; k++) {
              if(cluster_si(k)>j) {
                cluster_si(k) = cluster_si(k)-1;
              }
      }
      
      piparam = removeelement(piparam, j-1);
      lambda2param = removeelement(lambda2param, j-1);
      
      uvec j1(1);
      j1(0) = j-1;
      betaparam.shed_rows(j1);
      //Rprintf("n.rows.beta=%d,", betaparam.n_rows);
        
      sig2param = removeelement(sig2param, j-1);
      xiparam = removeelement(xiparam, j-1);
      
      J = J-1;
    }
    //(a)
    cluster_si(i) = 0; //% replace the membership value (the first obv i) with "0"!
    nj = histC( cluster_si(find(cluster_si>0)) ); // counts from each cluster without obv i
    vec probs(nj.size()+1, fill::zeros); // vector of length nj.size()+1 
    
    //(b)
    // xi = (1, X, Z)
    x_i(1) = X(i);
    x_i(2) = Z(i);
  
    for(int j=0; j<nj.size(); j++) {
      probs(j) = nj(j)/(n-1+alphaparam)*dlogsknorm(S(i), x_i, betaparam.row(j), sig2param(j), xiparam(j) )*
      //dbinom R::dbinom
      (x_i(2)*piparam(j)+(1-x_i(2))*(1-piparam(j)))*normpdf(x_i(1), X_bar, sqrt(lambda2param(j)));
      
      //Rprintf("n.logSN.beta=%.8f,", dlogsknorm(S(i), x_i, betaparam.row(j), sig2param(j), xiparam(j) ) );
    }
    
    //(c)
    j=nj.size();
    probs(j) = alphaparam/(n-1+alphaparam)*f0S(i)*f0X(i);
    
    //(d)
    newclust = rmultinomF(probs); //sample new cluster
    //Rprintf("newcl=%d", newclust); 
    
    cl_membership = cluster_si;
    cl_membership(i) = newclust;
    //Rprintf("cl_membership=%d", cl_membership(i));
    
    nj = histC(cl_membership); // cluster counts including observation i
  
    //(e)
    if( nj(cl_membership(i)-1)==1 ){
      Rprintf("new-cl");
      
      
      piparam = appendelement(piparam, R::rbeta(g0, h0));
      //Rprintf("piparam.size=%d,", piparam.size() );
      
      lambda2param = appendelement(lambda2param, rinvgamma(c0/2, d0/2));
      //Rprintf("lamb2.size=%d,", lambda2param.size() );
      
      betaparam.insert_rows(betaparam.n_rows, mvnrnd(m0, SIG_b0*varinf).t());
      //Rprintf("n.new.cl.beta=%d,", betaparam.n_rows);
      
      sig2param = appendelement(sig2param, rinvgamma(u0/2, v0/2));
      xiparam = appendelement(xiparam, R::rt(nu0));
      
      J = J+1;
    }
  }
  
  return List::create(_["cl_membership"]=cl_membership,
                      _["piparam"]=piparam,
                      _["lambda2param"]=lambda2param,
                      _["betaparam"]=betaparam,
                      _["sig2param"]=sig2param,
                      _["xiparam"]=xiparam);


  
}
  
