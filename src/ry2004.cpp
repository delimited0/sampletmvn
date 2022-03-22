#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include "tuvn.h"

// [[Rcpp::export]]
arma::mat ry2004_gibbs_iter_lg2015(int n, arma::vec alpha, arma::vec z, 
                                   arma::mat D, arma::vec Dz, 
                                   arma::vec lb, arma::vec ub) {
  
  int p = lb.n_elem;
  
  arma::mat samples(p, n);
  arma::vec Dj;
  
  for (int i = 0; i < n; i++) {
    
    for (int j = 0; j < p; j++) {
      
      Dj = D.col(j);
      
      arma::vec Dz_notj = Dz - (Dj * z(j));
      
      arma::vec a = (lb - Dz_notj) / Dj;
      arma::vec b = (ub - Dz_notj) / Dj;
      
      double lower_j = arma::max(a);
      double upper_j = arma::min(b);
      
      double zj_old = z(j);
      z(j) = rtuvn_lg2015(alpha(j), 1., lower_j, upper_j);
      
      Dz += Dj * (z(j) - zj_old);
    }
    
    samples.col(i) = z;
  }
  
  return(samples);
}

// [[Rcpp::export]]
arma::mat ry2004_gibbs_iter_be2017(int n, arma::vec alpha, arma::vec z, 
                                   arma::mat D, arma::vec Dz, 
                                   arma::vec lb, arma::vec ub) {
  
  int p = lb.n_elem;
  
  arma::mat samples(p, n);
  arma::vec Dj;
  
  for (int i = 0; i < n; i++) {
    
    for (int j = 0; j < p; j++) {
      
      Dj = D.col(j);
      
      arma::vec Dz_notj = Dz - (Dj * z(j));
      
      arma::vec a = (lb - Dz_notj) / Dj;
      arma::vec b = (ub - Dz_notj) / Dj;
      
      double lower_j = arma::max(a);
      double upper_j = arma::min(b);
      
      double zj_old = z(j);
      z(j) = rtuvn_be2017(alpha(j), 1., lower_j, upper_j);
      
      Dz += Dj * (z(j) - zj_old);
    }
    
    samples.col(i) = z;
  }
  
  return(samples);
}