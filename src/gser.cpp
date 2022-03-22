#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include "tuvn.h"

// [[Rcpp::export]]
arma::mat lg2015_gibbs_iter_lg2015(int n, arma::vec z, 
                                arma::mat R, arma::vec Rz, 
                                arma::vec a, arma::vec b) {
  
  int p = R.n_cols;
  
  double lower_pos;
  double upper_pos;
  double lower_neg;
  double upper_neg;
  
  arma::mat samples(p, n);
  arma::vec Rj;
  
  for (int i = 0; i < n; i++) {
    
    for (int j = 0; j < p; j++) {
      
      Rj = R.col(j);
      
      arma::vec Rz_notj = Rz - ( Rj * z(j) );
      
      arma::vec a_temp = a - Rz_notj;
      arma::vec b_temp = b - Rz_notj;
      
      arma::uvec pos = arma::find(Rj > 0);
      arma::uvec neg = arma::find(Rj < 0);
      
      if (pos.n_elem == 0) {
        lower_pos = -arma::datum::inf;
        upper_pos = arma::datum::inf;
      }
      else {
        lower_pos = arma::max(a_temp(pos) / Rj(pos));
        upper_pos = arma::min(b_temp(pos) / Rj(pos));
      }
      
      if (neg.n_elem == 0) {
        upper_neg = arma::datum::inf;
        lower_neg = -arma::datum::inf;
      } 
      else {
        upper_neg = arma::min(a_temp(neg) / Rj(neg));
        lower_neg = arma::max(b_temp(neg) / Rj(neg));
      }
      
      double lower_j = std::max(lower_pos, lower_neg);
      double upper_j = std::min(upper_pos, upper_neg);
      
      double zj_old = z(j);
      z(j) = rtuvn_lg2015(0., 1., lower_j, upper_j);
      
      Rz += R.col(j) * (z(j) - zj_old);
    }
   
    samples.col(i) = z;
  }
  
  return(samples);
}

// [[Rcpp::export]]
arma::mat lg2015_gibbs_iter_be2017(int n, arma::vec z, 
                                   arma::mat R, arma::vec Rz, 
                                   arma::vec a, arma::vec b) {
  
  int p = R.n_cols;
  
  double lower_pos;
  double upper_pos;
  double lower_neg;
  double upper_neg;
  
  arma::mat samples(p, n);
  arma::vec Rj;
  
  for (int i = 0; i < n; i++) {
    
    for (int j = 0; j < p; j++) {
      
      Rj = R.col(j);
      
      arma::vec Rz_notj = Rz - ( Rj * z(j) );
      
      arma::vec a_temp = a - Rz_notj;
      arma::vec b_temp = b - Rz_notj;
      
      arma::uvec pos = arma::find(Rj > 0);
      arma::uvec neg = arma::find(Rj < 0);
      
      if (pos.n_elem == 0) {
        lower_pos = -arma::datum::inf;
        upper_pos = arma::datum::inf;
      }
      else {
        lower_pos = arma::max(a_temp(pos) / Rj(pos));
        upper_pos = arma::min(b_temp(pos) / Rj(pos));
      }
      
      if (neg.n_elem == 0) {
        upper_neg = arma::datum::inf;
        lower_neg = -arma::datum::inf;
      } 
      else {
        upper_neg = arma::min(a_temp(neg) / Rj(neg));
        lower_neg = arma::max(b_temp(neg) / Rj(neg));
      }
      
      double lower_j = std::max(lower_pos, lower_neg);
      double upper_j = std::min(upper_pos, upper_neg);
      
      double zj_old = z(j);
      z(j) = rtuvn_be2017(0., 1., lower_j, upper_j);
      
      Rz += R.col(j) * (z(j) - zj_old);
    }
    
    samples.col(i) = z;
  }
  
  return(samples);
}
