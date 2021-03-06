#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include "tuvn.h"

// [[Rcpp::export]]
arma::mat ry2004(int n, arma::vec alpha, arma::vec z, 
                                   arma::mat D, arma::vec Dz, 
                                   arma::vec lb, arma::vec ub,
                                   Rcpp::String sampler) {
  
  std::unique_ptr<Tuvn> tuvn;
  if (sampler == "be2017") {
    tuvn = std::unique_ptr<Tuvn>(new BE2017Tuvn);
  }
  else {
    tuvn = std::unique_ptr<Tuvn>(new LG2015Tuvn);
  }
  
  int p = z.n_elem;
  
  arma::mat samples(p, n);
  arma::vec Dj;
  
  double lower_pos;
  double upper_pos;
  double lower_neg;
  double upper_neg;
  
  for (int i = 0; i < n; i++) {
    
    for (int j = 0; j < p; j++) {
      
      Dj = D.col(j);
      
      arma::vec Dz_notj = Dz - (Dj * z(j));
      
      arma::vec a = (lb - Dz_notj);
      arma::vec b = (ub - Dz_notj);
      
      arma::uvec pos = arma::find(Dj > 0.);
      arma::uvec neg = arma::find(Dj < 0.);
      
      if (pos.n_elem == 0) {
        lower_pos = -arma::datum::inf;
        upper_pos = arma::datum::inf;
      }
      else {
        lower_pos = arma::max(a(pos) / Dj(pos));
        upper_pos = arma::min(b(pos) / Dj(pos));
      }
      
      if (neg.n_elem == 0) {
        upper_neg = arma::datum::inf;
        lower_neg = -arma::datum::inf;
      } 
      else {
        upper_neg = arma::min(a(neg) / Dj(neg));
        lower_neg = arma::max(b(neg) / Dj(neg));
      }
      
      double lower_j = std::max(lower_pos, lower_neg);
      double upper_j = std::min(upper_pos, upper_neg);
      
      double zj_old = z(j);
      z(j) = tuvn->rtuvn(alpha(j), 1., lower_j, upper_j);
      
      Dz += Dj * (z(j) - zj_old);
    }
    
    samples.col(i) = z;
  }
  
  return(samples);
}