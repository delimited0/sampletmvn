#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat rsm(int n_samples, arma::vec mode, 
              arma::mat Sigma,
              arma::vec Prec_mode,
              arma::vec lb, arma::vec ub, arma::mat A) {
  
  arma::mat samples = arma::zeros(mode.n_elem, n_samples);
  
  for (int i = 0; i < n_samples; i++) {
    
    bool accepted = false;
    arma::vec proposal;
    
    while (!accepted) {
      
      proposal = arma::mvnrnd(mode, Sigma, 1);
      
      while (arma::any(A * proposal > ub) || arma::any(A * proposal < lb)) {
       
        proposal = arma::mvnrnd(mode, Sigma, 1);
        
      }
      
      double t = arma::as_scalar( 
        arma::exp( (mode - proposal).t() * Prec_mode ) 
      );
      
      double u = arma::randu();
      if (u <= t) accepted = true;
    }
    
    samples.col(i) = proposal;
  }
  
  return(samples);
}

// // [[Rcpp::export]]
// arma::mat rsm_axis(int n_samples, arma::vec mode, 
//                    arma::mat Sigma,
//                    arma::vec Prec_mode,
//                    arma::vec lb, arma::vec ub) {
//   
//   arma::mat samples = arma::zeros(mode.n_elem, n_samples);
//   
//   for (int i = 0; i < n_samples; i++) {
//     
//     bool accepted = false;
//     arma::vec proposal;
//     
//     while (!accepted) {
//       
//       proposal = arma::mvnrnd(mode, Sigma, 1);
//       
//       while (arma::any(proposal > ub) || arma::any(proposal < lb)) {
//         
//         proposal = arma::mvnrnd(mode, Sigma, 1);
//         
//       }
//       
//       double t = arma::as_scalar( 
//         arma::exp( (mode - proposal).t() * Prec_mode ) 
//       );
//       
//       double u = arma::randu();
//       if (u <= t) accepted = true;
//     }
//     
//     samples.col(i) = proposal;
//   }
//   
//   return(samples);
// }