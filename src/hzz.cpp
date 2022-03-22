#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat hamiltonian_zigzag(int n, arma::mat Prec, 
                             arma::mat A, arma::vec lb, arma::vec ub, 
                             arma::vec init, arma::vec p_init,
                             double T) {
  
  int d = init.n_elem;
  
  arma::vec x = init;
  arma::vec p = p_init;
  
  arma::mat samples(d, n);
  
  for (int i = 0; i < n; i++) {
    double tau = 0.0;
    arma::vec v = arma::sign(p);
    
    Rcpp::Rcout << "Pre-step state:" << std::endl;
    Rcpp::Rcout << "x: " << x << std::endl;
    Rcpp::Rcout << "v: " << v << std::endl;
    
    arma::vec Phix = Prec * x;
    arma::vec Phiv = Prec * x;
    
    int iter = 0;
    while (tau < T && iter < 300) {
      Rcpp::Rcout << "=======================" << std::endl;
      Rcpp::Rcout << "=======================" << std::endl;
      
      Rcpp::Rcout << "tau: " << tau << std::endl;
      
      arma::vec c = -p;
      
      // minimum positive root
      arma::vec t_grad(p);
      for (int j = 0; j < d; j++) {
        
        arma::vec coeffs = {.5*Phiv(j), Phix(j), c(j)};
        arma::cx_vec complex_roots = arma::roots(coeffs);
        arma::vec roots = real(complex_roots);
        arma::uvec ispositive = arma::find(roots > 0);
        
        if (ispositive.n_elem == 0) {
          t_grad(j) = arma::datum::inf;
        }
        else {
          t_grad(j) = arma::min(roots(ispositive));
        }
      }
      
      Rcpp::Rcout << "t_grad: " << t_grad << std::endl;
      
      double t_min_grad = arma::min(t_grad);
      arma::uword idx_grad = t_grad.index_min();
      
      Rcpp::Rcout << "tmin grad: " << t_min_grad << std::endl;
      
      // reflect on boundary
      arma::vec t_bdry = -x / v;
      arma::vec At_bdry = A * t_bdry;
      
      Rcpp::Rcout << "initial t_bdry: " << t_bdry << std::endl;
      
      arma::uvec out_of_bounds = arma::find(At_bdry < lb || At_bdry > ub);
      Rcpp::Rcout << "out_of_bounds: " << out_of_bounds << std::endl;
      
      t_bdry.elem(out_of_bounds) = arma::datum::inf * arma::ones(out_of_bounds.n_elem);
      
      Rcpp::Rcout << "t_bdry: " << t_bdry << std::endl;
      
      double t_min_bdry = arma::min(t_bdry);
      arma::uword idx_bdry = t_bdry.index_min();
      
      Rcpp::Rcout << "tmin bdry: " << t_min_bdry << std::endl;
      
      // next event time
      double t_min;
      arma::uword idx_min;
      if (t_min_grad < t_min_bdry) {
        t_min = t_min_grad;
        idx_min = idx_grad;
      }
      else {
        t_min = t_min_bdry;
        idx_min = idx_bdry;
      }
      
      Rcpp::Rcout << "t_min: " << t_min << std::endl;
      
      if (tau + t_min > T) {
        x += (T - tau) * v;
        p = p - (T - tau) * Phix - std::pow(T - tau, 2.) * Phiv / 2.;
        tau = T;
        
        Rcpp::Rcout << "Set tau to max time" << std::endl;
      }
      else {
        x += t_min * v;
        p = p - t_min * Phix - std::pow(t_min, 2.) * Phiv / 2.;
        double v_min = -v(idx_min);
        Phix += t_min * Phiv;
        Phiv += 2 * v_min * Prec.col(idx_min);
        tau += t_min;
        
        Rcpp::Rcout << "Added min to tau" << std::endl;
        
      }
      
      Rcpp::Rcout << "Updated state: " << std::endl;
      Rcpp::Rcout << "tau: " << tau << std::endl;
      Rcpp::Rcout << "p: " << p << std::endl;
      Rcpp::Rcout << "x: " << x << std::endl;  
      iter++;
    }
    
    samples.col(i) = x;
  }
  
  return(samples);
}