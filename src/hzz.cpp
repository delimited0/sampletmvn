#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

arma::vec zigzag(arma::vec x, arma::vec p, arma::mat & Prec, 
                 arma::mat & A, arma::vec & lb, arma::vec & ub, 
                 double T) {
  
    int d = x.n_elem;
    double tau = 0.0;
    arma::vec v = arma::sign(p);
    
    arma::vec Phix = Prec * x;
    arma::vec Phiv = Prec * v;
    
    int i = 0;
    while (tau < T) {
      
      arma::vec c = -p;
      
      // minimum positive root
      arma::vec t_grad(p);
      for (int j = 0; j < d; j++) {
        
        arma::vec coeffs = {.5*Phiv(j), Phix(j), c(j)};
        arma::cx_vec complex_roots = arma::roots(coeffs);
        arma::vec roots = real(complex_roots);
        arma::uvec ispositive = arma::find(roots > arma::datum::eps);
        
        if (ispositive.n_elem == 0) {
          t_grad(j) = arma::datum::inf;
        }
        else {
          t_grad(j) = arma::min(roots(ispositive));
        }
      }
      
      double t_min_grad = arma::min(t_grad);
      arma::uword idx_grad = t_grad.index_min();
      
      // reflect on boundary
      arma::vec t_bdry = -x / v;
      arma::vec At_bdry = A * t_bdry;
      
      arma::uvec out_of_bounds = arma::find(At_bdry < lb || At_bdry > ub);
      
      t_bdry.elem(out_of_bounds) = arma::datum::inf * arma::ones(out_of_bounds.n_elem);
      
      double t_min_bdry = arma::min(t_bdry);
      arma::uword idx_bdry = t_bdry.index_min();
      
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
      
      Rcpp::Rcout << "t_grad: " << t_grad << std::endl;
      Rcpp::Rcout << "t_bdry: " << t_bdry << std::endl;
      Rcpp::Rcout << "t_min: " << t_min << std::endl;
      
      if (tau + t_min > T) {
        // no further event occurred
        x += (T - tau) * v;
        p = p - (T - tau) * Phix - std::pow(T - tau, 2.) * Phiv / 2.;
        tau = T;
      }
      else {
        x += t_min * v;
        p = p - t_min * Phix - std::pow(t_min, 2.) * Phiv / 2.;
        // double v_min = -v(idx_min);
        v(idx_min) = -v(idx_min);
        Phix += t_min * Phiv;
        Phiv += 2. * v(idx_min) * Prec.col(idx_min);
        tau += t_min;
      }
      
      Rcpp::Rcout << "x: " << x << std::endl;
      i++;
    }
  
  return x;
}

// [[Rcpp::export]]
arma::mat hzz(int n, arma::mat Prec, 
              arma::mat A, arma::vec lb, arma::vec ub, 
              arma::vec init, double T) {
  
  int d = init.n_elem;
  arma::mat samples(d, n);
  
  arma::vec x = init;
  arma::vec p;
  
  for (int i = 0; i < n; i++) {
    p = arma::randg(d) - arma::randg(d);
    x = zigzag(x, p, Prec, A, lb, ub, T);
    samples.col(i) = x;
  }
  
  return(x);
}