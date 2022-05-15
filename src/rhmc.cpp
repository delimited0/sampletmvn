// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "rhmc.h"

// integrator
std::pair<arma::vec, arma::vec> Integrator::implicit_midpoint(arma::vec x, arma::vec v, Hamiltonian & ham) {
  
  double total_h = 0.0;
  std::pair<arma::vec, arma::vec> dz = ham.dH(x, v);
  arma::vec dx = dz.first;
  arma::vec dv = dz.second;
  
  while (total_h < .9 * this->traj_length) {
  
    double h = 
      std::min(
        std::min(this->max_relative_stepsize * ham.stepsize(x, dx),
                 this->max_stepsize), 
        this->traj_length - total_h
      );
    
    for (int i = 0; i < this->implicit_iter; i++) {
      dz = ham.dH(x + .5 * h * dx, v + .5 * h * dv);
      dx = dz.first;
      dv = dz.second;
      h = std::min(h, .5 * ham.stepsize(x, dx));
    }
    
    x = x + h * dx;
    v = v + h * dv;
    total_h = total_h + h;
  }
  
  return std::pair<arma::vec, arma::vec>(x, v);
}

arma::mat rhmc_simulator(int n, int burnin, arma::vec initial, Hamiltonian & ham, Integrator intg) {
  arma::mat samples(ham.d, n);
  arma::vec x = initial;
  std::pair<arma::vec, arma::vec> z;
  arma::vec v;
  
  for (int i = 0; i < (n + burnin); i++) {
    
    v = ham.generate(x);
        // Rcpp::Rcout << "z shape: " << z.n_elem << std::endl;
    z = intg.implicit_midpoint(x, v, ham);
    
    x = z.first;
    
    if (i >= burnin) 
      samples.col(i-burnin) = x;
  }
  
  return samples;
}


//' @param R upper cholesky factor of precision
//' @param lb lower bound with mean subtracted
//' @param ub upper bound with mean subtracted
// [[Rcpp::export]]
arma::mat rhmc(int n, arma::mat R, arma::vec lb, arma::vec ub, 
               int burnin, arma::vec initial, 
               double traj_length, double max_stepsize, double max_relative_stepsize,
               int implicit_iter) {
  
  TwoSidedBarrier barrier(lb, ub);
  PrecHamiltonian ham(R, barrier);
  Integrator intg(traj_length, max_stepsize, max_relative_stepsize, implicit_iter);
  
  arma::mat samples = rhmc_simulator(n, burnin, initial, ham, intg);
  
  return samples;
}