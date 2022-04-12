#include "hamiltonian.h"

// energy function
double Hamiltonian::H(const arma::vec & x, const arma::vec & v) {
  
  arma::vec H = this->barrier.hessian(x);
  arma::vec hess_inv = 1 / H;
  arma::vec g_inv_v = hess_inv % v;
  double e = .5 * arma::as_scalar(v.t() * g_inv_v);
  
  e = e + arma::sum(arma::log(H)) * .5 + this->f(x);
  return e;
}

// Hamiltonian dynamics
std::pair<arma::vec, arma::vec> Hamiltonian::dH(const arma::vec & x, const arma::vec & v) {
  
  arma::vec g_inv_v = (1 / this->barrier.hessian(x)) % v;
  
  arma::vec sigma = this->barrier.log_det_gradient(x);
  
  arma::vec dx = g_inv_v;
  arma::vec dv = 
    -this->df(x) + .5 * this->barrier.quadratic_form_gradient(x, dx) - .5 * sigma;
    
  std::pair<arma::vec, arma::vec> result(dx, dv);
  
  return result;
}

// sample momentum
arma::vec Hamiltonian::generate(const arma::vec & x) {
  
  arma::vec gv = arma::sqrt(this->barrier.hessian(x)) % arma::randn(this->d);
  
  return gv;
}

double Hamiltonian::stepsize(const arma::vec & x, const arma::vec & dx) {
  
  double t1 = this->barrier.stepsize(x, dx);
  double t2 = 1 / arma::max(arma::sqrt(this->barrier.hessian_norm(x, dx)));
  
  return std::min(t1, arma::min(t2));
}

double PrecHamiltonian::f(const arma::vec & x) {
  
  arma::vec Rx = this->R * x;
  
  return .5 * arma::dot(Rx, Rx);
}

arma::vec PrecHamiltonian::df(const arma::vec & x) {
  
   return this->R.t() * this->R * x;
}
