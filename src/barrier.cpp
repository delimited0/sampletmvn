#include "barrier.h"

bool TwoSidedBarrier::feasible(const arma::vec & x) {
  return arma::all((x > this->lb) && (x < this->ub));
}

arma::vec TwoSidedBarrier::hessian(const arma::vec & x) {
  arma::vec d = 
    1 / ((x - this->lb) % (x - this->lb)) +
    1 / ((this->ub - x) % (this->ub - x));
  return d;
}

arma::vec TwoSidedBarrier::tensor(const arma::vec & x) {
  arma::vec t = -2 * (1 / ((x - lb) % (x - lb) % (x - lb)) - 
    1 / ((ub - x) % (ub - x) % (ub - x)));
  return t;
}

arma::vec TwoSidedBarrier::quadratic_form_gradient(const arma::vec & x, const arma::mat & u) {
  arma::vec t = this->tensor(x);
  return arma::sum(arma::pow(u, 2), 1) % t;
}

arma::vec TwoSidedBarrier::log_det_gradient(const arma::vec & x) {
  arma::vec d = this->hessian(x);
  arma::vec t = this->tensor(x);
  return t / d;
}

double TwoSidedBarrier::stepsize(const arma::vec & x, const arma::vec & v) {
  double t;
  arma::uvec pi = arma::find(v > 0);
  if (pi.size() > 0)
    t = std::min(1e40, arma::min((this->ub(pi) - x(pi)) / v(pi)));
  else
    t = 1e40;
  arma::uvec ni = arma::find(v < 0);
  if (ni.size() > 0)
    t = std::min(t, arma::min((this->lb(ni) - x(ni)) / v(ni)));
  
  return(t);
}

arma::vec TwoSidedBarrier::hessian_norm(const arma::vec & x, const arma::vec & v) {
  arma::vec d = this->hessian(x);
  arma::vec u = (v % v) % d;
  return u;
}
