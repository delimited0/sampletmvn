#ifndef BARRIER_H
#define BARRIER_H

#include <RcppArmadillo.h>

class TwoSidedBarrier {
public:
  arma::vec lb;
  arma::vec ub;
  
  TwoSidedBarrier(arma::vec lb, arma::vec ub) : lb(lb), ub(ub) {
  }
  
  bool feasible(const arma::vec & x);
  arma::vec hessian(const arma::vec & x);
  arma::vec tensor(const arma::vec & x);
  arma::vec log_det_gradient(const arma::vec & x);
  arma::vec quadratic_form_gradient(const arma::vec & x, const arma::mat & u);
  double stepsize(const arma::vec & x, const arma::vec & v);
  arma::vec hessian_norm(const arma::vec & x, const arma::vec & v);
};

#endif