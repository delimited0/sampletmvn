#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

# include "barrier.h"

class Hamiltonian {
public: 
  int d;
  TwoSidedBarrier barrier;
  
  Hamiltonian(TwoSidedBarrier barrier) :
    d(barrier.lb.n_elem), barrier(barrier) {}
  
  virtual double f(const arma::vec & x) = 0;
  virtual arma::vec df(const arma::vec & x) = 0;
  double H(const arma::vec & x, const arma::vec & v);
  std::pair<arma::vec, arma::vec> dH(const arma::vec & x, const arma::vec & v);
  arma::vec generate(const arma::vec & x);
  double stepsize(const arma::vec & x, const arma::vec & dx);
};

class PrecHamiltonian : public Hamiltonian {
public:
  arma::mat R;  // cholesky factor of precision matrix
  
  PrecHamiltonian(arma::mat R, TwoSidedBarrier barrier) :
    Hamiltonian(barrier), R(R) { }
  
  double f(const arma::vec & x) override;
  arma::vec df(const arma::vec & x) override;
};

#endif