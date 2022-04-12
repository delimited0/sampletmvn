#ifndef RHMC_H
#define RHMC_H

#include "hamiltonian.h"

class Integrator {
public:
  Integrator(double traj_length, double max_stepsize, double max_relative_stepsize, int implicit_iter) :
    traj_length(traj_length), max_stepsize(max_stepsize),
    max_relative_stepsize(max_relative_stepsize), 
    implicit_iter(implicit_iter) {
  }  
  
  double traj_length;
  double max_relative_stepsize;
  double max_stepsize;
  int implicit_iter;
  
  std::pair<arma::vec, arma::vec> implicit_midpoint(arma::vec x, arma::vec v, 
                                                    Hamiltonian & ham);
};

#endif