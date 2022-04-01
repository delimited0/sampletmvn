#include "tuvn.h"
#include "RcppArmadillo.h"
#include <boost/math/distributions/normal.hpp>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]


/* Li and Ghosh 2015 */

const double a0 = 0.25696;

double unif_rej(double a, double b) {
  
  double u;
  double rho;
  double x;
  
  do {
    x = (b - a) * arma::randu() + a;
    u = arma::randu();
    
    if (a <= 0 && b >= 0) {
      rho = std::exp(-x*x / 2.);
    }
    if (a > 0) {
      rho = std::exp(-(x*x - a*a) / 2.);
    }
    if (b < 0) {
      rho = std::exp(-(x*x - b*b) / 2.);
    }
  } while (u > rho);
  
  return(x);
}

double norm_rej(double a, double b) {
  
  double x;
  
  do {
    x = arma::randn();
  } while (x < a || x > b);
  
  return(x);
}

double halfnorm_rej(double a, double b) {
  double x;
  do {
    x = arma::randn();
  } while (std::abs(x) < a || std::abs(x) > b);
  
  return(std::abs(x));
} 

double exp_rej(double a, double b) {
  
  double x;
  double u;
  double rho;
  double lambda = (a + std::sqrt(a*a + 4.) ) / 2.;
  
  std::random_device rd;
  std::mt19937 gen(rd());
  std::weibull_distribution<double> wd(1., 1. / lambda);
  
  do {
    x = wd(gen) + a;
    u = arma::randu();
    rho = std::exp(-std::pow((x - lambda), 2.) / 2);
  } while (u > rho || x >= b);
  
  return(x);
}

// [a, infty)
double rtuvn_case1(double a, double b) {
  double x;
  if (a < 0.) {
    x = norm_rej(a, b);
  }
  else {
    if (a < a0) {
      x = halfnorm_rej(a, b);
    }
    else {
      x = exp_rej(a, b);
    }
  }
  return(x);
}

// [a, b]
double rtuvn_case2(double a, double b) {
  double x;
  if (b > (std::sqrt(2*arma::datum::pi) + a)) {
    x = norm_rej(a, b);
  }
  else {
    x = unif_rej(a, b);
  }
  return(x);
}

// [a, b], a >= 0
double rtuvn_case3(double a, double b) {
  double x;
  if (a <= a0) {
    if (b <= (std::sqrt(arma::datum::pi / 2.) * std::exp(a*a / 2.) + a) ) {
      x = unif_rej(a, b);
    }
    else {
      x = halfnorm_rej(a, b);
    }
  }
  else {
    if (b <= (std::sqrt(2. * arma::datum::pi) + a) ) {
      x = unif_rej(a, b);
    }
    else {
      x = exp_rej(a, b);
    }
  }
  
  return(x);
}

/*
 * sample from univariate truncated normal
 * n number of samples
 * mean mean of the underlying univariate normal distribution.
 * sd standard deviation of the underlying univariate normal distribution.
 * a standardized lower bound
 * b standardized upper bound
 */
// [[Rcpp::export]]
arma::vec rtuvn(int n, double mean, double sd, double lower, double upper) {
  double x;
  
  double a = (lower - mean) / sd;
  double b = (upper - mean) / sd;
  
  arma::vec samples(n);
  
  for (int i = 0; i < n; i++) {
    
    if (std::isinf(a) || std::isinf(b)) {
      if (std::isinf(b)) {
        x = rtuvn_case1(a, b);
      }
      else {
        x = -rtuvn_case1(-b, -a);    
      }
    }
    else {
      if (a < 0 && b > 0) {
        x = rtuvn_case2(a, b);
      }
      if (a >= 0) {
        x = rtuvn_case3(a, b);
      }
      if (b <= 0) {
        x = -rtuvn_case3(-b, -a);
      }  
    }
    
    samples(i) = x;
  }
  
  samples = sd * samples + mean;
  
  return(samples);
}

double rtuvn_lg2015(double mean, double sd, double lower, double upper) {
  
  double x;
  
  double a = (lower - mean) / sd;
  double b = (upper - mean) / sd;
  
  if (std::isinf(a) || std::isinf(b)) {
    if (std::isinf(b)) {
      x = rtuvn_case1(a, b);
    }
    else {
      x = -rtuvn_case1(-b, -a);    
    }
  }
  else {
    if (a < 0 && b > 0) {
      x = rtuvn_case2(a, b);
    }
    if (a >= 0) {
      x = rtuvn_case3(a, b);
    }
    if (b <= 0) {
      x = -rtuvn_case3(-b, -a);
    }  
  }
  
  return sd * x + mean;
}

double LG2015Tuvn::rtuvn(const double mean, const double sd, 
                         const double lower, const double upper) {
  double x;
  
  double a = (lower - mean) / sd;
  double b = (upper - mean) / sd;
  
  if (std::isinf(a) || std::isinf(b)) {
    if (std::isinf(b)) {
      x = rtuvn_case1(a, b);
    }
    else {
      x = -rtuvn_case1(-b, -a);    
    }
  }
  else {
    if (a < 0 && b > 0) {
      x = rtuvn_case2(a, b);
    }
    if (a >= 0) {
      x = rtuvn_case3(a, b);
    }
    if (b <= 0) {
      x = -rtuvn_case3(-b, -a);
    }  
  }
  
  return sd * x + mean;
}

/* Botev and L'ecuyer  */
const double INV_THRESHOLD = 0.4;
const double TOL = 2.05;

double ntail(const double l, const double u) {
  double c = std::pow(l, 2.) / 2.;
  double q = std::expm1(c - std::pow(u, 2.) / 2.);
  double x;
  do {
    x = c - std::log(1 + arma::randu<double>() * q);
  } while (std::pow(arma::randu<double>(), 2.) * x > c);
  
  return std::sqrt(2. * x);
}

double trnd(const double l, const double u) {
  double x;
  do {
    x = arma::randn<double>();
  } while (x < l || x > u);
  
  return x;
}

double trninv(const double l, const double u) {
  boost::math::normal dist(0.0, 1.0);
  double pl = .5 * std::erfc(-l / std::sqrt(2.));
  double pu = .5 * std::erfc(-u / std::sqrt(2.));
  double x = quantile(dist, pl + (pu - pl) * arma::randu<double>());
  
  return(x);
}

double rtuvn_be2017(double mean, double sd, double lower, double upper) {
  double x;
  
  double l = (lower - mean) / sd;
  double u = (upper - mean) / sd;
  
  // case 1: a < l < u
  if (l > INV_THRESHOLD) {
    x = ntail(l, u);
  }
  
  // case 2: l < u < -a
  else if (u < -INV_THRESHOLD) {
    x = -ntail(-u, -l);
  }
  
  else {
    // case 3: abs(u - l) > tol
    if (std::abs(u - l) > TOL) {
      x = trnd(l, u);
    }
    // case 4: abs(u - l) < tol
    else {
      x = trninv(l, u);
    }
  }
  
  return(x);
}

double BE2017Tuvn::rtuvn(const double mean, const double sd,
                         const double lower, const double upper) {
  double x;

  double l = (lower - mean) / sd;
  double u = (upper - mean) / sd;

  // case 1: a < l < u
  if (l > INV_THRESHOLD) {
    x = ntail(l, u);
  }

  // case 2: l < u < -a
  else if (u < -INV_THRESHOLD) {
    x = -ntail(-u, -l);
  }

  else {
    // case 3: abs(u - l) > tol
    if (std::abs(u - l) > TOL) {
      x = trnd(l, u);
    }
    // case 4: abs(u - l) < tol
    else {
      x = trninv(l, u);
    }
  }

  return(x);
}
