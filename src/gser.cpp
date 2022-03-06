#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

/*
 *  univariate rejection samplers
 */

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

// [[Rcpp::export]]
arma::mat gibbs_mixed_rejection(int n, arma::vec z, 
                                arma::mat R, arma::vec Rz, arma::vec a, arma::vec b) {
  
  int p = R.n_cols;
  
  double lower_pos;
  double upper_pos;
  double lower_neg;
  double upper_neg;
  
  arma::mat samples(p, n);
  arma::vec Rj;
  
  for (int i = 0; i < n; i++) {
    
    for (int j = 0; j < p; j++) {
      
      Rj = R.col(j);
      
      arma::vec Rz_notj = Rz - ( Rj * z(j) );
      
      arma::vec a_temp = a - Rz_notj;
      arma::vec b_temp = b - Rz_notj;
      
      arma::uvec pos = arma::find(Rj > 0);
      arma::uvec neg = arma::find(Rj < 0);
      
      if (pos.n_elem == 0) {
        lower_pos = -arma::datum::inf;
        upper_pos = arma::datum::inf;
      }
      else {
        lower_pos = arma::max(a_temp(pos) / Rj(pos));
        upper_pos = arma::min(b_temp(pos) / Rj(pos));
      }
      
      if (neg.n_elem == 0) {
        upper_neg = arma::datum::inf;
        lower_neg = -arma::datum::inf;
      } 
      else {
        upper_neg = arma::min(a_temp(neg) / Rj(neg));
        lower_neg = arma::max(b_temp(neg) / Rj(neg));
      }
      
      double lower_j = std::max(lower_pos, lower_neg);
      double upper_j = std::min(upper_pos, upper_neg);
      
      double zj_old = z(j);
      z(j) = rtuvn(1, 0., 1., lower_j, upper_j)(0);
      
      Rz += R.col(j) * (z(j) - zj_old);
    }
   
    samples.col(i) = z;
  }
  
  return(samples);
}


