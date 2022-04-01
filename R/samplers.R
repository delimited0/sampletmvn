#' @export
sample_rsm = function(n, mu, Sigma, lb, ub, A = NULL) {
  
  d = length(mu)
  Prec = chol2inv(chol(Sigma))
  
  # find the mode of the truncated normal with quadratic programming
  if (is.null(A)) {
    Amat = diag(d)
  }
  else {
    Amat = A
  }
  
  # convert constraints to appropriate form for quadprog
  inf_idx = is.infinite(c(lb, ub))
  Amat = - t(rbind(Amat, -Amat)[!inf_idx, ])
  bvec = c(lb, -ub)[!inf_idx]
  dvec = rep(0, d)
  
  opt = quadprog::solve.QP(Prec, dvec, Amat, bvec)
  map = opt$solution
  
  Prec_mode = Prec %*% map
  
  if (!is.null(A)) {
    samples = rsm(n, map, Sigma, Prec_mode, lb, ub, A)  
  }
  else {
    samples = rsm_axis(n, map, Sigma, Prec_mode, lb, ub)
  }
  
  return(t(samples))
}

#' @param n number of random samples desired (sample size).
#' @param Mean mean vector of the underlying multivariate normal distribution.
#' @param Sigma positive definite covariance matrix of the underlying multivariate normal distribution.
#' @param D matrix or vector of coefficients of linear inequality constraints.
#' @param lower vector of lower bounds for truncation.
#' @param upper vector of upper bounds for truncation.
#' @param init initial value vector for Gibbs sampler (satisfying truncation), if \code{NULL} then determine automatically.
#' @param burn burn-in iterations discarded (default as \code{10}).
#' @param thin thinning lag (default as \code{1}).
#' 
#' @export
sample_gibbs_lg2015 <- 
  function(n, Mean, Sigma, lower, upper, D = diag(1, length(Mean)),
           init=NULL, burn=10, thin=1, tuvn_sampler = "lg2015") {
  
  if (length(Mean) == 1) {
    result <- rtuvn(n=n, mean=Mean, sd=c(Sigma), lower=lower, upper=upper)
  } 
  else {
    
    if ( any(lower >= upper)) stop("lower bound must be smaller than upper bound\n")
    
    bound.check <- 0
    
    if (!is.null(init)) {
      inits_test <- D %*% init
      lower.log <- inits_test >= lower + 1e-8  # small tol for get away from bound
      upper.log <- inits_test <= upper - 1e-8  # small tol for get away from bound
      bound.check <- prod(lower.log * upper.log)
      if (bound.check == 0) cat("initial is outside or too close from boundary, will be auto-corrected by ginv()!\n")
    } else if (bound.check == 0) {
      D.inv <- MASS::ginv(D)
      init <- D.inv%*%(lower + upper)/2
    }
    
    if ( any( c(burn, thin, n) %% 1 != 0) )  stop("burn, thin and n must be integer\n")
    if ( any( c(burn, thin, n -1) < 0) ) stop("burn, thin must be  non-negative interger, n must be positive integer\n")
    
    if (is.vector(D)) {
      Rtilde <- t(as.matrix(D))
      lower <- as.vector(lower)
      upper <- as.vector(upper)
    } else {
      Rtilde <- D
    }
    
    a <- lower - Rtilde%*%Mean
    b <- upper - Rtilde%*%Mean
    Sigma.chol <- t(chol(Sigma))
    R <- Rtilde %*% Sigma.chol
    
    p <- ncol(R) # number of parameters, i.e. length of beta vector
    #   m <- nrow(R) # number of constraints
    
    z <- solve(Sigma.chol, init - Mean)  # int is the initial value for the original problem
    Rz <- R %*% z
    
    # draw standardized samples
    total_samples = (thin+1)*n + burn
    samples = lg2015(total_samples, z, R, Rz, a, b, tuvn_sampler)
    
    
    final.ind <- 1:total_samples
    final.ind <- final.ind[(burn+1):length(final.ind)]
    final.ind <- seq(1,length(final.ind),by=thin+1) + thin + burn
    
    # unstandardize
    samples <- Sigma.chol %*% samples[, final.ind] + Mean
  }
  
  return(t(samples))
}

#' @param n number of random samples desired (sample size).
#' @param mu mean vector of the underlying multivariate normal distribution.
#' @param Sigma positive definite covariance matrix of the underlying multivariate normal distribution.
#' @param lb vector of lower bounds for truncation.
#' @param ub vector of upper bounds for truncation.
#' @param A matrix or vector of coefficients of linear inequality constraints.
#' @param init initial value vector for Gibbs sampler (satisfying truncation), if \code{NULL} then determine automatically.
#' @param burn burn-in iterations discarded (default as \code{10}).
#' @param thin thinning lag (default as \code{1}).
#' @export
sample_gibbs_ry2004 = function(n, mu, Sigma, lb, ub, A = diag(1, length(mu)),
                               init = NULL, burn = 10, thin = 1, 
                               tuvn_sampler = "lg2015") {

  if (is.null(init)) {
    A_inv <- MASS::ginv(A)
    init <- D_inv %*% (lb + ub) / 2  
  }
  
  L = t(chol(Sigma))
  D = A %*% L
  alpha = solve(L, mu)
  
  z = solve(L, init)
  Dz = D %*% z
  
  total_samples = (thin+1)*n + burn
  samples = ry2004(total_samples, alpha, z, D, Dz, lb, ub, tuvn_sampler)
  
  # if (tuvn_sampler == "lg2015")
  #   samples = ry2004_gibbs_iter_lg2015(total_samples, alpha, z, D, Dz, lb, ub)
  # else if (tuvn_sampler == "be2017")
  #   samples = ry2004_gibbs_iter_be2017(total_samples, alpha, z, D, Dz, lb, ub)
  # else 
  #   stop("Invalid truncated univariate normal sampler")
  
  final.ind <- 1:total_samples
  final.ind <- final.ind[(burn+1):length(final.ind)]
  final.ind <- seq(1, length(final.ind), by=thin+1) + thin + burn
  
  # unstandardize
  samples <- L %*% samples[, final.ind] + mu
  
  return(t(samples))
}

#' @export
sample_hamiltonian_zigzag = function(n, mu, Sigma, lb, ub, A, 
                                     intg_time, init = NULL, p_init = NULL) {
  d = length(mu)
  Prec = chol2inv(chol(Sigma))
  
  # p_init = rnorm(d)
  p_init = init
  
  samples = t(hamiltonian_zigzag(n, Prec, A, lb, ub, init, p_init, intg_time))
  return(samples)
}