#' @export
sample_rsm = function(n, mu, Sigma, lb, ub, A = NULL) {
  
  d = length(mu)
  Prec = chol2inv(chol(Sigma))
  
  inf_idx = is.infinite(c(lb, ub))
  Amat = - t(rbind(A, -A)[!inf_idx, ])
  bvec = c(lb, -ub)[!inf_idx]
  dvec = rep(0, d)
  
  opt = quadprog::solve.QP(Prec, dvec, Amat, bvec)
  map = opt$solution
  
  Prec_mode = Prec %*% map
  
  samples = rsm(n, map, Sigma, Prec_mode, lb, ub, A)
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
sample_gibbs_mixed_rejection <- 
  function(n, Mean, Sigma, D = diag(1, length(Mean)), lower, upper, 
           init=NULL, burn=10, thin=1) {
  
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
    samples = gibbs_mixed_rejection(total_samples, z, R, Rz, a, b)
    
    final.ind <- 1:total_samples
    final.ind <- final.ind[(burn+1):length(final.ind)]
    final.ind <- seq(1,length(final.ind),by=thin+1) + thin + burn
    
    samples <- Sigma.chol %*% samples[, final.ind] + Mean
  }
  
  return(t(samples))
}