#' @export
sample_rsm = function(n, mu, Sigma, lb, ub, A = NULL) {
  
  d = length(mu)
  
  Prec = chol2inv(chol(Sigma))
  
  # find the mode of the truncated normal with quadratic programming
  if (is.null(A)) A = diag(d)
  
  # convert constraints to appropriate form for quadprog
  inf_idx = is.infinite(c(lb, ub))
  Amat = - t(rbind(A, -A)[!inf_idx, ])
  bvec = c(lb - A %*% mu, -ub + A %*% mu)[!inf_idx]
  dvec = rep(0, d)
  
  opt = quadprog::solve.QP(Prec, dvec, Amat, bvec)
  map = opt$solution
  
  Prec_mode = Prec %*% map
  
  if (!is.null(A)) {
    samples = rsm(n, map, Sigma, Prec_mode, lb - A %*% mu, ub - A %*% mu, A)  
  }
  else {
    samples = rsm_axis(n, map, Sigma, Prec_mode, lb - A %*% mu, ub - A %*% mu)
  }
  
  return(t(samples) + matrix(rep(mu, n), nrow = n, byrow = TRUE))
}

#' @param n number of random samples desired (sample size).
#' @param mu mean vector of the underlying multivariate normal distribution.
#' @param Sigma positive definite covariance matrix of the underlying multivariate normal distribution.
#' @param A matrix or vector of coefficients of linear inequality constraints.
#' @param lower vector of lower bounds for truncation.
#' @param upper vector of upper bounds for truncation.
#' @param init initial value vector for Gibbs sampler (satisfying truncation), if \code{NULL} then determine automatically.
#' @param burn burn-in iterations discarded (default as \code{10}).
#' @param thin thinning lag (default as \code{1}).
#' 
#' @export
sample_gibbs_lg2015 <- 
  function(n, mu, Sigma, lower, upper, A = diag(length(mu)),
           init=NULL, burn=10, thin=0, tuvn_sampler = "lg2015") {
  
  if (length(mu) == 1) {
    result <- rtuvn(n=n, mean=mu, sd=c(Sigma), lower=lower, upper=upper)
  } 
  else {
    
    if ( any(lower >= upper)) stop("lower bound must be smaller than upper bound\n")
    
    bound.check <- 0
    
    if (!is.null(init)) {
      # check initial value validity
      inits_test <- A %*% init
      lower.log <- inits_test >= lower + 1e-8  # small tol for get away from bound
      upper.log <- inits_test <= upper - 1e-8  # small tol for get away from bound
      bound.check <- prod(lower.log * upper.log)
      if (bound.check == 0) 
        cat("initial is outside or too close from boundary, will be auto-corrected by ginv()!\n")
    } else if (bound.check == 0) {
      A.inv <- MASS::ginv(A)
      init <- A.inv%*%(lower + upper)/2
    }
    
    if ( any( c(burn, thin, n) %% 1 != 0) )  
      stop("burn, thin and n must be integer\n")
    if ( any( c(burn, thin, n -1) < 0) ) 
      stop("burn, thin must be  non-negative interger, n must be positive integer\n")
    
    if (is.vector(A)) {
      Rtilde <- t(as.matrix(A))
      lower <- as.vector(lower)
      upper <- as.vector(upper)
    } else {
      Rtilde <- A
    }
    
    a <- lower - Rtilde%*%mu
    b <- upper - Rtilde%*%mu
    Sigma.chol <- t(chol(Sigma))
    R <- Rtilde %*% Sigma.chol
    
    p <- ncol(R) # number of parameters, i.e. length of beta vector
    #   m <- nrow(R) # number of constraints
    
    z <- solve(Sigma.chol, init - mu)  # int is the initial value for the original problem
    Rz <- R %*% z
    
    # draw standardized samples
    total_samples = (thin+1)*n + burn
    samples = lg2015(total_samples, z, R, Rz, a, b, tuvn_sampler)
    
    
    final.ind <- 1:total_samples
    final.ind <- final.ind[(burn+1):length(final.ind)]
    final.ind <- seq(1,length(final.ind),by=thin+1) + thin + burn
    
    # unstandardize
    samples <- 
      Sigma.chol %*% samples[, final.ind, drop = FALSE] + 
      matrix(rep(mu, n), ncol = n, byrow = TRUE)
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
sample_gibbs_ry2004 = function(n, mu, Sigma, lb, ub, A = diag(length(mu)),
                               init = NULL, burn = 10, thin = 0, 
                               tuvn_sampler = "lg2015") {
  
  if (is.null(init)) {
    A_inv <- MASS::ginv(A)
    init <- A_inv %*% (lb + ub) / 2  
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
  samples <- L %*% samples[, final.ind, drop = FALSE]
  
  return(t(samples))
}

#' @export
sample_rhmc = function(n, mu, Sigma, lb, ub, 
                       initial = NULL, burnin = 0, 
                       traj_length = 2, max_stepsize = .1,
                       max_relative_stepsize = .2, implicit_iter = 1) {
  
  if (is.null(initial)) {
    initial = (lb + ub) / 2 - mu
  }    
  
  Prec = chol2inv(chol(Sigma))
  R = chol(Prec)
  
  samples = rhmc(n, R, lb - mu, ub - mu, burnin, initial - mu, 
                 traj_length, max_stepsize, max_relative_stepsize, implicit_iter)
  
  # un-center samples
  samples = samples + mu
  
  return(t(samples))
}
                       

#' @export
sample_hamiltonian_zigzag = function(n, mu, Sigma, lb, ub, A, 
                                     intg_time, init = NULL, p_init = NULL) {
  d = length(mu)
    Prec = chol2inv(chol(Sigma))
  
  # p_init = rnorm(d)
  p_init = init
  
  samples = t(hzz(n, Prec, A, lb, ub, init, intg_time))
  return(samples)
}