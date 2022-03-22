

# box test ----------------------------------------------------------------
n = 10000

d = 2
mu = rep(0, d)
Sigma = .5 * diag(d) + .5 * rep(1, d) %*% t(rep(1, d))
lb = rep(-2, d)
ub = rep(2, d)
A = diag(d)

rsm_samples = sampletmvn::sample_rsm(n, mu, Sigma, lb, ub, A)
plot(rsm_samples)

gibbs_mr_samples = sampletmvn::sample_gibbs_lg2015(n, mu, Sigma, A, lb, ub, 
                                                   init = c(1, 1), thin = 0,
                                                   tuvn_sampler = "lg2015")
plot(gibbs_mr_samples)

gibbs_mr_samples = sampletmvn::sample_gibbs_lg2015(n, mu, Sigma, A, lb, ub, 
                                                   init = c(1, 1), thin = 0,
                                                   tuvn_sampler = "be2017")
plot(gibbs_mr_samples)

gibbs_ry2004_samples = sampletmvn::sample_gibbs_ry2004(n, mu, Sigma, A, lb, ub,
                                                       init = c(1, 1), thin = 0,
                                                       tuvn_sampler = "lg2015")
plot(gibbs_ry2004_samples)

liness_samples = lincongauss::rtmvn(n, mu, Sigma, lb, ub, x_init = c(0, 0))
plot(liness_samples)

epess_samples = epmgpr::rtmvn(n, mu, Sigma, lb, ub)
plot(epess_samples)

rhmc_samples = truncmvnorm::rhmc_rtmvn(n, mu, Sigma, lb, ub, initial = c(0, 0))
plot(rhmc_samples)

hzz_samples = sampletmvn::sample_hamiltonian_zigzag(1, mu, Sigma, lb, ub, A, 
                                                    intg_time = 2, init = c(1, 1))

# positive orthant test --------------------------------------------------------
n = 1000

d = 2
mu = rep(0, d)
Sigma = .5 * diag(d) + .5 * rep(1, d) %*% t(rep(1, d))
lb = rep(0, d)
ub = rep(Inf, d)
A = diag(d)

rsm_samples = sampletmvn::sample_rsm(n, mu, Sigma, lb, ub, A)
plot(rsm_samples)

gibbs_mr_samples = sampletmvn::sample_gibbs_lg2015(n, mu, Sigma, A, lb, ub, 
                                                   init = c(1, 1), thin = 0)
plot(gibbs_mr_samples)

gibbs_ry2004_samples = sampletmvn::sample_gibbs_ry2004(n, mu, Sigma, A, lb, ub,
                                                       init = c(1, 1), thin = 0,
                                                       tuvn_sampler = "be2017")
plot(gibbs_ry2004_samples)

hzz_samples = sampletmvn::sample_hamiltonian_zigzag(10, mu, Sigma, lb, ub, A, 
                                                    intg_time = 1, 
                                                    init = c(1, 1), p_init = c(.1, .1))

# trapezoid test ----------------------------------------------------------
n = 1000
d = 2
muf <- function(d) rep(0, d)
Sigmaf <- function(d) diag(d)
lbf <- function(d) rep(-Inf, 2*d)
ubf <- function(d) c(0, rep(2, 2*d-1))
Af <- function(d) {
  lower_bounds <- -diag(d)
  upper_bounds <- diag(d)
  upper_bounds[1, ] <- c(2, 1, rep(0, d-2))
  A <- rbind(upper_bounds, lower_bounds)
  return(A)
}

mu = muf(d)
Sigma = Sigmaf(d)
lb = lbf(d)
ub = ubf(d)
A = Af(d)

rsm_samples = sampletmvn::sample_rsm(n, mu, Sigma, lb, ub, A)
plot(rsm_samples)

gibbs_mr_samples = sampletmvn::sample_gibbs_mixed_rejection(n, mu, Sigma, A, lb, ub, 
                                                            init = c(1, 1), thin = 0)
plot(gibbs_mr_samples)

