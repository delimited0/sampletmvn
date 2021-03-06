

# box test ----------------------------------------------------------------
n = 3000

d = 2
mu = rep(1, d)
Sigma = .5 * diag(d) + .5 * rep(1, d) %*% t(rep(1, d))
lb = rep(-2, d)
ub = rep(2, d)
A = diag(d)

par(mfrow = c(2, 4))

rsm_samples = sampletmvn::sample_rsm(n, mu, Sigma, lb, ub)
plot(rsm_samples, main = "RSM", pch = '+')

lg2015_samples = sampletmvn::sample_gibbs_lg2015(n, mu, Sigma,lb, ub, 
                                                   init = c(1, 1), thin = 0,
                                                   tuvn_sampler = "be2017")
plot(lg2015_samples, main = "LG2015", pch = '+')

ry2004_samples = sampletmvn::sample_gibbs_ry2004(n, mu, Sigma, lb, ub, A,
                                                 init = c(1, 1), thin = 0,
                                                 tuvn_sampler = "lg2015")
plot(ry2004_samples, main = "RY2004", pch = '+')

liness_samples = lincongauss::rtmvn(n, mu, Sigma, lb, ub, x_init = c(1, 1))
plot(liness_samples, main = "LINESS", pch = '+')

epess_samples = epmgpr::rtmvn(n, mu, Sigma, lb, ub, initial = c(1, 1))
plot(epess_samples, main = "EPESS", pch = "+")

met_samples = met::mvrandn(lb, ub, Sigma, n, mu)
plot(t(met_samples), main = "MET", pch = "+")

rhmc_samples = sampletmvn::sample_rhmc(n, mu, Sigma, lb, ub, initial = c(1, 1))
plot(rhmc_samples, main = "RHMC", pch = "+")

ghk_samples = tmvnsim::tmvnsim(n, d, lower = lb, upper = ub, means = mu, sigma = Sigma)
plot(ghk_samples$samp, main = "GHK", pch = "+")

hzz_samples = sampletmvn::sample_hamiltonian_zigzag(1, mu, Sigma, lb, ub, A, 
                                                    intg_time = 2, init = c(1, 1))


# positive orthant test --------------------------------------------------------
n = 1000

d = 2
mu = rep(2, d)
Sigma = .5 * diag(d) + .5 * rep(1, d) %*% t(rep(1, d))
lb = rep(0, d)
ub = rep(Inf, d)
A = diag(d)

par(mfrow = c(2, 3))

rsm_samples = sampletmvn::sample_rsm(n, mu, Sigma, lb, ub)
plot(rsm_samples, main = "RSM")

lg2015_samples = sampletmvn::sample_gibbs_lg2015(n, mu, Sigma, 
                                                 lower = lb, upper = ub, 
                                                 init = c(1, 1), thin = 0,
                                                 tuvn_sampler = "be2017")
plot(lg2015_samples, main = "LG2015")

ry2004_samples = sampletmvn::sample_gibbs_ry2004(n, mu, Sigma, 
                                                 lb, ub,
                                                 init = c(1, 1), thin = 0,
                                                 tuvn_sampler = "lg2015")
plot(ry2004_samples, main = "RY2004")

liness_samples = lincongauss::rtmvn(n, mu, Sigma, lb, ub, x_init = c(1, 1))
plot(liness_samples, main = "LINESS")

epess_samples = epmgpr::rtmvn(n, mu, Sigma, lb, ub, initial = c(1, 1))
plot(epess_samples, main = "EPESS")

met_samples = met::mvrandn(lb, ub, Sigma, n, mu)
plot(t(met_samples), main = "MET")

rhmc_samples = sampletmvn::sample_rhmc(n, mu, Sigma, lb, ub, initial = c(.1, .1))
plot(rhmc_samples, main = "RHMC")

hzz_samples = sampletmvn::sample_hamiltonian_zigzag(1, mu, Sigma, lb, ub, A, 
                                                    intg_time = 2, 
                                                    init = c(1, 1), p_init = c(.1, .1))

# trapezoid test ----------------------------------------------------------
n = 1000
d = 2
muf <- function(d) rep(1, d)
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

rsm_samples = sampletmvn::sample_rsm(n, mu, Sigma = Sigma, lb = lb, ub = ub,
                                     A = A)
plot(rsm_samples, main = "RSM")

lg2015_samples = sampletmvn::sample_gibbs_lg2015(n, mu, Sigma, lb, ub, A,
                                                 init = c(-1, -1), burn = 0,
                                                 tuvn_sampler = "be2017")
plot(lg2015_samples, main = "LG2015")

ry2004_samples = sampletmvn::sample_gibbs_ry2004(n, mu, Sigma, lb, ub, A, 
                                                 init = c(-1, -1), thin = 0,
                                                 tuvn_sampler = "lg2015")
plot(ry2004_samples, main = "RY2004")

liness_samples = lincongauss::rtmvn(n, mu, Sigma, lb, ub, A = A,
                                    x_init = c(1, 1))
plot(liness_samples, main = "LINESS")

epess_samples = epmgpr::rtmvn(n, mu, Sigma, lb, ub, A = A, initial = c(1, 1))
plot(epess_samples, main = "EPESS")


# Multinomial probit test -------------------------------------------------

max_polytope = function(y) {
  n = length(y)
  idx = which(y == 1)
  
  D = -1*diag(n-1)
  if (idx == 1)
    A = cbind(1, D)
  else if (idx == n)
    A = cbind(D, 1)
  else
    A = cbind(D[, 1:(idx-1)], rep(1, n-1), D[, idx:(n-1)])
  
  return(A)
}

y = c(0, 1, 0)
A = max_polytope(y)

lb = rep(0, 2)
ub = rep(Inf, 2)

n = 1000

mu = rep(0, 3)
Sigma = diag(3)

samples = sampletmvn::sample_gibbs_lg2015(n, mu, Sigma, lb, ub, A, tuvn_sampler = "lg2015", 
                                          init = 2 * y - 1)

sum(apply(samples, 1, function(s) {
  As = A %*% s 
  As <= ub & As >= lb
}))






