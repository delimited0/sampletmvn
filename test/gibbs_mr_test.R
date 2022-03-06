n = 1000
mu = 0
Sigma = 1

lb = -1
ub = Inf
samples = sampletmvn::rtuvn(n, mu, Sigma, lb, ub)
hist(samples)


lb = -Inf
ub = 1
samples = sampletmvn::rtuvn(n, mu, Sigma, lb, ub)
hist(samples)

lb = 2
ub = 3
samples = sampletmvn::rtuvn(n, mu, Sigma, lb, ub)
hist(samples)

lb = -2
ub = -1
samples = sampletmvn::rtuvn(n, mu, Sigma, lb, ub)
hist(samples)

pkg_samples = tmvmixnorm::rtuvn(n = n, lower = lb, upper = ub)
hist(pkg_samples)

lb = -2
ub = 0
samples = sampletmvn::rtuvn(n, mu, Sigma, lb, ub)
hist(samples)

pkg_samples = tmvmixnorm::rtuvn(n = n, lower = lb, upper = ub)
hist(pkg_samples)



