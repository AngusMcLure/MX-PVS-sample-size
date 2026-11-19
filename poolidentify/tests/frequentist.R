library(poolidentify)

# One binary observation per location gives the exact projected binomial interval.
one <- fit_prevalence(c(2, 8), 4, level = 0.9)
cp <- binom.test(8, 10, conf.level = 0.9)$conf.int
stopifnot(isTRUE(all.equal(unname(one$confidence), c(1 - (1 - cp[1])^(1/4), cp[2]))))

# K=1 with probability one is outside the two-pool binomial-mixture hull.
# Its nearest feasible count law is Bin(2, 1/2) in CDF distance.
bad <- fit_prevalence(c(0, 200, 0), 5, grid_size = 12)
stopifnot(max(abs(bad$fitted_pmf - c(1, 2, 1)/4)) < 1e-7,
          abs(bad$distance - 0.25) < 1e-7, diff(bad$identified) < 1e-7,
          bad$confidence[1] <= bad$identified[1], bad$confidence[2] >= bad$identified[2])
for (h in list(c(8, 0, 0), c(0, 0, 8))) {
  edge <- fit_prevalence(h, 10, grid_size = 10)
  truth <- as.numeric(h[3] > 0)
  stopifnot(all(abs(edge$identified - truth) < 1e-8),
            edge$confidence[1] <= truth, edge$confidence[2] >= truth)
}

# Exhaust all 45 histograms at eight locations. The true Q support is off-grid.
# Check the coverage of both sharp endpoints together, even on a coarse grid.
n <- 8
level <- 0.80
pi_true <- 0.76 * dbinom(0:2, 2, 0.107) + 0.24 * dbinom(0:2, 2, 0.893)
true_set <- identification_interval(pool_size = 5, probabilities = pi_true)
coverage <- 0
for (i in 0:n) {
  for (j in 0:(n-i)) {
    h <- c(i, j, n-i-j)
    fit <- fit_prevalence(h, 5, level = level, grid_size = 4)
    stopifnot(fit$confidence[1] <= fit$identified[1] + 1e-8,
              fit$confidence[2] >= fit$identified[2] - 1e-8)
    covered <- fit$confidence[1] <= true_set[1] && fit$confidence[2] >= true_set[2]
    coverage <- coverage + dmultinom(h, prob = pi_true) * covered
  }
}
stopifnot(coverage >= level - 1e-12)
