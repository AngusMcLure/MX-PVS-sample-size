library(poolidentify)

# The universal lower bound is attained with individual testing, irrespective
# of the sizes and number of the extra pools. Axes are sorted by pool size.
w <- 200
penalty <- sqrt(2 * log(40) / w)
h <- matrix(c(10, 15, 50, 5, 20, 100), nrow = 3)
single <- affine_pool_interval(c(8, 1, 1), w, histogram = h)
unpooled <- sum(rowSums(h) * c(1, .5, 0)) / w
stopifnot(all.equal(single$width_bound, penalty), single$structural_gap == 0,
          max(abs(single$interval - (unpooled + c(-.5, .5) * penalty))) < 1e-12,
          affine_pool_interval(c(8, 1), 1)$width_bound == 1)

# Independent dense evaluation uses the joint negative-count probabilities.
u <- seq(0, 1, length.out = 1001)
for (d in list(3, c(3, 7), c(7, 3, 3))) {
  fit <- affine_pool_interval(d, w, grid_size = 12)
  sizes <- sort(unique(d))
  v <- tabulate(match(d, sizes))
  e <- vapply(u, function(p) {
    cells <- Map(function(s, n) dbinom(0:n, n, (1-p)^s), sizes, v)
    law <- Reduce(function(a,b) as.vector(outer(a,b)), cells)
    sum(law * as.numeric(fit$score))
  }, numeric(1))
  stopifnot(min(u - e) >= fit$remainder[1] - 1e-12,
            max(u - e) <= fit$remainder[2] + 1e-12,
            fit$width_bound <= 1, fit$width_bound >= penalty - 1e-12)
}

# Analytic one-pool solution: independently optimise the score gap numerically.
s <- 5
expected <- optimize(function(score) {
  r <- u - score * (1 - (1-u)^s)
  diff(range(r)) + penalty * abs(score)
}, c(0, 2), tol = 1e-10)$objective
stopifnot(abs(affine_pool_interval(s, w)$width_bound - expected) < 1e-5)
stopifnot(inherits(try(affine_pool_interval(c(7, 3, 3), 200, histogram = t(h)),
                       silent = TRUE), "try-error"))
