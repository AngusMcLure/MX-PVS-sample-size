library(poolidentify)

# Compare screening and simulated-survey widths with the public inference
# entry points. The planning engine always takes negative-count histograms.
# This asymmetric, severely underdispersed sample needs the feasibility repair.
width <- poolidentify:::planning_width_engine(c(5, 5), 200, level = .91, grid_size = 16)
for (h in list(c(1, 198, 1), c(179, 18, 3), c(0, 0, 200))) {
  reference <- fit_prevalence(rev(h), 5, level = .91, grid_size = 16)
  stopifnot(abs(width(h) - diff(reference$confidence)) < 1e-9)
}
stopifnot(abs(width(c(1.00001, 197.99999, 1)) - width(c(1, 198, 1))) < 2e-6)

# An unsorted design with two sizes and multiplicities (2,1) ensures that
# reversing the count axes or ignoring one pool would change the answer.
width <- poolidentify:::planning_width_engine(c(9, 2, 2), 40, level = .9, grid_size = 14)
h <- array(c(0, 3, 1, 1, 8, 27), c(3, 2))
reference <- mixed_confidence(h, c(2, 9), level = .9, grid_size = 14)
stopifnot(abs(width(as.numeric(h)) - diff(reference$confidence)) < 1e-9)

# When every positive count is close to zero or one, integrate on the rare
# tail. Compare the truncated weighted sum to full enumeration, including
# both deterministic survey outcomes.
n <- 300
s <- 17
k <- 0:n
widths <- vapply(k, function(h) {
  ci <- stats::binom.test(h, n, conf.level = .9)$conf.int
  ci[2] + expm1(log1p(-ci[1]) / s)
}, numeric(1))
for (q in c(0, .004, .5, .996, 1)) {
  exact <- sum(stats::dbinom(k, n, q) * widths)
  near <- poolidentify:::planning_one_expected(c(q, 1-q), n, s, level = .9)
  stopifnot(abs(exact - near) < 8e-15)
}
