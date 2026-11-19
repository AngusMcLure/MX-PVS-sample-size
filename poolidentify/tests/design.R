library(poolidentify)

# Extremisers with the same individual mean and variance occupy opposite
# ends of the bias range. Calculate the two count laws directly.
p <- 0.025
rho <- 0.15
s <- 30
b <- p + rho * (1 - p)
a <- p * (1 - rho)
limits <- naive_bias_bounds(p, rho, s)
naive <- function(x, w) -expm1(log(sum(w * (1 - x)^s)) / s) - p
expected <- c(
  naive(c(0, b), c(1 - p/b, p/b)),
  naive(c(a, 1), c((1-p)/(1-a), (p-a)/(1-a))))
stopifnot(
  isTRUE(all.equal(unname(limits), expected)),
  all(naive_bias_bounds(p, rho, 1) == 0),
  all(naive_bias_bounds(p, 0, s) == 0),
  diff(naive_bias_bounds(p, rho, 2)) < 1e-14)
