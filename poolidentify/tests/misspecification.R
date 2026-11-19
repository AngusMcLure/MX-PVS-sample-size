library(poolidentify)

# Verify count probabilities against independent integration on the original
# prevalence or logit scale, including a broad distribution at low prevalence.
s <- 12
v <- 4
beta_law <- mixing_counts(s, v, "beta", c(0.4, 1/6))
ref_beta <- vapply(0:v, function(k) integrate(function(x) {
  dbinom(k, v, -expm1(s * log1p(-x))) * dbeta(x, 2, 3)
}, 0, 1, rel.tol = 1e-11)$value, numeric(1))
stopifnot(max(abs(beta_law$probabilities - ref_beta)) < 1e-11)
ln_law <- mixing_counts(s, v, "logitnormal", c(-8, 3.2))
ref_ln <- vapply(0:v, function(k) integrate(function(z) {
  dbinom(k, v, -expm1(s * plogis(8 - 3.2 * z, log.p = TRUE))) * dnorm(z)
}, -10, 10, rel.tol = 1e-11)$value, numeric(1))
stopifnot(max(abs(ln_law$probabilities - ref_ln)) < 2e-10)
ref_mean <- integrate(function(z) plogis(-8 + 3.2 * z) * dnorm(z),
                      -10, 10, rel.tol = 1e-11)$value
stopifnot(abs(ln_law$prevalence - ref_mean) < 2e-10)

# Two observable moments can be matched by both families at different means.
b <- mixing_counts(25, 2, "beta", c(0.01, 0.1))
bf <- working_projection(b$probabilities, 25, "beta")
lf <- working_projection(b$probabilities, 25, "logitnormal")
stopifnot(max(abs(bf$par - c(0.01, 0.1))) < 2e-6,
          abs(lf$prevalence - 0.01645824) < 2e-6,
          max(abs(lf$probabilities - b$probabilities)) < 2e-7)

# Independent enumeration of two-location Wilson coverage and finite bias.
prob <- c(0.35, 0.30, 0.35)
design <- expand.grid(first = 0:2, second = 0:2)
weights <- prob[design$first + 1] * prob[design$second + 1]
qhat <- (design$first + design$second)/4
z <- qnorm(.975)
width <- z * sqrt((qhat * (1 - qhat) + z^2/16)/4)
middle <- qhat + z^2/8
lower <- 1 - sqrt(1 - (middle - width) / (1 + z^2/4))
upper <- 1 - sqrt(pmax(0, 1 - (middle + width) / (1 + z^2/4)))
# The lower principal representation has nodes 0 and .7. Select an
# interior compatible prevalence for the coverage comparison.
p <- 0.35
performance <- working_performance(prob, p, 2, "naive", 2)
stopifnot(abs(performance$coverage - sum(weights[lower <= p & p <= upper])) < 1e-12,
          abs(performance$bias - (sum(weights * (1 - sqrt(1 - qhat))) - p)) < 1e-12)
stopifnot(inherits(try(working_performance(b$probabilities, 0.8, 25, "beta", 3),
                       silent = TRUE), "try-error"))

# A population with extremely frequent positives still retains negative-cell
# information about its naive target.
p <- 0.85
tiny <- mixing_counts(25, 2, "beta", c(p, 0))
stopifnot(tiny$probabilities[2] > 0, 1 - tiny$probabilities[3] == 0,
          abs(working_projection(tiny$probabilities, 25, "naive")$prevalence - p) < 1e-12)

set.seed(173)
for (family in c("beta", "logitnormal")) {
  good <- working_performance(b$probabilities, b$prevalence, 25, family, 500, replicates = 3)
  stopifnot(is.finite(good$bias), good$coverage >= 0, good$coverage <= 1,
            good$failure_fraction == 0)
  boundary <- working_performance(c(1, 0, 0), 0, 25, family, 5, replicates = 2)
  stopifnot(boundary$bias == 0, boundary$coverage == 0, boundary$failure_fraction == 1)
}
