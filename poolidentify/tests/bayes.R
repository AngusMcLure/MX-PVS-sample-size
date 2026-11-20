library(poolidentify)
set.seed(31603)

# With one pool there is a conjugate posterior even at an observed boundary.
fit <- bayes_prevalence(c(40, 0), 1, draws = 2500, chains = 2, within = c(1, 4))
stopifnot(abs(mean(fit$draws$c1) - 1 / 42) < 0.002,
          all(fit$draws$prevalence == fit$draws$c1),
          max(abs(fit$envelope - qbeta(c(.025, .975), 1, 41))) < .006)

# A tensor midpoint rule independently integrates the two-coordinate likelihood.
grid <- seq(.0025, .9975, by = .005)
g <- expand.grid(u = grid, r = grid)
var_q <- g$u * (1 - g$u) * g$r
mass <- ((1 - g$u)^2 + var_q)^24 * (2 * (g$u * (1 - g$u) - var_q))^4 *
  (g$u^2 + var_q)^2 * dbeta(g$r, 1, 3)
fit <- bayes_prevalence(c(24, 4, 2), 25, prior = "shrink", within = c(1, 4),
                        draws = 1500, burn = 100, chains = 3)
stopifnot(abs(mean(fit$draws$c1) - weighted.mean(g$u, mass)) < .009,
          abs(mean(fit$draws$c2) - weighted.mean(g$r, mass)) < .018,
          all(fit$draws$lower <= fit$draws$prevalence),
          all(fit$draws$prevalence <= fit$draws$upper))

# Three observed moments of uniform Q give equiprobable counts; Gauss and
# Simpson rules integrate the target at the lower and upper extremes.
z <- c(.5, 1/3, .5)
f <- function(x) -expm1(log1p(-x) / 8)
stopifnot(max(abs(poolidentify:::.pi_canonical_counts(z) - .25)) < 1e-13,
          max(abs(poolidentify:::.pi_canonical_ends(z, 8) -
                    c(mean(f(.5 + c(-1, 1) / sqrt(12))),
                      sum(c(1, 4, 1) * f(c(0, .5, 1))) / 6))) < 1e-12)

# 150-digit Jacobi quadrature reference. With nearly saturated pools the
# lower principal node of 1-Q is far below machine epsilon.
z <- c(1 - 1e-14, 1 - 1e-12, 1 - 1e-14, 1e-13)
stopifnot(max(abs(poolidentify:::.pi_canonical_ends(z, 50) -
                    c(.84151581989855329, .99999999999999001))) < 1e-12)
