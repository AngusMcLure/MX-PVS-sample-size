library(poolidentify)
set.seed(8176)
# A single mixed location: exactly one size-2 negative and one size-5 positive.
# Its likelihood at an atom X is g(X) = X^2(1-X^5). Dirichlet second moments
# and independence of atoms yield an analytic posterior mean of prevalence.
h <- matrix(0L, 2, 2)
h[2, 1] <- 1L
fit <- mixed_bayes_prevalence(h, c(2,5), concentration = .7,
                              draws = 1250, burn = 100, chains = 3)
M <- function(k) 2/(2+k)
g <- M(2) - M(7)
same <- (1 + .7 / fit$atoms) / 1.7
truth <- 1 - same * (M(3) - M(8)) / g - (1 - same) * M(1)
stopifnot(abs(fit$mean - truth) < .014,
          all(fit$draws$m2 >= fit$draws$m5 & fit$draws$m5 >= fit$draws$m7),
          all(fit$draws$prevalence >= 1-fit$draws$m2^(1/2) - 1e-12),
          all(fit$draws$prevalence <= 1-fit$draws$m2 + 1e-12))

# One unpooled individual per location. With atom_rate = 2 the prior
# on prevalence is Beta(1,2); 10 positives and 20 negatives update it exactly.
fit <- mixed_bayes_prevalence(c(10L,20L), 1, atom_rate = 2,
                              draws = 900, burn = 80, chains = 2)
stopifnot(fit$atoms == 1L, abs(fit$mean - 11/33) < .01,
          max(abs(fit$interval - qbeta(c(.025,.975), 11,22))) < .017)
