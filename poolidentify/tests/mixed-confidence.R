library(poolidentify)

# Construct two different populations that give exactly the same joint count law.
# Negative-count axes represent 0 and 1 negative pool for each size.
s <- c(2, 5)
p <- c(0.013, 0.141, 0.326, 0.673, 0.971)
kernel <- vapply(p, function(value) {
  negative <- (1 - value)^s
  as.vector(outer(c(1 - negative[1], negative[1]),
                  c(1 - negative[2], negative[2])))
}, numeric(4))
null <- qr.Q(qr(t(kernel)), complete = TRUE)[,5]
pos <- pmax(null, 0) / sum(pmax(null, 0))
neg <- pmax(-null, 0) / sum(pmax(-null, 0))
pmf <- drop(kernel %*% pos)
means <- range(sum(p * pos), sum(p * neg))
stopifnot(max(abs(kernel %*% neg - pmf)) < 1e-12, diff(means) > 1e-3)

# Count simultaneous inclusion of both observationally equivalent prevalences.
level <- 0.80
coverage <- 0
for (i in 0:4) for (j in 0:(4-i)) for (k in 0:(4-i-j)) {
  h <- matrix(c(i, j, k, 4-i-j-k), 2, 2)
  fit <- mixed_confidence(h, s, level, grid_size = 6)
  included <- fit$confidence[1] <= means[1] && fit$confidence[2] >= means[2]
  band_event <- all(fit$lower_probability <= pmf & pmf <= fit$upper_probability)
  stopifnot(!band_event || included,
            abs(sum(fit$fitted_pmf) - 1) < 1e-9)
  coverage <- coverage + dmultinom(as.numeric(h), prob = pmf) * included
}
stopifnot(coverage >= level - 1e-12)
# Axis permutation leaves the inference unchanged on a non-symmetric histogram.
h <- matrix(c(8, 12, 3, 77), 2)
a <- mixed_confidence(h, s, grid_size = 32)
b <- mixed_confidence(t(h), rev(s), grid_size = 32)
stopifnot(max(abs(a$confidence - b$confidence)) < 1e-6)
