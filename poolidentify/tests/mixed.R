library(poolidentify)

# Test all negative-count cells and repeated exponents for 2+2+4.
x <- c(0.12, 0.6, 1)
w <- c(0.1, 0.6, 0.3)
counts <- matrix(0, 3, 2)
for (i in seq_along(x))
  counts <- counts + w[i] * outer(dbinom(0:2, 2, x[i]^2), dbinom(0:1, 1, x[i]^4))
m <- joint_negative_moments(counts, c(2, 4))
stopifnot(identical(names(m), c("2", "4", "6", "8")),
          max(abs(m - drop(w %*% outer(x, c(2, 4, 6, 8), "^")))) < 1e-13)
est <- mixed_interval(m)
p <- 1 - sum(x*w)
stopifnot(est[1] <= p, est[2] >= p)
bad <- counts
bad[1, 1] <- bad[1, 1] + 0.01
bad[1, 2] <- bad[1, 2] - 0.01
stopifnot(inherits(try(joint_negative_moments(bad, c(2, 4)), silent = TRUE), "try-error"))

# From moments of X^2, both explicit Radau distributions attain sharp endpoints.
sharp <- c(1 - (.98 * sqrt(3/28) + .02), .75)
est <- mixed_interval(c(1/32, 1/8), c(4, 2))
stopifnot(est[1] <= sharp[1], est[2] >= sharp[2], max(abs(est - sharp)) < 3e-5)
stopifnot(identical(unname(mixed_interval(.25, 2)), c(.5, .75)),
          identical(unname(mixed_interval(c(.4, .22), c(1, 2))), c(.6, .6)))

# The 10+40 illustration uses two separately computed principal representations.
k <- c(10, 40, 50)
m <- exp(lbeta(.09, 8.91 + k) - lbeta(.09, 8.91))
coarse <- mixed_interval(m, k, cells = 64)
fine <- mixed_interval(m, k, cells = 1024)
stopifnot(fine[1] <= .0090021, fine[2] >= .0285132,
          max(abs(fine - c(.0090021, .0285132))) < 1.5e-5, diff(fine) < diff(coarse))
stopifnot(inherits(try(mixed_interval(c(.2, .4), c(2, 3)), silent = TRUE), "try-error"),
          inherits(try(mixed_interval(c(.75, .6, .59), 1:3), silent = TRUE), "try-error"))
