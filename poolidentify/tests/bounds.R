library(poolidentify)

near <- function(x, y, tol = 1e-8)
  isTRUE(all.equal(unname(x), unname(y), tolerance = tol, check.attributes = FALSE))
law <- function(v, q, w)
  drop(vapply(q, function(z) dbinom(0:v, v, z), numeric(v + 1L)) %*% w)

# Three simple endpoints: unpooled, all-negative and perfect within-site dependence.
stopifnot(near(identification_interval(0.4, 1), c(.4, .4)),
          near(identification_interval(0.36, 2), c(.2, .36)),
          near(identification_interval(pool_size = 12, probabilities = c(1,0,0)), c(0,0)),
          near(identification_interval(pool_size = 5, probabilities = c(.7,0,0,.3)), c(.3,.3)))

# Independent oracle: prescribe the principal representations themselves.
# Binomial laws then give known exact values for alternating Gauss/Radau/Lobatto rules.
for (v in 2:8) {
  odd <- v %% 2 == 1
  for (lower in c(TRUE, FALSE)) {
    endpoints <- if (lower && !odd) 0 else if (!lower && odd) c(0, 1) else if (!lower) 1
    n <- (v + 1 - length(endpoints)) %/% 2
    q <- c(endpoints, seq_len(n) / (n + 1))
    w <- seq_along(q) / sum(seq_along(q))
    prob <- law(v, q, w)
    truth <- sum(w * -expm1(log1p(-q) / 7))
    bounds <- identification_interval(probabilities = prob, pool_size = 7)
    stopifnot(near(bounds[if (lower) 1 else 2], truth),
              near(count_moments(prob), vapply(1:v, function(k) sum(w * q^k), 0.0)),
              near(identification_interval(count_moments(prob), 7), bounds, 1e-6))
  }
}

# Nested intervals for one fixed population and two design changes.
p <- c(.005, .04, .15, .35, .55, .82)
w <- c(.2,.3,.12,.22,.11,.05)
prev <- sum(w * p)
old <- c(0,1)
for (v in 1:8) {
  new <- identification_interval(pool_size=4, probabilities=law(v, 1-(1-p)^4, w))
  small <- identification_interval(pool_size=2, probabilities=law(v, 1-(1-p)^2, w))
  stopifnot(old[1] <= new[1]+1e-8, old[2] >= new[2]-1e-8,
            new[1] <= small[1]+1e-8, new[2] >= small[2]-1e-8,
            small[1] <= prev+1e-8, small[2] >= prev-1e-8)
  old <- new
}

# Nonnegative count probabilities alone do not establish feasibility.
for (args in list(list(moments=c(.75,.1), pool_size=2),
                  list(probabilities=c(0,.75,.25), pool_size=1),
                  list(probabilities=c(0,.5,0,.5,0), pool_size=3),
                  list(moments=.5, probabilities=c(.5,.5), pool_size=3))) {
  stopifnot(inherits(try(do.call(identification_interval, args), silent=TRUE), "try-error"))
}
