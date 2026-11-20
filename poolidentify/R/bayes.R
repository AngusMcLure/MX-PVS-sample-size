# Principal quadrature from canonical moments on [0,1]. Appending 0 or 1
# selects the lower or upper representation, respectively.
.pi_canonical_rule <- function(z, upper = FALSE) {
  v <- length(z)
  z <- c(z, as.integer(upper))
  zeta <- c(1, 1 - utils::head(z, -1L)) * z
  n <- if (upper) (v + 3L) %/% 2L else (v + 2L) %/% 2L
  # For even v both Radau representations use v/2+1 nodes.
  # The Jacobi matrix is B' B for upper bidiagonal B. Its singular values
  # resolve small positive eigenvalues without the cancellation of a
  # symmetric eigendecomposition near saturation.
  zz <- c(zeta, 0)
  odd <- 2L * seq_len(n) - 1L
  B <- diag(sqrt(zz[odd]), nrow = n)
  if (n > 1L) {
    i <- seq_len(n - 1L)
    B[cbind(i, i + 1L)] <- sqrt(zz[2L * i])
  }
  decomp <- svd(B, nu = 0L, nv = n)
  x <- pmin(1, decomp$d^2)
  if (upper) x[1L] <- 1
  if ((v %% 2L == 0L && !upper) || (v %% 2L == 1L && upper)) x[n] <- 0
  list(x = x, w = decomp$v[1L, ]^2)
}

.pi_canonical_counts <- function(z) {
  v <- length(z)
  u <- z[1L]
  if (v == 1L) return(c(1 - u, u))
  if (v == 2L) {
    cross <- u * (1 - u) * (1 - z[2L])
    return(c((1 - u) * (1 - u + u * z[2L]), 2 * cross,
             u * (u + (1 - u) * z[2L])))
  }
  if (u > 0.5) {
    odd <- seq.int(1L, v, 2L)
    z[odd] <- 1 - z[odd]
  }
  rule <- .pi_canonical_rule(z)
  counts <- if (u > 0.5) v:0 else 0:v
  drop(vapply(rule$x, function(x) stats::dbinom(counts, v, x),
              numeric(v + 1L)) %*% rule$w)
}

.pi_canonical_ends <- function(z, s) {
  u <- z[1L]
  fs <- function(q) -expm1(log1p(-q) / s)
  if (s == 1 || u == 0 || u == 1) return(c(u, u))
  v <- length(z)
  if (v == 1L) return(c(fs(u), u))
  if (v == 2L) {
    r <- z[2L]
    w1 <- u * r / (1 - u + u * r)
    return(c(u / (u + (1 - u) * r) *
               -expm1((log1p(-u) + log1p(-r)) / s),
             w1 + (1 - w1) * fs(u * (1 - r))))
  }
  # Work with 1-Q so a node very close to Q=1 need not be subtracted
  # from one before taking a fractional power. Reflection swaps the
  # Radau endpoints for even v and preserves the odd-v patterns.
  odd <- seq.int(1L, v, 2L)
  z[odd] <- 1 - z[odd]
  by_rule <- function(upper) {
    rule <- .pi_canonical_rule(z, upper)
    sum(rule$w * -expm1(log(rule$x) / s))
  }
  c(by_rule(v %% 2L == 0L), by_rule(v %% 2L == 1L))
}

.pi_bayes_diag <- function(x, chain) {
  series <- split(x, chain)
  half <- length(series[[1L]]) %/% 2L
  halves <- do.call(cbind, lapply(series, function(y)
    cbind(utils::head(y, half), utils::tail(y, half))))
  within <- mean(apply(halves, 2L, stats::var))
  between <- stats::var(colMeans(halves))
  rhat <- if (within == 0) {
    if (between > 0) Inf else NA_real_
  } else sqrt((half - 1) / half + between / within)
  block <- floor(sqrt(length(series[[1L]])))
  means <- unlist(lapply(series, function(y) {
    n <- length(y) %/% block
    colMeans(matrix(utils::head(y, n * block), nrow = block))
  }))
  c(split_Rhat = rhat, MCSE_mean = stats::sd(means) / sqrt(length(means)))
}

bayes_prevalence <- function(counts, pool_size, draws = 1000L, burn = 500L,
                             chains = 4L, prior = c("uniform", "shrink", "volume"),
                             shape = NULL, within = NULL, level = 0.95) {
  integer_scalar <- function(x, minimum)
    is.numeric(x) && length(x) == 1L && is.finite(x) && x >= minimum &&
      x <= .Machine$integer.max && x == floor(x)
  if (!is.numeric(counts) || !is.null(dim(counts)) || length(counts) < 2L ||
      any(!is.finite(counts)) || any(counts < 0 | counts != floor(counts)) ||
      sum(counts) == 0) stop("counts must be a nonempty location histogram for K = 0:v")
  if (!integer_scalar(pool_size, 1L) || !integer_scalar(draws, 8L) ||
      !integer_scalar(burn, 0L) || !integer_scalar(chains, 1L))
    stop("use positive integer pool_size and chains, draws >= 8, and burn >= 0")
  if (length(level) != 1L || !is.numeric(level) || !is.finite(level) ||
      level <= 0 || level >= 1) stop("level must be strictly between zero and one")
  v <- length(counts) - 1L
  prior <- match.arg(prior)
  if (is.null(shape)) {
    shape <- matrix(1, nrow = v, ncol = 2)
    if (prior == "shrink") {
      even <- which(seq_len(v) %% 2L == 0L)
      shape[even, 2L] <- 1 + even
    }
    if (prior == "volume") shape[] <- rep(v:1, 2L)
  } else {
    prior <- "custom"
  }
  if (!is.numeric(shape) || !identical(dim(shape), c(v, 2L)) ||
      any(!is.finite(shape)) || any(shape <= 0))
    stop("shape must be a positive finite v by 2 matrix of beta shapes")
  if (!is.null(within) && (!is.numeric(within) || length(within) != 2L ||
                            any(!is.finite(within)) || any(within <= 0)))
    stop("within must be NULL or two positive finite beta shapes")
  positive <- counts > 0
  ll <- function(z) {
    prob <- .pi_canonical_counts(z)[positive]
    sum(counts[positive] * log(prob))
  }
  chain_draws <- function(chain) {
    samples <- matrix(NA_real_, nrow = draws, ncol = v)
    if (v == 1L) {
      samples[, 1L] <- stats::rbeta(draws, shape[1L, 1L] + counts[2L],
                                     shape[1L, 2L] + counts[1L])
    } else {
      y <- stats::runif(v)
      z <- stats::qbeta(y, shape[, 1L], shape[, 2L])
      loglik <- ll(z)
      if (!is.finite(loglik)) stop("initial count likelihood underflow; check prior shapes and histogram")
      for (t in seq_len(burn + draws)) {
        for (k in seq_len(v)) {
          height <- loglik - stats::rexp(1L)
          left <- 0
          right <- 1
          repeat {
            proposal <- stats::runif(1L, left, right)
            z[k] <- stats::qbeta(proposal, shape[k, 1L], shape[k, 2L])
            proposed_ll <- ll(z)
            if (is.finite(proposed_ll) && proposed_ll >= height) break
            if (proposal < y[k]) left <- proposal else right <- proposal
            if (right <= left) stop("slice bracket exhausted at machine precision")
          }
          y[k] <- proposal
          loglik <- proposed_ll
        }
        if (t > burn) samples[t - burn, ] <- z
      }
    }
    endpoints <- t(vapply(seq_len(draws), function(i)
      .pi_canonical_ends(samples[i, ], pool_size), numeric(2)))
    colnames(samples) <- paste0("c", seq_len(v))
    result <- data.frame(chain = chain, samples, lower = endpoints[, 1L],
                         upper = endpoints[, 2L])
    if (!is.null(within)) {
      weight <- stats::rbeta(draws, within[1L], within[2L])
      result$prevalence <- result$lower + weight * (result$upper - result$lower)
    }
    result
  }
  posterior <- do.call(rbind, lapply(seq_len(chains), chain_draws))
  tail_prob <- (1 - level) / 2
  q <- function(x, p) unname(stats::quantile(x, p))
  env <- c(lower = q(posterior$lower, tail_prob),
           upper = q(posterior$upper, 1 - tail_prob))
  interval <- if (!is.null(within)) c(lower = q(posterior$prevalence, tail_prob),
                                     upper = q(posterior$prevalence, 1 - tail_prob))
  diagnostics <- t(vapply(posterior[-1L], .pi_bayes_diag, numeric(2L),
                          chain = posterior$chain))
  list(draws = posterior, envelope = env, interval = interval, level = level,
       diagnostics = diagnostics, prior = prior, shape = shape, within = within)
}
