# Conditional on the total number of infected individuals, their positions
# form a simple random sample. Exposing positions one at a time avoids the
# alternating sums of large beta moments in the usual inclusion-exclusion PMF.
.pi_occupancy <- function(s, v) {
  total <- s * v
  k <- 0:v
  a <- matrix(0, total + 1L, v + 1L)
  a[1L, 1L] <- 1
  for (h in seq_len(total)) {
    previous <- a[h, ]
    enter <- pmin(1, s * (v - k) / (total - h + 1))
    new <- previous * enter
    a[h + 1L, ] <- previous - new + c(0, new[seq_len(v)])
  }
  a
}

.pi_positive_integer <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0 && x == floor(x)
}

.pi_check_prob <- function(x) {
  if (!is.numeric(x) || length(x) < 2L || any(!is.finite(x)) ||
      any(x < 0) || abs(sum(x) - 1) > 1e-9)
    stop("Supply probabilities for counts 0 through v, summing to one.")
  unname(x / sum(x))
}

.pi_legendre <- function(n) {
  i <- seq_len(n - 1L)
  a <- matrix(0, n, n)
  a[cbind(i, i + 1L)] <- i / sqrt(4 * i^2 - 1)
  e <- eigen(a + t(a), symmetric = TRUE)
  list(x = e$values, w = 2 * e$vectors[1, ]^2)
}

# Make the design-dependent part once per fit or simulation.
.pi_count_engine <- function(s, v, family, quadrature) {
  if (!.pi_positive_integer(s) || !.pi_positive_integer(v))
    stop("pool_size and pools must be positive integers.")
  if (family == "beta") {
    occupancy <- .pi_occupancy(s, v)
    h <- 0:(s * v)
    logchoose <- lchoose(s * v, h)
    return(function(par) {
      p <- par[1]
      rho <- par[2]
      if (p == 0 || p == 1 || rho == 0) {
        prob <- stats::dbinom(0:v, v, -expm1(s * log1p(-p)))
      } else if (rho == 1) {
        prob <- c(1 - p, rep(0, v - 1L), p)
      } else {
        # Adjacent beta-binomial masses, normalized in log space.
        i <- seq_len(s * v)
        logmass <- logchoose + c(0, cumsum(log((i - 1) * rho + p * (1 - rho)) -
                             log((s * v - i) * rho + (1 - p) * (1 - rho))))
        mass <- exp(logmass - max(logmass))
        prob <- drop(crossprod(mass / sum(mass), occupancy))
      }
      list(probabilities = prob, prevalence = p)
    })
  }
  if (!.pi_positive_integer(quadrature) || quadrature < 16)
    stop("quadrature must be an integer of at least 16.")
  rule <- .pi_legendre(max(4, ceiling(quadrature / 8)))
  function(par) {
    mu <- par[1]
    sd <- par[2]
    if (sd == 0) {
      z <- 0
      w <- 1
    } else {
      # Resolve both the normal density and the transition of the binomial
      # kernel even at large logit standard deviations.
      cuts <- (seq(-40 - log(s * v), 40, by = 2) - mu) / sd
      cuts <- sort(unique(c(seq(-10, 10, by = 2), cuts[cuts > -10 & cuts < 10])))
      mid <- (cuts[-1L] + cuts[-length(cuts)]) / 2
      half <- diff(cuts) / 2
      z <- as.vector(outer(rule$x, half) + rep(mid, each = length(rule$x)))
      w <- as.vector(outer(rule$w, half)) * stats::dnorm(z)
      w <- w / sum(w)
    }
    eta <- mu + sd * z
    log_negative <- s * stats::plogis(-eta, log.p = TRUE)
    q <- -expm1(log_negative)
    negative <- exp(log_negative)
    prob <- vapply(0:v, function(k) {
      pmf <- ifelse(log_negative > -log(2), stats::dbinom(k, v, q),
                    stats::dbinom(v - k, v, negative))
      sum(w * pmf)
    }, numeric(1))
    list(probabilities = prob / sum(prob), prevalence = sum(w * stats::plogis(eta)))
  }
}

mixing_counts <- function(pool_size, pools, family = c("beta", "logitnormal"),
                          par, quadrature = 80L) {
  family <- match.arg(family)
  if (!is.numeric(par) || length(par) != 2L || any(!is.finite(par)) ||
      (family == "beta" && any(par < 0 | par > 1)) ||
      (family == "logitnormal" && par[2] < 0))
    stop("Invalid mixing-family parameters.")
  .pi_count_engine(pool_size, pools, family, quadrature)(unname(par))
}

.pi_working_fit <- function(prob, s, family, engine, start = NULL) {
  v <- length(prob) - 1L
  q <- sum((0:v) * prob) / v
  p0 <- -expm1(log1p(-q) / s)
  homogeneous <- stats::dbinom(0:v, v, q)
  if (family == "naive")
    return(list(prevalence = p0, probabilities = homogeneous, par = p0,
                status = "independence"))
  if (v < 2) stop("Estimating a two-parameter family requires at least two pools per location.")
  if (q == 0 || q == 1 || sum(prob[-c(1L, v + 1L)]) == 0)
    return(list(prevalence = q, probabilities = prob, par = c(NA_real_, NA_real_),
                status = "endpoint mixture"))
  positive <- prob > 0
  loss <- function(pmf) {
    if (any(pmf[positive] <= 0)) return(1e12)
    -sum(prob[positive] * log(pmf[positive]))
  }
  natural <- if (family == "beta") stats::plogis else function(x) c(x[1], exp(x[2]))
  objective <- function(x) loss(engine(natural(x))$probabilities)
  if (family == "beta") {
    lo <- c(-26, -18)
    hi <- c(26, 16)
    starts <- list(c(stats::qlogis(p0), -3), c(stats::qlogis(q), 0))
  } else {
    lo <- c(-120, -8)
    hi <- c(120, log(60))
    starts <- list(c(stats::qlogis(p0), -0.4), c(stats::qlogis(p0) - 2, 1))
  }
  starts <- c(if (!is.null(start)) list(start), starts)
  best <- list(prevalence = p0, probabilities = homogeneous,
               par = if (family == "beta") c(p0, 0) else c(stats::qlogis(p0), 0),
               status = "homogeneous")
  best_loss <- loss(homogeneous)
  success <- FALSE
  entropy <- -sum(prob[positive] * log(prob[positive]))
  for (x in starts) {
    fit <- tryCatch(stats::optim(pmin(hi, pmax(lo, x)), objective, method = "L-BFGS-B",
                                lower = lo, upper = hi,
                                control = list(factr = 1e4, ndeps = rep(2e-5, 2))),
                    error = function(e) NULL)
    if (is.null(fit) || fit$convergence != 0L) next
    success <- TRUE
    if (fit$value >= best_loss - 2e-12) next
    best_loss <- fit$value
    par <- natural(fit$par)
    best <- engine(par)
    best$par <- par
    best$internal <- fit$par
    best$status <- if (any(pmin(fit$par - lo, hi - fit$par) < 1e-4)) "search limit" else "interior"
    if (best_loss - entropy < 1e-13) break
  }
  if (!success) best$status <- "optimization failure"
  best
}

working_projection <- function(probabilities, pool_size,
                               family = c("naive", "beta", "logitnormal"), quadrature = 80L) {
  family <- match.arg(family)
  prob <- .pi_check_prob(probabilities)
  if (!.pi_positive_integer(pool_size)) stop("pool_size must be a positive integer.")
  engine <- if (family != "naive") .pi_count_engine(pool_size, length(prob) - 1L, family, quadrature)
  fit <- .pi_working_fit(prob, pool_size, family, engine)
  keep <- prob > 0
  fit$kl <- max(0, sum(prob[keep] * log(prob[keep] / fit$probabilities[keep])))
  fit
}

.pi_total_count <- function(prob, n) {
  size <- n * (length(prob) - 1L) + 1
  padded <- 2^ceiling(log2(size))
  spectrum <- stats::fft(c(prob, numeric(padded - length(prob))))
  mass <- Re(stats::fft(spectrum^n, inverse = TRUE))[seq_len(size)] / padded
  if (min(mass) < -1e-9) stop("Count convolution lost numerical precision.")
  mass <- pmax(0, mass)
  mass / sum(mass)
}

# Delta-method interval with model-based expected information per location.
.pi_model_interval <- function(fit, engine, family, n, z, v) {
  if (fit$status != "interior") return(c(NA_real_, NA_real_))
  natural <- if (family == "beta") stats::plogis else function(x) c(x[1], exp(x[2]))
  value <- function(x) {
    model <- engine(natural(x))
    c(model$probabilities, model$prevalence)
  }
  step <- 1e-4
  grad <- vapply(seq_len(2L), function(i) {
    d <- diag(2)[, i] * step
    (value(fit$internal + d) - value(fit$internal - d)) / (2 * step)
  }, numeric(v + 2L))
  pi <- fit$probabilities
  deriv <- grad[seq_along(pi), , drop = FALSE]
  positive <- pi > 0
  info <- crossprod(deriv[positive, , drop = FALSE] / sqrt(pi[positive]))
  if (rcond(info) < 1e-11) return(c(NA_real_, NA_real_))
  g <- grad[v + 2L, ]
  se <- sqrt(drop(crossprod(g, solve(info, g))) / n)
  pmin(1, pmax(0, fit$prevalence + c(-1, 1) * z * se))
}

working_performance <- function(probabilities, prevalence, pool_size, family = c("naive", "beta", "logitnormal"),
                                locations, level = 0.95, replicates = 200L, quadrature = 80L) {
  family <- match.arg(family)
  prob <- .pi_check_prob(probabilities)
  if (!is.numeric(prevalence) || length(prevalence) != 1L || !is.finite(prevalence) ||
      prevalence < 0 || prevalence > 1) stop("prevalence must be in [0, 1].")
  if (!.pi_positive_integer(locations) || !.pi_positive_integer(replicates) ||
      !.pi_positive_integer(pool_size)) stop("Invalid design or replication count.")
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("level must lie strictly between zero and one.")
  n <- locations
  s <- pool_size
  v <- length(prob) - 1L
  target <- working_projection(prob, s, family, quadrature)
  z <- stats::qnorm(1 - (1 - level)/2)
  if (family == "naive") {
    mass <- .pi_total_count(prob, n)
    total <- n * v
    qhat <- (0:total) / total
    point <- -expm1(log1p(-qhat)/s)
    # Invert the independent-pool score test to evaluate Wilson coverage.
    qnull <- -expm1(s * log1p(-prevalence))
    covered <- abs(qhat - qnull) <= z * sqrt(qnull * (1 - qnull) / total)
    coverage <- sum(mass[covered])
    return(list(target = target$prevalence, asymptotic_bias = target$prevalence - prevalence,
                bias = sum(mass * point) - prevalence, coverage = coverage,
                coverage_mcse = 0, failure_fraction = 0, method = "exact independent-pool Wilson"))
  }
  engine <- .pi_count_engine(s, v, family, quadrature)
  samples <- stats::rmultinom(replicates, n, prob)
  fitted <- rep(NA_real_, replicates)
  hit <- logical(replicates)
  available <- logical(replicates)
  for (b in seq_len(replicates)) {
    fit <- .pi_working_fit(samples[, b]/n, s, family, engine, target$internal)
    if (!fit$status %in% c("search limit", "optimization failure")) fitted[b] <- fit$prevalence
    ci <- .pi_model_interval(fit, engine, family, n, z, v)
    available[b] <- all(is.finite(ci))
    hit[b] <- available[b] && ci[1] <= prevalence && prevalence <= ci[2]
  }
  coverage <- mean(hit)
  list(target = target$prevalence, asymptotic_bias = target$prevalence - prevalence,
       projection_status = target$status, bias = if (all(is.finite(fitted))) mean(fitted) - prevalence else NA_real_,
       point_failure_fraction = mean(!is.finite(fitted)), coverage = coverage,
       coverage_mcse = sqrt(coverage * (1 - coverage) / replicates),
       failure_fraction = mean(!available), method = "Monte Carlo model-information Wald")
}
