# Empirical CDF projection uses one row per distinct count threshold.
# The last CDF ordinate is always one and adds no constraint.
location_cdf_grid <- function(prevalence, pools, pool_size) {
  negative <- exp(pool_size * log1p(-prevalence))
  matrix(stats::pbinom(rep(pools - seq_len(pools), times = length(prevalence)),
                       size = pools, prob = rep(negative, each = pools),
                       lower.tail = FALSE), nrow = pools)
}

confidence_grid <- function(n, pool_size, grid_size) {
  if (is.null(grid_size)) grid_size <- max(128, ceiling(4 * n^(1/4)))
  if (!is.numeric(grid_size) || length(grid_size) != 1L ||
      !is.finite(grid_size) || grid_size < 2 || grid_size != floor(grid_size)) {
    stop("grid_size must be an integer number of subdivisions of at least two.")
  }
  u <- seq(0, 1, length.out = grid_size + 1L)
  sort(unique(c(u, -expm1(log1p(-u) / pool_size))))
}

# Freeze the design matrices so that many future histograms can share them.
make_grid_projection <- function(design) {
  nodes <- ncol(design)
  v <- nrow(design)
  objective <- c(rep(0, nodes), 1)
  constraints <- t(rbind(cbind(design, -1), cbind(-design, -1), c(rep(1, nodes), 0)))
  directions <- c(rep("<=", 2L * v), "=")
  function(target) {
    rhs <- c(target, -target, 1)
    fit <- lpSolve::lp("min", objective, constraints, directions, rhs, transpose.constraints = FALSE)
    if (fit$status != 0L) {
      fit <- lpSolve::lp("min", objective, constraints, directions, rhs,
                         transpose.constraints = FALSE, scale = 0)
    }
    if (fit$status != 0L) stop("Could not project location counts onto feasible mixtures.")
    weights <- pmax(0, fit$solution[seq_len(nodes)])
    weights / sum(weights)
  }
}

make_outer_interval <- function(grid, kernel_lo, kernel_hi) {
  v <- nrow(kernel_lo)
  constraints <- t(rbind(kernel_hi, kernel_lo, 1))
  directions <- c(rep(">=", v), rep("<=", v), "=")
  left <- grid[-length(grid)]
  right <- grid[-1L]
  bound <- function(cost, lo, hi) {
    rhs <- c(lo, hi, 1)
    fit <- lpSolve::lp("min", cost, constraints, directions, rhs,
                       compute.sens = TRUE, transpose.constraints = FALSE)
    if (fit$status != 0L) {
      fit <- lpSolve::lp("min", cost, constraints, directions, rhs,
                         compute.sens = TRUE, transpose.constraints = FALSE, scale = 0)
    }
    if (fit$status != 0L) stop("Could not solve the prevalence envelope LP.")
    # Repair dual multipliers to certify the bound even for off-grid atoms.
    y <- pmax(0, fit$duals[seq_len(v)])
    z <- pmax(0, -fit$duals[v + seq_len(v)])
    intercept <- min(cost - drop(crossprod(y, kernel_hi)) + drop(crossprod(z, kernel_lo)))
    safety <- 64 * .Machine$double.eps * (v + 1) * (1 + sum(y) + sum(z))
    intercept + sum(y * lo) - sum(z * hi) - safety
  }
  function(lo, hi) {
    c(lower = max(0, bound(left, lo, hi)), upper = min(1, -bound(-right, lo, hi)))
  }
}

binomial_limits <- function(h, n, tail) {
  # Beta quantiles extend continuously to fractional planning histograms.
  lo <- numeric(length(h))
  hi <- rep(1, length(h))
  nonzero <- h > 0
  nonfull <- h < n
  lo[nonzero] <- stats::qbeta(tail, h[nonzero], n - h[nonzero] + 1)
  hi[nonfull] <- stats::qbeta(tail, h[nonfull] + 1, n - h[nonfull], lower.tail = FALSE)
  list(lo = lo, hi = hi)
}

equal_count_solver <- function(v, s, n, level, grid) {
  design <- location_cdf_grid(grid, v, s)
  project <- make_grid_projection(design)
  outer_interval <- make_outer_interval(grid, design[, -1L, drop = FALSE],
                                        design[, -ncol(design), drop = FALSE])
  dkw <- sqrt(log(2 / (1 - level)) / (2 * n))
  fit <- function(histogram) {
    cdf <- cumsum(histogram / n)[seq_len(v)]
    weights <- project(cdf)
    fitted_cdf <- drop(design %*% weights)
    distance <- max(abs(fitted_cdf - cdf))
    radius <- max(distance, dkw)
    confidence <- outer_interval(pmax(0, cdf - radius), pmin(1, cdf + radius))
    list(confidence = confidence, weights = weights, fitted_cdf = fitted_cdf,
         distance = distance, radius = radius)
  }
  list(fit = fit, outer_interval = outer_interval)
}

#' See man/fit_prevalence.Rd for definitions of the three interval summaries.
fit_prevalence <- function(histogram, pool_size, level = 0.95, grid_size = NULL) {
  if (!is.numeric(histogram) || !is.null(dim(histogram)) || length(histogram) < 2L ||
      any(!is.finite(histogram)) || any(histogram < 0) ||
      any(histogram != floor(histogram)) || sum(histogram) < 1 ||
      !is.finite(sum(histogram))) {
    stop("histogram must give nonnegative integer location counts for K=0,...,v.")
  }
  if (!is.numeric(pool_size) || length(pool_size) != 1L ||
      !is.finite(pool_size) || pool_size < 1 || pool_size != floor(pool_size)) {
    stop("pool_size must be a positive integer.")
  }
  if (!is.numeric(level) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1) {
    stop("level must lie strictly between zero and one.")
  }
  n <- sum(histogram)
  v <- length(histogram) - 1L
  grid <- confidence_grid(n, pool_size, grid_size)
  observed <- histogram / n
  if (v == 1L) {
    # Exact projection of a Clopper-Pearson interval through monotone endpoints.
    inverse <- function(q) -expm1(log1p(-q) / pool_size)
    endpoints <- c(lower = unname(inverse(observed[2])), upper = unname(observed[2]))
    ci <- binomial_limits(histogram[2], n, (1 - level) / 2)
    return(list(identified = endpoints, midpoint = unname(mean(endpoints)),
                confidence = c(lower = inverse(ci$lo), upper = ci$hi),
                fitted_pmf = observed, distance = 0, radius = NA_real_,
                grid_mesh = 0, identified_method = "sharp", level = level))
  }
  solver <- equal_count_solver(v, pool_size, n, level, grid)
  fit <- solver$fit(histogram)
  weights <- fit$weights
  active <- weights > 0
  negative <- exp(pool_size * log1p(-grid[active]))
  k <- 0:v
  pmfs <- matrix(stats::dbinom(rep(v - k, times = sum(active)), v,
                               rep(negative, each = v + 1L)), nrow = v + 1L)
  fitted <- drop(pmfs %*% weights[active])
  # A nonnegative polynomial of degree at most v, with simple roots at
  # support endpoints and double roots at interior atoms, exposes this face.
  support <- grid[active]
  index_twice <- sum(ifelse(support == 0 | support == 1, 1L, 2L))
  identified_method <- "sharp"
  if (index_twice <= v) {
    p <- sum(support * weights[active])
    identified <- c(lower = p, upper = p)
  } else {
    identified <- tryCatch(identification_interval(pool_size = pool_size, probabilities = fitted),
                           error = function(e) NULL)
    if (is.null(identified)) identified_method <- "grid outer bounds"

  }

  ci <- fit$confidence
  if (is.null(identified)) {
    # Near numerical saturation, a sharp quadrature can be unresolved. The
    # fitted CDF still defines a valid target for the off-grid envelope.
    rounding <- 64 * .Machine$double.eps * (v + 1)
    identified <- solver$outer_interval(pmax(0, fit$fitted_cdf - rounding),
                                        pmin(1, fit$fitted_cdf + rounding))
    identified[1L] <- max(identified[1L], ci[1L])
    identified[2L] <- min(identified[2L], ci[2L])
  }
  list(identified = identified, midpoint = unname(mean(identified)), confidence = ci,
       fitted_pmf = fitted, distance = fit$distance, radius = fit$radius,
       grid_mesh = max(diff(grid)), identified_method = identified_method, level = level)
}
