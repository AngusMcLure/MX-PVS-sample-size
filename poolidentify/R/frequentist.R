# Empirical CDF projection uses one row per distinct count threshold.
# The last CDF ordinate is always one and adds no constraint.
location_cdf_grid <- function(prevalence, pools, pool_size) {
  negative <- exp(pool_size * log1p(-prevalence))
  matrix(stats::pbinom(rep(pools - seq_len(pools), times = length(prevalence)),
                       size = pools, prob = rep(negative, each = pools),
                       lower.tail = FALSE), nrow = pools)
}

# A dual certificate remains a lower bound after repairing its multipliers.
# kernel_lo and kernel_hi bound each column throughout an entire prevalence cell.
kernel_envelope_min <- function(cost, kernel_lo, kernel_hi, lower, upper) {
  v <- length(lower)
  # For y,z >= 0 and arbitrary t, t + y'lower - z'upper is a
  # lower bound whenever t + y'kernel_hi - z'kernel_lo <= cost.
  constraints <- rbind(kernel_hi, kernel_lo, 1)
  directions <- c(rep(">=", v), rep("<=", v), "=")
  result <- lpSolve::lp("min", cost, constraints, directions, c(lower, upper, 1),
                        compute.sens = TRUE)
  if (result$status != 0L) {
    result <- lpSolve::lp("min", cost, constraints, directions, c(lower, upper, 1),
                          compute.sens = TRUE, scale = 0)
  }
  if (result$status != 0L) stop("Could not solve the prevalence envelope LP.")
  y <- pmax(0, result$duals[seq_len(v)])
  z <- pmax(0, -result$duals[v + seq_len(v)])
  intercept <- min(cost - drop(crossprod(y, kernel_hi)) +
                     drop(crossprod(z, kernel_lo)))
  safety <- 64 * .Machine$double.eps * (v + 1) * (1 + sum(y) + sum(z))
  intercept + sum(y * lower) - sum(z * upper) - safety
}

grid_mixture_weights <- function(design, target) {
  nodes <- ncol(design)
  v <- length(target)
  objective <- c(rep(0, nodes), 1)
  design_lp <- rbind(cbind(design, -1), cbind(-design, -1), c(rep(1, nodes), 0))
  directions <- c(rep("<=", 2L * v), "=")
  rhs <- c(target, -target, 1)
  projection <- lpSolve::lp("min", objective, design_lp, directions, rhs)
  if (projection$status != 0L) {
    projection <- lpSolve::lp("min", objective, design_lp, directions, rhs, scale = 0)
  }
  if (projection$status != 0L) stop("Could not project location counts onto feasible mixtures.")
  weights <- pmax(0, projection$solution[seq_len(nodes)])
  weights / sum(weights)
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
  if (is.null(grid_size)) grid_size <- max(128, ceiling(4 * n^(1/4)))
  if (!is.numeric(grid_size) || length(grid_size) != 1L ||
      !is.finite(grid_size) || grid_size < 2 || grid_size != floor(grid_size)) {
    stop("grid_size must be an integer number of subdivisions of at least two.")
  }
  observed <- histogram / n
  if (v == 1L) {
    # Exact projection of a Clopper-Pearson interval through monotone endpoints.
    inverse <- function(q) -expm1(log1p(-q) / pool_size)
    endpoints <- c(lower = unname(inverse(observed[2])), upper = unname(observed[2]))
    ci <- as.numeric(stats::binom.test(histogram[2], n, conf.level = level)$conf.int)
    return(list(identified = endpoints, midpoint = unname(mean(endpoints)),
                confidence = c(lower = inverse(ci[1]), upper = ci[2]),
                fitted_pmf = observed, distance = 0, radius = NA_real_,
                grid_mesh = 0, identified_method = "sharp", level = level))
  }
  u <- seq(0, 1, length.out = grid_size + 1L)
  grid <- sort(unique(c(u, -expm1(log1p(-u) / pool_size))))
  design <- location_cdf_grid(grid, v, pool_size)
  cdf <- cumsum(observed)[seq_len(v)]
  nodes <- length(grid)
  weights <- grid_mixture_weights(design, cdf)
  fitted_cdf <- drop(design %*% weights)
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

  distance <- max(abs(fitted_cdf - cdf))
  radius <- max(distance, sqrt(log(2 / (1 - level)) / (2 * n)))
  lower <- pmax(0, cdf - radius)
  upper <- pmin(1, cdf + radius)
  kernel_lo <- design[, -1L, drop = FALSE]
  kernel_hi <- design[, -nodes, drop = FALSE]
  left <- grid[-nodes]
  right <- grid[-1L]
  outer_interval <- function(lo, hi) {
    c(lower = max(0, kernel_envelope_min(left, kernel_lo, kernel_hi, lo, hi)),
      upper = min(1, -kernel_envelope_min(-right, kernel_lo, kernel_hi, lo, hi)))
  }
  ci <- outer_interval(lower, upper)
  if (is.null(identified)) {
    # Near numerical saturation, a sharp quadrature can be unresolved. The
    # fitted CDF still defines a valid target for the off-grid envelope.
    rounding <- 64 * .Machine$double.eps * (v + 1)
    identified <- outer_interval(pmax(0, fitted_cdf - rounding),
                                 pmin(1, fitted_cdf + rounding))
    identified[1L] <- max(identified[1L], ci[1L])
    identified[2L] <- min(identified[2L], ci[2L])
  }
  list(identified = identified, midpoint = unname(mean(identified)), confidence = ci,
       fitted_pmf = fitted, distance = distance, radius = radius,
       grid_mesh = max(diff(grid)), identified_method = identified_method, level = level)
}
