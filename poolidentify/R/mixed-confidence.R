mixed_count_grid <- function(shape, pool_sizes, grid) {
  nodes <- length(grid)
  cells <- nodes - 1L
  outcomes <- prod(shape)
  counts <- arrayInd(seq_len(outcomes), shape) - 1L
  design <- matrix(1, outcomes, nodes)
  kernel_lo <- kernel_hi <- matrix(1, outcomes, cells)
  for (a in seq_along(shape)) {
    v <- shape[a] - 1L
    k <- 0:v
    negative <- exp(pool_sizes[a] * log1p(-grid))
    factors <- matrix(stats::dbinom(rep(k, nodes), v, rep(negative, each = v + 1L)),
                      nrow = v + 1L)
    # As a function of the success probability, Bin(v, t)[k] is unimodal
    # with maximum at k/v, including the two boundary modes.
    mode <- -expm1(log(k / v) / pool_sizes[a])
    mode_cell <- pmin(cells, pmax(1L, findInterval(mode, grid)))
    lower <- pmin(factors[, -nodes, drop = FALSE], factors[, -1L, drop = FALSE])
    upper <- pmax(factors[, -nodes, drop = FALSE], factors[, -1L, drop = FALSE])
    upper[cbind(k + 1L, mode_cell)] <- stats::dbinom(k, v, k / v)
    select <- counts[, a] + 1L
    design <- design * factors[select, , drop = FALSE]
    kernel_lo <- kernel_lo * lower[select, , drop = FALSE]
    kernel_hi <- kernel_hi * upper[select, , drop = FALSE]
  }
  list(design = design, kernel_lo = kernel_lo, kernel_hi = kernel_hi)
}

mixed_count_solver <- function(shape, pool_sizes, n, level, grid) {
  kernels <- mixed_count_grid(shape, pool_sizes, grid)
  design <- kernels$design
  project <- make_grid_projection(design)
  outer_interval <- make_outer_interval(grid, kernels$kernel_lo, kernels$kernel_hi)
  tail <- (1 - level) / (2 * prod(shape))
  function(histogram) {
    h <- as.numeric(histogram)
    fit <- drop(design %*% project(h / n))
    # Bonferroni covers the multinomial cell probabilities simultaneously.
    band <- binomial_limits(h, n, tail)
    lo <- pmin(band$lo, fit)
    hi <- pmax(band$hi, fit)
    list(confidence = outer_interval(lo, hi), fitted_pmf = fit,
         lower_probability = lo, upper_probability = hi,
         distance = max(abs(fit - h / n)))
  }
}

mixed_confidence <- function(histogram, pool_sizes, level = 0.95, grid_size = NULL) {
  shape <- dim(histogram)
  if (is.null(shape)) shape <- length(histogram)
  if (!is.numeric(histogram) || any(shape < 2L) || any(!is.finite(histogram)) ||
      any(histogram < 0 | histogram != floor(histogram)) ||
      !is.finite(sum(histogram)) || sum(histogram) < 1) {
    stop("Supply a nonempty location histogram indexed by negative counts on each axis.")
  }
  if (!is.numeric(pool_sizes) || length(shape) != length(pool_sizes) ||
      any(!is.finite(pool_sizes)) || any(pool_sizes < 1 | pool_sizes != floor(pool_sizes))) {
    stop("Supply one positive integer pool size per axis of the histogram.")
  }
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) stop("level must lie strictly between zero and one.")
  n <- sum(histogram)
  grid <- confidence_grid(n, max(pool_sizes), grid_size)
  solve <- mixed_count_solver(shape, pool_sizes, n, level, grid)
  fit <- solve(histogram)
  as_hist <- function(x) array(x, shape, dimnames = dimnames(histogram))
  list(confidence = fit$confidence, fitted_pmf = as_hist(fit$fitted_pmf),
       lower_probability = as_hist(fit$lower_probability),
       upper_probability = as_hist(fit$upper_probability), distance = fit$distance,
       grid_mesh = max(diff(grid)), level = level)
}
