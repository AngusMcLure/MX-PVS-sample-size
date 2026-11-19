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
  if (is.null(grid_size)) grid_size <- max(128, ceiling(4 * n^(1/4)))
  if (!is.numeric(grid_size) || length(grid_size) != 1L || !is.finite(grid_size) ||
      grid_size < 2 || grid_size != floor(grid_size)) stop("grid_size must be an integer >= 2.")

  # Include nodes on the individual prevalence and largest-pool positivity scales.
  uniform <- seq(0, 1, length.out = grid_size + 1L)
  grid <- sort(unique(c(uniform, -expm1(log1p(-uniform) / max(pool_sizes)))))
  nodes <- length(grid)
  cells <- nodes - 1L
  outcomes <- length(histogram)
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

  h <- as.numeric(histogram)
  fit <- drop(design %*% grid_mixture_weights(design, h / n))
  # Bonferroni covers the joint multinomial vector without independent cells.
  tail <- (1 - level) / (2 * outcomes)
  lo <- rep(0, outcomes)
  hi <- rep(1, outcomes)
  nonzero <- h > 0
  nonfull <- h < n
  lo[nonzero] <- stats::qbeta(tail, h[nonzero], n - h[nonzero] + 1)
  hi[nonfull] <- stats::qbeta(tail, h[nonfull] + 1, n - h[nonfull], lower.tail = FALSE)
  lo <- pmin(lo, fit)
  hi <- pmax(hi, fit)
  ci <- c(lower = max(0, kernel_envelope_min(grid[-nodes], kernel_lo, kernel_hi, lo, hi)),
          upper = min(1, -kernel_envelope_min(-grid[-1L], kernel_lo, kernel_hi, lo, hi)))
  as_hist <- function(x) array(x, shape, dimnames = dimnames(histogram))
  list(confidence = ci, fitted_pmf = as_hist(fit), lower_probability = as_hist(lo),
       upper_probability = as_hist(hi), distance = max(abs(fit - h / n)),
       grid_mesh = max(diff(grid)), level = level)
}
