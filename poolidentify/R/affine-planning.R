# Fit one score per joint count cell. On the planning grid require
# E_u(z)-offset <= u <= E_u(z)-offset+gap, where 0 <= z <= span.
# Translation of all scores is absorbed by offset.
affine_grid_score <- function(kernel, grid, penalty) {
  cells <- nrow(kernel)
  nodes <- length(grid)
  index <- which(kernel != 0, arr.ind = TRUE)
  kernel_rows <- cbind(index[, 2L], index[, 1L], kernel[index])
  upper_rows <- kernel_rows
  upper_rows[, 1L] <- upper_rows[, 1L] + nodes
  i <- seq_len(cells)
  j <- seq_len(nodes)
  # Sparse triplets avoid building a cells-by-cells identity matrix.
  coefficients <- rbind(kernel_rows, upper_rows,
                        cbind(c(j, nodes + j), cells + 1L, -1),
                        cbind(nodes + j, cells + 2L, 1),
                        cbind(2L * nodes + i, i, 1),
                        cbind(2L * nodes + i, cells + 3L, -1))
  result <- lpSolve::lp("min", c(numeric(cells + 1L), 1, penalty),
                        dense.const = coefficients,
                        const.dir = c(rep("<=", nodes), rep(">=", nodes), rep("<=", cells)),
                        const.rhs = c(grid, grid, numeric(cells)), scale = 0)
  if (result$status != 0L || any(!is.finite(result$solution)))
    stop("Could not fit the design score on the prevalence grid.")
  z <- result$solution[i]
  z - mean(range(z))
}

affine_pool_interval <- function(design, locations, level = 0.95, grid_size = 64L,
                                 histogram = NULL) {
  design <- .planning_sizes(design)
  if (!.pi_positive_integer(locations)) stop("locations must be a positive integer.")
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) stop("level must lie strictly between zero and one.")
  if (!.pi_positive_integer(grid_size) || grid_size < 2)
    stop("grid_size must be an integer of at least two.")
  sizes <- sort(unique(design))
  pools <- tabulate(match(design, sizes))
  shape <- pools + 1L
  cells <- prod(shape)
  if (cells > .Machine$integer.max)
    stop("Joint count array is too large; use fewer distinct sizes or pools.")
  if (!is.null(histogram)) {
    dimensions <- dim(histogram)
    if (is.null(dimensions)) dimensions <- length(histogram)
    if (!is.numeric(histogram) || !identical(dimensions, shape) ||
        any(!is.finite(histogram)) || any(histogram < 0 | histogram != floor(histogram)) ||
        sum(histogram) != locations)
      stop("Supply a negative-count histogram with the design dimensions and total locations.")
  }
  penalty <- sqrt(2 * log(2 / (1 - level)) / locations)
  if (penalty >= 1) {
    score <- numeric(cells)
    remainder <- c(0, 1)
  } else if (sizes[1L] == 1L) {
    # The singleton class observes prevalence directly and achieves the
    # design-independent lower bound on this uniform Hoeffding width.
    score <- rep(1 - (0:pools[1L]) / pools[1L], length.out = cells)
    remainder <- c(0, 0)
  } else if (length(design) == 1L) {
    s <- as.double(design)
    slope <- min(1, penalty^(1/s - 1) / s)
    gap <- (1 - 1/s) * exp(-log(s * slope) / (s - 1))
    score <- c(slope, 0)
    remainder <- c(1 - slope - gap, 1 - slope)
  } else {
    grid <- confidence_grid(locations, max(sizes), grid_size)
    kernels <- mixed_count_grid(shape, sizes, grid)
    score <- affine_grid_score(kernels$design, grid, penalty)
    # Certify u - E_u(score) between nodes by two independent inequalities.
    # The conditional expectation is Bernstein of degree M = sum(design)
    # with coefficients in the score range (condition on the infection total).
    units <- sum(as.double(design))
    residual <- grid - drop(crossprod(score, kernels$design))
    error <- (diff(grid) * units) * (diff(grid) * (units - 1)) * diff(range(score)) / 4
    left <- seq_len(length(grid) - 1L)
    right <- left + 1L
    chord_lo <- pmin(residual[left], residual[right]) - error
    chord_hi <- pmax(residual[left], residual[right]) + error
    positive <- pmax(score, 0)
    negative <- pmin(score, 0)
    score_lo <- drop(crossprod(positive, kernels$kernel_lo) +
                     crossprod(negative, kernels$kernel_hi))
    score_hi <- drop(crossprod(positive, kernels$kernel_hi) +
                     crossprod(negative, kernels$kernel_lo))
    margin <- 64 * .Machine$double.eps * (cells + 1) * (1 + max(abs(score)))
    remainder <- c(min(pmax(chord_lo, grid[left] - score_hi)) - margin,
                   max(pmin(chord_hi, grid[right] - score_lo)) + margin)
    if (diff(remainder) + penalty * diff(range(score)) >= 1) {
      score <- numeric(cells)
      remainder <- c(0, 1)
    }
  }
  sampling <- penalty * diff(range(score))
  interval <- NULL
  if (!is.null(histogram)) {
    center <- sum((as.numeric(histogram) / locations) * score)
    interval <- stats::setNames(pmax(0, pmin(1, center + remainder + c(-.5, .5) * sampling)),
                                c("lower", "upper"))
  }
  list(width_bound = min(1, diff(remainder) + sampling), structural_gap = diff(remainder),
       sampling_allowance = sampling, score = array(score, dim = shape),
       remainder = stats::setNames(remainder, c("lower", "upper")), interval = interval)
}
