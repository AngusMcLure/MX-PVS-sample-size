# One tensor-product binomial factorial-moment transform per count dimension.
joint_negative_moments <- function(probabilities, pool_sizes, tolerance = 1e-9) {
  shape <- dim(probabilities)
  if (is.null(shape)) shape <- length(probabilities)
  if (!is.numeric(probabilities) || any(!is.finite(probabilities)) ||
      any(probabilities < 0) || any(shape < 2L) ||
      !is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance < 0 ||
      abs(sum(probabilities) - 1) > tolerance || sum(probabilities) == 0)
    stop("Supply a probability array with at least two cells on every axis.")
  if (!is.numeric(pool_sizes) || length(pool_sizes) != length(shape) ||
      any(!is.finite(pool_sizes)) || any(pool_sizes < 1) ||
      any(pool_sizes != floor(pool_sizes)))
    stop("Supply one positive integer pool size for each count axis.")
  values <- array(probabilities / sum(probabilities), dim = shape)
  for (j in seq_along(shape)) {
    sizes <- seq_len(shape[j]) - 1L
    transform <- outer(sizes, sizes, function(r, n) {
      exp(lchoose(n, r) - lchoose(max(sizes), r))
    })
    perm <- c(j, setdiff(seq_along(shape), j))
    rotated <- aperm(values, perm)
    dim(rotated) <- c(shape[j], length(values) / shape[j])
    transformed <- array(transform %*% rotated, dim = shape[perm])
    values <- aperm(transformed, order(perm))
  }
  powers <- drop((arrayInd(seq_along(values), shape) - 1L) %*% pool_sizes)
  levels <- sort(unique(powers[powers > 0]))
  if (any(!is.finite(levels))) stop("Pool-size exponents overflow.")
  ans <- vapply(levels, function(k) {
    matches <- values[powers == k]
    if (diff(range(matches)) > tolerance)
      stop("Repeated exponent ", k, " gives inconsistent population moments.")
    mean(matches)
  }, numeric(1))
  stats::setNames(ans, format(levels, scientific = FALSE, trim = TRUE))
}

mixed_interval <- function(moments, exponents = as.numeric(names(moments)),
                           cells = 512L) {
  if (!is.numeric(moments) || !length(moments) || any(!is.finite(moments)) ||
      any(moments < 0 | moments > 1))
    stop("Moments must be finite population expectations in [0,1].")
  if (!is.numeric(exponents) || length(exponents) != length(moments) ||
      any(!is.finite(exponents)) || any(exponents != floor(exponents)) ||
      any(exponents <= 0) || anyDuplicated(exponents))
    stop("Supply a distinct positive integer exponent for each moment.")
  if (!is.numeric(cells) || length(cells) != 1L || !is.finite(cells) ||
      cells < 2 || cells != floor(cells))
    stop("cells must be an integer of at least two.")
  idx <- order(exponents)
  k <- exponents[idx]
  m <- unname(moments[idx])
  # Necessary moment inequalities follow from x^b <= x^a and Jensen.
  if (length(k) > 1L && (any(diff(m) > 1e-12) ||
      any(diff(m^(1/k)) < -1e-12)))
    stop("Moments violate monotonicity or a power-mean inequality.")
  if (length(k) == 1L)
    return(c(lower = 1 - m^(1/k), upper = 1 - m))

  grid <- sin(seq(0, pi/2, length.out = cells + 1))^2
  basis <- outer(grid, k, "^")
  # An upper bound at either endpoint on the convex interpolation error
  # throughout each adjacent cell. The power k = 1 is exactly linear.
  curvature <- matrix(0, length(grid), length(k))
  higher <- k >= 2
  cell_error <- sweep(outer(grid[-1L], k[higher] - 2, "^"), 2L,
                      k[higher] * (k[higher] - 1) / 8, "*") * diff(grid)^2
  curvature[, higher] <- pmax(rbind(0, cell_error), rbind(cell_error, 0))
  if (!all(is.finite(curvature))) stop("Exponents are too large to certify.")
  # A majorant c0 + sum(cj*x^kj) uses free signed coefficients.
  # Only the positive part of cj can dip below the chord.
  design <- cbind(1, -1, basis - curvature, -basis)
  cost <- c(1, -1, m, -m)
  d <- length(k)
  upper_bound <- function(target) {
    # A growing coefficient budget also bounds the numerical problem at
    # boundary moments, where the exact dual can have unbounded minimisers.
    # Any fixed continuous envelope fits this budget for a fine enough grid.
    active <- unique(round(seq(1, length(grid), length.out = min(25, length(grid)))))
    best <- Inf
    repeat {
      fit <- lpSolve::lp("min", cost,
                         rbind(design[active, , drop = FALSE], rep(1, length(cost))),
                         c(rep(">=", length(active)), "<="),
                         c(target[active], max(100, cells)), scale = 0)
      if (fit$status != 0 || any(!is.finite(fit$solution))) {
        if (is.finite(best)) {
          warning("Refinement failed numerically; using the last verified envelope.")
          break
        }
        stop("Could not construct a numerical envelope (solver status ", fit$status, ").")
      }
      a <- fit$solution
      intercept <- a[1] - a[2]
      beta <- a[2L + seq_len(d)] - a[2L + d + seq_len(d)]
      # Independently shift the envelope until every cell satisfies its
      # sufficient inequality, including a roundoff allowance.
      gridded <- drop(intercept + basis %*% beta - curvature %*% pmax(beta, 0))
      scale <- 1 + abs(intercept) + sum(abs(beta) * (1 + apply(curvature, 2L, max)))
      margin <- 64 * .Machine$double.eps * (d + 2) * scale
      residual <- target - gridded
      bound <- intercept + sum(beta * m) + max(0, residual) + margin
      best <- min(best, bound)
      # Add the worst constraint from each remaining run of violations.
      residual[active] <- -Inf
      violations <- which(residual > 1e-10)
      if (!length(violations)) break
      runs <- cumsum(c(TRUE, diff(violations) > 1L))
      add <- vapply(split(violations, runs),
                    function(i) i[which.max(residual[i])], integer(1))
      active <- c(active, unname(add))
    }
    best
  }
  lo <- max(0, -upper_bound(grid - 1))
  hi <- min(1, upper_bound(1 - grid))
  if (!is.finite(lo) || !is.finite(hi) || lo > hi)
    stop("Inconsistent envelope bounds. Check moment feasibility and numerical scale.")
  if (k[1] == 1) return(c(lower = 1 - m[1], upper = 1 - m[1]))
  c(lower = lo, upper = hi)
}
