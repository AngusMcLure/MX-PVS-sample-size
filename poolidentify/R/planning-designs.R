.planning_sizes <- function(x) {
  if (!is.numeric(x) || !is.null(dim(x)) || !length(x) ||
      any(!is.finite(x)) || any(x < 1 | x > .Machine$integer.max | x != floor(x)))
    stop("Pool sizes must be positive integers.")
  as.integer(x)
}

candidate_pool_designs <- function(pool_sizes, max_pools, equal_only = FALSE,
                                   max_designs = 10000L) {
  sizes <- sort(unique(.planning_sizes(pool_sizes)))
  if (!.pi_positive_integer(max_pools) || max_pools > .Machine$integer.max)
    stop("max_pools must be a positive integer.")
  if (!is.logical(equal_only) || length(equal_only) != 1L || is.na(equal_only))
    stop("equal_only must be TRUE or FALSE.")
  if (!.pi_positive_integer(max_designs))
    stop("max_designs must be a positive integer.")
  n <- length(sizes)
  count <- if (equal_only || n == 1L) as.double(n) * max_pools else
    choose(n + as.double(max_pools), n) - 1
  if (!is.finite(count) || count > max_designs)
    stop("Candidate count exceeds max_designs. Narrow the search or increase the limit.")
  out <- vector("list", count)
  i <- 0L
  for (v in seq_len(max_pools)) {
    if (equal_only || n == 1L) {
      for (s in sizes) {
        i <- i + 1L
        out[[i]] <- rep.int(s, v)
      }
    } else {
      # Advance the last index that can increase; the suffix then repeats
      # its new value. Each nondecreasing tuple is visited once.
      indices <- rep.int(1L, v)
      repeat {
        i <- i + 1L
        out[[i]] <- sizes[indices]
        movable <- which(indices < n)
        if (!length(movable)) break
        last <- movable[length(movable)]
        indices[last:v] <- indices[last] + 1L
      }
    }
  }
  out
}

# Joint counts of negative pools, by distinct size in ascending order.
# Integrate the joint conditional distribution with a single common P per site.
planning_count_law <- function(design, prevalence, weights = rep(1, length(prevalence))) {
  design <- .planning_sizes(design)
  if (!is.numeric(prevalence) || !is.null(dim(prevalence)) || !length(prevalence) ||
      any(!is.finite(prevalence)) || any(prevalence < 0 | prevalence > 1))
    stop("Site prevalences must lie in [0, 1].")
  if (!is.numeric(weights) || !is.null(dim(weights)) ||
      length(weights) != length(prevalence) || any(!is.finite(weights)) ||
      any(weights < 0) || !any(weights > 0))
    stop("Supply finite, nonnegative site weights with positive total mass.")
  positive_mass <- weights > 0
  weights <- weights[positive_mass] / max(weights)
  weights <- weights / sum(weights)
  log_uninfected <- log1p(-prevalence[positive_mass])
  sizes <- sort(unique(design))
  pools <- tabulate(match(design, sizes))
  shape <- pools + 1L
  cells <- prod(shape)
  if (cells > .Machine$integer.max)
    stop("Joint count array is too large; use fewer distinct sizes or pools.")
  indices <- arrayInd(seq_len(cells), shape)
  mass <- numeric(cells)
  # Limit temporary matrices when many support points approximate a continuous law.
  block_size <- max(1, floor(100000 / cells))
  for (first in seq(1, length(weights), by = block_size)) {
    block <- first:min(first + block_size - 1, length(weights))
    joint <- matrix(1, cells, length(block))
    for (a in seq_along(sizes)) {
      log_neg <- sizes[a] * log_uninfected[block]
      # Evaluate the smaller binomial tail probability directly. Near p=0,
      # forming the positive probability by subtraction can erase it entirely.
      use_positive <- log_neg > -log(2)
      rate <- ifelse(use_positive, -expm1(log_neg), exp(log_neg))
      k <- rep.int(0:pools[a], length(block))
      k <- ifelse(rep(use_positive, each = shape[a]), pools[a] - k, k)
      factors <- matrix(stats::dbinom(k, pools[a], rep(rate, each = shape[a])),
                        nrow = shape[a])
      joint <- joint * factors[indices[, a], , drop = FALSE]
    }
    mass <- mass + drop(joint %*% weights[block])
  }
  array(mass / sum(mass), dim = shape)
}
