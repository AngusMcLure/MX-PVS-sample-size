.check_planning_grid <- function(cells, grid_size) {
  if (any(cells * (2 * as.double(grid_size) + 1) > 2e6))
    stop("Joint count grid is too large. Reduce distinct pool sizes or grid_size.")
}

# The first cell of a one-pool negative-count histogram is the number of
# positive pools. This vectorized version also permits fractional screening counts.
planning_one_width <- function(positives, locations, pool_size, level = 0.95) {
  ci <- binomial_limits(positives, locations, (1 - level) / 2)
  ci$hi + expm1(log1p(-ci$lo) / pool_size)
}

# O(sqrt(n)) significant binomial counts even near probability one half.
# Work in the rarer outcome to retain precision near either boundary.
planning_one_expected <- function(law, locations, pool_size, level = 0.95) {
  q <- min(law)
  positive_rare <- law[1L] <= law[2L]
  tail <- 4 * .Machine$double.eps
  bounds <- c(stats::qbinom(tail, locations, q),
              stats::qbinom(tail, locations, q, lower.tail = FALSE))
  counts <- seq(bounds[1L], bounds[2L])
  positives <- if (positive_rare) counts else locations - counts
  widths <- planning_one_width(positives, locations, pool_size, level)
  sum(stats::dbinom(counts, locations, q) * widths)
}

# Histograms here count NEGATIVE pools, in ascending size-class order.
# Evaluate only the CI width. The projection is still needed to make the
# confidence band nonempty on samples outside the mixture hull.
planning_width_engine <- function(design, locations, level = 0.95, grid_size = NULL) {
  sizes <- sort(unique(design))
  if (length(design) == 1L) {
    return(function(histogram) planning_one_width(histogram[1L], locations, sizes, level))
  }
  grid <- confidence_grid(locations, max(sizes), grid_size)
  if (length(sizes) == 1L) {
    solver <- equal_count_solver(length(design), sizes, locations, level, grid)
    return(function(histogram) unname(diff(solver$fit(rev(as.numeric(histogram)))$confidence)))
  }
  shape <- tabulate(match(design, sizes), nbins = length(sizes)) + 1L
  solve <- mixed_count_solver(shape, sizes, locations, level, grid)
  function(histogram) unname(diff(solve(histogram)$confidence))
}
