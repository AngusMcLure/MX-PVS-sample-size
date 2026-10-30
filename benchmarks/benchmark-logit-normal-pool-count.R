# Benchmark the FFT and direct-convolution implementations.
#
# Run from the repository root:
#   Rscript benchmarks/benchmark-logit-normal-pool-count.R
#
# This script performs numerical consistency checks before timing. It does not
# modify files or run automatically as part of package tests.

if (!requireNamespace("statmod", quietly = TRUE)) {
  stop("Install 'statmod' first: install.packages('statmod')")
}
if (!requireNamespace("bench", quietly = TRUE)) {
  stop("Install 'bench' first: install.packages('bench')")
}

source("R/logit-normal-pool-count-pmf.R")

# Keep the parameters fixed while increasing the number of sites. The direct
# method is deliberately benchmarked only at moderate sizes because its cost
# grows quadratically in the total number of pools.
v <- 10L
s <- 25L
mu <- qlogis(0.01)
sigma <- 0.8
n_quad <- 40L
site_counts <- c(5L, 10L, 20L, 40L)
tolerance <- 1e-9

benchmark_results <- vector("list", length(site_counts))
names(benchmark_results) <- paste0("w_", site_counts)

for (i in seq_along(site_counts)) {
  w <- site_counts[[i]]
  total_pools <- w * v

  message("Checking and benchmarking w = ", w,
          " (", total_pools, " total pools)")

  fft_pmf <- logit_normal_pool_count_pmf(
    w = w, v = v, s = s, mu = mu, sigma = sigma, n_quad = n_quad
  )
  direct_pmf <- logit_normal_pool_count_pmf_direct(
    w = w, v = v, s = s, mu = mu, sigma = sigma, n_quad = n_quad
  )

  difference <- max(abs(unname(fft_pmf) - unname(direct_pmf)))
  if (!isTRUE(all.equal(
    unname(fft_pmf), unname(direct_pmf), tolerance = tolerance
  ))) {
    stop(
      "PMF mismatch for w = ", w,
      "; maximum absolute difference = ", format(difference, scientific = TRUE)
    )
  }

  # Check that the CDF agrees with cumulative probabilities from the PMF.
  # The CDF API returns P(X < c), so integer cutoff c maps to cumsum(PMF)[c].
  cutoffs <- unique(pmax(1L, pmin(total_pools, c(1L, 3L, floor(total_pools / 10))))
  cdf_values <- logit_normal_pool_count_cdf(
    w = w, v = v, s = s, mu = mu, sigma = sigma,
    c = cutoffs, n_quad = n_quad
  )
  expected_cdf <- vapply(
    cutoffs,
    function(cutoff) sum(direct_pmf[seq_len(cutoff)]),
    numeric(1)
  )
  cdf_difference <- max(abs(unname(cdf_values) - expected_cdf))
  if (!isTRUE(all.equal(
    unname(cdf_values), expected_cdf, tolerance = tolerance
  ))) {
    stop(
      "CDF mismatch for w = ", w,
      "; maximum absolute difference = ",
      format(cdf_difference, scientific = TRUE)
    )
  }

  timing <- bench::mark(
    fft_pmf = logit_normal_pool_count_pmf(
      w = w, v = v, s = s, mu = mu, sigma = sigma, n_quad = n_quad
    ),
    direct_pmf = logit_normal_pool_count_pmf_direct(
      w = w, v = v, s = s, mu = mu, sigma = sigma, n_quad = n_quad
    ),
    truncated_cdf = logit_normal_pool_count_cdf(
      w = w, v = v, s = s, mu = mu, sigma = sigma,
      c = max(1L, floor(total_pools / 10)), n_quad = n_quad
    ),
    iterations = 3,
    check = FALSE,
    memory = FALSE
  )

  timing$w <- w
  timing$total_pools <- total_pools
  timing$max_abs_pmf_difference <- difference
  timing$max_abs_cdf_difference <- cdf_difference
  benchmark_results[[i]] <- timing
}

results <- do.call(rbind, benchmark_results)
print(results[, c(
  "expression", "w", "total_pools", "median", "itr/sec",
  "max_abs_pmf_difference", "max_abs_cdf_difference"
)])

message("All PMF and CDF consistency checks passed.")
