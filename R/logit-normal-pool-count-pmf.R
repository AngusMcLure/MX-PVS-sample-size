#' Distribution of positive pools under logit-normal site heterogeneity
#'
#' Calculates the probability mass function for the total number of positive
#' pools when pools are grouped within trapping locations and each location
#' has its own prevalence.
#'
#' @param w Number of independent trapping locations.
#' @param v Number of pools sampled at each location.
#' @param s Number of mosquitoes in each pool.
#' @param mu Mean of the location-specific logit prevalence.
#' @param sigma Standard deviation of the location-specific logit prevalence.
#' @param n_quad Number of Gauss-Hermite quadrature points.
#' @param x Optional integer count or vector of counts for which to return
#'   probabilities. If NULL, the complete PMF for counts 0, ..., w*v is
#'   returned.
#'
#' @details
#' At location j, logit(p_j) is normally distributed with mean mu and
#' standard deviation sigma. Conditional on p_j, the number of positive pools
#' at that location is binomial with size v and probability 1 - (1-p_j)^s.
#' The location-level PMF is integrated over the logit-normal distribution
#' using Gauss-Hermite quadrature from the statmod package.
#'
#' The total PMF is obtained by raising the discrete Fourier transform of the
#' location-level PMF to the power w, then applying the inverse transform.
#' Zero-padding to at least w*v + 1 elements ensures this gives linear rather
#' than circular convolution.
#'
#' @return If x is NULL, a numeric vector containing P(X = k), for
#'   k = 0, ..., w*v. Otherwise, a named numeric vector containing the
#'   requested probabilities.
#'
#' @examples
#' pmf <- logit_normal_pool_count_pmf(
#'   w = 30, v = 20, s = 25, mu = qlogis(0.01), sigma = 0.8
#' )
#' sum(pmf)
#'
#' logit_normal_pool_count_pmf(
#'   w = 30, v = 20, s = 25, mu = qlogis(0.01), sigma = 0.8,
#'   x = c(0, 5, 10)
#' )
#'
#' @export
logit_normal_pool_count_pmf <- function(
    w, v, s, mu, sigma,
    n_quad = 40L,
    x = NULL
) {
  check_scalar_integer <- function(value, name, minimum = 1) {
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value < minimum ||
        value != floor(value)) {
      stop(name, " must be a single integer >= ", minimum, ".",
           call. = FALSE)
    }
    invisible(as.integer(value))
  }

  w <- check_scalar_integer(w, "w")
  v <- check_scalar_integer(v, "v")
  s <- check_scalar_integer(s, "s")
  n_quad <- check_scalar_integer(n_quad, "n_quad", minimum = 2)

  if (!is.numeric(mu) || length(mu) != 1L || !is.finite(mu)) {
    stop("mu must be a single finite number.", call. = FALSE)
  }
  if (!is.numeric(sigma) || length(sigma) != 1L ||
      !is.finite(sigma) || sigma < 0) {
    stop("sigma must be a single finite number >= 0.", call. = FALSE)
  }

  max_count <- as.double(w) * v
  if (max_count > .Machine$integer.max - 1) {
    stop("w*v is too large to represent the PMF as an R vector.",
         call. = FALSE)
  }
  max_count <- as.integer(max_count)

  if (!is.null(x)) {
    if (!is.numeric(x) || any(!is.finite(x)) ||
        any(x != floor(x)) || any(x < 0 | x > max_count)) {
      stop("x must contain integer counts between 0 and w*v.",
           call. = FALSE)
    }
  }

  if (sigma == 0) {
    # With no site heterogeneity, all pools have the same positivity
    # probability and the total count is exactly binomial.
    p <- plogis(mu)
    q <- -expm1(s * log1p(-p))
    total_pmf <- dbinom(0:max_count, size = max_count, prob = q)

  } else {
    if (!requireNamespace("statmod", quietly = TRUE)) {
      stop(
        "Package 'statmod' is required for Gauss-Hermite quadrature. ",
        "Install it with install.packages('statmod').",
        call. = FALSE
      )
    }

    # Integrate with respect to a standard normal distribution.
    quadrature <- statmod::gauss.quad.prob(
      n = n_quad,
      dist = "normal"
    )

    local_prevalence <- plogis(mu + sigma * quadrature$nodes)

    # Probability that a pool is positive at each location.
    pool_positive <- -expm1(s * log1p(-local_prevalence))

    # Integrate the conditional Binomial(v, q_j) PMF over site effects.
    conditional_pmfs <- vapply(
      pool_positive,
      function(q_j) dbinom(0:v, size = v, prob = q_j),
      numeric(v + 1L)
    )
    site_pmf <- drop(conditional_pmfs %*% quadrature$weights)
    site_pmf <- site_pmf / sum(site_pmf)

    # The w-fold convolution is obtained by raising the DFT to power w.
    # Padding to length >= w*v + 1 prevents circular wrap-around.
    required_length <- max_count + 1
    n_fft <- 2^ceiling(log2(required_length))
    if (n_fft > .Machine$integer.max) {
      stop("The required FFT is too large for an R vector.",
           call. = FALSE)
    }
    n_fft <- as.integer(n_fft)

    padded_site_pmf <- c(
      site_pmf,
      numeric(n_fft - length(site_pmf))
    )
    transformed_pmf <- fft(padded_site_pmf)

    total_pmf <- Re(
      fft(transformed_pmf^w, inverse = TRUE)
    ) / n_fft
    total_pmf <- total_pmf[seq_len(max_count + 1L)]

    # FFT round-off can produce tiny negative values.
    if (min(total_pmf) < -1e-10) {
      warning(
        "FFT produced materially negative probabilities; ",
        "check numerical accuracy for these parameter values.",
        call. = FALSE
      )
    }
    total_pmf <- pmax(total_pmf, 0)

    total_mass <- sum(total_pmf)
    if (!is.finite(total_mass) || total_mass <= 0) {
      stop("Numerical failure when calculating the total PMF.",
           call. = FALSE)
    }
    total_pmf <- total_pmf / total_mass
  }

  names(total_pmf) <- as.character(0:max_count)

  if (is.null(x)) {
    return(total_pmf)
  }

  setNames(total_pmf[x + 1L], paste0("P(X=", x, ")"))
}
