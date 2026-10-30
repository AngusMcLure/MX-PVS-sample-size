# Gauss-Hermite quadrature for expectations under a standard normal
# distribution. The Golub-Welsch eigenvalue algorithm avoids an external
# package dependency for the clustered pool-count calculations.
.gauss_hermite_normal <- function(n) {
  if (length(n) != 1L || !is.finite(n) || n < 2 || n != as.integer(n)) {
    stop("n must be an integer greater than or equal to 2.", call. = FALSE)
  }
  n <- as.integer(n)

  jacobi <- matrix(0, nrow = n, ncol = n)
  off_diagonal <- sqrt(seq_len(n - 1L))
  jacobi[cbind(seq_len(n - 1L), 2:n)] <- off_diagonal
  jacobi[cbind(2:n, seq_len(n - 1L))] <- off_diagonal

  decomposition <- eigen(jacobi, symmetric = TRUE)
  ordering <- order(decomposition$values)
  weights <- decomposition$vectors[1L, ordering]^2
  weights <- weights / sum(weights)

  list(nodes = decomposition$values[ordering], weights = weights)
}


# Match the first two moments of a logit-normal distribution using only base R.
# For each sigma, solve for mu to match the requested mean, then solve for
# sigma to match the requested variance.
logit_normal_parameters_from_moments <- function(mean, variance) {
  if (length(mean) != 1L || !is.finite(mean) || mean <= 0 || mean >= 1) {
    stop("mean must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  if (length(variance) != 1L || !is.finite(variance) || variance < 0 ||
      variance >= mean * (1 - mean)) {
    stop("variance must be in [0, mean * (1 - mean)).", call. = FALSE)
  }
  if (variance == 0) {
    return(c(mu = qlogis(mean), sigma = 0))
  }

  normal_moments <- function(mu, sigma) {
    mean_value <- stats::integrate(
      function(z) plogis(mu + sigma * z) * stats::dnorm(z),
      lower = -10, upper = 10,
      subdivisions = 1000L, rel.tol = 1e-9, abs.tol = 1e-11
    )$value
    second_moment <- stats::integrate(
      function(z) plogis(mu + sigma * z)^2 * stats::dnorm(z),
      lower = -10, upper = 10,
      subdivisions = 1000L, rel.tol = 1e-9, abs.tol = 1e-11
    )$value
    c(mean = mean_value, variance = max(0, second_moment - mean_value^2))
  }

  mu_for_sigma <- function(sigma) {
    centre <- qlogis(mean)
    span <- 10 + 10 * sigma
    stats::uniroot(
      function(mu) normal_moments(mu, sigma)[["mean"]] - mean,
      interval = c(centre - span, centre + span),
      tol = 1e-10
    )$root
  }

  variance_for_sigma <- function(sigma) {
    mu <- mu_for_sigma(sigma)
    normal_moments(mu, sigma)[["variance"]]
  }

  upper <- 1
  upper_variance <- variance_for_sigma(upper)
  while (upper_variance < variance && upper < 64) {
    upper <- upper * 2
    upper_variance <- variance_for_sigma(upper)
  }
  if (upper_variance < variance) {
    stop("Could not find logit-normal parameters for the requested moments.",
         call. = FALSE)
  }

  sigma <- stats::uniroot(
    function(value) variance_for_sigma(value) - variance,
    interval = c(0, upper),
    tol = 1e-8
  )$root
  mu <- mu_for_sigma(sigma)

  c(mu = mu, sigma = sigma)
}


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
    quadrature <- .gauss_hermite_normal(n_quad)

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


#' Cumulative distribution function for the clustered pool count
#'
#' Calculates P(X < c) for the total number of positive pools under
#' logit-normal trapping-site heterogeneity.
#'
#' Unlike summing the FFT-derived PMF, this function calculates lower-tail
#' probabilities by repeated truncated convolution in probability space.
#' This avoids FFT cancellation and does not calculate probabilities above
#' the requested cutoff.
#'
#' @inheritParams logit_normal_pool_count_pmf
#' @param c Numeric cutoff or vector of cutoffs. Returns P(X < c); since X
#'   is integer-valued, this is P(X <= ceiling(c) - 1).
#'
#' @return A numeric vector of probabilities P(X < c), named by the supplied
#'   cutoff values.
#'
#' @examples
#' logit_normal_pool_count_cdf(
#'   w = 30, v = 20, s = 25, mu = qlogis(0.01), sigma = 0.8,
#'   c = c(5, 10, 15)
#' )
#'
#' @export
logit_normal_pool_count_cdf <- function(
    w, v, s, mu, sigma,
    c,
    n_quad = 40L
) {
  check_scalar_integer <- function(value, name, minimum = 1) {
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value < minimum ||
        value != floor(value)) {
      stop(name, " must be a single integer >= ", minimum, ".",
           call. = FALSE)
    }
    as.integer(value)
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
  if (!is.numeric(c) || length(c) == 0L || any(!is.finite(c))) {
    stop("c must be a non-empty vector of finite numeric cutoffs.",
         call. = FALSE)
  }

  max_count <- as.double(w) * v
  if (max_count > .Machine$integer.max - 1) {
    stop("w*v is too large to represent the distribution as an R vector.",
         call. = FALSE)
  }
  max_count <- as.integer(max_count)
  cutoffs <- ceiling(c) - 1

  # Boundary cases can be answered exactly without any distribution work.
  answer <- numeric(length(c))
  answer[cutoffs < 0] <- 0
  answer[cutoffs >= max_count] <- 1
  needs_calculation <- which(cutoffs >= 0 & cutoffs < max_count)

  if (length(needs_calculation) > 0L) {
    max_cutoff <- as.integer(max(cutoffs[needs_calculation]))

    if (sigma == 0) {
      p <- plogis(mu)
      q <- -expm1(s * log1p(-p))
      # pbinom is designed to evaluate binomial tails stably.
      calculated <- pbinom(
        cutoffs[needs_calculation],
        size = max_count,
        prob = q
      )
      answer[needs_calculation] <- calculated

    } else {
      if (!requireNamespace("statmod", quietly = TRUE)) {
        stop(
          "Package 'statmod' is required for Gauss-Hermite quadrature. ",
          "Install it with install.packages('statmod').",
          call. = FALSE
        )
      }

      quadrature <- statmod::gauss.quad.prob(
        n = n_quad,
        dist = "normal"
      )
      local_prevalence <- plogis(mu + sigma * quadrature$nodes)
      pool_positive <- -expm1(s * log1p(-local_prevalence))

      conditional_pmfs <- vapply(
        pool_positive,
        function(q_j) dbinom(0:v, size = v, prob = q_j),
        numeric(v + 1L)
      )
      site_pmf <- drop(conditional_pmfs %*% quadrature$weights)
      site_pmf <- site_pmf / sum(site_pmf)

      # Truncated convolution: retain only probabilities for counts
      # 0,...,max_cutoff at every step. All terms are non-negative, so this
      # avoids the cancellation that can affect FFT-based lower tails.
      accumulated <- 1
      for (site in seq_len(w)) {
        current_max <- length(accumulated) - 1L
        next_max <- min(max_cutoff, current_max + v)
        next_pmf <- numeric(next_max + 1L)

        for (k in 0:next_max) {
          j_min <- max(0L, k - current_max)
          j_max <- min(v, k)
          js <- j_min:j_max
          next_pmf[k + 1L] <- sum(
            accumulated[k - js + 1L] * site_pmf[js + 1L]
          )
        }

        accumulated <- next_pmf
      }

      # Use cumulative sums of the truncated, non-negative PMF.
      cumulative <- cumsum(accumulated)
      answer[needs_calculation] <- cumulative[
        cutoffs[needs_calculation] + 1L
      ]
    }
  }

  # Guard against tiny rounding excursions outside [0, 1].
  answer <- pmin(1, pmax(0, answer))
  names(answer) <- paste0("P(X<", c, ")")
  answer
}



#' Direct-convolution reference implementation of the clustered pool-count PMF
#'
#' Calculates the same distribution as [logit_normal_pool_count_pmf()] using
#' ordinary discrete convolution rather than FFTs. This implementation is
#' intended as a reference for numerical checks and benchmarking; it is much
#' slower for large distributions.
#'
#' @inheritParams logit_normal_pool_count_pmf
#' @return A numeric vector containing P(X = k), for k = 0, ..., w*v.
#' @keywords internal
logit_normal_pool_count_pmf_direct <- function(
    w, v, s, mu, sigma,
    n_quad = 40L
) {
  check_scalar_integer <- function(value, name, minimum = 1) {
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value < minimum ||
        value != floor(value)) {
      stop(name, " must be a single integer >= ", minimum, ".",
           call. = FALSE)
    }
    as.integer(value)
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

  if (sigma == 0) {
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

    quadrature <- .gauss_hermite_normal(n_quad)
    local_prevalence <- plogis(mu + sigma * quadrature$nodes)
    pool_positive <- -expm1(s * log1p(-local_prevalence))

    conditional_pmfs <- vapply(
      pool_positive,
      function(q_j) dbinom(0:v, size = v, prob = q_j),
      numeric(v + 1L)
    )
    site_pmf <- drop(conditional_pmfs %*% quadrature$weights)
    site_pmf <- site_pmf / sum(site_pmf)

    # Ordinary O(a*b) discrete convolution, implemented explicitly so this
    # reference path does not call an FFT-based convolution routine.
    total_pmf <- 1
    for (site in seq_len(w)) {
      current_max <- length(total_pmf) - 1L
      next_max <- current_max + v
      next_pmf <- numeric(next_max + 1L)

      for (k in 0:next_max) {
        j_min <- max(0L, k - current_max)
        j_max <- min(v, k)
        js <- j_min:j_max
        next_pmf[k + 1L] <- sum(
          total_pmf[k - js + 1L] * site_pmf[js + 1L]
        )
      }
      total_pmf <- next_pmf
    }

    total_mass <- sum(total_pmf)
    if (!is.finite(total_mass) || total_mass <= 0) {
      stop("Numerical failure when calculating the total PMF.",
           call. = FALSE)
    }
    total_pmf <- total_pmf / total_mass
  }

  names(total_pmf) <- as.character(0:max_count)
  total_pmf
}
