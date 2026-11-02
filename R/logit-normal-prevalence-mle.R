# Logit-normal pooled-data likelihood and exact exceedance probabilities.
#
# Location-specific prevalence follows logit(p_j) ~ N(mu, sigma^2).
# Each location contributes v independent pools conditional on p_j.

logit_normal_pool_count_probabilities <- function(
    mu, sigma, s, v, n_quad = 120L
) {
  stopifnot(
    length(mu) == 1L, is.finite(mu),
    length(sigma) == 1L, is.finite(sigma), sigma >= 0,
    length(s) == 1L, is.finite(s), s >= 1, s == as.integer(s),
    length(v) == 1L, is.finite(v), v >= 1, v == as.integer(v),
    length(n_quad) == 1L, is.finite(n_quad),
    n_quad >= 2, n_quad == as.integer(n_quad)
  )

  s <- as.integer(s)
  v <- as.integer(v)

  if (sigma == 0) {
    p <- plogis(mu)
    q <- -expm1(s * log1p(-p))
    return(dbinom(0:v, size = v, prob = q))
  }

  quadrature <- statmod::gauss.quad.prob(
    n = as.integer(n_quad),
    dist = "normal"
  )
  p <- plogis(mu + sigma * quadrature$nodes)
  q <- -expm1(s * log1p(-p))

  conditional <- vapply(
    q,
    function(q_j) dbinom(0:v, size = v, prob = q_j),
    numeric(v + 1L)
  )
  probabilities <- drop(conditional %*% quadrature$weights)
  probabilities / sum(probabilities)
}


# Derivative of each marginal location-count probability with respect to mu.
logit_normal_pool_count_score <- function(
    mu, sigma, s, v, n_quad = 120L
) {
  stopifnot(
    length(mu) == 1L, is.finite(mu),
    length(sigma) == 1L, is.finite(sigma), sigma >= 0,
    length(s) == 1L, is.finite(s), s >= 1, s == as.integer(s),
    length(v) == 1L, is.finite(v), v >= 1, v == as.integer(v),
    length(n_quad) == 1L, is.finite(n_quad),
    n_quad >= 2, n_quad == as.integer(n_quad)
  )

  s <- as.integer(s)
  v <- as.integer(v)

  if (sigma == 0) {
    p <- plogis(mu)
    q <- -expm1(s * log1p(-p))
    probabilities <- dbinom(0:v, size = v, prob = q)
    q_prime <- s * p * (1 - p)^s
    derivative <- probabilities *
      ((0:v) * q_prime / q - (v - 0:v) * s * p)
  } else {
    quadrature <- statmod::gauss.quad.prob(
      n = as.integer(n_quad),
      dist = "normal"
    )
    p <- plogis(mu + sigma * quadrature$nodes)
    q <- -expm1(s * log1p(-p))
    one_minus_q <- exp(s * log1p(-p))
    q_prime <- s * p * one_minus_q

    conditional <- vapply(
      seq_along(q),
      function(i) {
        k <- 0:v
        mass <- dbinom(k, size = v, prob = q[i])
        conditional_score <- k * q_prime[i] / q[i] -
          (v - k) * s * p[i]
        mass * conditional_score
      },
      numeric(v + 1L)
    )
    probabilities <- drop(
      vapply(
        q,
        function(q_j) dbinom(0:v, size = v, prob = q_j),
        numeric(v + 1L)
      ) %*% quadrature$weights
    )
    derivative <- drop(conditional %*% quadrature$weights)
  }

  probabilities <- probabilities / sum(probabilities)
  derivative <- derivative - probabilities * sum(derivative)
  derivative / probabilities
}


logit_normal_grand_mean <- function(mu, sigma, n_quad = 120L) {
  stopifnot(
    length(mu) == 1L, is.finite(mu),
    length(sigma) == 1L, is.finite(sigma), sigma >= 0
  )
  if (sigma == 0) return(plogis(mu))

  quadrature <- statmod::gauss.quad.prob(
    n = as.integer(n_quad),
    dist = "normal"
  )
  sum(quadrature$weights * plogis(mu + sigma * quadrature$nodes))
}


logit_normal_mu_from_grand_mean <- function(
    prevalence, sigma, n_quad = 120L
) {
  stopifnot(
    length(prevalence) == 1L, is.finite(prevalence),
    prevalence > 0, prevalence < 1,
    length(sigma) == 1L, is.finite(sigma), sigma >= 0
  )

  if (sigma == 0) return(qlogis(prevalence))

  stats::uniroot(
    function(mu) {
      logit_normal_grand_mean(mu, sigma, n_quad) - prevalence
    },
    interval = c(-40, 40),
    tol = 1e-10
  )$root
}


# Exact probability that the MLE of the marginal mean prevalence exceeds p_t.
#
# With sigma fixed, the MLE of the grand mean is the logit-normal marginal
# mean evaluated at the MLE of mu. The log-likelihood is concave in mu, so
# the MLE exceeds p_t exactly when the likelihood score at the mu corresponding
# to p_t is positive.
#
# The data consist of w independent location-level counts X_j in {0,...,v}.
# Their histogram has a multinomial distribution. For v=1 or v=2, the exact
# exceedance probability has an efficient binomial-tail calculation. For
# larger v, all multinomial histograms are enumerated, subject to a safeguard
# on the number of configurations.
prob_estimate_above_threshold_clustered <- function(
    mu, sigma, p_t, s, v, w,
    n_quad = 120L,
    max_configurations = 1e6
) {
  stopifnot(
    length(mu) == 1L, is.finite(mu),
    length(sigma) == 1L, is.finite(sigma), sigma >= 0,
    length(p_t) == 1L, is.finite(p_t), p_t >= 0, p_t < 1,
    length(s) == 1L, is.finite(s), s >= 1, s == as.integer(s),
    length(v) == 1L, is.finite(v), v >= 1, v == as.integer(v),
    length(w) == 1L, is.finite(w), w >= 1, w == as.integer(w),
    length(max_configurations) == 1L, is.finite(max_configurations),
    max_configurations >= 1
  )

  s <- as.integer(s)
  v <- as.integer(v)
  w <- as.integer(w)

  true_probabilities <- logit_normal_pool_count_probabilities(
    mu, sigma, s, v, n_quad
  )

  # A zero threshold is exceeded unless every pool is negative.
  if (p_t == 0) {
    return(1 - true_probabilities[1]^w)
  }

  # With no location heterogeneity, this is the ordinary pooled MLE.
  if (sigma == 0) {
    q <- -expm1(s * log1p(-plogis(mu)))
    q_t <- -expm1(s * log1p(-p_t))
    return(pbinom(floor(w * v * q_t), size = w * v,
                  prob = q, lower.tail = FALSE))
  }

  mu_threshold <- logit_normal_mu_from_grand_mean(
    p_t, sigma, n_quad
  )
  score <- logit_normal_pool_count_score(
    mu_threshold, sigma, s, v, n_quad
  )

  # The score contributions should increase with the number of positive
  # pools. Retain a general exact enumeration fallback if numerical
  # integration does not preserve that ordering.
  increasing_score <- all(diff(score) > 0)

  if (increasing_score && v == 1L) {
    cutoff <- floor(-w * score[1] / (score[2] - score[1]))
    return(pbinom(
      cutoff, size = w, prob = true_probabilities[2],
      lower.tail = FALSE
    ))
  }

  if (increasing_score && v == 2L) {
    # Let n1 be the number of locations with one positive pool. Conditional
    # on n1, the number n2 with two positive pools is binomial.
    a0 <- score[1]
    d1 <- score[2] - a0
    d2 <- score[3] - a0
    n1 <- 0:w
    remaining <- w - n1
    cutoff_n2 <- floor((-w * a0 - n1 * d1) / d2)
    conditional_probability <- true_probabilities[3] /
      (true_probabilities[1] + true_probabilities[3])
    conditional_tail <- if (conditional_probability <= 0) {
      as.numeric(cutoff_n2 < 0)
    } else if (conditional_probability >= 1) {
      as.numeric(cutoff_n2 < remaining)
    } else {
      pbinom(
        cutoff_n2, size = remaining,
        prob = conditional_probability, lower.tail = FALSE
      )
    }
    return(sum(
      dbinom(n1, size = w, prob = true_probabilities[2]) *
        conditional_tail
    ))
  }

  n_configurations <- exp(lchoose(w + v, v))
  if (!is.finite(n_configurations) ||
      n_configurations > max_configurations) {
    stop(
      "Exact clustered-MLE enumeration requires ",
      format(n_configurations, scientific = FALSE, digits = 4),
      " location-count histograms, exceeding the limit of ",
      format(max_configurations, scientific = FALSE), ". ",
      "Reduce the number of pools per trapping location or the number ",
      "of trapping locations.",
      call. = FALSE
    )
  }

  # Enumerate multinomial histograms recursively. Sequential conditional
  # binomial probabilities avoid factorial overflow and underflow.
  total_probability <- 0
  recurse <- function(category, remaining, log_probability, score_sum,
                      remaining_mass) {
    if (category == v) {
      last_probability <- true_probabilities[v + 1L] / remaining_mass
      if (score_sum + remaining * score[v + 1L] > 0) {
        total_probability <<- total_probability +
          exp(log_probability + dbinom(
            remaining, size = remaining,
            prob = last_probability, log = TRUE
          ))
      }
      return(invisible(NULL))
    }

    category_probability <- true_probabilities[category] /
      remaining_mass
    counts <- 0:remaining
    log_mass <- dbinom(
      counts, size = remaining,
      prob = category_probability, log = TRUE
    )
    for (i in seq_along(counts)) {
      count <- counts[i]
      if (!is.finite(log_mass[i])) next
      recurse(
        category + 1L,
        remaining - count,
        log_probability + log_mass[i],
        score_sum + count * score[category],
        remaining_mass - true_probabilities[category]
      )
    }
    invisible(NULL)
  }

  recurse(1L, w, 0, 0, 1)
  min(1, max(0, total_probability))
}
