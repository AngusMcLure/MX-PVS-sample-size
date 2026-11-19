# Bernstein probabilities are ordered by the number of *positive* pools.
.pi_count_law <- function(prob) {
  if (!is.numeric(prob) || !is.null(dim(prob)) || length(prob) < 2L ||
      any(!is.finite(prob)) || any(prob < 0) || abs(sum(prob) - 1) > 1e-10) {
    stop("Supply probabilities for counts 0 through v, summing to one.")
  }
  prob / sum(prob)
}

count_moments <- function(prob) {
  prob <- .pi_count_law(prob)
  v <- length(prob) - 1L
  vapply(seq_len(v), function(k) {
    r <- k:v
    sum(prob[r + 1L] * exp(lchoose(r, k) - lchoose(v, k)))
  }, numeric(1))
}

# E[q^a(1-q)^b] from the v-pool count law by sampling a+b
# distinct pools without replacement. No alternating sums are needed.
.pi_product <- function(prob, a, b) {
  v <- length(prob) - 1L
  r <- a:(v - b)
  log_prob <- lchoose(r, a) + lchoose(v - r, b) -
    lchoose(v, a + b) - lchoose(a + b, a)
  sum(prob[r + 1L] * exp(log_prob))
}

.pi_moments_to_law <- function(moments) {
  if (!is.numeric(moments) || !is.null(dim(moments)) || !length(moments) ||
      any(!is.finite(moments)) || any(moments < 0 | moments > 1)) {
    stop("Supply finite moments E[Q^k], k=1,...,v, in [0,1].")
  }
  v <- length(moments)
  differences <- c(1, moments)
  law <- numeric(v + 1L)
  for (b in 0:v) {
    law[v - b + 1L] <- choose(v, b) * differences[v - b + 1L]
    differences <- -diff(differences)
  }
  if (any(!is.finite(law)) || any(law < -1e-12) || abs(sum(law) - 1) > 1e-10) {
    stop("Infeasible or ill-conditioned raw moments. Use count probabilities if available.")
  }
  law <- pmax(law, 0)
  law / sum(law)
}

# The inverse pooling link has a cusp at q=1. Retain complementary
# coordinates supplied by the eigensolver when q is close to one.
.pi_unpool <- function(q, negative, pool_size) {
  log_negative <- ifelse(q < 0.5, log1p(-pmin(1, q)), log(negative))
  -expm1(log_negative / pool_size)
}

.pi_second_bounds <- function(prob, s) {
  second <- prob[3L]
  second_negative <- prob[1L]
  cross <- prob[2L] / 2
  q <- second + cross
  negative <- second_negative + cross
  # det(E[(1-Q,Q)'(1-Q,Q)]) equals Var(Q). This expression avoids
  # subtracting two numbers close to one at very high positivity.
  variance <- second_negative * second - cross^2
  tol <- 128 * .Machine$double.eps * max(second_negative * second, cross^2)
  if (variance < -tol) stop("Infeasible count law: negative mixing variance.")
  if (variance <= 0) {
    value <- .pi_unpool(q, negative, s)
    return(c(lower = value, upper = value))
  }
  if (cross == 0) return(c(lower = q, upper = q))
  left <- q * (q / second) * .pi_unpool(second / q, cross / q, s)
  right <- variance / second_negative +
    negative * (negative / second_negative) *
      .pi_unpool(cross / negative, second_negative / negative, s)
  c(lower = left, upper = right)
}

.pi_rule <- function(prob, fixed_zero, fixed_one, s) {
  v <- length(prob) - 1L
  n <- (v + 1L - fixed_zero - fixed_one) %/% 2L
  interior <- interior_negative <- numeric()
  if (n > 0) {
    # Gauss for the measure q^fixed_zero (1-q)^fixed_one dF.
    # Use the degree n-1 Bernstein basis (up to diagonal scaling).
    powers <- 0:(2L * n - 2L)
    a <- powers + fixed_zero
    b <- rev(powers) + fixed_one
    hankel_indices <- outer(seq_len(n), seq_len(n), "+") - 1L
    q_matrix <- matrix(mapply(function(i, j) .pi_product(prob, i + 1L, j),
                               a, b)[hankel_indices], n)
    negative_matrix <- matrix(mapply(function(i, j) .pi_product(prob, i, j + 1L),
                                      a, b)[hankel_indices], n)
    gram <- q_matrix + negative_matrix
    active <- diag(gram) > 0
    if (any(active)) {
      scales <- sqrt(diag(gram)[active])
      scaled <- function(x) x[active, active, drop = FALSE] / outer(scales, scales)
      factors <- eigen(scaled(gram), symmetric = TRUE)
      tol <- 128 * n * .Machine$double.eps * max(factors$values)
      if (min(factors$values) < -tol) stop("Infeasible moments: indefinite moment matrix.")
      keep <- factors$values > tol
      transform <- sweep(factors$vectors[, keep, drop = FALSE], 2L,
                         sqrt(factors$values[keep]), "/")
      nodes <- function(x) {
        operator <- crossprod(transform, scaled(x) %*% transform)
        vals <- eigen(operator, symmetric = TRUE, only.values = TRUE)$values
        if (min(vals) < -1e-8 || max(vals) > 1 + 1e-8)
          stop("Infeasible or ill-conditioned moments: nodes outside [0,1].")
        pmin(1, pmax(0, vals))
      }
      interior <- rev(nodes(q_matrix))
      interior_negative <- nodes(negative_matrix)
    }
  }
  q <- c(if (fixed_zero) 0, interior, if (fixed_one) 1)
  negative <- c(if (fixed_zero) 1, interior_negative, if (fixed_one) 0)
  # In a boundary solution an ostensibly interior node can coincide with
  # a prescribed endpoint. Use whichever coordinate resolves that end.
  distinct <- !duplicated(data.frame(q, negative))
  q <- q[distinct]
  negative <- negative[distinct]
  if (!length(q)) stop("No nodes available for these moments.")
  k <- 0:v
  basis <- vapply(seq_along(q), function(i) {
    if (q[i] <= 0.5) stats::dbinom(k, v, q[i])
    else stats::dbinom(v - k, v, negative[i])
  }, numeric(v + 1L))
  design <- qr(basis, tol = 64 * .Machine$double.eps)
  if (design$rank < ncol(basis)) {
    cols <- sort(design$pivot[seq_len(design$rank)])
    basis <- basis[, cols, drop = FALSE]
    q <- q[cols]
    negative <- negative[cols]
  }
  weights <- drop(qr.solve(basis, prob, tol = 64 * .Machine$double.eps))
  if (min(weights) < -1e-9 || max(abs(basis %*% weights - prob)) > 1e-10) {
    stop("Infeasible or ill-conditioned moments: quadrature cannot reproduce the count law.")
  }
  weights <- pmax(0, weights)
  list(value = sum(weights * .pi_unpool(q, negative, s)) / sum(weights),
       # Losing a node gives a boundary measure. Zero weights at a
       # prescribed endpoint alone only certify the boundary in exact arithmetic.
       rank_loss = length(interior) < n,
       boundary_at_one = fixed_one &&
         ((fixed_zero && weights[1L] <= 64 * .Machine$double.eps) ||
            (!fixed_zero && length(interior) > 0 &&
               min(interior) <= 64 * .Machine$double.eps)))
}

identification_interval <- function(moments = NULL, pool_size, probabilities = NULL) {
  if (length(pool_size) != 1L || !is.numeric(pool_size) || !is.finite(pool_size) ||
      pool_size < 1 || pool_size != floor(pool_size)) stop("pool_size must be a positive integer.")
  if (is.null(moments) == is.null(probabilities))
    stop("Supply exactly one of moments or probabilities.")
  prob <- if (is.null(moments)) .pi_count_law(probabilities) else .pi_moments_to_law(moments)
  v <- length(prob) - 1L
  if (v == 1L) {
    return(c(lower = .pi_unpool(prob[2L], prob[1L], pool_size), upper = prob[2L]))
  }
  if (all(prob[2:v] == 0)) return(c(lower = prob[v + 1L], upper = prob[v + 1L]))
  if (v == 2L) {
    bounds <- .pi_second_bounds(prob, pool_size)
  } else {
    odd <- v %% 2L == 1L
    lower <- .pi_rule(prob, !odd, FALSE, pool_size)
    upper <- .pi_rule(prob, odd, TRUE, pool_size)
    bounds <- c(lower = lower$value, upper = upper$value)
    # At a rank-deficient boundary the compatible law is unique. The
    # upper rule explicitly contains q=1, where the inverse link is singular.
    if (upper$rank_loss || lower$rank_loss || upper$boundary_at_one) {
      if (abs(diff(bounds)) > 1e-7)
        warning("Numerical boundary at q=1: endpoints sensitive to rounding; ",
                "returning the boundary-law value.")
      bounds[] <- bounds["upper"]
    }
  }
  if (pool_size == 1) bounds[] <- sum((0:v) * prob) / v
  if (bounds[1L] > bounds[2L] + 1e-9)
    stop("Ill-conditioned moments: quadrature endpoints crossed.")
  stats::setNames(sort(pmin(pmax(bounds, 0), 1)), c("lower", "upper"))
}
