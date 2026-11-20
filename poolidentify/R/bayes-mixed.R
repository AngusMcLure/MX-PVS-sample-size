# Use negative counts on every axis, as in joint_negative_moments().
mixed_bayes_prevalence <- function(histogram, pool_sizes, draws = 1000L,
                                   burn = 500L, chains = 4L, concentration = 1,
                                   atom_rate = 1, level = 0.95) {
  dimensions <- dim(histogram)
  if (is.null(dimensions)) dimensions <- length(histogram)
  if (!is.numeric(histogram) || any(dimensions < 2L) ||
      any(!is.finite(histogram)) || any(histogram < 0 | histogram != floor(histogram)) ||
      sum(histogram) < 1 || sum(histogram) > .Machine$integer.max)
    stop("Supply an integer location histogram indexed by negative counts with a positive total.")
  if (!is.numeric(pool_sizes) || length(pool_sizes) != length(dimensions) ||
      any(!is.finite(pool_sizes)) || any(pool_sizes < 1 | pool_sizes != floor(pool_sizes)))
    stop("Supply one positive integer pool size for each histogram axis.")
  good_number <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0
  good_count <- function(x, minimum) is.numeric(x) && length(x) == 1L &&
    is.finite(x) && x >= minimum && x == floor(x) && x <= .Machine$integer.max
  if (!good_count(draws, 8L) || !good_count(burn, 0L) || !good_count(chains, 1L) ||
      !good_number(concentration) || !good_number(atom_rate) ||
      !good_number(level) || level >= 1) stop("Invalid simulation, prior or interval settings.")
  cells <- arrayInd(seq_along(histogram), dimensions) - 1L
  exponents <- sort(unique(drop(cells %*% pool_sizes)))
  powers <- exponents[exponents > 0]
  # The lower principal representation of the augmented powers, including
  # x itself, represents every feasible (count law, prevalence) pair.
  atoms <- 1L + length(unique(c(1, powers))) %/% 2L
  used <- which(histogram > 0)
  freq <- as.numeric(histogram[used])
  negative <- cells[used, , drop = FALSE]
  positive <- sweep(-negative, 2L, dimensions - 1L, "+")
  small <- min(pool_sizes)
  ratios <- pool_sizes / small
  weighted_negative <- drop(negative %*% ratios)
  evaluate <- function(t, neg, pos) {
    # t = -log((1-P)^small); the multinomial allocations remove constants.
    -neg * t + sum(pos * log(-expm1(-ratios * t)))
  }
  chain_draws <- function(id) {
    rates <- stats::rexp(atoms, atom_rate)
    weights <- rep(1 / atoms, atoms)
    retained <- matrix(NA_real_, draws, 2L + length(powers))
    for (iteration in seq_len(burn + draws)) {
      lp_positive <- log(-expm1(-outer(ratios, rates)))
      log_mass <- sweep(positive %*% lp_positive - outer(weighted_negative, rates),
                        2L, log(weights), "+")
      allocations <- t(vapply(seq_along(freq), function(i) {
        mass <- log_mass[i, ]
        as.numeric(stats::rmultinom(1L, freq[i], exp(mass - max(mass))))
      }, numeric(atoms)))
      dim(allocations) <- c(length(freq), atoms)
      nk <- colSums(allocations)
      weights <- stats::rgamma(atoms, nk + concentration / atoms)
      weights <- weights / sum(weights)
      neg <- drop(crossprod(allocations, weighted_negative))
      pos <- crossprod(allocations, positive)
      for (k in seq_len(atoms)) {
        if (sum(pos[k, ]) == 0) {
          rates[k] <- stats::rexp(1L, atom_rate + neg[k])
        } else {
          # Under the exp(atom_rate) prior, exp(-atom_rate*t) is uniform.
          height <- evaluate(rates[k], neg[k], pos[k, ]) - stats::rexp(1L)
          current <- exp(-atom_rate * rates[k])
          left <- 0
          right <- 1
          repeat {
            proposal <- stats::runif(1L, left, right)
            proposed_rate <- -log(proposal) / atom_rate
            ll <- evaluate(proposed_rate, neg[k], pos[k, ])
            if (is.finite(ll) && ll >= height) break
            if (proposal < current) left <- proposal else right <- proposal
            if (right <= left) stop("Atom slice bracket exhausted at machine precision.")
          }
          rates[k] <- proposed_rate
        }
      }
      if (iteration > burn)
        retained[iteration - burn, ] <- c(sum(weights * -expm1(-rates / small)),
                                         1 / sum(weights^2),
                                         drop(exp(-outer(powers / small, rates)) %*% weights))
    }
    colnames(retained) <- c("prevalence", "effective_atoms", paste0("m", powers))
    data.frame(chain = id, retained)
  }
  posterior <- do.call(rbind, lapply(seq_len(chains), chain_draws))
  tail_prob <- (1 - level) / 2
  ci <- stats::quantile(posterior$prevalence, c(tail_prob, 1 - tail_prob), names = FALSE)
  list(draws = posterior, mean = mean(posterior$prevalence),
       interval = stats::setNames(ci, c("lower", "upper")), atoms = atoms,
       diagnostics = t(vapply(posterior[-1L], .pi_bayes_diag, numeric(2L),
                              chain = posterior$chain)),
       concentration = concentration, atom_rate = atom_rate, level = level)
}
