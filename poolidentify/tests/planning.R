library(poolidentify)

# Independent exact reference: observe one test per location and integrate the
# projected Clopper-Pearson width against the number of positive locations.
exact_width <- function(w, size, p, mass, level) {
  k <- 0:w
  q <- sum(mass * (1 - (1 - p)^size))
  widths <- vapply(k, function(x) {
    ci <- binom.test(x, w, conf.level = level)$conf.int
    ci[2] - (1 - (1 - ci[1])^(1/size))
  }, numeric(1))
  sum(dbinom(k, w, q) * widths)
}
designs <- list(individual = 1, pooled = 6, unaffordable = 101)
scenarios <- list(rare = list(prevalence = c(0, .15), weights = c(.9, .1)),
                  common = list(prevalence = .35, weights = 1))
args <- list(designs = designs, budget = 50, cost_location = .75,
             cost_pool = 1.25, cost_unit = .5, scenarios = scenarios,
             level = .9, grid_size = 12)
weights <- c(.8, .2)
weighted <- do.call(rank_pool_designs, c(args, list(scenario_weights = c(common = .2, rare = .8))))
worst <- do.call(rank_pool_designs, c(args, list(criterion = "worst")))
stopifnot(setequal(weighted$ranking$design, c("individual", "pooled")),
          identical(weighted$excluded$design, "unaffordable"),
          identical(weighted$best, as.integer(designs[[weighted$ranking$design[1]]])))
# This decimal budget buys three locations despite rounding in 0.3 / 0.1.
decimal <- rank_pool_designs(list(one = 1), .3, .1, 0, 0, scenarios = scenarios)
stopifnot(decimal$ranking$locations == 3,
          inherits(try(rank_pool_designs(list(one = 1), .3, 1, 0, 0), silent = TRUE), "try-error"))
for (name in c("individual", "pooled")) {
  s <- designs[[name]]
  w <- floor(50 / (.75 + 1.25 + .5*s))
  ref <- vapply(scenarios, function(x)
    exact_width(w, s, x$prevalence, x$weights, .9), numeric(1))
  row <- subset(weighted$ranking, design == name)
  stopifnot(row$locations == w, row$pools_per_location == 1,
            row$units_per_location == s, row$total_cost == w * (2 + .5*s),
            abs(row$width - sum(ref * weights)) < 1e-10,
            abs(subset(worst$ranking, design == name)$width - max(ref)) < 1e-10)
  details <- subset(weighted$by_scenario, design == name)
  stopifnot(max(abs(details$width - ref[details$scenario])) < 1e-10,
            all(details$mcse == 0))
}

# Exhaust all histograms for one size-two and one size-five pool per site.
# This distinguishes variation across surveys from a width evaluated at the
# average histogram, and checks the negative-count axes of a mixed design.
s <- c(2, 5)
scenario <- list(prevalence = c(.03, .45), weights = c(.7, .3))
pneg <- function(p) {
  q <- (1-p)^s
  as.vector(outer(c(1-q[1], q[1]), c(1-q[2], q[2])))
}
law <- drop(vapply(scenario$prevalence, pneg, numeric(4)) %*% scenario$weights)
n <- 4
first <- second <- 0
for (i in 0:n) for (j in 0:(n-i)) for (k in 0:(n-i-j)) {
  h <- c(i, j, k, n-i-j-k)
  width <- diff(mixed_confidence(matrix(h, 2), s, level = .8, grid_size = 8)$confidence)
  prob <- dmultinom(h, prob = law)
  first <- first + prob * width
  second <- second + prob * width^2
}
set.seed(719)
runs <- 128
ranked <- rank_pool_designs(list(mixed = s), budget = 44, cost_location = 2,
                            cost_pool = 1, cost_unit = 1,
                            scenarios = list(planning = scenario), level = .8,
                            replicates = runs, grid_size = 8)
truth_se <- sqrt((second - first^2) / runs)
stopifnot(ranked$ranking$locations == n, ranked$ranking$total_cost == 44,
          truth_se > 0,
          abs(ranked$ranking$width - first) < 6 * truth_se + 1e-8,
          abs(ranked$by_scenario$mcse - truth_se) < .65 * truth_se)

# Screening affects only the designs that need simulation. A one-pool design
# remains exactly evaluable even when another design has the smallest screen.
alternatives <- list(one = 1, two = c(2, 2), mixed = s, three = c(2, 2, 2))
short <- rank_pool_designs(alternatives, budget = 44, cost_location = 2,
                           cost_pool = 1, cost_unit = 1,
                           scenarios = list(planning = scenario),
                           shortlist = 1, replicates = 8, grid_size = 8)
r <- short$ranking
stopifnot(r$evaluated[r$design == "one"], sum(r$evaluated) == 2,
          all(is.na(r$width[!r$evaluated])), all(is.finite(r$screen_width)),
          !is.unsorted(r$width[r$evaluated]))

# No planning law: all designs are assessed deterministically, including
# those whose expected width would not make a scenario-based shortlist.
seed <- .Random.seed
uniform <- rank_pool_designs(alternatives, budget = 44, cost_location = 2,
                             cost_pool = 1, cost_unit = 1, grid_size = 8,
                             shortlist = 1)$ranking
stopifnot(identical(seed, .Random.seed), nrow(uniform) == 4,
          all(uniform$evaluated), all(uniform$width >= 0 & uniform$width <= 1),
          !is.unsorted(uniform$width))
