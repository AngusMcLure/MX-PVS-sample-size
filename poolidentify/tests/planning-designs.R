library(poolidentify)
generated <- candidate_pool_designs(c(10, 5, 5), 2)
stopifnot(identical(generated, list(5L, 10L, c(5L, 5L), c(5L, 10L), c(10L, 10L))),
          length(candidate_pool_designs(1:4, 4)) == choose(8, 4) - 1,
          length(candidate_pool_designs(1:4, 4, TRUE)) == 16L,
          inherits(try(candidate_pool_designs(1:20, 15), silent = TRUE), "try-error"))

# Mixed-size counts share one latent site prevalence. The observable factorial
# moments must reproduce each corresponding moment of 1-P, including cross terms.
law <- poolidentify:::planning_count_law
p <- c(0, .06, .4, 1)
mass <- c(2, 3, 4, 1)
joint <- law(c(10L, 5L, 5L), p, mass)
moments <- joint_negative_moments(joint, c(5L, 10L))
reference <- vapply(as.numeric(names(moments)), function(r)
  sum(mass * (1 - p)^r) / sum(mass), numeric(1))
stopifnot(identical(dim(joint), c(3L, 2L)),
          max(abs(moments - reference)) < 1e-13)

# With perfectly segregated locations all pools at a site agree.
seg <- law(c(1, 3, 3), c(0, 1), c(.8, .2))
stopifnot(seg[1, 1] == .2, seg[2, 3] == .8, sum(seg > 0) == 2)
# The rare count cell at tiny prevalence survives pooling and mixture evaluation.
tiny <- law(7L, 1e-20)
stopifnot(abs(tiny[1] / 7e-20 - 1) < 1e-13, abs(sum(tiny) - 1) < 1e-15)

# The list can fit its design-count limit and still require long pool vectors.
for (equal in c(FALSE, TRUE)) {
  designs <- candidate_pool_designs(c(2, 3), 4, equal_only = equal)
  entries <- sum(lengths(designs))
  stopifnot(identical(designs, candidate_pool_designs(c(2, 3), 4,
            equal_only = equal, max_entries = entries)),
            inherits(try(candidate_pool_designs(c(2, 3), 4, equal_only = equal,
            max_entries = entries - 1), silent = TRUE), "try-error"))
}
