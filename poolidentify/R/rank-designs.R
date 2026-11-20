# Average widths over independent simulated surveys, reusing results whenever
# the same histogram recurs (in particular for sparse or small surveys).
.planning_mc <- function(law, locations, replicates, width) {
  # Batch the draws to cap memory for the count matrix and the histogram keys.
  # Combine block sums of squared deviations to retain the sample variance.
  batch_limit <- max(1, floor(250000 / length(law)))
  n <- average <- sum_sq <- 0
  while (n < replicates) {
    batch <- min(replicates - n, batch_limit)
    counts <- stats::rmultinom(batch, locations, as.numeric(law))
    keys <- apply(counts, 2L, paste, collapse = ",")
    first <- which(!duplicated(keys))
    unique_widths <- vapply(first, function(i) width(counts[, i]), numeric(1))
    widths <- unique_widths[match(keys, keys[first])]
    block_mean <- mean(widths)
    delta <- block_mean - average
    sum_sq <- sum_sq + sum((widths - block_mean)^2) + delta^2 * (n / (n + batch)) * batch
    average <- average + delta * (batch / (n + batch))
    n <- n + batch
  }
  c(width = average, mcse = sqrt(sum_sq / (replicates - 1) / replicates))
}

.planning_inputs <- function(designs, budget, cost_location, cost_pool, cost_unit) {
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  costs <- list(cost_location, cost_pool, cost_unit)
  if (!scalar(budget) || budget <= 0 || !all(vapply(costs, scalar, logical(1))) ||
      any(unlist(costs) < 0) || all(unlist(costs) == 0))
    stop("Specify a positive budget and three nonnegative costs, at least one positive.")
  if (!is.list(designs) || !length(designs))
    stop("designs must be a nonempty list of vectors of individual pool sizes.")
  designs <- lapply(designs, function(x) sort(.planning_sizes(x)))
  keys <- vapply(designs, paste, collapse = "+", FUN.VALUE = character(1))
  if (anyDuplicated(keys)) stop("The same pool allocation appears more than once.")
  labels <- names(designs)
  if (is.null(labels)) labels <- rep("", length(keys))
  missing <- is.na(labels) | !nzchar(labels)
  labels[missing] <- keys[missing]
  if (anyDuplicated(labels)) stop("Design labels must be unique.")
  names(designs) <- labels
  units <- vapply(designs, function(x) sum(as.double(x)), numeric(1))
  pools <- lengths(designs)
  site_cost <- unname(cost_location + cost_pool * pools + cost_unit * units)
  locations <- floor(budget / site_cost)
  # Correct roundoff at a decimal budget boundary, e.g. 0.30 / 0.10.
  rounded <- round(budget / site_cost)
  close <- is.finite(rounded) & abs(rounded * site_cost - budget) <=
    4 * .Machine$double.eps * budget
  locations[which(close)] <- rounded[which(close)]
  if (any(locations > .Machine$integer.max))
    stop("The number of affordable locations exceeds the supported count range.")
  table <- data.frame(design = labels, pools_per_location = unname(pools),
                      units_per_location = unname(units), cost_per_location = site_cost,
                      locations = as.integer(locations), row.names = NULL)
  feasible <- table$locations >= 1L
  if (!any(feasible)) stop("No candidate design can afford one complete location.")
  excluded <- table[!feasible, c("design", "cost_per_location"), drop = FALSE]
  table <- table[feasible, , drop = FALSE]
  table$total_cost <- table$locations * table$cost_per_location
  list(designs = designs[feasible], table = table, excluded = excluded)
}

rank_pool_designs <- function(designs, budget, cost_location, cost_pool, cost_unit,
                              scenarios = NULL, criterion = c("weighted", "worst"),
                              scenario_weights = NULL, level = 0.95,
                              replicates = 64L, shortlist = 5L, grid_size = 64L) {
  criterion <- match.arg(criterion)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("level must lie strictly between zero and one.")
  if (!.pi_positive_integer(grid_size) || grid_size < 2L ||
      grid_size >= .Machine$integer.max)
    stop("grid_size must be an integer of at least two.")
  inputs <- .planning_inputs(designs, budget, cost_location, cost_pool, cost_unit)
  designs <- inputs$designs
  tab <- inputs$table

  # Score optimisation and confidence envelopes use dense joint-outcome
  # matrices. Check size before constructing an exponentially large count law.
  cells <- vapply(designs, function(d) prod(as.double(table(d)) + 1), numeric(1))
  if (!is.null(scenarios)) .check_planning_grid(cells[lengths(designs) > 1L], grid_size)
  tab$screen_width <- NA_real_
  tab$width <- NA_real_
  tab$mcse <- NA_real_
  tab$evaluated <- TRUE
  details <- NULL
  if (is.null(scenarios)) {
    if (!is.null(scenario_weights)) stop("scenario_weights requires planning scenarios.")
    for (i in seq_along(designs))
      tab$width[i] <- affine_pool_interval(designs[[i]], tab$locations[i], level = level,
                                           grid_size = grid_size)$width_bound
    objective <- "distribution-free width bound"
    tab$mcse <- 0
  } else {
    if (!is.list(scenarios) || !length(scenarios) ||
        any(!vapply(scenarios, function(x) is.list(x) && !is.null(x$prevalence), logical(1))))
      stop("Each scenario must be a list containing prevalence and optional weights.")
    if (!.pi_positive_integer(replicates) || replicates < 2 || replicates > .Machine$integer.max)
      stop("replicates must be an integer of at least two.")
    if (!(identical(shortlist, Inf) || .pi_positive_integer(shortlist)))
      stop("shortlist must be a positive integer or Inf.")
    ns <- length(scenarios)
    sn <- names(scenarios)
    if (is.null(sn)) sn <- paste0("scenario", seq_len(ns))
    if (anyNA(sn) || any(!nzchar(sn)) || anyDuplicated(sn))
      stop("Provide unique scenario names or leave all scenarios unnamed.")
    if (criterion == "worst" && !is.null(scenario_weights))
      stop("scenario_weights apply only to the weighted criterion.")
    if (is.null(scenario_weights)) scenario_weights <- rep(1, ns)
    if (!is.null(names(scenario_weights))) {
      if (!setequal(names(scenario_weights), sn) || anyDuplicated(names(scenario_weights)))
        stop("Named scenario_weights must match the scenario names.")
      scenario_weights <- scenario_weights[sn]
    }
    if (!is.numeric(scenario_weights) || length(scenario_weights) != ns ||
        any(!is.finite(scenario_weights)) || any(scenario_weights < 0) ||
        !any(scenario_weights > 0)) stop("Invalid scenario weights.")
    weights <- scenario_weights / max(scenario_weights)
    weights <- weights / sum(weights)
    aggregate_width <- if (criterion == "worst") max else function(x) sum(weights * x)
    exact <- tab$pools_per_location == 1L
    nd <- length(designs)
    screen <- estimate <- se <- matrix(NA_real_, nd, ns, dimnames = list(tab$design, sn))
    # Store only compact count laws across the screening step. A finalist's
    # confidence geometry is then shared by all scenarios and simulations.
    laws <- vector("list", nd)
    for (i in seq_len(nd)) {
      d <- designs[[i]]
      w <- tab$locations[i]
      laws[[i]] <- lapply(scenarios, function(sc) {
        mass <- sc$weights
        if (is.null(mass)) mass <- rep(1, length(sc$prevalence))
        planning_count_law(d, sc$prevalence, mass)
      })
      if (exact[i]) {
        screen[i, ] <- vapply(laws[[i]], planning_one_expected, numeric(1),
                              locations = w, pool_size = d[1], level = level)
        estimate[i, ] <- screen[i, ]
        se[i, ] <- 0
      } else {
        width <- planning_width_engine(d, w, level, grid_size)
        screen[i, ] <- vapply(laws[[i]], function(law) width(w * law), numeric(1))
      }
    }
    tab$screen_width <- apply(screen, 1, aggregate_width)
    multi <- which(!exact)
    ranked <- multi[order(tab$screen_width[multi])]
    finalists <- ranked[seq_len(min(length(multi), shortlist))]
    tab$evaluated <- exact | seq_len(nd) %in% finalists
    for (i in finalists) {
      width <- planning_width_engine(designs[[i]], tab$locations[i], level, grid_size)
      for (j in seq_len(ns)) {
        result <- .planning_mc(laws[[i]][[j]], tab$locations[i], replicates, width)
        estimate[i, j] <- result["width"]
        se[i, j] <- result["mcse"]
      }
    }
    for (i in which(tab$evaluated)) {
      tab$width[i] <- aggregate_width(estimate[i, ])
      if (criterion == "weighted") tab$mcse[i] <- sqrt(sum((weights * se[i, ])^2))
    }
    details <- data.frame(design = rep(tab$design, ns), scenario = rep(sn, each = nd),
                          screen_width = as.vector(screen), width = as.vector(estimate),
                          mcse = as.vector(se), stringsAsFactors = FALSE)
    objective <- paste(criterion, "scenario expected width")
  }
  tab <- tab[order(!tab$evaluated, tab$width, tab$screen_width, na.last = TRUE), , drop = FALSE]
  rownames(tab) <- NULL
  list(best = designs[[tab$design[1]]], ranking = tab, by_scenario = details,
       excluded = inputs$excluded, objective = objective, level = level, grid_size = grid_size)
}
