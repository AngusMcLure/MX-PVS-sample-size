\name{rank_pool_designs}
\alias{rank_pool_designs}
\title{Choose an affordable design by the width of a finite-survey prevalence interval}
\usage{
rank_pool_designs(designs, budget, cost_location, cost_pool, cost_unit,
                 scenarios = NULL, criterion = c("weighted", "worst"),
                 scenario_weights = NULL, level = 0.95,
                 replicates = 64L, shortlist = 5L, grid_size = 64L)
}
\arguments{
\item{designs}{List of candidate allocations. Each vector contains the size of
each pool at one location. For example, \code{c(5, 5, 20)} gives two pools of
five and one of twenty. List names label the results. Unnamed designs are
labelled by their sizes. See \code{\link{candidate_pool_designs}}.}
\item{budget}{Positive total budget, in the same currency as the costs.}
\item{cost_location}{Fixed cost of visiting one independent location.}
\item{cost_pool}{Cost of running one pool test.}
\item{cost_unit}{Cost of collecting one individual unit. Costs are nonnegative
and at least one must be positive.}
\item{scenarios}{Optional list of planning populations. Each is a list with
\code{prevalence} (possible local prevalences) and optional \code{weights}
(nonnegative relative location frequencies, default equal). An empirical
sample of local prevalences can be supplied directly in \code{prevalence}.
Without scenarios, compares distribution-free affine confidence intervals.}
\item{criterion}{With scenarios, minimise the weighted average or the maximum
of the scenario-specific expected widths. The latter ranges over the supplied
scenarios.}
\item{scenario_weights}{Relative frequencies or decision weights for the
weighted criterion, default equal. If named they must match the scenario names.}
\item{level}{Confidence level of the interval that will be reported.}
\item{replicates}{Integer number of independent surveys simulated per
multi-pool finalist per scenario, at least two. Use \code{set.seed} before
calling for reproducibility.}
\item{shortlist}{Evaluate this many multi-pool designs with the smallest
screening scores and every single-pool design. Set \code{Inf} to evaluate
all affordable candidates. With no scenarios all designs are evaluated.}
\item{grid_size}{Number of subdivisions for numerical bounds. Use the same
resolution during survey analysis. A finer grid can change a close ranking.}
}
\description{Compares repeating each proposed pool allocation at the largest
whole number of locations affordable under a fixed survey budget. Considers
uncertainty due to finite sampling and partial identification.}
\details{
The location cost is \eqn{c_L+v c_P+c_U\sum_i s_i} for pool sizes
\eqn{s_1,\ldots,s_v}. The budget fixes the number of independent locations
by integer division. Designs requiring more than the budget for one site
are excluded. Each site has the same fixed allocation, with perfect assays and
conditionally independent infection statuses given its local prevalence.

With planning scenarios, the observed location histogram has a multinomial
distribution. Width is the confidence upper endpoint minus the lower endpoint
for that survey. The analysis uses \code{\link{fit_prevalence}} for equal sizes
and \code{\link{mixed_confidence}} for mixed sizes. For one pool, expected
width is summed over binomial counts, omitting only tails of total probability
at most a few units of machine precision. For multiple pools, confidence width
at the expected (usually fractional) histogram screens candidates. It is a
heuristic and can misorder expected widths. Final widths average the actual
integer-histogram intervals over simulated surveys. This procedure never adds
the structural interval a second time. Screening and Monte Carlo error can
change the winner. To resolve near ties use \code{shortlist = Inf}, increase
\code{replicates} and check the grid resolution.

Without scenarios, the objective is the certified upper bound on the width of
the confidence interval returned by \code{\link{affine_pool_interval}}, for any
mixing law and any sample. Use that procedure on the subsequent survey with
the same design and settings. This deterministic bound and the scenario-based
expected widths refer to different interval procedures and objectives.
}
\value{A list with \code{best} (pool sizes for the first ranked design),
\code{ranking}, \code{by_scenario}, \code{excluded}, \code{objective},
\code{level} and \code{grid_size}. The ranking reports affordable locations,
allocation and cost, \code{width} (probability scale), \code{screen_width},
\code{mcse}, and \code{evaluated}. The lowest final width ranks first.
Skipped designs have missing final widths. For a weighted criterion, its
Monte Carlo standard error combines independent simulations. For the
maximum use the individual errors in \code{by_scenario}. Single-pool
expected widths and distribution-free bounds have zero simulation error.
The screen and scenario details apply only when scenarios are provided.
Multiply widths by 100 to express them in percentage points.}
\examples{
choices <- list(one = 5, two = c(5, 5), mixed = c(2, 5))
pop <- list(rare = list(prevalence = c(0, .1), weights = c(.9, .1)))
set.seed(31)
fit <- rank_pool_designs(choices, budget = 1000, cost_location = 10,
                         cost_pool = 4, cost_unit = .5, scenarios = pop,
                         replicates = 8, shortlist = Inf)
fit$ranking[, c("design", "locations", "total_cost", "width", "mcse")]
}
