# MX-PVS-sample-size
Sample size calculations for Post-Validation Surveillance for MX

## Access the document online

Read the latest version at https://angusmclure.shinyapps.io/mx-pvs-sample-size/

## Dependencies (if running locally)

Install the R packages required by the document, including PoolPoweR from GitHub:

```r
install.packages(c("shiny", "DT", "plotly", "statmod", "remotes"))
remotes::install_github("AngusMcLure/PoolPoweR", ref = "refactor")
```


## Theory papers and R package

[Partial identification of prevalence from clustered group tests](partial-identification.pdf)
studies what a survey can learn about mean individual prevalence when local
prevalence varies arbitrarily across locations. The [manuscript source](partial-identification.Rmd)
includes the proofs and generates its figures and finite-sample tables.
The [companion likelihood note](clustered-group-testing-theory.Rmd) develops
sufficiency, monotone likelihood ratios and fixed-ICC threshold calculations.

From the repository root, install the package and render the paper:

```r
install.packages(c("lpSolve", "rmarkdown", "knitr"))
install.packages("poolidentify", repos = NULL, type = "source")
rmarkdown::render("partial-identification.Rmd")  # PDF, requires XeLaTeX and Pandoc
# For an HTML copy:
rmarkdown::render("partial-identification.Rmd", "html_document")
```

Use a probability vector when the population count law is known, and
an integer histogram when counts come from a finite sample. Both vectors
start with locations having zero positive pools.

```r
library(poolidentify)
law <- mixing_counts(pool_size = 25, pools = 2, family = "beta", par = c(.01, .10))
identification_interval(probabilities = law$probabilities, pool_size = 25)
fit_prevalence(c(88, 10, 2), pool_size = 5)$confidence
set.seed(4216)
bayes_prevalence(c(24, 4, 2), pool_size = 25, prior = "shrink")$envelope
```

For a mixed-size design, `joint_negative_moments()` takes a joint array of
**negative** pool counts indexed by size class. `mixed_interval()` returns
numerical outer bounds over unrestricted mixing distributions. For a finite joint
histogram, use `mixed_confidence()` for simultaneous confidence limits or
`mixed_bayes_prevalence()` for posterior inference.

Working-model bias and coverage calculations use equal-sized pools. The design
planner accepts both equal and mixed sizes. To choose the number and sizes of pools, give each candidate as the allocation
at **one** location. The budget determines how many independent locations the
survey can afford. This example compares the expected width of a 95% interval
under two possible populations.

```r
designs <- candidate_pool_designs(c(2, 5, 10), max_pools = 2)
populations <- list(
  diffuse = list(prevalence = c(0, .06), weights = c(.5, .5)),
  concentrated = list(prevalence = c(0, .15), weights = c(.8, .2))
)
set.seed(42)
plan <- rank_pool_designs(designs, budget = 50000, cost_location = 20,
                          cost_pool = 5, cost_unit = .5, scenarios = populations)
head(plan$ranking[, c("design", "locations", "total_cost", "width", "mcse")])
```

Widths are on the prevalence scale. By default the planner screens multi-pool
designs and simulates the best five. Use `shortlist = Inf` to simulate every
multi-pool candidate. Single-pool expected widths are calculated directly.
With no planning population, the same function compares conservative width
guarantees for the package's affine interval procedure.

For a survey with a 50,000 budget and prevalence below 2%, a worked comparison
shows how a scientific upper bound on individual ICC changes the allocation
and the attainable 95% interval width. [Design comparison and interval lookup](documentation/low-prevalence-design.md)

When the survey is completed, use the confidence level and grid from the
ranking. For equal sizes, use `fit_prevalence(h, s, level = plan$level, grid_size = plan$grid_size)` with a positive-count histogram. For mixed
sizes, pass the joint negative-count array and its size classes to `mixed_confidence(h, sizes, level = plan$level, grid_size = plan$grid_size)`.
If you planned without scenarios, use `affine_pool_interval()` with the chosen allocation, number
of locations, level, grid, and observed negative-count histogram.

See the [package](poolidentify) for function documentation and the
[verification record](documentation/partial-identification-audit.md) for the scope and
independent checks. To run the package checks from the repository root:

```sh
R CMD build poolidentify
R CMD check poolidentify_0.1.0.tar.gz
```
