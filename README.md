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
See the [package](poolidentify) for function documentation and verification
examples. To run the package checks from the repository root:

```sh
R CMD build poolidentify
R CMD check poolidentify_0.1.0.tar.gz
```
