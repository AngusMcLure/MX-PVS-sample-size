# MX-PVS-sample-size
Sample size calculations for Post-Validation Surveillance for MX

## Dependencies

Install the R packages required by the document, including PoolPoweR from GitHub:

```r
install.packages(c("shiny", "DT", "plotly", "statmod", "remotes"))
remotes::install_github("AngusMcLure/PoolPoweR")
```

The clustered power calculation uses PoolPoweR's internal `mu_sigma_linknorm()`
helper to convert the requested mean and variance of site prevalence to the
logit-normal parameters. Because this helper is not exported, the document
calls it with `PoolPoweR:::`.
