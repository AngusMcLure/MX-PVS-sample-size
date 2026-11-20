\name{naive_bias_bounds}
\alias{naive_bias_bounds}
\title{Sharp bounds on asymptotic signed bias}
\usage{naive_bias_bounds(prevalence, icc, pool_size)}
\arguments{
\item{prevalence}{Mean individual prevalence in the location population.}
\item{icc}{Intraclass correlation for individual infection indicators.}
\item{pool_size}{Number of individuals per pool.}
}
\description{Extremises the probability limit of the homogeneous independent-pool fit over all mixing distributions with the specified mean and individual intraclass correlation.}
\value{Named lower and upper bounds for the probability limit minus the true prevalence, on the probability scale.}
\examples{
naive_bias_bounds(0.01, 0.1, 25)
}
