\name{fit_prevalence}
\alias{fit_prevalence}
\title{Estimate prevalence bounds and cover the entire identified set}
\usage{
fit_prevalence(histogram, pool_size, level = 0.95, grid_size = NULL)
}
\arguments{
\item{histogram}{Numbers of independent locations with zero through \eqn{v} positive pools, in order. Each location contributes the same \eqn{v} nonoverlapping pools of equal size. Include zero-frequency categories.}
\item{pool_size}{Number of individuals in each pool, a positive integer.}
\item{level}{Requested simultaneous coverage for the population identified interval.}
\item{grid_size}{Number of equal subdivisions on each of the prevalence and pool positivity scales. The method merges the two grids. The default grows with the number of locations. For more than one pool per location, a finer grid can improve the fitted law and tighten the outer confidence interval.}
}
\description{
Projects the empirical location-count distribution onto a grid of feasible
binomial mixtures in Kolmogorov distance. Reports sharp bounds for that
fitted distribution. Assumes a perfect assay and conditionally independent
individual outcomes given a location-specific prevalence drawn independently
from any common mixing law on the unit interval.
}
\details{
For multiple pools, a Dvoretzky--Kiefer--Wolfowitz band around the empirical
count CDF covers the true count CDF. The radius is enlarged if needed to include
the fitted distribution. Linear programs bound prevalence across all compatible
mixtures using the minimum and maximum conditional count CDF in each grid
cell. This includes mixing laws with support between nodes. Numerical grid
error can widen the confidence interval. With a single pool, the confidence
limits project an exact binomial interval through the analytic prevalence
bounds. The midpoint minimises worst-case absolute error given a fitted set.
When that set has positive population width, the data cannot yield an
asymptotically unbiased estimate of prevalence across all compatible laws.
If numerical precision prevents computing a sharp fitted interval, the
method instead calculates outer bounds on the fitted interval using the grid
cell program and labels the result in \code{identified_method}.}
\value{A list with \code{identified} (prevalence endpoints at the fitted count law,
with the calculation indicated by \code{identified_method}), \code{midpoint}, \code{confidence} (outer confidence limits
for the entire identified set), \code{fitted_pmf}, \code{distance} (maximum
absolute fitted CDF residual), \code{radius} (final CDF band half-width,
\code{NA} for the single-pool exact binomial calculation), \code{grid_mesh}
(maximum cell width in prevalence), \code{identified_method} (\code{"sharp"} or \code{"grid outer bounds"}), and \code{level}.}
\examples{
fit_prevalence(c(88, 10, 2), pool_size = 5)
}
