\name{identification_interval}
\alias{identification_interval}
\alias{count_moments}
\title{Sharp prevalence bounds for equal pool sizes}
\description{
For a known population law of the number of positive pools per location,
compute the range of mean individual prevalences compatible with that law.
Alternatively supply the raw moments of the location-specific pool
positivity probability \eqn{Q}.
}
\usage{
count_moments(prob)
identification_interval(moments = NULL, pool_size, probabilities = NULL)
}
\arguments{
  \item{prob, probabilities}{Numeric vector of population probabilities for
  counts \eqn{K=0,\ldots,v}, ordered from zero to all pools positive.}
  \item{moments}{Numeric vector \eqn{E(Q),E(Q^2),\ldots,E(Q^v)}.
  Omit the zeroth moment. Supply exactly one of this and \code{probabilities}.}
  \item{pool_size}{Positive integer number of individuals per pool.}
}
\value{
\code{count_moments} returns the vector of \eqn{v} factorial moments of
the counts, each divided by the corresponding factorial power of \eqn{v}.
These are the raw moments of \eqn{Q}. \code{identification_interval}
returns the two named numbers \code{lower} and \code{upper} on a 0--1 scale.
}
\details{
Conditional on \eqn{Q}, the number of positive pools has the binomial
distribution with parameters \eqn{v} and \eqn{Q}. The target is the expectation
of \eqn{1-(1-Q)^{1/s}} for pool size \eqn{s}.
At \eqn{v=1,2} the bounds are available in closed form.
At higher orders, Gaussian, Lobatto or Radau quadrature for the moment
functional gives the endpoint distributions. Their weights are checked by
reconstructing the input count law. The binomial basis avoids alternate
summation when count probabilities are supplied.
Infeasible input or ill-conditioned numerical quadrature may raise an
error. Empirical histograms need not have feasible population moments.
Prefer probabilities over raw moments for higher orders. Quadrature checks each cell on a relative scale. Near unit positivity,
the inverse link may amplify roundoff substantially. The routine stops
if complementary quadrature nodes have poorly resolved prevalence values.
After numerical rank loss, dropping one randomly selected pool gives
a lower-order enclosing interval. The routine returns this enclosure
only when its width is at most one ten-millionth on the prevalence scale,
and stops otherwise.
}
\examples{
# Q has mass 1/2 at 0.1 and at 0.7.
q <- c(0.1, 0.7)
pmf <- (dbinom(0:3, 3, q[1]) + dbinom(0:3, 3, q[2])) / 2
identification_interval(pool_size = 5, probabilities = pmf)
identification_interval(count_moments(pmf), pool_size = 5)
}
