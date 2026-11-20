\name{bayes_prevalence}
\alias{bayes_prevalence}
\title{Bayesian inference on feasible moments of pool positivity}
\usage{
bayes_prevalence(counts, pool_size, draws = 1000L, burn = 500L,
  chains = 4L, prior = c("uniform", "shrink", "volume"),
  shape = NULL, within = NULL, level = 0.95)
}
\arguments{
\item{counts}{Nonnegative integer location histogram in increasing order of
the number of positive pools, from zero to the common number of pools per location.}
\item{pool_size}{Number of individuals in every pool.}
\item{draws, burn, chains}{Number of retained iterations per chain (at least eight),
discarded sweeps and independent runs. The one-pool case uses direct beta sampling.}
\item{prior}{Independent canonical-coordinate priors: uniform coordinates,
shrinkage on even coordinates, or uniform Lebesgue volume in raw moments.}
\item{shape}{Optional matrix with one row per canonical coordinate and two positive beta shapes per row. Overrides \code{prior}.}
\item{within}{Optional pair of positive beta shapes for the relative position of prevalence
within each draw's sharp interval. NULL reports only the identified set.}
\item{level}{Posterior probability level strictly between zero and one.}
}
\description{Fit the multinomial likelihood for independent locations without
specifying a mixing family for location prevalence. Each canonical moment is
the fraction of the allowable range of the corresponding raw moment, conditional
on earlier moments.}
\details{For v pools, the volume prior uses Beta(v-j+1,v-j+1) on
coordinate j. The shrinkage prior uses Beta(1,1+j) on even coordinates j and
uniform distributions elsewhere. The even coordinates influence the
off-diagonal coefficients of the Jacobi matrix for quadrature. Shrinking them
toward zero favours lower effective quadrature rank. Every prior option has
full support on the closed feasible moment space, with no mass exactly on its
boundary. The sampler uses a coordinate slice update in prior probability
coordinates, with principal quadrature for the histogram likelihood.
Computations are intended for modest numbers of pools per location. Examine
the split-chain Rhat diagnostics and Monte Carlo standard errors and compare
longer independent runs, especially if cells are rare or priors concentrated.}
\value{A list. \code{draws} contains chain, canonical coordinates c1 through cv,
sharp lower and upper endpoints, and optional prevalence draws.
\code{envelope} uses opposite tail quantiles of the two endpoints and covers
the whole interval with at least the nominal posterior probability up to Monte
Carlo error. \code{interval} is the equal-tail prevalence credible interval
under \code{within}, or NULL. \code{diagnostics} gives ordinary split Rhat
(NA for a constant series) and batch-means Monte Carlo standard errors of
posterior means. \code{level}, \code{prior}, \code{shape}, and \code{within}
record the settings.}
\examples{
set.seed(431)
fit <- bayes_prevalence(c(21, 6, 3), 8, draws = 64, burn = 32, chains = 2)
fit$envelope
}
