\name{mixed_bayes_prevalence}
\alias{mixed_bayes_prevalence}
\title{Bayesian inference with several pool sizes per location}
\usage{
mixed_bayes_prevalence(histogram, pool_sizes, draws = 1000L,
  burn = 500L, chains = 4L, concentration = 1, atom_rate = 1, level = 0.95)
}
\arguments{
\item{histogram}{Integer location count array. Axis j indexes the number of
negative pools of size \code{pool_sizes[j]}, starting at zero. Each location has
the same pool allocation. A vector is accepted for a single axis.}
\item{pool_sizes}{Number of individuals per pool on each axis.}
\item{draws, burn, chains}{Retained sweeps per chain (at least eight),
discarded sweeps and independent runs.}
\item{concentration}{Positive total concentration of the symmetric Dirichlet
weight prior. Small values favour fewer effectively occupied components.}
\item{atom_rate}{Positive exponential rate for minus the logarithm of an atom's
negative probability in the smallest pool. Equivalently that negative
probability has a Beta(atom_rate,1) prior.}
\item{level}{Credibility level between zero and one.}
}
\description{Use atomic representations of the joint moments and prevalence.
The number of atoms includes every feasible observable law and mean in the
support of the induced prior.}
\details{If there are d distinct positive exponents in the mixed design augmented
by the exponent one, the lower principal representation needs at most
1+floor(d/2) atoms. Locations receive common exchangeable weights from a
symmetric Dirichlet prior. Changing the weight concentration or the atom prior
changes beliefs even conditional on known count probabilities. Label allocations
are sampled as multinomials per occupied histogram cell. Conditional updates
for atoms use conjugate exponential sampling when their pools are all negative
and a univariate slice update otherwise. Check mixing across chains and compare
hyperparameter settings. For a scientifically specified conditional prior
within a sharp interval for equal sizes see \code{bayes_prevalence}.
}
\value{A list with posterior \code{draws} (chain, prevalence, effective number
of atoms defined as the reciprocal of the sum of squared weights, and moments of
individual negativity labelled by exponent), posterior \code{mean}, equal-tail
\code{interval}, number of \code{atoms}, \code{diagnostics} (split Rhat and
batch-means Monte Carlo standard errors of means), and prior and level settings.
The interval is conditional on the specified prior.}
\examples{
set.seed(627)
h <- matrix(c(2, 3, 5, 20), 2, 2)
fit <- mixed_bayes_prevalence(h, c(4, 10), draws = 40, burn = 20, chains = 2)
fit$interval
}
