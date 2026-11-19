\name{working_models}
\alias{mixing_counts}
\alias{working_projection}
\alias{working_performance}
\title{Working mixing families, prevalence bias and interval coverage}
\description{Compute the count distribution for specified mixing parameters, fit working models to a count law, and evaluate bias and coverage under independent sampling of locations.}
\usage{
mixing_counts(pool_size, pools, family = c("beta", "logitnormal"),
              par, quadrature = 80L)
working_projection(probabilities, pool_size,
                   family = c("naive", "beta", "logitnormal"),
                   quadrature = 80L)
working_performance(probabilities, prevalence, pool_size,
                    family = c("naive", "beta", "logitnormal"),
                    locations, level = 0.95, replicates = 200L,
                    quadrature = 80L)
}
\arguments{
\item{pool_size}{Individuals per pool, a positive integer.}
\item{pools}{Pools per independent location, a positive integer.}
\item{family}{Working model for analysis or family used to generate counts.}
\item{par}{Beta mean and individual intraclass correlation in [0,1], or
logit-normal mean and standard deviation on the logit scale. Beta correlation
zero gives a common prevalence and one gives an all-or-nothing mixture.}
\item{quadrature}{Logit-normal integration resolution, integer of at least 16.
Each subinterval uses at least four Gauss--Legendre nodes. The default uses ten.
Increase for an accuracy check.}
\item{probabilities}{Population probabilities (or empirical proportions) for
zero through v positive pools at one location.}
\item{prevalence}{True mean prevalence under the data-generating law, in [0,1].
The supplied count law must be feasible and its sharp identification interval must contain this value.}
\item{locations}{Number of independent locations in each finite sample.}
\item{level}{Nominal confidence level strictly between zero and one.}
\item{replicates}{Number of Monte Carlo samples for beta and logit-normal coverage.
Set the random seed before calling for reproducibility.}
}
\details{
Beta calculations condition on the beta-binomial number of infected individuals.
An occupancy recursion distributes them uniformly among pool positions.
Logit-normal calculations integrate binomial probabilities against a standard
normal density, resolving the transition of pool positivity on the logit scale.

Independent-pool inference inverts average positivity. The family fit minimizes
cross-entropy of the location-count probabilities by multistart bounded numerical
optimization and checks the homogeneous limit. Both family parameters require
at least two pools per location. For beta fits each logit parameter is bounded
(the mean to [-26,26] and correlation to [-18,16]); logit-normal location lies
in [-120,120] and log standard deviation in [-8,log(60)]. A reported search
limit or optimization failure is not evidence of an interior optimum. Small
numerical KL divergence indicates agreement with the observable count law.

Naive finite-sample bias and coverage use the convolution of the count
distribution over all locations. Coverage is for the independent-pool Wilson
interval transformed to prevalence. Family results simulate multinomial
location histograms and refit the two-parameter likelihood for each sample.
Intervals use expected model information and the delta method, clipped to
[0,1]. A boundary fit or ill-conditioned information gives no regular Wald
interval. Coverage treats an unreported interval as an empty interval, and
failure fractions are reported explicitly. A nonzero point failure fraction
makes the corresponding unconditional finite-sample point bias unavailable.
These model-based intervals can have low coverage under heterogeneity or a
wrong mixing family. Sampling error of Monte Carlo coverage is reported;
zero empirical error when all outcomes agree does not establish zero true error.
}
\value{
\code{mixing_counts} returns \code{probabilities} and \code{prevalence}.
\code{working_projection} returns the fitted \code{prevalence}, cell
\code{probabilities}, natural \code{par}, fit \code{status} and
count-law Kullback--Leibler divergence \code{kl}. Regular two-parameter
fits also return optimization coordinates in \code{internal}.
\code{working_performance} returns the population \code{target},
\code{asymptotic_bias}, finite-sample \code{bias}, unconditional
\code{coverage}, its \code{coverage_mcse}, \code{failure_fraction} and
interval \code{method}. Family fits also include \code{projection_status}
and \code{point_failure_fraction}.
}
\examples{
truth <- mixing_counts(10, 2, "beta", c(0.05, 0.12))
working_projection(truth$probabilities, 10, "logitnormal")$prevalence
working_performance(truth$probabilities, truth$prevalence, 10, "naive",
                    locations = 50)$coverage
}
