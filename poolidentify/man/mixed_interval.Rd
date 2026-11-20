\name{mixed_interval}
\alias{mixed_interval}
\alias{joint_negative_moments}
\title{Prevalence Bounds from Several Pool Sizes}
\description{Compute continuous-support outer bounds for mean prevalence from sparse moments of the uninfected fraction at a location. The population joint count law can be converted into the required moments.}
\usage{
mixed_interval(moments, exponents = as.numeric(names(moments)), cells = 512L)

joint_negative_moments(probabilities, pool_sizes, tolerance = 1e-9)
}
\arguments{
\item{moments}{Numeric vector of feasible population moments of the uninfected fraction \eqn{X = 1-P}. Omit the constant moment.}
\item{exponents}{Distinct positive integers giving the power for each moment. By default read from the names of \code{moments}.}
\item{cells}{Number of intervals in the interpolation grid. Increase to refine the outer bounds.}
\item{probabilities}{Population probabilities as an array. Along each dimension the cells correspond to 0, 1, and so on up to the design maximum number of \emph{negative} pools in that size class. A numeric vector is sufficient for one class.}
\item{pool_sizes}{One positive integer pool size for each axis of the probability array.}
\item{tolerance}{Absolute tolerance on the probability sum and on disagreement between estimates of the same power from different tuples of counts.}
}
\value{\code{mixed_interval} returns a named numeric vector with lower and upper outer bounds for prevalence. \code{joint_negative_moments} returns a named numeric vector of moments, suitable for direct input to \code{mixed_interval}.}
\details{
The units are independent locations, with the same design at every location.
Conditional on prevalence, individuals are independent and testing is perfect.
Within a location, the product of binomial factorial moments of negative counts
identifies the moment whose exponent equals the sum of pool sizes in the selected
subset of pools. The conversion checks agreement when subsets lead to the same
exponent. A histogram from finite sampling can violate these equalities and other
feasibility restrictions.

The bounds are exact with a single supplied moment or when the first moment
is observed. Otherwise a linear program constructs two polynomial majorants, one
of prevalence and one of its negative. The solver uses a grid for interpolation
and checks the curvature allowance over entire cells. The grid does not restrict
the support of the mixing distribution. Floating-point rounding is allowed for
in a final check and constant shift of the majorant. With exact optimisation,
the outer bounds converge to the sharp bounds as the grid is refined.
They do not measure sampling uncertainty. The algorithm checks power-mean
inequalities and rejects inconsistent numerical envelopes. Passing these
checks does not establish that the input moments are feasible.
}
\examples{
# Two classes with one pool each. Independent locations have X uniform on [0,1].
m2 <- 1/3
m3 <- 1/4
m5 <- 1/6
counts <- matrix(c(1-m2-m3+m5, m2-m5, m3-m5, m5), 2, 2)
m <- joint_negative_moments(counts, c(2, 3))
mixed_interval(m)
}
