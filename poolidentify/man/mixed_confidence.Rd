\name{mixed_confidence}
\alias{mixed_confidence}
\title{Finite-location prevalence interval for several pool sizes}
\usage{
mixed_confidence(histogram, pool_sizes, level = 0.95, grid_size = NULL)
}
\arguments{
\item{histogram}{Joint histogram of independent locations. Along axis \eqn{a}, indices run from zero through \eqn{v_a} \emph{negative} pools of size \eqn{s_a} (the first index means zero negatives). Each location supplies the same pools. Include cells with zero observed frequency.}
\item{pool_sizes}{Positive integer pool sizes, one per histogram axis in the same order.}
\item{level}{Requested simultaneous coverage of the population identified set.}
\item{grid_size}{Number of subdivisions on each of the prevalence and largest-pool positivity scales. The default grows with sample size.}
}
\description{
Finds an outer confidence interval for prevalence under arbitrary mixing on the unit interval. Assumes independent locations, perfect assays, and conditionally independent individuals within locations.}
\details{
Each joint histogram cell yields an exact marginal binomial interval. Bonferroni adjustment covers all cells simultaneously. The intervals expand if needed to include a feasible joint probability array fitted by minimising the maximum cell-probability residual. Conditional joint probabilities are products of binomial kernels. Bounds on each factor over every prevalence grid cell include its mode when interior. Linear programming on the resulting products covers mixing laws supported anywhere within each grid cell. A coarse grid or a large number of joint outcomes can give conservative limits.}
\value{A list with \code{confidence} (lower and upper prevalence limits),
\code{fitted_pmf} (a feasible joint distribution), \code{lower_probability} and
\code{upper_probability} (simultaneous cell limits including the fit),
\code{distance} (maximum cell-probability residual), \code{grid_mesh}
(maximum prevalence cell width), and \code{level}. The fitted joint distribution
can be passed to \code{joint_negative_moments} and \code{mixed_interval} to
estimate population identification bounds.}
\examples{
# Each location has one size-two pool and one size-five pool.
fit <- mixed_confidence(matrix(c(8, 12, 3, 77), 2), c(2, 5))
fit$confidence
}
