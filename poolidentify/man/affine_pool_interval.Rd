\name{affine_pool_interval}
\alias{affine_pool_interval}
\title{Plan and apply a confidence interval with a uniform width guarantee}
\usage{
affine_pool_interval(design, locations, level = 0.95, grid_size = 64L,
                     histogram = NULL)
}
\arguments{
\item{design}{Integer sizes of all pools at a single location, repeated at every location.}
\item{locations}{Number of independent sampled locations.}
\item{level}{Confidence level for the whole prevalence identification interval.}
\item{grid_size}{At least two subdivisions on each of the prevalence and largest-pool positivity scales.}
\item{histogram}{Optional integer joint histogram of negative pool counts across locations. Axes correspond to distinct pool sizes in increasing order. Within each axis the first index represents zero negative pools. Counts must sum to \code{locations}.}
}
\description{Assigns a score to each possible location outcome before sampling. Returns a uniform upper bound on the width of the resulting confidence interval under arbitrary location prevalence distributions.}
\details{
The score is chosen using only the design, number of locations, confidence level and numerical resolution. For each location prevalence, its difference from the expected score is bounded on the full unit interval. Adding these bounds to the sample mean score accounts for structural uncertainty. A two-sided Hoeffding allowance for the independent locations accounts for sampling. This gives simultaneous coverage of all population prevalences compatible with the location count law, under perfect assays and conditionally independent infections. The confidence limits are clipped to the unit interval.

Multiple pooled sizes use a small linear program on grid nodes and a continuous-support correction between nodes. The result can be conservative. For one pool, or a design containing individual tests, the score is computed analytically. The structural gap is a uniform remainder bound for this score. The population-specific sharp identification interval may be narrower. Choose the design and grid before collecting the data.
}
\value{A list containing the uniform \code{width_bound}, its structural contribution \code{structural_gap}, and its full two-sided Hoeffding width contribution \code{sampling_allowance}. The \code{score} array and \code{remainder} limits specify the rule. The field \code{interval} contains the confidence interval for the supplied histogram, or NULL in planning mode. All widths use the probability scale.}
\examples{
affine_pool_interval(c(5, 10), 200)$width_bound
affine_pool_interval(c(5, 10), 200, histogram = matrix(c(6, 12, 2, 180), 2))$interval
}
