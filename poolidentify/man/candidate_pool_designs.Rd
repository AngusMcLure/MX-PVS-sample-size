\name{candidate_pool_designs}
\alias{candidate_pool_designs}
\title{Generate candidate pooling designs for each location}
\usage{
candidate_pool_designs(pool_sizes, max_pools, equal_only = FALSE,
                       max_designs = 10000L)
}
\arguments{
\item{pool_sizes}{Permitted positive integer numbers of individuals per pool. Duplicates are ignored.}
\item{max_pools}{Maximum total number of pools at one location.}
\item{equal_only}{If true, all pools within a candidate have the same size.}
\item{max_designs}{Maximum number of designs to generate. The function checks this limit before enumeration.}
}
\description{
Enumerates allocations of one through \code{max_pools} pools per location.
Each multiset of sizes appears once. With \eqn{b} distinct allowed sizes and
at most \eqn{v} pools, the unrestricted list has \eqn{{b+v \choose b}-1}
entries. Use a small allowed set or the equal-size option for a faster search.
}
\value{A list of integer vectors in order of increasing pool count, then
lexicographic order of pool sizes. Each element describes all the pools at one
location. For instance, \code{c(5L, 5L, 10L)} represents three pools
containing 20 units in total. A list of hand-selected designs has the same format.}
\examples{
candidate_pool_designs(c(5, 10), max_pools = 2)
candidate_pool_designs(c(5, 10), max_pools = 3, equal_only = TRUE)
}
