# Choosing a survey design below 2% prevalence

With 50,000 available, spend 100 per independent location, 10 per assay and 1 per individual. Assume the sampled locations are independent draws from a common population, the assays are perfect, and infections at a location are conditionally independent given its prevalence. The prevalence ceiling applies to the population mean. Individual locations may have much higher prevalence. ICC is the correlation of two individuals from the same location.

We compare the **worst expected width of a 95% prevalence interval**. For each possible population, average the full interval width over repeated surveys, then take the supremum over populations satisfying the scientific bounds. The intervals account for both sampling variation and the ambiguity that would remain with infinitely many locations. All widths below are percentage points. We use the closure of each strict bound when computing the supremum.

| Scientific constraints | Design at each location | Locations | Assays | Individuals | Spend | Worst expected width |
|---|---|---:|---:|---:|---:|---:|
| Mean ≤2%, unrestricted ICC | One individual test | 450 | 450 | 450 | 49,950 | 1.5808 |
| Mean ≤2%, ICC ≤0.2 | One pool of four | 438 | 438 | 1,752 | 49,932 | 1.3271 |
| Mean ≤2%, ICC ≤0.2 | Two pools of two | 403 | 806 | 1,612 | 49,972 | 1.3040–1.3096 |
| Mean ≤2%, ICC ≤0.2 | Two pools of three | 396 | 792 | 2,376 | 49,896 | 1.2956–1.3017 |
| Mean ≤2%, ICC ≤0.2 | Two pools of four | 390 | 780 | 3,120 | 49,920 | 1.3222–1.3284 |
| Mean ≤2%, ICC ≤0.2 | Two pools of five | 384 | 768 | 3,840 | 49,920 | 1.3602–1.3669 |
| **Mean ≤2%, ICC ≤0.2** | **Three pools of two** | **367** | **1,101** | **2,202** | **49,912** | **1.2876–1.2920** |

Single numbers are directly maximised expected widths. Ranges bracket the worst expected width for the specified procedure, rounded outwards. They are numerical bounds on a design criterion. The survey itself produces one confidence interval.

With the ICC ceiling, three disjoint pools of two gave the smallest verified worst expected width. This saves about 0.29 percentage points (18%) against the individual-testing benchmark. One pool of four is a practical alternative if running 1,101 assays is difficult. It needs 438 assays, and its worst expected interval is about 0.04 percentage points wider.

## Why the unrestricted answer is minimax

A population can have every individual in a location infected, or every individual uninfected. Even with many assays, that location then yields one Bernoulli observation of the population prevalence. No design can afford more than 450 sampled locations. Conversely, a single individual from each of 450 independent locations gives exactly 450 independent Bernoulli observations for **every** permitted mixing law. This design therefore attains the minimax information benchmark for the unrestricted class. Its reported width uses a central Clopper–Pearson interval with endpoints clipped to the scientific range. Its largest expected width occurs near **0.847% prevalence**, inside the planning range.

## How the ICC cap enters the interval

Write `p` for mean prevalence, `r` for the individual ICC ceiling and `s` for pool size. One pool tests positive with probability `q`. At any given `p`, the sharp range of `q` is

```text
b = r + (1 − r)p
p[1 − (1 − b)^s]/b  ≤  q  ≤  1 − (1 − p)^s.
```

The left endpoint comes from locations of prevalence zero or `b`. The right endpoint comes from constant prevalence. With one pool per location, the total number of positive pools is binomial. We invert an exact binomial interval through both inequalities, average its width over the binomial count and maximise over possible positivity. The derivative of the expected width certifies its maximum for the sizes evaluated. At larger pool sizes a lower bound from the population identification width excludes all remaining affordable sizes. One pool of four is the optimum over **all** one-pool sizes for this interval and the specified costs.

With two equal pools, two assays from the same location can be dependent. The two-pool calculation uses the total number of positive assays and accounts for their dependence. Let `Q = 1 − (1 − P)^s` at a location of prevalence `P`. The site's positive count has the same law as a mixture of a binomial pair with success probability `q` and a pair sharing a single Bernoulli outcome. The mixture weight is the pool-outcome ICC. [Derivation and implementation](../benchmarks/low-prevalence/two_pools.py)

For a population prevalence ceiling `U`, the individual ICC bound implies `E[P²] ≤ B E[P]` with `B = r + (1 − r)U`. Set `G(x) = 1 − (1 − x)^s`. Concavity gives

```text
G(x)[G(x) − G(B)] ≤ [G(B)/B] G′(B) x(x − B).
```

To see this, both `G(x)/x` and the secant slope between `x` and `B` decrease as `x` increases. Their product is at least its value at `B` to the left of `B` and at most that value to the right. Multiplying by `x(x − B)` gives the displayed inequality on both sides. Taking expectations gives `E[Q²] ≤ G(B)q` and bounds the pool-outcome ICC by `G(B)`.

At a fixed `q`, each site's count law is consequently a mixture of the zero-ICC law and the law at this upper bound. For each possible number of sites drawn from the upper-bound component, convolve their count distributions. The largest lower and upper tail probabilities across these component assignments give conservative finite-sample tests. Invert the two 2.5% tails over `q`, then project through the sharp mean bounds above. Stochastic monotonicity in `q` lets the inversion include entire grid cells. A high count that rejects the entire scientific range returns the boundary interval [2%, 2%]. Such a count also gives reason to revisit the scientific assumptions.

To bound the worst expected width, the lower bound uses attainable homogeneous populations. The upper bound takes the largest expected width across component assignments. Between grid points, a common-uniform coupling bounds the change by the largest width change between adjacent counts times the expected change in the total count. This accounts for discretisation between the 401 positivity grid points. The resulting guarantee applies to arbitrary mixing distributions satisfying the bounds. Picking a beta or logit-normal family is unnecessary.

### Three pools of two

Let `q = E[Q]`. When `q > 0`, write `a = E[Q²]/q` and `h = E[Q³]/q`. Conditional binomial sampling gives the following probabilities that at least one, at least two, or all three pools are positive at a location.

```text
q(3 − 3a + h),    q(3a − 2h),    qh.
```

From the preceding inequality, `q ≤ a ≤ G(B)`. Cauchy–Schwarz and `0 ≤ Q ≤ 1` also give `a² ≤ h ≤ a` and `(1 − q)(a − h) ≥ q(1 − a)²`. With `q = 0` all counts are zero. These restrictions give an outer region containing every feasible count law. We bound each of the three probabilities above and below in each small cell of this region. Adding independent locations preserves stochastic ordering, so convolving the two bounding site laws gives lower and upper tail bounds for the survey total. Inverting these bounds and projecting through the sharp one-pool mean inequalities gives an interval for prevalence. Finer subdivision near the tail cutoffs makes the calculation practical. [Three-pool calculation](../benchmarks/low-prevalence/three_pools.py)

The risk calculation uses the same ordering. The expected width equals the width at zero positives plus the sum, over counts, of each successive change in width times the probability of exceeding that count. Use an upper survival bound for increases and a lower one for decreases. Maximising over the cells gives the upper risk bound, while a homogeneous population supplies the lower one. Thus even the upper bound for three pools of two is below the lower bound for two pools of three for the interval procedures evaluated here.

## Comparison with mixed and more heavily replicated designs

We also examined a procedure that uses the numbers of locations with at least one positive, at least two positives and all pools positive. For mixed sizes it also uses the number with a positive among the smallest pools. Simultaneous binomial intervals for these events are projected onto the feasible prevalence range with a linear program. Simulations screened 34 equal- and mixed-size allocations. A refinement checked eight allocations against a wider set of planning populations. For example, one pool of five plus one of ten at 370 sites had mean width **1.440** (Monte Carlo SE 0.008) in a population with 0.9% prevalence and ICC 0.06. Two pools of eight at 367 sites gave **1.443** (SE 0.007) in that population. These are scenario-specific comparisons. Different interval procedures can change the ranking, and this search does not prove global optimality over all mixed or replicated designs. [Simulation and scenarios](../benchmarks/low-prevalence/mixed_pools.py)

For repeated individual tests, a separate check inverted count tests for a subset of feasible populations. Adding more populations can only enlarge this interval. With two through ten individual tests per site, that calculation already gave worst expected widths above 1.32 percentage points for the equal-tailed total-count interval. In particular, six individual tests at each of 301 locations reached at least 1.3356 points. [Repeated individual tests](../benchmarks/low-prevalence/singleton_repeats.py)

## Using the design after sampling

For the recommended ICC-constrained design, count positive pools among the 1,101 assays and look up that row. Zero positives give **0% to 0.282%**. Twenty positives give **0.464% to 1.760%**. The CSV endpoints are proportions. The lookup requires complete data from 367 independent locations with three disjoint pools of two per location and both scientific constraints. [Confidence interval lookup](../benchmarks/low-prevalence/intervals-367x3x2.csv)

```r
limits <- read.csv("benchmarks/low-prevalence/intervals-367x3x2.csv")
100 * subset(limits, positive_pools == 20, select = c(lower, upper))
```

A separate lookup applies to 396 locations with two pools of three at each. With that allocation, 20 positive pools among 792 assays give **0.458% to 1.755%**. [Two-pool confidence interval lookup](../benchmarks/low-prevalence/intervals-396x2x3.csv)

From the repository root, reproduce the comparisons with Python, NumPy and SciPy.

```sh
python benchmarks/low-prevalence/one_pool.py
python benchmarks/low-prevalence/two_pools.py 2 3 4 5
python benchmarks/low-prevalence/three_pools.py
python benchmarks/low-prevalence/check_three_pool_lookup.py
OPENBLAS_NUM_THREADS=1 python benchmarks/low-prevalence/mixed_pools.py refined
```

Generating the three-pool table and its risk bounds took about five and a half minutes on one CPU here. Once generated, analysis requires one row lookup. Direct polynomial convolution independently checked the full three-pool interval table across 101 prevalences, seven ICC values and zero-plus-hot, low-plus-one and beta populations. The smallest probability of covering the entire compatible mean range in these checks was 95.48%. The largest expected width observed was 1.287625 percentage points. [Independent check](../benchmarks/low-prevalence/check_three_pool_lookup.py)

For the two-pool comparison, direct polynomial convolution of the site count law reproduced the risk near the maximum and checked coverage at 101 prevalences and seven ICC values for homogeneous, two-point and beta populations. The smallest observed coverage in this check was 97.44%. The finite-sample coverage guarantee comes from tail inversion.

If either ceiling is only a planning belief, it can select populations for a sensitivity analysis. Using it to exclude values from the reported confidence interval requires it to hold in the surveyed population. For scale, if 2% is only a planning value, the ordinary 95% Clopper–Pearson interval for 450 independent individual tests has expected width 2.808 percentage points at 2% prevalence. Truncating this interval at 2% would require a hard prevalence ceiling.
