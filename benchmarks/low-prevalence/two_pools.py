"""Worst expected 95% width with two equal pools per site (NumPy and SciPy).

Example: python two_pools.py 2 3 4 5 --icc .2
Costs, prevalence cap and level match one_pool.py. The data summary is the
total number of positive pools. The bounds apply to arbitrary independent site
prevalences with the specified population mean and individual ICC constraints.
"""
import argparse
import numpy as np
from scipy.optimize import brentq

from one_pool import ALPHA, BUDGET, SITE, TEST, UNIT, U, inverse_pool, min_positivity


def positivity(p, size):
    if p == 1:
        return 1.0
    return -np.expm1(size * np.log1p(-p))


def component_pmfs(q, sites, pool_icc_cap, roots):
    """Row j is the total count PMF with j cap-ICC sites and sites-j independent sites."""
    # At fixed q the two-pool count law is a mixture of Bin(2,q) and
    # 2*Bernoulli(q). Interpolating its ICC therefore gives a mixture of the
    # site laws at zero and at the cap. Given the number j of cap sites, their
    # positions do not matter for the total count.
    independent = (1 - q + q * roots) ** 2
    cap = (1 - pool_icc_cap) * independent + pool_icc_cap * (1 - q + q * roots**2)

    def powers(pgf):
        result = np.empty((sites + 1, pgf.size), dtype=complex)
        result[0] = 1
        result[1:] = pgf
        return np.multiply.accumulate(result, axis=0, out=result)

    spectrum = powers(independent)[::-1] * powers(cap)
    pmfs = np.fft.irfft(spectrum, axis=1)[:, :2 * sites + 1]
    if pmfs.min() < -1e-10 or np.max(np.abs(pmfs.sum(axis=1) - 1)) > 1e-10:
        raise ArithmeticError("Count convolution lost precision.")
    pmfs = np.maximum(0, pmfs)
    return pmfs / pmfs.sum(axis=1, keepdims=True)


def design_intervals(size, individual_icc_cap, q_grid, components, sites):
    """Return one finite-sample confidence interval for each observable count."""
    acceptance = np.empty((len(q_grid), 2), dtype=int)
    tail = ALPHA / 2 - 1e-10  # round numerical decisions outward
    for i, q in enumerate(q_grid):
        pmfs = components(q)
        lower_tail = pmfs.cumsum(axis=1).max(axis=0)
        upper_tail = pmfs[:, ::-1].cumsum(axis=1).max(axis=0)[::-1]
        acceptance[i] = (np.flatnonzero(lower_tail >= tail)[0],
                         np.flatnonzero(upper_tail >= tail)[-1])

    # At fixed component assignment a larger q increases the total count
    # stochastically, by coupling each Bernoulli(q) with the same uniform.
    # For a cell [q_i,q_{i+1}] combine the two outward acceptance cutoffs.
    if np.any(np.diff(acceptance, axis=0) < 0):
        raise ArithmeticError("Acceptance cutoffs should increase with positivity.")
    upper_mean = np.empty_like(q_grid)
    min_q_at_cap = min_positivity(U, size, individual_icc_cap)
    for i, q in enumerate(q_grid):
        if q == 0 or individual_icc_cap == 0:
            upper_mean[i] = inverse_pool(q, size)
        elif q >= min_q_at_cap:
            upper_mean[i] = U
        else:
            upper_mean[i] = brentq(
                lambda p: min_positivity(p, size, individual_icc_cap) - q,
                0, U, xtol=1e-15)
    # If high counts reject every q in the capped parameter space, return the
    # singleton {U}. This retains coverage and explicitly defines its width.
    intervals = np.full((2 * sites + 1, 2), U)
    for count in range(len(intervals)):
        cells = np.flatnonzero((acceptance[:-1, 0] <= count)
                               & (acceptance[1:, 1] >= count))
        if cells.size:
            intervals[count] = (inverse_pool(q_grid[cells[0]], size),
                                upper_mean[cells[-1] + 1])
    return np.clip(intervals, 0, U)


def worst_width(size, individual_icc_cap=.2, grid_points=401, interval_path=None):
    """Return sites, spend, and lower and upper worst expected widths in pp."""
    if (not isinstance(size, (int, np.integer)) or size < 1
            or not 0 <= individual_icc_cap <= 1 or grid_points < 2):
        raise ValueError("Supply a positive pool size, an ICC in [0,1], and at least two nodes.")
    size = int(size)
    cost = SITE + 2 * (TEST + UNIT * size)
    sites = BUDGET // cost
    if sites == 0:
        raise ValueError("The design cannot afford one site.")
    # For G(p)=1-(1-p)^s and b=r+(1-r)U, E[P(P-b)]<=0. Both G(x)/x
    # and the secant slope of G between x and b decrease with x. Therefore
    # G(x)(G(x)-G(b)) <= [G(b)/b]*G'(b)*x*(x-b). Taking expectations
    # gives E[Q^2] <= G(b)E[Q], hence ICC_Q <= G(b). Singletons and r=0
    # retain the given individual ICC bound exactly.
    b = individual_icc_cap + (1 - individual_icc_cap) * U
    pool_icc_cap = (individual_icc_cap if size == 1 or individual_icc_cap == 0
                    else positivity(b, size))
    q_grid = np.linspace(0, positivity(U, size), grid_points)
    fft_length = 1 << (2 * sites).bit_length()
    roots = np.exp(-2j * np.pi * np.arange(fft_length // 2 + 1) / fft_length)

    def components(q):
        return component_pmfs(q, sites, pool_icc_cap, roots)

    intervals = design_intervals(size, individual_icc_cap, q_grid, components, sites)
    if interval_path is not None:
        # Prevalence limits use proportions; worst_width reports percentage points.
        counts = np.arange(intervals.shape[0])
        np.savetxt(interval_path, np.column_stack((counts, intervals)),
                   delimiter=",", fmt=["%d", "%.17g", "%.17g"],
                   header="positive_pools,lower,upper", comments="")
    widths = intervals[:, 1] - intervals[:, 0]
    lower = upper = 0.0
    for q in q_grid:
        # A second pass keeps memory independent of the number of q nodes.
        risks = np.sum(components(q) * widths, axis=1)
        lower = max(lower, risks[0])  # a homogeneous population attains this
        upper = max(upper, risks.max())  # any admissible risk averages rows

    # For any fixed pool ICC the common-uniform coupling changes the total by
    # 2*sites*|delta q| in expectation. Bound the change in risk from the
    # nearest q node by the largest difference between adjacent count widths.
    correction = np.max(np.abs(np.diff(widths))) * sites * (q_grid[1] - q_grid[0])
    return sites, sites * cost, 100 * max(0, lower - 1e-10), 100 * min(U, upper + correction + 1e-10)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sizes", nargs="*", type=int, default=[])
    parser.add_argument("--icc", type=float, default=.2)
    parser.add_argument("--grid-points", type=int, default=401)
    parser.add_argument("--write-intervals", metavar="CSV",
                        help="Save probability-scale confidence limits for a single design.")
    args = parser.parse_args()
    if not args.sizes:
        args.sizes = [3] if args.write_intervals else [1, 2, 3, 4, 5]
    if args.write_intervals and len(args.sizes) != 1:
        parser.error("--write-intervals requires a single design.")
    print("pool_size  sites  spend  worst_width_lower_pp  worst_width_upper_pp")
    for size in args.sizes:
        sites, spend, lower, upper = worst_width(size, args.icc, args.grid_points, args.write_intervals)
        print(f"{size:9d} {sites:6d} {spend:6d} {lower:20.6f} {upper:20.6f}", flush=True)
