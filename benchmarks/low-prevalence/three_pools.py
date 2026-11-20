"""Distribution-free 95% intervals for 367 sites, each with three pools of two.

The hard population bounds are mean prevalence <= .02 and individual ICC <= .2.
The only input to the interval is the total number of positive pools. Example
commands, run from this directory:

    python three_pools.py --lookup 12
    python three_pools.py --write-intervals intervals-367x3x2.csv
"""
import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.fft import irfft
from scipy.special import rel_entr
from scipy.stats import binom

from one_pool import ALPHA, BUDGET, SITE, TEST, UNIT, U, inverse_pool

SIZE, POOLS, ICC = 2, 3, .2
COST = SITE + POOLS * (TEST + UNIT * SIZE)
SITES = BUDGET // COST
Q_MAX = -np.expm1(SIZE * np.log1p(-U))
# E[Q^2]/E[Q] <= G(B) for B = ICC + (1-ICC)*U, as shown in two_pools.py.
A_MAX = -np.expm1(SIZE * np.log1p(-(ICC + (1 - ICC) * U)))
Q_CELLS, COARSE_CELLS, SUBDIVISIONS = 625, 12, 12
FFT_SIZE, COUNT_WINDOW = 256, 165
NUMERIC = 2e-10
TAIL = ALPHA / 2 - NUMERIC
TABLE = Path(__file__).with_name("intervals-367x3x2.csv")


@dataclass
class MomentCells:
    """Bounding count distributions and coordinates for a grid of moment cells."""
    low: np.ndarray
    high: np.ndarray
    row: np.ndarray
    col: np.ndarray
    a_upper: np.ndarray
    h_lower: np.ndarray

    def active(self, q_lower):
        # Hausdorff inequalities: q <= a, h <= a - q*(1-a)^2/(1-q).
        return ((self.a_upper >= q_lower - 1e-13) &
                (self.h_lower <= self.a_upper - q_lower *
                 (1 - self.a_upper)**2 / (1 - q_lower) + 1e-13))


def make_cells(n):
    """Enclose the count distribution on each feasible cell in (a,h).

    Write Q=1-(1-P)^2, q=E[Q], a=E[Q²]/q, h=E[Q³]/q. The three
    nonzero count survival probabilities are q*(3-3a+h, 3a-2h, h).
    Under Q-size-biasing, a is a mean and h a second moment, so a²<=h<=a.
    The degenerate case q=0 has no positive pools for any placeholder a,h.
    """
    nodes = np.linspace(0, A_MAX, n + 1)
    row, col = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    al, au = nodes[row], nodes[row + 1]
    hl, hu = nodes[col], nodes[col + 1]
    keep = (hl <= au) & (hu >= al**2)
    row, col, al, au, hl, hu = (v[keep] for v in (row, col, al, au, hl, hu))
    hl, hu = np.maximum(hl, al**2), np.minimum(hu, au)
    al, au = np.maximum(al, hl), np.minimum(au, np.sqrt(hu))
    low = np.column_stack([3 - 3*au + hl, 3*al - 2*hu, hl]).clip(0)
    high = np.column_stack([3 - 3*al + hu, 3*au - 2*hl, hu]).clip(0)

    # The FFT wraps counts >= FFT_SIZE. For each bounding law, Y/3 lies
    # in [0,1] with mean <= m. The iid Chernoff bound uses the Bernoulli
    # mgf envelope. Its tail allowance is absorbed in NUMERIC.
    m = Q_MAX * high.sum(axis=1).max() / POOLS
    x = FFT_SIZE / (POOLS * SITES)
    if not (m < x < 1 and np.exp(-SITES * (rel_entr(x, m) +
                                               rel_entr(1-x, 1-m))) < NUMERIC / 10):
        raise ArithmeticError("Increase the FFT length to control count aliasing.")
    roots = np.exp(-2j * np.pi * np.arange(FFT_SIZE // 2 + 1) / FFT_SIZE)

    def spectrum(survival):
        # Project the survival bounds to decreasing sequences. A valid
        # survival sequence between the two bounds remains between them.
        survival = np.minimum.accumulate(survival, axis=1)
        if Q_MAX * survival[:, 0].max() > 1:
            raise ArithmeticError("Invalid upper bounding count distribution.")
        return (roots - 1) * (survival[:, :1] + survival[:, 1:2] * roots +
                              survival[:, 2:3] * roots**2)

    return MomentCells(spectrum(low), spectrum(high), row, col, au, hl)


def count_cdf(q, spectrum):
    # The per-site PGF is 1+q*spectrum. Independence multiplies site PGFs.
    pmfs = irfft((1 + q * spectrum)**SITES, n=FFT_SIZE, axis=1)
    return pmfs[:, :COUNT_WINDOW].cumsum(axis=1).clip(0, 1)


def tighten(lower, cutoff, witness, coarse_cdf, active, fine_active, parent, fine, q):
    """Refine only the coarse cells that can move an equal-tailed cutoff."""
    while (cutoff < witness if lower else cutoff > witness):
        possible = (coarse_cdf[:, cutoff] >= TAIL if lower else
                    coarse_cdf[:, cutoff - 1] <= 1 - TAIL)
        rows = (possible & active)[parent] & fine_active
        if rows.any():
            tails = (count_cdf(q, fine.low[rows])[:, cutoff] if lower else
                     1 - count_cdf(q, fine.high[rows])[:, cutoff - 1])
            if tails.max() >= TAIL:
                break
        cutoff += 1 if lower else -1
    return cutoff


def acceptance_regions(grid, coarse, fine, parent, keep_risk):
    cuts = []
    cache = []
    for ql, qr in zip(grid[:-1], grid[1:]):
        active, fine_active = coarse.active(ql), fine.active(ql)
        # Stochastic bounds cover every q in [ql,qr] and every nuisance value
        # in the cell. Equal-tail tests may maximize the two tails separately.
        cdf_upper = count_cdf(ql, coarse.low)
        cdf_lower = count_cdf(qr, coarse.high)
        first = np.flatnonzero(cdf_upper[active].max(axis=0) >= TAIL)[0]
        last = np.flatnonzero(np.r_[1, 1-cdf_lower[active].min(axis=0)[:-1]] >= TAIL)[-1]
        if last == COUNT_WINDOW - 1:
            raise ArithmeticError("Increase the count window to include all accepted counts.")
        # The homogeneous location distribution is attainable. Its quantiles
        # stop the refinement once the count has an acceptance witness.
        first = tighten(True, first, binom.ppf(TAIL, POOLS*SITES, ql),
                        cdf_upper, active, fine_active, parent, fine, ql)
        last = tighten(False, last, binom.isf(TAIL, POOLS*SITES, qr),
                       cdf_lower, active, fine_active, parent, fine, qr)
        cuts.append((first, last))
        if keep_risk:
            cache.append((cdf_lower.astype("float32"), cdf_upper.astype("float32")))
    return np.array(cuts), cache


def project_intervals(grid, cuts):
    """Return count-specific confidence intervals on the probability scale."""
    intervals = np.full((POOLS * SITES + 1, 2), U)
    for count in range(COUNT_WINDOW):
        accepted = np.flatnonzero((cuts[:, 0] <= count) & (cuts[:, 1] >= count))
        if not accepted.size:
            continue  # use the singleton {U} if every parameter is rejected
        lower = inverse_pool(grid[accepted[0]], SIZE)
        q_upper = grid[accepted[-1] + 1]
        # For size two, min_positivity(p) = (2-ICC)*p - (1-ICC)*p².
        # Rationalizing its smaller quadratic root is stable near zero.
        upper = 2*q_upper / (2-ICC + np.sqrt((2-ICC)**2 - 4*(1-ICC)*q_upper))
        intervals[count] = (max(0, lower - 1e-14), min(U, upper + 1e-14))
    return intervals


def worst_width(grid, intervals, cache, coarse, fine, parent):
    """Enclose worst expected width over all admissible mixing distributions."""
    widths = np.diff(intervals[:COUNT_WINDOW], axis=1).ravel()
    differences = np.diff(widths)
    pos, neg = differences.clip(0), differences.clip(max=0)
    # Summation by parts: E[W(K)] = W(0) + sum_k dW_k Pr(K>k).
    # Use the appropriate stochastic envelope separately for each sign.
    # 2e-7 covers the float32 CDF cache and the bounded FFT alias/roundoff.
    def bound(cdf_lower, cdf_upper):
        return (widths[0] + np.einsum("ij,j->i", 1-cdf_lower[:, :-1], pos) +
                np.einsum("ij,j->i", 1-cdf_upper[:, :-1], neg) +
                2e-7 * np.abs(differences).sum())

    counts = np.arange(COUNT_WINDOW)
    lower = max(binom.pmf(counts, POOLS*SITES, q) @ widths for q in grid)
    upper = lower
    for i, (ql, qr) in enumerate(zip(grid[:-1], grid[1:])):
        outer = bound(*cache[i])
        could_raise = coarse.active(ql) & (outer > upper)
        rows = could_raise[parent] & fine.active(ql)
        if rows.any():
            upper = max(upper, bound(count_cdf(qr, fine.high[rows]),
                                     count_cdf(ql, fine.low[rows])).max())
    return 100 * max(0, lower - NUMERIC), 100 * min(U, upper + NUMERIC)


def calculate(verify=True):
    coarse = make_cells(COARSE_CELLS)
    fine = make_cells(COARSE_CELLS * SUBDIVISIONS)
    labels = np.full((COARSE_CELLS, COARSE_CELLS), -1)
    labels[coarse.row, coarse.col] = np.arange(coarse.row.size)
    parent = labels[fine.row // SUBDIVISIONS, fine.col // SUBDIVISIONS]
    if parent.min() < 0:
        raise ArithmeticError("A feasible fine cell is missing a parent.")
    grid = np.linspace(0, Q_MAX, Q_CELLS + 1)
    cuts, cache = acceptance_regions(grid, coarse, fine, parent, verify)
    intervals = project_intervals(grid, cuts)
    return intervals, worst_width(grid, intervals, cache, coarse, fine, parent) if verify else None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lookup", type=int, metavar="COUNT",
                        help="Read the 95%% prevalence interval from the distributed CSV.")
    parser.add_argument("--write-intervals", metavar="CSV",
                        help="Regenerate all 1,102 confidence intervals in the specified file.")
    parser.add_argument("--intervals-only", action="store_true",
                        help="Skip the expected-width risk calculation.")
    args = parser.parse_args()
    if args.lookup is not None:
        if not 0 <= args.lookup <= POOLS*SITES:
            parser.error(f"COUNT must be between 0 and {POOLS*SITES}.")
        if args.write_intervals or args.intervals_only:
            parser.error("--lookup reads the existing CSV; use it alone.")
        with TABLE.open(newline="") as handle:
            row = next(row for row in csv.DictReader(handle)
                       if int(row["positive_pools"]) == args.lookup)
        print(f"95% prevalence interval: {100*float(row['lower']):.6f}% to "
              f"{100*float(row['upper']):.6f}%")
    else:
        if args.intervals_only and not args.write_intervals:
            parser.error("--intervals-only requires --write-intervals CSV.")
        intervals, bounds = calculate(verify=not args.intervals_only)
        if args.write_intervals:
            counts = np.arange(intervals.shape[0])
            np.savetxt(args.write_intervals, np.column_stack((counts, intervals)),
                       delimiter=",", fmt=["%d", "%.17g", "%.17g"],
                       header="positive_pools,lower,upper", comments="")
        if bounds is not None:
            print(f"{SITES} locations, cost {COST*SITES}, worst expected width (pp) "
                  f"{bounds[0]:.6f} to {bounds[1]:.6f}")
