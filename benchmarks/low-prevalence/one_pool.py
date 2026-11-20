"""Worst expected 95% width for one pool per site. Requires NumPy and SciPy.

The two scientific bounds apply to the *population* mean and individual ICC.
Both endpoints of the projected Clopper-Pearson interval are clipped to [0,U].
"""
import numpy as np
from scipy.optimize import brentq
from scipy.stats import beta, binom

BUDGET, SITE, TEST, UNIT = 50_000, 100, 10, 1
U, ALPHA = 0.02, 0.05


def inverse_pool(q, size):
    return -np.expm1(np.log1p(-q) / size)


def min_positivity(p, size, r):
    b = r + (1 - r) * p
    if b == 1:
        return p
    return (p / b) * -np.expm1(size * np.log1p(-b))


def max_risk(size, r):
    sites = BUDGET // (SITE + TEST + UNIT * size)
    counts = np.arange(sites + 1)
    q_lower = np.r_[0, beta.ppf(ALPHA / 2, counts[1:], sites - counts[1:] + 1)]
    q_upper = np.r_[beta.isf(ALPHA / 2, counts[:-1] + 1, sites - counts[:-1]), 1]
    max_q = -np.expm1(size * np.log1p(-U))
    widths = np.array([
        U if q >= min_positivity(U, size, r) else
        brentq(lambda p: min_positivity(p, size, r) - q, 0, U, xtol=1e-15)
        for q in q_upper
    ]) - np.minimum(U, inverse_pool(q_lower, size))

    # The derivative of the expected width is a Bernstein polynomial with
    # coefficients proportional to adjacent differences of the interval widths.
    # A single sign change from + to - certifies at most one interior maximum.
    differences = np.diff(widths)
    keep = differences != 0
    signs = np.sign(differences[keep])
    if np.any(np.diff(signs) > 0):
        raise ValueError("The width sequence requires a general risk maximiser.")

    def derivative_sign(q):
        log_weights = binom.logpmf(counts[:-1][keep], sites - 1, q)
        if not np.any(np.isfinite(log_weights)):
            return np.sign(differences[keep][0 if q == 0 else -1])
        return np.exp(log_weights - max(log_weights)) @ differences[keep]

    if not np.any(keep) or derivative_sign(0) <= 0:
        q_star = 0
    elif derivative_sign(max_q) >= 0:
        q_star = max_q
    else:
        q_star = brentq(derivative_sign, 0, max_q, xtol=1e-14)
    return sites, q_star, binom.pmf(counts, sites, q_star) @ widths


def compare(r):
    candidates = []
    for size in range(1, (BUDGET - SITE - TEST) // UNIT + 1):
        # At q=min_positivity(U,size,r), two compatible means are U and
        # inverse_pool(q,size). Simultaneous coverage makes 0.95 times that
        # separation a lower bound on worst expected width. q <= U/b, so
        # this weaker bound increases with size and can exclude *all* larger pools.
        b = r + (1 - r) * U
        lower = (1 - ALPHA) * (U - inverse_pool(U / b, size))
        if candidates and lower > min(row[-1] for row in candidates):
            break
        candidates.append((size, *max_risk(size, r)))
    return sorted(candidates, key=lambda row: row[-1])


if __name__ == "__main__":
    for cap in (1., .2):
        print(f"ICC <= {cap}: pool size, sites, worst pool positivity, width (pp)")
        for size, sites, q, width in compare(cap)[:5]:
            print(f"{size:2d} {sites:4d} {q:.6f} {100*width:.6f}")
