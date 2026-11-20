"""Check the 367-location lookup by direct polynomial convolution (NumPy only)."""

from pathlib import Path

import numpy as np


def negative_moments(p, rho, shape):
    """E[(1-P)^s] for s=2,4,6 under a named, feasible site population."""
    powers = np.array([2, 4, 6])
    if p == 0 or rho == 0:
        return (1 - p) ** powers
    if shape == "zero-hot":
        hot = p + rho * (1 - p)
        return 1 - (p / hot) * (1 - (1 - hot) ** powers)
    if shape == "low-one":
        low = p * (1 - rho)
        return (1 - p) * (1 - low) ** (powers - 1)
    concentration = 1 / rho - 1
    i = np.arange(6)
    moments = np.cumprod(((1 - p) * concentration + i) / (concentration + i))
    return moments[powers - 1]


def survey_probabilities(moments):
    """Multiply independent site count polynomials, using ordinary convolution."""
    m2, m4, m6 = moments
    site = np.array([m6, 3 * (m4 - m6), 3 * (m2 - 2 * m4 + m6),
                     1 - 3 * m2 + 3 * m4 - m6])
    assert np.min(site) > -1e-12
    factor = np.maximum(site, 0)
    total = np.ones(1)
    sites = 367
    while sites:
        if sites % 2:
            total = np.convolve(total, factor)
        sites //= 2
        if sites:
            factor = np.convolve(factor, factor)
    assert len(total) == 1102 and abs(total.sum() - 1) < 1e-11
    return total


def check(path):
    table = np.genfromtxt(path, delimiter=",", names=True)
    counts, lower, upper = (table[name] for name in ("positive_pools", "lower", "upper"))
    assert np.array_equal(counts, np.arange(1102))
    assert np.all((0 <= lower) & (lower <= upper) & (upper <= .02))
    min_coverage, max_width = 1., 0.
    for p in np.linspace(0, .02, 101):
        for rho in [0, .025, .06, .10, .15, .195, .2]:
            for shape in ["zero-hot", "low-one", "beta"]:
                moments = negative_moments(p, rho, shape)
                probabilities = survey_probabilities(moments)
                q = 1 - moments[0]
                # With s=2, q is bounded between 1.8p-0.8p² and 2p-p².
                # Check simultaneous coverage of both compatible mean extremes.
                left = q / (1 + np.sqrt(1 - q))
                right = min(.02, 2 * q / (1.8 + np.sqrt(1.8**2 - 3.2 * q)))
                covered = (lower <= left + 1e-14) & (upper >= right - 1e-14)
                min_coverage = min(min_coverage, probabilities[covered].sum())
                max_width = max(max_width, probabilities @ (upper - lower))
    assert min_coverage >= .95 and max_width * 100 <= 1.2921
    print(f"Lowest full-set coverage: {min_coverage:.6%}")
    print(f"Largest expected width in checked populations: {100 * max_width:.6f} pp")


if __name__ == "__main__":
    check(Path(__file__).with_name("intervals-367x3x2.csv"))
