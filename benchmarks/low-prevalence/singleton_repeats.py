"""Lower bounds on worst expected width for repeated individual tests at a site.

The target interval inverts two 2.5% tests based on the total number positive.
At each hypothesised prevalence each tail is maximised over all mixing laws
with individual ICC at most 0.2. This script uses a subset of those laws
and prevalences, so the interval hull it constructs lies inside the target.
The comparison applies specifically to this equal-tailed count procedure.
Requires NumPy and SciPy.
"""

import argparse

import numpy as np
from scipy.stats import binom


ICC = np.array([0, .025, .06, .12, .17, .2])
PREVALENCE = np.linspace(0, .02, 401)


def site_laws(prevalence, replicates):
    """Conditional-binomial mixtures with the specified population mean and ICC."""
    counts = np.arange(replicates + 1)
    # Zero or a higher local prevalence, with mass chosen to give the mean.
    hot = prevalence + ICC * (1 - prevalence)
    hot_mass = np.divide(prevalence, hot, out=np.zeros_like(hot), where=hot > 0)
    zero_hot = hot_mass[:, None] * binom.pmf(counts, replicates, hot[:, None])
    zero_hot[:, 0] += 1 - hot_mass

    # A lower local prevalence or one, again matching the same two moments.
    low = prevalence * (1 - ICC)
    unit_mass = (prevalence - low) / (1 - low)
    low_one = (1 - unit_mass[:, None]) * binom.pmf(counts, replicates, low[:, None])
    low_one[:, -1] += unit_mass
    return np.vstack([zero_hot, low_one])


def lower_bound(replicates, prevalences=PREVALENCE):
    sites = 50_000 // (100 + 11 * replicates)
    total = sites * replicates
    # Omitting the nonnegative width contributions of larger counts is safe.
    counts = np.arange(min(125, total + 1))
    fft_size = 1 << total.bit_length()
    accepted = []
    for prevalence in prevalences:
        # Independent-site convolution. The transform contains the full
        # polynomial support so no high-count mass wraps to the left tail.
        site = site_laws(prevalence, replicates)
        transform = np.fft.rfft(site, n=fft_size, axis=1)
        survey = np.fft.irfft(transform ** sites, n=fft_size, axis=1)
        cdf = np.cumsum(survey[:, :len(counts)], axis=1)
        lower_tail = cdf.max(axis=0)
        upper_tail = 1 - np.r_[0., cdf.min(axis=0)[:-1]]
        # The tails can attain their maxima under different mixing laws.
        # The small inward tolerance avoids accepting a numerical near-tie.
        accepted.append((lower_tail > .025 + 1e-10) & (upper_tail > .025 + 1e-10))

    widths = np.zeros(len(counts))
    for k, keep in enumerate(np.asarray(accepted).T):
        if np.any(keep):
            widths[k] = np.ptp(prevalences[keep])
    # The homogeneous population is feasible for every prevalence. Its
    # exact binomial expectation supplies a lower bound on worst-case risk.
    risk = binom.pmf(counts, total, prevalences[:, None]) @ widths
    worst = risk.argmax()
    return sites, prevalences[worst], risk[worst]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("replicates", nargs="*", type=int, default=list(range(2, 11)))
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    if args.smoke:
        laws = site_laws(.013, 6)
        counts = np.arange(7)
        expected_second = .013**2 + ICC * .013 * (1 - .013)
        assert np.allclose(laws.sum(axis=1), 1)
        assert np.allclose(laws @ counts, 6 * .013)
        assert np.allclose(laws @ (counts * (counts - 1)) / 30, np.tile(expected_second, 2))
        sites, _, width = lower_bound(6, np.linspace(0, .02, 21))
        assert sites == 301 and 0 < width < .02
        print("Smoke check passed.")
    else:
        if any(v < 1 or 100 + 11 * v > 50_000 for v in args.replicates):
            parser.error("Choose a positive, affordable number of individual tests per location.")
        print("tests/site  sites  homogeneous prevalence  width lower bound (pp)")
        for v in args.replicates:
            sites, prevalence, width = lower_bound(v)
            print(f"{v:10d}  {sites:5d}  {prevalence:22.5f}  {100 * width:22.5f}")
