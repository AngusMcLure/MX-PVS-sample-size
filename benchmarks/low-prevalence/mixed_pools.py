"""Compare hard-capped 95% prevalence intervals for repeated pools.

Run with Python, NumPy and SciPy: mixed_pools.py initial|extended|refined|smoke.
For faster parallel runs set OPENBLAS_NUM_THREADS=1. Each simulation reports
the largest mean width among its evaluated scenarios, in percentage points,
and that scenario's Monte Carlo standard error. These are comparisons over
specified populations; they do not optimise over every feasible mixing law.
"""

import argparse
import warnings
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache, partial

import numpy as np
from scipy import optimize, stats


PREVALENCE_CAP = 0.02
ICC_CAP = 0.2
INITIAL = [
    (1, 1), (1,) * 6, (2, 2), (3, 3), (4, 4), (2,) * 4, (3,) * 4,
    (2,) * 6, (3,) * 6, (1, 4), (2, 4), (1, 8), (2, 6), (3, 6), (3, 3, 8),
]
EXTENDED = [
    (5, 5), (6, 6), (8, 8), (5, 8), (4, 8), (4, 10), (5, 10), (4, 12),
    (5, 5, 10), (4, 4, 8), (4, 4, 10), (4, 8, 12), (2, 6, 10), (2, 8, 12),
    (2, 4, 8, 16), (3, 6, 12), (3, 6, 15), (4,) * 4, (4, 4, 4, 8),
]
REFINED = [(8, 8), (10, 10), (12, 12), (16, 16), (5, 10), (3, 3, 8),
           (5, 5, 10), (4,) * 4]
SCENARIOS = [(p, shape) for p in (.003, .008, .012, .0199)
             for shape in ('hom', 'hot', 'tail')]
DENSE_SCENARIOS = [
    (p, shape) for p in (.001, .003, .005, .007, .009, .011, .013, .015, .017, .0199)
    for shape in ('hom', 'hot:.025', 'hot:.06', 'hot:.12', 'hot', 'tail:.06', 'tail:.12', 'tail')
]


def population(p, label):
    """Mass points and weights for homogeneous, zero-plus-hot, or unit-tail laws."""
    if label == 'hom':
        return np.array([p]), np.array([1.0])
    kind, _, suffix = label.partition(':')
    rho = float(suffix) if suffix else ICC_CAP
    if kind == 'hot':
        high = p + rho * (1 - p)
        return np.array([0, high]), np.array([1 - p / high, p / high])
    if kind != 'tail':
        raise ValueError(f'Unknown population: {label}')
    low = p * (1 - rho)
    weight = (p - low) / (1 - low)
    return np.array([low, 1.0]), np.array([1 - weight, weight])


def count_model(design):
    sizes = sorted(set(design))
    pools = [design.count(s) for s in sizes]
    counts = np.array(np.meshgrid(*[np.arange(n + 1) for n in pools], indexing='ij'))
    counts = counts.reshape(len(pools), -1)

    def joint(prevalences):
        """Joint conditional probability, outcome rows and site-prevalence columns."""
        pmf = np.ones((counts.shape[1], len(prevalences)))
        log_negative = np.log1p(-prevalences, where=prevalences < 1,
                                out=np.full(len(prevalences), -np.inf))
        for size, positives, n in zip(sizes, counts, pools):
            q = -np.expm1(size * log_negative)
            pmf *= stats.binom.pmf(positives[:, None], n, q[None, :])
        return pmf

    # These pre-specified events have monotone probability in local prevalence.
    totals = counts.sum(axis=0)
    events = np.array([totals >= k for k in np.unique([1, 2, len(design)])])
    if len(sizes) > 1:
        events = np.vstack((events, counts[0] > 0))
    return joint, events


def interval_width_solver(grid, kernels, sites):
    left, right = grid[:-1], grid[1:]
    # If mean <= cap, the ICC constraint implies E[P²-cP] <= 0. Use the
    # minimum of this quadratic over each prevalence cell for an outer bound.
    c = ICC_CAP + (1 - ICC_CAP) * PREVALENCE_CAP
    nearest = np.clip(c / 2, left, right)
    constraints = np.vstack([kernels[:, :-1], -kernels[:, 1:], left,
                             nearest * nearest - c * nearest])
    # Allocate half the total error to the any-positive event. The remaining
    # events share the other half. Exact binomial intervals cover simultaneously.
    n_events = len(kernels)
    allocations = np.r_[.5, np.repeat(.5 / (n_events - 1), n_events - 1)]
    tails = allocations * .05 / 2
    ix = np.arange(sites + 1)
    lower = np.array([np.r_[0, stats.beta.ppf(t, ix[1:], sites - ix[1:] + 1)]
                      for t in tails])
    upper = np.array([np.r_[stats.beta.isf(t, ix[:-1] + 1, sites - ix[:-1]), 1]
                      for t in tails])
    rows = np.arange(n_events)
    normalization = np.ones((1, len(left)))

    @lru_cache(maxsize=None)
    def width(event_counts):
        rhs = np.r_[upper[rows, event_counts], -lower[rows, event_counts],
                    PREVALENCE_CAP, 0.]
        endpoints = []
        for cost in (left, -right):
            with warnings.catch_warnings():
                # HiGHS accepts threads through the SciPy options pass-through.
                warnings.filterwarnings('ignore', 'Unrecognized options detected.*',
                                        optimize.OptimizeWarning)
                solution = optimize.linprog(cost, A_ub=constraints, b_ub=rhs,
                                            A_eq=normalization, b_eq=[1.0],
                                            bounds=(0, None), method='highs',
                                            options={'threads': 1})
            if not solution.success:
                return PREVALENCE_CAP
            # Repair dual feasibility and use the conservative objective.
            dual = np.minimum(solution.ineqlin.marginals, 0)
            endpoint = min(cost - dual @ constraints) + dual @ rhs
            endpoints.append(endpoint - 1e-10 * (1 + sum(abs(dual))))
        return max(0, min(PREVALENCE_CAP, -endpoints[1]) - max(0, endpoints[0]))

    return width


def compare(design, replicates, seed, grid_points=160, refine=False):
    rng = np.random.RandomState(seed + sum((i + 1) * size for i, size in enumerate(design)))
    sites = 50000 // (100 + 10 * len(design) + sum(design))
    grid = np.unique(np.r_[np.linspace(0, 1, grid_points + 1),
                           np.geomspace(.00003, .10, grid_points // 2)])
    joint, events = count_model(design)
    width = interval_width_solver(grid, events.astype(float) @ joint(grid), sites)

    scenarios = SCENARIOS
    if refine:
        # Rank 80 scenarios by width at rounded expected counts, then simulate
        # the first five and three reference populations on fresh random draws.
        screen = []
        for p, label in DENSE_SCENARIOS:
            local, weights = population(p, label)
            counts = np.rint((events @ joint(local) @ weights) * sites).astype(int)
            screen.append((width(tuple(counts)), p, label))
        selected = [(p, label) for _, p, label in sorted(screen, reverse=True)[:5]]
        scenarios = list(dict.fromkeys(selected + [(.012, 'tail'), (.012, 'hot'), (.008, 'hom')]))

    results = []
    for p, label in scenarios:
        local, weights = population(p, label)
        law = joint(local) @ weights
        survey = rng.multinomial(sites, law, size=replicates)
        widths = [width(tuple(counts)) for counts in survey @ events.T]
        results.append((np.mean(widths), np.std(widths, ddof=1) / np.sqrt(replicates), p, label))
    mean, se, p, label = max(results)
    return design, sites, mean * 100, se * 100, p, label


def smoke():
    joint, _ = count_model((2, 2))
    assert np.allclose(joint(np.array([0., 1.])), [[1, 0], [0, 0], [0, 1]])
    for label in ('hom', 'hot:.06', 'tail'):
        local, weights = population(.009, label)
        assert np.isclose(local @ weights, .009)
        assert np.isclose((joint(local) @ weights).sum(), 1)
    result = compare((2, 2), replicates=2, seed=4, grid_points=24)
    assert result[1] == 403 and 0 <= result[2] <= 100 * PREVALENCE_CAP
    print('Smoke checks passed.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('initial', 'extended', 'refined', 'smoke'))
    mode = parser.parse_args().mode
    if mode == 'smoke':
        smoke()
    else:
        runs = {'initial': (INITIAL, 48, 1901, 160), 'extended': (EXTENDED, 48, 3772, 160),
                'refined': (REFINED, 256, 674, 240)}
        designs, replicates, seed, points = runs[mode]
        worker = partial(compare, replicates=replicates, seed=seed, grid_points=points,
                         refine=mode == 'refined')
        print('design sites maximum_mean_width_pp mcse_pp prevalence population', flush=True)
        with ProcessPoolExecutor(max_workers=4) as executor:
            for design, sites, mean, se, p, label in executor.map(worker, designs):
                print(design, sites, f'{mean:.6f}', f'{se:.6f}', p, label, flush=True)
