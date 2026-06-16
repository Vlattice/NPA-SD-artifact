"""The vectorized per-variant NPA test must match the set-based definition.

The scan computes, per biallelic SNV, whether the child carries an allele absent
from both parents using boolean array operations. This checks that result against
the literal set difference over many random genotypes, including missing calls.
"""

import numpy as np


def _vectorized(child, father, mother):
    par0 = (father == 0).any(1) | (mother == 0).any(1)
    par1 = (father == 1).any(1) | (mother == 1).any(1)
    return (((child == 0).any(1) & ~par0)
            | ((child == 1).any(1) & ~par1))


def _set_based(child, father, mother):
    out = np.empty(len(child), dtype=bool)
    for i in range(len(child)):
        child_alleles = {a for a in child[i] if a >= 0}
        parental = {a for a in father[i] if a >= 0} | {a for a in mother[i] if a >= 0}
        out[i] = bool(child_alleles - parental)
    return out


def test_npa_matches_set_definition():
    rng = np.random.default_rng(0)
    n = 50_000
    g = rng.choice([-1, 0, 1], size=(n, 3, 2), p=[0.05, 0.6, 0.35])
    child, father, mother = g[:, 0], g[:, 1], g[:, 2]
    assert np.array_equal(
        _vectorized(child, father, mother),
        _set_based(child, father, mother),
    )
