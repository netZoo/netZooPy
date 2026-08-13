"""DRAGON cross-language parity test.

Runs the netZooPy DRAGON pipeline on shared inputs and asserts the outputs
match the snapshot in tests/dragon_parity/. The same inputs and gold values
live byte-identical in netZooR, where a mirror test asserts the same numbers.
If both tests pass, the two implementations agree to the documented tolerances.

Coverage: lambdas, shrunken covariance, precision matrix, partial correlation
matrix. Kappa / p-values are not covered because netZooR's
estimate_kappa_dragon and estimate_p_values_dragon are unimplemented stubs.
"""
import os
import numpy as np
import pytest
from netZooPy import dragon

FIXTURE_DIR = os.path.join(os.path.dirname(__file__), "dragon_parity")


def _load(name):
    return np.loadtxt(os.path.join(FIXTURE_DIR, name), delimiter=",")


def test_dragon_parity_lambdas():
    X1 = _load("X1.csv")
    X2 = _load("X2.csv")
    expected = np.loadtxt(os.path.join(FIXTURE_DIR, "lambdas.txt"))
    lambdas, _ = dragon.estimate_penalty_parameters_dragon(X1, X2)
    assert np.allclose(lambdas, expected, atol=1e-5, rtol=0), (
        "lambdas drift from gold: got=%s expected=%s" % (lambdas, expected)
    )


def test_dragon_parity_matrices():
    X1 = _load("X1.csv")
    X2 = _load("X2.csv")
    lambdas = tuple(np.loadtxt(os.path.join(FIXTURE_DIR, "lambdas.txt")))

    cov_gold  = _load("cov.csv")
    prec_gold = _load("prec.csv")
    ggm_gold  = _load("ggm.csv")

    cov  = dragon.get_shrunken_covariance_dragon(X1, X2, lambdas)
    prec, _ = dragon.get_precision_matrix_dragon(X1, X2, lambdas)
    ggm  = dragon.get_partial_correlation_dragon(X1, X2, lambdas)

    assert np.allclose(cov,  cov_gold,  atol=1e-5, rtol=0), (
        "cov max|diff|=%g" % np.abs(cov - cov_gold).max()
    )
    assert np.allclose(prec, prec_gold, atol=1e-5, rtol=0), (
        "prec max|diff|=%g" % np.abs(prec - prec_gold).max()
    )
    assert np.allclose(ggm,  ggm_gold,  atol=1e-5, rtol=0), (
        "ggm max|diff|=%g" % np.abs(ggm - ggm_gold).max()
    )
