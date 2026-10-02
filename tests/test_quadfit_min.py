"""Tests of the quadratic-regression minimiser on objectives with known minima.

Run: uv run pytest tests/test_quadfit_min.py
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from quadfit_min import QuadfitError, constrained_minimum, fit_quadratic, quad_design, quadfit_core  # noqa: E402

QUIET = dict(log=lambda m: None)


def make_quadratic(d, seed, correlated=True):
    """Positive-definite quadratic 0.5 (x - x0)^T H (x - x0) + f0 with a random correlation structure."""
    rng = np.random.default_rng(seed)
    sig = 10 ** rng.uniform(-3, 0, d)                 # parameter scales over three decades
    if correlated:
        M = rng.normal(size=(d, d))
        C = M @ M.T + d * np.eye(d)
        C /= np.sqrt(np.outer(np.diag(C), np.diag(C)))
        # strengthen one pair to a 0.9 correlation
        C[0, 1] = C[1, 0] = 0.9
        C = 0.5 * (C + C.T) + 0.05 * np.eye(d)
    else:
        C = np.eye(d)
    cov = C * np.outer(sig, sig)
    H = np.linalg.inv(cov)
    x0 = rng.normal(size=d) * sig
    f0 = 1234.5

    def f(xs):
        dx = np.atleast_2d(xs) - x0
        return f0 + 0.5 * np.einsum("ni,ij,nj->n", dx, H, dx)

    return f, x0, f0, cov, sig


def pseudo_noise(xs, amplitude, sig):
    """Deterministic, discontinuous-looking noise of given rms (like the CLASS likelihood noise)."""
    phase = (np.atleast_2d(xs) / sig) @ (1e4 * np.arange(1, xs.shape[-1] + 1) ** 0.5)
    return amplitude * np.sqrt(2) * np.sin(phase)


def wide_bounds(x0, sig):
    return x0 - 50 * sig, x0 + 50 * sig


def test_interior_minimum_exact():
    f, x0, f0, cov, sig = make_quadratic(6, 1)
    lo, hi = wide_bounds(x0, sig)
    start = x0 + 3 * sig                                # 3 sigma away in every parameter
    res = quadfit_core(f, list("abcdef"), lo, hi, 2.0 * cov, start, **QUIET)
    best = np.array([res["best"][n] for n in "abcdef"])
    assert res["converged"] and not res["active_bounds"]
    assert np.allclose(best, x0, atol=1e-6 * sig.max(), rtol=0) or np.max(np.abs(best - x0) / sig) < 1e-5
    assert abs(res["minuslogpost_fit"] - f0) < 1e-8
    assert abs(res["minuslogpost_direct"] - f0) < 1e-8


def test_correlated_poor_starting_covariance():
    """A diagonal starting covariance for a strongly correlated objective still converges."""
    f, x0, f0, cov, sig = make_quadratic(8, 2)
    lo, hi = wide_bounds(x0, sig)
    res = quadfit_core(f, [f"p{i}" for i in range(8)], lo, hi, np.diag(np.diag(cov)), x0 + 2 * sig,
                       max_iter=12, **QUIET)
    best = np.array(list(res["best"].values()))
    assert res["converged"]
    assert np.max(np.abs(best - x0) / sig) < 1e-4
    assert abs(res["minuslogpost_direct"] - f0) < 1e-6


def test_boundary_minimum():
    """The unconstrained minimum lies outside the box: the reported point must be on the bound,
    at the constrained minimum of the other parameters, and never outside."""
    d = 5
    f, x0, f0, cov, sig = make_quadratic(d, 3)
    names = [f"p{i}" for i in range(d)]
    lo, hi = wide_bounds(x0, sig)
    lo[0] = x0[0] + 1.5 * sig[0]                        # true minimum excluded by 1.5 sigma
    H = np.linalg.inv(cov)
    # analytic constrained minimum: x_0 = lo_0, the rest minimise the quadratic given x_0
    rest = np.arange(1, d)
    dx0 = lo[0] - x0[0]
    x_true = x0.copy()
    x_true[0] = lo[0]
    x_true[rest] = x0[rest] - np.linalg.solve(H[np.ix_(rest, rest)], H[rest, 0] * dx0)
    f_true = f(x_true)[0]
    start = x_true + np.array([2.0] + [1.0] * (d - 1)) * sig
    seen = []

    def f_checked(xs):
        seen.append(np.array(xs))
        assert np.all(xs > lo) and np.all(xs < hi), "evaluated a point outside the bounds"
        return f(xs)

    res = quadfit_core(f_checked, names, lo, hi, cov, start, **QUIET)
    best = np.array([res["best"][n] for n in names])
    assert res["converged"]
    assert res["active_bounds"] == {"p0": "lower"}
    assert np.all(best >= lo) and np.all(best <= hi)
    assert np.max(np.abs(best - x_true) / sig) < 1e-3
    assert abs(res["minuslogpost_fit"] - f_true) < 1e-5
    assert abs(res["minuslogpost_direct"] - f_true) < 1e-5
    for h in res["history"]:                             # every accepted centre respected the bounds
        assert h["acceptance"] > 0


def test_proposals_crossing_bounds_interior_minimum():
    """Minimum 0.3 sigma inside a bound: many proposals would cross it. The old minimiser clipped
    them and regressed on the unclipped coordinates; here the fit must still be exact."""
    d = 6
    f, x0, f0, cov, sig = make_quadratic(d, 4)
    names = [f"p{i}" for i in range(d)]
    lo, hi = wide_bounds(x0, sig)
    lo[2] = x0[2] - 0.3 * sig[2]
    hi[4] = x0[4] + 0.2 * sig[4]
    res = quadfit_core(f, names, lo, hi, cov, np.clip(x0 + 0.5 * sig, lo, hi), **QUIET)
    best = np.array([res["best"][n] for n in names])
    assert res["converged"] and not res["active_bounds"]
    assert res["final_acceptance"] < 0.9                 # the bounds really did cut the proposals
    assert np.max(np.abs(best - x0) / sig) < 1e-4
    assert abs(res["minuslogpost_fit"] - f0) < 1e-7


def test_clipping_would_bias():
    """Regression on unclipped coordinates of clipped points (the old behaviour) is biased;
    this documents the size of the bug the truncated sampling removes."""
    d = 4
    f, x0, f0, cov, sig = make_quadratic(d, 5, correlated=False)
    L = np.linalg.cholesky(cov)
    rng = np.random.default_rng(0)
    lo = x0 - 50 * sig
    lo[0] = x0[0] - 0.2 * sig[0]
    z = rng.normal(0, 0.5, size=(600, d))
    x = x0 + z @ L.T
    x_clipped = np.maximum(x, lo)
    c_bad, g_bad, A_bad, noise_bad, *_ = fit_quadratic(z, f(x_clipped))
    z_actual = np.linalg.solve(L, (x_clipped - x0).T).T
    with pytest.raises(QuadfitError):                    # clipped points are degenerate: caught by the rank check
        fit_quadratic(z_actual[x[:, 0] < lo[0]], f(x_clipped)[x[:, 0] < lo[0]])
    assert noise_bad > 1e-3                              # the mismatched regression does not fit its own points
    assert abs(A_bad[0, 0] - 1.0) > 0.05                 # and misestimates the curvature in the clipped direction


@pytest.mark.parametrize("amplitude", [0.03, 0.15])
def test_noisy_objective(amplitude):
    """Noise like the CLASS likelihoods (0.03 plik, 0.15 CamSpec), 15 parameters as in the CamSpec fit:
    the fitted minimum is within a few regression errors and far closer than single evaluations."""
    d = 15
    f, x0, f0, cov, sig = make_quadratic(d, 6)
    names = [f"p{i}" for i in range(d)]
    lo, hi = wide_bounds(x0, sig)

    def noisy(xs):
        xs = np.atleast_2d(xs)
        return f(xs) + pseudo_noise(xs, amplitude, sig)

    results = [quadfit_core(noisy, names, lo, hi, 1.5 * cov, x0 + 1.0 * sig, seed=seed, **QUIET) for seed in range(3)]
    fits = np.array([r["minuslogpost_fit"] for r in results])
    errs = np.array([r["minuslogpost_fit_err"] for r in results])
    assert all(r["converged"] for r in results)
    assert np.all(np.abs(fits - f0) < 5 * errs + 0.01)   # consistent with the truth
    assert np.all(errs < 0.3 * amplitude)                # the fit averages the noise down
    assert np.ptp(fits) < 4 * errs.mean() + 0.005        # repeatable across seeds, as the quoted error says
    for r in results:                                    # fitted rms recovers the noise level
        assert 0.7 * amplitude < r["fit_noise"] < 1.4 * amplitude
        best = np.array([r["best"][n] for n in names])
        dx = best - x0
        true_excess = 0.5 * dx @ np.linalg.inv(cov) @ dx
        assert true_excess < 0.05                        # the reported point is near the true minimum


def test_fitted_minimum_error_includes_curvature():
    """Monte Carlo check of the quoted error on the fitted minimum value."""
    d, amplitude = 3, 0.1
    f, x0, f0, cov, sig = make_quadratic(d, 7)
    L = np.linalg.cholesky(cov)
    lo, hi = wide_bounds(x0, sig)
    rng = np.random.default_rng(11)
    centre = x0 + L @ np.array([0.4, -0.3, 0.2])         # the fit is not centred on the minimum
    mins, errs = [], []
    for _ in range(300):
        z = rng.normal(0, 0.5, size=(60, d))
        y = f(centre + z @ L.T) + rng.normal(0, amplitude, 60)
        c, g, A, noise, cov_coef, *_ = fit_quadratic(z, y)
        step, _ = constrained_minimum(g, A, centre, L, lo, hi, 5.0)
        row, _ = quad_design(step[None, :])
        mins.append(c + g @ step + 0.5 * step @ A @ step)
        errs.append(np.sqrt(row[0] @ cov_coef @ row[0]))
    scatter, quoted = np.std(mins), np.mean(errs)
    assert 0.8 < quoted / scatter < 1.25
    # the old formula (c and g terms only, with 0.5 * step) underestimates the error
    assert abs(np.mean(mins) - f0) < 4 * scatter / np.sqrt(300) + 0.005


def test_failed_evaluations_raise():
    f, x0, f0, cov, sig = make_quadratic(4, 8)
    lo, hi = wide_bounds(x0, sig)

    def failing(xs):
        out = f(xs)
        out[::3] = np.inf                                # a third of the evaluations fail
        return out

    with pytest.raises(QuadfitError, match="evaluations failed"):
        quadfit_core(failing, list("abcd"), lo, hi, cov, x0, **QUIET)


def test_few_failures_tolerated_and_counted():
    f, x0, f0, cov, sig = make_quadratic(4, 9)
    lo, hi = wide_bounds(x0, sig)

    def failing(xs):
        out = f(xs)
        out[::40] = np.inf                               # 2.5% fail
        return out

    res = quadfit_core(failing, list("abcd"), lo, hi, cov, x0 + sig, **QUIET)
    assert res["converged"]
    assert 0 < res["history"][0]["failed_fraction"] < 0.05
    assert abs(res["minuslogpost_fit"] - f0) < 1e-7


def test_non_convergence_is_reported():
    """Too few iterations from far away: the result must say it did not converge."""
    f, x0, f0, cov, sig = make_quadratic(5, 10)
    lo, hi = wide_bounds(x0, sig)
    res = quadfit_core(f, list("abcde"), lo, hi, cov, x0 + 30 * sig, max_iter=0, step_cap=1.0,
                       clean_step_cap=1.0, **QUIET)
    assert not res["converged"]
    assert res["status_reasons"]


def test_degenerate_design_raises():
    z = np.zeros((50, 3))
    z[:, 0] = np.linspace(-1, 1, 50)                     # all points on a line
    with pytest.raises(QuadfitError, match="rank"):
        fit_quadratic(z, z[:, 0] ** 2)


def test_heavy_tailed_noise_from_far_start():
    """The case that broke the first production attempt: 15 parameters, CamSpec-like noise (0.1) with
    a few per cent of large outliers, started ~3 sigma from the minimum. Every seed must converge,
    without the centre random-walking, and agree with the truth within the quoted error."""
    d, amplitude = 15, 0.1
    f, x0, f0, cov, sig = make_quadratic(d, 6)
    names = [f"p{i}" for i in range(d)]
    lo, hi = wide_bounds(x0, sig)
    L = np.linalg.cholesky(cov)
    start = x0 + L @ (3.0 * np.ones(d) / np.sqrt(d))     # 3 sigma (Mahalanobis) away

    for seed in range(5):
        noise_rng = np.random.default_rng(100 + seed)

        def noisy(xs):
            xs = np.atleast_2d(xs)
            e = noise_rng.normal(0, amplitude, len(xs))
            out = noise_rng.random(len(xs)) < 0.03
            e[out] += noise_rng.choice([-1, 1], out.sum()) * noise_rng.uniform(0.5, 3.0, out.sum())
            return f(xs) + e

        res = quadfit_core(noisy, names, lo, hi, cov, start, seed=seed, **QUIET)
        best = np.array([res["best"][n] for n in names])
        dx = best - x0
        assert res["converged"], res["status_reasons"]
        assert len(res["history"]) <= 6
        assert 0.5 * dx @ np.linalg.inv(cov) @ dx < 0.1
        assert abs(res["minuslogpost_fit"] - f0) < 5 * res["minuslogpost_fit_err"] + 0.02
        assert 0.7 * amplitude < res["fit_noise"] < 1.5 * amplitude     # outliers did not inflate the noise
