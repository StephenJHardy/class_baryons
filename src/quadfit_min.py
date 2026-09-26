"""Noise-robust minimisation of -log(posterior) by iterated quadratic regression.

CLASS likelihoods have a small discontinuous noise (~0.03 in -log L, from
parameter-dependent sampling grids), which makes trust-region minimisers such as
BOBYQA stop early and scatter by ~0.5 in -log L between restarts. Near the
minimum the posterior is close to Gaussian, so instead we:
  1. work in whitened coordinates z = L^-1 (x - x_c), with L the Cholesky
     factor of a covariance (the M5 Fisher covmat, u fixed);
  2. evaluate -log P at N random points z ~ N(0, s^2 I) around the centre;
  3. fit f(z) = c + g.z + z.A.z / 2 by least squares (averaging the noise);
  4. move the centre to the fitted minimum (step capped), and iterate;
  5. on convergence, do a final larger fit and report the fitted minimum,
     its statistical uncertainty (from the regression) and the Hessian.
Evaluations run in parallel worker processes, each holding one Cobaya model.

Usage (library): quadfit_minimise(info_dict, covmat_file, x0) -> dict
"""

import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

_MODEL = None


def _init_worker(info_dict):
    global _MODEL
    from cobaya.model import get_model
    _MODEL = get_model(info_dict)


def _eval(point):
    try:
        return -float(_MODEL.logposterior(point).logpost)
    except Exception:  # noqa: BLE001 - CLASS failure: treat as very bad point
        return np.inf


def load_covmat(path, names):
    with open(path) as f:
        header = f.readline().lstrip("#").split()
    cov = np.loadtxt(path)
    idx = [header.index(n) for n in names]
    return cov[np.ix_(idx, idx)]


def quad_design(z):
    n, d = z.shape
    cols = [np.ones(n)] + [z[:, i] for i in range(d)]
    pairs = [(i, j) for i in range(d) for j in range(i, d)]
    cols += [z[:, i] * z[:, j] * (0.5 if i == j else 1.0) for i, j in pairs]
    return np.column_stack(cols), pairs


def fit_quadratic(z, f):
    X, pairs = quad_design(z)
    coef, *_ = np.linalg.lstsq(X, f, rcond=None)
    d = z.shape[1]
    c, g = coef[0], coef[1:1 + d]
    A = np.zeros((d, d))
    for k, (i, j) in enumerate(pairs):
        A[i, j] = A[j, i] = coef[1 + d + k]
    resid = f - X @ coef
    dof = max(len(f) - X.shape[1], 1)
    sigma2 = float(resid @ resid / dof)
    cov_coef = sigma2 * np.linalg.pinv(X.T @ X)
    return c, g, A, float(np.sqrt(sigma2)), cov_coef


def quadfit_minimise(info_dict, covmat_file, x0, scale=0.5, n_points=None, max_iter=8,
                     step_cap=2.0, tol_decrease=0.01, workers=4, threads_per_worker=4, seed=0, log=print,
                     covmat=None):
    from cobaya.model import get_model
    names = list(get_model(info_dict).parameterization.sampled_params())
    cov = np.array(covmat) if covmat is not None else load_covmat(covmat_file, names)
    L = np.linalg.cholesky(cov)
    d = len(names)
    n_points = n_points or 3 * (1 + d + d * (d + 1) // 2)
    rng = np.random.default_rng(seed)
    centre = np.array([x0[n] for n in names], dtype=float)
    os.environ["OMP_NUM_THREADS"] = str(threads_per_worker)
    history = []
    t0 = time.time()
    with ProcessPoolExecutor(workers, initializer=_init_worker, initargs=(info_dict,)) as pool:
        def evaluate(zs, c):
            pts = [dict(zip(names, c + L @ z)) for z in zs]
            return np.array(list(pool.map(_eval, pts)))

        for it in range(max_iter + 1):
            final = False
            zs = rng.normal(0, scale, size=(n_points, d))
            f = evaluate(zs, centre)
            ok = np.isfinite(f)
            c, g, A, noise, _ = fit_quadratic(zs[ok], f[ok])
            try:
                step = -np.linalg.solve(A, g)
                pos_def = bool(np.all(np.linalg.eigvalsh(A) > 0))
            except np.linalg.LinAlgError:
                step, pos_def = -g, False
            if not pos_def or np.linalg.norm(step) > step_cap:
                step = step / max(np.linalg.norm(step), 1e-12) * min(step_cap, np.linalg.norm(step))
            decrease = float(-0.5 * g @ step) if pos_def else np.inf
            history.append({"iteration": it, "min_f_sampled": float(f[ok].min()), "fit_c": float(c),
                            "predicted_decrease": decrease,
                            "step_norm": float(np.linalg.norm(step)), "fit_noise": noise, "pos_def": pos_def,
                            "n_failed": int((~ok).sum())})
            log(f"  iter {it}: fitted min at |step| = {np.linalg.norm(step):.3f} sigma, residual noise {noise:.3f}, "
                f"min sampled {f[ok].min():.3f}")
            centre = centre + L @ step
            if pos_def:
                # re-whiten with the fitted curvature: the true posterior covariance is L A^-1 L^T
                cov = L @ np.linalg.inv(A) @ L.T
                L = np.linalg.cholesky(0.5 * (cov + cov.T))
            if pos_def and decrease < tol_decrease:
                break
        # final, larger fit centred on the converged point
        zs = rng.normal(0, scale, size=(2 * n_points, d))
        f = evaluate(zs, centre)
        ok = np.isfinite(f)
        c, g, A, noise, cov_coef = fit_quadratic(zs[ok], f[ok])
        step = -np.linalg.solve(A, g)
        f_min = float(c + 0.5 * g @ step)
        # uncertainty of the fitted minimum value, propagated from the coefficient covariance
        grad = np.zeros(cov_coef.shape[0])
        grad[0] = 1.0
        grad[1:1 + d] = 0.5 * step
        f_min_err = float(np.sqrt(grad @ cov_coef @ grad))
        best = centre + L @ step
    cov_final = L @ np.linalg.inv(A) @ L.T
    return {"names": names, "best": dict(zip(names, best.tolist())), "minuslogpost": f_min,
            "covariance": cov_final.tolist(),
            "minuslogpost_err": f_min_err, "fit_noise": noise, "hessian_whitened": A.tolist(),
            "history": history, "n_evaluations": int(n_points * (len(history) + 2)), "runtime_s": time.time() - t0}
