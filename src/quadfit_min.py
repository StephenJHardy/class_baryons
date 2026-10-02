"""Noise-robust minimisation of -log(posterior) by iterated quadratic regression.

CLASS likelihoods have a small discontinuous noise (~0.03 in -log L for plik,
0.1-0.2 for CamSpec, from parameter-dependent sampling grids), which makes
trust-region minimisers such as BOBYQA stop early and scatter by ~0.5 in -log L
between restarts. Near the minimum the posterior is close to Gaussian, so we:
  1. work in whitened coordinates z = L^-1 (x - x_c), with L the Cholesky
     factor of a covariance;
  2. evaluate -log P at N random points z ~ N(0, s^2 I) around the centre,
     drawn so that every point lies inside the parameter bounds;
  3. fit f(z) = c + g.z + z.A.z / 2 by least squares (averaging the noise);
  4. move the centre to the minimum of the fitted quadratic inside the bounds
     and a trust region, and iterate;
  5. on convergence, do a final larger fit, report the fitted minimum and
     evaluate the objective directly at the reported point.

Revision of 2026-10 (minimiser audit). The first version clipped the sampled
points x = c + L z to the parameter bounds before evaluating them, but regressed
on the unclipped z, so whenever clipping occurred the regression coordinates did
not describe the evaluated points. Now:
  - proposals are drawn from the Gaussian truncated to the bounds (rejection
    sampling), so the regression always uses the coordinates that were evaluated;
  - steps minimise the fitted quadratic subject to the bounds and the trust
    region (boundary minima are found as such, and centres stay feasible);
  - failed evaluations and rank-deficient designs raise QuadfitError instead of
    producing an underdetermined fit;
  - the uncertainty of the fitted minimum is the regression prediction variance
    at the minimum, which includes the fitted curvature (by the envelope theorem
    d f_min / d beta = d q / d beta at the minimiser);
  - the objective is evaluated directly at the reported point and both values
    are returned;
  - a rejected final fit gives converged = False.
The regression error is a statistical error of the fit only. It does not
include optimisation error (a wrong basin, a non-quadratic posterior) or the
correlated part of the numerical noise; use independent seeds or starts to
assess those.

The core (quadfit_core) takes any callable objective, so it can be tested on
known functions (tests/test_quadfit_min.py). quadfit_minimise wraps it for a
Cobaya model evaluated in parallel worker processes.
"""

import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.optimize import LinearConstraint, NonlinearConstraint, minimize

OUTLIER = 30.0          # points more than this above the sampled minimum are not fitted
MIN_POINTS_PER_COEF = 1.5
ROBUST_CUT = 4.0        # residual outliers beyond this many robust sigma are dropped from the fit
REWHITEN_CLIP = (0.5, 2.0)   # fitted curvature may change the sampling variance by at most 2x per iteration


class QuadfitError(RuntimeError):
    """The fit could not be carried out reliably (failed evaluations, rank-deficient design...)."""


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


def fit_quadratic(z, f, outlier=OUTLIER):
    """Least-squares quadratic f(z) = c + g.z + z.A.z/2 through finite values f at points z.

    Returns c, g, A, the residual rms, the coefficient covariance, the number of points used
    and the number rejected as residual outliers.
    Raises QuadfitError if the design is rank deficient or has too few points.
    """
    finite = np.isfinite(f)
    keep = finite & (f < np.min(f[finite]) + outlier)
    z, f = z[keep], f[keep]
    X, pairs = quad_design(z)
    n_coef = X.shape[1]
    if len(f) < MIN_POINTS_PER_COEF * n_coef:
        raise QuadfitError(f"only {len(f)} usable points for {n_coef} coefficients")
    rank = np.linalg.matrix_rank(X)
    if rank < n_coef:
        raise QuadfitError(f"design matrix has rank {rank} < {n_coef} (duplicate or degenerate points)")
    coef, *_ = np.linalg.lstsq(X, f, rcond=None)
    # robust passes: the CLASS likelihood noise has heavy tails (a few evaluations per hundred lie
    # 4-5 sigma out), so drop points beyond ROBUST_CUT robust sigma (from the median absolute
    # residual) and refit, as long as enough points remain
    n_all = len(f)
    for _ in range(3):
        res = f - X @ coef
        sig = 1.4826 * np.median(np.abs(res - np.median(res))) * np.sqrt(len(f) / max(len(f) - n_coef, 1))
        good = np.abs(res) <= ROBUST_CUT * max(sig, 1e-12)
        if good.all() or good.sum() < MIN_POINTS_PER_COEF * n_coef:
            break
        z, f, X = z[good], f[good], X[good]
        if np.linalg.matrix_rank(X) < n_coef:
            raise QuadfitError("design matrix rank deficient after outlier rejection")
        coef, *_ = np.linalg.lstsq(X, f, rcond=None)
    d = z.shape[1]
    c, g = coef[0], coef[1:1 + d]
    A = np.zeros((d, d))
    for k, (i, j) in enumerate(pairs):
        A[i, j] = A[j, i] = coef[1 + d + k]
    resid = f - X @ coef
    sigma2 = float(resid @ resid / (len(f) - n_coef))
    cov_coef = sigma2 * np.linalg.inv(X.T @ X)
    return float(c), g, A, float(np.sqrt(sigma2)), cov_coef, int(len(f)), int(n_all - len(f))


def predict(z, c, g, A, cov_coef):
    """Value of the fitted quadratic at z and its regression standard error."""
    row, _ = quad_design(np.atleast_2d(z))
    value = c + g @ z + 0.5 * z @ A @ z
    return float(value), float(np.sqrt(row[0] @ cov_coef @ row[0]))


def sample_feasible(rng, centre, L, scale, n, lo, hi, max_factor=2000):
    """n points z ~ N(0, scale^2 I) restricted to lo < centre + L z < hi (rejection sampling)."""
    d = len(centre)
    zs, drawn = [], 0
    while len(zs) < n:
        z = rng.normal(0, scale, size=(max(n, 64), d))
        x = centre + z @ L.T
        ok = np.all((x > lo) & (x < hi), axis=1)
        zs.extend(z[ok])
        drawn += len(z)
        if drawn > max_factor * n:
            raise QuadfitError(f"feasible region too small: {len(zs)} of {drawn} proposals inside the bounds")
    zs = np.array(zs[:n])
    return zs, len(zs) / drawn if drawn else 0.0


def constrained_minimum(g, A, centre, L, lo, hi, radius):
    """Minimise g.z + z.A.z/2 subject to lo <= centre + L z <= hi and |z| <= radius.

    Returns the minimiser, and whether the trust-region constraint is active.
    """
    d = len(g)

    def q(z):
        return g @ z + 0.5 * z @ A @ z

    def dq(z):
        return g + A @ z

    cons = [LinearConstraint(L, lo - centre, hi - centre),
            NonlinearConstraint(lambda z: z @ z, 0.0, radius**2, jac=lambda z: 2 * z[None, :])]
    starts = [np.zeros(d)]
    try:   # the unconstrained Newton step, pulled into the feasible region, is a good second start
        newton = -np.linalg.solve(A, g)
        if np.all(np.isfinite(newton)):
            t = min(1.0, 0.99 * radius / max(np.linalg.norm(newton), 1e-300))
            x = centre + L @ (t * newton)
            if np.all((x >= lo) & (x <= hi)):
                starts.append(t * newton)
    except np.linalg.LinAlgError:
        pass
    best = None
    for z0 in starts:
        res = minimize(q, z0, jac=dq, method="SLSQP", constraints=cons,
                       options={"maxiter": 500, "ftol": 1e-12})
        z = res.x
        x = centre + L @ z
        feasible = np.all(x >= lo - 1e-9 * (hi - lo)) and np.all(x <= hi + 1e-9 * (hi - lo)) \
            and np.linalg.norm(z) <= radius * (1 + 1e-6)
        if feasible and (best is None or q(z) < q(best)):
            best = z
    if best is None:
        best = np.zeros(d)
    return best, bool(np.linalg.norm(best) > 0.98 * radius)


def active_bounds(x, lo, hi, names, tol=1e-4):
    span = hi - lo
    return {n: ("lower" if xi - l < tol * s else "upper")
            for n, xi, l, h, s in zip(names, x, lo, hi, span) if xi - l < tol * s or h - xi < tol * s}


def quadfit_core(evaluate, names, lo, hi, cov, x0, scale=0.5, n_points=None, max_iter=8, step_cap=2.0,
                 tol_decrease=0.01, seed=0, log=print, clean_step_cap=6.0, clean_noise=0.05,
                 max_failed=0.05, checkpoint=None):
    """Minimise evaluate(x) inside the box lo < x < hi.

    evaluate maps an array of points (n, d) to an array of n objective values (non-finite =
    failed evaluation). Returns a dict; see the module docstring. Raises QuadfitError if a
    fit cannot be made reliably.
    """
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    margin = 1e-9 * (hi - lo)
    lo_s, hi_s = lo + margin, hi - margin          # sampling region: strictly inside the box
    d = len(names)
    n_coef = 1 + d + d * (d + 1) // 2
    n_points = n_points or 3 * n_coef
    rng = np.random.default_rng(seed)
    centre = np.clip(np.asarray(x0, float), lo_s, hi_s)
    L = np.linalg.cholesky(np.asarray(cov, float))
    history, n_eval, t0 = [], 0, time.time()

    def sample_and_evaluate(n, centre, L, scale):
        nonlocal n_eval
        zs, acceptance = sample_feasible(rng, centre, L, scale, n, lo_s, hi_s)
        xs = centre + zs @ L.T
        f = np.asarray(evaluate(xs), float)
        n_eval += len(xs)
        failed = float((~np.isfinite(f)).mean())
        if failed > max_failed:
            raise QuadfitError(f"{100 * failed:.0f}% of evaluations failed (non-finite objective)")
        return zs, xs, f, acceptance, failed

    converged = False
    for it in range(max_iter + 1):
        # sample; if the outlier cut would leave the fit underdetermined (the scale is far too
        # wide for the posterior, e.g. a poor starting covariance), recentre and shrink
        for _attempt in range(8):
            zs, xs, f, acceptance, failed = sample_and_evaluate(n_points, centre, L, scale)
            ok = np.isfinite(f)
            n_kept = int(np.sum(f[ok] < np.min(f[ok]) + OUTLIER))
            if n_kept >= MIN_POINTS_PER_COEF * n_coef:
                break
            centre = xs[np.argmin(np.where(ok, f, np.inf))]      # an evaluated, feasible point
            scale /= 2
            log(f"  only {n_kept} of {len(f)} points within {OUTLIER:g} of the minimum: recentre, scale -> {scale:.3g}")
        else:
            raise QuadfitError("could not find a sampling scale with enough points near the minimum")
        c, g, A, noise, cov_coef, n_used, n_rejected = fit_quadratic(zs, f)
        eig = np.linalg.eigvalsh(A)
        pos_def = bool(eig[0] > 0)
        # trust region: a clean, positive-definite fit may take a longer step
        radius = clean_step_cap if (pos_def and noise < clean_noise) else step_cap
        step, at_trust = constrained_minimum(g, A, centre, L, lo, hi, radius)
        decrease = float(-(g @ step + 0.5 * step @ A @ step))
        # the predicted decrease is itself uncertain at the level of the regression error of the
        # fitted value, so a tolerance below that makes the centre random-walk and never converge
        # (seen at f_cl = 1, u = 1e-4 with CamSpec: decreases of 0.03-0.05 for several iterations)
        threshold = max(tol_decrease, 2.0 * float(np.sqrt(cov_coef[0, 0])))
        prev_noise = [h["fit_noise"] for h in history]
        anomalous = bool(prev_noise) and noise > 3 * float(np.median(prev_noise)) + 0.01
        new_centre = np.clip(centre + L @ step, lo_s, hi_s)
        history.append({"iteration": it, "min_f_sampled": float(f[ok].min()), "fit_c": c,
                        "predicted_decrease": decrease, "step_norm": float(np.linalg.norm(step)),
                        "fit_noise": noise, "pos_def": pos_def, "at_trust_region": at_trust,
                        "n_used": n_used, "n_rejected": n_rejected, "threshold": threshold,
                        "anomalous_noise": anomalous, "failed_fraction": failed, "acceptance": acceptance, "scale": scale,
                        "active_bounds": active_bounds(new_centre, lo, hi, names)})
        log(f"  iter {it}: step {np.linalg.norm(step):.3f} sigma, predicted decrease {decrease:.4f}, "
            f"noise {noise:.3f}, min sampled {f[ok].min():.3f}"
            + (f", at bounds: {sorted(history[-1]['active_bounds'])}" if history[-1]["active_bounds"] else ""))
        if anomalous:
            # a fit far noisier than the earlier ones is not trusted: keep the centre and the whitening
            log("    fit noise anomalous: centre and whitening kept, resampling")
            continue
        centre = new_centre
        if pos_def:
            # re-whiten with the fitted curvature (posterior covariance = L A^-1 L^T), but let sigma
            # change by at most 2x per iteration: noisy, weakly constrained directions would
            # otherwise blow the sampling scale up (seen at full-plik u = 2e-5)
            w, V = np.linalg.eigh(A)
            cov_new = L @ (V @ np.diag(1 / np.clip(w, *REWHITEN_CLIP)) @ V.T) @ L.T
            L = np.linalg.cholesky(0.5 * (cov_new + cov_new.T))
        if checkpoint:
            import json
            with open(checkpoint, "w") as fh:
                json.dump({"names": list(names), "centre": centre.tolist(), "covariance": (L @ L.T).tolist(),
                           "history": history}, fh)
        if pos_def and not at_trust and decrease < threshold:
            converged = True
            break

    # final, larger fit centred on the last centre; the reported minimum is the minimum of this
    # fit inside the bounds and within one sampling scale of the centre
    zs, xs, f, acceptance, failed = sample_and_evaluate(2 * n_points, centre, L, scale)
    ok = np.isfinite(f)
    c, g, A, noise, cov_coef, n_used, n_rejected = fit_quadratic(zs, f)
    eig = np.linalg.eigvalsh(A)
    step, at_trust = constrained_minimum(g, A, centre, L, lo, hi, 2.0)
    final_decrease = float(-(g @ step + 0.5 * step @ A @ step))
    iter_noise = float(np.median([h["fit_noise"] for h in history]))
    reasons = []
    if not converged:
        reasons.append("iterations did not converge")
    if eig[0] <= 0:
        reasons.append("final Hessian not positive definite")
    if noise > 3 * iter_noise + 0.01:
        reasons.append(f"final fit noise {noise:.3f} vs {iter_noise:.3f} in the iterations")
    if at_trust:
        reasons.append("final minimum at the trust-region boundary")
    final_threshold = max(tol_decrease, 2.0 * float(np.sqrt(cov_coef[0, 0])))
    if final_decrease > 5 * final_threshold:
        reasons.append(f"final fit predicts a further decrease of {final_decrease:.3f}")
    final_ok = not [r for r in reasons if r != "iterations did not converge"]
    if final_ok:
        best = np.clip(centre + L @ step, lo, hi)
        z_best = step
    else:
        log("  WARNING: final fit rejected (" + "; ".join(reasons) + "); reporting the centre, not converged")
        best, z_best = centre, np.zeros(d)
    f_fit, f_fit_err = predict(z_best, c, g, A, cov_coef)
    direct = float(np.asarray(evaluate(best[None, :]), float)[0])
    n_eval += 1
    i_low = int(np.argmin(np.where(ok, f, np.inf)))
    w, V = np.linalg.eigh(A)
    cov_final = L @ (V @ np.diag(1 / np.clip(w, 0.25, 4.0)) @ V.T) @ L.T if eig[0] > 0 else L @ L.T
    return {"names": list(names), "best": dict(zip(names, best.tolist())),
            "minuslogpost": f_fit, "minuslogpost_err": f_fit_err,           # regression estimate at `best`
            "minuslogpost_fit": f_fit, "minuslogpost_fit_err": f_fit_err,
            "minuslogpost_direct": direct,                                  # objective evaluated at `best`
            "lowest_sampled": {"minuslogpost": float(f[i_low]), "point": dict(zip(names, xs[i_low].tolist()))},
            "converged": bool(converged and final_ok), "final_fit_ok": bool(final_ok), "status_reasons": reasons,
            "active_bounds": active_bounds(best, lo, hi, names),
            "fit_noise": noise, "final_predicted_decrease": final_decrease, "final_step_norm": float(np.linalg.norm(step)),
            "final_n_used": n_used, "final_n_rejected": n_rejected, "final_failed_fraction": failed, "final_acceptance": acceptance,
            "hessian_whitened": A.tolist(), "hessian_min_eigenvalue": float(eig[0]),
            "covariance": (0.5 * (cov_final + cov_final.T)).tolist(),
            "bounds": {n: [float(a), float(b)] for n, a, b in zip(names, lo, hi)},
            "history": history, "n_evaluations": int(n_eval), "seed": seed, "scale": scale,
            "runtime_s": time.time() - t0}


# ---------------------------------------------------------------- Cobaya wrapper

_MODEL = None


def _init_worker(info_dict, class_build):
    global _MODEL
    if class_build:
        from likelihood_setup import select_class_build
        select_class_build(class_build)
    from cobaya.model import get_model
    _MODEL = get_model(info_dict)


def _eval(point):
    try:
        return -float(_MODEL.logposterior(point).logpost)
    except Exception:  # noqa: BLE001 - CLASS failure: a failed evaluation, counted by the caller
        return np.inf


def _eval_detail(point):
    """Objective at one point with its pieces, the derived parameters and the loaded CLASS module."""
    import sys
    out = {"classy_path": getattr(sys.modules.get("classy"), "__file__", None)}
    try:
        lp = _MODEL.logposterior(point)
        out.update({"minuslogpost": -float(lp.logpost), "minuslogprior": -float(np.sum(lp.logpriors)),
                    "chi2": {name: -2 * float(v) for name, v in zip(_MODEL.likelihood, lp.loglikes)},
                    "derived": dict(zip(_MODEL.parameterization.derived_params(), map(float, lp.derived)))})
    except Exception as err:  # noqa: BLE001
        out.update({"minuslogpost": np.inf, "error": repr(err)})
    return out


def quadfit_minimise(info_dict, covmat_file, x0, workers=None, threads_per_worker=None, covmat=None,
                     class_build=None, **kwargs):
    """Minimise -log(posterior) of the Cobaya model `info_dict`, starting at the dict x0.

    The covariance for the whitening comes from `covmat` (array, in sampled-parameter order)
    or `covmat_file`. Other keyword arguments go to quadfit_core.
    """
    from cobaya.model import get_model
    if class_build:
        from likelihood_setup import select_class_build
        select_class_build(class_build)
    probe = get_model(info_dict)
    names = list(probe.parameterization.sampled_params())
    bounds = np.array(probe.prior.bounds(confidence_for_unbounded=0.9999995))
    probe.close()
    cov = np.array(covmat) if covmat is not None else load_covmat(covmat_file, names)
    # QF_WORKERS x QF_THREADS should roughly equal the machine's vCPUs (default: 4 x 4 for 16 vCPUs)
    workers = workers or int(os.environ.get("QF_WORKERS", 4))
    threads_per_worker = threads_per_worker or int(os.environ.get("QF_THREADS", 4))
    os.environ["OMP_NUM_THREADS"] = str(threads_per_worker)
    with ProcessPoolExecutor(workers, initializer=_init_worker, initargs=(info_dict, class_build)) as pool:
        def evaluate(xs):
            return np.array(list(pool.map(_eval, [dict(zip(names, x)) for x in xs])))

        res = quadfit_core(evaluate, names, bounds[:, 0], bounds[:, 1], cov, [x0[n] for n in names], **kwargs)
        res["direct_detail"] = pool.submit(_eval_detail, res["best"]).result()
    return res
