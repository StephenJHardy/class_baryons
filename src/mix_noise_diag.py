"""Diagnose the regression residuals of the quadratic fit at one profile point.

Samples points around a given centre exactly as the minimiser does, fits the quadratic and
reports the distribution of residuals, in total and per likelihood, to show whether the
occasional large fit noise comes from a few outlying evaluations.

Usage: uv run python src/mix_noise_diag.py <f_cl> <kernel> <u> <old result.json for centre/cov> [n=408]
Output: results/mix/validation/noise_diag_f<f>_<kernel>_u<u>.json
"""

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from cosmology import ROOT
from likelihood_setup import info, select_class_build
import quadfit_min as q


def main(f_cl, kernel, u, start, n=408):
    select_class_build("patched")
    r = json.load(open(start))
    names = r["names"]
    centre = np.array([r["best"][k] for k in names])
    L = np.linalg.cholesky(np.array(r["covariance"]))
    inp = info("camspec_npipe", u=u, f_cl=f_cl, kernel=kernel if u > 0 else None)
    rng = np.random.default_rng(123)
    z = rng.normal(0, 0.5, size=(n, len(names)))
    xs = centre + z @ L.T
    os.environ["OMP_NUM_THREADS"] = "4"
    with ProcessPoolExecutor(4, initializer=q._init_worker, initargs=(inp, "patched")) as pool:
        det = list(pool.map(q._eval_detail, [dict(zip(names, x)) for x in xs]))
    f = np.array([d["minuslogpost"] for d in det])
    parts = {k: 0.5 * np.array([d["chi2"][k] for d in det]) for k in det[0]["chi2"]}
    out = {"f_cl": f_cl, "kernel": kernel, "u": u, "n": n, "z": z.tolist(), "minuslogpost": f.tolist(),
           "parts": {k: v.tolist() for k, v in parts.items()}, "residuals": {}}
    X, _ = q.quad_design(z)
    for name, y in [("total", f)] + list(parts.items()):
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        res = y - X @ coef
        dof = len(y) - X.shape[1]
        rms = np.sqrt(res @ res / dof)
        mad = 1.4826 * np.median(np.abs(res - np.median(res))) * np.sqrt(len(y) / dof)
        out["residuals"][name] = res.tolist()
        print(f"{name:40s} rms {rms:.3f}  robust sigma {mad:.3f}  max|res| {np.abs(res).max():.3f}  "
              f"n(|res| > 4 robust sigma) = {int(np.sum(np.abs(res) > 4 * mad))}")
    path = ROOT / "results" / "mix" / "validation" / f"noise_diag_f{f_cl:g}_{kernel}_u{u:g}.json"
    with open(path, "w") as fh:
        json.dump(out, fh)


if __name__ == "__main__":
    main(float(sys.argv[1]), sys.argv[2], float(sys.argv[3]), sys.argv[4], *(int(a) for a in sys.argv[5:6]))
