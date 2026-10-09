"""Numerical checks of the strong-coupling profile points (both kernels).

At the saved best fits of selected profile points, evaluate chi^2 with CLASS settings varied:
  default         the production settings
  tca_off_early   tight coupling switched off ten times earlier
  tca_first_order CLASS's first-order tight-coupling scheme
  precision_up    finer perturbation sampling and integration tolerance, l_linstep 10
Each is evaluated at the best fit and at 4 nearby points (0.3 sigma), and the variant minus
default difference is averaged, which suppresses the ~0.1 evaluation noise. The kernel difference
(isotropic - thomson at identical parameters) is evaluated the same way.

Usage: uv run python src/mix_strong_coupling_check.py
Output: results/mix/validation/strong_coupling_check.json
"""

import json
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from cosmology import ROOT
from likelihood_setup import info, select_class_build

OUT = ROOT / "results" / "mix"
POINTS = [(1.0, 2e-4), (0.1, 4e-3), (0.1, 6e-3), (0.03, 0.02), (0.03, 0.045), (0.01, 0.06), (0.01, 0.135), (0.01, 0.2)]
VARIANTS = {
    "default": {},
    "tca_off_early": {"tight_coupling_trigger_tau_c_over_tau_h": 0.0015, "tight_coupling_trigger_tau_c_over_tau_k": 0.001},
    "tca_first_order": {"tight_coupling_approximation": 2},
    "precision_up": {"l_linstep": 10, "perturbations_sampling_stepsize": 0.05, "tol_perturbations_integration": 1e-6},
}


def available(f_cl, u):
    return [k for k in ("thomson", "isotropic")
            if (OUT / "camspec_npipe" / f"f{f_cl:g}" / k / f"u{u:.3e}" / "result.json").exists()]


def job(args):
    f_cl, u, source, kernel, variant, points = args
    select_class_build("patched")
    from cobaya.model import get_model
    model = get_model(info("camspec_npipe", u=u, f_cl=f_cl, kernel=kernel, extra_classy=VARIANTS[variant]))
    vals = []
    for p in points:
        try:
            vals.append(-2 * float(model.logposterior(p).logpost))
        except Exception:  # noqa: BLE001
            vals.append(float("nan"))
    model.close()
    return {"f_cl": f_cl, "u": u, "points_from": source, "kernel": kernel, "variant": variant, "chi2": vals}


def main():
    rng = np.random.default_rng(5)
    jobs = []
    for f_cl, u in POINTS:
        for source in available(f_cl, u):                # best fit of this kernel's profile point
            r = json.load(open(OUT / "camspec_npipe" / f"f{f_cl:g}" / source / f"u{u:.3e}" / "result.json"))
            names = r["names"]
            L = np.linalg.cholesky(np.array(r["covariance"]))
            x = np.array([r["best"][n] for n in names])
            pts = [dict(zip(names, p)) for p in [x] + [x + L @ rng.normal(0, 0.3, len(names)) for _ in range(4)]]
            for variant in VARIANTS:
                jobs.append((f_cl, u, source, source, variant, pts))
            other = "isotropic" if source == "thomson" else "thomson"
            jobs.append((f_cl, u, source, other, "default", pts))
    os.environ["OMP_NUM_THREADS"] = "4"
    with ProcessPoolExecutor(4) as pool:
        rows = list(pool.map(job, jobs))
    key = {(r["f_cl"], r["u"], r["points_from"], r["kernel"], r["variant"]): np.array(r["chi2"]) for r in rows}
    summary = []
    for f_cl, u in POINTS:
        for source in available(f_cl, u):
            base = key[(f_cl, u, source, source, "default")]
            entry = {"f_cl": f_cl, "u": u, "kernel": source}
            for variant in VARIANTS:
                if variant == "default":
                    continue
                d = key[(f_cl, u, source, source, variant)] - base
                entry[variant] = {"mean": float(np.nanmean(d)), "sem": float(np.nanstd(d, ddof=1) / np.sqrt(len(d)))}
            other = "isotropic" if source == "thomson" else "thomson"
            d = key[(f_cl, u, source, other, "default")] - base
            entry[f"{other}_minus_{source}"] = {"mean": float(np.nanmean(d)), "sem": float(np.nanstd(d, ddof=1) / np.sqrt(len(d)))}
            summary.append(entry)
            print(f"f={f_cl:<5g} u={u:<7g} {source:9s} " + "  ".join(
                f"{k}: {v['mean']:+.3f}+-{v['sem']:.3f}" for k, v in entry.items() if isinstance(v, dict)), flush=True)
    path = OUT / "validation" / ("strong_coupling_check.json" if len(__import__("sys").argv) == 1
                                 else "strong_coupling_check_grid2d.json")
    with open(path, "w") as fh:
        json.dump({"variants": VARIANTS, "summary": summary, "rows": rows}, fh, indent=1)


if __name__ == "__main__":
    import sys
    select_class_build("patched")
    if len(sys.argv) > 1:     # e.g. 0.005:0.6 0.001:2 (kernels with a result at that point)
        POINTS[:] = [(float(a), float(b)) for a, b, *_ in (x.split(":") for x in sys.argv[1:])]
    main()
