"""Fixed-parameter scan of chi^2 against coupling for every clump fraction and kernel.

All other parameters are held at the PR4 CamSpec u = 0 best fit of the earlier profile, so
chi^2(u) - chi^2(0) here is an upper bound on the profile Delta chi^2 (the profile re-fits
the other parameters). It is used to (i) find where CLASS fails or its approximations
break down, (ii) choose the coupling grids for the profiles, and (iii) test the
tight-coupling treatment of both kernels over the full coupling range.

Usage: uv run python src/mix_direct_scan.py [variant ...]
Variants: default, tca_first_order, tca_off_early   (default: all)
Output: results/mix/validation/direct_scan.json
"""

import json
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from cosmology import ROOT
from likelihood_setup import KERNELS, class_build_identity, info, select_class_build

OUT = ROOT / "results" / "mix" / "validation"
LIKE = "camspec_npipe"
POINT = ROOT / "results" / "likelihood" / LIKE / "quadfit" / "u0.000e+00" / "result.json"
FRACTIONS = [1.0, 0.5, 0.1, 0.03, 0.01]
U_GRID = np.round(10 ** np.arange(-6.0, 1.01, 0.25), 12)
VARIANTS = {
    "default": {},
    # CLASS allows only its own first-order scheme or the default "compromise" scheme with idm_g
    "tca_first_order": {"tight_coupling_approximation": 2},
    # leave tight coupling ten times earlier (defaults 0.015 and 0.01)
    "tca_off_early": {"tight_coupling_trigger_tau_c_over_tau_h": 0.0015,
                      "tight_coupling_trigger_tau_c_over_tau_k": 0.001},
}


def scan(args):
    f_cl, kernel, variant = args
    select_class_build("patched")
    from cobaya.model import get_model
    point = json.load(open(POINT))["best"]
    inp = info(LIKE, u=1e-6, f_cl=f_cl, kernel=kernel, extra_classy=VARIANTS[variant])
    # make u a parameter so that one model serves the whole grid
    del inp["theory"]["classy"]["extra_args"]["u_idm_g"]
    inp["params"]["u_idm_g"] = {"prior": {"min": 0.0, "max": 1e3}}
    model = get_model(inp)
    rows = []
    for u in U_GRID:
        p = {n: point[n] for n in model.parameterization.sampled_params() if n != "u_idm_g"}
        p["u_idm_g"] = float(u)
        try:
            lp = model.logposterior(p)
            chi2 = {n: -2 * float(v) for n, v in zip(model.likelihood, lp.loglikes)}
            derived = dict(zip(model.parameterization.derived_params(), map(float, lp.derived)))
            rows.append({"u": float(u), "chi2_total": float(sum(chi2.values())), "chi2": chi2,
                         "sigma8": derived.get("sigma8"), "H0": derived.get("H0")})
        except Exception as err:  # noqa: BLE001 - a CLASS failure is a result here
            rows.append({"u": float(u), "error": repr(err)[:400]})
    model.close()
    return {"f_cl": f_cl, "kernel": kernel, "variant": variant, "rows": rows}


def main(variants):
    jobs = [(f, k, v) for v in variants for f in FRACTIONS for k in KERNELS]
    path = OUT / "direct_scan.json"
    done = json.load(open(path))["scans"] if path.exists() else []
    have = {(s["f_cl"], s["kernel"], s["variant"]) for s in done}
    jobs = [j for j in jobs if j not in have]
    import os
    os.environ["OMP_NUM_THREADS"] = "3"
    with ProcessPoolExecutor(5) as pool:
        for res in pool.map(scan, jobs):
            done.append(res)
            ok = [r for r in res["rows"] if "error" not in r]
            print(f"f = {res['f_cl']:<5g} {res['kernel']:9s} {res['variant']:16s}: {len(ok)}/{len(res['rows'])} ok, "
                  f"largest u ok = {max(r['u'] for r in ok) if ok else None}", flush=True)
            OUT.mkdir(parents=True, exist_ok=True)
            with open(path, "w") as fh:
                json.dump({"identity": class_build_identity("patched"), "u_grid": U_GRID.tolist(), "scans": done}, fh)


if __name__ == "__main__":
    select_class_build("patched")
    main(sys.argv[1:] or list(VARIANTS))
