"""Kernel difference at identical parameters, along the Thomson profiles.

For every completed Thomson profile point (f_cl, u > 0), evaluate the objective with both
kernels at the Thomson best-fit parameters. Because the parameters (and so CLASS's internal
grids) are identical, the numerical noise of the likelihood cancels in the difference, and by
the envelope theorem the difference equals the difference of the two profile minima to first
order in the kernel perturbation. This gives the kernel dependence of the profile far more
precisely than subtracting two independent minimisations.

Usage: uv run python src/mix_paired_kernel.py [likelihoods=camspec_npipe]
Output: results/mix/<likelihoods>/paired_kernel.json  {f_cl: [{u, thomson, isotropic, delta_chi2}]}
"""

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

from cosmology import ROOT
from likelihood_setup import KERNELS, info, select_class_build

OUT = ROOT / "results" / "mix"


def evaluate(args):
    likelihoods, f_cl, kernel, points = args
    select_class_build("patched")
    from cobaya.model import get_model
    inp = info(likelihoods, u=1e-6, f_cl=f_cl, kernel=kernel)
    del inp["theory"]["classy"]["extra_args"]["u_idm_g"]
    inp["params"]["u_idm_g"] = {"prior": {"min": 0.0, "max": 1e3}}     # one model for all u (flat: constant prior)
    model = get_model(inp)
    out = []
    for u, best in points:
        lp = model.logposterior({**best, "u_idm_g": u})
        out.append(-2 * float(sum(lp.loglikes)) - 2 * float(sum(lp.logpriors)))   # chi2-like, up to a constant
    model.close()
    return f_cl, kernel, out


def main(likelihoods="camspec_npipe"):
    jobs, index = [], {}
    for fdir in (OUT / likelihoods).glob("f*"):
        pts = []
        for path in sorted((fdir / "thomson").glob("u*/result.json")):
            r = json.load(open(path))
            if r.get("tag") or not r["converged"]:
                continue
            pts.append((r["u"], r["best"]))
        if pts:
            f_cl = float(fdir.name[1:])
            index[f_cl] = [u for u, _ in pts]
            jobs += [(likelihoods, f_cl, k, pts) for k in KERNELS]
    os.environ["OMP_NUM_THREADS"] = "4"
    res = {}
    with ProcessPoolExecutor(4) as pool:
        for f_cl, kernel, vals in pool.map(evaluate, jobs):
            res.setdefault(f_cl, {})[kernel] = vals
    out = {f"{f:g}": [{"u": u, "thomson": t, "isotropic": i, "delta_chi2": i - t}
                      for u, t, i in zip(index[f], res[f]["thomson"], res[f]["isotropic"])] for f in index}
    with open(OUT / likelihoods / "paired_kernel.json", "w") as fh:
        json.dump(out, fh, indent=1)
    for f, rows in out.items():
        print(f"f_cl = {f}: " + ", ".join(f"u={r['u']:g}: {r['delta_chi2']:+.3f}" for r in rows))


if __name__ == "__main__":
    select_class_build("patched")
    main(*sys.argv[1:2])
