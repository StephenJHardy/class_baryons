"""Profile likelihood in u by minimisation (M7/M8 first step; no MCMC).

For each fixed u (f_cl = 1), minimise -log L over the standard parameters and
(for the full plik likelihood) the Planck nuisance parameters, with Cobaya's
BOBYQA minimiser started from several points. Results go to
results/likelihood/<likelihoods>/u<value>/ .

Usage: uv run python src/profile_u.py <plik|plik_lite> <u> [<u> ...]
"""

import json
import sys
import time

from cobaya.run import run

from cosmology import ROOT
from likelihood_setup import info

COVMAT = ROOT / "results" / "likelihood" / "fisher_planck_like_u_fixed.covmat"   # src/fisher_covmat.py

OUT = ROOT / "results" / "likelihood"


def minimise(likelihoods, u, seed_runs=4):
    tag = f"u{u:.3e}"
    directory = OUT / likelihoods / tag
    directory.mkdir(parents=True, exist_ok=True)
    sampler = {"minimize": {"method": "bobyqa", "best_of": seed_runs, "confidence_for_unbounded": 0.9999995,
                            "covmat": str(COVMAT), "override_bobyqa": {"rhoend": 0.005, "maxfun": 4000, "objfun_has_noise": True}}}
    t0 = time.time()
    updated, sampler_out = run(info(likelihoods, u=u, sampler=sampler, output=str(directory / "min")),
                               force=True, no_mpi=True)
    products = sampler_out.products()
    best = products["minimum"].data.iloc[0].to_dict()
    result = {"likelihoods": likelihoods, "u": u, "minuslogpost": float(best["minuslogpost"]),
              "all_minima": {k: float(v[0]) for k, v in products.get("full_set_of_mins", {}).items()},
              "chi2": {k: float(v) for k, v in best.items() if k.startswith("chi2")},
              "params": {k: float(v) for k, v in best.items() if not k.startswith(("chi2", "minuslog", "weight"))},
              "runtime_s": time.time() - t0}
    with open(directory / "result.json", "w") as f:
        json.dump(result, f, indent=2, default=float)
    print(json.dumps(result, indent=1, default=float))
    return result




def quadfit_profile(likelihoods, us, start=None, seed=0):
    """Profile over u with the noise-robust quadratic-fit minimiser (src/quadfit_min.py).

    Each u point starts from the previous point's best fit and covariance, so the
    whitening improves along the profile.
    """
    from quadfit_min import quadfit_minimise
    start = start or json.load(open(OUT / likelihoods / "u0.000e+00" / "result.json"))["params"]
    cov = None
    for u in us:
        tag = f"u{u:.3e}"
        directory = OUT / likelihoods / "quadfit" / tag
        directory.mkdir(parents=True, exist_ok=True)
        print(f"u = {u:g}", flush=True)
        res = quadfit_minimise(info(likelihoods, u=u), str(COVMAT), start, covmat=cov, seed=seed,
                               log=lambda m: print(m, flush=True))
        res["u"] = u
        res["likelihoods"] = likelihoods
        with open(directory / "result.json", "w") as f:
            json.dump(res, f, indent=2, default=float)
        print(f"  -> -log P = {res['minuslogpost']:.3f} +- {res['minuslogpost_err']:.3f} "
              f"({res['n_evaluations']} evaluations, {res['runtime_s']:.0f} s)", flush=True)
        start, cov = res["best"], res["covariance"]


if __name__ == "__main__":
    if sys.argv[1] == "quadfit":
        quadfit_profile(sys.argv[2], [float(x) for x in sys.argv[3:]])
    else:
        for u in map(float, sys.argv[2:]):
            minimise(sys.argv[1], u)
