"""Coarse BOBYQA pre-minimisation to give quadfit a starting point near the minimum.

Used for a cold start with a new likelihood whose nuisance reference values are far
from its best fit (e.g. CamSpec NPIPE). BOBYQA's ~0.5 noise-limited precision is
irrelevant here; quadfit then polishes the result.

Usage: uv run python src/coarse_start.py <likelihoods> <u> [reference profile, default plik]
Writes results/likelihood/<likelihoods>/coarse_start_u<u>.json
"""

import json
import sys
import tempfile

import numpy as np
from cobaya.model import get_model
from cobaya.run import run

from likelihood_setup import PACKAGES, info
from profile_u import NUISANCE_COVMATS, OUT


def main(likelihoods, u, reference="plik"):
    probe = get_model(info(likelihoods, u=u))
    names = list(probe.parameterization.sampled_params())
    pinfo = probe.parameterization.sampled_params_info()
    ref0 = json.load(open(OUT / reference / "quadfit" / "u0.000e+00" / "result.json"))
    rn, rc = ref0["names"], np.array(ref0["covariance"])
    path = f"{PACKAGES}/data/planck_supp_data_and_covmats/covmats/{NUISANCE_COVMATS[likelihoods]}"
    pn = open(path).readline().lstrip("#").split()
    pc = np.loadtxt(path)
    cov = np.zeros((len(names), len(names)))
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if a in rn and b in rn:
                cov[i, j] = rc[rn.index(a), rn.index(b)]
            elif a in pn and b in pn and a not in rn and b not in rn:
                cov[i, j] = pc[pn.index(a), pn.index(b)]
    covfile = tempfile.NamedTemporaryFile("w", suffix=".covmat", delete=False)
    np.savetxt(covfile.name, cov, header=" ".join(names), comments="# ")
    inp = info(likelihoods, u=u, sampler={"minimize": {"method": "bobyqa", "best_of": 1, "covmat": covfile.name,
                                                        "override_bobyqa": {"rhoend": 0.05, "maxfun": 3000}}})
    for n in names:                                 # start from the reference profile's cosmology
        if n in ref0["best"] and n in inp["params"]:   # nuisance params live in the likelihood definitions
            inp["params"][n]["ref"] = ref0["best"][n]
    _, sampler = run(inp, force=True, no_mpi=True)
    best = sampler.products()["minimum"].data.iloc[0].to_dict()
    out = {"likelihoods": likelihoods, "u": u, "minuslogpost": float(best["minuslogpost"]),
           "best": {n: float(best[n]) for n in names}, "covariance": cov.tolist(), "names": names}
    OUT.joinpath(likelihoods).mkdir(parents=True, exist_ok=True)
    with open(OUT / likelihoods / f"coarse_start_u{u:.3e}.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"coarse minimum -log P = {out['minuslogpost']:.3f}")
    print({n: round(v, 4) for n, v in out["best"].items()})


if __name__ == "__main__":
    main(sys.argv[1], float(sys.argv[2]), *(sys.argv[3:4] or []))
