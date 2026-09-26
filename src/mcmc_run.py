"""Headline MCMC: Planck PR4 CamSpec + 2018 low-l TT/EE, f_cl = 1, flat prior on u.

Runs one resumable Cobaya chain (no MPI). Launch several with different chain
indices in parallel; convergence is then checked across chains with getdist
(src/mcmc_analyse.py) as well as by each chain's own R-1.

  - parameters: omega_b, omega_dm, 100 theta_s, n_s, ln 10^10 A_s, tau, u, + CamSpec nuisance
  - prior on u: flat in [0, 1e-3] (likelihood_setup.info(u_free=True))
  - proposal: covariance from the CamSpec profile at u = 0 (standard + nuisance
    parameters) with u added using the profile's curvature scale
  - fast-slow dragging (the nuisance parameters do not need a new CLASS run)

Usage: uv run python src/mcmc_run.py <chain index> [likelihoods=camspec_npipe]
Output: results/mcmc/<likelihoods>/chain_<i>.*   (Cobaya resumes if present)
"""

import json
import sys

import numpy as np
from cobaya.run import run

from cosmology import ROOT
from likelihood_setup import info

OUT = ROOT / "results" / "mcmc"
PROFILE = ROOT / "results" / "likelihood"
SIGMA_U_PROPOSAL = 5e-5   # from the u profiles (Delta chi^2 = 1 near u ~ 4-5e-5)


def proposal_covmat(likelihoods):
    prof = json.load(open(PROFILE / likelihoods / "quadfit" / "u0.000e+00" / "result.json"))
    names = prof["names"] + ["u_idm_g"]
    cov = np.zeros((len(names), len(names)))
    cov[:-1, :-1] = np.array(prof["covariance"])
    cov[-1, -1] = SIGMA_U_PROPOSAL**2
    path = OUT / likelihoods / "proposal.covmat"
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, cov, header=" ".join(names), comments="# ")
    return path, prof["best"]


def main(chain, likelihoods="camspec_npipe", max_samples=None):
    covmat, best = proposal_covmat(likelihoods)
    sampler = {"mcmc": {"covmat": str(covmat), "drag": True, "oversample_power": 0.4,
                        "learn_proposal": True, "Rminus1_stop": 0.01, "Rminus1_cl_stop": 0.2,
                        "max_tries": 10000, "seed": 1000 + chain}}
    if max_samples:
        sampler["mcmc"]["max_samples"] = max_samples
    inp = info(likelihoods, u_free=True, sampler=sampler,
               output=str(OUT / likelihoods / f"chain_{chain}"))
    # start near the profile's best fit. Nuisance parameters are defined inside the likelihood,
    # so copy their full definitions into the top-level block with the reference moved
    # (their default references can be far from the best fit, e.g. amp_143 = 10 vs 18.6).
    from cobaya.model import get_model
    probe = get_model(info(likelihoods, u_free=True))
    for name, pinfo in probe.parameterization.sampled_params_info().items():
        if name in best:
            definition = dict(inp["params"].get(name) or pinfo)
            scale = definition.get("proposal") or 1e-3
            definition["ref"] = {"dist": "norm", "loc": best[name], "scale": 0.5 * scale}
            inp["params"][name] = definition
    probe.close()
    run(inp, resume=True, no_mpi=True)


if __name__ == "__main__":
    import os
    main(int(sys.argv[1]), *(sys.argv[2:3] or []), max_samples=int(os.environ.get("MCMC_MAX_SAMPLES", 0)) or None)
