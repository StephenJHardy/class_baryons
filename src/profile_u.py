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
        if (directory / "result.json").exists():      # resumable (e.g. after spot preemption)
            res = json.load(open(directory / "result.json"))
            start, cov = res["best"], res["covariance"]
            print(f"u = {u:g}: already done, skipping", flush=True)
            continue
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



# Planck covmats use CosmoMC names for the plik nuisance parameters.
from likelihood_setup import PACKAGES  # noqa: E402

PLANCK_COVMAT = f"{PACKAGES}/data/planck_supp_data_and_covmats/covmats/base_plikHM_TTTEEE_lowE.covmat"
COSMOMC_TO_COBAYA = {
    "calPlanck": "A_planck", "acib217": "A_cib_217", "xi": "xi_sz_cib", "asz143": "A_sz", "aksz": "ksz_norm",
    "aps100": "ps_A_100_100", "aps143": "ps_A_143_143", "aps143217": "ps_A_143_217", "aps217": "ps_A_217_217",
    "kgal100": "gal545_A_100", "kgal143": "gal545_A_143", "kgal143217": "gal545_A_143_217", "kgal217": "gal545_A_217",
    "galfTE100": "galf_TE_A_100", "galfTE100143": "galf_TE_A_100_143", "galfTE100217": "galf_TE_A_100_217",
    "galfTE143": "galf_TE_A_143", "galfTE143217": "galf_TE_A_143_217", "galfTE217": "galf_TE_A_217",
    "cal0": "calib_100T", "cal2": "calib_217T",
}


PINNED_NUISANCE = {"xi_sz_cib": 0.0, "ksz_norm": 0.0}   # at their lower prior bound, where Planck's posteriors pile up


def full_plik_check(us, seed=0):
    """Profile points with the full plik likelihood (21 nuisance parameters).

    Covariance: cosmological block from the refined lite covariance at u = 0 (with A_planck),
    nuisance block from Planck's base_plikHM_TTTEEE_lowE covmat (renamed), no cross terms.
    Start: the lite best fit at the same u for cosmology, Cobaya reference values for nuisance.
    """
    import numpy as np
    from cobaya.model import get_model
    from quadfit_min import quadfit_minimise
    probe = get_model(info("plik", u=0.0, fixed_params=PINNED_NUISANCE))
    names = list(probe.parameterization.sampled_params())
    ref_info = probe.parameterization.sampled_params_info()
    with open(PLANCK_COVMAT) as f:
        planck_names = [COSMOMC_TO_COBAYA.get(n, n) for n in f.readline().lstrip("#").split()]
    planck_cov = np.loadtxt(PLANCK_COVMAT)
    lite0 = json.load(open(OUT / "plik_lite" / "quadfit" / "u0.000e+00" / "result.json"))
    lite_names, lite_cov = lite0["names"], np.array(lite0["covariance"])
    cov = np.zeros((len(names), len(names)))
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if a in lite_names and b in lite_names:
                cov[i, j] = lite_cov[lite_names.index(a), lite_names.index(b)]
            elif a in planck_names and b in planck_names and a not in lite_names and b not in lite_names:
                cov[i, j] = planck_cov[planck_names.index(a), planck_names.index(b)]
    previous = None
    for u in us:
        lite = json.load(open(OUT / "plik_lite" / "quadfit" / f"u{u:.3e}" / "result.json"))["best"]
        if previous is None:
            # cold start: lite cosmology + Cobaya reference values for the nuisance parameters
            start = {n: lite.get(n, ref_info[n]["ref"]["loc"] if isinstance(ref_info[n].get("ref"), dict)
                                 else ref_info[n].get("ref")) for n in names}
        else:
            # warm start: previous full-plik solution, with the cosmology moved as the lite profile moves
            lite_prev = json.load(open(OUT / "plik_lite" / "quadfit" / f"u{previous['u']:.3e}" / "result.json"))["best"]
            start = {n: previous["best"][n] + (lite[n] - lite_prev[n] if n in lite else 0.0) for n in names}
            cov = np.array(previous["covariance"])
        directory = OUT / "plik" / "quadfit" / f"u{u:.3e}"
        directory.mkdir(parents=True, exist_ok=True)
        if (directory / "result.json").exists():      # resumable
            previous = json.load(open(directory / "result.json"))
            print(f"u = {u:g}: already done, skipping", flush=True)
            continue
        if (directory / "checkpoint.json").exists():  # resume a preempted minimisation
            ck = json.load(open(directory / "checkpoint.json"))
            start, cov = dict(zip(ck["names"], ck["centre"])), np.array(ck["covariance"])
            print(f"u = {u:g}: resuming from checkpoint", flush=True)
        print(f"u = {u:g} (full plik, {len(names)} parameters)", flush=True)
        res = quadfit_minimise(info("plik", u=u, fixed_params=PINNED_NUISANCE), None, start, covmat=cov, seed=seed,
                               log=lambda m: print(m, flush=True), checkpoint=str(directory / "checkpoint.json"))
        res.update({"u": u, "likelihoods": "plik", "pinned_nuisance": PINNED_NUISANCE})
        previous = res
        with open(directory / "result.json", "w") as f:
            json.dump(res, f, indent=2, default=float)
        print(f"  -> -log P = {res['minuslogpost']:.3f} +- {res['minuslogpost_err']:.3f} "
              f"({res['n_evaluations']} evaluations, {res['runtime_s']:.0f} s)", flush=True)



def profile_generic(likelihoods, us, reference="plik", pinned=None, seed=0):
    """Full-likelihood profile for any likelihood set, seeded from an existing profile.

    Cosmological start and covariance come from the `reference` profile (same u where
    available, else its u = 0); nuisance parameters start at Cobaya's reference values
    with a diagonal covariance from their proposal widths. Points are chained (warm
    starts) and resumable, as in full_plik_check.
    """
    import numpy as np
    from cobaya.model import get_model
    from quadfit_min import quadfit_minimise
    probe = get_model(info(likelihoods, u=0.0, fixed_params=pinned))
    names = list(probe.parameterization.sampled_params())
    pinfo = probe.parameterization.sampled_params_info()
    ref0 = json.load(open(OUT / reference / "quadfit" / "u0.000e+00" / "result.json"))
    rn, rc = ref0["names"], np.array(ref0["covariance"])
    cov = np.zeros((len(names), len(names)))
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if a in rn and b in rn:
                cov[i, j] = rc[rn.index(a), rn.index(b)]
        if names[i] not in rn:
            width = pinfo[names[i]].get("proposal") or (pinfo[names[i]].get("ref") or {}).get("scale", 1.0)
            cov[i, i] = float(width) ** 2

    def ref_value(n):
        r = pinfo[n].get("ref")
        return r["loc"] if isinstance(r, dict) else r

    previous = None
    for u in us:
        tag = f"u{u:.3e}"
        directory = OUT / likelihoods / "quadfit" / tag
        directory.mkdir(parents=True, exist_ok=True)
        if (directory / "result.json").exists():
            previous = json.load(open(directory / "result.json"))
            print(f"u = {u:g}: already done, skipping", flush=True)
            continue
        seed_file = OUT / reference / "quadfit" / tag / "result.json"
        cosmo = json.load(open(seed_file))["best"] if seed_file.exists() else ref0["best"]
        if previous is None:
            start = {n: cosmo.get(n, ref_value(n)) for n in names}
        else:
            start = dict(previous["best"])
            prev_seed = OUT / reference / "quadfit" / f"u{previous['u']:.3e}" / "result.json"
            if seed_file.exists() and prev_seed.exists():
                prev_cosmo = json.load(open(prev_seed))["best"]
                for n in names:
                    if n in cosmo and n in prev_cosmo:
                        start[n] += cosmo[n] - prev_cosmo[n]
            cov = np.array(previous["covariance"])
        if (directory / "checkpoint.json").exists():
            ck = json.load(open(directory / "checkpoint.json"))
            start, cov = dict(zip(ck["names"], ck["centre"])), np.array(ck["covariance"])
            print(f"u = {u:g}: resuming from checkpoint", flush=True)
        print(f"u = {u:g} ({likelihoods}, {len(names)} parameters)", flush=True)
        res = quadfit_minimise(info(likelihoods, u=u, fixed_params=pinned), None, start, covmat=cov, seed=seed,
                               log=lambda m: print(m, flush=True), checkpoint=str(directory / "checkpoint.json"))
        res.update({"u": u, "likelihoods": likelihoods, "pinned": pinned or {}})
        with open(directory / "result.json", "w") as f:
            json.dump(res, f, indent=2, default=float)
        print(f"  -> -log P = {res['minuslogpost']:.3f} +- {res['minuslogpost_err']:.3f} "
              f"({res['n_evaluations']} evaluations, {res['runtime_s']:.0f} s)", flush=True)
        previous = res


if __name__ == "__main__":
    if sys.argv[1] == "quadfit":
        quadfit_profile(sys.argv[2], [float(x) for x in sys.argv[3:]])
    elif sys.argv[1] == "generic":
        profile_generic(sys.argv[2], [float(x) for x in sys.argv[3:]])
    elif sys.argv[1] == "fullplik":
        full_plik_check([float(x) for x in sys.argv[2:]])
    else:
        for u in map(float, sys.argv[2:]):
            minimise(sys.argv[1], u)
