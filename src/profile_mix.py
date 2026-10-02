"""Profile likelihoods for a mixture of collisionless CDM and photon-coupled clumps.

Model: omega_idm = f_cl * omega_dm, omega_cdm = (1 - f_cl) * omega_dm, with omega_dm fitted.
For each fixed clump fraction f_cl, scattering kernel and coupling u, minimise

    -log P = sum over likelihoods of (-log L) - log(prior)

over the six cosmological parameters and all nine CamSpec nuisance parameters (none is
pinned). The prior term is Planck's Gaussian priors on nuisance parameters plus constant
flat-prior normalisations; it is the same function at every (f_cl, kernel, u).

Results (one directory per point, resumable):
    results/mix/<likelihoods>/f<f>/u0/                      u = 0 (kernel-independent)
    results/mix/<likelihoods>/f<f>/<kernel>/u<u>[_<tag>]/   result.json | failure.json | checkpoint.json

Usage:
    uv run python src/profile_mix.py run <f_cl> <thomson|isotropic> <u> [<u> ...]
           [--seed N] [--tag NAME] [--start neighbour|lcdm|mirror] [--likelihoods camspec_npipe]
    uv run python src/profile_mix.py status

--start neighbour (default): warm start from the nearest completed point of the same fraction.
--start lcdm: independent start from the old LCDM (u = 0) solution and its covariance.
--start mirror: (isotropic) reuse the start point, covariance and seed of the Thomson fit at
    the same (f_cl, u), so that the two kernels are fitted on the same sample points and the
    numerical noise largely cancels in their difference.
Always uses the patched CLASS build with the kernel coefficients passed explicitly.
"""

import argparse
import json
import math
import time
import traceback

import numpy as np

from cosmology import ROOT
from likelihood_setup import KERNELS, class_build_identity, info, select_class_build
from quadfit_min import QuadfitError, quadfit_minimise

OUT = ROOT / "results" / "mix"
LCDM_START = ROOT / "results" / "likelihood" / "{likelihoods}" / "quadfit" / "u0.000e+00" / "result.json"
OBJECTIVE = ("-log(posterior) = sum of -log L over the likelihoods minus the log prior (Gaussian priors on "
             "CamSpec nuisance parameters; flat priors contribute constants)")


def point_dir(likelihoods, f_cl, kernel, u, tag=None):
    base = OUT / likelihoods / f"f{f_cl:g}"
    d = base / "u0" if u == 0 else base / kernel / f"u{u:.3e}"
    return d.with_name(f"{d.name}_{tag}") if tag else d


def completed(likelihoods, f_cl):
    """All untagged, converged results for this fraction: list of (u, kernel or None, result)."""
    out = []
    for path in (OUT / likelihoods / f"f{f_cl:g}").glob("**/result.json"):
        if "_" in path.parent.name:           # tagged repeat runs are not used as starts
            continue
        res = json.load(open(path))
        if res.get("converged"):
            out.append((res["u"], res.get("kernel"), res))
    return out


def choose_start(likelihoods, f_cl, kernel, u, mode):
    """Start point, covariance and seed, with a description of where they came from."""
    if mode == "mirror":
        res = json.load(open(point_dir(likelihoods, f_cl, "thomson", u) / "result.json"))
        return res["start"]["x0"], np.array(res["start"]["covariance"]), res["seed"], \
            f"mirror of thomson f={f_cl:g} u={u:g}"
    if mode == "neighbour":
        cands = completed(likelihoods, f_cl)
        if cands:
            def distance(c):
                cu, ck, _ = c
                du = abs(math.log10(max(cu, 1e-7)) - math.log10(max(u, 1e-7)))
                return du + (0.0 if ck in (kernel, None) else 0.05)
            cu, ck, res = min(cands, key=distance)
            return res["best"], np.array(res["covariance"]), None, f"neighbour f={f_cl:g} kernel={ck} u={cu:g}"
        # no point for this fraction yet: the u = 0 solution of the nearest fraction, if any
        others = [p for p in (OUT / likelihoods).glob("f*/u0/result.json")]
        others = [json.load(open(p)) for p in others]
        others = [r for r in others if r.get("converged")]
        if others:
            res = min(others, key=lambda r: abs(math.log10(max(r["f_cl"], 1e-9)) - math.log10(max(f_cl, 1e-9))))
            return res["best"], np.array(res["covariance"]), None, f"u=0 solution of f={res['f_cl']:g}"
    res = json.load(open(str(LCDM_START).format(likelihoods=likelihoods)))
    return res["best"], np.array(res["covariance"]), None, "earlier f=1, u=0 profile solution (LCDM-equivalent)"


def run_point(likelihoods, f_cl, kernel, u, seed=0, tag=None, start_mode="neighbour"):
    directory = point_dir(likelihoods, f_cl, kernel, u, tag)
    directory.mkdir(parents=True, exist_ok=True)
    label = f"f = {f_cl:g}, {kernel if u > 0 else 'u = 0'}, u = {u:g}" + (f" [{tag}]" if tag else "")
    if (directory / "result.json").exists():
        print(f"{label}: already done", flush=True)
        return json.load(open(directory / "result.json"))
    x0, cov, mirror_seed, origin = choose_start(likelihoods, f_cl, kernel, u, start_mode)
    if mirror_seed is not None:
        seed = mirror_seed
    if (directory / "checkpoint.json").exists() and start_mode != "mirror":
        ck = json.load(open(directory / "checkpoint.json"))
        x0, cov = dict(zip(ck["names"], ck["centre"])), np.array(ck["covariance"])
        origin += " (resumed from checkpoint)"
    inp = info(likelihoods, u=u, f_cl=f_cl, kernel=kernel if u > 0 else None)
    print(f"{label}: start from {origin}", flush=True)
    t0 = time.time()
    try:
        res = quadfit_minimise(inp, None, x0, covmat=cov, seed=seed, class_build="patched",
                               log=lambda m: print(m, flush=True), checkpoint=str(directory / "checkpoint.json"))
    except (QuadfitError, Exception) as err:  # noqa: BLE001 - record every failure; a failure is not an exclusion
        with open(directory / f"failure_{int(t0)}.json", "w") as fh:
            json.dump({"f_cl": f_cl, "kernel": kernel, "u": u, "seed": seed, "tag": tag, "start": origin,
                       "error": repr(err), "traceback": traceback.format_exc(),
                       "kind": "quadfit" if isinstance(err, QuadfitError) else "other"}, fh, indent=2)
        print(f"{label}: FAILED ({err!r})", flush=True)
        return None
    identity = class_build_identity("patched")
    loaded = res["direct_detail"].get("classy_path")
    if loaded != identity["classy_path"]:
        raise ImportError(f"workers loaded classy from {loaded}, expected {identity['classy_path']}")
    names = res["names"]
    res.update({"likelihoods": likelihoods, "f_cl": f_cl, "kernel": kernel if u > 0 else None, "u": u, "tag": tag,
                "objective": OBJECTIVE, "pinned_parameters": {},
                "kernel_coefficients": KERNELS[kernel] if u > 0 else None,
                "start": {"origin": origin, "mode": start_mode, "x0": {n: float(x0[n]) for n in names},
                          "covariance": np.asarray(cov).tolist()},
                "class": identity, "classy_extra_args": inp["theory"]["classy"]["extra_args"],
                "cobaya_input": json.loads(json.dumps(inp, default=str))})
    with open(directory / "result.json", "w") as fh:
        json.dump(res, fh, indent=2, default=float)
    print(f"{label}: -log P fit = {res['minuslogpost_fit']:.3f} +- {res['minuslogpost_fit_err']:.3f} (regression), "
          f"direct = {res['minuslogpost_direct']:.3f}, converged = {res['converged']}, "
          f"{res['n_evaluations']} evaluations, {res['runtime_s'] / 60:.1f} min", flush=True)
    return res


def status(likelihoods):
    for fdir in sorted((OUT / likelihoods).glob("f*"), key=lambda p: -float(p.name[1:])):
        rows = []
        for path in fdir.glob("**/result.json"):
            r = json.load(open(path))
            rows.append((r["u"], r["kernel"] or "-", r.get("tag") or "", r["minuslogpost_fit"], r["minuslogpost_fit_err"],
                         r["minuslogpost_direct"], r["converged"], sorted(r["active_bounds"])))
        fails = len(list(fdir.glob("**/failure_*.json")))
        print(f"{fdir.name}: {len(rows)} results, {fails} failures")
        for row in sorted(rows):
            print("   u = %-9.3g %-9s %-8s fit %.3f +- %.3f  direct %.3f  converged %s  bounds %s" % row)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=["run", "status"])
    ap.add_argument("f_cl", nargs="?", type=float)
    ap.add_argument("kernel", nargs="?", choices=list(KERNELS))
    ap.add_argument("u", nargs="*", type=float)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--tag")
    ap.add_argument("--start", default="neighbour", choices=["neighbour", "lcdm", "mirror"])
    ap.add_argument("--likelihoods", default="camspec_npipe")
    a = ap.parse_args()
    if a.command == "status":
        status(a.likelihoods)
    else:
        select_class_build("patched")
        for u in a.u:
            run_point(a.likelihoods, a.f_cl, a.kernel, u, seed=a.seed, tag=a.tag, start_mode=a.start)
