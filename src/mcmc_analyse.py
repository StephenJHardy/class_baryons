"""Analyse the headline MCMC chains written by src/mcmc_run.py.

  - loads every results/mcmc/<likelihoods>/chain_<i>.1.txt, discarding the first
    30% of each chain as burn-in (Cobaya's own R-1 stop uses the later part too)
  - Gelman-Rubin R-1 across chains for the mean and for the 95% bound on u
  - one-sided 95% (and 99%) upper limits on u, with Sigma/Q and sigma/M
  - marginalised constraints on the standard parameters, H0 and sigma8
  - figure: 1D posterior of u (with the profile limit marked) and a triangle plot

Usage: uv run python src/mcmc_analyse.py [likelihoods=camspec_npipe]
"""

import json
import sys

import matplotlib.pyplot as plt
import numpy as np
from getdist import MCSamples, loadMCSamples, plots

from cloud_mapping import sigma_over_m, sigma_surface_over_q
from cosmology import ROOT

OUT = ROOT / "results" / "mcmc"
BURN = 0.3
PARAMS = ["u_idm_g", "omega_b", "omega_dm", "theta_s_100", "n_s", "logA", "tau_reio", "H0", "sigma8"]


def upper_limit(samples, name, level):
    x = samples[name]
    w = samples.weights
    order = np.argsort(x)
    cdf = np.cumsum(w[order]) / w.sum()
    return float(np.interp(level, cdf, x[order]))


def load_combined(likelihoods="camspec_npipe"):
    """Per-chain samples (burn-in removed) and all chains combined into one MCSamples.

    Built from the arrays, because getCombinedSamplesWithSamples trips over Cobaya's
    dotted chi2__ parameter names."""
    chains = sorted((OUT / likelihoods).glob("chain_*.1.txt"))
    all_s = [loadMCSamples(str(c).removesuffix(".1.txt"), settings={"ignore_rows": BURN}) for c in chains]
    names = [p.name for p in all_s[0].paramNames.names]
    combined = MCSamples(samples=np.vstack([s.samples for s in all_s]),
                         weights=np.concatenate([s.weights for s in all_s]),
                         loglikes=np.concatenate([s.loglikes for s in all_s]),
                         names=names, labels=[p.label for p in all_s[0].paramNames.names],
                         ranges={"u_idm_g": [0, 1e-3]}, label=likelihoods, ignore_rows=0)
    return all_s, combined


def split_half_u95(all_s):
    """95% limit on u from the first and second halves of every chain (after burn-in)."""
    out = []
    for half in (0, 1):
        u, w = [], []
        for s in all_s:
            n = s.numrows // 2
            sl = slice(0, n) if half == 0 else slice(n, None)
            u.append(s["u_idm_g"][sl])
            w.append(s.weights[sl])
        u, w = np.concatenate(u), np.concatenate(w)
        order = np.argsort(u)
        out.append(float(np.interp(0.95, np.cumsum(w[order]) / w.sum(), u[order])))
    return out


def main(likelihoods="camspec_npipe"):
    root = OUT / likelihoods
    chains = sorted(root.glob("chain_*.1.txt"))
    per_chain = []
    for c in chains:
        s = loadMCSamples(str(c).removesuffix(".1.txt"), settings={"ignore_rows": BURN})
        per_chain.append({"chain": c.name, "n_accepted": int(s.numrows), "weight": float(s.weights.sum()),
                          "u_mean": float(s.mean("u_idm_g")), "u95": upper_limit(s, "u_idm_g", 0.95)})
    all_s, combined = load_combined(likelihoods)
    means = np.array([s.mean("u_idm_g") for s in all_s])
    variances = np.array([s.var("u_idm_g") for s in all_s])
    rminus1_mean = float(np.var(means, ddof=1) / np.mean(variances)) if len(all_s) > 1 else None
    u95s = np.array([p["u95"] for p in per_chain])
    u95 = upper_limit(combined, "u_idm_g", 0.95)
    u99 = upper_limit(combined, "u_idm_g", 0.99)
    summary = {
        "likelihoods": likelihoods, "burn_in_fraction": BURN, "chains": per_chain,
        "rminus1_u_mean": rminus1_mean,
        "u95_chain_scatter_over_mean": float(np.std(u95s, ddof=1) / np.mean(u95s)) if len(u95s) > 1 else None,
        "u95": u95, "u99": u99, "u95_split_halves": split_half_u95(all_s),
        "correlation_with_u": {p: float(combined.corr([combined.index["u_idm_g"], combined.index[p]])[0, 1])
                               for p in PARAMS[1:] if combined.paramNames.parWithName(p)},
        "Sigma_over_Q_min_95": sigma_surface_over_q(u95), "sigma_over_M_max_95": sigma_over_m(u95),
        "marginalised": {p: {"mean": float(combined.mean(p)), "std": float(combined.std(p))}
                         for p in PARAMS if combined.paramNames.parWithName(p)},
    }
    profile = ROOT / "results" / "likelihood" / likelihoods / "profile_summary.json"
    u_profile = json.load(open(profile))["limits"].get("u(dchi2=2.71)") if profile.exists() else None
    summary["profile_u95"] = u_profile
    with open(root / "mcmc_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=1))

    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    density = combined.get1DDensity("u_idm_g", boundaries=True)
    xs = np.linspace(0, max(4e-4, 1.5 * u99), 300)
    ax.plot(xs, density(xs) / density(xs).max(), "C0", label=f"MCMC posterior ({len(chains)} chains)")
    ax.axvline(u95, color="C0", ls="--", lw=0.8, label=f"95% upper limit {u95:.2e}")
    if u_profile:
        ax.axvline(u_profile, color="C3", ls=":", lw=0.8, label=f"profile Δχ² = 2.71: {u_profile:.2e}")
    ax.set(xlabel="u", ylabel="P(u) / P_max", xlim=(0, xs[-1]), ylim=(0, 1.05),
           title=f"Clump–photon coupling, f_cl = 1 ({likelihoods})")
    ax.legend(fontsize=8)
    fig.savefig(ROOT / "figures" / f"mcmc_u_{likelihoods}.png", dpi=150)

    g = plots.get_subplot_plotter(width_inch=9)
    g.triangle_plot(combined, [p for p in PARAMS if combined.paramNames.parWithName(p)], filled=True)
    g.export(str(ROOT / "figures" / f"mcmc_triangle_{likelihoods}.png"))


if __name__ == "__main__":
    main(*sys.argv[1:2])
