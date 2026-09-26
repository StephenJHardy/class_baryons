"""Analyse a u-profile computed by src/profile_u.py (quadfit mode).

Delta chi^2(u) = 2 [minuslogpost(u) - min_u minuslogpost]. Because u >= 0 is a
physical boundary and the question is a one-sided upper limit, the 95% limit
is taken where Delta chi^2 = 2.71 (one-sided, Wilks with a boundary at u = 0);
the two-sided 3.84 crossing is reported too. The M5 Fisher forecast
(Delta chi^2 = (u / sigma_u)^2 with sigma_u = 6.7e-5, Planck-like) is overlaid.

Usage: uv run python src/profile_analyse.py plik_lite
"""

import json
import sys

import matplotlib.pyplot as plt
import numpy as np

from cloud_mapping import sigma_over_m, sigma_surface_over_q
from cosmology import ROOT

OUT = ROOT / "results" / "likelihood"
FISHER_SIGMA_U = 6.7e-5          # M5, Planck-like, marginalised
PUBLISHED = {"Planck 2015 TT,TE,EE+lowTEB (Stadler & Boehm 2018)": 1.58e-4,
             "Planck 2015 TT+lowTEB (Stadler & Boehm 2018)": 2.25e-4}


def crossing(u, dchi2, level):
    order = np.argsort(u)
    u, d = np.asarray(u)[order], np.asarray(dchi2)[order]
    above = np.nonzero(d >= level)[0]
    if len(above) == 0 or above[0] == 0:
        return None
    i = above[0]
    return float(np.interp(level, [d[i - 1], d[i]], [u[i - 1], u[i]]))


def main(likelihoods):
    rows = []
    for path in sorted((OUT / likelihoods / "quadfit").glob("u*/result.json")):
        r = json.load(open(path))
        rows.append({"u": r["u"], "minuslogpost": r["minuslogpost"], "err": r["minuslogpost_err"],
                     "fit_noise": r["fit_noise"], "n_evaluations": r["n_evaluations"],
                     "iterations": len(r["history"]), "best": r["best"]})
    rows.sort(key=lambda r: r["u"])
    u = np.array([r["u"] for r in rows])
    m = np.array([r["minuslogpost"] for r in rows])
    err = np.array([r["err"] for r in rows])
    dchi2 = 2 * (m - m.min())
    # quadratic fit in u (Delta chi^2 = a (u - u0)^2 + b) as a smooth summary
    coef = np.polyfit(u, 2 * m, 2, w=1 / np.maximum(2 * err, 0.01))
    u_fit_min = -coef[1] / (2 * coef[0])
    sigma_fit = 1 / np.sqrt(coef[0])
    limits = {}
    for level in (2.71, 3.84):
        uu = crossing(u, dchi2, level)
        limits[f"u(dchi2={level})"] = uu
        if uu:
            limits[f"sigma_over_M(dchi2={level})"] = sigma_over_m(uu)
            limits[f"Sigma_over_Q_min(dchi2={level})"] = sigma_surface_over_q(uu)
    summary = {"likelihoods": likelihoods, "rows": rows, "delta_chi2": dchi2.tolist(),
               "quadratic_fit": {"u_min": float(u_fit_min), "sigma_u": float(sigma_fit)},
               "limits": limits, "fisher_sigma_u_M5": FISHER_SIGMA_U, "published": PUBLISHED}
    with open(OUT / likelihoods / "profile_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=float)

    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    ax.errorbar(u, dchi2, yerr=2 * err, fmt="o", label=f"profile ({likelihoods}), quadratic-fit minimiser")
    uu = np.linspace(0, u.max(), 200)
    ax.plot(uu, np.polyval(coef, uu) - 2 * m.min(), "C0-", lw=0.8, label=f"quadratic fit: σ(u) = {sigma_fit:.2g}")
    ax.plot(uu, (uu / FISHER_SIGMA_U) ** 2, "k--", lw=0.8, label=f"M5 Fisher forecast (σ = {FISHER_SIGMA_U:.1e})")
    for level, ls in ((2.71, ":"), (3.84, "-.")):
        ax.axhline(level, color="0.5", lw=0.6, ls=ls)
    for (label, val), color in zip(PUBLISHED.items(), ("C3", "C1")):
        ax.axvline(val, color=color, lw=0.8, label=f"published: {label}")
    ax.set(xlabel="u", ylabel="Δχ² (profile)", title="Profile likelihood in the clump–photon coupling (f_cl = 1)")
    ax.legend(fontsize=7)
    fig.savefig(ROOT / "figures" / f"profile_u_{likelihoods}.png", dpi=150)
    print(json.dumps({k: v for k, v in summary.items() if k != "rows"}, indent=1, default=float))
    for r, d in zip(rows, dchi2):
        print(f"u = {r['u']:.2e}: -logP = {r['minuslogpost']:.3f} +- {r['err']:.3f}, dchi2 = {d:.2f}, "
              f"{r['iterations']} iterations, {r['n_evaluations']} evals")


if __name__ == "__main__":
    main(sys.argv[1])
