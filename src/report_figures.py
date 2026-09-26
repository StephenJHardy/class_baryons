"""Figures for the results report (doc/report/).

  figures/report_profiles.pdf : Delta chi^2(u) profiles for Planck 2018 plik (full and
      lite) and PR4 CamSpec, with the 95% one-sided level and the published bounds
  figures/report_limits.pdf   : summary of 95% upper limits on u (and Sigma/Q)
  figures/report_mcmc.pdf     : MCMC marginal posterior of u, and posterior samples of u
      against the parameters it correlates with

Usage: uv run python src/report_figures.py
"""

import json

import matplotlib.pyplot as plt
import numpy as np

from cloud_mapping import sigma_surface_over_q, u_from_sigma_surface
from cosmology import ROOT

LIKE = ROOT / "results" / "likelihood"
FIG = ROOT / "figures"
PROFILES = [("plik", "Planck 2018 low-ℓ + plik TT,TE,EE", "C0", "o"),
            ("plik_lite", "Planck 2018 low-ℓ + plik-lite TT,TE,EE", "C2", "s"),
            ("camspec_npipe", "Planck 2018 low-ℓ + PR4 CamSpec TT,TE,EE", "C3", "D")]
PUBLISHED = [("Stadler & Bœhm 2018: Planck 2015 TT,TE,EE", 1.58e-4, "0.3", "--"),
             ("Stadler & Bœhm 2018: Planck 2015 TT", 2.25e-4, "0.3", ":"),
             ("Zhou et al. 2022: Planck 2018 (lite, τ fixed)", 1.55e-4, "C1", "--"),
             ("Zhou et al. 2022: Planck 2018 + lensing", 1.90e-4, "C1", ":")]
FISHER_SIGMA_U = 6.7e-5


def load(name):
    s = json.load(open(LIKE / name / "profile_summary.json"))
    u = np.array([r["u"] for r in s["rows"]])
    err = np.array([2 * max(r["err"], 0.025) for r in s["rows"]])   # >= fit noise floor
    return u, np.array(s["delta_chi2"]), err, s["limits"]


def sigma_axis(ax):
    top = ax.secondary_xaxis("top", functions=(lambda u: sigma_surface_over_q(np.maximum(u, 1e-12)),
                                                lambda s: u_from_sigma_surface(np.maximum(s, 1e-12))))
    return top


def profiles():
    fig, ax = plt.subplots(figsize=(7, 4.6), constrained_layout=True)
    for name, label, color, marker in PROFILES:
        u, d, err, lim = load(name)
        ax.errorbar(u * 1e4, d, yerr=err, fmt=marker + "-", color=color, ms=4, lw=0.9, capsize=2,
                    label=f"{label}: $u_{{95}}$ = {lim['u(dchi2=2.71)'] * 1e4:.2f}×10⁻⁴")
    uu = np.linspace(0, 4.2e-4, 200)
    ax.plot(uu * 1e4, (uu / FISHER_SIGMA_U) ** 2, color="0.6", lw=0.8, ls="-.",
            label="Fisher forecast (linear response at u = 0)")
    ax.axhline(2.71, color="k", lw=0.8)
    ax.text(4.1, 2.35, "95% one-sided (Δχ² = 2.71)", fontsize=8, ha="right", va="top")
    for label, val, color, ls in PUBLISHED:
        ax.axvline(val * 1e4, color=color, ls=ls, lw=1.0, label=label)
    ax.set(xlabel="u  [10⁻⁴]", ylabel="Δχ²  (profile)", xlim=(0, 4.2), ylim=(0, 17))
    ax.legend(fontsize=7.5, loc="upper left")
    fig.savefig(FIG / "report_profiles.pdf")
    fig.savefig(FIG / "report_profiles.png", dpi=150)


def mcmc(likelihoods="camspec_npipe"):
    from mcmc_analyse import load_combined
    _, comb = load_combined(likelihoods)
    summary = json.load(open(ROOT / "results" / "mcmc" / likelihoods / "mcmc_summary.json"))
    u95, u_prof = summary["u95"], load(likelihoods)[3]["u(dchi2=2.71)"]
    u = comb["u_idm_g"]
    rng = np.random.default_rng(1)
    pick = rng.choice(len(u), size=6000, p=comb.weights / comb.weights.sum())
    fig, axes = plt.subplots(2, 3, figsize=(10, 6.2), constrained_layout=True)
    ax = axes[0, 0]
    density = comb.get1DDensity("u_idm_g", boundaries=True)
    xs = np.linspace(0, 4e-4, 300)
    ax.plot(xs * 1e4, density(xs) / density(xs).max(), color="C3")
    ax.axvline(u95 * 1e4, color="C3", ls="--", lw=0.9, label=f"MCMC 95%: {u95 * 1e4:.2f}")
    ax.axvline(u_prof * 1e4, color="C3", ls=":", lw=0.9, label=f"profile 95%: {u_prof * 1e4:.2f}")
    ax.axvline(1.55, color="C1", ls="--", lw=0.9, label="Zhou+22: 1.55")
    ax.set(xlabel="u  [10⁻⁴]", ylabel="P(u) / P$_{max}$", xlim=(0, 4), ylim=(0, 1.05))
    ax.legend(fontsize=7.5, title="PR4 CamSpec, u₉₅ [10⁻⁴]", title_fontsize=7.5)
    for ax, (name, label) in zip(axes.flat[1:], [("sigma8", "σ₈"), ("theta_s_100", "100 θ$_s$"),
                                                 ("omega_dm", "ω$_{dm}$"), ("omega_b", "ω$_b$"), ("H0", "H₀")]):
        r = summary["correlation_with_u"][name]
        ax.scatter(u[pick] * 1e4, comb[name][pick], s=1.5, alpha=0.25, color="C0", rasterized=True)
        ax.axvline(u95 * 1e4, color="C3", ls="--", lw=0.9)
        ax.set(xlabel="u  [10⁻⁴]", ylabel=label, xlim=(0, 4))
        ax.text(0.97, 0.95, f"r = {r:+.2f}", transform=ax.transAxes, ha="right", va="top", fontsize=8)
    fig.savefig(FIG / "report_mcmc.pdf", dpi=200)
    fig.savefig(FIG / "report_mcmc.png", dpi=150)


def limits():
    rows = [("Fisher forecast, Planck-like (M5)", 1.1e-4, "0.5"),
            ("Profile: plik-lite", load("plik_lite")[3]["u(dchi2=2.71)"], "C2"),
            ("Profile: plik (full)", load("plik")[3]["u(dchi2=2.71)"], "C0"),
            ("Profile: PR4 CamSpec", load("camspec_npipe")[3]["u(dchi2=2.71)"], "C3"),
            ("MCMC (flat prior): PR4 CamSpec", json.load(open(ROOT / "results" / "mcmc" / "camspec_npipe"
                                                               / "mcmc_summary.json"))["u95"], "C4"),
            ("Published: Planck 2015 TT,TE,EE", 1.58e-4, "0.3"),
            ("Published: Planck 2015 TT", 2.25e-4, "0.3"),
            ("Published: Planck 2018 (lite, τ fixed)", 1.55e-4, "C1"),
            ("Published: Planck 2018 + lensing", 1.90e-4, "C1"),
            ("Forecast: SO + Planck (M5)", 1.2e-5, "0.5")]
    fig, ax = plt.subplots(figsize=(7, 4.3), constrained_layout=True)
    for i, (label, val, color) in enumerate(rows):
        y = len(rows) - 1 - i
        forecast = label.startswith(("Fisher", "Forecast"))
        ax.plot([val, 1e-2], [y, y], color=color, lw=6, alpha=0.35 if forecast else 0.8, solid_capstyle="butt")
        ax.plot(val, y, "|", color=color, ms=14, mew=2)
        ax.text(val * 0.88, y, f"{label}  ", ha="right", va="center", fontsize=8)
        m, e = f"{val:.2e}".split("e")
        ax.text(val * 1.1, y + 0.22, f"{float(m):.2f}×10$^{{{int(e)}}}$", va="bottom", fontsize=7.5)
    ax.set(xscale="log", xlim=(2e-7, 1e-3), ylim=(-0.7, len(rows) - 0.3), yticks=[],
           xlabel="u  (shaded: excluded at 95%; pale bars are forecasts)")
    top = sigma_axis(ax)
    top.set_xlabel("Σ/Q  [g cm⁻²]")
    fig.savefig(FIG / "report_limits.pdf")
    fig.savefig(FIG / "report_limits.png", dpi=150)


if __name__ == "__main__":
    profiles()
    limits()
    mcmc()
