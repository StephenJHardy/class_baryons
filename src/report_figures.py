"""Figures for the results report (doc/report/).

  figures/report_profiles.pdf : Delta chi^2(u) profiles for Planck 2018 plik (full and
      lite) and PR4 CamSpec, with the 95% one-sided level and the published bounds
  figures/report_limits.pdf   : summary of 95% upper limits on u (and Sigma/Q)

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
             ("Stadler & Bœhm 2018: Planck 2015 TT", 2.25e-4, "0.3", ":")]
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


def limits():
    rows = [("Fisher forecast, Planck-like (M5)", 1.1e-4, "0.5"),
            ("Profile: plik-lite", load("plik_lite")[3]["u(dchi2=2.71)"], "C2"),
            ("Profile: plik (full)", load("plik")[3]["u(dchi2=2.71)"], "C0"),
            ("Profile: PR4 CamSpec", load("camspec_npipe")[3]["u(dchi2=2.71)"], "C3"),
            ("Published: Planck 2015 TT,TE,EE", 1.58e-4, "0.3"),
            ("Published: Planck 2015 TT", 2.25e-4, "0.3"),
            ("Forecast: SO + Planck (M5)", 1.2e-5, "0.5")]
    fig, ax = plt.subplots(figsize=(7, 3.4), constrained_layout=True)
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
