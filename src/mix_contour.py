"""Map of the profile Delta chi^2 in the plane of clump fraction and clump surface density.

Each grid point is a full profile fit (src/profile_mix.py): all 15 parameters minimised at fixed
(f_cl, u). Delta chi^2 is relative to the common LCDM fit, so negative values fit better than
LCDM and the 2.71 / 3.84 contours mark couplings that worsen the fit by that much relative to
LCDM. This differs from the fixed-fraction thresholds (relative to each fraction's own minimum),
which are overplotted as markers. The surface is interpolated along lines of constant f_cl u (see surface())
and drawn only inside the sampled region.

Axes: f_cl (bottom) with omega_cl = f_cl * omega_dark (top), and Sigma/Q = 268/u g cm^-2.
Reference marks: Walker & Wardle (1998) clouds (M <~ 1e-3 Msun, R ~ 1-3 AU), and the present-day
masses of all stars and all cold gas as fractions of the cosmic baryons (approximate).

Usage: uv run python src/mix_contour.py
Output: figures/mix_contour_camspec_npipe.{png,pdf}
"""

import json

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap

from cloud_mapping import sigma_surface_over_q, u_from_sigma_surface
from cosmology import ROOT

LIKE = "camspec_npipe"
OUT = ROOT / "results" / "mix" / LIKE
FIG = ROOT / "figures"
INK, MUTED, SURFACE = "#1a1a19", "#6b6a63", "#fcfcfb"
OMEGA_DARK, OMEGA_B = 0.1199, 0.02213          # LCDM fit (f_cl = 0)
U_ZERO_PLOT = 2.68e-6                           # u = 0 points are drawn at Sigma/Q = 1e8 (the top edge)
LEVELS = [-4, -3, -2, -1, 0, 1, 2.71, 3.84, 6, 10, 20]
MSUN, AU = 1.989e33, 1.496e13


def load():
    lcdm = json.load(open(OUT / "f0" / "u0" / "result.json"))["minuslogpost_fit"]
    pts = {"thomson": [], "isotropic": []}
    for path in OUT.glob("f*/**/result.json"):
        r = json.load(open(path))
        if r.get("tag") or not r["converged"] or r["f_cl"] == 0:
            continue
        d = 2 * (r["minuslogpost_fit"] - lcdm)
        kernels = [r["kernel"]] if r["u"] > 0 else list(pts)       # u = 0 is shared by both kernels
        for k in kernels:
            pts[k].append((r["f_cl"], max(r["u"], U_ZERO_PLOT), d))
    return {k: np.array(sorted(v)) for k, v in pts.items()}


def cmap():
    # diverging: blue = better than LCDM, neutral grey at 0, orange = worse
    return LinearSegmentedColormap.from_list(
        "div", [(0.0, "#184f95"), (0.4, "#86b6ef"), (0.5, "#e6e5df"), (0.62, "#f4b393"), (1.0, "#a8361a")])


def thresholds():
    s = json.load(open(OUT / "summary.json"))["fractions"]
    out = {}
    for f, e in s.items():
        for k in ("thomson", "isotropic"):
            c = (e.get(k) or {}).get("crossings", {}).get("2.71")
            if c:
                out.setdefault(k, []).append((float(f), c["u"]))
    return out


def surface(data, logf, logsig):
    """Delta chi^2 on a regular (log f, log Sigma) grid.

    For each fitted fraction, Delta chi^2 is interpolated linearly in log(f_cl u) (u = 0 is taken
    as the value for all couplings below the smallest one fitted). Between fractions, it is
    interpolated linearly in log f_cl at fixed log(f_cl u), because the profiles line up in f_cl u.
    Points beyond the largest fitted coupling of either neighbouring fraction are left blank.
    """
    fr = np.unique(data[:, 0])
    curves = {}
    for f in fr:
        rows = data[data[:, 0] == f]
        rows = rows[np.argsort(rows[:, 1])]
        curves[f] = (np.log10(f * rows[:, 1]), rows[:, 2])
    lf = np.log10(fr)
    Z = np.full((len(logsig), len(logf)), np.nan)
    for j, x in enumerate(logf):
        if x < lf.min() - 1e-9 or x > lf.max() + 1e-9:
            continue
        i = min(np.searchsorted(lf, x), len(lf) - 1)
        lo, hi = (lf[i - 1], lf[i]) if i > 0 and lf[i] > x else (lf[i], lf[i])
        for k, y in enumerate(logsig):
            lfu = x + np.log10(u_from_sigma_surface(10 ** y))       # log(f u) at this (f, Sigma)
            vals = []
            for l in (lo, hi):
                cx, cy = curves[fr[np.argmin(abs(lf - l))]]
                vals.append(np.nan if lfu > cx.max() + 1e-9 else float(np.interp(lfu, cx, cy)))
            t = 0.0 if hi == lo else (x - lo) / (hi - lo)
            Z[k, j] = (1 - t) * vals[0] + t * vals[1]
    return Z


def panel(ax, data, kernel, thr, norm, cm):
    logf = np.linspace(-3, 0, 301)
    logsig = np.linspace(2, 8, 301)
    Z = surface(data, logf, logsig)
    cs = ax.contourf(logf, logsig, Z, levels=LEVELS, cmap=cm, norm=norm, extend="both")
    for level, ls, lw in ((0, "--", 0.8), (2.71, "-", 1.4), (3.84, "-", 0.9)):
        lines = ax.contour(logf, logsig, Z, levels=[level], colors=INK, linewidths=lw, linestyles=ls)
        ax.clabel(lines, fmt={level: f"{level:g}"}, fontsize=7, inline=True)
    f, u = data[:, 0], data[:, 1]
    fitted = u > U_ZERO_PLOT
    ax.plot(np.log10(f[fitted]), np.log10(sigma_surface_over_q(u[fitted])), "o", ms=2.2, color=INK, alpha=0.5,
            zorder=3, label="fitted grid points")
    mk = "v" if kernel == "isotropic" else "^"
    if kernel in thr:
        t = np.array(sorted(thr[kernel]))
        ax.plot(np.log10(t[:, 0]), np.log10(sigma_surface_over_q(t[:, 1])), mk, ms=6, mfc="none", mec=INK, mew=1.1,
                zorder=4, label="Δχ² = 2.71 above that fraction's own minimum")
    # Walker & Wardle (1998): M <~ 1e-3 Msun, R a few AU
    sig_lo, sig_hi = 1e-3 * MSUN / (np.pi * (3 * AU) ** 2), 1e-3 * MSUN / (np.pi * (1 * AU) ** 2)
    ax.axhspan(np.log10(sig_lo), np.log10(sig_hi), color="#4a3aa7", alpha=0.12, lw=0, zorder=1)
    ax.plot(0.0, np.log10(1e-3 * MSUN / (np.pi * (2 * AU) ** 2)), "*", ms=15, color="#4a3aa7", mec="white",
            mew=0.8, zorder=5, clip_on=False, label="Walker & Wardle (1998) clouds (10⁻³ M☉, R = 1–3 AU)")
    # present-day reference masses (fractions of the cosmic baryons; approximate)
    for share, text in ((0.06, "≈ mass in all stars"), (0.015, "≈ mass in all cold gas")):
        x = np.log10(share * OMEGA_B / OMEGA_DARK)
        ax.axvline(x, color=MUTED, lw=0.8, ls=":")
        ax.text(x, 7.9, f"{text} ", rotation=90, fontsize=7, color=MUTED, va="top", ha="right")
    ax.set_xlim(-3.05, 0.02)
    ax.set_ylim(2, 8)
    xt = [-3, -2, -1, 0]
    ax.set_xticks(xt, [f"{10**x:g}" for x in xt])
    ax.set_yticks(range(2, 9), [f"10$^{{{y}}}$" for y in range(2, 9)])
    ax.set_xlabel("clump fraction of the dark matter, f_cl", fontsize=9, color=INK)
    ax.set_title(f"{'Thomson' if kernel == 'thomson' else 'uniform / isotropic'} kernel", fontsize=10, color=INK)
    top = ax.secondary_xaxis("top", functions=(lambda x: x, lambda x: x))
    wt = [1e-4, 1e-3, 1e-2, 0.1]
    top.set_xticks([np.log10(w / OMEGA_DARK) for w in wt], [f"{w:g}" for w in wt])
    top.set_xlabel("ω_cl = f_cl ω_dark   (ordinary baryons: ω_b ≈ 0.022)", fontsize=9, color=INK)
    for a in (ax, top):
        a.tick_params(colors=MUTED, labelsize=8)
    return cs


def main():
    data = load()
    thr = thresholds()
    cm = cmap()
    norm = BoundaryNorm(LEVELS, cm.N, extend="both")
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.4), constrained_layout=True, sharey=True)
    fig.set_facecolor(SURFACE)
    for ax, kernel in zip(axes, ("isotropic", "thomson")):
        cs = panel(ax, data[kernel], kernel, thr, norm, cm)
    axes[0].set_ylabel("clump surface density Σ/Q  [g cm⁻²]", fontsize=9, color=INK)
    axes[0].legend(fontsize=7.5, frameon=False, loc="lower left")
    cb = fig.colorbar(cs, ax=axes, shrink=0.85, pad=0.01)
    cb.set_label("Δχ² relative to ΛCDM (profile over all other parameters)", fontsize=9, color=INK)
    cb.ax.tick_params(labelsize=8, colors=MUTED)
    fig.suptitle("PR4 CamSpec + 2018 low-ℓ: fit of a CDM + compact-clump mixture. "
                 "Blue fits better than ΛCDM, orange worse; contours are not confidence regions. Blank: not computed.",
                 fontsize=9, color=INK)
    fig.savefig(FIG / f"mix_contour_{LIKE}.png", dpi=170, facecolor=SURFACE)
    fig.savefig(FIG / f"mix_contour_{LIKE}.pdf", facecolor=SURFACE)


if __name__ == "__main__":
    main()
