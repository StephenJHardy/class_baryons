"""Analyse the mixed CDM + clump profiles written by src/profile_mix.py.

For every clump fraction f_cl and kernel:
  - the profile Delta chi^2(u) = 2 [min -log P (u) - min over the scanned u], without assuming
    that the minimum is at u = 0;
  - the same curve relative to a common LCDM reference (the f_cl = 0 fit);
  - nominal fixed-fraction profile thresholds: the first upward crossings of Delta chi^2 = 2.71
    and 3.84, by linear interpolation between bracketing points, with a numerical uncertainty
    from the per-point uncertainty and the local slope. No crossing is extrapolated.
The per-point numerical uncertainty is the regression error combined with the scatter of
repeated fits (tagged runs: other seeds and independent starts), pooled over all repeats.

Kernel comparison: the difference of the thresholds, and (more precisely) the direct
objective difference between kernels at identical parameters (src/mix_paired_kernel.py),
in which the CLASS numerical noise cancels.

Usage: uv run python src/profile_mix_analyse.py [likelihoods=camspec_npipe]
Output: results/mix/<likelihoods>/summary.json, summary.md and figures/mix_*.png
"""

import json
import sys

import matplotlib.pyplot as plt
import numpy as np

from cloud_mapping import sigma_surface_over_q
from cosmology import ROOT

OUT = ROOT / "results" / "mix"
FIG = ROOT / "figures"
LEVELS = (2.71, 3.84)
KERNEL_STYLE = {"thomson": dict(color="#2a78d6", marker="o", ls="-", label="Thomson"),
                "isotropic": dict(color="#eb6834", marker="s", ls="--", label="uniform / isotropic")}
INK, MUTED = "#1a1a19", "#6b6a63"
VALUE = "minuslogpost_fit"      # the regression estimate of the minimum; the direct value is tabulated too


def load(likelihoods):
    """{f: {"u0": result, kernel: {u: result}, "repeats": [(kernel, u, tag, result)]}}, and the LCDM fit."""
    data, lcdm = {}, None
    for path in sorted((OUT / likelihoods).glob("f*/**/result.json")):
        r = json.load(open(path))
        f = r["f_cl"]
        if f == 0:
            lcdm = r
            continue
        d = data.setdefault(f, {"u0": None, "thomson": {}, "isotropic": {}, "repeats": []})
        if r.get("tag"):
            d["repeats"].append((r["kernel"], r["u"], r["tag"], r))
        elif r["u"] == 0:
            d["u0"] = r
        else:
            d[r["kernel"]][r["u"]] = r
    if lcdm is not None:
        for path in (OUT / likelihoods / "f0").glob("**/result.json"):
            r = json.load(open(path))
            if r.get("tag"):
                data.setdefault(0.0, {"u0": lcdm, "thomson": {}, "isotropic": {}, "repeats": []})["repeats"].append(
                    (None, 0.0, r["tag"], r))
    return data, lcdm


def repeat_scatter(data):
    """Pooled rms scatter of the fitted minimum between repeated fits of the same point."""
    diffs, rows = [], []
    for f, d in data.items():
        groups = {}
        for kernel, u, tag, r in d["repeats"]:
            base = d["u0"] if u == 0 else d.get(kernel, {}).get(u)
            if base is None:
                continue
            groups.setdefault((kernel, u), [base]).append(r)
        for (kernel, u), runs in groups.items():
            vals = np.array([r[VALUE] for r in runs if r["converged"]])
            if len(vals) < 2:
                continue
            diffs.extend(vals - vals.mean())
            rows.append({"f_cl": f, "kernel": kernel, "u": u, "n": len(vals), "values": vals.tolist(),
                         "direct": [r["minuslogpost_direct"] for r in runs if r["converged"]],
                         "tags": [r.get("tag") or "main" for r in runs if r["converged"]],
                         "range": float(np.ptp(vals))})
    n_groups = len(rows)
    sigma = float(np.sqrt(np.sum(np.square(diffs)) / max(len(diffs) - n_groups, 1))) if diffs else None
    return sigma, rows


def crossing(u, d, err, level, i_min):
    """First upward crossing of `level` above the profile minimum; None if not bracketed."""
    for i in range(i_min, len(u) - 1):
        if d[i] < level <= d[i + 1]:
            slope = (d[i + 1] - d[i]) / (u[i + 1] - u[i])
            uc = u[i] + (level - d[i]) / slope
            sigma_d = float(np.hypot(max(err[i], err[i + 1]), err[i_min]))   # point and minimum both uncertain
            return {"u": float(uc), "u_err": float(sigma_d / slope), "bracket": [float(u[i]), float(u[i + 1])],
                    "Sigma_over_Q": sigma_surface_over_q(uc)}
    return None


def profile(d, kernel, lcdm_value, sigma_rep):
    pts = [d["u0"]] + [d[kernel][u] for u in sorted(d[kernel])] if d["u0"] else [d[kernel][u] for u in sorted(d[kernel])]
    ok = [r for r in pts if r["converged"]]
    bad = [r["u"] for r in pts if not r["converged"]]
    if not ok:
        return None
    u = np.array([r["u"] for r in ok])
    m = np.array([r[VALUE] for r in ok])
    direct = np.array([r["minuslogpost_direct"] for r in ok])
    err = np.array([np.hypot(r["minuslogpost_fit_err"], sigma_rep or 0.0) for r in ok])   # in -log P
    i_min = int(np.argmin(m))
    dchi2 = 2 * (m - m[i_min])
    out = {"u": u.tolist(), "minuslogpost_fit": m.tolist(), "minuslogpost_direct": direct.tolist(),
           "regression_err": [r["minuslogpost_fit_err"] for r in ok],
           "delta_chi2": dchi2.tolist(), "delta_chi2_err": (2 * np.sqrt(2) * err).tolist(),
           "delta_chi2_vs_lcdm": (2 * (m - lcdm_value)).tolist() if lcdm_value is not None else None,
           "u_at_minimum": float(u[i_min]),
           "minimum_below_u0": float(2 * (m[0] - m[i_min])) if u[0] == 0 else None,
           "not_converged_u": bad, "u_max_scanned": float(u.max()),
           "fit_noise": [r["fit_noise"] for r in ok],
           "active_bounds": {f"{r['u']:g}": r["active_bounds"] for r in ok if r["active_bounds"]},
           "best": {n: [r["best"][n] for r in ok] for n in ok[0]["names"]},
           "derived": {n: [r["direct_detail"]["derived"].get(n) for r in ok] for n in ("H0", "sigma8")},
           "chi2_parts": {n: [r["direct_detail"]["chi2"][n] for r in ok] for n in ok[0]["direct_detail"]["chi2"]}}
    out["crossings"] = {str(level): crossing(u, dchi2, 2 * err, level, i_min) for level in LEVELS}
    return out


def fmt_u(c):
    if c is None:
        return "not reached"
    digits = max(0, int(np.floor(np.log10(c["u"])) - np.floor(np.log10(max(c["u_err"], 1e-300)))))
    exp = int(np.floor(np.log10(c["u"])))
    return f"({c['u'] / 10**exp:.{min(digits, 2)}f} ± {c['u_err'] / 10**exp:.{min(digits, 2)}f})×10^{exp}"


def main(likelihoods="camspec_npipe"):
    data, lcdm = load(likelihoods)
    sigma_rep, repeats = repeat_scatter(data)
    lcdm_value = lcdm[VALUE] if lcdm else None
    fractions = sorted((f for f in data if f > 0), reverse=True)
    summary = {"likelihoods": likelihoods, "value_used": VALUE, "levels": LEVELS,
               "note": "nominal fixed-fraction profile thresholds; not joint (f_cl, u) confidence contours",
               "lcdm_reference": {"minuslogpost_fit": lcdm_value, "minuslogpost_direct": lcdm["minuslogpost_direct"],
                                  "best": lcdm["best"]} if lcdm else None,
               "repeat_scatter_minuslogpost": sigma_rep, "repeats": repeats, "fractions": {}}
    paired_path = OUT / likelihoods / "paired_kernel.json"
    paired = json.load(open(paired_path)) if paired_path.exists() else {}
    for f in fractions:
        entry = {k: profile(data[f], k, lcdm_value, sigma_rep) for k in KERNEL_STYLE if data[f][k]}
        if data[f]["u0"] is not None and lcdm_value is not None:
            entry["u0_minus_lcdm_chi2"] = 2 * (data[f]["u0"][VALUE] - lcdm_value)
        diff = {}
        for level in LEVELS:
            t = entry.get("thomson", {}) and entry["thomson"]["crossings"][str(level)]
            i = entry.get("isotropic", {}) and entry["isotropic"]["crossings"][str(level)]
            if t and i:
                frac = (i["u"] - t["u"]) / t["u"]
                diff[str(level)] = {"fractional_difference_iso_minus_thomson": frac,
                                    "err_independent": float(np.hypot(i["u_err"], t["u_err"]) / t["u"])}
        entry["kernel_threshold_difference"] = diff
        entry["paired_kernel"] = paired.get(f"{f:g}")
        summary["fractions"][f"{f:g}"] = entry
    (OUT / likelihoods).mkdir(parents=True, exist_ok=True)
    with open(OUT / likelihoods / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=1, default=float)

    # ---- table
    lines = ["| f_cl | kernel | u at min | u (Δχ²=2.71) | Σ/Q min [g/cm²] | u (Δχ²=3.84) | u scanned to | "
             "Δχ²(min) − ΛCDM |", "|---|---|---|---|---|---|---|---|"]
    for f in fractions:
        for k in KERNEL_STYLE:
            p = summary["fractions"][f"{f:g}"].get(k)
            if not p:
                continue
            c1, c2 = p["crossings"]["2.71"], p["crossings"]["3.84"]
            rel = min(p["delta_chi2_vs_lcdm"]) if p["delta_chi2_vs_lcdm"] else float("nan")
            lines.append(f"| {f:g} | {k} | {p['u_at_minimum']:g} | {fmt_u(c1)} | "
                         f"{c1['Sigma_over_Q']:.2g} | {fmt_u(c2)} | {p['u_max_scanned']:g} | {rel:+.2f} |"
                         if c1 else
                         f"| {f:g} | {k} | {p['u_at_minimum']:g} | not reached | – | {fmt_u(c2)} | "
                         f"{p['u_max_scanned']:g} | {rel:+.2f} |")
    table = "\n".join(lines)
    with open(OUT / likelihoods / "summary.md", "w") as fh:
        fh.write(table + "\n")
    print(table)
    print("repeat scatter in -log P:", sigma_rep)
    figures(summary, fractions, likelihoods)


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(colors=MUTED, labelsize=8)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.grid(True, color="#e6e5df", lw=0.5)
    ax.set_axisbelow(True)


def figures(summary, fractions, likelihoods):
    # 1. profiles: one panel per fraction, both kernels (small multiples; one y-scale meaning)
    n = len(fractions)
    fig, axes = plt.subplots(1, n, figsize=(3.1 * n, 3.4), constrained_layout=True, squeeze=False)
    for ax, f in zip(axes[0], fractions):
        style(ax)
        entry = summary["fractions"][f"{f:g}"]
        for k, st in KERNEL_STYLE.items():
            p = entry.get(k)
            if not p:
                continue
            ax.errorbar(np.array(p["u"]) * f * 1e4, p["delta_chi2"], yerr=p["delta_chi2_err"], color=st["color"],
                        marker=st["marker"], ls=st["ls"], lw=1.5, ms=4.5, capsize=2, label=st["label"])
        for level in LEVELS:
            ax.axhline(level, color=MUTED, lw=0.7, ls=":")
        ax.text(0.02, 2.71, " 2.71", color=MUTED, fontsize=7, va="bottom", transform=ax.get_yaxis_transform())
        ax.text(0.02, 3.84, " 3.84", color=MUTED, fontsize=7, va="bottom", transform=ax.get_yaxis_transform())
        ax.set_title(f"f_cl = {f:g}", fontsize=10, color=INK)
        ax.set_xlabel("f_cl · u  [10⁻⁴]", fontsize=9, color=INK)
        ax.set_ylim(-0.6, 9)
    axes[0][0].set_ylabel("Δχ² (profile, relative to its minimum)", fontsize=9, color=INK)
    axes[0][0].legend(fontsize=8, frameon=False, loc="upper left")
    fig.suptitle("Fixed-fraction profiles in the clump coupling, PR4 CamSpec + 2018 low-ℓ", fontsize=10, color=INK)
    fig.savefig(FIG / f"mix_profiles_{likelihoods}.png", dpi=160)
    fig.savefig(FIG / f"mix_profiles_{likelihoods}.pdf")

    # 2. boundary in (f_cl, u): nominal fixed-fraction thresholds
    fig, ax = plt.subplots(figsize=(6.2, 4.6), constrained_layout=True)
    style(ax)
    for k, st in KERNEL_STYLE.items():
        for level, alpha, lw in ((2.71, 1.0, 1.8), (3.84, 0.55, 1.2)):
            pts = [(f, summary["fractions"][f"{f:g}"][k]["crossings"][str(level)]) for f in fractions
                   if summary["fractions"][f"{f:g}"].get(k)]
            pts = [(f, c) for f, c in pts if c]
            if not pts:
                continue
            ff = np.array([f for f, _ in pts])
            uu = np.array([c["u"] for _, c in pts])
            ee = np.array([c["u_err"] for _, c in pts])
            ax.errorbar(ff, uu, yerr=ee, color=st["color"], marker=st["marker"], ls=st["ls"], lw=lw, ms=5,
                        alpha=alpha, capsize=2, label=f"{st['label']}, Δχ² = {level}")
    ref = summary["fractions"].get("1", {}).get("thomson")
    if ref and ref["crossings"]["2.71"]:
        ff = np.array([min(fractions), 1.0])
        ax.plot(ff, ref["crossings"]["2.71"]["u"] / ff, color=MUTED, lw=0.8, ls=":", label="u ∝ 1/f_cl (reference)")
    ax.set(xscale="log", yscale="log")
    ax.set_xlabel("clump fraction of the dark component, f_cl", fontsize=9, color=INK)
    ax.set_ylabel("coupling u_cl at the threshold", fontsize=9, color=INK)
    ax.set_title("Nominal fixed-fraction profile thresholds (not a joint confidence contour)\n"
                 "couplings above each curve exceed the stated Δχ² at that fraction", fontsize=9, color=INK)
    ax.legend(fontsize=8, frameon=False)
    sec = ax.secondary_yaxis("right", functions=(lambda u: sigma_surface_over_q(np.maximum(u, 1e-12)),
                                                 lambda s: 268.0 / np.maximum(s, 1e-12) * 1.0))
    sec.set_ylabel("Σ/Q  [g cm⁻²]", fontsize=9, color=INK)
    sec.tick_params(colors=MUTED, labelsize=8)
    fig.savefig(FIG / f"mix_boundary_{likelihoods}.png", dpi=160)
    fig.savefig(FIG / f"mix_boundary_{likelihoods}.pdf")

    # 3. fitted parameters along the Thomson profiles (small multiples, one parameter per panel)
    shown = [("sigma8", "σ₈", "derived"), ("theta_s_100", "100 θ_s", "best"), ("omega_dm", "ω_dark", "best"),
             ("omega_b", "ω_b", "best"), ("H0", "H₀", "derived"), ("n_s", "n_s", "best")]
    fig, axes = plt.subplots(2, 3, figsize=(9.5, 5.6), constrained_layout=True)
    greys = plt.get_cmap("Blues")(np.linspace(0.95, 0.35, len(fractions)))
    for ax, (name, label, where) in zip(axes.flat, shown):
        style(ax)
        for f, col in zip(fractions, greys):
            p = summary["fractions"][f"{f:g}"].get("thomson")
            if not p:
                continue
            ax.plot(np.array(p["u"]) * f * 1e4, p[where][name], color=col, marker="o", ms=3.5, lw=1.4,
                    label=f"f_cl = {f:g}")
        ax.set_xlabel("f_cl · u  [10⁻⁴]", fontsize=9, color=INK)
        ax.set_title(label, fontsize=10, color=INK)
    axes[0][0].legend(fontsize=7.5, frameon=False, title="Thomson kernel", title_fontsize=7.5)
    fig.suptitle("Best-fit parameters along the fixed-fraction profiles", fontsize=10, color=INK)
    fig.savefig(FIG / f"mix_parameters_{likelihoods}.png", dpi=160)
    fig.savefig(FIG / f"mix_parameters_{likelihoods}.pdf")


if __name__ == "__main__":
    main(*sys.argv[1:2])
