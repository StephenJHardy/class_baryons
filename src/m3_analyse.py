"""Analyse the M3 coupling scan: M4 mapping, diagnostics and figures.

Reads results/scans/m3/{fixed_h,fixed_theta}/ and writes
results/scans/m3/summary.json, results/scans/m3/summary.md and figures/m3_*.png.
"""

import json

import matplotlib.pyplot as plt
import numpy as np

from cloud_mapping import sigma_over_m, sigma_surface_over_q
from cosmology import ROOT, class_params, fix_h, load_config
from m3_coupling_scan import OUT, name_for
from rates import clump_optical_depth, crossing_redshift
from run_class import load
from small_scale import suppression

ELL_MAX = 2500
SHOW_U = [-5.0, -4.5, -4.0, -3.5, -3.0, -2.0, -1.0]
N_PEAKS = 7


def dl(ell, cl):
    return ell * (ell + 1) * cl / (2 * np.pi)


def cv_chi2(model, ref, key):
    """Ideal full-sky cosmic-variance Delta chi^2 for one spectrum, 2 <= l <= ELL_MAX."""
    ell = ref["ell"][2:]
    d = model[key][2:] - ref[key][2:]
    if key.startswith("te"):
        tt, ee = key.replace("te", "tt"), key.replace("te", "ee")
        var = (ref[key][2:] ** 2 + ref[tt][2:] * ref[ee][2:]) / (2 * ell + 1)
    else:
        var = 2 * ref[key][2:] ** 2 / (2 * ell + 1)
    return float(np.sum(d**2 / var))


def frac(model, ref, key, ell_min=2):
    sel = ref["ell"] >= ell_min
    if key.startswith("te"):
        tt, ee = key.replace("te", "tt"), key.replace("te", "ee")
        return (model[key][sel] - ref[key][sel]) / np.sqrt(ref[tt][sel] * ref[ee][sel])
    return (model[key][sel] - ref[key][sel]) / ref[key][sel]


def onset_ell(res, ell, level):
    idx = np.nonzero(np.abs(res) >= level)[0]
    return int(ell[idx[0]]) if len(idx) else None


def find_peaks(ell, d, n):
    """Parabola-refined positions/heights of the first n local maxima of D_l above l = 100.

    Peaks are detected independently in each model and matched by order; a
    window around the LCDM peaks fails at strong coupling (shifts > 100).
    """
    cand = [i for i in range(1, len(d) - 1) if d[i] > d[i - 1] and d[i] >= d[i + 1] and ell[i] > 100]
    cand = cand[:n]
    out = []
    for i in cand:
        y0, y1, y2 = d[i - 1], d[i], d[i + 1]
        shift = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
        out.append((ell[i] + shift, y1 - 0.25 * (y0 - y2) * shift))
    # Strongly damped spectra can have fewer than n peaks below l_max: pad with NaN.
    out += [(np.nan, np.nan)] * (n - len(out))
    return np.array(out)


def loaded_sound_horizon_ratio(config, h, z_rec):
    """r_s(z_rec) with idm loading the photon fluid, over the standard r_s."""
    from run_class import Class   # respects CLASS_BUILD; never import classy directly
    c = Class()
    c.set({k: v for k, v in fix_h(class_params(config, f_cl=1.0, u=0.0), h).items()
           if k not in ("output", "lensing", "l_max_scalars", "P_k_max_h/Mpc", "z_pk")})
    c.compute()
    bg = c.get_background()
    c.struct_cleanup(); c.empty()
    sel = bg["z"] >= z_rec
    tau = bg["conf. time [Mpc]"][sel]
    r_std = 3 * bg["(.)rho_b"][sel] / (4 * bg["(.)rho_g"][sel])
    r_load = 3 * (bg["(.)rho_b"][sel] + bg["(.)rho_idm"][sel]) / (4 * bg["(.)rho_g"][sel])
    rs = lambda r: np.trapezoid(1 / np.sqrt(3 * (1 + r)), tau)
    return rs(r_load) / rs(r_std)


def threshold_u(log10_u, values, level):
    """Smallest u (log-interpolated) at which values first reach level."""
    values = np.asarray(values)
    above = np.nonzero(values >= level)[0]
    if len(above) == 0 or above[0] == 0:
        return None
    i = above[0]
    lv = np.log10([values[i - 1], values[i]])
    return float(10 ** np.interp(np.log10(level), lv, [log10_u[i - 1], log10_u[i]]))


def main():
    config = load_config()
    status = json.load(open(OUT / "status.json"))
    ok = lambda label: [lu for lu in status["log10_u"] if status["status"][label][name_for(lu)] is None]
    lus = ok("fixed_h")
    ref, ref_meta = load("lcdm", OUT / "fixed_h")
    ref_theta, ref_theta_meta = load("lcdm", OUT / "fixed_theta")
    ell = ref["ell"][2:]
    ref_peaks = find_peaks(ref["ell"], dl(ref["ell"], ref["tt"]), N_PEAKS)
    ref_peaks_unl = find_peaks(ref["ell"], dl(ref["ell"], ref["tt_unlensed"]), N_PEAKS)

    rows, models = [], {}
    for lu in lus:
        m, meta = load(name_for(lu), OUT / "fixed_h")
        models[lu] = m
        u = 10**lu
        row = {"log10_u": lu, "u": u, "sigma_over_M_cm2_g": sigma_over_m(u),
               "Sigma_over_Q_g_cm2": sigma_surface_over_q(u)}
        for key in ("tt", "te", "ee", "pp", "tt_unlensed", "te_unlensed", "ee_unlensed"):
            res = frac(m, ref, key)
            row[f"max_{key}"] = float(np.max(np.abs(res[ell <= ELL_MAX])))
        for key in ("tt", "ee"):
            res = frac(m, ref, key)
            row[f"onset_ell_{key}_0.1pct"] = onset_ell(res, ell, 1e-3)
            row[f"onset_ell_{key}_1pct"] = onset_ell(res, ell, 1e-2)
            row[f"mean_{key}_l30_800"] = float(np.mean(res[(ell >= 30) & (ell <= 800)]))
            row[f"mean_{key}_l1500_2500"] = float(np.mean(res[(ell >= 1500) & (ell <= 2500)]))
        for key in ("tt", "te", "ee", "tt_unlensed", "ee_unlensed"):
            row[f"chi2cv_{key}"] = cv_chi2(m, ref, key)
        peaks = find_peaks(m["ell"], dl(m["ell"], m["tt"]), N_PEAKS)
        peaks_unl = find_peaks(m["ell"], dl(m["ell"], m["tt_unlensed"]), N_PEAKS)
        row["peak_dl"] = (peaks[:, 0] - ref_peaks[:, 0]).tolist()
        row["peak_height_ratio"] = (peaks[:, 1] / ref_peaks[:, 1]).tolist()
        row["peak_dl_unlensed"] = (peaks_unl[:, 0] - ref_peaks_unl[:, 0]).tolist()
        row["peak_height_ratio_unlensed"] = (peaks_unl[:, 1] / ref_peaks_unl[:, 1]).tolist()
        row["z_dec_drag"] = crossing_redshift(m["z_th"], m["gamma_cl_over_H"])
        row["z_dec_photon"] = crossing_redshift(m["z_th"], m["gamma_ph_over_H"])
        zrec = meta["derived"]["z_rec"]
        row["z_rec"] = zrec
        row["z_star"] = meta["derived"]["z_star"]
        # CLASS's thermodynamics table is ordered by increasing z.
        row["gamma_ph_over_H_at_zrec"] = float(np.interp(zrec, m["z_th"], m["gamma_ph_over_H"]))
        row["gamma_cl_over_H_at_zrec"] = float(np.interp(zrec, m["z_th"], m["gamma_cl_over_H"]))
        before = m["z_th"] > zrec
        row["clump_over_thomson_opacity_before_zrec"] = float(np.max(m["dmu_idm_g"][before] / m["dkappa_thomson"][before]))
        row["tau_clump_after_zrec"] = clump_optical_depth(m["z_th"], m["dmu_idm_g"], m["conformal_time"], zrec)
        row.update(suppression(m["k_h"], m["pk"], ref["pk"], meta["derived"]["Omega_m"]))
        row["sigma8"] = meta["derived"]["sigma8"]
        if lu in ok("fixed_theta"):
            _, tmeta = load(name_for(lu), OUT / "fixed_theta")
            row["h_at_fixed_theta"] = tmeta["derived"]["h"]
        rows.append(row)

    lu_arr = np.array(lus)
    thresholds = {}
    for key in ("tt", "te", "ee", "tt_unlensed", "ee_unlensed"):
        vals = [r[f"chi2cv_{key}"] for r in rows]
        thresholds[f"u_chi2cv_{key}=1"] = threshold_u(lu_arr, vals, 1.0)
    for key in ("tt", "ee", "pp", "tt_unlensed"):
        thresholds[f"u_max_{key}=0.1pct"] = threshold_u(lu_arr, [r[f"max_{key}"] for r in rows], 1e-3)
        thresholds[f"u_max_{key}=1pct"] = threshold_u(lu_arr, [r[f"max_{key}"] for r in rows], 1e-2)
    rs_ratio = loaded_sound_horizon_ratio(config, status["h_fixed"], ref_meta["derived"]["z_rec"])
    strongest = rows[-1]
    loaded = {
        "rs_loaded_over_rs_standard": rs_ratio,
        "expected_fractional_peak_shift_if_fully_loaded": 1 / rs_ratio - 1,
        "u_strongest": strongest["u"],
        "measured_mean_fractional_peak_shift_unlensed": float(np.nanmean(
            np.array(strongest["peak_dl_unlensed"]) / ref_peaks_unl[:, 0])),
        "measured_fractional_peak_shift_unlensed_per_peak": (
            np.array(strongest["peak_dl_unlensed"]) / ref_peaks_unl[:, 0]).tolist(),
    }
    summary = {"reference": "lcdm at fixed h (bit-identical to baseline_lcdm)",
               "reference_peaks_lensed": ref_peaks.tolist(),
               "reference_peaks_unlensed": ref_peaks_unl.tolist(),
               "clas_limit": {"max_log10_u_ok": max(lus),
                              "failures": {k: v for k, v in status["status"]["fixed_h"].items() if v}},
               "h_baseline": status["h_fixed"],
               "thresholds": thresholds, "strong_coupling": loaded, "rows": rows}
    with open(OUT / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=float)
    write_markdown(rows, thresholds, loaded)
    plots(models, ref, rows, ref_peaks)
    print(json.dumps({"thresholds": thresholds, "strong_coupling": loaded}, indent=2))


def fmt(x, spec=".2g"):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "—"
    return format(x, spec)


def write_markdown(rows, thresholds, loaded):
    lines = ["| log₁₀u | σ/M [cm²/g] | Σ/Q [g/cm²] | z_dec (drag) | Γ_γ/H at z_rec | clump/Thomson opacity (z>z_rec) | τ_clump (z<z_rec) | max ΔTT | max ΔTT unlensed | max ΔEE | max Δφφ | Δχ²_CV TT | Δχ²_CV EE | k_½ [h/Mpc] | M_hm [M☉/h] | h (fixed θ_s) |",
             "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in rows:
        lines.append("| " + " | ".join([
            f"{r['log10_u']:.3g}", fmt(r["sigma_over_M_cm2_g"]), fmt(r["Sigma_over_Q_g_cm2"]),
            fmt(r["z_dec_drag"], ".3g"), fmt(r["gamma_ph_over_H_at_zrec"]),
            fmt(r["clump_over_thomson_opacity_before_zrec"]), fmt(r["tau_clump_after_zrec"]),
            fmt(r["max_tt"]), fmt(r["max_tt_unlensed"]), fmt(r["max_ee"]), fmt(r["max_pp"]),
            fmt(r["chi2cv_tt"]), fmt(r["chi2cv_ee"]), fmt(r["k_half"], ".3g"), fmt(r["M_hm"]),
            fmt(r.get("h_at_fixed_theta"), ".4f")]) + " |")
    lines += ["", "Thresholds (log-interpolated u at which the quantity first reaches the level):", ""]
    lines += [f"- `{k}`: {fmt(v, '.3g')}" for k, v in thresholds.items()]
    lines += ["", f"Strong coupling: fully photon-loaded r_s ratio {loaded['rs_loaded_over_rs_standard']:.4f} "
              f"→ expected peak shift {100 * loaded['expected_fractional_peak_shift_if_fully_loaded']:.1f}%; "
              f"measured at u = {loaded['u_strongest']:.3g}: {100 * loaded['measured_mean_fractional_peak_shift_unlensed']:.2f}% (unlensed TT)."]
    (OUT / "summary.md").write_text("\n".join(lines) + "\n")


def plots(models, ref, rows, ref_peaks):
    figs = ROOT / "figures"
    show = [lu for lu in SHOW_U if lu in models]
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(show)))
    ell = ref["ell"]
    sel = (ell >= 2) & (ell <= ELL_MAX)

    for key, title, unl in (("tt", "TT", "tt_unlensed"), ("ee", "EE", "ee_unlensed"),
                            ("te", "TE", "te_unlensed"), ("pp", r"$\phi\phi$", None)):
        fig, (a1, a2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True, constrained_layout=True,
                                     gridspec_kw={"height_ratios": [2, 1.4]})
        scale = (ell * (ell + 1)) ** 2 / (2 * np.pi) if key == "pp" else dl(ell, 1.0)
        a1.plot(ell[sel], (ref[key] * scale)[sel], "k", lw=1.2, label="ΛCDM")
        for lu, c in zip(show, colors):
            m = models[lu]
            a1.plot(ell[sel], (m[key] * scale)[sel], color=c, lw=0.8, label=f"log₁₀u = {lu:g}")
            a2.plot(ell[sel], 100 * frac(m, ref, key)[ell[2:] <= ELL_MAX], color=c, lw=0.8)
            if unl:
                a2.plot(ell[sel], 100 * frac(m, ref, unl)[ell[2:] <= ELL_MAX], color=c, lw=0.8, ls="--")
        for lp in ref_peaks[:, 0]:
            for ax in (a1, a2):
                ax.axvline(lp, color="0.8", lw=0.6, zorder=0)
        a1.set_ylabel(r"$[L(L+1)]^2 C_L^{\phi\phi}/2\pi$" if key == "pp" else rf"$D_\ell^{{{title}}}$")
        a1.set_title(f"M3: {title} vs coupling (f_cl = 1, fixed H0); grey lines = ΛCDM TT peaks")
        a1.legend(fontsize=7, ncol=2)
        a2.set_ylabel("% residual" + ("\n(TE / √(TT EE))" if key == "te" else "") + ("\nsolid lensed, dashed unlensed" if unl else ""))
        a2.set_xlabel(r"$\ell$")
        a2.set_yscale("symlog", linthresh=0.1)
        a2.axhline(0, color="k", lw=0.5)
        if key == "pp":
            a1.set_xscale("log")
        fig.savefig(figs / f"m3_{key}.png", dpi=150)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    for lu, c in zip(show, colors):
        ax.semilogx(ref["k_h"], models[lu]["pk"] / ref["pk"], color=c, label=f"log₁₀u = {lu:g}")
    for level in (0.9, 0.5, 0.25):
        ax.axhline(level, color="0.7", lw=0.6, ls=":")
    ax.set(xlabel="k [h/Mpc]", ylabel=r"$T^2 = P/P_{\Lambda CDM}$", ylim=(0, 1.1),
           title="M3: matter power suppression (z = 0)")
    ax.legend(fontsize=7)
    fig.savefig(figs / "m3_pk.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    for lu, c in zip(show, colors):
        m = models[lu]
        ax.loglog(1 + m["z_th"], m["gamma_cl_over_H"], color=c, label=f"log₁₀u = {lu:g}")
        ax.loglog(1 + m["z_th"], m["gamma_ph_over_H"], color=c, ls="--")
    ax.axhline(1, color="k", lw=0.6)
    ax.axvline(1 + rows[0]["z_rec"], color="0.5", lw=0.6, ls=":")
    ax.set(xlabel="1 + z", ylabel="Γ / H", ylim=(1e-8, 1e8), xlim=(1, 1e9),
           title="M3: clump drag Γ_cl→γ/H (solid), photon opacity Γ_γ→cl/H (dashed); dotted: z_rec")
    ax.title.set_fontsize(9)
    ax.legend(fontsize=7)
    fig.savefig(figs / "m3_rates.png", dpi=150)
    plt.close(fig)

    lu = np.array([r["log10_u"] for r in rows])
    u = 10**lu
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    ax = axes[0, 0]
    for key, ls in (("tt", "-"), ("ee", "-"), ("te", "-"), ("pp", "-"), ("tt_unlensed", "--"), ("ee_unlensed", "--")):
        ax.loglog(u, [r[f"max_{key}"] for r in rows], ls, marker=".", label=key)
    ax.axhline(1e-3, color="0.7", lw=0.6, ls=":")
    ax.set(ylabel="max |ΔC_ℓ/C_ℓ|, 2 ≤ ℓ ≤ 2500", title="Size of the effect")
    ax.legend(fontsize=7)
    ax = axes[0, 1]
    for key, ls in (("tt", "-"), ("ee", "-"), ("te", "-"), ("tt_unlensed", "--"), ("ee_unlensed", "--")):
        ax.loglog(u, [max(r[f"chi2cv_{key}"], 1e-6) for r in rows], ls, marker=".", label=key)
    ax.axhline(1, color="k", lw=0.6)
    ax.set(ylabel="Δχ² (ideal full-sky cosmic variance)", title="Preliminary detectability (M5 refines)")
    ax.legend(fontsize=7)
    ax = axes[1, 0]
    ax.loglog(u, [r["z_dec_drag"] for r in rows], marker=".", label="z_dec: clump drag Γ_cl/H = 1")
    ax.loglog(u, [r["z_dec_photon"] for r in rows], marker=".", label="z_dec: photon opacity Γ_γ/H = 1")
    ax.loglog(u, [r["z_rec"] for r in rows], marker=".", label="z_rec (visibility peak)")
    ax.set(xlabel="u", ylabel="redshift", title="Decoupling epochs")
    ax.legend(fontsize=7)
    ax = axes[1, 1]
    for key, lab in (("k_10", "k_10 (T²=0.9)"), ("k_half", "k_½ (T²=0.5)"), ("k_hm", "k_hm (T=0.5)")):
        ax.loglog(u, [r[key] for r in rows], marker=".", label=lab)
    ax.set(xlabel="u", ylabel="k [h/Mpc]", title="Small-scale suppression (k ≤ 10 h/Mpc)")
    ax.legend(fontsize=7)
    for ax in axes.flat:
        top = ax.secondary_xaxis("top", functions=(lambda x: 268.0 / np.maximum(x, 1e-30),
                                                   lambda x: 268.0 / np.maximum(x, 1e-30)))
        top.set_xlabel("Σ/Q [g/cm²]", fontsize=8)
        ax.set_xlim(u.min(), u.max())
    fig.savefig(figs / "m3_summary.png", dpi=150)
    plt.close(fig)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    pcol = plt.cm.plasma(np.linspace(0, 0.85, N_PEAKS))
    for p in range(N_PEAKS):
        a1.semilogx(u, [r["peak_dl_unlensed"][p] for r in rows], color=pcol[p], label=f"peak {p + 1}")
        a1.semilogx(u, [r["peak_dl"][p] for r in rows], color=pcol[p], ls=":")
        a2.semilogx(u, [100 * (r["peak_height_ratio_unlensed"][p] - 1) for r in rows], color=pcol[p])
        a2.semilogx(u, [100 * (r["peak_height_ratio"][p] - 1) for r in rows], color=pcol[p], ls=":")
    a1.set(xlabel="u", ylabel="Δℓ of TT peak", title="Peak positions (solid unlensed, dotted lensed)")
    a1.set_yscale("symlog", linthresh=1)
    a2.set(xlabel="u", ylabel="% change in peak height", title="Peak heights (solid unlensed, dotted lensed)")
    a2.set_yscale("symlog", linthresh=0.1)
    a1.legend(fontsize=7)
    fig.savefig(figs / "m3_peaks.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
