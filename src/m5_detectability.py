"""M5: approximate detectability of the clump-photon coupling with realistic noise.

run      CLASS models to l = 3000:
         (a) the M3 coupling grid at fixed cosmology (fixed h, f_cl = 1);
         (b) Fisher derivative runs about u = 0 for {u, omega_b, omega_dm,
             100*theta_s, n_s, ln10^10 A_s, tau_reio}, two step sizes each.
analyse  Delta chi^2 against u for each noise configuration (src/m5_noise.py),
         and Fisher forecasts (conditional and marginalised sigma(u), 95% limits,
         degeneracies), compared with the published Planck bound u < 2.25e-4.

Usage: uv run python src/m5_detectability.py run | analyse
"""

import copy
import json
import sys

import matplotlib.pyplot as plt
import numpy as np

from cloud_mapping import sigma_over_m, sigma_surface_over_q
from cosmology import ROOT, class_params, fix_h, load_config
from m3_coupling_scan import LOG10_U, name_for
from m5_noise import configurations, delta_chi2, fisher
from run_class import load
from scan import run_grid

OUT = ROOT / "results" / "scans" / "m5"
ELL_MAX = 3000
PUBLISHED_U_MAX = 2.25e-4           # Stadler & Boehm 2018, 95%
U_STEPS = [3e-6, 1e-5, 3e-5]        # one-sided (u >= 0); response is linear in u here
PARAMS = {                          # parameter: (config key, step, fiducial taken from config)
    "omega_b": ("omega_b", 2e-4),
    "omega_dm": ("omega_dm", 2e-3),
    "100*theta_s": ("100*theta_s", 4e-4),
    "n_s": ("n_s", 4e-3),
    "ln10^{10}A_s": ("ln10^{10}A_s", 2e-2),
    "tau_reio": ("tau_reio", 5e-3),
}
FISHER_NAMES = ["u"] + list(PARAMS)
FIDUCIAL_U_STEP = 1e-5
GRID_LOG10_U = [lu for lu in LOG10_U if lu <= -1.0]


def m5_config():
    config = load_config()
    config["ell_max"] = ELL_MAX
    config["output"] = dict(config["output"], l_max_scalars=ELL_MAX + 500)
    return config


def fisher_params(config, u=0.0, shift=None):
    c = copy.deepcopy(config)
    if shift:
        key, delta = shift
        c["cosmology"][key] = c["cosmology"][key] + delta
    return class_params(c, f_cl=1.0, u=u)


def run():
    config = m5_config()
    h = load(config["name"])[1]["derived"]["h"]
    grid = [("lcdm", fix_h(class_params(config), h))]
    grid += [(name_for(lu), fix_h(class_params(config, f_cl=1.0, u=10**lu), h)) for lu in GRID_LOG10_U]
    status = run_grid(grid, config, OUT / "grid", with_rates=False)
    jobs = [("fid", fisher_params(config))]
    jobs += [(f"u_{s:g}", fisher_params(config, u=s)) for s in U_STEPS]
    for p, (key, step) in PARAMS.items():
        for frac in (1.0, 0.5):
            for sign in (+1, -1):
                jobs.append((f"{p}_{'p' if sign > 0 else 'm'}_{frac:g}", fisher_params(config, shift=(key, sign * frac * step))))
    status.update(run_grid(jobs, config, OUT / "fisher", with_rates=False))
    failed = {k: v for k, v in status.items() if v}
    print(f"{len(status) - len(failed)}/{len(status)} ok", failed or "")


KEYS = ("tt", "te", "ee", "pp")


def spectra(out, lensed=True):
    if lensed:
        return {k: out[k] for k in KEYS}
    return {"tt": out["tt_unlensed"], "te": out["te_unlensed"], "ee": out["ee_unlensed"], "pp": out["pp"]}


def derivatives(lensed=True, frac=1.0, u_step=FIDUCIAL_U_STEP):
    d = OUT / "fisher"
    fid = spectra(load("fid", d)[0], lensed)
    ders = [{k: (spectra(load(f"u_{u_step:g}", d)[0], lensed)[k] - fid[k]) / u_step for k in KEYS}]
    for p, (_, step) in PARAMS.items():
        plus = spectra(load(f"{p}_p_{frac:g}", d)[0], lensed)
        minus = spectra(load(f"{p}_m_{frac:g}", d)[0], lensed)
        ders.append({k: (plus[k] - minus[k]) / (2 * frac * step) for k in KEYS})
    return fid, ders


def forecast(fid, ders, cfg, lensing=False):
    f = fisher(fid, ders, cfg, "joint", lensing)
    cov = np.linalg.inv(f)
    sig_cond = 1 / np.sqrt(f[0, 0])
    sig_marg = np.sqrt(cov[0, 0])
    corr = {FISHER_NAMES[j]: float(cov[0, j] / np.sqrt(cov[0, 0] * cov[j, j])) for j in range(1, len(FISHER_NAMES))}
    return {"sigma_u_conditional": float(sig_cond), "sigma_u_marginalised": float(sig_marg),
            "u95_one_sided_1.645sigma": float(1.645 * sig_marg), "u95_half_gaussian_1.96sigma": float(1.96 * sig_marg),
            "corr_u_with": corr,
            "sigma_other_marginalised": {FISHER_NAMES[j]: float(np.sqrt(cov[j, j])) for j in range(1, len(FISHER_NAMES))}}


def u_at(lus, chi2, level):
    chi2 = np.asarray(chi2)
    i = np.nonzero(chi2 >= level)[0]
    if len(i) == 0 or i[0] == 0:
        return None
    i = i[0]
    return float(10 ** np.interp(np.log10(level), np.log10([chi2[i - 1], chi2[i]]), [lus[i - 1], lus[i]]))


def analyse():
    cfgs = configurations(ELL_MAX)
    ref = load("lcdm", OUT / "grid")[0]
    models = {lu: load(name_for(lu), OUT / "grid")[0] for lu in GRID_LOG10_U}
    result = {"delta_chi2": {}, "thresholds": {}, "fisher": {}, "fisher_stability": {}}

    # ---- Delta chi^2 at fixed cosmology ----
    variants = [("joint", "joint", False, True), ("tt", "tt", False, True), ("te", "te", False, True),
                ("ee", "ee", False, True), ("joint_unlensed", "joint", False, False),
                ("joint+pp_cv_envelope", "joint", True, True)]
    for cname, cfg in cfgs.items():
        result["delta_chi2"][cname] = {}
        result["thresholds"][cname] = {}
        for vname, which, lens_pp, lensed in variants:
            if lens_pp and not cfg.lensing:
                continue
            r = spectra(ref, lensed)
            vals = [delta_chi2(spectra(models[lu], lensed), r, cfg, which, lens_pp) for lu in GRID_LOG10_U]
            result["delta_chi2"][cname][vname] = vals
            result["thresholds"][cname][vname] = {"u(dchi2=1)": u_at(GRID_LOG10_U, vals, 1.0),
                                                  "u(dchi2=4)": u_at(GRID_LOG10_U, vals, 4.0)}

    # ---- Fisher ----
    for lensed in (True, False):
        fid, ders = derivatives(lensed=lensed)
        for cname, cfg in cfgs.items():
            key = f"{cname}{'' if lensed else ' [unlensed spectra]'}"
            result["fisher"][key] = forecast(fid, ders, cfg)
            if lensed and cfg.lensing:
                result["fisher"][f"{cname} + phiphi CV envelope"] = forecast(fid, ders, cfg, lensing=True)
    # stability: half steps, and different u steps
    for label, kw in (("half parameter steps", {"frac": 0.5}), ("u step 3e-6", {"u_step": 3e-6}), ("u step 3e-5", {"u_step": 3e-5})):
        fid, ders = derivatives(**kw)
        result["fisher_stability"][label] = {c: forecast(fid, ders, cfgs[c])["sigma_u_marginalised"]
                                             for c in ("planck_like", "so_baseline")}
    # tau calibration: Planck's large-scale polarization is systematics-limited, so the
    # noise-only low-l EE of the Planck-like set-up gives sigma(tau) ~2x tighter than
    # Planck 2018 (0.0075). Inflate the low-l EE noise to match and see what happens to sigma(u).
    fid, ders = derivatives()
    result["tau_calibration"] = []
    for factor in (1, 10, 100, 1000, 1e4):
        cfg = copy.deepcopy(cfgs["planck_like"])
        cfg.patches[0].n_ee = cfg.patches[0].n_ee * factor + (factor - 1) * 1e-20
        fc = forecast(fid, ders, cfg)
        result["tau_calibration"].append({"low_l_EE_noise_factor": factor,
                                          "sigma_tau": fc["sigma_other_marginalised"]["tau_reio"],
                                          "sigma_u_marginalised": fc["sigma_u_marginalised"]})
    for key, fc in result["fisher"].items():
        u95 = fc["u95_one_sided_1.645sigma"]
        fc["sigma_over_M_95_cm2_g"] = sigma_over_m(u95)
        fc["Sigma_over_Q_min_95_g_cm2"] = sigma_surface_over_q(u95)
        fc["ratio_to_published_planck_bound"] = u95 / PUBLISHED_U_MAX

    with open(OUT / "summary.json", "w") as f:
        json.dump(result, f, indent=2)
    plots(result, cfgs)
    print(json.dumps(result["thresholds"], indent=1))
    for key, fc in result["fisher"].items():
        top = sorted(fc["corr_u_with"].items(), key=lambda kv: -abs(kv[1]))[:3]
        print(f"{key:48s} sigma_cond {fc['sigma_u_conditional']:.2e} sigma_marg {fc['sigma_u_marginalised']:.2e} "
              f"u95(1.645) {fc['u95_one_sided_1.645sigma']:.2e}  ratio/published {fc['ratio_to_published_planck_bound']:.2f}  "
              f"top corr {[(k, round(v, 2)) for k, v in top]}")
    print(json.dumps(result["fisher_stability"], indent=1))
    print(json.dumps(result["tau_calibration"], indent=1))


def plots(result, cfgs):
    u = 10 ** np.array(GRID_LOG10_U)
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    styles = {"cv_fullsky": "k", "planck_like": "C0", "so_baseline": "C1", "so_baseline_white": "C2"}
    for cname in cfgs:
        ax.loglog(u, np.maximum(result["delta_chi2"][cname]["joint"], 1e-8), color=styles[cname], marker=".", label=f"{cfgs[cname].name}: TT+TE+EE")
        if "joint+pp_cv_envelope" in result["delta_chi2"][cname]:
            ax.loglog(u, result["delta_chi2"][cname]["joint+pp_cv_envelope"], color=styles[cname], ls="--",
                      label=f"{cfgs[cname].name}: + φφ (CV-only envelope)")
    ax.axhline(1, color="0.5", lw=0.6)
    ax.axhline(4, color="0.5", lw=0.6, ls=":")
    ax.axvline(PUBLISHED_U_MAX, color="C3", lw=0.8, label="published Planck bound (Stadler & Bœhm)")
    ax.set(xlabel="u", ylabel="Δχ² vs u = 0 (fixed cosmology)", ylim=(1e-4, 1e7),
           title="M5: detectability of the clump–photon coupling (f_cl = 1)")
    top = ax.secondary_xaxis("top", functions=(lambda x: 268.0 / np.maximum(x, 1e-30), lambda x: 268.0 / np.maximum(x, 1e-30)))
    top.set_xlabel("Σ/Q [g/cm²]")
    ax.legend(fontsize=7)
    fig.savefig(ROOT / "figures" / "m5_delta_chi2.png", dpi=150)
    plt.close(fig)

    # Fisher degeneracy: 1-sigma ellipses of u against each parameter, Planck-like and SO
    fid, ders = derivatives()
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
    for cname, color in (("planck_like", "C0"), ("so_baseline", "C1")):
        cov = np.linalg.inv(fisher(fid, ders, cfgs[cname]))
        for ax, j in zip(axes.flat, range(1, len(FISHER_NAMES))):
            sub = cov[np.ix_([j, 0], [j, 0])]
            vals, vecs = np.linalg.eigh(sub)
            t = np.linspace(0, 2 * np.pi, 200)
            ell = vecs @ (np.sqrt(np.maximum(vals, 0))[:, None] * np.array([np.cos(t), np.sin(t)]))
            ax.plot(ell[0], ell[1], color=color, label=cfgs[cname].name)
            ax.set(xlabel=f"Δ {FISHER_NAMES[j]}", ylabel="u")
            ax.ticklabel_format(axis="both", style="sci", scilimits=(-2, 2))
            ax.locator_params(axis="x", nbins=5)
            ax.axhline(0, color="0.7", lw=0.5); ax.axvline(0, color="0.7", lw=0.5)
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("M5: 1σ Fisher ellipses of u against the standard parameters (others marginalised)")
    fig.savefig(ROOT / "figures" / "m5_fisher_ellipses.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    {"run": run, "analyse": analyse}[sys.argv[1]]()
