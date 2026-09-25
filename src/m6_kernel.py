"""M6 (Experiment B, minimal): Thomson vs perfect-absorber clump kernel.

The patched CLASS (class_patches/idm_g_kernel.patch) adds two coefficients to
the clump part of the photon collision term:
  idm_g_quadrupole_coefficient  c_T  (1 = Thomson, 0 = isotropic absorber)
  idm_g_polarization_coefficient c_E  (1 = Thomson, 0 = isotropic absorber)

Usage (the CLASS build is chosen per process by CLASS_BUILD):
  CLASS_BUILD=stock   uv run python src/m6_kernel.py run
  CLASS_BUILD=patched uv run python src/m6_kernel.py run
  uv run python src/m6_kernel.py analyse
"""

import json
import sys

import matplotlib.pyplot as plt
import numpy as np

from cosmology import ROOT, class_params, fix_h, load_config
from diagnostics import max_abs, residuals
from m3_analyse import cv_chi2, frac
from run_class import CLASS_BUILD, load
from scan import run_grid

OUT = ROOT / "results" / "scans" / "m6"
TARGET_U = [1e-5, 1e-4, 2.25e-4, 1e-3]
EXTRA_U = [1e-6, 1e-2, 1e-1]
FINE_LOG10_U = np.round(np.arange(-5.5, -2.95, 0.1), 3)
KERNELS = {
    "thomson_default": None,                      # patched build, coefficients not passed
    "thomson_explicit": (1.0, 1.0),
    "absorber": (0.0, 0.0),
    "no_quadrupole_regen": (0.0, 1.0),            # c_T = 0 only
    "no_polarization_source": (1.0, 0.0),         # c_E = 0 only
}
TCA_VARIANTS = {   # absorber at u = 2.25e-4 with different tight-coupling handling
    "absorber_tca_first_order": {"tight_coupling_approximation": 2},
    "absorber_tca_early_off": {"tight_coupling_trigger_tau_c_over_tau_h": 0.0015,
                                "tight_coupling_trigger_tau_c_over_tau_k": 0.001},
}
KEYS = ("tt", "te", "ee", "pp", "tt_unlensed", "te_unlensed", "ee_unlensed")


def name(kernel, u):
    return f"{kernel}_u{np.log10(u):+.3f}"


def model(config, h, u, coeffs=None, extra=None):
    extra = dict(extra or {})
    if coeffs is not None:
        extra["idm_g_quadrupole_coefficient"], extra["idm_g_polarization_coefficient"] = coeffs
    return fix_h(class_params(config, f_cl=1.0, u=u, extra=extra), h)


def run():
    config = load_config()
    h = load(config["name"])[1]["derived"]["h"]
    us = TARGET_U + EXTRA_U
    jobs = [("lcdm", fix_h(class_params(config), h))]
    if CLASS_BUILD == "stock":
        jobs += [(name("stock", u), model(config, h, u)) for u in us]
        jobs += [(name("stock", 10**lu), model(config, h, 10**lu)) for lu in FINE_LOG10_U]
    else:
        for kernel, coeffs in KERNELS.items():
            jobs += [(name(kernel, u), model(config, h, u, coeffs)) for u in us]
        for kernel in ("thomson_default", "absorber"):
            jobs += [(name(kernel, 10**lu), model(config, h, 10**lu, KERNELS[kernel])) for lu in FINE_LOG10_U]
        for variant, extra in TCA_VARIANTS.items():
            jobs.append((name(variant, 2.25e-4), model(config, h, 2.25e-4, KERNELS["absorber"], extra)))
            jobs.append((name(variant.replace("absorber", "thomson"), 2.25e-4),
                         model(config, h, 2.25e-4, KERNELS["thomson_explicit"], extra)))
    status = run_grid(jobs, config, OUT / CLASS_BUILD)
    failed = {k: v for k, v in status.items() if v}
    print(f"{CLASS_BUILD}: {len(jobs) - len(failed)}/{len(jobs)} ok", failed or "")
    with open(OUT / f"status_{CLASS_BUILD}.json", "w") as f:
        json.dump(status, f, indent=2)


def get(build, key):
    return load(key, OUT / build)[0]


def summarise(a, b, ref):
    """Max |a - b| in units of ref (TE normalised by sqrt(TT EE)) and ideal-CV Delta chi^2."""
    out = {}
    for key in KEYS:
        out[f"max_{key}"] = float(np.max(np.abs(frac(a, ref, key) - frac(b, ref, key))))
    for key in ("tt", "te", "ee", "ee_unlensed"):
        out[f"chi2cv_{key}"] = cv_chi2(a, b, key)
    ell = ref["ell"][2:]
    for key in ("ee", "ee_unlensed", "te"):
        d = np.abs(frac(a, ref, key) - frac(b, ref, key))
        out[f"max_{key}_l<30"] = float(d[ell < 30].max())
        out[f"max_{key}_l>=30"] = float(d[ell >= 30].max())
    out["max_pk_dT2"] = float(np.max(np.abs(a["pk"] - b["pk"]) / ref["pk"]))
    return out


def analyse():
    ref = get("stock", "lcdm")
    result = {"validation": {}, "absorber_vs_thomson": {}, "decomposition": {}, "tca": {}, "smoothness": {}}

    # 1. Patched (Thomson coefficients) must reproduce stock.
    lc = get("patched", "lcdm")
    result["validation"]["lcdm"] = {**max_abs(residuals(lc, ref)),
                                    "bit_identical_tt": bool(np.array_equal(lc["tt"], ref["tt"]))}
    for u in TARGET_U + EXTRA_U:
        stock = get("stock", name("stock", u))
        for kernel in ("thomson_default", "thomson_explicit"):
            m = get("patched", name(kernel, u))
            res = max_abs(residuals(m, stock))
            res["bit_identical_tt"] = bool(np.array_equal(m["tt"], stock["tt"]))
            result["validation"][f"{kernel} u={u:g}"] = res

    # 2. Absorber vs Thomson, and the two single-coefficient variants.
    for u in TARGET_U + EXTRA_U:
        th = get("patched", name("thomson_explicit", u))
        result["absorber_vs_thomson"][f"u={u:g}"] = summarise(get("patched", name("absorber", u)), th, ref)
        for kernel in ("no_quadrupole_regen", "no_polarization_source"):
            result["decomposition"][f"{kernel} u={u:g}"] = summarise(get("patched", name(kernel, u)), th, ref)
        stock = get("stock", name("stock", u))
        result["absorber_vs_thomson"][f"u={u:g}"]["thomson_effect_vs_lcdm_chi2cv_tt"] = cv_chi2(th, ref, "tt")
        result["absorber_vs_thomson"][f"u={u:g}"]["thomson_effect_vs_lcdm_chi2cv_ee"] = cv_chi2(th, ref, "ee")
        ab = get("patched", name("absorber", u))
        for key in ("tt", "te", "ee"):
            # what matters for a bound: each kernel's distance from LCDM (the kernel
            # difference can add coherently to the signal, so the cross term counts)
            result["absorber_vs_thomson"][f"u={u:g}"][f"absorber_effect_vs_lcdm_chi2cv_{key}"] = cv_chi2(ab, ref, key)
            result["absorber_vs_thomson"][f"u={u:g}"][f"thomson_effect_vs_lcdm_chi2cv_{key}"] = cv_chi2(th, ref, key)

    # 3. Tight-coupling sensitivity at u = 2.25e-4.
    base = summarise(get("patched", name("absorber", 2.25e-4)), get("patched", name("thomson_explicit", 2.25e-4)), ref)
    result["tca"]["project settings"] = base
    for variant in TCA_VARIANTS:
        a = get("patched", name(variant, 2.25e-4))
        t = get("patched", name(variant.replace("absorber", "thomson"), 2.25e-4))
        result["tca"][variant] = summarise(a, t, ref)

    # 4. Glitch check on the kernel difference D(u) = r_absorber - r_thomson (vs LCDM).
    #    D is linear in u (it scales with the clump optical depth), so D/u should vary
    #    smoothly; a glitch in either run shows up as a jump in D/u between neighbours.
    #    Metric: max over l of |D_i/u_i - (D_{i-1}/u_{i-1} + D_{i+1}/u_{i+1})/2| * u_i, in
    #    the same units as C_l residuals, compared with max |D_i| itself.
    ell = ref["ell"][2:]
    for key in ("tt_unlensed", "ee_unlensed", "te_unlensed", "tt", "ee", "te"):
        us = 10**FINE_LOG10_U
        d = np.array([frac(get("patched", name("absorber", u)), ref, key)
                      - frac(get("patched", name("thomson_default", u)), ref, key) for u in us])
        n = d / us[:, None]
        jumps = [float(np.max(np.abs(n[i] - 0.5 * (n[i - 1] + n[i + 1]))[ell > 2]) * us[i])
                 for i in range(1, len(us) - 1)]
        sizes = [float(np.max(np.abs(d[i][ell > 2]))) for i in range(1, len(us) - 1)]
        result["smoothness"][key] = {"max_jump": max(jumps), "max_jump_over_difference": max(j / s for j, s in zip(jumps, sizes)),
                                     "n_jump_above_3e-4": int(sum(j > 3e-4 for j in jumps))}

    with open(OUT / "summary.json", "w") as f:
        json.dump(result, f, indent=2)
    plot(ref)
    print(json.dumps({k: v for k, v in result.items() if k != "validation"}, indent=1)[:6000])
    worst_val = max(max(v for k, v in r.items() if k != "bit_identical_tt") for r in result["validation"].values())
    print("validation: worst |patched Thomson - stock| =", worst_val,
          "bit-identical TT in", sum(r["bit_identical_tt"] for r in result["validation"].values()),
          "of", len(result["validation"]))


def plot(ref):
    ell = ref["ell"][2:]
    fig, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True, constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(0, 0.85, len(TARGET_U) + 1))
    for u, c in zip(TARGET_U + [1e-2], colors):
        th = get("patched", name("thomson_explicit", u))
        ab = get("patched", name("absorber", u))
        for ax, key in zip(axes, ("tt", "te", "ee")):
            d = frac(ab, ref, key) - frac(th, ref, key)
            ax.plot(ell, d, color=c, lw=0.8, label=f"u = {u:g}")
            d_unl = frac(ab, ref, key + "_unlensed") - frac(th, ref, key + "_unlensed")
            ax.plot(ell, d_unl, color=c, lw=0.8, ls="--")
    for ax, key in zip(axes, ("TT", "TE / √(TT EE)", "EE")):
        ax.set_xscale("log")
        ax.set_yscale("symlog", linthresh=1e-6)
        ax.axhline(0, color="k", lw=0.5)
        ax.set_ylabel(f"Δ{key}: absorber − Thomson\n(relative to ΛCDM)")
    axes[0].set_title("M6: absorber − Thomson kernel at equal u (solid lensed, dashed unlensed)", fontsize=10)
    axes[0].legend(fontsize=7)
    axes[2].set_xlabel(r"$\ell$")
    fig.savefig(ROOT / "figures" / "m6_absorber_vs_thomson.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    {"run": run, "analyse": analyse}[sys.argv[1]]()
