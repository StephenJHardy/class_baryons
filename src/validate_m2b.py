"""M2b: at fixed u, outputs must converge as m_idm grows (cold-mass limit).

P(k) convergence is judged by the change in the transfer function,
dT2 = (P - P_ref) / P_lcdm. At u = 1e-4, P(k) is suppressed by up to ~1e5
at k ~ 10 h/Mpc, and there dP/P is dominated by integrator noise (~1e-3
even at fixed mass). dP/P is still reported, together with that noise floor.
"""

import json

import matplotlib.pyplot as plt
import numpy as np

from cosmology import ROOT, class_params, load_config
from diagnostics import TOLERANCE, max_abs, residuals
from run_class import load, run, save

U = 1e-4
MASSES_EV = [1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e15, 1e18, 1e21, 1e24, 1e27, 1e30]
PLATEAU_CHECK_U = [1e-6, 1e-2]
PLATEAU_CHECK_MASSES_EV = [1e21, 1e24, 1e27]


def converged_metrics(model, ref_model, lcdm):
    res = max_abs(residuals(model, ref_model))
    res["dT2"] = float(np.max(np.abs((model["pk"] - ref_model["pk"]) / lcdm["pk"])))
    return res


def plateau_check(config):
    """Is the output bit-identical above ~1e24 eV at other couplings too?"""
    rows = []
    for u in PLATEAU_CHECK_U:
        outs = [run(class_params(config, f_cl=1.0, u=u, m_idm=m), config)[0]
                for m in PLATEAU_CHECK_MASSES_EV]
        for m, out in zip(PLATEAU_CHECK_MASSES_EV, outs):
            rows.append({"u_idm_g": u, "m_idm_eV": m,
                         "identical_to_heaviest": all(np.array_equal(out[k], outs[-1][k])
                                                      for k in ("tt", "te", "ee", "pp", "pk"))})
    return rows


def noise_floor(config, m):
    """Change from nudging the integrator tolerance at fixed mass."""
    a = run(class_params(config, f_cl=1.0, u=U, m_idm=m), config)[0]
    b = run(class_params(config, f_cl=1.0, u=U, m_idm=m,
                         extra={"tol_perturbations_integration": 0.9e-5}), config)[0]
    lcdm = load(config["name"])[0]
    return converged_metrics(b, a, lcdm)


def main():
    config = load_config()
    ref, _ = load(config["name"])
    runs, failures = {}, {}
    for m in MASSES_EV:
        try:
            out, meta = run(class_params(config, f_cl=1.0, u=U, m_idm=m), config)
        except Exception as e:  # classy raises CosmoSevereError / CosmoComputationError
            failures[f"{m:.0e}"] = str(e).strip().splitlines()[-1]
            print(f"m_idm = {m:.0e} eV failed: {failures[f'{m:.0e}']}")
            continue
        runs[m] = out
        save(f"m2b_u1e-4_m{m:.0e}", out, meta)
        print(f"m_idm = {m:.0e} eV: {meta['runtime_s']:.1f} s", flush=True)

    masses = sorted(runs)
    m_ref = masses[-1]
    table = []
    for m in masses:
        table.append({
            "m_idm_eV": m,
            "vs_heaviest": converged_metrics(runs[m], runs[m_ref], ref),
            "vs_lcdm": max_abs(residuals(runs[m], ref)),
        })
    judged = ("tt", "te", "ee", "pp", "dT2")
    converged = [row["m_idm_eV"] for row in table
                 if all(row["vs_heaviest"][k] <= TOLERANCE for k in judged)]
    # Converged means every mass from here up converges, not just this one.
    lightest = next(m for m in masses if all(x in converged for x in masses if x >= m))
    summary = {"u_idm_g": U, "tolerance": TOLERANCE, "judged_on": judged,
               "reference_mass_eV": m_ref, "lightest_converged_mass_eV": lightest,
               "project_default_m_idm_eV": config["idm"]["m_idm"],
               "noise_floor_at_reference_mass": noise_floor(config, m_ref),
               "plateau_check": plateau_check(config),
               "table": table, "failures": failures}
    with open(ROOT / "results" / "validation" / "m2b.json", "w") as f:
        json.dump(summary, f, indent=2)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for key in ("tt", "te", "ee", "pp", "dT2", "pk"):
        label = "pk (dP/P, noise-dominated tail)" if key == "pk" else key
        ax1.loglog(masses, [max(r["vs_heaviest"][key], 1e-12) for r in table], marker="o",
                   label=label, ls=":" if key == "pk" else "-")
    for key in ("tt", "te", "ee", "pp"):
        ax2.loglog(masses, [r["vs_lcdm"][key] for r in table], "o-", label=key)
    ax1.axhline(TOLERANCE, color="k", ls="--", lw=0.8)
    ax1.set(xlabel="m_idm [eV]", ylabel=f"max |Δ| vs m = {m_ref:.0e} eV",
            title=f"M2b: convergence in m_idm at u = {U:g}")
    ax2.set(xlabel="m_idm [eV]", ylabel="max |Δ| vs baseline_lcdm",
            title="Size of the u = 1e-4 effect itself")
    ax1.axvline(config["idm"]["m_idm"], color="k", lw=0.8, ls=":")
    ax1.legend(fontsize=7)
    ax2.legend(fontsize=7)
    fig.savefig(ROOT / "figures" / "m2b_mass_convergence.png", dpi=150)

    # Residual spectra vs the heaviest mass, for a few masses.
    fig, axes = plt.subplots(3, 1, figsize=(7, 8), constrained_layout=True)
    for m in masses[:-1:2]:
        res = residuals(runs[m], runs[m_ref])
        axes[0].plot(res["ell"], res["tt"], lw=0.8, label=f"{m:.0e} eV")
        axes[1].plot(res["ell"], res["ee"], lw=0.8)
        axes[2].plot(res["k_h"], (runs[m]["pk"] - runs[m_ref]["pk"]) / ref["pk"], lw=0.8)
    for ax, lab in zip(axes, ("ΔTT/TT", "ΔEE/EE", r"ΔP/P$_{\Lambda CDM}$")):
        ax.set_xscale("log"); ax.set_ylabel(lab)
        ax.axhspan(-TOLERANCE, TOLERANCE, color="0.9", zorder=0)
    axes[0].set_title(f"M2b: u = {U:g}, residuals vs m_idm = {m_ref:.0e} eV")
    axes[0].legend(fontsize=7)
    axes[1].set_xlabel(r"$\ell$"); axes[2].set_xlabel("k [h/Mpc]")
    fig.savefig(ROOT / "figures" / "m2b_residuals.png", dpi=150)

    for row in table:
        print(f"{row['m_idm_eV']:.0e}", {k: f"{v:.1e}" for k, v in row["vs_heaviest"].items()},
              "| vs LCDM", {k: f"{v:.1e}" for k, v in row["vs_lcdm"].items()})
    print("lightest converged mass:", f"{lightest:.0e} eV")
    print("noise floor:", summary["noise_floor_at_reference_mass"])
    print("plateau check:", summary["plateau_check"])


if __name__ == "__main__":
    main()
