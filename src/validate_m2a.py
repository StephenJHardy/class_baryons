"""M2a: f_cl = 1, u = 0 idm must reproduce baseline_lcdm to 1e-4."""

import json

import matplotlib.pyplot as plt
import numpy as np

from cosmology import ROOT, class_params, load_config
from diagnostics import TOLERANCE, max_abs, residuals
from run_class import load, run, save

LABELS = {"tt": "TT", "te": "TE / sqrt(TT EE)", "ee": "EE", "pp": r"$\phi\phi$"}


def plot(res, path):
    fig, axes = plt.subplots(5, 1, figsize=(7, 11), constrained_layout=True)
    for ax, key in zip(axes, ("tt", "te", "ee", "pp")):
        ax.plot(res["ell"], res[key], lw=0.8)
        ax.set_xscale("log")
        ax.set_ylabel(f"Δ{LABELS[key]}")
    axes[3].set_xlabel(r"$\ell$")
    axes[4].plot(res["k_h"], res["pk"], lw=0.8)
    axes[4].set_xscale("log")
    axes[4].set_xlabel(r"$k$ [h/Mpc]")
    axes[4].set_ylabel("ΔP(k)/P(k)")
    for ax in axes:
        ax.axhspan(-TOLERANCE, TOLERANCE, color="0.9", zorder=0)
        ax.axhline(0, color="0.5", lw=0.5)
    axes[0].set_title("M2a: idm (f_cl=1, u=0) vs baseline_lcdm; grey band = ±1e-4")
    fig.savefig(path, dpi=150)


def main():
    config = load_config()
    ref, ref_meta = load(config["name"])
    out, meta = run(class_params(config, f_cl=1.0, u=0.0), config)
    save("m2a_idm_u0", out, meta)
    res = residuals(out, ref)
    worst = max_abs(res)
    summary = {
        "tolerance": TOLERANCE,
        "max_abs_residual": worst,
        "pass": all(v <= TOLERANCE for v in worst.values()),
        "baseline_params_match": ref_meta["params"] == class_params(config),
    }
    (ROOT / "results" / "validation").mkdir(parents=True, exist_ok=True)
    with open(ROOT / "results" / "validation" / "m2a.json", "w") as f:
        json.dump(summary, f, indent=2)
    plot(res, ROOT / "figures" / "m2a_residuals.png")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
