"""M3 (Experiment A): scan the clump-photon coupling u at f_cl = 1.

Two passes over the same u grid:
  fixed_h:      h held at the baseline value, so every standard parameter is
                physically unchanged. Used for all physical diagnostics.
  fixed_theta:  100*theta_s held, as in baseline_lcdm. CLASS's theta_s is
                evaluated at z_rec, which clump opacity shifts, so h drifts;
                this pass records that drift.
Results go to results/scans/m3/<pass>/; analyse with src/m3_analyse.py.
"""

import json

import numpy as np

from cosmology import ROOT, class_params, fix_h, load_config
from run_class import load
from scan import run_grid

OUT = ROOT / "results" / "scans" / "m3"

LOG10_U = sorted(set(
    [-8.0, -7.0, -6.0]                                    # null checks
    + list(np.round(np.arange(-5.0, -2.0 + 1e-9, 0.125), 3))  # dense, CMB-sensitive
    + [-1.75, -1.5, -1.25, -1.0, -0.75, -0.5, -0.25]      # strongly coupled / find CLASS limit
))


def name_for(log10_u):
    return f"u{log10_u:+.3f}"


def main():
    config = load_config()
    h = load(config["name"])[1]["derived"]["h"]
    status = {}
    for label, build in (
        ("fixed_h", lambda lu: fix_h(class_params(config, f_cl=1.0, u=10**lu), h)),
        ("fixed_theta", lambda lu: class_params(config, f_cl=1.0, u=10**lu)),
    ):
        jobs = [(name_for(lu), build(lu)) for lu in LOG10_U]
        jobs.append(("lcdm", fix_h(class_params(config), h) if label == "fixed_h" else class_params(config)))
        status[label] = run_grid(jobs, config, OUT / label)
        failed = {k: v for k, v in status[label].items() if v}
        print(f"{label}: {len(jobs) - len(failed)}/{len(jobs)} ok")
        for k, v in failed.items():
            print(f"  {k} failed: {v[:160]}")
    with open(OUT / "status.json", "w") as f:
        json.dump({"log10_u": LOG10_U, "h_fixed": h, "status": status}, f, indent=2)


if __name__ == "__main__":
    main()
