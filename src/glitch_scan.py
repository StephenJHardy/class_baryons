"""Count sporadic C_l glitches across a fine u grid for candidate precision settings.

A glitch is a model whose residual (vs LCDM) departs from the average of its
two u-neighbours by more than GLITCH_LEVEL somewhere in 2 <= l <= 2500.
Physical responses are smooth in log u, so the second difference is small
for them (see the 2026-09-25 M3 experiment-log entry).

Usage: uv run python src/glitch_scan.py <setting> [<setting> ...]
"""

import json
import sys

import numpy as np

from cosmology import ROOT, class_params, fix_h, load_config
from m3_analyse import frac
from run_class import load
from scan import run_grid

GLITCH_LEVEL = 3e-4
LOG10_U = np.round(np.arange(-7.0, -3.99, 0.05), 3)
KEYS = ("tt_unlensed", "ee_unlensed", "te_unlensed", "tt", "ee", "te", "pp")
SETTINGS = {
    "project": {},
    "Nz_lin 100000": {"thermo_Nz_lin": 100000},
    "Nz_lin 120000": {"thermo_Nz_lin": 120000},
    "Nz_lin 160000": {"thermo_Nz_lin": 160000},
    "l_linstep 20": {"l_linstep": 20},
    "Nz_lin 100000 + l_linstep 20": {"thermo_Nz_lin": 100000, "l_linstep": 20},
    "l_linstep 10": {"l_linstep": 10},
    "l_linstep 20 + l_logstep 1.06": {"l_linstep": 20, "l_logstep": 1.06},
    "l_linstep 10 + l_logstep 1.03": {"l_linstep": 10, "l_logstep": 1.03},
}


def glitches(setting):
    config = load_config()
    config["precision"] = {**config["precision"], **SETTINGS[setting]}
    h = load(config["name"])[1]["derived"]["h"]
    out = ROOT / "results" / "scratch_glitch" / setting.replace(" ", "_").replace("+", "and")
    jobs = [(f"u{lu:+.3f}", fix_h(class_params(config, f_cl=1.0, u=10**lu), h)) for lu in LOG10_U]
    jobs.append(("lcdm", fix_h(class_params(config), h)))
    run_grid(jobs, config, out, workers=8)
    ref = load("lcdm", out)[0]
    res = {key: np.array([frac(load(f"u{lu:+.3f}", out)[0], ref, key) for lu in LOG10_U]) for key in KEYS}
    rows = []
    for i in range(1, len(LOG10_U) - 1):
        worst = {key: float(np.max(np.abs(r[i] - 0.5 * (r[i - 1] + r[i + 1])))) for key, r in res.items()}
        rows.append({"log10_u": float(LOG10_U[i]), "second_difference": worst})
    n_bad = sum(any(v > GLITCH_LEVEL for k, v in row["second_difference"].items() if k != "pp") for row in rows)
    return {"setting": setting, "extra_precision": SETTINGS[setting], "n_models": len(rows),
            "n_glitchy": n_bad, "rows": rows}


def main():
    results = [glitches(s) for s in sys.argv[1:]]
    for r in results:
        worst = max(max(v for k, v in row["second_difference"].items() if k != "pp") for row in r["rows"])
        print(f"{r['setting']:32s} glitchy {r['n_glitchy']:2d}/{r['n_models']}  worst CMB 2nd-diff {worst:.1e}")
    path = ROOT / "results" / "validation" / "m3_glitch_scan.json"
    old = json.load(open(path)) if path.exists() else {}
    old.update({r["setting"]: r for r in results})
    with open(path, "w") as f:
        json.dump(old, f, indent=2)


if __name__ == "__main__":
    main()
