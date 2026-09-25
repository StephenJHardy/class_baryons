"""Reproduce the M2a precision investigation (see the 2026-09-25 M2a log entry).

Each test compares an idm (u = 0) or perturbed-LCDM run against a reference
LCDM run at the same precision, starting from CLASS default precision
(the config's precision block is deliberately *not* used here).
"""

import json

import numpy as np

from cosmology import ROOT, class_params, load_config
from diagnostics import max_abs, residuals
from run_class import run

NZ4 = {"thermo_Nz_lin": 80000, "thermo_Nz_log": 20000}
EARLY = {"start_small_k_at_tau_c_over_tau_h": 1.5e-4, "start_large_k_at_tau_h_over_tau_k": 7.0e-3}


def summarise(model, ref):
    res = residuals(model, ref)
    out = max_abs(res)
    low = res["ell"] < 30
    out["low_ell_max"] = float(max(np.abs(res[k][low]).max() for k in ("tt", "te", "ee")))
    return out


def main():
    config = load_config()
    config["precision"] = {}
    lcdm = lambda extra=None: run(class_params(config, extra=extra), config)[0]
    idm = lambda extra=None, **kw: run(class_params(config, f_cl=1.0, u=0.0, extra=extra, **kw), config)[0]
    nlog_1e9 = int(5000 / np.log(5e6) * np.log(1e9))
    tests = {}

    ref = lcdm()
    tests["1_default: idm vs lcdm"] = summarise(idm(), ref)
    tests["2_default: lcdm(thermo_z_initial=1e9, rescaled Nz_log) vs lcdm"] = summarise(
        lcdm({"thermo_z_initial": 1e9, "thermo_Nz_log": nlog_1e9}), ref)

    ref4 = lcdm(NZ4)
    tests["3_Nz x4: lcdm(thermo_z_initial=1e9) vs lcdm"] = summarise(
        lcdm({**NZ4, "thermo_z_initial": 1e9, "thermo_Nz_log": int(4 * nlog_1e9)}), ref4)
    tests["4_Nz x4: idm vs lcdm"] = summarise(idm(NZ4), ref4)
    tests["5_Nz x4: idm(m_idm=1e15 eV) vs lcdm"] = summarise(idm(NZ4, m_idm=1e15), ref4)
    tests["6_Nz x4: idm(f_cl=0.1) vs lcdm"] = summarise(
        run(class_params(config, f_cl=0.1, u=0.0, extra=NZ4), config)[0], ref4)

    ref_final = lcdm({**NZ4, **EARLY})
    tests["7_Nz x4 + early start: idm vs lcdm"] = summarise(idm({**NZ4, **EARLY}), ref_final)

    path = ROOT / "results" / "validation" / "m2a_diagnosis.json"
    with open(path, "w") as f:
        json.dump(tests, f, indent=2)
    for name, row in tests.items():
        print(f"{name:66s}", {k: f"{v:.1e}" for k, v in row.items()})


if __name__ == "__main__":
    main()
