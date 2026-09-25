"""M12 cross-check: CLASS thermal-relic WDM half-mode scale vs the Viel et al. (2005) fit.

Thermal relic (one fermion family, deg_ncdm = 1) with omega_x = omega_dm:
    T_x / T_gamma = (4/11)^(1/3) (94.1 eV omega_x / m_x)^(1/3),
passed to CLASS as a second ncdm species (the first is the 0.06 eV neutrino).
m_ncdm and T_ncdm are given explicitly: giving omega_ncdm with m_ncdm instead
would renormalise the phase-space density at the wrong temperature.

Usage: uv run python src/m12_wdm_check.py
"""

import json
import time

import numpy as np

from cosmology import ROOT, class_params, fix_h, load_config
from m12_small_scale import K_MAX, K_PER_DECADE, N_K, OUT, wdm_k_hm
from run_class import Class, load
from small_scale import first_crossing

MASSES_KEV = [1.24, 3.5, 5.3, 6.5]


def pk(p, h):
    c = Class(); c.set(p); t0 = time.time(); c.compute()
    k = np.logspace(-3, np.log10(0.95 * K_MAX), N_K)
    P = np.array([c.pk_lin(kk * h, 0.0) for kk in k]) * h**3
    om = c.get_current_derived_parameters(["Omega_m"])["Omega_m"]
    c.struct_cleanup(); c.empty()
    return k, P, om, time.time() - t0


def main():
    config = load_config()
    h = load(config["name"])[1]["derived"]["h"]
    omega_x = config["cosmology"]["omega_dm"]
    base = fix_h(class_params(config), h)
    for key in ("lensing", "l_max_scalars"):
        base.pop(key)
    base.update({"output": "mPk", "P_k_max_h/Mpc": K_MAX, "k_per_decade_for_pk": K_PER_DECADE})
    k, p0, om0, _ = pk(base, h)
    rows = []
    for m in MASSES_KEV:
        t_x = (4 / 11) ** (1 / 3) * (94.1 * omega_x / (m * 1e3)) ** (1 / 3)
        p = dict(base)
        p.update({"omega_cdm": 0.0, "N_ncdm": 2, "m_ncdm": f"0.06, {m * 1e3}", "T_ncdm": f"0.71611, {t_x}"})
        k, pw, om, dt = pk(p, h)
        k_hm = first_crossing(k, pw / p0, 0.25)
        k_fit = wdm_k_hm(m, omega_x / h**2, h)
        rows.append({"m_WDM_keV": m, "T_x_over_T_gamma": t_x, "Omega_m_wdm_run": om, "Omega_m_lcdm": om0,
                     "k_hm_class": k_hm, "k_hm_viel_fit": k_fit, "ratio_class_over_fit": k_hm / k_fit,
                     "runtime_s": dt})
        print(f"m = {m} keV: k_hm CLASS {k_hm:.2f}, Viel fit {k_fit:.2f} h/Mpc (ratio {k_hm / k_fit:.3f}); "
              f"Omega_m {om:.5f} vs {om0:.5f}; {dt:.0f} s", flush=True)
    with open(OUT / "wdm_check.json", "w") as f:
        json.dump(rows, f, indent=2)


if __name__ == "__main__":
    main()
