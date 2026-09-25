"""M12 / E4: when does the discreteness of the clumps matter at the half-mode scale?

N clumps of mass M_cl in a smooth pressureless background add a Poisson
(isocurvature) power P_shot = 1/n_cl = M_cl / rho_cl (constant in k) at the
time they form, which then grows with the matter growth factor. For clumps
present from before matter-radiation equality, the isocurvature mode grows by
~ (3/2)(1 + z_eq) in amplitude by today (Meszaros, flat matter domination),
so P_poisson(z=0) ~ [(3/2)(1 + z_eq)]^2 M_cl / rho_cl.

The fluid (smooth-component) treatment used by CLASS is adequate at scale k
while P_poisson(z=0) << P_lin(k, z=0). This order-of-magnitude estimate gives
the clump mass at which the two are equal at the small-scale-limit half-mode
scale. It is a validity boundary for this analysis, not an observational bound.

Usage: uv run python src/m12_fluid_check.py
"""

import json

import numpy as np

from cosmology import ROOT, class_params, fix_h, load_config
from m12_small_scale import OUT
from rates import crossing_redshift, rates_from_class
from run_class import Class, load
from small_scale import RHO_CRIT_MSUN_H2_MPC3


def main():
    config = load_config()
    meta = load(config["name"])[1]
    h = meta["derived"]["h"]
    omega_m = meta["derived"]["Omega_m"]
    omega_dm = config["cosmology"]["omega_dm"] / h**2
    rho_cl = omega_dm * RHO_CRIT_MSUN_H2_MPC3            # (M_sun/h) / (Mpc/h)^3, f_cl = 1
    z_eq = 2.5e4 * omega_m * h**2 / (2.7255 / 2.7) ** 4 - 1   # rough; CLASS value is ~3400
    growth2 = (1.5 * (1 + z_eq)) ** 2
    with np.load(OUT / "lcdm.npz") as d:
        k, pk = d["k_h"], d["pk"]
    summary = json.load(open(OUT / "summary.json"))
    out = {"rho_cl_Msun_h_per_Mpc_h3": rho_cl, "z_eq_estimate": z_eq, "isocurvature_growth_power": growth2, "limits": {}}
    for label, b in summary["bounds_via_WDM_half_mode"].items():
        k_hm = b["k_hm_limit_h_Mpc"]
        p_lin = float(np.exp(np.interp(np.log(k_hm), np.log(k), np.log(pk))))
        m_equal = p_lin * rho_cl / growth2                     # M_sun/h where P_poisson = P_lin at k_hm
        out["limits"][label] = {"k_hm": k_hm, "P_lin_at_k_hm": p_lin, "M_cl_equal_power_Msun_h": m_equal,
                                "M_hm_Msun_h": b["M_hm_limit_Msun_h"]}
        print(f"{label}: k_hm = {k_hm:.1f} h/Mpc, P_lin = {p_lin:.3g} (Mpc/h)^3, "
              f"Poisson = linear power at M_cl ~ {m_equal:.2g} M_sun/h (M_hm = {b['M_hm_limit_Msun_h']:.2g})")
    # Epochs: horizon entry (k tau = 1) of the half-mode scale, and drag decoupling
    # (Gamma_cl->gamma = H) at the bound. The bound tests clumps that already exist,
    # with that Sigma, at z ~ z_enter.
    c = Class()
    p = fix_h(class_params(config), h)
    for key in ("lensing", "l_max_scalars", "P_k_max_h/Mpc", "z_pk"):
        p.pop(key)
    p["output"] = ""
    c.set(p); c.compute()
    bg = c.get_background(); c.struct_cleanup(); c.empty()
    for label, b in summary["bounds_via_WDM_half_mode"].items():
        k_mpc = out["limits"][label]["k_hm"] * h
        order = np.argsort(bg["conf. time [Mpc]"])
        z_enter = float(np.interp(1 / k_mpc, bg["conf. time [Mpc]"][order], bg["z"][order]))
        c = Class()
        pu = fix_h(class_params(config, f_cl=1.0, u=b["u_max"]), h)
        for key in ("lensing", "l_max_scalars", "P_k_max_h/Mpc", "z_pk"):
            pu.pop(key)
        pu["output"] = ""
        c.set(pu); c.compute()
        r = rates_from_class(c); c.struct_cleanup(); c.empty()
        z_dec = crossing_redshift(r["z_th"], r["gamma_cl_over_H"])
        out["limits"][label].update({"z_horizon_entry": z_enter, "z_dec_drag_at_u_max": z_dec,
                                     "T_photon_keV_at_z_enter": 2.7255 * 8.617e-8 * (1 + z_enter)})
        print(f"  {label}: k_hm enters horizon at z = {z_enter:.2g} (T_gamma = {2.7255 * 8.617e-8 * (1 + z_enter):.2g} keV); "
              f"drag decouples at z = {z_dec:.2g} for u = {b['u_max']:.2g}")
    # Late-time collisionality of clumps at the bound (geometric cross-sections):
    #  clump-clump: sigma/M = Q/Sigma, to compare with self-interacting DM scales ~ 1 cm^2/g;
    #  clump-ISM drag: stopping time t = Sigma / (rho_gas v) for a clump moving at v through gas.
    rho_gas = 1.67e-24          # g/cm^3 (1 proton per cm^3)
    v = 2.0e7                   # cm/s (200 km/s)
    year = 3.156e7
    for label, b in summary["bounds_via_WDM_half_mode"].items():
        sigma_min = b["Sigma_over_Q_min_g_cm2"]
        out["limits"][label].update({
            "clump_self_interaction_sigma_over_M_cm2_g_at_Q1": 1.0 / sigma_min,
            "ism_stopping_time_yr": sigma_min / (rho_gas * v) / year})
        print(f"  {label}: clump-clump sigma/M = {1 / sigma_min:.1g} cm^2/g; "
              f"stopping time in 1 cm^-3 gas at 200 km/s = {sigma_min / (rho_gas * v) / year:.1g} yr")
    with open(OUT / "fluid_check.json", "w") as f:
        json.dump(out, f, indent=2)


if __name__ == "__main__":
    main()
