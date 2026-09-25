"""M13 (Experiment F): energy exchange between opaque clumps and the CMB (estimates).

Clumps are opaque, grey absorbers (Q = 1) of surface density Sigma, all of the
dark matter (f_cl = 1), with photon absorption rate Gamma_gamma = (1+z) dmu_idm_g
(CLASS's opacity; for an absorber the absorption and momentum-transfer
cross-sections coincide, so u carries over). Composition: ionised H/He,
mean molecular weight mu_mol = 0.59, c_v = (3/2) k / (mu_mol m_p).

F1  Thermal locking and survival.
    t_th = Sigma c_v / (4 c a T^3): time for the clump surface to relax to T_gamma
    (absorbed - emitted power = pi R^2 c a (T_gamma^4 - T_cl^4) per clump).
    Lag of a clump cooling with the CMB: dT/T = H t_th.
    Gravitational binding of material at T_gamma: G M mu_mol m_p / R >= k T_gamma with
    M = pi R^2 Sigma  ->  M_min(z, Sigma) = (k T_gamma / (G mu_mol m_p))^2 / (pi Sigma).
    Density form (Jeans-like): a uniform clump of internal density rho bound at T_gamma needs
    M >= (4/3) pi rho (3 k T / (4 pi G rho mu_mol m_p))^(3/2); equivalently it survives only
    below the redshift where T_gamma = (4 pi G rho mu_mol m_p / (3 k)) (3 M / (4 pi rho))^(2/3).
    Skin heating: the Thomson-thick column heated by radiative diffusion in a Hubble
    time is Sigma_skin = sqrt(t_H rho c / kappa_T), with rho = 3 Sigma / (4 R).
F2  Energy exchange rate, as d(Q/rho_gamma)/dz:
    (a) passive, thermally locked blackbody clumps release their heat as they cool
        with the CMB: dQ/dt = rho_cl c_v H T_gamma;
    (b) a sustained offset delta = (T_cl - T_gamma)/T_gamma: dQ/dt = 4 Gamma_gamma delta rho_gamma;
    (c) an internal luminosity per mass L/M: dQ/dt = rho_cl L/M.
F3  mu = 1.401 int d(Q/rho_gamma)/dz J_mu dz,  y = (1/4) int d(Q/rho_gamma)/dz J_y dz
    with the Chluba (2016, MNRAS; arXiv:1603.02496, eqs 7, 8, 11, 'Method C') visibilities
        J_bb = exp(-(z/z_th)^(5/2)),  z_th = 1.98e6
        J_y  = [1 + ((1+z)/6e4)^2.58]^-1   (z >= z_rec; 0 otherwise)
        J_mu = J_bb [1 - exp(-((1+z)/5.8e4)^1.88)]
    compared with COBE/FIRAS |mu| < 9e-5, |y| < 1.5e-5 (95%; Fixsen et al. 1996, ApJ 473, 576).
    Also reported: the redshift above which clumps absorb and re-emit every photon at
    least once per Hubble time (Gamma_gamma > H), i.e. act as a thermaliser.

Usage: uv run python src/m13_distortions.py
"""

import json

import matplotlib.pyplot as plt
import numpy as np

from cloud_mapping import sigma_surface_over_q
from cosmology import ROOT, class_params, fix_h, load_config
from rates import analytic_dmu
from run_class import Class, load

OUT = ROOT / "results" / "scans" / "m13"
# cgs constants
K_B, M_P, G_N, C, A_RAD = 1.380649e-16, 1.67262e-24, 6.674e-8, 2.99792458e10, 7.5657e-15
MPC, MSUN, LSUN = 3.0857e24, 1.989e33, 3.828e33
MU_MOL = 0.59
C_V = 1.5 * K_B / (MU_MOL * M_P)
KAPPA_T = 0.35                      # cm^2/g, ionised H/He Thomson opacity
Z_TH = 1.98e6
FIRAS = {"mu": 9e-5, "y": 1.5e-5}

COUPLINGS = {                       # u values of interest (M3, M12)
    "CMB bound (Planck)": 2.25e-4,
    "MW satellites direct (Boehm+2014)": 5.5e-7,
    "small-scale via WDM half-mode (M12)": 6.1e-9,
}
CLUMP_MASSES_MSUN = [1e-6, 1e-3, 1.0]


def j_bb(z):
    return np.exp(-(z / Z_TH) ** 2.5)


def j_mu(z):
    return j_bb(z) * (1 - np.exp(-((1 + z) / 5.8e4) ** 1.88))


def j_y(z, z_rec):
    return np.where(z >= z_rec, 1 / (1 + ((1 + z) / 6e4) ** 2.58), 0.0)


def background(config, h):
    p = fix_h(class_params(config), h)
    for key in ("lensing", "l_max_scalars", "P_k_max_h/Mpc", "z_pk"):
        p.pop(key)
    p["output"] = ""
    c = Class(); c.set(p); c.compute()
    bg = c.get_background()
    z_rec = c.get_current_derived_parameters(["z_rec"])["z_rec"]
    c.struct_cleanup(); c.empty()
    z = np.geomspace(1e2, 1e8, 4000)
    order = np.argsort(bg["z"])
    f = lambda key: np.exp(np.interp(np.log(z), np.log(bg["z"][order] + 1e-30), np.log(bg[key][order])))
    to_cgs = 3 * C**2 / (8 * np.pi * G_N) / MPC**2   # CLASS rho [Mpc^-2] -> g/cm^3 (times c^2 -> erg/cm^3)
    return {"z": z, "H": f("H [1/Mpc]") * C / MPC,                  # 1/s
            "rho_gamma_erg": f("(.)rho_g") * to_cgs * C**2,         # erg/cm^3
            "rho_cl_g": f("(.)rho_cdm") * to_cgs,                   # g/cm^3 (f_cl = 1)
            "T": 2.7255 * (1 + z)}, z_rec


def integrate(z, dq_dz, weight):
    return float(np.trapezoid(dq_dz * weight, z))


def main():
    config = load_config()
    h = load(config["name"])[1]["derived"]["h"]
    bg, z_rec = background(config, h)
    z, H, rho_g, rho_cl, T = bg["z"], bg["H"], bg["rho_gamma_erg"], bg["rho_cl_g"], bg["T"]
    omega_cl = config["cosmology"]["omega_dm"]
    out = {"assumptions": {"mu_mol": MU_MOL, "c_v_erg_g_K": C_V, "kappa_T": KAPPA_T, "z_rec": z_rec,
                           "firas_95": FIRAS}, "F1": {}, "F2F3": {}}

    # ---- F1: thermal locking, gravitational binding, skin heating ----
    zs = [1e7, 1e6, 1e5, 4e4, 1e4]
    for label, u in COUPLINGS.items():
        sigma = sigma_surface_over_q(u)                      # Q = 1
        rows = []
        for zz in zs:
            i = np.argmin(np.abs(z - zz))
            t_th = sigma * C_V / (4 * C * A_RAD * T[i] ** 3)
            m_min = (K_B * T[i] / (G_N * MU_MOL * M_P)) ** 2 / (np.pi * sigma)
            row = {"z": float(z[i]), "T_gamma_eV": float(K_B * T[i] / 1.602e-12), "t_th_s": float(t_th),
                   "H_t_th (lag dT/T)": float(H[i] * t_th), "M_min_bound_Msun": float(m_min / MSUN),
                   "skin_fraction_per_Hubble_time": {}}
            for m in CLUMP_MASSES_MSUN:
                r = np.sqrt(m * MSUN / (np.pi * sigma))
                rho_int = 3 * sigma / (4 * r)
                skin = np.sqrt((1 / H[i]) * rho_int * C / KAPPA_T)
                row["skin_fraction_per_Hubble_time"][f"{m:g} Msun"] = float(min(skin / sigma, 1.0))
            rows.append(row)
        out["F1"][label] = {"u": u, "Sigma_over_Q_g_cm2": sigma, "rows": rows}

    # survival redshift for given mass and internal density (Jeans-like criterion)
    out["F1"]["survival_redshift"] = {}
    for rho_int in (0.1, 1.0, 10.0):
        for m in CLUMP_MASSES_MSUN:
            t_max = 4 * np.pi * G_N * rho_int * MU_MOL * M_P / (3 * K_B) * (3 * m * MSUN / (4 * np.pi * rho_int)) ** (2 / 3)
            out["F1"]["survival_redshift"][f"M={m:g} Msun, rho={rho_int:g} g/cm3"] = {
                "T_max_eV": t_max * K_B / 1.602e-12, "z_max": t_max / 2.7255 - 1}

    # ---- F2/F3 ----
    # (a) passive blackbody clumps: heat released as they cool with the CMB (u-independent)
    dq_dz_a = rho_cl * C_V * T / rho_g / (1 + z)            # d(Q/rho_gamma)/dz, positive = heating photons
    passive = {"mu": 1.401 * integrate(z, dq_dz_a, j_mu(z)), "y": 0.25 * integrate(z, dq_dz_a, j_y(z, z_rec))}
    # (c) internal luminosity per unit mass, per L_sun/M_sun, constant in time (u-independent)
    lm = LSUN / MSUN
    dq_dz_c = rho_cl * lm / (rho_g * H * (1 + z))
    per_lm = {"mu": 1.401 * integrate(z, dq_dz_c, j_mu(z)), "y": 0.25 * integrate(z, dq_dz_c, j_y(z, z_rec))}
    lum_bound = {k: FIRAS[k] / v for k, v in per_lm.items()}
    out["F2F3"]["passive_blackbody_clumps"] = passive
    out["F2F3"]["per_unit_L_over_M_in_Lsun_per_Msun"] = per_lm
    out["F2F3"]["FIRAS_bound_on_constant_L_over_M_Lsun_per_Msun"] = lum_bound
    # (b) sustained temperature offset delta, and thermalisation redshift, per coupling
    for label, u in COUPLINGS.items():
        gamma = analytic_dmu(z, omega_cl, u) * (1 + z) * C / MPC      # 1/s
        dq_dz_b = 4 * gamma / (H * (1 + z))                          # per unit delta
        a_mu = 1.401 * integrate(z, dq_dz_b, j_mu(z))
        a_y = 0.25 * integrate(z, dq_dz_b, j_y(z, z_rec))
        ratio = gamma / H
        z_therm = float(z[np.nonzero(ratio >= 1)[0][0]]) if np.any(ratio >= 1) else None
        out["F2F3"][label] = {"u": u, "mu_per_unit_delta": a_mu, "y_per_unit_delta": a_y,
                              "FIRAS_bound_on_delta_from_mu": FIRAS["mu"] / a_mu,
                              "FIRAS_bound_on_delta_from_y": FIRAS["y"] / a_y,
                              "z_where_Gamma_gamma_equals_H": z_therm}
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "summary.json", "w") as f:
        json.dump(out, f, indent=2, default=float)
    plots(bg, z_rec)
    print(json.dumps(out, indent=1, default=float))


def plots(bg, z_rec):
    z, T = bg["z"], bg["T"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    for label, u in COUPLINGS.items():
        sigma = sigma_surface_over_q(u)
        m_min = (K_B * T / (G_N * MU_MOL * M_P)) ** 2 / (np.pi * sigma) / MSUN
        a1.loglog(1 + z, m_min, label=f"Σ/Q = {sigma:.1e} g/cm² ({label})")
        a2.loglog(1 + z, bg["H"] * sigma * C_V / (4 * C * A_RAD * T**3), label=f"Σ/Q = {sigma:.1e}")
    for m in CLUMP_MASSES_MSUN:
        a1.axhline(m, color="0.6", ls=":", lw=0.8)
    a1.set(xlabel="1 + z", ylabel="M_min [M☉] for gravitational binding at T_γ",
           title="F1: minimum clump mass bound against thermal pressure at T_γ")
    a1.legend(fontsize=6.5)
    a2.axhline(1, color="k", lw=0.6)
    a2.set(xlabel="1 + z", ylabel="H t_th  (thermal lag δT/T)", title="F1: thermal locking of clump surfaces")
    a2.legend(fontsize=7)
    fig.savefig(ROOT / "figures" / "m13_thermal.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
