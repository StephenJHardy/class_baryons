"""Clump-photon interaction rates relative to expansion (plan Section 3).

CLASS stores dmu_idm_g = a n_cl sigma c, the conformal photon opacity due to
clumps [1/Mpc]. The physical rates are
    Gamma_gamma->cl = dmu_idm_g / a          (photon opacity)
    Gamma_cl->gamma = S_idm_g * dmu_idm_g / a (clump momentum drag),
with S_idm_g = 4 rho_gamma / (3 rho_idm), as in CLASS's perturbation equations.
"""

import numpy as np

from cloud_mapping import sigma_over_m

RHO_CRIT_H2_G_CM3 = 1.87834e-29   # 3 H0^2 / (8 pi G) for h = 1
MPC_CM = 3.085677581282e24        # CLASS _Mpc_over_m_ in cm


N_Z_SAVE = 2000   # CLASS's own table has ~1e5 points; resample for storage


def rates_from_class(cosmo):
    """Gamma/H for both rates, resampled onto N_Z_SAVE log-spaced redshifts."""
    th = cosmo.get_thermodynamics()
    if "dmu_idm_g" not in th:
        return {}
    bg = cosmo.get_background()
    z = th["z"]
    order = np.argsort(bg["z"])
    interp = lambda key: np.interp(z, bg["z"][order], bg[key][order])
    hubble = interp("H [1/Mpc]")
    s_idm_g = 4.0 / 3.0 * interp("(.)rho_g") / interp("(.)rho_idm")
    gamma_ph = th["dmu_idm_g"] * (1 + z)
    full = {
        "dmu_idm_g": th["dmu_idm_g"],
        "gamma_ph_over_H": gamma_ph / hubble,
        "gamma_cl_over_H": s_idm_g * gamma_ph / hubble,
        "visibility": th["g [Mpc^-1]"],
        "dkappa_thomson": th["kappa' [Mpc^-1]"],
        "conformal_time": th["conf. time [Mpc]"],
    }
    z_save = np.geomspace(1.0, 1.0 + z.max(), N_Z_SAVE) - 1.0
    out = {"z_th": z_save}
    for key, values in full.items():   # CLASS orders the table by increasing z
        out[key] = np.interp(np.log1p(z_save), np.log1p(z), values)
    return out


def analytic_dmu(z, omega_idm, u):
    """a n sigma c = (1+z)^2 rho_idm,0 (sigma/M) in 1/Mpc, for n_index_idm_g = 0."""
    return (1 + z) ** 2 * omega_idm * RHO_CRIT_H2_G_CM3 * sigma_over_m(u) * MPC_CM


def crossing_redshift(z, ratio):
    """Lowest redshift above which ratio >= 1 at all times (log interpolation).

    Returns nan if ratio < 1 everywhere, inf if ratio >= 1 everywhere.
    """
    order = np.argsort(z)
    z, ratio = z[order], ratio[order]
    above = ratio >= 1
    if not above.any():
        return np.nan
    if above.all():
        return np.inf
    i = np.nonzero(~above)[0][-1]   # last redshift where the rate is below H
    if i == len(z) - 1:
        return np.nan
    lz = np.interp(0.0, np.log([ratio[i], ratio[i + 1]]), np.log1p([z[i], z[i + 1]]))
    return float(np.expm1(lz))


def clump_optical_depth(z, dmu, conformal_time, z_max):
    """Optical depth to photons from clumps between z = 0 and z_max."""
    sel = z <= z_max
    order = np.argsort(conformal_time[sel])
    return float(np.trapezoid(dmu[sel][order], conformal_time[sel][order]))
