"""M12 (Experiment E): small-scale structure bound on the clump-photon coupling.

E2: linear P(k, z=0) to k = 300 h/Mpc for f_cl = 1, log10 u = -9 ... -3, giving
    T^2 = P/P_LCDM and the suppression scales k_10, k_1/2, k_hm, M_hm.
E3: approximate bound by half-mode matching to thermal-relic warm dark matter,
    using the Viel et al. (2005) WDM transfer function (hep-ph/0501562, eqs 6-7):
        T_WDM(k) = [1 + (alpha k)^(2 nu)]^(-5/nu),  nu = 1.12,
        alpha = 0.049 (m/keV)^-1.11 (Omega_x/0.25)^0.11 (h/0.7)^1.22  h^-1 Mpc.
    A clump model is taken to be as suppressed as a WDM model with the same
    half-mode wavenumber (T = 1/2). The verified WDM limits used are
        m_WDM > 6.5 keV  (95%, DES+PS1 Milky Way satellites; Nadler et al. 2021, PRL 126, 091101)
        m_WDM > 5.3 keV  (2 sigma, Lyman-alpha XQ-100 + HIRES/MIKE; Irsic et al. 2017, PRD 96, 023522)
        m_WDM > 3.5 keV  (same, allowing non-smooth IGM temperature evolution)
    and are compared with the only verified *direct* DM-photon small-scale bound,
        sigma_DM-gamma < 5.5e-9 sigma_T (m/GeV) (2 sigma, MW satellites;
        Boehm et al. 2014, MNRAS 445, L31, arXiv:1404.7012) -> u < 5.5e-7.

Usage: uv run python src/m12_small_scale.py run | analyse
"""

import json
import sys
import time

import matplotlib.pyplot as plt
import numpy as np

from cloud_mapping import sigma_over_m, sigma_surface_over_q
from cosmology import ROOT, class_params, fix_h, load_config
from run_class import CLASS_BUILD, Class, load
from small_scale import RHO_CRIT_MSUN_H2_MPC3, first_crossing

OUT = ROOT / "results" / "scans" / "m12"
LOG10_U = np.round(np.arange(-9.0, -2.99, 0.25), 3)
K_MAX = 300.0            # h/Mpc (P_k_max_h/Mpc); k_hm at u = 1e-9 is ~190 h/Mpc
K_PER_DECADE = 80        # k_hm converged to < 0.2% at 40/decade (default 10 gives ~5% errors)
N_K = 1200
WDM_NU = 1.12
WDM_LIMITS = {
    "MW satellites (DES+PS1), m > 6.5 keV, 95%": 6.5,
    "Lyman-alpha (XQ-100+HIRES/MIKE), m > 5.3 keV, 2 sigma": 5.3,
    "Lyman-alpha, non-smooth IGM T(z), m > 3.5 keV, 2 sigma": 3.5,
}
BOEHM_2014_U_MAX = 5.5e-9 * 100          # sigma < 5.5e-9 sigma_T (m/GeV)  ->  u = (sigma/sigma_T)(100 GeV/m)
BOEHM_2014_WDM_POINT = (2e-9 * 100, 1.24)  # their Fig. 1: u = 2e-7 has linear P(k) similar to 1.24 keV WDM


def params(config, h, u):
    p = fix_h(class_params(config, f_cl=1.0, u=u) if u > 0 else class_params(config), h)
    for key in ("lensing", "l_max_scalars"):
        p.pop(key)
    p.update({"output": "mPk", "P_k_max_h/Mpc": K_MAX, "k_per_decade_for_pk": K_PER_DECADE})
    return p


def compute(job):
    name, p, h = job
    c = Class()
    c.set(p)
    t0 = time.time()
    c.compute()
    k = np.logspace(-3, np.log10(0.95 * K_MAX), N_K)
    pk = np.array([c.pk_lin(kk * h, 0.0) for kk in k]) * h**3
    omega_m = c.Omega_m()
    c.struct_cleanup(); c.empty()
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez(OUT / f"{name}.npz", k_h=k, pk=pk)
    with open(OUT / f"{name}.json", "w") as f:
        json.dump({"params": p, "Omega_m": omega_m, "runtime_s": time.time() - t0,
                   "class_build": CLASS_BUILD}, f, indent=2, default=float)
    return name


def run():
    from concurrent.futures import ProcessPoolExecutor
    config = load_config()
    h = load(config["name"])[1]["derived"]["h"]
    jobs = [("lcdm", params(config, h, 0.0), h)]
    jobs += [(f"u{lu:+.3f}", params(config, h, 10**lu), h) for lu in LOG10_U]
    with ProcessPoolExecutor(4) as pool:
        for name in pool.map(compute, jobs):
            print(name, flush=True)


def wdm_alpha(m_kev, omega_x, h):
    return 0.049 * m_kev**-1.11 * (omega_x / 0.25) ** 0.11 * (h / 0.7) ** 1.22   # h^-1 Mpc


def wdm_t2(k, m_kev, omega_x, h):
    return (1 + (wdm_alpha(m_kev, omega_x, h) * k) ** (2 * WDM_NU)) ** (-10 / WDM_NU)


def wdm_k_hm(m_kev, omega_x, h):
    """k where T_WDM = 1/2 (h/Mpc)."""
    return (2 ** (WDM_NU / 5) - 1) ** (1 / (2 * WDM_NU)) / wdm_alpha(m_kev, omega_x, h)


def wdm_mass_from_k_hm(k_hm, omega_x, h):
    a_needed = (2 ** (WDM_NU / 5) - 1) ** (1 / (2 * WDM_NU)) / k_hm
    return (a_needed / (0.049 * (omega_x / 0.25) ** 0.11 * (h / 0.7) ** 1.22)) ** (-1 / 1.11)


def m_hm(k_hm, omega_m):
    return 4 * np.pi / 3 * omega_m * RHO_CRIT_MSUN_H2_MPC3 * (np.pi / k_hm) ** 3


def analyse():
    config = load_config()
    h = load(config["name"])[1]["derived"]["h"]
    omega_x = config["cosmology"]["omega_dm"] / h**2
    with np.load(OUT / "lcdm.npz") as d:
        k, p0 = d["k_h"], d["pk"]
    omega_m = json.load(open(OUT / "lcdm.json"))["Omega_m"]
    rows, t2s = [], {}
    for lu in LOG10_U:
        with np.load(OUT / f"u{lu:+.3f}.npz") as d:
            t2 = d["pk"] / p0
        t2s[lu] = t2
        row = {"log10_u": float(lu), "u": 10**lu, "sigma_over_M_cm2_g": sigma_over_m(10**lu),
               "Sigma_over_Q_g_cm2": sigma_surface_over_q(10**lu),
               "k_10": first_crossing(k, t2, 0.9), "k_half": first_crossing(k, t2, 0.5),
               "k_hm": first_crossing(k, t2, 0.25)}
        row["M_hm"] = m_hm(row["k_hm"], omega_m)
        row["m_WDM_equivalent_keV"] = wdm_mass_from_k_hm(row["k_hm"], omega_x, h)
        rows.append(row)

    lu = np.array([r["log10_u"] for r in rows])
    lk = np.log10([r["k_hm"] for r in rows])
    slope, intercept = np.polyfit(lu, lk, 1)
    fit_resid = float(np.max(np.abs(lk - (slope * lu + intercept))))

    def u_for_k_hm(k_target):   # interpolate log10 u(log10 k_hm); k_hm decreases with u
        return float(10 ** np.interp(np.log10(k_target), lk[::-1], lu[::-1]))

    bounds = {}
    for label, m in WDM_LIMITS.items():
        k_lim = wdm_k_hm(m, omega_x, h)
        u_max = u_for_k_hm(k_lim)
        bounds[label] = {"m_WDM_keV": m, "k_hm_limit_h_Mpc": float(k_lim), "M_hm_limit_Msun_h": float(m_hm(k_lim, omega_m)),
                         "u_max": u_max, "sigma_over_M_max_cm2_g": sigma_over_m(u_max),
                         "Sigma_over_Q_min_g_cm2": sigma_surface_over_q(u_max)}
    # Systematic: the same limits with k_hm from CLASS's own thermal-relic WDM
    # (src/m12_wdm_check.py) instead of the Viel fit, which is calibrated at k < 5 h/Mpc.
    check_path = OUT / "wdm_check.json"
    if check_path.exists():
        by_mass = {r["m_WDM_keV"]: r for r in json.load(open(check_path))}
        for label, b in bounds.items():
            if b["m_WDM_keV"] in by_mass:
                k_class = by_mass[b["m_WDM_keV"]]["k_hm_class"]
                b["k_hm_limit_class_wdm_h_Mpc"] = k_class
                b["u_max_class_wdm"] = u_for_k_hm(k_class)
                b["Sigma_over_Q_min_class_wdm_g_cm2"] = sigma_surface_over_q(b["u_max_class_wdm"])
    boehm = {"u_max": BOEHM_2014_U_MAX, "Sigma_over_Q_min_g_cm2": sigma_surface_over_q(BOEHM_2014_U_MAX),
             "our_WDM_equivalent_at_u_max_keV": wdm_mass_from_k_hm(10 ** (slope * np.log10(BOEHM_2014_U_MAX) + intercept), omega_x, h)}
    u_cal, m_cal = BOEHM_2014_WDM_POINT
    calibration = {"boehm_u": u_cal, "boehm_quoted_similar_WDM_keV": m_cal,
                   "our_half_mode_equivalent_keV": wdm_mass_from_k_hm(10 ** np.interp(np.log10(u_cal), lu, lk), omega_x, h)}
    cmb = {"published_planck_u_max": 2.25e-4, "Sigma_over_Q_min_g_cm2": sigma_surface_over_q(2.25e-4)}

    # shape comparison: best-matched WDM vs clump T^2 at the satellite bound
    m_sat = WDM_LIMITS["MW satellites (DES+PS1), m > 6.5 keV, 95%"]
    shape = {}
    lu_near = lu[np.argmin(np.abs(lu - np.log10(bounds[next(iter(bounds))]["u_max"])))]
    t2_clump = t2s[lu_near]
    m_match = wdm_mass_from_k_hm(first_crossing(k, t2_clump, 0.25), omega_x, h)
    t2_wdm = wdm_t2(k, m_match, omega_x, h)
    sel = k < 0.9 * K_MAX
    shape = {"log10_u": float(lu_near), "matched_m_WDM_keV": m_match,
             "k_10_clump": first_crossing(k, t2_clump, 0.9), "k_10_wdm": first_crossing(k, t2_wdm, 0.9),
             "k_half_clump": first_crossing(k, t2_clump, 0.5), "k_half_wdm": first_crossing(k, t2_wdm, 0.5),
             "T2_at_2k_hm_clump": float(np.interp(2 * first_crossing(k, t2_clump, 0.25), k[sel], t2_clump[sel])),
             "T2_at_2k_hm_wdm": float(np.interp(2 * first_crossing(k, t2_clump, 0.25), k[sel], t2_wdm[sel]))}

    summary = {"k_max_h_Mpc": K_MAX, "k_per_decade_for_pk": K_PER_DECADE, "Omega_x": omega_x, "h": h,
               "fit_log10_k_hm": {"slope": slope, "intercept": intercept, "max_residual_dex": fit_resid},
               "bounds_via_WDM_half_mode": bounds, "boehm_2014_direct": boehm,
               "boehm_2014_calibration": calibration, "cmb": cmb, "shape_check": shape, "rows": rows}
    with open(OUT / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=float)
    plots(k, t2s, rows, bounds, omega_x, h, shape)
    print(json.dumps({k_: v for k_, v in summary.items() if k_ != "rows"}, indent=1, default=float))


def plots(k, t2s, rows, bounds, omega_x, h, shape):
    figs = ROOT / "figures"
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    show = [-9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0]
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(show)))
    for lu, c in zip(show, colors):
        ax.semilogx(k, t2s[lu], color=c, label=f"log₁₀u = {lu:g}")
    m_match = shape["matched_m_WDM_keV"]
    ax.semilogx(k, wdm_t2(k, m_match, omega_x, h), "k--", lw=1,
                label=f"WDM {m_match:.1f} keV (same k_hm as log₁₀u = {shape['log10_u']:g})")
    ax.axhline(0.25, color="0.6", ls=":", lw=0.8)
    ax.set(xlabel="k [h/Mpc]", ylabel="T² = P/P_ΛCDM", ylim=(-0.02, 1.1), xlim=(0.01, 0.95 * K_MAX),
           title="M12: linear small-scale suppression from clump drag (f_cl = 1, z = 0)")
    ax.legend(fontsize=7)
    fig.savefig(figs / "m12_transfer.png", dpi=150)
    plt.close(fig)

    u = np.array([r["u"] for r in rows])
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    ax.loglog(u, [r["k_hm"] for r in rows], "o-", ms=3, label="k_hm (T = ½)")
    ax.loglog(u, [r["k_half"] for r in rows], "s-", ms=3, label="k_½ (T² = ½)")
    ax.loglog(u, [r["k_10"] for r in rows], "^-", ms=3, label="k_10 (T² = 0.9)")
    styles = iter(["-", "--", ":"])
    for label, b in bounds.items():
        ax.axhline(b["k_hm_limit_h_Mpc"], color="C3", ls=next(styles), lw=0.8, label=f"WDM k_hm: {label}")
    ax.axvline(BOEHM_2014_U_MAX, color="C4", lw=1, label="Boehm+2014 MW satellites (direct DM–γ): u < 5.5×10⁻⁷")
    ax.axvline(2.25e-4, color="C5", lw=1, label="Planck CMB (Stadler & Bœhm): u < 2.25×10⁻⁴")
    ax.set(xlabel="u", ylabel="k [h/Mpc]", title="M12: suppression scale vs coupling, and small-scale limits")
    top = ax.secondary_xaxis("top", functions=(lambda x: 268.0 / np.maximum(x, 1e-30), lambda x: 268.0 / np.maximum(x, 1e-30)))
    top.set_xlabel("Σ/Q [g/cm²]")
    ax.legend(fontsize=6.5, loc="lower left")
    fig.savefig(figs / "m12_bounds.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    {"run": run, "analyse": analyse}[sys.argv[1]]()
