"""Small-scale suppression diagnostics from P(k) (plan Section 6.3).

T^2(k) = P(k) / P_LCDM(k). Definitions used here:
  k_10:  first k where T^2 <= 0.9
  k_1/2: first k where T^2 <= 0.5
  k_hm:  first k where T <= 1/2 (T^2 <= 0.25): the conventional half-mode
  M_hm:  (4 pi / 3) rho_m (lambda_hm / 2)^3 with lambda_hm = 2 pi / k_hm [M_sun/h]
All k in h/Mpc. NaN means the threshold is not reached below k_max.
"""

import numpy as np

RHO_CRIT_MSUN_H2_MPC3 = 2.77536627e11


def first_crossing(k, t2, level):
    below = np.nonzero(t2 <= level)[0]
    if len(below) == 0:
        return np.nan
    i = below[0]
    if i == 0:
        return float(k[0])
    return float(np.exp(np.interp(level, [t2[i], t2[i - 1]], np.log([k[i], k[i - 1]]))))


def suppression(k_h, pk, pk_lcdm, omega_m):
    t2 = pk / pk_lcdm
    k_hm = first_crossing(k_h, t2, 0.25)
    m_hm = 4 * np.pi / 3 * omega_m * RHO_CRIT_MSUN_H2_MPC3 * (np.pi / k_hm) ** 3
    return {
        "k_10": first_crossing(k_h, t2, 0.9),
        "k_half": first_crossing(k_h, t2, 0.5),
        "k_hm": k_hm,
        "M_hm": float(m_hm),
    }
