"""Convert the CLASS coupling u to sigma/M and Sigma/Q (plan Section 2).

u = (sigma / sigma_T) (M / 100 GeV)^-1, where sigma is the momentum-transfer
cross-section. Constants match CLASS v3.4.0 (thermodynamics.h, background.h)
so that the mapping is consistent with what CLASS actually integrates.
"""

SIGMA_T_CM2 = 6.6524616e-29 * 1e4          # CLASS _sigma_ (m^2 -> cm^2)
EV_J = 1.602176487e-19                     # CLASS _eV_
C_M_S = 2.99792458e8                       # CLASS _c_
GRAM_PER_100GEV = 1e11 * EV_J / C_M_S**2 * 1e3
U_PER_SIGMA_OVER_M = GRAM_PER_100GEV / SIGMA_T_CM2   # u per (cm^2/g), ~268


def sigma_over_m(u):
    """Momentum-transfer cross-section per unit mass [cm^2/g]."""
    return u / U_PER_SIGMA_OVER_M


def u_from_sigma_over_m(sigma_m):
    return sigma_m * U_PER_SIGMA_OVER_M


def sigma_surface_over_q(u):
    """Surface density over momentum-transfer efficiency, Sigma/Q [g/cm^2]."""
    return 1.0 / sigma_over_m(u)


def u_from_sigma_surface(sigma_surface, q=1.0):
    return q * U_PER_SIGMA_OVER_M / sigma_surface
