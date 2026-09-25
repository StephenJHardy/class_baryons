"""M5: Gaussian noise models, detectability and Fisher matrices for TT/TE/EE (+ phiphi).

Noise spectra (white noise with a Gaussian beam, optionally 1/f):
    N_l = w^-1 exp[l(l+1) theta^2 / 8 ln 2] [1 + (l/l_knee)^alpha],   w^-1 = (Delta_T [uK rad])^2
Channels of one experiment are combined by inverse-variance weighting of N_l.
All C_l are in CLASS units (dimensionless, (dT/T)^2), so N_l is divided by T_cmb^2.

Sources (each checked against the paper itself; see the M5 experiment-log entry):
  Planck channels: Planck 2018 I (A&A 641, A1, 2020), Table 4:
      100/143/217 GHz beams 9.66/7.22/4.90 arcmin;
      temperature noise 1.29/0.55/0.78 uK deg; polarization noise 1.96/1.17/1.75 uK deg.
  Planck-like sky fractions and multipole ranges: the Planck forecast set-up of the
      Simons Observatory science-goals paper (Ade et al. 2019, JCAP 02, 056, Sec. on
      primordial parameters): T/E at 2 <= l <= 29 with f_sky = 0.8; TT/TE/EE at
      30 <= l <= 2500 with f_sky = 0.6; kappa-kappa at 8 <= L <= 400 with f_sky = 0.6.
  SO LAT: Ade et al. 2019, Table 1 (baseline): 93 GHz 2.2 arcmin 8.0 uK-arcmin,
      145 GHz 1.4 arcmin 10 uK-arcmin, f_sky = 0.4, polarization noise sqrt(2) higher;
      1/f model (their Eq. 1): polarization l_knee = 700, alpha = -1.4, N_red = N_white;
      temperature l_knee = 1000, alpha = -3.5 (N_red is quoted in uK^2 s and needs
      survey details not tabulated, so we APPROXIMATE N_red = N_white at l_knee).
"""

from dataclasses import dataclass, field

import numpy as np

T_CMB_UK = 2.7255e6
ARCMIN = np.pi / 180 / 60


def white_noise(delta_uk_arcmin, fwhm_arcmin, ell, l_knee=None, alpha=None):
    """N_l in CLASS units for one channel."""
    theta = fwhm_arcmin * ARCMIN
    n = (delta_uk_arcmin * ARCMIN) ** 2 * np.exp(ell * (ell + 1) * theta**2 / (8 * np.log(2)))
    if l_knee is not None:
        with np.errstate(divide="ignore"):
            n = n * (1 + (np.maximum(ell, 1) / l_knee) ** alpha)
    return n / T_CMB_UK**2


def combine(*noises):
    """Inverse-variance combination of channel noise spectra."""
    return 1.0 / np.sum([1.0 / n for n in noises], axis=0)


@dataclass
class Patch:
    """A sky region/multipole range with its own noise and sky fraction."""
    name: str
    f_sky: float
    ell_min: int
    ell_max: int
    n_tt: np.ndarray            # noise arrays indexed by ell (length >= ell_max + 1)
    n_ee: np.ndarray
    spectra: tuple = ("tt", "te", "ee")


@dataclass
class Config:
    name: str
    patches: list
    lensing: dict = field(default_factory=dict)   # {"f_sky", "L_min", "L_max", "n_pp" or None}
    notes: str = ""


def planck_channels(ell, one_over_f=False):
    beams = {100: 9.66, 143: 7.22, 217: 4.90}
    t_ukdeg = {100: 1.29, 143: 0.55, 217: 0.78}
    p_ukdeg = {100: 1.96, 143: 1.17, 217: 1.75}
    n_tt = combine(*[white_noise(60 * t_ukdeg[c], beams[c], ell) for c in beams])
    n_ee = combine(*[white_noise(60 * p_ukdeg[c], beams[c], ell) for c in beams])
    return n_tt, n_ee


def so_lat_channels(ell, one_over_f=True):
    beams = {93: 2.2, 145: 1.4}
    t_ukarcmin = {93: 8.0, 145: 10.0}
    kw_t = {"l_knee": 1000, "alpha": -3.5} if one_over_f else {}
    kw_p = {"l_knee": 700, "alpha": -1.4} if one_over_f else {}
    n_tt = combine(*[white_noise(t_ukarcmin[c], beams[c], ell, **kw_t) for c in beams])
    n_ee = combine(*[white_noise(np.sqrt(2) * t_ukarcmin[c], beams[c], ell, **kw_p) for c in beams])
    return n_tt, n_ee


def configurations(ell_max=3000):
    ell = np.arange(ell_max + 1, dtype=float)
    zero = np.zeros_like(ell)
    p_tt, p_ee = planck_channels(ell)
    s_tt, s_ee = so_lat_channels(ell, one_over_f=True)
    s_tt_w, s_ee_w = so_lat_channels(ell, one_over_f=False)
    planck_patches = [Patch("Planck low-l", 0.8, 2, 29, p_tt, p_ee),
                      Patch("Planck high-l", 0.6, 30, 2500, p_tt, p_ee)]

    def so_plus_planck(n_tt, n_ee, label):
        # SO LAT sky (0.4) is inside Planck's 0.6: combine noise there; Planck alone on the other 0.2.
        return [Patch("Planck low-l", 0.8, 2, 29, p_tt, p_ee),
                Patch(f"{label} + Planck (SO sky)", 0.4, 30, 3000, combine(n_tt, p_tt), combine(n_ee, p_ee)),
                Patch("Planck only (rest of Planck sky)", 0.2, 30, 2500, p_tt, p_ee)]

    return {
        "cv_fullsky": Config("cosmic variance, full sky, 2-2500", [Patch("CV", 1.0, 2, 2500, zero, zero)],
                             notes="reproduces the M3 statistic when cross-covariance is off"),
        "planck_like": Config("Planck-like", planck_patches,
                              lensing={"f_sky": 0.6, "L_min": 8, "L_max": 400, "n_pp": None},
                              notes="channels 100+143+217 GHz; SO-paper Planck set-up"),
        "so_baseline": Config("SO LAT baseline + Planck", so_plus_planck(s_tt, s_ee, "SO LAT"),
                              notes="SO LAT 93+145 GHz baseline with 1/f (T knee approximated)"),
        "so_baseline_white": Config("SO LAT baseline (white noise only) + Planck",
                                    so_plus_planck(s_tt_w, s_ee_w, "SO LAT white"),
                                    notes="same, without atmospheric 1/f noise"),
    }


def cov2(cl, patch, ell):
    """2x2 total covariance [[TT+N, TE], [TE, EE+N]] at each ell (shape (n, 2, 2))."""
    c = np.empty((len(ell), 2, 2))
    c[:, 0, 0] = cl["tt"][ell] + patch.n_tt[ell]
    c[:, 1, 1] = cl["ee"][ell] + patch.n_ee[ell]
    c[:, 0, 1] = c[:, 1, 0] = cl["te"][ell]
    return c


def delta2(d, ell):
    m = np.empty((len(ell), 2, 2))
    m[:, 0, 0] = d["tt"][ell]
    m[:, 1, 1] = d["ee"][ell]
    m[:, 0, 1] = m[:, 1, 0] = d["te"][ell]
    return m


def fisher_block(cl, derivs, patch, which="joint"):
    """Gaussian Fisher matrix from one patch.

    which = 'joint' (TT+TE+EE with full covariance), or 'tt', 'ee', 'te' alone.
    derivs: list of dicts of dC_l/dp (same keys as cl). With derivs = [cl_model - cl_ref]
    this gives Delta chi^2 for a small difference.
    """
    ell = np.arange(patch.ell_min, patch.ell_max + 1)
    pref = patch.f_sky * (2 * ell + 1)
    n = len(derivs)
    f = np.zeros((n, n))
    if which == "joint":
        cinv = np.linalg.inv(cov2(cl, patch, ell))
        m = [cinv @ delta2(d, ell) for d in derivs]
        for i in range(n):
            for j in range(i, n):
                tr = np.einsum("lab,lba->l", m[i], m[j])
                f[i, j] = f[j, i] = np.sum(0.5 * pref * tr)
        return f
    if which in ("tt", "ee"):
        var = 2 * (cl[which][ell] + (patch.n_tt if which == "tt" else patch.n_ee)[ell]) ** 2 / pref
    else:  # te
        var = (cl["te"][ell] ** 2 + (cl["tt"][ell] + patch.n_tt[ell]) * (cl["ee"][ell] + patch.n_ee[ell])) / pref
    for i in range(n):
        for j in range(i, n):
            f[i, j] = f[j, i] = np.sum(derivs[i][which][ell] * derivs[j][which][ell] / var)
    return f


def lensing_fisher(cl, derivs, lensing):
    """phiphi block. With n_pp None this is cosmic-variance only: an upper envelope."""
    if not lensing:
        return None
    ell = np.arange(lensing["L_min"], lensing["L_max"] + 1)
    n_pp = 0.0 if lensing["n_pp"] is None else lensing["n_pp"][ell]
    var = 2 * (cl["pp"][ell] + n_pp) ** 2 / (lensing["f_sky"] * (2 * ell + 1))
    n = len(derivs)
    f = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            f[i, j] = f[j, i] = np.sum(derivs[i]["pp"][ell] * derivs[j]["pp"][ell] / var)
    return f


def fisher(cl, derivs, config, which="joint", lensing=False):
    f = sum(fisher_block(cl, derivs, p, which) for p in config.patches)
    if lensing:
        lf = lensing_fisher(cl, derivs, config.lensing)
        if lf is not None:
            f = f + lf
    return f


def delta_chi2(model, ref, config, which="joint", lensing=False):
    d = {k: model[k] - ref[k] for k in ("tt", "te", "ee", "pp")}
    return float(fisher(ref, [d], config, which, lensing)[0, 0])
