"""Cobaya input for Planck 2018 likelihoods with the project's CLASS set-up.

The CLASS precision block and idm defaults come from configs/baseline.yaml, so
the likelihood runs use exactly the numerics validated in M2a-M5. The clump
model is f_cl = 1 by default (omega_cdm = 0, omega_idm = omega_dm), with
u_idm_g either fixed (profile points) or sampled/free.

Likelihood sets (Planck 2018; installed with cobaya-install into
/media/stephen/astro/class_baryons/cobaya_packages):
  "plik"      : planck_2018_lowl.TT + planck_2018_lowl.EE + planck_2018_highl_plik.TTTEEE
  "plik_lite" : planck_2018_lowl.TT + planck_2018_lowl.EE + planck_2018_highl_plik.TTTEEE_lite_native
"""

import os
import sys
from pathlib import Path

from cosmology import ROOT, load_config

# Scattering kernels of the clump-photon collision term (class_patches/idm_g_kernel.patch).
# Both are compared at equal momentum-transfer coupling u. "isotropic" is uniform re-emission
# in the clump rest frame with no polarization generation (the M6 absorber kernel).
KERNELS = {
    "thomson": {"idm_g_quadrupole_coefficient": 1.0, "idm_g_polarization_coefficient": 1.0},
    "isotropic": {"idm_g_quadrupole_coefficient": 0.0, "idm_g_polarization_coefficient": 0.0},
}
PATCHED_DIR = ROOT / "build" / "classy_patched"


def select_class_build(build):
    """Make `import classy` resolve to the stock or patched CLASS build, and verify it.

    Call before Cobaya creates a model (Cobaya's classy wrapper with path='global' imports
    whatever `classy` is first on sys.path). Returns the path of the loaded module. An
    earlier run that was meant to be patched silently imported stock CLASS; this raises instead.
    """
    if build not in ("stock", "patched"):
        raise ValueError(f"class build must be 'stock' or 'patched', not {build!r}")
    if build == "patched" and "classy" not in sys.modules:
        sys.path.insert(0, str(PATCHED_DIR))
    import classy
    path = Path(classy.__file__).resolve()
    patched = path.is_relative_to(PATCHED_DIR.resolve())
    if patched != (build == "patched"):
        raise ImportError(f"class build {build!r} requested but classy was loaded from {path}")
    return str(path)


def class_build_identity(build):
    """Version, module path and (for the patched build) the patch checksum, for the result files."""
    import hashlib
    import classy
    out = {"build": build, "classy_path": select_class_build(build), "classy_version": classy.__version__}
    if build == "patched":
        patch = ROOT / "class_patches" / "idm_g_kernel.patch"
        out["patch_sha256"] = hashlib.sha256(patch.read_bytes()).hexdigest()
    return out


# Cobaya packages: local astro disk by default; override with COBAYA_PACKAGES_PATH (e.g. on a cloud VM).
PACKAGES = os.environ.get("COBAYA_PACKAGES_PATH", "/media/stephen/astro/class_baryons/cobaya_packages")
LIKELIHOODS = {
    "plik": ["planck_2018_lowl.TT", "planck_2018_lowl.EE", "planck_2018_highl_plik.TTTEEE"],
    "plik_lite": ["planck_2018_lowl.TT", "planck_2018_lowl.EE", "planck_2018_highl_plik.TTTEEE_lite_native"],
    # Planck PR4 (NPIPE) high-l CamSpec (Rosenberg, Gratton & Efstathiou 2022) with the 2018 low-l likelihoods
    "camspec_npipe": ["planck_2018_lowl.TT", "planck_2018_lowl.EE", "planck_NPIPE_highl_CamSpec.TTTEEE"],
}
# Cosmological parameters: Planck-2018-like starting point, broad flat priors,
# proposal widths ~ Planck 2018 errors (used by the minimiser for scaling).
COSMO = {
    "omega_b": {"prior": {"min": 0.019, "max": 0.026}, "ref": 0.02237, "proposal": 0.00015, "latex": r"\omega_b"},
    "omega_dm": {"prior": {"min": 0.09, "max": 0.15}, "ref": 0.1200, "proposal": 0.0012, "latex": r"\omega_{\rm dm}"},
    "theta_s_100": {"prior": {"min": 1.03, "max": 1.05}, "ref": 1.0419, "proposal": 0.0003, "latex": r"100\theta_s"},
    "n_s": {"prior": {"min": 0.9, "max": 1.03}, "ref": 0.9649, "proposal": 0.004, "latex": "n_s"},
    "logA": {"prior": {"min": 2.8, "max": 3.3}, "ref": 3.044, "proposal": 0.014, "latex": r"\ln 10^{10}A_s", "drop": True},
    "tau_reio": {"prior": {"min": 0.01, "max": 0.12}, "ref": 0.0544, "proposal": 0.007, "latex": r"\tau"},
}


def info(likelihoods="plik", u=0.0, f_cl=1.0, u_free=False, extra_classy=None, sampler=None, output=None,
         fixed_params=None, kernel=None):
    """Cobaya input. With `kernel` ("thomson" or "isotropic") the kernel coefficients are passed
    explicitly, which needs the patched CLASS build (select_class_build("patched")); they are
    read by CLASS only when u > 0, so at u = 0 or f_cl = 0 they are omitted."""
    config = load_config()
    fixed = {k: v for k, v in config["cosmology"].items()
             if k not in ("omega_b", "omega_dm", "100*theta_s", "n_s", "ln10^{10}A_s", "tau_reio")}
    extra = dict(fixed)
    extra.update(config["precision"])
    extra.update({"non_linear": "none"})
    params = {k: dict(v) for k, v in COSMO.items()}
    # random starting points around the reference (used by the minimiser's restarts)
    for p in params.values():
        p["ref"] = {"dist": "norm", "loc": p["ref"], "scale": p["proposal"]}
    params["A_s"] = {"value": "lambda logA: 1e-10*np.exp(logA)", "latex": "A_s"}
    params["100*theta_s"] = {"value": "lambda theta_s_100: theta_s_100", "derived": False}
    params["omega_cdm"] = {"value": f"lambda omega_dm: {1.0 - f_cl!r}*omega_dm", "derived": False}
    if f_cl > 0:
        params["omega_idm"] = {"value": f"lambda omega_dm: {f_cl!r}*omega_dm", "derived": False}
        extra["m_idm"] = config["idm"]["m_idm"]
        if u_free:
            params["u_idm_g"] = {"prior": {"min": 0.0, "max": 1e-3}, "ref": 1e-5, "proposal": 3e-5, "latex": "u"}
            extra["n_index_idm_g"] = config["idm"]["n_index_idm_g"]
        else:
            extra["u_idm_g"] = u
            if u > 0:
                extra["n_index_idm_g"] = config["idm"]["n_index_idm_g"]
                if kernel is not None:
                    extra.update(KERNELS[kernel])
    params["omega_dm"]["drop"] = True
    params["theta_s_100"]["drop"] = True
    params["H0"] = {"latex": "H_0"}
    params["sigma8"] = {"latex": r"\sigma_8"}
    if extra_classy:
        extra.update(extra_classy)
    for name, value in (fixed_params or {}).items():   # e.g. pin boundary-hugging nuisance parameters
        params[name] = value
    out = {
        "theory": {"classy": {"extra_args": extra, "path": "global"}},
        "likelihood": {name: None for name in LIKELIHOODS[likelihoods]},
        "params": params,
        "packages_path": PACKAGES,
    }
    if sampler:
        out["sampler"] = sampler
    if output:
        out["output"] = output
    return out
