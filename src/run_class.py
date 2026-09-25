"""Run CLASS for one parameter dictionary and save its observables."""

import json
import time
from pathlib import Path

import numpy as np
from classy import Class

from cosmology import ROOT

SPECTRA_DIR = ROOT / "results" / "spectra"


def run(params, config):
    """Compute lensed CMB spectra, lensing potential and P(k, z=0).

    Returns a dict of numpy arrays. C_l are dimensionless (CLASS units);
    P(k) is in (Mpc/h)^3 on a fixed k grid in h/Mpc.
    """
    cosmo = Class()
    cosmo.set(params)
    t0 = time.time()
    cosmo.compute()
    runtime = time.time() - t0
    try:
        ell_max = config["ell_max"]
        cls = cosmo.lensed_cl(ell_max)
        raw = cosmo.raw_cl(ell_max)
        h = cosmo.h()
        k_h = np.logspace(np.log10(config["k_min_h"]), np.log10(config["k_max_h"]), config["n_k"])
        pk = np.array([cosmo.pk_lin(k * h, 0.0) for k in k_h]) * h**3
        derived = cosmo.get_current_derived_parameters(
            ["h", "100*theta_s", "Omega_m", "sigma8", "z_rec", "rs_rec", "age"])
        out = {
            "ell": cls["ell"],
            "tt": cls["tt"], "te": cls["te"], "ee": cls["ee"], "bb": cls["bb"],
            "pp": raw["pp"],
            "k_h": k_h, "pk": pk,
        }
    finally:
        cosmo.struct_cleanup()
        cosmo.empty()
    meta = {"params": params, "derived": derived, "runtime_s": runtime}
    return out, meta


def save(name, out, meta, directory=SPECTRA_DIR):
    directory.mkdir(parents=True, exist_ok=True)
    np.savez(directory / f"{name}.npz", **out)
    with open(directory / f"{name}.json", "w") as f:
        json.dump(meta, f, indent=2, default=float)


def load(name, directory=SPECTRA_DIR):
    with np.load(directory / f"{name}.npz") as data:
        out = dict(data)
    with open(directory / f"{name}.json") as f:
        meta = json.load(f)
    return out, meta
