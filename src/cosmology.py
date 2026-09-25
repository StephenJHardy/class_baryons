"""Build CLASS parameter dictionaries for the reference and clump models."""

from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "configs" / "baseline.yaml"


def load_config(path=DEFAULT_CONFIG):
    with open(path) as f:
        return yaml.safe_load(f)


def class_params(config, f_cl=0.0, u=0.0, m_idm=None, extra=None):
    """CLASS input for a model with a fraction f_cl of omega_dm in idm.

    f_cl = 0 gives pure LCDM with no idm species at all. The idm density is
    passed as omega_idm (not f_idm) so the bookkeeping stays explicit.
    """
    cosmo = dict(config["cosmology"])
    omega_dm = cosmo.pop("omega_dm")
    params = {**cosmo, **config["output"], **config["precision"]}
    params["omega_cdm"] = (1.0 - f_cl) * omega_dm
    if f_cl > 0:
        params["omega_idm"] = f_cl * omega_dm
        params["m_idm"] = config["idm"]["m_idm"] if m_idm is None else m_idm
        params["u_idm_g"] = u
        # CLASS only reads n_index_idm_g when u_idm_g > 0; passing it
        # otherwise makes classy fail with an unread-parameter error.
        if u > 0:
            params["n_index_idm_g"] = config["idm"]["n_index_idm_g"]
    if extra:
        params.update(extra)
    return params


def fix_h(params, h):
    """Replace the theta_s constraint by a fixed h (for fixed-H0 comparisons)."""
    params = {k: v for k, v in params.items() if k != "100*theta_s"}
    params["h"] = h
    return params
