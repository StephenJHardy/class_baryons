"""Residual metrics between two sets of saved observables."""

import numpy as np

TOLERANCE = 1e-4  # plan Section 5, M2a


def residuals(model, ref, ell_min=2):
    """Fractional differences of model vs ref.

    TT, EE, phiphi and P(k): (X - X_ref) / X_ref.
    TE: (TE - TE_ref) / sqrt(TT_ref EE_ref), avoiding TE zero-crossings.
    """
    sel = ref["ell"] >= ell_min
    res = {"ell": ref["ell"][sel], "k_h": ref["k_h"]}
    for key in ("tt", "ee", "pp"):
        res[key] = (model[key][sel] - ref[key][sel]) / ref[key][sel]
    res["te"] = (model["te"][sel] - ref["te"][sel]) / np.sqrt(ref["tt"][sel] * ref["ee"][sel])
    res["pk"] = (model["pk"] - ref["pk"]) / ref["pk"]
    return res


def max_abs(res):
    return {key: float(np.max(np.abs(res[key]))) for key in ("tt", "te", "ee", "pp", "pk")}
