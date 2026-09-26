"""Write a Cobaya covmat for the standard parameters from the M5 Fisher matrix.

The minimiser needs the posterior's scale and degeneracy directions in the
project's own parameterization (omega_dm, theta_s_100, logA), which Planck's
covmats (omega_cdm, theta_MC_100) do not match. For profile points u is fixed,
so the relevant covariance is the inverse of the Fisher block without u
(the conditional covariance of the six standard parameters). A_planck is added
with its prior width 0.0025.

Usage: uv run python src/fisher_covmat.py
"""

import numpy as np

from cosmology import ROOT
from m5_detectability import FISHER_NAMES, derivatives
from m5_noise import configurations, fisher

COBAYA_NAMES = {"omega_b": "omega_b", "omega_dm": "omega_dm", "100*theta_s": "theta_s_100",
                "n_s": "n_s", "ln10^{10}A_s": "logA", "tau_reio": "tau_reio"}
OUT = ROOT / "results" / "likelihood" / "fisher_planck_like_u_fixed.covmat"


def main():
    fid, ders = derivatives()
    f = fisher(fid, ders, configurations()["planck_like"])
    cov = np.linalg.inv(f[1:, 1:])                     # drop u: conditional on fixed u
    names = [COBAYA_NAMES[n] for n in FISHER_NAMES[1:]] + ["A_planck"]
    full = np.zeros((len(names), len(names)))
    full[:-1, :-1] = cov
    full[-1, -1] = 0.0025**2
    OUT.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(OUT, full, header=" ".join(names), comments="# ")
    print(OUT, "\n", dict(zip(names, np.sqrt(np.diag(full)))))


if __name__ == "__main__":
    main()
