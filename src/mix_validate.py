"""Direct objective evaluations that validate the mixed CDM + clump model set-up.

Evaluates -log(posterior) and the chi^2 of each likelihood at one fixed parameter point
(the PR4 CamSpec u = 0 best fit of the earlier profile) for a list of model configurations,
with one CLASS build per process:

  uv run python src/mix_validate.py stock   -> results/mix/validation/direct_stock.json
  uv run python src/mix_validate.py patched -> results/mix/validation/direct_patched.json
  uv run python src/mix_validate.py compare

Checks (see doc/experiment_log/2026-10-02_clump_fraction_kernels.md):
  1. patched Thomson (explicit coefficients) against stock CLASS;
  2. f_cl = 0 against LCDM (the same model: no idm species);
  3. several fractions at u = 0 against each other;
  4. the two kernels as u -> 0;
and records the loaded classy module for every evaluation.
"""

import json
import sys

from cosmology import ROOT
from likelihood_setup import class_build_identity, info, select_class_build

OUT = ROOT / "results" / "mix" / "validation"
LIKE = "camspec_npipe"
POINT = ROOT / "results" / "likelihood" / LIKE / "quadfit" / "u0.000e+00" / "result.json"

# (label, f_cl, u, kernel); kernel None = coefficients not passed (stock-compatible)
CONFIGS = [
    ("lcdm (f=0)", 0.0, 0.0, None),
    ("f=1 u=0", 1.0, 0.0, None),
    ("f=0.5 u=0", 0.5, 0.0, None),
    ("f=0.1 u=0", 0.1, 0.0, None),
    ("f=0.03 u=0", 0.03, 0.0, None),
    ("f=0.01 u=0", 0.01, 0.0, None),
    ("f=1 u=1e-4 default", 1.0, 1e-4, None),
    ("f=0.1 u=1e-3 default", 0.1, 1e-3, None),
    ("f=1 u=1e-9 default", 1.0, 1e-9, None),
]
PATCHED_ONLY = [
    ("f=1 u=1e-4 thomson", 1.0, 1e-4, "thomson"),
    ("f=1 u=1e-4 isotropic", 1.0, 1e-4, "isotropic"),
    ("f=0.1 u=1e-3 thomson", 0.1, 1e-3, "thomson"),
    ("f=0.1 u=1e-3 isotropic", 0.1, 1e-3, "isotropic"),
    ("f=1 u=1e-9 thomson", 1.0, 1e-9, "thomson"),
    ("f=1 u=1e-9 isotropic", 1.0, 1e-9, "isotropic"),
    ("f=0.01 u=1e-9 isotropic", 0.01, 1e-9, "isotropic"),
]


def evaluate(build):
    import sys as _sys
    from cobaya.model import get_model
    identity = class_build_identity(build)
    point = json.load(open(POINT))["best"]
    rows = []
    for label, f_cl, u, kernel in CONFIGS + (PATCHED_ONLY if build == "patched" else []):
        model = get_model(info(LIKE, u=u, f_cl=f_cl, kernel=kernel))
        lp = model.logposterior({n: point[n] for n in model.parameterization.sampled_params()})
        rows.append({"label": label, "f_cl": f_cl, "u": u, "kernel": kernel,
                     "minuslogpost": -float(lp.logpost),
                     "chi2": {n: -2 * float(v) for n, v in zip(model.likelihood, lp.loglikes)},
                     "derived": dict(zip(model.parameterization.derived_params(), map(float, lp.derived))),
                     "classy_path": _sys.modules["classy"].__file__,
                     "classy_extra_args": {k: v for k, v in model.theory["classy"].extra_args.items()
                                           if "idm" in k or k in ("non_linear",)}})
        print(f"{build:8s} {label:26s} -logP = {rows[-1]['minuslogpost']:.4f}", flush=True)
        model.close()
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / f"direct_{build}.json", "w") as fh:
        json.dump({"identity": identity, "point": point, "rows": rows}, fh, indent=2)


def compare():
    stock = {r["label"]: r for r in json.load(open(OUT / "direct_stock.json"))["rows"]}
    patched = {r["label"]: r for r in json.load(open(OUT / "direct_patched.json"))["rows"]}
    ref = stock["lcdm (f=0)"]["minuslogpost"]
    print(f"{'configuration':28s} {'stock':>12s} {'patched':>12s} {'patched-stock':>14s} {'patched-LCDM':>13s}")
    for label, p in patched.items():
        s = stock.get(label)
        print(f"{label:28s} {s['minuslogpost'] if s else float('nan'):12.4f} {p['minuslogpost']:12.4f} "
              f"{p['minuslogpost'] - s['minuslogpost'] if s else float('nan'):14.2e} {p['minuslogpost'] - ref:13.4f}")
    for a, b in (("f=1 u=1e-4 thomson", "f=1 u=1e-4 default"), ("f=0.1 u=1e-3 thomson", "f=0.1 u=1e-3 default"),
                 ("f=1 u=1e-4 isotropic", "f=1 u=1e-4 thomson"), ("f=0.1 u=1e-3 isotropic", "f=0.1 u=1e-3 thomson"),
                 ("f=1 u=1e-9 isotropic", "f=1 u=0"), ("f=1 u=1e-9 thomson", "f=1 u=0"),
                 ("f=0.01 u=1e-9 isotropic", "f=0.01 u=0")):
        print(f"patched: [{a}] - [{b}] = {patched[a]['minuslogpost'] - patched[b]['minuslogpost']:+.3e}")
    print("stock (default kernel) explicit-vs-stock: [patched thomson] - [stock default] at f=1, u=1e-4 = "
          f"{patched['f=1 u=1e-4 thomson']['minuslogpost'] - stock['f=1 u=1e-4 default']['minuslogpost']:+.3e}")
    print("classy modules:", {r["classy_path"] for r in stock.values()}, {r["classy_path"] for r in patched.values()})


if __name__ == "__main__":
    if sys.argv[1] == "compare":
        compare()
    else:
        select_class_build(sys.argv[1])
        evaluate(sys.argv[1])
