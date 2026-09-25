# 2026-09-25 — M0: environment

**Goal:** Milestone M0 (plan §15). Set up a reproducible environment with CLASS ≥ v3.2 and confirm that the idm–photon parameters are accepted.

## Configuration

| Item | Value |
|---|---|
| CLASS | `v3.4.0` (latest release tag on 2026-09-25), commit `64bbab707faf4de4779a9e04edd180fef18d98fa` |
| CLASS location | `class_public/` (gitignored), cloned by `scripts/setup_class.sh`, which checks the commit |
| Python | 3.12.3, managed by uv 0.12.17 (`pyproject.toml`, `uv.lock`) |
| Key packages | classy v3.4.0 (built from source), numpy 2.5.3, scipy 1.18.1, matplotlib 3.11.2, PyYAML 6.0.3 |
| Compiler | gcc 13.3.0 (Ubuntu 24.04) |
| Machine | 16 cores, 62 GB RAM |

## What was run

```bash
scripts/setup_class.sh     # clone at the pinned tag, check the commit, uv sync
```

Then a smoke test through classy: a default ΛCDM run, plus an idm run with `omega_cdm=0, omega_idm=0.12, m_idm=1e9, u_idm_g=1e-4, n_index_idm_g=0`.

## Outcome

- classy builds in under 30 s. A lensed-C_ℓ + P(k) run takes ~0.4 s at default precision.
- CLASS accepts all four parameters named in the plan: `omega_idm`, `m_idm`, `u_idm_g` and `n_index_idm_g`.
- Things found while reading the CLASS source (details in `IMPLEMENTATION_NOTES.md`):
  - `n_index_idm_g` is read only when `u_idm_g > 0`.
  - `m_idm < 1e6 eV` is rejected.
  - When idm is present, the thermodynamics start is moved to z = 10⁹.

## Follow-ups

- M1: freeze the baseline.
