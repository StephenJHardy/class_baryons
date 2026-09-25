# 2026-09-25 — M2b: cold-mass limit in m_idm

**Goal:** Milestone M2b. At fixed u = 10⁻⁴ (f_cl = 1), show that outputs converge as `m_idm` increases, then choose a project default inside the converged range.

## Configuration

- CLASS v3.4.0 @ `64bbab7`, `configs/baseline.yaml` with the M2a precision settings, and `n_index_idm_g = 0`.
- u is passed as `u_idm_g`, so σ/m stays fixed while m varies.
- Masses: 10⁶ (CLASS's minimum), 10⁷–10¹², 10¹⁵, 10¹⁸, 10²¹, 10²⁴, 10²⁷ and 10³⁰ eV.

## What was run

```bash
uv run python src/validate_m2b.py
```

Outputs: `results/spectra/m2b_u1e-4_m*.{npz,json}`, `results/validation/m2b.json`, `figures/m2b_mass_convergence.png` and `figures/m2b_residuals.png`.

## Outcome: PASS

- No failures or numerical warnings at any mass. Runtime ~1.2 s per model.
- **Convergence metric for P(k):** ΔT² ≡ (P − P_ref)/P_ΛCDM, not ΔP/P. At u = 10⁻⁴, P(k) is suppressed by up to 10⁵ at k ≈ 10 h/Mpc. In that tail, ΔP/P is dominated by integrator noise: at fixed m = 10³⁰ eV, nudging `tol_perturbations_integration` from 10⁻⁵ to 0.9×10⁻⁵ changes ΔP/P by 1.3×10⁻³. ΔP/P is still recorded in `m2b.json`. The change of metric is deliberate and stated here, not a quiet loosening.
- **Noise floor** at fixed mass (tolerance nudge): TT 1.4×10⁻⁵, TE 1.3×10⁻⁵, EE 9.5×10⁻⁶, φφ 6.5×10⁻⁷, ΔT² 2.1×10⁻⁵.
- **Residuals against m = 10³⁰ eV:**

| m_idm [eV] | TT | EE | φφ | ΔT² |
|---:|---:|---:|---:|---:|
| 10⁶ | 1.1×10⁻⁵ | 9.2×10⁻⁶ | **7.6×10⁻⁴** | **8.9×10⁻⁴** |
| 10⁷ | 6.1×10⁻⁶ | 3.6×10⁻⁶ | 7.7×10⁻⁵ | 7.1×10⁻⁵ |
| 10⁸ – 10¹⁸ | ≤1.2×10⁻⁵ | ≤2.7×10⁻⁶ | ≤8.3×10⁻⁶ | ≤2.9×10⁻⁵ (noise floor) |
| 10²¹ | 4.5×10⁻⁸ | 2.9×10⁻⁹ | 1.8×10⁻⁶ | 2.2×10⁻⁵ |
| ≥10²⁴ | 0 (bit-identical) | 0 | 0 | 0 |

- A genuine mass dependence is visible only below ~10⁸ eV, in φφ and P(k) at k ~ 2 h/Mpc: the thermal sound speed of light idm damping small scales. TT/TE/EE show no mass dependence above the noise floor anywhere in the range.
- **From 10²⁴ eV the output is bit-identical.** The terms proportional to T/m fall below double-precision resolution, so the mass drops out of CLASS exactly. The same was checked at u = 10⁻⁶ and u = 10⁻² (`plateau_check` in `m2b.json`).
- **Decision: project default `m_idm = 1e27 eV`.** It sits inside the bit-identical plateau, which is the physically correct limit for macroscopic clumps (their thermal velocity is effectively zero). It is recorded in `configs/baseline.yaml`.
- **Incidental observation (not an M3 result):** at u = 10⁻⁴, f_cl = 1 and fixed cosmology, the maximum deviations from ΛCDM are TT 1.2%, EE 1.0%, TE 0.6% (of √(TT·EE)) and φφ 54%. P(k) is almost completely suppressed at k ≈ 10 h/Mpc. The φφ deviation grows steadily with L: −0.4% at L = 100, −2.3% at 300, −16% at 1000 and −54% at 2500. This is consistent with lensing at high L seeing the suppressed small-scale P(k).

## Follow-ups

- M3 scan, with all models at `m_idm = 1e27 eV`.
- In M3, quantify how the φφ response depends on u (it may be a more sensitive probe than TT at high L).
