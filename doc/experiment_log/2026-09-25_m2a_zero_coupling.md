# 2026-09-25 — M2a: zero-coupling idm reproduces ΛCDM

**Goal:** Milestone M2a. With f_cl = 1 and u = 0 (all dark matter as `idm`, `omega_cdm = 0`), CLASS must reproduce `baseline_lcdm` to |ΔC_ℓ/C_ℓ| ≤ 10⁻⁴ for 2 ≤ ℓ ≤ 2500 (TE relative to √(C^TT C^EE)), with P(k) to the same tolerance for 10⁻⁴ ≤ k ≤ 10 h/Mpc.

## Configuration

CLASS v3.4.0 @ `64bbab7` and `configs/baseline.yaml`. The idm model is built by `class_params(config, f_cl=1, u=0)`.

## What was run

```bash
uv run python src/validate_m2a.py    # final validation -> results/validation/m2a.json, figures/m2a_residuals.png
uv run python src/diagnose_m2a.py    # reproduces the investigation -> results/validation/m2a_diagnosis.json
```

## Outcome

**At default CLASS precision M2a failed.** Max residuals: TT 5.2×10⁻⁴, TE 6.9×10⁻⁴, EE 1.9×10⁻³, φφ 7.4×10⁻⁴ and P(k) 1.0×10⁻⁴. The tolerance was not loosened. The failure was traced to two separate causes, neither of them physics:

### Cause 1: thermodynamics table resolution (the ℓ ≳ 30 part)

When idm is present, CLASS moves the start of the thermodynamics table from z = 5×10⁶ to z = 10⁹ (`thermo_z_initial_if_idm`). It also rescales `thermo_Nz_log`, so the table points land in different places. In **pure ΛCDM**, making the same change reproduces the whole discrepancy (TT 5.1×10⁻⁴, EE 2.1×10⁻³). So CLASS at default precision is only accurate to about 10⁻³ against trivial changes in sampling.

- `tol_thermo_integration`, `thermo_integration_stepsize` and `thermo_rate_smoothing_radius` barely help.
- Multiplying `thermo_Nz_lin` and `thermo_Nz_log` by 4 brings the ΛCDM self-consistency to ≤3.4×10⁻⁵, at no measurable runtime cost.
- Nz ×4 vs ×8: TT 1.2×10⁻⁵, EE 1.1×10⁻⁴, φφ 2.9×10⁻⁴. The ×4 table is adequate for comparisons made at matched precision.

### Cause 2: the idm initial velocity (the ℓ < 30 and P(k) part)

With Nz ×4, a residual remained (EE 5.5×10⁻⁴, TE 4.0×10⁻⁴ and TT 7×10⁻⁵, all at ℓ ≈ 9–11, plus P(k) 1.0×10⁻⁴ at k ≈ 1.3×10⁻³ h/Mpc). It:
- did not change with Nz, low-k sampling, reionization sampling or integrator tolerance;
- did not change with `m_idm` from 10⁶ to 10²⁰ eV, so it is not the idm sound speed;
- scaled linearly with f_cl (×0.1 at f_cl = 0.1);
- was not removed by the Newtonian gauge (it got larger there).

Comparing transfer functions showed θ_idm ≈ −1.3×10⁻⁹ at z = 900, nearly independent of k and decaying as 1/a.

The cause is in `source/perturbations.c:5929`: the adiabatic initial conditions set **θ_idm = θ_γ**, which is right for tightly coupled idm, while CDM gets θ = 0, which defines the synchronous gauge. For uncoupled idm this seeds a spurious decaying velocity mode. θ_γ,ini ∝ k⁴τ_ini³ and CLASS starts each mode at τ_ini ∝ 1/k, so after 1/τ decay the mode is independent of k. That matches what was measured.

**Test (without editing CLASS):** start the perturbations 10× earlier (`start_small_k_at_tau_c_over_tau_h` 1.5×10⁻³ → 1.5×10⁻⁴, `start_large_k_at_tau_h_over_tau_k` 0.07 → 0.007). The low-ℓ residual drops from 5.5×10⁻⁴ to 1.8×10⁻⁶. Starting 100× earlier gives no further gain, and makes P(k) at k ≳ 6 h/Mpc noisier (~10⁻⁴, which ΛCDM alone also shows there).

### Final result: PASS

With Nz ×4 and the 10× earlier start (both now in `configs/baseline.yaml`):

| TT | TE | EE | φφ | P(k) |
|---:|---:|---:|---:|---:|
| 1.9×10⁻⁵ | 1.2×10⁻⁵ | 8.3×10⁻⁶ | 1.9×10⁻⁶ | 4.1×10⁻⁵ |

The largest P(k) residual is at k ≈ 8 h/Mpc. Figure: `figures/m2a_residuals.png`. Runtime is unchanged (~1.2 s per model).

## Follow-ups

- The θ_idm = θ_γ initial condition also applies at small non-zero u, where idm is effectively uncoupled at early times. The earlier start deals with it there too. Worth remembering if Experiment B changes the initial conditions.
- P(k) at k ≳ 5 h/Mpc has ~10⁻⁴ numerical noise even in ΛCDM. The small-scale diagnostics of §6.3 may need higher k precision; check this in M3.
