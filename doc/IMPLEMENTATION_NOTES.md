# Implementation notes

*A living guide for anyone continuing this project: what was learned the hard way, organised by topic. Dated evidence lives in [`experiment_log/`](experiment_log/).*

## Reproducing the environment

```bash
scripts/setup_class.sh                 # clones CLASS v3.4.0 into class_public/, checks the commit, runs uv sync
uv run python src/run_baseline.py      # M1
uv run python src/validate_m2a.py      # M2a
uv run python src/validate_m2b.py      # M2b
uv run python src/diagnose_m2a.py      # reproduces the M2a precision investigation
```

- classy is installed from the local `class_public/` checkout through a uv path source. `uv sync` builds it. The C compiler needs OpenMP; gcc 13 works.
- Each model (lensed C_ℓ to ℓ = 2500, P(k) to 20 h/Mpc) takes ~1.2 s on 16 cores with the project precision settings.

## Code layout

| File | Purpose |
|---|---|
| `configs/baseline.yaml` | Frozen cosmology, output settings, **precision settings** and idm defaults. Any change needs a log entry. |
| `src/cosmology.py` | `class_params(config, f_cl, u, m_idm, extra)` builds the CLASS dictionary. f_cl = 0 means no idm species at all. |
| `src/run_class.py` | `run()` computes observables; `save()` and `load()` handle `results/spectra/<name>.{npz,json}`. |
| `src/diagnostics.py` | Residual metrics (TE normalised by √(TT·EE)) and the 10⁻⁴ tolerance. |

Scripts in `src/` import each other as top-level modules, so run them as `uv run python src/<script>.py`.

## CLASS idm: parameters and pitfalls

- **Use `omega_idm`, not `f_idm`.** With `f_idm`, CLASS derives densities from the CDM budget, and `f_idm > 0` with `Omega_cdm == 0` is an error. Passing `omega_cdm = (1-f)ω_dm` and `omega_idm = f·ω_dm` is explicit.
- **The CDM floor.** In the synchronous gauge CLASS raises `Omega0_cdm` to at least `Omega0_cdm_min_synchronous = 1e-10`, because CDM defines the gauge. That is negligible for us.
- **`n_index_idm_g` is read only when `u_idm_g > 0`.** classy fails with "did not read input parameter(s)" if you pass it at u = 0. `class_params` handles this.
- **`u_idm_g` vs `cross_idm_g`.** CLASS converts between them using `cross = u · σ_T · m_idm / 1e11` (m in eV). Always pass `u_idm_g`, so that σ/m stays fixed when `m_idm` changes.
- **`m_idm` has a floor of 10⁶ eV** (CLASS refuses smaller values). In the range 10⁶–10⁷ eV the idm thermal sound speed visibly damps φφ and P(k).
- **`m_idm ≥ 1e24 eV` makes the mass drop out exactly** (bit-identical output), because the T/m terms fall below double precision. The project default is 10²⁷ eV. CLASS's own default is 10⁹ eV.
- **The tight-coupling approximation.** When `u_idm_g > 0`, CLASS forces `tight_coupling_approximation = compromise_CLASS`, which is already the default. With idm–baryon coupling it would force `first_order_CLASS`.
- **The thermodynamics start moves** to `thermo_z_initial_if_idm = 1e9` whenever idm is present and not tightly coupled early. With early idm–photon tight coupling it moves to 100× the idm decoupling redshift instead. Both rescale `thermo_Nz_log`, so **the thermodynamics sampling depends on u**. The project precision settings make this harmless; see below.
- **The idm initial velocity is θ_idm = θ_γ** (`perturbations.c` ~l.5929, adiabatic initial conditions), not 0 as for CDM. That is right for tightly coupled idm. For uncoupled or weakly coupled idm it seeds a decaying velocity mode that is independent of k, and at default start times this leaves ~5×10⁻⁴ residuals in EE and TE at ℓ < 30 and 10⁻⁴ in P(k) at k ~ 10⁻³ h/Mpc. **Starting the perturbations 10× earlier suppresses it** (to ~2×10⁻⁶). Any change to the initial conditions, as Experiment B may make, should revisit this.

## Precision: what is needed for 10⁻⁴

Default CLASS v3.4.0 precision is only good to ~10⁻³ in C_ℓ for *relative* comparisons between models whose numerical set-up differs. In pure ΛCDM, just moving the thermodynamics table start changes EE by 2×10⁻³. The project settings (in `configs/baseline.yaml`):

| Parameter | Default | Project | Why |
|---|---|---|---|
| `thermo_Nz_lin` | 20000 | 80000 | thermodynamics table resolution (the dominant error) |
| `thermo_Nz_log` | 5000 | 20000 | same |
| `start_small_k_at_tau_c_over_tau_h` | 1.5e-3 | 1.5e-4 | earlier start suppresses the θ_idm initial-condition artefact |
| `start_large_k_at_tau_h_over_tau_k` | 0.07 | 7e-3 | same |

- These cost no measurable runtime.
- Things that did **not** help: `tol_thermo_integration`, `thermo_integration_stepsize`, `thermo_rate_smoothing_radius`, `reionization_sampling`, finer low-k sampling (`k_step_super`, `k_min_tau0`) and `tol_perturbations_integration`.
- Starting 100× earlier gains nothing over 10× and adds P(k) noise at k ≳ 6 h/Mpc.
- **Measured noise floor** (from nudging `tol_perturbations_integration` by 10%): C_ℓ ~1–2×10⁻⁵; ΔP/P_ΛCDM ~2×10⁻⁵.
- **P(k) at high k.** Even ΛCDM has ~10⁻⁴ noise at k ≳ 5 h/Mpc. For strongly suppressed models, ΔP/P in the damped tail (T² < 10⁻³) is noise-dominated at ~10⁻³. Judge P(k) changes with ΔT² = ΔP/P_ΛCDM. The §6.3 diagnostics (k_10%, k_½) sit where T² ≳ 0.1 and are unaffected. Anything that needs P(k) deep in the tail needs more k precision.

## Other CLASS gotchas

- **θ_s is not Planck's θ_MC.** CLASS `100*theta_s` ≈ 1.0419 for Planck-2018 cosmology; 1.04110 (θ_MC) gives h = 0.671.
- **P(k) from CLASS uses δ_m in the comoving gauge** (gauge-invariant) by default. Transfer functions at z_pk > z_rec need `matter_source_in_current_gauge = yes` or `output = mTk` only.
- **Newtonian and synchronous gauges disagree at the 10⁻³ level** at default precision (P(k) up to 5×10⁻³). Stay in the synchronous gauge, CLASS's default, for consistency.
- **classy rejects unread parameters.** When reusing a parameter dictionary for a different output (e.g. `mTk` only), remove `lensing` and `l_max_scalars`.

## Dead ends

- **Newtonian gauge as a check on the θ_idm initial condition.** It made the idm–ΛCDM residual larger, not smaller: the Newtonian-gauge transformation keeps the θ_γ offset, and ΛCDM itself is less self-consistent across gauges. The decisive test was the earlier start time.

## Open questions

- How the φφ response scales with u (at u = 10⁻⁴ it is −54% at L = 2500): M3.
- Whether the §6.3 small-scale diagnostics need more k precision: check in M3.
