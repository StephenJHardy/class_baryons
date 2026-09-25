# Implementation notes

*A living guide for anyone continuing this project: what was learned the hard way, organised by topic. Dated evidence lives in [`experiment_log/`](experiment_log/).*

## Reproducing the environment

```bash
scripts/setup_class.sh                 # clones CLASS v3.4.0 into class_public/, checks the commit, runs uv sync
uv run python src/run_baseline.py      # M1
uv run python src/validate_m2a.py      # M2a
uv run python src/validate_m2b.py      # M2b
uv run python src/diagnose_m2a.py      # reproduces the M2a precision investigation
uv run python src/m3_coupling_scan.py  # M3 grid (~70 s, parallel)
uv run python src/m3_analyse.py        # M3/M4 tables and figures
uv run python src/glitch_scan.py "project" "l_linstep 20"   # l-node glitch check (~2 min per setting)
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
| `src/cloud_mapping.py` | u ↔ σ/M ↔ Σ/Q, using CLASS's own constants (u ≈ 268 × σ/M in cm²/g). |
| `src/rates.py` | Γ_γ→cl/H and Γ_cl→γ/H from CLASS output, decoupling redshifts, clump optical depth. Rates are resampled to 2000 z points for storage. |
| `src/small_scale.py` | k_10, k_½, k_hm and M_hm from T²(k) = P/P_ΛCDM. |
| `src/scan.py` | `run_grid(jobs, ...)` runs models in parallel processes and records failures without stopping. |
| `src/cosmology.py: fix_h` | Swaps the θ_s constraint for a fixed h. |

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

## Patched CLASS (Experiment B kernel)

- **Files.** `class_patches/idm_g_kernel.patch` is the versioned patch. It is authored by `scripts/make_kernel_patch_edits.py`: rerun that script on a clean `class_public` copy and regenerate the patch with `git diff`. `scripts/setup_class.sh` applies the patch to `class_patched/` and installs it to `build/classy_patched/` (both gitignored).
- **Selecting the build.** Set `CLASS_BUILD=patched` in the environment. It is inherited by `scan.run_grid` worker processes. The default is stock.
- **Import rule.** Never `import classy` or `from classy import Class` anywhere except `src/run_class.py`. Import `Class` or `classy` *from* `run_class`. Once `classy` has been imported, `sys.path` changes have no effect: an early import silently gave stock CLASS in the first M6 run. `run_class` now refuses to load if `classy` is already imported, and checks the module path against `CLASS_BUILD`. Each saved run records `class_build` and `classy_path` in its JSON; check them.
- **New inputs** (read only when `u_idm_g > 0`): `idm_g_quadrupole_coefficient` (c_T) and `idm_g_polarization_coefficient` (c_E). Both default to 1 (Thomson); set both to 0 for a perfect absorber with isotropic re-emission. For what each modifies (hierarchy, line-of-sight sources, tight-coupling closure), see the M6 log and the patch comments.
- **Keeping the Thomson path bit-identical.** Write each modified term as the stock expression plus (regen − rate)·(source), set regen exactly equal to rate when a coefficient is 1, and call the closure helper only when a coefficient ≠ 1. Rearranging terms algebraically (−rate·x + regen·s) changed results by ~10⁻⁵ through rounding alone.
- **Not patched:** vector and tensor modes, and the 11/6 second-order tight-coupling correction (it keeps its Thomson form; the clump share of the rate is ≲10⁻⁴ while tight coupling is on).

## Interaction rates in CLASS

- `get_thermodynamics()` has columns `dmu_idm_g` (the conformal photon opacity from clumps, a·n·σ·c, in 1/Mpc), `ddmu_idm_g`, `T_idm [K]` and `c_idm^2`. The table is ordered by **increasing z**. It is easy to assume the opposite, and `np.interp` fails silently if you reverse it.
- Physical rates: Γ_γ→cl = dmu·(1+z), and Γ_cl→γ = S·dmu·(1+z) with S = 4ρ_γ/3ρ_idm (the `S_idm_g` of `perturbations.c`). Compare them with `H [1/Mpc]` from `get_background()`.
- `dmu_idm_g` matches the analytic (1+z)² ρ_idm,0 (σ/M) to 4×10⁻⁶ (constants in `rates.py`, `cloud_mapping.py`).
- **Clump opacity is included in the visibility function**: `g` uses dκ + dμ_idm and exp(−κ−μ_idm). So strong coupling moves z_rec.
- **With θ_s fixed, h drifts with u**, because θ_s is evaluated at the z_rec that clump opacity shifts. The drift is negligible for u ≤ 10⁻³ (3×10⁻⁵ in h), 0.2% at 10⁻², and 2.4% at 10⁻¹. Physical diagnostics hold h fixed instead.

## CLASS limits at strong coupling (f_cl = 1)

- Runs succeed up to u = 10^−0.875 ≈ 0.13.
- From u ≈ 0.18 they fail with `perturbations_vector_init: scalar initial conditions assume tight-coupling approximation`: the idm–photon trigger switches the approximation off before the initial time. This is independent of the project's start-time setting.
- At u ≈ 1 thermodynamics fails because clump opacity puts the visibility peak at z ≈ 470, below the hard-coded `_Z_REC_MIN_ = 500`.
- Going further would need changes to the CLASS source.

## Precision: what is needed for 10⁻⁴

Default CLASS v3.4.0 precision is only good to ~10⁻³ in C_ℓ for *relative* comparisons between models whose numerical set-up differs. In pure ΛCDM, just moving the thermodynamics table start changes EE by 2×10⁻³. The project settings (in `configs/baseline.yaml`):

| Parameter | Default | Project | Why |
|---|---|---|---|
| `thermo_Nz_lin` | 20000 | 80000 | thermodynamics table resolution (the dominant error) |
| `thermo_Nz_log` | 5000 | 20000 | same |
| `start_small_k_at_tau_c_over_tau_h` | 1.5e-3 | 1.5e-4 | earlier start suppresses the θ_idm initial-condition artefact |
| `start_large_k_at_tau_h_over_tau_k` | 0.07 | 7e-3 | same |
| `l_linstep` | 40 | 20 | removes sporadic ℓ-node glitches (below) |

- These cost no measurable runtime.
- Things that did **not** help: `tol_thermo_integration`, `thermo_integration_stepsize`, `thermo_rate_smoothing_radius`, `reionization_sampling`, finer low-k sampling (`k_step_super`, `k_min_tau0`) and `tol_perturbations_integration`.
- Starting 100× earlier gains nothing over 10× and adds P(k) noise at k ≳ 6 h/Mpc.
- **Measured noise floor** (from nudging `tol_perturbations_integration` by 10%): C_ℓ ~1–2×10⁻⁵; ΔP/P_ΛCDM ~2×10⁻⁵.
- **P(k) at high k.** Even ΛCDM has ~10⁻⁴ noise at k ≳ 5 h/Mpc. For strongly suppressed models, ΔP/P in the damped tail (T² < 10⁻³) is noise-dominated at ~10⁻³. Judge P(k) changes with ΔT² = ΔP/P_ΛCDM. The §6.3 diagnostics (k_10%, k_½) sit where T² ≳ 0.1 and are unaffected. Anything that needs P(k) deep in the tail needs more k precision.

### Sporadic ℓ-node glitches (found in M3)

- **The symptom.** With the default `l_linstep = 40`, about a third of idm–photon models, apparently at random in u, carry a deterministic error of ~1.5×10⁻³ in unlensed EE and ~5×10⁻⁴ in TT. It is largest at ℓ = 836 (sometimes 1111) and rings with period ~40 in ℓ, while neighbouring u values are clean. The glitchy models agree with each other, so the glitch is a distinct numerical state.
- **What doesn't matter:** integrator tolerance, k sampling, source time sampling, tight-coupling and radiation-streaming settings, the perturbation start time, the log-region thermodynamics table, the late-source cuts and Bessel sampling.
- **What does:** the linear thermodynamics table size, which only moves the glitch to other u, and above all the ℓ-node spacing.
- **The fix.** `l_linstep = 20` reduces glitchy models from 41/59 to 3/59; the three left are TE at ℓ = 2, which is irrelevant.
- **Status.** The root cause is **not understood**. `src/glitch_scan.py` measures it: the second difference in u of the residual curves, over a fine grid, flags glitches. Rerun it after any precision or CLASS-version change.
- **Lesson:** always scan u finely enough to see non-monotonic behaviour. A coarse grid mistook this glitch for an early physical onset (EE 0.1% at u ≈ 9×10⁻⁷ instead of 9×10⁻⁶).

## Other CLASS gotchas

- **θ_s is not Planck's θ_MC.** CLASS `100*theta_s` ≈ 1.0419 for Planck-2018 cosmology; 1.04110 (θ_MC) gives h = 0.671.
- **P(k) from CLASS uses δ_m in the comoving gauge** (gauge-invariant) by default. Transfer functions at z_pk > z_rec need `matter_source_in_current_gauge = yes` or `output = mTk` only.
- **Newtonian and synchronous gauges disagree at the 10⁻³ level** at default precision (P(k) up to 5×10⁻³). Stay in the synchronous gauge, CLASS's default, for consistency.
- **classy rejects unread parameters.** When reusing a parameter dictionary for a different output (e.g. `mTk` only), remove `lensing` and `l_max_scalars`.

## High-k P(k) (M12)

- For k up to a few hundred h/Mpc, set `P_k_max_h/Mpc` (300 is enough to u ~ 10⁻⁹) **and raise `k_per_decade_for_pk`**. At the default of 10, the half-mode scale of drag-suppressed models is off by ~5%, because T²(k) has damped oscillations. It is converged (≤ 0.2%) at 40; 80 is used.
- Use `output = mPk` alone for these runs. classy then rejects `lensing` and `l_max_scalars` as unread, and with `output = ''` it also rejects `P_k_max_h/Mpc` and `z_pk`; drop them from the dictionary.
- **CLASS thermal-relic WDM:** add it as a second ncdm species with `m_ncdm` (eV) and `T_ncdm = (4/11)^(1/3) (94.1 eV ω_x / m)^(1/3)`, and set `omega_cdm = 0`. Do **not** pass `omega_ncdm` together with `m_ncdm`: CLASS then renormalises the phase-space density at the given temperature, so the velocity distribution is wrong.
- The Viel et al. (2005) WDM fit underestimates CLASS's k_hm by 4–17% for 1.2–6.5 keV (it was fitted at k < 5 h/Mpc).

## Distortion and thermal estimates (M13)

- `src/m13_distortions.py` uses Chluba (2016) "Method C" visibility functions and FIRAS 95% limits; the constants and references are in the docstring. The background comes from CLASS, and cgs conversion of CLASS densities is ρ[g/cm³] = ρ_CLASS[Mpc⁻²]·3c²/(8πG)/Mpc².
- The survival criterion and skin-heating estimates are deliberately simple (uniform clump at T_γ, Thomson opacity). Refine them in Experiment D before quoting them beyond order of magnitude.

## Analysis gotchas

- **Peak finding.** Detect the peaks of D_ℓ independently in each model and match them by order. Tracking within a window around the ΛCDM peaks fails once shifts exceed the window (u ≳ 10⁻²). Strongly damped models can have fewer than 7 peaks below ℓ = 2500, so those are padded with NaN. Lensed peak 7 is poorly defined for u ≳ 5×10⁻⁴.
- **Storage.** Don't store CLASS's full thermodynamics table (~10⁵ points) per model: that made one scan 385 MB. Scratch scans go under `results/scratch_*` (gitignored).
- **Shell.** Never `pkill -f <pattern>` from a command whose own command line contains the pattern; it kills that shell too.

## Dead ends

- **Newtonian gauge as a check on the θ_idm initial condition.** It made the idm–ΛCDM residual larger, not smaller: the Newtonian-gauge transformation keeps the θ_γ offset, and ΛCDM itself is less self-consistent across gauges. The decisive test was the earlier start time.

- **Chasing the ℓ-node glitch** through tolerances, k and time sampling, approximation schemes, start time, thermodynamics tables and transfer cuts (all listed above). If you pick this up again, start from the ℓ-node values themselves: compare C_ℓ at the sampled nodes between a glitchy model (u = 10⁻⁶) and a clean one (u = 10^−6.1) before splining.

## Open questions

- **Clump survival (M13 F1).** Can any (M, ρ, z_form) both survive at T_γ and hide the extra baryons from the acoustic epoch? This is now the central question (Experiment D).

- The root cause of the ℓ-node glitch (above).
- The small-scale diagnostics: at CMB-relevant u, k_½ ≲ 1.5 h/Mpc, well inside k ≤ 10 h/Mpc where P(k) is accurate. For u ≲ 10⁻⁷ the suppression scale lies beyond k = 10 h/Mpc, which would need `P_k_max_h/Mpc` raised and the high-k precision checked.
- ~~Whether the low-ℓ EE deviation comes from clump scattering~~: **resolved in M6**. It is polarization generated by clump scattering (the c_E term), present at all ℓ and largest at ℓ ≈ 15–30.
