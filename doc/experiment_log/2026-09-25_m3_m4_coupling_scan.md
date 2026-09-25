# 2026-09-25 — M3 + M4: Experiment A coupling scan and physical mapping

**Goal:** Milestones M3 and M4 (plan §6.1–6.3, §3; Experiment A).
- Scan the clump–photon coupling u at f_cl = 1 with stock CLASS (Thomson-like kernel).
- Build physical intuition for how the spectra respond (the six questions of §6.2).
- Extract CLASS's interaction rates and decoupling redshifts.
- Map every model onto σ/M, Σ/Q, z_dec and the small-scale suppression scales.

M5 (a detectability statistic with noise, sky fraction and covariance) is **not** done here. The ideal cosmic-variance Δχ² below is the preliminary §6.4 diagnostic only.

## Configuration

- CLASS v3.4.0 @ `64bbab7`, with `configs/baseline.yaml` **plus one new precision setting, `l_linstep: 20`** (see "Numerical problem" below). That required re-running the baseline, M2a and M2b (details below).
- f_cl = 1, `m_idm = 1e27` eV, `n_index_idm_g = 0`, u passed as `u_idm_g`.
- **Two passes over the same grid:**
  - `fixed_h`: h = 0.6737666 (the baseline value). Every standard parameter is physically unchanged, and **all diagnostics use this pass**. Fixed-h ΛCDM is bit-identical to `baseline_lcdm`.
  - `fixed_theta`: 100θ_s = 1.04190, as in the baseline. CLASS evaluates θ_s at z_rec, and clump opacity shifts z_rec, so h drifts with u. This pass records that drift.
- **Grid:**
  - null checks: log₁₀u = −8, −7, −6;
  - dense: −5 to −2 in steps of 0.125;
  - strongly coupled: −1.75 to −0.25 in steps of 0.25.

## What was run

```bash
uv run python src/m3_coupling_scan.py   # 72 CLASS runs, ~70 s on 16 cores -> results/scans/m3/{fixed_h,fixed_theta}/
uv run python src/m3_analyse.py         # -> results/scans/m3/summary.{json,md}, figures/m3_*.png
uv run python src/glitch_scan.py "project" "l_linstep 20" ...   # numerical investigation -> results/validation/m3_glitch_scan.json
```

New modules:
- `cloud_mapping.py` (u ↔ σ/M ↔ Σ/Q, using CLASS's own constants);
- `rates.py` (Γ/H from CLASS output, decoupling redshifts, clump optical depth);
- `small_scale.py` (k_10, k_½, k_hm, M_hm);
- `scan.py` (runs a grid in parallel).

## Numerical problem found and fixed: sporadic C_ℓ glitches (ℓ-node sampling)

The first scan's max|ΔC_ℓ| was **not monotonic in u**: at u = 10⁻⁶, EE differed by 1.3×10⁻³, more than at 10⁻⁵. On a fine grid (log₁₀u from −7 to −4 in steps of 0.05), **about a third of the models** carried the same deterministic error: ~1.5×10⁻³ in unlensed EE and ~4–5×10⁻⁴ in TT, largest at ℓ = 836 (sometimes 1111), and ringing with a period of ~40 in ℓ. It switched on and off from one u to the next, independently of the size of u. Glitchy models agree with *each other* to 2×10⁻⁵, so the glitch is a well-defined alternative numerical state.

**Settings that did not change it at all:**
- `tol_perturbations_integration` down to 10⁻⁸;
- source time sampling (`perturbations_sampling_stepsize` 0.03);
- k sampling (`k_step_sub` 0.04 and 0.025, `k_max_tau0_over_l_max` 3, `q_linstep` 0.2);
- the tight-coupling scheme and its triggers, including the idm–photon trigger;
- switching off the radiation-streaming approximations;
- the perturbation start time (×0.1 to ×3);
- an identical log-region thermodynamics table for all models (`thermo_z_initial` 10¹³);
- the late-source time cut (`transfer_neglect_late_source` 3000, `neglect_CMB_sources_below_visibility` down to 10⁻⁷);
- Bessel sampling (`hyper_sampling_flat` 16).

**Things checked and found clean:** the thermodynamics output (x_e and g agree to ~10⁻⁶ between glitchy and clean models at all z), and `ddmu_idm_g` against its analytic form.

**What did change it:**
- the linear thermodynamics table size, `thermo_Nz_lin`, but only by moving the glitch to other u values (a 59-model count at 100000, 120000 and 160000 gave 33–45 glitchy models, the same as the project setting's 41);
- **the ℓ-node spacing.** With `l_linstep = 20` (default 40), glitchy models fall from 41/59 to 3/59, and the worst second-difference from 1.4×10⁻³ to 3.8×10⁻⁴. The three remaining flags are all TE at ℓ = 2 (cosmic variance ~50% there), so they are irrelevant. `l_linstep = 10` and a finer `l_logstep` gave no further improvement.

**Decision:** adopt `l_linstep: 20` project-wide. Runtime is 1.25 s per model, up from 1.17 s. **The root cause is not understood.** C_ℓ values at some ℓ nodes toggle between two states, and denser nodes suppress the effect. This is recorded as an open question in `IMPLEMENTATION_NOTES.md`.

**Re-validation with `l_linstep: 20`:**
- `baseline_lcdm` recomputed; derived parameters unchanged.
- M2a still **passes**: TT 2.6×10⁻⁵, TE 1.6×10⁻⁵, EE 8.2×10⁻⁶, φφ 1.9×10⁻⁶, P(k) 4.1×10⁻⁵.
- M2b is unchanged: converged from 10⁷ eV, noise floor TT 1.9×10⁻⁵.

**Effect on M3:** the glitch had created a false early onset. EE first reached 0.1% at u = 9×10⁻⁷ with the glitch and at 9.4×10⁻⁶ without it. The cosmic-variance thresholds barely moved.

## Analysis bugs caught (outputs discarded)

- **Γ_γ/H at z_rec was reported as 1.2×10⁻⁸ for u = 10⁻⁴** instead of 6.5×10⁻⁴. I reversed the thermodynamics table assuming it runs from high to low z, but CLASS orders it by *increasing* z. The plot, an order-of-magnitude estimate and the corrected table all agree now.
- **The TT peak tracker** searched ±40 in ℓ around the ΛCDM peaks and failed at strong coupling, giving peak shifts of +20% and then −84% at u = 0.1. It was replaced by independent detection, matched by order.
- **Stored thermodynamics arrays** at CLASS's full resolution (~1.1×10⁵ points) made the scan 385 MB. They are now resampled to 2000 log-spaced redshifts, giving 20 MB with no change to any derived number.

## Validation of the rates

CLASS's `dmu_idm_g` matches the analytic a·n_cl·σ·c = (1+z)² ρ_idm,0 (σ/M) to 4×10⁻⁶. The drag rate uses CLASS's own S_idm_g = 4ρ_γ/3ρ_idm.

## Outcome

### CLASS limits (stock v3.4.0, f_cl = 1)

- Runs succeed up to **u = 10⁻¹**.
- From **u ≥ 10^−0.75** they fail. The perturbation initial conditions assume the tight-coupling approximation, but the idm–photon trigger has already switched it off at the start time. The same happens with CLASS's default start time, so this is not caused by the project settings.
- At **u = 1** thermodynamics fails because clump opacity pushes the visibility peak (z_rec ≈ 470) below CLASS's hard-coded `_Z_REC_MIN_ = 500`.

### Physical mapping (M4), f_cl = 1, fixed cosmology

| | u | σ/M [cm²/g] | Σ/Q [g/cm²] | z_dec (drag) | max ΔTT | max ΔEE | max Δφφ | Δχ²_CV TT | σ₈ | k_½ [h/Mpc] | M_hm [M☉/h] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ideal-CV threshold | 1.1×10⁻⁵ | 4.1×10⁻⁸ | 2.4×10⁷ | 1.7×10⁵ | 0.15% | 0.12% | 13% | 1.0 | 0.805 | 1.4 | 1.5×10¹² |
| | 10⁻⁴ | 3.7×10⁻⁷ | 2.7×10⁶ | 5.6×10⁴ | 1.2% | 1.0% | 54% | 94 | 0.770 | 0.52 | 3.0×10¹³ |
| published Planck bound | 2.25×10⁻⁴ | 8.4×10⁻⁷ | 1.2×10⁶ | 3.8×10⁴ | 2.8% | 2.7% | 75% | 608 | 0.733 | 0.36 | 8.9×10¹³ |
| | 10⁻³ | 3.7×10⁻⁶ | 2.7×10⁵ | 1.8×10⁴ | 12% | 13% | 97% | 1.7×10⁴ | 0.606 | 0.20 | 5.9×10¹⁴ |

- The full table is in `results/scans/m3/summary.md`.
- "Max" means over 2 ≤ ℓ ≤ 2500, lensed.
- Δχ²_CV is ideal full-sky cosmic variance for one spectrum over 2 ≤ ℓ ≤ 2500 (§6.4). It is **not** an observational statement.
- **Thresholds** (u at which ideal-CV Δχ² = 1): TT 1.11×10⁻⁵, EE 1.07×10⁻⁵, TE 1.32×10⁻⁵; unlensed TT 2.0×10⁻⁵, unlensed EE 8.6×10⁻⁶.
- **Onsets:** the maximum deviation reaches 0.1% at u = 6.3×10⁻⁶ (TT), 9.4×10⁻⁶ (EE) and 6×10⁻⁸ (φφ).

### Mechanism: the photons barely see the clumps; the CMB responds through gravity

- **Photon-side opacity is negligible.** Before recombination, clump opacity is ≤ 1.1×10⁻⁴ of Thomson opacity at the Planck bound (5×10⁻⁵ at u = 10⁻⁴). After recombination the clump optical depth is τ_cl ≈ 1.0×10⁻³ (4.5×10⁻⁴ at u = 10⁻⁴), against τ_reio = 0.054. Γ_γ/H at z_rec is 1.5×10⁻³.
- **What matters is the drag on the clumps**, whose rate is larger by 4ρ_γ/3ρ_cl. Clumps are dragged until z_dec ≈ 5.6×10⁴ (u = 10⁻⁴), and modes that enter the horizon before then have their clump perturbations suppressed. Those are the k ≳ 0.2 h/Mpc modes where P(k) is cut off. The CMB then responds through the weaker gravitational potentials and through lensing.
- **Consequence:** at CMB-relevant couplings the observable signal is set by the ℓ = 1 momentum transfer σ_mt, which does not depend on the scattering kernel. One place where a photon-side contribution may show up is **low-ℓ EE**. EE departs from ΛCDM first at ℓ ≈ 7–17 for u ≳ 10^−4.5, with a low-ℓ change of ~0.4% at u = 10⁻⁴. The most plausible source is polarization generated by clump scattering after recombination (τ_cl ~ 4.5×10⁻⁴) through the Thomson kernel, but this has **not been isolated**: τ_cl builds up mostly just after recombination rather than at reionization, and the simple estimate 2τ_cl/τ_reio ≈ 1.7% overstates what is seen. If confirmed, it is exactly the kind of term Experiment B would change. It is small next to cosmic variance at ℓ < 20 and degenerate with τ.

### Answers to the six §6.2 questions

1. **Does increasing u mainly alter peak heights?** Yes, at CMB-relevant u. The damping tail is suppressed and the first two to three peaks are slightly enhanced. At u = 10⁻⁴: mean TT residual +0.12% for ℓ = 30–800 and −0.49% for ℓ = 1500–2500. Unlensed peak heights change by +0.09%, +0.22%, +0.01%, +0.12%, −0.32%, −0.15% and −0.79% for peaks 1–7.
2. **Phase shifts?** Small but systematic: peaks move to higher ℓ, more so at high ℓ. At fixed H0, u = 10⁻⁴ shifts the unlensed peaks 1–7 by Δℓ = 0.08, 0.18, 0.39, 0.51, 0.79, 0.84 and 1.18. At 10⁻³ the shifts are 0.7 to 8.7. With θ_s fixed, h moves by only 3×10⁻⁵ up to u = 10⁻³ (0.67377 → 0.67361). It moves by 0.2% at 10⁻² and 2.4% at 10⁻¹ (0.6574).
3. **Where does the effect first become significant?**
   - In lensed TT, at the highest multipoles: the 0.1% level first appears at ℓ ≈ 2456 (u = 10⁻⁵) and ℓ ≈ 877 (10^−4.5).
   - In EE, at low ℓ (see the mechanism section).
   - Lensing roughly doubles the TT signal at small u: max ΔTT is 0.14% lensed against 0.06% unlensed at u = 10⁻⁵. For EE, lensing partly *cancels* the signal: at u = 10⁻⁴ the unlensed Δχ² is 143 and the lensed one 82.
4. **Is the damping tail especially sensitive?** Yes. It carries the largest and first TT deviations, through both the unlensed damping-tail suppression and the change in lensing.
5. **When does matter-power suppression matter?** Suppression within k ≤ 10 h/Mpc starts at u ~ 10⁻⁷ (k_10 = 4.9 h/Mpc), far below CMB sensitivity. k_½ = 1.5, 0.52, 0.36 and 0.20 h/Mpc at u = 10⁻⁵, 10⁻⁴, 2.25×10⁻⁴ and 10⁻³. φφ is the most sensitive CMB observable: 0.16% at u = 10⁻⁷ and 1.5% at 10⁻⁶.
6. **Does strong coupling approach a photon-loaded fluid?** Partly, and CLASS cannot go further. If the clumps fully loaded the photon–baryon fluid, the sound horizon would shrink to 0.727 × standard, shifting the peaks by +37.5%. At u = 0.1, the largest runnable coupling, the unlensed peaks shift by +4.9% to +11.3% (mean +8.2%). The drag decouples at z_dec ≈ 2200, just before recombination, so the loading is only partial when the spectrum forms.

### Small-scale structure (diagnostic only, plan §6.3)

At the published CMB bound, M_hm ≈ 9×10¹³ M☉/h and σ₈ falls by 10%: structure below galaxy-cluster scale is strongly suppressed in linear theory. As plan §6.3 requires, this is **not** turned into a constraint here. It does suggest that a DM–photon small-scale analysis would be much more constraining than the CMB for this model.

## Follow-ups

- **Experiment B (for the user to decide):** since the photon-side kernel enters only at ≲10⁻⁴ of Thomson before recombination and through τ_cl ≲ 10⁻³ afterwards, B's expected effect on TT/TE/EE is tiny, apart from low-ℓ EE. Consider estimating it first (e.g. by scaling the clump-generated post-recombination polarization) before patching CLASS.
- **M5:** add noise, sky fraction, the TT/TE/EE covariance and a realistic φφ reconstruction noise. φφ dominates the ideal sensitivity, so its treatment will matter.
- **Open:** the root cause of the ℓ-node glitch.
- **Open:** a literature DM–photon small-scale bound, to put the M_hm numbers in context.
