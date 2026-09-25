# Results summary

*Current state of the project. Rewritten in place as results change; the history is in [`experiment_log/`](experiment_log/). Last updated 2026-09-25.*

## Headline

**Experiment A (Thomson-like clumps, stock CLASS), with fixed cosmology and no real data yet:**

- **With all dark matter in clumps (f_cl = 1), the CMB becomes sensitive at u ~ 10⁻⁵**, i.e. σ/M ~ 4×10⁻⁸ cm²/g or Σ/Q ~ 2×10⁷ g/cm², *if* measurements were limited only by full-sky cosmic variance. Real noise, sky coverage and refitting the other parameters will weaken this. The published Planck bound, u < 2.25×10⁻⁴ (Σ/Q ≳ 1.2×10⁶ g/cm²), is about 20× weaker, which is as expected; reproducing it is M7.
- **The CMB signal comes almost entirely from radiation *dragging the clumps* before z ~ 10⁴–10⁵, not from photons scattering off clumps.** At CMB-relevant couplings the clumps contribute ≲10⁻⁴ of the Thomson opacity before recombination and an optical depth τ_cl ≲ 10⁻³ after it. The observable effects are suppressed clump perturbations on small scales, weaker potentials (damping-tail TT/EE) and much less lensing.
- **Implication for Experiment B (open decision):** the drag depends only on the momentum-transfer cross-section, not on the angular or polarization kernel. So replacing the Thomson kernel with an opaque-clump kernel is expected to change TT/TE/EE only slightly, apart from a possible small low-ℓ EE term. See [the M3 log](experiment_log/2026-09-25_m3_m4_coupling_scan.md).
- **Small scales look far more constraining than the CMB (a diagnostic, not a constraint).** At the Planck bound the linear half-mode mass is ~9×10¹³ M☉/h and σ₈ falls by 10%.

## Status of milestones and validation gates

| Milestone | Status | Key number | Record |
|---|---|---|---|
| M0 environment | ✅ | CLASS v3.4.0 @ `64bbab7` | [log](experiment_log/2026-09-25_m0_environment.md) |
| M1 baseline | ✅ | h = 0.6738, σ₈ = 0.8107 | [log](experiment_log/2026-09-25_m1_baseline.md) |
| M2a u = 0 reproduces CDM | ✅ **pass** | max residual 4.1×10⁻⁵ (tolerance 10⁻⁴); re-passed with `l_linstep = 20` | [log](experiment_log/2026-09-25_m2a_zero_coupling.md), [M3 log](experiment_log/2026-09-25_m3_m4_coupling_scan.md) |
| M2b cold-mass limit | ✅ **pass** | converged from ~10⁷ eV; default m_idm = 10²⁷ eV | [log](experiment_log/2026-09-25_m2b_mass_convergence.md) |
| M2c reproduce published bound | ⏳ pending (at M7) | target u < 2.25×10⁻⁴ | — |
| M3 coupling scan (f_cl = 1) | ✅ | ideal-CV onset u ≈ 1.1×10⁻⁵ | [log](experiment_log/2026-09-25_m3_m4_coupling_scan.md) |
| M4 physical mapping | ✅ | table below | [log](experiment_log/2026-09-25_m3_m4_coupling_scan.md) |
| M5 – M11 | ⏳ not started | | |

| Experiment | Status |
|---|---|
| A: Thomson-like clumps | scan and mapping done; likelihood (M5, M7, M8) next |
| B: opaque-clump kernel | not started; M3 suggests its effect is small (decision needed) |
| C: partial fraction | not started |
| D: formation history | not started |

## Experiment A: coupling scan at f_cl = 1 (fixed H0, fixed cosmology)

| | u | σ/M [cm²/g] | Σ/Q [g/cm²] | z_dec (drag) | max ΔTT | max ΔEE | max Δφφ | Δχ²_CV (TT) | σ₈ | k_½ [h/Mpc] | M_hm [M☉/h] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ideal-CV threshold | 1.1×10⁻⁵ | 4.1×10⁻⁸ | 2.4×10⁷ | 1.7×10⁵ | 0.15% | 0.12% | 13% | 1 | 0.805 | 1.4 | 1.5×10¹² |
| | 10⁻⁴ | 3.7×10⁻⁷ | 2.7×10⁶ | 5.6×10⁴ | 1.2% | 1.0% | 54% | 94 | 0.770 | 0.52 | 3.0×10¹³ |
| Planck bound (literature) | 2.25×10⁻⁴ | 8.4×10⁻⁷ | 1.2×10⁶ | 3.8×10⁴ | 2.8% | 2.7% | 75% | 608 | 0.733 | 0.36 | 8.9×10¹³ |
| | 10⁻³ | 3.7×10⁻⁶ | 2.7×10⁵ | 1.8×10⁴ | 12% | 13% | 97% | 1.7×10⁴ | 0.606 | 0.20 | 5.9×10¹⁴ |

*Assumptions:*
- stock CLASS idm–photon model (Thomson-like kernel), f_cl = 1, all standard parameters fixed at the baseline (h fixed);
- "max" is over 2 ≤ ℓ ≤ 2500 (lensed);
- Δχ²_CV is ideal full-sky cosmic variance for a single spectrum, **not** a data constraint;
- M_hm is the linear half-mode mass (T = ½).

The full 33-point table is in [`results/scans/m3/summary.md`](../results/scans/m3/summary.md).

![M3 summary](../figures/m3_summary.png)

**How the spectra respond** (details in the log):
- The damping tail is suppressed and the low peaks slightly enhanced.
- The peaks shift slightly to higher ℓ (Δℓ ≈ 0.1–1.2 at u = 10⁻⁴).
- Lensing roughly doubles the TT signal at small u, while it partly cancels the EE signal.
- φφ is the most sensitive observable: 1% at u ≈ 6×10⁻⁷.
- With θ_s held fixed instead of H0, the inferred h hardly moves for u ≤ 10⁻³.

![TT](../figures/m3_tt.png)

**Interaction rates.** Clump drag stays coupled until z_dec (Γ_cl/H = 1), long before recombination for all CMB-relevant u. The photon opacity Γ_γ/H at recombination is ~10⁻³.

![rates](../figures/m3_rates.png)

**Strong coupling.** CLASS runs up to u = 0.1. At that coupling the TT peaks shift by +8% on average, compared with +37.5% for a fully photon-loaded fluid; the drag decouples at z ≈ 2200, just before recombination. Beyond u ≈ 0.13–0.18, stock CLASS fails (see the implementation notes).

## Validation results

- **M2a:** with all dark matter as zero-coupling idm, the spectra match ΛCDM to TT 2.6×10⁻⁵, TE 1.6×10⁻⁵, EE 8.2×10⁻⁶, φφ 1.9×10⁻⁶ and P(k) 4.1×10⁻⁵. This needs three non-default CLASS precision settings: a 4× thermodynamics table, a 10× earlier perturbation start and `l_linstep = 20` (see the implementation notes).
- **M2b:** outputs are independent of `m_idm` from ~10⁷ eV, and bit-identical from 10²⁴ eV.

![M2a residuals](../figures/m2a_residuals.png)
