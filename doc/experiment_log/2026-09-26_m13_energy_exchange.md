# 2026-09-26 — M13: Experiment F, energy exchange and spectral distortions (estimates)

**Goal:** Milestone M13 (plan revision 2, §11b), as estimates only.
- **F1:** how tightly clumps are thermally locked to the CMB, and what that implies.
- **F2:** the energy-exchange budget.
- **F3:** μ and y distortions compared with FIRAS.
- **F4:** interpretation, and a go/no-go on a full Boltzmann treatment.

## Configuration and method

- Script: `src/m13_distortions.py` → `results/scans/m13/summary.json` and `figures/m13_thermal.png`.
- Background (H, ρ_γ, ρ_cl = ρ_dm) from stock CLASS at the baseline cosmology (fixed h). Photon absorption rate Γ_γ = (1+z)·`dmu_idm_g` from the validated analytic form (`rates.analytic_dmu`).
- **Clump model:**
  - opaque grey absorber (Q = 1, so absorption and momentum-transfer cross-sections coincide and u carries over);
  - f_cl = 1;
  - ionised H/He, mean molecular weight 0.59, so c_v = (3/2)k/(0.59 m_p) = 2.1×10⁸ erg g⁻¹ K⁻¹;
  - Thomson opacity κ_T = 0.35 cm²/g.
- **Couplings examined:**
  - CMB bound u = 2.25×10⁻⁴ (Σ/Q = 1.2×10⁶ g/cm²);
  - direct satellite bound u = 5.5×10⁻⁷ (Σ/Q = 4.9×10⁸);
  - M12 small-scale bound u = 6.1×10⁻⁹ (Σ/Q = 4.4×10¹⁰).

**Verified references:**
- COBE/FIRAS: |μ| < 9×10⁻⁵ and |y| < 1.5×10⁻⁵ (95%). Fixsen et al. 1996, ApJ 473, 576; checked against the abstract.
- Distortion visibility functions from Chluba 2016, MNRAS, arXiv:1603.02496, eqs. 7, 8, 11 ("Method C"), checked in the full text:
  - y = ¼∫ d(Q/ρ_γ)/dz · J_y dz and μ = 1.401 ∫ d(Q/ρ_γ)/dz · J_μ dz;
  - J_bb = exp[−(z/z_th)^(5/2)], z_th = 1.98×10⁶;
  - J_y = [1 + ((1+z)/6×10⁴)^2.58]⁻¹ (z ≥ z_rec);
  - J_μ = J_bb·{1 − exp[−((1+z)/5.8×10⁴)^1.88]}.

## F1 — Thermal locking and its consequence: survival

**Locking.** A clump in an isotropic blackbody field absorbs πR²·c·aT_γ⁴ and emits πR²·c·aT_cl⁴. Linearising, its surface relaxes to T_γ on

  t_th = Σ c_v / (4 c a T_γ³).

A clump following the cooling CMB lags by δT/T = H·t_th, which is **≪ 1 at every relevant epoch and Σ**. The largest value, at z = 10⁴ and Σ = 4.4×10¹⁰ g/cm², is 1.2×10⁻³. For the CMB-bound Σ it is ≤ 3×10⁻⁸ at z ≥ 10⁴. **Clump surfaces are held at T_γ.**

**Heating through.** Radiative diffusion heats a Thomson-thick column Σ_skin = √(t_H ρ c/κ_T) per Hubble time, with ρ = 3Σ/4R. For planetary-mass clumps (10⁻⁶–10⁻³ M☉) the whole clump is heated within a Hubble time at z ≲ 10⁵–10⁶ for every Σ considered. At z ~ 10⁷ and the small-scale-bound Σ, 1–8% of the column is heated per Hubble time. (Order of magnitude: it ignores the opacity of any cold neutral interior.)

**Survival.** Material at T_γ is gravitationally bound only if G M μ m_p/R ≳ k T_γ.
- With M = πR²Σ this gives **M_min(z, Σ) = (k T_γ/(G μ m_p))² / (πΣ)**:

| z (T_γ) | M_min at Σ = 1.2×10⁶ (CMB bound) | at 4.9×10⁸ (Boehm+2014) | at 4.4×10¹⁰ (M12) |
|---|---:|---:|---:|
| 10⁷ (2.3 keV) | 4×10⁵ M☉ | 1×10³ M☉ | 12 M☉ |
| 10⁶ (235 eV) | 4×10³ M☉ | 11 M☉ | 0.12 M☉ |
| 10⁵ (23 eV) | 44 M☉ | 0.11 M☉ | 1.2×10⁻³ M☉ |
| 4×10⁴ (9.4 eV) | 7 M☉ | 0.017 M☉ | 1.9×10⁻⁴ M☉ |
| 10⁴ (2.3 eV) | 0.44 M☉ | 1.1×10⁻³ M☉ | 1.2×10⁻⁵ M☉ |

- Equivalently, a uniform clump of internal density ρ is bound only below the redshift where T_γ = (4πGρμm_p/3k)·(3M/4πρ)^(2/3):

| M | ρ = 0.1 g/cm³ | ρ = 1 g/cm³ | ρ = 10 g/cm³ |
|---|---:|---:|---:|
| 10⁻⁶ M☉ (~Earth) | z < 210 | z < 450 | z < 960 |
| 10⁻³ M☉ (~Jupiter) | z < 2.1×10⁴ | z < 4.5×10⁴ | z < 9.6×10⁴ |
| 1 M☉ | z < 2.1×10⁶ | z < 4.5×10⁶ | z < 9.6×10⁶ |

**Consequence.** For the drag bounds to apply, the clumps must already exist, with that Σ, when the relevant modes are affected. That is z ~ 4×10⁴ for the CMB bound (drag decoupling at u = 2.25×10⁻⁴, M3) and z ~ 10⁷ for the small-scale bound (M12). More fundamentally, for extra baryons to stay out of the photon–baryon fluid at all, they must be clumped **before recombination**. Planetary-mass clumps at ordinary densities cannot be bound at those epochs: Earth-mass clumps only after z ~ 200–1000, which is at or after recombination. Clumps at the M12 epoch (z ~ 10⁷) must be stellar-mass or heavier (M ≳ 10 M☉ at Σ ~ 4×10¹⁰ g/cm²), or of degenerate density. That, in turn, is in tension with the fluid approximation (M_cl ≪ 3 M☉, M12 E4).

## F2 — Energy-exchange budget

Three cases, each written as d(Q/ρ_γ)/dz:
- **(a) Passive, thermally locked blackbody clumps.** They are bound, so they don't cool adiabatically; they follow T_γ by radiating their heat to the photons: dQ/dt = ρ_cl c_v H T_γ. This is **independent of u**, limited by the clumps' heat capacity. Per e-fold it is ρ_cl c_v/(aT_γ³) ~ 10⁻⁹ of the photon energy.
- **(b) A sustained temperature offset** δ = (T_cl − T_γ)/T_γ, maintained by some internal process: dQ/dt = 4Γ_γ δ ρ_γ, proportional to u·δ.
- **(c) An internal luminosity** per mass L/M, all emitted into the photons: dQ/dt = ρ_cl L/M, **independent of u**.

## F3 — Distortions

| Case | μ | y | FIRAS implication |
|---|---:|---:|---|
| (a) passive blackbody clumps | +1.6×10⁻⁸ | +3.1×10⁻⁹ | 6000× and 5000× below the limits; comparable to standard-ΛCDM distortions (~10⁻⁸) |
| (c) per 1 L☉/M☉ of constant internal luminosity | 1.1×10⁻⁹ | 1.8×10⁻⁸ | L/M < 850 L☉/M☉ (y); < 8.4×10⁴ L☉/M☉ (μ) |
| (b) per unit δ, CMB-bound u = 2.25×10⁻⁴ | 25.8 | 0.19 | \|δ\| < 3.5×10⁻⁶ during the μ era |
| (b) per unit δ, u = 5.5×10⁻⁷ | 0.063 | 4.8×10⁻⁴ | \|δ\| < 1.4×10⁻³ |
| (b) per unit δ, u = 6.1×10⁻⁹ | 7×10⁻⁴ | 5×10⁻⁶ | \|δ\| < 0.13 |

**Thermalisation.** At the CMB-bound coupling, clumps absorb and re-emit every photon more than once per Hubble time (Γ_γ > H) above z ≈ 3.7×10⁵, inside the μ era. Blackbody clumps would then partially *erase* μ distortions from other sources, such as the ΛCDM acoustic-damping μ ~ 10⁻⁸. That is below FIRAS sensitivity, but a possible signature for a future spectrometer. For the small-scale-allowed couplings (u ≲ 5×10⁻⁷), Γ_γ < H throughout z < 10⁸.

## F4 — Interpretation and go/no-go

1. **Spectral distortions do not constrain passive opaque clumps.** Energy exchange is limited by the clumps' own heat capacity (~10⁻⁹ of ρ_γ per e-fold), giving μ ~ 2×10⁻⁸ and y ~ 3×10⁻⁹, orders of magnitude below FIRAS. As expected, a blackbody clump at T_γ escapes.
2. **FIRAS constrains only clump microphysics:** internal luminosity (L/M ≲ 10³ L☉/M☉, weak) and sustained temperature offsets, whose limit scales as 1/u (|δ| ≲ 3.5×10⁻⁶ at the CMB-bound coupling). Future spectrometers could see the thermalisation (μ suppression) at large u.
3. **The significant result of Experiment F is F1, not F3.** Thermal locking is fast, so clumps are heated to T_γ. Gravitational binding at T_γ then requires clump masses far above planetary for the epochs where the CMB and small-scale constraints test the clumps. Taken at face value, planetary-mass baryonic clumps cannot be the dark matter during the acoustic epoch, whatever their photon coupling. That questions the premise of Experiments A–E for the planetary-mass case, and points directly to Experiment D (formation history) together with a clump-survival criterion.
4. **Go/no-go:** **no-go** on a full Boltzmann treatment of energy exchange for passive clumps (effects ≲10⁻⁸). **Go** on making clump survival (F1) part of Experiment D.

## Caveats

- These are order-of-magnitude estimates. The binding criterion treats the clump as uniform material at T_γ; a clump with a cold, opaque core and a heated skin would survive longer, losing mass by ablation of its surface layers. The skin-heating estimate uses Thomson opacity only.
- Non-thermal support is ignored (e.g. degeneracy pressure, which requires densities far beyond planetary), as are magnetic fields and rotation.
- The distortion visibilities are the "Method C" approximation (10–20% accuracy on the μ/y split).
- Case (a) assumes the whole clump stays locked. If only the skin is locked, the heat released, and the distortion, are smaller.

## Follow-ups

- Proposed plan change (in `RESULTS.md`): promote Experiment D, built around a survival/formation criterion. When can clumps of (M, ρ, Σ) exist, and does any combination both survive and hide the extra baryons from the acoustic epoch?
