# 2026-09-26 — M12: Experiment E, small-scale structure

**Goal:** Milestone M12 (plan revision 2, §11a). Turn the drag-induced cutoff in the linear matter power spectrum into an (approximate) bound on u, and so on Σ, and compare it with the CMB bound.
- **E1:** find and verify the literature.
- **E2:** extend P(k) to k ~ 100 h/Mpc and compute the half-mode scales.
- **E3:** get an approximate bound by warm-dark-matter half-mode matching.
- **E4:** write up the caveats.

## Configuration

- Stock CLASS v3.4.0 @ `64bbab7`, `configs/baseline.yaml`, f_cl = 1, m_idm = 10²⁷ eV, fixed h = 0.6737666, Thomson kernel (M6 showed the kernel doesn't matter for P(k)).
- P(k) only (`output = mPk`), `P_k_max_h/Mpc = 300`, `k_per_decade_for_pk = 80`; P(k, z=0) evaluated at 1200 log-spaced points, 10⁻³–285 h/Mpc.
- Grid: log₁₀u from −9 to −3 in steps of 0.25 (25 models), plus ΛCDM.

## What was run

```bash
OMP_NUM_THREADS=4 uv run python src/m12_small_scale.py run    # 26 models, ~6.5 min
uv run python src/m12_wdm_check.py                            # CLASS thermal-relic WDM, 4 masses, ~1 min
uv run python src/m12_small_scale.py analyse                  # -> results/scans/m12/summary.json, figures/m12_*.png
uv run python src/m12_fluid_check.py                          # E4 estimates -> results/scans/m12/fluid_check.json
```

## E1 — Literature (each source checked against the paper itself)

| Source | What it constrains | Result used | Status |
|---|---|---|---|
| Boehm, Schewtschenko, Wilkinson, Baugh & Pascoli 2014, MNRAS 445, L31 (arXiv:1404.7012v2) | **DM–photon** elastic scattering, from N-body counts of Milky Way satellites (classical + SDSS) | σ_DM−γ < 5.5×10⁻⁹ σ_Th (m_DM/GeV), 2σ, most conservative host-mass bin (2.3–2.7)×10¹² M☉ → **u < 5.5×10⁻⁷**. Also: σ = 2×10⁻⁹ σ_Th (m/GeV), i.e. u = 2×10⁻⁷, gives linear P(k) "similar to" 1.24 keV WDM. | ✅ verified in the full text (§4 and Fig. 1). The only direct DM–photon small-scale bound found. |
| Nadler et al. (DES Collaboration) 2021, PRL 126, 091101 (arXiv:2008.00022) | Thermal-relic WDM, Milky Way satellites (DES + PS1) | m_WDM > 6.5 keV, 95% | ✅ verified (abstract) |
| Iršič et al. 2017, PRD 96, 023522 (arXiv:1702.01764) | Thermal-relic WDM, Lyman-α (XQ-100 + HIRES/MIKE) | m_WDM > 5.3 keV (2σ); 3.5 keV with non-smooth IGM temperature evolution | ✅ verified (abstract) |
| Viel, Lesgourgues, Haehnelt, Matarrese & Riotto 2005, PRD 71, 063534 (hep-ph/0501562) | WDM transfer-function fit | T(k) = [1 + (αk)^(2ν)]^(−5/ν), ν = 1.12, α = 0.049 (m/keV)^(−1.11) (Ω_x/0.25)^(0.11) (h/0.7)^(1.22) h⁻¹Mpc (eqs. 6–7; fitted for k < 5 h/Mpc) | ✅ verified in the full text |

**Checked and not used**, because they aren't DM–photon analyses:
- Escudero et al. 2015 (arXiv:1505.06735): DM–neutrino.
- Hooper et al. 2022, "One likelihood to bind them all" (arXiv:2206.08188): feebly interacting DM, DM–baryon, mixed WDM.
- Schewtschenko et al. 2016 (arXiv:1512.06774): DM–radiation satellite *structure*; its abstract quotes no cross-section bound.

A 2026 review (Nadler, Rogers & Drlica-Wagner, arXiv:2607.28564) was found, but its abstract doesn't show whether it covers DM–photon scattering. **Not verified** and not used; worth reading for M11.

## E2 — P(k) to 300 h/Mpc

**Precision.**
- k_max doesn't matter: 300 and 600 h/Mpc give the same k_hm to 0.01%.
- The **k sampling does matter**. With CLASS's default `k_per_decade_for_pk = 10`, k_hm is 4.6% low at u = 10⁻⁸, because T²(k) has damped acoustic oscillations from the drag. At 40 per decade it is converged: 80 and 160 change k_hm by ≤ 0.2% at u = 10⁻⁹, 10⁻⁸ and 10⁻⁶. `tol_perturbations_integration` 10⁻⁶ changes nothing (≤ 0.02%).
- Production uses 80 per decade (~10 s per model).

**Result.** The suppression scale follows k_hm ≈ 10^(−2.05) u^(−0.477) h/Mpc (maximum residual 0.045 dex over u = 10⁻⁹–10⁻³). Examples:

| u | k_10 | k_½ | k_hm [h/Mpc] | M_hm [M☉/h] | equivalent m_WDM (Viel fit) |
|---:|---:|---:|---:|---:|---:|
| 10⁻⁹ | 51 | 133 | 186 | 1.8×10⁶ | 14.5 keV |
| 10⁻⁸ | 16.4 | 42.6 | 59.8 | 5.3×10⁷ | 5.2 keV |
| 10⁻⁷ | 5.3 | 13.8 | 19.3 | 1.6×10⁹ | 1.9 keV |
| 10⁻⁶ | 1.75 | 4.49 | 6.31 | 4.5×10¹⁰ | 0.69 keV |

The shape is WDM-like, with damped oscillations: a secondary bump in T² of ~0.07 at ~3k_hm (see `figures/m12_transfer.png`).

## E3 — Approximate bound by half-mode matching

A clump model is taken to be as suppressed as the thermal-relic WDM model with the same half-mode wavenumber (T = ½). Primary: the Viel fit.

| Small-scale limit | WDM k_hm [h/Mpc] (Viel / CLASS WDM) | **u_max** (Viel / CLASS) | σ/M max [cm²/g] | **Σ/Q min [g/cm²]** (Viel / CLASS) |
|---|---:|---:|---:|---:|
| MW satellites, m > 6.5 keV (95%) | 76.5 / 89.7 | **6.1×10⁻⁹ / 4.4×10⁻⁹** | 2.3×10⁻¹¹ | **4.4×10¹⁰ / 6.1×10¹⁰** |
| Lyman-α, m > 5.3 keV (2σ) | 61.0 / 70.3 | 9.6×10⁻⁹ / 7.2×10⁻⁹ | 3.6×10⁻¹¹ | 2.8×10¹⁰ / 3.7×10¹⁰ |
| Lyman-α non-smooth T(z), m > 3.5 keV (2σ) | 38.5 / 42.9 | 2.5×10⁻⁸ / 2.0×10⁻⁸ | 9.1×10⁻¹¹ | 1.1×10¹⁰ / 1.4×10¹⁰ |
| *Direct DM–γ satellites (Boehm+2014)* | — | *5.5×10⁻⁷* | *2.1×10⁻⁹* | *4.9×10⁸* |
| *CMB (Planck; Stadler & Bœhm)* | — | *2.25×10⁻⁴* | *8.4×10⁻⁷* | *1.2×10⁶* |

**Cross-checks and systematics:**
- **Calibration against Boehm+2014.** Their u = 2×10⁻⁷ "similar to 1.24 keV WDM" maps to 1.39 keV in our half-mode matching, a 12% agreement (they used a Planck 2013 cosmology). At their own bound (u = 5.5×10⁻⁷) the half-mode-equivalent WDM mass is 0.92 keV. That is consistent with satellite WDM limits of that era being ~1–2 keV. Today's 6.5 keV limit gives a ~100× tighter bound on u.
- **Viel fit vs CLASS WDM.** CLASS's own thermal-relic WDM (m_ncdm and T_ncdm set for ω_x = 0.12; Ω_m matches ΛCDM to 3×10⁻⁴) has k_hm higher than the fit by 4%, 11%, 15% and 17% at 1.24, 3.5, 5.3 and 6.5 keV. The fit is calibrated at k < 5 h/Mpc. Using CLASS instead tightens the u bounds by 20–40%. Which is appropriate depends on the transfer function behind each published limit; that has not been checked.
- **Shape.** At matched k_hm (6.7 keV, log₁₀u = −8.25), the clump suppression starts earlier (k_10 = 21.8 against 24.3 h/Mpc for WDM) but keeps more power beyond the half-mode scale (T² = 0.011 against 0.005 at 2k_hm). Matching at k_10 instead of k_hm would move u_max by ~25%.
- Overall, **treat the small-scale bounds as uncertain by a factor ~2 in u**. This is still a ~10⁴-fold improvement on the CMB.

## E4 — Caveats

1. **Epoch.** The half-mode scale at the satellite limit (k = 76.5 h/Mpc) enters the horizon at **z ≈ 2.4×10⁷** (T_γ ≈ 5.6 keV), and at the bound the drag decouples at z ≈ 7.1×10⁶. The Lyman-α limits give z ≈ 1.2–1.9×10⁷. **The bound therefore applies only if clumps with that Σ already exist at z ~ 10⁷**: after BBN, but long before any conventional structure formation. If clumps form or compact later, the drag history differs and this bound doesn't apply as stated. That is Experiment D's question, and it now determines whether the small-scale bound is relevant.
2. **Fluid approximation.** The discreteness of the clumps adds Poisson (isocurvature) power P ≈ [(3/2)(1 + z_eq)]² M_cl/ρ_cl by today. It equals the linear ΛCDM power at the limit's k_hm for M_cl ≈ 2.6 M☉/h (satellites), 4.8 (Lyman-α) and 17 (conservative Lyman-α). Planetary-mass clumps (≲10⁻³ M☉) are therefore in the fluid regime by more than three orders of magnitude. Stellar-mass clumps would not be, and their Poisson power would itself alter small-scale structure. The hierarchy M_cl ≪ M_hm (≈2.5×10⁷ M☉/h) holds comfortably.
3. **Collisionless after decoupling.**
   - Clump–clump: geometric cross-sections give σ/M = Q/Σ ≈ 2×10⁻¹¹ cm²/g at the bound, ~10¹¹ times below self-interacting-DM scales (~1 cm²/g).
   - Clump–ISM: the stopping time in 1 cm⁻³ gas at 200 km/s is Σ/(ρ_gas v) ≈ 4×10¹⁹ yr.
   - So the collisionless-DM assumption behind the WDM simulations holds. Evaporation, tidal disruption and cloud survival are **not** addressed.
4. **Linear matching only.** No simulation of this model was run; the WDM limits come from simulations or emulators with WDM initial conditions. Clump models keep somewhat more small-scale power (the oscillation bump), which Boehm+2014 also noted.
5. **f_cl < 1.** Not computed here. With only part of the dark matter suppressed, T² saturates at ~(1 − f_cl)², and the WDM mapping no longer applies; mixed-DM limits would be needed (part of Experiment C).

## Outcome

- **The small-scale structure bound is ~10⁴ times tighter in u than the CMB:** u ≲ 5×10⁻⁹, i.e. **Σ/Q ≳ 5×10¹⁰ g/cm²** (Milky Way satellites via WDM half-mode matching, factor ~2 systematic).
- The one direct DM–photon small-scale bound (Boehm+2014) gives u < 5.5×10⁻⁷ (Σ/Q > 4.9×10⁸ g/cm²), from an older satellite sample.
- **Interpretation:** for clumps making up all the dark matter and present by z ~ 10⁷, the minimum surface density is ~5×10¹⁰ g/cm², not ~10⁶ g/cm² (the CMB bound). Whether that applies depends on the clumps' formation history.

## Follow-ups

- Experiment D (formation history) is now the critical question for E: when must the clumps exist?
- Experiment C should use mixed-DM small-scale limits for f_cl < 1.
- Check which WDM transfer function underlies the Nadler+2021 and Iršič+2017 limits (Viel fit or Boltzmann code); that resolves most of the 20–40% systematic.
