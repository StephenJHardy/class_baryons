# 2026-09-26 — Likelihood set-up (M7/M8 first step) and a correction to the M5 comparison

**Goal:** Set up real-data Planck likelihoods with the project's CLASS configuration, check which data the published bound used, and start profiling the likelihood in u by minimisation (the user asked for minimisation rather than MCMC as a first step, pending a decision on compute).

## Correction to M5: which data the published bound used

Reading Stadler & Bœhm 2018 (JCAP 10, 009; arXiv:1802.06589v3) in full:
- **u < 2.25×10⁻⁴** (the "most conservative" bound, 95% CL) is from **Planck 2015 TT + lowTEB** (Table II). That is temperature only at high ℓ, from the 2015 release.
- With Planck 2015 lensing it becomes 1.81×10⁻⁴.
- **With high-ℓ polarization (Planck 2015 TT,TE,EE + lowTEB, Table III) it is u < 1.58×10⁻⁴**, and 1.49×10⁻⁴ with lensing.
- Their analysis used MontePython with a modified CLASS v2.6 and sampled the six standard parameters plus u. The prior on u isn't stated explicitly in the text; 10⁴u is reported linearly, suggesting a flat prior. **Not verified.**

The M5 log compared the Planck-like forecast (TT,TE,EE with 2018 noise; u₉₅ = 1.1×10⁻⁴ at 1.645σ, 1.3×10⁻⁴ at 1.96σ) with 2.25×10⁻⁴ and called it a "factor ~2" difference. **The like-for-like comparison is with 1.58×10⁻⁴** (TT,TE,EE + low ℓ). The forecast is 0.70–0.82 of that, i.e. consistent to about 20–30%, which is as expected given 2015 against 2018-like noise and the idealisations of the Fisher forecast. `RESULTS.md` is updated accordingly.

For M7: Cobaya 3.6 no longer ships Planck 2015 likelihoods. Reproducing the 2015 TT+lowTEB number exactly would need the 2015 clik files from the Planck Legacy Archive, installed by hand. The 2018 likelihoods are the natural target for new results.

## Likelihood installation

- **Cobaya 3.6.2**, added to the uv environment, with packages at `/media/stephen/astro/class_baryons/cobaya_packages` (3.3 TB free on that disk).
- **Installed:**
  - `planck_2018_lowl.TT` and `planck_2018_lowl.EE`, both native Python;
  - `planck_2018_highl_plik.TTTEEE_lite_native`;
  - `planck_2018_highl_plik.TTTEEE`, the full plik run through **clipy** (a Python reimplementation of clik). Cobaya installs it with `pip`, so the uv environment needed `pip` added as a dev dependency, and `astropy` added.
- **clipy self-test** on Planck's reference spectra: log L = −1172.47, the Planck 2018 V Table 20 value (χ² = 2344.94), to 4×10⁻⁶.
- **Set-up:** `src/likelihood_setup.py` builds the Cobaya input from `configs/baseline.yaml`, so the M2–M5 precision settings are used unchanged.
  - Clump model: f_cl = 1, with ω_cdm = 0 and ω_idm = ω_dm via Cobaya lambdas.
  - Sampled parameters: ω_b, ω_dm, 100θ_s (CLASS θ_s), n_s, ln 10¹⁰A_s and τ, plus the Planck nuisance parameters (1 for lite, 21 for full plik).
- **Checks at the baseline cosmology** (nuisance parameters at reference values):
  - low-ℓ TT χ² = 23.46 (Planck 2018 V: 23.25 at their best fit);
  - low-ℓ EE −2 log L = 396.11 (Planck: 2 × 198.02);
  - f_cl = 1, u = 0 matches ΛCDM to Δχ² = 0.02.
- **Speed:** one likelihood evaluation takes 1.1 s on 16 threads, lite or full.

## Minimisation: the numerical-noise problem

- **Cobaya's BOBYQA minimiser** (lite likelihood, ΛCDM point) scatters by Δ(−log L) = 0.6–1.2 between restarts. That is far too much for a profile, which needs Δχ² precision of ~0.1. Two causes were found:
  1. **Scaling.** Planck's covmats use θ_MC and ω_cdm, not the project's θ_s and ω_dm, so the automatic covmat didn't match. A covmat in the project's parameterization was built from the M5 Fisher matrix, u-fixed block (`src/fisher_covmat.py`). The scatter remained.
  2. **Noise.** Perturbing the parameters by 10⁻³ σ raises −log L by 0.025–0.032 every time. The CLASS likelihood has small discrete jumps, presumably from parameter-dependent sampling grids. They defeat BOBYQA's local quadratic model, and its `objfun_has_noise` mode with restarts was too slow to be practical (stopped).
- **Response:** `src/quadfit_min.py` minimises by **iterated quadratic regression** in Fisher-whitened coordinates. It evaluates ~100 points within ~½σ of the current centre in parallel, fits a full quadratic by least squares (averaging the jumps), moves to the fitted minimum and repeats. It reports the fitted minimum with a regression uncertainty, and the Hessian.

## Quadratic-fit minimiser: test

- **ΛCDM point (u = 0, lite):** the fitted minimum was stable at −log P = 480.376–480.382 over iterations 3–6. That is 0.16 below BOBYQA's best restart, with residual noise in the regression of 0.015.
- **Scaling:** the whitened Hessian eigenvalues spanned 0.04–12, so the Fisher covmat is off by up to 5× in some directions (τ–A_s among them, as M5 anticipated). The position therefore wandered along flat directions while the *value* stayed pinned.
- **Fix:** the minimiser now re-whitens with the fitted curvature after each iteration and stops on the predicted decrease (< 0.01 in −log P).
- **Result:** each profile point then converged in 2 iterations (432 evaluations, ~5.8 min on 4 workers × 4 threads), with a regression uncertainty of ±0.002–0.003 in −log P.

## First result: profile likelihood in u (Planck 2018 lowl TT + lowl EE + plik_lite TT,TE,EE; f_cl = 1)

```bash
uv run python src/profile_u.py quadfit plik_lite 0 2e-5 5e-5 1e-4 1.5e-4 2e-4 3e-4 4e-4   # ~47 min
uv run python src/profile_analyse.py plik_lite   # -> results/likelihood/plik_lite/profile_summary.json, figures/profile_u_plik_lite.png
```

| u | −log P | Δχ² | H0 | σ₈ | χ² plik_lite |
|---:|---:|---:|---:|---:|---:|
| 0 | 480.382 ± 0.003 | 0 | 67.16 | 0.812 | 584.38 |
| 2×10⁻⁵ | 480.661 | 0.56 | 67.16 | 0.802 | 584.88 |
| 5×10⁻⁵ | 481.060 | 1.36 | 67.15 | 0.790 | 585.61 |
| 10⁻⁴ | 481.763 | 2.76 | 67.15 | 0.773 | 587.00 |
| 1.5×10⁻⁴ | 482.533 | 4.30 | 67.16 | 0.759 | 588.60 |
| 2×10⁻⁴ | 483.400 | 6.04 | 67.18 | 0.746 | 590.46 |
| 3×10⁻⁴ | 485.471 | 10.18 | 67.22 | 0.724 | 594.82 |
| 4×10⁻⁴ | 488.053 | 15.34 | 67.29 | 0.705 | 600.27 |

- **Validation at u = 0** against Planck 2018 VI Table 2 (TT,TE,EE+lowE): ω_b = 0.02232 (0.02236 ± 0.00015), ω_dm = 0.1204 (0.1202 ± 0.0014), n_s = 0.9635 (0.9649 ± 0.0044), τ = 0.0539 (0.0544), H0 = 67.16 (67.27 ± 0.60). All within 0.3σ.
- **Limits:** the one-sided 95% limit (Δχ² = 2.71) is **u < 9.8×10⁻⁵ (σ/M < 3.7×10⁻⁷ cm²/g, Σ/Q > 2.7×10⁶ g/cm²)**. The Δχ² = 3.84 crossing is at u = 1.35×10⁻⁴ (Σ/Q > 2.0×10⁶).
- **Comparison:**
  - The M5 Fisher forecast crosses 2.71 at 1.1×10⁻⁴.
  - The published Bayesian bound from Planck 2015 TT,TE,EE+lowTEB is 1.58×10⁻⁴.
  - The profile is somewhat tighter than the published number; the datasets (2018 against 2015) and methods (profile against posterior) differ.
- **Shape:**
  - Δχ² rises almost linearly from u = 0 (≈ 0.27 per 10⁻⁵ up to 10⁻⁴), then more slowly than the Fisher parabola.
  - A single parabola fitted in u has its minimum at an unphysical u ≈ −2.7×10⁻⁴ with curvature width σ ≈ 1.6×10⁻⁴.
  - So the data lean mildly against the direction in which u moves the spectra. The rise comes almost entirely from the high-ℓ likelihood (plik_lite), with the low-ℓ terms nearly flat.
  - *Hypothesis, not tested:* u lowers the lensing amplitude (σ₈ falls from 0.81 to 0.77 at u = 10⁻⁴), and Planck's temperature spectrum is known to prefer somewhat more lensing smoothing than ΛCDM predicts. That would produce such a one-sided pull. A test is to repeat the profile with the lensing amplitude A_L free, or with unlensed-equivalent information.
- **Degeneracies along the profile:** ω_b, ω_dm and θ_s drift upward and n_s and A_s slightly, as the M5 Fisher analysis predicted. H0 barely changes.

## Caveats and next steps

- **Lite likelihood.** It pre-marginalises Planck's foreground and calibration nuisance parameters. The full plik likelihood (21 nuisance parameters) has not yet been profiled; with the quadratic-fit method it costs ~1,200 evaluations per iteration (~1–1.5 h per u point here). A check at u = 0 and u = 10⁻⁴ would show whether lite and full agree at the crossing.
- **No lensing reconstruction.** Adding it would likely tighten the bound (Stadler & Bœhm find ~6% with 2015 data).
- **Minimisation only (no MCMC),** as requested pending a decision on compute.
