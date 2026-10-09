# 2026-10-02 to 10-05: CDM + clump mixtures, two scattering kernels, and a minimiser audit

Branch `clump-fraction-kernels`. Profile likelihoods only; no MCMC.

## Question and model

**Question.** Do the CMB anisotropies exclude a significant baryonic dark component if it consists of existing, compact, pressureless, optically thick clumps? The clumps' existence is assumed. Formation, survival, BBN and small-scale structure are outside this task.

**Model.**
- **Two dark components:** ω_cl = f_cl·ω_dark and ω_CDM = (1 − f_cl)·ω_dark.
  - ω_dark is fitted (`omega_dm`).
  - The diffuse baryon density ω_b is a separate fitted parameter.
  - CDM has no photon coupling. The clumps have coupling u = `u_idm_g`.
- **Kernels:** compared at equal u, with patched CLASS (`class_patches/idm_g_kernel.patch`):
  - **Thomson:** `idm_g_quadrupole_coefficient = 1`, `idm_g_polarization_coefficient = 1`. Standard angular redistribution and polarization generation.
  - **Uniform / isotropic:** both coefficients 0. Uniform re-emission in the clump frame, with no polarization generation.
  - The patch was re-read for this task. Both coefficients enter only the regeneration and source terms (ℓ = 2 temperature P0 term, polarization hierarchy, line-of-sight Π and E sources), with damping at the full rate. The tight-coupling closure is generalised exactly (M6 log).
- **Likelihoods:** planck_2018_lowl.TT + planck_2018_lowl.EE + planck_NPIPE_highl_CamSpec.TTTEEE, the same as the headline analysis.
- **Theory settings:** unchanged from `configs/baseline.yaml`:
  - thermodynamics table 80000/20000;
  - early perturbation start;
  - `l_linstep 20`;
  - `non_linear: none` (linear spectra; the lensing potential is linear);
  - m_idm = 10²⁷ eV, n_index_idm_g = 0;
  - one massive ν of 0.06 eV, N_ur = 2.0328.
- **Objective profiled:** −log P = Σ(−log L) − log(prior). The prior term contains CamSpec's Gaussian nuisance priors and constant flat-prior normalisations. It is the same at every point; Δχ² = 2Δ(−log P).
- **Free parameters (15):** ω_b, ω_dark, 100θ_s, n_s, ln(10¹⁰A_s), τ, and all 9 CamSpec nuisance parameters (A_planck, amp_143, amp_217, amp_143x217, n_143, n_217, n_143x217, calTE, calEE). **None is pinned.** The two plik parameters pinned in earlier work (`xi_sz_cib`, `ksz_norm`) do not exist in CamSpec.

## 1. Minimiser audit and fix (`src/quadfit_min.py`)

**The bug.** The evaluator built x = c + Lz, clipped x to the prior bounds, and evaluated the clipped point, but the regression used the unclipped z.

**Fixes:**
- Proposals are drawn from the Gaussian truncated to the bounds (rejection sampling), so the regression always uses the evaluated coordinates.
- Steps minimise the fitted quadratic subject to the bounds and a trust region (SLSQP), so boundary minima are found as such and centres stay feasible.
- Failed evaluations above 5% and rank-deficient designs raise `QuadfitError`, which is recorded as a failure file.
- The objective is evaluated directly at the reported point, and both the fitted and direct values are saved.
- The uncertainty of the fitted minimum is the regression prediction variance at the minimiser. This includes the curvature terms (envelope theorem).
- A rejected final fit gives `converged = False`. It never falls back silently.

**Second problem, found in the first production attempt.** At f_cl = 1, u = 10⁻⁴ the centre random-walked for 8 iterations without converging.
- Diagnosis (`src/mix_noise_diag.py`, 408 evaluations): the noise is entirely in CamSpec high-ℓ.
  - rms 0.10, with a robust σ of 0.085.
  - Heavy tails: a few evaluations per hundred lie 4–5σ out.
  - The convergence tolerance (0.01) was below the regression error of the predicted decrease (0.03–0.05).
- Fixes:
  - robust refits that drop residuals beyond 4 robust σ;
  - a convergence threshold of max(0.01, 2 × SE of the fitted value);
  - an anomalously noisy iteration is resampled without moving or re-whitening;
  - re-whitening may change the sampling variance by at most 2× per iteration (previously 4×).
- The old minimiser had the same weakness. This probably explains the two old CamSpec points (u = 1.5×10⁻⁴ and 2×10⁻⁴) whose final fits were rejected.

**Tests** (`tests/test_quadfit_min.py`, 13, all pass): interior minimum; a strongly correlated objective with a poor diagonal starting covariance; a boundary minimum, checked against the analytic constrained minimum, with no evaluation outside the bounds; proposals crossing bounds near an interior minimum; a demonstration that the old clipped regression is biased; noisy objectives at the plik (0.03) and CamSpec (0.15) noise levels with 15 parameters; heavy-tailed noise from a 3σ start over 5 seeds; Monte-Carlo validation of the fitted-minimum error; failed evaluations; non-convergence reporting; and a degenerate design.

**Effect on the old f_cl = 1 result.**

| u | Old Δχ² | New Δχ² (relative to own u = 0) |
|---|---|---|
| 5×10⁻⁵ | 1.05 | 1.15 |
| 10⁻⁴ | 2.18 | 2.42 |
| 1.5×10⁻⁴ | 3.74 | 3.76 |
| 2×10⁻⁴ | 5.27 | 5.35 |

- **Nominal 2.71 crossing:** 1.17×10⁻⁴ before, **(1.12 ± 0.02)×10⁻⁴** now, a −4% change.
- **Cause:** not clipping. The old best fits were ≥5σ from every bound, so clipping rarely fired. The cause is the old fits' sensitivity to outlying evaluations; the u = 10⁻⁴ point had fit noise 0.23.
- **Repeatability of the new fits** (other seeds, plus an independent start from ΛCDM):

  | Point | Repeat values (−log P) | Spread |
  |---|---|---|
  | u = 0 | 5468.877, 5468.888 | |
  | u = 10⁻⁴ | 5470.086, 5470.089, 5470.070 | |
  | u = 1.5×10⁻⁴ | 5470.759, 5470.761 | |
  | pooled | | 0.008 rms in −log P |

  This is consistent with the regression errors (0.01–0.015), and the −4% shift is larger than the scatter.
- **Downstream:** the old profile limit (1.17×10⁻⁴) and the headline MCMC (1.69×10⁻⁴, which does not use the minimiser) are **not** invalidated. The profile limit should be quoted as 1.12×10⁻⁴ from now on.

## 2. Validation before production (`src/mix_validate.py`, `src/mix_direct_scan.py`)

1. **Build identity.**
   - `select_class_build()` puts `build/classy_patched` first on sys.path and raises if `classy` resolves elsewhere.
   - Cobaya is given `path: global`.
   - Every result records the classy path, version and patch SHA-256.
   - Each worker's loaded module is checked against the expected one after the fit.
   - Stock runs load `.venv/.../classy`; patched runs load `build/classy_patched/classy`.
2. **Patched Thomson = stock:** −log P is identical to the printed precision (difference 0.0) at f = 1, u = 10⁻⁴ and f = 0.1, u = 10⁻³, with explicit and default coefficients.
3. **f_cl = 0 against ΛCDM:** the same model, since no idm species is created.
4. **Fractions at u = 0** (f = 1, 0.5, 0.1, 0.03, 0.01) at a common point: within ±0.009 of ΛCDM in −log P, which is numerical noise. The fitted u = 0 minima agree with ΛCDM to within 0.02.
5. **Both kernels at u → 0:** at u = 10⁻⁹, Thomson and isotropic differ from u = 0 by ≤ 0.004.
6. **Fixed-parameter scan:** u from 10⁻⁶ to 10, all fractions, both kernels, at the old ΛCDM best fit.
   - CLASS runs without failure up to u = 0.1 (f = 1), 0.18 (0.5), 1 (0.1), 3.2 (0.03) and 10 (0.01), far beyond the couplings that matter.
   - This scan set the coarse grids.
7. **Tight coupling:** leaving tight coupling 10× earlier changes χ² by ≤ 0.05 over each fraction's relevant range. CLASS's first-order scheme gives ≤ 0.43 (f = 0.01 at u = 0.018), and ≤ 0.06 in the kernel difference.
8. **Strong-coupling points, rechecked at the profile best fits** (`src/mix_strong_coupling_check.py`). At f = 1 (2×10⁻⁴), 0.1 (4×10⁻³, 6×10⁻³), 0.03 (0.02, 0.045) and 0.01 (0.06, 0.135, 0.2), for both kernels, averaged over the best fit and 4 nearby points:
   - **TCA off early:** ≤ 0.04.
   - **First-order TCA:** −0.15 to +0.10.
   - **Precision up** (`l_linstep 10`, sampling step 0.05, tolerance 10⁻⁶): −0.12 to −0.26. That is almost entirely the same offset as at ΛCDM (−0.154 ± 0.010), so profile differences shift by ≲ 0.1.
   - **None of these moves the kernel difference**, which is several χ² at the strongest couplings. That difference is physical, not an artefact of the patch's approximations.

## 3. Production

**Commands:**
```bash
scripts/run_detached.sh results/mix/logs/production.log scripts/run_mix_production.sh   # coarse grids, repeats
scripts/run_detached.sh results/mix/logs/extend.log scripts/run_mix_extend.sh           # extended grids
uv run python src/profile_mix_analyse.py        # -> results/mix/camspec_npipe/summary.{json,md}, figures/mix_*.png
uv run python src/mix_strong_coupling_check.py  # -> results/mix/validation/strong_coupling_check.json
uv run python src/profile_mix.py status         # table of every point
```

**Results layout.** One directory per point: `results/mix/camspec_npipe/f<f>/<kernel>/u<u>[_<tag>]/result.json`. Each file holds:
- the complete Cobaya input and CLASS extra args;
- the build identity;
- the start point, its origin and the covariance;
- best-fit parameters, fitted and direct −log P, and the per-likelihood χ² and derived parameters at the best point;
- convergence history, active bounds and status reasons.

Failures go to `failure_*.json`.

**Totals:**
- 86 fits converged: the coarse and extended grids, both kernels, plus 5 repeats.
- One fit did not converge: f = 0.03 isotropic, u = 3×10⁻³, first attempt. It is kept as `u3.000e-03_notconverged1` and was refitted with a new seed.
- No CLASS failures inside the scanned ranges.
- Runs took 15–64 min each, about 40 h of local compute in total; no cloud resources were used.

**Grids:**
- Chosen from the fixed-parameter scan, not from 1/f_cl.
- Extended upward wherever the coarse profile was still falling or had not crossed.
- Starting points: neighbour warm starts.
- Isotropic fits are "mirrored": the same start, covariance and seed as Thomson, so the sample points coincide and the noise largely cancels in the kernel difference.

## 4. Results

### Nominal fixed-fraction profile thresholds

These are the PR4 CamSpec + 2018 low-ℓ thresholds.
- Δχ² is measured from each fraction's own profile minimum over the scanned u, which is not assumed to be at u = 0.
- The uncertainties come from the regression errors plus the pooled repeat scatter, propagated through the local slope.
- These are not joint (f_cl, u) confidence regions, and no frequentist coverage is claimed.

| f_cl | Kernel | u at min | u (Δχ² = 2.71) | Σ/Q min [g cm⁻²] | u (Δχ² = 3.84) | Scanned to | Δχ²(min) − ΛCDM |
|---|---|---|---|---|---|---|---|
| 1 | Thomson | 0 | (1.12 ± 0.02)×10⁻⁴ | 2.4×10⁶ | (1.52 ± 0.01)×10⁻⁴ | 5×10⁻⁴ | +0.02 |
| 1 | isotropic | 0 | (1.16 ± 0.02)×10⁻⁴ | 2.3×10⁶ | (1.56 ± 0.02)×10⁻⁴ | 5×10⁻⁴ | +0.02 |
| 0.5 | Thomson | 0 | (2.66 ± 0.04)×10⁻⁴ | 1.0×10⁶ | (3.66 ± 0.03)×10⁻⁴ | 8×10⁻⁴ | +0.03 |
| 0.5 | isotropic | 0 | (2.74 ± 0.04)×10⁻⁴ | 9.8×10⁵ | (3.76 ± 0.03)×10⁻⁴ | 8×10⁻⁴ | +0.03 |
| 0.1 | Thomson | 0 | (4.06 ± 0.04)×10⁻³ | 6.6×10⁴ | (5.06 ± 0.04)×10⁻³ | 9×10⁻³ | +0.02 |
| 0.1 | isotropic | 0 | (4.01 ± 0.03)×10⁻³ | 6.7×10⁴ | (4.87 ± 0.03)×10⁻³ | 9×10⁻³ | +0.02 |
| 0.03 | Thomson | 0.02 | (3.79 ± 0.02)×10⁻² | 7.1×10³ | (4.30 ± 0.02)×10⁻² | 0.045 | −1.91 |
| 0.03 | isotropic | 0.0135 | (3.77 ± 0.02)×10⁻² | 7.1×10³ | (4.28 ± 0.02)×10⁻² | 0.045 | −1.83 |
| 0.01 | Thomson | 0.06 | (1.48 ± 0.01)×10⁻¹ | 1.8×10³ | (1.72 ± 0.01)×10⁻¹ | 0.2 | −2.56 |
| 0.01 | isotropic | 0.06 | **not reached** | – | not reached | 0.2 | −2.88 |

- Σ/Q = 268/u g cm⁻²; Q remains explicit, since a kernel does not fix Q.
- The quoted u uncertainties are numerical only and assume linear interpolation between grid points. The grids are coarse (steps of 1.3–1.5 in u), and interpolation adds a few per cent of systematic error that is not included. Quote two significant figures.
- **f_cl = 0.01 isotropic:** the profile reaches Δχ² = 1.8 above its minimum at the largest scanned coupling, u = 0.2, and is still rising. No upper limit is claimed. The reliable calculation continues beyond u = 0.2 (CLASS ran to u = 10 in the fixed scan), but the grid stops there.

### How much the CDM weakens the coupling constraint

At the nominal 2.71 thresholds, f_cl·u_thr is:

| f_cl | f_cl·u_thr | Relative to f_cl = 1 |
|---|---|---|
| 1 | 1.1×10⁻⁴ | 1 |
| 0.5 | 1.3×10⁻⁴ | ×1.2 |
| 0.1 | 4.1×10⁻⁴ | ×3.6 |
| 0.03 | 1.1×10⁻³ | ×10 |
| 0.01 | 1.5×10⁻³ (Thomson) | ×13 |

- The threshold on u therefore grows **faster than 1/f_cl**: by 2.4× from f = 1 to 0.5, 36× to f = 0.1, 340× to f = 0.03, and 1300× to f = 0.01.
- In surface density, the minimum Σ/Q falls from 2.4×10⁶ g cm⁻² (f = 1) to about 1×10⁶ (0.5), 7×10⁴ (0.1), 7×10³ (0.03) and 2×10³ g cm⁻² (0.01, Thomson only).
- **Why faster than 1/f_cl:**
  - The drag suppresses only the clump perturbations. The CDM keeps growing and pulls the clumps along after decoupling, so σ₈ falls less at a given f_cl·u. For example, at f_cl·u = 10⁻⁴, σ₈ = 0.772 at f = 0.5 against 0.768 at f = 1.
  - Below f ≈ 0.1 the observable effect stops being dominated by the σ₈/lensing suppression (next subsection).

### The dip below ΛCDM at small fractions

For f_cl ≤ 0.03 the profile first falls below ΛCDM:
- by Δχ² = −1.9 at f = 0.03 (u ≈ 0.014–0.02);
- by −2.6 (Thomson) or −2.9 (isotropic) at f = 0.01 (u ≈ 0.06).

It then rises steeply. The decomposition at the minimum is:
- **Low-ℓ TT:** −1.4 to −2.0.
- **CamSpec high-ℓ:** 0 to −1.7.
- **Low-ℓ EE:** about 0.

Along the profile:
- n_s rises steadily, from 0.961 to 0.97–0.99;
- ln A_s and θ_s rise;
- τ, ω_b and σ₈ barely move (σ₈ ≈ 0.80–0.82 at f = 0.01).

**Interpretation.**
- A small, strongly coupled clump component modifies the damping tail, so the fit can raise n_s. A higher n_s lowers the large-scale power, which the known low-ℓ TT deficit of Planck favours.
- The rise beyond the minimum is the high-ℓ likelihood objecting (up to +4.5 to +5.2 at the largest couplings).
- **This is not evidence for clumps:**
  - the gain is Δχ² ≈ 2–3 for two extra parameters (f_cl, u);
  - it was found after looking;
  - it is driven by the same low-ℓ feature that other one-parameter extensions (running n_s, A_L) use.
- It is a degeneracy direction of the model, and the reason these profiles have non-zero minima.

### Kernel dependence

- **Thresholds:** isotropic minus Thomson at the 2.71 threshold is **+3% (f = 1), +3% (0.5), −1% (0.1) and −1% (0.03)**. The independent-fit errors on these differences are 1–3%, so the kernel changes the thresholds by ≲ 3% wherever both crossings are reliable.
- **Profiles:** at f ≥ 0.1 the profiles differ by ≤ 0.2 in Δχ² below the thresholds. Isotropic is almost always slightly lower, from the absence of clump-generated polarization (M6).
- **Strong coupling at small f:** the kernel matters.
  - At f = 0.01 the Thomson profile turns up much faster: Δχ² relative to ΛCDM at u = 0.2 is +2.6 for Thomson against −1.1 for isotropic, which is why the isotropic profile has no threshold in range.
  - At identical parameters the kernel difference is several χ² there (strong-coupling check above). It comes from CamSpec high-ℓ, where Thomson clump scattering generates polarization.
  - It is not moved by the TCA or precision variants.
  - **So the kernel does not change the conclusion for f_cl ≥ 0.03, but at f_cl ≈ 0.01 the upper limit on u exists only for the Thomson kernel within the scanned range.**

### Fitted parameters

ω_dark does not depend on f_cl itself:
- At u = 0 it is 0.1199 for every fraction.
- Along each profile it rises slightly, by about 1σ at most. At f ≥ 0.5 it tracks f_cl·u, and so do θ_s and ω_b.

At f ≥ 0.5 the coupling acts through σ₈ (from 0.806 to 0.74 at the largest u). H₀ stays at 67.1 throughout (`figures/mix_parameters_camspec_npipe.png`).

## Figures

- `figures/mix_profiles_camspec_npipe.png`: profile Δχ² against f_cl·u, one panel per fraction, both kernels.
- `figures/mix_vs_lcdm_camspec_npipe.png`: all fractions against the common ΛCDM reference.
- `figures/mix_boundary_camspec_npipe.png`: nominal fixed-fraction thresholds in (f_cl, u), with the Σ/Q axis.
- `figures/mix_parameters_camspec_npipe.png`: best-fit parameters along the Thomson profiles.

## Limitations

- The grids are coarse, and crossings are linearly interpolated.
- f = 0.01 isotropic has no threshold up to u = 0.2.
- Fixed-fraction profiles only; there is no joint (f_cl, u) analysis.
- Single likelihood combination, no lensing reconstruction.
- Linear theory (`non_linear: none`).
- The strong-coupling regime (f ≤ 0.03, u ≳ 0.01) is where post-recombination clump scattering becomes comparable to Thomson scattering on the residual free electrons. There the kernel is a real model choice.

## 5. Two-dimensional grid in (f_cl, u) (2026-10-06 to 10-09)

**Purpose:** map the shape of the small-fraction valley without any prior (`scripts/run_mix_grid2d.sh`).

**Grid:**
- New fractions 0.05, 0.02, 0.005, 0.003 and 0.001, with couplings at f_cl·u = 2, 4, 6, 9 and 13.5 ×10⁻⁴ plus u = 0. f_cl = 0.01 isotropic was extended to u = 0.3 and 0.45.
- Isotropic first, then Thomson on the same couplings. For f_cl ≤ 0.005 the isotropic grid continues to f_cl·u = 3×10⁻³; the Thomson grid stops at 1.35×10⁻³.

**Fits:** 63 fits, all converged except f_cl = 0.01 isotropic at u = 0.45, which is excluded and not needed because u = 0.3 brackets that threshold.

**Numerical checks at the extreme isotropic points** (f_cl = 0.005 at u = 0.4 and 0.6; 0.003 at 0.67 and 1; 0.001 at 0.9 and 2; `results/mix/validation/strong_coupling_check_grid2d.json`):
- TCA off early: ≤ 0.02.
- First-order TCA: −0.05 to −0.12.
- Precision up: −0.11 to −0.15, the same offset as at ΛCDM (−0.154).
- The Thomson-minus-isotropic difference at the isotropic best fits is 3–39 in χ², unchanged by these variants. It is physical.

**Results** (Δχ² of the profile minimum relative to the ΛCDM fit, and the nominal fixed-fraction threshold):

| f_cl | Kernel | u at min | Min − ΛCDM | u at Δχ² = 2.71 above own min | Scanned to |
|---|---|---|---|---|---|
| 0.05 | Thomson | 0.008 | −0.49 | 1.9×10⁻² | 0.027 |
| 0.05 | isotropic | 0.008 | −0.27 | 1.8×10⁻² | 0.027 |
| 0.02 | Thomson | 0.03 | −2.61 | 6.2×10⁻² | 0.0675 |
| 0.02 | isotropic | 0.03 | −2.64 | 6.6×10⁻² | 0.0675 |
| 0.01 | isotropic | 0.06 | −2.88 | 0.28 (Σ/Q ≈ 9.6×10² g cm⁻²) | 0.3 |
| 0.005 | Thomson | 0.08 | −2.08 | not reached | 0.27 |
| 0.005 | isotropic | 0.4 | −2.75 | not reached | 0.6 |
| 0.003 | Thomson | 0.133 | −1.53 | not reached | 0.45 |
| 0.003 | isotropic | 0.67 | −2.48 | not reached | 1 |
| 0.001 | Thomson | 0.2 | −0.69 | not reached | 1.35 |
| 0.001 | isotropic | 0.9 | −1.17 | not reached | 3 |

**Shape of the valley.**
- **Isotropic:** a broad region at Δχ² ≈ −2.5 to −2.9 from f_cl ≈ 0.02 to 0.003. Its best coupling moves to lower Σ/Q as f_cl falls (Σ/Q ≈ 9×10³ at 0.02, ~4×10² at 0.003), and it becomes shallower below f_cl ≈ 0.003.
- **Thomson:** the valley is narrower and shallower below f_cl = 0.01, because clump-generated polarization penalises strong coupling.

**Limits.**
- For f_cl ≤ 0.005 neither kernel reaches a threshold within the scanned range. That range runs down to Σ/Q ≈ 1.3×10² (isotropic, f_cl = 0.001) and ≈ 2×10² (Thomson, f_cl = 0.001). No limit is claimed there.
- Isotropic f_cl = 0.01 now has a threshold at u = 0.28 (it was "not reached" at u ≤ 0.2).

**Significance of the valley** (relative to ΛCDM, Δχ² up to 2.9):
- **Two fitted parameters:** p ≈ 0.23, about 1.2σ.
- **If u were fixed in advance (one parameter, boundary at f_cl = 0):** p ≈ 0.044, about 1.7σ. That is an upper bound, since u was in fact searched.
- The look-elsewhere effect over the valley lowers it further.
- Not evidence for clumps. It is the degeneracy with the Planck low-ℓ deficit described above.

**Walker & Wardle clouds** (M ≲ 10⁻³ M☉, R ≈ 1–3 AU, so Σ ≈ 3×10²–3×10³ g cm⁻²; Walker & Wardle 1998, ApJ 498, L125):
- With the isotropic kernel, this band lies in the region fitting as well as or slightly better than ΛCDM for f_cl ≈ 0.001–0.01.
- It is excluded (beyond the fixed-fraction threshold) only for f_cl ≳ 0.01–0.02.
- With Thomson, f_cl = 0.01 is excluded below Σ/Q ≈ 1.8×10³. At smaller fractions the band reaches about the ΛCDM level within the computed range.
- This assumes the clouds exist through the acoustic epoch; formation, survival and BBN are not addressed.

**Figure:** `figures/mix_contour_camspec_npipe.png` (`src/mix_contour.py`).
- Δχ² relative to ΛCDM in (f_cl, Σ/Q) for both kernels, with ω_cl on the top axis.
- Marked: the fixed-fraction thresholds, the Walker–Wardle band and the present-day masses in stars and cold gas (approximate shares of the cosmic baryons: 6% and 1.5%).
- Interpolation runs along lines of constant f_cl·u between fitted fractions. The u = 0 fits are set to Δχ² = 0 (they equal ΛCDM to within ±0.03 numerically). Regions beyond the fitted couplings are blank.

**Next:** an MCMC with free f_cl and u (log priors) to give the posterior shape, planned after this grid, needs the cloud VM and the user's go-ahead.
