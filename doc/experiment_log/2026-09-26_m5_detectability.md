# 2026-09-26 — M5: approximate detectability with realistic noise, and Fisher forecasts

**Goal:** Milestone M5 (plan §6.4 and §15, revision 3).
- Go beyond the ideal cosmic-variance statistic of M3: add instrumental noise, sky fractions, multipole ranges and the full TT/TE/EE covariance for a Planck-like and a future (Simons Observatory) configuration.
- Forecast the bound on u when the six standard parameters are free (Fisher), and compare the Planck-like forecast with the published Planck bound u < 2.25×10⁻⁴.

## Configuration

- Stock CLASS v3.4.0 @ `64bbab7`, `configs/baseline.yaml`, f_cl = 1, Thomson kernel, m_idm = 10²⁷ eV.
- Spectra computed to ℓ = 3000 (`l_max_scalars = 3500`) for the SO configuration.
- **Coupling grid:** the M3 u values up to 10⁻¹, at fixed h.
- **Fisher runs:** about u = 0 (f_cl = 1, `omega_cdm = 0`) in {u, ω_b, ω_dm, 100θ_s, n_s, ln 10¹⁰A_s, τ}.
  - Central differences with steps 2×10⁻⁴, 2×10⁻³, 4×10⁻⁴, 4×10⁻³, 0.02 and 5×10⁻³, and half of each.
  - One-sided steps in u (u ≥ 0) of 3×10⁻⁶, 10⁻⁵ (default) and 3×10⁻⁵.

## What was run

```bash
uv run python src/m5_detectability.py run       # 61 CLASS runs, ~1 min
uv run python src/m5_detectability.py analyse   # -> results/scans/m5/summary.json, figures/m5_*.png
```

New module: `src/m5_noise.py`. It provides noise spectra, a `Config` made of sky patches, the Gaussian Fisher matrix and Δχ².

## Noise model and sources (each checked against the paper itself)

- **Noise spectrum:** N_ℓ = (Δ_T)² exp[ℓ(ℓ+1)θ²/8 ln 2]·[1 + (ℓ/ℓ_knee)^α] per channel, with channels combined by inverse variance. It is converted to CLASS units by dividing by T_cmb².
- **Planck channels:** Planck 2018 I (A&A 641, A1), Table 4.
  - 100/143/217 GHz beams: 9.66′, 7.22′, 4.90′.
  - Temperature noise: 1.29, 0.55, 0.78 µK·deg.
  - Polarization noise: 1.96, 1.17, 1.75 µK·deg.
  - Combined, this is ≈ 36 µK·arcmin-equivalent at ℓ = 1000 including the beam. N_TT/C_TT is 0.017 at ℓ = 1000 and 1.6 at ℓ = 2000; N_EE/C_EE = 1.9 at ℓ = 1000.
- **Planck-like sky fractions and ranges:** the Planck forecast set-up of the SO science-goals paper (Ade et al. 2019, JCAP 02, 056), which that paper reports as reproducing the Planck 2018 errors approximately.
  - T/E at 2 ≤ ℓ ≤ 29 with f_sky = 0.8;
  - TT/TE/EE at 30 ≤ ℓ ≤ 2500 with f_sky = 0.6;
  - κκ at 8 ≤ L ≤ 400 with f_sky = 0.6.
  - For reference, the actual Planck 2018 likelihood (Planck 2018 V, A&A 641, A5) uses the T66/T57/T47 and P70/P50/P41 masks, TT to ℓ = 2508 and TE/EE to ℓ = 1996.
- **SO LAT:** Ade et al. 2019, Table 1, baseline.
  - 93 GHz: 2.2′, 8.0 µK·arcmin; 145 GHz: 1.4′, 10 µK·arcmin; f_sky = 0.4; polarization noise √2 higher.
  - 1/f model (their Eq. 1): polarization ℓ_knee = 700, α = −1.4, N_red = N_white.
  - Temperature: ℓ_knee = 1000, α = −3.5. The paper gives N_red in µK²·s, which needs survey details it doesn't tabulate, so **we approximate N_red = N_white at ℓ_knee**. A white-noise-only case is also shown.
  - The SO configuration is "SO + Planck":
    - Planck low ℓ (2–29, f_sky = 0.8);
    - SO and Planck noise combined on the SO sky (f_sky = 0.4, ℓ = 30–3000);
    - Planck alone on the rest of the Planck sky (f_sky = 0.2, ℓ = 30–2500).
- **Covariance:** per multipole and sky patch, Δχ² and the Fisher matrix use the 2×2 covariance [[TT + N_T, TE], [TE, EE + N_E]] (the full Gaussian TT/TE/EE covariance) via (2ℓ+1) f_sky/2 · Tr[C⁻¹ΔC C⁻¹ΔC]. Single-spectrum versions use Var = 2(C+N)²/((2ℓ+1)f_sky), and (TE² + (TT+N)(EE+N))/((2ℓ+1)f_sky) for TE.
- **Check:** the cosmic-variance full-sky configuration reproduces the M3 single-spectrum Δχ² exactly (TT 94.24, EE 82.45, TE 45.73 at u = 10⁻⁴).

## Lensing

- **φφ reconstruction is excluded from the baseline numbers.** I have not verified a Planck lensing-reconstruction noise curve, so any number would be unsourced.
- **Upper envelope instead:** φφ is added with cosmic variance only over 8 ≤ L ≤ 400 at f_sky = 0.6 (the SO paper's Planck range). This ignores reconstruction noise and the overlap with lensed-spectrum information, so it can only overstate what φφ adds.
- **Lensed-spectrum information** is always included (lensed TT/TE/EE), and its effect is shown by rerunning with unlensed spectra.

## Results

### Δχ² at fixed cosmology (TT+TE+EE, relative to u = 0)

u at which Δχ² reaches 4:

| Configuration | joint | TT | TE | EE | joint, unlensed spectra | joint + φφ envelope |
|---|---:|---:|---:|---:|---:|---:|
| cosmic variance, full sky, ℓ ≤ 2500 | 1.3×10⁻⁵ | 2.3×10⁻⁵ | 2.7×10⁻⁵ | 2.2×10⁻⁵ | 1.2×10⁻⁵ | — |
| **Planck-like** | **6.8×10⁻⁵** | 8.2×10⁻⁵ | 1.5×10⁻⁴ | 2.2×10⁻⁴ | 6.6×10⁻⁵ | 3.1×10⁻⁵ |
| SO LAT baseline + Planck | 9.2×10⁻⁶ | 9.7×10⁻⁶ | 5.1×10⁻⁵ | 5.1×10⁻⁵ | 2.3×10⁻⁵ | — |
| same, white noise only | 9.1×10⁻⁶ | 9.7×10⁻⁶ | 4.9×10⁻⁵ | 4.8×10⁻⁵ | 2.2×10⁻⁵ | — |

Figure: `figures/m5_delta_chi2.png`.
- For Planck, TT carries most of the sensitivity; EE is noise-limited.
- At small u, SO + Planck exceeds the "cosmic variance, ℓ ≤ 2500" curve. The reason is that it extends to ℓ = 3000, and the extra damping-tail multipoles carry lensing information.
- For SO, the lensed spectra are worth 2.5× in u at fixed cosmology (9.2×10⁻⁶ lensed against 2.3×10⁻⁵ unlensed).
- The atmospheric 1/f noise has almost no effect (< 3%), because Planck covers the low multipoles.

### Fisher forecasts (fiducial u = 0; the six standard parameters free)

| Configuration | σ(u) conditional | **σ(u) marginalised** | u₉₅ (1.645σ) | u₉₅ (1.96σ, half-Gaussian) | Σ/Q min at u₉₅ [g/cm²] | u₉₅ / published |
|---|---:|---:|---:|---:|---:|---:|
| cosmic variance, full sky | 6.2×10⁻⁶ | 1.0×10⁻⁵ | 1.7×10⁻⁵ | 2.0×10⁻⁵ | 1.6×10⁷ | 0.07 |
| **Planck-like** | 3.1×10⁻⁵ | **6.7×10⁻⁵** | **1.1×10⁻⁴** | 1.3×10⁻⁴ | 2.4×10⁶ | **0.49** |
| Planck-like + φφ CV envelope | 1.4×10⁻⁵ | 3.4×10⁻⁵ | 5.6×10⁻⁵ | 6.6×10⁻⁵ | 4.8×10⁶ | 0.25 |
| SO LAT baseline + Planck | 4.6×10⁻⁶ | 7.4×10⁻⁶ | 1.2×10⁻⁵ | 1.4×10⁻⁵ | 2.2×10⁷ | 0.05 |
| Planck-like, unlensed spectra | 3.1×10⁻⁵ | 9.0×10⁻⁵ | 1.5×10⁻⁴ | 1.8×10⁻⁴ | 1.8×10⁶ | 0.65 |
| SO + Planck, unlensed spectra | 1.1×10⁻⁵ | 2.7×10⁻⁵ | 4.5×10⁻⁵ | 5.3×10⁻⁵ | 6.0×10⁶ | 0.20 |

- **Consistency:** the fixed-cosmology Δχ² = 1 point for Planck-like (3.3×10⁻⁵, fixed h) agrees with the conditional σ(u) (3.1×10⁻⁵, fixed θ_s) to 7%.
- **Degeneracies** (correlation of u with the others, Planck-like):
  - 100θ_s 0.77, ω_b 0.45, ω_dm 0.29; n_s, A_s and τ ≤ 0.1.
  - Marginalising over the standard parameters doubles σ(u).
  - Adding φφ breaks the θ_s degeneracy (correlation drops to 0.49), but u then correlates with ω_dm (0.78), since both change the lensing amplitude.
  - For SO + Planck the correlations are all ≤ 0.4.
  - Figure: `figures/m5_fisher_ellipses.png`.
- **Lensing information in the spectra** helps break the degeneracies: removing it (unlensed spectra) worsens the marginalised σ(u) by 1.35× for Planck-like, 2.0× for cosmic variance and 3.7× for SO.
- **Stability:**
  - Halving all parameter steps changes σ(u) by < 0.1%.
  - Changing the u step from 10⁻⁵ to 3×10⁻⁶ or 3×10⁻⁵ changes it by −2%/+3% (Planck-like) and −13%/+21% (SO). The response to u is mildly sublinear, which matters only at the SO sensitivity; quote the SO forecast as ±20%.
- **Validation of the Fisher set-up against Planck 2018** (Planck 2018 VI, A&A 641, A6, Table 2, TT,TE,EE+lowE, checked in the paper). Forecast (Planck-like) against published σ:
  - ω_b 0.00014 vs 0.00015; ω_c 0.0012 vs 0.0014; n_s 0.0032 vs 0.0044; 100θ 0.00041 (θ_s) vs 0.00031 (θ_MC).
  - ln 10¹⁰A_s 0.0069 vs 0.016 and τ 0.0034 vs ~0.0075.
  - The core parameters are reproduced within ~25%. τ and A_s are 2× too tight, because noise-only low-ℓ EE ignores the systematics that limit Planck's large-scale polarization.
  - **Test:** inflating the low-ℓ EE noise until σ(τ) = 0.0064–0.015 changes σ(u) by < 4% (u barely correlates with τ). So this idealisation doesn't affect the u forecast.

## Comparison with the published bound

The Planck-like forecast, u₉₅ = 1.1×10⁻⁴ (1.3×10⁻⁴ half-Gaussian), is **a factor ~2 tighter** than the published u < 2.25×10⁻⁴ (Stadler & Bœhm 2018). That is within the factor ~2–3 expected, so it does not trigger an investigation. Plausible contributors, **not quantified here**:
- the Fisher forecast is idealised: no foreground or calibration nuisance parameters, Gaussian likelihood, single f_sky, TE/EE to ℓ = 2500 against 1996 in Planck's likelihood;
- the published analysis used a different data release and likelihood, not checked here;
- its prior on u affects a one-sided bound; not checked here.

M7 (reproducing the bound with the real likelihood) will settle this.

## Outcome

- **Planck-like:** at fixed cosmology the coupling becomes detectable (Δχ² = 4) at u ≈ 7×10⁻⁵. With the standard parameters free, the forecast 95% limit is **u < 1.1×10⁻⁴ (Σ/Q > 2.4×10⁶ g/cm²)**, consistent with the published bound to a factor ~2. TT dominates, and the main degeneracy is with θ_s and ω_b.
- **SO LAT + Planck:** the forecast limit is **u < 1.2×10⁻⁵ (Σ/Q > 2.2×10⁷ g/cm²)**, about 10–20× tighter than Planck. The gain comes largely from lensing information in the damping tail. Atmospheric noise barely matters.
- **φφ reconstruction** could at most halve the Planck-like σ(u) (a cosmic-variance envelope).
- **Even the future CMB bound stays ~10³ times weaker in u than the small-scale bound (M12),** which, however, assumes clumps already in place by z ~ 10⁷.

## Follow-ups

- M7: reproduce the published bound with a real likelihood (Cobaya), which also resolves the factor-2 difference.
- For a better φφ treatment, find and verify a published Planck lensing-reconstruction noise curve, or compute a quadratic-estimator N⁽⁰⁾.
- In M8/M9, use profile likelihoods; the Fisher forecast gives the expected scale and the most degenerate directions (θ_s, ω_b, ω_dm).
