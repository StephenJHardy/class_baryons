# 2026-09-26 — M6: Experiment B, minimal kernel check

**Goal:** Milestone M6 (plan revision 2, §7). Patch CLASS so that the clump part of the photon collision term can use a kernel other than Thomson. Validate the patch against stock CLASS. Compare a perfect absorber (isotropic re-emission in the clump frame) with Thomson at equal momentum transfer u, isolate the low-ℓ EE term seen in M3, and reach a verdict.

## Configuration

- CLASS v3.4.0 @ `64bbab7` plus `class_patches/idm_g_kernel.patch`.
  - Authored by `scripts/make_kernel_patch_edits.py`, applied by `scripts/setup_class.sh` to a separate tree `class_patched/`.
  - Installed to `build/classy_patched/` and selected with `CLASS_BUILD=patched`. The stock build is untouched and remains the default.
- `configs/baseline.yaml` (including `l_linstep: 20`); f_cl = 1, m_idm = 10²⁷ eV, n_index_idm_g = 0, fixed h = 0.6737666.

### What the patch does

It adds two input parameters, read only when `u_idm_g > 0`:

| Parameter | Symbol | Thomson (default) | Absorber |
|---|---|---|---|
| `idm_g_quadrupole_coefficient` | c_T | 1 | 0 |
| `idm_g_polarization_coefficient` | c_E | 1 | 0 |

Let a = κ′ (Thomson rate) and b = μ′ (clump rate).
- **Photon hierarchy.** Damping of every moment with ℓ ≥ 2 keeps the full rate (a + b). The regeneration or source terms use (a + c·b): the P0 term in the temperature ℓ = 2 equation uses c_T, and the P0 terms in the polarization hierarchy use c_E. The ℓ = 1 momentum exchange is unchanged, since the kernel doesn't affect it.
- **Line-of-sight sources.** The Π term in temperature (t2) and the E source (p) use g·(a + c·b)/(a + b) in place of the total visibility g.
- **Tight coupling.** CLASS's Thomson closed forms (shear = 16/45·τ₂·θ with τ₂ = 1/(a + b), and P0 = (5/8)·s₂·shear) are generalised exactly. Solving the quasi-static ℓ = 2 temperature and polarization equations gives:
  - P0 = 2·s₂·shear/(8 − 4.8r), with r = (a + c_E·b)/(a + b);
  - shear = (4/15)·θ/D, with D = (a + b) − 0.8·(a + c_T·b)/(8 − 4.8r);
  - so τ₂ = 0.75/D, whose analytic time derivative feeds the second-order terms.
  All of these reduce exactly to CLASS's expressions when c_T = c_E = 1.
- **Bit-identity for Thomson.** Every modified term is written as the stock expression plus (regen − rate)·(source), and regen is set to exactly the stock rate when the coefficient is 1. The closure helper is called only when a coefficient differs from 1. So the Thomson path is **bit-identical** to stock.
- **Not patched:** vector and tensor modes (the project doesn't compute them). The 11/6 second-order tight-coupling correction keeps its Thomson form; the clump contribution is ≲10⁻⁴ of the rate while tight coupling is on.

## What was run

```bash
scripts/setup_class.sh                                      # clones, patches, builds both
CLASS_BUILD=stock   uv run python src/m6_kernel.py run      # 34 models
CLASS_BUILD=patched uv run python src/m6_kernel.py run      # 92 models, ~2 min
uv run python src/m6_kernel.py analyse                      # -> results/scans/m6/summary.json, figures/m6_absorber_vs_thomson.png
OMP_NUM_THREADS=2 CLASS_BUILD=patched uv run python src/glitch_scan.py project
```

The patched set comprises:
- five kernel variants (Thomson with default coefficients, Thomson with explicit coefficients, absorber, c_T = 0 only, c_E = 0 only) at u = 10⁻⁶, 10⁻⁵, 10⁻⁴, 2.25×10⁻⁴, 10⁻³, 10⁻² and 10⁻¹;
- a fine grid (log₁₀u from −5.5 to −3 in steps of 0.1) for Thomson and the absorber;
- two tight-coupling variants at the Planck bound.

## A bug caught in the pipeline (outputs discarded)

The first patched run silently used the **stock** classy.
- **Symptom:** every model passing the new coefficients failed with "did not read input parameter(s)", while the default-coefficient "patched" models ran.
- **Cause:** `m6_kernel.py` imports `m3_analyse`, which imported `classy` at module level *before* `run_class` could put `build/classy_patched` on `sys.path`. Python reused the cached stock module, and each run's recorded `classy_path` pointed at `.venv`.
- **Fix:**
  - `run_class` now refuses to load if `classy` was already imported, and checks that the loaded module lives where `CLASS_BUILD` says it should.
  - `m3_analyse` and `run_baseline` get `Class` and `classy` from `run_class` instead of importing them directly.
- The first set of results was deleted and rerun. Every result below was checked to come from the intended build: `classy_path` is `build/...` for patched runs and `.venv/...` for stock runs.

## Outcome

### Validation: PASS (bit-identical)

The patched build with Thomson coefficients (defaults, or passed explicitly) reproduces stock CLASS **bit for bit** in all 15 comparisons: ΛCDM plus 7 couplings × 2. An earlier draft of the patch, which rewrote −rate·(x − s) as −rate·x + regen·s, differed from stock by ≤ 9×10⁻⁶ purely through floating-point rounding. It was replaced by the bit-identical form.

The ℓ-node glitch scan on the patched build (59 models, log₁₀u from −7 to −4) gives 3/59 flagged, worst 3.8×10⁻⁴. That is identical to stock: all three are TE at ℓ = 2.

The *kernel difference* itself is smooth in u. On the fine grid, the absorber − Thomson difference D scales linearly with u, and the largest jump in D/u between neighbours is ≤ 3.4×10⁻⁵ in any spectrum (noise level), 1.3% of D for EE.

### Absorber vs Thomson at equal u

| u | max ΔTT | max ΔTE | max ΔEE (ℓ<30 / ℓ≥30) | max Δφφ | Δχ²_CV kernel difference (TT / TE / EE) | Δχ²_CV of Thomson vs ΛCDM (TT / EE) |
|---:|---:|---:|---:|---:|---:|---:|
| 10⁻⁵ | 3.9×10⁻⁶ | 1.7×10⁻⁴ | 6.4×10⁻⁴ / 5.4×10⁻⁴ | 5×10⁻⁷ | 9×10⁻⁶ / 3×10⁻⁴ / 2.6×10⁻³ | 0.81 / 0.88 |
| 10⁻⁴ | 4.3×10⁻⁵ | 1.7×10⁻³ | 6.4×10⁻³ / 5.4×10⁻³ | 4×10⁻⁷ | 1.6×10⁻⁴ / 0.032 / 0.27 | 94 / 82 |
| 2.25×10⁻⁴ | 4.6×10⁻⁵ | 3.9×10⁻³ | 1.4×10⁻² / 1.2×10⁻² | 2×10⁻⁷ | 7.5×10⁻⁴ / 0.16 / 1.3 | 608 / 509 |
| 10⁻³ | 1.8×10⁻⁴ | 1.7×10⁻² | 6.5×10⁻² / 5.5×10⁻² | 1×10⁻⁶ | 0.020 / 3.1 / 26 | 1.7×10⁴ / 1.5×10⁴ |

Residuals are relative to ΛCDM, with TE relative to √(TT·EE); "max" is over 2 ≤ ℓ ≤ 2500, lensed.

- **TT, φφ and P(k) are kernel-independent.** Differences are ≤ 5×10⁻⁵ in TT at the Planck bound, at the noise floor in φφ, and ΔT² ≤ 1.2×10⁻⁴ in P(k). The drag, and hence everything gravitational, is unchanged, as expected.
- **EE (and TE) carry a real kernel dependence, linear in u.** At the Planck bound the absorber has less EE than Thomson by 1.4% at ℓ ≈ 15–30, 0.9% at ℓ = 50, 0.4% at ℓ = 100, 0.27% at ℓ = 300 and 0.06% at ℓ ~ 1200. It is negative everywhere, since the absorber creates no polarization.
- **Decomposition:** switching off only the polarization source (c_E = 0) reproduces the whole EE difference (6.4×10⁻³ at u = 10⁻⁴). Switching off only quadrupole regeneration (c_T = 0) gives 1.4×10⁻⁵. So **the kernel-dependent signal is polarization generated by clump scattering.** This confirms the M3 hypothesis for the low-ℓ EE deviation and extends it to all ℓ. The relative size grows where the clumps' share of the visibility function grows: in the tail after recombination, as x_e falls.
- **Effect on a bound.** What matters is each kernel's distance from ΛCDM, including the cross term with the signal. Absorber against Thomson, the ideal-CV Δχ² relative to ΛCDM changes by:
  - EE: −2.8% to +3.6% depending on u;
  - TE: ≤ 3%;
  - TT: ≤ 0.5%;
  - TT + TE + EE combined: ≤ 1.6%.
  Since Δχ² ∝ u², that shifts the ideal-CV threshold on u by **< 1%**.
- **Tight-coupling handling doesn't matter.** At the Planck bound, the first-order scheme and switching tight coupling off early both reproduce the absorber − Thomson difference to within 2×10⁻⁵ (TT) and 2×10⁻⁴ (EE).

### Verdict

The Thomson-kernel CMB constraint carries over to opaque clumps. Replacing the Thomson kernel with a perfect absorber at equal momentum transfer changes the ideal-CV detectability of u by less than 1%. The only kernel-dependent observable is polarization generated by the clumps themselves: a few-percent part of EE at CMB-relevant u, degenerate in practice with τ at low ℓ. The kernel matters for Σ mainly through Q: Σ_min = Q·268/u_max, with Q = 1 (absorber or specular) or 13/9 (Lambertian). As the plan specifies, the full B1–B3 programme is **not** triggered.

## Follow-ups

- C and the likelihood work (M7, M8) can use the Thomson kernel. The kernel enters the Σ mapping only through Q.
- A realistic-noise version of the kernel cross term belongs to M5. With Planck-level EE noise it will be smaller than the ideal-CV numbers above.
