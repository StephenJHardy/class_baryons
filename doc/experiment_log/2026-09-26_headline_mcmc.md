# 2026-09-26: Headline MCMC (PR4 CamSpec + 2018 low-ℓ, f_cl = 1, flat prior on u)

## Set-up

- **Likelihoods:** planck_2018_lowl.TT + planck_2018_lowl.EE + planck_NPIPE_highl_CamSpec.TTTEEE.
- **Parameters:** 6 cosmological + u + 9 CamSpec nuisance. u has a flat prior on [0, 10⁻³].
- **Sampler:** Cobaya mcmc with dragging (fast nuisance parameters, oversample_power 0.4) and learn_proposal. The initial proposal is the covariance from the CamSpec profile at u = 0, with σ_u = 5×10⁻⁵; starting points are near the profile best fit (`src/mcmc_run.py`).
- **Compute:** 24 independent chains without MPI on the spot VM (c3d-highcpu-180): 12 chains at 15 threads and 12 at 7 threads. CLASS ran at ~0.86 evaluations/s per chain. The run took about 6.5 hours of VM time, and the VM was stopped afterwards.
- **Watchdog:** `scripts/vm_mcmc_watchdog.sh` (restarts after preemption; chains resume).
  - The first launch did nothing. Its running-chain check (`ps | grep`) matched the remote `bash -c` command line, which contains the chain command text, so no chain was ever started.
  - Fixed by excluding shell processes (commit fdb13fd). See IMPLEMENTATION_NOTES.

## Convergence

- **Stopping rule:** Cobaya's per-chain R−1 (split batches) was still 0.3–2 when we stopped. It would take many more hours to reach its 0.01 target, so we judged convergence across chains:

| check | accepted (after 30% burn-in) | R−1 of mean(u) across chains | u95 |
|---|---:|---:|---:|
| +1 h | 5,100 | 0.32 | 1.8×10⁻⁴ |
| +2.5 h | 16,300 | 0.052 | 1.75×10⁻⁴ |
| +3.5 h | 18,400 | 0.040 | 1.73×10⁻⁴ |
| +5 h | 26,500 | 0.022 | 1.70×10⁻⁴ |
| final | 29,900 | 0.021 | 1.69×10⁻⁴ |

- **Split-half test:** u95 from the first and second halves of every chain gives 1.694 and 1.685×10⁻⁴ (< 1% apart).
- **Chain storage:** the chains (22 MB) are archived at `/media/stephen/astro/class_baryons/mcmc/` and are not committed. The summary is in `results/mcmc/camspec_npipe/mcmc_summary.json`.

## Result

- **u < 1.69×10⁻⁴ (95%), < 2.33×10⁻⁴ (99%)**, i.e. σ/M < 6.3×10⁻⁷ cm²/g and **Σ/Q > 1.6×10⁶ g/cm²** (95%).
- **Shape:** the posterior peaks at u = 0 and falls monotonically. Mean u = 6.0×10⁻⁵, s.d. 5.4×10⁻⁵.
- **Other parameters:** H0 = 67.15 ± 0.51, n_s = 0.9616 ± 0.0041, τ = 0.052 ± 0.008. σ8 = 0.784 ± 0.019 is pulled down by u.
- **Correlations with u:** σ8 −0.92, 100θ_s +0.73, ω_dm +0.28, ω_b +0.16, H0 −0.03, n_s −0.02.
- **Figures:** `figures/report_mcmc.png` (posterior plus scatter against correlated parameters), `figures/mcmc_u_camspec_npipe.png`, `figures/mcmc_triangle_camspec_npipe.png`.

## Profile against marginal

- **The two limits:** the profile 95% limit for the same likelihood is 1.17×10⁻⁴; the MCMC limit is 45% higher.
- **Cause:** Δχ²(u) is nearly linear, so the likelihood is roughly exponential in u. For an exponential, the 95% probability point lies at Δχ² = 2 ln 20 = 5.99, which is 2.2 times the u of the Δχ² = 2.71 crossing. For a quadratic Δχ² the ratio would be 1.19.
- **Check:** integrating exp(−Δχ²/2) from the profile over the flat prior gives 1.96×10⁻⁴ (CamSpec) and 1.82×10⁻⁴ (plik). The MCMC's 1.69×10⁻⁴ lies between the profile crossing and this value, so there is no inflating volume effect.
- **Published limits:** the MCMC agrees with the published Bayesian limits (Zhou+22, Planck 2018: 1.55×10⁻⁴; Stadler & Bœhm 2018, Planck 2015: 1.58×10⁻⁴). The earlier gap between our profile limit and the literature was therefore mostly profile against marginal.

## Interpretation

- **σ8:** the anticorrelation of u with σ8 (−0.92) shows directly that the constraint comes from suppressed structure, seen through CMB lensing.
- **θ_s:** the correlation with θ_s (+0.73) matches the M5 Fisher degeneracy. It is probably a small phase shift of the peaks from the weaker potentials.
- **H0:** u is uncorrelated with H0, so the coupling does not help the Hubble tension.
