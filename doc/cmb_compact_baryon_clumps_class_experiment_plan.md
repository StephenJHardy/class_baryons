# CMB Constraints on Compact Weakly Photon-Coupled Matter

## Experimental Plan Using CLASS

### Research question

This project asks:

> **Given a fraction \(f_{\rm cl}\) of the would-be dark matter in compact, pressureless objects, what photon momentum-transfer cross-section per unit mass can that population have while remaining compatible with the observed CMB?**

and, more sharply:

> **Does an optically thick macroscopic baryonic object have the same CMB collision operator as the particle-scattering model conventionally used to derive dark-matter–photon bounds? If not, how does the CMB constraint on its surface density \(\Sigma\) change when the correct macro-object radiative-transfer physics is used?**

The project deliberately does **not** model cloud formation, Big Bang nucleosynthesis, or detailed internal cloud hydrodynamics. It asks whether compact clumps already present before the CMB acoustic epoch can behave sufficiently like cold dark matter to avoid altering the observed anisotropies.

The public **CLASS** Boltzmann code already contains an interacting dark matter species (`idm`) with photon coupling (`u_idm_g` / `cross_idm_g`). That implementation follows the Stadler & Bœhm elastic DM–photon scattering formalism, i.e. a **Thomson-like angular collision operator**. It is therefore treated here as a **reference approximation** (a momentum-drag proxy), not as the physical model of an opaque cloud.

### Project structure

| Experiment | Question | CLASS modification? |
|---|---|---|
| **A** — Thomson-like clumps | Reproduce the established interacting-DM bound at \(f_{\rm cl}=1\) and translate it into \(\Sigma\). | No |
| **B** — Opaque-clump collision operator | Does replacing the particle-scattering kernel with a macro-object kernel materially change TT/TE/EE, and hence the \(\Sigma\) bound? | Yes (photon hierarchy) |
| **C** — Partial clumped fraction | Map the allowed \((f_{\rm cl},\Sigma/Q)\) space, with whichever kernel(s) B shows to matter. | Only if B requires it |
| **D** — Formation history | \(\Sigma(z)\), \(f_{\rm cl}(z)\), \(z_{\rm form}\). Only if warranted. | Yes |

Reproducing a Planck \(u\lesssim10^{-4}\) bound is established science and serves as validation. **Experiment B is the scientifically novel core.**

---

## 1. Physical model

Start from a standard \(\Lambda\)CDM cosmology, but split the usual dark matter density into two pressureless components:

\[
\omega_{\rm dm}
=
\omega_{\rm cdm}
+
\omega_{\rm cl},
\]

where:

- \(\omega_{\rm cdm}\) is perfectly collisionless cold dark matter;
- \(\omega_{\rm cl}\) is the compact-clump component;
- \(\omega_b\) is ordinary diffuse baryonic matter participating in standard recombination.

Define the compact-clump fraction

\[
f_{\rm cl}
\equiv
\frac{\omega_{\rm cl}}
{\omega_{\rm cdm}+\omega_{\rm cl}},
\qquad
\omega_{\rm cl}=f_{\rm cl}\omega_{\rm dm},
\qquad
\omega_{\rm cdm}=(1-f_{\rm cl})\omega_{\rm dm}.
\]

Use `omega_idm = omega_cl` explicitly in CLASS rather than `f_idm`, so that the bookkeeping remains transparent.

### Assumptions

The compact component is assumed to satisfy

\[
w_{\rm cl}=0,
\qquad
c_{s,\rm cl}\simeq 0,
\]

and the clumps are assumed to:

1. already exist at the earliest redshifts relevant to the CLASS calculation;
2. have a constant comoving abundance;
3. behave as pressureless matter;
4. couple to photons through an effective momentum-transfer cross-section;
5. be numerous enough to be treated as a **continuous fluid**.

On assumption 5: CLASS evolves a smooth component. For planetary-mass clouds the number of objects per cosmological perturbation volume is astronomically large, so the fluid approximation is safe. Sufficiently massive, rare objects would introduce Poisson (shot-noise) perturbations not represented by the standard equations; such masses are outside the scope of this plan.

These assumptions define a **phenomenological CMB experiment**, not a cloud-formation theory.

---

## 2. Interaction parameter and physical mapping

CLASS parameterises the photon interaction with

\[
u_{\gamma{\rm cl}}
=
\frac{\sigma_{\gamma{\rm cl}}}{\sigma_T}
\left(
\frac{M_{\rm cl}}{100\,{\rm GeV}}
\right)^{-1}.
\]

Because the Thomson kernel is forward–backward symmetric, \(\sigma\) here is equivalently the **momentum-transfer** cross-section. That is the quantity that must be carried across to any other kernel (Section 5).

Numerically, since \(\sigma_T/(100\,{\rm GeV}) = 3.73\times10^{-3}\ {\rm cm^2\,g^{-1}}\),

\[
u
\simeq
268
\left[
\frac{\sigma/M}
{\mathrm{cm^2\,g^{-1}}}
\right],
\qquad
\boxed{
\frac{\sigma}{M}
\simeq
3.73\times10^{-3}u
\ {\rm cm^2\,g^{-1}}
}
\]

For an optically thick spherical clump,

\[
\sigma_{\rm mt}\simeq Q\pi R^2,
\]

where \(Q\) is the momentum-transfer efficiency of the surface/radiative law. In the geometric-optics limit:

| Surface law | \(Q\) |
|---|---:|
| Perfect absorber, isotropic re-emission in clump frame | \(1\) |
| Specular reflection from a sphere | \(1\) |
| Diffuse (Lambertian) reflection from a sphere | \(13/9\) |

\(Q\) is therefore **not** a free fudge factor once the kernel is specified in Experiment B.

With surface density \(\Sigma\equiv M/\pi R^2\),

\[
\frac{\sigma_{\rm mt}}{M}
\simeq
\frac{Q}{\Sigma},
\qquad
\boxed{
u\simeq
\frac{268\,Q}
{\Sigma/({\rm g\,cm^{-2}})}
}
\]

Every result computed in \(u\) is also reported as \(\sigma/M\) and \(\Sigma/Q\).

| \(u\) | \(\sigma/M\) [cm\(^2\)/g] | \(\Sigma/Q\) [g/cm\(^2\)] |
|---:|---:|---:|
| \(10^{-2}\) | \(3.7\times10^{-5}\) | \(2.7\times10^4\) |
| \(10^{-3}\) | \(3.7\times10^{-6}\) | \(2.7\times10^5\) |
| \(2.25\times10^{-4}\) | \(8.4\times10^{-7}\) | \(1.2\times10^6\) |
| \(10^{-4}\) | \(3.7\times10^{-7}\) | \(2.7\times10^6\) |
| \(10^{-5}\) | \(3.7\times10^{-8}\) | \(2.7\times10^7\) |

### Reference literature bound

Stadler & Bœhm obtain, at 95% confidence,

\[
\sigma_{\rm DM-\gamma}
<
2.25\times10^{-6}\sigma_T
\left(\frac{m_{\rm DM}}{\rm GeV}\right)
\quad\Longrightarrow\quad
u < 2.25\times10^{-4},
\]

which maps to

\[
\frac{\sigma}{M}\lesssim 8.4\times10^{-7}\ {\rm cm^2\,g^{-1}},
\qquad
\boxed{\Sigma/Q \gtrsim 1.2\times10^6\ {\rm g\,cm^{-2}}.}
\]

Wilkinson et al. show the characteristic suppression of small-scale CMB and matter power around \(u\sim10^{-4}\). Before using this bound as a validation target, record exactly which data set (Planck release, likelihood combination) and which prior on \(u\) it was derived with; see Section 9.

---

## 3. Interaction rates and decoupling

There are two distinct rates, and they answer different questions:

\[
\Gamma_{\gamma\rightarrow{\rm cl}}
\sim n_{\rm cl}\sigma c
= \frac{\rho_{\rm cl}}{M}\sigma c
\propto f_{\rm cl}\,u
\]

is the **photon opacity** due to clumps, while

\[
\Gamma_{{\rm cl}\rightarrow\gamma}
\sim
\frac{4\rho_\gamma}{3\rho_{\rm cl}}
\Gamma_{\gamma\rightarrow{\rm cl}}
\propto u
\]

is the **clump momentum-drag rate**.

The question "when do the clumps cease to be dynamically dragged by radiation?" is answered by

\[
\Gamma_{{\rm cl}\rightarrow\gamma}(z_{\rm dec}) \sim H(z_{\rm dec}).
\]

Extract both rates from CLASS's own thermodynamics output (the idm–photon interaction rate and its \(4\rho_\gamma/3\rho_{\rm idm}\)-scaled counterpart) and validate them against the expressions above, rather than inventing an external diagnostic. Plot both \(\Gamma/H\) against redshift for representative couplings, and report \(z_{\rm dec}\) for every model:

> the clumps must dynamically decouple from the radiation field before approximately \(z_{\rm dec}\).

The different scalings (\(f_{\rm cl}u\) vs \(u\)) motivate the degeneracy test in Experiment C.

---

## 4. Reference cosmology and code version

### CLASS version

Pin a modern CLASS release (v3.2 or later) that exposes the unified interacting-DM parameters `omega_idm`, `f_idm`, `m_idm`, `u_idm_g` / `cross_idm_g`, and `n_index_idm_g`. Record the exact git commit.

Always set explicitly:

- `n_index_idm_g = 0`: constant (temperature-independent) cross-section, which is correct for geometric cross-sections. This is also the CLASS default, but state it anyway.
- `m_idm`: see validation M2b. CLASS defaults to \(10^9\) eV = 1 GeV.

### Baseline

Freeze a single Planck-like baseline, e.g.

\[
\omega_b \simeq 0.0224,
\qquad
\omega_{\rm dm} \simeq 0.120,
\]

with fixed \(H_0\) or \(\theta_s\), \(n_s\), \(A_s\), \(\tau\), neutrino assumptions, and helium/recombination settings.

Generate \(C_\ell^{TT}, C_\ell^{TE}, C_\ell^{EE}, C_L^{\phi\phi}\) and \(P(k,z=0)\), and save the exact parameter dictionary as `baseline_lcdm`.

---

## 5. Validation (before any new science)

### M2a — Zero coupling reproduces CDM

\[
\omega_{\rm cdm}=0,
\qquad
\omega_{\rm idm}=\omega_{\rm dm},
\qquad
u_{\rm idm\_g}=0
\quad(f_{\rm cl}=1,\ u=0).
\]

Compare against `baseline_lcdm` for TT, TE, EE, lensing and \(P(k)\):

\[
\Delta_\ell
=
\frac{
C_\ell({\rm IDM},u=0)
-
C_\ell(\Lambda{\rm CDM})
}{
C_\ell(\Lambda{\rm CDM})
},
\qquad
|\Delta_\ell|\lesssim 10^{-4}
\]

over the scientifically relevant range. For TE near zero-crossings, use absolute differences scaled by \(\sqrt{C^{TT}C^{EE}}\).

### M2b — Cold-mass limit in `m_idm`

At fixed non-zero \(u\) (e.g. \(10^{-4}\)), scan

\[
m_{\rm idm} = 1\,{\rm GeV},\ 10^3\,{\rm GeV},\ 10^6\,{\rm GeV},\ \ldots
\]

until all outputs become numerically invariant. We do not need \(m_{\rm idm}\) to equal a planet's literal mass. We need to show that the thermal/sound-speed contribution of the CLASS particle (\(c_s^2\sim T_{\rm idm}/m_{\rm idm}\)) has become irrelevant. Fix `m_idm` at a value safely inside the invariant regime for all later runs, and watch for numerical problems at extreme masses.

### M2c — Reproduce the published bound

Reproduce the Stadler & Bœhm \(u\lesssim2.25\times10^{-4}\) bound using **their data set and their prior**. This validates the whole likelihood pipeline end to end. It depends on the likelihood infrastructure (M7) and is scheduled there, but no novel result is claimed until it passes.

---

## 6. Experiment A: Thomson-like clumps (stock CLASS)

### 6.1 One-dimensional coupling scan

Set \(f_{\rm cl}=1\). Concentrate the grid where the CMB is sensitive:

\[
\log_{10}u \in [-5,-2]\ \text{in steps of } 0.125,
\]

plus a few null-check points at \(u\sim10^{-8}\), and a few larger values to illustrate the strongly coupled regime.

For each model record \(u\), \(\sigma/M\), \(\Sigma/Q\), \(z_{\rm dec}\).

### 6.2 Diagnostic plots

For TT, plot \(D_\ell^{TT}=\ell(\ell+1)C_\ell^{TT}/2\pi\) for several \(u\), with the fractional residual

\[
100\times\frac{C_\ell(u)-C_\ell(0)}{C_\ell(0)}
\]

below it. Repeat for EE, TE (with care at zero-crossings), lensing and \(P(k)\). Mark the first several acoustic peaks.

Questions to answer:

1. Does increasing \(u\) primarily alter peak heights?
2. Are there measurable phase shifts?
3. At which multipoles does the effect first become significant?
4. Is the damping tail especially sensitive?
5. At what coupling does matter-power suppression become important?
6. Does the strongly coupled limit approach the expected behaviour of a photon-loaded matter component?

### 6.3 Small-scale suppression diagnostic

For **every** CLASS run (all experiments), derive quantitative small-scale suppression measures from the transfer function \(T^2(k)=P(k)/P_{\Lambda{\rm CDM}}(k)\):

- \(k_{10\%}\): where \(T^2\) first drops by 10%;
- \(k_{1/2}\): where \(T^2=1/2\);
- the corresponding half-mode mass \(M_{1/2}\).

For CMB-allowed \(u\lesssim10^{-4}\), the suppression is expected to move to roughly \(k\gtrsim0.2\,h\,{\rm Mpc}^{-1}\) (Wilkinson et al.). This costs essentially nothing. **Converting it into an observational constraint is a separate stage.** Do not import small-scale bounds into this plan until a paper is identified that analyses **DM–photon** scattering specifically with those data. DM–proton (e.g. Maamari et al.) and DM–dark-radiation Lyman-\(\alpha\) bounds are not interchangeable with this model. One candidate to check is Bœhm et al. (2014), on Milky Way satellites and DM–radiation interactions; it has not yet been verified.

### 6.4 First-pass detectability statistic

Before the real likelihood, use a cosmic-variance diagnostic:

\[
\Delta\chi^2_{\rm CV}
\simeq
\sum_{\ell}
\frac{
\left(
C_\ell^{\rm model}-C_\ell^{u=0}
\right)^2
}{
{\rm Var}(C_\ell)
},
\qquad
{\rm Var}(C_\ell)
\simeq
\frac{2}{2\ell+1}C_\ell^2
\ \text{(full-sky TT)}.
\]

Extend it to TT/TE/EE covariance, noise, sky fraction and realistic \(\ell\) ranges. This is a **diagnostic**, not an exclusion.

---

## 7. Experiment B: opaque-clump collision operator

### Motivation

Stock `idm_g` uses the **Thomson** angular collision operator: scattered radiation partly regenerates the photon quadrupole (the \(\Pi\) term), and the Thomson polarization kernel sources E-modes. An optically thick macroscopic object exchanges momentum with radiation through a different microphysical law, and **the Thomson kernels used by stock CLASS need not apply to it**. TE and EE therefore require dedicated treatment.

### Structure of the problem

For any kernel, the dipole (\(\ell=1\)) momentum exchange is fixed by \(\sigma_{\rm mt}\) and the clump velocity. That is what \(u\) encodes. For ideal instantaneous thermalisation the monopole is energy-conserving. The kernels differ in:

- the **\(\ell\ge2\) redistribution**, i.e. how much anisotropy is regenerated rather than erased (Thomson: partial quadrupole regeneration; ideal isotropic re-emission in the clump frame: none);
- the **polarization source**. Thomson has a known one. Ideal isotropic thermal re-emission plausibly gives none, but reflection, limb effects, anisotropic incident radiation and non-instantaneous thermalisation can each change this. State the polarization law explicitly for each kernel rather than assuming it.

### Candidate kernels

1. **B1 — perfect absorber, isotropic thermal re-emission in the clump rest frame** (\(Q=1\)).
2. **B2 — diffuse (Lambertian) surface reflection** (\(Q=13/9\)).
3. **B3 — specular reflection from a sphere** (\(Q=1\)).

For each, derive the Legendre moments of the angular redistribution function and the polarization source, and express them as coefficients multiplying the \(\ell=2\) and polarization source terms in the photon hierarchy. Thomson corresponds to specific values of these coefficients, so stock CLASS is recovered as a special case.

### Implementation

- Add kernel coefficients (quadrupole-regeneration factor, polarization-source factor) to the idm–photon collision terms in CLASS's perturbation module. Default them to the Thomson values.
- The idm–photon **tight-coupling approximation** also embeds the Thomson kernel (e.g. through the shear/quadrupole closure). Either rederive it for the new kernels or disable it, and verify the difference.
- Keep \(u\) defined as the **momentum-transfer** cross-section per unit mass, so that A and B are compared at equal drag. Map to \(\Sigma\) using the kernel-appropriate \(Q\).
- Validation: with Thomson coefficients, the modified code must reproduce stock CLASS bit-for-bit (or to numerical precision).

### Outputs

- \(\Delta C_\ell^{TT,TE,EE}\) between kernels at equal \(u\).
- The shift in the \(f_{\rm cl}=1\) upper limit on \(u\), and hence on \(\Sigma\), for each kernel.
- A verdict: is the Thomson-kernel bound an adequate proxy for opaque clumps, or not?

---

## 8. Experiment C: partial clumped fraction

Vary both parameters, e.g.

\[
f_{\rm cl} = 0.01,\ 0.03,\ 0.1,\ 0.3,\ 0.5,\ 0.8,\ 1,
\qquad
\log_{10}u \in [-5,-1]
\]

(extending to larger \(u\) at small \(f_{\rm cl}\)). Run with the Thomson kernel and with any Experiment B kernel that differs materially.

For each point compute TT, TE, EE, lensing, \(P(k)\), the small-scale diagnostics and \(z_{\rm dec}\).

### Degeneracy test

Photon opacity scales as \(f_{\rm cl}u\), while the clump drag scales as \(u\) (Section 3). Plot the constraints in three coordinate systems:

\[
(f_{\rm cl},u),
\qquad
(f_{\rm cl},f_{\rm cl}u),
\qquad
(f_{\rm cl},\Sigma/Q).
\]

If contours collapse approximately in \(f_{\rm cl}u\), photon opacity dominates the observable effect. Where they do not collapse, clump dynamics matters. Test this rather than assuming it; either way the result says something about the physics.

---

## 9. Statistical inference

### Fixed cosmology

With the six standard parameters fixed, evaluate the CMB likelihood with **CLASS + Cobaya** and report

\[
\Delta\chi^2(u,f_{\rm cl})
=
-2\left[\ln L(u,f_{\rm cl})-\ln L(u=0,f_{\rm cl})\right],
\]

i.e. relative to the nested no-interaction model, not an unrelated global best fit. This answers: how strongly does the measured CMB object to the interaction if nothing else changes?

### Refitted cosmology: profile likelihood as primary result

The fixed-cosmology test overstates detectability, because standard parameters can absorb part of the effect. The science question is a **one-sided upper limit** with \(u=0\) as the nested \(\Lambda\)CDM boundary. Therefore:

\[
\boxed{\text{primary result: profile likelihood in }u\text{ (and }f_{\rm cl}\text{), marginalising nothing, optimising over the standard parameters}}
\]

Profiled parameters: \(\{\omega_b,\omega_{\rm dm},\theta_s,n_s,A_s,\tau\}\) plus nuisance parameters.

### Bayesian cross-check

Run posterior sampling as a secondary check, with **linear or Jeffreys-type priors on \(u\)**. A log-flat prior on \(u\) gives an upper bound that depends strongly on the arbitrary lower cutoff (Diacoumis & Wong, specifically for DM–photon scattering). Use a log-flat prior **only** where needed to reproduce a published analysis for M2c.

---

## 10. Main result

Confidence regions in

\[
\boxed{f_{\rm cl}\ \text{vs}\ \sigma_{\rm mt}/M}
\qquad\text{and}\qquad
\boxed{f_{\rm cl}\ \text{vs}\ \Sigma/Q},
\]

for the Thomson kernel and for each physically motivated opaque-clump kernel.

The headline numbers are the \(f_{\rm cl}=1\) limits:

> Could all of the conventional dark-matter component be compact baryonic clumps, and if so, what minimum surface density \(\Sigma_{\min}\) do they need, **for the correct collision operator**?

The expected Thomson-kernel answer is \(\Sigma/Q\gtrsim1.2\times10^6\ {\rm g\,cm^{-2}}\). The novel content is how, and whether, that changes for opaque clumps, and how it scales with \(f_{\rm cl}\).

### Control

Include explicitly \(f_{\rm cl}=1,\ u=0\). After optimisation it must be observationally equivalent to CDM. This demonstrates that the primary CMB is sensitive to density, pressure, sound speed, anisotropic stress, interactions and perturbation evolution, **not** to whether collisionless pressureless matter is baryonic. The baryonic nature of a compact object becomes visible only through additional microphysics, and Experiment B is exactly a test of such microphysics.

---

## 11. Experiment D: formation history (only if warranted)

Replace constant coupling with a redshift-dependent one:

\[
u(z)=
\begin{cases}
u_{\rm early}, & z>z_{\rm form},\\
u_{\rm cl}, & z<z_{\rm form},
\end{cases}
\]

or more physically \(u(z)=268\,Q/\Sigma(z)\) for a collapsing proto-cloud, possibly with \(f_{\rm cl}(z)\). The model becomes \((f_{\rm cl},\Sigma_{\rm final},z_{\rm form})\). This requires further CLASS modification and should only be attempted after A–C are complete.

---

## 12. Explicitly excluded

The project does **not** establish that baryonic dark matter is viable. It does not address:

- Big Bang nucleosynthesis;
- the origin of baryon inhomogeneities;
- cloud formation, cooling, ionisation or evaporation;
- internal radiative transfer beyond the surface kernel of Experiment B;
- clump–clump collisions;
- Poisson/shot-noise effects of rare massive objects;
- microlensing or other gravitational-lensing limits;
- Galactic dynamics;
- spectral distortions and energy injection;
- late-time gas-cloud observational constraints;
- observational small-scale-structure constraints (Section 6.3 produces diagnostics only).

The conclusion should remain narrow:

> **Assuming a pressureless compact component already exists before the relevant CMB epoch and exchanges momentum with photons with cross-section \(\sigma_{\rm mt}/M\) through a specified collision operator, CMB anisotropy data permit or exclude a region of \((f_{\rm cl},\sigma_{\rm mt}/M)\), equivalently \((f_{\rm cl},\Sigma/Q)\).**

---

## 13. Outputs and documentation (requirement)

The work must be **consumable by someone who was not present**: a reader wanting the answers, a reviewer checking how they were obtained, or a future implementer continuing the project. These three documents are required deliverables, not optional write-ups. A milestone is not complete until they are updated.

**Format.** Documents are Markdown in `doc/` (extended with LaTeX maths, tables and embedded figures from `figures/`). Markdown is the canonical record because it diffs cleanly in git and renders on GitHub. Notebooks in `notebooks/` are welcome for presenting and exploring results. Anything a notebook establishes must also be recorded in the Markdown documents, and every notebook must re-run from a clean kernel using only committed code and saved results.

### 13.1 Results summary — `doc/RESULTS.md`

A short, current statement of what is known. It should be readable in ten minutes.

- Headline numbers with uncertainties and confidence level. The key ones are the \(f_{\rm cl}=1\) limits on \(u\), \(\sigma_{\rm mt}/M\) and \(\Sigma/Q\), per kernel.
- Principal figures: spectra residuals, \((f_{\rm cl},\Sigma/Q)\) contours, and the \(\Gamma/H\) and \(z_{\rm dec}\) plots.
- For each result: the assumptions it rests on, the kernel, the data set and the statistical method (profile or posterior, and which prior).
- Status of each experiment and validation gate (A–D, M2a–c): passed, pending or failed.
- Rewritten in place as results change. The history lives in the experimental record.

### 13.2 Experimental record — `doc/experiment_log/`

An append-only lab notebook, with one dated entry (`YYYY-MM-DD_<short-name>.md`) per experiment or significant run batch. Each entry records:

- **Goal:** the question or milestone being addressed.
- **Configuration:** CLASS version and git commit, any patch applied, config file(s), the parameter dictionary or a pointer to the saved one, and the Python environment.
- **What was run:** grid, number of models, compute time, and the command or script used.
- **Outputs:** paths under `results/` and `figures/`.
- **Outcome:** what was found, including null and failed results and runs that were discarded, with the reason.
- **Follow-ups:** what this changes or motivates next.

Entries are never rewritten after the fact. Corrections go in a new entry that links back to the old one. Every number in `RESULTS.md` must trace to a log entry.

### 13.3 Implementation notes — `doc/IMPLEMENTATION_NOTES.md`

A living guide for future implementers, organised by topic rather than by date, covering what was learned the hard way. For example:

- CLASS parameter behaviour and pitfalls: `m_idm` invariance thresholds, `n_index_idm_g`, `omega_idm` vs `f_idm` bookkeeping, precision settings required for the \(10^{-4}\) agreement.
- Numerical issues: stiffness at large \(u\), tight-coupling behaviour, and extreme-mass problems.
- Experiment B: the kernel derivations (or a pointer to them), which CLASS source locations were patched and why, and how the Thomson-limit regression test works.
- Likelihood and sampling: Cobaya configuration, the reproduction of the published bound, prior choices, and profile-likelihood technique.
- Dead ends and approaches that did not work, with the reason.
- Open questions and recommended next steps.

---

## 14. Suggested repository structure

```text
class_baryons/
├── README.md
├── environment.yml
├── doc/
│   ├── cmb_compact_baryon_clumps_class_experiment_plan.md
│   ├── RESULTS.md                 # current results summary
│   ├── IMPLEMENTATION_NOTES.md    # lessons for future implementers
│   └── experiment_log/            # append-only dated run records
│       └── YYYY-MM-DD_<name>.md
├── configs/
│   ├── baseline.yaml
│   └── scan.yaml
├── src/
│   ├── cosmology.py
│   ├── run_class.py
│   ├── cloud_mapping.py
│   ├── rates.py
│   ├── small_scale.py
│   ├── scan.py
│   └── diagnostics.py
├── class_patches/          # Experiment B kernel modifications
├── notebooks/
│   ├── 01_validation.ipynb
│   ├── 02_expA_coupling_scan.ipynb
│   ├── 03_expB_kernels.ipynb
│   ├── 04_expC_fraction_scan.ipynb
│   └── 05_likelihood.ipynb
├── results/
│   ├── spectra/
│   ├── scans/
│   ├── profiles/
│   └── chains/
└── figures/
```

CLASS execution logic lives in normal Python scripts, not only in notebooks.

---

## 15. Milestones

**Definition of done (every milestone):** add an experiment-log entry, update `RESULTS.md` if any result or gate status changed, and add anything learned to `IMPLEMENTATION_NOTES.md` (Section 13).

### M0 — Environment
- Install a pinned CLASS (≥ v3.2) / `classy`; record the git commit.
- Verify the standard example; freeze the Python environment.
- Create `RESULTS.md`, `IMPLEMENTATION_NOTES.md` and `experiment_log/` with their templates.

### M1 — Baseline
- Run the reference \(\Lambda\)CDM model; save spectra, \(P(k)\) and the parameter dictionary.

### M2 — Validation
- **M2a:** \(u=0\) IDM reproduces CDM.
- **M2b:** `m_idm` scan reaches the cold-mass limit; fix `m_idm`. Set `n_index_idm_g = 0` explicitly.
- (**M2c** is scheduled at M7.)

### M3 — Experiment A scan
- \(f_{\rm cl}=1\), dense grid in \(10^{-5}\lesssim u\lesssim10^{-2}\) plus null checks.
- Diagnostic plots, CLASS-native \(\Gamma/H\) rates, \(z_{\rm dec}\), small-scale diagnostics.

### M4 — Physical mapping
- Convert to \(\sigma/M\), \(\Sigma/Q\), \(z_{\rm dec}\), \(k_{1/2}\), \(M_{1/2}\).

### M5 — Approximate detectability
- Cosmic-variance / covariance \(\Delta\chi^2\) relative to \(u=0\).

### M6 — Experiment B kernels
- Derive the B1–B3 kernel coefficients (angular redistribution + polarization).
- Patch CLASS, including the tight-coupling treatment; verify that Thomson coefficients reproduce stock CLASS.
- Compare TT/TE/EE between kernels at equal momentum-transfer \(u\).

### M7 — Likelihood pipeline and published-bound reproduction
- Integrate CLASS with Cobaya.
- **M2c:** reproduce the Stadler & Bœhm \(u\lesssim2.25\times10^{-4}\) bound with matching data and prior. Gate: no novel claim until this passes.

### M8 — Profile likelihoods
- Fixed-cosmology \(\Delta\chi^2\) relative to \(u=0\).
- Profile likelihood in \(u\) at \(f_{\rm cl}=1\), per kernel.
- Bayesian cross-check with linear/Jeffreys priors.

### M9 — Experiment C
- \((f_{\rm cl},u)\) grid and profiles; degeneracy test in \(f_{\rm cl}u\).

### M10 — Principal result
- Confidence regions in \((f_{\rm cl},\sigma_{\rm mt}/M)\) and \((f_{\rm cl},\Sigma/Q)\) per kernel; \(f_{\rm cl}=1\) limits on \(\Sigma\).

### M11 — Decision point
Assess whether to proceed to Experiment D (redshift-dependent coupling, formation), recombination modifications, BBN, observational small-scale constraints for DM–photon coupling, or other astrophysical constraints.

---

## Summary

The CMB constrains how matter interacts, not what it is made of. Stock CLASS encodes that interaction with an elementary-particle (Thomson-like) collision operator. This project:

1. reproduces the established DM–photon bound and translates it into a minimum clump surface density (\(\Sigma/Q\gtrsim10^6\ {\rm g\,cm^{-2}}\) expected);
2. **tests whether an optically thick macroscopic object has the same CMB collision operator, and how the \(\Sigma\) bound changes if not**;
3. maps the allowed \((f_{\rm cl},\Sigma/Q)\) space and tests the \(f_{\rm cl}u\) degeneracy;
4. reports decoupling redshifts and small-scale suppression scales throughout;
5. bases its main limits on profile likelihoods, avoiding prior-driven upper bounds.

Only after this should the project move on to cloud formation, early-universe baryon inhomogeneity, BBN, recombination microphysics and other astrophysical constraints.
