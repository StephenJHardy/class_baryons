# CMB Constraints on Compact Weakly Photon-Coupled Matter

## Experimental Plan Using CLASS

> **Revision 2 (2026-09-25, after M3/M4).**
> - **The finding:** the M3 coupling scan showed that at CMB-allowed couplings the CMB responds to radiation *dragging the clumps*, which depends only on the momentum-transfer cross-section. Photons scattering off clumps is negligible: clump opacity is ≲10⁻⁴ of Thomson before recombination, and τ_cl ≲ 10⁻³ after.
> - **Experiment B:** the collision-operator question it was built around is therefore expected to have a near-null answer. B is cut down to a minimal confirmation check (§7).
> - **New experiments:** the novel direction moves to physics that *does* distinguish opaque baryonic clumps: small-scale structure (Experiment E, §11a) and energy exchange / spectral distortions (Experiment F, §11b). The formation history (Experiment D) is unchanged.
> - **Details:** `doc/RESULTS.md` and `doc/experiment_log/2026-09-25_m3_m4_coupling_scan.md`.

### Research question

This project asks:

> **Given a fraction \(f_{\rm cl}\) of the would-be dark matter in compact, pressureless objects, what photon momentum-transfer cross-section per unit mass can that population have while remaining compatible with the observed CMB?**

and, more sharply:

> **Does an optically thick macroscopic baryonic object have the same CMB collision operator as the particle-scattering model conventionally used to derive dark-matter–photon bounds? If not, how does the CMB constraint on its surface density \(\Sigma\) change when the correct macro-object radiative-transfer physics is used?**

*(Revision 2: M3 indicates the answer is "no, but it barely matters, because the constraint is set by momentum transfer alone". Experiment B confirms this minimally.)* The questions that now carry the novelty are:

> **Which observables other than the CMB anisotropies constrain opaque baryonic clumps more tightly, specifically small-scale structure (the drag-induced cutoff in P(k)) and CMB spectral distortions (absorption and re-emission exchange energy, unlike elastic Thomson scattering)? What minimum \(\Sigma\) do they require?**

The project deliberately does **not** model cloud formation, Big Bang nucleosynthesis, or detailed internal cloud hydrodynamics. It asks whether compact clumps already present before the CMB acoustic epoch can behave sufficiently like cold dark matter to avoid altering the observed anisotropies.

The public **CLASS** Boltzmann code already contains an interacting dark matter species (`idm`) with photon coupling (`u_idm_g` / `cross_idm_g`). That implementation follows the Stadler & Bœhm elastic DM–photon scattering formalism, i.e. a **Thomson-like angular collision operator**. It is therefore treated here as a **reference approximation** (a momentum-drag proxy), not as the physical model of an opaque cloud.

### Project structure

| Experiment | Question | CLASS modification? |
|---|---|---|
| **A** — Thomson-like clumps | Reproduce the established interacting-DM bound at \(f_{\rm cl}=1\) and translate it into \(\Sigma\). | No |
| **B** — Opaque-clump collision operator (minimal check) | Confirm that swapping the Thomson kernel for an opaque-clump kernel leaves TT/TE/EE unchanged at equal momentum transfer, and isolate the possible low-ℓ EE term. | Yes (small patch: two coefficients) |
| **C** — Partial clumped fraction | Map the allowed \((f_{\rm cl},\Sigma/Q)\) space. | No (unless B surprises) |
| **D** — Formation history | \(\Sigma(z)\), \(f_{\rm cl}(z)\), \(z_{\rm form}\). Only if warranted. | Yes |
| **E** — Small-scale structure | Turn the drag-induced P(k) cutoff into a bound on u (and Σ) and compare it with the CMB bound. | No (higher k precision only) |
| **F** — Energy exchange and spectral distortions | Do absorbing/re-emitting clumps distort the CMB spectrum during the μ and y eras, and what Σ does FIRAS (and future missions) allow? | Estimate first; possibly CLASS's distortions module |

Reproducing a Planck \(u\lesssim10^{-4}\) bound is established science and serves as validation. **Revision 2: the scientifically novel content is now expected from Experiments E and F.** B is a short confirmation that the Thomson-kernel CMB bound carries over to opaque clumps.

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

plus a few null-check points at \(u\sim10^{-8}\), and a few larger values to illustrate the strongly coupled regime. *(M3: stock CLASS v3.4.0 runs only up to \(u\approx0.13\) at \(f_{\rm cl}=1\).)*

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

## 7. Experiment B: opaque-clump collision operator (minimal confirmation)

### Motivation

Stock `idm_g` uses the **Thomson** angular collision operator. Scattered radiation partly regenerates the photon quadrupole (the \(\Pi\) term), and the Thomson polarization kernel sources E-modes. An optically thick macroscopic object exchanges momentum with radiation through a different microphysical law, so **the Thomson kernels used by stock CLASS need not apply to it**.

### Why this is expected to matter little (M3)

- **The drag is kernel-independent.** To first order in velocity and anisotropy, the force on any scatterer is the photon momentum flux times \(\sigma_{\rm mt}=\int(1-\cos\theta)\,d\sigma\). So the clump drag, which is what produces the CMB signal, depends only on \(u\).
- **The kernel only enters the photon-side \(\ell\ge2\) and polarization terms.** Those are weighted by clump opacity, which is ≲10⁻⁴ of Thomson before recombination and adds up to only τ_cl ≲ 10⁻³ afterwards.
- **One exception to test.** EE departs from ΛCDM first at ℓ ≈ 7–17. This could be clump-generated polarization after recombination, and that would depend on the kernel.

### Scope

Two kernels only:
- **Thomson (stock).** Quadrupole-regeneration and polarization-source coefficients = 1.
- **B1, a perfect absorber with isotropic thermal re-emission in the clump rest frame (\(Q=1\)).** Clump opacity damps every photon moment with \(\ell\ge2\) and the polarization hierarchy, with no \(\Pi\) regeneration and no polarization source. Both coefficients are 0.

B1 is the extreme case: it removes every photon-side term the Thomson kernel adds. If B1 matches Thomson, any intermediate kernel will too.
- B2 (Lambertian, \(Q=13/9\)) and B3 (specular, \(Q=1\)) enter only through their \(Q\) values in the \(\Sigma\) mapping.
- The kernels also differ in their polarization law: ideal re-emission plausibly gives none, but reflection, limb effects and non-instantaneous thermalisation can each change that. Derive the full B2/B3 coefficients only if B1 differs materially from Thomson.

### Implementation

- In CLASS's perturbation module, add the two coefficients (quadrupole regeneration and polarization source) to the idm–photon collision terms. Default them to the Thomson values.
- The idm–photon **tight-coupling approximation** also embeds the Thomson kernel, through the shear/quadrupole closure. Disable it for B1, or rederive it, and check the difference with Thomson coefficients.
- Keep \(u\) defined as the momentum-transfer cross-section per unit mass, so that Thomson and B1 are compared at equal drag.
- Keep the patch as a versioned diff in `class_patches/`, applied by the setup script to a separate build, so stock CLASS stays available.
- **Validation:** with Thomson coefficients, the patched code must reproduce stock CLASS to numerical precision, using the project precision settings and the `src/glitch_scan.py` smoothness check.

### Outputs

- \(\Delta C_\ell^{TT,TE,EE,\phi\phi}\) between B1 and Thomson at \(u = 10^{-5}, 10^{-4}, 2.25\times10^{-4}, 10^{-3}\), with ideal-CV \(\Delta\chi^2\).
- An isolated low-ℓ EE contribution from clump scattering: the B1 − Thomson difference measures it directly.
- A verdict. **Expected:** differences ≲10⁻⁴ outside low-ℓ EE, so the Thomson-kernel bound applies to opaque clumps with \(\Sigma_{\min}=Q\cdot268/u_{\max}\) and the kernel-appropriate \(Q\). If instead differences exceed the detectability threshold, restore the full B1–B3 programme and carry the kernels through C and the likelihood.

---

## 8. Experiment C: partial clumped fraction

Vary both parameters, e.g.

\[
f_{\rm cl} = 0.01,\ 0.03,\ 0.1,\ 0.3,\ 0.5,\ 0.8,\ 1,
\qquad
\log_{10}u \in [-5,-1]
\]

(extending to larger \(u\) at small \(f_{\rm cl}\)). Run with the Thomson kernel, plus any Experiment B kernel that turns out to differ materially.

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

*(Revision 2: M3 shows the drag dominates at f_cl = 1. Collapse in \(f_{\rm cl}u\) is therefore **not** expected. More likely the constraint behaves like \(f_{\rm cl}\times g(u)\): the drag rate depends only on \(u\), while the resulting loss of gravitational clustering scales with the clumped fraction.)*

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

for the Thomson kernel (plus B1 only if Experiment B finds a material difference). Compare them with the bounds from Experiments E and F in the same coordinates.

The headline numbers are the \(f_{\rm cl}=1\) limits:

> Could all of the conventional dark-matter component be compact baryonic clumps, and if so, what minimum surface density \(\Sigma_{\min}\) do they need, **for the correct collision operator**?

The expected Thomson-kernel answer is \(\Sigma/Q\gtrsim1.2\times10^6\ {\rm g\,cm^{-2}}\). Revision 2: the novel content is
- confirming that this carries over to opaque clumps (B);
- how it scales with \(f_{\rm cl}\) (C);
- above all, whether small-scale structure (E) or spectral distortions (F) require a much larger \(\Sigma_{\min}\) than the CMB anisotropies do. M3's linear half-mode mass of ~10¹⁴ M☉/h at the CMB bound suggests E will.

### Control

Include explicitly \(f_{\rm cl}=1,\ u=0\). After optimisation it must be observationally equivalent to CDM. This demonstrates that the primary CMB is sensitive to density, pressure, sound speed, anisotropic stress, interactions and perturbation evolution, **not** to whether collisionless pressureless matter is baryonic. The baryonic nature of a compact object becomes visible only through additional microphysics. Experiments B (angular kernel) and F (energy exchange) test such microphysics directly.

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

*(Revision 2: D gains importance if E or F turn out to be binding. A cloud that forms or compacts late avoids early drag and early energy exchange, so \(z_{\rm form}\) could relax both.)*

---

## 11a. Experiment E: small-scale structure

### Motivation

- The drag suppresses clump perturbations on scales that enter the horizon before \(z_{\rm dec}\).
- At the published CMB bound, M3 finds a linear half-mode mass of \(M_{\rm hm}\approx9\times10^{13}\,M_\odot/h\) and a 10% drop in σ₈. Structure below galaxy-cluster scale is strongly suppressed.
- Galaxy-scale structure exists, so small-scale observables (the Lyman-α forest, Milky Way satellite counts, strong-lensing substructure) are likely to bound \(u\), and hence \(\Sigma\), orders of magnitude more tightly than the CMB.

### Steps

1. **E1 — Literature.** Find analyses that constrain **DM–photon** scattering with small-scale data (a Bœhm et al. 2014 Milky Way satellite analysis is a candidate to verify) and convert their bounds to \(u\). Do not substitute DM–baryon or DM–dark-radiation bounds.
2. **E2 — Transfer-function mapping.** Extend P(k) to \(k\sim10^2\ h/{\rm Mpc}\): raise `P_k_max_h/Mpc` and re-check the high-k precision (M2b found ~10⁻³ noise in strongly damped tails). Compute \(T^2(k)\), \(k_{\rm hm}\) and \(M_{\rm hm}\) down to \(u\sim10^{-9}\).
3. **E3 — Approximate bound.** Match \(k_{\rm hm}\) to thermal-relic warm dark matter and apply published WDM half-mode bounds (Lyman-α, satellites, lensing) as an **approximate** constraint on \(u\). State the approximation clearly: DM–photon transfer functions have damped oscillations and are not identical in shape to WDM.
4. **E4 — Clump-specific checks.**
   - The fluid approximation must hold at the half-mode scale (\(M_{\rm hm}/M_{\rm cl}\gg1\)).
   - After decoupling the clumps behave as CDM, but they are baryons: note, without yet modelling them, any late-time effects (collisions, ram pressure) that would invalidate the CDM-like treatment on galactic scales.

### Outputs

- \(u_{\max}\) and \(\Sigma_{\min}/Q\) from small scales at \(f_{\rm cl}=1\), with the approximation stated.
- Their dependence on \(f_{\rm cl}\), run on the Experiment C grid; mixed CDM + suppressed-component transfer functions weaken the bound at small \(f_{\rm cl}\).
- A comparison figure: CMB versus small-scale \(\Sigma_{\min}\) against \(f_{\rm cl}\).

---

## 11b. Experiment F: energy exchange and spectral distortions

### Motivation

- **Thomson scattering is elastic** (energy-conserving to \(O(v^2)\)). An opaque clump **absorbs** photons and re-emits them at its own temperature, with its own emissivity spectrum.
- **At early times photons interact with clumps many times per Hubble time.** M3 finds \(\Gamma_{\gamma\to{\rm cl}}/H=1\) at \(z\approx8\times10^5\) for \(u=10^{-4}\), inside the μ era (\(5\times10^4\lesssim z\lesssim2\times10^6\)).
- **So any departure of the re-emitted spectrum from a blackbody at \(T_\gamma\) produces a spectral distortion.** Such departures include a clump temperature different from \(T_\gamma\), non-grey emissivity, or energy sources inside the clumps. FIRAS limits (\(|\mu|\lesssim9\times10^{-5}\), \(|y|\lesssim1.5\times10^{-5}\)) are strong, and future missions would improve them by orders of magnitude.
- This channel has no analogue for elementary-particle dark matter.

### Steps

1. **F1 — Thermal equilibrium.** Estimate how tightly a clump's (surface) temperature is locked to \(T_\gamma\). Compare the radiative exchange time (heat capacity over \(\sigma_{\rm SB}T^4\) × area) with the Hubble time and the adiabatic-cooling mismatch between gas and radiation. Identify the regimes where \(T_{\rm cl}\neq T_\gamma\).
2. **F2 — Energy-exchange budget.** Derive the fractional energy exchange \(\Delta\rho_\gamma/\rho_\gamma\) per Hubble time as a function of \(u\), \(T_{\rm cl}-T_\gamma\) and emissivity. Include the kinematic (clump bulk-velocity) contribution, which is \(O(\tau_{\rm cl}v^2)\) and y-type.
3. **F3 — Distortion estimate.** Convert F2 into \(\mu\) and \(y\) with the standard distortion visibility functions. Where needed, use CLASS's spectral-distortions module with a custom heating rate. Compare with FIRAS and with a representative future mission.
4. **F4 — Interpretation.** Which cloud properties (temperature offset, emissivity, internal heating) would be excluded, and does a pure blackbody clump at \(T_\gamma\) escape entirely? The latter is expected, which makes F a constraint on clump **microphysics** rather than on \(u\) alone.

### Outputs

- \(\mu(u,\ldots)\) and \(y(u,\ldots)\) estimates, and the region of cloud parameter space that FIRAS excludes.
- A go/no-go on a full treatment. A full treatment would need a Boltzmann-level energy-exchange term and should not be started before F1–F3.

---

## 12. Explicitly excluded

The project does **not** establish that baryonic dark matter is viable. It does not address:

- Big Bang nucleosynthesis;
- the origin of baryon inhomogeneities;
- cloud formation, cooling, ionisation or evaporation;
- internal radiative transfer beyond the surface kernel of Experiment B and the thermal estimates of Experiment F;
- clump–clump collisions;
- Poisson/shot-noise effects of rare massive objects;
- microlensing or other gravitational-lensing limits;
- Galactic dynamics;
- a full Boltzmann treatment of energy exchange (Experiment F is estimates first);
- late-time gas-cloud observational constraints;
- a full small-scale likelihood analysis (Experiment E uses literature bounds and half-mode matching).

The conclusion should remain narrow:

> **Assuming a pressureless compact component already exists before the relevant CMB epoch and exchanges momentum with photons with cross-section \(\sigma_{\rm mt}/M\) through a specified collision operator, CMB anisotropy data permit or exclude a region of \((f_{\rm cl},\sigma_{\rm mt}/M)\), equivalently \((f_{\rm cl},\Sigma/Q)\). Small-scale structure (approximately) and spectral distortions (by estimate) add further, clearly labelled, constraints on the same space and on clump microphysics.**

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

**Status and execution order (revision 2).** M0–M4 are done. Milestone IDs are kept stable because the logs refer to them. The execution order is now:

**M6 (B, minimal) → M12 (E) → M13 (F) → M5 → M7 → M8 → M9 (C) → M10 → M11**

The cheap B confirmation and the E/F estimates come before the heavy likelihood work, because they may change which constraint is the headline.

### M0 — Environment ✅
- Install a pinned CLASS (≥ v3.2) / `classy`; record the git commit.
- Verify the standard example; freeze the Python environment.
- Create `RESULTS.md`, `IMPLEMENTATION_NOTES.md` and `experiment_log/` with their templates.

### M1 — Baseline ✅
- Run the reference \(\Lambda\)CDM model; save spectra, \(P(k)\) and the parameter dictionary.

### M2 — Validation ✅ (M2a, M2b)
- **M2a:** \(u=0\) IDM reproduces CDM.
- **M2b:** `m_idm` scan reaches the cold-mass limit; fix `m_idm`. Set `n_index_idm_g = 0` explicitly.
- (**M2c** is scheduled at M7.)

### M3 — Experiment A scan ✅
- \(f_{\rm cl}=1\), dense grid in \(10^{-5}\lesssim u\lesssim10^{-2}\) plus null checks.
- Diagnostic plots, CLASS-native \(\Gamma/H\) rates, \(z_{\rm dec}\), small-scale diagnostics.

### M4 — Physical mapping ✅
- Convert to \(\sigma/M\), \(\Sigma/Q\), \(z_{\rm dec}\), \(k_{1/2}\), \(M_{1/2}\).

### M5 — Approximate detectability
- Cosmic-variance / covariance \(\Delta\chi^2\) relative to \(u=0\).

### M6 — Experiment B, minimal check (next)
- Patch CLASS with two kernel coefficients, and handle tight coupling for B1. With Thomson values the patch must reproduce stock CLASS.
- Compare B1 with Thomson at equal \(u\); isolate the low-ℓ EE contribution.
- Verdict. Expand to the full B1–B3 programme only if B1 differs materially.

### M7 — Likelihood pipeline and published-bound reproduction
- Integrate CLASS with Cobaya.
- **M2c:** reproduce the Stadler & Bœhm \(u\lesssim2.25\times10^{-4}\) bound with matching data and prior. Gate: no novel claim until this passes.

### M8 — Profile likelihoods
- Fixed-cosmology \(\Delta\chi^2\) relative to \(u=0\).
- Profile likelihood in \(u\) at \(f_{\rm cl}=1\) (Thomson kernel; B1 only if M6 requires it).
- Bayesian cross-check with linear/Jeffreys priors.

### M9 — Experiment C
- \((f_{\rm cl},u)\) grid and profiles; degeneracy test in \(f_{\rm cl}u\).

### M10 — Principal result
- Confidence regions in \((f_{\rm cl},\sigma_{\rm mt}/M)\) and \((f_{\rm cl},\Sigma/Q)\); \(f_{\rm cl}=1\) limits on \(\Sigma\).
- Overlay the small-scale (E) and spectral-distortion (F) constraints.

### M11 — Decision point
Assess whether to proceed to Experiment D (redshift-dependent coupling, formation), a full Boltzmann treatment of energy exchange (F), a full small-scale analysis (E), recombination modifications, BBN, or other astrophysical constraints.

### M12 — Experiment E: small-scale structure
- E1: literature DM–photon small-scale bounds, converted to \(u\).
- E2: P(k) to \(k\sim10^2\ h/{\rm Mpc}\), with a precision check; \(k_{\rm hm}\), \(M_{\rm hm}\) down to \(u\sim10^{-9}\).
- E3: approximate bound via WDM half-mode matching, approximation stated.
- E4: fluid-approximation and baryonic-clump caveats.

### M13 — Experiment F: energy exchange and spectral distortions
- F1: clump thermal locking to \(T_\gamma\).
- F2: energy-exchange budget against \(u\), \(T_{\rm cl}-T_\gamma\) and emissivity.
- F3: μ and y estimates against FIRAS and a future mission.
- F4: interpretation, and a go/no-go on a full treatment.

---

## Summary

The CMB constrains how matter interacts, not what it is made of. Stock CLASS encodes that interaction with an elementary-particle (Thomson-like) collision operator. This project:

1. reproduces the established DM–photon bound and translates it into a minimum clump surface density (\(\Sigma/Q\gtrsim10^6\ {\rm g\,cm^{-2}}\) expected);
2. confirms, with a minimal kernel swap, that the bound depends on the momentum-transfer cross-section alone and so carries over to opaque clumps (M3 shows the CMB signal comes from drag, not photon scattering);
3. **asks which observables distinguish opaque baryonic clumps and constrain them more tightly: small-scale structure (E) and spectral distortions from absorption and re-emission (F)**;
4. maps the allowed \((f_{\rm cl},\Sigma/Q)\) space and tests the \(f_{\rm cl}u\) degeneracy;
5. reports decoupling redshifts and small-scale suppression scales throughout;
6. bases its main CMB limits on profile likelihoods, avoiding prior-driven upper bounds.

Only after this should the project move on to cloud formation, early-universe baryon inhomogeneity, BBN, recombination microphysics and other astrophysical constraints.
