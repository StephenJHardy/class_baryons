"""Apply the Experiment B kernel edits to a CLASS v3.4.0 source tree.

This is the *authoring* tool for class_patches/idm_g_kernel.patch: run it
on a clean copy of CLASS, check the build, then regenerate the patch with
`git diff > ../class_patches/idm_g_kernel.patch`. Users apply the .patch
file, not this script (see scripts/setup_class.sh).

Physics: CLASS's idm-photon interaction uses the Thomson angular kernel.
Two coefficients generalise the clump part of the photon collision term:
  c_T (idm_g_quadrupole_coefficient): fraction of Thomson quadrupole
      regeneration (the P0 term in the temperature l=2 equation and the
      Pi term in the temperature line-of-sight source);
  c_E (idm_g_polarization_coefficient): fraction of the Thomson
      polarization source (P0 terms in the polarization hierarchy and the
      E-mode line-of-sight source).
c_T = c_E = 1 is stock CLASS (Thomson); c_T = c_E = 0 is a perfect absorber
with isotropic re-emission in the clump frame. Damping of all photon
moments with l >= 2 and the l = 1 momentum exchange are unchanged.
"""

import sys
from pathlib import Path

root = Path(sys.argv[1])


def edit(rel, old, new, count=1):
    path = root / rel
    text = path.read_text()
    found = text.count(old)
    if found != count:
        raise SystemExit(f"{rel}: expected {count} occurrence(s), found {found}:\n{old}")
    path.write_text(text.replace(old, new))


# --- thermodynamics.h: new parameters -------------------------------------------------
edit("include/thermodynamics.h",
     "  int n_index_idm_g;     /**< temperature dependence of the interactions between dark matter and photons */",
     "  int n_index_idm_g;     /**< temperature dependence of the interactions between dark matter and photons */\n"
     "  double idm_g_quad_coeff; /**< [class_baryons] clump quadrupole-regeneration coefficient c_T (1 = Thomson, 0 = isotropic absorber) */\n"
     "  double idm_g_pol_coeff;  /**< [class_baryons] clump polarization-source coefficient c_E (1 = Thomson, 0 = isotropic absorber) */")

# --- input.c: defaults and reading ----------------------------------------------------
edit("source/input.c",
     "  pth->n_index_idm_g = 0;\n",
     "  pth->n_index_idm_g = 0;\n  pth->idm_g_quad_coeff = 1.;\n  pth->idm_g_pol_coeff = 1.;\n")
edit("source/input.c",
     "    if (pth->u_idm_g > 0) {\n      class_read_double(\"n_index_idm_g\",pth->n_index_idm_g);\n    }",
     "    if (pth->u_idm_g > 0) {\n      class_read_double(\"n_index_idm_g\",pth->n_index_idm_g);\n"
     "      /* [class_baryons] angular/polarization kernel of the clump-photon collision term */\n"
     "      class_read_double(\"idm_g_quadrupole_coefficient\",pth->idm_g_quad_coeff);\n"
     "      class_read_double(\"idm_g_polarization_coefficient\",pth->idm_g_pol_coeff);\n"
     "    }")

# --- perturbations.c: helper for the tight-coupling closure ---------------------------
edit("source/perturbations.c",
     '#include "perturbations.h"\n',
     '#include "perturbations.h"\n\n'
     "/**\n"
     " * [class_baryons] Tight-coupling closure for a photon collision term made of\n"
     " * Thomson scattering (rate a = kappa') plus clumps (rate b = dmu_idm_g) whose\n"
     " * kernel has quadrupole-regeneration coefficient cT and polarization-source\n"
     " * coefficient cE. Solving the quasi-static l=2 temperature and polarization\n"
     " * equations gives P0 = 2 s_2 shear / (8 - 4.8 r) with r = (a + cE b)/(a + b),\n"
     " * and shear = (4/15)(theta + metric_shear) / D with\n"
     " * D = (a + b) - 0.8 (a + cT b)/(8 - 4.8 r). Returns tau2 = 0.75/D (so that\n"
     " * shear = 16/45 tau2 (...), reducing to 1/(a+b) for Thomson), its conformal\n"
     " * time derivative given da, db, and p_over_shear = 2/(8 - 4.8 r).\n"
     " */\n"
     "static void idm_g_kernel_closure(double a, double b, double da, double db,\n"
     "                                 double cT, double cE,\n"
     "                                 double * tau2, double * dtau2, double * p_over_shear) {\n"
     "  double s = a + b, ds = da + db;\n"
     "  double r = (a + cE*b)/s;\n"
     "  double dr = ((da + cE*db)*s - (a + cE*b)*ds)/s/s;\n"
     "  double m = 8. - 4.8*r, dm = -4.8*dr;\n"
     "  double n = a + cT*b, dn = da + cT*db;\n"
     "  double D = s - 0.8*n/m;\n"
     "  double dD = ds - 0.8*(dn*m - n*dm)/m/m;\n"
     "  *tau2 = 0.75/D;\n"
     "  *dtau2 = -0.75*dD/D/D;\n"
     "  *p_over_shear = 2./m;\n"
     "}\n")

# first-order TCA shear used in the Einstein equations (two gauges)
edit("source/perturbations.c",
     "          shear_g = 16./45./(ppw->pvecthermo[pth->index_th_dkappa] + ppw->pvecthermo[pth->index_th_dmu_idm_g])*(y[ppw->pv->index_pt_theta_g]+k2*ppw->pvecmetric[ppw->index_mt_alpha]);",
     "          double tau2_kc = 1./(ppw->pvecthermo[pth->index_th_dkappa] + ppw->pvecthermo[pth->index_th_dmu_idm_g]), dtau2_kc, pos_kc;\n"
     "          if ((pth->idm_g_quad_coeff != 1.) || (pth->idm_g_pol_coeff != 1.))\n"
     "            idm_g_kernel_closure(ppw->pvecthermo[pth->index_th_dkappa], ppw->pvecthermo[pth->index_th_dmu_idm_g], 0., 0.,\n"
     "                                 pth->idm_g_quad_coeff, pth->idm_g_pol_coeff, &tau2_kc, &dtau2_kc, &pos_kc);\n"
     "          shear_g = 16./45.*tau2_kc*(y[ppw->pv->index_pt_theta_g]+k2*ppw->pvecmetric[ppw->index_mt_alpha]);")
edit("source/perturbations.c",
     "          shear_g = 16./45./(ppw->pvecthermo[pth->index_th_dkappa] + ppw->pvecthermo[pth->index_th_dmu_idm_g])*y[ppw->pv->index_pt_theta_g];",
     "          double tau2_kc = 1./(ppw->pvecthermo[pth->index_th_dkappa] + ppw->pvecthermo[pth->index_th_dmu_idm_g]), dtau2_kc, pos_kc;\n"
     "          if ((pth->idm_g_quad_coeff != 1.) || (pth->idm_g_pol_coeff != 1.))\n"
     "            idm_g_kernel_closure(ppw->pvecthermo[pth->index_th_dkappa], ppw->pvecthermo[pth->index_th_dmu_idm_g], 0., 0.,\n"
     "                                 pth->idm_g_quad_coeff, pth->idm_g_pol_coeff, &tau2_kc, &dtau2_kc, &pos_kc);\n"
     "          shear_g = 16./45.*tau2_kc*y[ppw->pv->index_pt_theta_g];")

# tau_2 and its derivative in the TCA slip/shear function
edit("source/perturbations.c",
     "    tau_2_idm_g = 1./(pvecthermo[pth->index_th_dkappa] + dmu_idm_g);\n"
     "    dtau_2_idm_g = - (pvecthermo[pth->index_th_ddkappa] + ddmu_idm_g)*tau_2_idm_g*tau_2_idm_g;",
     "    tau_2_idm_g = 1./(pvecthermo[pth->index_th_dkappa] + dmu_idm_g);\n"
     "    dtau_2_idm_g = - (pvecthermo[pth->index_th_ddkappa] + ddmu_idm_g)*tau_2_idm_g*tau_2_idm_g;\n"
     "    if ((pth->idm_g_quad_coeff != 1.) || (pth->idm_g_pol_coeff != 1.)) {\n"
     "      double pos_kc;\n"
     "      idm_g_kernel_closure(pvecthermo[pth->index_th_dkappa], dmu_idm_g,\n"
     "                           pvecthermo[pth->index_th_ddkappa], ddmu_idm_g,\n"
     "                           pth->idm_g_quad_coeff, pth->idm_g_pol_coeff,\n"
     "                           &tau_2_idm_g, &dtau_2_idm_g, &pos_kc);\n"
     "    }")

# P during TCA in the source function
edit("source/perturbations.c",
     "        P = 5.* ppw->s_l[2] * ppw->tca_shear_g/8.; /* (2.5+0.5+2)shear_g/8 */",
     "        {\n"
     "          P = 5.* ppw->s_l[2] * ppw->tca_shear_g/8.; /* (2.5+0.5+2)shear_g/8 */\n"
     "          if ((pth->has_idm_g == _TRUE_) && ((pth->idm_g_quad_coeff != 1.) || (pth->idm_g_pol_coeff != 1.))) {\n"
     "            double tau2_kc, dtau2_kc, pos_kc;\n"
     "            idm_g_kernel_closure(pvecthermo[pth->index_th_dkappa], pvecthermo[pth->index_th_dmu_idm_g], 0., 0.,\n"
     "                                 pth->idm_g_quad_coeff, pth->idm_g_pol_coeff, &tau2_kc, &dtau2_kc, &pos_kc);\n"
     "            P = pos_kc * ppw->s_l[2] * ppw->tca_shear_g;\n"
     "          }\n"
     "        }")

# line-of-sight sources: Pi term (t2) in both gauges, and E source
edit("source/perturbations.c",
     "        _set_source_(ppt->index_tp_t2) = ppt->switch_pol * g * P;\n",
     "        _set_source_(ppt->index_tp_t2) = ppt->switch_pol * g_quad_kc * P;\n")
edit("source/perturbations.c",
     "        _set_source_(ppt->index_tp_t2) =\n          ppt->switch_pol * g * P;\n",
     "        _set_source_(ppt->index_tp_t2) =\n          ppt->switch_pol * g_quad_kc * P;\n")
edit("source/perturbations.c",
     "      _set_source_(ppt->index_tp_p) = _SQRT6_ * g * P;\n",
     "      _set_source_(ppt->index_tp_p) = _SQRT6_ * g_pol_kc * P;\n")
edit("source/perturbations.c",
     "  if (pth->has_idm_g == _TRUE_) {\n"
     "    dmu_idm_g = pvecthermo[pth->index_th_dmu_idm_g];\n"
     "    ddmu_idm_g = pvecthermo[pth->index_th_ddmu_idm_g];\n"
     "    exp_mu_idm_g = pvecthermo[pth->index_th_exp_mu_idm_g];\n"
     "  }\n",
     "  /* [class_baryons] visibility weighted by the kernel coefficients: the clump part of\n"
     "     g = (kappa' + mu') e^{-kappa-mu} contributes c_T (Pi term) or c_E (E source) */\n"
     "  double g_quad_kc = g, g_pol_kc = g;\n"
     "  if (pth->has_idm_g == _TRUE_) {\n"
     "    dmu_idm_g = pvecthermo[pth->index_th_dmu_idm_g];\n"
     "    ddmu_idm_g = pvecthermo[pth->index_th_ddmu_idm_g];\n"
     "    exp_mu_idm_g = pvecthermo[pth->index_th_exp_mu_idm_g];\n"
     "    if ((pth->idm_g_quad_coeff != 1.) || (pth->idm_g_pol_coeff != 1.)) {\n"
     "      double rate_kc = pvecthermo[pth->index_th_dkappa] + dmu_idm_g;\n"
     "      g_quad_kc = g * (pvecthermo[pth->index_th_dkappa] + pth->idm_g_quad_coeff*dmu_idm_g)/rate_kc;\n"
     "      g_pol_kc = g * (pvecthermo[pth->index_th_dkappa] + pth->idm_g_pol_coeff*dmu_idm_g)/rate_kc;\n"
     "    }\n"
     "  }\n")

# photon hierarchy (full equations, scalar modes)
edit("source/perturbations.c",
     "  double dmu_idm_g = 0., photon_scattering_rate;",
     "  double dmu_idm_g = 0., photon_scattering_rate;\n"
     "  double photon_regen_T, photon_regen_E; /* [class_baryons] kernel-weighted source rates */")
edit("source/perturbations.c",
     "  photon_scattering_rate = pvecthermo[pth->index_th_dkappa];\n",
     "  photon_scattering_rate = pvecthermo[pth->index_th_dkappa];\n"
     "  photon_regen_T = pvecthermo[pth->index_th_dkappa];\n"
     "  photon_regen_E = pvecthermo[pth->index_th_dkappa];\n")
edit("source/perturbations.c",
     "      photon_scattering_rate += pvecthermo[pth->index_th_dmu_idm_g];\n",
     "      photon_scattering_rate += pvecthermo[pth->index_th_dmu_idm_g];\n"
     "      photon_regen_T += pth->idm_g_quad_coeff * pvecthermo[pth->index_th_dmu_idm_g];\n"
     "      photon_regen_E += pth->idm_g_pol_coeff * pvecthermo[pth->index_th_dmu_idm_g];\n")
edit("source/perturbations.c",
     "  photon_regen_T = pvecthermo[pth->index_th_dkappa];\n"
     "  photon_regen_E = pvecthermo[pth->index_th_dkappa];\n",
     "  photon_regen_T = pvecthermo[pth->index_th_dkappa];\n"
     "  photon_regen_E = pvecthermo[pth->index_th_dkappa];\n"
     "  /* (the regeneration rates are reset to photon_scattering_rate below when the kernel is Thomson,\n"
     "     so that stock CLASS is reproduced bit for bit) */\n")
edit("source/perturbations.c",
     "            -photon_scattering_rate*(y[pv->index_pt_shear_g]-2./5./s_l[2]*P0);",
     "            -photon_scattering_rate*(y[pv->index_pt_shear_g]-2./5./s_l[2]*P0)\n"
     "            +(photon_regen_T-photon_scattering_rate)*2./5./s_l[2]*P0;",
     count=2)
edit("source/perturbations.c",
     "            -photon_scattering_rate*(y[pv->index_pt_pol0_g]-4.*P0);",
     "            -photon_scattering_rate*(y[pv->index_pt_pol0_g]-4.*P0)\n"
     "            +(photon_regen_E-photon_scattering_rate)*4.*P0;")
edit("source/perturbations.c",
     "            -photon_scattering_rate*(y[pv->index_pt_pol2_g]-4./5.*P0);",
     "            -photon_scattering_rate*(y[pv->index_pt_pol2_g]-4./5.*P0)\n"
     "            +(photon_regen_E-photon_scattering_rate)*4./5.*P0;")
edit("source/perturbations.c",
     "            -photon_scattering_rate*(y[pv->index_pt_E2] + _SQRT6_*P0);",
     "            -photon_scattering_rate*(y[pv->index_pt_E2] + _SQRT6_*P0)\n"
     "            -(photon_regen_E-photon_scattering_rate)*_SQRT6_*P0;")
edit("source/perturbations.c",
     "      photon_regen_E += pth->idm_g_pol_coeff * pvecthermo[pth->index_th_dmu_idm_g];\n",
     "      photon_regen_E += pth->idm_g_pol_coeff * pvecthermo[pth->index_th_dmu_idm_g];\n"
     "      if (pth->idm_g_quad_coeff == 1.) photon_regen_T = photon_scattering_rate;\n"
     "      if (pth->idm_g_pol_coeff == 1.) photon_regen_E = photon_scattering_rate;\n")
print("edits applied")
