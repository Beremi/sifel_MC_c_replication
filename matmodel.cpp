#include "matmodel.h"
#include "vector.h"
#include "matrix.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
  This file is a one-material-point rewrite of the Matlab pair

    matlab_export/constitutive_problem.m
    matlab_export/stiffness_matrix.m

  The Matlab code is vectorized over all elements. Here the same logic is evaluated
  for one material point directly inside the original-template API functions.
*/

/**
  The function reads material model parameters and stress return setup from the opened text.

  @param[in] in - pointer to the opened text file.

  @retval 0 - on success
  @retval 1 - in the case of an error
*/
long matmodel::read(FILE *in)
{
  matmodel_params tmp;
  double young, poisson, cohesion, phi;

  if (in == NULL)
    return 1;

  if (fscanf(in, "%lf %lf %lf %lf", &young, &poisson, &cohesion, &phi) != 4)
    return 1;

  // Precompute the elastic constants and Mohr-Coulomb parameters used in S and DS.
  tmp.young = young;
  tmp.poisson = poisson;
  tmp.c = cohesion;
  tmp.phi = phi;
  tmp.shear = young/(2.0*(1.0 + poisson));
  tmp.bulk = young/(3.0*(1.0 - 2.0*poisson));
  tmp.lame = tmp.bulk - 2.0*tmp.shear/3.0;
  tmp.sin_phi = sin(phi);
  tmp.cos_phi = cos(phi);
  tmp.c_bar = 2.0*cohesion*tmp.cos_phi;

  if (tmp.young <= 0.0)
    return 1;
  if ((tmp.poisson <= -1.0) || (tmp.poisson >= 0.5))
    return 1;
  if (tmp.c < 0.0)
    return 1;
  if ((tmp.phi <= 1.0e-14) || (tmp.phi >= 0.5*M_PI - 1.0e-14))
    return 1;

  par = tmp;
  return 0;
}

/**
  The function prints material model parameters and stress return setup to the opened text
  file.

  @param[in] out - pointer to the opened text output file.

  @return The function does not return anything.
*/
void matmodel::print(FILE *out)
{
  if (out == NULL)
    return;

  fprintf(out, "MC perfect plastic associative 2D model\n");
  fprintf(out, "  E      = %.15g\n", par.young);
  fprintf(out, "  nu     = %.15g\n", par.poisson);
  fprintf(out, "  c      = %.15g\n", par.c);
  fprintf(out, "  phi    = %.15g rad\n", par.phi);
  fprintf(out, "  G      = %.15g\n", par.shear);
  fprintf(out, "  K      = %.15g\n", par.bulk);
  fprintf(out, "  lambda = %.15g\n", par.lame);
  fprintf(out, "  other  = %d components\n", MATMODEL_NCOMP_OTHER);
  fprintf(out, "  eqother= %d components\n", MATMODEL_NCOMP_EQOTHER);
  fprintf(out, "  buffer = epsp(4), return_type(1), eig(3), proj(12), hess_red(27), sigma_princ(3)\n");
}

/**
  The function computes actual stresses with respect to the attained strains and state
  variables.

  @param[in] strain - array of actual strain, components are ordered as follows
                      eps_x, eps_y, gamma_xy, eps_z - for the plane strain problem
                      eps_x, eps_y, eps_z, gamma_yz, gamma_xz, gamma_xy - for the space
                      stress problem

  @param[in] eqstatev - packed point buffer, where the first 4 entries contain the
                        incoming plastic strain history; remaining entries are ignored
                        by nlstresses()

  @param[out] stress - array of the resulting stress components, it must be computed,
                       components are ordered as follows
                       sig_x, sig_y, tau_xy, sig_z - for the plane strain problem
                       sig_x, sig_y, sig_z, tau_yz, tau_xz, tau_xy - for the space
                       stress problem

  @param[out] statev - array of the resulting packed point buffer for the actual strains

  @return The function does not return anything, the results are stored in the arrays
          stress and statev passed in as arguments.
*/
void matmodel::nlstresses(const vector &strain,
                          const vector &eqstatev,
                          vector &stress,
                          vector &statev)
{
  /*
    Scalar material-point version of Matlab:

      [S,eig_1,eig_2,eig_3,Eig_1,Eig_2,Eig_3,EIG_1,EIG_2,EIG_3,...
        sigma_1,sigma_2,sigma_3] = constitutive_problem(E_new, Ep_prev)

    plus the plastic-strain recovery used later in matlab_export/newton.m:

      Ep = -Inv_ELAST*S;
      Ep(1:3,:) = Ep(1:3,:) + E;
  */
  vector E(MATMODEL_NCOMP_STRAIN), Ep_prev(4), E_trial(MATMODEL_NCOMP_STRAIN);
  vector Ep(4), eps_el(MATMODEL_NCOMP_STRAIN), eig(3), sigma(3), S(MATMODEL_NCOMP_STRESS);
  vector Eig0_1(MATMODEL_NCOMP_STRAIN), Eig0_2(MATMODEL_NCOMP_STRAIN), Eig0_3(MATMODEL_NCOMP_STRAIN);
  vector Eig_1(MATMODEL_NCOMP_STRAIN), Eig_2(MATMODEL_NCOMP_STRAIN), Eig_3(MATMODEL_NCOMP_STRAIN);
  matrix metric(MATMODEL_NCOMP_STRAIN, MATMODEL_NCOMP_STRAIN);
  matrix EIG0_1(MATMODEL_NCOMP_STRAIN, MATMODEL_NCOMP_STRAIN);
  matrix EIG0_2(MATMODEL_NCOMP_STRAIN, MATMODEL_NCOMP_STRAIN);
  matrix EIG0_3(MATMODEL_NCOMP_STRAIN, MATMODEL_NCOMP_STRAIN);
  matrix EIG_1(MATMODEL_NCOMP_STRAIN, MATMODEL_NCOMP_STRAIN);
  matrix EIG_2(MATMODEL_NCOMP_STRAIN, MATMODEL_NCOMP_STRAIN);
  matrix EIG_3(MATMODEL_NCOMP_STRAIN, MATMODEL_NCOMP_STRAIN);
  double eig0_1, eig0_2, eig0_3;
  double trace_E, f_tr, gamma_sl, gamma_sr, gamma_la, gamma_ra;
  double denom_s, denom_l, denom_r;
  double lambda_s, lambda_l, lambda_r;
  int return_type = MATMODEL_RETURN_ELASTIC;
  long idx;

  // Resize the output arrays before filling this material point.
  reallocv(MATMODEL_NCOMP_STRESS, stress);
  reallocv(MATMODEL_NCOMP_OTHER, statev);
  nullv(E);
  nullv(Ep_prev);
  nullv(E_trial);
  nullv(Ep);
  nullv(eps_el);
  nullv(eig);
  nullv(sigma);
  nullv(S);
  nullv(Eig0_1);
  nullv(Eig0_2);
  nullv(Eig0_3);
  nullv(Eig_1);
  nullv(Eig_2);
  nullv(Eig_3);
  nullm(metric);
  nullm(EIG0_1);
  nullm(EIG0_2);
  nullm(EIG0_3);
  nullm(EIG_1);
  nullm(EIG_2);
  nullm(EIG_3);

  // Read the current strain and only the leading plastic-strain entries from the fixed buffer layout.
  copyv(strain, E);
  copyv(&eqstatev[MATMODEL_IO_EP_XX], Ep_prev);

  /*
    Matlab:
      E_trial = -Ep_prev;
      E_trial(1:3,:) = E_trial(1:3,:) + E_new;
  */
  addmultv(E, 1.0, Ep_prev, -1.0, E_trial);

  // This is the reduced identity metric used in the eigenvalue Hessian formula.
  metric(0,0) = 1.0;
  metric(1,1) = 1.0;
  metric(2,2) = 0.5;

  /*
    Matlab block:

      I1 = E_trial(1,:) + E_trial(2,:);
      I2 = sqrt((E_trial(1,:)-E_trial(2,:)).^2 + E_trial(3,:).^2);
      eig0_1 = (I1+I2)/2;
      eig0_2 = (I1-I2)/2;
      eig0_3 = E_trial(4,:);

      Eig0_1 = ...
      Eig0_2 = ...
      Eig0_3 = ...

      EIG0_1 = ...
      EIG0_2 = -EIG0_1;
      EIG0_3 = 0
  */
  {
    const double ex = E_trial[0];
    const double ey = E_trial[1];
    const double g = E_trial[2];
    const double d = ex - ey;
    const double i1 = ex + ey;
    const double i2 = sqrt(d*d + g*g);
    const double tol = 1.0e-14;

    eig0_1 = 0.5*(i1 + i2);
    eig0_2 = 0.5*(i1 - i2);
    eig0_3 = E_trial[3];

    // Matlab uses a special fallback when the first two trial eigenvalues coincide.
    if (i2 <= tol)
    {
      Eig0_1[0] = 1.0;
      Eig0_1[1] = 1.0;
    }
    else
    {
      // Compute the in-plane eigenprojections of the trial strain.
      Eig0_1[0] = (ex - eig0_2)/i2;
      Eig0_1[1] = (ey - eig0_2)/i2;
      Eig0_1[2] = g/(2.0*i2);

      Eig0_2[0] = 1.0 - Eig0_1[0];
      Eig0_2[1] = 1.0 - Eig0_1[1];
      Eig0_2[2] = -Eig0_1[2];

      // Assemble the reduced Hessians of the in-plane trial eigenvalues.
      for (long a=0; a<MATMODEL_NCOMP_STRAIN; a++)
      {
        for (long b=0; b<MATMODEL_NCOMP_STRAIN; b++)
        {
          EIG0_1(a,b) = (metric(a,b) - Eig0_1[a]*Eig0_1[b] - Eig0_2[a]*Eig0_2[b])/i2;
          EIG0_2(a,b) = -EIG0_1(a,b);
        }
      }

      // The zz direction is decoupled, so its Hessian rows and columns stay zero.
      for (long a=0; a<MATMODEL_NCOMP_STRAIN; a++)
      {
        EIG0_1(3,a) = 0.0;
        EIG0_1(a,3) = 0.0;
        EIG0_2(3,a) = 0.0;
        EIG0_2(a,3) = 0.0;
      }
    }

    Eig0_3[3] = 1.0;
  }

  /*
    Matlab:
      test2 = (eig0_1>=eig0_3)&(eig0_3>eig0_2);
      test3 = (eig0_3>eig0_1);
  */
  // Start from the natural order and swap when the zz mode crosses the in-plane modes.
  eig[0] = eig0_1;
  eig[1] = eig0_2;
  eig[2] = eig0_3;
  copyv(Eig0_1, Eig_1);
  copyv(Eig0_2, Eig_2);
  copyv(Eig0_3, Eig_3);
  copym(EIG0_1, EIG_1);
  copym(EIG0_2, EIG_2);
  copym(EIG0_3, EIG_3);

  if ((eig0_1 >= eig0_3) && (eig0_3 > eig0_2))
  {
    eig[1] = eig0_3;
    eig[2] = eig0_2;
    copyv(Eig0_3, Eig_2);
    copyv(Eig0_2, Eig_3);
    copym(EIG0_3, EIG_2);
    copym(EIG0_2, EIG_3);
  }
  if (eig0_3 > eig0_1)
  {
    eig[0] = eig0_3;
    eig[1] = eig0_1;
    eig[2] = eig0_2;
    copyv(Eig0_3, Eig_1);
    copyv(Eig0_1, Eig_2);
    copyv(Eig0_2, Eig_3);
    copym(EIG0_3, EIG_1);
    copym(EIG0_1, EIG_2);
    copym(EIG0_2, EIG_3);
  }

  /*
    Matlab:
      trace_E = eig_1 + eig_2 + eig_3;
      f_tr    = ...
      gamma_sl = ...
      gamma_sr = ...
      gamma_la = ...
      gamma_ra = ...

      denom_s = ...
      denom_l = ...
      denom_r = ...

      lambda_s = ...
      lambda_l = ...
      lambda_r = ...
  */
  // Evaluate the Mohr-Coulomb trial function and the return candidates.
  trace_E = eig[0] + eig[1] + eig[2];
  f_tr = 2.0*par.shear*((1.0 + par.sin_phi)*eig[0] - (1.0 - par.sin_phi)*eig[2]) + 2.0*par.lame*par.sin_phi*trace_E - par.c_bar;
  gamma_sl = (eig[0] - eig[1])/(1.0 + par.sin_phi);
  gamma_sr = (eig[1] - eig[2])/(1.0 - par.sin_phi);
  gamma_la = (eig[0] + eig[1] - 2.0*eig[2])/(3.0 - par.sin_phi);
  gamma_ra = (2.0*eig[0] - eig[1] - eig[2])/(3.0 + par.sin_phi);

  denom_s = 4.0*par.lame*par.sin_phi*par.sin_phi + 2.0*par.shear*(1.0 + par.sin_phi)*(1.0 + par.sin_phi) + 2.0*par.shear*(1.0 - par.sin_phi)*(1.0 - par.sin_phi);
  denom_l = 4.0*par.lame*par.sin_phi*par.sin_phi + par.shear*(1.0 + par.sin_phi)*(1.0 + par.sin_phi) + 2.0*par.shear*(1.0 - par.sin_phi)*(1.0 - par.sin_phi);
  denom_r = 4.0*par.lame*par.sin_phi*par.sin_phi + 2.0*par.shear*(1.0 + par.sin_phi)*(1.0 + par.sin_phi) + par.shear*(1.0 - par.sin_phi)*(1.0 - par.sin_phi);

  lambda_s = f_tr/denom_s;
  lambda_l = (par.shear*((1.0 + par.sin_phi)*(eig[0] + eig[1]) - 2.0*(1.0 - par.sin_phi)*eig[2]) + 2.0*par.lame*par.sin_phi*trace_E - par.c_bar)/denom_l;
  lambda_r = (par.shear*(2.0*(1.0 + par.sin_phi)*eig[0] - (1.0 - par.sin_phi)*(eig[1] + eig[2])) + 2.0*par.lame*par.sin_phi*trace_E - par.c_bar)/denom_r;

  /*
    Matlab return classification:

      test_el = (f_tr<=0)
      test_s  = ...
      test_l  = ...
      test_r  = ...
      test_a  = ~(test_el|test_s|test_l|test_r)
  */
  // Choose the return branch and compute the principal stresses after return mapping.
  if (f_tr <= 0.0)
  {
    // The trial state is admissible, so the response stays elastic.
    return_type = MATMODEL_RETURN_ELASTIC;
    sigma[0] = par.lame*trace_E + 2.0*par.shear*eig[0];
    sigma[1] = par.lame*trace_E + 2.0*par.shear*eig[1];
    sigma[2] = par.lame*trace_E + 2.0*par.shear*eig[2];
  }
  else if (lambda_s <= ((gamma_sl < gamma_sr) ? gamma_sl : gamma_sr))
  {
    // Smooth-face return updates all three principal stresses separately.
    return_type = MATMODEL_RETURN_SMOOTH;
    sigma[0] = par.lame*trace_E + 2.0*par.shear*eig[0] - lambda_s*(2.0*par.lame*par.sin_phi + 2.0*par.shear*(1.0 + par.sin_phi));
    sigma[1] = par.lame*trace_E + 2.0*par.shear*eig[1] - lambda_s*(2.0*par.lame*par.sin_phi);
    sigma[2] = par.lame*trace_E + 2.0*par.shear*eig[2] - lambda_s*(2.0*par.lame*par.sin_phi - 2.0*par.shear*(1.0 - par.sin_phi));
  }
  else if ((gamma_sl < gamma_sr) && (lambda_l >= gamma_sl) && (lambda_l <= gamma_la))
  {
    // Left-edge return collapses the first two principal stresses to one value.
    return_type = MATMODEL_RETURN_LEFT_EDGE;
    sigma[0] = par.lame*trace_E + par.shear*(eig[0] + eig[1]) - lambda_l*(2.0*par.lame*par.sin_phi + par.shear*(1.0 + par.sin_phi));
    sigma[1] = sigma[0];
    sigma[2] = par.lame*trace_E + 2.0*par.shear*eig[2] - lambda_l*(2.0*par.lame*par.sin_phi - 2.0*par.shear*(1.0 - par.sin_phi));
  }
  else if ((gamma_sl > gamma_sr) && (lambda_r >= gamma_sr) && (lambda_r <= gamma_ra))
  {
    // Right-edge return collapses the last two principal stresses to one value.
    return_type = MATMODEL_RETURN_RIGHT_EDGE;
    sigma[0] = par.lame*trace_E + 2.0*par.shear*eig[0] - lambda_r*(2.0*par.lame*par.sin_phi + 2.0*par.shear*(1.0 + par.sin_phi));
    sigma[2] = par.lame*trace_E + par.shear*(eig[1] + eig[2]) - lambda_r*(2.0*par.lame*par.sin_phi - par.shear*(1.0 - par.sin_phi));
    sigma[1] = sigma[2];
  }
  else
  {
    // Apex return gives a hydrostatic principal stress state.
    return_type = MATMODEL_RETURN_APEX;
    sigma[0] = par.c_bar/(2.0*par.sin_phi);
    sigma[1] = sigma[0];
    sigma[2] = sigma[0];
  }

  /*
    Matlab stress reconstruction:

      elastic: S = ELAST*E_trial
      smooth:  S = sigma_1*Eig_1 + sigma_2*Eig_2 + sigma_3*Eig_3
      left:    S = sigma_1*(Eig_1+Eig_2) + sigma_3*Eig_3
      right:   S = sigma_1*Eig_1 + sigma_3*(Eig_2+Eig_3)
      apex:    S = iota*sigma_1
  */
  // Rebuild the full stress vector from the principal stresses and projections.
  nullv(S);
  addmultv(S, Eig_1, sigma[0]);
  addmultv(S, Eig_2, sigma[1]);
  addmultv(S, Eig_3, sigma[2]);

  /*
    Matlab computes the updated plastic strain in newton.m after the constitutive call:

      Ep = -Inv_ELAST*S;
      Ep(1:3,:) = Ep(1:3,:) + E;
  */
  // Convert the stress vector back to elastic strain and recover the plastic part.
  eps_el[0] = (S[0] - par.poisson*(S[1] + S[3]))/par.young;
  eps_el[1] = (S[1] - par.poisson*(S[0] + S[3]))/par.young;
  eps_el[2] = S[2]/par.shear;
  eps_el[3] = (S[3] - par.poisson*(S[0] + S[1]))/par.young;
  addmultv(E, 1.0, eps_el, -1.0, Ep);

  /*
    statev is the packed per-point payload used by the subsequent tangent call:

      epsp, return_type, eig_*, Eig_*, EIG_*, sigma_*
  */
  // Start with a clean point buffer before packing the current constitutive payload.
  nullv(statev);

  // Store the updated plastic strain for the next stress update.
  copyv(Ep, &statev[MATMODEL_IO_EP_XX]);

  // Store the chosen return type and the ordered trial eigenvalues.
  statev[MATMODEL_IO_RETURN_TYPE] = static_cast<double>(return_type);
  statev[MATMODEL_IO_EIG_1] = eig[0];
  statev[MATMODEL_IO_EIG_2] = eig[1];
  statev[MATMODEL_IO_EIG_3] = eig[2];

  // Store the first derivatives of the ordered trial eigenvalues.
  copyv(Eig_1, &statev[MATMODEL_IO_PROJ_1]);
  copyv(Eig_2, &statev[MATMODEL_IO_PROJ_2]);
  copyv(Eig_3, &statev[MATMODEL_IO_PROJ_3]);

  // Pack the reduced 3x3 Hessian blocks in Matlab column-major order.
  idx = MATMODEL_IO_HESS_1;
  for (long j=0; j<3; j++)
  {
    for (long i=0; i<3; i++)
      statev[idx++] = EIG_1(i,j);
  }
  for (long j=0; j<3; j++)
  {
    for (long i=0; i<3; i++)
      statev[idx++] = EIG_2(i,j);
  }
  for (long j=0; j<3; j++)
  {
    for (long i=0; i<3; i++)
      statev[idx++] = EIG_3(i,j);
  }

  // Store the returned principal stresses and the reconstructed stress vector.
  copyv(sigma, &statev[MATMODEL_IO_SIGMA_1]);
  copyv(S, stress);
}

/**
  The function computes material stiffness matrix with respect to the attained strains
  and state variables.

  @param[in] strain - array of actual strain, components are ordered as follows
                      eps_x, eps_y, gamma_xy, eps_z - for the plane strain problem
                      eps_x, eps_y, eps_z, gamma_yz, gamma_xz, gamma_xy - for the space
                      stress problem

  @param[in] eqstatev - packed per-point buffer produced by nlstresses(); this buffer
                        contains the same auxiliary data that Matlab passes from
                        constitutive_problem.m to stiffness_matrix.m

  @param[in] stress - array of the resulting stress components, it must be computed,
                      components are ordered as follows
                      sig_x, sig_y, tau_xy, sig_z - for the plane strain problem
                      sig_x, sig_y, sig_z, tau_yz, tau_xz, tau_xy - for the space
                      stress problem

  @param[out] d - the resulting material stiffness matrix, it must be calculated.

  @return The function does not return anything, the resulting matrix is stored in the
          matrix d passed in as an argument.
*/
void matmodel::stiffmat(const vector &strain,
                        const vector &eqstatev,
                        const vector &stress,
                        matrix &d)
{
  /*
    Scalar material-point version of Matlab:

      Sderiv = stiffness_matrix(eig_1,eig_2,eig_3,Eig_1,Eig_2,Eig_3,...
                                EIG_1,EIG_2,EIG_3,sigma_1,sigma_2,sigma_3)

    In this SIFEL-style rewrite, eqstatev is the packed point buffer slice that carries
    the auxiliary constitutive data for the current point. stiffmat() therefore reads
    eig_*, Eig_*, EIG_*, sigma_* directly from eqstatev and does not recompute the
    return-map classification.
  */
  vector sigma(3), iota(4), Eig_1(MATMODEL_NCOMP_STRAIN), Eig_2(MATMODEL_NCOMP_STRAIN), Eig_3(MATMODEL_NCOMP_STRAIN), Eig_12(MATMODEL_NCOMP_STRAIN), Eig_23(MATMODEL_NCOMP_STRAIN), Eig_6(MATMODEL_NCOMP_STRAIN);
  matrix EIG_1(MATMODEL_NCOMP_STRAIN, MATMODEL_NCOMP_STRAIN), EIG_2(MATMODEL_NCOMP_STRAIN, MATMODEL_NCOMP_STRAIN), EIG_3(MATMODEL_NCOMP_STRAIN, MATMODEL_NCOMP_STRAIN), EIG_12(MATMODEL_NCOMP_STRAIN, MATMODEL_NCOMP_STRAIN), EIG_23(MATMODEL_NCOMP_STRAIN, MATMODEL_NCOMP_STRAIN);
  matrix zero(MATMODEL_NCOMP_STRAIN, MATMODEL_NCOMP_STRAIN), outer(MATMODEL_NCOMP_STRAIN, MATMODEL_NCOMP_STRAIN);
  double eig_1, eig_2, eig_3;
  double denom_s, denom_l, denom_r;
  int return_type = MATMODEL_RETURN_ELASTIC;
  long idx;
  (void)strain;
  (void)stress;

  // Allocate the material tangent and the local work arrays used in the assembly.
  reallocm(MATMODEL_NCOMP_STRAIN, MATMODEL_NCOMP_STRAIN, d);
  nullm(d);
  nullv(sigma);
  nullv(iota);
  nullv(Eig_1);
  nullv(Eig_2);
  nullv(Eig_3);
  nullv(Eig_12);
  nullv(Eig_23);
  nullv(Eig_6);
  nullm(EIG_1);
  nullm(EIG_2);
  nullm(EIG_3);
  nullm(EIG_12);
  nullm(EIG_23);
  nullm(zero);
  nullm(outer);

  // Read the stored ordered trial eigenvalues from their fixed positions in the point buffer.
  eig_1 = eqstatev[MATMODEL_IO_EIG_1];
  eig_2 = eqstatev[MATMODEL_IO_EIG_2];
  eig_3 = eqstatev[MATMODEL_IO_EIG_3];
  (void)eig_1;
  (void)eig_2;
  (void)eig_3;

  // Read the stored return branch directly from the packed point buffer.
  return_type = static_cast<int>(eqstatev[MATMODEL_IO_RETURN_TYPE] + 0.5);

  // Unpack the stored eigenprojections of the ordered trial strain.
  copyv(&eqstatev[MATMODEL_IO_PROJ_1], Eig_1);
  copyv(&eqstatev[MATMODEL_IO_PROJ_2], Eig_2);
  copyv(&eqstatev[MATMODEL_IO_PROJ_3], Eig_3);

  // Rebuild the full 4x4 Hessian blocks from the packed reduced 3x3 Matlab layout.
  idx = MATMODEL_IO_HESS_1;
  for (long j=0; j<3; j++)
  {
    for (long i=0; i<3; i++)
    {
      EIG_1(i,j) = eqstatev[idx];
      idx++;
    }
  }
  for (long j=0; j<3; j++)
  {
    for (long i=0; i<3; i++)
    {
      EIG_2(i,j) = eqstatev[idx];
      idx++;
    }
  }
  for (long j=0; j<3; j++)
  {
    for (long i=0; i<3; i++)
    {
      EIG_3(i,j) = eqstatev[idx];
      idx++;
    }
  }

  // Read the returned principal stresses from their fixed positions in the point buffer.
  sigma[0] = eqstatev[MATMODEL_IO_SIGMA_1];
  sigma[1] = eqstatev[MATMODEL_IO_SIGMA_2];
  sigma[2] = eqstatev[MATMODEL_IO_SIGMA_3];

  // Only the branch denominators depend on material constants alone and are rebuilt here.
  denom_s = 4.0*par.lame*par.sin_phi*par.sin_phi + 2.0*par.shear*(1.0 + par.sin_phi)*(1.0 + par.sin_phi) + 2.0*par.shear*(1.0 - par.sin_phi)*(1.0 - par.sin_phi);
  denom_l = 4.0*par.lame*par.sin_phi*par.sin_phi + par.shear*(1.0 + par.sin_phi)*(1.0 + par.sin_phi) + 2.0*par.shear*(1.0 - par.sin_phi)*(1.0 - par.sin_phi);
  denom_r = 4.0*par.lame*par.sin_phi*par.sin_phi + 2.0*par.shear*(1.0 + par.sin_phi)*(1.0 + par.sin_phi) + par.shear*(1.0 - par.sin_phi)*(1.0 - par.sin_phi);

  // The plane-strain unit tensor enters the volumetric part of every tangent branch.
  iota[0] = 1.0;
  iota[1] = 1.0;
  iota[3] = 1.0;

  // Assemble the 4x4 SIFEL tangent from the already packed Matlab-style state.
  switch (return_type)
  {
    case MATMODEL_RETURN_ELASTIC:
      // Elastic tangent is the constant isotropic stiffness matrix.
      d(0,0) = par.lame + 2.0*par.shear;
      d(0,1) = par.lame;
      d(0,3) = par.lame;
      d(1,0) = par.lame;
      d(1,1) = par.lame + 2.0*par.shear;
      d(1,3) = par.lame;
      d(2,2) = par.shear;
      d(3,0) = par.lame;
      d(3,1) = par.lame;
      d(3,3) = par.lame + 2.0*par.shear;
      break;

    case MATMODEL_RETURN_SMOOTH:
      // Smooth-face return adds the Hessian terms and the consistency correction.
      addmultm(zero, 0.0, EIG_1, sigma[0], d);
      addmultm(zero, 0.0, EIG_2, sigma[1], d);
      addmultm(zero, 0.0, EIG_3, sigma[2], d);
      vxv(iota, iota, outer);
      addmultm(zero, 0.0, outer, par.lame, d);
      vxv(Eig_1, Eig_1, outer);
      addmultm(zero, 0.0, outer, 2.0*par.shear, d);
      vxv(Eig_2, Eig_2, outer);
      addmultm(zero, 0.0, outer, 2.0*par.shear, d);
      vxv(Eig_3, Eig_3, outer);
      addmultm(zero, 0.0, outer, 2.0*par.shear, d);

      // Eig_6 is Matlab's consistency-direction vector for the smooth return.
      addmultv(Eig_1, 2.0*par.shear*(1.0 + par.sin_phi), Eig_3, -2.0*par.shear*(1.0 - par.sin_phi), Eig_6);
      addmultv(Eig_6, iota, 2.0*par.lame*par.sin_phi);
      vxv(Eig_6, Eig_6, outer);
      addmultm(zero, 0.0, outer, -1.0/denom_s, d);
      break;

    case MATMODEL_RETURN_LEFT_EDGE:
      // Left-edge return uses the merged first-two-eigenvalue branch.
      addmultv(Eig_1, 1.0, Eig_2, 1.0, Eig_12);
      addm(EIG_1, EIG_2, EIG_12);
      addmultm(zero, 0.0, EIG_12, sigma[0], d);
      addmultm(zero, 0.0, EIG_3, sigma[2], d);
      vxv(iota, iota, outer);
      addmultm(zero, 0.0, outer, par.lame, d);
      vxv(Eig_12, Eig_12, outer);
      addmultm(zero, 0.0, outer, par.shear, d);
      vxv(Eig_3, Eig_3, outer);
      addmultm(zero, 0.0, outer, 2.0*par.shear, d);

      // Eig_6 is Matlab's consistency-direction vector for the left-edge return.
      addmultv(Eig_12, par.shear*(1.0 + par.sin_phi), Eig_3, -2.0*par.shear*(1.0 - par.sin_phi), Eig_6);
      addmultv(Eig_6, iota, 2.0*par.lame*par.sin_phi);
      vxv(Eig_6, Eig_6, outer);
      addmultm(zero, 0.0, outer, -1.0/denom_l, d);
      break;

    case MATMODEL_RETURN_RIGHT_EDGE:
      // Right-edge return uses the merged last-two-eigenvalue branch.
      addmultv(Eig_2, 1.0, Eig_3, 1.0, Eig_23);
      addm(EIG_2, EIG_3, EIG_23);
      addmultm(zero, 0.0, EIG_1, sigma[0], d);
      addmultm(zero, 0.0, EIG_23, sigma[2], d);
      vxv(iota, iota, outer);
      addmultm(zero, 0.0, outer, par.lame, d);
      vxv(Eig_1, Eig_1, outer);
      addmultm(zero, 0.0, outer, 2.0*par.shear, d);
      vxv(Eig_23, Eig_23, outer);
      addmultm(zero, 0.0, outer, par.shear, d);

      // Eig_6 is Matlab's consistency-direction vector for the right-edge return.
      addmultv(Eig_1, 2.0*par.shear*(1.0 + par.sin_phi), Eig_23, -par.shear*(1.0 - par.sin_phi), Eig_6);
      addmultv(Eig_6, iota, 2.0*par.lame*par.sin_phi);
      vxv(Eig_6, Eig_6, outer);
      addmultm(zero, 0.0, outer, -1.0/denom_r, d);
      break;

    case MATMODEL_RETURN_APEX:
    default:
      // Apex tangent is zero in the Matlab formulation.
      break;
  }
}
