#include "mohrc3d_ugn.h"
#include "vector.h"
#include "matrix.h"
#include "iotools.h"
#include "global.h"
#include "mechmat.h"

#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
  This file is a one-material-point rewrite of the Matlab pair

    matlab_export_3d/constitutive_problem.m
    matlab_export_3d/stiffness_matrix.m

  The Matlab code is vectorized over all elements. Here the same logic is evaluated
  for one material point directly inside the original-template API functions.

  The 3D components are kept in the SIFEL ordering throughout:
  [xx, yy, zz, yz, xz, xy].
*/

/**
  The function computes the ordered trial eigenvalues together with the tensor
  arrays needed later by the 3D spectral projection formulas.

  @param[in] E_trial - trial strain vector in the SIFEL 3D ordering
  @param[out] E_tr - trial strain in stress notation
  @param[out] E_square - square of the trial strain in stress notation
  @param[out] DER_E_square - derivative of the squared trial strain tensor
  @param[out] eig - ordered trial eigenvalues

  @return The function does not return anything, the results are stored in the
          output arguments.
*/
void mohrc3d_ugn::compute_ordered_trial_spectral(const vector &E_trial,
                                                 vector &E_tr,
                                                 vector &E_square,
                                                 matrix &DER_E_square,
                                                 vector &eig) const
{
  double Q_raw, Q, R, theta0, theta;
  double I1, I2, I3;

  // Clear the spectral work arrays before filling the current trial state.
  nullv(E_tr);
  nullv(E_square);
  nullm(DER_E_square);
  nullv(eig);

  // Trial strain in stress notation and its square.
  E_tr[0] = E_trial[0];
  E_tr[1] = E_trial[1];
  E_tr[2] = E_trial[2];
  E_tr[3] = 0.5 * E_trial[3];
  E_tr[4] = 0.5 * E_trial[4];
  E_tr[5] = 0.5 * E_trial[5];

  E_square[0] = E_tr[0]*E_tr[0] + E_tr[5]*E_tr[5] + E_tr[4]*E_tr[4];
  E_square[1] = E_tr[1]*E_tr[1] + E_tr[5]*E_tr[5] + E_tr[3]*E_tr[3];
  E_square[2] = E_tr[2]*E_tr[2] + E_tr[4]*E_tr[4] + E_tr[3]*E_tr[3];
  E_square[3] = E_tr[5]*E_tr[4] + E_tr[1]*E_tr[3] + E_tr[2]*E_tr[3];
  E_square[4] = E_tr[0]*E_tr[4] + E_tr[5]*E_tr[3] + E_tr[2]*E_tr[4];
  E_square[5] = E_tr[0]*E_tr[5] + E_tr[1]*E_tr[5] + E_tr[4]*E_tr[3];

  // Derivative of the square of the trial strain tensor.
  DER_E_square(0,0) = 2.0 * E_tr[0];
  DER_E_square(4,0) = E_tr[4];
  DER_E_square(5,0) = E_tr[5];

  DER_E_square(1,1) = 2.0 * E_tr[1];
  DER_E_square(3,1) = E_tr[3];
  DER_E_square(5,1) = E_tr[5];

  DER_E_square(2,2) = 2.0 * E_tr[2];
  DER_E_square(3,2) = E_tr[3];
  DER_E_square(4,2) = E_tr[4];

  DER_E_square(1,3) = E_tr[3];
  DER_E_square(2,3) = E_tr[3];
  DER_E_square(3,3) = 0.5 * (E_tr[1] + E_tr[2]);
  DER_E_square(4,3) = 0.5 * E_tr[5];
  DER_E_square(5,3) = 0.5 * E_tr[4];

  DER_E_square(0,4) = E_tr[4];
  DER_E_square(2,4) = E_tr[4];
  DER_E_square(3,4) = 0.5 * E_tr[5];
  DER_E_square(4,4) = 0.5 * (E_tr[0] + E_tr[2]);
  DER_E_square(5,4) = 0.5 * E_tr[3];

  DER_E_square(0,5) = E_tr[5];
  DER_E_square(1,5) = E_tr[5];
  DER_E_square(3,5) = 0.5 * E_tr[4];
  DER_E_square(4,5) = 0.5 * E_tr[3];
  DER_E_square(5,5) = 0.5 * (E_tr[0] + E_tr[1]);

  // Ordered eigenvalues of the trial strain.
  I1 = E_tr[0] + E_tr[1] + E_tr[2];
  I2 = E_tr[0]*E_tr[1] + E_tr[0]*E_tr[2] + E_tr[1]*E_tr[2] -
       E_tr[3]*E_tr[3] - E_tr[4]*E_tr[4] - E_tr[5]*E_tr[5];
  I3 = E_tr[0]*E_tr[1]*E_tr[2] - E_tr[0]*E_tr[3]*E_tr[3] -
       E_tr[1]*E_tr[4]*E_tr[4] - E_tr[2]*E_tr[5]*E_tr[5] +
       2.0*E_tr[3]*E_tr[4]*E_tr[5];

  Q_raw = (I1*I1 - 3.0*I2) / 9.0;
  Q = (Q_raw > 0.0) ? Q_raw : 0.0;
  R = (-2.0*I1*I1*I1 + 9.0*I1*I2 - 27.0*I3) / 54.0;
  theta0 = 0.0;

  if (Q > 0.0)
    theta0 = R / sqrt(Q*Q*Q);
  if (theta0 < -1.0)
    theta0 = -1.0;
  if (theta0 > 1.0)
    theta0 = 1.0;

  theta = acos(theta0) / 3.0;
  eig[0] = -2.0*sqrt(Q)*cos(theta + 2.0*M_PI/3.0) + I1/3.0;
  eig[1] = -2.0*sqrt(Q)*cos(theta - 2.0*M_PI/3.0) + I1/3.0;
  eig[2] = -2.0*sqrt(Q)*cos(theta) + I1/3.0;
}



/**
  The function computes the return denominators used by the smooth, edge, and
  apex branches of the 3D Mohr-Coulomb return mapping.

  @param[out] denom_s - denominator for the smooth-face return
  @param[out] denom_l - denominator for the left-edge return
  @param[out] denom_r - denominator for the right-edge return
  @param[out] denom_a - denominator for the apex return

  @return The function does not return anything, the results are stored in the
          output arguments.
*/
void mohrc3d_ugn::compute_return_denominators(double &denom_s,
                                              double &denom_l,
                                              double &denom_r,
                                              double &denom_a) const
{
  // These branch constants depend only on the material parameters.
  denom_s = 4.0*lame*sin_phi*sin_phi + 4.0*shear*(1.0 + sin_phi*sin_phi);
  denom_l = 4.0*lame*sin_phi*sin_phi +
            shear*(1.0 + sin_phi)*(1.0 + sin_phi) +
            2.0*shear*(1.0 - sin_phi)*(1.0 - sin_phi);
  denom_r = 4.0*lame*sin_phi*sin_phi +
            2.0*shear*(1.0 + sin_phi)*(1.0 + sin_phi) +
            shear*(1.0 - sin_phi)*(1.0 - sin_phi);
  denom_a = 4.0*bulk*sin_phi*sin_phi;
}



/**
  The function computes the returned principal stresses for the already known
  return branch.

  @param[in] eig - ordered trial eigenvalues
  @param[in] trace_E - trace of the ordered trial eigenvalues
  @param[in] return_type - selected return branch
  @param[in] lambda_s - smooth-face plastic multiplier
  @param[in] lambda_l - left-edge plastic multiplier
  @param[in] lambda_r - right-edge plastic multiplier
  @param[out] sigma - returned principal stresses

  @return The function does not return anything, the results are stored in the
          output arguments.
*/
void mohrc3d_ugn::compute_returned_principal_stresses(const vector &eig,
                                                      double trace_E,
                                                      int return_type,
                                                      double lambda_s,
                                                      double lambda_l,
                                                      double lambda_r,
                                                      vector &sigma) const
{
  // Recompute the returned principal stresses on the already selected branch.
  switch (return_type){
    case RET_ELASTIC:
      sigma[0] = lame*trace_E + 2.0*shear*eig[0];
      sigma[1] = lame*trace_E + 2.0*shear*eig[1];
      sigma[2] = lame*trace_E + 2.0*shear*eig[2];
      break;

    case RET_SMOOTH:
      sigma[0] = lame*trace_E + 2.0*shear*eig[0] -
        lambda_s*(2.0*lame*sin_phi + 2.0*shear*(1.0 + sin_phi));
      sigma[1] = lame*trace_E + 2.0*shear*eig[1] -
        lambda_s*(2.0*lame*sin_phi);
      sigma[2] = lame*trace_E + 2.0*shear*eig[2] -
        lambda_s*(2.0*lame*sin_phi - 2.0*shear*(1.0 - sin_phi));
      break;

    case RET_LEFT_EDGE:
      sigma[0] = lame*trace_E + shear*(eig[0] + eig[1]) -
        lambda_l*(2.0*lame*sin_phi + shear*(1.0 + sin_phi));
      sigma[1] = sigma[0];
      sigma[2] = lame*trace_E + 2.0*shear*eig[2] -
        lambda_l*(2.0*lame*sin_phi - 2.0*shear*(1.0 - sin_phi));
      break;

    case RET_RIGHT_EDGE:
      sigma[0] = lame*trace_E + 2.0*shear*eig[0] -
        lambda_r*(2.0*lame*sin_phi + 2.0*shear*(1.0 + sin_phi));
      sigma[1] = lame*trace_E + shear*(eig[1] + eig[2]) -
        lambda_r*(2.0*lame*sin_phi - shear*(1.0 - sin_phi));
      sigma[2] = sigma[1];
      break;

    case RET_APEX:
    default:
      sigma[0] = c_bar / (2.0 * sin_phi);
      sigma[1] = sigma[0];
      sigma[2] = sigma[0];
      break;
  }
}



/**
  The function reads material model parameters and stress return setup from the
  opened text.

  @param[in] in - pointer to the opened text file.

  @retval 0 - on success
  @retval 1 - in the case of an error
*/
long mohrc3d_ugn::read(XFILE *in)
{
  double young_read, poisson_read, cohesion_read, phi_read;
  double shear_read, bulk_read, lame_read, sin_phi_read, cos_phi_read, c_bar_read;

  if (in == NULL)
    return 1;

  if (xfscanf(in, "%lf %lf %lf %lf", &young_read, &poisson_read, &cohesion_read, &phi_read) != 4)
    return 1;

  // Precompute the elastic constants and Mohr-Coulomb parameters used in S and DS.
  shear_read = young_read / (2.0 * (1.0 + poisson_read));
  bulk_read = young_read / (3.0 * (1.0 - 2.0 * poisson_read));
  lame_read = bulk_read - 2.0 * shear_read / 3.0;
  sin_phi_read = sin(phi_read);
  cos_phi_read = cos(phi_read);
  c_bar_read = 2.0 * cohesion_read * cos_phi_read;

  if (young_read <= 0.0)
    return 1;
  if ((poisson_read <= -1.0) || (poisson_read >= 0.5))
    return 1;
  if (cohesion_read < 0.0)
    return 1;
  if ((phi_read <= 1.0e-14) || (phi_read >= 0.5 * M_PI - 1.0e-14))
    return 1;

  young = young_read;
  poisson = poisson_read;
  c = cohesion_read;
  phi = phi_read;
  shear = shear_read;
  bulk = bulk_read;
  lame = lame_read;
  sin_phi = sin_phi_read;
  cos_phi = cos_phi_read;
  c_bar = c_bar_read;
  return 0;
}



/**
  The function prints material model parameters and stress return setup to the
  opened text file.

  @param[in] out - pointer to the opened text output file.

  @return The function does not return anything.
*/
void mohrc3d_ugn::print(FILE *out)
{
  if (out == NULL)
    return;

  fprintf(out, "%le %le %le %le ", young, poisson, c, phi);

  /*
    fprintf(out, "MC perfect plastic associative 3D model\n");
    fprintf(out, "  E      = %.15g\n", young);
    fprintf(out, "  nu     = %.15g\n", poisson);
    fprintf(out, "  c      = %.15g\n", c);
    fprintf(out, "  phi    = %.15g rad\n", phi);
    fprintf(out, "  G      = %.15g\n", shear);
    fprintf(out, "  K      = %.15g\n", bulk);
    fprintf(out, "  lambda = %.15g\n", lame);
    fprintf(out, "  other  = %d components\n", NCOMP_OTHER);
    fprintf(out, "  eqother= %d components\n", NCOMP_EQOTHER);
    fprintf(out, "  other  = epsp(6), epsp_prev(6), return_type(1)\n");
    fprintf(out, "  eqother= epsp(6)\n");
  */
}



/**
  Computes the actual stresses based on the attained strains and state
  variables.

  @param[in] ipp Integration point number in the mechmat integration-point
                 array.
  @param[in] ido Index shift for the given material in the ipp other array.
                 It determines where the state variables of the given
                 material model begin.

  @return This function does not return a value. The computed results are
          stored in the stress and other arrays of the given integration
          point.
*/
void mohrc3d_ugn::nlstresses(long ipp, long ido)
{
  long ncompstr = Mm->ip[ipp].ncompstr; // the number of stress/strain components at the given integration point
  vector eps(ASTCKVEC(ncompstr));
  vector sig(ASTCKVEC(ncompstr));
  vector other;
  vector eqother(ASTCKVEC(NCOMP_EQOTHER));
  const long lcid = 0;
  Mm->givestrain(lcid, ipp, eps); // take actual strains from the ipp-th integration point
  Mm->givestress(lcid, ipp, sig); // take actual stresses from the ipp-th integration point
  other.makerefv(NCOMP_OTHER, Mm->ip[ipp].other+ido); // make reference vector to the dedicated part of the state variable array at the ipp-th point
  eqother.makerefv(NCOMP_EQOTHER, Mm->ip[ipp].eqother+ido); // make reference vector to the dedicated part of the equilibrium state variable array at the ipp-th point

  nlstresses(eps, eqother, sig, other);

  Mm->storestress(lcid, ipp, sig);
}



/**
  Computes the actual stresses based on the attained nonlocal strains and state
  variables.

  @param[in] ipp Integration point number in the mechmat integration-point array.
  @param[in] ido Index shift for the given material in the ipp other array.
                 It determines where the state variables of the given material
                 model begin.

  @return This function does not return a value. The computed results are stored
          in the stress and other arrays of the given integration point.
*/
void mohrc3d_ugn::nonloc_nlstresses (long ipp, long ido)
{
  long ncompstr = Mm->ip[ipp].ncompstr; // the number of stress/strain components at the given integration point
  vector eps(ASTCKVEC(ncompstr));
  vector sig(ASTCKVEC(ncompstr));
  vector other;
  vector eqother;
  const long lcid = 0;
  Mm->givestrain(lcid, ipp, eps); // take actual strains from the ipp-th integration point
  Mm->givestress(lcid, ipp, sig); // take actual stresses from the ipp-th integration point
  other.makerefv(NCOMP_OTHER, Mm->ip[ipp].other+ido);  // make reference vector to the dedicated part of the state variable array at the ipp-th point
  eqother.makerefv(NCOMP_EQOTHER, Mm->ip[ipp].nonloc); // make reference vector to the averaged plastic strain array at the ipp-th point

  nlstresses(eps, eqother, sig, other);

  Mm->storestress(lcid, ipp, sig);
}



/**
  The function computes actual stresses with respect to the attained strains and
  state variables.

  @param[in] strain - array of actual strain, components are ordered as follows
                      eps_x, eps_y, eps_z, gamma_yz, gamma_xz, gamma_xy - for the
                      space stress problem

  @param[in] eqstatev - point buffer from the last converged state; nlstresses()
                        reads the plastic-strain history from its first 6 entries

  @param[out] stress - array of the resulting stress components, it must be
                       computed, components are ordered as follows
                       sig_x, sig_y, sig_z, tau_yz, tau_xz, tau_xy - for the
                       space stress problem

  @param[out] statev - array of the resulting current-point buffer for the actual
                       strains

  @return The function does not return anything, the results are stored in the
          arrays stress and statev passed in as arguments.
*/
void mohrc3d_ugn::nlstresses(const vector &strain,
                             const vector &eqstatev,
                             vector &stress,
                             vector &statev)
{
  /*
    Scalar material-point version of Matlab:

    [S,eig_1,eig_2,eig_3,Eig_1,Eig_2,Eig_3,EIG_1,EIG_2,EIG_3,...
    sigma_1,sigma_2,sigma_3,return_type] = constitutive_problem(E_new, Ep_prev)

    plus the plastic-strain recovery used later in matlab_export_3d/newton.m:

    Ep = -Inv_ELAST*S;
    Ep(1:6,:) = Ep(1:6,:) + E;
  */
  vector E(ASTCKVEC(NCOMP_STRAIN)), Ep_prev(ASTCKVEC(NCOMP_EQOTHER)), E_trial(ASTCKVEC(NCOMP_STRAIN));
  vector Ep(ASTCKVEC(NCOMP_EQOTHER)), eps_el(ASTCKVEC(NCOMP_STRAIN)), eig(ASTCKVEC(3)), sigma(ASTCKVEC(3)), S(ASTCKVEC(NCOMP_STRESS));
  vector E_tr(ASTCKVEC(NCOMP_STRAIN)), E_square(ASTCKVEC(NCOMP_STRAIN)), iota(ASTCKVEC(NCOMP_STRAIN));
  vector Eig_1(ASTCKVEC(NCOMP_STRAIN)), Eig_2(ASTCKVEC(NCOMP_STRAIN)), Eig_3(ASTCKVEC(NCOMP_STRAIN));
  vector Eig_12(ASTCKVEC(NCOMP_STRAIN)), Eig_23(ASTCKVEC(NCOMP_STRAIN));
  matrix DER_E_square(ASTCKMAT(NCOMP_STRAIN, NCOMP_STRAIN));
  double trace_E, f_tr, gamma_sl, gamma_sr, gamma_la, gamma_ra;
  double denom_s, denom_l, denom_r, denom_a;
  double lambda_s, lambda_l, lambda_r;
  int return_type = RET_ELASTIC;

  // Resize the point outputs before filling the current material response.
  reallocv(RSTCKVEC(NCOMP_STRESS, stress));
  reallocv(RSTCKVEC(NCOMP_OTHER, statev));
  nullv(E);
  nullv(Ep_prev);
  nullv(E_trial);
  nullv(Ep);
  nullv(eps_el);
  nullv(eig);
  nullv(sigma);
  nullv(S);
  nullv(E_tr);
  nullv(E_square);
  nullv(iota);
  nullv(Eig_1);
  nullv(Eig_2);
  nullv(Eig_3);
  nullv(Eig_12);
  nullv(Eig_23);
  nullm(DER_E_square);
  nullv(stress);
  nullv(statev);

  // The space-stress unit tensor enters the volumetric part of every branch.
  iota[0] = 1.0;
  iota[1] = 1.0;
  iota[2] = 1.0;

  // Read the actual strain and the converged plastic-strain history.
  copyv(strain, E);
  copyv(&eqstatev[O_EP_XX], Ep_prev);

  // Matlab forms the trial strain by subtracting the old plastic strain.
  addmultv(E, 1.0, Ep_prev, -1.0, E_trial);

  // Rebuild the ordered trial eigenvalues and the arrays used by the projection formulas.
  compute_ordered_trial_spectral(E_trial, E_tr, E_square, DER_E_square, eig);

  // Evaluate the Mohr-Coulomb trial function and the return candidates.
  trace_E = eig[0] + eig[1] + eig[2];
  f_tr = 2.0 * shear * ((1.0 + sin_phi) * eig[0] - (1.0 - sin_phi) * eig[2]) +
    2.0 * lame * sin_phi * trace_E - c_bar;
  gamma_sl = (eig[0] - eig[1]) / (1.0 + sin_phi);
  gamma_sr = (eig[1] - eig[2]) / (1.0 - sin_phi);
  gamma_la = (eig[0] + eig[1] - 2.0 * eig[2]) / (3.0 - sin_phi);
  gamma_ra = (2.0 * eig[0] - eig[1] - eig[2]) / (3.0 + sin_phi);

  // These denominators are the same branch constants as in the Matlab code.
  compute_return_denominators(denom_s, denom_l, denom_r, denom_a);
  (void)denom_a;

  // Compute the plastic multipliers for the smooth and edge return candidates.
  lambda_s = f_tr / denom_s;
  lambda_l = (shear * ((1.0 + sin_phi) * (eig[0] + eig[1]) - 2.0 * (1.0 - sin_phi) * eig[2]) +
              2.0 * lame * sin_phi * trace_E - c_bar) / denom_l;
  lambda_r = (shear * (2.0 * (1.0 + sin_phi) * eig[0] - (1.0 - sin_phi) * (eig[1] + eig[2])) +
              2.0 * lame * sin_phi * trace_E - c_bar) / denom_r;

  // Choose the return branch from the Matlab admissibility tests.
  if (f_tr <= 0.0)
    return_type = RET_ELASTIC;
  else if (lambda_s <= ((gamma_sl < gamma_sr) ? gamma_sl : gamma_sr))
    return_type = RET_SMOOTH;
  else if ((gamma_sl < gamma_sr) && (lambda_l >= gamma_sl) && (lambda_l <= gamma_la))
    return_type = RET_LEFT_EDGE;
  else if ((gamma_sl > gamma_sr) && (lambda_r >= gamma_sr) && (lambda_r <= gamma_ra))
    return_type = RET_RIGHT_EDGE;
  else
    return_type = RET_APEX;

  // Compute the returned principal stresses for the selected branch.
  compute_returned_principal_stresses(eig, trace_E, return_type, lambda_s, lambda_l, lambda_r, sigma);

  // Rebuild the full stress vector from the returned principal stresses.
  nullv(S);
  switch (return_type){
    case RET_ELASTIC:
      // Matlab: S(:,test_el)=ELAST*E_trial(:,test_el)
      for (long i=0; i<NCOMP_STRESS; i++)
        S[i] = lame * trace_E * iota[i] + 2.0 * shear * E_tr[i];
      break;

    case RET_SMOOTH:
      // Matlab smooth branch:
      // S(:,test_s)=sigma_1*Eig_1+sigma_2*Eig_2+sigma_3*Eig_3
      addmultv(E_square, 1.0, E_tr, -(eig[1] + eig[2]), Eig_1);
      addmultv(Eig_1, iota, eig[1] * eig[2]);
      cmulv(1.0 / ((eig[0] - eig[1]) * (eig[0] - eig[2])), Eig_1);
      addmultv(E_square, 1.0, E_tr, -(eig[0] + eig[2]), Eig_2);
      addmultv(Eig_2, iota, eig[0] * eig[2]);
      cmulv(1.0 / ((eig[1] - eig[0]) * (eig[1] - eig[2])), Eig_2);
      addmultv(E_square, 1.0, E_tr, -(eig[0] + eig[1]), Eig_3);
      addmultv(Eig_3, iota, eig[0] * eig[1]);
      cmulv(1.0 / ((eig[2] - eig[0]) * (eig[2] - eig[1])), Eig_3);
      addmultv(S, Eig_1, sigma[0]);
      addmultv(S, Eig_2, sigma[1]);
      addmultv(S, Eig_3, sigma[2]);
      break;

    case RET_LEFT_EDGE:
      // Matlab left-edge branch:
      // S(:,test_l)=sigma_1*Eig_12+sigma_3*Eig_3
      addmultv(E_square, 1.0, E_tr, -(eig[0] + eig[1]), Eig_3);
      addmultv(Eig_3, iota, eig[0] * eig[1]);
      cmulv(1.0 / ((eig[2] - eig[0]) * (eig[2] - eig[1])), Eig_3);
      addmultv(iota, 1.0, Eig_3, -1.0, Eig_12);
      addmultv(S, Eig_12, sigma[0]);
      addmultv(S, Eig_3, sigma[2]);
      break;

    case RET_RIGHT_EDGE:
      // Matlab right-edge branch:
      // S(:,test_r)=sigma_1*Eig_1+sigma_3*Eig_23
      addmultv(E_square, 1.0, E_tr, -(eig[1] + eig[2]), Eig_1);
      addmultv(Eig_1, iota, eig[1] * eig[2]);
      cmulv(1.0 / ((eig[0] - eig[1]) * (eig[0] - eig[2])), Eig_1);
      addmultv(iota, 1.0, Eig_1, -1.0, Eig_23);
      addmultv(S, Eig_1, sigma[0]);
      addmultv(S, Eig_23, sigma[2]);
      break;

    case RET_APEX:
    default:
      // Matlab: S(:,test_a)=iota*sigma_1
      addmultv(S, iota, sigma[0]);
      break;
  }

  // Convert stress back to elastic strain and recover the updated plastic strain.
  eps_el[0] = (S[0] - poisson * (S[1] + S[2])) / young;
  eps_el[1] = (S[1] - poisson * (S[0] + S[2])) / young;
  eps_el[2] = (S[2] - poisson * (S[0] + S[1])) / young;
  eps_el[3] = S[3] / shear;
  eps_el[4] = S[4] / shear;
  eps_el[5] = S[5] / shear;
  addmultv(E, 1.0, eps_el, -1.0, Ep);

  // Pack only the persistent history and the minimal tangent helper payload.
  copyv(Ep, &statev[O_EP_XX]);
  copyv(Ep_prev, &statev[O_EP_PREV_XX]);
  statev[O_RET] = static_cast<double>(return_type);

  // Return the stress vector to the caller.
  copyv(S, stress);
}



/**
  Computes the material stiffness %matrix at the given integration point.

  @param[out] d   Resulting material stiffness %matrix. The %matrix must be
                  allocated with dimensions ncompstr x ncompstr.
  @param[in]  ipp Integration point number in the mechmat integration-point
                  array.
  @param[in]  ido Index shift for the given material in the ipp other array.
                  It determines where the state variables of the given
                  material model begin.

  @return The function returns the material stiffness %matrix in the argument d.
*/
void mohrc3d_ugn::matstiff(matrix &d, long ipp, long ido)
{
  long ncompstr = Mm->ip[ipp].ncompstr; // the number of stress/strain components at the given integration point
  vector eps(ASTCKVEC(ncompstr));
  vector sig(ASTCKVEC(ncompstr));
  vector other;
  const long lcid = 0;

  Mm->givestrain(lcid, ipp, eps); // take actual strains from the ipp-th integration point
  Mm->givestress(lcid, ipp, sig); // take actual stresses from the ipp-th integration point
  other.makerefv(NCOMP_OTHER, Mm->ip[ipp].other+ido); // make reference vector to the dedicated part of the state variable array at the ipp-th point

  stiffmat(eps, other, sig, d);
}



/**
  The function computes material stiffness matrix with respect to the attained
  strains and state variables.

  @param[in] strain - array of actual strain, components are ordered as follows
                      eps_x, eps_y, eps_z, gamma_yz, gamma_xz, gamma_xy - for the
                      space stress problem

  @param[in] statev - current-point buffer produced by nlstresses(); this buffer
                      stores the updated plastic strain, the previous plastic
                      strain, and the selected return type

  @param[in] stress - array of the resulting stress components, it must be
                      computed, components are ordered as follows
                      sig_x, sig_y, sig_z, tau_yz, tau_xz, tau_xy - for the
                      space stress problem

  @param[out] d - the resulting material stiffness matrix, it must be calculated.

  @return The function does not return anything, the resulting matrix is stored
          in the matrix d passed in as an argument.
*/
void mohrc3d_ugn::stiffmat(const vector &strain,
                           const vector &statev,
                           const vector &stress,
                           matrix &d)
{
  /*
    Scalar material-point version of Matlab stiffness_matrix.m with minimal
    recomputation. The tangent call reuses the stored return type and rebuilds
    the trial spectral data from strain and epsp_prev.
  */
  vector Ep_prev(ASTCKVEC(NCOMP_EQOTHER)), E_trial(ASTCKVEC(NCOMP_STRAIN)), eig(ASTCKVEC(3)), sigma(ASTCKVEC(3));
  vector E_tr(ASTCKVEC(NCOMP_STRAIN)), E_square(ASTCKVEC(NCOMP_STRAIN)), iota(ASTCKVEC(NCOMP_STRAIN));
  vector Eig_1(ASTCKVEC(NCOMP_STRAIN)), Eig_2(ASTCKVEC(NCOMP_STRAIN)), Eig_3(ASTCKVEC(NCOMP_STRAIN));
  vector Eig_12(ASTCKVEC(NCOMP_STRAIN)), Eig_23(ASTCKVEC(NCOMP_STRAIN)), Eig_6(ASTCKVEC(NCOMP_STRAIN));
  matrix DER_E_square(ASTCKMAT(NCOMP_STRAIN, NCOMP_STRAIN));
  matrix EIG_1(ASTCKMAT(NCOMP_STRAIN, NCOMP_STRAIN));
  matrix EIG_2(ASTCKMAT(NCOMP_STRAIN, NCOMP_STRAIN));
  matrix EIG_3(ASTCKMAT(NCOMP_STRAIN, NCOMP_STRAIN));
  matrix IDENT(ASTCKMAT(NCOMP_STRAIN, NCOMP_STRAIN));
  matrix zero(ASTCKMAT(NCOMP_STRAIN, NCOMP_STRAIN));
  matrix outer(ASTCKMAT(NCOMP_STRAIN, NCOMP_STRAIN));
  matrix outer_1(ASTCKMAT(NCOMP_STRAIN, NCOMP_STRAIN));
  matrix outer_2(ASTCKMAT(NCOMP_STRAIN, NCOMP_STRAIN));
  matrix outer_3(ASTCKMAT(NCOMP_STRAIN, NCOMP_STRAIN));
  matrix outer_4(ASTCKMAT(NCOMP_STRAIN, NCOMP_STRAIN));
  matrix outer_5(ASTCKMAT(NCOMP_STRAIN, NCOMP_STRAIN));
  matrix outer_6(ASTCKMAT(NCOMP_STRAIN, NCOMP_STRAIN));
  double trace_E, f_tr, lambda_s, lambda_l, lambda_r;
  double denom_s, denom_l, denom_r, denom_a;
  double denom_1, denom_2, denom_3;
  int return_type = RET_ELASTIC;
  (void)stress; // ???!!!

  // Allocate and clear the tangent and its local work arrays.
  reallocm(NCOMP_STRAIN, NCOMP_STRAIN, d);
  nullm(d);
  nullv(Ep_prev);
  nullv(E_trial);
  nullv(eig);
  nullv(sigma);
  nullv(E_tr);
  nullv(E_square);
  nullv(iota);
  nullv(Eig_1);
  nullv(Eig_2);
  nullv(Eig_3);
  nullv(Eig_12);
  nullv(Eig_23);
  nullv(Eig_6);
  nullm(DER_E_square);
  nullm(EIG_1);
  nullm(EIG_2);
  nullm(EIG_3);
  nullm(IDENT);
  nullm(zero);
  nullm(outer);
  nullm(outer_1);
  nullm(outer_2);
  nullm(outer_3);
  nullm(outer_4);
  nullm(outer_5);
  nullm(outer_6);

  // Read the stored return branch from the current point buffer.
  return_type = static_cast<int>(statev[O_RET] + 0.5);

  // The space-stress unit tensor enters the volumetric part of every tangent branch.
  iota[0] = 1.0;
  iota[1] = 1.0;
  iota[2] = 1.0;

  // This is the 3D identity metric used in the Matlab projection Hessian formulas.
  IDENT(0,0) = 1.0;
  IDENT(1,1) = 1.0;
  IDENT(2,2) = 1.0;
  IDENT(3,3) = 0.5;
  IDENT(4,4) = 0.5;
  IDENT(5,5) = 0.5;

  // Rebuild the trial strain from the actual strain and the previous plastic strain.
  copyv(&statev[O_EP_PREV_XX], Ep_prev);
  addmultv(strain, 1.0, Ep_prev, -1.0, E_trial);

  // Recompute the ordered trial spectral data needed by stiffness_matrix.m.
  compute_ordered_trial_spectral(E_trial, E_tr, E_square, DER_E_square, eig);

  // Rebuild the branch denominators and the branch-specific plastic multipliers.
  trace_E = eig[0] + eig[1] + eig[2];
  compute_return_denominators(denom_s, denom_l, denom_r, denom_a);
  (void)denom_a;
  f_tr = 2.0 * shear * ((1.0 + sin_phi) * eig[0] - (1.0 - sin_phi) * eig[2]) +
    2.0 * lame * sin_phi * trace_E - c_bar;
  lambda_s = f_tr / denom_s;
  lambda_l = (shear * ((1.0 + sin_phi) * (eig[0] + eig[1]) - 2.0 * (1.0 - sin_phi) * eig[2]) +
              2.0 * lame * sin_phi * trace_E - c_bar) / denom_l;
  lambda_r = (shear * (2.0 * (1.0 + sin_phi) * eig[0] - (1.0 - sin_phi) * (eig[1] + eig[2])) +
              2.0 * lame * sin_phi * trace_E - c_bar) / denom_r;

  // Recompute the returned principal stresses on the already selected branch.
  compute_returned_principal_stresses(eig, trace_E, return_type, lambda_s, lambda_l, lambda_r, sigma);

  // Assemble the 6x6 SIFEL tangent from the rebuilt Matlab spectral data.
  switch (return_type){
    case RET_ELASTIC:
      // Matlab: Sderiv(:,test_el)=Elast(:)
      d(0, 0) = lame + 2.0 * shear;
      d(0, 1) = lame;
      d(0, 2) = lame;
      d(1, 0) = lame;
      d(1, 1) = lame + 2.0 * shear;
      d(1, 2) = lame;
      d(2, 0) = lame;
      d(2, 1) = lame;
      d(2, 2) = lame + 2.0 * shear;
      d(3, 3) = shear;
      d(4, 4) = shear;
      d(5, 5) = shear;
      break;

    case RET_SMOOTH:
      // Matlab smooth branch:
      // Sderiv(:,test_s)=mat1_s+mat2_s+mat3_s+mat4_s+mat5_s-mat6_s
      denom_1 = (eig[0] - eig[1]) * (eig[0] - eig[2]);
      denom_2 = (eig[1] - eig[0]) * (eig[1] - eig[2]);
      denom_3 = (eig[2] - eig[0]) * (eig[2] - eig[1]);
      addmultv(E_square, 1.0, E_tr, -(eig[1] + eig[2]), Eig_1);
      addmultv(Eig_1, iota, eig[1] * eig[2]);
      cmulv(1.0 / denom_1, Eig_1);
      addmultv(E_square, 1.0, E_tr, -(eig[0] + eig[2]), Eig_2);
      addmultv(Eig_2, iota, eig[0] * eig[2]);
      cmulv(1.0 / denom_2, Eig_2);
      addmultv(E_square, 1.0, E_tr, -(eig[0] + eig[1]), Eig_3);
      addmultv(Eig_3, iota, eig[0] * eig[1]);
      cmulv(1.0 / denom_3, Eig_3);

      vxv(Eig_1, Eig_1, outer_1);
      vxv(Eig_2, Eig_2, outer_2);
      vxv(Eig_3, Eig_3, outer_3);
      copym(DER_E_square, EIG_1);
      addmultm(zero, 0.0, IDENT, -(eig[1] + eig[2]), EIG_1);
      addmultm(zero, 0.0, outer_1, -(2.0 * eig[0] - eig[1] - eig[2]), EIG_1);
      addmultm(zero, 0.0, outer_2, -(eig[1] - eig[2]), EIG_1);
      addmultm(zero, 0.0, outer_3, eig[1] - eig[2], EIG_1);
      cmulm(1.0 / denom_1, EIG_1);
      copym(DER_E_square, EIG_2);
      addmultm(zero, 0.0, IDENT, -(eig[0] + eig[2]), EIG_2);
      addmultm(zero, 0.0, outer_2, -(2.0 * eig[1] - eig[0] - eig[2]), EIG_2);
      addmultm(zero, 0.0, outer_1, -(eig[0] - eig[2]), EIG_2);
      addmultm(zero, 0.0, outer_3, eig[0] - eig[2], EIG_2);
      cmulm(1.0 / denom_2, EIG_2);
      copym(DER_E_square, EIG_3);
      addmultm(zero, 0.0, IDENT, -(eig[0] + eig[1]), EIG_3);
      addmultm(zero, 0.0, outer_3, -(2.0 * eig[2] - eig[0] - eig[1]), EIG_3);
      addmultm(zero, 0.0, outer_1, -(eig[0] - eig[1]), EIG_3);
      addmultm(zero, 0.0, outer_2, eig[0] - eig[1], EIG_3);
      cmulm(1.0 / denom_3, EIG_3);

      addmultm(zero, 0.0, EIG_1, sigma[0], d);
      addmultm(zero, 0.0, EIG_2, sigma[1], d);
      addmultm(zero, 0.0, EIG_3, sigma[2], d);
      vxv(iota, iota, outer);
      addmultm(zero, 0.0, outer, lame, d);
      vxv(Eig_1, Eig_1, outer);
      addmultm(zero, 0.0, outer, 2.0 * shear, d);
      vxv(Eig_2, Eig_2, outer);
      addmultm(zero, 0.0, outer, 2.0 * shear, d);
      vxv(Eig_3, Eig_3, outer);
      addmultm(zero, 0.0, outer, 2.0 * shear, d);

      // Matlab: Eig_6=2*shear*((1+sin_phi)*Eig_1-(1-sin_phi)*Eig_3)+2*lame*sin_phi*[1;1;1;0;0;0]
      addmultv(Eig_1, 2.0 * shear * (1.0 + sin_phi), Eig_3, -2.0 * shear * (1.0 - sin_phi), Eig_6);
      addmultv(Eig_6, iota, 2.0 * lame * sin_phi);
      vxv(Eig_6, Eig_6, outer);
      addmultm(zero, 0.0, outer, -1.0 / denom_s, d);
      break;

    case RET_LEFT_EDGE:
      // Matlab left-edge branch:
      // Sderiv(:,test_l)=mat1_l+mat2_l+mat3_l+mat5_l-mat6_l
      denom_3 = (eig[2] - eig[0]) * (eig[2] - eig[1]);
      addmultv(E_square, 1.0, E_tr, -(eig[0] + eig[1]), Eig_3);
      addmultv(Eig_3, iota, eig[0] * eig[1]);
      cmulv(1.0 / denom_3, Eig_3);
      addmultv(iota, 1.0, Eig_3, -1.0, Eig_12);

      vxv(Eig_12, Eig_12, outer_1);
      vxv(Eig_3, Eig_3, outer_2);
      vxv(Eig_12, E_tr, outer_3);
      vxv(E_tr, Eig_12, outer_4);
      vxv(Eig_12, Eig_3, outer_5);
      vxv(Eig_3, Eig_12, outer_6);
      copym(DER_E_square, EIG_3);
      addmultm(zero, 0.0, IDENT, -(eig[0] + eig[1]), EIG_3);
      addmultm(zero, 0.0, outer_4, -1.0, EIG_3);
      addmultm(zero, 0.0, outer_3, -1.0, EIG_3);
      addmultm(zero, 0.0, outer_1, eig[0] + eig[1], EIG_3);
      addmultm(zero, 0.0, outer_2, eig[0] + eig[1] - 2.0 * eig[2], EIG_3);
      addmultm(zero, 0.0, outer_5, eig[2], EIG_3);
      addmultm(zero, 0.0, outer_6, eig[2], EIG_3);
      cmulm(1.0 / denom_3, EIG_3);

      addmultm(zero, 0.0, EIG_3, sigma[2] - sigma[0], d);
      vxv(iota, iota, outer);
      addmultm(zero, 0.0, outer, lame, d);
      vxv(Eig_12, Eig_12, outer);
      addmultm(zero, 0.0, outer, shear, d);
      vxv(Eig_3, Eig_3, outer);
      addmultm(zero, 0.0, outer, 2.0 * shear, d);

      // Matlab: Eig_12=iota-Eig_3 and Eig_6=shear*((1+sin_phi)*Eig_12-2*(1-sin_phi)*Eig_3)+...
      addmultv(Eig_12, shear * (1.0 + sin_phi), Eig_3, -2.0 * shear * (1.0 - sin_phi), Eig_6);
      addmultv(Eig_6, iota, 2.0 * lame * sin_phi);
      vxv(Eig_6, Eig_6, outer);
      addmultm(zero, 0.0, outer, -1.0 / denom_l, d);
      break;

    case RET_RIGHT_EDGE:
      // Matlab right-edge branch:
      // Sderiv(:,test_r)=mat1_r+mat2_r+mat3_r+mat5_r-mat6_r
      denom_1 = (eig[0] - eig[1]) * (eig[0] - eig[2]);
      addmultv(E_square, 1.0, E_tr, -(eig[1] + eig[2]), Eig_1);
      addmultv(Eig_1, iota, eig[1] * eig[2]);
      cmulv(1.0 / denom_1, Eig_1);
      addmultv(iota, 1.0, Eig_1, -1.0, Eig_23);

      vxv(Eig_1, Eig_1, outer_1);
      vxv(Eig_23, Eig_23, outer_2);
      vxv(Eig_23, E_tr, outer_3);
      vxv(E_tr, Eig_23, outer_4);
      vxv(Eig_23, Eig_1, outer_5);
      vxv(Eig_1, Eig_23, outer_6);
      copym(DER_E_square, EIG_1);
      addmultm(zero, 0.0, IDENT, -(eig[1] + eig[2]), EIG_1);
      addmultm(zero, 0.0, outer_4, -1.0, EIG_1);
      addmultm(zero, 0.0, outer_3, -1.0, EIG_1);
      addmultm(zero, 0.0, outer_2, eig[1] + eig[2], EIG_1);
      addmultm(zero, 0.0, outer_1, eig[1] + eig[2] - 2.0 * eig[0], EIG_1);
      addmultm(zero, 0.0, outer_5, eig[0], EIG_1);
      addmultm(zero, 0.0, outer_6, eig[0], EIG_1);
      cmulm(1.0 / denom_1, EIG_1);

      addmultm(zero, 0.0, EIG_1, sigma[0] - sigma[2], d);
      vxv(iota, iota, outer);
      addmultm(zero, 0.0, outer, lame, d);
      vxv(Eig_1, Eig_1, outer);
      addmultm(zero, 0.0, outer, 2.0 * shear, d);
      vxv(Eig_23, Eig_23, outer);
      addmultm(zero, 0.0, outer, shear, d);

      // Matlab: Eig_23=iota-Eig_1 and Eig_6=shear*(2*(1+sin_phi)*Eig_1-(1-sin_phi)*Eig_23)+...
      addmultv(Eig_1, 2.0 * shear * (1.0 + sin_phi), Eig_23, -shear * (1.0 - sin_phi), Eig_6);
      addmultv(Eig_6, iota, 2.0 * lame * sin_phi);
      vxv(Eig_6, Eig_6, outer);
      addmultm(zero, 0.0, outer, -1.0 / denom_r, d);
      break;

    case RET_APEX:
    default:
      // Matlab: Sderiv(:,test_a)=zeros(36,nt_a)
      break;
  }
}



/**
  The function updates the state variable array eqother from the actual
  converged iteration values.

  @param[in]  ipp - integration point identfier,
  @param[in]  im  - index of the given material in the material model
                    description array at the given integration point,
  @param[in]  ido - index shift of the state variable array other.
  @return The function does not return anything, the resulting converged state
          of the given integration point other array is stored in the array
          eqother.
*/
void mohrc3d_ugn::updateval(long ipp, long im, long ido)
{
  vector statev;
  vector eqstatev;

  makerefv(statev, Mm->ip[ipp].other+ido, NCOMP_OTHER);

  long idoeq = Mm->givencompeqother(ipp, 0) - Mm->givencompeqother(ipp, im);
  makerefv(eqstatev, Mm->ip[ipp].eqother+idoeq, NCOMP_EQOTHER);

  updateval(statev, eqstatev);
}



/**
  Function returns irreversible plastic strains.

  @param ipp   - integration point number in the mechmat ip array.
  @param ido   - index of the first internal variable for given material in the
                 ipp other array
  @param epsp  - %vector of irreversible strains

  @return The function returns %vector of irreversible strains in the parameter
          epsp.
*/
void mohrc3d_ugn::giveirrstrains (long ipp, long ido, vector &epsp)
{
  long i;
  for (i=0;i<epsp.n;i++)
    epsp[i] = Mm->ip[ipp].other[ido+i];
}



/**
  The function updates the point buffer from the actual converged iteration
  values.

  @param[in] statev - current-point buffer produced by nlstresses()

  @param[out] eqstatev - point buffer for the last converged state

  @return The function does not return anything, the resulting converged state
          is stored in the array eqstatev passed in as an argument.
*/
void mohrc3d_ugn::updateval(const vector &statev, vector &eqstatev)
{
  // Persist only the updated plastic strain between converged steps.
  reallocv(NCOMP_EQOTHER, eqstatev);
  copyv(&statev[O_EP_XX], eqstatev);
}
