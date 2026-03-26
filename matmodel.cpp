#include "matmodel.h"
#include "vector.h"
#include "matrix.h"

#include <math.h>
#include <stdio.h>

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
  The function computes the ordered trial eigenvalues together with their first and
  second derivatives with respect to the plane-strain Voigt components.

  @param[in] E_trial - trial strain vector
  @param[out] eig - ordered trial eigenvalues
  @param[out] Eig_1 - first derivative of the first ordered eigenvalue
  @param[out] Eig_2 - first derivative of the second ordered eigenvalue
  @param[out] Eig_3 - first derivative of the third ordered eigenvalue
  @param[out] EIG_1 - second derivative of the first ordered eigenvalue
  @param[out] EIG_2 - second derivative of the second ordered eigenvalue
  @param[out] EIG_3 - second derivative of the third ordered eigenvalue

  @return The function does not return anything, the results are stored in the output
          arguments.
*/
void matmodel::compute_ordered_trial_spectral(const vector &E_trial,
                                              vector &eig,
                                              vector &Eig_1,
                                              vector &Eig_2,
                                              vector &Eig_3,
                                              matrix &EIG_1,
                                              matrix &EIG_2,
                                              matrix &EIG_3) const
{
    vector Eig0_1(NCOMP_STRAIN), Eig0_2(NCOMP_STRAIN), Eig0_3(NCOMP_STRAIN);
    matrix metric(NCOMP_STRAIN, NCOMP_STRAIN);
    matrix EIG0_1(NCOMP_STRAIN, NCOMP_STRAIN);
    matrix EIG0_2(NCOMP_STRAIN, NCOMP_STRAIN);
    matrix EIG0_3(NCOMP_STRAIN, NCOMP_STRAIN);
    double eig0_1, eig0_2, eig0_3;

    // Clear the spectral work arrays before filling the current trial state.
    nullv(Eig0_1);
    nullv(Eig0_2);
    nullv(Eig0_3);
    nullm(metric);
    nullm(EIG0_1);
    nullm(EIG0_2);
    nullm(EIG0_3);
    nullv(eig);
    nullv(Eig_1);
    nullv(Eig_2);
    nullv(Eig_3);
    nullm(EIG_1);
    nullm(EIG_2);
    nullm(EIG_3);

    // This is the reduced identity metric used in the Matlab Hessian formula.
    metric(0,0) = 1.0;
    metric(1,1) = 1.0;
    metric(2,2) = 0.5;

    {
        // Split the in-plane trial strain into invariants used by the spectral formulas.
        const double ex = E_trial[0];
        const double ey = E_trial[1];
        const double g = E_trial[2];
        const double d = ex - ey;
        const double i1 = ex + ey;
        const double i2 = sqrt(d*d + g*g);
        const double tol = 1.0e-14;

        // Compute the unsorted trial eigenvalues exactly as in constitutive_problem.m.
        eig0_1 = 0.5*(i1 + i2);
        eig0_2 = 0.5*(i1 - i2);
        eig0_3 = E_trial[3];

        // Matlab uses a special fallback when the in-plane eigenvalues coincide.
        if (i2 <= tol)
        {
            Eig0_1[0] = 1.0;
            Eig0_1[1] = 1.0;
        }
        else
        {
            // Build the first derivatives of the in-plane trial eigenvalues.
            Eig0_1[0] = (ex - eig0_2)/i2;
            Eig0_1[1] = (ey - eig0_2)/i2;
            Eig0_1[2] = g/(2.0*i2);

            Eig0_2[0] = 1.0 - Eig0_1[0];
            Eig0_2[1] = 1.0 - Eig0_1[1];
            Eig0_2[2] = -Eig0_1[2];

            // Assemble the Hessians of the in-plane trial eigenvalues.
            for (long a=0; a<NCOMP_STRAIN; a++)
            {
                for (long b=0; b<NCOMP_STRAIN; b++)
                {
                    EIG0_1(a,b) = (metric(a,b) - Eig0_1[a]*Eig0_1[b] - Eig0_2[a]*Eig0_2[b])/i2;
                    EIG0_2(a,b) = -EIG0_1(a,b);
                }
            }

            // The zz direction is decoupled from the in-plane spectral problem.
            for (long a=0; a<NCOMP_STRAIN; a++)
            {
                EIG0_1(3,a) = 0.0;
                EIG0_1(a,3) = 0.0;
                EIG0_2(3,a) = 0.0;
                EIG0_2(a,3) = 0.0;
            }
        }

        // The third principal mode is the zz direction.
        Eig0_3[3] = 1.0;
    }

    // Start from the natural order before applying the Matlab reordering rules.
    eig[0] = eig0_1;
    eig[1] = eig0_2;
    eig[2] = eig0_3;
    copyv(Eig0_1, Eig_1);
    copyv(Eig0_2, Eig_2);
    copyv(Eig0_3, Eig_3);
    copym(EIG0_1, EIG_1);
    copym(EIG0_2, EIG_2);
    copym(EIG0_3, EIG_3);

    // Swap the zz mode into the Matlab descending order when needed.
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
}

/**
  The function computes the three return denominators used by the smooth and edge
  branches of the Mohr-Coulomb return mapping.

  @param[out] denom_s - denominator for the smooth-face return
  @param[out] denom_l - denominator for the left-edge return
  @param[out] denom_r - denominator for the right-edge return

  @return The function does not return anything, the results are stored in the output
          arguments.
*/
void matmodel::compute_return_denominators(double &denom_s,
                                           double &denom_l,
                                           double &denom_r) const
{
    // These branch constants depend only on the material parameters.
    denom_s = 4.0 * lame * sin_phi * sin_phi +
              2.0 * shear * (1.0 + sin_phi) * (1.0 + sin_phi) +
              2.0 * shear * (1.0 - sin_phi) * (1.0 - sin_phi);
    denom_l = 4.0 * lame * sin_phi * sin_phi +
              shear * (1.0 + sin_phi) * (1.0 + sin_phi) +
              2.0 * shear * (1.0 - sin_phi) * (1.0 - sin_phi);
    denom_r = 4.0 * lame * sin_phi * sin_phi +
              2.0 * shear * (1.0 + sin_phi) * (1.0 + sin_phi) +
              shear * (1.0 - sin_phi) * (1.0 - sin_phi);
}

/**
  The function computes the returned principal stresses for the already known return
  branch.

  @param[in] eig - ordered trial eigenvalues
  @param[in] trace_E - trace of the ordered trial eigenvalues
  @param[in] return_type - selected return branch
  @param[in] lambda_s - smooth-face plastic multiplier
  @param[in] lambda_l - left-edge plastic multiplier
  @param[in] lambda_r - right-edge plastic multiplier
  @param[out] sigma - returned principal stresses

  @return The function does not return anything, the results are stored in the output
          arguments.
*/
void matmodel::compute_returned_principal_stresses(const vector &eig,
                                                   double trace_E,
                                                   int return_type,
                                                   double lambda_s,
                                                   double lambda_l,
                                                   double lambda_r,
                                                   vector &sigma) const
{
    // Recompute the returned principal stresses on the already selected branch.
    switch (return_type)
    {
    case RET_ELASTIC:
        sigma[0] = lame * trace_E + 2.0 * shear * eig[0];
        sigma[1] = lame * trace_E + 2.0 * shear * eig[1];
        sigma[2] = lame * trace_E + 2.0 * shear * eig[2];
        break;

    case RET_SMOOTH:
        sigma[0] = lame * trace_E + 2.0 * shear * eig[0] -
                   lambda_s * (2.0 * lame * sin_phi + 2.0 * shear * (1.0 + sin_phi));
        sigma[1] = lame * trace_E + 2.0 * shear * eig[1] -
                   lambda_s * (2.0 * lame * sin_phi);
        sigma[2] = lame * trace_E + 2.0 * shear * eig[2] -
                   lambda_s * (2.0 * lame * sin_phi - 2.0 * shear * (1.0 - sin_phi));
        break;

    case RET_LEFT_EDGE:
        sigma[0] = lame * trace_E + shear * (eig[0] + eig[1]) -
                   lambda_l * (2.0 * lame * sin_phi + shear * (1.0 + sin_phi));
        sigma[1] = sigma[0];
        sigma[2] = lame * trace_E + 2.0 * shear * eig[2] -
                   lambda_l * (2.0 * lame * sin_phi - 2.0 * shear * (1.0 - sin_phi));
        break;

    case RET_RIGHT_EDGE:
        sigma[0] = lame * trace_E + 2.0 * shear * eig[0] -
                   lambda_r * (2.0 * lame * sin_phi + 2.0 * shear * (1.0 + sin_phi));
        sigma[2] = lame * trace_E + shear * (eig[1] + eig[2]) -
                   lambda_r * (2.0 * lame * sin_phi - shear * (1.0 - sin_phi));
        sigma[1] = sigma[2];
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
  The function reads material model parameters and stress return setup from the opened text.

  @param[in] in - pointer to the opened text file.

  @retval 0 - on success
  @retval 1 - in the case of an error
*/
long matmodel::read(FILE *in)
{
    double young_read, poisson_read, cohesion_read, phi_read;
    double shear_read, bulk_read, lame_read, sin_phi_read, cos_phi_read, c_bar_read;

    if (in == NULL)
        return 1;

    if (fscanf(in, "%lf %lf %lf %lf", &young_read, &poisson_read, &cohesion_read, &phi_read) != 4)
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
    fprintf(out, "  E      = %.15g\n", young);
    fprintf(out, "  nu     = %.15g\n", poisson);
    fprintf(out, "  c      = %.15g\n", c);
    fprintf(out, "  phi    = %.15g rad\n", phi);
    fprintf(out, "  G      = %.15g\n", shear);
    fprintf(out, "  K      = %.15g\n", bulk);
    fprintf(out, "  lambda = %.15g\n", lame);
    fprintf(out, "  other  = %d components\n", NCOMP_OTHER);
    fprintf(out, "  eqother= %d components\n", NCOMP_EQOTHER);
    fprintf(out, "  other  = epsp(4), epsp_prev(4), return_type(1)\n");
    fprintf(out, "  eqother= epsp(4)\n");
}

/**
  The function computes actual stresses with respect to the attained strains and state
  variables.

  @param[in] strain - array of actual strain, components are ordered as follows
                      eps_x, eps_y, gamma_xy, eps_z - for the plane strain problem
                      eps_x, eps_y, eps_z, gamma_yz, gamma_xz, gamma_xy - for the space
                      stress problem

  @param[in] eqstatev - point buffer from the last converged state; nlstresses() reads
                        the plastic-strain history from its first 4 entries

  @param[out] stress - array of the resulting stress components, it must be computed,
                       components are ordered as follows
                       sig_x, sig_y, tau_xy, sig_z - for the plane strain problem
                       sig_x, sig_y, sig_z, tau_yz, tau_xz, tau_xy - for the space
                       stress problem

  @param[out] statev - array of the resulting current-point buffer for the actual strains

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
    vector E(NCOMP_STRAIN), Ep_prev(4), E_trial(NCOMP_STRAIN);
    vector Ep(4), eps_el(NCOMP_STRAIN), eig(3), sigma(3), S(NCOMP_STRESS);
    vector Eig_1(NCOMP_STRAIN), Eig_2(NCOMP_STRAIN), Eig_3(NCOMP_STRAIN);
    matrix EIG_1(NCOMP_STRAIN, NCOMP_STRAIN);
    matrix EIG_2(NCOMP_STRAIN, NCOMP_STRAIN);
    matrix EIG_3(NCOMP_STRAIN, NCOMP_STRAIN);
    double trace_E, f_tr, gamma_sl, gamma_sr, gamma_la, gamma_ra;
    double denom_s, denom_l, denom_r;
    double lambda_s, lambda_l, lambda_r;
    int return_type = RET_ELASTIC;

    // Resize the point outputs before filling the current material response.
    reallocv(NCOMP_STRESS, stress);
    reallocv(NCOMP_OTHER, statev);
    nullv(E);
    nullv(Ep_prev);
    nullv(E_trial);
    nullv(Ep);
    nullv(eps_el);
    nullv(eig);
    nullv(sigma);
    nullv(S);
    nullv(Eig_1);
    nullv(Eig_2);
    nullv(Eig_3);
    nullm(EIG_1);
    nullm(EIG_2);
    nullm(EIG_3);

    // Read the actual strain and the converged plastic-strain history.
    copyv(strain, E);
    copyv(&eqstatev[O_EP_XX], Ep_prev);

    // Matlab forms the trial strain by subtracting the old plastic strain.
    addmultv(E, 1.0, Ep_prev, -1.0, E_trial);

    // Rebuild the ordered trial eigenvalues and their derivatives from the trial strain.
    compute_ordered_trial_spectral(E_trial, eig, Eig_1, Eig_2, Eig_3, EIG_1, EIG_2, EIG_3);

    // Evaluate the Mohr-Coulomb trial function and the return candidates.
    trace_E = eig[0] + eig[1] + eig[2];
    f_tr = 2.0 * shear * ((1.0 + sin_phi) * eig[0] - (1.0 - sin_phi) * eig[2]) +
           2.0 * lame * sin_phi * trace_E - c_bar;
    gamma_sl = (eig[0] - eig[1]) / (1.0 + sin_phi);
    gamma_sr = (eig[1] - eig[2]) / (1.0 - sin_phi);
    gamma_la = (eig[0] + eig[1] - 2.0 * eig[2]) / (3.0 - sin_phi);
    gamma_ra = (2.0 * eig[0] - eig[1] - eig[2]) / (3.0 + sin_phi);

    // These denominators are the same branch constants as in the Matlab code.
    compute_return_denominators(denom_s, denom_l, denom_r);

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
    addmultv(S, Eig_1, sigma[0]);
    addmultv(S, Eig_2, sigma[1]);
    addmultv(S, Eig_3, sigma[2]);

    // Convert stress back to elastic strain and recover the updated plastic strain.
    eps_el[0] = (S[0] - poisson * (S[1] + S[3])) / young;
    eps_el[1] = (S[1] - poisson * (S[0] + S[3])) / young;
    eps_el[2] = S[2] / shear;
    eps_el[3] = (S[3] - poisson * (S[0] + S[1])) / young;
    addmultv(E, 1.0, eps_el, -1.0, Ep);

    // Pack only the persistent history and the minimal tangent helper payload.
    nullv(statev);
    copyv(Ep, &statev[O_EP_XX]);
    copyv(Ep_prev, &statev[O_EP_PREV_XX]);
    statev[O_RET] = static_cast<double>(return_type);

    // Return the stress vector to the caller.
    copyv(S, stress);
}

/**
  The function computes material stiffness matrix with respect to the attained strains
  and state variables.

  @param[in] strain - array of actual strain, components are ordered as follows
                      eps_x, eps_y, gamma_xy, eps_z - for the plane strain problem
                      eps_x, eps_y, eps_z, gamma_yz, gamma_xz, gamma_xy - for the space
                      stress problem

  @param[in] statev - current-point buffer produced by nlstresses(); this buffer stores
                      the updated plastic strain, the previous plastic strain, and the
                      selected return type

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
                        const vector &statev,
                        const vector &stress,
                        matrix &d)
{
    /*
      Scalar material-point version of Matlab stiffness_matrix.m with minimal
      recomputation. The tangent call reuses the stored return type and rebuilds the
      trial spectral data from strain and epsp_prev.
    */
    vector Ep_prev(4), E_trial(NCOMP_STRAIN), eig(3), sigma(3);
    vector iota(4);
    vector Eig_1(NCOMP_STRAIN), Eig_2(NCOMP_STRAIN), Eig_3(NCOMP_STRAIN);
    vector Eig_12(NCOMP_STRAIN), Eig_23(NCOMP_STRAIN), Eig_6(NCOMP_STRAIN);
    matrix EIG_1(NCOMP_STRAIN, NCOMP_STRAIN);
    matrix EIG_2(NCOMP_STRAIN, NCOMP_STRAIN);
    matrix EIG_3(NCOMP_STRAIN, NCOMP_STRAIN);
    matrix EIG_12(NCOMP_STRAIN, NCOMP_STRAIN);
    matrix EIG_23(NCOMP_STRAIN, NCOMP_STRAIN);
    matrix zero(NCOMP_STRAIN, NCOMP_STRAIN);
    matrix outer(NCOMP_STRAIN, NCOMP_STRAIN);
    double trace_E, f_tr, lambda_s, lambda_l, lambda_r;
    double denom_s, denom_l, denom_r;
    int return_type = RET_ELASTIC;
    (void)stress;

    // Allocate and clear the tangent and its local work arrays.
    reallocm(NCOMP_STRAIN, NCOMP_STRAIN, d);
    nullm(d);
    nullv(Ep_prev);
    nullv(E_trial);
    nullv(eig);
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

    // Read the stored return branch from the current point buffer.
    return_type = static_cast<int>(statev[O_RET] + 0.5);

    // Rebuild the trial strain from the actual strain and the previous plastic strain.
    copyv(&statev[O_EP_PREV_XX], Ep_prev);
    addmultv(strain, 1.0, Ep_prev, -1.0, E_trial);

    // Recompute the ordered trial spectral data needed by stiffness_matrix.m.
    compute_ordered_trial_spectral(E_trial, eig, Eig_1, Eig_2, Eig_3, EIG_1, EIG_2, EIG_3);

    // Rebuild the branch denominators and the branch-specific plastic multipliers.
    trace_E = eig[0] + eig[1] + eig[2];
    compute_return_denominators(denom_s, denom_l, denom_r);
    f_tr = 2.0 * shear * ((1.0 + sin_phi) * eig[0] - (1.0 - sin_phi) * eig[2]) +
           2.0 * lame * sin_phi * trace_E - c_bar;
    lambda_s = f_tr / denom_s;
    lambda_l = (shear * ((1.0 + sin_phi) * (eig[0] + eig[1]) - 2.0 * (1.0 - sin_phi) * eig[2]) +
                2.0 * lame * sin_phi * trace_E - c_bar) / denom_l;
    lambda_r = (shear * (2.0 * (1.0 + sin_phi) * eig[0] - (1.0 - sin_phi) * (eig[1] + eig[2])) +
                2.0 * lame * sin_phi * trace_E - c_bar) / denom_r;

    // Recompute the returned principal stresses on the already selected branch.
    compute_returned_principal_stresses(eig, trace_E, return_type, lambda_s, lambda_l, lambda_r, sigma);

    // The plane-strain unit tensor enters the volumetric part of every tangent branch.
    iota[0] = 1.0;
    iota[1] = 1.0;
    iota[3] = 1.0;

    // Assemble the 4x4 SIFEL tangent from the rebuilt Matlab spectral data.
    switch (return_type)
    {
    case RET_ELASTIC:
        // Matlab: Sderiv(:,test_el)=Elast(:)
        d(0, 0) = lame + 2.0 * shear;
        d(0, 1) = lame;
        d(0, 3) = lame;
        d(1, 0) = lame;
        d(1, 1) = lame + 2.0 * shear;
        d(1, 3) = lame;
        d(2, 2) = shear;
        d(3, 0) = lame;
        d(3, 1) = lame;
        d(3, 3) = lame + 2.0 * shear;
        break;

    case RET_SMOOTH:
        // Matlab smooth branch:
        // Sderiv(:,test_s)=mat1_s+mat2_s+mat3_s+mat4_s+mat5_s-mat6_s
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

        // Matlab: Eig_6=2*shear*((1+sin_phi)*Eig_1-(1-sin_phi)*Eig_3)+2*lame*sin_phi*[1;1;0]
        addmultv(Eig_1, 2.0 * shear * (1.0 + sin_phi), Eig_3, -2.0 * shear * (1.0 - sin_phi), Eig_6);
        addmultv(Eig_6, iota, 2.0 * lame * sin_phi);
        vxv(Eig_6, Eig_6, outer);
        addmultm(zero, 0.0, outer, -1.0 / denom_s, d);
        break;

    case RET_LEFT_EDGE:
        // Matlab left-edge branch:
        // Sderiv(:,test_l)=mat1_l+mat2_l+mat3_l+mat5_l-mat6_l
        addmultv(Eig_1, 1.0, Eig_2, 1.0, Eig_12);
        addm(EIG_1, EIG_2, EIG_12);
        addmultm(zero, 0.0, EIG_12, sigma[0], d);
        addmultm(zero, 0.0, EIG_3, sigma[2], d);
        vxv(iota, iota, outer);
        addmultm(zero, 0.0, outer, lame, d);
        vxv(Eig_12, Eig_12, outer);
        addmultm(zero, 0.0, outer, shear, d);
        vxv(Eig_3, Eig_3, outer);
        addmultm(zero, 0.0, outer, 2.0 * shear, d);

        // Matlab: Eig_12=Eig_1+Eig_2 and Eig_6=shear*((1+sin_phi)*Eig_12-2*(1-sin_phi)*Eig_3)+...
        addmultv(Eig_12, shear * (1.0 + sin_phi), Eig_3, -2.0 * shear * (1.0 - sin_phi), Eig_6);
        addmultv(Eig_6, iota, 2.0 * lame * sin_phi);
        vxv(Eig_6, Eig_6, outer);
        addmultm(zero, 0.0, outer, -1.0 / denom_l, d);
        break;

    case RET_RIGHT_EDGE:
        // Matlab right-edge branch:
        // Sderiv(:,test_r)=mat1_r+mat2_r+mat3_r+mat5_r-mat6_r
        addmultv(Eig_2, 1.0, Eig_3, 1.0, Eig_23);
        addm(EIG_2, EIG_3, EIG_23);
        addmultm(zero, 0.0, EIG_1, sigma[0], d);
        addmultm(zero, 0.0, EIG_23, sigma[2], d);
        vxv(iota, iota, outer);
        addmultm(zero, 0.0, outer, lame, d);
        vxv(Eig_1, Eig_1, outer);
        addmultm(zero, 0.0, outer, 2.0 * shear, d);
        vxv(Eig_23, Eig_23, outer);
        addmultm(zero, 0.0, outer, shear, d);

        // Matlab: Eig_23=Eig_2+Eig_3 and Eig_6=shear*(2*(1+sin_phi)*Eig_1-(1-sin_phi)*Eig_23)+...
        addmultv(Eig_1, 2.0 * shear * (1.0 + sin_phi), Eig_23, -shear * (1.0 - sin_phi), Eig_6);
        addmultv(Eig_6, iota, 2.0 * lame * sin_phi);
        vxv(Eig_6, Eig_6, outer);
        addmultm(zero, 0.0, outer, -1.0 / denom_r, d);
        break;

    case RET_APEX:
    default:
        // Matlab: Sderiv(:,test_a)=zeros(9,nt_a)
        break;
    }
}

/**
  The function updates the point buffer from the actual converged iteration values.

  @param[in] statev - current-point buffer produced by nlstresses()

  @param[out] eqstatev - point buffer for the last converged state

  @return The function does not return anything, the resulting converged state is stored
          in the array eqstatev passed in as an argument.
*/
void matmodel::updateval(const vector &statev, vector &eqstatev)
{
    // Persist only the updated plastic strain between converged steps.
    reallocv(NCOMP_EQOTHER, eqstatev);
    copyv(&statev[O_EP_XX], eqstatev);
}
