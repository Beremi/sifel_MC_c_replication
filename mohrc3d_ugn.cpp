#include "mohrc3d_ugn.h"
#include "vector.h"
#include "matrix.h"
#include "iotools.h"
#include "global.h"
#include "mechmat.h"

#include <algorithm>
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace {

const int SIFEL_TO_INTERNAL[6] = {0, 1, 2, 5, 3, 4};
const int INTERNAL_TO_SIFEL[6] = {0, 1, 2, 4, 5, 3};
const double IDENT_DIAG[6] = {1.0, 1.0, 1.0, 0.5, 0.5, 0.5};
const double IOTA[6] = {1.0, 1.0, 1.0, 0.0, 0.0, 0.0};

void zero6(double a[6])
{
  for (int i=0; i<6; i++)
    a[i] = 0.0;
}

void zero66(double a[6][6])
{
  for (int i=0; i<6; i++)
  {
    for (int j=0; j<6; j++)
      a[i][j] = 0.0;
  }
}

void add_outer(double a[6][6],
               const double u[6],
               const double v[6],
               double coeff)
{
  for (int i=0; i<6; i++)
  {
    for (int j=0; j<6; j++)
      a[i][j] += coeff * u[i] * v[j];
  }
}

void form_outer(const double u[6], const double v[6], double a[6][6])
{
  for (int i=0; i<6; i++)
  {
    for (int j=0; j<6; j++)
      a[i][j] = u[i] * v[j];
  }
}

void add_scaled_matrix(double dst[6][6],
                       const double src[6][6],
                       double coeff)
{
  for (int i=0; i<6; i++)
  {
    for (int j=0; j<6; j++)
      dst[i][j] += coeff * src[i][j];
  }
}

void add_iota_outer(double a[6][6], double coeff)
{
  for (int i=0; i<6; i++)
  {
    for (int j=0; j<6; j++)
      a[i][j] += coeff * IOTA[i] * IOTA[j];
  }
}

void elastic_tangent(double lame, double shear, double d[6][6])
{
  zero66(d);
  add_iota_outer(d, lame);
  for (int i=0; i<6; i++)
    d[i][i] += 2.0 * shear * IDENT_DIAG[i];
}

void array_sifel_to_internal(const double sifel[6], double internal[6])
{
  for (int i=0; i<6; i++)
    internal[i] = sifel[SIFEL_TO_INTERNAL[i]];
}

void array_internal_to_sifel(const double internal[6], double sifel[6])
{
  for (int i=0; i<6; i++)
    sifel[i] = internal[INTERNAL_TO_SIFEL[i]];
}

void matrix_internal_to_sifel(const double internal[6][6], matrix &sifel)
{
  reallocm(6, 6, sifel);
  for (int i=0; i<6; i++)
  {
    for (int j=0; j<6; j++)
      sifel(i,j) = internal[INTERNAL_TO_SIFEL[i]][INTERNAL_TO_SIFEL[j]];
  }
}

void vector_to_array6(const vector &src, double dst[6])
{
  for (int i=0; i<6; i++)
    dst[i] = src[i];
}

void internal_trial_from_sifel(const vector &strain,
                               const double ep_prev_sifel[6],
                               double e_trial_internal[6])
{
  double strain_sifel[6], ep_prev_internal[6], strain_internal[6];

  vector_to_array6(strain, strain_sifel);
  array_sifel_to_internal(strain_sifel, strain_internal);
  array_sifel_to_internal(ep_prev_sifel, ep_prev_internal);

  for (int i=0; i<6; i++)
    e_trial_internal[i] = strain_internal[i] - ep_prev_internal[i];
}

void compute_e_tr_and_square(const double e_trial[6],
                             double e_tr[6],
                             double e_square[6])
{
  for (int i=0; i<6; i++)
    e_tr[i] = IDENT_DIAG[i] * e_trial[i];

  e_square[0] = e_tr[0]*e_tr[0] + e_tr[3]*e_tr[3] + e_tr[5]*e_tr[5];
  e_square[1] = e_tr[1]*e_tr[1] + e_tr[3]*e_tr[3] + e_tr[4]*e_tr[4];
  e_square[2] = e_tr[2]*e_tr[2] + e_tr[4]*e_tr[4] + e_tr[5]*e_tr[5];
  e_square[3] = e_tr[0]*e_tr[3] + e_tr[1]*e_tr[3] + e_tr[4]*e_tr[5];
  e_square[4] = e_tr[3]*e_tr[5] + e_tr[1]*e_tr[4] + e_tr[2]*e_tr[4];
  e_square[5] = e_tr[0]*e_tr[5] + e_tr[3]*e_tr[4] + e_tr[2]*e_tr[5];
}

void compute_der_e_square(const double e_tr[6], double der[6][6])
{
  zero66(der);

  der[0][0] = 2.0 * e_tr[0];
  der[3][0] = e_tr[3];
  der[5][0] = e_tr[5];

  der[1][1] = 2.0 * e_tr[1];
  der[3][1] = e_tr[3];
  der[4][1] = e_tr[4];

  der[2][2] = 2.0 * e_tr[2];
  der[4][2] = e_tr[4];
  der[5][2] = e_tr[5];

  der[0][3] = e_tr[3];
  der[1][3] = e_tr[3];
  der[3][3] = 0.5 * (e_tr[0] + e_tr[1]);
  der[4][3] = 0.5 * e_tr[5];
  der[5][3] = 0.5 * e_tr[4];

  der[1][4] = e_tr[4];
  der[2][4] = e_tr[4];
  der[3][4] = 0.5 * e_tr[5];
  der[4][4] = 0.5 * (e_tr[1] + e_tr[2]);
  der[5][4] = 0.5 * e_tr[3];

  der[0][5] = e_tr[5];
  der[2][5] = e_tr[5];
  der[3][5] = 0.5 * e_tr[4];
  der[4][5] = 0.5 * e_tr[3];
  der[5][5] = 0.5 * (e_tr[0] + e_tr[2]);
}

void compute_eigenvalues(const double e_tr[6], double eig[3])
{
  const double I1 = e_tr[0] + e_tr[1] + e_tr[2];
  const double I2 = e_tr[0]*e_tr[1] + e_tr[0]*e_tr[2] + e_tr[1]*e_tr[2] -
                    e_tr[3]*e_tr[3] - e_tr[4]*e_tr[4] - e_tr[5]*e_tr[5];
  const double I3 = e_tr[0]*e_tr[1]*e_tr[2] - e_tr[2]*e_tr[3]*e_tr[3] -
                    e_tr[1]*e_tr[5]*e_tr[5] - e_tr[0]*e_tr[4]*e_tr[4] +
                    2.0*e_tr[3]*e_tr[4]*e_tr[5];
  const double Q_raw = (I1*I1 - 3.0*I2) / 9.0;
  const double Q = (Q_raw > 0.0) ? Q_raw : 0.0;
  const double R = (-2.0*I1*I1*I1 + 9.0*I1*I2 - 27.0*I3) / 54.0;
  double theta0 = 0.0;

  if (Q > 0.0)
    theta0 = R / sqrt(Q*Q*Q);
  theta0 = std::max(-1.0, std::min(1.0, theta0));

  const double theta = acos(theta0) / 3.0;
  eig[0] = -2.0*sqrt(Q)*cos(theta + 2.0*M_PI/3.0) + I1/3.0;
  eig[1] = -2.0*sqrt(Q)*cos(theta - 2.0*M_PI/3.0) + I1/3.0;
  eig[2] = -2.0*sqrt(Q)*cos(theta) + I1/3.0;
}

void compute_projection(const double e_square[6],
                        const double e_tr[6],
                        double eig_i,
                        double eig_j,
                        double eig_k,
                        double proj[6])
{
  const double denom = (eig_i - eig_j) * (eig_i - eig_k);

  for (int i=0; i<6; i++)
    proj[i] = (e_square[i] - (eig_j + eig_k) * e_tr[i] + IOTA[i] * eig_j * eig_k) / denom;
}

void complement_projection(const double proj[6], double complement[6])
{
  for (int i=0; i<6; i++)
    complement[i] = IOTA[i] - proj[i];
}

void compute_smooth_projection_derivative(const double der_e_square[6][6],
                                          const double proj_i[6],
                                          const double proj_j[6],
                                          const double proj_k[6],
                                          double eig_i,
                                          double eig_j,
                                          double eig_k,
                                          double hess[6][6])
{
  const double denom = (eig_i - eig_j) * (eig_i - eig_k);
  double outer_i[6][6], outer_j[6][6], outer_k[6][6];

  form_outer(proj_i, proj_i, outer_i);
  form_outer(proj_j, proj_j, outer_j);
  form_outer(proj_k, proj_k, outer_k);

  for (int r=0; r<6; r++)
  {
    for (int c=0; c<6; c++)
    {
      double v = der_e_square[r][c];
      if (r == c)
        v -= IDENT_DIAG[r] * (eig_j + eig_k);
      v -= (2.0*eig_i - eig_j - eig_k) * outer_i[r][c];
      v -= (eig_j - eig_k) * (outer_j[r][c] - outer_k[r][c]);
      hess[r][c] = v / denom;
    }
  }
}

void compute_left_edge_projection_derivative(const double der_e_square[6][6],
                                             const double e_tr[6],
                                             const double eig12[6],
                                             const double eig3[6],
                                             double eig_1,
                                             double eig_2,
                                             double eig_3,
                                             double hess3[6][6])
{
  const double denom = (eig_3 - eig_1) * (eig_3 - eig_2);
  double e12_x_e12[6][6], e3_x_e3[6][6], e12_x_etr[6][6], etr_x_e12[6][6];
  double e12_x_e3[6][6], e3_x_e12[6][6];

  form_outer(eig12, eig12, e12_x_e12);
  form_outer(eig3, eig3, e3_x_e3);
  form_outer(eig12, e_tr, e12_x_etr);
  form_outer(e_tr, eig12, etr_x_e12);
  form_outer(eig12, eig3, e12_x_e3);
  form_outer(eig3, eig12, e3_x_e12);

  for (int r=0; r<6; r++)
  {
    for (int c=0; c<6; c++)
    {
      double v = der_e_square[r][c];
      if (r == c)
        v -= IDENT_DIAG[r] * (eig_1 + eig_2);
      v -= etr_x_e12[r][c] + e12_x_etr[r][c];
      v += (eig_1 + eig_2) * e12_x_e12[r][c];
      v += (eig_1 + eig_2 - 2.0*eig_3) * e3_x_e3[r][c];
      v += eig_3 * (e12_x_e3[r][c] + e3_x_e12[r][c]);
      hess3[r][c] = v / denom;
    }
  }
}

void compute_right_edge_projection_derivative(const double der_e_square[6][6],
                                              const double e_tr[6],
                                              const double eig1[6],
                                              const double eig23[6],
                                              double eig_1,
                                              double eig_2,
                                              double eig_3,
                                              double hess1[6][6])
{
  const double denom = (eig_1 - eig_2) * (eig_1 - eig_3);
  double e1_x_e1[6][6], e23_x_e23[6][6], e23_x_etr[6][6], etr_x_e23[6][6];
  double e23_x_e1[6][6], e1_x_e23[6][6];

  form_outer(eig1, eig1, e1_x_e1);
  form_outer(eig23, eig23, e23_x_e23);
  form_outer(eig23, e_tr, e23_x_etr);
  form_outer(e_tr, eig23, etr_x_e23);
  form_outer(eig23, eig1, e23_x_e1);
  form_outer(eig1, eig23, e1_x_e23);

  for (int r=0; r<6; r++)
  {
    for (int c=0; c<6; c++)
    {
      double v = der_e_square[r][c];
      if (r == c)
        v -= IDENT_DIAG[r] * (eig_2 + eig_3);
      v -= etr_x_e23[r][c] + e23_x_etr[r][c];
      v += (eig_2 + eig_3) * e23_x_e23[r][c];
      v += (eig_2 + eig_3 - 2.0*eig_1) * e1_x_e1[r][c];
      v += eig_1 * (e23_x_e1[r][c] + e1_x_e23[r][c]);
      hess1[r][c] = v / denom;
    }
  }
}

void stress_to_elastic_strain_sifel(const double stress_sifel[6],
                                    double young,
                                    double poisson,
                                    double shear,
                                    double eps_el_sifel[6])
{
  eps_el_sifel[0] = (stress_sifel[0] - poisson * (stress_sifel[1] + stress_sifel[2])) / young;
  eps_el_sifel[1] = (stress_sifel[1] - poisson * (stress_sifel[0] + stress_sifel[2])) / young;
  eps_el_sifel[2] = (stress_sifel[2] - poisson * (stress_sifel[0] + stress_sifel[1])) / young;
  eps_el_sifel[3] = stress_sifel[3] / shear;
  eps_el_sifel[4] = stress_sifel[4] / shear;
  eps_el_sifel[5] = stress_sifel[5] / shear;
}

} // namespace

void mohrc3d_ugn::compute_response_internal(const double strain_internal[6],
                                            int forced_return_type,
                                            double stress_internal[6],
                                            double tangent_internal[6][6],
                                            int &return_type,
                                            double eig[3],
                                            double sigma[3]) const
{
  double e_tr[6], e_square[6], der_e_square[6][6];
  double eig1[6], eig2[6], eig3[6], eig12[6], eig23[6];
  double hess1[6][6], hess2[6][6], hess3[6][6];
  double d_phi[6], d_psi[6];
  double lambda_s, lambda_l, lambda_r, lambda_a;

  zero6(stress_internal);
  zero66(tangent_internal);
  zero6(eig1);
  zero6(eig2);
  zero6(eig3);
  zero6(eig12);
  zero6(eig23);
  zero66(hess1);
  zero66(hess2);
  zero66(hess3);
  zero6(d_phi);
  zero6(d_psi);

  compute_e_tr_and_square(strain_internal, e_tr, e_square);
  compute_der_e_square(e_tr, der_e_square);
  compute_eigenvalues(e_tr, eig);

  const double I1 = eig[0] + eig[1] + eig[2];
  const double f_tr = 2.0*shear*((1.0 + sin_phi)*eig[0] - (1.0 - sin_phi)*eig[2]) +
                      2.0*lame*sin_phi*I1 - c_bar;
  const double gamma_sl = (eig[0] - eig[1]) / (1.0 + sin_phi);
  const double gamma_sr = (eig[1] - eig[2]) / (1.0 - sin_phi);
  const double gamma_la = (eig[0] + eig[1] - 2.0*eig[2]) / (3.0 - sin_phi);
  const double gamma_ra = (2.0*eig[0] - eig[1] - eig[2]) / (3.0 + sin_phi);
  const double denom_s = 4.0*lame*sin_phi*sin_phi + 4.0*shear*(1.0 + sin_phi*sin_phi);
  const double denom_l = 4.0*lame*sin_phi*sin_phi +
                         shear*(1.0 + sin_phi)*(1.0 + sin_phi) +
                         2.0*shear*(1.0 - sin_phi)*(1.0 - sin_phi);
  const double denom_r = 4.0*lame*sin_phi*sin_phi +
                         2.0*shear*(1.0 + sin_phi)*(1.0 + sin_phi) +
                         shear*(1.0 - sin_phi)*(1.0 - sin_phi);
  const double denom_a = 4.0*bulk*sin_phi*sin_phi;

  lambda_s = f_tr / denom_s;
  lambda_l = (shear*((1.0 + sin_phi)*(eig[0] + eig[1]) - 2.0*(1.0 - sin_phi)*eig[2]) +
              2.0*lame*sin_phi*I1 - c_bar) / denom_l;
  lambda_r = (shear*(2.0*(1.0 + sin_phi)*eig[0] - (1.0 - sin_phi)*(eig[1] + eig[2])) +
              2.0*lame*sin_phi*I1 - c_bar) / denom_r;
  lambda_a = (2.0*bulk*sin_phi*I1 - c_bar) / denom_a;
  (void)lambda_a;

  if (forced_return_type >= 0)
    return_type = forced_return_type;
  else if (f_tr <= 0.0)
    return_type = RET_ELASTIC;
  else if (lambda_s <= std::min(gamma_sl, gamma_sr))
    return_type = RET_SMOOTH;
  else if ((gamma_sl < gamma_sr) && (lambda_l >= gamma_sl) && (lambda_l <= gamma_la))
    return_type = RET_LEFT_EDGE;
  else if ((gamma_sl > gamma_sr) && (lambda_r >= gamma_sr) && (lambda_r <= gamma_ra))
    return_type = RET_RIGHT_EDGE;
  else
    return_type = RET_APEX;

  sigma[0] = 0.0;
  sigma[1] = 0.0;
  sigma[2] = 0.0;

  switch (return_type){
    case RET_ELASTIC:
      sigma[0] = lame*I1 + 2.0*shear*eig[0];
      sigma[1] = lame*I1 + 2.0*shear*eig[1];
      sigma[2] = lame*I1 + 2.0*shear*eig[2];
      for (int i=0; i<6; i++)
        stress_internal[i] = lame*I1*IOTA[i] + 2.0*shear*e_tr[i];
      elastic_tangent(lame, shear, tangent_internal);
      break;

    case RET_SMOOTH:
      compute_projection(e_square, e_tr, eig[0], eig[1], eig[2], eig1);
      compute_projection(e_square, e_tr, eig[1], eig[0], eig[2], eig2);
      compute_projection(e_square, e_tr, eig[2], eig[0], eig[1], eig3);
      sigma[0] = lame*I1 + 2.0*shear*eig[0] -
                 lambda_s*(2.0*lame*sin_phi + 2.0*shear*(1.0 + sin_phi));
      sigma[1] = lame*I1 + 2.0*shear*eig[1] -
                 lambda_s*(2.0*lame*sin_phi);
      sigma[2] = lame*I1 + 2.0*shear*eig[2] -
                 lambda_s*(2.0*lame*sin_phi - 2.0*shear*(1.0 - sin_phi));
      for (int i=0; i<6; i++)
        stress_internal[i] = sigma[0]*eig1[i] + sigma[1]*eig2[i] + sigma[2]*eig3[i];

      compute_smooth_projection_derivative(der_e_square, eig1, eig2, eig3, eig[0], eig[1], eig[2], hess1);
      compute_smooth_projection_derivative(der_e_square, eig2, eig1, eig3, eig[1], eig[0], eig[2], hess2);
      compute_smooth_projection_derivative(der_e_square, eig3, eig1, eig2, eig[2], eig[0], eig[1], hess3);
      add_scaled_matrix(tangent_internal, hess1, sigma[0]);
      add_scaled_matrix(tangent_internal, hess2, sigma[1]);
      add_scaled_matrix(tangent_internal, hess3, sigma[2]);
      add_iota_outer(tangent_internal, lame);
      add_outer(tangent_internal, eig1, eig1, 2.0*shear);
      add_outer(tangent_internal, eig2, eig2, 2.0*shear);
      add_outer(tangent_internal, eig3, eig3, 2.0*shear);
      for (int i=0; i<6; i++)
      {
        d_phi[i] = 2.0*shear*((1.0 + sin_phi)*eig1[i] - (1.0 - sin_phi)*eig3[i]) +
                   2.0*lame*sin_phi*IOTA[i];
        d_psi[i] = d_phi[i];
      }
      add_outer(tangent_internal, d_psi, d_phi, -1.0 / denom_s);
      break;

    case RET_LEFT_EDGE:
      compute_projection(e_square, e_tr, eig[2], eig[0], eig[1], eig3);
      complement_projection(eig3, eig12);
      sigma[0] = lame*I1 + shear*(eig[0] + eig[1]) -
                 lambda_l*(2.0*lame*sin_phi + shear*(1.0 + sin_phi));
      sigma[1] = sigma[0];
      sigma[2] = lame*I1 + 2.0*shear*eig[2] -
                 lambda_l*(2.0*lame*sin_phi - 2.0*shear*(1.0 - sin_phi));
      for (int i=0; i<6; i++)
        stress_internal[i] = sigma[0]*eig12[i] + sigma[2]*eig3[i];

      compute_left_edge_projection_derivative(der_e_square, e_tr, eig12, eig3, eig[0], eig[1], eig[2], hess3);
      add_scaled_matrix(tangent_internal, hess3, sigma[2] - sigma[0]);
      add_iota_outer(tangent_internal, lame);
      add_outer(tangent_internal, eig12, eig12, shear);
      add_outer(tangent_internal, eig3, eig3, 2.0*shear);
      for (int i=0; i<6; i++)
      {
        d_phi[i] = shear*((1.0 + sin_phi)*eig12[i] - 2.0*(1.0 - sin_phi)*eig3[i]) +
                   2.0*lame*sin_phi*IOTA[i];
        d_psi[i] = d_phi[i];
      }
      add_outer(tangent_internal, d_psi, d_phi, -1.0 / denom_l);
      break;

    case RET_RIGHT_EDGE:
      compute_projection(e_square, e_tr, eig[0], eig[1], eig[2], eig1);
      complement_projection(eig1, eig23);
      sigma[0] = lame*I1 + 2.0*shear*eig[0] -
                 lambda_r*(2.0*lame*sin_phi + 2.0*shear*(1.0 + sin_phi));
      sigma[1] = lame*I1 + shear*(eig[1] + eig[2]) -
                 lambda_r*(2.0*lame*sin_phi - shear*(1.0 - sin_phi));
      sigma[2] = sigma[1];
      for (int i=0; i<6; i++)
        stress_internal[i] = sigma[0]*eig1[i] + sigma[2]*eig23[i];

      compute_right_edge_projection_derivative(der_e_square, e_tr, eig1, eig23, eig[0], eig[1], eig[2], hess1);
      add_scaled_matrix(tangent_internal, hess1, sigma[0] - sigma[2]);
      add_iota_outer(tangent_internal, lame);
      add_outer(tangent_internal, eig1, eig1, 2.0*shear);
      add_outer(tangent_internal, eig23, eig23, shear);
      for (int i=0; i<6; i++)
      {
        d_phi[i] = shear*(2.0*(1.0 + sin_phi)*eig1[i] - (1.0 - sin_phi)*eig23[i]) +
                   2.0*lame*sin_phi*IOTA[i];
        d_psi[i] = d_phi[i];
      }
      add_outer(tangent_internal, d_psi, d_phi, -1.0 / denom_r);
      break;

    case RET_APEX:
    default:
      sigma[0] = c_bar / (2.0 * sin_phi);
      sigma[1] = sigma[0];
      sigma[2] = sigma[0];
      for (int i=0; i<6; i++)
        stress_internal[i] = IOTA[i] * sigma[0];
      zero66(tangent_internal);
      break;
  }
}

long mohrc3d_ugn::read(XFILE *in)
{
  double young_read, poisson_read, cohesion_read, phi_read;

  if (in == NULL)
    return 1;

  if (xfscanf(in, "%lf %lf %lf %lf", &young_read, &poisson_read, &cohesion_read, &phi_read) != 4)
    return 1;

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
  shear = young / (2.0 * (1.0 + poisson));
  bulk = young / (3.0 * (1.0 - 2.0 * poisson));
  lame = bulk - 2.0 * shear / 3.0;
  sin_phi = sin(phi);
  cos_phi = cos(phi);
  c_bar = 2.0 * c * cos_phi;
  return 0;
}

void mohrc3d_ugn::print(FILE *out)
{
  if (out == NULL)
    return;

  fprintf(out, "%le %le %le %le ", young, poisson, c, phi);
}

void mohrc3d_ugn::nlstresses(long ipp, long ido)
{
  long ncompstr = Mm->ip[ipp].ncompstr;
  vector eps(ASTCKVEC(ncompstr));
  vector sig(ASTCKVEC(ncompstr));
  vector other;
  vector eqother(ASTCKVEC(NCOMP_EQOTHER));
  const long lcid = 0;

  Mm->givestrain(lcid, ipp, eps);
  Mm->givestress(lcid, ipp, sig);
  other.makerefv(NCOMP_OTHER, Mm->ip[ipp].other+ido);
  eqother.makerefv(NCOMP_EQOTHER, Mm->ip[ipp].eqother+ido);

  nlstresses(eps, eqother, sig, other);

  Mm->storestress(lcid, ipp, sig);
}

void mohrc3d_ugn::nonloc_nlstresses(long ipp, long ido)
{
  long ncompstr = Mm->ip[ipp].ncompstr;
  vector eps(ASTCKVEC(ncompstr));
  vector sig(ASTCKVEC(ncompstr));
  vector other;
  vector eqother;
  const long lcid = 0;

  Mm->givestrain(lcid, ipp, eps);
  Mm->givestress(lcid, ipp, sig);
  other.makerefv(NCOMP_OTHER, Mm->ip[ipp].other+ido);
  eqother.makerefv(NCOMP_EQOTHER, Mm->ip[ipp].nonloc);

  nlstresses(eps, eqother, sig, other);

  Mm->storestress(lcid, ipp, sig);
}

void mohrc3d_ugn::nlstresses(const vector &strain,
                             const vector &eqstatev,
                             vector &stress,
                             vector &statev)
{
  double strain_sifel[6], ep_prev_sifel[6], e_trial_internal[6];
  double stress_internal[6], stress_sifel[6], tangent_internal[6][6];
  double eps_el_sifel[6], epsp_sifel[6], eig[3], sigma[3];
  int return_type = RET_ELASTIC;

  reallocv(NCOMP_STRESS, stress);
  reallocv(NCOMP_OTHER, statev);
  nullv(stress);
  nullv(statev);

  vector_to_array6(strain, strain_sifel);
  for (int i=0; i<6; i++)
    ep_prev_sifel[i] = eqstatev[i];

  internal_trial_from_sifel(strain, ep_prev_sifel, e_trial_internal);
  compute_response_internal(e_trial_internal, -1, stress_internal, tangent_internal, return_type, eig, sigma);
  (void)tangent_internal;
  (void)eig;
  (void)sigma;

  array_internal_to_sifel(stress_internal, stress_sifel);
  stress_to_elastic_strain_sifel(stress_sifel, young, poisson, shear, eps_el_sifel);

  for (int i=0; i<6; i++)
  {
    epsp_sifel[i] = strain_sifel[i] - eps_el_sifel[i];
    statev[O_EP_XX + i] = epsp_sifel[i];
    statev[O_EP_PREV_XX + i] = ep_prev_sifel[i];
    stress[i] = stress_sifel[i];
  }

  statev[O_RET] = static_cast<double>(return_type);
}

void mohrc3d_ugn::matstiff(matrix &d, long ipp, long ido)
{
  long ncompstr = Mm->ip[ipp].ncompstr;
  vector eps(ASTCKVEC(ncompstr));
  vector sig(ASTCKVEC(ncompstr));
  vector other;
  const long lcid = 0;

  Mm->givestrain(lcid, ipp, eps);
  Mm->givestress(lcid, ipp, sig);
  other.makerefv(NCOMP_OTHER, Mm->ip[ipp].other+ido);

  stiffmat(eps, other, sig, d);
}

void mohrc3d_ugn::stiffmat(const vector &strain,
                           const vector &statev,
                           const vector &stress,
                           matrix &d)
{
  double ep_prev_sifel[6], e_trial_internal[6];
  double stress_internal[6], tangent_internal[6][6], eig[3], sigma[3];
  int return_type = static_cast<int>(statev[O_RET] + 0.5);
  (void)stress;

  for (int i=0; i<6; i++)
    ep_prev_sifel[i] = statev[O_EP_PREV_XX + i];

  internal_trial_from_sifel(strain, ep_prev_sifel, e_trial_internal);
  compute_response_internal(e_trial_internal, return_type, stress_internal, tangent_internal, return_type, eig, sigma);
  (void)stress_internal;
  (void)eig;
  (void)sigma;

  matrix_internal_to_sifel(tangent_internal, d);
}

void mohrc3d_ugn::updateval(long ipp, long im, long ido)
{
  vector statev;
  vector eqstatev;

  makerefv(statev, Mm->ip[ipp].other+ido, NCOMP_OTHER);

  long idoeq = Mm->givencompeqother(ipp, 0) - Mm->givencompeqother(ipp, im);
  makerefv(eqstatev, Mm->ip[ipp].eqother+idoeq, NCOMP_EQOTHER);

  updateval(statev, eqstatev);
}

void mohrc3d_ugn::giveirrstrains(long ipp, long ido, vector &epsp)
{
  for (long i=0; i<epsp.n; i++)
    epsp[i] = Mm->ip[ipp].other[ido+i];
}

void mohrc3d_ugn::updateval(const vector &statev, vector &eqstatev)
{
  reallocv(NCOMP_EQOTHER, eqstatev);
  for (int i=0; i<6; i++)
    eqstatev[i] = statev[O_EP_XX + i];
}
