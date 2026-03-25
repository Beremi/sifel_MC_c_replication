#include "matmodel.h"
#include "vector.h"
#include "matrix.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace {

struct matmodel_cache
{
  double epsp[MATMODEL_NCOMP_EQOTHER];
  int return_type;

  double eig[3];
  double proj[3][MATMODEL_NCOMP_STRAIN];
  double hess[3][MATMODEL_NCOMP_STRAIN][MATMODEL_NCOMP_STRAIN];
  double sig_princ[3];
  double stress[MATMODEL_NCOMP_STRESS];
};

void init_params(matmodel_params *par,
                 double young,
                 double poisson,
                 double cohesion,
                 double phi)
{
  if (par == NULL)
    return;

  par->young = young;
  par->poisson = poisson;
  par->c = cohesion;
  par->phi = phi;
  par->shear = young/(2.0*(1.0 + poisson));
  par->bulk = young/(3.0*(1.0 - 2.0*poisson));
  par->lame = par->bulk - 2.0*par->shear/3.0;
  par->sin_phi = sin(phi);
  par->cos_phi = cos(phi);
  par->c_bar = 2.0*cohesion*par->cos_phi;
}

int validate_params(const matmodel_params *par)
{
  if (par == NULL)
    return 0;
  if (par->young <= 0.0)
    return 0;
  if ((par->poisson <= -1.0) || (par->poisson >= 0.5))
    return 0;
  if (par->c < 0.0)
    return 0;
  if ((par->phi <= 1.0e-14) || (par->phi >= 0.5*M_PI - 1.0e-14))
    return 0;
  return 1;
}

void compliance_times_stress(const matmodel_params *par,
                             const double stress[4],
                             double eps_el[4])
{
  eps_el[0] = (stress[0] - par->poisson*(stress[1] + stress[3]))/par->young;
  eps_el[1] = (stress[1] - par->poisson*(stress[0] + stress[3]))/par->young;
  eps_el[2] = stress[2]/par->shear;
  eps_el[3] = (stress[3] - par->poisson*(stress[0] + stress[1]))/par->young;
}

void compute_inplane_spectral(const double trial[4],
                              double eig[3],
                              double proj[3][4],
                              double hess[3][4][4])
{
  static const double metric[4][4] = {
    {1.0, 0.0, 0.0, 0.0},
    {0.0, 1.0, 0.0, 0.0},
    {0.0, 0.0, 0.5, 0.0},
    {0.0, 0.0, 0.0, 0.0}
  };

  const double ex = trial[0];
  const double ey = trial[1];
  const double g = trial[2];
  const double d = ex - ey;
  const double i1 = ex + ey;
  const double i2 = sqrt(d*d + g*g);
  const double tol = 1.0e-14;

  eig[0] = 0.5*(i1 + i2);
  eig[1] = 0.5*(i1 - i2);
  eig[2] = trial[3];

  memset(proj, 0, 3*4*sizeof(double));
  memset(hess, 0, 3*4*4*sizeof(double));

  if (i2 <= tol)
  {
    proj[0][0] = 1.0;
    proj[0][1] = 1.0;
    proj[1][0] = 0.0;
    proj[1][1] = 0.0;
  }
  else
  {
    proj[0][0] = (ex - eig[1])/i2;
    proj[0][1] = (ey - eig[1])/i2;
    proj[0][2] = g/(2.0*i2);

    proj[1][0] = 1.0 - proj[0][0];
    proj[1][1] = 1.0 - proj[0][1];
    proj[1][2] = -proj[0][2];

    for (long a=0; a<4; a++)
    {
      for (long b=0; b<4; b++)
      {
        hess[0][a][b] = (metric[a][b] - proj[0][a]*proj[0][b] - proj[1][a]*proj[1][b])/i2;
        hess[1][a][b] = -hess[0][a][b];
      }
    }

    for (long a=0; a<4; a++)
    {
      hess[0][3][a] = 0.0;
      hess[0][a][3] = 0.0;
      hess[1][3][a] = 0.0;
      hess[1][a][3] = 0.0;
    }
  }

  proj[2][3] = 1.0;
}

void copy_mode(double dst_proj[4],
               double dst_hess[4][4],
               double *dst_eig,
               const double src_proj[4],
               const double src_hess[4][4],
               double src_eig)
{
  copyv(src_proj, dst_proj, 4);
  copyv(&src_hess[0][0], &dst_hess[0][0], 16);
  *dst_eig = src_eig;
}

void elastic_matrix(const matmodel_params *par, matrix &d)
{
  reallocm(4, 4, d);
  nullm(d);
  if (par == NULL)
    return;

  d(0,0) = par->lame + 2.0*par->shear;
  d(0,1) = par->lame;
  d(0,3) = par->lame;

  d(1,0) = par->lame;
  d(1,1) = par->lame + 2.0*par->shear;
  d(1,3) = par->lame;

  d(2,2) = par->shear;

  d(3,0) = par->lame;
  d(3,1) = par->lame;
  d(3,3) = par->lame + 2.0*par->shear;
}

void compute_response(const matmodel_params *par,
                      const double strain[4],
                      const double eqother[4],
                      matmodel_cache *cache)
{
  double trial[4], epsp_prev[4], eps_el[4];
  double eig0[3], proj0[3][4], hess0[3][4][4];
  double trace_E, f_tr, gamma_sl, gamma_sr, gamma_la, gamma_ra;
  double denom_s, denom_l, denom_r;
  double lambda_s, lambda_l, lambda_r;

  if ((par == NULL) || (cache == NULL) || (strain == NULL))
    return;

  memset(cache, 0, sizeof(*cache));
  nullv(epsp_prev, 4);
  if (eqother != NULL)
    copyv(eqother, epsp_prev, 4);

  for (long i=0; i<4; i++)
    trial[i] = strain[i] - epsp_prev[i];

  compute_inplane_spectral(trial, eig0, proj0, hess0);

  copy_mode(cache->proj[0], cache->hess[0], &cache->eig[0], proj0[0], hess0[0], eig0[0]);
  copy_mode(cache->proj[1], cache->hess[1], &cache->eig[1], proj0[1], hess0[1], eig0[1]);
  copy_mode(cache->proj[2], cache->hess[2], &cache->eig[2], proj0[2], hess0[2], eig0[2]);

  if ((eig0[0] >= eig0[2]) && (eig0[2] > eig0[1]))
  {
    copy_mode(cache->proj[1], cache->hess[1], &cache->eig[1], proj0[2], hess0[2], eig0[2]);
    copy_mode(cache->proj[2], cache->hess[2], &cache->eig[2], proj0[1], hess0[1], eig0[1]);
  }
  if (eig0[2] > eig0[0])
  {
    copy_mode(cache->proj[0], cache->hess[0], &cache->eig[0], proj0[2], hess0[2], eig0[2]);
    copy_mode(cache->proj[1], cache->hess[1], &cache->eig[1], proj0[0], hess0[0], eig0[0]);
    copy_mode(cache->proj[2], cache->hess[2], &cache->eig[2], proj0[1], hess0[1], eig0[1]);
  }

  trace_E = cache->eig[0] + cache->eig[1] + cache->eig[2];
  f_tr = 2.0*par->shear*((1.0 + par->sin_phi)*cache->eig[0] - (1.0 - par->sin_phi)*cache->eig[2])
       + 2.0*par->lame*par->sin_phi*trace_E - par->c_bar;

  gamma_sl = (cache->eig[0] - cache->eig[1])/(1.0 + par->sin_phi);
  gamma_sr = (cache->eig[1] - cache->eig[2])/(1.0 - par->sin_phi);
  gamma_la = (cache->eig[0] + cache->eig[1] - 2.0*cache->eig[2])/(3.0 - par->sin_phi);
  gamma_ra = (2.0*cache->eig[0] - cache->eig[1] - cache->eig[2])/(3.0 + par->sin_phi);

  denom_s = 4.0*par->lame*par->sin_phi*par->sin_phi
          + 2.0*par->shear*(1.0 + par->sin_phi)*(1.0 + par->sin_phi)
          + 2.0*par->shear*(1.0 - par->sin_phi)*(1.0 - par->sin_phi);

  denom_l = 4.0*par->lame*par->sin_phi*par->sin_phi
          + par->shear*(1.0 + par->sin_phi)*(1.0 + par->sin_phi)
          + 2.0*par->shear*(1.0 - par->sin_phi)*(1.0 - par->sin_phi);

  denom_r = 4.0*par->lame*par->sin_phi*par->sin_phi
          + 2.0*par->shear*(1.0 + par->sin_phi)*(1.0 + par->sin_phi)
          + par->shear*(1.0 - par->sin_phi)*(1.0 - par->sin_phi);

  lambda_s = f_tr/denom_s;
  lambda_l = (par->shear*((1.0 + par->sin_phi)*(cache->eig[0] + cache->eig[1]) - 2.0*(1.0 - par->sin_phi)*cache->eig[2])
             + 2.0*par->lame*par->sin_phi*trace_E - par->c_bar)/denom_l;
  lambda_r = (par->shear*(2.0*(1.0 + par->sin_phi)*cache->eig[0] - (1.0 - par->sin_phi)*(cache->eig[1] + cache->eig[2]))
             + 2.0*par->lame*par->sin_phi*trace_E - par->c_bar)/denom_r;

  if (f_tr <= 0.0)
  {
    cache->return_type = MATMODEL_RETURN_ELASTIC;
    cache->sig_princ[0] = par->lame*trace_E + 2.0*par->shear*cache->eig[0];
    cache->sig_princ[1] = par->lame*trace_E + 2.0*par->shear*cache->eig[1];
    cache->sig_princ[2] = par->lame*trace_E + 2.0*par->shear*cache->eig[2];
  }
  else if (lambda_s <= ((gamma_sl < gamma_sr) ? gamma_sl : gamma_sr))
  {
    cache->return_type = MATMODEL_RETURN_SMOOTH;
    cache->sig_princ[0] = par->lame*trace_E + 2.0*par->shear*cache->eig[0]
                        - lambda_s*(2.0*par->lame*par->sin_phi + 2.0*par->shear*(1.0 + par->sin_phi));
    cache->sig_princ[1] = par->lame*trace_E + 2.0*par->shear*cache->eig[1]
                        - lambda_s*(2.0*par->lame*par->sin_phi);
    cache->sig_princ[2] = par->lame*trace_E + 2.0*par->shear*cache->eig[2]
                        - lambda_s*(2.0*par->lame*par->sin_phi - 2.0*par->shear*(1.0 - par->sin_phi));
  }
  else if ((gamma_sl < gamma_sr) && (lambda_l >= gamma_sl) && (lambda_l <= gamma_la))
  {
    cache->return_type = MATMODEL_RETURN_LEFT_EDGE;
    cache->sig_princ[0] = par->lame*trace_E + par->shear*(cache->eig[0] + cache->eig[1])
                        - lambda_l*(2.0*par->lame*par->sin_phi + par->shear*(1.0 + par->sin_phi));
    cache->sig_princ[1] = cache->sig_princ[0];
    cache->sig_princ[2] = par->lame*trace_E + 2.0*par->shear*cache->eig[2]
                        - lambda_l*(2.0*par->lame*par->sin_phi - 2.0*par->shear*(1.0 - par->sin_phi));
  }
  else if ((gamma_sl > gamma_sr) && (lambda_r >= gamma_sr) && (lambda_r <= gamma_ra))
  {
    cache->return_type = MATMODEL_RETURN_RIGHT_EDGE;
    cache->sig_princ[0] = par->lame*trace_E + 2.0*par->shear*cache->eig[0]
                        - lambda_r*(2.0*par->lame*par->sin_phi + 2.0*par->shear*(1.0 + par->sin_phi));
    cache->sig_princ[2] = par->lame*trace_E + par->shear*(cache->eig[1] + cache->eig[2])
                        - lambda_r*(2.0*par->lame*par->sin_phi - par->shear*(1.0 - par->sin_phi));
    cache->sig_princ[1] = cache->sig_princ[2];
  }
  else
  {
    cache->return_type = MATMODEL_RETURN_APEX;
    cache->sig_princ[0] = par->c_bar/(2.0*par->sin_phi);
    cache->sig_princ[1] = cache->sig_princ[0];
    cache->sig_princ[2] = cache->sig_princ[0];
  }

  vector stress_v;
  vector proj_v0;
  vector proj_v1;
  vector proj_v2;

  stress_v.makerefv(MATMODEL_NCOMP_STRESS, cache->stress);
  proj_v0.makerefv(MATMODEL_NCOMP_STRAIN, cache->proj[0]);
  proj_v1.makerefv(MATMODEL_NCOMP_STRAIN, cache->proj[1]);
  proj_v2.makerefv(MATMODEL_NCOMP_STRAIN, cache->proj[2]);

  nullv(stress_v);
  addmultv(stress_v, proj_v0, cache->sig_princ[0]);
  addmultv(stress_v, proj_v1, cache->sig_princ[1]);
  addmultv(stress_v, proj_v2, cache->sig_princ[2]);

  compliance_times_stress(par, cache->stress, eps_el);
  for (long i=0; i<4; i++)
    cache->epsp[i] = strain[i] - eps_el[i];
}

void compute_tangent(const matmodel_params *par,
                     const matmodel_cache *cache,
                     matrix &d)
{
  double denom_s, denom_l, denom_r;
  vector iota(4), aux(4), proj0(4), proj1(4), proj2(4), proj12(4), proj23(4);
  matrix zero(4,4), outer(4,4), hess0(4,4), hess1(4,4), hess2(4,4), hess12(4,4), hess23(4,4);

  reallocm(4, 4, d);
  nullm(d);
  if ((par == NULL) || (cache == NULL))
    return;

  nullm(zero);
  fillv(0.0, iota);
  iota[0] = 1.0;
  iota[1] = 1.0;
  iota[3] = 1.0;

  copyv(cache->proj[0], proj0);
  copyv(cache->proj[1], proj1);
  copyv(cache->proj[2], proj2);
  copym(&cache->hess[0][0][0], hess0);
  copym(&cache->hess[1][0][0], hess1);
  copym(&cache->hess[2][0][0], hess2);

  denom_s = 4.0*par->lame*par->sin_phi*par->sin_phi
          + 2.0*par->shear*(1.0 + par->sin_phi)*(1.0 + par->sin_phi)
          + 2.0*par->shear*(1.0 - par->sin_phi)*(1.0 - par->sin_phi);

  denom_l = 4.0*par->lame*par->sin_phi*par->sin_phi
          + par->shear*(1.0 + par->sin_phi)*(1.0 + par->sin_phi)
          + 2.0*par->shear*(1.0 - par->sin_phi)*(1.0 - par->sin_phi);

  denom_r = 4.0*par->lame*par->sin_phi*par->sin_phi
          + 2.0*par->shear*(1.0 + par->sin_phi)*(1.0 + par->sin_phi)
          + par->shear*(1.0 - par->sin_phi)*(1.0 - par->sin_phi);

  switch (cache->return_type)
  {
    case MATMODEL_RETURN_ELASTIC:
      elastic_matrix(par, d);
      break;

    case MATMODEL_RETURN_SMOOTH:
      addmultm(zero, 0.0, hess0, cache->sig_princ[0], d);
      addmultm(zero, 0.0, hess1, cache->sig_princ[1], d);
      addmultm(zero, 0.0, hess2, cache->sig_princ[2], d);
      vxv(iota, iota, outer);
      addmultm(zero, 0.0, outer, par->lame, d);
      vxv(proj0, proj0, outer);
      addmultm(zero, 0.0, outer, 2.0*par->shear, d);
      vxv(proj1, proj1, outer);
      addmultm(zero, 0.0, outer, 2.0*par->shear, d);
      vxv(proj2, proj2, outer);
      addmultm(zero, 0.0, outer, 2.0*par->shear, d);
      addmultv(proj0, 2.0*par->shear*(1.0 + par->sin_phi),
               proj2, -2.0*par->shear*(1.0 - par->sin_phi), aux);
      addmultv(aux, iota, 2.0*par->lame*par->sin_phi);
      vxv(aux, aux, outer);
      addmultm(zero, 0.0, outer, -1.0/denom_s, d);
      break;

    case MATMODEL_RETURN_LEFT_EDGE:
      addmultv(proj0, 1.0, proj1, 1.0, proj12);
      addm(hess0, hess1, hess12);
      addmultm(zero, 0.0, hess12, cache->sig_princ[0], d);
      addmultm(zero, 0.0, hess2, cache->sig_princ[2], d);
      vxv(iota, iota, outer);
      addmultm(zero, 0.0, outer, par->lame, d);
      vxv(proj12, proj12, outer);
      addmultm(zero, 0.0, outer, par->shear, d);
      vxv(proj2, proj2, outer);
      addmultm(zero, 0.0, outer, 2.0*par->shear, d);
      addmultv(proj12, par->shear*(1.0 + par->sin_phi),
               proj2, -2.0*par->shear*(1.0 - par->sin_phi), aux);
      addmultv(aux, iota, 2.0*par->lame*par->sin_phi);
      vxv(aux, aux, outer);
      addmultm(zero, 0.0, outer, -1.0/denom_l, d);
      break;

    case MATMODEL_RETURN_RIGHT_EDGE:
      addmultv(proj1, 1.0, proj2, 1.0, proj23);
      addm(hess1, hess2, hess23);
      addmultm(zero, 0.0, hess0, cache->sig_princ[0], d);
      addmultm(zero, 0.0, hess23, cache->sig_princ[2], d);
      vxv(iota, iota, outer);
      addmultm(zero, 0.0, outer, par->lame, d);
      vxv(proj0, proj0, outer);
      addmultm(zero, 0.0, outer, 2.0*par->shear, d);
      vxv(proj23, proj23, outer);
      addmultm(zero, 0.0, outer, par->shear, d);
      addmultv(proj0, 2.0*par->shear*(1.0 + par->sin_phi),
               proj23, -par->shear*(1.0 - par->sin_phi), aux);
      addmultv(aux, iota, 2.0*par->lame*par->sin_phi);
      vxv(aux, aux, outer);
      addmultm(zero, 0.0, outer, -1.0/denom_r, d);
      break;

    case MATMODEL_RETURN_APEX:
    default:
      break;
  }
}

void pack_other(const matmodel_cache *cache, double other[MATMODEL_NCOMP_OTHER])
{
  long idx;

  if ((cache == NULL) || (other == NULL))
    return;

  for (long i=0; i<MATMODEL_NCOMP_OTHER; i++)
    other[i] = 0.0;

  other[MATMODEL_IO_EP_XX] = cache->epsp[0];
  other[MATMODEL_IO_EP_YY] = cache->epsp[1];
  other[MATMODEL_IO_GP_XY] = cache->epsp[2];
  other[MATMODEL_IO_EP_ZZ] = cache->epsp[3];
  other[MATMODEL_IO_RETURN_TYPE] = static_cast<double>(cache->return_type);
  other[MATMODEL_IO_EIG_1] = cache->eig[0];
  other[MATMODEL_IO_EIG_2] = cache->eig[1];
  other[MATMODEL_IO_EIG_3] = cache->eig[2];

  for (long i=0; i<4; i++)
  {
    other[MATMODEL_IO_PROJ_1 + i] = cache->proj[0][i];
    other[MATMODEL_IO_PROJ_2 + i] = cache->proj[1][i];
    other[MATMODEL_IO_PROJ_3 + i] = cache->proj[2][i];
  }

  idx = MATMODEL_IO_HESS_1;
  for (long k=0; k<3; k++)
  {
    for (long j=0; j<4; j++)
    {
      for (long i=0; i<4; i++)
        other[idx++] = cache->hess[k][i][j];
    }
  }

  other[MATMODEL_IO_SIGMA_1] = cache->sig_princ[0];
  other[MATMODEL_IO_SIGMA_2] = cache->sig_princ[1];
  other[MATMODEL_IO_SIGMA_3] = cache->sig_princ[2];
}

void unpack_other(const double other[MATMODEL_NCOMP_OTHER], matmodel_cache *cache)
{
  long idx;

  if ((cache == NULL) || (other == NULL))
    return;

  memset(cache, 0, sizeof(*cache));

  cache->epsp[0] = other[MATMODEL_IO_EP_XX];
  cache->epsp[1] = other[MATMODEL_IO_EP_YY];
  cache->epsp[2] = other[MATMODEL_IO_GP_XY];
  cache->epsp[3] = other[MATMODEL_IO_EP_ZZ];
  cache->return_type = static_cast<int>(other[MATMODEL_IO_RETURN_TYPE] + 0.5);
  cache->eig[0] = other[MATMODEL_IO_EIG_1];
  cache->eig[1] = other[MATMODEL_IO_EIG_2];
  cache->eig[2] = other[MATMODEL_IO_EIG_3];

  for (long i=0; i<4; i++)
  {
    cache->proj[0][i] = other[MATMODEL_IO_PROJ_1 + i];
    cache->proj[1][i] = other[MATMODEL_IO_PROJ_2 + i];
    cache->proj[2][i] = other[MATMODEL_IO_PROJ_3 + i];
  }

  idx = MATMODEL_IO_HESS_1;
  for (long k=0; k<3; k++)
  {
    for (long j=0; j<4; j++)
    {
      for (long i=0; i<4; i++)
        cache->hess[k][i][j] = other[idx++];
    }
  }

  cache->sig_princ[0] = other[MATMODEL_IO_SIGMA_1];
  cache->sig_princ[1] = other[MATMODEL_IO_SIGMA_2];
  cache->sig_princ[2] = other[MATMODEL_IO_SIGMA_3];
}

} // namespace

matmodel::matmodel()
{
  init_params(&par, 0.0, 0.25, 0.0, 0.52359877559829887308);
  plane_strain_only = 1;
  memset(cached_strain, 0, sizeof(cached_strain));
  memset(cached_eqother, 0, sizeof(cached_eqother));
  memset(cached_stress, 0, sizeof(cached_stress));
  memset(cached_other, 0, sizeof(cached_other));
  cached_response_valid = 0;
}

long matmodel::read(FILE *in)
{
  double young, poisson, cohesion, phi;

  if (in == NULL)
    return 1;

  if (fscanf(in, "%lf %lf %lf %lf", &young, &poisson, &cohesion, &phi) != 4)
    return 1;

  init_params(&par, young, poisson, cohesion, phi);
  return (validate_params(&par) ? 0 : 1);
}

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
  fprintf(out, "  buffer = epsp(4), return_type(1), eig(3), proj(12), hess(48), sigma_princ(3)\n");
}

void matmodel::cache_response(const double strain[4],
                              const double eqother[4],
                              const double stress[4],
                              const double statev[MATMODEL_NCOMP_OTHER])
{
  copyv(strain, cached_strain, MATMODEL_NCOMP_STRAIN);
  copyv(eqother, cached_eqother, MATMODEL_NCOMP_EQOTHER);
  copyv(stress, cached_stress, MATMODEL_NCOMP_STRESS);
  copyv(statev, cached_other, MATMODEL_NCOMP_OTHER);
  cached_response_valid = 1;
}

int matmodel::cached_response_matches(const vector &strain,
                                      const vector &eqstatev,
                                      const vector &stress) const
{
  const double tol = 1.0e-14;
  double value = 0.0;

  if (!cached_response_valid)
    return 0;

  for (long i=0; i<MATMODEL_NCOMP_STRAIN; i++)
  {
    value = (i < strain.n) ? strain[i] : 0.0;
    if (fabs(value - cached_strain[i]) > tol)
      return 0;
  }

  for (long i=0; i<MATMODEL_NCOMP_EQOTHER; i++)
  {
    value = (i < eqstatev.n) ? eqstatev[i] : 0.0;
    if (fabs(value - cached_eqother[i]) > tol)
      return 0;
  }

  for (long i=0; i<MATMODEL_NCOMP_STRESS; i++)
  {
    value = (i < stress.n) ? stress[i] : 0.0;
    if (fabs(value - cached_stress[i]) > tol)
      return 0;
  }

  return 1;
}

void matmodel::fill_response(const vector &strain,
                             const vector &eqstatev,
                             vector &stress,
                             vector &statev) const
{
  const long nstrain = (strain.n < MATMODEL_NCOMP_STRAIN) ? strain.n : static_cast<long>(MATMODEL_NCOMP_STRAIN);
  const long neqother = (eqstatev.n < MATMODEL_NCOMP_EQOTHER) ? eqstatev.n : static_cast<long>(MATMODEL_NCOMP_EQOTHER);
  double e[4] = {0.0, 0.0, 0.0, 0.0};
  double eqp[4] = {0.0, 0.0, 0.0, 0.0};
  double other[MATMODEL_NCOMP_OTHER];
  matmodel_cache cache;
  vector e_v;
  vector eqp_v;

  reallocv(MATMODEL_NCOMP_STRESS, stress);
  reallocv(MATMODEL_NCOMP_OTHER, statev);
  e_v.makerefv(MATMODEL_NCOMP_STRAIN, e);
  eqp_v.makerefv(MATMODEL_NCOMP_EQOTHER, eqp);

  if (nstrain > 0)
    rcopyv(strain, 0, e_v, 0, nstrain);
  if (neqother > 0)
    rcopyv(eqstatev, 0, eqp_v, 0, neqother);

  compute_response(&par, e, eqp, &cache);
  pack_other(&cache, other);

  copyv(cache.stress, stress);
  copyv(other, statev);

  const_cast<matmodel *>(this)->cache_response(e, eqp, cache.stress, other);
}

void matmodel::nlstresses(const vector &strain,
                          const vector &eqstatev,
                          vector &stress,
                          vector &statev)
{
  fill_response(strain, eqstatev, stress, statev);
}

void matmodel::stiffmat_from_statev(const vector &statev, matrix &d)
{
  const long nstatev = (statev.n < MATMODEL_NCOMP_OTHER) ? statev.n : static_cast<long>(MATMODEL_NCOMP_OTHER);
  double other[MATMODEL_NCOMP_OTHER];
  matmodel_cache cache;
  vector other_v;

  nullv(other, MATMODEL_NCOMP_OTHER);
  other_v.makerefv(MATMODEL_NCOMP_OTHER, other);
  if (nstatev > 0)
    rcopyv(statev, 0, other_v, 0, nstatev);

  unpack_other(other, &cache);
  compute_tangent(&par, &cache, d);
}

void matmodel::stiffmat(const vector &strain,
                        const vector &eqstatev,
                        const vector &stress,
                        matrix &d)
{
  vector statev;
  vector hstress;

  if (cached_response_matches(strain, eqstatev, stress))
  {
    reallocv(MATMODEL_NCOMP_OTHER, statev);
    copyv(cached_other, statev);
  }
  else
  {
    fill_response(strain, eqstatev, hstress, statev);
  }

  stiffmat_from_statev(statev, d);
}

void matmodel::updateval(const vector &statev, vector &eqstatev)
{
  const long nstatev = (statev.n < MATMODEL_NCOMP_EQOTHER) ? statev.n : static_cast<long>(MATMODEL_NCOMP_EQOTHER);
  reallocv(MATMODEL_NCOMP_EQOTHER, eqstatev);
  nullv(eqstatev);
  if (nstatev > 0)
    rcopyv(statev, 0, eqstatev, 0, nstatev);
}
