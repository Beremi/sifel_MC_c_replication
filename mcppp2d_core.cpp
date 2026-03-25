#include "mcppp2d_core.h"

#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void mcppp2d_zero4(double v[4])
{
  long i;
  for (i=0; i<4; i++)
    v[i] = 0.0;
}

static void mcppp2d_zero44(double a[4][4])
{
  long i, j;
  for (i=0; i<4; i++)
    for (j=0; j<4; j++)
      a[i][j] = 0.0;
}

static void mcppp2d_copy4(const double src[4], double dst[4])
{
  long i;
  for (i=0; i<4; i++)
    dst[i] = src[i];
}

static void mcppp2d_copy44(const double src[4][4], double dst[4][4])
{
  long i, j;
  for (i=0; i<4; i++)
    for (j=0; j<4; j++)
      dst[i][j] = src[i][j];
}


static void mcppp2d_add44(double dst[4][4], const double src[4][4], double scale)
{
  long i, j;
  for (i=0; i<4; i++)
    for (j=0; j<4; j++)
      dst[i][j] += scale*src[i][j];
}

static void mcppp2d_add_outer4(double dst[4][4], const double a[4], const double b[4], double scale)
{
  long i, j;
  for (i=0; i<4; i++)
    for (j=0; j<4; j++)
      dst[i][j] += scale*a[i]*b[j];
}

static void mcppp2d_sub_outer4(double dst[4][4], const double a[4], const double b[4], double scale)
{
  long i, j;
  for (i=0; i<4; i++)
    for (j=0; j<4; j++)
      dst[i][j] -= scale*a[i]*b[j];
}

static void mcppp2d_add4(double dst[4], const double src[4], double scale)
{
  long i;
  for (i=0; i<4; i++)
    dst[i] += scale*src[i];
}

static void mcppp2d_compliance_times_stress(const mcppp2d_params *par,
                                             const double stress[4],
                                             double eps_el[4])
{
  eps_el[0] = (stress[0] - par->poisson*(stress[1] + stress[3]))/par->young;
  eps_el[1] = (stress[1] - par->poisson*(stress[0] + stress[3]))/par->young;
  eps_el[2] = stress[2]/par->shear;
  eps_el[3] = (stress[3] - par->poisson*(stress[0] + stress[1]))/par->young;
}

static void mcppp2d_compute_inplane_spectral(const double trial[4],
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
  const double g  = trial[2];
  const double d  = ex - ey;
  const double i1 = ex + ey;
  const double i2 = sqrt(d*d + g*g);
  const double tol = 1.0e-14;

  long a, b;

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
    proj[0][3] = 0.0;

    proj[1][0] = 1.0 - proj[0][0];
    proj[1][1] = 1.0 - proj[0][1];
    proj[1][2] =     - proj[0][2];
    proj[1][3] = 0.0;

    for (a=0; a<4; a++)
    {
      for (b=0; b<4; b++)
      {
        hess[0][a][b] = (metric[a][b] - proj[0][a]*proj[0][b] - proj[1][a]*proj[1][b])/i2;
        hess[1][a][b] = -hess[0][a][b];
      }
    }

    for (a=0; a<4; a++)
    {
      hess[0][3][a] = 0.0;
      hess[0][a][3] = 0.0;
      hess[1][3][a] = 0.0;
      hess[1][a][3] = 0.0;
    }
  }

  proj[2][0] = 0.0;
  proj[2][1] = 0.0;
  proj[2][2] = 0.0;
  proj[2][3] = 1.0;
}

static void mcppp2d_copy_mode(double dst_proj[4],
                              double dst_hess[4][4],
                              double *dst_eig,
                              const double src_proj[4],
                              const double src_hess[4][4],
                              double src_eig)
{
  mcppp2d_copy4(src_proj, dst_proj);
  mcppp2d_copy44(src_hess, dst_hess);
  *dst_eig = src_eig;
}

void mcppp2d_init_params(mcppp2d_params *par,
                         double young,
                         double poisson,
                         double cohesion,
                         double phi)
{
  if (par == NULL)
    return;

  par->young    = young;
  par->poisson  = poisson;
  par->c        = cohesion;
  par->phi      = phi;
  par->shear    = young/(2.0*(1.0 + poisson));
  par->bulk     = young/(3.0*(1.0 - 2.0*poisson));
  par->lame     = par->bulk - 2.0*par->shear/3.0;
  par->sin_phi  = sin(phi);
  par->cos_phi  = cos(phi);
  par->c_bar    = 2.0*cohesion*par->cos_phi;
}

int mcppp2d_validate_params(const mcppp2d_params *par)
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

const char *mcppp2d_return_name(int return_type)
{
  switch (return_type)
  {
    case MCPPP2D_RETURN_ELASTIC:    return "elastic";
    case MCPPP2D_RETURN_SMOOTH:     return "smooth";
    case MCPPP2D_RETURN_LEFT_EDGE:  return "left_edge";
    case MCPPP2D_RETURN_RIGHT_EDGE: return "right_edge";
    case MCPPP2D_RETURN_APEX:       return "apex";
    default:                        return "unknown";
  }
}

void mcppp2d_elastic_matrix(const mcppp2d_params *par, double d[4][4])
{
  mcppp2d_zero44(d);
  if (par == NULL)
    return;

  d[0][0] = par->lame + 2.0*par->shear;
  d[0][1] = par->lame;
  d[0][3] = par->lame;

  d[1][0] = par->lame;
  d[1][1] = par->lame + 2.0*par->shear;
  d[1][3] = par->lame;

  d[2][2] = par->shear;

  d[3][0] = par->lame;
  d[3][1] = par->lame;
  d[3][3] = par->lame + 2.0*par->shear;
}

void mcppp2d_compute_response(const mcppp2d_params *par,
                              const double strain[4],
                              const double eqother[4],
                              mcppp2d_cache *cache)
{
  double trial[4], epsp_prev[4], eps_el[4];
  double eig0[3], proj0[3][4], hess0[3][4][4];
  double trace_E, f_tr, gamma_sl, gamma_sr, gamma_la, gamma_ra;
  double denom_s, denom_l, denom_r;
  double lambda_s, lambda_l, lambda_r;
  long i;

  if ((par == NULL) || (cache == NULL) || (strain == NULL))
    return;

  memset(cache, 0, sizeof(*cache));
  mcppp2d_zero4(epsp_prev);
  if (eqother != NULL)
    mcppp2d_copy4(eqother, epsp_prev);

  for (i=0; i<4; i++)
    trial[i] = strain[i] - epsp_prev[i];

  mcppp2d_compute_inplane_spectral(trial, eig0, proj0, hess0);

  mcppp2d_copy_mode(cache->proj[0], cache->hess[0], &cache->eig[0], proj0[0], hess0[0], eig0[0]);
  mcppp2d_copy_mode(cache->proj[1], cache->hess[1], &cache->eig[1], proj0[1], hess0[1], eig0[1]);
  mcppp2d_copy_mode(cache->proj[2], cache->hess[2], &cache->eig[2], proj0[2], hess0[2], eig0[2]);

  /* exact Matlab ordering */
  if ((eig0[0] >= eig0[2]) && (eig0[2] > eig0[1]))
  {
    mcppp2d_copy_mode(cache->proj[1], cache->hess[1], &cache->eig[1], proj0[2], hess0[2], eig0[2]);
    mcppp2d_copy_mode(cache->proj[2], cache->hess[2], &cache->eig[2], proj0[1], hess0[1], eig0[1]);
  }
  if (eig0[2] > eig0[0])
  {
    mcppp2d_copy_mode(cache->proj[0], cache->hess[0], &cache->eig[0], proj0[2], hess0[2], eig0[2]);
    mcppp2d_copy_mode(cache->proj[1], cache->hess[1], &cache->eig[1], proj0[0], hess0[0], eig0[0]);
    mcppp2d_copy_mode(cache->proj[2], cache->hess[2], &cache->eig[2], proj0[1], hess0[1], eig0[1]);
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
          +       par->shear*(1.0 + par->sin_phi)*(1.0 + par->sin_phi)
          + 2.0*par->shear*(1.0 - par->sin_phi)*(1.0 - par->sin_phi);

  denom_r = 4.0*par->lame*par->sin_phi*par->sin_phi
          + 2.0*par->shear*(1.0 + par->sin_phi)*(1.0 + par->sin_phi)
          +       par->shear*(1.0 - par->sin_phi)*(1.0 - par->sin_phi);

  lambda_s = f_tr/denom_s;
  lambda_l = (par->shear*((1.0 + par->sin_phi)*(cache->eig[0] + cache->eig[1]) - 2.0*(1.0 - par->sin_phi)*cache->eig[2])
             + 2.0*par->lame*par->sin_phi*trace_E - par->c_bar)/denom_l;
  lambda_r = (par->shear*(2.0*(1.0 + par->sin_phi)*cache->eig[0] - (1.0 - par->sin_phi)*(cache->eig[1] + cache->eig[2]))
             + 2.0*par->lame*par->sin_phi*trace_E - par->c_bar)/denom_r;

  if (f_tr <= 0.0)
  {
    cache->return_type = MCPPP2D_RETURN_ELASTIC;
    cache->sig_princ[0] = par->lame*trace_E + 2.0*par->shear*cache->eig[0];
    cache->sig_princ[1] = par->lame*trace_E + 2.0*par->shear*cache->eig[1];
    cache->sig_princ[2] = par->lame*trace_E + 2.0*par->shear*cache->eig[2];
  }
  else if (lambda_s <= ((gamma_sl < gamma_sr) ? gamma_sl : gamma_sr))
  {
    cache->return_type = MCPPP2D_RETURN_SMOOTH;
    cache->sig_princ[0] = par->lame*trace_E + 2.0*par->shear*cache->eig[0]
                        - lambda_s*(2.0*par->lame*par->sin_phi + 2.0*par->shear*(1.0 + par->sin_phi));
    cache->sig_princ[1] = par->lame*trace_E + 2.0*par->shear*cache->eig[1]
                        - lambda_s*(2.0*par->lame*par->sin_phi);
    cache->sig_princ[2] = par->lame*trace_E + 2.0*par->shear*cache->eig[2]
                        - lambda_s*(2.0*par->lame*par->sin_phi - 2.0*par->shear*(1.0 - par->sin_phi));
  }
  else if ((gamma_sl < gamma_sr) && (lambda_l >= gamma_sl) && (lambda_l <= gamma_la))
  {
    cache->return_type = MCPPP2D_RETURN_LEFT_EDGE;
    cache->sig_princ[0] = par->lame*trace_E + par->shear*(cache->eig[0] + cache->eig[1])
                        - lambda_l*(2.0*par->lame*par->sin_phi + par->shear*(1.0 + par->sin_phi));
    cache->sig_princ[1] = cache->sig_princ[0];
    cache->sig_princ[2] = par->lame*trace_E + 2.0*par->shear*cache->eig[2]
                        - lambda_l*(2.0*par->lame*par->sin_phi - 2.0*par->shear*(1.0 - par->sin_phi));
  }
  else if ((gamma_sl > gamma_sr) && (lambda_r >= gamma_sr) && (lambda_r <= gamma_ra))
  {
    cache->return_type = MCPPP2D_RETURN_RIGHT_EDGE;
    cache->sig_princ[0] = par->lame*trace_E + 2.0*par->shear*cache->eig[0]
                        - lambda_r*(2.0*par->lame*par->sin_phi + 2.0*par->shear*(1.0 + par->sin_phi));
    cache->sig_princ[2] = par->lame*trace_E + par->shear*(cache->eig[1] + cache->eig[2])
                        - lambda_r*(2.0*par->lame*par->sin_phi - par->shear*(1.0 - par->sin_phi));
    cache->sig_princ[1] = cache->sig_princ[2];
  }
  else
  {
    cache->return_type = MCPPP2D_RETURN_APEX;
    cache->sig_princ[0] = par->c_bar/(2.0*par->sin_phi);
    cache->sig_princ[1] = cache->sig_princ[0];
    cache->sig_princ[2] = cache->sig_princ[0];
  }

  mcppp2d_zero4(cache->stress);
  mcppp2d_add4(cache->stress, cache->proj[0], cache->sig_princ[0]);
  mcppp2d_add4(cache->stress, cache->proj[1], cache->sig_princ[1]);
  mcppp2d_add4(cache->stress, cache->proj[2], cache->sig_princ[2]);

  mcppp2d_compliance_times_stress(par, cache->stress, eps_el);
  for (i=0; i<4; i++)
    cache->epsp[i] = strain[i] - eps_el[i];
}

void mcppp2d_compute_tangent(const mcppp2d_params *par,
                             const mcppp2d_cache *cache,
                             double d[4][4])
{
  static const double iota[4] = {1.0, 1.0, 0.0, 1.0};
  double denom_s, denom_l, denom_r;
  double aux[4], proj12[4], proj23[4];
  double hess12[4][4], hess23[4][4];
  long i, j;

  mcppp2d_zero44(d);
  if ((par == NULL) || (cache == NULL))
    return;

  denom_s = 4.0*par->lame*par->sin_phi*par->sin_phi
          + 2.0*par->shear*(1.0 + par->sin_phi)*(1.0 + par->sin_phi)
          + 2.0*par->shear*(1.0 - par->sin_phi)*(1.0 - par->sin_phi);

  denom_l = 4.0*par->lame*par->sin_phi*par->sin_phi
          +       par->shear*(1.0 + par->sin_phi)*(1.0 + par->sin_phi)
          + 2.0*par->shear*(1.0 - par->sin_phi)*(1.0 - par->sin_phi);

  denom_r = 4.0*par->lame*par->sin_phi*par->sin_phi
          + 2.0*par->shear*(1.0 + par->sin_phi)*(1.0 + par->sin_phi)
          +       par->shear*(1.0 - par->sin_phi)*(1.0 - par->sin_phi);

  switch (cache->return_type)
  {
    case MCPPP2D_RETURN_ELASTIC:
      mcppp2d_elastic_matrix(par, d);
      break;

    case MCPPP2D_RETURN_SMOOTH:
      mcppp2d_add44(d, cache->hess[0], cache->sig_princ[0]);
      mcppp2d_add44(d, cache->hess[1], cache->sig_princ[1]);
      mcppp2d_add44(d, cache->hess[2], cache->sig_princ[2]);
      mcppp2d_add_outer4(d, iota, iota, par->lame);
      mcppp2d_add_outer4(d, cache->proj[0], cache->proj[0], 2.0*par->shear);
      mcppp2d_add_outer4(d, cache->proj[1], cache->proj[1], 2.0*par->shear);
      mcppp2d_add_outer4(d, cache->proj[2], cache->proj[2], 2.0*par->shear);
      for (i=0; i<4; i++)
        aux[i] = 2.0*par->shear*((1.0 + par->sin_phi)*cache->proj[0][i] - (1.0 - par->sin_phi)*cache->proj[2][i])
               + 2.0*par->lame*par->sin_phi*iota[i];
      mcppp2d_sub_outer4(d, aux, aux, 1.0/denom_s);
      break;

    case MCPPP2D_RETURN_LEFT_EDGE:
      for (i=0; i<4; i++)
        proj12[i] = cache->proj[0][i] + cache->proj[1][i];
      for (i=0; i<4; i++)
        for (j=0; j<4; j++)
          hess12[i][j] = cache->hess[0][i][j] + cache->hess[1][i][j];
      mcppp2d_add44(d, hess12,         cache->sig_princ[0]);
      mcppp2d_add44(d, cache->hess[2], cache->sig_princ[2]);
      mcppp2d_add_outer4(d, iota, iota, par->lame);
      mcppp2d_add_outer4(d, proj12, proj12, par->shear);
      mcppp2d_add_outer4(d, cache->proj[2], cache->proj[2], 2.0*par->shear);
      for (i=0; i<4; i++)
        aux[i] = par->shear*((1.0 + par->sin_phi)*proj12[i] - 2.0*(1.0 - par->sin_phi)*cache->proj[2][i])
               + 2.0*par->lame*par->sin_phi*iota[i];
      mcppp2d_sub_outer4(d, aux, aux, 1.0/denom_l);
      break;

    case MCPPP2D_RETURN_RIGHT_EDGE:
      for (i=0; i<4; i++)
        proj23[i] = cache->proj[1][i] + cache->proj[2][i];
      for (i=0; i<4; i++)
        for (j=0; j<4; j++)
          hess23[i][j] = cache->hess[1][i][j] + cache->hess[2][i][j];
      mcppp2d_add44(d, cache->hess[0], cache->sig_princ[0]);
      mcppp2d_add44(d, hess23,         cache->sig_princ[2]);
      mcppp2d_add_outer4(d, iota, iota, par->lame);
      mcppp2d_add_outer4(d, cache->proj[0], cache->proj[0], 2.0*par->shear);
      mcppp2d_add_outer4(d, proj23, proj23, par->shear);
      for (i=0; i<4; i++)
        aux[i] = par->shear*(2.0*(1.0 + par->sin_phi)*cache->proj[0][i] - (1.0 - par->sin_phi)*proj23[i])
               + 2.0*par->lame*par->sin_phi*iota[i];
      mcppp2d_sub_outer4(d, aux, aux, 1.0/denom_r);
      break;

    case MCPPP2D_RETURN_APEX:
    default:
      /* zero tangent */
      break;
  }
}

void mcppp2d_pack_eqother(const mcppp2d_cache *cache, double eqother[4])
{
  if ((cache == NULL) || (eqother == NULL))
    return;
  mcppp2d_copy4(cache->epsp, eqother);
}

void mcppp2d_pack_other(const mcppp2d_cache *cache, double other[71])
{
  long i, j, k, idx;

  if ((cache == NULL) || (other == NULL))
    return;

  for (i=0; i<71; i++)
    other[i] = 0.0;

  other[MCPPP2D_IO_EP_XX] = cache->epsp[0];
  other[MCPPP2D_IO_EP_YY] = cache->epsp[1];
  other[MCPPP2D_IO_GP_XY] = cache->epsp[2];
  other[MCPPP2D_IO_EP_ZZ] = cache->epsp[3];

  other[MCPPP2D_IO_RETURN_TYPE] = (double)cache->return_type;

  other[MCPPP2D_IO_EIG_1] = cache->eig[0];
  other[MCPPP2D_IO_EIG_2] = cache->eig[1];
  other[MCPPP2D_IO_EIG_3] = cache->eig[2];

  for (i=0; i<4; i++)
  {
    other[MCPPP2D_IO_PROJ_1 + i] = cache->proj[0][i];
    other[MCPPP2D_IO_PROJ_2 + i] = cache->proj[1][i];
    other[MCPPP2D_IO_PROJ_3 + i] = cache->proj[2][i];
  }

  idx = MCPPP2D_IO_HESS_1;
  for (k=0; k<3; k++)
  {
    for (j=0; j<4; j++)
      for (i=0; i<4; i++)
        other[idx++] = cache->hess[k][i][j];
  }

  other[MCPPP2D_IO_SIGMA_1] = cache->sig_princ[0];
  other[MCPPP2D_IO_SIGMA_2] = cache->sig_princ[1];
  other[MCPPP2D_IO_SIGMA_3] = cache->sig_princ[2];
}

void mcppp2d_unpack_other(const double other[71], mcppp2d_cache *cache)
{
  long i, j, k, idx;

  if ((cache == NULL) || (other == NULL))
    return;

  memset(cache, 0, sizeof(*cache));

  cache->epsp[0] = other[MCPPP2D_IO_EP_XX];
  cache->epsp[1] = other[MCPPP2D_IO_EP_YY];
  cache->epsp[2] = other[MCPPP2D_IO_GP_XY];
  cache->epsp[3] = other[MCPPP2D_IO_EP_ZZ];

  cache->return_type = (int)(other[MCPPP2D_IO_RETURN_TYPE] + 0.5);

  cache->eig[0] = other[MCPPP2D_IO_EIG_1];
  cache->eig[1] = other[MCPPP2D_IO_EIG_2];
  cache->eig[2] = other[MCPPP2D_IO_EIG_3];

  for (i=0; i<4; i++)
  {
    cache->proj[0][i] = other[MCPPP2D_IO_PROJ_1 + i];
    cache->proj[1][i] = other[MCPPP2D_IO_PROJ_2 + i];
    cache->proj[2][i] = other[MCPPP2D_IO_PROJ_3 + i];
  }

  idx = MCPPP2D_IO_HESS_1;
  for (k=0; k<3; k++)
  {
    for (j=0; j<4; j++)
      for (i=0; i<4; i++)
        cache->hess[k][i][j] = other[idx++];
  }

  cache->sig_princ[0] = other[MCPPP2D_IO_SIGMA_1];
  cache->sig_princ[1] = other[MCPPP2D_IO_SIGMA_2];
  cache->sig_princ[2] = other[MCPPP2D_IO_SIGMA_3];
}
