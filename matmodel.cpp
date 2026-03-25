#include "matmodel.h"
#include "vector.h"
#include "matrix.h"

#include <stdio.h>

matmodel::matmodel()
{
  mcppp2d_init_params(&par, 0.0, 0.25, 0.0, 0.52359877559829887308);
  plane_strain_only = 1;
}

long matmodel::read(FILE *in)
{
  double young, poisson, cohesion, phi;

  if (in == NULL)
    return 1;

  /*
     Expected input order in this template implementation:
       young  poisson  cohesion  phi
     where phi is in radians.
  */
  if (fscanf(in, "%lf %lf %lf %lf", &young, &poisson, &cohesion, &phi) != 4)
    return 1;

  mcppp2d_init_params(&par, young, poisson, cohesion, phi);
  return (mcppp2d_validate_params(&par) ? 0 : 1);
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
  fprintf(out, "  other  = %d components\n", MCPPP2D_NCOMP_OTHER);
  fprintf(out, "  eqother= %d components\n", MCPPP2D_NCOMP_EQOTHER);
  fprintf(out, "  buffer = epsp(4), return_type(1), eig(3), proj(12), hess(48), sigma_princ(3)\n");
}

void matmodel::write_matrix_to_statev(const double a[4][4], long start, vector &statev) const
{
  long i, j, idx = start;
  for (j=0; j<4; j++)
  {
    for (i=0; i<4; i++)
      statev[idx++] = a[i][j];
  }
}

void matmodel::write_vector_to_statev(const double a[4], long start, vector &statev) const
{
  long i;
  for (i=0; i<4; i++)
    statev[start+i] = a[i];
}

void matmodel::write_tangent_to_matrix(const double a[4][4], matrix &d) const
{
  long i, j;
  reallocm(4, 4, d);
  for (i=0; i<4; i++)
    for (j=0; j<4; j++)
      d(i,j) = a[i][j];
}

void matmodel::fill_response(const vector &strain,
                             const vector &eqstatev,
                             vector &stress,
                             vector &statev) const
{
  double e[4]   = {0.0, 0.0, 0.0, 0.0};
  double eqp[4] = {0.0, 0.0, 0.0, 0.0};
  double other[71];
  mcppp2d_cache cache;
  long i;

  reallocv(4, stress);
  reallocv(MCPPP2D_NCOMP_OTHER, statev);

  for (i=0; (i<4) && (i<strain.n); i++)
    e[i] = strain[i];
  for (i=0; (i<4) && (i<eqstatev.n); i++)
    eqp[i] = eqstatev[i];

  mcppp2d_compute_response(&par, e, eqp, &cache);
  mcppp2d_pack_other(&cache, other);

  for (i=0; i<4; i++)
    stress[i] = cache.stress[i];
  for (i=0; i<MCPPP2D_NCOMP_OTHER; i++)
    statev[i] = other[i];
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
  double dd[4][4];
  double other[71];
  mcppp2d_cache cache;
  long i;

  for (i=0; i<MCPPP2D_NCOMP_OTHER; i++)
    other[i] = (i < statev.n) ? statev[i] : 0.0;

  mcppp2d_unpack_other(other, &cache);
  mcppp2d_compute_tangent(&par, &cache, dd);
  write_tangent_to_matrix(dd, d);
}

void matmodel::stiffmat(const vector &strain,
                        const vector &eqstatev,
                        const vector &stress,
                        matrix &d)
{
  vector statev;
  vector hstress;
  (void)stress;

  /*
     The supplied template does not pass the current statev/other buffer into
     stiffmat(). Therefore the cache is recomputed here. In the real SIFEL
     material class, matstiff() should read Mm->ip[ipp].other and call
     stiffmat_from_statev() instead.
  */
  fill_response(strain, eqstatev, hstress, statev);
  stiffmat_from_statev(statev, d);
}

void matmodel::updateval(const vector &statev, vector &eqstatev)
{
  long i;
  reallocv(MCPPP2D_NCOMP_EQOTHER, eqstatev);
  for (i=0; i<MCPPP2D_NCOMP_EQOTHER; i++)
    eqstatev[i] = (i < statev.n) ? statev[i] : 0.0;
}
