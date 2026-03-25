#ifndef MCPPP2D_CORE_H
#define MCPPP2D_CORE_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
  2D associative, perfectly plastic Mohr-Coulomb model in the 4-component
  Voigt notation used by the supplied SIFEL template:

    strain = [ eps_xx, eps_yy, gamma_xy, eps_zz ]
    stress = [ sig_xx, sig_yy, tau_xy,   sig_zz ]

  The model assumes a 3D isotropic elastic law with only the in-plane shear
  gamma_xy active. This is the natural constitutive setting for plane strain.
*/

enum
{
  MCPPP2D_NCOMP_STRAIN  = 4,
  MCPPP2D_NCOMP_STRESS  = 4,
  MCPPP2D_NCOMP_EQOTHER = 4,
  MCPPP2D_NCOMP_OTHER   = 71
};

enum mcppp2d_return_type
{
  MCPPP2D_RETURN_ELASTIC    = 0,
  MCPPP2D_RETURN_SMOOTH     = 1,
  MCPPP2D_RETURN_LEFT_EDGE  = 2,
  MCPPP2D_RETURN_RIGHT_EDGE = 3,
  MCPPP2D_RETURN_APEX       = 4
};

/*
  Layout of the SIFEL "other" buffer.

  Hessians are flattened in column-major order to stay close to the original
  Matlab vectorization, i.e.

    idx = row + 4*col

  where row,col are in <0,3> for the 4x4 Hessian block.
*/
enum mcppp2d_other_index
{
  MCPPP2D_IO_EP_XX = 0,
  MCPPP2D_IO_EP_YY = 1,
  MCPPP2D_IO_GP_XY = 2,
  MCPPP2D_IO_EP_ZZ = 3,

  MCPPP2D_IO_RETURN_TYPE = 4,

  MCPPP2D_IO_EIG_1 = 5,
  MCPPP2D_IO_EIG_2 = 6,
  MCPPP2D_IO_EIG_3 = 7,

  MCPPP2D_IO_PROJ_1 = 8,   /* 8..11  */
  MCPPP2D_IO_PROJ_2 = 12,  /* 12..15 */
  MCPPP2D_IO_PROJ_3 = 16,  /* 16..19 */

  MCPPP2D_IO_HESS_1 = 20,  /* 20..35 */
  MCPPP2D_IO_HESS_2 = 36,  /* 36..51 */
  MCPPP2D_IO_HESS_3 = 52,  /* 52..67 */

  MCPPP2D_IO_SIGMA_1 = 68,
  MCPPP2D_IO_SIGMA_2 = 69,
  MCPPP2D_IO_SIGMA_3 = 70
};

typedef struct mcppp2d_params
{
  double young;
  double poisson;
  double c;
  double phi;      /* radians */

  double shear;
  double bulk;
  double lame;
  double sin_phi;
  double cos_phi;
  double c_bar;
} mcppp2d_params;

typedef struct mcppp2d_cache
{
  double epsp[MCPPP2D_NCOMP_EQOTHER];
  int    return_type;

  double eig[3];
  double proj[3][MCPPP2D_NCOMP_STRAIN];
  double hess[3][MCPPP2D_NCOMP_STRAIN][MCPPP2D_NCOMP_STRAIN];
  double sig_princ[3];
  double stress[MCPPP2D_NCOMP_STRESS];
} mcppp2d_cache;

void mcppp2d_init_params(mcppp2d_params *par,
                         double young,
                         double poisson,
                         double cohesion,
                         double phi);

int mcppp2d_validate_params(const mcppp2d_params *par);

const char *mcppp2d_return_name(int return_type);

void mcppp2d_elastic_matrix(const mcppp2d_params *par,
                            double d[MCPPP2D_NCOMP_STRAIN][MCPPP2D_NCOMP_STRAIN]);

void mcppp2d_compute_response(const mcppp2d_params *par,
                              const double strain[MCPPP2D_NCOMP_STRAIN],
                              const double eqother[MCPPP2D_NCOMP_EQOTHER],
                              mcppp2d_cache *cache);

void mcppp2d_compute_tangent(const mcppp2d_params *par,
                             const mcppp2d_cache *cache,
                             double d[MCPPP2D_NCOMP_STRAIN][MCPPP2D_NCOMP_STRAIN]);

void mcppp2d_pack_eqother(const mcppp2d_cache *cache,
                          double eqother[MCPPP2D_NCOMP_EQOTHER]);

void mcppp2d_pack_other(const mcppp2d_cache *cache,
                        double other[MCPPP2D_NCOMP_OTHER]);

void mcppp2d_unpack_other(const double other[MCPPP2D_NCOMP_OTHER],
                          mcppp2d_cache *cache);

#ifdef __cplusplus
}
#endif

#endif
