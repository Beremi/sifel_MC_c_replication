#ifndef MATMODEL_H
#define MATMODEL_H

#include <stdio.h>
#include <string.h>

struct XFILE;
struct vector;
struct matrix;

enum
{
  // Plane-strain Voigt layout used both by the Matlab prototype and by this rewrite:
  // [xx, yy, xy, zz] for strains and [xx, yy, xy, zz] for stresses.
  MATMODEL_NCOMP_STRAIN  = 4,
  MATMODEL_NCOMP_STRESS  = 4,

  // The point buffer stores the current plastic strain plus the auxiliary quantities
  // that Matlab passes from constitutive_problem.m to stiffness_matrix.m.
  MATMODEL_NCOMP_OTHER   = 50,

  // eqstatev and statev share the same per-point packed layout.
  MATMODEL_NCOMP_EQOTHER = MATMODEL_NCOMP_OTHER
};

enum matmodel_return_type
{
  MATMODEL_RETURN_ELASTIC    = 0,
  MATMODEL_RETURN_SMOOTH     = 1,
  MATMODEL_RETURN_LEFT_EDGE  = 2,
  MATMODEL_RETURN_RIGHT_EDGE = 3,
  MATMODEL_RETURN_APEX       = 4
};

enum matmodel_other_index
{
  // Current plastic strain epsp in Voigt order.
  MATMODEL_IO_EP_XX = 0,
  MATMODEL_IO_EP_YY = 1,
  MATMODEL_IO_GP_XY = 2,
  MATMODEL_IO_EP_ZZ = 3,

  // Return type stored as double so it can travel through the generic state array.
  MATMODEL_IO_RETURN_TYPE = 4,

  // Ordered trial eigenvalues eig_1, eig_2, eig_3.
  MATMODEL_IO_EIG_1 = 5,
  MATMODEL_IO_EIG_2 = 6,
  MATMODEL_IO_EIG_3 = 7,

  // First derivatives of the ordered trial eigenvalues. These correspond to the Matlab
  // arrays Eig_1, Eig_2, Eig_3 returned by constitutive_problem.m.
  MATMODEL_IO_PROJ_1 = 8,
  MATMODEL_IO_PROJ_2 = 12,
  MATMODEL_IO_PROJ_3 = 16,

  // Second derivatives of the ordered trial eigenvalues. Matlab stores these as the
  // reduced 9-component arrays EIG_1..EIG_3 in column-major order:
  // [11,21,31,12,22,32,13,23,33].
  MATMODEL_IO_HESS_1 = 20,
  MATMODEL_IO_HESS_2 = 29,
  MATMODEL_IO_HESS_3 = 38,

  // Ordered principal stresses sigma_1, sigma_2, sigma_3 after return mapping.
  MATMODEL_IO_SIGMA_1 = 47,
  MATMODEL_IO_SIGMA_2 = 48,
  MATMODEL_IO_SIGMA_3 = 49
};

struct matmodel_params
{
  double young;
  double poisson;
  double c;
  double phi;

  double shear;
  double bulk;
  double lame;
  double sin_phi;
  double cos_phi;
  double c_bar;
};

class matmodel
{
 public:
  long read(FILE *in);
  void print(FILE *out);

  void nlstresses(const vector &strain, const vector &eqstatev, vector &stress, vector &statev);
  void stiffmat(const vector &strain, const vector &eqstatev, const vector &stress, matrix &d);

  // data members represent material parameters and stress return algorithm setup
  matmodel_params par;
};

#endif
