#ifndef MATMODEL_H
#define MATMODEL_H

#include <stdio.h>

struct XFILE;
struct vector;
struct matrix;

enum
{
  MATMODEL_NCOMP_STRAIN  = 4,
  MATMODEL_NCOMP_STRESS  = 4,
  MATMODEL_NCOMP_EQOTHER = 4,
  MATMODEL_NCOMP_OTHER   = 71
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
  MATMODEL_IO_EP_XX = 0,
  MATMODEL_IO_EP_YY = 1,
  MATMODEL_IO_GP_XY = 2,
  MATMODEL_IO_EP_ZZ = 3,

  MATMODEL_IO_RETURN_TYPE = 4,

  MATMODEL_IO_EIG_1 = 5,
  MATMODEL_IO_EIG_2 = 6,
  MATMODEL_IO_EIG_3 = 7,

  MATMODEL_IO_PROJ_1 = 8,
  MATMODEL_IO_PROJ_2 = 12,
  MATMODEL_IO_PROJ_3 = 16,

  MATMODEL_IO_HESS_1 = 20,
  MATMODEL_IO_HESS_2 = 36,
  MATMODEL_IO_HESS_3 = 52,

  MATMODEL_IO_SIGMA_1 = 68,
  MATMODEL_IO_SIGMA_2 = 69,
  MATMODEL_IO_SIGMA_3 = 70
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
  matmodel();

  long read(FILE *in);
  void print(FILE *out);

  void nlstresses(const vector &strain, const vector &eqstatev, vector &stress, vector &statev);
  void stiffmat(const vector &strain, const vector &eqstatev, const vector &stress, matrix &d);
  void stiffmat_from_statev(const vector &statev, matrix &d);
  void updateval(const vector &statev, vector &eqstatev);

  long give_num_of_statev() const;
  long give_num_of_eqstatev() const;

  matmodel_params par;
  long plane_strain_only;

 protected:
  void write_tangent_to_matrix(const double a[4][4], matrix &d) const;
  void fill_response(const vector &strain, const vector &eqstatev, vector &stress, vector &statev) const;
};

inline long matmodel::give_num_of_statev() const
{
  return MATMODEL_NCOMP_OTHER;
}

inline long matmodel::give_num_of_eqstatev() const
{
  return MATMODEL_NCOMP_EQOTHER;
}

#endif
