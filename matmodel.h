#ifndef MATMODEL_H
#define MATMODEL_H

#include <stdio.h>
#include <string.h>
#include "vector.h"

struct XFILE;
struct matrix;

enum
{
  MATMODEL_NCOMP_STRAIN  = 4,
  MATMODEL_NCOMP_STRESS  = 4,
  MATMODEL_NCOMP_EQOTHER = 4,
  MATMODEL_NCOMP_OTHER   = 68
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

  MATMODEL_IO_PROJ_1 = 5,
  MATMODEL_IO_PROJ_2 = 9,
  MATMODEL_IO_PROJ_3 = 13,

  MATMODEL_IO_HESS_1 = 17,
  MATMODEL_IO_HESS_2 = 33,
  MATMODEL_IO_HESS_3 = 49,

  MATMODEL_IO_SIGMA_1 = 65,
  MATMODEL_IO_SIGMA_2 = 66,
  MATMODEL_IO_SIGMA_3 = 67
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
  void updateval(const vector &statev, vector &eqstatev);

  matmodel_params par;

 protected:
  void fill_response(const vector &strain, const vector &eqstatev, vector &stress, vector &statev);

  vector cached_strain;
  vector cached_eqother;
  vector cached_other;
  long cached_response_valid;

 private:
  matmodel(const matmodel &) = delete;
  matmodel &operator=(const matmodel &) = delete;
};

#endif
