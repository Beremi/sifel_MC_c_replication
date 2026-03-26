#ifndef MATMODEL_H
#define MATMODEL_H

#include <stdio.h>
#include <string.h>

struct XFILE;
struct vector;
struct matrix;

class matmodel
{
 public:
  enum
  {
    // Plane-strain Voigt layout used both by the Matlab prototype and by this rewrite:
    // [xx, yy, xy, zz] for strains and [xx, yy, xy, zz] for stresses.
    NCOMP_STRAIN  = 4,
    NCOMP_STRESS  = 4,

    // Current-iteration point buffer:
    // epsp(4), epsp_prev(4), return_type(1)
    NCOMP_OTHER   = 9,

    // Converged point buffer:
    // epsp(4)
    NCOMP_EQOTHER = 4
  };

  enum return_type_code
  {
    RET_ELASTIC    = 0,
    RET_SMOOTH     = 1,
    RET_LEFT_EDGE  = 2,
    RET_RIGHT_EDGE = 3,
    RET_APEX       = 4
  };

  enum other_offset
  {
    // Numeric offsets in the packed point arrays statev / eqstatev.
    // Current plastic strain epsp in Voigt order.
    O_EP_XX = 0,
    O_EP_YY = 1,
    O_GP_XY = 2,
    O_EP_ZZ = 3,

    // Previous plastic strain copied from eqstatev for tangent recomputation.
    O_EP_PREV_XX = 4,
    O_EP_PREV_YY = 5,
    O_GP_PREV_XY = 6,
    O_EP_PREV_ZZ = 7,

    // Return type stored as double so it can travel through the generic state array.
    O_RET = 8
  };

  long read(FILE *in);
  void print(FILE *out);

  void nlstresses(const vector &strain, const vector &eqstatev, vector &stress, vector &statev);
  void stiffmat(const vector &strain, const vector &statev, const vector &stress, matrix &d);
  void updateval(const vector &statev, vector &eqstatev);

  // Material parameters and precomputed constants for the return mapping.
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

 protected:
  void compute_ordered_trial_spectral(const vector &E_trial,
                                      vector &eig,
                                      vector &Eig_1,
                                      vector &Eig_2,
                                      vector &Eig_3,
                                      matrix &EIG_1,
                                      matrix &EIG_2,
                                      matrix &EIG_3) const;
  void compute_return_denominators(double &denom_s,
                                   double &denom_l,
                                   double &denom_r) const;
  void compute_returned_principal_stresses(const vector &eig,
                                           double trace_E,
                                           int return_type,
                                           double lambda_s,
                                           double lambda_l,
                                           double lambda_r,
                                           vector &sigma) const;
};

#endif
