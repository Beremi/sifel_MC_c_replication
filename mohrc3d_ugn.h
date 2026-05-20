#ifndef MOHRC3D_UGN_H
#define MOHRC3D_UGN_H

#include <stdio.h>
#include <string.h>

struct XFILE;
struct vector;
struct matrix;

class mohrc3d_ugn
{
 public:
  enum
  {
    // SIFEL 3D Voigt order:
    // [xx, yy, zz, yz, xz, xy] for strains and stresses.
    NCOMP_STRAIN  = 6,
    NCOMP_STRESS  = 6,

    // Current-iteration point buffer:
    // epsp(6), epsp_prev(6), return_type(1)
    NCOMP_OTHER   = 13,

    // Converged point buffer:
    // epsp(6)
    NCOMP_EQOTHER = 6
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
    O_EP_XX = 0,
    O_EP_YY = 1,
    O_EP_ZZ = 2,
    O_GP_YZ = 3,
    O_GP_XZ = 4,
    O_GP_XY = 5,

    O_EP_PREV_XX = 6,
    O_EP_PREV_YY = 7,
    O_EP_PREV_ZZ = 8,
    O_GP_PREV_YZ = 9,
    O_GP_PREV_XZ = 10,
    O_GP_PREV_XY = 11,

    O_RET = 12
  };

  long read(XFILE *in);
  void print(FILE *out);
  void nlstresses(long ipp, long ido);
  void nonloc_nlstresses(long ipp, long ido);
  void nlstresses(const vector &strain, const vector &eqstatev, vector &stress, vector &statev);
  void stiffmat(const vector &strain, const vector &statev, const vector &stress, matrix &d);
  void matstiff(matrix &d, long ipp, long ido);
  void updateval(const vector &statev, vector &eqstatev);
  void updateval(long ipp, long im, long ido);
  void giveirrstrains(long ipp, long ido, vector &epsp);

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
  void compute_response_internal(const double strain_internal[6],
                                 int forced_return_type,
                                 double stress_internal[6],
                                 double tangent_internal[6][6],
                                 int &return_type,
                                 double eig[3],
                                 double sigma[3]) const;
};

#endif
