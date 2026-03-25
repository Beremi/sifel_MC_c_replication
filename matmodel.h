#ifndef MATMODEL_H
#define MATMODEL_H

#include <stdio.h>
#include "mcppp2d_core.h"

struct XFILE;
struct vector;
struct matrix;

class matmodel
{
 public:
  matmodel();

  /// function reads material model parameters from the given opened text file
  long read(FILE *in);

  /// function prints material model parameters and buffer layout summary
  void print(FILE *out);

  /// function computes actual stresses with respect to the attained strains and state variables
  void nlstresses(const vector &strain, const vector &eqstatev, vector &stress, vector &statev);

  /// function computes material stiffness matrix; in this template variant it recomputes the cache
  void stiffmat(const vector &strain, const vector &eqstatev, const vector &stress, matrix &d);

  /// function computes material stiffness matrix directly from the cached "other" state variables
  void stiffmat_from_statev(const vector &statev, matrix &d);

  /// function updates converged state variables (copy current plastic strains to eqstatev)
  void updateval(const vector &statev, vector &eqstatev);

  /// number of cached variables in the current state buffer (SIFEL other)
  long give_num_of_statev() const;

  /// number of converged variables in the equilibrium state buffer (SIFEL eqother)
  long give_num_of_eqstatev() const;

  mcppp2d_params par;

  /// only plane strain / 2D-in-3D stress state is supported by this model
  long plane_strain_only;

 protected:
  void write_matrix_to_statev(const double a[4][4], long start, vector &statev) const;
  void write_vector_to_statev(const double a[4], long start, vector &statev) const;
  void write_tangent_to_matrix(const double a[4][4], matrix &d) const;
  void fill_response(const vector &strain, const vector &eqstatev, vector &stress, vector &statev) const;
};

inline long matmodel::give_num_of_statev() const
{
  return MCPPP2D_NCOMP_OTHER;
}

inline long matmodel::give_num_of_eqstatev() const
{
  return MCPPP2D_NCOMP_EQOTHER;
}

#endif
