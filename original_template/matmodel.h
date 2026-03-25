#ifndef MATMODEL_H
#define MATMODEL_H

#include <stdio.h>



struct XFILE;
struct vector;
struct matrix;

class matmodel
{
 public:
  /// function reads material model parameters and stress return algorithm setup from the given opend text file
  long read(FILE *in);

  /// function prints material model parameters and stress return setup
  void print(FILE *out);

  /// function computes actual stresses with respect to the attained strains and state variables
  void nlstresses(const vector &strain, const vector &eqstatev, vector &stress, vector &statev);

  /// function computes material stiffness matrix with respect to the attained strains and state variables
  void stiffmat(const vector &strain, const vector &eqstatev, const vector &stress, matrix &d);

  // data members represent material parameters and stress return algorithm setup
};


#endif
