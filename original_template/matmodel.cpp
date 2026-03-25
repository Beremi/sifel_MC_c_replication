#include "matmodel.h"
#include "vector.h"
#include "matrix.h"

#include <stdio.h>

/*
 Example of vector type usage:

  vector a(n); // defines variable a type of vector allocated to the dimension n

  reallocv(n, a) - (re)allocates vector a to the dimensions n
  a.n - represents the number of vector components
  a[i] - returns i-th vector component
  a(i) - returns i-th vector component  
  scprd(a,b) - returns vector scalar product of a.b
  normv(a) - returns the vector norm, ||a||
  cmulv (c, a) - scales vector a by factor c
  fillv(c, a) - fills all vector components by value c
  nullv(a) - zeroes all components of vector a

 Example of matrix type usage:

  matrix a(nr,nc);  // defines variable a type of matrix allocated to the dimension nr x nc

  a.m represents the number of matrix rows in a
  a.n represents the number of matrix columns in a

  a[i][j] - returns matrix element at the i-th row and j-th column
  a(i,j)  - returns matrix element at the i-th row and j-th column

  reallocm(nr, nc, a) - (re)allocates matrix a to the dimensions nr x nc,
                        where nr is the number of rows and nc is the number of columns
  mxm(a,b,c) - computes matrix multiplication c = a*b
  mtxm(a,b,c) - computes matrix multiplication c = a^T*b (matrix a is transposed)
  mxmt(a,b,c) - computes matrix multiplication c = a*b^T (matrix b is transposed)
  mxv(a,v,u) - computes multiplication of matrix a by vector v  u = a*v
  normm(a) - returns the Euklidean matrix norm ||a||
  cmulm (c, a) - scales matrix a by factor c
  fillm(c, a) - fills all matrix components by value c
  nullm(a) - zeroes all components of matrix a
*/



/**
  The function reads material model parameters and stress return setup from the opened text.

  @param[in] in - pointer to the opened text file.

  @retval 0 - on success
  @retval 1 - in the case of an error
*/
long matmodel::read(FILE *in)
{
  // reading of data members representing the material model parameters from the text file
  // will be written by Tomas Koudelka

  return 0;
}



/**
  The function prints material model parameters and stress return setupto the opened text file.

  @param[in] out - pointer to the opened text output file.

*/  
void matmodel::print(FILE *out)
{
  // printing of data members representing the material model parameters and stress return algorithm setup
  // to the text file
  // will be written by Tomas Koudelka
}



/**
  The function computes actual stresses with respect to the attained strains and state variables.

  @param[in] strain - array of actual strain, components are ordered as follows
                      eps_x, eps_y, gamma_xy, eps_z - for the plain strain/plain stress problem
                      eps_x, eps_y, eps_z, gamma_yz, gamma_xz, gamma_xy - for the space stress problem

  @param[in] eqstatev - array of state variables from the last equilibrium state
                        state variables order is arbitrary, the author of the model defines it

  @param[out] stress - array of the resulting stress components, it must be computed,
                       components are ordered as follows
                       sig_x, sig_y, tau_xy, sig_z - for the plain strain/plain stress problem
                       sig_x, sig_y, sig_z, tau_yz, tau_xz, tau_xy - for the space stress problem  

  @param[in] statev - array of the resulting state variables calculated for the actual strains
                      state variables order is arbitrary, the author of the model defines it,
                      but it must be the same as in eqstatev.

  @return The function does not return anything, the results are stored in the arrays stress and
          statev passed in as arguments.
*/
void matmodel::nlstresses(const vector &strain, const vector &eqstatev, vector &stress, vector &statev)
{
  // will be written by UGN (export from Matlab)
}



/**
  The function computes material stiffness matrix with respect to the attained strains and
  state variables.

  @param[in] strain - array of actual strain, components are ordered as follows
                      eps_x, eps_y, gamma_xy, eps_z - for the plain strain/plain stress problem
                      eps_x, eps_y, eps_z, gamma_yz, gamma_xz, gamma_xy - for the space stress problem

  @param[in] eqstatev - array of state variables from the last equilibrium state
                        state variables order is arbitrary, the author of the model defines it

  @param[in] stress - array of the resulting stress components, it must be computed,
                      components are ordered as follows
                      sig_x, sig_y, tau_xy, sig_z - for the plain strain/plain stress problem
                      sig_x, sig_y, sig_z, tau_yz, tau_xz, tau_xy - for the space stress problem

  @param[out] d -  the resulting material stiffness matrix, it must be calculated.

  
  @return The function does not return anything, the resulting matrix is stored in the
          matrix d passed in as arguments.
*/
void matmodel::stiffmat(const vector &strain, const vector &eqstatev, const vector &stress, matrix &d)
{
  // will be written by UGN (export from Matlab)
}
