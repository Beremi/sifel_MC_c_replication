#ifndef DEF_MATRIX
#define DEF_MATRIX
#include <assert.h>
#include <string.h>
#include "vector.h"

//#define DEBUG_MATRIX

struct matstat
{
  unsigned dealloc:1;
  unsigned symm:1;
  unsigned transp:1;
  matstat() {dealloc=0; symm=0; transp=0;};
};



struct matrix
/** This file declares struct matrix, which implements
    %matrix with %elements type of double.
    There are also declarations of the functions for the %matrix and %vector computing */
{
  long m;    ///< number of rows
  long n;    ///< number of columns
  double *a; ///< pointer to onedimensional array with matrix elements stored in the rows
  long size; ///< real length of array a (due to reallocation of memory)
  matstat stat; ///< bit array of vector status flags (deallocation, symmetric, transpose, etc)   

  inline matrix() {m = n = size = 0L; a = NULL;}; ///< default constructor
  matrix(long m, long n);               ///< allocating constructor
  matrix(const matrix &mat);            ///< copy constructor
  matrix& operator=(const matrix &mat); ///< assignment operator

  /**
    The constructor assigns preallocated memory to the member array a.
   
    @param m is the number of rows
    @param n is the number of columns
    @param dealloc is the flag for deallocation control (0=no deallocation ,1=release memory by delete operator)
    @param ptr is the pointer to the allocated memory which will be used for storage of m*n elements
   
    created  4.5.2015, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  */
  inline matrix(long m, long n, unsigned dealloc, double *ptr)
  {
    matrix::m = m;
    matrix::n = n;
    size = m*n;
    a = ptr;
    stat.dealloc = dealloc;

    memset (a, 0, size*sizeof(*a));

   #ifdef DEBUG_MATRIX
    Acm++;
    Ama += size;
    if (Ama > Ammax)
      Ammax = Ama;
   #endif
  };
  
  /**
   The operator enables access to the rows of elements of the member array a.

    @param i is the number of row

   created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
   modified  9.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  */
  inline double* operator [] (long i)
  {
    #ifdef DEBUG_MATRIX
     if (i >= m)
       print_err("matrix row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, m-1);

     assert(i < m);
    #endif
    return (&a[i * n]);
  };


  /**
   The operator enables CONSTANT access to the rows of elements of the member array a.

    @param i is the number of row

   created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
   modified  9.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  */
  inline const double* operator [] (long i) const
  {
    #ifdef DEBUG_MATRIX
     if (i >= m)
       print_err("matrix row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, m-1);

     assert(i < m);
    #endif
    return (&a[i * n]);
  };


  /**
    The operator enables access to the elements of the member array a.

    @param i is the index of row
    @param j is the index of column

    created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
    modified  9.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  */
  inline double& operator () (long i, long j)
  {
    #ifdef DEBUG_MATRIX
     if (i >= m)
       print_err("matrix row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, m-1);
     if (j >= n)
       print_err("matrix column index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, j, n-1);

     assert((i < m) && (j < n));
    #endif
    return (a[i * n + j]);
  };

  /**
    The operator enables CONSTANT access to the elements of the member array a.

    @param i is the index of row
    @param j is the index of column

    created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
    modified  9.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  */
  inline const double& operator () (long i, long j) const
  {
    #ifdef DEBUG_MATRIX
     if (i >= m)
       print_err("matrix row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, m-1);
     if (j >= n)
       print_err("matrix column index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, j, n-1);

     assert((i < m) && (j < n));
    #endif
    return (a[i * n + j]);
  };

  /// creates %matrix which contains reference to %matrix src
  void makerefm(matrix &src);

  /// creates %matrix which contains reference to array ptr
  void makerefm(long mp, long np, double *ptr);
  
  inline ~matrix() ///< destructor
  {
   #ifdef DEBUG_MATRIX
    if (a != NULL)
    {
      Acm--;
      Ama -= size;
    }
   #endif

    if (stat.dealloc)
      delete [] a;
  };
};


/** 
  Struct imatrix implements %matrix with %elements type of long. 
*/
struct imatrix
{
  long m;    ///< number of rows
  long n;    ///< number of columns
  long *a;   ///< pointer to onedimensional array with matrix elements stored in the rows
  long size; ///< real length of array a (due to reallocation of memory)
  matstat stat; ///< bit array of vector status flags (deallocation, symmetric, transpose, etc)   
  
  inline imatrix() {m = n = 0L; a = NULL;}; ///< default constructor
  imatrix(long m, long n);          ///< allocating constructor
  imatrix(const imatrix &mat);      ///< copy constructor

  /**
    The constructor assigns preallocated memory to the member array a.
   
    @param m is the number of rows
    @param n is the number of columns
    @param dealloc is the flag for deallocation control (0=no deallocation ,1=release memory by delete operator)
    @param ptr is the pointer to the allocated memory which will be used for storage of m*n elements
    
    created  4.5.2015, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  */
  inline imatrix(long m, long n, unsigned dealloc, long *ptr)
  {
    imatrix::m = m;
    imatrix::n = n;
    size = m*n;
    a = ptr;
    stat.dealloc = dealloc;

    memset (a, 0, size*sizeof(*a));

   #ifdef DEBUG_MATRIX
    Acm++;
    Ama+=size;
    if (Ama > Ammax)
      Ammax = Ama;
   #endif
  };
  

  /**
   The operator enables access to the rows of elements of the member array a.

    @param i is the number of row

   created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
   modified  9.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  */
  inline long* operator [] (long i)
  {
    #ifdef DEBUG_MATRIX
     if (i >= m)
       print_err("imatrix row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, m-1);

     assert(i < m);
    #endif
    return (&a[i * n]);
  };

  /**
   The operator enables access to the elements of the member array a.

    @param i is the index of row
    @param j is the index of column

   created  29.5.2000, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
   modified  9.5.2001, Tomas Koudelka, koudelka@cml.fsv.cvut.cz
  */
  inline long& operator () (long i, long j)
  {
    #ifdef DEBUG_MATRIX
     if (i >= m)
       print_err("imatrix row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, m-1);
     if (j >= n)
       print_err("imatrix column index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, j, n-1);

     assert((i < m) && (j < n));
    #endif
    return (a[i * n + j]);
  };

  inline ~imatrix() ///< destructor
  {
   #ifdef DEBUG_MATRIX
    if (a != NULL)
    {
      Acm--;
      Aima -= size;
    }
   #endif

    if (stat.dealloc)
      delete [] a;
  };
};


#ifdef DEBUG_MATRIX
 long give_acm();
 long give_ama();
 long give_ammax();
#endif

/// allocates %matrix to the dimensions m x n
long allocm(long m, long n, matrix &mat);
/// allocates imatrix to the dimensions m x n
long allocm(long m, long n, imatrix &mat);
/// allocates double array to the dimensions m x n
long allocm(long m, long n, double **&mat);
/// allocates and initializes identity %matrix of dimension n x n
long allocim(long n, matrix &mat);
/// allocates and initializes identity %matrix of dimension n x n with help of preallocated memory
long allocim(long n, matrix &mat, unsigned dealloc, double *ptr);

/// reallocates matrix to the dimensions m x n
long reallocm(long m, long n, matrix &mat);
/// reallocates matrix to the dimensions m x n with help of preallocated memory
long reallocm(long m, long n, matrix &mat, unsigned dealloc, double *ptr);
/// reallocates imatrix to the dimensions m x n
long reallocm(long m, long n, imatrix &mat);
/// reallocates imatrix to the dimensions m x n with help of preallocated memory
long reallocm(long m, long n, imatrix &mat, unsigned dealloc, long *ptr);

/// creates %matrix which contains reference to %matrix src
long makerefm(matrix &ref, matrix &src);
/// creates %matrix which contains reference to array ptr
long makerefm(matrix &ref, double *ptr, long m, long n);

/// copies contents of matrix from the src to dest
long copym(const matrix &src, matrix &dest);

/// copies contents of matrix from the src to dest
long copym(const matrix &src, double *dest);

/// copies contents of matrix from the src to dest
long copym(const double *src, matrix &dest);

/// zeros contents of %matrix
long nullm(matrix &mat);

/// fills contents of %matrix with given value c
long fillm(double c, matrix &mat);
/// fills i-th row of %matrix with given value c
long fillrow(double c, long i, matrix &mat);
/// fills i-th column of %matrix with given value c
long fillcol(double c, long i, matrix &mat);
/// sets contents of mat to identity %matrix
long identm(matrix &mat);

/// deallocates %matrix
long destrm(matrix &mat);
/// deallocates %matrix given by 2 dimensional array
long destrm(double **&mat, long m);
/// deallocates integer %matrix given by 2 dimensional array
long destrm(long **&mat, long m);

/// adds 2 matrices C := A+B
long addm(const matrix &a, const matrix &b, matrix &c);

/// substracts 2 matrices C := A-B
long subm(const matrix &a, const matrix &b, matrix &c);

/// adds 2 matrices multplied by constant to the given matrix C: C := C + ac*A + bc*B
long addmultm(const matrix &a, double ac, const matrix &b, double bc, matrix &c);

/// performs tensor product of two vectors c_{ij} = a_{i}. b{j}
void tensprd (vector &a, vector &b, matrix &c);

/// multiplies two matrices C := A.B
long mxm(const matrix &a, const matrix &b, matrix &c);
/// multiplies two matrices C := A.B given by arrays
void mxm(const double *a, const double *b, double *c, long l, long m, long n);

/// multiplies two matrices C := A.B^T
long mxmt(const matrix &a, const matrix &b, matrix &c);
/// multiplies two matrices C := A.B^T given by arrays
void mxmt(const double *a, const double *b, double *c, long l, long m, long n);

/// multiplies two matrices C := A^T.B
long mtxm(const matrix &a, const matrix &b, matrix &c);
/// multiplies two matrices C := A^T.B given by arrays
void mtxm(const double *a, const double *b, double *c, long l, long m, long n);

/// multiplies transposed %matrix a from right by %matrix b, C := A^T . B, matrices are stored by columns
void mtxmccr(const double *a, const double *b, double *c, long l, long m, long n);

/// multiplies %matrix by constant B := c.A
long cmulm(const double c, const matrix &a, matrix &b);
/// multiplies %matrix by constant A := c.A
long cmulm(double c, matrix &a);

/// transposes %matrix A_t := A^T
long tranm(const matrix &a, matrix &at);
/// transposes %matrix A := A^T
long tranm(matrix &a);

/// multiplies %matrix a from left by row %vector u  v := u.A
long vxm(const vector &u, const matrix &a, vector &v);
/// multiplies transposed %matrix a from left by row %vector u  v := u.A^T
long vxmt(const vector &u, const matrix &a, vector &v);

/// multiplies %matrix a from right by %vector u  v := A.u
long mxv(const matrix &a, const vector &u, vector &v);
/// multiplies %matrix a given by array from right by %vector u given by array  v := A.u
void mxv(const double *a, const double *b, double *c, long m, long n);

/// multiplies %matrix a stored by columns from the right by %vector b   c := A.b 
void mxvc(const double *a, const double *b, double *c, long m, long n);
/// computes i-th component of the %vector v given by product of %matrix a and %vector u  v := A.u
long mixv(const matrix &a, const vector &u, double &vi, long i);

/// multiplies transposed %matrix a from right by %vector u  v := A^T.u
long mtxv(const matrix &a, const vector &u, vector &v);
/// multiplies transposed %matrix a given by array from right by %vector u given by array  v := A^T.u
void mtxv(const double *a, const double *b, double *c, long m, long n);
/// multiplies transposed %matrix a stored by columns from the right by %vector b   c := A^T.b 
void mtvc(const double *a, const double *b, double *c, long m, long n);

/// multiplies column vector u with row vector, i.e. it performs tensor product, A := u * v
long vxv(const vector &u, const vector &v, matrix &a);

/// multiplies row vector %v with symetric %matrix m and column %vector v, i.e. answer := (v^T).A.v
long vxmxv (const vector &v, const matrix &m, double &answer);
///  multiplies row vector %v given by array with symetric %matrix m given by array and column %vector v, i.e. answer := (v^T).A.v
long vxmxv (const double *v, const double *m, long dim, double &answer);


/// performs Gauss elimination on the %matrix a, result is stored in the %matrix b
long gause(const matrix &a, matrix &b,double zero);
/// solves system of linear algebraic equations by the Gauss elimination (A.x = r), resulting %matrix is b, resulting %vector is stored in r
long gause(const matrix &a, matrix &b, vector &r,double zero);
/// solves system of linear algebraic equations by the Gauss elimination (A.x = r), resulting %matrix is b, resulting %vector is stored in r
long gause(const matrix &a, matrix &b, vector &r, vector &sol,double zero);

long gemp (matrix &a,vector &x,vector &y,double zero,long pivot);

/// inverts %matrix a, i.e. computes A^{-1}
long invm(const matrix &a, matrix &b,double zero);

/// computes determinant of the %matrix a, i.e. det(A)  
long detm(const matrix &a, double &det);
void det2x2 (const matrix &a,double &det);
void det3x3 (const matrix &a,double &det);

/// prints the %matrix to the file with the given width and precision
long printm(const matrix &a, FILE *out = stdout, int prec = 3, int width = 11);
/// prints the %matrix to the file
long printm(FILE *out, const matrix &a);

/// swaps %matrix columns i and j
long mswapc(matrix &mat, long i, long j);
/// swaps %matrix rows i and j
long mswapr(matrix &mat, long i, long j);

/// performs condesation of selected dofs of the given stiffness %matrix
long condense_matrix(matrix &sm, ivector &cu);
/// performs condesation of selected dofs of the given load %vector nf according to corresponding stiffness %matrix sm
long condense_vector(matrix &sm, vector &nf, ivector &cu);

/// computes %matrix product K += B^T . D . B. jac (only upper triangle components are evaluated)
void bdbj (double *k,const double *b,const double *d,double jac,long m,long n);
/// computes %matrix product D += A^T . B . C . jac
void bdbjac (matrix &d,const matrix &a,const matrix &b,const matrix &c,double jac);

/// computes %matrix product M += N^T . N . jac
void nnj (double *m,const double *n,double jac,long mm,long nn);
/// computes %matrix product D += A^T . B . jac
void nnjac (matrix &d,const matrix &a,const matrix &b,double jac);


/// transforms %vector l given in the local coordinate system to %vector g in global coordinate system (g = T . l)
void lgvectortransf (vector &g,const vector &l,const matrix &tmat);
/// transforms %vector g given in the global coordinate system to %vector l in local coordinate system (l = T^T . g)
void glvectortransf (const vector &g,vector &l,const matrix &tmat);
/// transforms tensor(%matrix) a given in the global coordinate system to the local one (A_l := T^T . A_g . T)
void glmatrixtransf (matrix &a,const matrix &tmat);
/// transforms tensor(%matrix) a given in the local coordinate system to the global one (A_g := T . A_l . T^T)
void lgmatrixtransf (matrix &a,const matrix &tmat);

/// transforms block %matrix a given in the local coordinate system to the global one (A_g := T . A_l . T^T) (number of blocks is given by dimension of T)
void lgmatrixtransfblock (matrix &a, const matrix &t);
/// transforms block %matrix a given in the global coordinate system to the local one (A_l := T^T . A_g . T) (number of blocks is given by dimension of T)
void glmatrixtransfblock (matrix &a, const matrix &t);
/// transforms block %vector v given in the global coordinate system to the local one (v_l := T^T . v_g) (number of blocks is given by dimension of T)
void glvectortransfblock (vector &v, const matrix &t);
/// transforms block %vector v given in the local coordinate system to the global one (v_g := T^T . v_l) (number of blocks is given by dimension of T)
void lgvectortransfblock (vector &v, const matrix &t);

/// localizes components of local(element) %vector to the global(problem) %vector
void locglob (double *gv, const double *lv,const long *cn,long n);
/// localizes components of local(element) %vector to the global(problem) %vector
void globloc (const double *gv, double *lv, const long *cn,long n);

///
void locvecmat(double *mat, const double *vect, const long *cn, long ci, long m, long n);
///  localizes components of local %matrix lm to the global one gm
void mat_localize (matrix &gm, matrix &lm, long *rcn, long *ccn);

/// computes determinant of %matrix which equals to double area of triangle given by coordinate vectors x and y
double det2d (const double *x,const double *y);
/// computes determinant of %matrix which equals to six-multiple of tetrahedron volume which is given coordinate vectors x,y,z
double det3d (const double *x,const double *y,const double *z);

/// computes eigenvalues(principal values) of the given %matrix v by Jacobi method
void princ_val (const matrix &v, vector &pval, matrix &pvect, long ni, double err, double zero, long n, long normalize);

/// solves system of algebraic equations by Gauss method with pivot selection
void gemp (double *a,double *x,double *y,long n,long m,double limit,long pivot);

/// solves system of algebraic equations by LU decomposition
void lu_full (double *a,double *x,double *y,long n,double zero,long tc);

///  assembles %matrix of least square problem for extrapolation of nodal values from integration point values
void matassem_lsm (double *lsm,vector &natcoord);

///  assembles right hand side %vector of least square problem for extrapolation of nodal values from integration point values
void rhsassem_lsm (double *rhs,vector &natcoord,vector &values);

/// solves system of equation given by least square problem for extrapolation of nodal values from integration point values
void solve_lsm (double *lsm,double *lhs,double *rhs,double zero,long n,long m);

/// extracts squared %matrix b(ncomp,ncomp) from %matrix a
void extractm (matrix &a,matrix &b,long fi, long ncomp);
/// extracts submatrix b(nr,nc) from %matrix a
void extractblock (double *a,double *b,long n,long nr,long nc,long fri,long fci);
/// extracts submatrix b from %matrix a
void extractblock (matrix &a, matrix &b, long fri,long fci);
/// stores submatrix b in the %matrix a
void storeblock (matrix &a, matrix &b, long fri,long fci);
/// extracts the i-th row from the %matrix m and stores it to %vector dest
long extractrow(const matrix &m, long i, vector &dest);
/// extracts the j-th column from %matrix m and stores it to %vector dest
long extractcol(const matrix &m, long j, vector &dest);

///  function diagonalizes a %matrix mat
void diagonalization (matrix &mat);

/// function computes normalized matrix
double normalize(matrix &mat);

/// function computes Euclidean %matrix norm
double norm(matrix &a);

///  power method for determination of the largest eigenvalue of a matrix
double power_method (matrix &a,double err,long ni,double zero);

void jacobi_rot (double *a,double *evec,double *eval,long n,long ni,long &ani,double limit, long normalize);
#endif
