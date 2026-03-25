#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifdef DEBUG_MATRIX
 static unsigned long Acm; ///< alocation counter for matrix class
 static unsigned long Ammax; ///< peak of alocation memory of all matrix instances
 static unsigned long Ama; ///< actual allocated memory of all matrix instances

 static unsigned long Acim; ///< alocation counter for imatrix class
 static unsigned long Aimmax; ///< peak of alocation memory of all imatrix instances
 static unsigned long Aima; ///< actual allocated memory of all imatrix instances
#endif 



/**
   The constructor allocates memory form the heap to the member array a.
   
   @param[in] m is the number of rows
   @param[in] n is the number of columns
   
   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
matrix :: matrix(long m, long n)
{
  matrix::m = m;
  matrix::n = n;
  size = m*n;
  a = new double[size];
  if (a)
  {
    stat.dealloc = 1;

    memset (a, 0, size*sizeof(*a));

    #ifdef DEBUG_MATRIX
    Acm++;
    Ama += size;
    if (Ama > Ammax)
      Ammax = Ama;
   #endif
  }
  else
    print_err("cannot allocate memory for matrix a(%ld,%ld)", __FILE__, __LINE__, __func__, m, n);
}



/**
   The copy constructor creates copy of the object given by the constant reference parameter.
   
   @param[in] mat is reference to the %matrix object which should be copied to the given instance
   
   @author Created  16.10.2003 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
matrix :: matrix (const matrix &mat)
{
  long asize = mat.m*mat.n; // mat.size may be greater than actual matrix size, so actual size must be calculated

  m = mat.m;
  n = mat.n;
  if (asize > size)
  {  
    if (a && stat.dealloc)
      delete [] a;

    size = asize;
    a = new double[size];
    if (a)
    {
      stat.dealloc = 1;

      print_warning("Copy constructor of matrix is called.\n"
                    "Please check your code and make sure you want use a copy constructor.\n"
                    "Consider to use copym function rather", 
                    __FILE__, __LINE__, __func__);

      #ifdef DEBUG_MATRIX
      Acm++;
      Ama += size;
      if (Ama > Ammax)
        Ammax = Ama;
      #endif
    }
    else
      print_err("cannot allocate memory for matrix a(%ld,%ld)", __FILE__, __LINE__, __func__, mat.m, mat.n);
  }

  memcpy(a, mat.a, asize*sizeof(*mat.a));
}



/**
  Copy assignment operator for %matrix class.

  @param mat[in] - %matrix whose values to be assigned to the given instance.

  @author Created by Tomas Koudelka, 05.2023
*/
matrix& matrix::operator=(const matrix &mat)
{
  if (this != &mat){ // check for self-assignment
    long asize = mat.m*mat.n; // mat.size may be greater than actual matrix size, so actual size must be calculated
    if (asize > size){
      if (a && stat.dealloc)  delete [] a;
      
      size = asize;
      a = new double[size];
      if (a){
        stat.dealloc = 1;
        
        print_warning("Copy assignment operator of matrix is called.\n"
                      "Please check your code and make sure you want use an assignment operator.\n"
                      "Consider to use copym function rather", 
                      __FILE__, __LINE__, __func__);

        #ifdef DEBUG_MATRIX
         Acm++;
         Ama += size;
         if (Ama > Ammax)
           Ammax = Ama;
        #endif
      }
      else
        print_err("cannot allocate memory for matrix a(%ld,%ld)", __FILE__, __LINE__, __func__, mat.m, mat.n);
    }
    memcpy(a, mat.a, asize*sizeof(*mat.a));
  }
  return *this;
}



/**
  The function creates %matrix referencing to %matrix src.
  Deallocation of the allocated memory is left to src matrix. 

  @param[in] src - the referenced matrix.
  
  @return 

  @date Created by TKo, 18.9.2024
  @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz
*/
void matrix::makerefm(matrix &src)
{
  if (a && stat.dealloc)
  {
    delete [] a;
    #ifdef DEBUG_MATRIX
    Acm--;
    Ama -= mat.size;
    #endif
  }

  size = src.m*src.n;
  m = src.m;
  n = src.n;
  a = src.a;
  stat.dealloc = 0;
}



/**
  The function creates %matrix which contains reference to array ptr.

  @param[in] mp  - the number of rows stored in ptr
  @param[in] np  - the number of columns stored in ptr
  @param[in] ptr - the pointer to the allocated array of mp x np components which will be referenced by the given %matrix.
                   It is supposed that ptr is the pointer to the array of matrix components stored by rows.

  @date Created by TKo, 18.9.2024
  @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz
*/
void matrix::makerefm(long mp, long np, double *ptr)
{
  if (a && stat.dealloc)
  {
    delete [] a;
    #ifdef DEBUG_MATRIX
    Acm--;
    Ama -= mat.size;
    #endif
  }

  size = mp*np;
  m = mp;
  n = np;
  a = ptr;
  stat.dealloc = 0;
}



/**
   The constructor allocates memory form the heap to the member array a.
   
   @param[in] m is the number of rows
   @param[in] n is the number of columns
   
   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
imatrix :: imatrix(long m, long n)
{
  imatrix::m = m;
  imatrix::n = n;
  size = m*n;
  a = new long[size];
  if (a)
  {
    stat.dealloc = 1;

    memset (a, 0, size*sizeof(*a));

    #ifdef DEBUG_MATRIX
    Acim++;
    Aima += size;
    if (Aima > Aimmax)
      Aimmax = Aima;
    #endif
  }
  else
    print_err("cannot allocate memory for imatrix a(%ld,%ld)", __FILE__, __LINE__, __func__, m, n);
}



/**
   The copy constructor creates copy of the object given by the constant reference parameter.
   
   @param[in] mat is reference to the object of the %imatrix which should be copied
   
   @author Created  16.10.2003 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
imatrix :: imatrix (const imatrix &mat)
{
  long  asize = mat.m*mat.n;

  m = mat.m;
  n = mat.n;
  if (asize > size)
  {
    size = asize;

    if (a && stat.dealloc)
      delete [] a;

    a = new long[size];
    if (a)
    {
      stat.dealloc = 1;

      print_warning("Copy constructor of imatrix is called.\n"
                    "Please check your code and make sure you want use copy constructor.\n"
                    "Consider the usage of copym function.",
                    __FILE__, __LINE__, __func__);

      #ifdef DEBUG_MATRIX
      Acim++;
      Aima += size;
      if (Aima > Aimmax)
        Aimmax = Aima;
      #endif
    }
    else
      print_err("cannot allocate memory for imatrix a(%ld,%ld)", __FILE__, __LINE__, __func__, mat.m, mat.n);
  }

  memcpy(a, mat.a, asize*sizeof(*mat.a));
}



#ifdef DEBUG_MATRIX
 long give_acm()
 {
   return Acm;
 }

 long give_ama()
 {
   return Ama;
 }

 long give_ammax()
 {
   return Ammax;
 }
#endif



/**
   The function allocates memory from heap to the mat's member array a.
   
   @param[in] m is the number of rows
   @param[in] n is the number of columns
   @param[out] mat is the structure for allocated matrix
   
   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long allocm(long m, long n, matrix &mat)
{
  mat.m = m;
  mat.n = n;
  mat.size = m*n;
  mat.a = new double[mat.size];

  if (mat.a)
  {
    memset (mat.a, 0, mat.size*sizeof(*mat.a));
    mat.stat.dealloc = 1;

    #ifdef DEBUG_MATRIX
    Acm++;
    Ama += mat.size;
    if (Ama > Ammax)
      Ammax = Ama;
    #endif
    return (0);
  }
  else
    print_err("cannot allocate memory for matrix a(%ld,%ld)", __FILE__, __LINE__, __func__, mat.m, mat.n);

  return (1);
}



/**
   The function allocates memory from heap to the mat's member array a.
   
   @param[in] m is the number of rows
   @param[in] n is the number of columns
   @param[out] mat is the structure for allocated imatrix
   
   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long allocm(long m, long n, imatrix &mat)
{
  mat.m = m;
  mat.n = n;
  mat.size = m*n;
  mat.a = new long[mat.size];

  if (mat.a)
  {
    memset (mat.a,0,m*n*sizeof(long));
    mat.stat.dealloc = 1;

    #ifdef DEBUG_MATRIX
    Acim++;
    Aima += mat.size;
    if (Aima > Aimmax)
      Aimmax = Aima;
    #endif
    return (0);
  }
  else
    print_err("cannot allocate memory for imatrix a(%ld,%ld)", __FILE__, __LINE__, __func__, mat.m, mat.n);

  return (1);
}



/**
   The function allocates memory from heap to the mat's member array a.
   
   @param[in] m is the number of rows
   @param[in] n is the number of columns
   @param[out] mat is the allocated double array
   
   @author Created  20.11.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
long allocm(long m, long n, double **&mat)
{
  long i;
  
  mat = new double* [m];
  if (mat){
    for (i=0; i<m; i++)
      mat[i] = new double [n];
    return (0);
  }
  else
    print_err("cannot allocate memory for array=matrix a(%ld,%ld)", __FILE__, __LINE__, __func__, m, n);

  return (1);
}



/**
   The function allocates memory from heap to the mat's member array a and
   it is initialized to be identity %matrix.
   
   @param[in] m is the number of rows
   @param[in] n is the number of columns
   @param[out] mat is the structure for allocated %matrix
   
   @author Created  27.11.2014, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long allocim(long n, matrix &mat)
{
  mat.m = n;
  mat.n = n;
  mat.size = n*n;
  mat.a = new double[mat.size];
  if (mat.a)
  {
    mat.stat.dealloc = 1;
    memset (mat.a, 0, mat.size*sizeof(*mat.a)); 
    for (long i=0; i<n; i++)
      mat.a[i*n+i] = 1.0;

     #ifdef DEBUG_MATRIX
     Acm++;
     Ama+=mat.size;
     if (Ama > Ammax)
       Ammax = Ama;
     #endif
     return (0);
  }
  else
    print_err("cannot allocate memory for matrix a(%ld,%ld)", __FILE__, __LINE__, __func__, mat.m, mat.n);

  return (1);
}



/**
   The function uses preallocated memory referenced by ptr for the creation square matrix (nxn) and
   initialization to be identity %matrix.
   
   @param[in] n is the number of rows and columns
   @param[out] mat is the structure for allocated %matrix
   @param[in] dealloc - flag for deallocation of provided memory (0 = no deallocation, 1 = deallocation in destructor)
   @param[in] ptr - pointer to the preallocated memory for the %matrix components
   
   @author Created  27.11.2014, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long allocim(long n, matrix &mat, unsigned dealloc, double *ptr)
{
  if (ptr)
  {
    mat.m = n;
    mat.n = n;
    mat.a = ptr;
    mat.size = n*n;
    mat.stat.dealloc = dealloc;
    memset (mat.a, 0, mat.size*sizeof(*mat.a));
    for (long i=0; i<n; i++)
      mat(i,i) = 1.0;

    #ifdef DEBUG_MATRIX
    Acm++;
    Ama+=mat.size;
    if (Ama > Ammax)
      Ammax = Ama;
    #endif
    return (0);
  }
  else
    print_err("cannot allocate memory for matrix a(%ld,%ld)", __FILE__, __LINE__, __func__, mat.m, mat.n);

  return (1);
}



/**
   The function reallocates memory for the member array a.
   If the memory requirements are less or equal to the 
   allocated one, the memory is not reallocated and the matrix dimensions are changed only. 
   
   @param[in] m is the number of rows
   @param[in] n is the number of columns
   @param[out] mat is the structure for allocated matrix
   
   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long reallocm(long m, long n, matrix &mat)
{
  long asize = m*n;
  
  if (mat.size >= asize)
  {
    mat.m = m;
    mat.n = n;
    memset (mat.a, 0, asize*sizeof(*mat.a));
    return (0);
  }
  else
  {
    if (mat.a && mat.stat.dealloc)
      delete [] mat.a;

    mat.m = m;
    mat.n = n;
    mat.size = asize;
    mat.a = new double[asize];
    if (mat.a)
    {
      mat.stat.dealloc = 1;
      memset (mat.a, 0, asize*sizeof(*mat.a));

      #ifdef DEBUG_MATRIX
      Ama += asize;
      if (Ama > Ammax)
        Ammax = Ama;
      #endif
      return (0);
    }
    else
      print_err("cannot allocate memory for matrix a(%ld,%ld)", __FILE__, __LINE__, __func__, m, n);
  }

  return (1);
}



/**
   The function reallocates memory for the %matrix's member array a.
   If the new %matrix size is greater than previously allocated one, 
   array a is deleted and new memory is allocated with help of pointer ptr which
   references to  preallocated memory for m*n element %matrix. Additionally in this case,
   the dealloc flag specifies whether the ptr was allocated on the stack (dealloc=0) or 
   heap (dealloc=1).
   If the new matrix size is less or equal, then the memory 
   allocated for the array a is left and the numbers of %matrix rows and columns
   are changed. In such case, arguments ptr and dealloc are not used and NULL or 
   zero values should be passed. Components of the array a are set to zero.
   
   @param[in] m   is the new number of the %matrix rows
   @param[in] n   is the new number of the %matrix columns
   @param[out] mat is the structure for reallocated vector
   @param[in] dealloc is the flag for deallocation control (0=no deallocation ,1=release memory by delete operator)
   @param[in] ptr is the pointer to the allocated memory which will be used for storage of n %vector elements
   
   @retval 0 - on success
   @retval 1 - if fails allocating memory
   
   @author Created  5.2015 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long reallocm(long m, long n, matrix &mat, unsigned dealloc, double *ptr)
{
  long asize = m*n;

  if (mat.size >= asize)
  {
    mat.m = m;
    mat.n = n;
    memset (mat.a, 0, asize*sizeof(*mat.a));

    return (0);
  }
  else
  {
    if (ptr || (asize == 0))
    {

      if (mat.a && mat.stat.dealloc)
        delete [] mat.a;

      mat.m = m;
      mat.n = n;
      mat.size = asize;
      mat.a = ptr;
      mat.stat.dealloc = dealloc;
      memset (mat.a, 0, asize*sizeof(*mat.a));

      #ifdef DEBUG_MATRIX
      Aima += asize;
      if (Aima > Aimmax)
        Aimmax = Aima;
      #endif

      return (0);
    }
    else
      print_err("cannot allocate memory for matrix a(%ld,%ld)", __FILE__, __LINE__, __func__, m, n);
  }

  return (1);
}



/**
   The function reallocates memory for the member array a.
   If the memory requirements are less or equal to the 
   allocated one, the memory is not reallocated and the matrix dimensions are changed only. 
   
   @param[in] m is the number of rows
   @param[in] n is the number of columns
   @param[out] mat is the structure for allocated imatrix
   
   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long reallocm(long m, long n, imatrix &mat)
{
  long asize = m*n;

  if (mat.size >= asize)
  {
    mat.m = m;
    mat.n = n;
    memset (mat.a, 0, asize*sizeof(*mat.a));

    return (0);
  }
  else
  {
    if (mat.a && mat.stat.dealloc)
      delete [] mat.a;

    mat.m = m;
    mat.n = n;
    mat.size = asize;
    mat.a = new long[asize];
    if (mat.a)
    {
      memset (mat.a, 0, asize*sizeof(*mat.a));
      mat.stat.dealloc = 1;

      #ifdef DEBUG_MATRIX
      Aima += asize;
      if (Aima > Aimmax)
        Aimmax = Aima;
      #endif

      return (0);
    }
    else
      print_err("cannot allocate memory for imatrix a(%ld,%ld)", __FILE__, __LINE__, __func__, m, n);
  }
  return (1);
}



/**
   The function reallocates memory for the %imatrix's member array a.
   If the new %imatrix size is greater than previously allocated one, 
   array a is deleted and new memory is allocated with help of pointer ptr which
   references to  preallocated memory for m*n element %imatrix. Additionally in this case,
   the dealloc flag specifies whether the ptr was allocated on the stack (dealloc=0) or 
   heap (dealloc=1).
   If the new imatrix size is less or equal, then the memory 
   allocated for the array a is left and the numbers of %matrix rows and columns
   are changed. In such case, arguments ptr and dealloc are not used and NULL or 
   zero values should be passed. Components of the array a are set to zero.
   
   @param[in] m   is the new number of the %imatrix rows
   @param[in] n   is the new number of the %imatrix columns
   @param[out] mat is the structure for reallocated vector
   @param[in] dealloc is the flag for deallocation control (0=no deallocation ,1=release memory by delete operator)
   @param[in] ptr is the pointer to the allocated memory which will be used for storage of n %vector elements
   
   @retval 0 - on success
   @retval 1 - if fails allocating memory
   
   @author Created  5.2015 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long reallocm(long m, long n, imatrix &mat, unsigned dealloc, long *ptr)
{
  long asize = m*n;

  if (mat.size >= asize)
  {
    mat.m = m;
    mat.n = n;
    memset (mat.a, 0, asize*sizeof(*mat.a));
    return (0);
  }
  else
  {
    if (ptr || (asize==0))
    {

      if (mat.a && mat.stat.dealloc)
        delete [] mat.a;

      mat.m = m;
      mat.n = n;
      mat.size = asize;
      mat.a = ptr;
      mat.stat.dealloc = dealloc;
      memset (mat.a, 0, asize*sizeof(*mat.a));

      #ifdef DEBUG_MATRIX    
      Aima += asize;
      if (Aima > Aimmax)
        Aimmax = Aima;
      #endif

      return (0);
    }
    else
      print_err("cannot allocate memory for imatrix a(%ld,%ld)", __FILE__, __LINE__, __func__, m, n);
  }

  return (1);
}



/**
   The function initializes pointer for the matrix's member array a by the parameter
   ptr and thus creates reference to the ptr. Deallocation flag of ref %matrix is set to 0.
   Component array a of the %matrix ref is deallocated if the stat.dealloc flag is set on.
   
   @param[out] ref is the structure of initialized %matrix by ptr pointer
   @param[in] ptr is the pointer to the allocated array of mxn components which will be referenced by ref %matrix.
              It is supposed that ptr is the pointer to the array of matrix components stored by rows.
   @param[in] m   is the number of ref/ptr rows
   @param[in] n   is the number of the ref/ptr columns
   
   @retval 0 - on success
   
   @author Created  3.2019 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long makerefm(matrix &ref, double *ptr, long m, long n)
{

  if (ref.a && ref.stat.dealloc)
  {
    delete [] ref.a;
    #ifdef DEBUG_MATRIX
    Acm--;
    Ama -= mat.size;
    #endif
  }

  ref.size = m*n;
  ref.m = m;
  ref.n = n;
  ref.a = ptr;
  ref.stat.dealloc = 0;

  return (0);
}



/**
   The function creates reference %matrix ref to the source %matrix src.
   The ref %matrix is initialized by the content of src %matrix but only SHALOW
   copy of data members is performed except of deallocation flag of ref %matrix which is set to 0.
   Component array a of the %matrix ref is deallocated if the stat.dealloc flag is set on.
   
   @param[out] ref is the structure for %matrix used for storage of new refrence to src
   @param[in] src structure with the source %matrix which will be referenced by ref 
   
   @retval 0 - on success
   
   @author Created  3.2019 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long makerefv(matrix &ref, matrix &src)
{

  if (ref.a && ref.stat.dealloc)
  {
    delete [] ref.a;
    #ifdef DEBUG_MATRIX
    Acm--;
    Ama -= mat.size;
    #endif
  }

  ref.size = src.size;
  ref.m = src.m;
  ref.n = src.n;
  ref.a = src.a;
  ref.stat.dealloc = 0;

  return (0);
}



/**
   The function copies %matrix given by src to dest.
   
   @param[in] src is the structure of source %matrix to copy
   @param[out] dest is the structure of destination %matrix to which will be copied contents of src
   
   @b Requests :
   dest has to be setuped dimensions and allocated memory array for elements which is enough
   large to hold all contents of the %matrix src.
   
   @retval 0 : on succes
   @retval 1 : in case incompatibility sizes of the src and dest matrices
   
   @author created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long copym(const matrix &src, matrix &dest)
{
  if ((src.m == dest.m) && (src.n == dest.n))
  {
    memcpy(dest.a, src.a, sizeof(*src.a)*src.m*src.n);
    return (0);
  }
  else
    print_err("cannot copy matrix - incompatible size of matrices\n"
              " src(%ld,%ld) X dest(%ld,%ld)", __FILE__, __LINE__, __func__, src.m, src.n, dest.m, dest.n);

  return (1);
}



/**
   The function copies matrix given by src to dest.
   
   @param[in] src is the structure of source %matrix to copy
   @param[out] dest is pointer to the allocated array for destination %matrix to which will be copied contents of src
   
   @b Requests :
   dest has to be setuped dimensions and allocated memory array for elements which is enough
   large to hold all contents of the matrix src.
   
   @retval 0 : on succes
   
   @author Created  5.3.2019, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long copym(const matrix &src, double *dest)
{
  memcpy(dest, src.a, sizeof(*src.a)*src.m*src.n);
  return (0);
}



/**
   The function copies matrix given by src to dest.
   
   @param[in] src is pointer to the array with source %matrix, 
   @param[out] dest is the structure of destination %matrix to which will be copied contents of src
   
   @b Requests :
   src has to be the same dimensions as src %matrix
   
   @retval 0 : on succes
   
   @author Created  22.5.2024, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long copym(const double *src, matrix &dest)
{
  memcpy(dest.a, src, sizeof(*src)*dest.m*dest.n);
  return (0);
}



/**
   The function fills memory of mat's member array a with zeros.
   
   @param[out] mat is the given matrix for the zeroing.
   
   @return always zero
   
   @author Created  22.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long nullm(matrix &mat)
{
  memset(mat.a, 0, sizeof(*mat.a)*mat.m*mat.n);
  return (0);
}



/**
   The function fills memory of mat's member array a with constant c.
   
   @param[in]  c is constant, which will be used for filling entire matrix
   @param[out] mat is the structure for the filled matrix
   
   @return always zero
   
   @author Created  22.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long fillm(double c, matrix &mat)
{
  long i, j;
  for (i = 0; i < mat.m; i++)
    for (j = 0; j < mat.n; mat[i][j] = c, j++);
  return (0);
}



/**
 The function fills mat's row i with the constant c.

  @param[in] c   is constant, which will be used for filling entire matrix,
  @param[in] i   is index of given row,
  @param[out] mat is the given matrix for the filling.

  @return always zero

  @author Created  22.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long fillrow(double c, long i, matrix &mat)
{
  long j;

  if ((i < mat.m) && (i >= 0))
  {
    for (j=0; j < mat.n; j++)
      mat[i][j] = c;

    return (0);
  }
  else
    print_err("matrix row index %ld is out of range <0,%ld>", 
              __FILE__, __LINE__, __func__, i, mat.m-1);

  return (1);
}



/**
 The function fills mat's column i with the constant c.

  @param[in] c   is constant, which will be used for filling entire matrix,
  @param[in] i   is index of given column,
  @param[out] mat is the given matrix for the filling.

  @return always zero

  @author Created  22.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long fillcol(double c, long i, matrix &mat)
{
  long j;

  if (( i < mat.n) && (i >= 0))
  {
    for (j=0; j < mat.m; j++)
      mat[j][i] = c;

    return (0);
  }
  else
    print_err("matrix column index %ld is out of range <0,%ld>", 
              __FILE__, __LINE__, __func__, i, mat.n-1);

  return (1);
}



/**
   The function sets the given %matrix to the identity one.
   
   @param[out] mat is the given matrix which becomes identity one.
   
   @return always zero
   
   @author Created  11.12.2014, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long identm(matrix &mat)
{
  long i;

  memset(mat.a, 0, mat.m*mat.n*sizeof(*mat.a));
  for (i = 0; i < mat.m; i++)
    mat(i,i) = 1.0;

  return (0);
}



/**
   The function deallocates memory occupied by mat's member array a.
   
   @param[out] mat is the given deallocated matrix
   
   @return always zero
   
   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long destrm(matrix &mat)
{
 #ifdef DEBUG_MATRIX
  if (mat.a != NULL)
  {
    Acm--;
    Ama -= mat.size;
  }
 #endif

  if (mat.stat.dealloc)
    delete [] mat.a;

  mat.m = 0;
  mat.n = 0;
  mat.a = NULL;
  mat.size = 0;
  return (0);
}



/**
   The function deallocates memory occupied by double array mat.
   
   @param[in] m - the number of %matrix rows,
   @param[out] mat - the given deallocated array.
   
   @return always zero
   
   @author Created  5.5.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
long destrm(double **&mat, long m)
{
  for (long i=0; i<m; i++)
    delete [] mat[i];
  
  delete [] mat;
  mat = NULL;
  return (0);
}



/**
   The function deallocates memory occupied by long array mat.
   
   @param[in] m - the number of %matrix rows,
   @param[out] mat - the given deallocated array.
   
   @return always zero
   
   @author Created  5.5.2003, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
long destrm(long **&mat, long m)
{
  for (long i=0; i<m; i++)
    delete [] mat[i];
  
  delete [] mat;
  mat = NULL;
  return (0);
}



/**
   The function adds matrix given by a to matrix given by b, the result is stored in c.
   
   @param[in] a is the structure of the first added matrix
   @param[in] b is the structure of the second added matrix
   @param[out] c is the structure of the result matrix
   
   @b Requests :
   a, b and c have to be same dimensions, c has to have allocated memory array for elements
   which is enough large to hold contents of the result.
   
   @retval 0 : on succes.
   @retval 1 : in case incompatibility sizes of a, b and c matrices.
   
   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  5.9.2016, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long addm(const matrix &a, const matrix &b, matrix &c)
{
  long i,j;
  
  if ((a.m == b.m) && (a.n == b.n) &&
      (c.m == b.m) && (c.n == b.n))
  {
    for (i = 0; i < a.m; i++)
      for (j = 0; j < a.n; j++)
        c(i,j) = a(i,j) + b(i,j);

    return(0);
  }
  else
    print_err("cannot add matrices due to their incompatible dimensions\n"
              " a(%ld,%ld) X b(%ld,%ld) X c(%ld,%ld)", __FILE__, __LINE__, __func__, a.m, a.n, b.m, b.n, c.m, c.n);

  return (1);
}



/**
   The function subtracts matrix given by b from matrix given by a, the result is stored in the c
   
   @param[in] a is the structure of the matrix, from which is subtracted
   @param[in] b is the structure of subtracted matrix
   @param[out] c is the structure of the result matrix
   
   @b Requests :
   a, b and c have to be same dimensions, c has to have allocated memory array for elements
   which is enough large to hold contents of the result.
   
   @retval 0 : on succes
   @retval 1 : in case incompatibility sizes of a, b and c matrices
   
   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long subm(const matrix &a, const matrix &b, matrix &c)
{
  long i,j;
  
  if ((a.m == b.m) && (a.n == b.n) &&
      (b.m == c.m) && (b.n == c.n))
  {
    for (i = 0; i < a.m; i++)
      for (j = 0; j < a.n; j++)
        c(i,j) = a(i,j) - b(i,j);

    return(0);
  }
  else
    print_err("cannot subtract matrices due to their incompatible dimensions\n"
              " a(%ld,%ld) X b(%ld,%ld) X c(%ld,%ld)", __FILE__, __LINE__, __func__, a.m, a.n, b.m, b.n, c.m, c.n);

  return (1);
}



/**
  The function adds 2 matrices multplied by constants to the matrix C: C := C + ac*A + bc*B.

  @param[in] a - the first %matrix adder,
  @param[in] ac -  a scaling factor of the first %matrix adder,
  @param[in] b - the second %matrix adder,
  @param[in] bc -  a scaling factor of the second %matrix adder,
  @param[in,out] c - the resulting sum %matrix

  @retval 0 : on succes
  @retval 1 : in case incompatibility sizes of a, b and c matrices
   
  @author Created  02.2024, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long addmultm(const matrix &a, double ac, const matrix &b, double bc, matrix &c)
{
  if ((a.m == b.m) && (a.n == b.n) &&
      (c.m == b.m) && (c.n == b.n))
  {
    for (long i = 0; i < a.m; i++)
      for (long j = 0; j < a.n; j++)
        c(i,j) += ac*a(i,j) + bc*b(i,j);

    return(0);
  }
  else
    print_err("cannot add matrices due to their incompatible dimensions\n"
              " a(%ld,%ld) X b(%ld,%ld) X c(%ld,%ld)", __FILE__, __LINE__, __func__, a.m, a.n, b.m, b.n, c.m, c.n);

  return (1);  
}



/**
  The function performs tensor product of two vectors, i.e. it performs C[m,n] = a[m] x b[n].
   
   @param[in] a is the structure of the first vector
   @param[in] b is the structure of the second vector
   @param[out] c is the structure of the result matrix
   
   @b Requests :
   c has to have allocated memory array for elements
   which is enough large to hold contents of the result.
   
   @retval 0 : on succes
   @retval 1 : in case incompatibility sizes of a, b and c matrices
   
   @author Created by ?
*/
long tensprd (const vector &a, const vector &b, matrix &c)
{
  long i,j;
  
  if ((a.n == c.m) || (b.n == c.n))
  {
    for (i=0;i<a.n;i++)
      for (j=0;j<b.n;j++)
        c(i,j)=a[i]*b[j];
  
    return(0);
  }
  else
    print_err("cannot make tensor product of vectors due to incompatible dimensions\n"
              " a(%ld) X b(%ld) X c(%ld,%ld)", __FILE__, __LINE__, __func__, a.n, b.n, c.m, c.n);

  return (1);
}



/**
   The function multiplies matrix given by a from right with matrix given by b,
   the result is stored in c, i.e. it performs C[m,p] := A[m,n].B[n,p] 
   
   @param[in] a is the structure of the matrix, from which is multiplied
   @param[in] b is the structure of the multiplicating matrix
   @param[out] c is the structure of the result matrix
   
   @b Requests :
   a, b and c have to be following dimensions
   a (m,n),  b (n,p),  c (m,p)
   c has to have allocated memory array for elements
   which is enough large to hold contents of the result.
   
   @retval 0 : on succes
   @retval 1 : in case incompatibility sizes of a, b and c matrices
   
   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long mxm(const matrix &a, const matrix &b, matrix &c)
{
  long i,j,k;
  double aux;
  
  if ((a.n == b.m) &&
      (c.m == a.m) && (c.n == b.n))
  {
    for (i=0; i < a.m; i++)
    {
      for (j=0; j < b.n; j++)
      {
        aux = 0.0;
        for (k=0; k < a.n; k++)
          aux += a(i,k) * b(k,j);
        c(i,j) = aux;
      }
    }

    return(0);
  }
  else
    print_err("cannot multiply matrices due to their incompatible dimensions\n"
              " a(%ld) X b(%ld) X c(%ld,%ld)", __FILE__, __LINE__, __func__, a.m, a.n, b.m, b.n, c.m, c.n);

  return (1);
}



/**
   The function multiplies matrix given by 1D double array a from right
   with matrix given by 1D double array b, the result is stored in 1D double array c.
   A[l,m] . B[m,n] = C[l,n]
   
   @param[in] a is the 1D double array of first matrix, from which is multiplied
   @param[in] b is the 1D double array of the second=multiplicating matrix
   @param[out] c is the 1D double array of the result matrix
   
   @b Requests :
   a, b and c have to be following dimensions
   a (l*m),  b (m*n),  c (l*n)
   
   @retval 0 : always
   
   @author Created  19.2.1997 by JK
*/
void mxm(const double *a, const double *b, double *c, long l, long m, long n)
{
  long i,j,k,ac,acb,acu,acl;
  double s;

  acl=0;  ac=0;
  for (i=0;i<l;i++){
    acu=acl+m;
    for (j=0;j<n;j++){
      s=0.0;
      acb=j;
      for (k=acl;k<acu;k++){
        s+=a[k]*b[acb];
        acb+=n;
      }
      c[ac]=s;
      ac++;
    }
    acl+=m;
  }
}



/**
   The function multiplies matrix given by a from right with transposed matrix given by b,
   the result is stored in c, i.e. it performs C[m,p] := A[m,n].B[p,n]^T
   
   @param[in] a is the structure of the matrix, which will be multiplied form right
   @param[in] b is the structure of the multiplicating matrix which will be transposed
   @param[out] c is the structure of the resulting matrix
   
   @b Requests :
   a, b and c have to be following dimensions
   a (m,n),  b (p,n),  c (m,p)
   c has to have allocated memory array for elements
   which is enough large to hold contents of the result.
   
   @retval 0 : on succes
   @retval 1 : in case incompatibility sizes of a, b and c matrices
   
   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long mxmt(const matrix &a, const matrix &b, matrix &c)
{
  long i,j,k;
  double aux;
  
  if ((a.n == b.n) &&
      (c.m == a.m) && (c.n == b.m))
  {
    for (i=0; i < a.m; i++)
    {
      for (j=0; j < b.m; j++)
      {
        aux = 0.0;
        for (k=0; k < a.n; k++)
          aux += a(i,k) * b(j,k);
        c(i,j) = aux;
      }
    }

    return(0);
  }
  else
    print_err("cannot multiply matrices due to their incompatible dimensions\n"
              " a(%ld,%ld) X b(%ld,%ld) X c(%ld,%ld)", __FILE__, __LINE__, __func__, a.m, a.n, b.m, b.n, c.m, c.n);

  return (1);
}



/**
   The function multiplies matrix given by 1D double array a from right
   with transposed matrix given by 1D double array b, the result is stored in c
   A[l,m] . B[m,n]^T = C[l,n]
      
   @param[in] a is the 1D double array of first matrix, from which is multiplied
   @param[in] b is the 1D double array of the second=multiplicating matrix
   @param[out] c is the 1D double array of the result matrix
   
   @b Requests :
   a, b and c have to be following dimensions
   a (l*m),  b (m*n),  c (l*n)
   
   @retval 0 : always
   
   @author Created  9.11.1999 by JK
*/
void mxmt(const double *a, const double *b, double *c, long l, long m, long n)
{
  long i,j,k,ii,kk,lk,uk;
  double s;
  
  ii=0;
  for (i=0;i<l;i++){
    lk=i*m;  uk=lk+m;
    for (j=0;j<n;j++){
      s=0.0;  kk=j*m;
      for (k=lk;k<uk;k++){
        s+=a[k]*b[kk];  kk++;
      }
      c[ii]=s;  ii++;
    }
  }
}



/**
   The function multiplies transposed matrix given by a from right with matrix given by b,
   the result is stored in c, i.e. it performs C[m,p] := A[n,m]^T . B[n,p]
   
   @param[in] a is the structure of the matrix, which will be transposed and multiplied from the right
   @param[in] b is the structure of the multiplicating matrix
   @param[out] c is the structure of the result matrix
   
   @b Requests :
   a, b and c have to be following dimensions
   a (n,m),  b (n,p),  c (m,p)
   c has to have allocated memory array for elements
   which is enough large to hold contents of the result.
   
   @retval 0 : on succes
   @retval 1 : in case incompatibility sizes of a, b and c matrices
   
   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long mtxm(const matrix &a, const matrix &b, matrix &c)
{
  long i,j,k;
  double aux;
  
  if ((a.m == b.m) &&
      (c.m == a.n) && (c.n == b.n))
  {
    for (i=0; i < a.n; i++)
    {
      for (j=0; j < b.n; j++)
      {
        aux = 0.0;
        for (k=0; k < a.m; k++)
          aux += a(k,i) * b(k,j);
        c(i,j) = aux;
      }
    }

    return(0);
  }
  else
    print_err("cannot multiply matrices due to their incompatible dimensions\n"
              " a(%ld,%ld) X b(%ld,%ld) X c(%ld,%ld)", __FILE__, __LINE__, __func__, a.m, a.n, b.m, b.n, c.m, c.n);

  return (1);
}



/**
   The function multiplies transposed matrix given by 1D double array a from right
   with matrix given by 1D double array b, the result is stored in c
   A[l,m]^T . B[m,n] = C[l,n]
      
   @param[in] a is the 1D double array of first matrix, from which is multiplied
   @param[in] b is the 1D double array of the second=multiplicating matrix
   @param[out] c is the 1D double array of the result matrix
   
   @b Requests :
   a, b and c have to be following dimensions
   a (l*m),  b (m*n),  c (l*n)
   
   @retval 0 : always
   
   @author Created  13.4.1998 by JK
*/
void mtxm (const double *a, const double *b, double *c, long l, long m, long n)
{
  long i,j,k,aca,acb,acc;
  double s;
  
  acc=0;
  for (i=0;i<m;i++){
    for (j=0;j<n;j++){
      s=0.0;  aca=i;  acb=j;
      for (k=0;k<l;k++){
        s+=a[aca]*b[acb];
        aca+=m;  acb+=n;
      }
      c[acc]=s;  acc++;
    }
  }
}



/**
  The function computes %matrix product C := A^T . B,
  all matrices are stored by columns.

  @param[in] a - array with %matrix A(m,n)
  @param[in] b - array with %matrix B(m,p)
  @param[out] c - array with %matrix C(n,p)
  @param[in] l -
  @param[in] m -
  @param[in] n - 

  @author Created by JK, 7.12.1998
*/
void mtxmccr (const double *a, const double *b, double *c, long l, long m, long n)
{
  long i,j,k,ii,lk,uk,ac;
  double s;
  
  ac=0;
  for (i=0;i<m;i++){
    for (j=0;j<n;j++){
      s=0.0;  lk=i*l;  uk=lk+l;  ii=j*l;
      for (k=lk;k<uk;k++){
        s+=a[k]*b[ii];  ii++;
      }
      c[ac]=s;  ac++;
    }
  }
}



/**
   The function multiplies %matrix given by a by constant c,
   the result is stored in b, i.e. it performs B[m,n] := c A[m,n]
   
   @param[in] c is the real number
   @param[in] a is the given multiplied %matrix
   @param[out] b is the resulting %matrix
   
   @b Requests :
   a and b have to be same dimensions, b has to have allocated memory array for elements
   which is enough large to hold contents of the result.
   
   @retval 0 : on succes
   @retval 1 : in case incompatibility sizes of a and b matrices
   
   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long cmulm(const double c, const matrix &a, matrix &b)
{
  long i,j;
  
  if ((a.m == b.m) || (a.n == b.n))
  {
    for (i=0; i < a.m; i++)
    {
      for (j=0; j < a.n; j++)
        b(i,j) = c * a(i,j);
    }

    return (0);
  }
  else
    print_err("incompatible dimensions of  matrices - a(%ld,%ld) X b(%ld,%ld)",
              __FILE__, __LINE__, __func__, a.m, a.n, b.m, b.n);

  return (1);
}   



/**
   The function multiplies %matrix given by a with constant c,
   the result is stored in a, i.e. it performs A[m,n] := c A[m,n]
   
   @param[in] c is the real number
   @param[out] a is the given multiplied %matrix
   
   @return always zero
   
   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long cmulm(double c, matrix &a)
{
  long i,j;
  
  for (i=0; i < a.m; i++)
  {
    for (j=0; j < a.n; j++)
      a(i,j) *= c;
  }
  
  return (0);
}



/**
   The function transposes %matrix given by a,
   the result is stored in at, i.e. it performs A_t[n,m] := A^T[m,n]
   
   @param[in] a  is the %matrix being transposed,
   @param[out] at is the resulting %matrix.
   
   @b Requests :
   For a(m,n), dimensions of at have to be (n,m), at has to have allocated memory array for elements
   which is enough large to hold contents of the result.
   
   @retval 0 : on succes
   @retval 1 : in case incompatibility sizes of a and at matrices
   
   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long tranm(const matrix &a, matrix &at)
{
  long i,j;
  
  if ((a.m == at.n) && (a.n == at.m))
  {
    for (i = 0; i < at.m; i++)
    {
      for (j = 0; j < at.n; j++)
        at(i,j) = a(j,i);
    }

    return(0);
  }
  else
    print_err("cannot transpose matrix - incompatible dimensions of  matrices\n"
              " a(%ld,%ld) X at(%ld,%ld)", __FILE__, __LINE__, __func__, a.m, a.n, at.m, at.n);

  return (1);
}



/**
   The function transposes %matrix given by a(m,n),
   the result is stored in a whose dimensions are changed to (n,m).
   
   @param[in,out] a  - the given matrix being transposed
   
   @retval always zero
   
   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long tranm(matrix &a)
{
  long i,j;
  matrix temp(a.n, a.m);

  for (i = 0; i < temp.m; i++)
  {
    for (j = 0; j < temp.n; j++)
      temp(i,j) = a(j,i);
  }
  memcpy(a.a, temp.a, a.m*a.n*sizeof(*a.a));
  long t = a.n;
  a.n = a.m;
  a.m = t;    
  return(0);
}



/**
   The function multiplies row %vector given by u with %matrix given by a,
   the result is stored in row %vector v, i.e. it performs v[n]^T := u[m]^T.A[m,n]
   
   @param[in] u - the first factor, i.e. row %vector
   @param[in] a - the second factor, i.e. %matrix
   @param[out] v - the resulting row %vector
   
   @b Requests :
   u ,a and v have to be following dimensions
   u(m),  a(m,n),  v(n)
   v has to have allocated memory array for elements
   which is enough large to hold contents of the result.
   
   @retval 0 : on succes
   @retval 1 : in case incompatibility sizes of the matrix a
   
   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long vxm(const vector &u, const matrix &a, vector &v)
{
  long i,j;
  double aux;
  
  if ((u.n == a.m) && (v.n == a.n))
  {
    for (i=0; i < v.n; i++)
    {
      aux = 0.0;
      for (j=0; j < a.m; j++)
        aux += u[j] * a(j,i);
      v[i] = aux;
    }

    return(0);
  }
  else
    print_err("cannot multiply vector by matrix due to their incompatible dimensions\n"
              " u(%ld) X a(%ld,%ld) X v(%ld)", __FILE__, __LINE__, __func__, u.n, a.m, a.n, v.n);

  return (1);
}



/**
   The function multiplies row vector given by u with transposed matrix given by a,
   the result is stored in row vector v, i.e. it performs v[m]^T := u[n]^T.A[m,n]^T
   
   @param[in] u - the first factor, i.e. a row %vector,
   @param[in] a - the second factor, i.e. %matrix to be transposed,
   @param[out] v - the resulting row vector.
   
   @b Requests :
   u ,a and v have to be following dimensions
   u(n),  a(m,n),  v(m)
   v has to have allocated memory array for elements
   which is enough large to hold contents of the result.
   
   @retval 0 : on succes
   @retval 1 : in case incompatibility sizes of the matrix a
   
   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long vxmt(const vector &u, const matrix &a, vector &v)
{
  long i,j;
  double aux;
  
  if ((u.n == a.n) && (v.n == a.m))
  {
    for (i=0; i < a.m; i++)
    {
      aux = 0.0;
      for (j=0; j < a.n; j++)
        aux += u[j] * a(i,j);
      v[i] = aux;
    }

    return(0);
  }
  else
    print_err("cannot multiply vector by matrix due to their incompatible dimensions\n"
              "u(%ld) X a(%ld,%ld) X v(%ld)", __FILE__, __LINE__, __func__, u.n, a.m, a.n, v.n);

  return (1);
}



/**
   The function multiplies %matrix given by a with column %vector u,
   the result is stored in column %vector v, i.e. it performs v[m] := A[m,n].u[n]
   
   @param[in] a - the first factor, i.e. a %matrix,
   @param[in] u - the second factor, i.e. a column %vector,
   @param[out] v is the resulting column vector.
   
   @b Requests :
   u ,a and v have to be following dimensions
   v(m),  a(m,n),  u(n)
   v has to have allocated memory array for the elements
   which is enough large to hold contents of the result.
   
   @retval 0 : on succes
   @retval 1 : in case incompatibility sizes of the %matrix a
   
   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long mxv(const matrix &a, const vector &u, vector &v)
{
  long i,j;
  double aux;
  
  if ((a.n == u.n) && (v.n == a.m))
  {
    for (i=0; i < a.m; i++)
    {
      aux = 0.0;
      for (j=0; j < a.n; j++)
        aux += u[j] * a(i,j);
      v[i] = aux;
    }

    return(0);
  }
  else
    print_err("cannot multiply matrix by vector due to their incompatible dimensions\n"
              " a(%ld,%ld) X u(%ld) X v(%ld)",__FILE__, __LINE__, __func__, a.m, a.n, u.n, v.n);

  return (1);
}



/**
  The function multiplies matrix stored in 1D array a with vector b stored in 1D array,
  i.e. it performs   c[m] :=  A[m,n].b[n]

  @param[in] a - array with %matrix A components stored by rows
  @param[in] b - array with vector b components
  @param[out] c - result vector
  @param[in] m - number of rows in matrix A, dimension of c
  @param[in] n - number of columns in matrix A, dimension of b

  @author Created 19.2.1997 by Jaroslav Kruis
*/
void mxv (const double *a, const double *b, double *c, long m, long n)
{
  long i,j,aca;
  double s;

  aca=0;
  for (i=0;i<m;i++){
    s=0.0;
    for (j=0;j<n;j++){
      s+=a[aca]*b[j];
      aca++;
    }
    c[i]=s;
  }
}



/**
  The function multiplies %matrix a stored by columns in the 
  array by the %vector b stored in array. The resulting %vector
  is stored in the array c. (c := A.b)

  @param[in] a - array with %matrix components stored by columns
  @param[in] b - array with components of %vector b
  @param[in] m - the number of %matrix rows
  @param[in] n - the number of %matrix/%vector columns 
  @param[out] c - array of the resulting %vector
   
  @author Created 15.12.1998 by Jaroslav Kruis
*/
void mxvc (const double *a, const double *b, double *c, long m, long n)
{
  long i,j,k;
  double s;
  
  for (i=0;i<m;i++){
    s=0.0;  k=i;
    for (j=0;j<n;j++){
      s+=a[k]*b[j];  k+=m;
    }
    c[i]=s;
  }
}



/** 
  The function computes i-th component of the %vector v given by product 
  of %matrix a and %vector u  v := A.u
   
  @param[in] a - structure with the multiplied %matrix,
  @param[in] u - structure with the column %vector,
  @param[out] vi - i-th component of the resulting %vector v,
  @param[in] i - index of required %vector component.
   
   @b Requests :
   u ,a have to be following dimensions
   a (m,n),  u (n)
   
   @retval 0 : on succes
   @retval 1 : in case incompatibility sizes of the %matrix a and %vector u.
   
   @author Created  17.12.2014, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long mixv(const matrix &a, const vector &u, double &vi, long i)
{
  long j;
  double aux;

  if (a.n == u.n)
  {
    aux = 0.0;
    for (j=0; j<a.n; j++)
      aux += a(i,j)*u[j];
    vi = aux;

    return 0;
  }
  else
    print_err("cannot multiply matrix by vector due to their incompatible dimensions\n"
              " a(%ld,%ld) X u(%ld)",__FILE__, __LINE__, __func__, a.m, a.n, u.n);

  return (1);
}



/**
   The function multiplies transposed %matrix given by a with column %vector u,
   the result is stored in column %vector v, i.e. it performs v[n] := A[m,n]^T.u[m]
   
   @param[in] a is the first factor, i.e. %matrix to be transposed and multiplied,
   @param[in] u is the second factor, i.e. column %vector,
   @param[out] v is the resulting column %vector.
   
   @b Requests :
   u ,a and v have to be following dimensions
   u(m),  a(m,n),  v(n)
   v has to have allocated memory array for the elements
   which is enough large to hold contents of the result.
   
   @retval 0 : on succes
   @retval 1 : in case incompatibility sizes of the matrix a
   
   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long mtxv(const matrix &a, const vector &u, vector &v)
{
  long i,j;
  double aux;
  
  if ((a.m == u.n) && (v.n == a.n))
  {
    for (i=0; i < v.n; i++)
    {
      aux = 0.0;
      for (j=0; j < a.m; j++)
        aux += u[j] * a(j,i);
      v[i] = aux;
    }

    return(0);
  }
  else
    print_err("cannot multiply matrix by vector due to their incompatible dimensions\n"
              " a(%ld,%ld) X u(%ld) X v(%ld)", __FILE__, __LINE__, __func__, a.m, a.n, u.n, v.n);

  return (1);
}



/**
   The function multiplies transposed %matrix by %vector from right, i.e. c:= A^T.b.
   
   @param[in] a(m,n) - %matrix A given by array,
   @param[in] b(m,1) - %vector b given by array,
   @param[out] c(n,1) - the resulting %vector c given by array.
   @param[in] m - number of %matrix rows
   @param[in] n - number of %matrix columns, number of %vector components
   
   @author Created by JK, 19.2.1997
*/
void mtxv(const double *a, const double *b, double *c, long m, long n)
{
  long i,j,aca;
  double s;
  
  for (i=0;i<n;i++){
    s=0.0;  aca=i;
    for (j=0;j<m;j++){
      s+=a[aca]*b[j];
      aca+=n;
    }
    c[i]=s;
  }
}



/**
   The function multiplies transposed %matrix a by the %vector from right, i.e. c:= A^T.b,
   The %matrix a is stored by columns.

   @param[in] a - array with %matrix a(m,n) stored by columns
   @param[in] b - array with %vector b(m)
   @param[in] m - number of %matrix rows / %vector b components
   @param[in] n - number of %matrix columns / %vector c components
   @param[out] c - array with %vector c(n)

   @author Created 7.12.1998 by JK
*/
void mtxvc(const double *a, const double *b, double *c, long m, long n)
{
  long i,j,k;
  double s;
  
  for (i=0;i<n;i++){
    s=0.0;  k=i*m;
    for (j=0;j<m;j++){
      s+=a[k]*b[j];  k++;
    }
    c[i]=s;
  }
}



/**
   The function multiplies column %vector given by u with row %vector given by v,
   the result is stored in %matrix a, i.e. a[m,n] := u[m] x v[n]^T
   
   @param[in] u - the first factor, i.e. a column %vector,
   @param[in] v - the second factor, i.e. a row %vector,
   @param[out] a - the resulting %matrix.
   
   @b Requests :
   u, a and v have to be following dimensions
   u(m),  v(n), a(m,n)
   a has to have allocated memory array for elements
   which is enough large to hold contents of the result.
   
   @retval 0 : on succes
   @retval 1 : in case incompatibility sizes of the %matrix a
   
   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long vxv(const vector &u, const vector &v, matrix &a)
{
  long i,j;
  
  if ((u.n == a.m) && (v.n == a.n))
  {
    for (i=0; i < a.m; i++)
    {
      for (j=0; j < a.n; j++)
        a(i,j)=u[i]*v[j];
    }

    return(0);
  }
  else
    print_err("cannot perform tensor product of vectors due to their incompatible dimensions\n"
              " u(%ld) X v(%ld) X a(%ld,%ld)", __FILE__, __LINE__, __func__, u.n, v.n, a.m, a.n);

  return (1);
}



/**
   The function multiplies row %vector given by "v" with symmetric %matrix given by "m" and
   with column %vector given by "v" -> v.m.v , the result is stored in "answer".

   @param[in] v is the given %vector to be multiplied,
   @param[in] m is the given %matrix to be multiplied,
   @param[out] answer is the resulting real number.
   
   @b Requests :
   "v" and "m" have to be following dimensions: v (n),  m (n,n)

   @retval 0 : on succes
   @retval 1 : in case incompatibility sizes of the matrix a
   
   @author Created  25.10.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
long vxmxv (const vector &v, const matrix &m, double &answer)
{
  long i,j;
  double aux;
  
  if ((v.n == m.n) && (v.n == m.m)) 
  {
    aux = 0.0;

    for (i=0;i<v.n-1;i++)
      for (j=i+1;j<v.n;j++)
        aux += 2*v[i]*m(i,j)*v[j];

    for (i=0;i<v.n;i++)
      aux += v[i]*m(i,i)*v[i];

    answer = aux;
    return(0);
  }
  else
    print_err("cannot vxmxv cannot be performed due to incompatible dimensions\n"
              "v(%ld) X m(%ld,%ld)", __FILE__, __LINE__, __func__, v.n, m.m, m.n);

  return (1);
}



/**
   The function multiplies row vector given by "v" with symmetric %matrix given by "m"
   with column %vector given by "v" -> v.m.v , the result is stored in "answer".

   @param[in] v is pointer on double array,
   @param[in] m is pointer on double array,
   @param[in] dim is dimensions of "v" and "m",
   @param[out] answer is the resulting real number.
   
   @b Requests :
   "v" and "m" have to be following dimensions: v (dim),  m (dim,dim)
   
   @retval 0 : always
    
   @author Created  25.10.2002, Ladislav Svoboda, termit@cml.fsv.cvut.cz
*/
long vxmxv (const double *v, const double *m, long dim, double &answer)
{
  long i,j;
  
  answer = 0.0;

  for (i=0;i<dim-1;i++)
    for (j=i+1;j<dim;j++)
      answer += 2*v[i]*m[i*dim+j]*v[j];

  for (i=0;i<dim;i++)
    answer += v[i]*m[i*dim+i]*v[i];
  
  return(0);
}


/**
   The function performs Gauss elimination on the %matrix a,
   the result is stored in the %matrix b.
   
   @param[in] a - elimanted %matrix,
   @param[out] b - the resulting %matrix,
   @param[in] zero - computer zero.
   
   @b Requests :
   a has to have the same number of the rows and columns
   b has to have same dimensions as the matrix a
   b has to have allocated memory array for the elements
     which is enough large to hold contents of the result.
   
   @retval 0 : on succes
   @retval 1 : in case incompatibility sizes of a and b matrices
   @retval 2 : in case dimension of the %matrix a  <  2
   @retval 3 : in case some diagonal elements is to small or zero
               the value which is assumed as zero is driven by macro
               USREPS, which is defined in the matrix.h
               
   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long gause(const matrix &a, matrix &b,double zero)
{
  long i,j,k;

  if ((a.m != b.m) || (a.n != b.n) || (a.m != a.n))
  {
    print_err("Gauss elimination - incompatible dimension of matrices"
              "a(%ld,%ld) X b(%ld,%ld)", __FILE__, __LINE__, __func__, a.m, a.n, b.m, b.n);
    return (1);
  }
  if (a.m < 2)
  {
    print_err("matrix has less than 2 rows", __FILE__, __LINE__, __func__);
    return (2);
  }

  for (i = 0; i < a.n; i++)
    b(0,i) = a(0,i);  // copying of the first row of the a to the b

  for (i = 0; i < a.m-1; i++)
  {
    for (j = i+1; j < a.m; j++)
    {
      // testing for the zero on the diagonal
      if (fabs(b(i,i)) <= zero)
      {
        print_err("diagonal element %ld is too small or zero (%le)", __FILE__, __LINE__, __func__, i+1, b(i,i));
        return (3);
      }
      if ((i == 0) && (fabs(a(j,i)) <= zero))
      // testing for the zero under the diagonal in the current row in case of the first step
      {
        // zero is hit so we have to copy this row from the a to the b to
        // fill empty b's rows
        for (k = 0; k < b.n; k++)
          b(j,k) = a(j,k);
        continue; // we may continue on the next step
      }
      if ((i != 0) && (fabs(b(j,i)) <= zero))
      // testing for the zero under the diagonal in the current row in case of nonfirst step
        continue; // we may continue on the next step
      
      double r;
      // solving row ratio coef.
      if (i == 0)
      // in case of the first step we must take values from a, because we suppose b is empty
        r = a(i,i)/a(j,i);
      else
        r = b(i,i)/b(j,i);

      b(j, i) = 0.0; // element under diagonal is right zero
      for (k = i+1; k < a.n; k++)
      {
        if (i == 0)
        // eliminating of the row in case of the first step
          b(j,k) = a(i,k) - a(j,k) * r;
        else
        // eliminating of the row in case of the nonfirst step
          b(j,k) = b(i,k) - b(j,k) * r;
      }
    }
  }
  return (0);
}



/**
   The function solves the system of equtions by the Gauss elimination. A.x = r
   The resulting eliminated %matrix is stored b, solution %vector x is stored in 
   the %vector r.
   
   @param[in]  a - the system equation %matrix,
   @param[out] b - the resulting factorised %matrix,
   @param[in,out]  r - the righthand side %vector/ resulting %vector,
   @param[in] zero - computer zero.
   
   @b Requests :
   a has to have the same number of the rows and columns
   b has to have same dimensions as the %matrix a
   b has to have allocated memory array for the elements
     which is enough large to hold contents of the result.
   r has to have same dimension as number of rows in the %matrix a
   
   @retval 0 : on succes
   @retval 1 : in case incompatibility sizes of a and b matrices
   @retval 2 : in case dimension of the %matrix a  <  2
   @retval 3 : in case some diagonal elements is to small or zero
               the value which is assumed as zero is driven by macro
               USREPS, which is defined in the matrix.h
               
   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long gause(const matrix &a, matrix &b, vector &r,double zero)
{
  long i,j,k;

  if ((a.m != b.m) || (a.n != b.n) || (r.n != a.m) || (a.m != a.n))
  {
    print_err("Gauss elimination - incompatible dimension of matrices or vector"
              "a(%ld,%ld) X b(%ld,%ld) X r(%ld)", __FILE__, __LINE__, __func__, a.m, a.n, b.m, b.n, r.n);
    return (1);
  }
  if (a.m < 2)
  {
    print_err("matrix has less than 2 rows", __FILE__, __LINE__, __func__);
    return (2);
  }

  for (i = 0; i < a.n; i++)
    b(0,i) = a(0,i);  // copying of the first row of the a to the b

  for (i = 0; i < a.m-1; i++)
  {
    for (j = i+1; j < a.m; j++)
    {
      if (fabs(b(i,i)) <= zero)
      // testing for the zero on the diagonal
      {
        print_err("diagonal element %ld is too small or zero (%le)", __FILE__, __LINE__, __func__, i+1, b(i,i));
        return (3);
      }
      if ((i == 0) && (fabs(a(j,i)) <= zero))
      // testing for the zero under the diagonal in the current row in case of the first step
      {
        // zero is hit so we have to copy this row from the a to the b to
        // fill empty b's rows
        for (k = 0; k < b.n; k++)
          b(j,k) = a(j,k);
        continue; // we may continue on the next step
      }
      if ((i != 0) && (fabs(b(j,i)) <= zero))
      // testing for the zero under the diagonal in the current row in case of nonfirst step
        continue; // we may continue on the next step
      
      double q;
      // solving row ratio coef.
      if (i == 0)
      // in case of the first step we must take values from a, because we suppose b is empty
        q = a(i,i)/a(j,i);
      else
        q = b(i,i)/b(j,i);
      
      b(j, i) = 0.0; // element under diagonal is right zero
      r[j] = r[i] - r[j] * q; // eliminating righthand side
      for (k = i+1; k < a.n; k++)
      {
        if (i == 0)
        // eliminating of the row in case of the first step
          b(j,k) = a(i,k) - a(j,k) * q;
        else
        // eliminating of the row in case of the nonfirst step
          b(j,k) = b(i,k) - b(j,k) * q;
      }
    }
  }

  if (fabs(b(b.m-1,b.n-1)) < zero)
  {
    print_err("diagonal element %ld is too small or zero (%le)", __FILE__, __LINE__, __func__, b.m, b(b.m-1,b.n-1));
    return (3);
  }
  // computing solution, the solution is stored in the vector r
  r[r.n-1] = r[r.n-1] / b(b.m-1, b.n-1);
  for (i = b.m-2; i >= 0; i--)
  {
    for (j = b.n-1; j > i; j--)
      r[i] -= b(i,j) * r[j];
    r[i] /= b(i,i);
  }
  return (0);
}



/**
   The function solves system of linear algebraic equtions by the Gauss elimination method.
   The system of equation is given by %matrix a, the resulting eliminated %matrix is 
   stored in b, solution is stored in the %vector sol.
   
   @param[in] a - system %matrix,
   @param[out] b - resulting eliminated %matrix,
   @param[in] r - righthand side %vector,
   @param[out] sol - solution %vector,
   @param[in] zero - computer zero.
   
   @b Requests :
   a has to have the same number of the rows and columns
   b has to have same dimensions as the matrix A
   b has to have allocated memory array for the elements
     which is enough large to hold contents of the result.
   r has to have same dimension as number of rows in the %matrix a
   sol has to have same dimension as the %vector r
   sol has to have allocated memory array for the elements
       which is enough large to hold contents of the result.
       
   @retval 0 : on succes
   @retval 1 : in case incompatibility sizes of a and b matrices
   @retval 2 : in case dimension of the matrix a  <  2
   @retval 3 : in case some diagonal elements is to small or zero
               the value which is assumed as zero is driven by macro
               USREPS, which is defined in the matrix.h

   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long gause(const matrix &a, matrix &b, vector &r, vector &sol,double zero)
{
  long i,j,k;

  if ((a.m != b.m) || (a.n != b.n) || (r.n != a.m) ||
      (sol.n != a.m) || (a.m != a.n))
  {
    print_err("Gauss elimination - incompatible dimensions of matrices or vectors"
              "a(%ld,%ld) X b(%ld,%ld) X r(%ld) X sol(%ld)", __FILE__, __LINE__, __func__, a.m, a.n, b.m, b.n, r.n, sol.n);
    return (1);
  }
  if (a.m < 2)
  {
    print_err("matrix a has less than 2 rows", __FILE__, __LINE__, __func__);
    return (2);
  }

  for (i = 0; i < a.n; i++)
    b(0,i) = a(0,i);  // copying of the first row of the a to the b

  for (i = 0; i < a.m-1; i++)
  {
    for (j = i+1; j < a.m; j++)
    {
      if (fabs(b(i,i)) <= zero)
      // testing for the zero on the diagonal
      {
        print_err("diagonal element %ld is too small or zero (%le)", __FILE__, __LINE__, __func__, i+1, b(i,i));
        return (4);
      }
      if ((i == 0) && (fabs(a(j,i)) <= zero))
      // testing for the zero under the diagonal in the current row in case of the first step
      {
        // zero is hit so we have to copy this row from the a to the b to
        // fill empty b's rows
        for (k = 0; k < b.n; k++)
          b(j,k) = a(j,k);
        continue; // we may continue on the next step
      }
      if ((i != 0) && (fabs(b(j,i)) <= zero))
      // testing for the zero under the diagonal in the current row in case of nonfirst step
        continue; // we may continue on the next step

      double q;
      // solving row ratio coef.
      if (i == 0)
      // in case of the first step we must take values from a, because we suppose b is empty
        q = a(i,i)/a(j,i);
      else
        q = b(i,i)/b(j,i);

      b(j, i) = 0.0; // element under diagonal is right zero
      r[j] = r[i] - r[j] * q; // eliminating righthand side
      for (k = i+1; k < a.n; k++)
      {
        if (i == 0)
        // eliminating of the row in case of the first step
          b(j,k) = a(i,k) - a(j,k) * q;
        else
        // eliminating of the row in case of the nonfirst step
          b(j,k) = b(i,k) - b(j,k) * q;
      }
    }
  }

  if (fabs(b(b.m-1,b.n-1)) < zero)
  {
    print_err("diagonal element %ld is too small or zero (%le)", __FILE__, __LINE__, __func__, b.m, b(b.m-1,b.n-1));
    return (3);
  }
  // computing solution, the solution is stored in the vector r
  sol[sol.n-1] = r[r.n-1] / b(b.m-1, b.n-1);
  for (i = b.m-2; i >= 0; i--)
  {
    sol[i] = r[i];
    for (j = b.n-1; j > i; j--)
      sol[i] -= b(i,j) * sol[j];
    sol[i] /= b(i,i);
  }
  return (0);
}



/**
  The function solves equation system by Gauss elimination method.
  The %matrix is stored as dense %matrix, x and y are also dense matrices, 
  stored by rows and right hand side vectors and vectors of unknowns are columns!

  @param[in,out] a - %matrix of the system of equations,
  @param[out] x - %vector of unknowns,
  @param[in,out] y - %vector of right hand side,
  @param[in] zero - computer zero,
  @param[in] pivot - type of pivoting.
   
  pivot=1 - pivot is searched only when there is zero on the diagonal
  pivot=2 - pivot is searched every time
   
  @author JK, 24. 8. 2016
*/
long gemp (matrix &a,vector &x,vector &y,double zero,long pivot)
{
  long i,j,k,ii,jj,ir,ic,stop,n;
  long *order;
  double f,s;
  
  n=a.m;
  order = new long [n];
  stop=0;
  
  //  initialization of order array
  for (i=0;i<n;i++){
    order[i]=i;
  }
  
  for (k=0;k<n-1;k++){
    ir=k;  ic=k;
    
    if (pivot==1){
      //  pivot is searched only when A(k,k)=0
      if (fabs(a[k][k])<zero){
        s=0.0;
        for (i=k;i<n;i++){
          for (j=k;j<n;j++){
            f=fabs(a[i][j]);
            if (f>s){
              s=f;  ir=i;  ic=j;
            }
          }
        }
        if (s<zero){
          //print_err("singular matrix is detected", __FILE__, __LINE__, __func__);
          stop=1;
          break;
        }
      }
    }
    
    if (pivot==2){
      //  pivot is searched in every step
      s=0.0;
      for (i=k;i<n;i++){
        for (j=k;j<n;j++){
          f=fabs(a[k][k]);
          if (f>s){
            s=f;  ir=i;  ic=j;
          }
        }
      }
      if (s<zero){
        //print_err("singular matrix is detected", __FILE__, __LINE__, __func__);
        stop=1;
        break;
      }
    }
    
    //  rows exchange
    if (ir!=k){
      for (i=0;i<n;i++){
        s=a[k][i];
        a[k][i]=a[ir][i];
        a[ir][i]=s;
      }

      s=y[k];
      y[k]=y[ir];
      y[ir]=s;
    }

    //  column exchange
    if (ic!=k){
      i=order[k];
      order[k]=order[ic];
      order[ic]=i;
      
      for (i=0;i<n;i++){
        s=a[i][k];
        a[i][k]=a[i][ic];
        a[i][ic]=s;
      }                         
    }
    
    
    //  elimination
    
    for (i=k+1;i<n;i++){
      s=a[i][k]/a[k][k];

      //  modification of matrix A
      for (j=k;j<n;j++){
        a[i][j]-=s*a[k][j];
      }

      //  modification of right hand side
      y[i]-=s*y[k];
    }
    
    
  }
  
  if (fabs(a[n-1][n-1])<zero){
    //print_err("singular matrix is detected", __FILE__, __LINE__, __func__);
    stop=1;
  }
  

  //  back substitution
  if (stop==0){
    for (k=n-1;k>-1;k--){
      s=0.0;
      for (j=n-1;j>k;j--){
        s+=a[k][j]*x[j];
      }
      x[k]=(y[k]-s)/a[k][k];
    }
    
    //  reordering of unknowns into original positions
    for (i=0;i<n;i++){
      if (order[i]!=i){
        for (j=i;j<n;j++){
          if (order[j]==i){
            jj=j;  break;
          }
        }
        
        ii=order[i];
        order[i]=order[jj];
        order[jj]=ii;
        
        s=x[i];
        x[i]=x[jj];
        x[jj]=s;
      }
    }
  }
  
  delete [] order;
  return stop;
}




/**
   The function inverts %matrix A by the Gauss elimination, i.e. B := A^{-1}.
   The result is stored in the matrix b, contents of the a is not changed.
   
   @param[in] a - %matrix to be inverted
   @param[out] b - resulting inverse %matrix
   
   @b Requests :
   a has to have the same number of the rows and columns
   b has to have same dimensions as the %matrix a
   b has to have allocated memory array for the elements
     which is enough large to hold contents of the result.
     
   @retval 0 : on succes
   @retval 1 : in case incompatibility sizes of a and b matrices
   @retval 2 : in case dimension of the %matrix a  <  2
   
   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long invm(const matrix &a, matrix &b,double zero)
{
  long i,j,k;
  double r;
  //long *order1, *order2;
  
  if ((a.m != b.m) || (a.n != b.n) || (a.m != a.n))
    {
      print_err("cannot invert matrix due to incompatible dimensions of inverted and resulting matrices"
                "a(%ld,%ld) X b(%ld,%ld)", __FILE__, __LINE__, __func__, a.m, a.n, b.m, b.n);
      return (1);
    }
  if (a.m < 2)
    {
      print_err("matrix has less than 2 rows", __FILE__, __LINE__, __func__);
      return (2);
    }
  
  //order1 = new long [a.m];
  //order2 = new long [a.m];
  
  
  //  initialization of order array
  //for (i=0;i<a.m;i++){
  //order1[i]=i;
  //order2[i]=i;
  //}
  
  matrix c(a.m, a.n); // new temporary variable
  // copying contents of a to the c, changing contents of the b to elementary matrix
  for (i = 0; i < b.m; i++)
    for (j = 0; j < b.n; j++)
      {
	c(i,j) = a(i,j);
	i == j ? b(i,j) = 1.0 : b(i,j) = 0.0;
      }
  
  //downward elimination
  for (i = 0; i < c.m-1; i++)
    {
      /* //      if (pivot == 3)
      //      {
      p = fabs(c(i,i));  ir = i;
      for(k=i+1; k< c.m; k++)
      {
      if (fabs(c(k,i)) > p)
      {
      p = fabs(c(k,i));
      ir = k;
      }
      }
      if (ir != i)//  row exchange
      {
      mswapr(c, i, ir);
      mswapr(b, i, ir);
      ii =  order1[i];
      order1[i]=order1[ir];
      order1[ii]=i;
      }
      //      }
      */
      
      for (j = i+1; j < c.m; j++)
	{
	  if (fabs(c(i,i)) < zero)
	    // testing for the zero on the diagonal
	    {
              print_err("diagonal element %ld is too small or zero (%le)", __FILE__, __LINE__, __func__, i+1, c(i,i));
	      return (3);
	    }
	  if (fabs(c(j,i)) <= zero)
	    // testing for the zero under the diagonal in the current row
	    continue; // we may continue on the next step
	  
	  
	  // solving row ratio coef.
	  r = c(j,i)/c(i,i);
	  //      r = c(i,i)/c(j,i);
	  c(j,i) = 0.0;  // element under diagonal is right zero
	  for (k = 0; k < c.n; k++)
	    {
	      if (k > i)
		// in case c matrix we may elimante from position under diagonal of current row
		c(j,k) = c(j,k) - c(i,k) * r;
	      b(j,k) = b(j,k) - b(i,k) * r; // but on the elementary matrix we must do same operation on the whole row
	      //          c(j,k) = c(i,k) - c(j,k) * r;
	      //        b(j,k) = b(i,k) - b(j,k) * r; // but on the elementary matrix we must do same operation on the whole row
	    }
	}
    }
  
  // upward elimination
  for (i = c.m-1; i > 0; i--)
    {
      /* //      if (pivot == 3)
      //      {
      p = fabs(c(i,i));  ir = i;
      for(k=i-1; k>=0; k--)
      {
      if (fabs(c(k,i)) > p)
      {
      p = fabs(c(k,i));
      ir = k;
      }
      }
      if (ir != i)//  row exchange
      {
      mswapr(c, i, ir);
      mswapr(b, i, ir);
      ii =  order2[i];
      order2[i]=order2[ir];
      order2[ii]=i;
      }
      
      //      }
      */
      
      for (j = i-1; j >= 0; j--)
	{
	  if (fabs(c(i,i)) < zero)
	    // testing for the zero on the diagonal
	    {
              print_err("diagonal element %ld is too small or zero (%le)", __FILE__, __LINE__, __func__, i+1, c(i,i));
	      return (3);
	    }
	  if (fabs(c(j,i)) <= zero)
	    // testing for the zero under the diagonal in the current row
	    continue; // we may continue on the next step
	  
	  // solving row ratio coef.
	  r = c(j,i)/c(i,i);
	  //      r = c(i,i)/c(j,i);
	  c(j,i) = 0.0;  // element under diagonal is right zero
	  for (k = c.n-1; k >= 0; k--)
	    {
	      if (k < i)
		// in case c matrix we may elimante from position under diagonal of current row
		c(j,k) = c(j,k) - c(i,k) * r;
	      b(j,k) = b(j,k) - b(i,k) * r; // but on the elementary matrix we must do same operation on the whole row
	      //          c(j,k) = c(i,k) - c(j,k) * r;
	      //        b(j,k) = b(i,k) - b(j,k) * r; // but on the elementary matrix we must do same operation on the whole row
	    }
	}
    }
  
  // finally we must divide rows in the elementary matrix b by the numbers on the a diagonal
  for (i = 0; i < c.m; i++)
    for (j = 0; j < b.n; j++)
      {
	if (fabs(c(i,i)) < zero)
	  // testing for the zero on the diagonal
	  {
            print_err("diagonal element %ld is too small or zero (%le)", __FILE__, __LINE__, __func__, i+1, c(i,i));
	    return (3);
	  }
	b(i,j) /= c(i,i);
      }
  
  /* //  reordering into original positions
     for (i = 0; i < a.m; i++){
     if (order2[i]!=i){
     ir=order2[i];
     mswapr(b, ir, i);
     }
     }
     
     for (i = 0; i < a.m; i++){
     if (order1[i]!=i){
     ir=order1[i];
     mswapr(b, ir, i);
     }
     }  
     
     
     delete [] order1;
     delete [] order2;
  */
  
  return (0);
}



/**
   The function solves determinant of the %matrix a by the Gauss elimination,
   the result is stored in the det, contents of the a is not changed
   
   @param[in] a    - the given %matrix, whose determinant will being computed,
   @param[out] det - the resuting determinent value.
   
   @b Requests :
   a has to have the same number of the rows and columns
   
   @retval 0 : on succes
   @retval 1 : in case incompatibility sizes of the A matrix
   
   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long detm(const matrix &a, double &det)
{
  long i;
  double zero=1.0e-20;

  if (a.m != a.n)
  {
    print_err("wrong dimensions of matrix a(%ld,%ld)", __FILE__, __LINE__, __func__, a.m, a.n);
    return (1);
  }
  
  matrix b(a.m, a.n); // new temporary matrix
  if (gause(a, b, zero) == 3)
  // eliminating of the A matrix, result is stored in the B
  // in case zero on the diagonal, determinatnt is zero too
  {
    det = 0.0;
    return (0);
  }
  // determinant is solved as the product of the diagonal elements
  det = b(0,0);
  for (i = 1; i < b.m; i++)
    det *= b(i,i);
  return (0);
}



/**
  The function computes determinant of a %matrix 2x2.
   
  @param[in] a - the given %matrix,
  @param[out] det - the resulting value of the determinant of the %matrix.
   
  @author JK, 24. 8. 2016
*/
void det2x2 (const matrix &a,double &det)
{
  det = a[0][0]*a[1][1]-a[0][1]*a[1][0];
}



/**
  The function computes determinant of a %matrix 3x3.
   
  @param[in] a - the given %matrix,
  @param[out] det - the resulting value of determinant of the %matrix.
   
  @author JK, 24. 8. 2016
*/
void det3x3 (const matrix &a,double &det)
{
  det  = a[0][0]*a[1][1]*a[2][2] + a[0][1]*a[1][2]*a[2][0] + a[0][2]*a[1][0]*a[2][1];
  det += 0.0 - a[0][2]*a[1][1]*a[2][0] - a[1][2]*a[2][1]*a[0][0] - a[2][2]*a[0][1]*a[1][0];
}



/**
   The function prints out the contents of the %matrix a to the stdout file,
   in precision 3 digits, width of the number is 11
   optionally, can be specified file to print, precision and the field width.
   
   @param[in] a     is the structure of the %matrix which is printed.
   
   @b Optionally :
   @param out   is the structure with opened file for %matrix output (default is stdout)
   @param prec  is desired precision of the %matrix elements  (default is 3)
   @param width is desired width of the %matrix elements on the output (default is 11)
   
   @return always zero.
   
   @author Created  29.5.2000, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   @author modified  9.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long printm(const matrix &a, FILE *out, int prec, int width)
{
  long i,j;

  for (i = 0; i < a.m; i++)
  {
    for (j = 0; j < a.n; j++)
      fprintf(out, "%*.*e ", width, prec, a(i,j));
    fprintf(out, "\n");
  }
  return (0);
}



/**
   The function prints out the contents of the %matrix a to the file given by out.
   
   @param[in] out   is the structure with opened file for %matrix output,
   @param[in] a     is the structure of the %matrix which is printed.
   
   @return always zero.
   
   @author Created  10.12.2014, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long printm(FILE *out, const matrix &a)
{
  long i,j;

  for (i = 0; i < a.m; i++)
  {
    for (j = 0; j < a.n; j++)
      fprintf(out, "% e ", a(i,j));
    fprintf(out, "\n");
  }
  return (0);
}



/**
  The function swaps matrix's columns i and j.

  @param[in,out] mat is the structure of the matrix,
  @param[in] i   number of the first column which will be swapped,
  @param[in] j   number of the second column which will be swapped.

  @return always zero

  @author Created  22.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long mswapc(matrix &mat, long i, long j)
{
  long k;
  double s;
  for (k = 0; k < mat.m; k++)
  {
    s = mat[k][i];
    mat[k][i] = mat[k][j];
    mat[k][j] = s;
  }
  return (0);
}



/**
 The function swaps matrix's rows i and j.

  @param[in,out] mat - the given %matrix,
  @param[in] i - number of the first row which will be swapped,
  @param[in] j - number of the second row which will be swapped.

  @return always zero

  @author Created  22.5.2001, Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long mswapr(matrix &mat, long i, long j)
{
  long k;
  double s;
  for (k = 0; k < mat.n; k++)
  {
    s = mat[i][k];
    mat[i][k] = mat[j][k];
    mat[j][k] = s;
  }
  return (0);
}



/**
   The function performs a condensation of selected dofs in the stiffness %matrix sm.
   
   @param[in,out] sm - stiffness %matrix,
   @param[in] cu - %vector with number one on the place of condensed dofs.
   
   @author Created  23.5.2012, Jan Zitny, zitny.ja@gmail.com
*/
long condense_matrix(matrix &sm, ivector &cu)
{
  //allocates an array of specific places of condensed variables
  long n = sm.m;
  long x,i,j,z,v=0;
  ivector q;
  matrix  auxm(n,n);

  for(i=0;i<n;i++)
  {
    if(cu[i])
      v++;
  }
  reallocv(v, q);
  
//  for(i=0,j=0;i<n,j<v;i++)
  for(i=0,j=0;i<n && j<v;i++)
  {
    if(cu[i])
    {
      q[j]=i;
      j++;
    }
  }
  
  
  
  for(z=0;z<v;z++)
  {
    x = q[z];
    
    if(sm[x][x]==0.0)
      return 1;      //returns 1 if entered condensation was nonsense
    
    for(i=0;i<n;i++)
    {
      for(j=0;j<n;j++)
        auxm[i][j] = (sm[i][j]-((sm[i][x]*sm[x][j])/sm[x][x]));
    }
    copym(auxm, sm);
  }

  return 0;  //returns 1 if succeed
}



/**
   The function performs a condensation of selected dofs of the load %vector 
   according to corresponding stiffness %matrix.

   @param[in] sm - stiffness %matrix,
   @param[in,out] f - load %vector,
   @param[in] cu - is a %vector with number one on the position of condensed dofs.
   
   @author Created  23.5.2012, Jan Zitny, zitny.ja@gmail.com
*/
long condense_vector(matrix &sm, vector &nf, ivector &cu)
{
  // allocates an array of specific positions of condensed variables
  long n = sm.m;
  long x,i,j,z,v=0;
  ivector q;
  
  for(i=0;i<n;i++)
  {
    if(cu[i])
      v++;
  }
  
  reallocv(v, q);
  
//  for(i=0,j=0;i<n,j<v;i++)
  for(i=0,j=0;i<n && j<v;i++)
  {
    if(cu[i])
    {
      q[j]=i;
      j++;
    }
  }
  
  for(z=0;z<v;z++)
  {
    x = q[z];
    
    if(sm[x][x]==0.0)
    {
      return 1;      //  returns 1 if entered condensation was nonsense
    }
    
    for(i=0;i<n;i++)
      nf[i] -= ((nf[x]*sm[i][x])/sm[x][x]);
  }
  
  return 0;  // returns 1 if succeed
}



/**
  The function computes %matrix multiplication in the form B^T . D . B . jac
  and the result is added to the %matrix k. It is supposed that the resulting 
  %matrix k is symmetric and therefore only upper triangle components are 
  evaluated while the remaining are assigned only.

  @param[out] k - array of resulting %matrix K(n,n) = B^T . D . B . jac,
  @param[in] b - array of %matrix B(m,n),
  @param[in] d - array of %matrix D(m,m),
  @param[in] jac - jacobian (constant),
  @param[in] m - number of rows of %matrix B, number of rows/columns of %matrix D,
  @param[in] n - number of columns of %matrix B.

  @author Created by JK 12.7.1996
*/
void bdbj (double *k, const double *b, const double *d,double jac,long m,long n)
{
  long i,j,l,ii,jj,mm,ll,ul;
  double s;

  /*  auxiliary matrix  (am=b.d)  stored due to better memory management in vector structure */
  vector am(m*n);
  
  /*  product of matrix multiplication b^T.d=am  */
  ll=0;
  for (i=0;i<n;i++){
    for (j=0;j<m;j++){
      s=0.0;  ii=i;  jj=j;
      for (l=0;l<m;l++){
        s+=b[ii]*d[jj];
        ii+=n;  jj+=m;
      }
      am[ll]=s;  ll++;
    }
  }
  
  /*  product of matrix multiplication am.b=k  */
  for (i=0;i<n;i++){
    ii=i*n+i;  mm=ii;
    s=0.0;  ll=i*m;  ul=ll+m;  jj=i;
    for (l=ll;l<ul;l++){
      s+=am[l]*b[jj];
      jj+=n;
    }
    k[ii]+=s*jac;
    ii++;  mm+=n;
    
    for (j=i+1;j<n;j++){
      s=0.0;  jj=j;
      for (l=ll;l<ul;l++){
        s+=am[l]*b[jj];
        jj+=n;
      }
      k[ii]+=s*jac;  k[mm]+=s*jac;
      ii++;  mm+=n;
    }
  }
}



/**
  The function computes %matrix multiplication in the form A^T . B . C . jac
  and the result is added to the %matrix d. 

  @param[out] d - resulting %matrix D(n,n) += A^T . B . C . jac,
  @param[in] a - %matrix A(m,n),
  @param[in] b - %matrix B(m,m),
  @param[in] b - %matrix C(m,n),
  @param[in] jac - jacobian (constant).

  @author Created by JK 12.7.1996
   
  @date 10.5.2002
*/
void bdbjac (matrix &d, const matrix &a, const matrix &b, const matrix &c, double jac)
{
  matrix x(a.n, b.n), y(d.m, d.n);
  mtxm (a,b,x);
  mxm (x,c,y);
  cmulm (jac,y);
  addm (y,d,d);
}



/**  
  The function evaluates %matrix product N^T . N . jac,
  the result is added to the %matrix m.

  @param[out] m - array of resulting %matrix M += N^T . N. jac,
  @param[in] n - array of %matrix N,
  @param[in] jac - jacobian (constant),
  @param[in] mm - number of rows of %matrix n,
  @param[in] nn - number of columns of %matrix n.

  @author Created by JK 9.7.2001
*/
void nnj (double *m, const double *n, double jac, long mm, long nn)
{
  long i,j,k,ii,jj,kk,ll;
  double s;
  
  for (i=0;i<nn;i++){
    kk=i*nn+i;  ll=kk;
    s=0.0;  ii=i;  jj=i;
    for (k=0;k<mm;k++){
      s+=n[ii]*n[jj];
      ii+=nn;  jj+=nn;
    }
    m[kk]+=s*jac;
    kk++;  ll+=nn;
    
    for (j=i+1;j<nn;j++){
      s=0.0;  ii=i;  jj=j;
      for (k=0;k<mm;k++){
        s+=n[ii]*n[jj];
        ii+=nn;  jj+=nn;
      }
      m[kk]+=s*jac;  m[ll]+=s*jac;
      kk++;  ll+=nn;
    }
  }
}



/**
  The function evaluates %matrix product A^T . B . jac,
  the result is added to the %matrix d.

  @param[out] d - resulting %matrix D += A^T . B . jac,
  @param[in] a - %matrix A,
  @param[in] b - %matrix B,
  @param[in] jac - jacobian (constant).

   
  @author Created by JK 10.5.2002
*/
void nnjac (matrix &d,const matrix &a,const matrix &b,double jac)
{
  matrix x(a.n,b.n);

  mtxm (a,b,x);
  cmulm (jac,x);
  addm (x,d,d);
}



/**
  The function transforms %vector of element from local element 
  coordinate system to global problem coordinate system
  g = T l

  The transformation %matrix has to be in the form x_g = T x_l
   
  @param[out] g - %vector in the global problem coordinate system
  @param[in] l - %vector in the local element coordinate system
  @param[in] tmat - transformation %matrix T
   
  @author Created by JK, 28.11.2006
*/
void lgvectortransf (vector &g,const vector &l,const matrix &tmat)
{
  mxv (tmat,l,g);
}



/**
  The function transforms %vector of element from global problem coordinate system
  to local element coordinate system
  l = T^T g

  The transformation %matrix T has to be in the form x_g = T x_l
   
  @param[in] g - %vector in the global problem coordinate system
  @param[out] l - %vector in the local element coordinate system
  @param[in] tmat - transformation %matrix T
   
  @author Created by JK, 28.11.2006
*/
void glvectortransf (const vector &g,vector &l,const matrix &tmat)
{
  mtxv (tmat,g,l);
}



/**
  The function transforms %matrix a of element from global problem coordinate system
  to local element coordinate system  A_l = T^T A_g T.

  The transformation %matrix T has to be in the form x_g = T x_l,
  the transormed %matrix is stored in the %matrix a.
   
  @param[in,out] a - %matrix of element A_g (input)/resulting transformed %matrix A_l(output)
  @param[in] tmat - transformation %matrix T
   
  @author Created by JK, 28.11.2006
*/
void glmatrixtransf (matrix &a,const matrix &tmat)
{
  long n;
  n=a.m;
  matrix am(n,n);
  
  mtxm (tmat,a,am);
  mxm (am,tmat,a);
}



/**
  The function transforms %matrix a of element from local element coordinate system
  to global problem coordinate system A_g = T A_l T^T.

  The transformation %matrix T has to be in the form x_g = T x_l,
  the transformed %matrix is stored in the %matrix a.
   
  @param[in,out] a - %matrix of element A_l(at input)/resulting transformed %matrix A_g (at output),
  @param[in] tmat - transformation %matrix T,
   
  @author Created by JK, 28.11.2006
*/
void lgmatrixtransf (matrix &a,const matrix &tmat)
{
  long n;
  n=a.m;
  matrix am(n,n);
  
  mxm (tmat,a,am);
  mxmt (am,tmat,a);
}



/**
  The function transforms %matrix expressed in local coordinate
  system to %matrix expressed in global coordinate system  A_g = T . A_l . T^T

  The dimension (dim) of the problem is determined from the size of
  the transformation %matrix T. The %matrix A may contain several blocks 
  with size (dim,dim)

  @param a - local %matrix (at input), global %matrix (at output),
  @param t - transformation %matrix T.
   
  @author Created by JK,
*/
void lgmatrixtransfblock (matrix &a, const matrix &t)
{
  long i,j,k,ii,jj;
  double s;
  //  problem dimension, 2D or 3D
  long dim=t.n;
  matrix b(dim, dim);

  for (ii=0;ii<a.n;ii+=dim){
    for (jj=0;jj<a.n;jj+=dim){
      for (i=0;i<dim;i++){
	for (j=0;j<dim;j++){
	  s=0.0;
	  for (k=0;k<dim;k++){
	    s+=a[ii+i][jj+k]*t[j][k];
	  }
	  b[i][j]=s;
	}
      }
      for (i=0;i<dim;i++){
	for (j=0;j<dim;j++){
	  s=0.0;
	  for (k=0;k<dim;k++){
	    s+=t[i][k]*b[k][j];
	  }
	  a[ii+i][jj+j]=s;
	}
      }
    }
  }
}



/**
  The function transforms %matrix expressed in global coordinate
  system to %matrix expressed in local coordinate system  A_l = T^T . A_g . T.
  Thus the transformation %matrix T is defined as v_g = T * v_l, where v_l is the 
  %vector in the given local coordinate system and v_g is the %vector in the global 
  coordinate system.

  The dimension (dim) of the problem is determined from the size of
  the transformation %matrix T. The %matrix A may contain several blocks 
  with size (dim,dim)

  @param[in,out] a - local %matrix (at input), global %matrix (at output),
  @param t - transformation %matrix T.
   
  @author Created by TKo according to JK, 09.2023
*/
void glmatrixtransfblock (matrix &a, const matrix &t)
{
  long i,j,k,ii,jj;
  double s;
  //  problem dimension, 2D or 3D
  long dim=t.n;
  matrix b(dim, dim);

  for (ii=0;ii<a.n;ii+=dim){
    for (jj=0;jj<a.n;jj+=dim){
      for (i=0;i<dim;i++){
	for (j=0;j<dim;j++){
	  s=0.0;
	  for (k=0;k<dim;k++){
	    s+=a[ii+i][jj+k]*t[k][j];
	  }
	  b[i][j]=s;
	}
      }
      for (i=0;i<dim;i++){
	for (j=0;j<dim;j++){
	  s=0.0;
	  for (k=0;k<dim;k++){
	    s+=t[k][i]*b[k][j];
	  }
	  a[ii+i][jj+j]=s;
	}
      }
    }
  }
}



/**
  The function transforms %vector v expressed in global coordinate
  system to %vector expressed in local coordinate system  v_l = T^T . v_g

  The dimension (dim) of the problem is determined from the size of
  the transformation %matrix T and %vector v may contain several 
  blocks with size (dim).

  @param[in,out] v - global %vector v_g (at input), local %vector v_l (at output),
  @param[in] t - transformation %matrix T.
   
  @author Created by JK,
*/
void glvectortransfblock (vector &v, const matrix &t)
{
  /*
  long i,i1,i2;
  double c1,c2;
  for (i=0;i<v.n;i=i+3){
    i1=i+1;
    i2=i+2;
    c1   =t[0][0]*v[i]+t[1][0]*v[i1]+t[2][0]*v[i2];
    c2   =t[0][1]*v[i]+t[1][1]*v[i1]+t[2][1]*v[i2];
    v[i2]=t[0][2]*v[i]+t[1][2]*v[i1]+t[2][2]*v[i2];
    v[i] =c1;
    v[i1]=c2;
  }
  */
  
  long i,j,ii;
  double s;
  //  problem dimension, 2D or 3D
  long dim=t.n;
  vector b(dim);
  
  
  for (ii=0;ii<v.n;ii=ii+dim){
    for (i=0;i<dim;i++){
      s=0.0;
      for (j=0;j<dim;j++){
	s+=t[j][i]*v[ii+j];
      }
      b[i]=s;
    }
    for (i=0;i<dim;i++){
      v[ii+i]=b[i];
    }
  }
  
}



/**
  The function transforms %vector v expressed in local coordinate
  system to %vector expressed in global coordinate system  v_g = T * v_l

  The dimension (dim) of the problem is determined from the size of
  the transformation %matrix T and  %vector v may contain several blocks 
  with size (dim)

  @param[in,out] v - local %vector v_l (at input), global %vector v_g (at output),
  @param[in] t - transformation %matrix T.
   
  @author Created by JK,
*/
void lgvectortransfblock (vector &v, const matrix &t)
{
  /*
  long i,i1,i2;
  double c1,c2;
  for (i=0;i<v.n;i=i+3){
      i1=i+1;
      i2=i+2;
      c1   =t[0][0]*v[i]+t[0][1]*v[i1]+t[0][2]*v[i2];
      c2   =t[1][0]*v[i]+t[1][1]*v[i1]+t[1][2]*v[i2];
      v[i2]=t[2][0]*v[i]+t[2][1]*v[i1]+t[2][2]*v[i2];
      v[i] =c1;
      v[i1]=c2;
  }
  */
  
  long i,j,ii;
  double s;
  //  problem dimension, 2D or 3D
  long dim=t.n;
  vector b(dim);
  
  for (ii=0;ii<v.n;ii=ii+dim){
    for (i=0;i<dim;i++){
      s=0.0;
      for (j=0;j<dim;j++){
	s+=t[i][j]*v[ii+j];
      }
      b[i]=s;
    }
    for (i=0;i<dim;i++){
      v[ii+i]=b[i];
    }
  }
  
}



/**
  The function localizes components of local %vector lv to the global %vector gv.
   
  @param[in,out]  gv - array of global %vector (%vector of the whole problem/structure),
  @param[in]  lv - array of local %vector (%vector of one element),
  @param[in]  cn - array containing code numbers,
  @param[in]  n - number of components in the local %vector.
   
  @author Created by JK  25.7.2001
*/
void locglob (double *gv,const double *lv,const long *cn,long n)
{
  long i,j;
  
  for (i=0;i<n;i++){
    j=cn[i];
    if (j<1)  continue;
    else{
      j--;
      //#ifdef INC_OPENMP      
      //#pragma omp atomic update
      //#endif
      gv[j]+=lv[i];
    }
  }
}



/**
  The function localizes components of global %vector gv to the local %vector lv.
   
  @param[in,out] gv - array of global %vector (%vector of the whole problem/structure),
  @param[in] lv - array of local %vector (%vector of one element),
  @param[in] cn - array containing code numbers,
  @param[in] n - number of components in the local %vector.

  @author Created by JK
*/
void globloc (const double *gv, double *lv, const long *cn, long n)
{
  long i,j;
  
  for (i=0;i<n;i++){
    j=cn[i];
    if (j<1)  continue;
    else{
      j--;
      lv[i]=gv[j];
    }
  }
}



void locvecmat(double *mat, const double *vect, const long *cn, long ci, long m, long n)
{
  long i;
  
  for (i = 0; i < n; i++){
    if (cn[i] < 1)
      continue;
    mat[(cn[i]-1)*m+ci] += vect[i];
  }
}



/**
  The function localizes components of local %matrix lm to the global %matrix gm.
   
  @param[in,out] gv - global %matrix (%matrix of the whole problem/structure),
  @param[in] lv - local %matrix (%vmatrix of one element),
  @param[in] rcn - array containing code numbers for %matrix rows,
  @param[in] ccn - array containing code numbers for %matrix columns.
   
  @author Created by JK  25.7.2001
*/
void mat_localize (matrix &gm, matrix &lm, long *rcn, long *ccn)
{
  long i,j,m,n;
  
  m=lm.m;  n=lm.n;
  
  for (i=0;i<m;i++){
    if (rcn[i]<=0)  continue; //if (rcn[i]==0)  continue; //changed by TKr 26/08/2010
    for (j=0;j<n;j++){
      if (ccn[j]<=0)  continue; //if (ccn[j]==0)  continue; //changed by TKr 26/08/2010
      gm[rcn[i]-1][ccn[j]-1]+=lm[i][j];
    }
  }
  
}



/**
  The function computes determinant of the %matrix

  | 1  x1  y1 |
  | 1  x2  y2 |
  | 1  x3  y3 |

  which is equal to double area of the triangle given by 
  vectors of vertex coordinates.

  @param[in] x - array of x-coordinates of the triangle verteces,
  @param[in] y - array of y-coordinates of the triangle verteces.

  @author Created by JK, 24.8.2001
*/
double det2d (const double *x, const double *y)
{
  double d;
  d = x[0]*(y[1]-y[2])+ x[1]*(y[2]-y[0]) + x[2]*(y[0]-y[1]);
  return d;
}



/**
  The function computes determinant of the %matrix

  | 1 x1 y1 z1 |
  | 1 x2 y2 z2 |
  | 1 x3 y3 z3 |
  | 1 x4 y4 z4 |

  which is equal to six-multiple of the volume of the 
  tetrahedron whose vertices are given by coordinate 
  %vectors x,y,z.

  @param[in] x - %vector of x-coordinates of tetrahedron vertices,
  @param[in] y - %vector of y-coordinates of tetrahedron vertices,
  @param[in] z - %vector of z-coordinates of tetrahedron vertices.

  @author Created by JK, 24.8.2001, 
  @author Modified by JK 23.9.2011
*/
double det3d (const double *x,const double *y,const double *z)
{
  double d;
  /*
  d  = x[0]*(y[1]*(z[3]-z[2]) + y[2]*(z[1]-z[3]) + y[3]*(z[2]-z[1]));
  d += x[1]*(y[0]*(z[2]-z[3]) + y[2]*(z[3]-z[0]) + y[3]*(z[0]-z[2]));
  d += x[2]*(y[0]*(z[3]-z[1]) + y[1]*(z[0]-z[3]) + y[3]*(z[1]-z[0]));
  d += x[3]*(y[0]*(z[1]-z[2]) + y[1]*(z[2]-z[0]) + y[2]*(z[0]-z[1]));
*/
  /*
  d  = x[1]*(y[2]*z[3]-y[3]*z[2]) - x[2]*(y[1]*z[3]-y[3]*z[1]) + x[3]*(y[1]*z[2]-y[2]*z[1]);
  d -= x[0]*(y[2]*z[3]-y[3]*z[2]) - x[2]*(y[0]*z[3]-y[3]-z[0]) + x[3]*(y[0]*z[2]-y[2]*z[0]);
  d += x[0]*(y[1]*z[3]-y[3]*z[1]) - x[1]*(y[0]*z[3]-y[3]*z[0]) + x[3]*(y[0]*z[1]-y[1]*z[0]);
  d -= x[0]*(y[1]*z[2]-y[2]*z[1]) - x[1]*(y[0]*z[2]-y[2]*z[0]) + x[2]*(y[0]*z[1]-y[1]*z[0]);
  */
  double a[3],b[3],c[3];

  a[0]=x[0]-x[3];
  a[1]=y[0]-y[3];
  a[2]=z[0]-z[3];

  b[0]=x[1]-x[3];
  b[1]=y[1]-y[3];
  b[2]=z[1]-z[3];

  c[0]=x[2]-x[3];
  c[1]=y[2]-y[3];
  c[2]=z[2]-z[3];
  
  d = a[0]*b[1]*c[2] + b[0]*c[1]*a[2] + c[0]*a[1]*b[2] - c[0]*b[1]*a[2] - c[1]*b[2]*a[0] - c[2]*b[0]*a[1];

  return d;
}



/**
  Function computes principal values and vectors of
  tensors of 2D and 3D problems .

  @param[in] v - matrix of second order tensor,
  @param[out] pval - %vector of sorted principal values (from min to max),
  @param[out] pvect - %matrix of principal direction vectors stored in columns,
  @param[in] ni - number of iterations in jacobi method,
  @param[in] err - required error of Jacobi method,
  @param[in] zero - computer zero for explicit computation method,
  @param[in] n - order of problem (only 2 for 2D problem and 3 for 3D problem are supported),
  @param[in] normalize - switch for normalizing of matrix in Jacobi method.

  @return eigenvectors create columns of the matrix pvect and eigenvalues are stored and sorted in pval vector

  @author Created 29.8.2001
  @author Modified 17.8.2007 TKo (added normalizing in Jacobi method)
*/
void princ_val (const matrix &v,vector &pval,matrix &pvect,long ni,double err,double zero,long n, long normalize)
{
  long ani;
  
  if (n==2){
    double i1,i2,d,s;
    //  first invariant
    i1=v[0][0]+v[1][1];
    //  second invariant
    i2=v[0][0]*v[1][1]-v[0][1]*v[1][0];
    //  discriminant
    d=i1*i1-4.0*i2;
    if (d<0.0)  print_err("negative discriminant (%le)",__FILE__,__LINE__,__func__,d);
    pval[0]=(i1+sqrt(d))/2.0;
    pval[1]=(i1-sqrt(d))/2.0;
    
    s=v[0][0]-pval[0];
    if (fabs(s)<zero){
      pvect[0][0]=1.0;
      pvect[1][0]=0.0;
    }
    else{
      pvect[0][0]=-v[0][1]/s;
      pvect[1][0]=1.0;
    }
    
    s=v[0][0]-pval[1];
    if (fabs(s)<zero){
      pvect[0][1]=1.0;
      pvect[1][1]=0.0;
    }
    else{
      pvect[0][1]=-v[0][1]/s;
      pvect[1][1]=1.0;
    }
    
    //  normalization
    s=pvect[0][0]*pvect[0][0]+pvect[1][0]*pvect[1][0];
    s=sqrt(s);
    pvect[0][0]/=s;  pvect[1][0]/=s;
    
    s=pvect[0][1]*pvect[0][1]+pvect[1][1]*pvect[1][1];
    s=sqrt(s);
    pvect[0][1]/=s;  pvect[1][1]/=s;

    //  sorting
    if (pval[0]>pval[1]){
      s=pval[0];
      pval[0]=pval[1];
      pval[1]=s;
      
      s=pvect[0][0];  pvect[0][0]=pvect[0][1];  pvect[0][1]=s;
      s=pvect[1][0];  pvect[1][0]=pvect[1][1];  pvect[1][1]=s;
    }
  }

  if (n==3){
    double a[9],b[9];
    
    a[0]=v[0][0];  a[1]=v[0][1];  a[2]=v[0][2];
    a[3]=v[1][0];  a[4]=v[1][1];  a[5]=v[1][2];
    a[6]=v[2][0];  a[7]=v[2][1];  a[8]=v[2][2];

    jacobi_rot (a,b,pval.a,n,ni,ani,err, normalize);
    
    pvect[0][0]=b[0];  pvect[0][1]=b[1];  pvect[0][2]=b[2];
    pvect[1][0]=b[3];  pvect[1][1]=b[4];  pvect[1][2]=b[5];
    pvect[2][0]=b[6];  pvect[2][1]=b[7];  pvect[2][2]=b[8];
    
    double temp;
    if (pval[0] > pval[1])
      {
        temp = pval[0];
        pval[0] = pval[1];
        pval[1] = temp;
        mswapc(pvect,0,1);
      }
    if (pval[0] > pval[2])
      {
        temp = pval[0];
        pval[0] = pval[2];
        pval[2] = temp;
        mswapc(pvect,0,2);
      }
    if (pval[1] > pval[2])
      {
        temp = pval[1];
        pval[1] = pval[2];
        pval[2] = temp;
        mswapc(pvect,1,2);
      }
  }
}



/**
  The function performes Gaussian elimination on %matrix a.
  The %matrix is stored as dense %matrix, x and y are also dense matrices, stored by rows.
  The right hand side vectors and vectors of unknowns are columns!

  @param[in,out] a - array containing matrix of the system,
  @param[out] x - array containing vectors of unknowns,
  @param[in,out] y - array containing right hand sides,
  @param[in] n - order of matrix,
  @param[in] m - number of right hand sides,
  @param[in] limit - computer zero,
  @param[in] pivot=1 - pivot is searched only when there is zero on the diagonal,
             pivot=2 - pivot is searched every time.

  @author Created by JK, 23.7.2001
*/
void gemp (double *a,double *x,double *y,long n,long m,double limit,long pivot)
{
  long i,j,k,ii,jj,kk,ir,ic,li,ui,*order;
  double f,s;
  
  order = new long [n];
  
  //  initialization of order array
  for (i=0;i<n;i++){
    order[i]=i;
  }
  
  for (k=0;k<n-1;k++){
    ir=k;  ic=k;

    if (pivot==1){
      //  pivot is searched only when A(k,k)=0
      if (fabs(a[k*n+k])<limit){
        s=0.0;  ii=k*n+k;
        for (i=k;i<n;i++){
          for (j=k;j<n;j++){
            f=fabs(a[ii]);
            if (f>s){
              s=f;  ir=i;  ic=j;
            }
            ii++;
          }
          ii+=k;
        }
        if (s<limit){
          print_err("singular matrix (%le)", __FILE__, __LINE__, __func__,s);
          break;
        }
      }
    }
    
    if (pivot==2){
      //  pivot is searched in every step
      s=0.0;  ii=k*n+k;
      for (i=k;i<n;i++){
        for (j=k;j<n;j++){
          f=fabs(a[ii]);
          if (f>s){
            s=f;  ir=i;  ic=j;
          }
          ii++;
        }
        ii+=k;
      }
      if (s<limit){
        print_err("singular matrix (%le)", __FILE__, __LINE__, __func__,s);
        break;
      }
    }
    
    //  rows exchange
    if (ir!=k){
      li=k*n+k;  ui=(k+1)*n;  j=ir*n+k;
      for (i=li;i<ui;i++){
        s=a[i];
        a[i]=a[j];
        a[j]=s;
        j++;
      }

      li=k*m;  ui=li+m;  j=ir*m;
      for (i=li;i<ui;i++){
        s=y[i];
        y[i]=y[j];
        y[j]=s;
        j++;
      }
    }

    //  column exchange
    if (ic!=k){
      i=order[k];
      order[k]=order[ic];
      order[ic]=i;
      
      j=k;  ii=ic;
      for (i=0;i<n;i++){
        s=a[j];
        a[j]=a[ii];
        a[ii]=s;
        j+=n;  ii+=n;
      }                         
    }
    
    
    //  elimination
    
    for (i=k+1;i<n;i++){
      ii=i*n+k;  kk=k*n+k;
      s=a[ii]/a[kk];

      //  modification of matrix A
      for (j=k;j<n;j++){
        a[ii]-=s*a[kk];
        ii++;  kk++;
      }

      //  modification of right hand sides
      ii=i*m;  kk=k*m;
      for (j=0;j<m;j++){
        y[ii]-=s*y[kk];
        ii++;  kk++;
      }
    }
    
    
  }
  
  if (fabs(a[(n-1)*n+n-1])<limit){
    print_err("singular matrix (%le)", __FILE__, __LINE__, __func__,a[(n-1)*n+n-1]);
  }
  

  //  back substitution
  
  for (k=n-1;k>-1;k--){
    f=a[k*n+k];  kk=k*m;
    for (i=0;i<m;i++){
      s=0.0;  ii=k*n+n-1;  jj=(n-1)*m+i;
      for (j=n-1;j>k;j--){
        s+=a[ii]*x[jj];
        ii--;  jj-=m;
      }
      x[kk]=(y[kk]-s)/f;
      kk++;
    }
  }
  
  //  reordering of unknowns into original positions
  for (i=0;i<n;i++){
    if (order[i]!=i){
      for (j=i;j<n;j++){
        if (order[j]==i){
          jj=j;  break;
        }
      }
      
      ii=order[i];
      order[i]=order[jj];
      order[jj]=ii;
      
      ii=i*m;  kk=jj*m;
      for (j=0;j<m;j++){
        s=x[ii];
        x[ii]=x[kk];
        x[kk]=s;
        ii++;  kk++;
      }
    }
  }

  delete [] order;
}



/**
  The function decomposes %matrix A into L.U form.
   
  @param[in,out] a - the given %matrix,
  @param[out] x - left hand side %vector,
  @param[in,out] y - right hand side %vector,
  @param[in] n - number of rows/columns of %matrix,
  @param[in] zero - computer zero,
  @param[in] tc - type of computation.
   
  tc = 1 - decomposition of matrix and back substitution
  tc = 2 - decomposition of matrix
  tc = 3 - back substitution
   
  @author JK, 6.1.2000
*/
void lu_full (double *a,double *x,double *y,long n,double zero,long tc)
{
  long i,j,k;
  double s;
  char emsg[200];
  
  if (tc==1 || tc==2){
    for (i=0;i<n;i++){
      for (j=0;j<i;j++){
        s=0.0;
        for (k=0;k<j;k++){
          s+=a[i*n+k]*a[k*n+j];
        }
        a[i*n+j]-=s;
      }
      
      s=0.0;
      for (k=0;k<i;k++){
        s+=a[i*n+k]*a[k*n+i];
      }
      a[i*n+i]-=s;
      if (fabs(a[i*n+i])<zero){
        sprintf(emsg, "zero diagonal entry (%e) in LU decomposition in matrix row %ld", a[i*n+i], i+1);
        print_err(emsg, __FILE__, __LINE__, __func__);
      }
      
      for (j=i+1;j<n;j++){
        s=0.0;
        for (k=0;k<i;k++){
          s+=a[i*n+k]*a[k*n+j];
        }
        a[i*n+j]=(a[i*n+j]-s)/a[i*n+i];
      }
    }
  }
  
  if (tc==1 || tc==3){
    for (i=0;i<n;i++){
      s=0.0;
      for (j=0;j<i;j++){
        s+=a[i*n+j]*y[j];
      }
      y[i]=(y[i]-s)/a[i*n+i];
    }
    for (i=n-1;i>-1;i--){
      s=0.0;
      for (j=i+1;j<n;j++){
        s+=a[i*n+j]*x[j];
      }
      x[i]=y[i]-s;
    }
  }
  
}



/**
  The function assembles %matrix of least square problem which is used in extrapolation 
  of values from integration points to nodes.
   
  @param[out] lsm - array containing %matrix,
  @param[in] natcoord - %vector containing natural coordinates of integration points.
   
  @author Created by JK, 10.5.2002
*/
void matassem_lsm (double *lsm,vector &natcoord)
{
  if (natcoord.n==1){
    lsm[0*2+0]++;  lsm[0*2+1]+=natcoord[0];
                   lsm[1*2+1]+=natcoord[0]*natcoord[0];
  }
  if (natcoord.n==2){
    lsm[0*3+0]++;  lsm[0*3+1]+=natcoord[0];              lsm[0*3+2]+=natcoord[1];
                   lsm[1*3+1]+=natcoord[0]*natcoord[0];  lsm[1*3+2]+=natcoord[1]*natcoord[0];
                                                         lsm[2*3+2]+=natcoord[1]*natcoord[1];
  }
  if (natcoord.n==3){
    lsm[0*4+0]++;  lsm[0*4+1]+=natcoord[0];              lsm[0*4+2]+=natcoord[1];              lsm[0*4+3]+=natcoord[2];
                   lsm[1*4+1]+=natcoord[0]*natcoord[0];  lsm[1*4+2]+=natcoord[1]*natcoord[0];  lsm[1*4+3]+=natcoord[2]*natcoord[0];
                                                         lsm[2*4+2]+=natcoord[1]*natcoord[1];  lsm[2*4+3]+=natcoord[2]*natcoord[1];
                                                                                               lsm[3*4+3]+=natcoord[2]*natcoord[2];
  }
}



/**
  The function assembles right hand side of least square problems.
  It is used in extrapolation of values from integration points to nodes.
   
  @param[out] rhs - array containing %matrix
  @param[in] natcoord - %vector of natural coordinates of integration points
  @param[in] values - %vector of values in integration point
   
  @author Created by JK, 10.5.2002
*/
void rhsassem_lsm (double *rhs,vector &natcoord,vector &values)
{
  long i;
  if (natcoord.n==1){
    for (i=0;i<values.n;i++){
      rhs[i*2+0]+=values[i];
      rhs[i*2+1]+=values[i]*natcoord[0];
    }
  }
  if (natcoord.n==2){
    for (i=0;i<values.n;i++){
      rhs[i*3+0]+=values[i];
      rhs[i*3+1]+=values[i]*natcoord[0];
      rhs[i*3+2]+=values[i]*natcoord[1];
    }
  }
  if (natcoord.n==3){
    for (i=0;i<values.n;i++){
      rhs[i*4+0]+=values[i];
      rhs[i*4+1]+=values[i]*natcoord[0];
      rhs[i*4+2]+=values[i]*natcoord[1];
      rhs[i*4+3]+=values[i]*natcoord[2];
    }
  }
}



/**
  The function solves least square problems.
  It is used in extrapolation of values from integration points to nodes
   
  @param[in,out] lsm - array of %matrix of least square problem,
  @param[out] lhs - array of left hand side (contains unknown coefficients),
  @param[in,out] rhs - array of right hand side,
  @param[in] zero - computer zero,
  @param[in] n - number of rows/columns of %matrix,
  @param[in] m - number of values.
   
  @author Created by JK, 10.5.2002
*/
void solve_lsm (double *lsm, double *lhs, double *rhs, double zero, long n, long m)
{
  long i,j,k;

  if (n==3){
    lsm[1*3+0]=lsm[0*3+1];
    lsm[2*3+0]=lsm[0*3+2];  lsm[2*3+1]=lsm[1*3+2];
  }
  if (n==4){
    lsm[1*4+0]=lsm[0*4+1];
    lsm[2*4+0]=lsm[0*4+2];  lsm[2*4+1]=lsm[1*4+2];
    lsm[3*4+0]=lsm[0*4+3];  lsm[3*4+1]=lsm[1*4+3];  lsm[3*4+2]=lsm[2*4+3];
  }
  
  if (lsm[0]<3.0){
    k=0;
    for (i=0;i<m;i++){
      lhs[k]=rhs[k]/lsm[0];
      k++;
      for (j=1;j<n;j++){
        lhs[k]=0.0;
        k++;
      }
    }
  }
  else{
    lu_full (lsm,lhs,rhs,n,zero,2);
    for (i=0;i<m;i++){
      lu_full (lsm,lhs+i*n,rhs+i*n,n,zero,3);
    }
  }
  
}




/**
  The function extracts a squared matrix with ncomp components from %matrix b
  and put them into %matrix a.

  @param[in] a - %matrix containing extracted components,
  @param[out] b - %matrix from where components are extracted,
  @param[in] fi - first index,
  @param[in] ncomp - number of components.

  @author Created by JK 25.2.2004
*/
void extractm (matrix &a,matrix &b,long fi, long ncomp)
{
  long i,j;
  
  for (i=0;i<ncomp;i++){
    for (j=0;j<ncomp;j++){
      a[i][j]=b[i+fi][j+fi];
    }
  }
}



/**
  The function picks up submatrix b from the original %matrix a.
   
  @param[in] a - array of original %matrix of type A(n,n)
  @param[out] b - array of selected submatrix B(nr,nc)
  @param[in] n - number of columns in %matrix A
  @param[in] nr - number of rows in the %matrix B
  @param[in] nc - number of columns in the %matrix B
  @param[in] fri, first row index of submatrix in the original %matrix
  @param[in] fci - first column index of submatrix in the original %matrix
   
  @author Created by JK, 9.12.2004
*/
void extractblock (double *a,double *b,long n,long nr,long nc,long fri,long fci)
{
  long i,j;
  
  for (i=0;i<nr;i++){
    for (j=0;j<nc;j++){
      b[i*nc+j]=a[(fri+i)*n+(fci+j)];
    }
  }
}



/**
  The function picks up submatrix b from the original %matrix a.
   
  @param[in] a - original %matrix of type A(nra,nca),
  @param[out] b - selected submatrix B(nrb,ncb),
  @param[in] fri - first row index of submatrix in the original %matrix,
  @param[in] fci - first column index of submatrix in the original %matrix.

  @pre Both matrices must be allocated, fri+nrb-1 <= nra, fci+ncb-1 <= nca
   
  @author Created by TKo, 02.2020
*/
void extractblock (matrix &a, matrix &b, long fri,long fci)
{
  long i,j;
  
  for (i=0; i<b.m; i++){
    for (j=0; j<b.n; j++){
      b(i,j) = a(fri+i, fci+j);
    }
  }
}



/**
  The function stores submatrix b to the destination %matrix a.
   
  @param[in] a - destination %matrix of type A(nra,nca) where the submatrix B will be stored,
  @param[out] b - submatrix B(nrb,ncb) to be stored,
  @param[in] fri - first row index of submatrix in the destination %matrix,
  @param[in] fci - first column index of submatrix in the destination %matrix.

  @pre Both matrices must be allocated, fri+nrb-1 <= nra, fci+ncb-1 <= nca
   
  @author Created by TKo, 02.2020
*/
void storeblock (matrix &a, matrix &b, long fri,long fci)
{
  long i,j;
  
  for (i=0; i<b.m; i++){
    for (j=0; j<b.n; j++){
      a(fri+i, fci+j) = b(i,j);
    }
  }
}



/**
  The function extracts the i-th row of the %matrix m and stores it to %vector dest.
  @param[in] m - the given %matrix whose i-th row will be extracted,
  @param[in] i - index of the extracted row,
  @aparam[out] dest - destination %vector where the extracted row will be stored to

  @return The function returns extracted row in the argument dest.
  @retval 0 - on success,
  @retval 1 - in the case of incompatible dimensions of %matrix m and %vector dest.

  @author Created by TKo, 02.2024
*/
long extractrow(const matrix &m, long i, vector &dest)
{
  if (m.n != dest.n)
  {
    print_err("cannot extract %ld-th matrix row - incompatible dimensions of the destination vector\n"
              "m(%ld,%ld) X dest(%ld)", __FILE__, __LINE__, __func__, m.m, m.n, dest.n);
    return 1;
  }
  for (long j=0; j<m.n; j++)
    dest[j] = m(i, j);

  return 0;
}



/**
  The function extracts the j-th column of the %matrix m and stores it to %vector dest.
  @param[in] m - the given %matrix whose j-th column will be extracted,
  @param[in] i - index of the extracted column,
  @aparam[out] dest - destination %vector where the extracted column will be stored to

  @return The function returns extracted column in the argument dest.
  @retval 0 - on success,
  @retval 1 - in the case of incompatible dimensions of %matrix m and %vector dest.

  @author Created by TKo, 02.2024
*/
long extractcol(const matrix &m, long j, vector &dest)
{
  if (m.n != dest.n)
  {
    print_err("cannot extract %ld-th matrix row - incompatible dimensions of the destination vector\n"
              "m(%ld,%ld) X dest(%ld)", __FILE__, __LINE__, __func__, m.m, m.n, dest.n);
    return 1;
  }
  for (long i=0; i<m.m; i++)
    dest[i] = m(i, j);

  return 0;
}



/**
   function diagonalizes a %matrix m
   
   JK, 14. 6. 2019
*/
void diagonalization (matrix &mat)
{
  long i,j;
  double s;
  
  if (mat.m != mat.n)
    print_err("matrix is not a square matrix (%ld,%ld)", __FILE__, __LINE__, __func__, mat.m, mat.n);
  
  for (i=0;i<mat.m;i++){
    s=0.0;
    for (j=0;j<mat.n;j++){
      s+=mat[i][j];
      mat[i][j]=0.0;
    }
    mat[i][i]=s;
  }
  
}



/**
  The function normalizes components of the given %matrix by the maximum 
  value from the absolute values of its components.

  @param[in,out] mat - matrix to be normalized

  @return The function returns factor of the matrix component normalization, i.e. 
          reciprocal of the maximum matrix component.
  @author Created by TKo, 07.2020
 */
double normalize(matrix &mat)
{
  long i;
  double max = 0.0;
  
  for (i=0; i<mat.m*mat.n; i++){
    if (fabs(mat.a[i]) > max)  max = fabs(mat.a[i]);
  }
  if (max == 0.0)
    return 1.0;
  
  cmulm(1.0/max, mat);

  return 1.0/max;
}



/**
  The function computes Euclidean norm of given %matrix components.

  @param a - given matrix for the norm computation


  @retval Value of Euclidean %matrix norm.
*/
double norm(matrix &a)
{
  double ret = 0.0;
  long i, j;

  for(i=0; i<a.m; i++){
    for(j=0; j<a.n; j++)
      ret += a(i,j)*a(i,j);
  }
  return sqrt(ret);
}



/**
   Function computes the largest eigenvalue of a square %matrix A
   
   @param err - required error
   @param ni - maximum number of iterations
   @param zero - computer zero
   
   24. 3. 2021, JK
*/
double power_method (matrix &a,double err,long /*ni*/,double zero)
{
  long i,j;
  double rho,rhoprev,norv;
  
  if (a.n != a.m){
    print_err("matrix is not a square matrix (%ld,%ld)", __FILE__, __LINE__, __func__, a.m, a.n);
    abort ();
  }
  
  vector u(a.n),v(a.n);
  fillv (1.0,v);
  
  rhoprev=0.0;
  
  for (i=0;i<a.n;i++){

    // norm of the vector v
    norv = normv (v);
    
    if (norv<zero){
      print_err("norm of the approximation of the eigenvector is zero", __FILE__, __LINE__, __func__);
      break;
    }
    
    for (j=0;j<a.n;j++){
      u[j] = v[j]/norv;
    }

    //  v is a new approximation of eigenvector
    mxv (a,u,v);
    
    //  Rayleigh quotient
    scprd (u,v,rho);
    
    if (rho-rhoprev<err)
      break;
    
  }
  return rho;
}



/**
   function computes all eigenvalues and eigenvectors of the matrix
   matrix is stored as dense matrix and will be overwritten!
   
   eigenvectors are columns of the matrix evec but matrix is stored
   in usual way, by rows.
   
   @param a - array containing matrix
   @param evec - array containing eigenvectors
   @param eval - array containing eigenvalues
   @param n - order of matrix a (number of rows or columns)
   @param ni - maximum number of iterations
   @param ani - number of performed iterations
   @param limit - maximum acceptable absolute value of offdiagonal element
   @param normalize - switch for normalizing of matrix (on>0, off=0)
   
   !!!! nutno predelat trideni vlastnich vektoru
   
   23.7.2001
   updated 17.8.2007 by TKo. - normalizing of matrix
*/
void jacobi_rot (double *a,double *evec,double *eval,long n,long ni,long &ani,double limit, long normalize)
{
  long i, j, k, l, ii, jj, en, nn;
  double c, s, q, r, t, ai, aj, maxa, av;
  
  if (normalize){
    maxa = 0.0; 
    for (i=0; i<n*n; i++){ // find maximal nonzero element of matrix
      if (a[i] == 0.0)
        continue;
      av = fabs(a[i]);
      if (av > maxa)
        maxa = av;
    }
    if (maxa != 0.0){
      for (i=0; i<n*n; i++) // normalizing
        a[i] /= maxa;
    }
  }
  en=(n*n-n)/2;
  //  initial values
  k=0;
  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      evec[k]=0.0;  k++;
    }
  }
  for (i=0;i<n;i++){
    evec[i*n+i]=1.0;
  }

  //  main iteration loop
  for (k=0;k<ni;k++){
    
    nn=0;
    for (i=0;i<n;i++){
      for (j=i+1;j<n;j++){
	if (fabs(a[i*n+j])<limit){
	  nn++;  continue;
	}
	else{
	  q=(a[j*n+j]-a[i*n+i])/2.0/a[i*n+j];
	  r=sqrt(1.0+q*q);
	  if (q>=0.0)  t=-q+r;
	  else         t=-q-r;
	  r=sqrt(1.0+t*t);
	  c=1.0/r;
	  s=t*c;
	  
	  //  A.J
	  ii=i;  jj=j;
	  for (l=0;l<n;l++){
	    ai=a[ii];
	    aj=a[jj];
	    a[ii]=ai*c-aj*s;
	    a[jj]=ai*s+aj*c;
	    ii+=n;  jj+=n;
	  }

	  //  J.A.J
	  ii=i*n;  jj=j*n;
	  for (l=0;l<n;l++){
	    ai=a[ii];
	    aj=a[jj];
	    a[ii]=ai*c-aj*s;
	    a[jj]=ai*s+aj*c;
	    ii++;  jj++;
	  }

	  //  Q=Q.J
	  ii=i;  jj=j;
	  for (l=0;l<n;l++){
	    ai=evec[ii];
	    aj=evec[jj];
	    evec[ii]=ai*c-aj*s;
	    evec[jj]=ai*s+aj*c;
	    ii+=n;  jj+=n;
	  }

	}
      }
    }
    if (nn==en){
      break;
    }
  }

  //  number of performed iterations
  ani=k;

  for (i=0;i<n;i++){
    eval[i]=a[i*n+i];
  }
  if ((normalize) && (maxa != 0.0))
  {
    for (i=0;i<n;i++)
      eval[i]*=maxa;
  }
  /*
  //  sorting eigenvalues and eigenvectors according to magnitudes of eigenvalues
  for (i=n-1;i>=0;i--){
    c=0.0;
    for (j=i;j>=0;j--){
      s=fabs(a[j*n+j]);
      if (s>c){
	c=s;  k=j;
      }
    }
    eval[i]=c;
    a[k*n+k]=a[i*n+i];

    ii=i;  jj=k;
    for (j=0;j<n;j++){
      t=evec[ii];
      evec[ii]=evec[jj];
      evec[jj]=t;
      ii+=n;  jj+=n;
    }
  }
  */
}
