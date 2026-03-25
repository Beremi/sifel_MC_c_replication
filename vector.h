#ifndef DEF_VECTOR
#define DEF_VECTOR
#include <stdio.h>
#include <assert.h>
#include <string.h>

//#define DEBUG_VECTOR



/**
  Structure with bit array of flags for handling with vectors.

  @date Created 4.5.2015 by TKo.

  @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz
*/
struct vecstat
{
  unsigned dealloc:1;  ///< flag for vector memory deallocation (0=no deallocation, 1=release memory by delete)
  vecstat() {dealloc=0;};
};



/**
  This struct implements %vector with real elements type of double.
  There are also declarations of the functions for the %vector computing.

  @date Created  29.5.2000 by TKo.

  @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz
*/
struct vector
{
  long n;       ///< the number of %vector entries (components)
  double *a;    ///< pointer to memory ,where are stored elements of %vector
  long size;    ///< real length of array a (due to memory reallocation)
  vecstat stat; ///< bit array of vector status flags (deallocation)

  inline vector() { n = size = 0L; a = NULL;}; ///< default constructor
  vector(long n);                              ///< allocating constructor
  vector(const vector &v);                     ///< copy constructor

  /**
    The constructor assigns preallocated memory to the member array a.
   
    @param[in] n - the number of %vector elements,
    @param[in] dealloc - the flag for deallocation control (0=no deallocation ,1=release memory by delete operator),
    @param[in] ptr - the pointer to the allocated memory which will be used for storage of n elements.
   
    @date Created  4.5.2015 by TKo.

    @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz
  */
  inline vector (long n, unsigned dealloc, double *ptr)
  {
    vector::n = size = n;
    a = ptr;
    stat.dealloc = dealloc;

    memset (a, 0, n*sizeof(*a));

   #ifdef DEBUG_VECTOR
    Acv++;
    Ava += size;
    if (Ava > Avmax)
      Avmax = Ava;
   #endif
  };


  /**
    The function operator enables access to the elements of the member array a.
     
    @param[in] i - the number of the desired %vector %element.

    @retval The function returns the i-th element of the given %vector.
     
    @date Created  29.5.2000 by TKo,
    @date Modified  9.5.2001 by TKo.

    @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz
  */
  inline double& operator () (long i)
  {
    #ifdef DEBUG_VECTOR
    if ((i >= n) || (i < 0))
       fprintf(stderr, "vector row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, n-1);

     assert(i < n);
    #endif
    return (a[i]);
  };
  
  /**
    The function operator enables CONSTANT access to the elements of the member array a.
     
    @param[in] i - the number of the required %vector %element.

    @retval The function returns the i-th element of the given %vector.
     
    @date Created  29.5.2000 by TKo, 
    @date modified  9.5.2001, TKo.

    @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz
  */
  inline double& operator () (long i) const
  {
    #ifdef DEBUG_VECTOR
    if ((i >= n) || (i < 0))
       fprintf(stderr, "vector row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, n-1);

     assert(i < n);
    #endif
    return (a[i]);
  };
  
  /**
    The array index operator enables access to the elements of the member array a.
     
    @param[in] i - the number of the required %vector %element.

    @retval The function returns the i-th %element of the given %vector.
     
    @date Created  29.5.2000 by TKo,
    @date modified  9.5.2001 by TKo.

    @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz
  */
  inline double& operator [] (long i)
  {
    #ifdef DEBUG_VECTOR
     if  ((i >= n) || (i < 0))
       fprintf(stderr, "vector row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, n-1);

     assert(i < n);
    #endif
    return (a[i]);
  };
  
  /**
    The array index operator enables CONSTANT access to the elements of the member array a
     
    @param[in] i - the number of the required %vector %element.

    @retval The function returns the i-th %element of the given %vector.
     
    @date Created  29.5.2000 by TKo,
    @date modified  9.5.2001 by TKo.

    @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz
  */
  inline double& operator [] (long i) const
  {
    #ifdef DEBUG_VECTOR
     if  ((i >= n) || (i < 0))
       fprintf(stderr, "vector row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, n-1);

     assert(i < n);
    #endif
    return (a[i]);
  };

  inline ~vector() ///< destructor
  {
   #ifdef DEBUG_VECTOR
    if (a != NULL)
    {
      Acv--;
      Ava -= size;
    }
   #endif

    if (stat.dealloc)
      delete [] a;
  };
  /// reallocates the given %vector to the required number of components
  long reallocv(long n);
  /// reallocates %vector with help of preallocated memory to the required number of components
  long reallocv(long n, unsigned dealloc, double *ptr);
  /// makes the given %vector to referenece to the %vector src
  long makerefv(const vector &src);
  /// makes the given %vector to referenece to the array ptr
  long makerefv(long nc, double *ptr);
};



/**
  This struct implements %vector with integer elements type of long.

  @date Created 29.5.2000 by TKo.

  @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz
*/
struct ivector
{
  long n;   ///< number of %vector elements
  long *a;  ///< pointer to memory ,where are stored elements of %vector
  long size; ///< real length of array a (due to memory reallocation)
  vecstat stat; ///< bit array of vector status flags (deallocation)
  
  inline ivector() { n = size = 0L; a = NULL;}; ///< default constructor
  ivector(long n);                      ///< allocating constructor
  ivector(const ivector &v);            ///< copy constructor
  ivector& operator=(const ivector &v); ///< copy assignment operator

  /**
    The constructor assigns preallocated memory to the member array a.
   
    @param[in] n - the number of %vector elements,
    @param[in] dealloc - the flag for deallocation control (0=no deallocation ,1=release memory by delete operator),
    @param[in] ptr - the pointer to the allocated memory which will be used for storage of n elements.
   
    @date Created 4.5.2015 by TKo.

    @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz
  */
  inline ivector (long n, unsigned dealloc, long *ptr)
  {
    ivector::n = size = n;
    a = ptr;
    stat.dealloc = dealloc;

    memset (a, 0, n*sizeof(*a));
  
   #ifdef DEBUG_VECTOR
    Acv++;
    Ava += size;
    if (Ava > Avmax)
      Avmax = Ava;
   #endif
  };
  
  /**
    The function operator allows for access to the elements of the member array a.
     
    @param[in] i - the number of the required %vector %element.

    @retval The function returns the i-th element of the given %vector.
     
    @date Created  29.5.2000 by TKo,
    @date modified  9.5.2001 by TKo.

    @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz
  */
  inline long& operator () (long i)
  {
    #ifdef DEBUG_VECTOR
     if  ((i >= n) || (i < 0))
       fprintf(stderr, "ivector row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, n-1);

     assert(i < n);
    #endif
    return (a[i]);
  };

  /**
    The function operator enables CONSTANT access to the elements of the member array a.
     
    @param[in] i - the number of the required %vector %element.

    @retval The function returns the i-th element of the given %vector.
     
    @date Created  29.5.2000 by TKo, 
    @date modified  9.5.2001 by TKo.

    @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz
  */
  inline long& operator () (long i) const
  {
    #ifdef DEBUG_VECTOR
     if  ((i >= n) || (i < 0))
       fprintf(stderr, "ivector row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, n-1);

     assert(i < n);
    #endif
    return (a[i]);
  };

  /**
    The array index operator enables access to the elements of the member array a.
   
    @param[in] i - the number of the required %vector %element.

    @retval The function returns the i-th %element of the given %vector.
   
    @date Created  29.5.2000 by TKo,
    @date modified  9.5.2001 by TKo.
 
    @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz
  */
  inline long& operator [] (long i)
  {
    #ifdef DEBUG_VECTOR
     if  ((i >= n) || (i < 0))
       fprintf(stderr, "ivector row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, n-1);

     assert(i < n);
    #endif
    return (a[i]);
  };

  /**
    The array index operator enables CONSTANT access to the elements of the member array a.
   
    @param[in] i is the number of the required %vector %element.

    @retval The function returns the i-th %element of the given %vector.
   
    @date Created  29.5.2000 by TKo,
    @date modified  9.5.2001 by TKo.
   
    @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz
  */
  inline long& operator [] (long i) const
  {
    #ifdef DEBUG_VECTOR
     if  ((i >= n) || (i < 0))
       fprintf(stderr, "ivector row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, n-1);

     assert(i < n);
    #endif
    return (a[i]);
  };

  /// makes the given %vector to referenece to the %vector src
  long makerefv(const ivector &src);
  /// makes the given %vector to referenece to the array ptr
  long makerefv(long nc, long *ptr);
  /** reallocates the given ivector to the given size n and fills all components with a sequence
      starting with the initial value init incremented by incr */
  long realloc(long nc, long init, long incr);
  /// fills the given %vector components starting from fi with a sequence with initial value init incremented by incr
  void fillseq(long fi, long nc, long init, long incr);

  inline ~ivector() ///< destructor
  {
   #ifdef DEBUG_VECTOR
    if (a != NULL)
    {
      Aciv--;
      Aiva -= size;
    }
   #endif

    if (stat.dealloc)
      delete [] a;
  };
};



/**
  This struct implements %vector with logical elements type of bool.

  @date Created  29.5.2000 by TKo.

  @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz
*/
struct bvector
{
  long n;   ///< number of %vector elements
  bool *a;  ///< pointer to memory ,where are stored elements of %vector
  long size; ///< real length of array a (due to memory reallocation)
  vecstat stat; ///< bit array of vector status flags (deallocation)
  
  inline bvector() { n = size = 0L; a = NULL;}; ///< default constructor
  bvector(long n);                      ///< allocating constructor
  bvector(const bvector &v);            ///< copy constructor
  bvector& operator=(const bvector &v); ///< copy assignment operator


  /**
    The constructor assigns preallocated memory to the member array a.
   
    @param[in] n - the number of %vector elements,
    @param[in] dealloc - the flag for deallocation control (0=no deallocation ,1=release memory by delete operator),
    @param[in] ptr - the pointer to the allocated memory which will be used for storage of n elements.
   
    @date Created  4.5.2024 by TKo.

    @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz
  */
  inline bvector (long n, unsigned dealloc, bool *ptr)
  {
    bvector::n = size = n;
    a = ptr;
    stat.dealloc = dealloc;

    memset (a, 0, n*sizeof(*a));
  
   #ifdef DEBUG_VECTOR
    Acv++;
    Ava += size;
    if (Ava > Avmax)
      Avmax = Ava;
   #endif
  };
  
  /**
    The function operator enables access to the elements of the member array a.
     
    @param[in] i - the number of the required %vector element.

    @retval The function returns the i-th element of the given %vector.
     
    @date Created  29.5.2000 by TKo,
    @date modified  9.5.2001 by TKo.

    @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz
  */
  inline bool& operator () (long i)
  {
    #ifdef DEBUG_VECTOR
     if  ((i >= n) || (i < 0))
       fprintf(stderr, "bvector row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, n-1);

     assert(i < n);
    #endif
    return (a[i]);
  };

  /**
    The function operator enables CONSTANT access to the elements of the member array a.
     
    @param[in] i - the number of the required %vector element.

    @retval The function returns the i-th element of the given %vector.
     
    @date Created  29.5.2000 by TKo,
    @date modified  9.5.2001 by TKo.

    @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz
  */
  inline bool& operator () (long i) const
  {
    #ifdef DEBUG_VECTOR
     if  ((i >= n) || (i < 0))
       fprintf(stderr, "bvector row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, n-1);

     assert(i < n);
    #endif
    return (a[i]);
  };

  /**
    The array index operator enables access to the elements of the member array a.
   
    @param[in] i - the number of the required %vector %element.

    @retval The function returns the i-th %element of the given %vector.
   
    @date Created  29.5.2000 by TKo,
    @date modified  9.5.2001 by TKo.

    @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz
  */
  inline bool& operator [] (long i)
  {
    #ifdef DEBUG_VECTOR
     if  ((i >= n) || (i < 0))
       fprintf(stderr, "bvector row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, n-1);

     assert(i < n);
    #endif
    return (a[i]);
  };

  /**
    The array index operator enables CONSTANT access to the elements of the member array a
   
    @param[in] i - the number of the required %vector %element.

    @retval The function returns the i-th %element of the given %vector.
   
    @date Created  29.5.2000 by TKo,
    @date modified  9.5.2001 by TKo.
  */
  inline bool& operator [] (long i) const
  {
    #ifdef DEBUG_VECTOR
     if  ((i >= n) || (i < 0))
       fprintf(stderr, "bvector row index %ld is out of range <0,%ld>", __FILE__, __LINE__, __func__, i, n-1);

     assert(i < n);
    #endif
    return (a[i]);
  };

  inline ~bvector() ///< destructor
  {
   #ifdef DEBUG_VECTOR
    if (a != NULL)
    {
      Aciv--;
      Aiva -= size;
    }
   #endif

    if (stat.dealloc)
      delete [] a;
  };
};



#ifdef DEBUG_VECTOR
 long give_acv();
 long give_avmax();
#endif

/// allocates %vector
long allocv(long n, vector &vec);
/// allocates %ivector
long allocv(long n, ivector &vec);
/// creates %vector which contains reference to %vector src
long makerefv(vector &ref, vector &src);
/// creates %vector which contains reference to array ptr
long makerefv(vector &ref, double *ptr, long n);
/// creates %vector which contains reference to %vector src
long makerefv(ivector &ref, ivector &src);
/// creates %vector which contains reference to array ptr
long makerefv(ivector &ref, long *ptr, long n);
/// reallocates %vector
long reallocv(long n, vector &vec);
/// reallocates %vector with help of preallocated memory
long reallocv(long n, vector &vec, unsigned dealloc, double *ptr);
/// reallocates %ivector
long reallocv(long n, ivector &vec);
/// reallocates %ivector with help of preallocated memory
long reallocv(long n, ivector &vec, unsigned dealloc, long *ptr);
/// reallocates %bvector
long reallocv(long n, bvector &vec);
/// reallocates %ivector with help of preallocated memory
long reallocv(long n, bvector &vec, unsigned dealloc, long *ptr);

/// copies contents of %vector
long copyv(const vector &src, vector &dest);
/// copies contents of %ivector
long copyv(const ivector &src, ivector &dest);
/// copies contents of %bvector
long copyv(const bvector &src, bvector &dest);
/// copies contents of array
long copyv(const long *src,long *dest,long n);
/// copies contents of array
long copyv(const double *src, double *dest, long n);
/// copies contents of array
long copyv(const bool *src, bool *dest, long n);
/// copies contents of array to %vector
long copyv(const double *src, vector &dest);
/// copies contents of array to %ivector
long copyv(const long *src, ivector &dest);
/// copies contents of array to %bvector
long copyv(const bool *src, bvector &dest);
/// copies contents of %vector to array
long copyv(const vector &src, double *dest);
/// copies contents of %ivector to array
long copyv(const ivector &src, long *dest);
/// copies contents of %bvector to array
long copyv(const bvector &src, bool *dest);

/// swaps pointers to component arrays of two vectors, i.e. swaps their content
long swapv(vector &a, vector &b);

/// copies range of components between vectors
long rcopyv(const vector &src, long src_fi, vector &dest, long dest_fi, long n);
/// copies range of components between vectors given by arrays
long rcopyv(const double *src, long src_fi, long src_n, double *dest, long dest_fi, long dest_n, long n);

/// copies contents of %vector multiplied by scalar
long copymultv(const vector &src, vector &dest, double c);
/// copies contents of %ivector multiplied by scalar
long copymultv(const ivector &src, ivector &dest, long c);
/// copies contents of array multiplied by scalar
long copymultv(const double *src, double *dest, double c, long n);

/// fills contents of %vector with value c
long fillv(double c, vector &vec);
/// fills contents of %ivector with value c
long fillv(long c, ivector &vec);
/// fills contents of %bvector with value c
long fillv(bool c, bvector &vec);
/// fills contents of integer %vector given by array with value c
void fillv(long c, long *vec,long n);
/// fills contents of %vector given by array with value c
void fillv(double c, double *vec,long n);

/// fills array with zeros
long nullv (double *a, long n);
/// fills integer array with zeros
long nullv (long *a, long n);
/// fills %vector with zeros
long nullv (vector &a);
/// fills %ivector with zeros
long nullv (ivector &a);
/// fills %bvector with false values
long nullv (bvector &a);

/// changes sign of %vector components
void chsgnv(vector &v);
/// changes sign of array components
void chsgnv(double *a, long n);

/// deallocates %vector
long destrv(vector &vec);
/// deallocates %ivector
long destrv(ivector &vec);

/// adds 2 vectors, c = a + b
long addv(const vector &a, const vector &b, vector &c);
/// adds vectors, a += b
long addv(vector &a, const vector &b);
/// computes expression c = a + bc*b with vectors a, b and c and real constant bc
void addmultv(const vector &a, const vector &b, double bc, vector &c);
/// computes expression a += bc*b with vectors a and b and real constant bc
void addmultv(vector &a, const vector &b, double bc);

/// adds 2 ivectors c = a + b
long addv(const ivector &a, const ivector &b, ivector &c);

/// adds 2 double arrays a += b
void addv(double *a, const double *b, long n);

/// adds 2 double arrays c = a + b
void addv(const double *a, const double *b, double *c, long n);

/// computes expression a += bc*b with double arrays a and b and real constant bc
void addmultv (double *a, const double *b, double bc, long n);
/// computes expression c = a + bc*b with double arrays a, b and c and real constant bc
void addmultv (const double *a, const double *b, double bc, double *c, long n);

/// computes expression c = ac*a + bc*b with vectors a, b and c and real constants ac and bc
void addmultv (const vector &a, double ac, const vector &b, double bc, vector &c);
/// computes expression a = ac*a + bc*b with double arrays a and b and real constants ac and bc
void addmultv (double *a, double ac, const double *b, double bc, long n);
/// computes expression c = ac*a + bc*b with double arrays a, b and c and real constants ac and bc
void addmultv (const double *a, double ac, const double *b, double bc, double *c, long n);

/// subtracts 2 vectors, c = a - b
long subv(const vector &a, const vector &b, vector &c);
/// subtracts 2 vectors, a -= b
long subv(vector &a, const vector &b);
/// subtracts 2 ivectors, c = a-b
long subv(const ivector &a, const ivector &b, ivector &c);

/// subtracts 2 arrays, a -= b
void subv(double *a, const double *b, long n);
/// subtracts 2 arrays, result is stored to the another array, c = a-b
void subv (const double *a, const double *b, double *c, long n);

/// performs cross product of 2 vectors, c = a x b
long crprd(const vector &a, const vector &b, vector &c);
/// performs cross product of 2 ivectors, c = a x b
long crprd(const ivector &a, const ivector &b, ivector &c);

/// performs scalar product of 2 vectors, scalar = a . b
long scprd(const vector &a, const vector &b, double &scalar);
/// returns scalar product a . b
double scprd(const vector &a, const vector &b);
/// performs scalar product of 2 ivectors, scalar = a . b
long scprd(const ivector &a, const ivector &b, long &scalar);
/// returns scalar product of 2 arrays, a . b
double scprd(const double *a, const double *b, long n);
/// returns scalar product of 2 vectors, a . b computed on several threads 
double scprd_t(const vector &a, const vector &b, long nt);
/// returns scalar product of two arrays a . b computed on several threads 
double scprd_t(const double *a, const double *b, long n, long nt);
/// performs scalar product of 2 arrays, a . b
double ss (double *a,double *b,long n);

/// multiplies %vector by constant, v = c u
long cmulv(double c, const vector &u, vector &v);
/// multiplies %ivector by constant of type double, v = c u
long cmulv(double c, const ivector &u, vector &v);
/// multiplies %ivector by constant of type long, v = c u
long cmulv(long c, const ivector &u, ivector &v);

/// multiplies %vector by constant, a *= c
long cmulv(double c, vector &a);
/// multiplies %ivector by constant of type long, a *= c
long cmulv(long c, ivector &a);
/// multiplies double array by constant, a *= c
void cmulv(double c, double *a, long n);
/// multiplies double array by constant, u = c a
void cmulv(double c, double *a, double *u, long n);

/// computes norm of %vector, |a|
double normv(const vector &a);
/// computes norm of %vector, |a|
double normv(const double *a, long n);
/// computes norm of %ivector, |a|
double normv(const ivector &a);

/// computes selective norm of %vector
double snormv(const vector &a, long fi, long nc);
/// computes selective norm of %vector
double snormv(const double *a, long n, long fi, long nc);

/// normalizes %vector components so that after function call |a| = 1
void normalize(vector &a);
/// normalizes %vector components and stores them into b, |b| = 1
void normalize(const vector &a, vector &b);
/// normalizes %vector components so that after function call |a| = 1
void normalize(double *a, long n);
/// normalizes %vector components and stores them into b, |b| = 1
void normalize(double *a, long n, double *b);

/// returns maximum of vector components
//double maxcompv(const vector &a);
//double maxcompv(const double *a, long n);

/// cosine angle of 2 vectors, cos = a . b/(|a||b|)
long cosav(const vector &a, const vector &b, double &cos);
/// cosine angle of 2 ivectors, cos = a . b/(|a||b|)
long cosav(const ivector &a, const ivector &b, double &cos);

/// absolute value of vector components
void absv(vector &a);
/// absolute value of vector components
void absv(ivector &a);

/// prints contents of %vector to the file out with given precision and field width
long printv(const vector &u, FILE *out = stdout, int prec = 3, int width = 11);
/// prints contents of %vector given by array to the file out with given precision and field width
long printv(const double *u, long n, FILE *out = stdout, int prec = 3, int width = 11);
/// prints %vector u to the file out
long printv(FILE *out, vector &u);
/// prints contents of %ivector
long printv(const ivector &u, FILE *out = stdout, int width = 11);
/// prints contents of %bvector
long printv(const bvector &u, FILE *out = stdout, int width = 11);
/// prints contents of integer %vector given by array
long printv(const long *u, long n, FILE *out = stdout, int width = 11);

/// reads contents of %vector from string
long readv(char *in, vector &u);
/// reads contents of integer %vector from string
long readv(char *in, ivector &u);
/// reads contents of integer %bvector from string
long readv(char *in, bvector &u);

/// prints complete record of %vector u to the file out with given precision and field width
long printva(const vector &u, FILE *out = stdout, int prec = 3, int width = 11);
/// prints complete record of %vector u to the file out
long printva(FILE *out, const vector &u);


/// function extracts ncomp components (from index fi) from %vector b
void extract (vector &a, const vector &b, long fi, long ncomp);

/// extracts only positive components from %vector a and stores them to %vector b
long extractposv (const vector &a, vector &b);
/// extracts only negative components from %vector a and stores them to %vector b
long extractnegv (const vector &a, vector &b);

/// Shell sort of ivector
void shell_sort(ivector &v);

/// Shell sort for long integer array
void shell_sort(long *v, long n);

/// function prints error message to standard error output device, format specifiers are allowed in emsg
void print_err(const char *emsg, const char *errfile, int errln, const char *errfunc, ...);

/// function prints warning message to standard error output device, format specifiers are allowed in wmsg
void print_warning(const char *wmsg, const char *warnfile, int warnln, const char *warnfunc, ...);
#endif
