#include "vector.h"
#define STCKGLOBVAR //global variables for the control of maximum stack size will be defined in this file
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>


#ifdef DEBUG_VECTOR
 static unsigned long Acv; ///< alocation counter for vector class
 static unsigned long Avmax; ///< peak of alocation memory of all vector instances
 static unsigned long Ava; ///< actual allocated memory of all vector instances

 static unsigned long Aciv; ///< alocation counter for ivector class
 static unsigned long Aivmax; ///< peak of alocation memory of all ivector instances
 static unsigned long Aiva; ///< actual allocated memory of all ivector instances

 static unsigned long Acbv; ///< alocation counter for bvector class
 static unsigned long Abvmax; ///< peak of alocation memory of all bvector instances
 static unsigned long Abva; ///< actual allocated memory of all bvector instances
#endif

/**
   The constructor allocates memory from the heap to the member array a.
   
   @param n is the number of %vector elements
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
vector :: vector (long n)
{
  vector::n = n;
  size = n;
  a = new double[size];
  if (a == NULL)
    print_err("cannot allocate memory for vector v(%ld)", __FILE__, __LINE__, __func__, n);
  stat.dealloc = 1;

  memset (a, 0, n*sizeof(*a));

  #ifdef DEBUG_VECTOR
   Acv++;
   Ava += size;
   if (Ava > Avmax)
     Avmax = Ava;
  #endif
}



/**
   The copy constructor creates copy of the object given by the constant reference parameter.
   
   @param v is reference to the object of the vector which should be copied
   
   created  16.10.2003 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
vector :: vector (const vector &v)
{
  n = v.n;
  if (n > size)
  {
    if (a && stat.dealloc)
      delete [] a;

    size = n;
    a = new double[size];
    if (a == NULL)
      print_err("cannot allocate memory for vector v(%ld)", __FILE__, __LINE__, __func__, n);
    stat.dealloc = 1;

    print_warning("Copy constructor of vector is called.\n"
                  "Please check your code and make sure you want use copy constructor."
                  "Consider usage of copyv function rather\n", __FILE__, __LINE__, __func__);

   #ifdef DEBUG_VECTOR
    Acv++;
    Ava += size;
    if (Ava > Avmax)
      Avmax = Ava;
   #endif
  }
  memcpy(a, v.a, n*sizeof(*v.a));
}



/** 
  The function reallocates the given %vector to the required number of components nc.

  @param[in] nc - the new number of %vector components.

  @retval 0 - on success
  @retval 1 - in the case of an error

  @date Created by TKo, 12.9.2024.
  @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz
*/
long vector::reallocv(long nc)
{
  if (size >= nc){
    n = nc;
    memset(a, 0, n*sizeof(*a));
  }
  else{
    if (a && stat.dealloc){
      delete [] a;
      #ifdef DEBUG_VECTOR
       Acv--;
       Ava -= size;
      #endif
    }

    size = n = nc;
    a = new double[size];
    if (a == NULL){
      print_err("cannot reallocate memory for vector v(%ld)", __FILE__, __LINE__, __func__, n);
      return (1);
    }
    stat.dealloc = 1;

    memset(a, 0, n*sizeof(*a));

   #ifdef DEBUG_VECTOR
    Acv++;
    Ava += size;
    if (Ava > Avmax)
      Avmax = Ava;
   #endif
  }
  return (0);
}



/** 
  The function reallocates the given %vector with the help of preallocated memory 
  to the required number of components nc.

  @param[in] nc - the new number of %vector components.
  @param[in] dealloc - flag for memory deallocation, i.e. it determines whether the memory 
                       referenced by ptr should be deallocated in the given %vector instance or 
                       the deallocation is managed outside the %vector class.
  @param[in] ptr - pointer to the preallocated memory array for %vector components.

  @retval 0 - on success
  @retval 1 - in the case of an error

  @date Created by TKo, 12.9.2024.
  @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz
*/
long vector::reallocv(long nc, unsigned dealloc, double *ptr)
{
  if (size >= nc)
  {
    n = nc;
    memset (a, 0, n*sizeof(*a));
    return 0;
  }
  else
  {
    if ((ptr == NULL) && n)
    {
      print_err("cannot reallocate memory for vector v(%ld)", __FILE__, __LINE__, __func__, n);
      return (1);
    }

    if (a && stat.dealloc)
    {
      delete [] a;
      #ifdef DEBUG_VECTOR
       Acv--;
       Ava -= size;
      #endif
    }
    
    size = n = nc;
    a = ptr;
    stat.dealloc = dealloc;

    memset (a, 0, n*sizeof(*a));

   #ifdef DEBUG_VECTOR
    Acv++;
    Ava += size;
    if (Ava > Avmax)
      Avmax = Ava;
   #endif
  }
  return (0);
}



/**
  The function makes the given %vector to referenece to the %vector src.

  @param[in] src - the source %vector to be referenced.

  @retval 0 - on success

  @date Created by TKo, 12.9.2024.
  @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz  
*/
long vector::makerefv(const vector &src)
{

  if (a && stat.dealloc)
  {
    delete [] a;
    #ifdef DEBUG_VECTOR
     Acv--;
     Ava -= size;
    #endif
  }

  size = src.size;
  n = src.n;
  a = src.a;
  stat.dealloc = 0;

  return (0);
}



/**
  The function makes the given %vector to referenece to the %vector src.

  @param[in] nc - the new number of %vector components.
  @param[in] ptr - pointer to the preallocated memory array for %vector components.

  @retval 0 - on success

  @date Created by TKo, 12.9.2024.
  @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz  
*/
long vector::makerefv(long nc, double *ptr)
{
  if (a && stat.dealloc){
    delete [] a;
    #ifdef DEBUG_VECTOR
     Acv--;
     Ava -= size;
    #endif
  }

  size = n = nc;
  a = ptr;
  stat.dealloc = 0;

  return (0);
}



/**
   The constructor allocates memory form the heap to the member array a.

   @param n is the number of ivector elements
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
ivector :: ivector (long n)
{
  ivector::n = n;
  size = n;
  a = new long[size];
  if (a == NULL)
    print_err("cannot allocate memory for ivector v(%ld)", __FILE__, __LINE__, __func__, n);
  stat.dealloc = 1;

  memset (a, 0, n*sizeof(*a));

 #ifdef DEBUG_VECTOR
  Aciv++;
  Aiva += size;
  if (Aiva > Aivmax)
    Aivmax = Aiva;
 #endif
}



/**
   The copy constructor creates copy of the object given by the constant reference parameter.
   
   @param v is reference to the object of the ivector which should be copied
   
   created  16.10.2003 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
ivector :: ivector (const ivector &v)
{
  n = v.n;
  if (n > size)
  {
    if (a && stat.dealloc)
      delete [] a;

    size = n;
    a = new long[size];
    if (a == NULL)
      print_err("cannot allocate memory for ivector v(%ld)", __FILE__, __LINE__, __func__, n);
    stat.dealloc = 1;

    print_warning("Copy constructor of ivector is called.\n"
                  "Please check your code and make sure you want use the copy constructor."
                  "Consider usage of copyv function rather\n", __FILE__, __LINE__, __func__);

   #ifdef DEBUG_VECTOR
    Aciv++;
    Aiva += size;
    if (Aiva > Aivmax)
      Aivmax = Aiva;
   #endif
  }
  memcpy(a, v.a, n*sizeof(*v.a));
}



/**
  The function makes the given %vector to referenece to the %vector src.

  @param[in] src - the source %vector to be referenced.

  @retval 0 - on success

  @date Created by TKo, 12.9.2024.
  @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz  
*/
long ivector::makerefv(const ivector &src)
{

  if (a && stat.dealloc)
  {
    delete [] a;
    #ifdef DEBUG_VECTOR
     Aciv--;
     Aiva -= ref.size;
    #endif
  }

  size = src.size;
  n = src.n;
  a = src.a;
  stat.dealloc = 0;

  return (0);
}



/**
  The function makes the given %vector to referenece to the %vector src.

  @param[in] nc - the new number of %vector components.
  @param[in] ptr - pointer to the preallocated memory array for %vector components.

  @retval 0 - on success

  @date Created by TKo, 12.9.2024.
  @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz  
*/
long ivector::makerefv(long nc, long *ptr)
{
  if (a && stat.dealloc){
    delete [] a;
    #ifdef DEBUG_VECTOR
     Aciv--;
     Aiva -= ref.size;
    #endif
  }

  size = n = nc;
  a = ptr;
  stat.dealloc = 0;

  return (0);
}




/** 
  The function reallocates an %ivector instance to the given size n and fills all components with a 
  sequence starting with the initial value init and incremented by incr 

  @param[in] nc - the number of components, which the %ivector will be reallocated to,
  @param[in] init - the initial value of the fi-th component,
  @param[in] incr - the increment of each successive component after the fi-th one.

  @date Created by TKo, 03.10.2024.
  @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz
*/
long ivector :: realloc(long nc, long init, long incr)
{
  long i, val;
  
  if (size >= nc)
  {
    n = nc;
    for (i=0, val=init; i<nc; i++){
      a[i] = val;
      val += incr;
    }
    return 0;
  } 
  else
  {
    if (a && stat.dealloc)
    {
      delete [] a;
      #ifdef DEBUG_VECTOR
       Aciv--;
       Aiva -= size;
      #endif
    }
   
    n = nc;
    size = nc;
    a = new long[size];
    if (a == NULL)
    {
      print_err("cannot reallocate memory for ivector v(%ld)", __FILE__, __LINE__, __func__, nc);
      return (1);
    }
    stat.dealloc = 1;

    for (i=0, val=init; i<nc; i++){
      a[i] = val;
      val += incr;
    }

    #ifdef DEBUG_VECTOR
     Aciv++;
     Aiva += size;
     if (Aiva > Aivmax)
       Aivmax = Aiva;
    #endif
  }
  return (0);  
}



/** 
  The function fills nc components of an %ivector instance starting from the fi-th component with a sequence 
  that starts with the initial value init incremented by incr.

  @param[in] fi - the index of the first component where the filling begins,
  @param[in] nc - the number of filled components,
  @param[in] init - the initial value of the fi-th component,
  @param[in] incr - the increment of each successive component after the fi-th one.

  @date Created by TKo, 03.10.2024.
  @author Tomas Koudelka (TKo), tomas.koudelka@fsv.cvut.cz
*/
void ivector::fillseq(long fi, long nc, long init, long incr)
{
  long i, val;
  
  for (i=fi, val=init; i<fi+nc && i<n; i++){
    a[i] = val;
    val += incr;
  }
}



/**
   The copy assignment operator copies the object given by the constant reference argumnet v to the given instance.
   
   @param v[in] is reference to the object of the %vector which should be copied to the given instance
   
   created  05.2023 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
ivector& ivector::operator=(const ivector &v)
{
  if (this != &v){ // check for self-assignmnet
    n = v.n;
    if (n > size){
      if (a && stat.dealloc)  delete [] a;
        
      size = n;
      a = new long[size];
      if (a == NULL)
        print_err("cannot allocate memory for vector v(%ld)", __FILE__, __LINE__, __func__, n);
      stat.dealloc = 1;
      
      print_warning("Copy assignment operator of ivector is called.\n"
                    "Please check your code and make sure you want use assignment operator."
                    "Consider usage of copyv function rather\n", __FILE__, __LINE__, __func__);

      #ifdef DEBUG_VECTOR
       Aciv++;
       Aiva += size;
       if (Aiva > Aivmax)
         Aivmax = Aiva;
      #endif
    }
    memcpy(a, v.a, n*sizeof(*v.a));
  }
  return *this;
}



/**
   The constructor allocates memory form the heap to the member array a.

   @param n is the number of bvector elements
   
   created  16.4.2024 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
bvector :: bvector (long n)
{
  bvector::n = n;
  size = n;
  a = new bool[size];
  if (a == NULL)
    print_err("cannot allocate memory for bvector v(%ld)", __FILE__, __LINE__, __func__, n);
  stat.dealloc = 1;

  memset (a, 0, n*sizeof(*a));

 #ifdef DEBUG_VECTOR
  Acbv++;
  Abva += size;
  if (Abva > Abvmax)
    Abvmax = Abva;
 #endif
}



/**
   The copy constructor creates copy of the object given by the constant reference parameter.
   
   @param v is reference to the object of the bvector which should be copied
   
   created  16.10.2003 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
bvector :: bvector (const bvector &v)
{
  n = v.n;
  if (n > size)
  {
    if (a && stat.dealloc)
      delete [] a;

    size = n;
    a = new bool[size];
    if (a == NULL)
      print_err("cannot allocate memory for bvector v(%ld)", __FILE__, __LINE__, __func__, n);
    stat.dealloc = 1;

    print_warning("Copy constructor of bvector is called.\n"
                  "Please check your code and make sure you want use the copy constructor."
                  "Consider usage of copyv function rather\n", __FILE__, __LINE__, __func__);

   #ifdef DEBUG_VECTOR
    Acbv++;
    Abva += size;
    if (Abva > Abvmax)
      Abvmax = Abva;
   #endif
  }
  memcpy(a, v.a, n*sizeof(*v.a));
}




/**
   The copy assignment operator copies the object given by the constant reference argumnet v to the given instance.
   
   @param v[in] is reference to the object of the %vector which should be copied to the given instance
   
   created  05.2023 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
bvector& bvector::operator=(const bvector &v)
{
  if (this != &v){ // check for self-assignmnet
    n = v.n;
    if (n > size){
      if (a && stat.dealloc)  delete [] a;
        
      size = n;
      a = new bool[size];
      if (a == NULL)
        print_err("cannot allocate memory for vector v(%ld)", __FILE__, __LINE__, __func__, n);
      stat.dealloc = 1;
      
      print_warning("Copy assignment operator of bvector is called.\n"
                    "Please check your code and make sure you want use assignment operator."
                    "Consider usage of copyv function rather\n", __FILE__, __LINE__, __func__);

      #ifdef DEBUG_VECTOR
       Acbv++;
       Abva += size;
       if (Abva > Abvmax)
         Abvmax = Abva;
      #endif
    }
    memcpy(a, v.a, n*sizeof(*v.a));
  }
  return *this;
}



#ifdef DEBUG_VECTOR
 long give_acv()
 {
   return Acv;
 }
 long give_avmax()
 {
   return Avmax;
 }
#endif



/**
   The function allocates memory from heap to the %vector's member array a.
   
   @param n   is the number of the vector elements
   @param vec is the structure for allocated vector
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long allocv(long n, vector &vec)
{
  vec.n = n;
  vec.size = n;  
  vec.a = new double[vec.size];

  if (vec.a)
  {
    vec.stat.dealloc = 1;

    memset(vec.a, 0, n*sizeof(*vec.a));

   #ifdef DEBUG_VECTOR
    Acv++;
    Ava += vec.size;
    if (Ava > Avmax)
      Avmax = Ava;
   #endif
    return (0);
  }
  else
  {
    print_err("cannot allocate memory for vector v(%ld)", __FILE__, __LINE__, __func__, n);
    return (1);
  }
  return 0;
}



/**
   The function reallocates memory for the vector's member array a.
   If the new vector size is greater than previously allocated one, 
   array a is deleted and new memory is allocted from the heap. 
   If the new vector size is less or equal, then the memory 
   allocated for the array a is left and the number of vector components 
   is changed. Components of the array a are set to zero.
   
   
   @param n   is the new number of the %vector components
   @param vec is the structure for reallocated %vector
   
   @retval 0 - on success
   @retval 1 - if fails allocating memory
   
   created  12.2012 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long reallocv(long n, vector &vec)
{
  if (vec.size >= n)
  {
    vec.n = n;
    memset (vec.a, 0, n*sizeof(*vec.a));
  }
  else
  {
    if (vec.a && vec.stat.dealloc)
    {
      delete [] vec.a;
      #ifdef DEBUG_VECTOR
       Acv--;
       Ava -= vec.size;
      #endif
    }

    vec.n = n;
    vec.size = n;
    vec.a = new double[vec.size];
    if (vec.a == NULL)
    {
      print_err("cannot reallocate memory for vector v(%ld)", __FILE__, __LINE__, __func__, n);
      return (1);
    }
    vec.stat.dealloc = 1;

    memset(vec.a, 0, n*sizeof(*vec.a));

   #ifdef DEBUG_VECTOR
    Acv++;
    Ava += vec.size;
    if (Ava > Avmax)
      Avmax = Ava;
   #endif
  }
  return (0);
}



/**
   The function reallocates memory for the vector's member array a.
   If the new vector size is greater than previously allocated one, 
   array a is deleted and new memory is allocated with help of pointer ptr which
   references to  preallocated memory for n element %vector. Additionally in this case,
   the dealloc flag specifies whether the ptr was allocated on the stack (dealloc=0) or 
   heap (dealloc=1).
   If the new vector size is less or equal, then the memory 
   allocated for the array a is left and the number of vector components 
   is changed. In such case, arguments ptr and dealloc are not used 
   and NULL or zero values should be passed.
   Components of the array a are set to zero.
   
   @param n   is the new number of the %vector components
   @param vec is the structure for reallocated vector
   @param dealloc is the flag for deallocation control (0=no deallocation ,1=release memory by delete operator)
   @param ptr is the pointer to the allocated memory which will be used for storage of n %vector elements
   
   @retval 0 - on success
   @retval 1 - if fails allocating memory
   
   created  5.2015 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long reallocv(long n, vector &vec, unsigned dealloc, double *ptr)
{
  if (vec.size >= n)
  {
    vec.n = n;
    memset (vec.a, 0, n*sizeof(*vec.a));
    return 0;
  }
  else
  {
    if ((ptr == NULL) && n)
    {
      print_err("cannot reallocate memory for vector v(%ld)", __FILE__, __LINE__, __func__, n);
      return (1);
    }

    if (vec.a && vec.stat.dealloc)
    {
      delete [] vec.a;
      #ifdef DEBUG_VECTOR
       Acv--;
       Ava -= vec.size;
      #endif
    }
    
    vec.n = n;
    vec.size = n;
    vec.a = ptr;
    vec.stat.dealloc = dealloc;

    memset (vec.a, 0, n*sizeof(*vec.a));

   #ifdef DEBUG_VECTOR
    Acv++;
    Ava += vec.size;
    if (Ava > Avmax)
      Avmax = Ava;
   #endif
  }
  return (0);
}



/**
   The function allocates memory from heap to the ivector's member array a.
   
   @param n   is the number of the ivector elements
   @param vec is the structure for allocated ivector
   
   @retval 0 - on success
   @retval 1 - if fails allocating memory
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long allocv(long n, ivector &vec)
{
  vec.n = n;
  vec.size = n;
  vec.a = new long[vec.size];
  if (vec.a)
  {
    vec.stat.dealloc = 1;

    memset (vec.a, 0, n*sizeof(*vec.a));

   #ifdef DEBUG_VECTOR
    Aciv++;
    Aiva += vec.size;
    if (Aiva > Aivmax)
      Aivmax = Aiva;
   #endif
    return 0;
  }
  else
  {
    print_err("cannot allocate memory for ivector v(%ld)", __FILE__, __LINE__, __func__, n);
    return (1);
  }

  return (0);
}



/**
   The function initializes pointer for the vector's member array a by the parameter
   ptr and thus creates reference to the ptr. Deallocation flag of ref %vector is set to 0.
   Component array a of the %vector ref is deallocated if the stat.dealloc flag is set on.
   
   @param ref is the structure of initialized %vector by ptr pointer
   @param ptr is the pointer to the allocated array of n components which will be referenced by ref %vector.
   @param n   is the number of the ref/ptr components
   
   @retval 0 - on success
   
   created  5.2016 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long makerefv(vector &ref, double *ptr, long n)
{

  if (ref.a && ref.stat.dealloc)
  {
    delete [] ref.a;
    #ifdef DEBUG_VECTOR
     Acv--;
     Ava -= ref.size;
    #endif
  }

  ref.size = ref.n = n;
  ref.a = ptr;
  ref.stat.dealloc = 0;

  return (0);
}



/**
   The function creates reference %vector ref to the source %vector src.
   The ref %vector is initialized by the content of src %vector but only SHALOW
   copy of data members is performed except of deallocation flag of ref %vector which is set to 0.
   Component array a of the %vector ref is deallocated if the stat.dealloc flag is set on.
   
   @param ref is the structure for %vector used for storage of new refrence to src
   @param src structure with the source %vector which will be referenced by ref 
   
   @retval 0 - on success
   
   created  5.2015 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long makerefv(vector &ref, vector &src)
{

  if (ref.a && ref.stat.dealloc)
  {
    delete [] ref.a;
    #ifdef DEBUG_VECTOR
     Acv--;
     Ava -= ref.size;
    #endif
  }

  ref.size = src.size;
  ref.n = src.n;
  ref.a = src.a;
  ref.stat.dealloc = 0;

  return (0);
}



/**
   The function initializes pointer for the ivector's member array a by the parameter
   ptr and thus creates reference to the ptr. Deallocation flag of ref %vector is set to 0.
   Component array a of the %ivector ref is deallocated if the stat.dealloc flag is set on.
   
   @param ref is the structure of initialized %ivector by ptr pointer
   @param ptr is the pointer to the allocated array of n components which will be referenced by ref %ivector.
   @param n   is the number of the ref/ptr components
   
   @retval 0 - on success
   
   created  1.2020 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long makerefv(ivector &ref, long *ptr, long n)
{

  if (ref.a && ref.stat.dealloc)
  {
    delete [] ref.a;
    #ifdef DEBUG_VECTOR
     Aciv--;
     Aiva -= ref.size;
    #endif
  }

  ref.size = ref.n = n;
  ref.a = ptr;
  ref.stat.dealloc = 0;

  return (0);
}



/**
   The function creates reference %ivector ref to the source %ivector src.
   The ref %ivector is initialized by the content of src %ivector but only SHALOW
   copy of data members is performed except of deallocation flag of ref %ivector which is set to 0.
   Component array a of the %ivector ref is deallocated if the stat.dealloc flag is set on.
   
   @param ref is the structure for %ivector used for storage of new refrence to src
   @param src structure with the source %ivector which will be referenced by ref 
   
   @retval 0 - on success
   
   created  1.2020 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long makerefv(ivector &ref, ivector &src)
{

  if (ref.a && ref.stat.dealloc)
  {
    delete [] ref.a;
    #ifdef DEBUG_VECTOR
     Aciv--;
     Aiva -= ref.size;
    #endif
  }

  ref.size = src.size;
  ref.n = src.n;
  ref.a = src.a;
  ref.stat.dealloc = 0;

  return (0);
}



/**
   The function reallocates memory for the ivector's member array a.
   If the new ivector size is greater than previously allocated one, 
   array a is deleted and new memory is allocted from the heap. 
   If the new ivector size is less or equal, then the memory 
   allocated for the array a is left and the number of %ivector components 
   is changed. Components of the array a are set to zero.
   
   
   @param n   is the new number of the ivector components
   @param vec is the structure for reallocated ivector
   
   @retval 0 - on success
   @retval 1 - if fails allocating memory
   
   created  12.2012 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long reallocv(long n, ivector &vec)
{
  if (vec.size >= n)
  {
    vec.n = n;
    memset (vec.a, 0, n*sizeof(*vec.a));
    return 0;
  } 
  else
  {
    if (vec.a && vec.stat.dealloc)
    {
      delete [] vec.a;
      #ifdef DEBUG_VECTOR
       Aciv--;
       Aiva -= vec.size;
      #endif
    }
   
    vec.n = n;
    vec.size = n;
    vec.a = new long[vec.size];
    if (vec.a == NULL)
    {
      print_err("cannot reallocate memory for ivector v(%ld)", __FILE__, __LINE__, __func__, n);
      return (1);
    }
    vec.stat.dealloc = 1;

    memset (vec.a, 0, n*sizeof(*vec.a));

    #ifdef DEBUG_VECTOR
     Aciv++;
     Aiva += vec.size;
     if (Aiva > Aivmax)
       Aivmax = Aiva;
    #endif
  }
  return (0);
}



/**
   The function reallocates memory for the ivector's member array a.
   If the new ivector size is greater than previously allocated one, 
   array a is deleted and new memory is allocated with help of pointer ptr which
   references to  preallocated memory for n element %ivector. Additionally in this case,
   the dealloc flag specifies whether the ptr was allocated on the stack (dealloc=0) or 
   heap (dealloc=1).
   If the new %ivector size is less or equal to the actual one, then the memory 
   allocated for the array a is left and the number of vector components 
   is changed. In such case, arguments ptr and dealloc are not used and NULL or zero values 
   should be passed. Components of the array a are set to zero.
   
   @param n   is the new number of the ivector components
   @param vec is the structure for reallocated %vector
   @param dealloc is the flag for deallocation control (0=no deallocation ,1=release memory by delete operator)
   @param ptr is the pointer to the allocated memory which will be used for storage of n %vector elements
   
   @retval 0 - on success
   @retval 1 - if fails allocating memory
   
   created  5.2015 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long reallocv(long n, ivector &vec, unsigned dealloc, long *ptr)
{
  if (vec.size >= n)
  {
    vec.n = n;
    memset (vec.a, 0, n*sizeof(*vec.a));
    return 0;
  }
  else
  {
    if ((ptr == NULL) && n)
    {
      print_err("cannot reallocate memory for ivector v(%ld)", __FILE__, __LINE__, __func__, n);
      return (1);
    }

    if (vec.a && vec.stat.dealloc)
    {
      delete [] vec.a;
      #ifdef DEBUG_VECTOR
       Aciv--;
       Aiva -= vec.size;
      #endif
    }
    
    vec.n = n;
    vec.size = n;
    vec.a = ptr;
    vec.stat.dealloc = dealloc;

    memset (vec.a, 0, n*sizeof(*vec.a));

   #ifdef DEBUG_VECTOR
    Aciv++;
    Aiva += vec.size;
    if (Aiva > Aivmax)
      Aivmax = Aiva;
   #endif
  }
  return (0);
}



/**
   The function reallocates memory for the bvector's member array a.
   If the new bvector size is greater than previously allocated one, 
   array a is deleted and new memory is allocted from the heap. 
   If the new bvector size is less or equal, then the memory 
   allocated for the array a is left and the number of %bvector components 
   is changed. Components of the array a are set to zero.
   
   
   @param n   is the new number of the bvector components
   @param vec is the structure for reallocated bvector
   
   @retval 0 - on success
   @retval 1 - if fails allocating memory
   
   created  12.2012 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long reallocv(long n, bvector &vec)
{
  if (vec.size >= n)
  {
    vec.n = n;
    memset (vec.a, 0, n*sizeof(*vec.a));
    return 0;
  } 
  else
  {
    if (vec.a && vec.stat.dealloc)
    {
      delete [] vec.a;
      #ifdef DEBUG_VECTOR
       Acbv--;
       Abva -= vec.size;
      #endif
    }
   
    vec.n = n;
    vec.size = n;
    vec.a = new bool[vec.size];
    if (vec.a == NULL)
    {
      print_err("cannot reallocate memory for bvector v(%ld)", __FILE__, __LINE__, __func__, n);
      return (1);
    }
    vec.stat.dealloc = 1;

    memset (vec.a, 0, n*sizeof(*vec.a));

    #ifdef DEBUG_VECTOR
     Acbv++;
     Abva += vec.size;
     if (Abva > Abvmax)
       Abvmax = Abva;
    #endif
  }
  return (0);
}



/**
   The function reallocates memory for the bvector's member array a.
   If the new bvector size is greater than previously allocated one, 
   array a is deleted and new memory is allocated with help of pointer ptr which
   references to  preallocated memory for n element %bvector. Additionally in this case,
   the dealloc flag specifies whether the ptr was allocated on the stack (dealloc=0) or 
   heap (dealloc=1).
   If the new %bvector size is less or equal to the actual one, then the memory 
   allocated for the array a is left and the number of vector components 
   is changed. In such case, arguments ptr and dealloc are not used and NULL or zero values 
   should be passed. Components of the array a are set to zero.
   
   @param n   is the new number of the bvector components
   @param vec is the structure for reallocated %vector
   @param dealloc is the flag for deallocation control (0=no deallocation ,1=release memory by delete operator)
   @param ptr is the pointer to the allocated memory which will be used for storage of n %vector elements
   
   @retval 0 - on success
   @retval 1 - if fails allocating memory
   
   created  5.2015 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long reallocv(long n, bvector &vec, unsigned dealloc, bool *ptr)
{
  if (vec.size >= n)
  {
    vec.n = n;
    memset (vec.a, 0, n*sizeof(*vec.a));
    return 0;
  }
  else
  {
    if ((ptr == NULL) && n)
    {
      print_err("cannot reallocate memory for bvector v(%ld)", __FILE__, __LINE__, __func__, n);
      return (1);
    }

    if (vec.a && vec.stat.dealloc)
    {
      delete [] vec.a;
      #ifdef DEBUG_VECTOR
       Acbv--;
       Abva -= vec.size;
      #endif
    }
    
    vec.n = n;
    vec.size = n;
    vec.a = ptr;
    vec.stat.dealloc = dealloc;

    memset (vec.a, 0, n*sizeof(*vec.a));

   #ifdef DEBUG_VECTOR
    Acbv++;
    Abva += vec.size;
    if (Abva > Abvmax)
      Abvmax = Abva;
   #endif
  }
  return (0);
}



/**
   The function copies vector given by src to dest.
   
   @param src  is the structure of source vector to copy
   @param dest is the structure of destination vector to which will be copied contents of src
   
   @b Requests :
   dest has to be setuped dimensions and allocated array for elements
   which is enough large to hold all contents of the vector src.

   @retval 0 : on succes
   @retval 1 : in case incompatibility sizes of the src and dest vectors
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long copyv(const vector &src, vector &dest)
{
  if (src.n != dest.n)
  {
    print_err("cannot copy vectors - incompatible dimensions of vectors\n"
              "src(%ld) X dest(%ld)", __FILE__, __LINE__, __func__, src.n, dest.n);
    return (1);
  }
  memcpy(dest.a, src.a, sizeof(*src.a)*src.n);
  return (0);
}



/**
   The function copies %vector given by src to dest.
   
   @param src  is the structure of source %ivector to copy
   @param dest is the structure of destination %ivector to which will be copied contents of src
   
   @b Requests :
   dest has to be setuped dimensions and allocated array for elements
   which is enough large to hold all contents of the %ivector src.
   
   @retval 0 on succes
   @retval 1 in case incompatibility sizes of the src and dest vectors
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long copyv(const ivector &src, ivector &dest)
{
  if (src.n != dest.n)
  {
    print_err("cannot copy ivector - incompatible dimensions of ivectors\n"
              "src(%ld) X dest(%ld)", __FILE__, __LINE__, __func__, src.n, dest.n);
    return (1);
  }
  memcpy(dest.a, src.a, sizeof(*src.a)*src.n);
  return (0);
}



/**
   The function copies %vector given by src to dest.
   
   @param src  is the structure of source %bvector to copy
   @param dest is the structure of destination %bvector to which will be copied contents of src
   
   @b Requests :
   dest has to be setuped dimensions and allocated array for elements
   which is enough large to hold all contents of the %bvector src.
   
   @retval 0 on succes
   @retval 1 in case incompatibility sizes of the src and dest vectors
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long copyv(const bvector &src, bvector &dest)
{
  if (src.n != dest.n)
  {
    print_err("cannot copy bvector - incompatible dimensions of bvectors\n"
              "src(%ld) X dest(%ld)", __FILE__, __LINE__, __func__, src.n, dest.n);
    return (1);
  }
  memcpy(dest.a, src.a, sizeof(*src.a)*src.n);
  return (0);
}



/**
   The function copies long array given by src to dest.
   
   @param src  is the double array to copy
   @param dest is the double array to which will be copied contents of src
   
   @retval always zero
   
   created  3.11.2001 by JK
*/
long copyv(const long *src, long *dest, long n)
{
  memcpy(dest, src, sizeof(*src)*n);
  return (0);
}



/**
   The function copies double array given by src to dest.
   
   @param src  is the double array to copy
   @param dest is the double array to which will be copied contents of src
   
   @retval always zero
   
   created  3.11.2001 by JK
*/
long copyv (const double *src, double *dest, long n)
{
  memcpy(dest, src, sizeof(*src)*n);
  return (0);
}



/**
   The function copies long array given by src to dest.
   
   @param src  is the double array to copy
   @param dest is the double array to which will be copied contents of src
   
   @retval always zero
   
   created  3.11.2001 by JK
*/
long copyv(const bool *src, bool *dest, long n)
{
  memcpy(dest, src, sizeof(*src)*n);
  return (0);
}



/**
   The function copies array given by src to %vector dest.
   
   @param src  is source array to copy
   @param dest is the structure of destination vector which will be copied contents of src to
   
   @b Requests :
   dest has to be setuped dimensions and allocated array for required elements.
   Number of copied elements is given by the dest dimension

   @retval 0 : on succes
   
   created  4.10.2007 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
 */
long copyv(const double *src, vector &dest)
{
  memcpy(dest.a, src, sizeof(*dest.a)*dest.n);
  return (0);
}



/**
   The function copies integer array given by src to %vector dest.
   
   @param src  is source array to copy
   @param dest is the structure of destination vector which will be copied contents of src to
   
   @b Requests :
   dest has to be setuped dimensions and allocated array for required elements.
   Number of copied elements is given by the dest dimension

   @retval 0 : on succes
   
   created  4.10.2007 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
 */
long copyv(const long *src, ivector &dest)
{
  memcpy(dest.a, src, sizeof(*dest.a)*dest.n);
  return (0);
}



/**
   The function copies integer array given by src to %vector dest.
   
   @param src  is source array to copy
   @param dest is the structure of destination vector which will be copied contents of src to
   
   @b Requests :
   dest has to be setuped dimensions and allocated array for required elements.
   Number of copied elements is given by the dest dimension

   @retval 0 : on succes
   
   created  4.10.2007 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
 */
long copyv(const bool *src, bvector &dest)
{
  memcpy(dest.a, src, sizeof(*dest.a)*dest.n);
  return (0);
}



/**
   The function copies vector given by src to array dest.
   
   @param src  is the structure of source vector to copy
   @param dest is the destination array which will be copied contents of src to
   
   @b Requests :
   dest has to have allocated memory which is enough large to hold all content
   of the vector src. Number of copied elements is given by the src dimension

   @retval 0 : on succes
   
   created  4.10.2007 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
 */
long copyv(const vector &src, double *dest)
{
  memcpy(dest, src.a, sizeof(*src.a)*src.n);
  return (0);
}



/**
   The function copies ivector given by src to array dest.
   
   @param src  is the structure of source ivector to copy
   @param dest is the destination array which will be copied contents of src to
   
   @b Requests :
   dest has to have allocated memory which is enough large to hold all content
   of the ivector src. Number of copied elements is given by the src dimension

   @retval 0 : on succes
   
   created  21.05.2021 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
 */
long copyv(const ivector &src, long *dest)
{
  memcpy(dest, src.a, sizeof(*src.a)*src.n);
  return (0);
}



/**
   The function copies bvector given by src to array dest.
   
   @param src  is the structure of source bvector to copy
   @param dest is the destination array which will be copied contents of src to
   
   @b Requests :
   dest has to have allocated memory which is enough large to hold all content
   of the bvector src. Number of copied elements is given by the src dimension

   @retval 0 : on succes
   
   created  21.05.2021 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
 */
long copyv(const bvector &src, bool *dest)
{
  memcpy(dest, src.a, sizeof(*src.a)*src.n);
  return (0);
}



/**
  The function swaps pointers to component arrays of two vectors, i.e. it swaps 
  their content. In case of vectors allocated on the stack, the user must be careful
  because the function also changes lifetime of pointers with respect to the stack frame 
  in which they were allocated originally. Generally, two vectors whose component
  arrays were stack allocated can be swapped safely only if they were allocated within 
  the same stack frame (i.e function or satement block).

  @param a - the first %vector
  @param b - the second %vector

  @retval 0 - on success
  @retval 1 - in the case of incompatible dimensions of vectors

  Created by Tomas Koudelka, 9.10.2015
*/
long swapv(vector &a, vector &b)
{
  double *tmp;
  long aux;
  unsigned auxd;

  if (a.n != b.n)
  {
    print_err("cannot swap vectors - incompatible dimensions of vectors\n"
              "a(%ld) X b(%ld)", __FILE__, __LINE__, __func__, a.n, b.n);
    return (1);
  }
  tmp = a.a;
  a.a = b.a;
  b.a = tmp;  
  aux = a.size;
  a.size = b.size;
  b.size = aux;
  auxd = a.stat.dealloc;
  a.stat.dealloc = b.stat.dealloc;
  b.stat.dealloc = auxd;

  return 0;
}



/**
  The function copies range of components from vector %src to %vector dest.
  
  @param src - source %vector
  @param src_fi  - the first index of the component to be copied in the source %vector (in C++ notation)
  @param dest  - destination %vector
  @param dest_fi  - the first index from which the copied components will be stored in the destination %vector (in C++ notation)
  @param n - number of components to be copied

  @retval 0 : on succes
  @retval 1 : in case incompatibility sizes of the src and dest vectors
   
  Created  21.4.2015 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long rcopyv(const vector &src, long src_fi, vector &dest, long dest_fi, long n)
{
  if ((src_fi < 0) || (dest_fi < 0))
  {
    print_err("cannot copy vectors - range indeces are negative\n" 
              "src_fi=%ld, dest_fi=%ld", __FILE__, __LINE__, __func__, src_fi, dest_fi);
    return (1);
  }
  if (((src_fi + n) > src.n) || ((dest_fi + n) > dest.n))
  {
    print_err("cannot copy vectors - incompatible dimensions of vectors and required index range\n" 
              "src(%ld), dest(%ld), src_fi=%ld, dest_fi=%ld, n=%ld", __FILE__, __LINE__, __func__, src.n, dest.n, src_fi, dest_fi, n);
    return (1);
  }
  memcpy(dest.a+dest_fi, src.a+src_fi, sizeof(*src.a)*n);
  return 0;
}



/**
  The function copies range of components from vector %src to %vector dest. Vectors are given by pointers 
  to arrays of components.
  
  @param src - pointer to the array of source %vector components
  @param src_fi  - the first index of the component to be copied in the source %vector (in C++ notation)
  @param src_n  - the total number of components of the source %vector
  @param dest  - pointer to the array of destination %vector components
  @param dest_fi  - the first index from which the copied components will be stored in the destination %vector (in C++ notation)
  @param dest_n  - the total number of components of the destination %vector
  @param n - number of components to be copied

  @retval 0 : on succes
  @retval 1 : in case incompatibility sizes of the src and dest vectors
   
  Created  04.2016 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long rcopyv(const double *src, long src_fi, long src_n, double *dest, long dest_fi, long dest_n, long n)
{
  if ((src_fi < 0) || (dest_fi < 0))
  {
    print_err("cannot copy vectors - range indeces are negative\n" 
              "src_fi=%ld, dest_fi=%ld", __FILE__, __LINE__, __func__, src_fi, dest_fi);
    return (1);
  }
  if (((src_fi + n) > src_n) || ((dest_fi + n) > dest_n))
  {
    print_err("cannot copy vectors - incompatible dimensions of vectors and required index range\n" 
              "src(%ld), dest(%ld), src_fi=%ld, dest_fi=%ld, n=%ld", __FILE__, __LINE__, __func__, src_n, dest_n, src_fi, dest_fi, n);
    return (1);
  }
  memcpy(dest+dest_fi, src+src_fi, sizeof(*src)*n);
  return 0;
}



/**
   The function copies vector given by src multiplied by c to dest.
   
   @param src  - structure with source %vector to copy
   @param dest - structure with destination %vector to which will be copied contents of the src
   @param c    - scalar multiplier
   
   @b Requests :
   dest has to have dimensions and allocated array for %vector components
   which is large enough to hold all contents of the %vector src.

   @retval 0 : on succes
   @retval 1 : in case incompatibility sizes of the src and dest vectors
   
   created  08.2011 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long copymultv(const vector &src, vector &dest, double c)
{
  if (src.n != dest.n)
  {
    print_err("cannot copy vectors - incompatible dimensions of vectors\n"
              "src(%ld) X dest(%ld)", __FILE__, __LINE__, __func__, src.n, dest.n);
    return (1);
  }
  long i;  
  for (i=0; i<src.n; i++)
    dest[i] = c*src[i];
  return (0);
}



/**
   The function copies ivector given by src multiplied by scalar c to dest.
   
   @param src  - structure with the source ivector to copy
   @param dest - structure with the destination %ivector to which will be copied the contents of the src
   @param c    - scalar multiplier
   
   @b Requests :
   dest has to have dimensions and allocated array for %vector components
   which is large enough to hold all contents of the %ivector src.
   
   @retval 0 on succes
   @retval 1 in case incompatibility sizes of the src and dest vectors
   
   created  08.2011 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long copymultv(const ivector &src, ivector &dest, long c)
{
  if (src.n != dest.n)
  {
    print_err("cannot copy ivector - incompatible dimensions of ivectors\n" 
              "src(%ld) X dest(%ld)", __FILE__, __LINE__, __func__, src.n, dest.n);
    return (1);
  }
  long i;  
  for (i=0; i<src.n; i++)
    dest[i] = c*src[i];
  return (0);
}



/**
   The function copies double array given by the src multiplied by scalar c to the dest.
   
   @param[in]  src  - double array to copy,
   @param[out] dest - double array to which will be copied contents of src,
   @param[in]  c    - scalar multiplier,
   @param[in]  n    - length of arrays src and dest.
   
   @retval always zero
   
   created  08.2011 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long copymultv(const double *src, double *dest, double c, long n)
{
  long i;  
  for (i=0; i<n; i++)
    dest[i] = c*src[i];
  return (0);
}



/**
   The function fills vec's member array a with value c.
   
   @param[in]     c   - a value used for the %vector filling,
   @param[in,out] vec - a %vector to be filled.
   
   @return always zero
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long fillv(double c, vector &vec)
{
  for (long i = 0; i < vec.n; vec(i) = c, i++);
  return (0);
}



/**
   The function fills vec's member array a with value c.
   
   @param[in]     c   - a value used for the array filling,
   @param[in,out] vec - %ivector to be filled.
   
   @retval always zero
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long fillv(long c, ivector &vec)
{
  for (long i = 0; i < vec.n; vec(i) = c, i++);
  return (0);
}



/**
   The function fills vec's member array a with value c
   
   @param[in]     c   - a value used for the filling vec,
   @param[in,out] vec - %bvector to be filled.
   
   @retval always zero
   
   created  16.4.2024 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long fillv(bool c, bvector &vec)
{
  for (long i = 0; i < vec.n; vec(i) = c, i++);
  return (0);
}



/**
   The function fills array vec with value c
   
   @param[in] c       - a value, which will be used for filling memory,
   @param[in,out] vec - array to be filled,
   @param[in] n       - dimension of array vec.
   
   @retval always zero
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void fillv(long c, long *vec, long n)
{
  long i;
  
  for (i=0;i<n;i++)
    vec[i]=c;
}



/**
   The function fills vec's member array a with value c
   
   @param[in]     c   - value, which will be used for filling memory
   @param[in,out] vec - the filled %ivector
   
   @retval always zero
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void fillv(double c, double *vec, long n)
{
  long i;
  
  for (i=0;i<n;i++)
    vec[i]=c;
}



/**
   The function fills double array a with zeros
   
   @param[in,out] a is array which should be filled
   
   @retval always zero
*/
long nullv (double *a, long n)
{
  memset (a,0,n*sizeof(*a));
  return (0);
}



/**
   The function fills long integer array a with zeros
   
   @param[in,out] a is array which should be filled
   
   @retval always zero
*/
long nullv (long *a, long n)
{
  memset (a,0,n*sizeof(*a));
  return (0);
}



/**
   The function fills components of %vector a with zeros.
   
   @param[in,out] a is the vector structure which should be filled
   
   @retval always zero
*/
long nullv (vector &a)
{
  memset (a.a,0,a.n*sizeof(*a.a));
  return (0);
}



/**
   The function fills components of %ivector a with zeros.
   
   @param[in,out] a is the ivector structure which should be filled
   
   @retval always zero
*/
long nullv (ivector &a)
{
  memset (a.a,0,a.n*sizeof(*a.a));
  return (0);
}



/**
   The function fills components of %ivector a with false value.
   
   @param[in,out] a is the bvector structure which should be filled
   
   @retval always zero
*/
long nullv (bvector &a)
{
  memset (a.a,0,a.n*sizeof(*a.a));
  return (0);
}



/**
  The function changes sign of all components of the %vector v.

  @param[in,out] v - is the structure with vector whose components sign should be changed

  @return The function does not return anything, it changes the content of v.

  Created by Tomas Koudelka, 11.2015
*/
void chsgnv(vector &v)
{
  long i;
  for (i=0; i<v.n; i++)
    v(i) = -v(i);
}



/**
  The function changes sign of all components of the %vector v.

  @param[in,out] a - is the array whose components sign should be changed
  @param[in]     n - the number of array components

  @return The function does not return anything, it changes the content of a.

  Created by Tomas Koudelka, 11.2015
*/
void chsgnv(double *a, long n)
{
  long i;
  for (i=0; i<n; i++)
    a[i] = -a[i];
}



/**
   The function deallocates memory occupied by %vector vec
   
   @param[in] vec is the structure for the deallocated %vector
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long destrv(vector &vec)
{
  if (vec.a && vec.stat.dealloc)
  {
    delete [] vec.a;
    #ifdef DEBUG_VECTOR
     Acv--;
     Ava -= vec.size;
    #endif
  }

  vec.n = 0;
  vec.size = 0;
  vec.a = NULL;
  vec.stat.dealloc = 0;

  return (0);
}



/**
   The function deallocates memory occupied by %ivector vec
   
   @param[in] vec is the structure for the deallocated %ivector
   
   @retval 0
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long destrv(ivector &vec)
{
  if (vec.a && vec.stat.dealloc)
  {
    delete [] vec.a;
    #ifdef DEBUG_VECTOR
     Aciv--;
     Aiva -= vec.size;
    #endif
  }

  vec.n = 0;
  vec.size = 0;
  vec.a = NULL;
  vec.stat.dealloc = 0;

  return (0);
}



/**
   The function adds vector given by a to %vector given by b, the result is stored in c
   
   @param[in]  a - the first summand %vector
   @param[in]  b - the second summand %vector
   @param[out] c - the resulting %vector
   
   @b Requests :
   a, b and c have to be same dimensions, c has to have allocated array for elements
   which is enough large to hold contents of the result.
   
   @retval 0 : on succes
   @retval 1 : in case incompatibility sizes of a, b and c vectors
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long addv(const vector &a, const vector &b, vector &c)
{
  if ((a.n != b.n) || (b.n != c.n))
  {
    print_err("cannot add vectors due to their incompatible dimensions\n"
              "a(%ld) X b(%ld) X c(%ld)", __FILE__, __LINE__, __func__, a.n, b.n, c.n);
    return (1);
  }
  for (long i = 0; i < a.n; i++)
    c.a[i] = a.a[i] + b.a[i];

  return(0);
}



/**
   The function adds %vector given by argument b to the vector given by argument a, i.e. a += b.
   
   @param[in,out] a is the structure of the first vector which will be added to.
   @param[in]     b is the structure of the second added vector.
   
   @b Requests :
   a and b have to be same dimensions
   
   @retval 0 : on succes
   @retval 1 : in case incompatibility sizes of a and b.   
   Created  by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz, 09.2023
*/
long addv(vector &a, const vector &b)
{
  if (a.n != b.n)
  {
    print_err("cannot add vectors due to their incompatible dimensions\n"
              "a(%ld) X b(%ld)", __FILE__, __LINE__, __func__, a.n, b.n);
    return (1);
  }
  for (long i = 0; i < a.n; i++)
    a.a[i] += b.a[i];

  return(0);
}



/**
   The function adds vector given by a to %vector given by b, the result is stored in c.
   
   @param[in] a is the structure of the first added vector
   @param[in] b is the structure of the second added vector
   @param[out] c is the structure of the result vector
   
   @b Requests :
   a, b and c must be the same dimensions, c must have allocated array for elements
   which is large enough to hold contents of the result.
   
   @return The function returns the result in the parameter c.

   created  08.2011 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void addmultv (const vector &a, const vector &b, double bc, vector &c)
{
  if ((a.n != b.n) || (b.n != c.n))
  {
    print_err("cannot add vectors due to their incompatible dimensions\n"
              "a(%ld) X b(%ld) X c(%ld)", __FILE__, __LINE__, __func__, a.n, b.n, c.n);
    return;
  }
  for (long i=0;i<a.n;i++)
    c[i] = a[i]+b[i]*bc;
}



/**
   The function adds vector given by a to %vector given by b which is multiplied by bc, the result is stored in a.
   
   @param[in,out] a is the structure of the first added vector
   @param[in]     b is the structure of the second added vector
   
   @b Requests :
   a and b must have the same dimensions, c must have allocated array for elements
   which is large enough to hold contents of the result.
   
   @return The function returns the result in the parameter a.

   created  08.2011 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void addmultv (vector &a, const vector &b, double bc)
{
  if (a.n != b.n)
  {
    print_err("cannot add vectors due to their incompatible dimensions\n"
              "a(%ld) X b(%ld)", __FILE__, __LINE__, __func__, a.n, b.n);
    return;
  }
  for (long i=0;i<a.n;i++)
    a[i] += b[i]*bc;
}



/**
  The function adds 2 vectors a and b, that are multiplied by the independent constants ac and bc and 
  result is stored in the %vector c:
  c_i = ac*a_i + bc*b_i
  
  @param[in]  a  - the first %vector added
  @param[in]  ac - constant for the scaling of %vector a
  @param[in]  b  - the second %vector added
  @param[in]  bc - constant for the scaling of %vector b
  @param[out] c  - resulting %vector (output)

  @return The function stores result of ac*a_i + bc*b_i in the argument c.

  Created by Tomas Koudelka, 4.12.2015
*/
void addmultv (const vector &a, double ac, const vector &b, double bc, vector &c)
{
  long i;

  if ((a.n == b.n) && (b.n == c.n))
  {
    for(i=0; i<c.n; i++)
      c(i) = ac*a(i) + bc*b(i);
  }
  else
  {
    print_err("incompatible dimensions of vectors: a(%ld), b(%ld), c(%ld)", __FILE__, __LINE__, __func__, a.n, b.n, c.n);
    abort();
  }
}



/**
  The function adds 2 arrays a and b, that are multiplied by the independent constants ac and bc and 
  result is stored in the array a:
  a[i] = ac*a[i] + bc*b[i]
  All arrays must be allocated and they must have length at least n.
  
  @param[in,out] a  - pointer to the first array added (input/output)
  @param[in] ac - constant for the scaling of array a
  @param[in] b  - pointer to the second array added
  @param[in] bc - constant for the scaling of array b
  @param[in] n  - number of components of arrays a and b

  @return The function stores result of ac*a[i] + bc*b[i] in the argument a.

  Created by Tomas Koudelka, 4.12.2015
*/
void addmultv (double *a, double ac, const double *b, double bc, long n)
{
  long i;

  for(i=0; i<n; i++)
    a[i] = ac*a[i] + bc*b[i];
}



/**
  The function adds 2 arrays a and b, that are multiplied by the independent constants ac and bc and 
  result is stored in the array c:
  c[i] = ac*a[i] + bc*b[i]
  All arrays must be allocated and they must have length at least n.
  
  @param a  - pointer to the first array added
  @param ac - constant for the scaling of array a
  @param b  - pointer to the second array added
  @param bc - constant for the scaling of array b
  @param c  - pointer to the resulting array (output)
  @param n  - number of components of arrays a, b and c

  @return The function stores result of ac*a[i] + bc*b[i] in the argument c.

  Created by Tomas Koudelka, 4.12.2015
*/
void addmultv (const double *a, double ac, const double *b, double bc, double *c, long n)
{
  long i;

  for(i=0; i<n; i++)
    c[i] = ac*a[i] + bc*b[i];
}



/**
   The function adds ivector given by a to ivector given by b, the result is stored in c
   
   @param a is the structure of the first added ivector
   @param b is the structure of the second added ivector
   @param c is the structure of the result ivector
   
   @b Requests :
   a, b and c should have same dimensions, c has to have allocated array for elements
   which is enough large to hold contents of the result.
   
   @retval 0 on succes
   @retval 1 in case incompatibility sizes of a, b and c vectors
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long addv(const ivector &a, const ivector &b, ivector &c)
{
  if ((a.n != b.n) || (b.n != c.n))
  {
    print_err("cannot add ivectors due to their incompatible dimensions\n"
              "a(%ld) X b(%ld) X c(%ld)", __FILE__, __LINE__, __func__, a.n, b.n, c.n);
    return (1);
  }
  for (long i = 0; i < a.n; i++)
    c.a[i] = a.a[i] + b.a[i];

  return(0);
}



/**
   The function adds double array b to double array given by a, the result is stored in a
   
   @param a is result array
   @param b is added array
   @param n is the number of vector components
   
   @retval always zero
   
   created  by JK
*/
void addv (double *a, const double *b, long n)
{
  for (long i=0;i<n;i++)
    a[i]+=b[i];
}



/**
   The function adds two double arrays a and b, the result is stored in c
   
   @param a is added array
   @param b is added array
   @param c is result array
   @param n is the number of vector components
   
   created  by JK
*/
void addv (const double *a, const double *b, double *c, long n)
{
  for (long i=0;i<n;i++)
    c[i] = a[i]+b[i];
}



/**
   The function adds double array b multiplied by scalar c to double array given by a, the result is stored in a
   
   @param a is result array
   @param b is added array
   @param bc is scalar
   @param n is the number of vector components
   
   @return The result is stored in the array a.
   
   created  by JK
*/
void addmultv (double *a, const double *b, double bc, long n)
{
  for (long i=0;i<n;i++)
    a[i]+=b[i]*bc;
}



/**
   The function adds double array b multiplied by scalar bc to double array given by a, the result is stored in c
   
   @param a is result array
   @param b is added array
   @param bc is scalar
   @param c is scalar
   @param n is the number of vector components
   
   @return The result is stored in the array c.

   created  by TKo
*/
void addmultv (const double *a, const double *b, double bc, double* c, long n)
{
  for (long i=0;i<n;i++)
    c[i] = a[i]+b[i]*bc;
}



/**
   The function subtracts vector given by b from %vector given by a, the result is stored in the c
   
   @param a is the structure of the %vector, from which is subtracted
   @param b is the structure of subtracted %vector
   @param c is the structure of the result %vector
   
   @b Requests :
   a, b and c have to be same dimensions, c has to have allocated array for elements
   which is enough large to hold contents of the result.
   
   @retval 0 on succes
   @retval 1 in case incompatibility sizes of a, b and c vectors
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long subv(const vector &a, const vector &b, vector &c)
{
  if ((a.n != b.n) || (b.n != c.n))
  {
    print_err("cannot subtract vectors due to their incompatible dimensions\n"
              "a(%ld) X b(%ld) X c(%ld)", __FILE__, __LINE__, __func__, a.n, b.n, c.n);
    return (1);
  }
  for (long i = 0; i < c.n; i++)
    c.a[i] = a.a[i] - b.a[i];

  return(0);
}



/**
   The function subtracts vector given by b from %vector given by a, the result is stored in the a
   
   @param a is the structure of the %vector, from which is subtracted
   @param b is the structure of subtracted %vector
   
   @b Requests :
   a and b have to have same dimensions
   
   @retval 0 on succes
   @retval 1 in case incompatibility sizes of a and b vectors
   
   created  08.2011 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long subv(vector &a, const vector &b)
{
  if (a.n != b.n)
  {
    print_err("cannot subtract vectors due to their incompatible dimensions\n"
              "a(%ld) X b(%ld)", __FILE__, __LINE__, __func__, a.n, b.n);
    return (1);
  }
  for (long i = 0; i < a.n; i++)
    a.a[i] -= b.a[i];

  return(0);
}



/**
   The function subtracts %ivector given by b from %ivector given by a, the result is stored in the c
   
   @param a is the structure of the %ivector, from which is subtracted
   @param b is the structure of subtracted %ivector
   @param c is the structure of the result %ivector
   
   @b Requests :
   a, b and c should have same dimensions, c has to have allocated array for elements
   which is enough large to hold contents of the result.
   
   @retval 0 on succes
   @retval 1 in case incompatibility sizes of a, b and c ivectors
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long subv(const ivector &a, const ivector &b, ivector &c)
{
  if ((a.n != b.n) || (b.n != c.n))
  {
    print_err("cannot subtract vectors due to their incompatible dimensions\n"
              "a(%ld) X b(%ld) X c(%ld)", __FILE__, __LINE__, __func__, a.n, b.n, c.n);
    return (1);
  }
  for (long i = 0; i < c.n; i++)
    c.a[i] = a.a[i] - b.a[i];
      
  return(0);
}



/**
   The function subtracts double array given by b from double array given by a, the result is stored in the a
   
   @param a is the array, from which is subtracted
   @param b is the subtracted array
   @param n is number of members of arrays
   
   @retval always zero
   
   created by JK
*/
void subv (double *a, const double *b, long n)
{
  for (long i=0;i<n;i++)
    a[i]-=b[i];
}



/**
  The function subtracts double array given by b from double array by a, the result is stored 
  in the array c.

  @param a - double array of minuend
  @param b - double array of subtrahent
  @param c - double array of resulting difference
  @param n - length of arrays a, b and c

  @return The result is returned in the parameter c.

  Created by JK
*/
void subv (const double *a, const double *b, double *c, long n)
{
  long i;
  
  for (i=0;i<n;i++){
    c[i] = a[i]-b[i];
  }
}



/**
   The function performs vector product, multiplies %vector given by a
   from right with %vector given by b, the result is stored in c
   
   @param a is the structure of the %vector, from which is multiplied
   @param b is the structure of the multiplicating %vector
   @param c is the structure of the result %vector
   
   @b Requests :
   a, b and c have to be same dimension equal 3
   c has to have allocated array for elements
   which is enough large to hold contents of the result.
   
   @retval 0 on succes
   @retval 1 in case incompatibility sizes of a, b and c vectors
   @retval 2 in case a and b vectors have dimension inequal 3
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long crprd(const vector &a, const vector &b, vector &c)
{
  if ((a.n != b.n) || (b.n != c.n))
  {
    print_err("cannot perform cross product of vectors due to their incompatible dimensions\n"
              "a(%ld) X b(%ld) X c(%ld)", __FILE__, __LINE__, __func__, a.n, b.n, c.n);
    return (1);
  }
  if ((a.n != 3) || (b.n != 3))
  {
    print_err("cannot perform cross product of vectors due to wrong dimensions (!= 3)\n"
              "a(%ld) X b(%ld) X c(%ld)", __FILE__, __LINE__, __func__, a.n, b.n, c.n);
    return (2);
  }
  c.a[0] = a.a[1] * b.a[2] - a.a[2] * b.a[1];
  c.a[1] = a.a[2] * b.a[0] - a.a[0] * b.a[2];
  c.a[2] = a.a[0] * b.a[1] - a.a[1] * b.a[0];
  
  return (0);
}



/**
   The function performs vector product, multiplies %ivector given by a
   from right with %ivector given by b, the result is stored in c
   
   @param a is the structure of the %ivector, from which is multiplied
   @param b is the structure of the multiplicating %ivector
   @param c is the structure of the result %ivector
   
   @b Requests :
   a, b and c should have same dimension equal 3
   c has to have allocated array for elements
   which is enough large to hold contents of the result.
   
   
   @retval 0 on succes
   @retval 1 in case incompatibility sizes of a, b and c vectors
   @retval 2 in case a and b vectors have dimension inequal 3
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long crprd(const ivector &a, const ivector &b, ivector &c)
{
  if ((a.n != b.n) || (b.n != c.n))
  {
    print_err("cannot perform cross product of ivectors due to their incompatible dimensions\n"
              "a(%ld) X b(%ld) X c(%ld)", __FILE__, __LINE__, __func__, a.n, b.n, c.n);
    return (1);
  }
  if ((a.n != 3) || (b.n != 3))
  {
    print_err("cannot perform cross product of ivectors due to wrong dimensions (!= 3)\n"
              "a(%ld) X b(%ld) X c(%ld)", __FILE__, __LINE__, __func__, a.n, b.n, c.n);
    return (2);
  }
  c.a[0] = a.a[1] * b.a[2] - a.a[2] * b.a[1];
  c.a[1] = a.a[2] * b.a[0] - a.a[0] * b.a[2];
  c.a[2] = a.a[0] * b.a[1] - a.a[1] * b.a[0];
  
  return (0);
}



/**
   The function performs scalar product of the %vector given by a and %vector given by b,
   the result is stored in the c.

  @param a is the structure of the scalar multiplied %vector
  @param b is the structure of the scalar multiplied %vector
  @param scalar is the real number type of double by which the result is passed

  @b Requests :
  a and b have to be same dimension
  
  @retval 0 : on succes
  @retval 1 : in case incompatibility sizes of a and b
  
  created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
  modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long scprd(const vector &a, const vector &b, double &scalar)
{
  if (a.n != b.n)
  {
    print_err("cannot perform scalar product of vectors due to their incompatible dimensions\n"
              "a(%ld) X b(%ld)", __FILE__, __LINE__, __func__, a.n, b.n);
    return (1);
  }
  scalar = 0.0;
  for (long i=0; i < a.n; i++)
    scalar += a(i) * b(i);
  return (0);
}



/**
   The function computes scalar product of the %vector given by a and %vector given by b, 
   i.e. it returns a.b.

  @param[in] a is the structure of the %vector with left factor
  @param[in] b is the structure of the scalar multiplied %vector

  @b Requests :
  a and b have to be same dimension
  
  @returns The function returns the scalar product of the vectors a and b.
  
  created  03.2024 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
double scprd(const vector &a, const vector &b)
{
  if (a.n != b.n)
  {
    print_err("cannot perform scalar product of vectors due to their incompatible dimensions\n"
              "a(%ld) X b(%ld)", __FILE__, __LINE__, __func__, a.n, b.n);
    abort();
  }
  double scalar = 0.0;
  for (long i=0; i < a.n; i++)
    scalar += a(i) * b(i);
  return (scalar);
}



/**
  The function computes the scalar product of the %vector given by a and %vector given by b.

  @param a - onedimensional array of the first multiplied %vector
  @param b - onedimensional array of the second multiplied %vector
  @param n - number of components in arrays a and b.

  @b Requests :
  a and b must have the same dimension
  
  @return The function returns the scalar product of vectors a and b.

  created  12.2011 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
double scprd(const double *a, const double *b, long n)
{
  double scalar = 0.0;
  for (long i=0; i < n; i++)
    scalar += a[i] * b[i];
  return (scalar);
}



#ifdef INC_OPENMP  
/**
   The function computes scalar product of the %vector given by a and %vector given by b on several threads. 
   It returns a.b. 

  @param[in] a is the structure of the %vector with left factor
  @param[in] b is the structure of the scalar multiplied %vector
  @param[in] nt - number of threads used

  @b Requests :
  a and b have to be same dimension
  
  @returns The function returns the scalar product of the vectors a and b.
  
  created  03.2024 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
double scprd_t(const vector &a, const vector &b, long nt)
{
  if (a.n != b.n)
  {
    print_err("cannot perform scalar product of vectors due to their incompatible dimensions\n"
              "a(%ld) X b(%ld)", __FILE__, __LINE__, __func__, a.n, b.n);
    abort();
  }
  double scalar = 0.0;
  #pragma omp parallel num_threads(nt)
  {
   #pragma omp for reduction(+:scalar)
    for (long i=0; i < a.n; i++)
      scalar += a(i) * b(i);
  }
  return (scalar);
}



/**
  The function computes the scalar product of the %vector given by a and %vector given by b
  on several threads.

  @param a[in]  - onedimensional array of the first multiplied %vector
  @param b[in]  - onedimensional array of the second multiplied %vector
  @param n[in]  - number of components in arrays a and b.
  @param nt[in] - number of threads used

  @b Requests :
  a and b must have the same dimension
  
  @return The function returns the scalar product of vectors a and b.

  created  11.2022 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
double scprd_t(const double *a, const double *b, long n, long nt)
{
  double scalar = 0.0;
  #pragma omp parallel num_threads(nt)
  {
   #pragma omp for reduction(+:scalar)
   for (long i=0; i < n; i++)
     scalar += a[i] * b[i];
  }
  return (scalar);
}
#endif



/**
   The function performs scalar product of the %ivector given by a and %ivector given by b,
   the result is stored in the c
   
   @param a is the structure of the scalar multiplied %ivector
   @param b is the structure of the scalar multiplied %ivector
   @param scalar is the integer number type of long by which the result is passed
   
   @b Requests :
   a and b should have same dimension
   
   @retval 0 on succes
   @retval 1 in case incompatibility sizes of a and b
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long scprd(const ivector &a, const ivector &b, long &scalar)
{
  if (a.n != b.n)
  {
    print_err("cannot perform scalar product of ivectors due to their incompatible dimensions\n"
              "a(%ld) X b(%ld)", __FILE__, __LINE__, __func__, a.n, b.n);
    return (1);
  }
  scalar = 0;
  for (long i=0; i < a.n; i++)
    scalar += a.a[i] * b.a[i];
  return (0);
}



/**
   The function performs scalar product of the double array given by a and double array given by b,
   the result is stored in the c
   
   @param a is the first array
   @param b is the second array
   @param n is number of members of arrays
   
   created  by JK
*/
double ss (double *a,double *b,long n)
{
  double s = 0.0;
  for (long i=0;i<n;i++)
    s+=a[i]*b[i];
  
  return s;
}



/**
   The function multiplies %vector given by u by real constant c,
   the result is stored in the double %vector v, i.e. v = c.u
   
   @param c is the real number type double
   @param u is the structure of the %vector, which is multiplied
   @param v is the structure of the result %vector
   
   @b Requests :
   u and v have to be same dimensions, u has to have allocated array for elements
   which is enough large to hold contents of the result.
   
   @retval 0 : on succes
   @retval 1 : in case incompatibility sizes of u and v vectors
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long cmulv(double c, const vector &u, vector &v)
{
  if (u.n != v.n)
  {
    print_err("cannot multiply vector by constant due to their incompatible dimensions\n"
              "u(%ld) X v(%ld)", __FILE__, __LINE__, __func__, u.n, v.n);
    return (1);
  }
  for (long i=0; i < u.n; i++)
    v.a[i] = c * u.a[i];

  return (0);
}



/**
   The function multiplies %ivector given by u by real constant c,
   the result is stored in the double %ivector v, i.e. v = c.u
   
   @param c is the real number type double
   @param u is the structure of the %ivector, which is multiplied
   @param v is the structure of the result %vector
   
   @b Requests :
   u and v should have same dimensions, u has to have allocated array for elements
   which is enough large to hold contents of the result.
   
   @retval 0 on succes
   @retval 1 in case incompatibility sizes of u and v vectors
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long cmulv(double c, const ivector &u, vector &v)
{
  if (u.n != v.n)
  {
    print_err("cannot multiply ivector by constant due to their incompatible dimensions\n"
              "u(%ld) X v(%ld)", __FILE__, __LINE__, __func__, u.n, v.n);
    return (1);
  }
  for (long i=0; i < u.n; i++)
    v.a[i] = c * u.a[i];
  
  return (0);
}



/**
   The function multiplies %ivector given by u by integer constant c,
   the result is stored in the %vector v, i.e. v = c.u
   
   @param c is the integer number type long
   @param u is the structure of the %ivector, which is multiplied
   @param v is the structure of the result %ivector
   
   @b Requests :
   u and v should have same dimensions, u has to have allocated array for elements
   which is enough large to hold contents of the result.
   
   @retval 0 on succes
   @retval 1 in case incompatibility sizes of u and v %vector
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long cmulv(long c, const ivector &u, ivector &v)
{
  if (u.n != v.n)
  {
    print_err("cannot multiply ivector by constant due to their incompatible dimensions\n"
              "u(%ld) X v(%ld)", __FILE__, __LINE__, __func__, u.n, v.n);
    return (1);
  }
  for (long i=0; i < u.n; i++)
    v.a[i] = c * u.a[i];
  
  return (0);
}



/**
   The function multiplies %vector given by a with constant c,
   the result is stored in a
   
   @param c is the real number
   @param a is the structure of the %vector, which is multiplied
   
   @return always zero
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long cmulv(double c, vector &a)
{
  for (long i=0; i < a.n; i++)
    a.a[i] *= c;
  
  return (0);
}



/**
   The function multiplies %ivector given by a with integer constant c,
   the result is stored in a
   
   @param c is the integer number
   @param a is the structure of the %ivector, which is multiplied
   
   @retval 0 always
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long cmulv(long c, ivector &a)
{
  for (long i=0; i < a.n; i++)
    a.a[i] *= c;
  
  return (0);
}



/**
   The function multiplies double array given by a with constant c,
   the result is stored in a
   
   @param c is the real number
   @param a is the array, which is multiplied
   @param n is number of members of arrays
   
   @return always zero
   
   created by JK
*/
void cmulv(double c,double *a,long n)
{
  for (long i=0; i < n; i++)
    a[i] *= c;
  
}



/**
   The function multiplies double array given by a with constant c,
   the result is stored in u.
   
   @param c is the real number
   @param a is the array, which is multiplied
   @param u is the array, where the result is stored
   @param n is number of members of arrays
   
   @return always zero
   
   created by JK
*/
void cmulv(double c,double *a, double *u, long n)
{
  for (long i=0; i < n; i++)
    u[i] = c*a[i];
}



/**
   The function computes norm of %vector given by a.
   
   @param a - structure of %vector components 
   
   @return The function returns computed norm of the %vector.
   
   created  02.2010 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
double normv(const vector &a)
{
  long i;
  double norm = 0.0;
  
  for (i=0; i < a.n; i++)
    norm += a.a[i] * a.a[i];
  return sqrt(norm);
}



/*
  The function computes norm of the %vector a.

  @param a - array with vector values
  @param n - number of %vector components (i.e. size of the array a)
  
  @return The function returns computed norm of the %vector.

  created  02.2010 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
double normv(const double *a,long n)
{
  long i;
  double norm = 0.0;
  
  for (i=0;i<n;i++){
    norm += a[i]*a[i];
  }
  return sqrt(norm);
}



/**
   The function computes norm of the %ivector given by a,
   the result is stored in the parameter norm
   
   @param a is the structure of the ivector, which size is requested
   
   @return Returns computed size of %vector norm.
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
double normv(const ivector &a)
{
  long i;
  double norm = 0.0;
  for (i=0; i < a.n; i++)
    norm += a.a[i] * a.a[i];
  norm = sqrt(norm);
  
  return(norm);
}



/**
   The function computes selective norm of %vector given by a. Only range of nc %vector
   components starting from component fi is considered in the norm computation.
   
   @param a - structure of %vector components
   @param fi - index of the first %vector component considered in the computation (in C/C++ notation)
   @param nc - the number of %vector components considered in the computation
   
   @return The function returns computed selective norm of the %vector.
   
   created  02.2016 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
double snormv(const vector &a, long fi, long nc)
{
  double norm = 0.0;
  long li = fi+nc;

  if ((fi >= a.n) || (li > a.n)) 
  {
    print_err("invalid first index or number of components (tnc=%ld; fi=%ld; nc=%ld)", 
              __FILE__, __LINE__, __func__, a.n, fi, nc); 
    abort();
  }
  for (long i=fi; i < li; i++)
    norm += a.a[i] * a.a[i];
  return sqrt(norm);
}



/*
   The function computes selective norm of %vector given by a. Only range of nc %vector
   components starting from component fi is considered in the norm computation.
   
   @param a - array with vector values
   @param n - the total number of %vector components (i.e. size of the array a)
   @param fi - index of the first %vector component considered in the computation (in C/C++ notation)
   @param nc - the number of %vector components considered in the computation
   
   @return The function returns computed selective norm of the %vector.

  created  02.2016 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
double snormv(const double *a,long n, long fi, long nc)
{
  long i;
  double norm = 0.0;
  long li = fi+nc;

  if ((fi >= n) || (li > n)) 
  {
    print_err("invalid first index or number of components (tnc=%ld; fi=%ld; nc=%ld)", 
              __FILE__, __LINE__, __func__, n, fi, nc); 
    abort();
  }

  for (i=fi;i<fi+nc;i++){
    norm += a[i]*a[i];
  }
  return sqrt(norm);
}



/**
   The function normalizes components of %vector given by a, |a| = 1 after function call.
   
   @param[in/out] a - structure of %vector components 
   
   @return The function normalizes %vector in the argument.
   
   created  06.2020 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void normalize(vector &a)
{
  long i;
  double norm = 0.0;
  for (i=0; i < a.n; i++)
    norm += a.a[i] * a.a[i];
  norm = 1.0/sqrt(norm);
  for (i=0; i < a.n; i++)
    a.a[i] *= norm;
}



/*
  The function normalizes components of the %vector a, |a| = 1 after function call.

  @param a[in/out] - array with %vector component values
  @param n[in] - number of %vector components (i.e. dimension of the array a)
  
  @return The function normalizes %vector in the argument.

  created  06.2020 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void normalize(double *a,long n)
{
  long i;
  double norm = 0.0;
  for (i=0; i < n; i++)
    norm += a[i] * a[i];
  norm = 1.0/sqrt(norm);
  for (i=0; i < n; i++)
    a[i] *= norm;
}



/**
   The function normalizes components of %vector given by a, |a| = 1 after function call.
   
   @param[in/out] a - structure of %vector components 
   
   @return The function normalizes %vector in the argument.
   
   created  06.2020 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void normalize(const vector &a, vector &b)
{
  long i;
  double norm = 0.0;
  for (i=0; i < a.n; i++)
    norm += a.a[i] * a.a[i];
  norm = 1.0/sqrt(norm);
  for (i=0; i < a.n; i++)
    b.a[i] = norm*a.a[i];
}



/*
  The function normalizes components of the %vector a, |a| = 1 after function call.

  @param a[in] - array with %vector component values to be normalized
  @param n[in] - number of %vector components (i.e. dimension of the array a)
  @param b[out] - array with normalized %vector components
  
  @return The function normalizes %vector in the argument.

  created  06.2020 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void normalize(double *a, long n, double *b)
{
  long i;
  double s = 0.0;
  for (i=0; i < n; i++)
    s += a[i] * a[i];
  s = 1.0/sqrt(s);
  for (i=0; i < n; i++)
    b[i] = s*a[i];
}



/**
  Returns maximum of the %vector components.
  @param a - given %vector
*/
/*long maxcompv(const vector &a)
  {
  long ret = 0;
  
  if (a.n < 1)
  return ret;
  
  double max = a[0];
  for(long i=0; i<a.n; i++){
  if (a[i] > ret){
  max = a[i];
  ret = i;
  }
  }
  return ret;
  }
*/



/**
  Returns maximum of the %vector components.
  @param a - given %vector
  @param n - length of %vector a
  
*/
/* long maxcompv(const double *a, long n)
   {
   long ret = 0;
   
   if (n < 1)
   return ret;
   
   double max = a[0];
   
   for(long i=0; i<n; i++){
   if (a[i] > max){
   max = a[i];
   ret = i;
   }
   }
   return ret;
   }
*/


/**
   The function computes cosine of the angle between %vector given by a and b
   the result is stored in the cos
   
   @param a is the structure of the first %vector
   @param b is the structure of the second %vector
   @param cos is real number, type double, where result is stored
   
   @b Requests :
   vectors a nd b should have the same dimension
   
   @retval 0 on succes
   @retval 1 in case
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long cosav(const vector &a, const vector &b, double &cos)
{
  if (a.n != b.n)
  {
    print_err("cannot determine cos angle of vectors due to their incompatible dimensions\n"
              "a(%ld) X b(%ld)", __FILE__, __LINE__, __func__, a.n, b.n);
    return (1);
  }
  
  double scalar, norm_a, norm_b;
  scprd(a, b, scalar);
  norm_a = normv(a);
  norm_b = normv(b);
  cos = scalar/(norm_a * norm_b);
  
  return (0);
}



/**
   The function computes cosine of the angle between %ivector given by a and b
   the result is stored in the cos
   
   @param a is the structure of the first %ivector
   @param b is the structure of the second %ivector
   @param cos is real number, type double, where result is stored
   
   @b Requests :
   vectors a nd b should have the same dimension

   @retval 0 on succes
   @retval 1 in case
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long cosav(const ivector &a, const ivector &b, double &cos)
{
  if (a.n != b.n)
  {
    print_err("cannot determine cos angle of ivectors due to their incompatible dimensions\n"
              "a(%ld) X b(%ld)", __FILE__, __LINE__, __func__, a.n, b.n);
    return (1);
  }
  
  long scalar;
  double norm_a, norm_b;
  scprd(a, b, scalar);
  norm_a = normv(a);
  norm_b = normv(b);
  cos = scalar/(norm_a * norm_b);
  
  return (0);
}



/**
   The function changes components of %vector a to their absolute value
   
   @param a is the structure of the first %ivector
   
   created  7.2008 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void absv(vector &a)
{
  long i;
  for (i=0; i<a.n; i++)
    a[i] = fabs(a[i]);
}



/**
   The function changes components of %vector a to their absolute value
   
   @param a is the structure of the first %ivector
   
   created  7.2008 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void absv(ivector &a)
{
  long i;
  for (i=0; i<a.n; i++)
    a[i] = labs(a[i]);
}



/**
   The function prints out the contents of the %vector u to the out file,
   Optionally, the file, precision and field width can be specified file.
   
   @param u is the structure of the %vector which is printed
   
   @b Optionally :
   @param out   is the structure with opened file for %vector output (default is stdout)
   @param prec  is desired precision of the %vector elements on the output (default is 3)
   @param width is desired width of the %vector elements on the output (default is 11)
   
   @retval 0 : on succes
   @retval 1 : in case NULL pointer of the out parameter
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long printv(const vector &u, FILE *out, int prec, int width)
{
  if (out == NULL)
  {
    print_err("cannot print vector - parameter for output is a NULL pointer", 
              __FILE__, __LINE__, __func__);
    return (1);
  }
  for (long i = 0; i < u.n; i++)
    fprintf(out, "% *.*e\n", width, prec, u.a[i]);
  return (0);
}



/**
   The function prints out the contents of the %vector u given by array to the out file.
   Optionally, file, precision and field width can be specified
   
   @param u - pointer to array of double values which will be printed
   @param n - size of array u
   
   @b Optionally :
   @param out   is the structure with opened file for %vector output (default is stdout)
   @param prec  is desired precision of the %vector elements on the output (default is 3)
   @param width is desired width of the %vector elements on the output (default is 11)
   
   @retval 0 : on succes
   @retval 1 : in case NULL pointer of the out parameter
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long printv(const double *u, long n, FILE *out, int prec, int width)
{
  if (out == NULL)
  {
    print_err("cannot print vector - parameter for output is a NULL pointer", 
              __FILE__, __LINE__, __func__);
    return (1);
  }
  for (long i = 0; i < n; i++)
    fprintf(out, "% *.*e\n", width, prec, u[i]);
  return (0);
}



/**
   The function prints contents of the %vector u into the file out,
   
   @param out -  pointer to the structure with opened output file
   @param u   -  %vector which contents will be printed
   
   @b Requirements :
   Vector u has to be allocated.
   Number of read values is given by the %vector n member.
   
   @retval 0 : on succes
   @retval 1 : in case NULL pointer of the in parameter
   @retval 2 : in case error reading number from the file
   
   TKr, 11/02/2013 accroding to readv(XFILE *in, vector u)
*/
long printv(FILE *out, vector &u)
{
  for (long i = 0; i < u.n; i++)
  {
    fprintf(out, "  % le", u[i]);
  }
  return (0);
}



/**
   The function prints out the contents of the %ivector u to the file out.
   Optionally, the file and field width can be specified.
   
   @param u is the structure of the %ivector which is printed
   
   @b Optionally :
   @param out   is the structure with opened file for %ivector output (default is stdout)
   @param width is desired width of the %ivector elements on the output (default is 11)
   
   @retval 0 on succes
   @retval 1 in case NULL pointer of the out parameter
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long printv(const ivector &u, FILE *out, int width)
{
  if (out == NULL)
  {
    print_err("cannot print ivector - parameter for output is a NULL pointer", 
              __FILE__, __LINE__, __func__);
    return (1);
  }
  for (long i = 0; i < u.n; i++)
    fprintf(out, "% *ld\n", width, u.a[i]);
  return (0);
}



/**
   The function prints out the contents of the %bvector u to the file out.
   Optionally, the file and field width can be specified.
   
   @param u is the structure of the %bvector which is printed
   
   @b Optionally :
   @param out   is the structure with opened file for %bvector output (default is stdout)
   @param width is desired width of the %bvector elements on the output (default is 11)
   
   @retval 0 on succes
   @retval 1 in case NULL pointer of the out parameter
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long printv(const bvector &u, FILE *out, int width)
{
  if (out == NULL)
  {
    print_err("cannot print bvector - parameter for output is a NULL pointer", 
              __FILE__, __LINE__, __func__);
    return (1);
  }
  for (long i = 0; i < u.n; i++)
    fprintf(out, "% *d\n", width, u.a[i]);
  return (0);
}



/**
   The function prints out the contents of the %vector u (integer values are assumed) 
   given by array. Optionally, the file and field width can be specified.
   
   @param u - pointer to array of integer values which will be printed
   @param n - size of array u
   
   @b Optionally :
   @param out   is the structure with opened file for %ivector output (default is stdout)
   @param width is desired width of the %ivector elements on the output (default is 11)
   
   @retval 0 on succes
   @retval 1 in case NULL pointer of the out parameter
   
   created  29.5.2000 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified  9.5.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long printv(const long *u, long n, FILE *out, int width)
{
  if (out == NULL)
  {
    print_err("cannot print vector - parameter for output is a NULL pointer", 
              __FILE__, __LINE__, __func__);
    return (1);
  }
  for (long i = 0; i < n; i++)
    fprintf(out, "%*ld\n", width, u[i]);
  return (0);
}



/**
   The function reads components of the %vector u from the string in,

   @param in    pointer to the input string
   @param u     vector which components will be read
   
   @b Requirements :
   Vector u have to be allocated.
   Number of read values is given by the vector n member.
   
   @retval 0 : on succes
   @retval 1 : in case NULL pointer of the in parameter
   @retval 2 : in case error reading number from the file
   
   created  6.9.2001 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
   modified 16.4.2024 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long readv(char *in, vector &u)
{
  if (in == NULL)
  {
    print_err("cannot read vector - parameter for input string is a NULL pointer", 
              __FILE__, __LINE__, __func__);
    return (1);
  }
  int br, tbr = 0;
  for (long i = 0; i < u.n; i++)
  {
    if (sscanf(in+tbr, "%le%n", &u.a[i], &br) != 1)
      return(2);
    tbr += br;
  }
  return (0);
}



/**
   The function reads components of the integer %vector u from the string in,

   @param in    pointer to the input string
   @param u     integer %vector which components will be read
   
   @b Requirements :
   Vector u have to be allocated.
   Number of read values is given by the ivector n member.
   
   @retval 0 : on succes
   @retval 1 : in case NULL pointer of the in parameter
   @retval 2 : in case error reading number from the file
   
   created  30.4.2023 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long readv(char *in, ivector &u)
{
  if (in == NULL)
  {
    print_err("cannot read ivector - parameter for input string is a NULL pointer", 
              __FILE__, __LINE__, __func__);
    return (1);
  }
  int br, tbr = 0;
  for (long i = 0; i < u.n; i++)
  {
    if (sscanf(in+tbr, "%ld%n", &u.a[i], &br) != 1)
      return(2);
    tbr += br;
  }
  return (0);
}



/**
   The function reads components of the integer %vector u from the string in,

   @param in    pointer to the input string
   @param u     integer %vector which components will be read
   
   @b Requirements :
   Vector u have to be allocated.
   Number of read values is given by the bvector n member.
   
   @retval 0 : on succes
   @retval 1 : in case NULL pointer of the in parameter
   @retval 2 : in case error reading number from the file
   
   created  16.4.2024 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long readv(char *in, bvector &u)
{
  if (in == NULL)
  {
    print_err("cannot read bvector - parameter for input string is a NULL pointer", 
              __FILE__, __LINE__, __func__);
    return (1);
  }
  int br, tbr=0, aux;
  char auxs[6];
  for (long i = 0; i < u.n; i++)
  {
    if (sscanf(in+tbr, "%5s%n", auxs, &br) == 1){
      tbr += br;
      if (sscanf(auxs, "%d", &aux) == 1){
        u.a[i] = aux;
        continue;
      }
      if (strcmp(auxs, "true")){
        u.a[i] = true;
        continue;
      }
      if (strcmp(auxs, "false")){
        u.a[i] = false;
        continue;
      }
      return(2);
    }
  }
  return (0);
}



/**
  The function prints out the complete record of the %vector u to the out file.
  Optionally, file, precision and field width can be specified
   
  @param[in] u       - a %vector whose record will be written to the out
   
  @b Optionally :
  @param[in] out - a structure with opened file for %vector output (default is stdout)
  @param[in] prec    - desired precision of the %vector elements on the output (default is 3)
  @param[in] width   - desired width of the %vector elements on the output (default is 11)
   
  @retval 0 : on succes
  @retval 1 : in case NULL pointer of the out parameter
   
  created  18.4.2024 by Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
long printva(const vector &u, FILE *out, int prec, int width)
{
  if (out == NULL){
    print_err("cannot print vector - parameter for output is a NULL pointer", 
              __FILE__, __LINE__, __func__);
    return (1);
  }
  fprintf(out, "%ld\n", u.n);
  for (long i = 0; i < u.n; i++)
    fprintf(out, "% *.*e\n", width, prec, u[i]);
  return (0);
}



/**
  The function prints complete record of the %vector u into the file out,
   
  @param out -  pointer to the structure with opened output file
  @param u   -  %vector which contents will be printed
  
  @b Requirements :
  Vector u has to be allocated.
   
  @retval 0 : on succes
   
  Created by TKo, 18.04.2024
*/
long printva(FILE *out, const vector &u)
{
  fprintf(out, "%ld\n", u.n);
  for (long i = 0; i < u.n; i++)
    fprintf(out, "%le\n", u[i]);

  return (0);
}



/**
   function extracts ncomp components from %vector b
   and put them into %vector a

   @param a - %vector containing extracted components
   @param b - %vector from where components are extracted
   @param fi - first index
   @param ncomp - number of components

   30.11.2002
*/
void extract (vector &a, const vector &b, long fi, long ncomp)
{
  long i,j;

  j=fi;
  for (i=0;i<ncomp;i++){
    a[i]=b[j];
    j++;
  }
}



/**
   Function extracts positive components from vector a
   and puts them into vector b in corresponding position as in the vector a.
   If component is negative the zero is stored.

   @param a - %vector containing compared components
   @param b - %vector where extracted components are stored

   3.3.2003
*/
long extractposv (const vector &a, vector &b)
{
  long i;
  if (a.n != b.n)
  {
    print_err("cannot extract positive components - incompatible dimensions of vectors\n"
              "a(%ld) X b(%ld)", __FILE__, __LINE__, __func__, a.n, b.n);
    return (1);
  }
  for (i = 0; i < a.n; i++)
  {
    if (a[i] > 0.0)
      b[i] = a[i];
    else
      b[i] = 0.0;
  }
  return 0;
}



/**
   Function extracts negative components from vector a
   and puts them into vector b in corresponding position as in the vector a.
   If component is positive the zero is stored.

   @param a - %vector containing compared components
   @param b - %vector where extracted components are stored

   25.2.2004
*/
long extractnegv (vector &a,vector &b)
{
  long i;
  if (a.n != b.n)
  {
    print_err("cannot extract negative components - incompatible dimensions of vectors\n"
              "a(%ld) X b(%ld)", __FILE__, __LINE__, __func__, a.n, b.n);
    return (1);
  }
  for (i = 0; i < a.n; i++)
  {
    if (a[i] < 0.0)
      b[i] = a[i];
    else
      b[i] = 0.0;
  }
  return 0;
}



/**
  The function sorts ivector v by Shell sort algorithm in ascending order.

  @param v - vector of sorted numbers

  @return The function does not return anything but changes order of v elements.

  Created by Tomas Koudelka, 2.7.2013
*/
void shell_sort(ivector &v)
{
  long flag = 1, d = v.n, i, temp;

  while(flag || (d>1))   // boolean flag (true when not equal to 0)
  {
    flag = 0;           // reset flag to 0 to check for future swaps
    d = (d+1)/2;
    for (i=0; i<(v.n-d); i++)
    {
      if (v[i+d] < v[i])
      {
        temp   = v[i+d]; // swap items at positions i+d and d
        v[i+d] = v[i];
        v[i]   = temp;
        flag   = 1;        // indicate that a swap has occurred
      }
    }
  }
  return;
}



/**
  The function sorts array of integer numbers v by Shell sort algorithm in ascending order.

  @param v - array of sorted numbers
  @param n - length of array v

  @return The function does not return anything but changes order of v elements.

  Created by Tomas Koudelka, 3.7.2013
*/
void shell_sort(long *v, long n)
{
  long flag = 1, d = n, i, temp;

  while(flag || (d>1))   // boolean flag (true when not equal to 0)
  {
    flag = 0;           // reset flag to 0 to check for future swaps
    d = (d+1)/2;
    for (i=0; i<(n-d); i++)
    {
      if (v[i+d] < v[i])
      {
        temp   = v[i+d]; // swap items at positions i+d and d
        v[i+d] = v[i];
        v[i]   = temp;
        flag   = 1;        // indicate that a swap has occurred
      }
    }
  }
  return;
}



/**
  Function writes error on standard error output device, format specifiers are allowed.

  Parameters:
  @param emsg - string with error message and format specifiers, one terminating character \n will automatically added
  @param errfile - string with source file name where the error was generated
  @param errln - line number in source file  where the error was generated
  @param errfunc - string with function name where the error was generated
  @param ...    - other parameters used in format string emsg passed to the vfprintf function

  @return  Function does not return anything.
  
  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void print_err(const char *emsg, const char *errfile, int errln, const char *errfunc, ...)
{
  va_list params;
  fflush(stdout); // force standard output to be written
  fflush(stderr); // force standard error output to be written
  fprintf(stderr, "\n\nError: ");
  va_start(params, errfunc);
  vfprintf(stderr, emsg, params);
  va_end(params);
  fprintf(stderr, "\n");  
  fprintf(stderr, "detected in source file %s, line %d, function %s\n", errfile, errln, errfunc);
  //stack_trace(stderr, 2);
  fflush(stderr);
  fflush(stdout); // force standard output to be written
}



/**
  Function writes warning on standard error output device, format specifiers are allowed in wmsg.

  Parameters:
  @param wmsg - string with warning message, one terminating character \n will automatically added
  @param warnfile - string with source file name where the error was generated
  @param warnln - line number in source file  where the warning was generated
  @param warnfunc - string with function name where the warning was generated
  @param ...    - other parameters used in format string wmsg passed to the vfprintf function

  @return Function does not return anything.
  
  created  12.2006,  Tomas Koudelka, tomas.koudelka@fsv.cvut.cz
*/
void print_warning(const char *wmsg, const char *warnfile, int warnln, const char *warnfunc, ...)
{
  va_list params;
  fflush(stdout); // force standard output to be written
  fflush(stderr); // force standard error output to be written
  fprintf(stderr, "\n\nWarning: ");
  va_start(params, warnfunc);
  vfprintf(stderr, wmsg, params);
  va_end(params);
  fprintf(stderr,"\n");
  fprintf(stderr, "detected in file %s, line %d, function %s\n", warnfile, warnln, warnfunc);
  //stack_trace(stderr, 2);
  fflush(stderr);
}
