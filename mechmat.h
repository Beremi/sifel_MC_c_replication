#ifndef SIFEL_COMPAT_MECHMAT_H
#define SIFEL_COMPAT_MECHMAT_H

#include "vector.h"

struct mechmat_ip
{
  long ncompstr;
  double *other;
  double *eqother;
  double *nonloc;

  mechmat_ip() : ncompstr(0), other(NULL), eqother(NULL), nonloc(NULL) {}
};

class mechmat
{
 public:
  mechmat_ip *ip;

  mechmat() : ip(NULL) {}

  void givestrain(long, long, vector &) const {}
  void givestress(long, long, vector &) const {}
  void storestress(long, long, const vector &) {}
  long givencompeqother(long, long) const { return 0; }
};

#endif
