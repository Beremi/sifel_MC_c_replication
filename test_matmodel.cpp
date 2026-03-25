#include "matmodel.h"
#include "vector.h"
#include "matrix.h"

#include <cmath>
#include <cstdio>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace {

struct TestCase
{
  const char *name;
  double strain[4];
  int expected_return;
};

void fill_vector4(const double src[4], vector &dst)
{
  reallocv(4, dst);
  for (long i=0; i<4; i++)
    dst[i] = src[i];
}

double norm44(const matrix &a)
{
  double s = 0.0;
  for (long i=0; i<4; i++)
  {
    for (long j=0; j<4; j++)
      s += a(i,j)*a(i,j);
  }
  return std::sqrt(s);
}

void diff44(const matrix &a, const matrix &b, matrix &out)
{
  reallocm(4, 4, out);
  for (long i=0; i<4; i++)
  {
    for (long j=0; j<4; j++)
      out(i,j) = a(i,j) - b(i,j);
  }
}

void compute_response(matmodel &mm,
                      const double strain_data[4],
                      const double eq_data[4],
                      vector &stress,
                      vector &statev)
{
  vector strain;
  vector eqstatev;

  fill_vector4(strain_data, strain);
  fill_vector4(eq_data, eqstatev);
  mm.nlstresses(strain, eqstatev, stress, statev);
}

void finite_difference_tangent(matmodel &mm,
                               const double strain[4],
                               const double eqother[4],
                               matrix &dnum)
{
  const double h = 1.0e-7;
  reallocm(4, 4, dnum);

  for (long j=0; j<4; j++)
  {
    double ep[4], em[4];
    vector stress_p, stress_m;
    vector statev_p, statev_m;

    for (long i=0; i<4; i++)
    {
      ep[i] = strain[i];
      em[i] = strain[i];
    }
    ep[j] += h;
    em[j] -= h;

    compute_response(mm, ep, eqother, stress_p, statev_p);
    compute_response(mm, em, eqother, stress_m, statev_m);

    for (long i=0; i<4; i++)
      dnum(i,j) = (stress_p[i] - stress_m[i])/(2.0*h);
  }
}

int read_params(matmodel &mm)
{
  FILE *in = tmpfile();
  if (in == NULL)
    return 1;

  std::fprintf(in, "%.17g %.17g %.17g %.17g\n",
               20000.0, 0.49, 50.0, M_PI/9.0);
  std::rewind(in);
  const long ret = mm.read(in);
  std::fclose(in);
  return static_cast<int>(ret);
}

} // namespace

int main()
{
  matmodel mm;
  if (read_params(mm) != 0)
  {
    std::fprintf(stderr, "Invalid material parameters.\n");
    return 1;
  }

  const TestCase tests[] = {
    {"elastic",    { 1.0e-5,  0.0,    0.0,   0.0  }, MATMODEL_RETURN_ELASTIC},
    {"smooth",     {-5.0e-3,  0.0,    0.0,   5.0e-3}, MATMODEL_RETURN_SMOOTH},
    {"left_edge",  {-5.0e-3,  5.0e-3, 0.0,   5.0e-3}, MATMODEL_RETURN_LEFT_EDGE},
    {"right_edge", {-5.0e-3, -5.0e-3, 0.0,   1.0e-2}, MATMODEL_RETURN_RIGHT_EDGE},
    {"apex",       { 0.0,     0.0,    0.0,   5.0e-3}, MATMODEL_RETURN_APEX}
  };

  const double eqother[4] = {0.0, 0.0, 0.0, 0.0};
  bool ok = true;

  for (size_t it=0; it<sizeof(tests)/sizeof(tests[0]); it++)
  {
    vector stress;
    vector statev;
    vector eqnext;
    matrix d, dnum, derr;
    double relerr;

    compute_response(mm, tests[it].strain, eqother, stress, statev);
    mm.stiffmat_from_statev(statev, d);
    finite_difference_tangent(mm, tests[it].strain, eqother, dnum);
    diff44(d, dnum, derr);
    mm.updateval(statev, eqnext);

    relerr = norm44(derr)/(norm44(dnum) + 1.0e-14);

    std::printf("[%s]\n", tests[it].name);
    std::printf("  return type : %.0f\n", statev[MATMODEL_IO_RETURN_TYPE]);
    std::printf("  stress      : [% .10e, % .10e, % .10e, % .10e]\n",
                stress[0], stress[1], stress[2], stress[3]);
    std::printf("  eps_p       : [% .10e, % .10e, % .10e, % .10e]\n",
                statev[0], statev[1], statev[2], statev[3]);
    std::printf("  D rel. err. : %.6e\n", relerr);

    if (static_cast<int>(statev[MATMODEL_IO_RETURN_TYPE] + 0.5) != tests[it].expected_return)
    {
      std::printf("  ERROR: unexpected return type.\n");
      ok = false;
    }

    if (relerr > 1.0e-6)
    {
      std::printf("  ERROR: tangent check failed.\n");
      ok = false;
    }

    for (long i=0; i<4; i++)
    {
      if (std::fabs(eqnext[i] - statev[i]) > 1.0e-14)
      {
        std::printf("  ERROR: updateval failed at epsp[%ld].\n", i);
        ok = false;
      }
    }

    std::printf("\n");
  }

  if (!ok)
  {
    std::printf("TEST STATUS: FAIL\n");
    return 1;
  }

  std::printf("TEST STATUS: PASS\n");
  return 0;
}
