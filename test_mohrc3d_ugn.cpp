#include "iotools.h"
#include "mohrc3d_ugn.h"
#include "vector.h"
#include "matrix.h"

#include <cmath>
#include <cstdio>
#include <cstddef>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace {

struct TestCase
{
  const char *name;
  double strain[6];
  int expected_return;
};

void fill_vector6(const double src[6], vector &dst)
{
  reallocv(6, dst);
  for (long i=0; i<6; i++)
    dst[i] = src[i];
}

void fill_state_buffer(const double history[6], vector &dst)
{
  reallocv(mohrc3d_ugn::NCOMP_EQOTHER, dst);
  nullv(dst);
  for (long i=0; i<6; i++)
    dst[i] = history[i];
}

double norm66(const matrix &a)
{
  double s = 0.0;
  for (long i=0; i<6; i++)
  {
    for (long j=0; j<6; j++)
      s += a(i,j)*a(i,j);
  }
  return std::sqrt(s);
}

void diff66(const matrix &a, const matrix &b, matrix &out)
{
  reallocm(6, 6, out);
  for (long i=0; i<6; i++)
  {
    for (long j=0; j<6; j++)
      out(i,j) = a(i,j) - b(i,j);
  }
}

void compute_response(mohrc3d_ugn &mm,
                      const double strain_data[6],
                      const double eq_data[6],
                      vector &stress,
                      vector &statev)
{
  vector strain;
  vector eqstatev;

  fill_vector6(strain_data, strain);
  fill_state_buffer(eq_data, eqstatev);
  mm.nlstresses(strain, eqstatev, stress, statev);
}

void finite_difference_tangent(mohrc3d_ugn &mm,
                               const double strain[6],
                               const double eqother[6],
                               matrix &dnum)
{
  const double h = 1.0e-8;
  reallocm(6, 6, dnum);

  for (long j=0; j<6; j++)
  {
    double ep[6], em[6];
    vector stress_p, stress_m;
    vector statev_p, statev_m;

    for (long i=0; i<6; i++)
    {
      ep[i] = strain[i];
      em[i] = strain[i];
    }
    ep[j] += h;
    em[j] -= h;

    compute_response(mm, ep, eqother, stress_p, statev_p);
    compute_response(mm, em, eqother, stress_m, statev_m);

    for (long i=0; i<6; i++)
      dnum(i,j) = (stress_p[i] - stress_m[i])/(2.0*h);
  }
}

int read_params(mohrc3d_ugn &mm)
{
  FILE *in = tmpfile();
  if (in == NULL)
    return 1;

  std::fprintf(in, "%.17g %.17g %.17g %.17g\n",
               20000.0, 0.49, 50.0, M_PI/9.0);
  std::rewind(in);
  XFILE xin = {in};
  const long ret = mm.read(&xin);
  std::fclose(in);
  return static_cast<int>(ret);
}

} // namespace

int main()
{
  mohrc3d_ugn mm;
  if (read_params(mm) != 0)
  {
    std::fprintf(stderr, "Invalid 3D material parameters.\n");
    return 1;
  }

  const TestCase tests[] = {
    {"elastic",    { 1.0e-5,  0.0,     0.0,    0.0, 0.0, 0.0}, mohrc3d_ugn::RET_ELASTIC},
    {"smooth",     {-5.0e-3,  0.0,     5.0e-3, 0.0, 0.0, 0.0}, mohrc3d_ugn::RET_SMOOTH},
    {"left_edge",  {-5.0e-3,  5.0e-3,  5.0e-3, 0.0, 0.0, 0.0}, mohrc3d_ugn::RET_LEFT_EDGE},
    {"right_edge", {-5.0e-3, -5.0e-3,  1.0e-2, 0.0, 0.0, 0.0}, mohrc3d_ugn::RET_RIGHT_EDGE},
    {"apex",       { 0.0,     0.0,     5.0e-3, 0.0, 0.0, 0.0}, mohrc3d_ugn::RET_APEX}
  };

  const double eqother[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  bool ok = true;

  for (size_t it=0; it<sizeof(tests)/sizeof(tests[0]); it++)
  {
    vector strain;
    vector eqstatev;
    vector stress;
    vector statev;
    matrix d, dnum, derr;

    fill_vector6(tests[it].strain, strain);
    fill_state_buffer(eqother, eqstatev);
    mm.nlstresses(strain, eqstatev, stress, statev);
    mm.stiffmat(strain, statev, stress, d);
    mm.updateval(statev, eqstatev);
    finite_difference_tangent(mm, tests[it].strain, eqother, dnum);
    diff66(d, dnum, derr);

    const double relerr = norm66(derr)/(norm66(dnum) + 1.0e-14);

    std::printf("[%s]\n", tests[it].name);
    std::printf("  return type : %.0f\n", statev[mohrc3d_ugn::O_RET]);
    std::printf("  stress      : [% .10e, % .10e, % .10e, % .10e, % .10e, % .10e]\n",
                stress[0], stress[1], stress[2], stress[3], stress[4], stress[5]);
    std::printf("  eps_p       : [% .10e, % .10e, % .10e, % .10e, % .10e, % .10e]\n",
                statev[0], statev[1], statev[2], statev[3], statev[4], statev[5]);
    std::printf("  D rel. err. : %.6e\n", relerr);

    if (stress.n != mohrc3d_ugn::NCOMP_STRESS)
    {
      std::printf("  ERROR: unexpected stress size.\n");
      ok = false;
    }

    if (statev.n != mohrc3d_ugn::NCOMP_OTHER)
    {
      std::printf("  ERROR: unexpected statev size.\n");
      ok = false;
    }

    if (static_cast<int>(statev[mohrc3d_ugn::O_RET] + 0.5) != tests[it].expected_return)
    {
      std::printf("  ERROR: unexpected return type.\n");
      ok = false;
    }

    for (long i=0; i<6; i++)
    {
      if (std::fabs(statev[mohrc3d_ugn::O_EP_PREV_XX + i] - eqother[i]) > 1.0e-12)
      {
        std::printf("  ERROR: previous plastic strain was not copied to statev.\n");
        ok = false;
        break;
      }
    }

    for (long i=0; i<mohrc3d_ugn::NCOMP_EQOTHER; i++)
    {
      if (std::fabs(eqstatev[i] - statev[i]) > 1.0e-12)
      {
        std::printf("  ERROR: statev -> eqstatev copy mismatch.\n");
        ok = false;
        break;
      }
    }

    if (relerr > 1.0e-5)
    {
      std::printf("  ERROR: tangent check failed.\n");
      ok = false;
    }

    std::printf("\n");
  }

  if (!ok)
  {
    std::printf("MOHRC3D_UGN TEST STATUS: FAIL\n");
    return 1;
  }

  std::printf("MOHRC3D_UGN TEST STATUS: PASS\n");
  return 0;
}
