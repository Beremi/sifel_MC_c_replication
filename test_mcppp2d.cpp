#include "mcppp2d_core.h"

#include <cmath>
#include <cstdio>
#include <cstring>

namespace {

struct TestCase
{
  const char *name;
  double strain[4];
  int expected_return;
};

void copy4(const double src[4], double dst[4])
{
  for (int i=0; i<4; i++)
    dst[i] = src[i];
}

double norm44(const double a[4][4])
{
  double s = 0.0;
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      s += a[i][j]*a[i][j];
  return std::sqrt(s);
}

void diff44(const double a[4][4], const double b[4][4], double out[4][4])
{
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      out[i][j] = a[i][j] - b[i][j];
}

void finite_difference_tangent(const mcppp2d_params &par,
                               const double strain[4],
                               const double eqother[4],
                               double dnum[4][4])
{
  const double h = 1.0e-7;
  for (int j=0; j<4; j++)
  {
    double ep[4], em[4];
    mcppp2d_cache cp, cm;

    copy4(strain, ep);
    copy4(strain, em);
    ep[j] += h;
    em[j] -= h;

    mcppp2d_compute_response(&par, ep, eqother, &cp);
    mcppp2d_compute_response(&par, em, eqother, &cm);

    for (int i=0; i<4; i++)
      dnum[i][j] = (cp.stress[i] - cm.stress[i])/(2.0*h);
  }
}

} // namespace

int main()
{
  mcppp2d_params par;
  mcppp2d_init_params(&par, 20000.0, 0.49, 50.0, M_PI/9.0);

  if (!mcppp2d_validate_params(&par))
  {
    std::fprintf(stderr, "Invalid material parameters.\n");
    return 1;
  }

  const TestCase tests[] = {
    {"elastic",    { 1.0e-5,  0.0,    0.0,   0.0  }, MCPPP2D_RETURN_ELASTIC},
    {"smooth",     {-5.0e-3,  0.0,    0.0,   5.0e-3}, MCPPP2D_RETURN_SMOOTH},
    {"left_edge",  {-5.0e-3,  5.0e-3, 0.0,   5.0e-3}, MCPPP2D_RETURN_LEFT_EDGE},
    {"right_edge", {-5.0e-3, -5.0e-3, 0.0,   1.0e-2}, MCPPP2D_RETURN_RIGHT_EDGE},
    {"apex",       { 0.0,     0.0,    0.0,   5.0e-3}, MCPPP2D_RETURN_APEX}
  };

  const double eqother[4] = {0.0, 0.0, 0.0, 0.0};
  bool ok = true;

  for (size_t it=0; it<sizeof(tests)/sizeof(tests[0]); it++)
  {
    mcppp2d_cache cache, decoded;
    double d[4][4], dnum[4][4], derr[4][4];
    double other[71], eqnext[4];
    double relerr;

    mcppp2d_compute_response(&par, tests[it].strain, eqother, &cache);
    mcppp2d_compute_tangent(&par, &cache, d);
    finite_difference_tangent(par, tests[it].strain, eqother, dnum);

    mcppp2d_pack_other(&cache, other);
    mcppp2d_unpack_other(other, &decoded);
    mcppp2d_pack_eqother(&cache, eqnext);

    diff44(d, dnum, derr);
    relerr = norm44(derr)/(norm44(dnum) + 1.0e-14);

    std::printf("[%s]\n", tests[it].name);
    std::printf("  return type : %s\n", mcppp2d_return_name(cache.return_type));
    std::printf("  stress      : [% .10e, % .10e, % .10e, % .10e]\n",
                cache.stress[0], cache.stress[1], cache.stress[2], cache.stress[3]);
    std::printf("  eps_p       : [% .10e, % .10e, % .10e, % .10e]\n",
                cache.epsp[0], cache.epsp[1], cache.epsp[2], cache.epsp[3]);
    std::printf("  D rel. err. : %.6e\n", relerr);

    if (cache.return_type != tests[it].expected_return)
    {
      std::printf("  ERROR: expected return type %s\n",
                  mcppp2d_return_name(tests[it].expected_return));
      ok = false;
    }

    if (relerr > 1.0e-6)
    {
      std::printf("  ERROR: tangent check failed.\n");
      ok = false;
    }

    for (int i=0; i<4; i++)
    {
      if (std::fabs(cache.epsp[i] - decoded.epsp[i]) > 1.0e-14)
      {
        std::printf("  ERROR: other-buffer roundtrip failed at epsp[%d].\n", i);
        ok = false;
      }
      if (std::fabs(cache.epsp[i] - eqnext[i]) > 1.0e-14)
      {
        std::printf("  ERROR: eqother-buffer packing failed at epsp[%d].\n", i);
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
