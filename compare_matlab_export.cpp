#include "matmodel.h"
#include "vector.h"
#include "matrix.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

struct ExportMeta
{
  long nt;
  double young;
  double poisson;
  double cohesion;
  double phi;
};

static int read_meta(const char *path, ExportMeta &meta)
{
  FILE *in = std::fopen(path, "r");
  char key[128];
  double value;

  if (in == NULL)
    return 1;

  meta.nt = 0;
  meta.young = 0.0;
  meta.poisson = 0.0;
  meta.cohesion = 0.0;
  meta.phi = 0.0;

  while (std::fscanf(in, "%127s %lf", key, &value) == 2)
  {
    if (std::strcmp(key, "nt") == 0)
      meta.nt = static_cast<long>(value + 0.5);
    else if (std::strcmp(key, "young") == 0)
      meta.young = value;
    else if (std::strcmp(key, "poisson") == 0)
      meta.poisson = value;
    else if (std::strcmp(key, "cohesion") == 0)
      meta.cohesion = value;
    else if (std::strcmp(key, "phi") == 0)
      meta.phi = value;
  }

  std::fclose(in);
  return (meta.nt > 0) ? 0 : 1;
}

static int read_ascii_values(const char *path, long count, std::vector<double> &out)
{
  FILE *in = std::fopen(path, "r");

  if (in == NULL)
    return 1;

  out.resize(count);
  for (long i=0; i<count; i++)
  {
    if (std::fscanf(in, "%lf", &out[i]) != 1)
    {
      std::fclose(in);
      return 1;
    }
  }

  std::fclose(in);
  return 0;
}

static int read_params(matmodel &mm, const ExportMeta &meta)
{
  FILE *in = std::tmpfile();

  if (in == NULL)
    return 1;

  std::fprintf(in, "%.17g %.17g %.17g %.17g\n",
               meta.young, meta.poisson, meta.cohesion, meta.phi);
  std::rewind(in);
  const long ret = mm.read(in);
  std::fclose(in);
  return static_cast<int>(ret);
}

static void fill_eqstatev_from_ep_prev(const std::vector<double> &ep_prev,
                                       long nt,
                                       long ip,
                                       vector &buf)
{
  reallocv(matmodel::NCOMP_EQOTHER, buf);
  nullv(buf);
  for (long i=0; i<4; i++)
    buf[i] = ep_prev[i*nt + ip];
}

static void fill_strain_from_matlab(const std::vector<double> &E_final,
                                    long nt,
                                    long ip,
                                    vector &strain)
{
  reallocv(matmodel::NCOMP_STRAIN, strain);
  nullv(strain);
  strain[0] = E_final[ip];
  strain[1] = E_final[nt + ip];
  strain[2] = E_final[2*nt + ip];

  // Matlab constitutive_problem.m receives only the in-plane strain components.
  strain[3] = 0.0;
}

static double update_max_abs(double current, double a, double b)
{
  const double diff = std::fabs(a - b);
  return (diff > current) ? diff : current;
}

int main(int argc, char **argv)
{
  const std::string root = (argc > 1) ? argv[1] : "replication_output";
  const std::string matlab_dir = root + "/matlab";
  const std::string meta_path = root + "/meta.txt";
  const double tol_buffer = 1.0e-10;
  const double tol_stress = 1.0e-10;
  const double tol_tangent = 1.0e-8;

  ExportMeta meta;
  std::vector<double> E_final;
  std::vector<double> Ep_prev;
  std::vector<double> Ep_final_matlab;
  std::vector<double> S_matlab;
  std::vector<double> return_type_matlab;
  std::vector<double> Sderiv_matlab;
  matmodel mm;
  vector strain;
  vector eqstatev;
  vector statev;
  vector eqstatev_new;
  vector stress;
  matrix d;
  double max_epsp_diff = 0.0;
  double max_epsp_prev_copy_diff = 0.0;
  double max_return_diff = 0.0;
  double max_stress_diff = 0.0;
  double max_tangent_diff = 0.0;
  double max_update_copy_diff = 0.0;

  if (read_meta(meta_path.c_str(), meta) != 0)
  {
    std::fprintf(stderr, "Cannot read metadata from %s\n", meta_path.c_str());
    return 1;
  }

  if (read_params(mm, meta) != 0)
  {
    std::fprintf(stderr, "Invalid material parameters in export metadata.\n");
    return 1;
  }

  if (read_ascii_values((matlab_dir + "/E_final.txt").c_str(), 3*meta.nt, E_final) != 0 ||
      read_ascii_values((matlab_dir + "/Ep_prev.txt").c_str(), 4*meta.nt, Ep_prev) != 0 ||
      read_ascii_values((matlab_dir + "/Ep_final_matlab.txt").c_str(), 4*meta.nt, Ep_final_matlab) != 0 ||
      read_ascii_values((matlab_dir + "/S_matlab.txt").c_str(), 4*meta.nt, S_matlab) != 0 ||
      read_ascii_values((matlab_dir + "/return_type_matlab.txt").c_str(), meta.nt, return_type_matlab) != 0 ||
      read_ascii_values((matlab_dir + "/Sderiv_matlab_reduced.txt").c_str(), 9*meta.nt, Sderiv_matlab) != 0)
  {
    std::fprintf(stderr, "Cannot read Matlab export files from %s\n", matlab_dir.c_str());
    return 1;
  }

  for (long ip=0; ip<meta.nt; ip++)
  {
    fill_strain_from_matlab(E_final, meta.nt, ip, strain);
    fill_eqstatev_from_ep_prev(Ep_prev, meta.nt, ip, eqstatev);

    mm.nlstresses(strain, eqstatev, stress, statev);
    mm.stiffmat(strain, statev, stress, d);
    mm.updateval(statev, eqstatev_new);

    for (long i=0; i<4; i++)
    {
      max_epsp_diff = update_max_abs(max_epsp_diff, statev[i], Ep_final_matlab[i*meta.nt + ip]);
      max_epsp_prev_copy_diff = update_max_abs(max_epsp_prev_copy_diff,
                                               statev[matmodel::O_EP_PREV_XX + i],
                                               Ep_prev[i*meta.nt + ip]);
      max_stress_diff = update_max_abs(max_stress_diff, stress[i], S_matlab[i*meta.nt + ip]);
      max_update_copy_diff = update_max_abs(max_update_copy_diff, eqstatev_new[i], statev[i]);
    }

    max_return_diff = update_max_abs(max_return_diff,
                                     statev[matmodel::O_RET],
                                     return_type_matlab[ip]);

    for (long j=0; j<3; j++)
    {
      for (long i=0; i<3; i++)
      {
        const long idx = j*3 + i;
        max_tangent_diff = update_max_abs(max_tangent_diff, d(i,j), Sderiv_matlab[idx*meta.nt + ip]);
      }
    }
  }

  std::printf("Matlab export comparison\n");
  std::printf("  nt                    : %ld\n", meta.nt);
  std::printf("  max epsp diff         : %.6e\n", max_epsp_diff);
  std::printf("  max epsp_prev copy    : %.6e\n", max_epsp_prev_copy_diff);
  std::printf("  max return diff       : %.6e\n", max_return_diff);
  std::printf("  max stress diff       : %.6e\n", max_stress_diff);
  std::printf("  max tangent diff      : %.6e\n", max_tangent_diff);
  std::printf("  max update copy diff  : %.6e\n", max_update_copy_diff);

  if ((max_epsp_diff > tol_buffer) ||
      (max_epsp_prev_copy_diff > tol_buffer) ||
      (max_return_diff > 0.5) ||
      (max_stress_diff > tol_stress) ||
      (max_tangent_diff > tol_tangent) ||
      (max_update_copy_diff > tol_buffer))
  {
    std::printf("COMPARE STATUS: FAIL\n");
    return 1;
  }

  std::printf("COMPARE STATUS: PASS\n");
  return 0;
}
