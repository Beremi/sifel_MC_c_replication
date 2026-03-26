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

static void fill_buffer_from_ep_prev(const std::vector<double> &ep_prev,
                                     long nt,
                                     long ip,
                                     vector &buf)
{
  reallocv(MATMODEL_NCOMP_EQOTHER, buf);
  nullv(buf);
  for (long i=0; i<4; i++)
    buf[i] = ep_prev[i*nt + ip];
}

static void fill_strain_from_matlab(const std::vector<double> &E_final,
                                    long nt,
                                    long ip,
                                    vector &strain)
{
  reallocv(MATMODEL_NCOMP_STRAIN, strain);
  nullv(strain);
  strain[0] = E_final[ip];
  strain[1] = E_final[nt + ip];
  strain[2] = E_final[2*nt + ip];

  // Matlab constitutive_problem.m receives only the in-plane strain components. For
  // export comparison, keep the SIFEL zz component at zero so the packed payload is
  // compared against the same snapshot.
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
  std::vector<double> eig_matlab;
  std::vector<double> sigma_matlab;
  std::vector<double> return_type_matlab;
  std::vector<double> proj1_matlab;
  std::vector<double> proj2_matlab;
  std::vector<double> proj3_matlab;
  std::vector<double> hess1_matlab;
  std::vector<double> hess2_matlab;
  std::vector<double> hess3_matlab;
  std::vector<double> Sderiv_matlab;
  matmodel mm;
  vector strain;
  vector eqstatev;
  vector statev;
  vector stress;
  matrix d;
  double max_epsp_diff = 0.0;
  double max_return_diff = 0.0;
  double max_eig_diff = 0.0;
  double max_proj_diff = 0.0;
  double max_hess_diff = 0.0;
  double max_sigma_diff = 0.0;
  double max_stress_diff = 0.0;
  double max_tangent_diff = 0.0;

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
      read_ascii_values((matlab_dir + "/eig_matlab.txt").c_str(), 3*meta.nt, eig_matlab) != 0 ||
      read_ascii_values((matlab_dir + "/sigma_matlab.txt").c_str(), 3*meta.nt, sigma_matlab) != 0 ||
      read_ascii_values((matlab_dir + "/return_type_matlab.txt").c_str(), meta.nt, return_type_matlab) != 0 ||
      read_ascii_values((matlab_dir + "/proj1_matlab.txt").c_str(), 4*meta.nt, proj1_matlab) != 0 ||
      read_ascii_values((matlab_dir + "/proj2_matlab.txt").c_str(), 4*meta.nt, proj2_matlab) != 0 ||
      read_ascii_values((matlab_dir + "/proj3_matlab.txt").c_str(), 4*meta.nt, proj3_matlab) != 0 ||
      read_ascii_values((matlab_dir + "/hess1_matlab_reduced.txt").c_str(), 9*meta.nt, hess1_matlab) != 0 ||
      read_ascii_values((matlab_dir + "/hess2_matlab_reduced.txt").c_str(), 9*meta.nt, hess2_matlab) != 0 ||
      read_ascii_values((matlab_dir + "/hess3_matlab_reduced.txt").c_str(), 9*meta.nt, hess3_matlab) != 0 ||
      read_ascii_values((matlab_dir + "/Sderiv_matlab_reduced.txt").c_str(), 9*meta.nt, Sderiv_matlab) != 0)
  {
    std::fprintf(stderr, "Cannot read Matlab export files from %s\n", matlab_dir.c_str());
    return 1;
  }

  for (long ip=0; ip<meta.nt; ip++)
  {
    fill_strain_from_matlab(E_final, meta.nt, ip, strain);
    fill_buffer_from_ep_prev(Ep_prev, meta.nt, ip, eqstatev);

    mm.nlstresses(strain, eqstatev, stress, statev);

    for (long i=0; i<4; i++)
      max_epsp_diff = update_max_abs(max_epsp_diff, statev[i], Ep_final_matlab[i*meta.nt + ip]);

    max_return_diff = update_max_abs(max_return_diff,
                                     statev[MATMODEL_IO_RETURN_TYPE],
                                     return_type_matlab[ip]);

    max_eig_diff = update_max_abs(max_eig_diff, statev[MATMODEL_IO_EIG_1], eig_matlab[ip]);
    max_eig_diff = update_max_abs(max_eig_diff, statev[MATMODEL_IO_EIG_2], eig_matlab[meta.nt + ip]);
    max_eig_diff = update_max_abs(max_eig_diff, statev[MATMODEL_IO_EIG_3], eig_matlab[2*meta.nt + ip]);

    for (long i=0; i<4; i++)
    {
      max_proj_diff = update_max_abs(max_proj_diff, statev[MATMODEL_IO_PROJ_1 + i], proj1_matlab[i*meta.nt + ip]);
      max_proj_diff = update_max_abs(max_proj_diff, statev[MATMODEL_IO_PROJ_2 + i], proj2_matlab[i*meta.nt + ip]);
      max_proj_diff = update_max_abs(max_proj_diff, statev[MATMODEL_IO_PROJ_3 + i], proj3_matlab[i*meta.nt + ip]);
      max_stress_diff = update_max_abs(max_stress_diff, stress[i], S_matlab[i*meta.nt + ip]);
    }

    for (long i=0; i<9; i++)
    {
      max_hess_diff = update_max_abs(max_hess_diff, statev[MATMODEL_IO_HESS_1 + i], hess1_matlab[i*meta.nt + ip]);
      max_hess_diff = update_max_abs(max_hess_diff, statev[MATMODEL_IO_HESS_2 + i], hess2_matlab[i*meta.nt + ip]);
      max_hess_diff = update_max_abs(max_hess_diff, statev[MATMODEL_IO_HESS_3 + i], hess3_matlab[i*meta.nt + ip]);
    }

    max_sigma_diff = update_max_abs(max_sigma_diff, statev[MATMODEL_IO_SIGMA_1], sigma_matlab[ip]);
    max_sigma_diff = update_max_abs(max_sigma_diff, statev[MATMODEL_IO_SIGMA_2], sigma_matlab[meta.nt + ip]);
    max_sigma_diff = update_max_abs(max_sigma_diff, statev[MATMODEL_IO_SIGMA_3], sigma_matlab[2*meta.nt + ip]);

    copyv(statev, eqstatev);
    mm.stiffmat(strain, eqstatev, stress, d);

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
  std::printf("  nt               : %ld\n", meta.nt);
  std::printf("  max epsp diff    : %.6e\n", max_epsp_diff);
  std::printf("  max return diff  : %.6e\n", max_return_diff);
  std::printf("  max eig diff     : %.6e\n", max_eig_diff);
  std::printf("  max proj diff    : %.6e\n", max_proj_diff);
  std::printf("  max hess diff    : %.6e\n", max_hess_diff);
  std::printf("  max sigma diff   : %.6e\n", max_sigma_diff);
  std::printf("  max stress diff  : %.6e\n", max_stress_diff);
  std::printf("  max tangent diff : %.6e\n", max_tangent_diff);

  if ((max_epsp_diff > tol_buffer) ||
      (max_return_diff > 0.5) ||
      (max_eig_diff > tol_buffer) ||
      (max_proj_diff > tol_buffer) ||
      (max_hess_diff > tol_buffer) ||
      (max_sigma_diff > tol_buffer) ||
      (max_stress_diff > tol_stress) ||
      (max_tangent_diff > tol_tangent))
  {
    std::printf("COMPARE STATUS: FAIL\n");
    return 1;
  }

  std::printf("COMPARE STATUS: PASS\n");
  return 0;
}
