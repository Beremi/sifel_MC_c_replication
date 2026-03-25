#include "mcppp2d_core.h"

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

namespace {

struct Matrix
{
  Matrix() : rows(0), cols(0) {}
  Matrix(size_t r, size_t c) : rows(r), cols(c), data(r*c, 0.0) {}

  double &operator()(size_t r, size_t c)
  {
    return data[r*cols + c];
  }

  double operator()(size_t r, size_t c) const
  {
    return data[r*cols + c];
  }

  size_t rows;
  size_t cols;
  std::vector<double> data;
};

struct Stats
{
  Stats()
  : count(0),
    max_abs(0.0),
    sum_abs(0.0),
    sum_sq_diff(0.0),
    sum_sq_ref(0.0),
    worst_flat(0),
    worst_cpp(0.0),
    worst_ref(0.0)
  {}

  size_t count;
  double max_abs;
  double sum_abs;
  double sum_sq_diff;
  double sum_sq_ref;
  size_t worst_flat;
  double worst_cpp;
  double worst_ref;
};

struct QuantitySummary
{
  std::string name;
  Stats stats;
  size_t rows;
  size_t cols;
};

std::string join_path(const std::string &lhs, const std::string &rhs)
{
  if (lhs.empty())
    return rhs;
  if (lhs[lhs.size()-1] == '/')
    return lhs + rhs;
  return lhs + "/" + rhs;
}

void ensure_directory(const std::string &path)
{
  if (path.empty())
    return;

  if (::mkdir(path.c_str(), 0777) == 0)
    return;

  if (errno == EEXIST)
    return;

  throw std::runtime_error("Failed to create directory: " + path);
}

std::map<std::string, double> load_metadata(const std::string &path)
{
  std::ifstream in(path.c_str());
  std::map<std::string, double> meta;
  std::string key;
  double value = 0.0;

  if (!in)
    throw std::runtime_error("Cannot open metadata file: " + path);

  while (in >> key >> value)
    meta[key] = value;

  return meta;
}

double require_meta(const std::map<std::string, double> &meta, const std::string &key)
{
  std::map<std::string, double>::const_iterator it = meta.find(key);
  if (it == meta.end())
    throw std::runtime_error("Missing metadata key: " + key);
  return it->second;
}

Matrix load_matrix(const std::string &path, size_t rows, size_t cols)
{
  std::ifstream in(path.c_str());
  Matrix out(rows, cols);
  size_t idx = 0;

  if (!in)
    throw std::runtime_error("Cannot open matrix file: " + path);

  while ((idx < out.data.size()) && (in >> out.data[idx]))
    idx++;

  if (idx != out.data.size())
  {
    std::ostringstream oss;
    oss << "Unexpected matrix size in " << path
        << ": expected " << out.data.size()
        << " values, got " << idx;
    throw std::runtime_error(oss.str());
  }

  double dummy = 0.0;
  if (in >> dummy)
    throw std::runtime_error("Matrix file contains extra data: " + path);

  return out;
}

void write_matrix(const std::string &path, const Matrix &m)
{
  std::ofstream out(path.c_str());

  if (!out)
    throw std::runtime_error("Cannot write matrix file: " + path);

  out << std::scientific << std::setprecision(17);
  for (size_t i=0; i<m.rows; i++)
  {
    for (size_t j=0; j<m.cols; j++)
    {
      if (j > 0)
        out << ' ';
      out << m(i, j);
    }
    out << '\n';
  }
}

void update_stats(Stats &stats, double cpp_value, double ref_value, size_t flat_index)
{
  const double diff = cpp_value - ref_value;
  const double abs_diff = std::fabs(diff);

  stats.count++;
  stats.sum_abs += abs_diff;
  stats.sum_sq_diff += diff*diff;
  stats.sum_sq_ref += ref_value*ref_value;

  if (abs_diff > stats.max_abs)
  {
    stats.max_abs = abs_diff;
    stats.worst_flat = flat_index;
    stats.worst_cpp = cpp_value;
    stats.worst_ref = ref_value;
  }
}

Stats compare_matrices(const Matrix &cpp_matrix, const Matrix &ref_matrix)
{
  Stats stats;

  if ((cpp_matrix.rows != ref_matrix.rows) || (cpp_matrix.cols != ref_matrix.cols))
    throw std::runtime_error("Matrix size mismatch in comparison.");

  for (size_t idx=0; idx<cpp_matrix.data.size(); idx++)
    update_stats(stats, cpp_matrix.data[idx], ref_matrix.data[idx], idx);

  return stats;
}

double rel_l2(const Stats &stats)
{
  return std::sqrt(stats.sum_sq_diff)/(std::sqrt(stats.sum_sq_ref) + 1.0e-30);
}

double rms_abs(const Stats &stats)
{
  if (stats.count == 0)
    return 0.0;
  return std::sqrt(stats.sum_sq_diff/static_cast<double>(stats.count));
}

double mean_abs(const Stats &stats)
{
  if (stats.count == 0)
    return 0.0;
  return stats.sum_abs/static_cast<double>(stats.count);
}

void write_stats_row(std::ofstream &out, const QuantitySummary &summary)
{
  const size_t worst_row = summary.stats.worst_flat/summary.cols;
  const size_t worst_col = summary.stats.worst_flat%summary.cols;

  out << "| " << summary.name
      << " | " << summary.rows << "x" << summary.cols
      << " | " << std::scientific << std::setprecision(6) << summary.stats.max_abs
      << " | " << mean_abs(summary.stats)
      << " | " << rms_abs(summary.stats)
      << " | " << rel_l2(summary.stats)
      << " | row " << (worst_row + 1) << ", elem " << (worst_col + 1)
      << " |\n";
}

void fill_hess_reduced_column_major(const double hess[4][4], double out[9])
{
  size_t idx = 0;
  for (size_t col=0; col<3; col++)
  {
    for (size_t row=0; row<3; row++)
      out[idx++] = hess[row][col];
  }
}

void fill_full_column_major(const double a[4][4], double out[16])
{
  size_t idx = 0;
  for (size_t col=0; col<4; col++)
  {
    for (size_t row=0; row<4; row++)
      out[idx++] = a[row][col];
  }
}

const char *return_type_name(int value)
{
  switch (value)
  {
    case MCPPP2D_RETURN_ELASTIC:    return "elastic";
    case MCPPP2D_RETURN_SMOOTH:     return "smooth";
    case MCPPP2D_RETURN_LEFT_EDGE:  return "left_edge";
    case MCPPP2D_RETURN_RIGHT_EDGE: return "right_edge";
    case MCPPP2D_RETURN_APEX:       return "apex";
    default:                        return "unknown";
  }
}

} // namespace

int main(int argc, char **argv)
{
  try
  {
    const std::string root = (argc > 1) ? argv[1] : "replication_output";
    const std::string matlab_dir = join_path(root, "matlab");
    const std::string cpp_dir = join_path(root, "cpp");

    ensure_directory(root);
    ensure_directory(cpp_dir);

    const std::map<std::string, double> meta = load_metadata(join_path(root, "meta.txt"));
    const size_t nt = static_cast<size_t>(require_meta(meta, "nt") + 0.5);
    const size_t nv = static_cast<size_t>(require_meta(meta, "nv") + 0.5);

    Matrix E_final = load_matrix(join_path(matlab_dir, "E_final.txt"), 3, nt);
    Matrix Ep_prev = load_matrix(join_path(matlab_dir, "Ep_prev.txt"), 4, nt);
    Matrix Ep_final_matlab = load_matrix(join_path(matlab_dir, "Ep_final_matlab.txt"), 4, nt);
    Matrix S_matlab = load_matrix(join_path(matlab_dir, "S_matlab.txt"), 4, nt);
    Matrix eig_matlab = load_matrix(join_path(matlab_dir, "eig_matlab.txt"), 3, nt);
    Matrix sigma_matlab = load_matrix(join_path(matlab_dir, "sigma_matlab.txt"), 3, nt);
    Matrix return_type_matlab = load_matrix(join_path(matlab_dir, "return_type_matlab.txt"), 1, nt);
    Matrix proj1_matlab = load_matrix(join_path(matlab_dir, "proj1_matlab.txt"), 4, nt);
    Matrix proj2_matlab = load_matrix(join_path(matlab_dir, "proj2_matlab.txt"), 4, nt);
    Matrix proj3_matlab = load_matrix(join_path(matlab_dir, "proj3_matlab.txt"), 4, nt);
    Matrix hess1_matlab = load_matrix(join_path(matlab_dir, "hess1_matlab_reduced.txt"), 9, nt);
    Matrix hess2_matlab = load_matrix(join_path(matlab_dir, "hess2_matlab_reduced.txt"), 9, nt);
    Matrix hess3_matlab = load_matrix(join_path(matlab_dir, "hess3_matlab_reduced.txt"), 9, nt);
    Matrix Sderiv_matlab = load_matrix(join_path(matlab_dir, "Sderiv_matlab_reduced.txt"), 9, nt);

    Matrix S_cpp(4, nt);
    Matrix Ep_final_cpp(4, nt);
    Matrix eig_cpp(3, nt);
    Matrix sigma_cpp(3, nt);
    Matrix return_type_cpp(1, nt);
    Matrix proj1_cpp(4, nt);
    Matrix proj2_cpp(4, nt);
    Matrix proj3_cpp(4, nt);
    Matrix hess1_cpp_reduced(9, nt);
    Matrix hess2_cpp_reduced(9, nt);
    Matrix hess3_cpp_reduced(9, nt);
    Matrix hess1_cpp_full(16, nt);
    Matrix hess2_cpp_full(16, nt);
    Matrix hess3_cpp_full(16, nt);
    Matrix Sderiv_cpp_reduced(9, nt);
    Matrix D_cpp_full(16, nt);

    mcppp2d_params par;
    mcppp2d_init_params(&par,
                        require_meta(meta, "young"),
                        require_meta(meta, "poisson"),
                        require_meta(meta, "cohesion"),
                        require_meta(meta, "phi"));

    if (!mcppp2d_validate_params(&par))
      throw std::runtime_error("Invalid material parameters loaded from metadata.");

    std::vector<size_t> matlab_counts(5, 0);
    std::vector<size_t> cpp_counts(5, 0);
    std::vector<size_t> mismatch_elements;

    for (size_t elem=0; elem<nt; elem++)
    {
      double strain[4] = {E_final(0, elem), E_final(1, elem), E_final(2, elem), 0.0};
      double eqother[4] = {Ep_prev(0, elem), Ep_prev(1, elem), Ep_prev(2, elem), Ep_prev(3, elem)};
      double tangent[4][4];
      double reduced[9];
      double full16[16];
      mcppp2d_cache cache;
      int matlab_type = static_cast<int>(return_type_matlab(0, elem) + 0.5);

      mcppp2d_compute_response(&par, strain, eqother, &cache);
      mcppp2d_compute_tangent(&par, &cache, tangent);

      for (size_t i=0; i<4; i++)
      {
        S_cpp(i, elem) = cache.stress[i];
        Ep_final_cpp(i, elem) = cache.epsp[i];
        proj1_cpp(i, elem) = cache.proj[0][i];
        proj2_cpp(i, elem) = cache.proj[1][i];
        proj3_cpp(i, elem) = cache.proj[2][i];
      }

      for (size_t i=0; i<3; i++)
      {
        eig_cpp(i, elem) = cache.eig[i];
        sigma_cpp(i, elem) = cache.sig_princ[i];
      }

      fill_hess_reduced_column_major(cache.hess[0], reduced);
      for (size_t i=0; i<9; i++)
        hess1_cpp_reduced(i, elem) = reduced[i];
      fill_hess_reduced_column_major(cache.hess[1], reduced);
      for (size_t i=0; i<9; i++)
        hess2_cpp_reduced(i, elem) = reduced[i];
      fill_hess_reduced_column_major(cache.hess[2], reduced);
      for (size_t i=0; i<9; i++)
        hess3_cpp_reduced(i, elem) = reduced[i];

      fill_full_column_major(cache.hess[0], full16);
      for (size_t i=0; i<16; i++)
        hess1_cpp_full(i, elem) = full16[i];
      fill_full_column_major(cache.hess[1], full16);
      for (size_t i=0; i<16; i++)
        hess2_cpp_full(i, elem) = full16[i];
      fill_full_column_major(cache.hess[2], full16);
      for (size_t i=0; i<16; i++)
        hess3_cpp_full(i, elem) = full16[i];

      fill_hess_reduced_column_major(tangent, reduced);
      for (size_t i=0; i<9; i++)
        Sderiv_cpp_reduced(i, elem) = reduced[i];

      fill_full_column_major(tangent, full16);
      for (size_t i=0; i<16; i++)
        D_cpp_full(i, elem) = full16[i];

      return_type_cpp(0, elem) = static_cast<double>(cache.return_type);

      if ((matlab_type >= 0) && (matlab_type < 5))
        matlab_counts[matlab_type]++;
      if ((cache.return_type >= 0) && (cache.return_type < 5))
        cpp_counts[cache.return_type]++;
      if (cache.return_type != matlab_type)
        mismatch_elements.push_back(elem + 1);
    }

    write_matrix(join_path(cpp_dir, "S_cpp.txt"), S_cpp);
    write_matrix(join_path(cpp_dir, "Ep_final_cpp.txt"), Ep_final_cpp);
    write_matrix(join_path(cpp_dir, "eig_cpp.txt"), eig_cpp);
    write_matrix(join_path(cpp_dir, "sigma_cpp.txt"), sigma_cpp);
    write_matrix(join_path(cpp_dir, "return_type_cpp.txt"), return_type_cpp);
    write_matrix(join_path(cpp_dir, "proj1_cpp.txt"), proj1_cpp);
    write_matrix(join_path(cpp_dir, "proj2_cpp.txt"), proj2_cpp);
    write_matrix(join_path(cpp_dir, "proj3_cpp.txt"), proj3_cpp);
    write_matrix(join_path(cpp_dir, "hess1_cpp_reduced.txt"), hess1_cpp_reduced);
    write_matrix(join_path(cpp_dir, "hess2_cpp_reduced.txt"), hess2_cpp_reduced);
    write_matrix(join_path(cpp_dir, "hess3_cpp_reduced.txt"), hess3_cpp_reduced);
    write_matrix(join_path(cpp_dir, "hess1_cpp_full.txt"), hess1_cpp_full);
    write_matrix(join_path(cpp_dir, "hess2_cpp_full.txt"), hess2_cpp_full);
    write_matrix(join_path(cpp_dir, "hess3_cpp_full.txt"), hess3_cpp_full);
    write_matrix(join_path(cpp_dir, "Sderiv_cpp_reduced.txt"), Sderiv_cpp_reduced);
    write_matrix(join_path(cpp_dir, "D_cpp_full.txt"), D_cpp_full);

    std::vector<QuantitySummary> summaries;
    summaries.push_back(QuantitySummary{"stress S", compare_matrices(S_cpp, S_matlab), 4, nt});
    summaries.push_back(QuantitySummary{"plastic strain Ep", compare_matrices(Ep_final_cpp, Ep_final_matlab), 4, nt});
    summaries.push_back(QuantitySummary{"ordered eig", compare_matrices(eig_cpp, eig_matlab), 3, nt});
    summaries.push_back(QuantitySummary{"principal stress sigma", compare_matrices(sigma_cpp, sigma_matlab), 3, nt});
    summaries.push_back(QuantitySummary{"projection Eig_1", compare_matrices(proj1_cpp, proj1_matlab), 4, nt});
    summaries.push_back(QuantitySummary{"projection Eig_2", compare_matrices(proj2_cpp, proj2_matlab), 4, nt});
    summaries.push_back(QuantitySummary{"projection Eig_3", compare_matrices(proj3_cpp, proj3_matlab), 4, nt});
    summaries.push_back(QuantitySummary{"Hessian EIG_1 reduced", compare_matrices(hess1_cpp_reduced, hess1_matlab), 9, nt});
    summaries.push_back(QuantitySummary{"Hessian EIG_2 reduced", compare_matrices(hess2_cpp_reduced, hess2_matlab), 9, nt});
    summaries.push_back(QuantitySummary{"Hessian EIG_3 reduced", compare_matrices(hess3_cpp_reduced, hess3_matlab), 9, nt});
    summaries.push_back(QuantitySummary{"tangent Sderiv reduced", compare_matrices(Sderiv_cpp_reduced, Sderiv_matlab), 9, nt});

    std::ofstream report(join_path(root, "report.md").c_str());
    if (!report)
      throw std::runtime_error("Cannot write report file.");

    report << "# Matlab to C++ Replication Report\n\n";
    report << "## Run Summary\n\n";
    report << "- Final accepted load-step count: " << static_cast<long>(require_meta(meta, "accepted_steps") + 0.5) << "\n";
    report << "- Final printed step index: " << static_cast<long>(require_meta(meta, "final_printed_step") + 0.5) << "\n";
    report << "- Final gravity load factor zeta: " << std::scientific << std::setprecision(17) << require_meta(meta, "final_zeta") << "\n";
    report << "- Final settlement at point A: " << require_meta(meta, "final_settlement") << "\n";
    report << "- Final Newton iteration count: " << static_cast<long>(require_meta(meta, "last_iter") + 0.5) << "\n";
    report << "- Mesh size: nv = " << nv << ", nt = " << nt << "\n";
    report << "- Exit reason code: " << static_cast<long>(require_meta(meta, "exit_reason_code") + 0.5)
           << " (1 = too_large_settlement, 2 = too_small_load_increment, 3 = max_steps)\n\n";

    report << "## Return Types\n\n";
    report << "| Return type | Matlab count | C++ count |\n";
    report << "| --- | ---: | ---: |\n";
    for (int value=0; value<5; value++)
    {
      report << "| " << return_type_name(value)
             << " | " << matlab_counts[value]
             << " | " << cpp_counts[value]
             << " |\n";
    }
    report << "\n";
    report << "- Return-type mismatches: " << mismatch_elements.size() << "\n";
    if (!mismatch_elements.empty())
    {
      report << "- Mismatching element indices:";
      for (size_t i=0; i<mismatch_elements.size(); i++)
      {
        report << ' ' << mismatch_elements[i];
        if (i + 1 >= 20)
        {
          report << " ...";
          break;
        }
      }
      report << "\n";
    }
    report << "\n";

    report << "## Quantity Comparison\n\n";
    report << "| Quantity | Shape | Max abs diff | Mean abs diff | RMS abs diff | Rel L2 diff | Worst location |\n";
    report << "| --- | --- | ---: | ---: | ---: | ---: | --- |\n";
    for (size_t i=0; i<summaries.size(); i++)
      write_stats_row(report, summaries[i]);
    report << "\n";

    report << "## Exported Data\n\n";
    report << "- Matlab inputs and reference outputs: `" << matlab_dir << "`\n";
    report << "- C++ replay outputs: `" << cpp_dir << "`\n";
    report << "- Full 4x4 C++ Hessians are exported as `hess*_cpp_full.txt`.\n";
    report << "- Full 4x4 C++ consistent tangent is exported as `D_cpp_full.txt`.\n";
    report << "- Matlab only provides reduced 3x3 Hessians and reduced 3x3 tangent, so the numerical comparison uses those reduced blocks.\n";

    report.close();

    std::cout << "Replication report written to " << join_path(root, "report.md") << "\n";
    std::cout << "Return-type mismatches: " << mismatch_elements.size() << "\n";
    for (size_t i=0; i<summaries.size(); i++)
    {
      std::cout << summaries[i].name
                << ": max_abs=" << std::scientific << std::setprecision(6) << summaries[i].stats.max_abs
                << ", rel_l2=" << rel_l2(summaries[i].stats)
                << "\n";
    }
  }
  catch (const std::exception &ex)
  {
    std::cerr << "ERROR: " << ex.what() << "\n";
    return 1;
  }

  return 0;
}
