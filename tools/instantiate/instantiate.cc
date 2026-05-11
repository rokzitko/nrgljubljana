// Initial data instantiation tool
// First production slice: in-process Wilson-chain generation.

#include <cstdlib>
#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <sstream>
#include <system_error>
#include <utility>
#include <vector>

#include <diag.hpp>
#include <tridiag.hpp>
#include <read-input.hpp>
#include <traits.hpp>

#include "matrix_evaluator.hpp"
#include "nrgchain.hpp"

namespace {

struct Options {
  std::string param_filename = "param";
  std::string template_dir = "template";
  bool wilson_only = false;
  bool diag_seed_only = false;
  bool generate_temporaries = false;
};

void usage(std::ostream &out = std::cout) {
  out << "Usage: instantiate [--wilson-only | --diag-seed-only] [--generate-temporaries] [--param FILE] [--template-dir DIR]\n";
}

Options parse_options(const int argc, char *argv[]) {
  Options options;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "-h" || arg == "--help") {
      usage();
      std::exit(EXIT_SUCCESS);
    }
    if (arg == "--wilson-only") {
      options.wilson_only = true;
      continue;
    }
    if (arg == "--diag-seed-only") {
      options.diag_seed_only = true;
      continue;
    }
    if (arg == "--generate-temporaries") {
      options.generate_temporaries = true;
      continue;
    }
    if (arg == "--param") {
      if (i + 1 >= argc) throw std::invalid_argument("--param requires a filename.");
      options.param_filename = argv[++i];
      continue;
    }
    if (arg == "--template-dir") {
      if (i + 1 >= argc) throw std::invalid_argument("--template-dir requires a directory.");
      options.template_dir = argv[++i];
      continue;
    }
    throw std::invalid_argument("Unknown argument: " + arg);
  }
  return options;
}

std::string trim(std::string value) {
  const auto first = std::find_if_not(value.begin(), value.end(), [](const unsigned char ch) { return std::isspace(ch); });
  const auto last = std::find_if_not(value.rbegin(), value.rend(), [](const unsigned char ch) { return std::isspace(ch); }).base();
  if (first >= last) return {};
  return std::string(first, last);
}

std::string read_text_file(const std::filesystem::path &filename) {
  std::ifstream in(filename);
  if (!in) throw std::runtime_error("Can't open " + filename.string() + " for reading.");
  std::ostringstream text;
  text << in.rdbuf();
  return text.str();
}

std::filesystem::path resolve_template_file(const std::filesystem::path &template_dir, const std::filesystem::path &filename) {
  if (filename.is_absolute()) return filename;
  const auto in_template_dir = template_dir / filename;

  std::error_code ec;
  if (std::filesystem::exists(in_template_dir, ec)) return in_template_dir;
  ec.clear();
  if (std::filesystem::exists(filename, ec)) return filename;
  return in_template_dir;
}

struct ParamSections {
  std::map<std::string, std::vector<std::string>> lines;
  std::map<std::string, std::map<std::string, std::string>> values;
};

ParamSections parse_param_sections(const std::string &filename) {
  std::ifstream in(filename);
  if (!in) throw std::runtime_error("Can't open " + filename + " for reading.");

  ParamSections sections;
  std::string section;
  std::string line;
  while (std::getline(in, line)) {
    const auto stripped = trim(line);
    if (stripped.size() >= 2 && stripped.front() == '[' && stripped.back() == ']') {
      section = stripped.substr(1, stripped.size() - 2);
    }
    if (section.empty()) continue;
    sections.lines[section].push_back(line);

    if (stripped.empty() || stripped.front() == '#' || stripped.front() == '[') continue;
    const auto eq = stripped.find('=');
    if (eq == std::string::npos) continue;
    sections.values[section][trim(stripped.substr(0, eq))] = trim(stripped.substr(eq + 1));
  }
  return sections;
}

void write_param_section_files(const ParamSections &sections, const std::string &param_filename) {
  for (const auto &[section, lines] : sections.lines) {
    std::ofstream out(param_filename + "." + section);
    if (!out) throw std::runtime_error("Can't open " + param_filename + "." + section + " for writing.");
    for (const auto &line : lines) out << line << '\n';
  }
}

void remove_file_noexcept(const std::filesystem::path &filename) {
  std::error_code ec;
  std::filesystem::remove(filename, ec);
}

void remove_prefixed_files_noexcept(const std::string &prefix) {
  std::error_code ec;
  for (std::filesystem::directory_iterator it(".", ec), end; !ec && it != end; it.increment(ec)) {
    const auto filename = it->path().filename().string();
    if (filename.rfind(prefix + ".", 0) != 0) continue;

    std::error_code type_ec;
    if (it->is_directory(type_ec)) continue;
    remove_file_noexcept(it->path());
  }
}

void cleanup_instantiation_temporaries(const std::string &param_filename, const ParamSections &sections) {
  const std::vector<std::string> filenames = {"xi.dat",      "zeta.dat",     "theta.dat", "xi1.dat",     "zeta1.dat",
                                             "theta1.dat",  "val",          "vec",       "mat",         "mat.res",
                                             "ham",         "param.param",  "param.dmft", "param.extra", "solverlog",
                                             "solverlogneg"};
  for (const auto &filename : filenames) remove_file_noexcept(filename);
  for (const auto &[section, lines] : sections.lines) {
    (void)lines;
    remove_file_noexcept(param_filename + "." + section);
  }
  remove_prefixed_files_noexcept("vec");
  remove_prefixed_files_noexcept("ham");
}

class TemporaryCleanupGuard {
 private:
  bool enabled = false;
  std::string param_filename;
  ParamSections sections;

 public:
  TemporaryCleanupGuard(const bool enabled_, std::string param_filename_, ParamSections sections_)
      : enabled(enabled_), param_filename(std::move(param_filename_)), sections(std::move(sections_)) {}
  TemporaryCleanupGuard(const TemporaryCleanupGuard &) = delete;
  TemporaryCleanupGuard &operator=(const TemporaryCleanupGuard &) = delete;
  ~TemporaryCleanupGuard() {
    if (enabled) cleanup_instantiation_temporaries(param_filename, sections);
  }
};

std::map<std::string, double> numeric_variables(const std::map<std::string, std::string> &values) {
  std::map<std::string, double> variables;
  for (const auto &[key, value] : values) {
    char *end = nullptr;
    const auto parsed = std::strtod(value.c_str(), &end);
    if (end != value.c_str() && trim(end).empty()) variables[key] = parsed;
  }
  return variables;
}

void write_values(const std::string &filename, const std::vector<double> &values) {
  std::ofstream out(filename);
  if (!out) throw std::runtime_error("Can't open " + filename + " for writing.");
  out << std::setprecision(16);
  for (const auto value : values) out << value << '\n';
}

void write_scalar(const std::string &filename, const double value) {
  std::ofstream out(filename);
  if (!out) throw std::runtime_error("Can't open " + filename + " for writing.");
  out << std::setprecision(16) << value << '\n';
}

void write_wilson_channel(const NRG::Tools::NrgChain::WilsonChannel &channel, const size_t index) {
  const auto suffix = std::to_string(index) + ".dat";
  write_values("xi" + suffix, channel.xi);
  write_values("zeta" + suffix, channel.zeta);
  write_scalar("theta" + suffix, channel.theta);
}

void save_text_matrix(const std::string &filename, const NRG::Matrix_traits<double> &matrix) {
  std::ofstream out(filename);
  if (!out) throw std::runtime_error("Can't open " + filename + " for writing.");
  out << std::setprecision(18);
  for (Eigen::Index row = 0; row < matrix.rows(); ++row) {
    for (Eigen::Index col = 0; col < matrix.cols(); ++col) {
      if (col != 0) out << ' ';
      out << matrix(row, col);
    }
    out << '\n';
  }
}

void save_binary_matrix(const std::string &filename, const NRG::Matrix_traits<double> &matrix) {
  if (matrix.rows() > static_cast<Eigen::Index>(std::numeric_limits<unsigned int>::max()) ||
      matrix.cols() > static_cast<Eigen::Index>(std::numeric_limits<unsigned int>::max()))
    throw std::runtime_error("Matrix too large for legacy binary format.");
  std::ofstream out(filename, std::ios::binary);
  if (!out) throw std::runtime_error("Can't open " + filename + " for writing.");
  const auto rows = static_cast<unsigned int>(matrix.rows());
  const auto cols = static_cast<unsigned int>(matrix.cols());
  out.write(reinterpret_cast<const char *>(&rows), sizeof(rows));
  out.write(reinterpret_cast<const char *>(&cols), sizeof(cols));
  for (Eigen::Index row = 0; row < matrix.rows(); ++row)
    for (Eigen::Index col = 0; col < matrix.cols(); ++col) {
      const auto value = matrix(row, col);
      out.write(reinterpret_cast<const char *>(&value), sizeof(value));
    }
}

void save_values(const std::string &filename, const std::vector<double> &values) {
  std::ofstream out(filename);
  if (!out) throw std::runtime_error("Can't open " + filename + " for writing.");
  out << std::setprecision(18);
  for (const auto value : values) out << value << '\n';
}

std::vector<std::string> fields(const std::string &line) {
  std::istringstream in(line);
  std::vector<std::string> result;
  std::string field;
  while (in >> field) result.push_back(field);
  return result;
}

class DataTemplateReader {
 private:
  std::ifstream in;
  std::ostream *comments_out = nullptr;
  double scale_factor = -1.0;

  static std::optional<double> scale_from_comment(const std::string &line) {
    std::istringstream input(line);
    std::string hash, keyword;
    double scale = 0.0;
    input >> hash >> keyword >> scale;
    if (hash == "#" && keyword == "SCALE") return 1.0 / scale;
    return std::nullopt;
  }

 public:
  explicit DataTemplateReader(const std::filesystem::path &filename, std::ostream *comments_out_ = nullptr)
      : in(filename), comments_out(comments_out_) {
    if (!in) throw std::runtime_error("Can't open " + filename.string() + " for reading.");
  }

  std::optional<std::string> next_data_line() {
    std::string line;
    while (std::getline(in, line)) {
      const auto stripped = trim(line);
      if (stripped.empty()) continue;
      if (stripped.front() == '#') {
        if (const auto scale = scale_from_comment(stripped)) scale_factor = *scale;
        if (comments_out != nullptr) *comments_out << line << '\n';
        continue;
      }
      return stripped;
    }
    return std::nullopt;
  }

  double factor() const {
    if (scale_factor <= 0.0) throw std::runtime_error("No # SCALE line found in template/data.in.");
    return scale_factor;
  }
};

struct TemplateHeader {
  size_t channels = 0;
  size_t chain_sites = 0;
  size_t subspaces = 0;
};

struct SeedSubspace {
  std::string qn_line;
  NRG::Invar qn;
  size_t dimension = 0;
  std::vector<double> eigenvalues;
  NRG::Matrix_traits<double> eigenvectors;
};

struct SeedData {
  double factor = 1.0;
  double smallest = std::numeric_limits<double>::max();
  double ground_energy = 0.0;
  std::vector<SeedSubspace> subspaces;
  std::map<NRG::Invar, size_t> dimensions;
  std::map<NRG::Invar, NRG::Matrix_traits<double>> eigenvectors;
};

struct MatrixElementData {
  std::string qn_line;
  NRG::Invar bra;
  NRG::Invar ket;
  NRG::Matrix_traits<double> matrix;
};

struct MatrixBlockData {
  std::string header_line;
  char block_type = '\0';
  std::vector<MatrixElementData> elements;
};

struct TemplateTailData {
  bool has_ground_energy = false;
  double ground_energy = 0.0;
  std::vector<MatrixBlockData> matrix_blocks;
};

std::string legacy_qn_name(const NRG::Invar &qn) {
  std::ostringstream out;
  qn.insertor(out, ".");
  return out.str();
}

NRG::Invar parse_invar_line(const std::string &line, const std::string &context) {
  std::istringstream input(line);
  NRG::Invar qn;
  try {
    input >> qn;
  } catch (const std::exception &e) {
    throw std::runtime_error("Failed reading quantum numbers for " + context + ": " + e.what());
  }
  std::string extra;
  if (input >> extra) throw std::runtime_error("Too many quantum-number fields for " + context + ": " + line);
  return qn;
}

std::pair<NRG::Invar, NRG::Invar> parse_invar_pair_line(const std::string &line, const std::string &context) {
  std::istringstream input(line);
  NRG::Invar qn1;
  NRG::Invar qn2;
  try {
    input >> qn1 >> qn2;
  } catch (const std::exception &e) {
    throw std::runtime_error("Failed reading quantum-number pair for " + context + ": " + e.what());
  }
  std::string extra;
  if (input >> extra) throw std::runtime_error("Too many quantum-number fields for " + context + ": " + line);
  return {qn1, qn2};
}

size_t parse_size_t_value(const std::string &text, const std::string &context) {
  const auto stripped = trim(text);
  if (stripped.empty() || stripped.front() == '-' || stripped.front() == '+')
    throw std::runtime_error("Invalid unsigned integer for " + context + ": " + text);
  errno = 0;
  char *end = nullptr;
  const auto value = std::strtoull(stripped.c_str(), &end, 10);
  if (errno == ERANGE || end == stripped.c_str() || !trim(end).empty() || value > std::numeric_limits<size_t>::max())
    throw std::runtime_error("Invalid unsigned integer for " + context + ": " + text);
  return static_cast<size_t>(value);
}

double parse_double_value(const std::string &text, const std::string &context) {
  char *end = nullptr;
  const auto value = std::strtod(text.c_str(), &end);
  if (end == text.c_str() || !trim(end).empty() || !std::isfinite(value))
    throw std::runtime_error("Invalid number for " + context + ": " + text);
  return value;
}

std::string param_value(const ParamSections &sections, const std::string &section, const std::string &key) {
  const auto section_it = sections.values.find(section);
  if (section_it == sections.values.end()) throw std::runtime_error("Missing [" + section + "] section in param.");
  const auto value_it = section_it->second.find(key);
  if (value_it == section_it->second.end()) throw std::runtime_error("Missing parameter " + key + " in [" + section + "].");
  return value_it->second;
}

bool optional_bool_param(const ParamSections &sections, const std::string &section, const std::string &key, const bool fallback) {
  const auto section_it = sections.values.find(section);
  if (section_it == sections.values.end()) return fallback;
  const auto value_it = section_it->second.find(key);
  if (value_it == section_it->second.end()) return fallback;
  const auto value = value_it->second;
  if (value == "true") return true;
  if (value == "false") return false;
  throw std::runtime_error("Parameter " + key + " must be true or false.");
}

TemplateHeader read_template_header(DataTemplateReader &data_in) {
  const auto header_line = data_in.next_data_line();
  if (!header_line) throw std::runtime_error("Missing data.in header.");
  std::istringstream header_stream(*header_line);
  TemplateHeader header;
  header_stream >> header.channels >> header.chain_sites >> header.subspaces;
  if (!header_stream) throw std::runtime_error("Invalid data.in header.");
  return header;
}

NRG::Matrix_traits<double> parse_inline_matrix(const std::vector<std::string> &rows, const size_t expected_cols,
                                               const std::string &context) {
  NRG::Matrix_traits<double> matrix(static_cast<Eigen::Index>(rows.size()), static_cast<Eigen::Index>(expected_cols));
  for (size_t row = 0; row < rows.size(); ++row) {
    const auto tokens = fields(rows[row]);
    if (tokens.size() != expected_cols)
      throw std::runtime_error("Inline matrix dimension mismatch for " + context + ".");
    for (size_t col = 0; col < expected_cols; ++col)
      matrix(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col)) = parse_double_value(tokens[col], context);
  }
  return matrix;
}

void write_matrix(std::ostream &out, const NRG::Matrix_traits<double> &matrix) {
  out << std::setprecision(18);
  for (Eigen::Index row = 0; row < matrix.rows(); ++row) {
    for (Eigen::Index col = 0; col < matrix.cols(); ++col) {
      if (col != 0) out << ' ';
      const auto value = matrix(row, col);
      out << (std::abs(value) > 1e-14 ? value : 0.0);
    }
    out << '\n';
  }
}

SeedData read_seed_data(DataTemplateReader &data_in, const size_t subspaces, NRG::Spawn::MatrixEvaluator &evaluator,
                        const std::filesystem::path &template_dir, const bool write_legacy_outputs) {
  SeedData seed;

  for (size_t subspace = 0; subspace < subspaces; ++subspace) {
    const auto qn_line = data_in.next_data_line();
    const auto size_line = data_in.next_data_line();
    const auto diag_line = data_in.next_data_line();
    if (!qn_line || !size_line || !diag_line) throw std::runtime_error("Unexpected end of data.in while reading subspaces.");
    if (subspace == 0) seed.factor = data_in.factor();

    const auto qn = parse_invar_line(*qn_line, "seed subspace");
    const auto qn_name = legacy_qn_name(qn);
    const auto expected_size = parse_size_t_value(*size_line, "subspace size for " + qn_name);
    if (expected_size == 0) throw std::runtime_error("Invalid subspace size for " + qn_name + ".");

    std::istringstream diag_stream(*diag_line);
    std::string keyword, matrix_file;
    diag_stream >> keyword >> matrix_file;
    if (keyword != "DIAG" || matrix_file.empty()) throw std::runtime_error("Expected DIAG line for " + qn_name + ".");

    auto matrix = evaluator.evaluate_matrix(read_text_file(resolve_template_file(template_dir, matrix_file)), matrix_file);
    if (matrix.rows() != static_cast<Eigen::Index>(expected_size) || matrix.cols() != static_cast<Eigen::Index>(expected_size))
      throw std::runtime_error("Matrix dimension mismatch for " + qn_name + ".");

    if (write_legacy_outputs) {
      save_text_matrix("ham." + qn_name, matrix);
      save_text_matrix("ham", matrix);
    }

    auto raw = NRG::diagonalise_dsyev(matrix);
    for (auto &value : raw.val) value *= seed.factor;
    if (!raw.val.empty()) seed.smallest = std::min(seed.smallest, raw.val.front());

    if (write_legacy_outputs) {
      save_values("val", raw.val);
      save_binary_matrix("vec", raw.vec);
      save_binary_matrix("vec." + qn_name, raw.vec);
    }

    if (seed.dimensions.contains(qn)) throw std::runtime_error("Duplicate seed subspace " + qn_name + ".");
    seed.dimensions[qn] = expected_size;
    seed.eigenvectors[qn] = raw.vec;
    seed.subspaces.push_back(SeedSubspace{*qn_line, qn, expected_size, std::move(raw.val), std::move(raw.vec)});
  }

  if (seed.subspaces.empty() || seed.smallest == std::numeric_limits<double>::max())
    throw std::runtime_error("No seed subspaces found in data.in.");
  seed.ground_energy = seed.smallest / seed.factor;
  return seed;
}

void write_seed_energy_block(std::ostream &out, const SeedData &seed) {
  out << std::setprecision(18);
  for (const auto &subspace : seed.subspaces) {
    out << subspace.qn_line << '\n';
    out << subspace.dimension << '\n';
    for (size_t i = 0; i < subspace.eigenvalues.size(); ++i) {
      if (i != 0) out << ' ';
      out << subspace.eigenvalues[i] - seed.smallest;
    }
    out << '\n';
  }
}

NRG::Matrix_traits<double> read_operator_matrix(DataTemplateReader &data_in, const std::string &first_line,
                                                const size_t rows, const size_t cols,
                                                NRG::Spawn::MatrixEvaluator &evaluator,
                                                const std::filesystem::path &template_dir,
                                                const std::string &context) {
  if (!first_line.empty() && first_line.front() == 'o') {
    auto matrix = evaluator.evaluate_matrix(read_text_file(resolve_template_file(template_dir, first_line)), first_line);
    if (matrix.rows() != static_cast<Eigen::Index>(rows) || matrix.cols() != static_cast<Eigen::Index>(cols))
      throw std::runtime_error("Matrix dimension mismatch for " + context + ".");
    return matrix;
  }

  std::vector<std::string> matrix_rows;
  matrix_rows.push_back(first_line);
  for (size_t row = 1; row < rows; ++row) {
    const auto line = data_in.next_data_line();
    if (!line) throw std::runtime_error("Unexpected end of data.in while reading matrix for " + context + ".");
    matrix_rows.push_back(*line);
  }
  return parse_inline_matrix(matrix_rows, cols, context);
}

MatrixBlockData read_matrix_block(const std::string &header_line, DataTemplateReader &data_in, const SeedData &seed,
                                  NRG::Spawn::MatrixEvaluator &evaluator, const std::filesystem::path &template_dir) {
  const auto count_line = data_in.next_data_line();
  if (!count_line) throw std::runtime_error("Unexpected end of data.in while reading matrix-element count.");
  const auto count = parse_size_t_value(*count_line, "matrix-element count");

  MatrixBlockData block;
  block.header_line = header_line;
  block.block_type = header_line.front();
  block.elements.reserve(count);

  for (size_t i = 0; i < count; ++i) {
    const auto qn_line = data_in.next_data_line();
    if (!qn_line) throw std::runtime_error("Unexpected end of data.in while reading matrix-element quantum numbers.");
    const auto [qn1, qn2] = parse_invar_pair_line(*qn_line, "matrix-element block");
    const auto qn1_name = legacy_qn_name(qn1);
    const auto qn2_name = legacy_qn_name(qn2);
    const auto dim1_it = seed.dimensions.find(qn1);
    const auto dim2_it = seed.dimensions.find(qn2);
    if (dim1_it == seed.dimensions.end() || dim2_it == seed.dimensions.end())
      throw std::runtime_error("Matrix-element block references unknown subspace " + qn1_name + " -> " + qn2_name + ".");
    const auto vec1_it = seed.eigenvectors.find(qn1);
    const auto vec2_it = seed.eigenvectors.find(qn2);
    if (vec1_it == seed.eigenvectors.end() || vec2_it == seed.eigenvectors.end())
      throw std::runtime_error("Missing eigenvectors for subspace " + qn1_name + " -> " + qn2_name + ".");
    const auto dim1 = dim1_it->second;
    const auto dim2 = dim2_it->second;

    const auto first_matrix_line = data_in.next_data_line();
    if (!first_matrix_line) throw std::runtime_error("Unexpected end of data.in while reading matrix for " + qn1_name + " -> " + qn2_name + ".");
    const auto matrix = read_operator_matrix(data_in, *first_matrix_line, dim1, dim2, evaluator, template_dir, qn1_name + " -> " + qn2_name);
    NRG::Matrix_traits<double> transformed = vec1_it->second * matrix * vec2_it->second.transpose();
    block.elements.push_back(MatrixElementData{*qn_line, qn1, qn2, std::move(transformed)});
  }

  return block;
}

void write_matrix_block(std::ostream &out, const MatrixBlockData &block) {
  out << block.header_line << '\n';
  out << block.elements.size() << '\n';
  for (const auto &element : block.elements) {
    out << element.qn_line << '\n';
    write_matrix(out, element.matrix);
  }
}

bool is_operator_block(const char block_type) {
  return block_type == 's' || block_type == 'p' || block_type == 'd' || block_type == 'v' || block_type == 't' ||
         block_type == 'o' || block_type == 'g' || block_type == 'q';
}

TemplateTailData process_template_tail(DataTemplateReader &data_in, std::ostream &out, const SeedData &seed,
                                       NRG::Spawn::MatrixEvaluator &evaluator, const std::filesystem::path &template_dir,
                                       const size_t expected_wilson_blocks) {
  TemplateTailData tail;
  size_t wilson_blocks = 0;
  while (const auto line = data_in.next_data_line()) {
    const auto block_type = line->front();
    if (block_type == 'f') {
      if (++wilson_blocks > expected_wilson_blocks)
        throw std::runtime_error("More Wilson-chain operator blocks than expected for the selected symmetry in data.in.");
      auto block = read_matrix_block(*line, data_in, seed, evaluator, template_dir);
      write_matrix_block(out, block);
      tail.matrix_blocks.push_back(std::move(block));
      continue;
    }

    if (*line == "e") {
      const auto discarded = data_in.next_data_line();
      if (!discarded) throw std::runtime_error("Missing placeholder GS energy after e block in data.in.");
      out << "e\n" << std::setprecision(18) << seed.ground_energy << '\n';
      tail.has_ground_energy = true;
      tail.ground_energy = seed.ground_energy;
      continue;
    }

    if (block_type == 'z' || block_type == 'Z' || block_type == 'T') break;

    if (is_operator_block(block_type)) {
      auto block = read_matrix_block(*line, data_in, seed, evaluator, template_dir);
      write_matrix_block(out, block);
      tail.matrix_blocks.push_back(std::move(block));
      continue;
    }

    throw std::runtime_error("Unsupported data.in block: " + *line);
  }

  if (wilson_blocks != expected_wilson_blocks)
    throw std::runtime_error("Wilson-chain operator block count does not match the selected symmetry.");
  if (!tail.has_ground_energy) throw std::runtime_error("Missing e block in template/data.in.");
  return tail;
}

void write_coefficient_table(std::ostream &out, const std::vector<double> &values, const size_t max_index) {
  if (values.size() <= max_index) throw std::runtime_error("Wilson coefficient table is shorter than Nmax.");
  out << max_index << '\n';
  out << std::setprecision(16);
  for (size_t i = 0; i <= max_index; ++i) out << values[i] << '\n';
}

void write_z_coefficients(std::ostream &out, const NRG::Tools::NrgChain::WilsonData &wilson, const size_t coefchannels,
                          const size_t max_index) {
  if (wilson.channels.size() != coefchannels)
    throw std::runtime_error("Only one coefficient table per physical channel is supported in this instantiate slice.");
  out << "z\n";
  for (const auto &channel : wilson.channels) write_coefficient_table(out, channel.xi, max_index);
  for (const auto &channel : wilson.channels) write_coefficient_table(out, channel.zeta, max_index);
}

NRG::Spawn::MatrixEvaluator make_matrix_evaluator(const ParamSections &sections, const NRG::Tools::NrgChain::WilsonData &wilson) {
  NRG::Spawn::MatrixEvalContext context;
  if (const auto it = sections.values.find("extra"); it != sections.values.end()) context.variables = numeric_variables(it->second);
  context.wilson = wilson;
  return NRG::Spawn::MatrixEvaluator(std::move(context));
}

class ScopedCoutRedirect {
 private:
  std::ostringstream sink;
  std::streambuf *old_buffer = nullptr;

 public:
  ScopedCoutRedirect() : old_buffer(std::cout.rdbuf(sink.rdbuf())) {}
  ScopedCoutRedirect(const ScopedCoutRedirect &) = delete;
  ScopedCoutRedirect &operator=(const ScopedCoutRedirect &) = delete;
  ~ScopedCoutRedirect() { std::cout.rdbuf(old_buffer); }
};

std::unique_ptr<NRG::Params> make_instantiate_params(const std::string &param_filename) {
  ScopedCoutRedirect quiet;
  auto params = std::make_unique<NRG::Params>(param_filename, "param", std::make_unique<NRG::Workdir>(), true, true);
  params->silent = true;
  return params;
}

void initialize_template_symmetry(NRG::Params &params, const TemplateHeader &header) {
  if (params.symtype.value().empty()) throw std::runtime_error("Parameter symtype must be set before parsing template/data.in.");
  ScopedCoutRedirect quiet;
  const auto symmetry = NRG::set_symmetry<double>(params, params.symtype.value(), header.channels);
  (void)symmetry;
}

void validate_generated_data(const std::string &param_filename, const std::string &data_text) {
  std::istringstream data_in(data_text);
  ScopedCoutRedirect quiet;
  auto params = make_instantiate_params(param_filename);
  const NRG::InputData<double> input(*params, data_in);
}

void run_diag_seed_only(const Options &options) {
  auto sections = parse_param_sections(options.param_filename);
  write_param_section_files(sections, options.param_filename);
  auto params = make_instantiate_params(options.param_filename);

  const auto wilson = NRG::Tools::NrgChain::calculate_from_file(options.param_filename);
  if (wilson.channels.size() != 1)
    throw std::runtime_error("Only single-channel Wilson generation is supported in this instantiate slice.");
  write_wilson_channel(wilson.channels.front(), 1);

  auto evaluator = make_matrix_evaluator(sections, wilson);

  const auto template_dir = std::filesystem::path(options.template_dir);
  DataTemplateReader data_in(resolve_template_file(template_dir, "data.in"));
  const auto header = read_template_header(data_in);
  if (header.channels != 1) throw std::runtime_error("Only single-channel data.in templates are supported in this slice.");
  initialize_template_symmetry(*params, header);
  const auto seed = read_seed_data(data_in, header.subspaces, evaluator, template_dir, true);
  std::cout << "E_gs=" << std::setprecision(18) << seed.ground_energy << '\n';
}

void run_full_instantiation(const Options &options) {
  auto sections = parse_param_sections(options.param_filename);
  TemporaryCleanupGuard cleanup(!options.generate_temporaries, options.param_filename, sections);
  if (options.generate_temporaries) write_param_section_files(sections, options.param_filename);
  auto params = make_instantiate_params(options.param_filename);

  const auto wilson = NRG::Tools::NrgChain::calculate_from_file(options.param_filename);
  if (wilson.channels.size() != 1)
    throw std::runtime_error("Only single-channel Wilson generation is supported in this instantiate slice.");
  if (options.generate_temporaries) write_wilson_channel(wilson.channels.front(), 1);

  const auto nmax = parse_size_t_value(param_value(sections, "param", "Nmax"), "Nmax");
  const auto polarized = optional_bool_param(sections, "param", "polarized", false);
  if (polarized) throw std::runtime_error("Spin-polarized coefficient tables are not supported in this instantiate slice.");

  auto evaluator = make_matrix_evaluator(sections, wilson);

  std::ostringstream data_buffer;

  const auto template_dir = std::filesystem::path(options.template_dir);
  DataTemplateReader data_in(resolve_template_file(template_dir, "data.in"), &data_buffer);
  const auto header = read_template_header(data_in);
  if (header.channels != 1) throw std::runtime_error("Only single-channel data.in templates are supported in this slice.");
  initialize_template_symmetry(*params, header);
  data_buffer << header.channels << ' ' << nmax << ' ' << header.subspaces << '\n';

  const auto seed = read_seed_data(data_in, header.subspaces, evaluator, template_dir, options.generate_temporaries);
  write_seed_energy_block(data_buffer, seed);
  const auto tail = process_template_tail(data_in, data_buffer, seed, evaluator, template_dir, params->channels * params->perchannel);
  if (tail.ground_energy != seed.ground_energy) throw std::runtime_error("Internal ground-energy mismatch while instantiating data.");
  write_z_coefficients(data_buffer, wilson, params->coefchannels, nmax);

  const auto data_text = data_buffer.str();
  validate_generated_data(options.param_filename, data_text);

  std::ofstream data_out("data");
  if (!data_out) throw std::runtime_error("Can't open data for writing.");
  data_out << data_text;

  std::cout << "E_gs=" << std::setprecision(18) << seed.ground_energy << '\n';
}

void run_wilson_only(const Options &options) {
  const auto wilson = NRG::Tools::NrgChain::calculate_from_file(options.param_filename);
  if (wilson.channels.size() != 1)
    throw std::runtime_error("Only single-channel Wilson generation is supported in this instantiate slice.");
  write_wilson_channel(wilson.channels.front(), 1);
}

} // namespace

int main(int argc, char *argv[]) {
  try {
    const auto options = parse_options(argc, argv);
    if (options.wilson_only && options.diag_seed_only)
      throw std::runtime_error("Use only one of --wilson-only and --diag-seed-only.");
    if (options.wilson_only) {
      run_wilson_only(options);
    } else if (options.diag_seed_only) {
      run_diag_seed_only(options);
    } else {
      run_full_instantiation(options);
    }
    return EXIT_SUCCESS;
  } catch (const std::exception &e) {
    std::cerr << "instantiate: error: " << e.what() << '\n';
    return EXIT_FAILURE;
  }
}
