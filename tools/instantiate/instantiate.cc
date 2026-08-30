// Initial data instantiation tool
// First production slice: in-process Wilson-chain generation.

#include <cstdlib>
#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cstddef>
#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ios>
#include <istream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <sstream>
#include <system_error>
#include <streambuf>
#include <utility>
#include <vector>

#include <diag.hpp>
#include <tridiag.hpp>
#include <read-input.hpp>
#include <traits.hpp>
#include <workdir.hpp>
#include <parse_bool.hpp>

#include <common/version.hpp>

#include "../common/diagnostics.hpp"
#include "../common/parser.hpp"
#include "../common/tabulated_density.hpp"
#include "matrix_evaluator.hpp"
#include "nrgchain.hpp"

namespace {

using std::size_t;

struct Options {
  std::string param_filename = "param";
  std::string template_dir = "template";
  bool wilson_only = false;
  bool diag_seed_only = false;
  bool generate_temporaries = false;
  int verbosity = 0;
};

void usage(std::ostream &out = std::cout) {
  out << "Usage: instantiate [options]\n"
      << "  --wilson-only          generate Wilson-chain files only\n"
      << "  --diag-seed-only       generate and diagonalize the seed only\n"
      << "  --generate-temporaries publish legacy temporary files\n"
      << "  --param FILE           read parameters from FILE (default: param)\n"
      << "  --template-dir DIR     find template inputs under DIR (default: template)\n"
      << "  -v                     show resolved configuration on standard error\n"
      << "  -vv                    increase verbosity further\n"
      << "  -V, --version          show project version\n"
      << "  -h, --help             show this help\n";
}

Options parse_options(const int argc, char *argv[]) {
  Options options;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "-h" || arg == "--help") {
      usage();
      std::exit(EXIT_SUCCESS);
    }
    if (arg.size() >= 2 && arg.front() == '-' &&
        std::all_of(arg.begin() + 1, arg.end(), [](const char ch) { return ch == 'v'; })) {
      options.verbosity += static_cast<int>(arg.size() - 1);
      continue;
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

using ArtifactManifest = std::map<std::filesystem::path, std::filesystem::path>;

void close_output_checked(std::ofstream &out, const std::filesystem::path &filename) {
  out.close();
  if (!out) throw std::runtime_error("Error writing " + filename.string() + ".");
}

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
    close_output_checked(out, param_filename + "." + section);
  }
}

void register_artifact(ArtifactManifest &artifacts, const std::filesystem::path &staged,
                       const std::filesystem::path &destination) {
  artifacts[std::filesystem::absolute(destination)] = staged;
}

void write_staged_param_section_files(const ParamSections &sections, const std::string &param_filename,
                                      const std::filesystem::path &staging_dir, ArtifactManifest &artifacts) {
  size_t index = 0;
  for (const auto &[section, lines] : sections.lines) {
    const auto staged = staging_dir / ("param-section-" + std::to_string(index++));
    std::ofstream out(staged);
    if (!out) throw std::runtime_error("Can't open " + staged.string() + " for writing.");
    for (const auto &line : lines) out << line << '\n';
    close_output_checked(out, staged);
    register_artifact(artifacts, staged, param_filename + "." + section);
  }
}

void register_existing_artifact(ArtifactManifest &artifacts, const std::filesystem::path &staging_dir,
                                const std::filesystem::path &filename) {
  const auto staged = staging_dir / filename;
  std::error_code ec;
  const auto exists = std::filesystem::exists(staged, ec);
  if (ec) throw std::runtime_error("Can't inspect " + staged.string() + ": " + ec.message());
  if (exists) register_artifact(artifacts, staged, filename);
}

void publish_file(const std::filesystem::path &staged, const std::filesystem::path &destination) {
  std::error_code ec;
  std::filesystem::rename(staged, destination, ec);
  if (!ec) return;
  if (ec != std::errc::cross_device_link)
    throw std::runtime_error("Can't publish " + destination.string() + " from " + staged.string() + ": " + ec.message());

  NRG::Workdir local_staging(destination.parent_path().string(), true);
  const auto local_file = std::filesystem::path(local_staging.get()) / "artifact";
  ec.clear();
  std::filesystem::copy_file(staged, local_file, ec);
  if (ec)
    throw std::runtime_error("Can't stage " + destination.string() + " from " + staged.string() + ": " + ec.message());
  std::filesystem::rename(local_file, destination, ec);
  if (ec)
    throw std::runtime_error("Can't publish " + destination.string() + " from " + local_file.string() + ": " + ec.message());
}

void write_text_file_checked(const std::filesystem::path &filename, const std::string &text) {
  std::ofstream out(filename);
  if (!out) throw std::runtime_error("Can't open " + filename.string() + " for writing.");
  out << text;
  close_output_checked(out, filename);
}

std::map<std::string, double> numeric_variables(const std::map<std::string, std::string> &values) {
  std::map<std::string, double> variables;
  for (const auto &[key, value] : values) {
    char *end = nullptr;
    const auto parsed = std::strtod(value.c_str(), &end);
    if (end != value.c_str() && trim(end).empty()) variables[key] = parsed;
  }
  return variables;
}

void write_values(const std::filesystem::path &filename, const std::vector<double> &values) {
  std::ofstream out(filename);
  if (!out) throw std::runtime_error("Can't open " + filename.string() + " for writing.");
  out << std::setprecision(16);
  for (const auto value : values) out << value << '\n';
  close_output_checked(out, filename);
}

void write_scalar(const std::filesystem::path &filename, const double value) {
  std::ofstream out(filename);
  if (!out) throw std::runtime_error("Can't open " + filename.string() + " for writing.");
  out << std::setprecision(16) << value << '\n';
  close_output_checked(out, filename);
}

void write_wilson_channel(const NRG::Tools::NrgChain::WilsonChannel &channel, const size_t index,
                          const std::filesystem::path &output_dir = ".") {
  const auto suffix = std::to_string(index) + ".dat";
  write_values(output_dir / ("xi" + suffix), channel.xi);
  write_values(output_dir / ("zeta" + suffix), channel.zeta);
  write_scalar(output_dir / ("theta" + suffix), channel.theta);
}

void save_text_matrix(const std::filesystem::path &filename, const NRG::Matrix_traits<double> &matrix) {
  std::ofstream out(filename);
  if (!out) throw std::runtime_error("Can't open " + filename.string() + " for writing.");
  out << std::setprecision(18);
  for (Eigen::Index row = 0; row < matrix.rows(); ++row) {
    for (Eigen::Index col = 0; col < matrix.cols(); ++col) {
      if (col != 0) out << ' ';
      out << matrix(row, col);
    }
    out << '\n';
  }
  close_output_checked(out, filename);
}

void save_binary_matrix(const std::filesystem::path &filename, const NRG::Matrix_traits<double> &matrix) {
  if (matrix.rows() > static_cast<Eigen::Index>(std::numeric_limits<unsigned int>::max()) ||
      matrix.cols() > static_cast<Eigen::Index>(std::numeric_limits<unsigned int>::max()))
    throw std::runtime_error("Matrix too large for legacy binary format.");
  std::ofstream out(filename, std::ios::binary);
  if (!out) throw std::runtime_error("Can't open " + filename.string() + " for writing.");
  const auto rows = static_cast<unsigned int>(matrix.rows());
  const auto cols = static_cast<unsigned int>(matrix.cols());
  out.write(reinterpret_cast<const char *>(&rows), sizeof(rows));
  out.write(reinterpret_cast<const char *>(&cols), sizeof(cols));
  for (Eigen::Index row = 0; row < matrix.rows(); ++row)
    for (Eigen::Index col = 0; col < matrix.cols(); ++col) {
      const auto value = matrix(row, col);
      out.write(reinterpret_cast<const char *>(&value), sizeof(value));
    }
  close_output_checked(out, filename);
}

void save_values(const std::filesystem::path &filename, const std::vector<double> &values) {
  std::ofstream out(filename);
  if (!out) throw std::runtime_error("Can't open " + filename.string() + " for writing.");
  out << std::setprecision(18);
  for (const auto value : values) out << value << '\n';
  close_output_checked(out, filename);
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

TemplateHeader read_template_header(DataTemplateReader &data_in) {
  const auto header_line = data_in.next_data_line();
  if (!header_line) throw std::runtime_error("Missing data.in header.");
  std::istringstream header_stream(*header_line);
  TemplateHeader header;
  header_stream >> header.channels >> header.chain_sites >> header.subspaces;
  if (!header_stream) throw std::runtime_error("Invalid data.in header.");
  return header;
}

using NrgChainTableMode = NRG::Tools::NrgChain::TableMode;

struct ResolvedNrgChainConfiguration {
  NrgChainTableMode mode;
  double Lambda;
  double z;
  double bandrescale;
  double xmax;
  unsigned int Nmax;
  unsigned int mMAX;
  unsigned int preccpp;
  bool adapt;
  bool rescalexi;
  std::string band;
  std::string dos;
  NRG::Tools::InterpolationMethod density_interpolation;
  bool requested_tables_save;
  bool requested_tables_load;
  bool requested_tridiagonalize;
  bool tables_save;
  bool tables_load;
  bool tridiagonalize;
};

struct NrgChainInvocation {
  std::map<std::string, std::string> parameters;
  ResolvedNrgChainConfiguration configuration;
};

std::map<std::string, std::string> read_nrgchain_parameters(const std::string &filename) {
  std::ifstream input(filename);
  if (!input) throw std::runtime_error("Can't open " + filename + " for reading.");
  std::map<std::string, std::string> parameters;
  if (NRG::Tools::find_block(input, "param")) NRG::Tools::parse_key_value_block(input, parameters);
  return parameters;
}

ResolvedNrgChainConfiguration resolve_nrgchain_configuration(
  const std::map<std::string, std::string> &parameters, const NrgChainTableMode mode) {
  auto text = [&parameters](const std::string_view key) -> const std::string * {
    const auto it = parameters.find(std::string(key));
    return it == parameters.end() ? nullptr : &it->second;
  };
  auto number = [&text](const std::string_view key, const double default_value) {
    const auto *value = text(key);
    return value ? NRG::Tools::parse_parameter_double(*value, key) : default_value;
  };
  auto integer = [&text](const std::string_view key, const unsigned int default_value) {
    const auto *value = text(key);
    if (!value) return default_value;
    const auto parsed = NRG::Tools::parse_parameter_int(*value, key);
    if (parsed < 0) throw std::invalid_argument(std::string(key) + " must be nonnegative.");
    return static_cast<unsigned int>(parsed);
  };
  auto boolean = [&text](const std::string_view key, const bool default_value) {
    const auto *value = text(key);
    return value ? NRG::parse_bool(*value) : default_value;
  };
  auto string = [&text](const std::string_view key, std::string default_value) {
    const auto *value = text(key);
    return value ? *value : std::move(default_value);
  };

  ResolvedNrgChainConfiguration config{
    .mode = mode,
    .Lambda = number("Lambda", 2.0),
    .z = number("z", 1.0),
    .bandrescale = number("bandrescale", 1.0),
    .xmax = number("xmax", 30.0),
    .Nmax = integer("Nmax", 0),
    .mMAX = 0,
    .preccpp = integer("preccpp", 2000),
    .adapt = boolean("adapt", false),
    .rescalexi = boolean("rescalexi", false),
    .band = string("band", "adapt"),
    .dos = string("dos", "Delta.dat"),
    .density_interpolation = NRG::Tools::parse_density_interpolation_method(
      string("density_interpolation", "linear")),
    .requested_tables_save = boolean("nrgchain_tables_save", false),
    .requested_tables_load = boolean("nrgchain_tables_load", false),
    .requested_tridiagonalize = parameters.contains("nrgchain_tridiag")
                                      ? boolean("nrgchain_tridiag", true)
                                      : boolean("nrgchains_tridiag", true),
    .tables_save = false,
    .tables_load = false,
    .tridiagonalize = false,
  };
  if (!(config.Lambda > 1.0)) throw std::invalid_argument("Lambda must be greater than 1.");
  if (!(config.z > 0.0 && config.z <= 1.0)) throw std::invalid_argument("z must satisfy 0 < z <= 1.");
  if (!(config.bandrescale > 0.0)) throw std::invalid_argument("bandrescale must be positive.");
  if (!(config.xmax >= 1.0)) throw std::invalid_argument("xmax must be greater than or equal to 1.");
  if (config.Nmax > static_cast<unsigned int>(std::numeric_limits<int>::max() / 2))
    throw std::invalid_argument("Nmax is too large to derive mMAX.");
  config.mMAX = integer("mMAX", 2U * config.Nmax);
  if (config.mMAX == 0) throw std::invalid_argument("mMAX must be greater than 0.");
  if (config.preccpp <= 10) throw std::invalid_argument("preccpp must be greater than 10.");
  config.tables_save = config.requested_tables_save;
  config.tables_load = config.requested_tables_load;
  config.tridiagonalize = config.requested_tridiagonalize;
  if (mode == NrgChainTableMode::SaveOnly) {
    config.tables_load = false;
    config.tables_save = true;
    config.tridiagonalize = false;
  } else if (mode == NrgChainTableMode::LoadAndTridiagonalize) {
    config.tables_load = true;
    config.tables_save = false;
    config.tridiagonalize = true;
  }
  if (config.tables_load && config.tables_save)
    throw std::invalid_argument("nrgchain_tables_load and nrgchain_tables_save cannot both be true.");
  return config;
}

NrgChainInvocation make_nrgchain_invocation(const std::string &filename, const NrgChainTableMode mode) {
  // Keep reporting and core execution on the same parser output and mode.
  auto parameters = read_nrgchain_parameters(filename);
  auto configuration = resolve_nrgchain_configuration(parameters, mode);
  if (!configuration.tridiagonalize)
    throw std::invalid_argument("instantiate requires nrgchain_tridiag=true.");
  return {std::move(parameters), std::move(configuration)};
}

std::string_view nrgchain_mode_name(const NrgChainTableMode mode) {
  switch (mode) {
    case NrgChainTableMode::Calculate: return "parameter-file";
    case NrgChainTableMode::SaveOnly: return "save-only";
    case NrgChainTableMode::LoadAndTridiagonalize: return "load-and-tridiagonalize";
  }
  throw std::logic_error("Unknown nrgchain table mode.");
}

void report_configuration(const Options &options, const NrgChainInvocation &invocation) {
  if (options.verbosity == 0) return;

  const auto &values = invocation.parameters;
  const auto &config = invocation.configuration;
  NRG::Tools::ConfigurationReport report("instantiate");
  report.value("verbosity", options.verbosity);
  if (options.wilson_only)
    report.value("mode", "wilson-only");
  else if (options.diag_seed_only)
    report.value("mode", "diag-seed-only");
  else
    report.resolved("mode", "full", "no restricted-mode option");
  report.value("generate_temporaries", options.generate_temporaries);
  report.value("param.file", options.param_filename);
  report.value("nrgchain.mode", nrgchain_mode_name(config.mode));

  if (!options.wilson_only) {
    report.value("template.directory", options.template_dir);
    const auto data_filename = resolve_template_file(options.template_dir, "data.in");
    report.resolved("template.data", data_filename.string(), "template search");
    DataTemplateReader data_in(data_filename);
    const auto header = read_template_header(data_in);
    report.value("template.channels", header.channels);
    report.value("template.chain_sites", header.chain_sites);
    report.value("template.subspaces", header.subspaces);
    for (int seed_header_line = 0; seed_header_line < 3; ++seed_header_line)
      if (!data_in.next_data_line()) throw std::runtime_error("Template contains incomplete seed data after its header.");
    report.resolved("template.scale_factor", data_in.factor(), "template # SCALE value");
    report.value("matrix.output_precision", 18);
    report.value("matrix.chop_tolerance", 1e-14);
  }

  auto parameter = [&values, &report](const std::string_view name, const std::string_view key, const auto &value,
                                     const std::string_view reason = "nrgchain default") {
    if (values.contains(std::string(key)))
      report.value(name, value);
    else
      report.resolved(name, value, reason);
  };
  auto requested_tridiagonalization = [&](const std::string_view name) {
    if (values.contains("nrgchain_tridiag")) {
      report.value(name, config.requested_tridiagonalize);
      report.value("nrgchain.tridiagonalize.parameter", "nrgchain_tridiag");
      if (values.contains("nrgchains_tridiag"))
        report.value("nrgchain.nrgchains_tridiag.ignored", values.at("nrgchains_tridiag"));
    } else if (values.contains("nrgchains_tridiag")) {
      report.resolved(name, config.requested_tridiagonalize, "legacy nrgchains_tridiag parameter");
      report.value("nrgchain.tridiagonalize.parameter", "nrgchains_tridiag (legacy)");
    } else {
      report.resolved(name, config.requested_tridiagonalize, "nrgchain default");
      report.resolved("nrgchain.tridiagonalize.parameter", "<default>", "no tridiagonalization parameter");
    }
  };

  parameter("nrgchain.Lambda", "Lambda", config.Lambda);
  parameter("nrgchain.z", "z", config.z);
  parameter("nrgchain.bandrescale", "bandrescale", config.bandrescale);
  parameter("nrgchain.xmax", "xmax", config.xmax);
  parameter("nrgchain.Nmax", "Nmax", config.Nmax);
  parameter("nrgchain.mMAX", "mMAX", config.mMAX, "2 * Nmax");
  parameter("nrgchain.preccpp", "preccpp", config.preccpp);
  report.value("nrgchain.output_precision", 16);
  parameter("nrgchain.adapt", "adapt", config.adapt);
  parameter("nrgchain.rescalexi", "rescalexi", config.rescalexi);
  parameter("nrgchain.band", "band", config.band);
  parameter("nrgchain.dos", "dos", config.dos);
  const auto density_method = NRG::Tools::interpolation_method_name(config.density_interpolation);
  if (config.tables_load) {
    report.resolved("nrgchain.density_interpolation", "inactive", "coefficient tables are loaded");
    report.resolved("nrgchain.density_integration", "inactive", "coefficient tables are loaded");
  } else if (config.band == "flat") {
    report.resolved("nrgchain.density_interpolation", density_method, "flat density is interpolation-independent");
    report.value("nrgchain.density_integration", "exact interpolant primitive");
  } else {
    parameter("nrgchain.density_interpolation", "density_interpolation", density_method);
    report.value("nrgchain.density_integration", "exact interpolant primitive");
  }

  if (config.mode == NrgChainTableMode::Calculate) {
    parameter("nrgchain.tables_save", "nrgchain_tables_save", config.tables_save);
    parameter("nrgchain.tables_load", "nrgchain_tables_load", config.tables_load);
    requested_tridiagonalization("nrgchain.tridiagonalize");
  } else {
    parameter("nrgchain.tables_save.requested", "nrgchain_tables_save", config.requested_tables_save);
    parameter("nrgchain.tables_load.requested", "nrgchain_tables_load", config.requested_tables_load);
    requested_tridiagonalization("nrgchain.tridiagonalize.requested");
    const auto reason = config.mode == NrgChainTableMode::SaveOnly ? "save-only mode override"
                                                                   : "load-and-tridiagonalize mode override";
    report.resolved("nrgchain.tables_save", config.tables_save, reason);
    report.resolved("nrgchain.tables_load", config.tables_load, reason);
    report.resolved("nrgchain.tridiagonalize", config.tridiagonalize, reason);
  }

  if (config.tables_load) {
    report.resolved("nrgchain.density_source", "coefficient-tables", "tables_load=true");
    report.resolved("nrgchain.dos.active", false, "coefficient tables are loaded");
    report.value("nrgchain.theta_input", "theta.dat");
  } else if (config.band == "flat") {
    report.value("nrgchain.density_source", "flat");
    report.resolved("nrgchain.dos.active", false, "flat band");
  } else {
    report.value("nrgchain.density_source", "dos-file");
    report.value("nrgchain.dos.active", true);
    report.value("nrgchain.f_positive", "FSOL.dat");
    report.value("nrgchain.f_negative", "FSOLNEG.dat");
    if (config.adapt) {
      report.value("nrgchain.g_positive", "GSOL.dat");
      report.value("nrgchain.g_negative", "GSOLNEG.dat");
    }
  }

  if (!options.wilson_only) {
    report.value("output.stdout", "E_gs");
  } else {
    report.resolved("output.stdout", "none", "wilson-only mode");
  }
  if (!options.wilson_only && !options.diag_seed_only) {
    report.value("output.data", "data");
  } else {
    report.resolved("output.data", "inactive", "restricted mode");
  }
  if (options.wilson_only || options.diag_seed_only || options.generate_temporaries) {
    report.value("output.wilson_theta", "theta1.dat");
    report.value("output.wilson_xi", "xi1.dat");
    report.value("output.wilson_zeta", "zeta1.dat");
  } else {
    report.resolved("output.wilson_files", "temporary", "full mode without --generate-temporaries");
  }
  if (options.diag_seed_only || options.generate_temporaries) {
    const auto sections = parse_param_sections(options.param_filename);
    std::string section_outputs;
    for (const auto &[section, lines] : sections.lines) {
      (void)lines;
      if (!section_outputs.empty()) section_outputs += ',';
      section_outputs += options.param_filename + "." + section;
    }
    report.value("output.parameter_sections", section_outputs.empty() ? "none" : section_outputs);
  }
  if (options.diag_seed_only) {
    report.value("output.seed_artifacts", "ham*,val,vec*");
  } else if (options.generate_temporaries) {
    report.value("output.legacy_artifacts", "parameter sections, Wilson tables, ham*, val, vec*");
  } else {
    report.resolved("output.legacy_artifacts", "inactive", "--generate-temporaries not selected");
  }
  report.write(std::cerr);
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
                        const std::filesystem::path &template_dir, const bool write_legacy_outputs,
                        const std::filesystem::path &output_dir = ".") {
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
      save_text_matrix(output_dir / ("ham." + qn_name), matrix);
      save_text_matrix(output_dir / "ham", matrix);
    }

    auto raw = NRG::diagonalise_dsyev(matrix);
    for (auto &value : raw.val) value *= seed.factor;
    if (!raw.val.empty()) seed.smallest = std::min(seed.smallest, raw.val.front());

    if (write_legacy_outputs) {
      save_values(output_dir / "val", raw.val);
      save_binary_matrix(output_dir / "vec", raw.vec);
      save_binary_matrix(output_dir / ("vec." + qn_name), raw.vec);
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
    const auto intermediate = NRG::matrix_prod<double>(matrix, NRG::trans(vec2_it->second));
    NRG::Matrix_traits<double> transformed = NRG::matrix_prod<double>(vec1_it->second, intermediate);
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

std::unique_ptr<NRG::Params> make_instantiate_params(const std::string &param_filename,
                                                     const std::filesystem::path &workdir_parent = ".") {
  ScopedCoutRedirect quiet;
  auto workdir = std::make_unique<NRG::Workdir>(workdir_parent.string(), true);
  auto params = std::make_unique<NRG::Params>(param_filename, "param", std::move(workdir), true, true);
  params->silent = true;
  return params;
}

void initialize_template_symmetry(NRG::Params &params, const TemplateHeader &header) {
  if (params.symtype.value().empty()) throw std::runtime_error("Parameter symtype must be set before parsing template/data.in.");
  ScopedCoutRedirect quiet;
  const auto symmetry = NRG::set_symmetry<double>(params, params.symtype.value(), header.channels);
  (void)symmetry;
}

void validate_generated_data(const std::string &param_filename, const std::string &data_text,
                             const std::filesystem::path &workdir_parent = ".") {
  std::istringstream data_in(data_text);
  ScopedCoutRedirect quiet;
  auto params = make_instantiate_params(param_filename, workdir_parent);
  const NRG::InputData<double> input(*params, data_in);
}

void register_full_legacy_artifacts(ArtifactManifest &artifacts, const std::filesystem::path &staging_dir,
                                    const SeedData &seed) {
  const std::vector<std::string> filenames = {
      "xi.dat",     "zeta.dat",  "theta.dat", "de_pos.dat", "de_neg.dat", "du_pos.dat", "du_neg.dat",
      "xi1.dat",    "zeta1.dat", "theta1.dat", "ham",        "val",        "vec",
  };
  for (const auto &filename : filenames) register_existing_artifact(artifacts, staging_dir, filename);
  for (const auto &subspace : seed.subspaces) {
    const auto qn_name = legacy_qn_name(subspace.qn);
    register_existing_artifact(artifacts, staging_dir, "ham." + qn_name);
    register_existing_artifact(artifacts, staging_dir, "vec." + qn_name);
  }
}

void run_diag_seed_only(const Options &options) {
  auto sections = parse_param_sections(options.param_filename);
  const auto nrgchain = make_nrgchain_invocation(options.param_filename, NrgChainTableMode::Calculate);
  report_configuration(options, nrgchain);
  write_param_section_files(sections, options.param_filename);
  auto params = make_instantiate_params(options.param_filename);

  const auto wilson = NRG::Tools::NrgChain::calculate_from_params(nrgchain.parameters, nrgchain.configuration.mode);
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
  const auto nrgchain = make_nrgchain_invocation(options.param_filename, NrgChainTableMode::Calculate);
  report_configuration(options, nrgchain);
  const auto data_destination = std::filesystem::absolute("data");
  NRG::Workdir staging_workdir(data_destination.parent_path().string(), true);
  const auto staging_dir = std::filesystem::path(staging_workdir.get());
  ArtifactManifest legacy_artifacts;
  if (options.generate_temporaries)
    write_staged_param_section_files(sections, options.param_filename, staging_dir, legacy_artifacts);
  auto params = make_instantiate_params(options.param_filename, staging_dir);

  const auto wilson = NRG::Tools::NrgChain::calculate_from_params(
      nrgchain.parameters, nrgchain.configuration.mode, staging_dir);
  if (wilson.channels.size() != 1)
    throw std::runtime_error("Only single-channel Wilson generation is supported in this instantiate slice.");
  if (options.generate_temporaries) write_wilson_channel(wilson.channels.front(), 1, staging_dir);

  const auto nmax = parse_size_t_value(param_value(sections, "param", "Nmax"), "Nmax");
  const auto polarized = params->polarized.value();
  if (polarized) throw std::runtime_error("Spin-polarized coefficient tables are not supported in this instantiate slice.");

  auto evaluator = make_matrix_evaluator(sections, wilson);

  std::ostringstream data_buffer;

  const auto template_dir = std::filesystem::path(options.template_dir);
  DataTemplateReader data_in(resolve_template_file(template_dir, "data.in"), &data_buffer);
  const auto header = read_template_header(data_in);
  if (header.channels != 1) throw std::runtime_error("Only single-channel data.in templates are supported in this slice.");
  initialize_template_symmetry(*params, header);
  data_buffer << header.channels << ' ' << nmax << ' ' << header.subspaces << '\n';

  const auto seed = read_seed_data(data_in, header.subspaces, evaluator, template_dir, options.generate_temporaries, staging_dir);
  write_seed_energy_block(data_buffer, seed);
  const auto tail = process_template_tail(data_in, data_buffer, seed, evaluator, template_dir, params->channels * params->perchannel);
  if (tail.ground_energy != seed.ground_energy) throw std::runtime_error("Internal ground-energy mismatch while instantiating data.");
  write_z_coefficients(data_buffer, wilson, params->coefchannels, nmax);

  const auto data_text = data_buffer.str();
  validate_generated_data(options.param_filename, data_text, staging_dir);

  const auto staged_data = staging_dir / "data";
  write_text_file_checked(staged_data, data_text);
  if (options.generate_temporaries) {
    register_full_legacy_artifacts(legacy_artifacts, staging_dir, seed);
    for (const auto &[destination, staged] : legacy_artifacts) publish_file(staged, destination);
  }
  publish_file(staged_data, data_destination);

  std::cout << "E_gs=" << std::setprecision(18) << seed.ground_energy << '\n';
}

void run_wilson_only(const Options &options) {
  const auto nrgchain = make_nrgchain_invocation(options.param_filename, NrgChainTableMode::Calculate);
  report_configuration(options, nrgchain);
  const auto wilson = NRG::Tools::NrgChain::calculate_from_params(nrgchain.parameters, nrgchain.configuration.mode);
  if (wilson.channels.size() != 1)
    throw std::runtime_error("Only single-channel Wilson generation is supported in this instantiate slice.");
  write_wilson_channel(wilson.channels.front(), 1);
}

} // namespace

int main(int argc, char *argv[]) {
  if (NRG::Tools::report_version_if_requested(argc, argv, "instantiate")) return EXIT_SUCCESS;
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
