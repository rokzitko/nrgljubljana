#include <nrgljubljana/solver_preparation.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

namespace fs = std::filesystem;
using Entries = std::map<std::string, std::string>;
unsigned assertions = 0;

void check(bool condition, const std::string &message) {
  ++assertions;
  if (!condition) throw std::runtime_error(message);
}

template <typename Function>
void rejects_call(Function function, const std::string &detail) {
  try {
    function();
  } catch (const std::runtime_error &error) {
    check(std::string(error.what()).find(detail) != std::string::npos,
          "Expected error containing '" + detail + "', got '" + error.what() + "'");
    return;
  }
  throw std::runtime_error("Expected rejection containing: " + detail);
}

// Fixture shorthand only. Production callers select [param] and detect duplicate
// keys/sections before constructing the map; neither is recoverable from a map.
Entries settings(const std::string &text) {
  Entries result;
  std::istringstream input(text);
  for (std::string line; std::getline(input, line);) {
    const auto equals = line.find('=');
    check(equals != std::string::npos, "Invalid parameter fixture: " + line);
    check(result.emplace(line.substr(0, equals), line.substr(equals + 1)).second, "Duplicate fixture key: " + line);
  }
  return result;
}

void resolves(const std::string &text, unsigned expected) {
  const auto parameters = settings(text);
  const auto before = parameters;
  const auto actual = NRG::Workflow::resolve_chain_length(parameters);
  check(actual == expected, "Expected Nmax=" + std::to_string(expected) + ", got " + std::to_string(actual) + " for:\n" + text);
  check(NRG::Workflow::resolve_chain_length(parameters) == actual && parameters == before, "Resolution changed input/state");
}

void rejects_chain(const std::string &text, const std::string &detail) {
  const auto parameters = settings(text);
  const auto before = parameters;
  rejects_call([&] { NRG::Workflow::resolve_chain_length(parameters); }, detail);
  check(parameters == before, "Failed resolution changed input");
}

std::string number(double value) {
  std::ostringstream text;
  text.imbue(std::locale::classic());
  text << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
  return text.str();
}

const std::string y_scale = "Lambda=4\nz=1\nbandrescale=1\ndiscretization=Y\n";
const std::string ratio_cutoff = "T=0.3125\nTmin_ratio=0.5\n" + y_scale;
const std::string missing = "Specify Nmax, a positive Tmin, or positive T and Tmin_ratio";

void chain_fixtures() {
  // Scalar results from the DMFT native and Perl instantiate tests. Preserve all
  // valid results; only ignored-value evaluation with explicit Tmin is corrected.
  resolves("Nmax=7\n", 7);
  for (const auto &[z, expected] : {std::pair{"0.25", 54u}, {"0.50", 54u}, {"0.75", 53u}, {"1.00", 53u}})
    resolves("Tmin=1e-7\nLambda=2\nz=" + std::string(z) + "\nbandrescale=10\ndiscretization=Z\n", expected);
  resolves("Tmin=0.15625\n" + y_scale, 3);
  resolves("Tmin=0.34756349141427556\nLambda=1.1\nz=0.1\nbandrescale=1\ndiscretization=Y\n", 24);
  resolves("Nmax=999\n" + ratio_cutoff, 3);
  resolves("T=1\nTmin_ratio=1\nTmin=0.15625\n" + y_scale, 3);
  resolves("Nmax=7\nT=1\nTmin_ratio=0\n", 7);
  rejects_chain("Nmax=7\nTmin=1e-3\n", "Specify either Nmax or Tmin, not both!");
  rejects_chain("Tmin=invalid\n", "Tmin must be a finite number");
  rejects_chain("Nmax=7\nT=1e-9999\nTmin_ratio=1\n", "Numeric underflow in T");
  rejects_chain("Nmax=999\n", "range 1..998");
  rejects_chain("Tmin=1e-3\nLambda=1\n", "Lambda must be greater than 1");
  rejects_chain("Tmin=2\n" + y_scale, "resolves to Nmax=0");
  rejects_chain("Tmin_ratio=1\n", missing);
  resolves("Nmax=7\nunknown=model_expression\n", 7);
}

void chain_precedence() {
  // Contract: mutual exclusion first, then Nmax syntax, then explicit Tmin OR
  // the paired ratio inputs. An overridden integer needs no representable range.
  for (const auto *nmax : {"0", "-0", "-7", "999", "1000", "999999999999999999999999999999999999999",
                           "-999999999999999999999999999999999999999"})
    resolves("Nmax=" + std::string(nmax) + '\n' + ratio_cutoff, 3);
  for (const auto *tmin : {"0", "-0", "-1"}) {
    rejects_chain("Tmin=" + std::string(tmin) + '\n' + ratio_cutoff, missing);
    rejects_chain("Tmin=" + std::string(tmin) + '\n', missing);
    rejects_chain("Tmin=" + std::string(tmin) + "\nT=invalid\nTmin_ratio=invalid\n", missing);
  }
  for (const auto *inactive : {"0", "-0", "-1", "0e-9999", "-0.000e+9999"}) {
    resolves("Nmax=7\nT=" + std::string(inactive) + "\nTmin_ratio=1\n", 7);
    resolves("Nmax=7\nT=1\nTmin_ratio=" + std::string(inactive) + '\n', 7);
    rejects_chain("T=" + std::string(inactive) + "\nTmin_ratio=1\n", missing);
    rejects_chain("T=1\nTmin_ratio=" + std::string(inactive) + '\n', missing);
  }
  for (const auto *unpaired : {"T=invalid\n", "Tmin_ratio=invalid\n"}) {
    resolves("Nmax=7\n" + std::string(unpaired), 7);
    resolves("Tmin=0.15625\n" + std::string(unpaired) + y_scale, 3);
    rejects_chain(unpaired, missing);
  }
  rejects_chain("", missing);
  resolves("Nmax=7\nLambda=invalid\nz=invalid\nbandrescale=invalid\ndiscretization=invalid\n", 7);
  rejects_chain("Nmax=invalid\nTmin=invalid\n", "Specify either Nmax or Tmin, not both!");
  rejects_chain("Nmax=7\nTmin=0\n" + ratio_cutoff, "Specify either Nmax or Tmin, not both!");
  rejects_chain("Nmax=invalid\nT=invalid\nTmin_ratio=invalid\n", "Nmax must be an integer");
  rejects_chain("Nmax=7\nT=0\nTmin_ratio=invalid\n", "Tmin_ratio must be a finite number");
  rejects_chain("Nmax=7\nT=invalid\nTmin_ratio=0\n", "T must be a finite number");

  // Reconciled nrginit precedence: ignored ratio inputs must not be parsed or
  // multiplied, even if malformed, nonfinite, underflowing or overflowing.
  for (const auto *ignored : {"T=invalid\nTmin_ratio=0\n", "T=0\nTmin_ratio=invalid\n",
                              "T=nan\nTmin_ratio=inf\n", "T=1e-9999\nTmin_ratio=1e9999\n",
                              "T=1e308\nTmin_ratio=1e308\n", "T=1e-300\nTmin_ratio=1e-300\n"}) {
    resolves("Tmin=0.15625\n" + std::string(ignored) + y_scale, 3);
    rejects_chain("Tmin=invalid\n" + std::string(ignored), "Tmin must be a finite number");
  }
}

void chain_validation() {
  resolves("Nmax=1\n", 1);
  resolves("Nmax=998\n", 998);
  resolves("Nmax=000000000000000000000000000000000000000007\n", 7);
  resolves("Nmax=000998\n", 998);
  resolves("Nmax= \t7 \t\r\n", 7);
  for (const auto *nmax : {"", "+7", "7.0", "7e0", "0x7", "7junk", "7 # comment", "7 0", "--7", "-", "nan", "inf"}) {
    rejects_chain("Nmax=" + std::string(nmax) + '\n', "Nmax must be an integer:");
    rejects_chain("Nmax=" + std::string(nmax) + '\n' + ratio_cutoff, "Nmax must be an integer:");
  }
  for (const auto *nmax : {"0", "0000", "-0", "-0007", "-7", "999", "000999", "1000",
                           "999999999999999999999999999999999999999"})
    rejects_chain("Nmax=" + std::string(nmax) + '\n', "range 1..998");

  resolves("Tmin=+.15625E+0\n" + y_scale, 3);
  resolves("Tmin=15625.e-5\n" + y_scale, 3);
  resolves("Tmin= \t.15625\t\r\n" + y_scale, 3);
  for (const auto *name : {"Tmin", "T", "Tmin_ratio", "Lambda", "z", "bandrescale"}) {
    const std::string key = name;
    const std::string rest = key == "Tmin" ? "" : key == "T" ? "Nmax=7\nTmin_ratio=1\n" :
                             key == "Tmin_ratio" ? "Nmax=7\nT=1\n" : "Tmin=1e-3\n";
    for (const auto *value : {"", "nan", "NaN", "Inf", "-inf", "0x1", "1/2", "1 0", "1e", ".", "1 # comment",
                              "EPS", "ZZZ", "1e9999", "-1e9999"})
      rejects_chain(rest + key + '=' + value + '\n', key + " must be a finite number");
    for (const auto *value : {"1e-9999", "-1e-9999"})
      rejects_chain(rest + key + '=' + value + '\n', "Numeric underflow in " + key);
  }
  for (const auto *lambda : {"1", "0", "-2"})
    rejects_chain("Tmin=1e-3\nLambda=" + std::string(lambda) + '\n', "Lambda must be greater than 1");
  for (const auto *z : {"0", "-1", "1.0000000000000002"})
    rejects_chain("Tmin=1e-3\nz=" + std::string(z) + '\n', "z must satisfy 0 < z <= 1");
  for (const auto *band : {"0", "-1"})
    rejects_chain("Tmin=1e-3\nbandrescale=" + std::string(band) + '\n', "bandrescale must be positive");
  for (const auto *scheme : {"", "y", "c", "z", "X", "1"})
    rejects_chain("Tmin=1e-3\ndiscretization=" + std::string(scheme) + '\n', "Unknown discretization");

  for (const auto *product : {"T=1e308\nTmin_ratio=1e308\n", "T=1e-300\nTmin_ratio=1e-300\n"}) {
    rejects_chain(product, "T*Tmin_ratio must be a finite representable number");
    rejects_chain("Nmax=7\n" + std::string(product), "T*Tmin_ratio must be a finite representable number");
  }
  rejects_chain("Tmin=1\nLambda=1e308\nz=0.5\nbandrescale=1e308\n", "Wilson scale is not finite");
}

void chain_scales() {
  // Defaults and first-character scheme selection are intentional scalar legacy
  // behavior; validation of unused Lambda/z/bandrescale/scheme is not performed.
  resolves("Tmin=0.15625\n", 5);
  for (const auto *scheme : {"C", "Z", "Canything", "Zanything"})
    resolves("Tmin=0.15\nLambda=4\ndiscretization=" + std::string(scheme) + '\n', 2);
  for (const auto *scheme : {"Y", "Yanything"})
    resolves("Tmin=0.15\nLambda=4\ndiscretization=" + std::string(scheme) + '\n', 3);
  resolves("Tmin=1.5625\nLambda=4\nz=1\nbandrescale=10\ndiscretization=Y\n", 3);
  resolves("Tmin=0.15625\nLambda=4\nz=0.5\nbandrescale=1\ndiscretization=Y\n", 4);

  // Exact powers of two pin both sides of inclusive cutoffs, including 1/998/999.
  // Do not replace the loop with logarithmic inversion or reassociate products.
  for (const unsigned n : {1u, 3u, 25u, 998u}) {
    const double cutoff = std::ldexp(0.625, -static_cast<int>(n - 1));
    resolves("Tmin=" + number(cutoff) + '\n' + y_scale, n);
    resolves("Tmin=" + number(std::nextafter(cutoff, 0.0)) + '\n' + y_scale, n);
    const auto above = "Tmin=" + number(std::nextafter(cutoff, 1.0)) + '\n' + y_scale;
    if (n == 1) rejects_chain(above, "resolves to Nmax=0");
    else resolves(above, n - 1);
  }
  const double cutoff999 = std::ldexp(0.625, -998);
  rejects_chain("Tmin=" + number(cutoff999) + '\n' + y_scale, "Nmax >= 999");
  resolves("Tmin=" + number(std::nextafter(cutoff999, 1.0)) + '\n' + y_scale, 998);
  rejects_chain("Tmin=1e-300\n", "Nmax >= 999");

  // Representable subnormals are valid; nonzero scalar inputs rounded to zero
  // are not. Scale underflow itself can legitimately end the counting loop.
  const auto tiny = number(std::numeric_limits<double>::denorm_min());
  resolves("Tmin=" + tiny + "\nLambda=4\nz=1\nbandrescale=" + tiny + "\ndiscretization=Y\n", 1);
  resolves("T=" + tiny + "\nTmin_ratio=1\nLambda=4\nz=1\nbandrescale=" + tiny + "\ndiscretization=Y\n", 1);
  resolves("Tmin=1e-3\nz=" + tiny + '\n', 21);
  rejects_chain("Tmin=" + tiny + "\nLambda=16\nz=1\nbandrescale=" + tiny + "\ndiscretization=Z\n", "resolves to Nmax=0");
}

void preparation_parameters(const std::string &text, const std::string &error = "") {
  const auto parameters = settings(text);
  const auto before = parameters;
  const auto out_flags = std::cout.flags(), err_flags = std::cerr.flags();
  const auto out_precision = std::cout.precision(), err_precision = std::cerr.precision();
  std::ostringstream output, diagnostics;
  auto *out = std::cout.rdbuf(output.rdbuf());
  auto *err = std::cerr.rdbuf(diagnostics.rdbuf());
  try {
    for (unsigned repeat = 0; repeat < 2; ++repeat) {
      if (error.empty()) NRG::Workflow::validate_preparation_parameters(parameters);
      else rejects_call([&] { NRG::Workflow::validate_preparation_parameters(parameters); }, error);
    }
  } catch (...) {
    std::cout.rdbuf(out);
    std::cerr.rdbuf(err);
    throw;
  }
  std::cout.rdbuf(out);
  std::cerr.rdbuf(err);
  check(parameters == before, "Parameter validation changed input");
  check(output.str().empty() && diagnostics.str().empty(), "Parameter validation emitted diagnostics");
  check(std::cout.flags() == out_flags && std::cerr.flags() == err_flags &&
        std::cout.precision() == out_precision && std::cerr.precision() == err_precision,
        "Parameter validation changed stream formatting");
}

void combined_preparation() {
  // Combined discretize(integral=true) + full-instantiate static contract, not
  // the chain-only contract (which still permits xmax=1 and rescalexi=true).
  // These checks use the actual initializer resolver and adapt's shared parser.
  const std::string base = "Nmax=7\n";
  preparation_parameters(base);
  preparation_parameters(base + "xmax=1\n", "xmax must be greater than 1");
  preparation_parameters(base + "xmax=" + number(std::nextafter(1.0, 2.0)) + '\n');
  preparation_parameters(base + "xfine=1\n", "xfine must be greater than 1");
  preparation_parameters(base + "xfine=" + number(std::nextafter(1.0, 2.0)) + '\n');
  for (const auto *value : {"true", "1", "YES"})
    preparation_parameters(base + "rescalexi=" + value + '\n', "rescalexi=true is not supported for full instantiation");
  for (const auto *value : {"false", "0", "NO"}) preparation_parameters(base + "rescalexi=" + value + '\n');

  // Effective Nmax is authoritative here, unlike resolve_chain_length's input
  // map. A rendered twist must not be resolved again against stale cutoff keys.
  preparation_parameters("Nmax=1\n");
  preparation_parameters("Nmax=998\n");
  preparation_parameters(base + "Tmin=invalid\nT=invalid\nTmin_ratio=invalid\n");
  preparation_parameters("Nmax=999\n" + ratio_cutoff, "resolved value in [1,998]");
  preparation_parameters("Nmax=0\nmMAX=1\n", "resolved value in [1,998]");
  preparation_parameters("Nmax=-1\n", "Nmax must be nonnegative");
  preparation_parameters("Nmax=7.0\n", "Nmax must be an integer");
  preparation_parameters("Tmin=0.15625\nT=invalid\nTmin_ratio=invalid\n" + y_scale);
  preparation_parameters(ratio_cutoff);
  preparation_parameters("", missing);
  preparation_parameters(base + "z=ZZZ\n", "z must be a finite representable number");
  for (const auto *z : {"0.25", "0.50", "0.75", "1"}) preparation_parameters(base + "z=" + z + '\n');
  for (const auto *z : {"0", "-1", "1.0000000000000002"})
    preparation_parameters(base + "z=" + z + '\n', "z must satisfy 0 < z <= 1");

  for (const auto *key : {"Lambda", "z", "bandrescale", "xmax", "xfine", "outputstep", "dx_fine", "dx_fast",
                          "allowed_error", "max_abs", "secant_eps", "secant_factor", "boundary"}) {
    for (const auto *value : {"nan", "inf", "invalid", "1e9999", "1e-9999"})
      preparation_parameters(base + key + '=' + value + '\n', std::string(key) + " must be a finite representable number");
  }
  for (const auto *lambda : {"1", "0", "-2"})
    preparation_parameters(base + "Lambda=" + lambda + '\n', "Lambda must be greater than 1");
  for (const auto *key : {"bandrescale", "dx_fine", "dx_fast", "allowed_error", "max_abs", "secant_eps", "secant_factor"}) {
    for (const auto *value : {"0", "-1"})
      preparation_parameters(base + key + '=' + value + '\n', key);
  }
  preparation_parameters(base + "outputstep=0\n", "outputstep must be finite and in the range (0, 1]");
  preparation_parameters(base + "outputstep=1\n");
  preparation_parameters(base + "outputstep=" + number(std::nextafter(1.0, 2.0)) + '\n', "outputstep must be finite");
  for (const auto *key : {"max_subdiv", "secant_max_iter", "mMAX", "preccpp"})
    preparation_parameters(base + key + "=1.5\n", std::string(key) + " must be an integer");
  preparation_parameters(base + "max_subdiv=-1\n", "max_subdiv must be nonnegative");
  preparation_parameters(base + "max_subdiv=0\n");
  preparation_parameters(base + "secant_max_iter=0\n", "secant_max_iter must be positive");
  preparation_parameters(base + "secant_max_iter=1\n");

  // Pin the strict sub-step boundary and the inclusive default integral tolerance
  // boundary. CQUAD validation is pure; this API must not install a GSL handler.
  const double epsilon = std::numeric_limits<double>::epsilon();
  const double step = std::ldexp(epsilon, 10);
  preparation_parameters(base + "dx_fine=" + number(step) + '\n', "sub-step below machine precision");
  preparation_parameters(base + "dx_fine=" + number(std::nextafter(step, 1.0)) + '\n');
  preparation_parameters(base + "allowed_error=" + number(epsilon) + '\n');
  preparation_parameters(base + "allowed_error=" + number(std::nextafter(epsilon, 0.0)) + '\n', "Relative CQUAD tolerance is too small");

  preparation_parameters(base + "mMAX=0\n", "mMAX must be greater than 0");
  preparation_parameters("Nmax=1\nmMAX=1\n");
  preparation_parameters(base + "preccpp=10\n", "preccpp must be greater than 10");
  preparation_parameters(base + "preccpp=11\n");
  preparation_parameters(base + "tridiag_method=rkpw\npreccpp=0\n");
  preparation_parameters(base + "tridiag_method=rkpw\npreccpp=-1\n", "preccpp must be nonnegative");
  preparation_parameters(base + "tridiag_method=unknown\n", "Unknown tridiag_method");
  preparation_parameters(base + "tri=cpp\n", "tri=cpp is not supported");
  preparation_parameters(base + "nrgchain_tridiag=false\n", "instantiate requires nrgchain_tridiag=true");
  preparation_parameters(base + "nrgchains_tridiag=false\n", "instantiate requires nrgchain_tridiag=true");
  preparation_parameters(base + "nrgchain_tridiag=true\nnrgchains_tridiag=ignored\n");
  preparation_parameters(base + "nrgchain_tridiag=false\nnrgchains_tridiag=true\n", "instantiate requires nrgchain_tridiag=true");
  preparation_parameters(base + "nrgchain_tables_load=true\n", "does not support nrgchain_tables_load=true");
  preparation_parameters(base + "nrgchain_tables_save=true\n");
  preparation_parameters(base + "nrgchain_tables_load=true\nnrgchain_tables_save=true\n", "cannot both be true");
  for (const auto *key : {"adapt", "hardgap", "rescalexi", "nrgchain_tables_save", "nrgchain_tables_load",
                          "nrgchain_tridiag", "nrgchains_tridiag"})
    preparation_parameters(base + key + "=invalid\n", "Invalid boolean value");

  // Leave DMFT policy to the caller. Adaptive meshes remain legal except for the
  // native nonzero-hardgap/flat-band combinations; no fdm or QS checks live here.
  preparation_parameters(base + "adapt=true\nfdm=false\nsymtype=app-owned\n");
  preparation_parameters(base + "adapt=true\nloadg=invalid\n", "Invalid boolean value");
  preparation_parameters(base + "adapt=false\nloadg=ignored\n");
  preparation_parameters(base + "adapt=true\nloadg=true\ndos=/not/read/by/preflight\n");
  preparation_parameters(base + "band=flat\nadapt=true\n", "analytic flat-band generation require adapt=false");
  preparation_parameters(base + "band=flat\nadapt=false\n");
  preparation_parameters(base + "hardgap=true\nboundary=0\nadapt=true\n");
  preparation_parameters(base + "hardgap=true\nboundary=0.5\nadapt=true\n", "Nonzero hardgap supports only a fixed mesh");
  preparation_parameters(base + "hardgap=true\nboundary=0.5\nf_method=ignored\n");
  preparation_parameters(base + "hardgap=true\nboundary=" + number(std::nextafter(1.0, 0.0)) + '\n');
  for (const auto *boundary : {"-0.1", "1"})
    preparation_parameters(base + "hardgap=true\nboundary=" + boundary + '\n', "boundary must be in [0,1)");
  preparation_parameters(base + "hardgap=false\nboundary=2\n");
  for (const auto *method : {"linear", "steffen"}) preparation_parameters(base + "density_interpolation=" + method + '\n');
  preparation_parameters(base + "density_interpolation=cspline\n", "density_interpolation must be either 'linear' or 'steffen'");
}

Entries minimal() {
  return {{"data.in", "# Input file for NRG Ljubljana\n# symtype  QS\n#!9\n1 0 1\n"
                     "# SCALE 1\n0 1\n1\nDIAG ham_0.1\nf 0 0\n0\ne\n0\n"},
          {"ham_0.1", "{{0}}\n"}};
}

void replace(Entries &files, const std::string &from, const std::string &to, const std::string &path = "data.in") {
  auto &text = files.at(path);
  const auto pos = text.find(from);
  check(pos != std::string::npos, "Invalid test substitution: " + from);
  text.replace(pos, from.size(), to);
}

void accepts(const Entries &files, unsigned subspaces = 1) {
  const auto before = files;
  for (unsigned i = 0; i < 2; ++i) {
    const auto info = NRG::Workflow::inspect_template_bundle(files);
    check(info.symmetry == "QS" && info.channels == 1 && info.subspaces == subspaces, "Incorrect template metadata");
    check(files == before, "Inspection changed input bytes");
  }
}

void rejects_template(const Entries &files, const std::string &description) {
  const auto before = files;
  try {
    rejects_call([&] { NRG::Workflow::inspect_template_bundle(files); }, "NRG template: ");
  } catch (const std::exception &error) {
    throw std::runtime_error(description + ": " + error.what());
  }
  check(files == before, "Failed inspection changed input bytes");
}

void template_structure() {
  // The native version/symmetry header is intentionally byte-sensitive; the body
  // follows DataTemplateReader. This is not a permissive Mathematica file loader.
  const auto valid = minimal();
  accepts(valid);
  rejects_template({}, "missing data.in");
  const auto &text = valid.at("data.in");
  for (std::size_t end = 0; end < text.size(); ++end) {
    if (end && text[end - 1] != '\n') continue;
    auto files = valid;
    files["data.in"] = text.substr(0, end);
    rejects_template(files, "truncation at byte " + std::to_string(end));
  }
  auto no_newline = valid;
  no_newline["data.in"].pop_back();
  accepts(no_newline);
  for (const auto &[from, to] : std::vector<std::pair<std::string, std::string>>{
         {"#!9", "#!8"}, {"#!9", "#!"}, {"# symtype  QS", "# symtype  QSZ"},
         {"# symtype  QS", "# symtype  QS "}, {"# symtype  QS", "# symtype  QS\r"},
         {"#!9", "#!9\r\n# another header comment"}, {"#!9", "#!9 \n# another header comment"},
         {"# Input file", "# COMPLEX Input file"}, {"# Input file", "\n# Input file"},
         {"# Input file", " # Input file"}, {"1 0 1", "2 0 1"}, {"1 0 1", "1 0 0"},
         {"1 0 1", "1 -1 1"}, {"1 0 1", "1 0 1 extra"}, {"1 0 1", "1 0 +1"},
         {"1 0 1", "1 0 184467440737095516160"}, {"1 0 1", "1 0 4294967296"}, {"1 0 1", "1 0 2"},
         {"# SCALE 1\n", ""}, {"# SCALE 1", "# SCALE 0"}, {"# SCALE 1", "# SCALE -1"},
         {"# SCALE 1", "# SCALE nan"}, {"# SCALE 1", "# SCALE 1e-320"}, {"# SCALE 1", "# SCALE 1e9999"},
         {"# SCALE 1", "# SCALE"}, {"# SCALE 1", "# SCALE 1e-9999"},
         {"0 1\n1\nDIAG", "0\n1\nDIAG"}, {"0 1\n1\nDIAG", "0 1 2\n1\nDIAG"},
         {"0 1\n1\nDIAG", "0 0\n1\nDIAG"}, {"0 1\n1\nDIAG", "0 -1\n1\nDIAG"},
         {"0 1\n1\nDIAG", "2147483648 1\n1\nDIAG"},
         {"0 1\n1\nDIAG", "0 1\n0\nDIAG"}, {"0 1\n1\nDIAG", "0 1\n-1\nDIAG"},
         {"DIAG ham_0.1", "0"}, {"DIAG ham_0.1", "DIAG"}, {"DIAG ham_0.1", "DIAG ham_0.1 extra"},
         {"DIAG ham_0.1", "DIAG ham_0.1 # comment"}, {"f 0 0\n0\n", ""},
         {"f 0 0", "f"}, {"f 0 0", "f 0 0 extra"}, {"f 0 0\n0", "f 0 0\n-1"},
         {"f 0 0\n0", "f 0 0\n1"}, {"e\n0\n", ""}, {"e\n0", "e"},
         {"e\n0", "e\n0\nf 0 0\n0"}, {"f 0 0\n0\ne\n0", "e\n0\nf 0 0\n0"},
         {"e\n0", "e\n0\nvname\n0"}, {"e\n0", "e\n0\nunknown\n0"}}) {
    auto files = valid;
    replace(files, from, to);
    rejects_template(files, from + " -> " + to);
  }
  auto duplicate = valid;
  replace(duplicate, "1 0 1", "1 0 2");
  replace(duplicate, "f 0 0", "0 1\n1\nDIAG ham_0.1\nf 0 0");
  rejects_template(duplicate, "duplicate subspace");
  auto nul = valid;
  nul["data.in"].insert(0, 1, '\0');
  rejects_template(nul, "embedded NUL");

  // CRLF is harmless only in these positions. Placeholder energy is discarded;
  // Wilson labels are unsigned but ignored by the native Opch reader.
  for (const auto &[from, to] : std::vector<std::pair<std::string, std::string>>{
         {"# Input file for NRG Ljubljana\n", "# Input file for NRG Ljubljana\r\n"},
         {"#!9\n", "#!9\r\n"}, {"1 0 1\n", "1 0 1\r\n"},
         {"# SCALE 1\n", "\n  # SCALE 1\r\n\t\r\n"},
         {"DIAG ham_0.1\n", "  DIAG ham_0.1\r\n"}, {"f 0 0", "f 7 8"},
         {"e\n0", "e\nunused-placeholder"}, {"e\n0", "e\n0\ne\nunused-placeholder"}}) {
    auto files = valid;
    replace(files, from, to);
    accepts(files);
  }
  auto optional = valid;
  for (const char type : std::string("spdtogq")) optional["data.in"] += std::string(1, type) + "name\n0\n";
  accepts(optional);
  for (const char marker : std::string("zZT")) {
    auto files = optional;
    files["data.in"] += std::string(1, marker) + "\nignored old coefficient tables\n!not-executed\n";
    accepts(files);
    files = valid;
    replace(files, "e\n0", std::string(1, marker) + "\ne\n0");
    rejects_template(files, "e block only in discarded suffix");
  }
}

void template_dependencies() {
  const auto valid = minimal();
  for (const std::string name : {"/etc/passwd", "../ham_0.1", "sub/../ham_0.1", "../template/ham_0.1",
                                 "missing", "ham_0.1/", "ham_0.1/."}) {
    auto files = valid;
    replace(files, "DIAG ham_0.1", "DIAG " + name);
    rejects_template(files, "unsafe or absent DIAG reference " + name);
  }
  auto files = valid;
  files.erase("ham_0.1");
  rejects_template(files, "missing DIAG matrix");
  for (const auto &[name, canonical] : std::vector<std::pair<std::string, std::string>>{
         {"./ham_0.1", "ham_0.1"}, {"sub/./ham", "sub/ham"}, {"sub//ham", "sub/ham"}}) {
    files = valid;
    replace(files, "DIAG ham_0.1", "DIAG " + name);
    files[canonical] = "{{0}}";
    accepts(files);
  }
  files = valid;
  replace(files, "DIAG ham_0.1", "DIAG sub/ham");
  files["sub/ham"] = "{{0}}";
  files["sub"] = "not a directory";
  rejects_template(files, "file used as matrix parent");

  // There is no special 'template/' namespace in this API. The app must strip
  // its import prefix; a real nested directory named template is still legal.
  Entries prefixed;
  for (const auto &[name, bytes] : valid) prefixed["template/" + name] = bytes;
  rejects_template(prefixed, "application import keys are not root-relative");
  files = valid;
  files.erase("ham_0.1");
  files["template/ham_0.1"] = "{{0}}";
  rejects_template(files, "no implicit prefix lookup");
  replace(files, "DIAG ham_0.1", "DIAG template/ham_0.1");
  accepts(files);

  for (const std::string key : {"", "/absolute", "../outside", "./alias", "a//b", "a/../b", "a/", "a/."}) {
    files = valid;
    files[key] = "unused";
    rejects_template(files, "noncanonical bundle key");
  }
  files = valid;
  files[std::string("a\0b", 3)] = "unused";
  rejects_template(files, "NUL in bundle key");

  files = valid;
  files["data.in"] += "soperator\n1\n0 1 0 1\nop.matrix\n";
  rejects_template(files, "missing operator matrix");
  files["template/op.matrix"] = "{{0}}";
  rejects_template(files, "operator matrix in wrong namespace");
  files["op.matrix"] = "{{0}}";
  accepts(files);
  for (const std::string name : {"ops/../../ham_0.1", "op/../ham_0.1", "/etc/passwd", "../ham_0.1"}) {
    auto unsafe = files;
    replace(unsafe, "op.matrix\n", name + '\n');
    rejects_template(unsafe, "unsafe operator matrix " + name);
  }
  replace(files, "op.matrix\n", "operators/nested matrix\n");
  files["operators/nested matrix"] = "{{0}}";
  accepts(files);
  replace(files, "operators/nested matrix\n", "./operators/nested matrix\n");
  rejects_template(files, "operator filenames must begin with o, not a guessed extension");
}

void template_matrices() {
  const auto valid = minimal();
  for (const std::string matrix : {"0", "{0}", "{{0}}", "# comment\n{{0}} # tail\n", "{{eps1 + U1}}",
        "{{Sqrt[2]*gammaPolCh[1] + coefzeta(1, 0)}}", "{{-(-1.25e-3) / +0x1p2}}",
        "{ # row\n Conjugate[eps1] + exp(log(2)) / Pi }", "{{undefined / 0}}", "sqrt(2) + coefxi[(1), +0]"}) {
    auto files = valid;
    files["ham_0.1"] = matrix;
    accepts(files);
  }
  for (const std::string matrix : {"", "# no matrix", "{}", "{{}}", "{0", "{{0}", "{{0}} junk", "{{0,}}",
        "{{0},{0,1}}", "{{0,1}}", "{{{0}}}", "{{(0]}}", "{{0 +}}", "{{0 0}}", "{{1e999}}",
        "{{coefxi[1,]}}", "{{(0,1)}}", "parse \"outside\"", "{{parse}}", "!touch outside\n{{0}}",
        "{{0}}\n!touch outside", "{{0}}; Get[\"outside\"]", "Sqrt[]", "Sqrt[1,2]", "coefxi[1]",
        "coefxi[1,0,0]", "coefdelta[1,0]", "coefkappa[1,0]", "coefV[1,1]", "unknown[1]", "I[1]"}) {
    auto files = valid;
    files["ham_0.1"] = matrix;
    rejects_template(files, "invalid or unsupported matrix " + matrix);
  }
  auto nul = valid;
  nul["ham_0.1"] = std::string("{{0}}\0suffix", 12);
  rejects_template(nul, "matrix NUL");
  auto deep = valid;
  deep["ham_0.1"] = std::string(4096, '(') + '0' + std::string(4096, ')');
  accepts(deep);
  deep["unreferenced"] = "!not-executed\nparse outside";
  accepts(deep);

  Entries rectangular{{"data.in", "# symtype QS\n#!9\n1 0 2\n# SCALE 1\n"
                                  "0 1\n2\nDIAG ham2\n1 2\n1\nDIAG ham1\n"
                                  "f 0 0\n1\n0 1 1 2\n1\n2\ne\n0\n"
                                  "dA_d\n1\n1 2 0 1\nop.row\n"},
                      {"ham2", "{{eps1, Sqrt[2]}, {Sqrt[2], U1}}"},
                      {"ham1", "0"}, {"op.row", "{0, coefxi[1, 0]}"}};
  accepts(rectangular, 2);
  for (const auto &[from, to] : std::vector<std::pair<std::string, std::string>>{
         {"0 1 1 2", "0 1 1 3"}, {"0 1 1 2", "0 1 1"}, {"0 1 1 2", "0 1 1 2 3"},
         {"1\n2\ne", "1 2\n2\ne"}, {"1\n2\ne", "1\ne"}, {"1\n2\ne", "nan\n2\ne"},
         {"1\n2\ne", "eps1\n2\ne"}, {"op.row\n", ""}, {"op.row\n", "0 1\n2\n"},
         {"dA_d\n1", "dA_d\n2"}}) {
    auto files = rectangular;
    replace(files, from, to);
    rejects_template(files, "invalid operator block " + to);
  }
  for (const std::string matrix : {"{{0}, {1}}", "{{0,1},{2,3}}", "{{0,1}", "{{0}}"}) {
    auto files = rectangular;
    files["op.row"] = matrix;
    rejects_template(files, "invalid external operator matrix " + matrix);
  }
  auto inline_row = rectangular;
  replace(inline_row, "op.row", "0 1");
  inline_row.erase("op.row");
  accepts(inline_row, 2);
}

void real_template(const fs::path &root) {
  // Optional read-only smoke test. Pass the template directory itself, e.g.
  // test/mixchain+nrgrun/templates/QS or the DMFT application's template/.
  // No scratch directory, copies, shell calls, initialization or diagonalization.
  check(fs::is_directory(fs::symlink_status(root)), "Expected a nonsymlink template directory");
  Entries files;
  for (const auto &entry : fs::recursive_directory_iterator(root)) {
    const auto status = entry.symlink_status();
    check(fs::is_directory(status) || fs::is_regular_file(status), "Unsafe real template entry");
    if (!fs::is_regular_file(status)) continue;
    std::ifstream input(entry.path(), std::ios::binary);
    check(bool(input), "Cannot read real template entry");
    files.emplace(entry.path().lexically_relative(root).generic_string(), std::string(std::istreambuf_iterator<char>(input), {}));
    check(!input.bad(), "Error reading real template entry");
  }
  const auto info = NRG::Workflow::inspect_template_bundle(files);
  accepts(files, info.subspaces);
  std::cout << "Read-only template: symmetry=" << info.symmetry << " channels=" << info.channels
            << " subspaces=" << info.subspaces << " files=" << files.size() << '\n';
}

} // namespace

int main(int argc, char **argv) {
  try {
    check(argc <= 2, "Usage: workflow_preparation [read-only-template-root]");
    const auto cwd = fs::current_path();
    chain_fixtures();
    chain_precedence();
    chain_validation();
    chain_scales();
    combined_preparation();
    template_structure();
    template_dependencies();
    template_matrices();
    if (argc == 2) real_template(fs::path(argv[1]));
    check(fs::current_path() == cwd, "Preparation changed CWD");
    std::cout << "Solver preparation tests passed (" << assertions << " checks).\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "Solver preparation tests failed: " << error.what() << '\n';
    return 1;
  }
}
