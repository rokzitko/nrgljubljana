#ifndef NRGLJUBLJANA_SOLVER_PREPARATION_HPP
#define NRGLJUBLJANA_SOLVER_PREPARATION_HPP

#include <map>
#include <string>

namespace NRG::Workflow {

// In-memory preparation only: these functions never modify their input, perform I/O,
// change CWD or solver state, evaluate model expressions, or diagonalize.
// Invalid/unsupported input throws std::runtime_error. No diagnostics are printed.
// Callers own parameter-section selection and filesystem/symlink policy.

// param_section contains only [param] entries, not a flattened parameter file.
// Numeric values are scalar decimal tokens (outer whitespace is ignored); expressions,
// hexadecimal numbers, nonfinite numbers and nonzero values rounded to zero are
// rejected. Representable subnormals are accepted. Unknown keys are ignored.
//
// Nmax and explicit Tmin are mutually exclusive, even when inactive/invalid.
// Explicit Tmin supersedes T and Tmin_ratio WITHOUT parsing either ignored value.
// Otherwise both T and Tmin_ratio, if present, must be finite scalars; a positive
// pair supplies a finite positive product which overrides Nmax. An unpaired value
// is ignored. Nmax syntax is checked even when overridden: -?[0-9]+, without '+'.
// Only the effective Nmax is range-checked, and must lie in [1,998]. Nonpositive
// explicit Tmin is an error, not permission to fall back to the ratio cutoff.
//
// Positive cutoffs use Lambda=2, z=1, bandrescale=1, discretization=Z by default.
// Require Lambda>1, 0<z<=1, bandrescale>0. The scheme's first character is C/Y/Z.
// Count scales >= cutoff, preserving the scalar nrginit/instantiate operation
// order (not logarithmic inversion): sc0 * pow(Lambda,-n/2) * bandrescale, n>=0.
// sc0 = prefactor * pow(Lambda,1-z); prefactor is (1+1/Lambda)/2 for Y and
// (1-1/Lambda)/log(Lambda) for C/Z. A nonfinite scale or result outside [1,998] fails.
// Scale parameters are not evaluated when an explicit Nmax remains effective.
unsigned resolve_chain_length(const std::map<std::string, std::string> &param_section);

// Static parameter checks for the combined discretize(integral=true) + full
// instantiate workflow. Supply parsed/trimmed [param] entries for EACH rendered
// twist; z must be numeric, not a placeholder (missing z uses the native default 1).
// A supplied Nmax is already effective and must be in [1,998]: Tmin/T/Tmin_ratio
// are then ignored, including their mutual exclusion with Nmax. If Nmax is absent,
// resolve_chain_length supplies it in a local copy, without changing the caller.
// Other numeric/bool parsing and defaults are the native adapt/initializer ones.
//
// Checked keys: Lambda, z, bandrescale, discretization (C/Y/Z), xmax (>1 for discretization), xfine,
// outputstep, dx_fine, dx_fast, allowed_error (including integral CQUAD minimum),
// max_subdiv (including the dx_fine sub-step bound), max_abs, secant_eps,
// secant_factor, secant_max_iter, adapt, hardgap, boundary, density_interpolation,
// Nmax, mMAX, preccpp, tridiag_method, rescalexi, nrgchain_tables_save,
// nrgchain_tables_load, nrgchain_tridiag (legacy nrgchains_tridiag fallback), tri.
// band selects the native flat-band/adapt compatibility check; dos is not opened.
// loadg's boolean syntax is checked only when adapt=true. f_method is ignored
// because integral=true overrides it; CLI-only integration options are not keys.
// Reject rescalexi=true, disabled chain tridiagonalization, table loading and
// tri=cpp. Keep the native hardgap/flat-band restrictions, but impose no general
// adapt=false, fdm, symmetry or other DMFT-specific policy; these remain caller-owned.
// Success is not validation of density/table data, expressions or numerical
// convergence. Throws std::runtime_error without logging or changing GSL handlers.
void validate_preparation_parameters(const std::map<std::string, std::string> &rendered_param_section);

struct TemplateInfo {
  std::string symmetry;
  unsigned channels;
  unsigned subspaces;
};

// files maps canonical '/'-separated paths RELATIVE TO THE TEMPLATE ROOT to bytes:
// "data.in", "ham_0.1", "operators/op.matrix", never an implicit "template/" prefix.
// Keys must be nonempty relative file paths without '.', '..', repeated/trailing
// separators or NUL. References may contain './' or repeated separators, but not
// '..', absolute paths or file-valued parents. There is no CWD/root fallback.
//
// Supported profile: real version-9, one-channel QS. Inspect the native header,
// positive seed dimensions/scale, unique QS subspaces, DIAG dependencies, one
// leading f block, an e placeholder, and s/p/d/t/o/g/q operator blocks. Check
// referenced matrix shapes and scalar syntax, not numerical values of expressions.
// Scalar/vector/matrix forms, + - * /, grouping, identifiers, and the initializer's
// Sqrt/sqrt/exp/log/Conjugate/gammaPolCh/coefxi/coefzeta calls are supported.
// Unknown functions, complex profiles, parse directives and shell escapes fail.
// Matrix numeric lexing follows the initializer's strtod grammar and assumes the
// C numeric locale; this API does not change the process locale.
// External operators follow the initializer's first-character 'o' filename rule.
// At z/Z/T, stop: the initializer discards and regenerates that entire suffix.
// Unreferenced file contents are ignored. Success is structural metadata, not a
// guarantee of evaluability, Hermiticity, coefficient availability or solver success.
TemplateInfo inspect_template_bundle(const std::map<std::string, std::string> &files);

} // namespace NRG::Workflow

#endif
