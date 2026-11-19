// spectrum_names.hpp - Names of the spectrum files

#ifndef _spectrum_names_hpp_
#define _spectrum_names_hpp_

#include <string>

namespace NRG {

// The suffix of a doublet spectrum under a symmetry type with a field: -u or -d for spin +1 or -1, none for 0
inline std::string spectrum_spin_suffix(const int spin) { return spin == 0 ? "" : (spin == 1 ? "-u" : "-d"); }

// The correlator of two operators, with the spin suffix
inline std::string spectrum_correlator(const std::string &op1, const std::string &op2, const int spin = 0) {
  return op1 + "-" + op2 + spectrum_spin_suffix(spin);
}

// The file of a spectrum without its extension (.dat for the broadened spectrum, .bin for the bins): prefix of the
// spectrum list (spec, corr, ...), algorithm (FT, DMNRG, CFS, FDM, ...), correlator
inline std::string spectrum_file_stem(const std::string &prefix, const std::string &algorithm, const std::string &correlator) {
  return prefix + "_" + algorithm + "_dens_" + correlator;
}

} // namespace NRG

#endif
