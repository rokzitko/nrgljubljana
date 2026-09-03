#ifndef _operators_hpp_
#define _operators_hpp_

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <ios>
#include <istream>
#include <map>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <ostream>
#include <string>
#include <string_view>
#include <stdexcept>
#include <system_error>
#include <vector>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/range/adaptor/map.hpp>
#include <range/v3/all.hpp>

#include "invar.hpp"
#include "misc.hpp"
#include "traits.hpp"
#include "eigen.hpp"
#include "params.hpp"
#include "h5.hpp"
#include "numerics.hpp"// read_matrix

namespace NRG {

inline constexpr std::string_view density_checkpoint_magic = "NRGLjubljana-density\n";
inline constexpr std::uint32_t density_checkpoint_version = 1;

inline void validate_density_checkpoint_magic(std::istream &input, const std::string &filename) {
  std::array<char, density_checkpoint_magic.size()> magic{};
  input.read(magic.data(), static_cast<std::streamsize>(magic.size()));
  if (!input || !std::equal(magic.begin(), magic.end(), density_checkpoint_magic.begin()))
    throw std::runtime_error(fmt::format("Unsupported density checkpoint format in {}", filename));
}

template<scalar S>
void validate_density_checkpoint_metadata(boost::archive::binary_iarchive &ia, const size_t N,
                                          const std::string &prefix, const Params &P, const std::string &filename) {
  const auto version = read_one<std::uint32_t>(ia);
  if (version != density_checkpoint_version)
    throw std::runtime_error(fmt::format("Unsupported density checkpoint format in {}", filename));
  const auto stored_N = read_one<size_t>(ia);
  const auto stored_prefix = read_one<std::string>(ia);
  const auto stored_Nmax = read_one<size_t>(ia);
  const auto stored_Ninit = read_one<size_t>(ia);
  const auto stored_symtype = read_one<std::string>(ia);
  const auto stored_complex = read_one<bool>(ia);
  const auto stored_param_fingerprint = read_one<std::uint64_t>(ia);
  const auto stored_data_fingerprint = read_one<std::uint64_t>(ia);
  const auto stored_unitary_fingerprint = read_one<std::uint64_t>(ia);
  if (stored_N != N || stored_prefix != prefix || stored_Nmax != P.Nmax || stored_Ninit != P.Ninit.value() ||
      stored_symtype != std::string(P.symtype) || stored_complex != is_complex<S>::value ||
      stored_param_fingerprint != P.checkpoint_param_fingerprint || stored_data_fingerprint != P.checkpoint_data_fingerprint)
    throw std::runtime_error(fmt::format("Density checkpoint {} is incompatible with the current calculation", filename));
  if (stored_unitary_fingerprint != P.checkpoint_unitary_fingerprint)
    throw std::runtime_error(fmt::format("Density checkpoint {} belongs to a different unitary chain", filename));
}

template<scalar S, typename t_matel = matel_traits<S>, typename Matrix = Matrix_traits<S>>
class MatrixElements : public std::map<Twoinvar, Matrix> {
 public:
   MatrixElements() = default;
   MatrixElements(std::istream &fdata, const DiagInfo<S> &diag) {
     const auto nf = read_one<size_t>(fdata); // Number of I1 x I2 combinations
     for ([[maybe_unused]] const auto i : range0(nf)) {
       const auto I1 = read_one<Invar>(fdata);
       const auto I2 = read_one<Invar>(fdata);
       if (const auto it1 = diag.find(I1), it2 = diag.find(I2); it1 != diag.end() && it2 != diag.end())
         (*this)[{I1, I2}] = read_matrix<t_matel>(fdata, it1->second.getnrstored(), it2->second.getnrstored());
       else
         throw std::runtime_error("Corrupted input file.");
     }
     my_assert(this->size() == nf);
   }
   // We trim the matrices containing the irreducible matrix elements of the operators to the sizes that are actually
   // required in the next iterations. This saves memory and leads to better cache usage in recalc_general()
   // recalculations. Note: this is only needed for strategy=all; copying is avoided for strategy=kept.
   void trim(const DiagInfo<S> &diag) {
     for (auto &[II, mat] : *this) {
       const auto &[I1, I2] = II;
       const auto &[dim1, dim2] = diag.dims(I1, I2); // Target matrix dimensions
       mat = trim_matrix(mat, dim1, dim2);
     }
   }
   std::ostream &insertor(std::ostream &os) const {
     for (const auto &[II, mat] : *this)
       os << "----" << II << "----" << std::endl << mat << std::endl;
     return os;
   }
   friend std::ostream &operator<<(std::ostream &os, const MatrixElements<S> &m) { return m.insertor(os); }
   friend void dump_diagonal_op(const std::string &name, const MatrixElements<S> &m, const size_t max_nr, std::ostream &F = std::cout) {
     F << "Diagonal matrix elements of operator " << name << std::endl;
     for (const auto &[II, mat] : m) {
       const auto & [I1, I2] = II;
       if (I1 == I2) {
         F << I1 << ": ";
         dump_diagonal_matrix(mat, max_nr, F);
       }
     }
   }
   friend void dump_diagonal_op_I_n(const std::string &name, const MatrixElements<S> &m, const Invar I, const size_t n, std::ostream &F, const Params &P) {
     for (const auto &[II, mat] : m) {
       const auto & [I1, I2] = II;
       if (I == I1 && I == I2)
         F << "<" << name << ">=" << prec_output(mat(n, n), P) << " ";
     }
   }
   void h5save(H5Easy::File &fd, const std::string &name) const {
     for (const auto &[II, mat] : *this) {
       const auto &[I1, I2] = II;
       h5_dump_matrix(fd, name + "/" + I1.name() + "/" + I2.name() + "/matrix", mat);
     }
   }
};

template<scalar S, typename Matrix = Matrix_traits<S>>
class DensMatElements : public std::map<Invar, Matrix> {
 public:
   template <typename MF> auto trace(MF mult) const {
     return ranges::accumulate(*this, 0.0, {},
                               [mult](const auto &z) { const auto &[I, mat] = z; return mult(I) * trace_real(mat); });
   }
   void save(const size_t N, const Params &P, const std::string &prefix) const {
     const auto fn = P.workdir->rhofn(N, prefix);
     const auto tmp_fn = fn + ".tmp";
     struct temp_file_guard {
       std::string filename;
       bool active = true;
       ~temp_file_guard() {
         if (active) {
           std::error_code ec;
           std::filesystem::remove(filename, ec);
         }
       }
     } guard{tmp_fn};
     std::ofstream MATRIXF(tmp_fn, std::ios::binary | std::ios::out);
     if (!MATRIXF) throw std::runtime_error(fmt::format("Can't open file {} for writing.", fn));
     MATRIXF.write(density_checkpoint_magic.data(), static_cast<std::streamsize>(density_checkpoint_magic.size()));
     if (!MATRIXF) throw std::runtime_error(fmt::format("Error writing {}", fn));
     {
       boost::archive::binary_oarchive oa(MATRIXF);
       oa << density_checkpoint_version
          << N
          << prefix
          << P.Nmax
          << P.Ninit.value()
          << std::string(P.symtype)
          << is_complex<S>::value
          << P.checkpoint_param_fingerprint
          << P.checkpoint_data_fingerprint
          << P.checkpoint_unitary_fingerprint;
       oa << this->size();
       for (const auto &[I, mat] : *this) {
         oa << I;
         NRG::save(oa, mat);
         if (MATRIXF.bad()) throw std::runtime_error(fmt::format("Error writing {}", fn));
       }
     }
     MATRIXF.flush();
     if (!MATRIXF) throw std::runtime_error(fmt::format("Error writing {}", fn));
     MATRIXF.close();
     if (!MATRIXF) throw std::runtime_error(fmt::format("Error writing {}", fn));
     std::filesystem::rename(tmp_fn, fn);
     guard.active = false;
   }
   void load(const size_t N, const Params &P, const std::string &prefix, const bool remove_files) {
     const auto fn = P.workdir->rhofn(N, prefix);
     std::ifstream MATRIXF(fn, std::ios::binary | std::ios::in);
     if (!MATRIXF) throw std::runtime_error(fmt::format("Can't open file {} for reading", fn));
     DensMatElements<S> loaded;
     validate_density_checkpoint_magic(MATRIXF, fn);
     {
       boost::archive::binary_iarchive ia(MATRIXF);
       validate_density_checkpoint_metadata<S>(ia, N, prefix, P, fn);
       const auto nr = read_one<size_t>(ia);
       for ([[maybe_unused]] const auto cnt : range0(nr)) {
         const auto inv = read_one<Invar>(ia);
         loaded[inv] = NRG::load<S>(ia);
         if (MATRIXF.bad()) throw std::runtime_error(fmt::format("Error reading {}", fn));
       }
       if (loaded.size() != nr) throw std::runtime_error(fmt::format("Duplicate subspaces in density checkpoint {}", fn));
     }
     MATRIXF.close();
     if (!MATRIXF) throw std::runtime_error(fmt::format("Error reading {}", fn));
     if (remove_files) {
       if (NRG::remove(fn)) throw std::runtime_error(fmt::format("Error removing {}", fn));
     }
     this->swap(loaded);
   }
   static bool compatible_file(const size_t N, const Params &P, const std::string &prefix) {
     try {
       DensMatElements<S> loaded;
       loaded.load(N, P, prefix, false);
       return true;
     } catch (const std::exception &) {
       return false;
     }
   }
};

// Map of operator matrices
template <scalar S>
struct CustomOp : public std::map<std::string, MatrixElements<S>> {
  void trim(const DiagInfo<S> &diag) {
    for (auto &op : *this | boost::adaptors::map_values) op.trim(diag);
  }
  void h5save(H5Easy::File &fd, const std::string &name) const {
    for (const auto &[n, op] : *this) op.h5save(fd, name + "/" + n);
  }
};

// Vector containing irreducible matrix elements of f operators.
template<scalar S>
using OpchChannel = std::vector<MatrixElements<S>>;

// Each channel contains P.perchannel OpchChannel matrices.
template<scalar S>
class Opch : public std::vector<OpchChannel<S>> {
 public:
   Opch() = default;
   explicit Opch(const Params &P) {
      this->resize(P.channels);
      for (auto &oc: *this) {
        oc.resize(P.perchannel);
        for (auto &o: oc)
          o.clear(); // set all matrix elements to zero
      }
   }
   explicit Opch(std::istream &fdata, const DiagInfo<S> &diag, const Params &P) {
     skip_comments(fdata);
     this->resize(P.channels);
     for (auto &oc : *this) {
       oc.resize(P.perchannel);
       for (auto &o : oc) {
         [[maybe_unused]] const auto ch = read_one<char>(fdata);
         [[maybe_unused]] const auto iread = read_one<size_t>(fdata);
         [[maybe_unused]] const auto jread = read_one<size_t>(fdata);
         my_assert(ch == 'f');
         o = MatrixElements<S>(fdata, diag);
       }
     }
   }
   void dump(std::ostream &F = std::cout) {
     F << std::endl;
     for (const auto &&[i, ch] : *this | ranges::views::enumerate)
       for (const auto &&[j, mat] : ch | ranges::views::enumerate)
         F << fmt::format("<f> dump, i={} j={}\n", i, j) << mat << std::endl;
     F << std::endl;
   }
};

// Object of class Operators cotains full information about matrix representations when entering stage N of the NRG
// iteration.
template<scalar S>
class Operators {
 public:
   Opch<S> opch;     // f operators (channels)
   CustomOp<S> ops;  // singlet operators (even parity)
   CustomOp<S> opsp; // singlet operators (odd parity)
   CustomOp<S> opsg; // singlet operators [global op]
   CustomOp<S> opd;  // doublet operators (spectral functions)
   CustomOp<S> opt;  // triplet operators (dynamical spin susceptibility)
   CustomOp<S> opq;  // quadruplet operators (spectral functions for J=3/2)
   CustomOp<S> opot; // orbital triplet operators

   // Dump diagonal matrix element values for max_nr states in each invariant subspace
   void dump_diagonal(const size_t max_nr, std::ostream &F = std::cout) const {
     if (max_nr) {
       for (const auto &[name, m] : ops)  dump_diagonal_op(name, m, max_nr, F);
       for (const auto &[name, m] : opsg) dump_diagonal_op(name, m, max_nr, F);
     }
   }
   // Dump the diagonal matrix element values for state number n in invariant subspace I
   void dump_diagonal_I_n(const Invar I, const size_t n, std::ostream &F, const Params &P) const {
     for (const auto &[name, m] : ops)  dump_diagonal_op_I_n(name, m, I, n, F, P);
     for (const auto &[name, m] : opsg) dump_diagonal_op_I_n(name, m, I, n, F, P);
     F << std::endl;
   }
   void trim_matrices(const DiagInfo<S> &diag) {
     ops.trim(diag);
     opsp.trim(diag);
     opsg.trim(diag);
     opd.trim(diag);
     opt.trim(diag);
     opq.trim(diag);
     opot.trim(diag);
   }
   void h5save(H5Easy::File &fd, const std::string &name) const {
     ops.h5save(fd, name + "/s");
     opsp.h5save(fd, name + "/sp");
     opsg.h5save(fd, name + "/sg");
     opd.h5save(fd, name + "/d");
     opt.h5save(fd, name + "/t");
     opq.h5save(fd, name + "/q");
     opot.h5save(fd, name + "/ot");
   }
};

} // namespace

#endif
