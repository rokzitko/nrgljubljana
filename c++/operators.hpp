#ifndef _operators_hpp_
#define _operators_hpp_

#include <map>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include "invar.hpp"
#include "traits.hpp"
#include "eigen.hpp"
#include "numerics.hpp" // read_matrix
#include "params.hpp"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/range/adaptor/map.hpp>
#include <range/v3/all.hpp>

namespace NRG {

template<typename S>
class MatrixElements : public std::map<Twoinvar, typename traits<S>::Matrix> {
 public:
   MatrixElements() = default;
   MatrixElements(std::ifstream &fdata, const DiagInfo<S> &diag) {
     size_t nf = 0; // Number of I1 x I2 combinations
     fdata >> nf;
     for (const auto i : range0(nf)) {
       Invar I1, I2;
       fdata >> I1 >> I2;
       if (const auto it1 = diag.find(I1), it2 = diag.find(I2); it1 != diag.end() && it2 != diag.end())
         read_matrix(fdata, (*this)[{I1, I2}], it1->second.getnrstored(), it2->second.getnrstored());
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
       // Current matrix dimensions
       const auto size1 = mat.size1();
       const auto size2 = mat.size2();
       if (size1 == 0 || size2 == 0) continue;
       const auto &[nr1, nr2] = diag.dims(I1, I2); // Target matrix dimensions
       my_assert(nr1 <= size1 && nr2 <= size2);
       if (nr1 == size1 && nr2 == size2)           // Trimming not necessary!!
         continue;
       ublas::matrix_range<typename traits<S>::Matrix> m2(mat, ublas::range(0, nr1), ublas::range(0, nr2));
       typename traits<S>::Matrix m2new = m2;
       mat.swap(m2new);
     }
   }
   std::ostream &insertor(std::ostream &os) const {
     for (const auto &[II, mat] : *this)
       os << "----" << II << "----" << std::endl << mat << std::endl;
     return os;
   }
   friend std::ostream &operator<<(std::ostream &os, const MatrixElements<S> &m) { return m.insertor(os); }
   friend void dump_diagonal_op(const std::string &name, const MatrixElements<S> &m, const size_t max_nr, std::ostream &F) {
     F << "Diagonal matrix elements of operator " << name << std::endl;
     for (const auto &[II, mat] : m) {
       const auto & [I1, I2] = II;
       if (I1 == I2) {
         F << I1 << ": ";
         dump_diagonal_matrix(mat, max_nr, F);
       }
     }
   }
};

template<typename S>
class DensMatElements : public std::map<Invar, typename traits<S>::Matrix> {
 public:
   template <typename MF>
     auto trace(MF mult) const {
       return ranges::accumulate(*this, 0.0, [mult](double acc, const auto z) { const auto &[I, mat] = z; 
         return acc + mult(I) * trace_real(mat); });
     }
   void save(const size_t N, const Params &P, const std::string &prefix) const {
     const auto fn = P.workdir.rhofn(N, prefix);
     std::ofstream MATRIXF(fn, std::ios::binary | std::ios::out);
     if (!MATRIXF) throw std::runtime_error(fmt::format("Can't open file {} for writing.", fn));
     boost::archive::binary_oarchive oa(MATRIXF);
     oa << this->size();
     for (const auto &[I, mat] : *this) {
       oa << I;
       NRG::save(oa, mat);
       if (MATRIXF.bad()) throw std::runtime_error(fmt::format("Error writing {}", fn));  // Check each time
     }
     MATRIXF.close();
   }
   void load(const size_t N, const Params &P, const std::string &prefix, const bool remove_files) {
     const auto fn = P.workdir.rhofn(N, prefix);
     std::ifstream MATRIXF(fn, std::ios::binary | std::ios::in);
     if (!MATRIXF) throw std::runtime_error(fmt::format("Can't open file {} for reading", fn));
     boost::archive::binary_iarchive ia(MATRIXF);
     size_t nr = 0;
     ia >> nr;
     for (const auto cnt : range0(nr)) {
       Invar inv;
       ia >> inv;
       NRG::load(ia, (*this)[inv]);
       if (MATRIXF.bad()) throw std::runtime_error(fmt::format("Error reading {}", fn));  // Check each time
     }
     MATRIXF.close();
     if (remove_files)
	if (NRG::remove(fn)) throw std::runtime_error(fmt::format("Error removing {}", fn));
   }
};

// Map of operator matrices
template<typename S>
struct CustomOp : public std::map<std::string, MatrixElements<S>> {
   void trim(const DiagInfo<S> &diag) {
     for (auto &op : *this | boost::adaptors::map_values) op.trim(diag);
   }
};

// Vector containing irreducible matrix elements of f operators.
template<typename S>
using OpchChannel = std::vector<MatrixElements<S>>;

// Each channel contains P.perchannel OpchChannel matrices.
template<typename S>
class Opch : public std::vector<OpchChannel<S>> {
 public:
   Opch() = default;
   explicit Opch(const size_t nrch) { this->resize(nrch); }
   Opch(std::ifstream &fdata, const DiagInfo<S> &diag, const Params &P) {
     this->resize(P.channels);
     for (const auto i : range0(size_t(P.channels))) {
       (*this)[i] = OpchChannel<S>(P.perchannel);
       for (const auto j : range0(size_t(P.perchannel))) {
         char ch = 0;
         size_t iread = 0, jread = 0;
         fdata >> ch >> iread >> jread;
         my_assert(ch == 'f' && i == iread && j == jread);
         (*this)[i][j] = MatrixElements<S>(fdata, diag);
       }
     }
   }
   void dump() {
     std::cout << std::endl;
     for (const auto &&[i, ch] : *this | ranges::views::enumerate)
       for (const auto &&[j, mat] : ch | ranges::views::enumerate)
         std::cout << fmt::format("<f> dump, i={} j={}\n", i, j) << mat << std::endl;
     std::cout << std::endl;
   }
};

// Object of class IterInfo cotains full information about matrix representations when entering stage N of the NRG
// iteration.
template<typename S>
class IterInfo {
 public:
   Opch<S> opch;     // f operators (channels)
   CustomOp<S> ops;  // singlet operators (even parity)
   CustomOp<S> opsp; // singlet operators (odd parity)
   CustomOp<S> opsg; // singlet operators [global op]
   CustomOp<S> opd;  // doublet operators (spectral functions)
   CustomOp<S> opt;  // triplet operators (dynamical spin susceptibility)
   CustomOp<S> opq;  // quadruplet operators (spectral functions for J=3/2)
   CustomOp<S> opot; // orbital triplet operators

   void dump_diagonal(const size_t max_nr, std::ostream &F = std::cout) const {
     if (max_nr) {
       for (const auto &[name, m] : ops)  dump_diagonal_op(name, m, max_nr, F);
       for (const auto &[name, m] : opsg) dump_diagonal_op(name, m, max_nr, F);
     }
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
};

} // namespace

#endif
