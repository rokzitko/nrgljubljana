#ifndef _output_hpp_
#define _output_hpp_

#include <memory>
#include <ostream>
#include <iomanip>
#include <string>
#include <map>
#include <list>
#include <utility>
#include <vector>

#include <range/v3/all.hpp>
#include <boost/range/adaptor/map.hpp>

#include "traits.hpp"
#include "params.hpp"
#include "io.hpp" // formatted_output
#include "misc.hpp" // range1
#include "step.hpp"
#include "stats.hpp"
#include "symmetry.hpp"
#include "h5.hpp"
#include "eigen.hpp"

namespace NRG {

// Formatted output of the computed expectation values
template<scalar S, typename t_expv = expv_traits<S>>
class ExpvOutput {
 private:
   std::ofstream F;                     // output stream
   std::map<std::string, t_expv> &m;    // reference to the name->value mapping (e.g. from class Stats)
   const std::list<std::string> fields; // list of fields to be output (may be a subset of the fields actually present in m)
   const Params &P;
   void field_numbers() {     // Consecutive numbers for the columns
     F << '#' << formatted_output(1, P) << ' ';
     for (const auto ctr : range1(fields.size())) F << formatted_output(1 + ctr, P) << ' ';
     F << std::endl;
   }
   // Label and field names. Label is the first column (typically the temperature).
   void field_names(const std::string &labelname = "T") {
     F << '#' << formatted_output(labelname, P) << ' ';
     std::transform(fields.cbegin(), fields.cend(), std::ostream_iterator<std::string>(F, " "), [this](const auto op) { return formatted_output(op, P); });
     F << std::endl;
   }
 public:
   // Output the current values for the label and for all the fields
   void field_values(const double labelvalue, const bool cout_dump = true) {
     F << ' ' << formatted_output(labelvalue, P) << ' ';
     std::transform(fields.cbegin(), fields.cend(), std::ostream_iterator<std::string>(F, " "), [this](const auto op) { return formatted_output(m[op], P); }); // NOTE: only real part stored
     F << std::endl;
     if (cout_dump)
       for (const auto &op: fields)
         color_print(P.pretty_out, fmt::emphasis::bold | fg(fmt::color::red), "<{}>={}\n", op, formatted_output(m[op], P)); // NOTE: real and imaginary part shown
   }
   ExpvOutput(const std::string &fn, std::map<std::string, t_expv> &m,
              std::list<std::string> fields, const Params &P) : m(m), fields(std::move(fields)), P(P) {
     F.open(fn);
     field_numbers();
     field_names();
   }
};

// scaled = true -> output scaled energies (i.e. do not multiply by the rescale factor). note: dumpEscale is not applied here.
template<typename T, scalar S>
inline auto scaled_energy(const T e, const Step &step, const Stats<S> &stats, const Params &P) {
   return (e * (P.dumpscaled ? 1.0 : step.scale()) + (P.dumpabs ? stats.total_energy : 0.0))/P.dumpEscale;
}

// Store eigenvalue & quantum numbers information (RG flow diagrams)
class Annotated {
 private:
   std::ofstream F;
   const Params &P;
 public:
   explicit Annotated(const Params &P) : P(P) {}
   template<scalar S, typename MF>
   void dump(const Step &step, const DiagInfo<S> &diag, const Stats<S> &stats,
             MF mult, const std::string &filename = "annotated.dat") {
     if (!P.dumpannotated) return;
     if (!F.is_open()) { // open output file
       F.open(filename);
       F << std::setprecision(P.dumpprecision);
     }
     std::vector<std::pair<double, Invar>> seznam;
     for (const auto &[I, eig] : diag)
       for (const auto e : eig.values.all_rel_zero())
         seznam.emplace_back(e, I);
     ranges::sort(seznam);
     size_t len = std::min<size_t>(seznam.size(), P.dumpannotated); // non-const
     // If states are clustered, we dump the full cluster
     while (len < seznam.size()-1 && my_fcmp(seznam[len].first, seznam[len-1].first, P.grouptol) == 0) len++;
     const auto scale = [&step, &stats, this](auto x) { return scaled_energy(x, step, stats, P); };
     if (P.dumpgroups) {
       // Group by degeneracies
       for (size_t i = 0; i < len;) { // i increased in the while loop below
         const auto [e0, I0] = seznam[i];
         F << scale(e0);
         std::vector<std::string> QNstrings;
         size_t total_degeneracy = 0; // Total number of levels (incl multiplicity)
         while (i < len && my_fcmp(seznam[i].first, e0, P.grouptol) == 0) {
           const auto [e, I] = seznam[i];
           QNstrings.push_back(to_string(I));
           total_degeneracy += mult(I);
           i++;
         }
         ranges::sort(QNstrings);
         for (const auto &j : QNstrings) F << " (" << j << ")";
         F << " [" << total_degeneracy << "]" << std::endl;
       }
     } else {
       seznam.resize(len); // truncate!
       for (const auto &[e, I] : seznam)
         F << scale(e) << " " << I << std::endl;
     }
     F << std::endl; // Consecutive iterations are separated by an empty line
   }
};

template<scalar S>
auto singlet_operators_for_expv_evaluation(const Operators<S> &operators)
{
  std::list<std::string> ops;
  for (const auto &name : operators.ops  | boost::adaptors::map_keys) ops.push_back(name);
  for (const auto &name : operators.opsg | boost::adaptors::map_keys) ops.push_back(name);
  return ops;
}

// Handle all output
template<scalar S>
struct Output {
  const RUNTYPE runtype;
  const Params &P;
  Annotated annotated;
  std::ofstream Fenergies;  // all energies (different file for NRG and for DMNRG)
  std::ofstream Fstates;    // all states
  std::ofstream Freport;    // select state energies and diagonal matrix elements of singlet operators
  std::unique_ptr<ExpvOutput<S>> custom;
  std::unique_ptr<ExpvOutput<S>> customfdm;
  std::unique_ptr<H5Easy::File> h5raw;
  Output(const RUNTYPE &runtype, const Operators<S> &operators, Stats<S> &stats, const Params &P,
         const std::string filename_energies= "energies.nrg"s,
         const std::string filename_states = "states.nrg"s,
         const std::string filename_report = "report.nrg"s,
         const std::string filename_custom = "custom",
         const std::string filename_customfdm = "customfdm")
    : runtype(runtype), P(P), annotated(P)
    {
      if (runtype == RUNTYPE::NRG) {
        if (P.dumpenergies) Fenergies.open(filename_energies);
        if (P.dumpstates) Fstates.open(filename_states);
        if (P.reportdiagonal) Freport.open(filename_report);
      }
      const auto ops = singlet_operators_for_expv_evaluation(operators);
      if (runtype == RUNTYPE::NRG)
        custom = std::make_unique<ExpvOutput<S>>(filename_custom, stats.expv, ops, P);
      else if (runtype == RUNTYPE::DMNRG && P.fdmexpv)
        customfdm = std::make_unique<ExpvOutput<S>>(filename_customfdm, stats.fdmexpv, ops, P);
      if (P.h5raw) {
        const auto filename_h5 = runtype == RUNTYPE::NRG ? "raw.h5" : "raw-dm.h5";
        h5raw = std::make_unique<H5Easy::File>(filename_h5, H5Easy::File::Overwrite);
        P.h5save(*h5raw);
      }
    }
  // Dump eigenvalues from the diagonalisation to a file.
  void dump_energies(const int N, const DiagInfo<S> &diag, const double rescaled_by = 1.0) {
    if (!Fenergies) return;
    Fenergies << std::endl << "===== Iteration number: " << N << std::endl;
    dump_all_energies(diag, rescaled_by, Fenergies, P);
  }
  void dump_states(const int N, const DiagInfo<S> &diag, const double rescaled_by = 1.0) {
    if (!Fstates) return;
    Fstates << std::endl << "===== Iteration number: " << N << std::endl;
    dump_all_states(diag, rescaled_by, Fstates, P);
  }
  void reportdiagonal(const Step &step,
                      Stats<S> &stats,
                      const DiagInfo<S> &diag,
                      const Operators<S> &operators,
                      const Params &P) {
    if (!Freport) return;
    if (P.reportdiagonallast && !step.last()) return;
    Freport << "=== Report N=" << step.ndx() << std::endl;
    for (auto &[I, eig] : diag) {
      const size_t nmax = std::min(size_t(P.reportdiagonal), size_t(eig.getnrkept()));
      if (nmax) {
        Freport << "Sector I=" << I << std::endl;
        for (size_t n = 0; n < nmax; n++) {
          Freport << "I=" << I << " n=" << n << " E=" << scaled_energy(eig.values.rel_zero(n), step, stats, P) << " ";
          operators.dump_diagonal_I_n(I, n, Freport);
        }
      }
    }
    if (!step.last())
      Freport << std::endl;
  }};
} // namespace

#endif
