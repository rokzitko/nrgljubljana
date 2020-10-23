#ifndef _output_h_
#define _output_h_

#include <memory>
#include <ostream>
#include <iomanip>
#include <string>
#include <map>
#include <list>
#include <vector>
#include "traits.h"
#include "params.h"
#include "io.h" // formatted_output
#include "misc.h" // range1
#include "step.h"
#include "stats.h"
#include "symmetry.h"
#include <range/v3/all.hpp>
#include <boost/range/adaptor/map.hpp>

// Formatted output of the computed expectation values
template<typename S>
class ExpvOutput {
 private:
   using t_expv = typename traits<S>::t_expv;
   std::ofstream F;                     // output stream
   std::map<std::string, t_expv> &m;    // reference to the name->value mapping
   const std::list<std::string> fields; // list of fields to be output (may be a subset of the fields actually present in m)
   const Params &P;
   void field_numbers() {     // Consecutive numbers for the columns
     F << '#' << formatted_output(1, P) << ' ';
     for (const auto ctr : range1(fields.size())) F << formatted_output(1 + ctr, P) << ' ';
     F << std::endl;
   }
   // Label and field names. Label is the first column (typically the temperature).
   void field_names(const std::string labelname = "T") {
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
         fmt::print(fmt::emphasis::bold | fg(fmt::color::red), "<{}>={}\n", op, to_string(m[op])); // NOTE: real and imaginary part shown
   }
   ExpvOutput(const std::string &fn, std::map<std::string, t_expv> &m_, 
                   const std::list<std::string> &fields_, const Params &P_) : m(m_), fields(fields_), P(P_) {
     F.open(fn);
     field_numbers();
     field_names();
   }
};

// Store eigenvalue & quantum numbers information (RG flow diagrams)
class Annotated {
 private:
   std::ofstream F;
   // scaled = true -> output scaled energies (i.e. do not multiply by the rescale factor)
   template<typename S>
     inline auto scaled_energy(typename traits<S>::t_eigen e, const Step &step, const Stats<S> &stats,
                               bool scaled = true, bool absolute = false) {
     return e * (scaled ? 1.0 : step.scale()) + (absolute ? stats.total_energy : 0.0);
   }
   const Params &P;
 public:
   explicit Annotated(const Params &P) : P(P) {}
   template<typename S> void dump(const Step &step, const DiagInfo<S> &diag, const Stats<S> &stats, 
                                  std::shared_ptr<Symmetry<S>> Sym, const std::string filename = "annotated.dat") {
     if (!P.dumpannotated) return;
     if (!F.is_open()) { // open output file
       F.open(filename);
       F << std::setprecision(P.dumpprecision);
     }
     std::vector<std::pair<double, Invar>> seznam;
     for (const auto &[I, eig] : diag)
       for (const auto e : eig.value_zero)
         seznam.emplace_back(e, I);
     ranges::sort(seznam);
     size_t len = std::min<size_t>(seznam.size(), P.dumpannotated); // non-const
     // If states are clustered, we dump the full cluster
     while (len < seznam.size()-1 && my_fcmp(seznam[len].first, seznam[len-1].first, P.grouptol) == 0) len++;
     auto scale = [&step, &stats, this](auto x) { return scaled_energy(x, step, stats, P.dumpscaled, P.dumpabs); };
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
           total_degeneracy += Sym->mult(I);
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

// Handle all output
template<typename S>
struct Output {
  const RUNTYPE runtype;
  const Params &P;
  Annotated annotated;
  std::ofstream Fenergies;  // all energies (different file for NRG and for DMNRG)
  std::unique_ptr<ExpvOutput<S>> custom;
  std::unique_ptr<ExpvOutput<S>> customfdm;
  Output(const RUNTYPE &runtype, const IterInfo<S> &iterinfo, Stats<S> &stats, const Params &P,
              const std::string filename_energies= "energies.nrg"s,
              const std::string filename_custom = "custom", 
              const std::string filename_customfdm = "customfdm")
    : runtype(runtype), P(P), annotated(P) {
      // We dump all energies to separate files for NRG and DM-NRG runs. This is a very convenient way to check if both
      // runs produce the same results.
      if (P.dumpenergies && runtype == RUNTYPE::NRG) Fenergies.open(filename_energies);
      std::list<std::string> ops;
      for (const auto &name : iterinfo.ops  | boost::adaptors::map_keys) ops.push_back(name);
      for (const auto &name : iterinfo.opsg | boost::adaptors::map_keys) ops.push_back(name);
      if (runtype == RUNTYPE::NRG)
        custom = std::make_unique<ExpvOutput<S>>(filename_custom, stats.expv, ops, P);
      else if (runtype == RUNTYPE::DMNRG && P.fdmexpv) 
        customfdm = std::make_unique<ExpvOutput<S>>(filename_customfdm, stats.fdmexpv, ops, P);
    }
  // Dump all energies in diag to a file
  void dump_all_energies(const DiagInfo<S> &diag, const int N) {
    if (!Fenergies) return;
    Fenergies << std::endl << "===== Iteration number: " << N << std::endl;
    diag.dump_value_zero(Fenergies);
  }
};

#endif
