#ifndef _oprecalc_hpp_
#define _oprecalc_hpp_

#include <string>
#include <set>
#include <memory>
#include "time_mem.hpp"
#include "operators.hpp"
#include "symmetry.hpp"
#include "step.hpp"
#include "params.hpp"
#include "algo.hpp"
#include "stats.hpp"

namespace NRG {

template<typename S>
class Oprecalc {
 private:
   RUNTYPE runtype;
   std::shared_ptr<Symmetry<S>> Sym;
   MemTime &mt; // ref!
   const Params &P;
 public:
   // Operators required to calculate expectation values and spectral densities
   struct Ops : public std::set<std::pair<std::string,std::string>> {
     void report(std::ostream &F = std::cout) {
       F << std::endl << "Computing the following operators:" << std::endl;
       for (const auto &[type, name]: *this) fmt::print("{} {}\n", name, type);
     }
     // Singlet operators are always recomputed in the first NRG run, so that we can calculate the expectation values.
     [[nodiscard]] bool do_s(const std::string &name, const Params &P, const Step &step) {
       if (step.nrg()) return true;                                          // for computing <O>
       if (step.dmnrg() && P.fdmexpv && step.N() <= P.fdmexpvn) return true; // for computing <O> using FDM algorithm
       return this->count({"s", name});
     }
     [[nodiscard]] bool do_g(const std::string &name, const Params &P, const Step &step) {
       if (step.nrg()) return true;                                          // for computing <O>
       if (step.dmnrg() && P.fdmexpv && step.N() <= P.fdmexpvn) return true; // for computing <O> using FDM algorithm
       return this->count({"g", name});
     }
   };
   Ops ops;

   // Spectral densities
   struct SL : public speclist<S> {
     void calc(const Step &step, const DiagInfo<S> &diag, DensMatElements<S> &rho, DensMatElements<S> &rhoFDM,
               const Stats<S> &stats, MemTime &mt, const Params &P) {
       const auto section_timing = mt.time_it("spec");
       for (auto &i : *this) i.calc(step, diag, rho, rhoFDM, stats);
     }
   };
   SL sl;
 
   // Wrapper routine for recalculations
   template <typename RecalcFnc>
     MatrixElements<S> recalc(const std::string &name, const MatrixElements<S> &mold, RecalcFnc recalc_fnc, const std::string &tip,
                              const Step &step, const DiagInfo<S> &diag, const SubspaceStructure &substruct) {
       nrglog('0', "Recalculate " << tip << " " << name);
       auto mnew = recalc_fnc(diag, substruct, mold);
       if (tip == "g") Sym->recalc_global(step, diag, substruct, name, mnew);
       return mnew;
     }

   template <typename ... Args>
     MatrixElements<S> recalc_or_clear(const bool selected, Args&& ... args) {
       return selected ? recalc(std::forward<Args>(args)...) : MatrixElements<S>();
     }

   // Recalculate operator matrix representations
   void recalculate_operators(Operators<S> &a, const Step &step, const DiagInfo<S> &diag, const SubspaceStructure &substruct) {
       const auto section_timing = mt.time_it("recalc");
       for (auto &[name, m] : a.ops)
         m = recalc_or_clear(ops.do_s(name, P, step), name, m, [this](const auto &... pr) { return Sym->recalc_singlet(pr..., 1);  }, "s", step, diag, substruct);
       for (auto &[name, m] : a.opsp)
         m = recalc_or_clear(ops.count({"p", name}),  name, m, [this](const auto &... pr) { return Sym->recalc_singlet(pr..., -1); }, "p", step, diag, substruct);
       for (auto &[name, m] : a.opsg)
         m = recalc_or_clear(ops.do_g(name, P, step), name, m, [this](const auto &... pr) { return Sym->recalc_singlet(pr...,  1); }, "g", step, diag, substruct);
       for (auto &[name, m] : a.opd)
         m = recalc_or_clear(ops.count({"d", name}),  name, m, [this](const auto &... pr) { return Sym->recalc_doublet(pr...);     }, "d", step, diag, substruct);
       for (auto &[name, m] : a.opt)
         m = recalc_or_clear(ops.count({"t", name}),  name, m, [this](const auto &... pr) { return Sym->recalc_triplet(pr...);     }, "t", step, diag, substruct);
       for (auto &[name, m] : a.opot)
         m = recalc_or_clear(ops.count({"ot", name}), name, m, [this](const auto &... pr) { return Sym->recalc_orb_triplet(pr...); }, "ot", step, diag, substruct);
       for (auto &[name, m] : a.opq)
         m = recalc_or_clear(ops.count({"q", name}),  name, m, [this](const auto &... pr) { return Sym->recalc_quadruplet(pr...);  }, "q", step, diag, substruct);
     }

   // Establish the data structures for storing spectral information [and prepare output files].
   template<typename A, typename M>
     [[nodiscard]] bool prepare_spec_algo(std::string prefix, const Params &P, FactorFnc ff, CheckFnc cf, M && op1, M && op2, int spin,
                            std::string name, const gf_type gt) {
       BaseSpectrum<S> spec(std::forward<M>(op1), std::forward<M>(op2), spin, std::make_shared<A>(name, prefix, gt, P), ff, cf);
       sl.push_back(spec);
       return true; // recalculation of operators required
     }

   template<typename ... Args>
     [[nodiscard]] bool prepare_spec(std::string prefix, Args && ... args) {
       bool b = false;
       if (prefix == "gt") {
         if (runtype == RUNTYPE::NRG) b |= prepare_spec_algo<Algo_GT<S,0>>(prefix, P, std::forward<Args>(args)...);
         return b;
       }
       if (prefix == "i1t") {
         if (runtype == RUNTYPE::NRG) b |= prepare_spec_algo<Algo_GT<S,1>>(prefix, P, std::forward<Args>(args)...);
         return b;
       }
       if (prefix == "i2t") {
         if (runtype == RUNTYPE::NRG) b |= prepare_spec_algo<Algo_GT<S,2>>(prefix, P, std::forward<Args>(args)...);
         return b;
       }
       if (prefix == "chit") {
         if (runtype == RUNTYPE::NRG) b |= prepare_spec_algo<Algo_CHIT<S>>(prefix, P, std::forward<Args>(args)...);
         return b;
       }
       // If we did not return from this funciton by this point, what we are computing is the spectral function. There are
       // several possibilities in this case, all of which may be enabled at the same time.
       if (runtype == RUNTYPE::NRG) {
         if (P.finite)     b |= prepare_spec_algo<Algo_FT<S>>    (prefix, P, std::forward<Args>(args)...);
         if (P.finitemats) b |= prepare_spec_algo<Algo_FTmats<S>>(prefix, P, std::forward<Args>(args)...);
       }
       if (runtype == RUNTYPE::DMNRG) {
         if (P.dmnrg)     b |= prepare_spec_algo<Algo_DMNRG<S>>(prefix, P, std::forward<Args>(args)...);
         if (P.dmnrgmats) b |= prepare_spec_algo<Algo_DMNRGmats<S>>(prefix, P, std::forward<Args>(args)...);
         if (P.cfs)       b |= prepare_spec_algo<Algo_CFS<S>>(prefix, P, std::forward<Args>(args)...);
         if (P.cfsgt)     b |= prepare_spec_algo<Algo_CFSgt<S>>(prefix, P, std::forward<Args>(args)...);
         if (P.cfsls)     b |= prepare_spec_algo<Algo_CFSls<S>>(prefix, P, std::forward<Args>(args)...);
         if (P.fdm)       b |= prepare_spec_algo<Algo_FDM<S>>(prefix, P, std::forward<Args>(args)...);
         if (P.fdmgt)     b |= prepare_spec_algo<Algo_FDMgt<S>>(prefix, P, std::forward<Args>(args)...);
         if (P.fdmls)     b |= prepare_spec_algo<Algo_FDMls<S>>(prefix, P, std::forward<Args>(args)...);
         if (P.fdmmats)   b |= prepare_spec_algo<Algo_FDMmats<S>>(prefix, P, std::forward<Args>(args)...);
       }
       return b;
     }

   // Construct the suffix of the filename for spectral density files: 'A_?-A_?'.
   // If SPIN == 1 or SPIN == -1, '-u' or '-d' is appended to the string.
   [[nodiscard]] auto sdname(const std::string &a, const std::string &b, const int spin) {
     return a + "-" + b + (spin == 0 ? "" : (spin == 1 ? "-u" : "-d"));
   }

   void loopover(const CustomOp<S> &set1, const CustomOp<S> &set2,
                 const string_token &stringtoken, FactorFnc ff, CheckFnc cf,
                 const std::string &prefix,
                 const std::string &type1, const std::string &type2, const gf_type gt, const int spin) {
    for (const auto &[name1, op1] : set1) {
      for (const auto &[name2, op2] : set2) {
        if (const auto name = sdname(name1, name2, spin); stringtoken.find(name)) {
          if (prepare_spec(prefix, ff, cf, op1, op2, spin, name, gt)) {
            ops.insert({type1, name1});
            ops.insert({type2, name2});
          }
        }
      }
    }
  }

  // Reset lists of operators which need to be iterated
  Oprecalc(const RUNTYPE &runtype, const Operators<S> &a, std::shared_ptr<Symmetry<S>> Sym, MemTime &mt, const Params &P) : 
    runtype(runtype), Sym(Sym), mt(mt), P(P) 
  {
    std::cout << std::endl << "Computing the following spectra:" << std::endl;
    // Correlators (singlet operators of all kinds)
    string_token sts(P.specs);
    loopover(a.ops,  a.ops,  sts, Sym->CorrelatorFactorFnc(), Sym->TrivialCheckSpinFnc(), "corr", "s", "s", gf_type::bosonic, 0);
    loopover(a.opsp, a.opsp, sts, Sym->CorrelatorFactorFnc(), Sym->TrivialCheckSpinFnc(), "corr", "p", "p", gf_type::bosonic, 0);
    loopover(a.opsg, a.opsg, sts, Sym->CorrelatorFactorFnc(), Sym->TrivialCheckSpinFnc(), "corr", "g", "g", gf_type::bosonic, 0);
    loopover(a.ops,  a.opsg, sts, Sym->CorrelatorFactorFnc(), Sym->TrivialCheckSpinFnc(), "corr", "s", "g", gf_type::bosonic, 0);
    loopover(a.opsg, a.ops,  sts, Sym->CorrelatorFactorFnc(), Sym->TrivialCheckSpinFnc(), "corr", "g", "s", gf_type::bosonic, 0);
    // Global susceptibilities (global singlet operators)
    string_token stchit(P.specchit);
    loopover(a.ops,  a.ops,  stchit, Sym->CorrelatorFactorFnc(), Sym->TrivialCheckSpinFnc(), "chit", "s", "s", gf_type::bosonic, 0);
    loopover(a.ops,  a.opsg, stchit, Sym->CorrelatorFactorFnc(), Sym->TrivialCheckSpinFnc(), "chit", "s", "g", gf_type::bosonic, 0);
    loopover(a.opsg, a.ops,  stchit, Sym->CorrelatorFactorFnc(), Sym->TrivialCheckSpinFnc(), "chit", "g", "s", gf_type::bosonic, 0);
    loopover(a.opsg, a.opsg, stchit, Sym->CorrelatorFactorFnc(), Sym->TrivialCheckSpinFnc(), "chit", "g", "g", gf_type::bosonic, 0);
    // Dynamic spin susceptibilities (triplet operators)
    string_token stt(P.spect);
    loopover(a.opt, a.opt, stt, Sym->SpinSuscFactorFnc(), Sym->TrivialCheckSpinFnc(),  "spin", "t", "t", gf_type::bosonic, 0);
    string_token stot(P.specot);
    loopover(a.opot, a.opot, stot, Sym->OrbSuscFactorFnc(), Sym->TrivialCheckSpinFnc(), "orbspin", "ot", "ot", gf_type::bosonic, 0);
    const auto varmin = Sym->isfield() ? -1 : 0;
    const auto varmax = Sym->isfield() ? +1 : 0;
    // Spectral functions (doublet operators)
    string_token std(P.specd);
    for (int SPIN = varmin; SPIN <= varmax; SPIN += 2)
      loopover(a.opd, a.opd, std, Sym->SpecdensFactorFnc(), Sym->SpecdensCheckSpinFnc(), "spec", "d", "d", gf_type::fermionic, SPIN);
    string_token stgt(P.specgt);
    for (int SPIN = varmin; SPIN <= varmax; SPIN += 2)
      loopover(a.opd, a.opd, stgt, Sym->SpecdensFactorFnc(), Sym->SpecdensCheckSpinFnc(), "gt", "d", "d", gf_type::fermionic, SPIN);
    string_token sti1t(P.speci1t);
    for (int SPIN = varmin; SPIN <= varmax; SPIN += 2)
      loopover(a.opd, a.opd, sti1t, Sym->SpecdensFactorFnc(), Sym->SpecdensCheckSpinFnc(),"i1t", "d", "d", gf_type::fermionic, SPIN);
    string_token sti2t(P.speci2t);
    for (int SPIN = varmin; SPIN <= varmax; SPIN += 2)
      loopover(a.opd, a.opd, sti2t, Sym->SpecdensFactorFnc(), Sym->SpecdensCheckSpinFnc(), "i2t", "d", "d", gf_type::fermionic, SPIN);
    // Spectral functions (quadruplet operators)
    string_token stq(P.specq);
    loopover(a.opq, a.opq, stq, Sym->SpecdensquadFactorFnc(), Sym->TrivialCheckSpinFnc(),  "specq", "q", "q", gf_type::fermionic, 0);
    ops.report();
  }
};

// Recalculate irreducible matrix elements for Wilson chains.
template<typename S>
void recalc_irreducible(const Step &step, const DiagInfo<S> &diag, const SubspaceStructure &substruct, Opch<S> &opch, 
                        const Symmetry<S> *Sym, MemTime &mt, const Params &P) {
  const auto section_timing = mt.time_it("recalc f");
  if (!P.substeps) {
    opch = Sym->recalc_irreduc(step, diag, substruct);
  } else {
    const auto [N, M] = step.NM();
    for (const auto i: range0(size_t(P.channels)))
      if (i == M) {
        opch[i] = Sym->recalc_irreduc_substeps(step, diag, substruct, i);
      } else {
        for (const auto j: range0(size_t(P.perchannel)))
          opch[i][j] = Sym->recalc_doublet(diag, substruct, opch[i][j]);
      }
  }
}

} // namespace

#endif
