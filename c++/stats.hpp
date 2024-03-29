#ifndef _stats_hpp_
#define _stats_hpp_

#include <string>
#include <map>
#include <vector>
#include "constants.hpp"
#include "traits.hpp"
#include "outfield.hpp"
#include "mp.hpp"
#include "step.hpp"

namespace NRG {

// Structure for storing various statistical quantities calculated during the iteration
template<scalar S, typename t_eigen = eigen_traits<S>, typename t_expv  = expv_traits<S>>
class Stats {
 private:
   const Params &P;
 public:
   t_eigen Egs{};

   // ** Thermodynamic quantities
   double Z{};
   double Zft{};   // grand-canonical partition function (at shell n)
   double Zgt{};   // grand-canonical partition function for computing G(T)
   double Zchit{}; // grand-canonical partition function for computing chi(T)
   TD td;

   //  ** Expectation values
   std::map<std::string, t_expv> expv;    // expectation values of custom operators
   std::map<std::string, t_expv> fdmexpv; // Expectation values computed using the FDM algorithm

   // ** Energies
   // "total_energy" is the total energy of the ground state at the current iteration. This is the sum of all the
   // zero state energies (eigenvalue shifts converted to absolute energies) for all the iteration steps so far.
   t_eigen total_energy{};
   // GS_energy is the energy of the ground states in absolute units. It is equal to the value of the variable
   // "total_energy" at the end of the iteration.
   t_eigen GS_energy{};
   std::vector<double> rel_Egs;        // Values of 'Egs' for all NRG steps.
   std::vector<double> abs_Egs;        // Values of 'Egs' (multiplied by the scale, i.e. in absolute scale) for all NRG steps.
   std::vector<double> energy_offsets; // Values of "total_energy" for all NRG steps.

   // ** Containers related to the FDM-NRG approach
   // Consult A. Weichselbaum, J. von Delft, PRL 99, 076402 (2007).
   vmpf ZnDG;                    // Z_n^D=\sum_s^D exp(-beta E^n_s), sum over **discarded** states at shell n
   vmpf ZnDN;                    // Z'_n^D=Z_n^D exp(beta E^n_0)=\sum_s^D exp[-beta(E^n_s-E^n_0)]
   std::vector<double> ZnDNd;    //
   std::vector<double> wn;       // Weights w_n. They sum to 1.
   std::vector<double> wnfactor; // wn/ZnDG
   double ZZG{};                 // grand-canonical partition function with energies referred to the ground state energy
   double Z_fdm{};               // grand-canonical partition function (full-shell) at temperature T
   double F_fdm{};               // free-energy at temperature T
   double E_fdm{};               // energy at temperature T
   double C_fdm{};               // heat capacity at temperature T
   double S_fdm{};               // entropy at temperature T
   TD_FDM td_fdm;

  explicit Stats(const Params &_P, const std::vector<std::string> &td_fields, const double GS_energy_0,
                 const std::string &filename_td = "td"s, const std::string &filename_tdfdm = "tdfdm"s) :
     P(_P), td(P, filename_td), total_energy(GS_energy_0), rel_Egs(MAX_NDX), abs_Egs(MAX_NDX), energy_offsets(MAX_NDX),
     ZnDG(MAX_NDX), ZnDN(MAX_NDX), ZnDNd(MAX_NDX), wn(MAX_NDX), wnfactor(MAX_NDX), td_fdm(P, filename_tdfdm) {
       td.allfields.add(td_fields, 1);
     }

   void update(const Step &step) {
     total_energy += Egs * step.scale(); // stats.Egs has already been initialized
     std::cout << "Total energy=" << HIGHPREC(total_energy) << "  Egs=" << HIGHPREC(Egs) << std::endl;
     rel_Egs[step.ndx()] = Egs;
     abs_Egs[step.ndx()] = Egs * step.scale();
     energy_offsets[step.ndx()] = total_energy;
   }

   // Called after first run
   void h5save_nrg(H5Easy::File &fd) const {
     h5_dump_scalar(fd, "stats/GS_energy", GS_energy);
     for (Step step(P); !step.end(); step.next()) {
       const auto ndx = step.ndx(); // 0-based!
       const auto prefix = "stats/" + std::to_string(ndx+1); // 1-based!
       h5_dump_scalar(fd, prefix + "/energyscale", step.energyscale());
       h5_dump_scalar(fd, prefix + "/scale", step.scale());
       h5_dump_scalar(fd, prefix + "/Teff", step.Teff());
       h5_dump_scalar(fd, prefix + "/rel_Egs", rel_Egs[ndx]);
       h5_dump_scalar(fd, prefix + "/abs_Egs", rel_Egs[ndx]);
       h5_dump_scalar(fd, prefix + "/energy_offset", energy_offsets[ndx]);
     }
     for (const auto &[name, value]: expv)
       h5_dump_scalar(fd, "expv/" + name, value);
     h5_dump_scalar(fd, "total_energy", total_energy);
   }

   // Called after second run
   void h5save_dmnrg(H5Easy::File &fd) const {
     if (P.fdm) {
       h5_dump_scalar(fd, "fdm/Z", Z_fdm);
       h5_dump_scalar(fd, "fdm/F", F_fdm);
       h5_dump_scalar(fd, "fdm/E", E_fdm);
       h5_dump_scalar(fd, "fdm/C", C_fdm);
       h5_dump_scalar(fd, "fdm/S", S_fdm);
       if (P.fdmexpv)
         for (const auto &[name, value]: fdmexpv)
           h5_dump_scalar(fd, "fdm/expv/" + name, value);
     }
   }
};

} // namespace

#endif
