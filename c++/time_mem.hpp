#ifndef _time_mem_hpp_
#define _time_mem_hpp_

#include <algorithm>
#include <utility>
#include <chrono>
#include <map>
#include <string>
#include "io.hpp"
#include "portabil.hpp"
#include "basicio.hpp"

#define FMT_HEADER_ONLY
#include <fmt/format.h>

namespace NRG {

// Warning: not thread safe!
class Timing {
 private:
   using tp = std::chrono::steady_clock::time_point;
   using dp = std::chrono::duration<double>;
   const tp start_time;         // time when Timing object constructed
   tp timer;                    // start time for timing sections
   bool running{};                // currently timing a section
   std::map<std::string, dp> t; // accumulators
 public:
   Timing() : start_time(std::chrono::steady_clock::now()) {}
   [[nodiscard]] tp now() const  {
     return std::chrono::steady_clock::now();
   }
   void start() {
     my_assert(!running);
     running = true;
     timer   = now();
   }
   dp stop() {
     my_assert(running);
     running = false;
     tp end  = now();
     return end - timer;
   }
   void add(const std::string &timer) {
     t[timer] += stop();
   }
   [[nodiscard]] dp total() const {
     const tp end_time = now();
     return end_time - start_time;
   }
   [[nodiscard]] auto total_in_seconds() const {
     return total().count();
   }
   void report(const double non_negligible_ratio = 0.01) const {
     const auto T_WIDTH  = 12;
     const auto t_all = total();
     std::cout << std::endl << "Timing report" << std::endl;
     std::cout << std::setw(T_WIDTH) << "All"
       << ": " << prec3(t_all.count()) << " s" << std::endl;
     dp t_sum{};
     for (const auto &[name, val] : t) {
       if (val/t_all > non_negligible_ratio) {
         std::cout << std::setw(T_WIDTH) << name << ": " << prec3(val.count()) << " s" << std::endl;
         t_sum += val;
       }
     }
     std::cout << std::setw(T_WIDTH) << "Other" << ": " << prec3((t_all-t_sum).count()) << " s" << std::endl;
   }
};

// Higher-level timing code: time a section for as long as the object is in scope.
class TimeScope {
 private:
   Timing &timer;
   const std::string timer_name;
 public:
   TimeScope(Timing &timer, const std::string &timer_name) : timer(timer), timer_name(timer_name) { timer.start(); }
   TimeScope(const TimeScope &) = delete;
   TimeScope(TimeScope &&) = delete;
   TimeScope & operator=(const TimeScope &) = delete;
   TimeScope & operator=(TimeScope &&) = delete;
   ~TimeScope() { timer.add(timer_name); }
};

class MemoryStats {
 private:
   mutable long peakusage{};
 public:
   auto used() const {
     const auto memused = memoryused();
     peakusage          = std::max<long>(peakusage, memused);
     return memused;
   }
   void report() const {
#ifdef HAS_MEMORY_USAGE
     fmt::print("\nPeak usage: {} MB\n", peakusage / 1024); // NOLINT
#endif
   }
};

class MemTime {
 private:
   MemoryStats ms;
   Timing tm;
 public:
   void brief_report() const {
#ifdef HAS_MEMORY_USAGE
     std::cout << "Memory used: " << long(ms.used() / 1024) << " MB "; // NOLINT
#endif
     std::cout << "Time elapsed: " << prec3(tm.total_in_seconds()) << " s" << std::endl;
   }
   void report() const {
     ms.report();
     tm.report();
   }
   auto time_it(const std::string &name) {
     return TimeScope(tm, name); // measures time when this object exists in a given scope (assign to variable!)
   }
};

} // namespace

#endif
