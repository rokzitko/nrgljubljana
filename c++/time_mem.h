#ifndef _time_mem_h_
#define _time_mem_h_

#include <utility>
#include <chrono>

// Returns a string with a floating value in fixed (non-exponential) format with 
// N digits of precision after the decimal point.
string prec(double x, int N)
{
  ostringstream s;
  s << std::fixed << std::setprecision(N) << x;
  return s.str();
}

string prec3(double x)
{
  return prec(x, 3);
}
   
namespace time_mem {

// Warning: not thread safe!
class Timing {
  private:
    using tp = chrono::steady_clock::time_point;
    using dp = chrono::duration<double>;
    tp start_time;     // time when Timing object constructed
    tp timer;          // start time for timing sections
    bool running;      // currently timing a section
    map<string, dp> t; // accumulators
  public:
  tp now() noexcept {
    return chrono::steady_clock::now();
  }
  Timing() {
    start_time = now();
    running    = false;
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
  void add(string timer) {
    t[timer] += stop();
  }
  dp total() {
    tp end_time = now();
    return end_time - start_time;
  }
  double total_in_seconds() {
    return total().count();
  }
  void report() {
    const int T_WIDTH  = 12;
    const dp t_all = total();
    cout << endl;
    cout << "Timing report [" << myrank() << "]" << endl;
    cout << "=============" << endl;
    cout << setw(T_WIDTH) << "All"
         << ": " << prec3(t_all.count()) << " s" << endl;
    dp t_sum;
    for (const auto &i : t) {
      // Only show those that contribute more than 1% of the total time!
      if (i.second/t_all > 0.01) {
        cout << setw(T_WIDTH) << i.first << ": " << prec3(i.second.count()) << " s" << endl;
        if (i.first[0] != '*') t_sum += i.second;
      }
    }
    cout << setw(T_WIDTH) << "Other"
         << ": " << prec3((t_all-t_sum).count()) << " s" << endl;
  }
};

// Higher-level timing code: time a section for as long as the object
// is in scope.
class TimeScope {
  private:
  Timing &timer;
  string timer_name;
  public:
  TimeScope(Timing &_timer, string _timer_name) : timer(_timer), timer_name(std::move(_timer_name)) { timer.start(); }
  ~TimeScope() { timer.add(timer_name); }
};

#define TIME(timer_name) time_mem::TimeScope timer(time_mem::tm, timer_name)
#define TIME_SECTION(timer_name, section)                                                                                                            \
  {                                                                                                                                                  \
    TIME(timer_name);                                                                                                                                \
    section                                                                                                                                          \
  }

// Stores maximal memory usage at various breakpoints.
// This is useful for estimating memory requirements at various
// points of the execution path.

class MemoryStats {
  private:
  map<string, int> maxvals;
  int peakusage;
  public:
  MemoryStats() { peakusage = 0; }
  int used() {
    const int memused = memoryused();
    peakusage         = max(peakusage, memused);
    return memused;
  }
  // Sample memory usage at an arbitrarily named "breakpoint".
  int check(string breakpoint) {
    const int memused = used();
    maxvals[breakpoint] = max(maxvals[breakpoint], memused);
    return memused;
  }
  void report(ostream &F = cout) const {
#ifdef HAS_MEMORY_USAGE
    const int MS_WIDTH = 12;
    F << endl;
    F << "Memory usage report [" << myrank() << "]" << endl;
    F << "===================" << endl;
    int topusage = 0; // top usage recorded by check()
    for (const auto &i : maxvals) topusage = max(topusage, i.second);
    if (topusage != 0) {
      for (const auto &i : maxvals) F << setw(MS_WIDTH) << i.first << ": " << i.second << " kB" << endl;
    }
    my_assert(topusage <= peakusage);
    F << endl << "Peak usage: " << peakusage << " kB" << endl;
#endif
  }
};

} // namespace time_mem
#endif
