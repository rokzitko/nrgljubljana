#ifndef _time_mem_h_
#define _time_mem_h_

// Warning: not thread safe!
class Timing {
private:
  double all, timer;
  bool running;
  map<string, double> t;
  bool reportonexit;
public:  
  Timing(bool _reportonexit = true) : reportonexit(_reportonexit) {
     all = gettime();
     timer = gettime();
     running = false;
  }
  void start(void) {
    my_assert(!running);
    running = true;
    timer = gettime();
  }
  double stop(void) {
    my_assert(running);
    running = false;
    double end = gettime();
    return end-timer;
  }
  void add(string timer) {
    if (!t.count(timer)) 
       t[timer] = 0.0; // initialize
    t[timer] += stop();
  }
  double value(string timer) const {
    if (!t.count(timer)) 
       return 0.0;
    return t.find(timer)->second;
  }
  double total() {
     double end = gettime();
     return end-all;
  }   
  void report(void) {
    const int T_WIDTH = 12;
    const double t_all = total();
    cout << endl;
    cout << "Timing report [" << myrank() << "]" << endl;
    cout << "=============" << endl;
    cout << setw(T_WIDTH) << "All" << ": " << long(t_all) << " s" << endl;
    double t_sum = 0.0;
    for (const auto &i : t) {
       // Only show those that contribute more than 1% of the total time!
       if ( (i.second/t_all) > 0.01) {
	  cout << setw(T_WIDTH) << i.first << ": " << long(i.second) << " s" << endl;
	  if (i.first[0] != '*') 
	     t_sum += i.second;
       }
    }
     cout << setw(T_WIDTH) << "Other" << ": " << long(t_all-t_sum) << " s" << endl;
  }
  ~Timing() { if (reportonexit) report(); }
};

// Higher-level timing code: time a section for as long as the object
// is in scope.
class TimeScope 
{
 private:
   Timing & timer;
   string timer_name;
 public:
   TimeScope(Timing & _timer, const string &_timer_name) :
    timer(_timer), timer_name(_timer_name) {
       timer.start();
    }
   ~TimeScope() {
      timer.add(timer_name);
   }
};

#define TIME(timer_name) TimeScope timer(t, timer_name)
#define TIME_SECTION(timer_name, section)  { TIME(timer_name); section }

// Stores maximal memory usage at various breakpoints.
// This is useful for estimating memory requirements at various
// points of the execution path.

class MemoryStats {
private:
  map<string, int> maxvals;
  bool reportonexit;
  int peakusage;
public:
   MemoryStats(bool rpexit = true) : reportonexit(rpexit) {
    peakusage = 0;
  }
  // Intermediate level routine.
  int used() {
     const int memused = memoryused();
     peakusage = max(peakusage, memused);
     return memused; 
  }
  // Sample memory usage at an arbitrarily named "breakpoint".
  int check(string breakpoint) {
    const int memused = used();
    if (maxvals.count(breakpoint) == 0) 
      maxvals[breakpoint] = 0;
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
     for (const auto &i : maxvals)
	topusage = max(topusage, i.second);
     if (topusage != 0) {
	for (const auto &i : maxvals)
	   F << setw(MS_WIDTH) << i.first << ": " << i.second << " kB" << endl;
     }
     my_assert(topusage <= peakusage);
     F << endl << "Peak usage: " << peakusage << " kB" << endl;
#endif
  }   
  ~MemoryStats() { if (reportonexit) report(); }
};

#endif
