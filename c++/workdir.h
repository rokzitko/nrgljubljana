#ifndef _workdir_h_
#define _workdir_h_

#include <memory>
#include <string>
using namespace std::string_literals;
#include <cstring> // strncpy
#include <cstdlib> // mkdtemp, getenv
#include "portabil.h" // remove(std::string)

inline const auto default_workdir{"."s};

class Workdir {
 private:
   std::string workdir {};
   bool remove_at_exit {true}; // XXX: tie to P.removefiles?

 public:
   Workdir(const std::string &dir) {
     const auto workdir_template = dir + "/XXXXXX";
     const auto len = workdir_template.length()+1;
     auto x = std::make_unique<char[]>(len);
     strncpy(x.get(), workdir_template.c_str(), len);
     if (char *w = mkdtemp(x.get())) // create a unique directory
       workdir = w;
     else
       workdir = default_workdir;
     std::cout << "workdir=" << workdir << std::endl << std::endl;
   }
   std::string rhofn(const std::string &fn, const int N) const {  // density matrix files
     return workdir + "/" + fn + std::to_string(N); 
   }
   std::string unitaryfn(size_t N, const std::string &filename = "unitary"s) const { // eigenstates files
     return workdir + "/" + filename + std::to_string(N); 
   }
   void remove() {
     if (workdir != "")
       ::remove(workdir);
   }
   ~Workdir() {
     if (remove_at_exit)
       remove();
   }
};

inline auto set_workdir(const std::string &dir_) {
  std::string dir = default_workdir;
  if (const char *env_w = std::getenv("NRG_WORKDIR")) dir = env_w;
  if (!dir_.empty()) dir = dir_;
  return Workdir(dir);
}

#endif
