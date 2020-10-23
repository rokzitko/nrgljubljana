#ifndef _workdir_h_
#define _workdir_h_

#include <memory>
#include <string>
using namespace std::string_literals;
#include <optional>
#include <cstring> // strncpy
#include <cstdlib> // mkdtemp, getenv
#include "portabil.h" // remove(std::string)
#include <cstdio> // C remove()

inline const auto default_workdir{"."s};

// create a unique directory
inline auto dtemp(const std::string &path)
{
  const auto workdir_template = path + "/XXXXXX";
  const auto len = workdir_template.length()+1;
  auto x = std::make_unique<char[]>(len);
  strncpy(x.get(), workdir_template.c_str(), len);
  char *w = mkdtemp(x.get());
  return w ? std::optional<std::string>(w) : std::nullopt;
}

inline int remove(const std::string &filename) { return remove(filename.c_str()); }

class Workdir {
 private:
   const std::string workdir {};
   bool remove_at_exit {true}; // XXX: tie to P.removefiles?
 public:
   Workdir(const std::string &dir) : workdir(dtemp(dir).value_or(default_workdir)) {
     std::cout << "workdir=" << workdir << std::endl << std::endl;
   }
   std::string rhofn(const int N, const std::string &filename) const {  // density matrix files
     return workdir + "/" + filename + std::to_string(N);
   }
   std::string unitaryfn(const int N, const std::string &filename = "unitary"s) const { // eigenstates files
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
