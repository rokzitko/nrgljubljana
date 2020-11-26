#ifndef _workdir_hpp_
#define _workdir_hpp_

#include <memory>
#include <string>
#include <optional>
#include <cstring> // strncpy
#include <cstdlib> // mkdtemp, getenv
#include "portabil.hpp" // remove(std::string)
#include <cstdio> // C remove()

namespace NRG {

using namespace std::string_literals;

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

inline int remove(const std::string &filename) { return std::remove(filename.c_str()); }

class Workdir {
 private:
   const std::string workdir {};
   bool remove_at_exit {true}; // XXX: tie to P.removefiles?
 public:
   explicit Workdir(const std::string &dir) : workdir(dtemp(dir).value_or(default_workdir)) {
     std::cout << "workdir=" << workdir << std::endl << std::endl;
   }
   Workdir(const Workdir &) = delete;
   Workdir(Workdir &&) = delete;
   Workdir & operator=(const Workdir &) = delete;
   Workdir & operator=(Workdir &&) = delete;
   [[nodiscard]] auto rhofn(const size_t N, const std::string &filename) const {  // density matrix files
     return workdir + "/" + filename + std::to_string(N);
   }
   [[nodiscard]] auto unitaryfn(const size_t N, const std::string &filename = "unitary"s) const { // eigenstates files
     return workdir + "/" + filename + std::to_string(N);
   }
   void remove_workdir() {
     if (workdir != "") NRG::remove(workdir);
   }
  ~Workdir() {
    if (remove_at_exit) remove_workdir();
  }
};

inline auto set_workdir(const std::string &dir_) {
  std::string dir = default_workdir;
  if (const char *env_w = std::getenv("NRG_WORKDIR")) dir = env_w;
  if (!dir_.empty()) dir = dir_;
  return Workdir(dir);
}

} // namespace

#endif
