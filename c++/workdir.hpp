#ifndef _workdir_hpp_
#define _workdir_hpp_

#include <memory>
#include <string>
#include <optional>
#include <stdexcept>
#include <filesystem>
#include <iostream>
#include <cstring> // strncpy
#include <cstdlib> // mkdtemp, getenv
#include "portabil.hpp" // remove(std::string)
#include <cstdio> // C remove()

namespace NRG {

using namespace std::string_literals;

inline const auto default_workdir{"."s};

// create a unique directory
inline auto dtemp(const std::string &path, const std::string &pattern = "/XXXXXX"s)
{
  const auto workdir_template = path + pattern;
  const auto len = workdir_template.length()+1;
  auto x = std::make_unique<char[]>(len);
  strncpy(x.get(), workdir_template.c_str(), len);
  char *w = mkdtemp(x.get());
  return w ? std::optional<std::string>(w) : std::nullopt;
}

// Note: This will remove a directory only if it is empty!
inline int remove(const std::string &filename) { return std::remove(filename.c_str()); }

class Workdir {
 private:
   const std::string workdir {};
   bool remove_at_exit {true}; // XXX: tie to P.removefiles?
  public:
    explicit Workdir(const std::string &dir, const bool quiet = false) : workdir([&dir]() {
      const auto temp = dtemp(dir);
      if (!temp) throw std::runtime_error("Failed to create temporary workdir in " + dir);
      return temp.value();
    }()) {
      if (!quiet) std::cout << "workdir=" << workdir << std::endl << std::endl;
    }
   explicit Workdir() : Workdir(default_workdir, true) {} // defaulted version (for testing purposes)
   Workdir(const Workdir &) = delete;
   Workdir(Workdir &&) = delete;
   Workdir & operator=(const Workdir &) = delete;
   Workdir & operator=(Workdir &&) = delete;
   [[nodiscard]] auto get() const { return workdir; }
   [[nodiscard]] auto rhofn(const size_t N, const std::string &filename) const {  // density matrix files
     return workdir + "/" + filename + std::to_string(N);
   }
    [[nodiscard]] auto unitaryfn(const size_t N, const std::string &filename = "unitary"s) const { // eigenstates files
      return workdir + "/" + filename + std::to_string(N);
    }
     bool remove_workdir() {
       if (workdir == "") return true;
       std::error_code ec;
       std::filesystem::remove_all(workdir, ec);
       if (ec) {
         std::cerr << "Failed to remove workdir " << workdir << ": " << ec.message() << std::endl;
         return false;
       }
       return true;
     }
   ~Workdir() {
     if (remove_at_exit) static_cast<void>(remove_workdir());
   }
};

inline auto set_workdir(const std::string &dir_) {
  std::string dir = default_workdir;
  if (const char *env_w = std::getenv("NRG_WORKDIR")) dir = env_w;
  if (!dir_.empty()) dir = dir_;
  return std::make_unique<Workdir>(dir);
}

} // namespace

#endif
