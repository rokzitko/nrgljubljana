#ifndef _workdir_hpp_
#define _workdir_hpp_

#include <cstddef>
#include <memory>
#include <ostream>
#include <string>
#include <optional>
#include <stdexcept>
#include <filesystem>
#include <iostream>
#include <system_error>
#include <cstring> // strncpy
#include <cstdlib> // mkdtemp, getenv
#include <cerrno>
#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>
#include "portabil.hpp" // remove(std::string)
#include <cstdio> // C remove()

namespace NRG {

using namespace std::string_literals;

inline const auto default_workdir{"."s};

enum class WorkdirMode { unique_temporary, persistent_exact };

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
   const bool remove_at_exit {true};
   const int lock_fd {-1};

   static auto create(const std::string &dir, const WorkdirMode mode) {
     if (mode == WorkdirMode::unique_temporary) {
       const auto temp = dtemp(dir);
       if (!temp) throw std::runtime_error("Failed to create temporary workdir in " + dir);
       return temp.value();
     }

     if (dir.empty()) throw std::runtime_error("Persistent exact workdir path must not be empty");
     std::error_code ec;
     static_cast<void>(std::filesystem::create_directories(dir, ec));
     if (ec || !std::filesystem::is_directory(dir, ec)) {
       const auto reason = ec ? ": " + ec.message() : "";
       throw std::runtime_error("Failed to create or open persistent exact workdir " + dir + reason);
     }
     return dir;
   }

   static auto acquire_lock(const std::string &dir, const WorkdirMode mode) {
     if (mode != WorkdirMode::persistent_exact) return -1;
     const auto filename = dir + "/.nrg.lock";
     const int fd = open(filename.c_str(), O_CREAT | O_RDWR | O_CLOEXEC, 0600);
     if (fd == -1) throw std::runtime_error("Failed to open checkpoint lock " + filename + ": " + std::strerror(errno));
     if (flock(fd, LOCK_EX | LOCK_NB) == -1) {
       const auto reason = std::string(std::strerror(errno));
       close(fd);
       throw std::runtime_error("Checkpoint directory is already in use: " + dir + ": " + reason);
     }
     return fd;
   }

  public:
   Workdir(const std::string &dir, const WorkdirMode mode, const bool quiet = false) :
     workdir(create(dir, mode)), remove_at_exit(mode == WorkdirMode::unique_temporary), lock_fd(acquire_lock(workdir, mode)) {
     if (!quiet) std::cout << "workdir=" << workdir << std::endl << std::endl;
   }
   explicit Workdir(const std::string &dir, const bool quiet = false) : Workdir(dir, WorkdirMode::unique_temporary, quiet) {}
   explicit Workdir() : Workdir(default_workdir, true) {} // defaulted version (for testing purposes)
   Workdir(const Workdir &) = delete;
   Workdir(Workdir &&) = delete;
   Workdir & operator=(const Workdir &) = delete;
   Workdir & operator=(Workdir &&) = delete;
   [[nodiscard]] auto get() const { return workdir; }
   [[nodiscard]] bool persistent() const noexcept { return !remove_at_exit; }
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
     if (lock_fd != -1) close(lock_fd);
   }
};

inline auto set_workdir(const std::string &dir_, const WorkdirMode mode) {
  if (mode == WorkdirMode::persistent_exact) return std::make_unique<Workdir>(dir_, mode);

  std::string dir = default_workdir;
  if (const char *env_w = std::getenv("NRG_WORKDIR")) dir = env_w;
  if (!dir_.empty()) dir = dir_;
  return std::make_unique<Workdir>(dir, mode);
}

inline auto set_workdir(const std::string &dir_) { return set_workdir(dir_, WorkdirMode::unique_temporary); }

} // namespace

#endif
