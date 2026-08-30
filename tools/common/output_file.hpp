#ifndef _tools_common_output_file_hpp_
#define _tools_common_output_file_hpp_

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <system_error>

namespace NRG::Tools {

inline auto normalized_file_path(const std::filesystem::path &input) {
  std::error_code error;
  const auto normalized = std::filesystem::weakly_canonical(std::filesystem::absolute(input), error);
  if (!error) return normalized;
  if (error != std::errc::no_such_file_or_directory)
    throw std::filesystem::filesystem_error("Unable to normalize file path", input, error);
  return std::filesystem::weakly_canonical(std::filesystem::absolute(input).parent_path()) / input.filename();
}

inline auto ascii_lowercase(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](const unsigned char character) { return static_cast<char>(std::tolower(character)); });
  return value;
}

inline auto files_refer_to_same_location(const std::filesystem::path &first, const std::filesystem::path &second) {
  if (std::filesystem::exists(first) && std::filesystem::exists(second)) {
    std::error_code error;
    if (std::filesystem::equivalent(first, second, error)) return true;
    if (error) throw std::filesystem::filesystem_error("Unable to compare file paths", first, second, error);
  }

  const auto normalized_first = normalized_file_path(first);
  const auto normalized_second = normalized_file_path(second);
  if (normalized_first == normalized_second) return true;
  return normalized_first.parent_path() == normalized_second.parent_path()
         && ascii_lowercase(normalized_first.filename().string())
              == ascii_lowercase(normalized_second.filename().string());
}

inline void finish_output(std::ostream &output, const std::string &destination) {
  try {
    output.flush();
  } catch (const std::ios_base::failure &error) {
    throw std::runtime_error(destination + ": output write or flush failed: " + error.what());
  }
  if (!output) throw std::runtime_error(destination + ": output write or flush failed.");
}

inline void write_output_file(const std::filesystem::path &destination, const std::string &contents) {
  if (destination == "/dev/stdout") {
    std::cout << contents;
    finish_output(std::cout, "<stdout>");
    return;
  }

  std::ofstream output(destination);
  if (!output) throw std::runtime_error("Can't open " + destination.string() + " for writing.");
  output << contents;
  finish_output(output, destination.string());
  try {
    output.close();
  } catch (const std::ios_base::failure &error) {
    throw std::runtime_error(destination.string() + ": output close failed: " + error.what());
  }
  if (!output) throw std::runtime_error(destination.string() + ": output close failed.");
}

} // namespace NRG::Tools

#endif
