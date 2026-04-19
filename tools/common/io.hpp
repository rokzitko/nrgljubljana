#ifndef _tools_common_io_hpp_
#define _tools_common_io_hpp_

#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

namespace NRG::Tools {

inline void open_input(std::ifstream &stream, const std::string &filename) {
  stream.open(filename.c_str());
  if (!stream)
    throw std::runtime_error("Can't open " + filename + " for reading.");
}

inline void open_output(std::ofstream &stream, const std::string &filename, const int precision = 16) {
  stream.open(filename.c_str());
  if (!stream)
    throw std::runtime_error("Can't open " + filename + " for writing.");
  stream << std::setprecision(precision);
}

template<typename Stream>
std::string next_data_line(Stream &stream) {
  std::string line;
  while (stream) {
    std::getline(stream, line);
    if (!stream) break;
    if (line.empty()) continue;
    if (line[0] == '#') continue;
    return line;
  }
  return "";
}

} // namespace NRG::Tools

#endif
