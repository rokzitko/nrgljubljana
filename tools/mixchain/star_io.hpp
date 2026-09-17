// Channel-mixing discretization for NRG
// ** Reading and writing of the star file

#ifndef _mixchain_star_io_hpp_
#define _mixchain_star_io_hpp_

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <ios>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

#include "../common/io.hpp"
#include "../common/tabulated.hpp"
#include "blocks.hpp"
#include "star.hpp"
#include "types.hpp"

namespace NRG::MixChain {

// THE STAR FILE
//
// A text file describing the star Hamiltonian
//
//   H = sum_k E_k c_k^dag c_k + sum_k sum_i ( v_{k,i} d_i^dag c_k + h.c. ),
//
// with one row per bath level k. It is written by the 's' mode and read by the 'l' mode, and it is the only thing
// that passes between them.
//
// Lines beginning with '#' are comments. Exactly one of them is read back: the line carrying 'channels=', whose
// whitespace-separated key=value pairs are the header. Unknown keys are ignored, which is what makes the diagnostic
// comments and the column legend harmless. Everything else is a data row.
//
// Header keys:
//
//   channels     the dimension N of Gamma, and the number of components of every coupling vector
//   mMAX         the largest interval index; the file holds 2*channels*(mMAX+1) rows
//   z            the twist parameter of the mesh
//   Lambda       the discretization parameter
//   bandrescale  the band rescaling that was applied to Gamma when it was read; the energies below are in the
//                rescaled band, whose edge is 1
//   complex      1 if the coupling vectors are complex, 0 if they are real. It fixes the number of columns.
//   untabulated  optional: 'from,to', the part of the band the mesh reaches where the input is not tabulated,
//                from < |omega| < to in the rescaled band, between the accumulation point of the mesh and the
//                innermost tabulated frequency, the widest over the frequency branches and blocks; or 'none' when
//                the mesh reaches no such part, as when it accumulates at a gap edge. There the star follows the
//                constant continuation of the input, and the chain stage reports from which site the chain samples
//                it. Absent means unknown.
//
// Blocks line, written right after the header when Gamma was discretized in more than one block (see blocks.hpp):
//
//   # blocks= {1,3} {2}
//
// Each group is one block, with the channels numbered from 1 as in the input file names. Without the line the star is
// a single block of all channels, which is what a Gamma that does not split, or split_blocks=false, gives. The chain
// stage maps each block onto a chain of its own. On loading, every level must couple only to channels of the block
// its branch label belongs to.
//
// Columns of a data row:
//
//   m      the interval index, 0..mMAX. The interval is [eps(z+m+2), eps(z+m+1)] on that branch's mesh, so a larger
//          m lies closer to the Fermi level
//   sign   the frequency branch the level came from: '+' for omega>0 and '-' for omega<0
//   a      the eigenvalue branch of Gamma the level came from, 0..channels-1. This is a label from the branch
//          tracking, not a channel: a branch is an eigenvector direction of Gamma, which in general points across
//          several channels and rotates with omega. With blocks, the branches of the first block are numbered first,
//          then those of the next, so the first block owns the labels 0..size-1
//   E      the representative energy of that interval and branch, carrying the sign of its frequency branch
//   v_i    the coupling vector in the channel basis: component i is the amplitude between impurity orbital i and
//          this bath level. For complex=1 each component is written as a pair, Re_v_i Im_v_i.
//
//          v = sqrt(w) u, with w the weight of that branch over that interval and u its eigenvector at E, in the
//          normalization of the input Gamma. Neither the weight nor the vector is renormalized, so
//          sum_k v_k v_k^dag is the hybridisation weight Theta of the star as it was built.
//
// m, sign and a are not used to build the chain, which sees only the pairs (E, v). They are written so that the file
// can be read by a person and checked on loading: the row count, the ranges of m and a, and the agreement of the
// sign column with the sign of E.
//
// The diagnostics of the star stage are written as comments and are not read back: they are properties of Gamma
// rather than of the star, so a star loaded from a file has them empty. With several blocks, each diagnostic line
// starts with the block it belongs to, as in "# block {1,3}: max_cquad_error=...".

inline constexpr auto star_default_filename = "star.dat";

struct StarHeader {
  int channels{};
  unsigned int mMAX{};
  double z{};
  double Lambda{};
  double bandrescale{1.0};
  bool complex_data{};
  double untabulated_from{};
  double untabulated_to{};
  bool untabulated_known{}; // false if the file does not record it
  Blocks blocks;            // empty if the file has no blocks line
};

namespace detail {

inline auto read_lines(const std::string &filename) {
  std::ifstream F;
  NRG::Tools::open_input(F, filename);
  std::vector<std::string> lines;
  std::string line;
  while (std::getline(F, line)) lines.push_back(line);
  if (F.bad()) throw std::runtime_error("Error reading " + filename + ".");
  return lines;
}

inline auto is_comment(const std::string &line) {
  const auto first = line.find_first_not_of(" \t");
  return first != std::string::npos && line[first] == '#';
}

// The number of columns of a data row: m, sign, a, E and then the components of the coupling vector.
inline auto star_columns(const int channels, const bool complex_data) {
  return static_cast<std::size_t>(4 + channels * (complex_data ? 2 : 1));
}

// An integer field, which save_star() writes as a plain decimal integer. Anything else, such as "1.5" or "1e3", is
// rejected rather than truncated.
inline int parse_integer(const std::string &text, const std::string &what, const std::string &filename) {
  int value        = 0;
  const auto *last = text.data() + text.size();
  const auto [end, error] = std::from_chars(text.data(), last, value);
  if (error != std::errc{} || end != last)
    throw std::runtime_error(filename + ": " + what + " must be an integer, not '" + text + "'.");
  return value;
}

inline void parse_header_line(const std::string &line, StarHeader &header, const std::string &filename) {
  const auto fields = NRG::Tools::split_fields(line);
  bool have_channels = false, have_mmax = false, have_lambda = false, have_z = false, have_complex = false;
  for (const auto &field : fields) {
    const auto separator = field.find('=');
    if (separator == std::string::npos) continue;
    const auto key   = field.substr(0, separator);
    const auto value = field.substr(separator + 1);
    const auto number = [&] { return NRG::Tools::parse_tabulated_double(value); };
    if (key == "channels") {
      header.channels = parse_integer(value, "channels", filename);
      have_channels   = true;
    } else if (key == "mMAX") {
      // Checked before the conversion to unsigned, where a negative value would wrap around.
      const auto mMAX = parse_integer(value, "mMAX", filename);
      if (mMAX < 1) throw std::runtime_error(filename + ": mMAX must be greater than 0.");
      header.mMAX = static_cast<unsigned int>(mMAX);
      have_mmax   = true;
    } else if (key == "z") {
      header.z = number();
      have_z   = true;
    } else if (key == "Lambda") {
      header.Lambda = number();
      have_lambda   = true;
    } else if (key == "bandrescale") {
      header.bandrescale = number();
    } else if (key == "untabulated") {
      header.untabulated_known = true;
      if (value != "none") {
        const auto comma = value.find(',');
        if (comma == std::string::npos)
          throw std::runtime_error(filename + ": untabulated must be 'from,to' or 'none', not '" + value + "'.");
        header.untabulated_from = NRG::Tools::parse_tabulated_double(value.substr(0, comma));
        header.untabulated_to   = NRG::Tools::parse_tabulated_double(value.substr(comma + 1));
        if (!(header.untabulated_to >= header.untabulated_from && header.untabulated_from >= 0.0))
          throw std::runtime_error(filename + ": untabulated must satisfy 0 <= from <= to, not '" + value + "'.");
      }
    } else if (key == "complex") {
      const auto flag = parse_integer(value, "complex", filename);
      if (flag != 0 && flag != 1) throw std::runtime_error(filename + ": complex must be 0 or 1.");
      header.complex_data = flag == 1;
      have_complex        = true;
    }
    // Unknown keys belong to the diagnostic comments and are ignored.
  }
  if (!(have_channels && have_mmax && have_z && have_lambda && have_complex))
    throw std::runtime_error(filename + ": the star header must give channels, mMAX, z, Lambda and complex.");
  if (header.channels < 1 || header.channels > max_channels)
    throw std::runtime_error(filename + ": channels must be between 1 and " + std::to_string(max_channels) + ".");
  if (!(header.Lambda > 1.0)) throw std::runtime_error(filename + ": Lambda must be greater than 1.");
  if (!(header.z > 0.0 && header.z <= 1.0)) throw std::runtime_error(filename + ": z must be in (0,1].");
  if (!(std::isfinite(header.bandrescale) && header.bandrescale > 0.0))
    throw std::runtime_error(filename + ": bandrescale must be a positive finite number.");
}

} // namespace detail

namespace detail {

inline constexpr auto blocks_key = "blocks=";

// The text after "blocks=" if the line is the blocks line, which is recognised by that key as its first field.
inline std::optional<std::string> blocks_value(const std::string &line) {
  const auto start = line.find_first_not_of(" \t", line.find('#') + 1);
  if (start == std::string::npos || line.compare(start, std::string(blocks_key).size(), blocks_key) != 0)
    return std::nullopt;
  return line.substr(start + std::string(blocks_key).size());
}

} // namespace detail

// The header alone, with the blocks line if there is one. Call this first: 'complex' decides the scalar type the star
// must be loaded with.
inline auto read_star_header(const std::string &filename) {
  std::optional<std::string> header_line, blocks_line;
  for (const auto &line : detail::read_lines(filename)) {
    if (!detail::is_comment(line)) continue;
    if (const auto value = detail::blocks_value(line)) {
      if (blocks_line) throw std::runtime_error(filename + ": more than one blocks line.");
      blocks_line = value;
    } else if (!header_line && line.find("channels=") != std::string::npos) {
      header_line = line;
    }
  }
  if (!header_line) throw std::runtime_error(filename + ": no star header found.");

  StarHeader header;
  detail::parse_header_line(*header_line, header, filename);
  if (blocks_line) {
    try {
      header.blocks = parse_blocks(*blocks_line, header.channels);
    } catch (const std::invalid_argument &error) { throw std::runtime_error(filename + ": " + error.what()); }
  }
  return header;
}

inline auto star_is_complex(const std::string &filename) { return read_star_header(filename).complex_data; }

template<typename S> void save_star(const Star<S> &star, std::ostream &out) {
  out << std::setprecision(18);
  out << "# mixchain star" << std::endl;
  out << "# channels=" << star.channels << " mMAX=" << star.mMAX << " z=" << star.z << " Lambda=" << star.Lambda
      << " bandrescale=" << star.bandrescale << " complex=" << (is_complex_v<S> ? 1 : 0);
  if (star.untabulated_known) {
    out << " untabulated=";
    if (star.untabulated_to > star.untabulated_from)
      out << star.untabulated_from << "," << star.untabulated_to;
    else
      out << "none";
  }
  out << std::endl;
  if (star.blocks.size() > 1) out << "# " << detail::blocks_key << " " << blocks_name(star.blocks) << std::endl;
  for (std::size_t b = 0; b < star.diagnostics.size(); b++) {
    const auto &diagnostics = star.diagnostics[b];
    // With a single block the lines carry no prefix.
    const auto block = star.blocks.size() > 1 ? "block " + blocks_name({star.blocks[b]}) + ": " : std::string();
    out << "# " << block << "max_interval_deviation=" << diagnostics.max_interval_deviation
        << " at_omega=" << diagnostics.max_interval_omega << std::endl;
    out << "# " << block << "max_cquad_error=" << diagnostics.max_cquad_error << std::endl;
    for (const auto &[name, crossings] : {std::pair{"crossings_pos", &diagnostics.crossings_pos},
                                          std::pair{"crossings_neg", &diagnostics.crossings_neg}}) {
      out << "# " << block << name << "=" << crossings->size();
      for (const auto omega : *crossings) out << " " << omega;
      out << std::endl;
    }
  }
  out << "# m sign a E";
  for (int i = 1; i <= star.channels; i++) {
    if constexpr (is_complex_v<S>)
      out << " Re_v" << i << " Im_v" << i;
    else
      out << " v" << i;
  }
  out << std::endl;

  for (const auto &level : star.levels) {
    out << level.m << " " << (level.sign == Sign::POS ? '+' : '-') << " " << level.branch << " " << level.energy;
    for (int i = 0; i < star.channels; i++) {
      const auto component = level.coupling(i);
      if constexpr (is_complex_v<S>)
        out << " " << std::real(component) << " " << std::imag(component);
      else
        out << " " << component;
    }
    out << std::endl;
  }
}

template<typename S> void save_star(const Star<S> &star, const std::string &filename) {
  std::ofstream F;
  NRG::Tools::open_output(F, filename, 18);
  save_star(star, F);
  F.close();
  if (!F) throw std::runtime_error("Error writing " + filename + ".");
}

// Load a star written by save_star(). Theta is recomputed from the rows in file order, so it is identical to the
// value the star was built with; theta_exact and the diagnostics are not stored and are left empty.
template<typename S> auto load_star(const std::string &filename) {
  const auto header = read_star_header(filename);
  if (header.complex_data != is_complex_v<S>)
    throw std::logic_error(filename + ": the star was written in "
                           + (header.complex_data ? "complex" : "real") + " arithmetic.");

  const auto dimension = static_cast<Eigen::Index>(header.channels);
  Star<S> star;
  star.channels    = header.channels;
  star.mMAX        = header.mMAX;
  star.z           = header.z;
  star.Lambda      = header.Lambda;
  star.bandrescale     = header.bandrescale;
  star.untabulated_from  = header.untabulated_from;
  star.untabulated_to    = header.untabulated_to;
  star.untabulated_known = header.untabulated_known;
  star.blocks          = header.blocks;
  if (star.blocks.empty()) {
    star.blocks.emplace_back(static_cast<std::size_t>(header.channels));
    std::iota(star.blocks.front().begin(), star.blocks.front().end(), 0);
  }
  // The block that owns each branch label: the first block has the labels 0..size-1, the next one the following.
  std::vector<std::size_t> block_of_branch;
  for (std::size_t b = 0; b < star.blocks.size(); b++)
    block_of_branch.insert(block_of_branch.end(), star.blocks[b].size(), b);
  star.theta       = Matrix<S>::Zero(dimension, dimension);
  star.theta_exact = Matrix<S>::Zero(dimension, dimension);

  const auto columns = detail::star_columns(header.channels, header.complex_data);
  std::size_t number = 0;
  for (const auto &line : detail::read_lines(filename)) {
    if (detail::is_comment(line)) continue;
    const auto fields = NRG::Tools::split_fields(line);
    if (fields.empty()) continue;
    number++;
    if (fields.size() != columns)
      throw std::runtime_error(filename + ": row " + std::to_string(number) + " has " + std::to_string(fields.size())
                               + " columns instead of " + std::to_string(columns) + ".");

    StarLevel<S> level;
    level.m      = detail::parse_integer(fields[0], "row " + std::to_string(number) + ": the interval index", filename);
    level.branch = detail::parse_integer(fields[2], "row " + std::to_string(number) + ": the branch index", filename);
    level.energy = NRG::Tools::parse_tabulated_double(fields[3]);
    if (fields[1] == "+")
      level.sign = Sign::POS;
    else if (fields[1] == "-")
      level.sign = Sign::NEG;
    else
      throw std::runtime_error(filename + ": row " + std::to_string(number) + " has the sign '" + fields[1] + "'.");

    if (level.m < 0 || static_cast<unsigned int>(level.m) > header.mMAX)
      throw std::runtime_error(filename + ": row " + std::to_string(number) + " has the interval index "
                               + std::to_string(level.m) + ".");
    if (level.branch < 0 || level.branch >= header.channels)
      throw std::runtime_error(filename + ": row " + std::to_string(number) + " has the branch index "
                               + std::to_string(level.branch) + ".");
    if (!(std::isfinite(level.energy) && level.energy != 0.0))
      throw std::runtime_error(filename + ": row " + std::to_string(number) + " has a vanishing or non-finite "
                               "energy.");
    if ((level.energy > 0.0) != (level.sign == Sign::POS))
      throw std::runtime_error(filename + ": row " + std::to_string(number)
                               + " has an energy whose sign contradicts its frequency branch.");

    level.coupling = Vector<S>::Zero(dimension);
    for (int i = 0; i < header.channels; i++) {
      const auto first = static_cast<std::size_t>(4 + (header.complex_data ? 2 * i : i));
      const auto real_part = NRG::Tools::parse_tabulated_double(fields[first]);
      const auto imaginary_part =
        header.complex_data ? NRG::Tools::parse_tabulated_double(fields[first + 1]) : 0.0;
      level.coupling(i) = make_scalar<S>(real_part, imaginary_part);
    }
    const auto &block = star.blocks[block_of_branch[static_cast<std::size_t>(level.branch)]];
    for (int i = 0; i < header.channels; i++) {
      if (level.coupling(i) == S(0) || std::find(block.begin(), block.end(), i) != block.end()) continue;
      throw std::runtime_error(filename + ": row " + std::to_string(number) + " (branch " + std::to_string(level.branch)
                               + ", block " + blocks_name({block}) + ") couples to channel " + std::to_string(i + 1)
                               + ", which is outside its block.");
    }
    star.theta += level.coupling * level.coupling.adjoint();
    star.levels.push_back(std::move(level));
  }

  const auto expected = 2 * static_cast<std::size_t>(header.channels) * (header.mMAX + 1);
  if (star.levels.size() != expected)
    throw std::runtime_error(filename + ": expected " + std::to_string(expected) + " levels, found "
                             + std::to_string(star.levels.size()) + ".");
  return star;
}

} // namespace NRG::MixChain

#endif
