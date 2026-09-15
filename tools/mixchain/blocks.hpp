// Channel-mixing discretization for NRG
// ** Blocks: sets of channels that Gamma couples

#ifndef _mixchain_blocks_hpp_
#define _mixchain_blocks_hpp_

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include "load.hpp"
#include "types.hpp"

namespace NRG::MixChain {

// BLOCKS
//
// If Gamma_ij vanishes exactly at every frequency for all i and j in different sets of channels, Gamma is block
// diagonal up to a permutation of the channels, and each block is an independent problem: it can be discretized on a
// mesh of its own and mapped onto a chain of its own. A diagonal Gamma is the extreme case, with one block per channel.
//
// Only exact zeros separate blocks. A coupling that is small but nonzero is part of the input.

using Block  = std::vector<int>;   // channel indices, 0-based and ascending
using Blocks = std::vector<Block>; // ordered by the first channel of each block

namespace detail {

// The connected components of the graph on 0..channels-1 with an edge wherever linked(i, j), for i < j, is true.
template<typename Linked> Blocks connected_blocks(const int channels, Linked &&linked) {
  std::vector<int> parent(static_cast<std::size_t>(channels));
  std::iota(parent.begin(), parent.end(), 0);
  const auto root = [&parent](int i) {
    while (parent[static_cast<std::size_t>(i)] != i) i = parent[static_cast<std::size_t>(i)];
    return i;
  };
  for (int i = 0; i < channels; i++)
    for (int j = i + 1; j < channels; j++)
      if (linked(i, j)) {
        const auto a = root(i);
        const auto b = root(j);
        if (a != b) parent[static_cast<std::size_t>(std::max(a, b))] = std::min(a, b);
      }

  // Channels are visited in ascending order, so each block is ascending and the blocks come ordered by their first
  // channel.
  Blocks blocks;
  std::vector<int> index_of_root(static_cast<std::size_t>(channels), -1);
  for (int i = 0; i < channels; i++) {
    auto &index = index_of_root[static_cast<std::size_t>(root(i))];
    if (index < 0) {
      index = static_cast<int>(blocks.size());
      blocks.emplace_back();
    }
    blocks[static_cast<std::size_t>(index)].push_back(i);
  }
  return blocks;
}

} // namespace detail

// The blocks of Gamma: channels i and j are in the same block if Gamma_ij is nonzero at some node of either frequency
// branch, or if a chain of such elements connects them.
template<typename S> Blocks gamma_blocks(const GammaInput<S> &input) {
  const auto nonzero_somewhere = [](const GammaBranch<S> &branch, const int i, const int j) {
    return std::any_of(branch.gamma.begin(), branch.gamma.end(),
                       [i, j](const Matrix<S> &gamma) { return gamma(i, j) != S(0); });
  };
  return detail::connected_blocks(input.channels, [&](const int i, const int j) {
    return nonzero_somewhere(input.pos, i, j) || nonzero_somewhere(input.neg, i, j);
  });
}

// Gamma restricted to the channels of one block, on the same frequency grids.
template<typename S> GammaInput<S> restrict_input(const GammaInput<S> &input, const Block &block) {
  const auto size     = static_cast<Eigen::Index>(block.size());
  const auto restrict = [&](const GammaBranch<S> &branch) {
    GammaBranch<S> result;
    result.omega     = branch.omega;
    result.innermost = branch.innermost;
    result.gamma.reserve(branch.size());
    for (const auto &gamma : branch.gamma) {
      Matrix<S> part(size, size);
      for (Eigen::Index i = 0; i < size; i++)
        for (Eigen::Index j = 0; j < size; j++)
          part(i, j) = gamma(block[static_cast<std::size_t>(i)], block[static_cast<std::size_t>(j)]);
      result.gamma.push_back(part);
    }
    return result;
  };
  GammaInput<S> result;
  result.channels = static_cast<int>(block.size());
  result.pos      = restrict(input.pos);
  result.neg      = restrict(input.neg);
  return result;
}

// "{1,3} {2,4}": 1-based, as in the names of the input files.
inline std::string blocks_name(const Blocks &blocks) {
  std::string text;
  for (const auto &block : blocks) {
    if (!text.empty()) text += ' ';
    text += '{';
    for (std::size_t k = 0; k < block.size(); k++) text += (k ? "," : "") + std::to_string(block[k] + 1);
    text += '}';
  }
  return text;
}

// The inverse of blocks_name(), for a Gamma with the given number of channels. The blocks must partition the channels,
// each channel appearing exactly once; they are returned in the canonical order of connected_blocks().
inline Blocks parse_blocks(const std::string &text, const int channels) {
  const auto fail = [&text](const std::string &why) {
    throw std::invalid_argument("Invalid block list '" + text + "': " + why + ".");
  };
  const auto skip_space = [&text](std::size_t position) {
    while (position < text.size() && std::isspace(static_cast<unsigned char>(text[position]))) position++;
    return position;
  };

  Blocks blocks;
  std::vector<bool> seen(static_cast<std::size_t>(channels), false);
  auto position = skip_space(0);
  while (position < text.size()) {
    if (text[position] != '{') fail("expected '{'");
    position = skip_space(position + 1);
    Block block;
    while (true) {
      const auto start = position;
      while (position < text.size() && std::isdigit(static_cast<unsigned char>(text[position]))) position++;
      if (position == start) fail("expected a channel number");
      if (position - start > 3) fail("channel " + text.substr(start, position - start) + " is out of range");
      const auto channel = std::stoi(text.substr(start, position - start));
      if (channel < 1 || channel > channels)
        fail("channel " + std::to_string(channel) + " is not in 1.." + std::to_string(channels));
      if (seen[static_cast<std::size_t>(channel - 1)]) fail("channel " + std::to_string(channel) + " appears twice");
      seen[static_cast<std::size_t>(channel - 1)] = true;
      block.push_back(channel - 1);
      position = skip_space(position);
      if (position < text.size() && text[position] == ',') {
        position = skip_space(position + 1);
        continue;
      }
      if (position < text.size() && text[position] == '}') break;
      fail("expected ',' or '}'");
    }
    std::sort(block.begin(), block.end());
    blocks.push_back(block);
    position = skip_space(position + 1);
  }
  for (int i = 0; i < channels; i++)
    if (!seen[static_cast<std::size_t>(i)]) fail("channel " + std::to_string(i + 1) + " is missing");
  std::sort(blocks.begin(), blocks.end(), [](const Block &a, const Block &b) { return a.front() < b.front(); });
  return blocks;
}

} // namespace NRG::MixChain

#endif
