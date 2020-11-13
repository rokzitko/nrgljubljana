// Discretization ODE solver for NRG
// ** Parsing of the parameter file

#ifndef _adapt_parser_hpp_
#define _adapt_parser_hpp_

#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>

namespace NRG::Adapt {

// Locate block [name] in a file stream. Returns true if succeessful.
inline bool find_block(std::ifstream &F, const std::string &s) {
  std::string target = "[" + s + "]";
  F.clear();
  F.seekg(0, std::ios::beg);
  while (F) {
    std::string line;
    std::getline(F, line);
    if (F && target.compare(line) == 0) { break; }
  }
  return bool(F); // True if found.
}

class Params : public std::map<std::string, std::string> {
 private:
   // Parse a block of "keyword=value" lines.
   void parse_block(std::ifstream &F) {
     while (F) {
       std::string line;
       std::getline(F, line);
       if (!F) { break; }
       if (line[0] == '[') // new block, we're done!
         break;
       if (line.length() == 0) // skip empty lines
         continue;
       if (line[0] == '#') // skip comment lines
         continue;
       const auto pos_eq = line.find_first_of('=');
       if (pos_eq == std::string::npos) // not found
         continue;
       const auto keyword = line.substr(0, pos_eq);
       const auto value   = line.substr(pos_eq + 1);
       this->insert_or_assign(keyword, value);
     }
   }
 public:
   Params(const std::string &filename) {
     std::ifstream F;
     safe_open(F, filename);
     if (find_block(F, "param")) { parse_block(F); }
   }
   // Return a parameter of type double, use default value if not found.
   auto P(const std::string &keyword, const double def) const {
     if (this->count(keyword) == 0) return def;
     return atof(this->at(keyword));
   }
   auto Pint(const std::string &keyword, const int def) const {
     if (this->count(keyword) == 0) return def;
     return atoi(this->at(keyword));
   }
   auto Pstr(const std::string &keyword, const std::string &def) const {
     if (this->count(keyword) == 0) return def;
     return this->at(keyword);
   }
   auto Pbool(const std::string &keyword, const bool def) const {
     if (this->count(keyword) == 0) return def;
     return this->at(keyword) == "true";     
   }
};

} // namespace

#endif
