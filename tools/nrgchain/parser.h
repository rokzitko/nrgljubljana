// Discretization ODE solver for NRG
//
// ** Parsing of the parameter file
//
// Rok Zitko, zitko@theorie.physik.uni-goettingen.de, Dec 2008
// $Id: parser.h,v 1.1 2009/03/20 09:53:41 rok Exp $

map<string, string> params;

// Locate block [name] in a file stream. Returns true if succeessful.
bool find_block(ifstream &F, const string &s)
{
  string target = "[" + s + "]";
  F.clear();
  F.seekg(0, ios::beg);
  while (F) {
    string line;
    getline(F, line);
    if (F && target.compare(line) == 0) {
      break;
    }
  }
  return bool(F); // True if found.
}

// Return a parameter of type double, use default value if not found.
double P(const string &keyword, double def)
{
  if (params.count(keyword) == 0)
    return def;

  return atof(params[keyword]);
}

int Pint(const string &keyword, int def)
{
  if (params.count(keyword) == 0)
    return def;

  return atoi(params[keyword]);
}

string Pstr(const string &keyword, string def)
{
  if (params.count(keyword) == 0)
    return def;

  return params[keyword];
}

bool Pbool(const string &keyword, bool def)
{
  if (params.count(keyword) == 0)
    return def;

  return params[keyword] == "true";
}

// Parse a block of "keyword=value" lines.
void parse_block(ifstream &F)
{
  while (F) {
    string line;
    getline(F, line);
    if (!F) {
      break;
    }

    if (line[0] == '[') // new block, we're done!
      break;

    if (line.length() == 0) // skip empty lines
      continue;

    if (line[0] == '#') // skip comment lines
      continue;

    string::size_type pos_eq = line.find_first_of('=');
    if (pos_eq == string::npos) // not found
      continue;

    string keyword = line.substr(0, pos_eq);
    string value = line.substr(pos_eq+1);

    params[keyword] = value;
  }
}

void parser(const string &filename)
{
  ifstream F;
  safe_open(F, filename);

  if (find_block(F, "param")) {
    parse_block(F);
  }
}
