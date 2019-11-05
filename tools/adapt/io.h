// Discretization ODE solver for NRG
//
// ** Input/output code
//
// Rok Zitko, zitko@theorie.physik.uni-goettingen.de, Dec 2008
// $Id: io.h,v 1.1 2009/03/20 09:53:41 rok Exp $

inline double atof(const string &s)
{
  return atof(s.c_str());
}

inline int atoi(const string &s)
{
  return atoi(s.c_str());
}

void safe_open(ifstream &F, const string &filename)
{
  F.open(filename.c_str());
  if (!F) {
    cerr << "Can't open " << filename << " for reading." << endl;
    exit(1);
  }
}

const int PREC = 16;

void safe_open(ofstream &F, const string &filename)
{
  F.open(filename.c_str());
  if (!F) {
    cerr << "Can't open " << filename << " for writing." << endl;
    exit(1);
  }
  F << setprecision(PREC);
}

// Get next line from stream F, skipping empty lines and comments.
string getnextline(ifstream &F)
{
  string line;
  while (F) {
    getline(F, line);
    if (!F) // bail out
      break;

    if (line.length() == 0) // skip empty lines
      continue;

    if (line[0] == '#') // skip comment lines
      continue;

    return line;
  }
  return ""; // error
}
