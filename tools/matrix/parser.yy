/*
   Parser/calculator tool for Mathematica matrices
   Part of "NRG Ljubljana"
   Rok Zitko, rok.zitko@ijs.si
*/

%{
#include <cstddef>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ios>
#include <istream>
#include <ostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cassert>
#include <unistd.h>
#include <string>
#include <vector>

#include "matrix.h"

using namespace std;

struct symtab symtab[NSYMS];

void dump_vector(struct vec *dvec);
void dump_matrix(struct mat *dmat);
void free_vector(struct vec *dvec);
void free_matrix(struct mat *dmat);
void clear_symtab();

// channel index - 0 based
// Wilson chain site index - 0 based
// V - 0 based (variables are V11, V12, etc.)

double gammapolch(int);
double coefxi(int, int);
double coefzeta(int, int);
double coefdelta(int, int);
double coefkappa(int, int);
double coefV(int, int); // Nambu indexes; channel index not implemented yet
void yyerror(const char *);
int yylex();

bool numberedch = false; // number suffix?
bool sc = false;
int nrchannels = 1;

ostream & OUT = cout;
string prefix = "";

bool verbose = false;
bool veryverbose = false;
%}

%union {
	double dval;
	int ival;
	char *str;
	struct vec *dvec;
	struct mat *dmat;
	struct symtab *symp;
}
%token <str> STRING
%token <symp> NAME
%token <dval> NUMBER
%token <ival> INTEGER
%left '-' '+'
%left '*' '/'
%left ','
%nonassoc UMINUS

%token PARSE EXIT
%token GAMMAPOLCH COEFZETA COEFXI COEFDELTA COEFKAPPA COEFV

%type <dval> expression
%type <dvec> expressionlist
%type <dvec> vector
%type <dmat> vectorlist
%type <dmat> matrix
%%
statement_list:	statement '\n'
	|	statement_list statement '\n'
	;

statement:	NAME '=' expression	{
                   if (veryverbose) {
		     cerr << "Defining " << $1->name << "=" << $3 << endl;
		   }
                   $1->value = $3; 
		}
	|	expression		{ OUT << prefix << $1 << endl; }
	|       vector                  { OUT << prefix; dump_vector($1); free_vector($1); }
	|       matrix                  { OUT << prefix; dump_matrix($1); free_matrix($1); }
        |       EXIT {
	           /* We are done! */
		   exit(1);
		}
        |          /* null statement (empty line) */
	;

expression:	expression '+' expression { $$ = $1 + $3; }
	|	expression '-' expression { $$ = $1 - $3; }
	|	expression '*' expression { $$ = $1 * $3; }
	|	expression '/' expression
				{	if($3 == 0.0)
						yyerror("divide by zero");
					else
						$$ = $1 / $3;
				}
	|	'-' expression %prec UMINUS	{ $$ = -$2; }
	|	'(' expression ')'	{ $$ = $2; }
	|	NUMBER
	|       INTEGER                 { $$ = $1; }
	|	NAME			{ $$ = $1->value; }
	|	NAME '(' expression ')'	{
			if($1->funcptr)
				$$ = ($1->funcptr)($3);
			else {
				cerr << $1->name << " not a function" << endl;
				$$ = 0.0;
			}
		}
	|	NAME '[' expression ']'	{
			if($1->funcptr)
				$$ = ($1->funcptr)($3);
			else {
			        cerr << $1->name << " not a function" << endl;
				$$ = 0.0;
			}
		}
	;

expression:   GAMMAPOLCH '[' INTEGER ']' { $$ = gammapolch($3); }
        |     COEFXI     '[' INTEGER ',' INTEGER ']' { $$ = coefxi($3, $5); } 
	|     COEFZETA         '[' INTEGER ',' INTEGER ']' { $$ = coefzeta($3, $5); }
	|     COEFDELTA      '[' INTEGER ',' INTEGER ']' { $$ = coefdelta($3, $5); }
	|     COEFKAPPA      '[' INTEGER ',' INTEGER ']' { $$ = coefkappa($3, $5); }
	|     COEFV            '[' INTEGER ',' INTEGER ']' { $$ = coefV($3, $5); }
	;

expressionlist: expression  {
                $$ = new(struct vec);
		$$->val = $1;
    		$$->next = 0;
            }
  | expression ',' expressionlist {
                struct vec *new_node = new(struct vec);
		new_node->val = $1;
		new_node->next = $3;
		$$ = new_node;
	    }
  ;

vector: '{' '}'    { $$ = 0; }
  | '{' expressionlist '}' { $$ = $2; }
  ;
  
vectorlist: vector {
               $$ = new(struct mat);
	       $$->vec = $1;
	       $$->next = 0;
         }
  | vector ',' vectorlist {
                struct mat *new_node = new(struct mat);
		new_node->vec = $1;
		new_node->next = $3;
		$$ = new_node;
	 }
  ;
  
matrix: '{' '{' '}' '}' { $$ = 0; }
  | '{' vectorlist '}' { $$ = $2; }
  ;

%%
#include <common/version.hpp>

#include "../common/diagnostics.hpp"

/* look up a symbol table entry, add if not present */
struct symtab * symlook(char *s)
{
	struct symtab *sp;
	
	for(sp = symtab; sp < &symtab[NSYMS]; sp++) {
		/* is it already here? */
		if(sp->name && !strcmp(sp->name, s))
			return sp;
		
		/* is it free */
		if(!sp->name) {
			sp->name = strdup(s);
			return sp;
		}
		/* otherwise continue to next */
	}
	yyerror("Too many symbols");
	exit(1);	/* cannot continue */
} /* symlook */

void addfunc(const char *name, double (*func)(double))
{
   struct symtab *sp = symlook(const_cast<char *>(name));
   sp->funcptr = func;
}

void usage()
{
  std::cout << "Usage: matrix [options] <file1> <file2> ...\n"
            << "  -h             show this help\n"
            << "  -v             show resolved configuration and verbose diagnostics on standard error\n"
            << "  -vv            also show detailed parser diagnostics\n"
            << "  -V, --version  show project version\n"
            << "  -c channels    use numbered coefficient files for this many channels\n"
            << "  -s             accepted for compatibility; every coefficient file that exists is now read\n"
            << "  -p | -P        omit or add '= ' to output lines" << std::endl;
}

void parse_param(int argc, char *argv[])
{
  int c;
  while ((c = getopt(argc, argv, "hc:vpPs")) != -1) {
  switch (c) {
    case 'h':
      usage();
      exit(EXIT_SUCCESS);
  
    case 'v':
       if (verbose) veryverbose = true;
       verbose = true;
       break;
  
    case 'c':
       nrchannels = atoi(optarg);
       numberedch = true;
       break;
	
    case 's': // accepted for compatibility: the tables are now loaded by what is on disk, not by this switch
      sc = true;
      break;
  
      case 'p':
        prefix = "";
	break;
	
      case 'P':
        prefix = "= ";
	break;
	
      default:
        abort();
    }
  }
}

vector<double> theta;
vector<vector<double>> xi;
vector<vector<double>> zeta;
vector<vector<double>> delta;
vector<vector<double>> kappa;
vector<vector<vector<double>>> V;

// Which optional tables were found, per channel. xi and zeta are not here: every consumer reads them, so a missing
// one is an error while loading. The rest are required only where an expression uses them, which is what lets a
// template that has no pairing terms run without scdelta and sckappa, and one that reads coefV run without theta.
vector<char> has_theta;
vector<char> has_delta;
vector<char> has_kappa;
vector<char> has_V;

string coeffile(const string &name, int ch)
{
  return name + (numberedch ? to_string(ch) : "") + ".dat";
}

// The table an expression asks for must have been loaded. Called from the accessors rather than from the loader, so
// that the message names both the file and the symbol that wanted it.
void require_table(bool loaded, const string &filename, const char *symbol)
{
  if (!loaded) {
    cerr << "matrix: " << filename << " is needed by " << symbol << " but was not found." << endl;
    exit(1);
  }
}

bool load_vector(string filename, vector<double> &v)
{
  ifstream F(filename);
  if (!F) return false;
  while (F.good()) {
    double x;
    F >> x;
    if (!F.fail()) {
      v.push_back(x);
      if (verbose) {
        cerr << filename << " " << x << endl;
      }
    }
  }
  F.close();
  return true;
}

bool load_scalar(string filename, double &value)
{
  ifstream F(filename);
  if (!F) return false;
  F >> value;
  if (F.fail()) {
    cerr << "Can't read a number from " << filename << "." << endl;
    exit(1);
  }
  F.close();
  if (verbose) {
    cerr << filename << " " << value << endl;
  }
  return true;
}

void load_required_vector(string filename, vector<double> &v)
{
  if (!load_vector(filename, v)) {
    cerr << "Can't open " << filename << " for reading." << endl;
    exit(1);
  }
}

// Every coefficient table of every channel, whichever of them exist. Which ones a run needs depends on the
// expressions it evaluates, not on a command line switch: a template with no Z block reads coefV and no pairing
// tables, while a superconducting one reads scdelta and sckappa as well.
void load_discretization()
{
  for (int ch = 1; ch <= nrchannels ; ch++) {
    if (verbose)
      cerr << "Channel " << ch << endl;
    const string suffix = (numberedch ? to_string(ch) : "") + ".dat";

    load_required_vector("xi" + suffix, xi[ch-1]);
    load_required_vector("zeta" + suffix, zeta[ch-1]);

    has_theta[ch-1] = load_scalar("theta" + suffix, theta[ch-1]);
    if (has_theta[ch-1]) assert(theta[ch-1] >= 0);

    has_delta[ch-1] = load_vector("scdelta" + suffix, delta[ch-1]);
    has_kappa[ch-1] = load_vector("sckappa" + suffix, kappa[ch-1]);

    // The Nambu structure is all four elements or none of them: a partial set is a staging mistake rather than a
    // run that does not need V.
    V[ch-1].resize(2);
    int found = 0;
    for (int i = 1; i <= 2; i++) {
      V[ch-1][i-1].resize(2);
      for (int j = 1; j <= 2; j++)
        if (load_scalar("V" + to_string(i) + to_string(j) + suffix, V[ch-1][i-1][j-1])) found++;
    }
    if (found != 0 && found != 4) {
      cerr << "matrix: channel " << ch << " has " << found << " of the four V{i}{j}" << suffix
           << " files; write all of them or none." << endl;
      exit(1);
    }
    has_V[ch-1] = (found == 4);
  }
}

// Global variables with filenames
char **filelist;
int remaining;
int current = 0;
FILE *file;
extern FILE *yyin;

extern "C" {
int yywrap(void)
{
    if (current != 0 && file) {
      fclose(file);
      file = nullptr;
    }

    if (current == remaining) {
      // we're done!
      return 1;
    }

    char *filename = filelist[current];
    
    if (verbose) {
      cerr << "Parsing " << filename << endl;
    }
    file = fopen(filename, "r");
    if (!file) {
      cerr << "Can't open " << filename << " for reading." << endl;
      exit(1);
    }
    yyin = file;
    current++;
    return 0;
}
}

void report_configuration()
{
  if (!verbose) return;
  NRG::Tools::ConfigurationReport report("matrix");
  report.value("verbosity", veryverbose ? 2 : 1);
  // -s no longer selects which tables are read; it is reported only so that a caller passing it sees that it was seen.
  report.value("superconducting_switch", sc);
  report.value("channels", nrchannels);
  report.value("numbered_coefficient_files", numberedch);
  report.value("output.prefix", prefix.empty() ? "none" : "equals");
  report.value("output.precision", 18);
  report.resolved("input.mode", remaining == 0 ? "stdin" : "files", "positional arguments");
  report.value("input.count", remaining);
  if (remaining == 0) {
    report.resolved("input.0", "<stdin>", "no input files");
  } else {
    for (int i = 0; i < remaining; ++i) report.value("input." + to_string(i + 1), filelist[i]);
  }

  for (int ch = 1; ch <= nrchannels; ++ch) {
    const auto base = "channel." + to_string(ch) + ".";
    const auto suffix = (numberedch ? to_string(ch) : "") + ".dat";
    report.value(base + "xi.file", "xi" + suffix);
    report.value(base + "xi.count", xi[ch - 1].size());
    if (!xi[ch - 1].empty()) report.resolved(base + "xi.max_index", xi[ch - 1].size() - 1, "loaded coefficient file");
    report.value(base + "zeta.file", "zeta" + suffix);
    report.value(base + "zeta.count", zeta[ch - 1].size());
    if (!zeta[ch - 1].empty()) report.resolved(base + "zeta.max_index", zeta[ch - 1].size() - 1, "loaded coefficient file");
    // Each optional table is reported as it was found, so the report says what this run has to work with.
    report.value(base + "theta.file", "theta" + suffix);
    if (has_theta[ch - 1])
      report.value(base + "theta", theta[ch - 1]);
    else
      report.value(base + "theta", "not found");
    report.value(base + "delta.file", "scdelta" + suffix);
    report.value(base + "delta.count", has_delta[ch - 1] ? to_string(delta[ch - 1].size()) : "not found");
    report.value(base + "kappa.file", "sckappa" + suffix);
    report.value(base + "kappa.count", has_kappa[ch - 1] ? to_string(kappa[ch - 1].size()) : "not found");
    for (int i = 1; i <= 2; ++i)
      for (int j = 1; j <= 2; ++j) {
        const auto name = base + "V" + to_string(i) + to_string(j);
        report.value(name + ".file", "V" + to_string(i) + to_string(j) + suffix);
        if (has_V[ch - 1])
          report.value(name, V[ch - 1][i - 1][j - 1]);
        else
          report.value(name, "not found");
      }
  }
  report.write(cerr);
}

int main(int argc, char *argv[])
{
 if (NRG::Tools::report_version_if_requested(argc, argv, "matrix")) return EXIT_SUCCESS;
 cout << setprecision(18);

 parse_param(argc, argv);
 remaining = argc-optind; // arguments left = filenames!
 filelist = argv+optind; // copy the pointer
 
 theta.resize(nrchannels);
 zeta.resize(nrchannels);
 xi.resize(nrchannels);
 delta.resize(nrchannels);
 kappa.resize(nrchannels);
 V.resize(nrchannels);
 has_theta.assign(nrchannels, 0);
 has_delta.assign(nrchannels, 0);
 has_kappa.assign(nrchannels, 0);
 has_V.assign(nrchannels, 0);

 load_discretization();

 report_configuration();

 addfunc((const char *)"Sqrt", sqrt);
 addfunc((const char *)"sqrt", sqrt);
 addfunc((const char *)"exp", exp);
 addfunc((const char *)"log", log);

 if (remaining != 0) {
   // Read from files!
   yywrap();
 }

 yyparse();

 if (verbose) {
  cerr << "DONE!" << endl;
 }

  clear_symtab();

 return 0;
}

void dump_vector(struct vec *dvec)
{
  struct vec *ptr = dvec;
  
  while (ptr) {
    OUT << ptr->val << (ptr->next != 0 ? " " : "");
    ptr = ptr->next;
  }
  OUT << endl;
}

void dump_matrix(struct mat *dmat)
{
  struct mat *ptr = dmat;
  
  while (ptr) {
    dump_vector(ptr->vec);
    if (ptr->next != 0) { 
      OUT << prefix;
    }
    ptr = ptr->next;
  }
}

void free_vector(struct vec *dvec)
{
  while (dvec) {
    auto *next = dvec->next;
    delete dvec;
    dvec = next;
  }
}

void free_matrix(struct mat *dmat)
{
  while (dmat) {
    free_vector(dmat->vec);
    auto *next = dmat->next;
    delete dmat;
    dmat = next;
  }
}

void clear_symtab()
{
  for (auto &entry : symtab) {
    free(entry.name);
    entry.name = nullptr;
    entry.funcptr = nullptr;
    entry.value = 0.0;
  }
}

struct vec * duplicate_vector(struct vec *dvec)
{
  struct vec *ptr = dvec;
  struct vec *head = 0;
  struct vec *previous = 0;
  
  while (ptr) {
    struct vec *new_node = new(struct vec);
    if (!head) {
      head = new_node;
    }
    if (previous) {
      previous->next = new_node;
    }
    previous = new_node;
    new_node->val = ptr->val;
    new_node->next = 0;

    ptr = ptr->next;
  }
  
  return head;
}

double gammapolch(int ch)
{
  // In initial.m, gammaPolCh[] is defined as sqrt(theta/pi * Gamma) !!
  assert(1 <= ch && ch <= nrchannels);
  require_table(has_theta[ch-1], coeffile("theta", ch), "gammaPolCh");
  return sqrt(theta[ch-1]/M_PI);
}

// channel nr. has offset 1, i is 0 based
double coefzeta(int ch, int i)
{
  assert(1 <= ch && ch <= nrchannels);
  assert(i < zeta[ch-1].size());
  return zeta[ch-1][i];
}

double coefxi(int ch, int i)
{
  assert(1 <= ch && ch <= nrchannels);
  assert(i < xi[ch-1].size());
  return xi[ch-1][i];
}

double coefdelta(int ch, int i)
{
  assert(1 <= ch && ch <= nrchannels);
  require_table(has_delta[ch-1], coeffile("scdelta", ch), "coefdelta");
  assert(i < delta[ch-1].size());
  return delta[ch-1][i];
}

double coefkappa(int ch, int i)
{
  assert(1 <= ch && ch <= nrchannels);
  require_table(has_kappa[ch-1], coeffile("sckappa", ch), "coefkappa");
  assert(i < kappa[ch-1].size());
  return kappa[ch-1][i];
}

// channel number has offset 1
// i,j have offset 1 (unlike in zeta/xi/... which are 0 based)
double coefV(int i, int j)
{
  const int ch = 1;
  require_table(has_V[ch-1], coeffile("V11", ch), "coefV");
  assert(1 <= i && i <= V[ch-1].size());
  assert(1 <= j && j <= V[ch-1][i-1].size());
  return V[ch-1][i-1][j-1];
}

void yyerror(const char *error)
{
  cerr << error << endl;
}
