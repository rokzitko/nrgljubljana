/*
   Parser/calculator tool for Mathematica matrices
   Part of "NRG Ljubljana"
   Rok Zitko, rok.zitko@ijs.si, 2009-2025
*/

%{
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cassert>
#include <unistd.h>
#include <vector>

#include "matrix.h"

using namespace std;

struct symtab symtab[NSYMS];

void dump_vector(struct vec *dvec);
void dump_matrix(struct mat *dmat);

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
	|       vector                  { OUT << prefix; dump_vector($1); }
	|       matrix                  { OUT << prefix; dump_matrix($1); }
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
   char *dup = strdup(name);
   struct symtab *sp = symlook(dup);
   sp->funcptr = func;
}

void usage()
{
  std::cout << "Usage: matrix [-h] [-vV] [-c channels] [-p | -P] <file1> <file2> ..." << std::endl;
}

void parse_param(int argc, char *argv[])
{
  char c;
  while (c = getopt(argc, argv, "hc:vpPs"), c != -1) {
  switch (c) {
    case 'h':
      usage();
      exit(EXIT_SUCCESS);
  
    case 'v':
       verbose = true;
       break;
  
    case 'c':
       nrchannels = atoi(optarg);
       numberedch = true;
       break;
	
    case 's': // superconducting case with Nambu structure 
      sc = true;
      break;
  
      case 'V':
        veryverbose = true;
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

void load_vector(string filename, vector<double> &v)
{
  ifstream F(filename);
  if (!F) {
    cerr << "Can't open " << filename << " for reading." << endl;
    exit(1);
  }
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
}

void load_discretization()
{
  for (int ch = 1; ch <= nrchannels ; ch++) {
    if (verbose)
      cerr << "Channel " << ch << endl;
    string suffix = (numberedch ? to_string(ch) : "") + ".dat";
    string fntheta = "theta" + suffix;
    ifstream THETA(fntheta);
    if (!THETA) {
      cerr << "Can't open " << fntheta << endl;
      exit(1);
    }
    THETA >> theta[ch-1];
    THETA.close();
    if (verbose) {
      cerr << "theta[" << ch << "]=" << theta[ch-1] << endl;
    }
    assert(theta[ch-1] >= 0);
  
    string fnxi = "xi" + suffix;
    load_vector(fnxi, xi[ch-1]);
    string fnzeta = "zeta" + suffix;
    load_vector(fnzeta, zeta[ch-1]);
  }
}

void load_discretization_sc()
{
  for (int ch = 1; ch <= nrchannels ; ch++) {
    if (verbose)
      cerr << "Channel " << ch << endl;
    string suffix = (numberedch ? to_string(ch) : "") + ".dat";
    V[ch-1].resize(2);
    for (int i = 1; i <= 2; i++) {
      V[ch-1][i-1].resize(2);
      for (int j = 1; j <= 2; j++) {
        string fnV = "V" + to_string(i) + to_string(j) + suffix;
        ifstream FV(fnV);
        if (!FV) {
          cerr << "Can't open " << fnV << endl;
          exit(1);
        }
        FV >> V[ch-1][i-1][j-1];
        FV.close();
        if (verbose) {
          cerr << "V[" << ch << "](" << i << ", " << j << ")=" << V[ch-1][i-1][j-1] << endl;
        }
      }
    }  
  
    string fnxi = "xi" + suffix;
    load_vector(fnxi, xi[ch-1]);
    string fnzeta = "zeta" + suffix;
    load_vector(fnzeta, zeta[ch-1]);
    string fndelta = "scdelta" + suffix;
    load_vector(fndelta, delta[ch-1]);
    string fnkappa = "sckappa" + suffix;
    load_vector(fnkappa, kappa[ch-1]);
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
    if (current != 0) {
      fclose(file);
    }
    
    char *filename = filelist[current];
    if (current == remaining) {
      // we're done!
      return 1;
    }
    
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

int main(int argc, char *argv[])
{
 cout << setprecision(18);

 parse_param(argc, argv);
 
 theta.resize(nrchannels);
 zeta.resize(nrchannels);
 xi.resize(nrchannels);
 delta.resize(nrchannels);
 kappa.resize(nrchannels);
 V.resize(nrchannels);
 
 if (!sc) {
   load_discretization();
 } else {
   load_discretization_sc();
 }

 addfunc((const char *)"Sqrt", sqrt);
 addfunc((const char *)"sqrt", sqrt);
 addfunc((const char *)"exp", exp);
 addfunc((const char *)"log", log);

 remaining = argc-optind; // arguments left = filenames!
 filelist = argv+optind; // copy the pointer
 
 if (remaining != 0) {
   // Read from files!
   yywrap();
 }

 yyparse();

 if (verbose) {
  cerr << "DONE!" << endl;
 }

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

struct vec * duplicate_vector(struct vec *dvec)
{
  struct vec *ptr = dvec;
  struct vec *previous = 0;
  
  while (ptr) {
    struct vec *new_node = new(struct vec);
    if (previous) {
      previous->next = new_node;
    }
    previous = new_node;
    new_node->val = ptr->val;
    new_node->next = 0;

    ptr = ptr->next;
  }
  
  return 0;
}

double gammapolch(int ch)
{
  // In initial.m, gammaPolCh[] is defined as sqrt(theta/pi * Gamma) !!
  assert(1 <= ch && ch <= nrchannels);
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
  assert(i < delta[ch-1].size());
  return delta[ch-1][i];
}

double coefkappa(int ch, int i)
{
  assert(1 <= ch && ch <= nrchannels);
  assert(i < kappa[ch-1].size());
  return kappa[ch-1][i];
}

// channel number has offset 1
// i,j have offset 1 (unlike in zeta/xi/... which are 0 based)
double coefV(int i, int j)
{
  const int ch = 1;
  assert(1 <= i && i <= V[ch-1].size());
  assert(1 <= j && j <= V[ch-1][i-1].size());
  return V[ch-1][i-1][j-1];
}

void yyerror(const char *error)
{
  cerr << error << endl;
}

