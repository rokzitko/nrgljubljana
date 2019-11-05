/*
   Parser/calculator tool for Mathematica matrices
   Part of "NRG Ljubljana"
   Rok Zitko, rok.zitko@ijs.si, June 2009
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
double gammapolch(int);
double coefxi(int, int);
double coefzeta(int, int);
void yyerror(const char *);
int yylex();

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
%token GAMMAPOLCH COEFZETA COEFXI

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
	|     COEFZETA   '[' INTEGER ',' INTEGER ']' { $$ = coefzeta($3, $5); }
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

void parse_param(int argc, char *argv[])
{
  char c;
  
  while ((c = getopt(argc, argv, "vpP")) != -1) {
    switch (c) {
      case 'v':
        verbose = true;
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

int discretization_loaded = 0;

double theta = -1;
vector<double> xi;
vector<double> zeta;

void load_vector(const char *filename, vector<double> &v)
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
  ifstream THETA("theta.dat");
  if (!THETA) {
    cerr << "Can't open theta.dat" << endl;
    exit(1);
  }
  THETA >> theta;
  THETA.close();
  if (verbose) {
    cerr << "theta=" << theta << endl;
  }
  assert(theta >= 0);
  
  load_vector("xi.dat", xi);
  load_vector("zeta.dat", zeta);

  discretization_loaded = 1;
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
 
 load_discretization();

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

double gammapolch(int i)
{
  assert(discretization_loaded);
  // In initial.m, gammaPolCh[] is defined as sqrt(theta/pi * Gamma) !!
  return sqrt(theta/M_PI);
}

// channel nr. has offset 1
double coefzeta(int ch, int i)
{
  assert(discretization_loaded);
  return zeta[i];
}

double coefxi(int ch, int i)
{
  assert(discretization_loaded);
  return xi[i];
}

void yyerror(const char *error)
{
  cerr << error << endl;
}

