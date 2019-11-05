%{
#include <iostream>
#include <cstdlib>
#include "parser.hh"
#include "matrix.h"
#include <math.h>

using namespace std;

extern bool verbose;
extern bool veryverbose;

extern "C" int yywrap(void);

void parse_file(const char *);

int count_braces = 0;

%}

%x incl
%x sys

ALPHANAME ([A-Za-z][_A-Za-z0-9]*)
TEXT ([A-Za-z][^\n]*)

%%
^#[^\n]*  ; /* Ignore comment lines */

gammaPolCh { return GAMMAPOLCH; }
coefxi     { return COEFXI; }
coefzeta   { return COEFZETA; }
parse      { BEGIN(incl); }

<incl>[ \t]* ;
<incl>[^ \t\n]+ {
	 parse_file(yytext);
	 BEGIN(INITIAL);
       }

^! { BEGIN(sys); }

<sys>[ \t]* ;
<sys>[^ \t\n]+ {
         system(yytext);
         BEGIN(INITIAL);
       }

[0-9]+ { 
         yylval.ival = atoi(yytext);
	 return INTEGER;
       }

([1-9]e[-+]?[0-9]+) {
  yylval.dval = atof(yytext);
  return NUMBER;
  }

([0-9]+|([0-9]*\.[0-9]*)([eE][-+]?[0-9]+)?) {
		yylval.dval = atof(yytext);
		return NUMBER;
	}


[ \t]	;		 /* ignore white space */

\[{ALPHANAME}\] ;          /* ignore sections (param) */

{ALPHANAME}={TEXT} { 
                     if (veryverbose) {
                       cerr << yytext << endl;
		     } 
		}

{ALPHANAME}	{	/* return symbol pointer */
		struct symtab *sp = symlook(yytext);

		yylval.symp = sp;
		return NAME;
	}

<<EOF>> {
           if (verbose) {
             cerr << "EOF" << endl;
	   }
           yypop_buffer_state();

	   if (!YY_CURRENT_BUFFER) {
	     yyterminate();
	     return EXIT;
 	   }
	}

"$"	{ return 0; }

"{"     { count_braces++;
          return yytext[0];
	}
	
"}"     { if (count_braces == 0) {
             cerr << "Brace mismatch." << endl;
	     exit(1);
	  }
          count_braces--;
          return yytext[0];
	}

\n	{
          if (count_braces > 0) {
	    // ignore line breaks when lists span multiple lines
	  } else {
	    return yytext[0];
	  }
	}
.	return yytext[0];
%%

// The following function is specific to flex lexer.
// (Called from parser, see matrix.y *)
void parse_file(const char *filename)
{
   printf("parsing \"%s\"\n", filename);
   FILE *F = fopen(filename, "r");
   if (!F) {
     printf("Can't open %s for reading.\n", filename);
     exit(1);
   }
   YY_BUFFER_STATE buffer = yy_new_buffer(F, YY_BUF_SIZE);
   if (!buffer) {
     printf("yy_new_buffer() failed.\n");
     exit(1);
   }
   yypush_buffer_state(buffer);
}
