/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.0.4"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* Copy the first part of user declarations.  */
#line 7 "parser.yy" /* yacc.c:339  */

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

#line 111 "parser.cc" /* yacc.c:339  */

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* In a future release of Bison, this section will be replaced
   by #include "parser.hh".  */
#ifndef YY_YY_PARSER_HH_INCLUDED
# define YY_YY_PARSER_HH_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    STRING = 258,
    NAME = 259,
    NUMBER = 260,
    INTEGER = 261,
    UMINUS = 262,
    PARSE = 263,
    EXIT = 264,
    GAMMAPOLCH = 265,
    COEFZETA = 266,
    COEFXI = 267,
    COEFDELTA = 268,
    COEFKAPPA = 269,
    COEFV = 270
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 52 "parser.yy" /* yacc.c:355  */

	double dval;
	int ival;
	char *str;
	struct vec *dvec;
	struct mat *dmat;
	struct symtab *symp;

#line 176 "parser.cc" /* yacc.c:355  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_PARSER_HH_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 193 "parser.cc" /* yacc.c:358  */

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef YY_ATTRIBUTE
# if (defined __GNUC__                                               \
      && (2 < __GNUC__ || (__GNUC__ == 2 && 96 <= __GNUC_MINOR__)))  \
     || defined __SUNPRO_C && 0x5110 <= __SUNPRO_C
#  define YY_ATTRIBUTE(Spec) __attribute__(Spec)
# else
#  define YY_ATTRIBUTE(Spec) /* empty */
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# define YY_ATTRIBUTE_PURE   YY_ATTRIBUTE ((__pure__))
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# define YY_ATTRIBUTE_UNUSED YY_ATTRIBUTE ((__unused__))
#endif

#if !defined _Noreturn \
     && (!defined __STDC_VERSION__ || __STDC_VERSION__ < 201112)
# if defined _MSC_VER && 1200 <= _MSC_VER
#  define _Noreturn __declspec (noreturn)
# else
#  define _Noreturn YY_ATTRIBUTE ((__noreturn__))
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN \
    _Pragma ("GCC diagnostic push") \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")\
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif


#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYSIZE_T yynewbytes;                                            \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / sizeof (*yyptr);                          \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  37
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   158

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  29
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  8
/* YYNRULES -- Number of rules.  */
#define YYNRULES  34
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  86

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   270

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      21,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      23,    24,     9,     8,    11,     7,     2,    10,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    22,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    25,     2,    26,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    27,     2,    28,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,    12,    13,    14,    15,    16,    17,    18,    19,
      20
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint8 yyrline[] =
{
       0,    78,    78,    79,    82,    88,    89,    90,    91,    95,
      98,    99,   100,   101,   107,   108,   109,   110,   111,   112,
     120,   130,   131,   132,   133,   134,   135,   138,   143,   151,
     152,   155,   160,   168,   169
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "STRING", "NAME", "NUMBER", "INTEGER",
  "'-'", "'+'", "'*'", "'/'", "','", "UMINUS", "PARSE", "EXIT",
  "GAMMAPOLCH", "COEFZETA", "COEFXI", "COEFDELTA", "COEFKAPPA", "COEFV",
  "'\\n'", "'='", "'('", "')'", "'['", "']'", "'{'", "'}'", "$accept",
  "statement_list", "statement", "expression", "expressionlist", "vector",
  "vectorlist", "matrix", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,    45,    43,    42,
      47,    44,   262,   263,   264,   265,   266,   267,   268,   269,
     270,    10,    61,    40,    41,    91,    93,   123,   125
};
# endif

#define YYPACT_NINF -19

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-19)))

#define YYTABLE_NINF -1

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     110,    15,   -19,   -19,   134,   -19,   -18,   -11,   -10,    17,
      32,    33,   134,    45,    16,    20,    92,   -19,   -19,   134,
     134,   134,   -17,   -19,    53,    60,    61,    63,    65,    77,
      -6,    70,   -19,    71,    56,    80,    64,   -19,    82,   -19,
     134,   134,   134,   134,    92,    46,     3,    78,   100,   101,
     108,   109,   111,   -19,    93,   134,   -19,    96,   -19,   -19,
      18,    18,   -19,   -19,   -19,   -19,   -19,   125,   126,   128,
     129,   130,   -19,   -19,    90,   -19,   116,   117,   118,   119,
     120,   -19,   -19,   -19,   -19,   -19
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       9,    18,    16,    17,     0,     8,     0,     0,     0,     0,
       0,     0,     0,     0,     9,     0,     5,     6,     7,     0,
       0,     0,    18,    14,     0,     0,     0,     0,     0,     0,
       0,     0,    29,    27,     0,    31,     0,     1,     0,     2,
       0,     0,     0,     0,     4,     0,     0,     0,     0,     0,
       0,     0,     0,    15,    29,     0,    30,     0,    34,     3,
      11,    10,    12,    13,    19,    20,    21,     0,     0,     0,
       0,     0,    33,    28,     0,    32,     0,     0,     0,     0,
       0,    23,    22,    24,    25,    26
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
     -19,   -19,   133,     5,   103,   -13,    91,   -19
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,    14,    15,    33,    34,    17,    36,    18
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_uint8 yytable[] =
{
      35,    40,    41,    42,    43,    16,    20,    24,    21,    23,
      40,    41,    42,    43,    25,    26,    37,    30,    53,    16,
       1,     2,     3,     4,    44,    45,    46,    42,    43,    65,
       5,     6,     7,     8,     9,    10,    11,    19,    20,    12,
      21,    39,    27,    13,    35,    60,    61,    62,    63,    22,
       2,     3,     4,    40,    41,    42,    43,    28,    29,    47,
       6,     7,     8,     9,    10,    11,    48,    49,    12,    50,
      64,    51,    31,    32,    22,     2,     3,     4,    40,    41,
      42,    43,    55,    52,    56,     6,     7,     8,     9,    10,
      11,    57,    58,    12,    22,     2,     3,     4,    54,    40,
      41,    42,    43,    59,    66,     6,     7,     8,     9,    10,
      11,    67,    68,    12,     1,     2,     3,     4,    32,    69,
      70,    72,    71,    74,     5,     6,     7,     8,     9,    10,
      11,    76,    77,    12,    78,    79,    80,    13,    22,     2,
       3,     4,    81,    82,    83,    84,    85,    38,    75,     6,
       7,     8,     9,    10,    11,     0,     0,    12,    73
};

static const yytype_int8 yycheck[] =
{
      13,     7,     8,     9,    10,     0,    23,    25,    25,     4,
       7,     8,     9,    10,    25,    25,     0,    12,    24,    14,
       4,     5,     6,     7,    19,    20,    21,     9,    10,    26,
      14,    15,    16,    17,    18,    19,    20,    22,    23,    23,
      25,    21,    25,    27,    57,    40,    41,    42,    43,     4,
       5,     6,     7,     7,     8,     9,    10,    25,    25,     6,
      15,    16,    17,    18,    19,    20,     6,     6,    23,     6,
      24,     6,    27,    28,     4,     5,     6,     7,     7,     8,
       9,    10,    11,     6,    28,    15,    16,    17,    18,    19,
      20,    11,    28,    23,     4,     5,     6,     7,    28,     7,
       8,     9,    10,    21,    26,    15,    16,    17,    18,    19,
      20,    11,    11,    23,     4,     5,     6,     7,    28,    11,
      11,    28,    11,    27,    14,    15,    16,    17,    18,    19,
      20,     6,     6,    23,     6,     6,     6,    27,     4,     5,
       6,     7,    26,    26,    26,    26,    26,    14,    57,    15,
      16,    17,    18,    19,    20,    -1,    -1,    23,    55
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     4,     5,     6,     7,    14,    15,    16,    17,    18,
      19,    20,    23,    27,    30,    31,    32,    34,    36,    22,
      23,    25,     4,    32,    25,    25,    25,    25,    25,    25,
      32,    27,    28,    32,    33,    34,    35,     0,    31,    21,
       7,     8,     9,    10,    32,    32,    32,     6,     6,     6,
       6,     6,     6,    24,    28,    11,    28,    11,    28,    21,
      32,    32,    32,    32,    24,    26,    26,    11,    11,    11,
      11,    11,    28,    33,    27,    35,     6,     6,     6,     6,
       6,    26,    26,    26,    26,    26
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    29,    30,    30,    31,    31,    31,    31,    31,    31,
      32,    32,    32,    32,    32,    32,    32,    32,    32,    32,
      32,    32,    32,    32,    32,    32,    32,    33,    33,    34,
      34,    35,    35,    36,    36
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     2,     3,     3,     1,     1,     1,     1,     0,
       3,     3,     3,     3,     2,     3,     1,     1,     1,     4,
       4,     4,     6,     6,     6,     6,     6,     1,     3,     2,
       3,     1,     3,     4,     3
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;                                                  \
    }                                                           \
while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256



/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)

/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  YYUSE (yytype);
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, int yyrule)
{
  unsigned long int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[yyssp[yyi + 1 - yynrhs]],
                       &(yyvsp[(yyi + 1) - (yynrhs)])
                                              );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
yystrlen (const char *yystr)
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
yystpcpy (char *yydest, const char *yysrc)
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
        switch (*++yyp)
          {
          case '\'':
          case ',':
            goto do_not_strip_quotes;

          case '\\':
            if (*++yyp != '\\')
              goto do_not_strip_quotes;
            /* Fall through.  */
          default:
            if (yyres)
              yyres[yyn] = *yyp;
            yyn++;
            break;

          case '"':
            if (yyres)
              yyres[yyn] = '\0';
            return yyn;
          }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
{
  YYUSE (yyvaluep);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

int
yyparse (void)
{
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        YYSTYPE *yyvs1 = yyvs;
        yytype_int16 *yyss1 = yyss;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * sizeof (*yyssp),
                    &yyvs1, yysize * sizeof (*yyvsp),
                    &yystacksize);

        yyss = yyss1;
        yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yytype_int16 *yyss1 = yyss;
        union yyalloc *yyptr =
          (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
                  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yylex ();
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 4:
#line 82 "parser.yy" /* yacc.c:1646  */
    {
                   if (veryverbose) {
		     cerr << "Defining " << (yyvsp[-2].symp)->name << "=" << (yyvsp[0].dval) << endl;
		   }
                   (yyvsp[-2].symp)->value = (yyvsp[0].dval); 
		}
#line 1334 "parser.cc" /* yacc.c:1646  */
    break;

  case 5:
#line 88 "parser.yy" /* yacc.c:1646  */
    { OUT << prefix << (yyvsp[0].dval) << endl; }
#line 1340 "parser.cc" /* yacc.c:1646  */
    break;

  case 6:
#line 89 "parser.yy" /* yacc.c:1646  */
    { OUT << prefix; dump_vector((yyvsp[0].dvec)); }
#line 1346 "parser.cc" /* yacc.c:1646  */
    break;

  case 7:
#line 90 "parser.yy" /* yacc.c:1646  */
    { OUT << prefix; dump_matrix((yyvsp[0].dmat)); }
#line 1352 "parser.cc" /* yacc.c:1646  */
    break;

  case 8:
#line 91 "parser.yy" /* yacc.c:1646  */
    {
	           /* We are done! */
		   exit(1);
		}
#line 1361 "parser.cc" /* yacc.c:1646  */
    break;

  case 10:
#line 98 "parser.yy" /* yacc.c:1646  */
    { (yyval.dval) = (yyvsp[-2].dval) + (yyvsp[0].dval); }
#line 1367 "parser.cc" /* yacc.c:1646  */
    break;

  case 11:
#line 99 "parser.yy" /* yacc.c:1646  */
    { (yyval.dval) = (yyvsp[-2].dval) - (yyvsp[0].dval); }
#line 1373 "parser.cc" /* yacc.c:1646  */
    break;

  case 12:
#line 100 "parser.yy" /* yacc.c:1646  */
    { (yyval.dval) = (yyvsp[-2].dval) * (yyvsp[0].dval); }
#line 1379 "parser.cc" /* yacc.c:1646  */
    break;

  case 13:
#line 102 "parser.yy" /* yacc.c:1646  */
    {	if((yyvsp[0].dval) == 0.0)
						yyerror("divide by zero");
					else
						(yyval.dval) = (yyvsp[-2].dval) / (yyvsp[0].dval);
				}
#line 1389 "parser.cc" /* yacc.c:1646  */
    break;

  case 14:
#line 107 "parser.yy" /* yacc.c:1646  */
    { (yyval.dval) = -(yyvsp[0].dval); }
#line 1395 "parser.cc" /* yacc.c:1646  */
    break;

  case 15:
#line 108 "parser.yy" /* yacc.c:1646  */
    { (yyval.dval) = (yyvsp[-1].dval); }
#line 1401 "parser.cc" /* yacc.c:1646  */
    break;

  case 17:
#line 110 "parser.yy" /* yacc.c:1646  */
    { (yyval.dval) = (yyvsp[0].ival); }
#line 1407 "parser.cc" /* yacc.c:1646  */
    break;

  case 18:
#line 111 "parser.yy" /* yacc.c:1646  */
    { (yyval.dval) = (yyvsp[0].symp)->value; }
#line 1413 "parser.cc" /* yacc.c:1646  */
    break;

  case 19:
#line 112 "parser.yy" /* yacc.c:1646  */
    {
			if((yyvsp[-3].symp)->funcptr)
				(yyval.dval) = ((yyvsp[-3].symp)->funcptr)((yyvsp[-1].dval));
			else {
				cerr << (yyvsp[-3].symp)->name << " not a function" << endl;
				(yyval.dval) = 0.0;
			}
		}
#line 1426 "parser.cc" /* yacc.c:1646  */
    break;

  case 20:
#line 120 "parser.yy" /* yacc.c:1646  */
    {
			if((yyvsp[-3].symp)->funcptr)
				(yyval.dval) = ((yyvsp[-3].symp)->funcptr)((yyvsp[-1].dval));
			else {
			        cerr << (yyvsp[-3].symp)->name << " not a function" << endl;
				(yyval.dval) = 0.0;
			}
		}
#line 1439 "parser.cc" /* yacc.c:1646  */
    break;

  case 21:
#line 130 "parser.yy" /* yacc.c:1646  */
    { (yyval.dval) = gammapolch((yyvsp[-1].ival)); }
#line 1445 "parser.cc" /* yacc.c:1646  */
    break;

  case 22:
#line 131 "parser.yy" /* yacc.c:1646  */
    { (yyval.dval) = coefxi((yyvsp[-3].ival), (yyvsp[-1].ival)); }
#line 1451 "parser.cc" /* yacc.c:1646  */
    break;

  case 23:
#line 132 "parser.yy" /* yacc.c:1646  */
    { (yyval.dval) = coefzeta((yyvsp[-3].ival), (yyvsp[-1].ival)); }
#line 1457 "parser.cc" /* yacc.c:1646  */
    break;

  case 24:
#line 133 "parser.yy" /* yacc.c:1646  */
    { (yyval.dval) = coefdelta((yyvsp[-3].ival), (yyvsp[-1].ival)); }
#line 1463 "parser.cc" /* yacc.c:1646  */
    break;

  case 25:
#line 134 "parser.yy" /* yacc.c:1646  */
    { (yyval.dval) = coefkappa((yyvsp[-3].ival), (yyvsp[-1].ival)); }
#line 1469 "parser.cc" /* yacc.c:1646  */
    break;

  case 26:
#line 135 "parser.yy" /* yacc.c:1646  */
    { (yyval.dval) = coefV((yyvsp[-3].ival), (yyvsp[-1].ival)); }
#line 1475 "parser.cc" /* yacc.c:1646  */
    break;

  case 27:
#line 138 "parser.yy" /* yacc.c:1646  */
    {
                (yyval.dvec) = new(struct vec);
		(yyval.dvec)->val = (yyvsp[0].dval);
    		(yyval.dvec)->next = 0;
            }
#line 1485 "parser.cc" /* yacc.c:1646  */
    break;

  case 28:
#line 143 "parser.yy" /* yacc.c:1646  */
    {
                struct vec *new_node = new(struct vec);
		new_node->val = (yyvsp[-2].dval);
		new_node->next = (yyvsp[0].dvec);
		(yyval.dvec) = new_node;
	    }
#line 1496 "parser.cc" /* yacc.c:1646  */
    break;

  case 29:
#line 151 "parser.yy" /* yacc.c:1646  */
    { (yyval.dvec) = 0; }
#line 1502 "parser.cc" /* yacc.c:1646  */
    break;

  case 30:
#line 152 "parser.yy" /* yacc.c:1646  */
    { (yyval.dvec) = (yyvsp[-1].dvec); }
#line 1508 "parser.cc" /* yacc.c:1646  */
    break;

  case 31:
#line 155 "parser.yy" /* yacc.c:1646  */
    {
               (yyval.dmat) = new(struct mat);
	       (yyval.dmat)->vec = (yyvsp[0].dvec);
	       (yyval.dmat)->next = 0;
         }
#line 1518 "parser.cc" /* yacc.c:1646  */
    break;

  case 32:
#line 160 "parser.yy" /* yacc.c:1646  */
    {
                struct mat *new_node = new(struct mat);
		new_node->vec = (yyvsp[-2].dvec);
		new_node->next = (yyvsp[0].dmat);
		(yyval.dmat) = new_node;
	 }
#line 1529 "parser.cc" /* yacc.c:1646  */
    break;

  case 33:
#line 168 "parser.yy" /* yacc.c:1646  */
    { (yyval.dmat) = 0; }
#line 1535 "parser.cc" /* yacc.c:1646  */
    break;

  case 34:
#line 169 "parser.yy" /* yacc.c:1646  */
    { (yyval.dmat) = (yyvsp[-1].dmat); }
#line 1541 "parser.cc" /* yacc.c:1646  */
    break;


#line 1545 "parser.cc" /* yacc.c:1646  */
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYTERROR;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;


      yydestruct ("Error: popping",
                  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  return yyresult;
}
#line 172 "parser.yy" /* yacc.c:1906  */

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

