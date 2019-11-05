#define NSYMS 200	/* maximum number of symbols */

struct symtab {
	char *name;
	double (*funcptr)(double);
	double value;
};

struct symtab *symlook(char *);

struct vec 
{
   double val;
   struct vec *next;
};

struct mat
{
   struct vec *vec;
   struct mat *next;
};
