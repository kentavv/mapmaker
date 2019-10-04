/******************************************************************************

  ####    #####  #####   #          #    #####           #    #
 #          #    #    #  #          #    #    #          #    #
  ####      #    #    #  #          #    #####           ######
      #     #    #####   #          #    #    #   ###    #    #
 #    #     #    #   #   #          #    #    #   ###    #    #
  ####      #    #    #  ######     #    #####    ###    #    #

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

/***** String and parsing functions for the helpers library *****/

/****************************************************************************
Unfortunately, much of the C-library of string functions is not very
portable: Most functions are in general the same, although they differ
in a number of minor anoying ways. Of particular annayance are strlen()
and the strn... functions (strncmp(), strncpy(), and strncat()).

The functions strchr(), strrchr(), strspn(), strcspn(), strbrk(), and
strtok() do not exist in all C implementations. Also, the string <->
number conversion routines, including atof(), strtod(), strto(),
atoi(), atol() and strtoul() may have compatibility problems. Some of
the helper functions provide similar functionalities while preserving
compatibility: use them instead!

The functions toupper() etc. mentioned in K&R seem to be portable. 
However, they are macros in some implementations! Thus, don't try
toupper(ptr++), as ptr may get incremented more than once! Note that
contrary to what K&R says, it is entirely reasonable to assume ASCII 
codes are in use. 

The only library functions known to be portable are strcat(),
strcpy(), and strcmp(). Otherwise, you should use the functions
provided here. Note that strcat() and strcpy() are inherently 
dangerous, and they do not provide any bounds checking, and may 
result in a crahing program. nstrcat() etc. are much preferred!
****************************************************************************/

/* Replacements for the C-Library... */

/* These will return TRUE if they don't have to truncate, FALSE otherwise. */
/* THIS IS NOT IMPLEMENTED YET! */
void nstrcpy();	  /* args: char *to, *from; int num; copy <=num chars */
void nstrcat();	  /* args: char *to, *from; int num; append <=num chars */
void maxstrcat(); /* args: char *to, *from; int max; length kept <=max chars */

/* To avoid the ANSI size_t idiocy... */

#define nstrcmp(s1,s2,max_chars) strncmp(s1,s2,((int)max_chars))
#define len(str) ((int) strlen(str)) 

/* Other useful stuff... */

char *mkstrcpy(); /* args: char *s; returns an allocated copy */ 
char *ptr_to();   /* args: constant character string - return a ptr to it */
void strdel();    /* args: char *s; int num; deletes num chars at s */ 
#define endof(str) ((str)+len(str))

int strfinder();  /* args: char *str, c; get index of c in str or NOT_FOUND */ 
#define NOT_FOUND (-1) 
#define strin(str,chr) (strfinder(str,chr)!=NOT_FOUND)

void strins();    /* args: char *a, *b; insert string b at a */
void nstrins();	  /* args: char *a, *b; int num; insert at most num chars */
void maxstrins(); /* args: char *a, *b; int max; length kept <= max */

bool nullstr();      /* args: char *s; TRUE if s==NULL or is all whitespace */
extern char *null_string; /* set to "" */
#define streq(s1,s2)    (!strcmp(s1,s2)) 
#define nstreq(s1,s2,n) (!nstrcmp(s1,s2,n)) 

extern char Cw, Ct;
#define white(chr) ((Cw=(chr))==' ' || Cw=='\t' || Cw=='\n')
#define trash(chr) (((Ct=(chr))!='\0') && (Ct<' ' || Ct>'~') && !white(Ct))

bool nmatches(); /* args: char *s, *t; int num; */ 
/* returns TRUE if chars in the 1st token in s match those in the token in 
   "template" t, and if the token in s is at least num chars long. For 
   example: "land" matches the template "lander", but not "lampshade". 
   NOTE: s and t should be despace()ed and filter()ed, but not necessarily 
   lowercase()ed (the match is case INSENSITIVE anyway). */
#define matches(s,t) nmatches(s,t,1)

int xstreq(); 	/* currently broken? */
#define istrlen len /* THIS IS AN OBSOLETE NAME- DON'T USE IT */

/****************************************************************************
Various string crunching routines: despace() changes all globs of
whitespace to single spaces and does away with leading and trailing
whitespace entirely. filter() removes all non-printing ASCII
characters, and lowercase() converts uppercase letters to lowercase().
crunch() invokes despace, filter and lowercase.  filter_nonspaces() is
like filter, except that in addition tabs and newlines are converted
to spaces. uppercase() is the obvious opposite of lowercase().
truncstr() limits the length of a string to some number of characters
(not including the trailing '\0') pad_to_len() adds spaces to the end
of the string until it is a particular length, while append_spaces()
simply adds the requested number of spaces to the string. All
side-effect their argument str, and return a pointer to it for yucks.
****************************************************************************/

char *despace(); 	   /* args: char *str; side-effected */ 
char *lowercase();	   /* args: char *str; side-effected */ 
char *uppercase();	   /* args: char *str; side-effected */
char *_filter();		   /* args: char *str; side-effected */ 
char *filter_nonspaces();  /* args: char *str; side-effected */ 
char *crunch(); 	   /* args: char *str; despace(_filter(lowercase())) */

char *truncstr();	/* args: char *str; int max_chars; str side-effected */
char *pad_to_len();     /* args: char *str; int max_chars; adds spaces */
char *append_spaces();  /* args: char *str; int num_spaces; also adds spaces */

/****************************************************************************
Each of the token-parsing functions work as follows: 

bool itoken(), ltoken(), rtoken();
args: char **p_str; <value type> default_value, *result; 

bool stoken();
args: char **p_str; char *default_value, *result; 

bool nstoken(), maxstoken();
args: char **p_str; char *default_value,*result; int num_chars;

bool stokenof();
args: char **p_str; char *default_value, *result; char *parsable_chars; 

bool nstokenof(), maxstokenof();
args: char **p_str; char *default_value, *result; int num_chars; 
      char *parsable_chars;

If one of these succeeds: *p_str points to the delimiting character which
follows the token (which may be '\0'), *result is set, and TRUE is returned.

If no token is avail: *p_str points to the '\0' at the end of the
string. In this case, if a default is given then TRUE is returned and
*result is set.  Otherwise, FALSE is returned, and *result is
undefined.

If the token is bad: *p_str points to the beginning of the token (so
that nullstr(*p_str)==FALSE, see bad_token() below). If a default is
available, then *result is set to it, otherwise *result is undefined.
FALSE is always returned.

Note that for itoken(), rtoken(), ltoken(), stoken(), and stokenof(),
the length of a token is limited to TOKLEN chars. (Thus, stoken() etc.
should be passed a pointer to a string of at least TOKLEN+1 chars to
hold the result.) Longer tokens are truncated (which always makes
numbers 'bad').

For nstoken(), the user may specify the length of the result string to
use instead of TOKLEN. For maxstoken() the length of the token must be
<= num_chars characters, otherwise the result string is truncated,
*p_str is left untouched, and FALSE is returned.

stokenof(), nstokenof(), and maxstokenof() specify the legal
characters which may comprise the token. Other characters (excluding
the self_delimiting, described below) cause FALSE to be returned, and
no token is still parsed from the string. ANYCHAR (really NULL,
defined below for parse_char()) may be used to indicate that for
parsable_chars to indicate that any character is OK.
****************************************************************************/

#define TOKLEN 40 
int itoken();	 /* int token */ 
int ltoken();	 /* long int token */ 
int rtoken();	 /* real token */

int stok(); /* INTERNAL USE ONLY! */ 
int stoken(); /* args: p,def,val; does stok(p,def,val,TOKLEN,TRUE,NULL) */
#define nstoken(p_str,def,val,num) stok(p_str,def,val,num,TRUE,NULL) 
#define maxstoken(p_str,def,val,num) stok(p_str,def,val,num,FALSE,NULL)
#define stokenof(p_str,def,val,chrs) stok(p_str,def,val,TOKLEN,TRUE,chrs) 
#define nstokenof(p_str,def,val,num,chrs)   stok(p_str,def,val,num,TRUE,chrs)
#define maxstokenof(p_str,def,val,num,chrs) stok(p_str,def,val,num,FALSE,chrs)

/* Possible default values */ 
#define sREQUIRED NULL 
#define iREQUIRED (-32768) 
#define lREQUIRED -1073741823L 
#define rREQUIRED ((real)-1.2345e31)

/* To decipher FALSE responses of the token parsers */ 
#define no_token(p_str)  (**p_str=='\0') 
#define bad_token(p_str) (**p_str!='\0')

/*** Usually, tokens are separated by whitespace. Certain characters
however, which are listed in the self_delimiting string, will always
be parsed as separate tokens, whether surrounded by whitespace or not.
For example: the string "53 - (14+2)" with self_delimiting equal to
"()-+" will parse into tokens "53","-","(","14","+","2", and ")". ***/
extern char *self_delimiting;

/*** Other useful stuff for parsing ***/
int count_tokens();  /* args: char *str; */
bool is_a_token();   /* args: char *str; must be despaced, or from stoken() */
bool split_string(); /* args: char *str, **rest, divider; rest side-effected */


/****************************************************************************
The range functions work as follows:

   args: <value type> *value, low, high; returns bool;

   If *value is in the range [low,high] (inclusive), then TRUE is
returned.  Otherwise, *value is set to the apprpriate limit (low or
high) and FALSE is returned.
****************************************************************************/

int irange();	 /* integer range */ 
int lrange(); 	 /* long int range */ 
int rrange();	 /* real range */

/****************************************************************************
Parse_char parsing is similar to (and compatible with) the token
parsing routines shown above.  

bool parse_char(char **p_str,*parsable_chars; int skip_whitespace; char *c;)

First, if skip_white is TRUE, *p_str is incremented until a non-
whitespace char or the '\0' is encountered.  If **p_str is '\0', FALSE
is returned and *c is set to '\0'.  If **p_str is in the
parsable_chars string, or if parsable_chars is NULL, then *c=**p_str
and *p_str is incremented.  If **p_str is not in the parsable_chars
string, then c=**p_str, *p_str is NOT incremented, and FALSE is
returned.

parse_whitespace() moves the ptr along until a non-whitespace character
is encountered.
****************************************************************************/

bool parse_char();	  /* args shown above */
void parse_whitespace();  /* args: char **p_str; *p_str is side-effected */

/* Arguments to parse_char() */ 
#define ANYCHAR    NULL 
#define SKIPWHITE  TRUE 
#define NOSKIP     FALSE

/* To decipher FALSE responses of parse_char */ 
#define no_char(p_str) (**p_str=='\0') 
#define bad_char(p_str) (**p_str!='\0' && !white(**p_str)) 
#define white_char(p_str) (white(**p_str))


/***********************************************************************
The global pool of strings for printing things into. These strings are
"allocated" from a reusable global pool of strings, which when
exhausted begins reusing strings previously allocated. Thus, strings
gotten using get_temp_str() should be considered VERY temporary
storage: use mkstrcpy() etc. to make more permanent storage. The same
caution applies to any routines which return strings "allocated" by
this function, including pr(), prn(), prd(), and others. Be sure to
doccument this in your functions which use these or which call
get_temp_str() directly.
***********************************************************************/

char *get_temp_string(); /* returns the next available string for bashing */

#define NUM_TEMP_STRINGS 50 
#define TEMP_STRING_LEN 500

/************************************************************************
rs() etc.: The funky real number printing routines...

rs(), meaning "real to string", takes a format number (like one given
to sprintf) and a data number (both reals) and returns a string
containing a human readable form of the number. Unlike sprintf(), if
the number can't fit properly, decimal places are thrown away. If it
still can't fit, the string is filled with asterisks. Thus, unlike
sprintf()ed strings, the returned string will ALWAYS have the
specified length. The format number must be of the form n.m, where n
and m are SINGLE DIGITS, m<=n-2, and m.n>0.0. n specifies the printed
string length, and m specifies the desired number of decimal places.
If the number to be printed is negative, then the minus sign will take
one space of the string's length. Note that decimal places which are
truncated are simply cut out, not rounded out!

rsn(), with 'n' for "negative", is like rs() except that if the number
is positive, a leading space is printed where the minus sign will
appear for negative numbers. This way columnar output will have have
the leading digits (rather than leading digit OR minus sign) line up,
and positive and negative numbers will both be printed to the same
precision.

rsd(), with 'd' for "decimal", is also essentially the same as rs() except
that it adds spaces to the front rather than rear for shorter numbers,
in order to make columns of numbers line up on their decimal points.

Each of these routines return a string "allocated" (appropriated) using
get_temp_string(). Thus, heed the warnings above. 
************************************************************************/

char *rs();   /* args: real format, num_to_print; */
char *rsn();  /* args: real format, num_to_print; */
char *rsd();  /* args: real format, num_to_print; */

/* other output formating stuff... */ 
char *binary(); /* args: int num_to_print, num_bits; char *str; */

/* macro char *maynl();  args: int chars; 
   returns a ptr to a string containing only the \n character if the number of 
   chars is too many to fit on a screen line, and otherwise returns a null 
   (e.g. zero length) string. For example:

   sprintf(str,"successfully loaded file: %s%s\n",maynl(len(name)+26),name); 

   maynls(str,chars) is the same as maynl(len(str)+chars) */

#define maynl(chars) ((chars)>LINE ? "\n" : "")
#define maynls(str,chars) ((len(str)+chars)>LINE ? "\n" : "")

/* macro char *maybe_s(num); returns a string containing "s" if num!=1,
   or a string of "" if num==1. Useful for making words maybe plural. */
#define maybe_s(n) ((n)!=1 ? "s" : "")
/* macro char *maybe_sp(num); returns a string containing "" if num!=1,
   or a string of " " if num==1. Useful with maybe_s for lining things up */
#define maybe_sp(n) ((n)!=1 ? "" : " ")

void str_init();







