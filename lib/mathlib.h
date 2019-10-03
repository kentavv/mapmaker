/******************************************************************************

 #    #    ##     #####  #    #  #          #    #####           #    #
 ##  ##   #  #      #    #    #  #          #    #    #          #    #
 # ## #  #    #     #    ######  #          #    #####           ######
 #    #  ######     #    #    #  #          #    #    #   ###    #    #
 #    #  #    #     #    #    #  #          #    #    #   ###    #    #
 #    #  #    #     #    #    #  ######     #    #####    ###    #    #

******************************************************************************/

/* Numeric stuff for the helpers library */

/********** THE FOLLOWING MATH FUNCTIONS AWRE KNOWN TO BE PORTABLE. **********

THESE MUST BE CALLED WITH REAL (THAT IS, DOUBLE) ARGS:
log10(); exp(); pow(); log(); floor(); ceil(); fabs(); sqrt(); 
sin(); cos(); tan(); asin(); acos(); atan(); sinh(); cosh(); tanh(); atan2();

INTEGER FUNCTIONS:  abs(); 

DO NOT USE ANY OTHER MATH FUNCTIONS WHICH EXIST ON PARTICULAR SYSTEMS, 
INCLUDING:

hypot(); acosh(); asinh(); atanh(); cabs(); fmod(); ldexp(); frexp();
ldiv(); div(); labs(); modf(); bessel and gamma functions; many other
BSD UNIX math library functions; and particular random number
generator functions (use randnum() below).
******************************************************************************/

/* Other useful ones... */
#define exp10(x) pow(10.0,x)
#define eq(x,y)  (rabs((x)-(y))<VERY_SMALL)
real sq(); /* arg: real x; returns x squared */

/* Portable functions to extract the real and fractional portions of real 
   numbers. For example:
       rint(4.2)=4.0; frac(4.2)=0.2; 
       rint(-4.2)=-4.0; frac(-4.2)=-0.2;
   Thus, it is always true that x==rint(x)+frac(x). Note that this is not the
   BSD Unix rint() function! */

#define rint(x)	   (x>=0.0 ? floor(x) : ceil(x)) 
#define frac(x)    (x>=0.0 ? (x-floor(x)) : (x-ceil(x)))
#define sign(x)    (x>=0.0 ? 1.0 : -1.0)

/* NOTE: all of the args here are ints, even if a long is returned */
long lpow();  /* arg: int x,i; returns x to the i, for i<=31*/
int  ipow();  /* arg: int x,i; returns x to the i, for small x,i */
long lpow2(); /* arg: int i;   returns 2 to the i, for i<=31*/
int  ipow2(); /* arg: int i;   returns 2 to the i, for i<=15*/


/* Remember, its a no-no to change any of these defs! VAX is the smallest! */
#define REAL(x)    ((real)(x))
#define MAXREAL    ((real) 1e30)
#define MINREAL    ((real) 1e-30)
#define VERY_BIG   MAXREAL
#define VERY_SMALL MINREAL
#define pi (3.14159265358978)
#define OBSCURE_REAL  (-1.23456789e+31)
#define OBSCURE_REAL2 (-1.23456789e+30)

#define INT(x)   ((int)(x))
#define MAXINT   32767
#define MININT   (-MAXINT)
#define NONINT   (-32768)
#define MAXLONG  2147483647l
#define MINLONG  (-MAXLONG)
#define INTBITS  15
#define LONGBITS 31
#define enbool(i) (((int)(i))!=0 ? TRUE : FALSE)
#define unbool(x) ((x) ? 1 : 0)


/* The min() and max() macros only really work right if the arguments
   are very simple (eg variables or constants). However, they should work on 
   any atomic type (ints, shorts, longs, doubles, chars, and even identical 
   pointer types) If the arguments are functions or expressions, they may fail 
   [for ex: max(10,i++) does the wrong thing], or they may take too long 
   [for ex: max(log(x),log(y)) takes longer that rmaxf(log(x),log(y)); */

#define max(a,b) ((a)>(b) ? (a) : (b))
#define min(a,b) ((a)<(b) ? (a) : (b))

/* The following functions exist so that rmaxf(f(x),f(y)) works, where f() 
   side-effects something, or where f() takes lots of time. Hard code
   a conditional expression instead if speed is really critical... */

real rmaxf();   /* args: real r,s; returns largest of r and s */
real rminf();   /* args: real r,s; returns smallest of r and s */
int  imaxf();   /* args: int i,j;  returns largest of i and j */
int  iminf();   /* args: int i,j;  returns smallest of i and j */
long lmaxf();   /* args: long i,j; returns largest of i and j */
long lminf();   /* args: long i,j; returns smallest of i and j */


real randnum(); /* no args; number from 0...(1-epsilon) */
#define coin_flip() (randnum()>0.5)
#define biased_flip(p_heads) (randnum()<(p_heads))

void do_seedrand();
#define seedrand(x) do_seedrand((long)(x))
/* effectively: void seedrand(x); long x; x may be one of: */
#define RANDOM    0l
#define NONRANDOM 1l

#define math_init() { seedrand(NONRANDOM); } /* called by lib_init() */


