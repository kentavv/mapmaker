#ifndef _MATHLIB_H_
#define _MATHLIB_H_

/******************************************************************************

 #    #    ##     #####  #    #  #          #    #####           #    #
 ##  ##   #  #      #    #    #  #          #    #    #          #    #
 # ## #  #    #     #    ######  #          #    #####           ######
 #    #  ######     #    #    #  #          #    #    #   ###    #    #
 #    #  #    #     #    #    #  #          #    #    #   ###    #    #
 #    #  #    #     #    #    #  ######     #    #####    ###    #    #

******************************************************************************/

/* Numeric stuff for the helpers library */

///********** THE FOLLOWING MATH FUNCTIONS AWRE KNOWN TO BE PORTABLE. **********
//
//THESE MUST BE CALLED WITH REAL (THAT IS, DOUBLE) ARGS:
//log10(); exp(); pow(); log(); floor(); ceil(); fabs(); sqrt();
//sin(); cos(); tan(); asin(); acos(); atan(); sinh(); cosh(); tanh(); atan2();
//
//INTEGER FUNCTIONS:  abs();
//
//DO NOT USE ANY OTHER MATH FUNCTIONS WHICH EXIST ON PARTICULAR SYSTEMS,
//INCLUDING:
//
//hypot(); acosh(); asinh(); atanh(); cabs(); fmod(); ldexp(); frexp();
//ldiv(); div(); labs(); modf(); bessel and gamma functions; many other
//BSD UNIX math library functions; and particular random number
//generator functions (use randnum() below).
//******************************************************************************/

/* Other useful ones... */
#define exp10(x) pow(10.0,x)
#define eq(x, y)  (rabs((x)-(y))<VERY_SMALL)


/***** SUPPORT FOR SOME GENERALLY USEFUL MATH FUNCTIONS *****/

real sq(real r);

real rmaxf(real r, real s);

real rminf(real r, real s);

long lpow2(int i);

int ipow(int x, int i);

int ichoose(int n, int k);

int imaxf(int r, int s);

int iminf(int r, int s);

int icomp(const void *x, const void *y);

int lcomp(const void *x, const void *y);

int rcomp(const void *x, const void *y);

int scomp(const void *x, const void *y);

int inv_icomp(const void *x, const void *y);

int inv_rcomp(const void *x, const void *y);

real rmean(real *data, int length);

void mat_invert(real **m, int size, real **m_inverse);

void eliminate(int row, int col, int row_to_sub, real **m, int n_cols);

void array_times_matrix(const real *a, real **b, int rows, int columns, real *c);



//real sq(); /* arg: real x; returns x squared */

/* Portable functions to extract the real and fractional portions of real 
   numbers. For example:
       rint(4.2)=4.0; frac(4.2)=0.2; 
       rint(-4.2)=-4.0; frac(-4.2)=-0.2;
   Thus, it is always true that x==rint(x)+frac(x). Note that this is not the
   BSD Unix rint() function! */

/* NOTE: all of the args here are ints, even if a long is returned */
//long lpow();  /* arg: int x,i; returns x to the i, for i<=31*/
//int  ipow();  /* arg: int x,i; returns x to the i, for small x,i */
//long lpow2(); /* arg: int i;   returns 2 to the i, for i<=31*/
//int  ipow2(); /* arg: int i;   returns 2 to the i, for i<=15*/


/* Remember, its a no-no to change any of these defs! VAX is the smallest! */
#define REAL(x)    ((real)(x))
#define MAXREAL    ((real) 1e30)
#define MINREAL    ((real) 1e-30)
#define VERY_BIG   MAXREAL
#define VERY_SMALL MINREAL
#define OBSCURE_REAL  (-1.23456789e+31)
#define OBSCURE_REAL2 (-1.23456789e+30)

#define MAXINT   32767
#define INTBITS  15
#define LONGBITS 31


/* The min() and max() macros only really work right if the arguments
   are very simple (eg variables or constants). However, they should work on 
   any atomic type (ints, shorts, longs, doubles, chars, and even identical 
   pointer types) If the arguments are functions or expressions, they may fail 
   [for ex: max(10,i++) does the wrong thing], or they may take too long 
   [for ex: max(log(x),log(y)) takes longer that rmaxf(log(x),log(y)); */

#define max(a, b) ((a)>(b) ? (a) : (b))
#define min(a, b) ((a)<(b) ? (a) : (b))

/* The following functions exist so that rmaxf(f(x),f(y)) works, where f() 
   side-effects something, or where f() takes lots of time. Hard code
   a conditional expression instead if speed is really critical... */

//real rmaxf();   /* args: real r,s; returns largest of r and s */
//real rminf();   /* args: real r,s; returns smallest of r and s */
//int  imaxf();   /* args: int i,j;  returns largest of i and j */
//int  iminf();   /* args: int i,j;  returns smallest of i and j */
//long lmaxf();   /* args: long i,j; returns largest of i and j */
//long lminf();   /* args: long i,j; returns smallest of i and j */


//real randnum(); /* no args; number from 0...(1-epsilon) */
#define coin_flip() (randnum()>0.5)

//void do_seedrand();
#define seedrand(x) do_seedrand((long)(x))
/* effectively: void seedrand(x); long x; x may be one of: */
#define RANDOM    0l
#define NONRANDOM 1l

#define math_init() { seedrand(NONRANDOM); } /* called by lib_init() */

/***** OTHER USEFUL OPERATIONS FOR REAL NUMBER ARRAYS *****/

//real rmean();	/* real *data; int array_len; returns mean */
//real rmaxin();	/* real *data; int array_len; returns largest */
//real rminin();	/* real *data; int array_len; returns smallest */
//real rmedian();	/* real *data; int array_len; returns median */
//real rmiddle();	/* real *data; int array_len; returns middle entry */
//void rcopy(); 	/* real *to, *from; int array_len; copies array */

///* This one needs work */
//bool rhistogram();	/* real *data; int array_len, min_num_buckets;... */
//			/* real foo, bar; (unused) - returns T or F */
//
///***** AND MATRICIES... *****/
//
//void mat_invert(); /* args: real **m; int size; real **m_inverse; */
///* Invert square matrix m by Gauss' method, side-effecting the
//   (already allocated!) matrix m_inverse. m_inverse should be 2*size
//   columns (2nd index) by size rows (first index) - its left square will be
//   left with the result (ignore the right side). */
//
//void mat_mult(); /* args: real **m, **m2; int size; real **m_times_m2; */
///* Multiply square matricies m and m2, side-effecting the (already allocated!)
//   matrix m_times_m2, which should be the same size as m and m2 */
//
//void array_times_matrix(); /* real *a, **m; int rows, columns; real *result; */
///* Multiply array a (length=rows) times matrix b (indicies=[row][column]),
//   side effecting the (already allocated!) array c (length=columns). */


#endif
