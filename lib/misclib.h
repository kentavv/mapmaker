/******************************************************************************

 #    #     #     ####    ####   #          #    #####           #    #
 ##  ##     #    #       #    #  #          #    #    #          #    #
 # ## #     #     ####   #       #          #    #####           ######
 #    #     #         #  #       #          #    #    #   ###    #    #
 #    #     #    #    #  #    #  #          #    #    #   ###    #    #
 #    #     #     ####    ####   ######     #    #####    ###    #    #

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

/***** GETTING THE TIME (code is in syscode.c) *****/
/* DO NOT use the time() system function alone, it's not portable! */

char *time_string(); /* no args; returns NULL if it fails */
real usertime(); /* args: bool do_reset; */
/* returns time in seconds, or -1.0 if it fails: do_reset makes it work like
   a stopwatch (it starts counting when lib_init is executed) */

/***** SUBPROCESSES (code is in syscode.c) *****/

/* These both return FALSE on failure, TRUE if OK */
bool shell_command(); /* arg: char *command; defined in system.h */
bool subshell();      /* no args. returns T/F */


/***** GET/SET DIRECTORY (code is in syscode.c) *****/

/* These both also return FALSE on failure, TRUE if OK */
bool get_directory();     /* args: char *str; str>=PATH_LENGTH+1 chars */
bool change_directory();  /* args: char *str; str is directory name */


/***** SORT OPERATIONS FOR SIMPLE ARRAYS (code is in mathlib.c) *****/

/* effective declarations: 
   void isort(); args: int  *data; int array_len; 
   void lsort(); args: long *data; int array_len; 
   void rsort(); args: real *data; int array_len; 
   void ssort(); args: real *data; int array_len; STRINGS, NOT SHORTS!
   void inv_isort(); args: int *data; int array_len; ascending order 
   void inv_rsort(); well, figure it out!       */


#define isort(p,n) qsort(QSORT_CAST(p),(QSORT_LENGTH)n,sizeof(int),icomp); 
#define lsort(p,n) qsort(QSORT_CAST(p),(QSORT_LENGTH)n,sizeof(long),lcomp); 
#define rsort(p,n) qsort(QSORT_CAST(p),(QSORT_LENGTH)n,sizeof(real),rcomp); 
#define ssort(p,n) qsort(QSORT_CAST(p),(QSORT_LENGTH)n,sizeof(char*),scomp); 
#define inv_isort(p,n) qsort(QSORT_CAST(p),(QSORT_LENGTH)n,sizeof(int),inv_icomp);
#define inv_rsort(p,n) qsort(QSORT_CAST(p),(QSORT_LENGTH)n,sizeof(real),inv_rcomp);
 
int icomp(), lcomp(), scomp(), rcomp(), inv_icomp(), inv_rcomp();

/* For arrays of pointers to things (eg: structs) use:
   void psort(); args: <type> *data; int array_len; <type>; int comp_func();

   Comp_func() will be passed two ptrs to ptrs to <type>s: you should
   derefernce these carefully (as shown below) for portability. For
   Example, here we sort an array of ptrs to FOOs on the field "key":
   Note that FOOs could have variable sizes (only the pointers to them
   are getting shuffled).

   \* Make and sort an array of pointers to FOOs *\ 
   typedef struct { int key; ... } FOO; 
   FOO **data;                             
   array(data, length, FOO*);
   for (i=0; i<length; i++) single(data[i], FOO);  \* or use parray() *\
   ...
   int compare_foos(p1,p2)
   PSORT_COMPARE_TYPE(FOO) **p1, **p2;     \* ptrs to ptrs to FOOs *\  
   {
       int key1, key2;
       key1=(*(FOO**)p1)->key;  
       key2=(*(FOO**)p2)->key;

       if (key1 < key2) return(-1); 
       else if (key1==key2) return(0); 
       else return(1);
   }
   ...
   psort(data, 100, FOO, compare_foos);

You may want to look in system.h to see how these work in order to write 
your own sorting routines. */

#define psort(p,n,t,f) qsort(QSORT_CAST(p),(QSORT_LENGTH)n,sizeof(t*),f)
#define PSORT_COMPARE_TYPE(t) QSORT_COMPARE_PTR_TO(t)
 

/***** OTHER USEFUL OPERATIONS FOR REAL NUMBER ARRAYS *****/

real rmean();	/* real *data; int array_len; returns mean */
real rmaxin();	/* real *data; int array_len; returns largest */
real rminin();	/* real *data; int array_len; returns smallest */
real rmedian();	/* real *data; int array_len; returns median */	
real rmiddle();	/* real *data; int array_len; returns middle entry */
void rcopy(); 	/* real *to, *from; int array_len; copies array */

/* This one needs work */
bool rhistogram();	/* real *data; int array_len, min_num_buckets;... */
			/* real foo, bar; (unused) - returns T or F */

/***** AND MATRICIES... *****/ 

void mat_invert(); /* args: real **m; int size; real **m_inverse; */
/* Invert square matrix m by Gauss' method, side-effecting the
   (already allocated!) matrix m_inverse. m_inverse should be 2*size
   columns (2nd index) by size rows (first index) - its left square will be
   left with the result (ignore the right side). */

void mat_mult(); /* args: real **m, **m2; int size; real **m_times_m2; */
/* Multiply square matricies m and m2, side-effecting the (already allocated!)
   matrix m_times_m2, which should be the same size as m and m2 */

void array_times_matrix(); /* real *a, **m; int rows, columns; real *result; */
/* Multiply array a (length=rows) times matrix b (indicies=[row][column]),
   side effecting the (already allocated!) array c (length=columns). */


/***** SUPPORT FOR SOME GENERALLY USEFUL MATH FUNCTIONS *****/

typedef struct {
    int entries;
    real start, increment, mean;
    real *deviation, *prob; /* [entries], prob is cumulative */
} DISTRIBUTION;

DISTRIBUTION *make_distribution(); 
/* args: real increment,limit,mean,(*function)();
   Makes a probability distribution for d= mean-limit to mean+limit by 
   increment, where function(d) is the probability of deviation d from mean. */

DISTRIBUTION *make_normal_distribution();
/* args: real mu, sigma, increment, limit; 
   Makes a normal probability distribution with sigma about mean mu, where
   increment and limit are both expressed as fractions of sigma (eg std
   deviations, unlike the case in make_distribution()) */

real pick_from_distribution();  /* args: DISTRIBUTION *dist; real *prob;
  Return a randomly chosen deviation d from the distribution, and optionally
  set *prob (if non-NULL) to d's CUMULATIVE probability. */

