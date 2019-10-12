#ifndef _MISCLIB_H_
#define _MISCLIB_H_

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

#include "syscode.h"

/***** GETTING THE TIME (code is in syscode.c) *****/
/* DO NOT use the time() system function alone, it's not portable! */

//char *time_string(); /* no args; returns NULL if it fails */
//real usertime(); /* args: bool do_reset; */
///* returns time in seconds, or -1.0 if it fails: do_reset makes it work like
//   a stopwatch (it starts counting when lib_init is executed) */

/***** SUBPROCESSES (code is in syscode.c) *****/

///* These both return FALSE on failure, TRUE if OK */
//bool shell_command(); /* arg: char *command; defined in system.h */
//bool subshell();      /* no args. returns T/F */


/***** GET/SET DIRECTORY (code is in syscode.c) *****/

///* These both also return FALSE on failure, TRUE if OK */
//bool get_directory();     /* args: char *str; str>=PATH_LENGTH+1 chars */
//bool change_directory();  /* args: char *str; str is directory name */


/***** SORT OPERATIONS FOR SIMPLE ARRAYS (code is in mathlib.c) *****/

/* effective declarations: 
   void isort(); args: int  *data; int array_len; 
   void lsort(); args: long *data; int array_len; 
   void rsort(); args: real *data; int array_len; 
   void ssort(); args: real *data; int array_len; STRINGS, NOT SHORTS!
   void inv_isort(); args: int *data; int array_len; ascending order 
   void inv_rsort(); well, figure it out!       */

//int icomp(), lcomp(), scomp(), rcomp(), inv_icomp(), inv_rcomp();

#define isort(p,n) qsort(p,n,sizeof(int),icomp);
#define lsort(p,n) qsort(p,n,sizeof(long),lcomp);
#define rsort(p,n) qsort(p,n,sizeof(real),rcomp);
#define ssort(p,n) qsort(p,n,sizeof(char*),scomp);
#define inv_isort(p,n) qsort(p,n,sizeof(int),inv_icomp);
#define inv_rsort(p,n) qsort(p,n,sizeof(real),inv_rcomp);
 
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
 
#include "mathlib.h"

#endif
