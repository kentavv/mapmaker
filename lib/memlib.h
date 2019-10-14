#ifndef _MEMLIB_H_
#define _MEMLIB_H_

/******************************************************************************

 #    #  ######  #    #  #          #    #####           #    #
 ##  ##  #       ##  ##  #          #    #    #          #    #
 # ## #  #####   # ## #  #          #    #####           ######
 #    #  #       #    #  #          #    #    #   ###    #    #
 #    #  #       #    #  #          #    #    #   ###    #    #
 #    #  ######  #    #  ######     #    #####    ###    #    #

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

/***************************************************************************
The following are good specs for writing memory allocating/freeing routines:

     FILL IN LATER

Do not use any of the C-library mem... functions, including memcmp(), 
memcpy() and so on. It is also strongly recomended that you avoid 
directly using malloc(), calloc(), free() or any other C memory allocation
routines, for many reasons. Use those defined in here instead. 
***************************************************************************/

typedef real REAL2[2];
typedef real REAL2x2[2][2];
typedef real REAL4[4];

extern bool verbose_mem;    /* may be set by user code */
extern int yy, zz;        /* INTERNAL USE ONLY! */

/***************************************************************************
   Syntax for allocating and freeing things:	
   All array dimensions must be ints (not unsigned, or size_t, or...),
   and arrays are limited to 32K bytes for portability! (Matrixes, being 
   arrays of pointers to arrays, each being up to 32K, thus may be huge).
   We also limit each xalloc() call to 64K bytes, and thus the effective 
   portable maximum for most array dimensions is 6K-8K (real), 16K 
   (long or pointer), 32K (short or char).
   
	typedef struct {...} foo;
	foo *z;
	single(z, foo);
	unsingle(z, foo);
	
	int *x; 
	array(x, 10, int);
	unarray(x, int);

	double **y;
	matrix(y, 10, 20, double);	* i,j dimensions (eg y[i][j])
	unmatrix(y, 10, int); 		* The i dimension must be given

	foo **z;		* parray: each z[i] is a pointer to a foo
	parray(z, 10, foo);	* parray allocates the pointer array, and
	unparray(z, foo);	* the FOO's	

   Note that this last construct is identical in function to:

	array(z, 10, (foo*));		
	for (i=0; i<10; i++) single(z[i], foo);

These routines have the following nice behaviors. If a allocation error is 
encountered mid-execution, then the entire array or matrix is free()ed and
the NOMEMORY message is sent. The free routines all ignore all NULL pointers,
including those pointing to the pointer array of a matrix, and those in the
pointer array. 

It is strongly recomended that all allocating and freeing routines be written 
to similar specs.
***************************************************************************/

/***** allocate *****/

#define xL(i) ((size_t)(i))
#define xS(c) sizeof(c)
#define REALLY_1 ((size_t)(-32768))

#define single(var, cell) \
{ var=NULL; if ((var=(cell*)xalloc(REALLY_1,xS(cell)))==NULL) send(NOMEMORY); }

#define array(var, i, cell) \
{ var=NULL; if ((var=(cell*)xalloc(xL(i),xS(cell)))==NULL) send(NOMEMORY); }

#define matrix(var, i, j, cell) \
{ var=NULL; if ((var=(cell**)xalloc(xL(i),xS(cell*)))==NULL) send(NOMEMORY); \
    for (zz=0; zz<i; zz++) { \
        if ((var[zz]=(cell*)xalloc(xL(j),xS(cell)))==NULL) { \
        for (yy=0; yy<zz; yy++) unarray(var[yy], cell);  \
        unarray(var, (cell*)); send(NOMEMORY); \
    } \
    } \
}

#define parray(var, num, cell) matrix(var,num,REALLY_1,cell)

/***** free *****/

#define unarray(ptr, cell) {if(ptr!=NULL) free(ptr); ptr=NULL;}

#define unsingle(ptr, cell) unarray(ptr,cell)

#define unmatrix(ptr, i, cell) \
{ if (ptr!=NULL) for(zz=0; zz<i; zz++) unarray(ptr[zz],cell); \
  unarray(ptr,(cell*)); }

#define unparray(ptr, num, cell) unmatrix(ptr,num,cell)

/***************************************************************************
Other useful allocators go here...
***************************************************************************/

void *xalloc(size_t num, size_t cell_sizeof); /* INTERNAL USE ONLY! */

char ***alloc_char_3d_matrix(int i, int j, int k);

void free_char_3d_matrix(char ***p, int i, int j);

void mem_init(void);

#endif