/******************************************************************************

 #    #  ######  #    #  #          #    #####            ####
 ##  ##  #       ##  ##  #          #    #    #          #    #
 # ## #  #####   # ## #  #          #    #####           #
 #    #  #       #    #  #          #    #    #   ###    #
 #    #  #       #    #  #          #    #    #   ###    #    #
 #    #  ######  #    #  ######     #    #####    ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

//#define INC_IO
//#define INC_STR
//#define INC_MSG
//#define INC_MEM
//#define INC_HELP_DEFS
#include "system.h"

int verbose_mem;       /* may be set by user */
long n_alloced;

int yy, zz;

#define ALLOC_MSG \
    "alloc:\t%ldx%ld\tbytes @ %8lxH. total bytes alloced=%ld\n"
#define ALLOC_FAIL_MSG \
    "FAILED:\t%ldx%ld\tbytes              total bytes alloced=%ld\n"

void *xalloc(size_t num, size_t cell_sizeof) {
    void *p;

    if (num == REALLY_1) num = 1;
    else if (num <= 0) send(CRASH);
    else num++; /* paranoia */

//	long chars=num*cell_sizeof;
    /* Note that this test will be inadequate if word-aligning
       causes things to take up more space than their sizeof()! */
    /* if (chars > 65472L) send(CRASH);  64K-64 bytes max KLUDGE */

    if ((p = calloc(num, cell_sizeof)) == NULL) {
        if (TRUE || verbose_mem) { /* KLUDGE */
            sprintf(ps_, ALLOC_FAIL_MSG, num, cell_sizeof, n_alloced);
            print(ps_);
        }
        NOMEMORY_num_cells = (int) num;
        NOMEMORY_cell_size = (int) cell_sizeof;
        return (NULL);
    } else { /* all is OK */
        if (verbose_mem) {
            sprintf(ps_, ALLOC_MSG, num, cell_sizeof, ((long) p), n_alloced);
            print(ps_);
        }
        n_alloced += ((long) num) * ((long) cell_sizeof);
        return (p);
    }
}


real ***
alloc_real_3d_matrix(int i, int j, int k) {
    int z;
    real ***p;

    run {
            array(p, i, real**);
            for (z = 0; z < i; z++) p[z] = NULL;
            for (z = 0; z < i; z++) matrix(p[z], j, k, real);
        } except_when(NOMEMORY) {
        if (p != NULL) {
            for (z = 0; z < i; z++) if (p[z] != NULL) unmatrix(p[z], j, real);
            unarray(p, (real * *));
            relay;
        }
    }
    return (p);
}


void
free_real_3d_matrix(real ***p, int i, int j) {
    int z;

    if (p == NULL) return;
    for (z = 0; z < i; z++) unmatrix(p[z], j, real);
    unarray(p, real**);
}

char ***
alloc_char_3d_matrix(int i, int j, int k) {
    int z;
    char ***p;

    run {
            array(p, i, char**);
            for (z = 0; z < i; z++) p[z] = NULL;
            for (z = 0; z < i; z++) matrix(p[z], j, k, char);
        } except_when(NOMEMORY) {
        if (p != NULL) {
            for (z = 0; z < i; z++) if (p[z] != NULL) unmatrix(p[z], j, char);
            unarray(p, (char**));
            relay;
        }
    }
    return (p);
}


void
free_char_3d_matrix(char ***p, int i, int j) {
    int z;

    if (p == NULL) return;
    for (z = 0; z < i; z++) unmatrix(p[z], j, char);
    unarray(p, char**);
}


void
pmat_r(char *n, real **x, int a, int b) {
    int i, j;
    print(n);
    print(":\n");
    print("       [i][0]     [i][1]     [i][2]     ");
    print("[i][3]     [i][4]     [i][5]\n");
    for (i = 0; i < a; i++) {
        sprintf(ps, "%2d %2d ", i, 0);
        print(ps);
        for (j = 0; j < b; j++) {
            sprintf(ps, "%10.2le ", x[i][j]);
            print(ps);
            if (j % 6 == 5 && j != b - 1) {
                sprintf(ps, "\n   %2d ", j);
                print(ps);
            }
        }
        nl();
    }
}


void
pary_r(char *n, real *x, int a) {
    int i;

    print(n);
    print(":\n    [i]\n");
    for (i = 0; i < a; i++) {
        sprintf(ps, "%2d %10.2le\n", i, x[i]);
        print(ps);
    }
}


void
pary_r2x2(char *n, REAL2x2 *x, int a) {
    int i;

    print(n);
    print(":\n");
    print("    [0][0]     [0][1]     [1][0]     [1][1]\n");
    for (i = 0; i < a; i++) {
        sprintf(ps, "%2d %10.2le %10.2le %10.2le %10.2le\n",
                i, x[i][0][0], x[i][0][1], x[i][1][0], x[i][1][1]);
        print(ps);
    }
}


void
pary_r4(char *n, REAL4 *x, int a) {
    int i;
    print(n);
    print(":\n");
    print("    [0]        [1]        [2]        [3]\n");
    for (i = 0; i < a; i++) {
        sprintf(ps, "%2d %10.2le %10.2le %10.2le %10.2le\n",
                i, x[i][0], x[i][1], x[i][2], x[i][3]);
        print(ps);
    }
}


void
mem_init(void) {
    n_alloced = 0l;
    verbose_mem = FALSE;
}

