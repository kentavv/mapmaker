/******************************************************************************

 #    #    ##     #####  #    #  #          #    #####            ####
 ##  ##   #  #      #    #    #  #          #    #    #          #    #
 # ## #  #    #     #    ######  #          #    #####           #
 #    #  ######     #    #    #  #          #    #    #   ###    #
 #    #  #    #     #    #    #  #          #    #    #   ###    #    #
 #    #  #    #     #    #    #  ######     #    #####    ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#include "system.h"

real sq(real r) {
    return (r * r);
}

real rmaxf(real r, real s) {
    return (r >= s ? r : s);
}

real rminf(real r, real s) {
    return (r <= s ? r : s);
}

long lpow2(int i) {
    /* Return 2 to the i-th power for 0<=i<=LONGBITS (31) */
    long result;
    int j;

    if (i < 0 || i > LONGBITS) send(CRASH);
    for (result = (long) 1, j = 0; j < i; j++) result <<= 1;
    return (result);
}

int ipow(int x, int i) {
    /* Return 2 to the i-th power 0<=i<=INTBITS (15) */
    int result;
    int j;

    if (i < 0 || i > INTBITS) send(CRASH);
    for (result = 1, j = 0; j < i; j++) result *= x;
    return (result);
}

int ichoose(int n, int k) {
    /* Best algorithm for small k */
    long a, b, Ln, Lk, L1, Lx;

    if (k > n || n < 1 || k < 1) send(CRASH);
    Ln = (long) n;
    Lk = (long) k;
    Lx = (long) (n - k + 1);
    L1 = (long) 1;

    for (a = L1; Ln >= Lx; Ln--) a *= Ln; /* is product of n...n-k+1, or n!/(n-k)! */
    for (b = L1; Lk > L1; Lk--) b *= Lk;
    return ((int) (a / b));
}

int imaxf(int r, int s) {
    return (r >= s ? r : s);
}

int iminf(int r, int s) {
    return (r <= s ? r : s);
}


/*** Compare Functions for Array Sort Routines ***/

int icomp(const void *x, const void *y) {
    return ((*(int *) y) - (*(int *) x));
}

int lcomp(const void *x, const void *y) {
    long a, b;
    a = *(long *) x;
    b = *(long *) y;
    if (a < b) return (-1); else if (a == b) return (0); else return (1);
}

int rcomp(const void *x, const void *y) {
    real a, b;
    a = *(real *) x;
    b = *(real *) y;
    if (a < b) return (-1); else if (a == b) return (0); else return (1);
}

int scomp(const void *x, const void *y) {
    /* does this work? */
    return ((int) strcmp(*(char **) x, *(char **) y));
}

int inv_icomp(const void *x, const void *y) {
    return ((*(int *) x) - (*(int *) y));
}

int inv_rcomp(const void *x, const void *y) {
    real a, b;
    a = *(real *) x;
    b = *(real *) y;
    if (a > b) return (-1); else if (a == b) return (0); else return (1);
}


real rmean(real *data, int length) {
    int i;
    real total;

    if (length == 0) return (0.0);
    for (i = 0, total = 0.0; i < length; i++) total += data[i];
    return (total / REAL(length));
}


void mat_invert(real **m, int size, real **m_inverse) {
    /* Invert square matrix by Gauss' method */
    /* m_inverse should be 2*size columns (2nd index) by size rows (first index)
       its left square will be left with the result (ignore the right side) */
    int row, col, c, twice_size;
    real value;

    twice_size = 2 * size;
    for (row = 0; row < size; row++)
        for (col = 0; col < size; col++) {
            m_inverse[row][col] = m[row][col];
            m_inverse[row][col + size] = (row == col ? 1.0 : 0.0);
        }

    for (col = 0; col < size; col++) {
        if (m_inverse[col][col] == 0.0) {
            for (row = 0; row < size && m_inverse[row][col] == 0.0; row++);
            if (row == size) send(CRASH);
            for (c = 0; c < twice_size; c++)
                m_inverse[col][c] += m_inverse[row][c];
        }
        for (row = 0; row < size; row++)
            if (row != col) eliminate(row, col, col, m_inverse, twice_size);
    }

    for (row = 0; row < size; row++) {
        value = m_inverse[row][row];
        for (col = 0; col < twice_size; col++)
            m_inverse[row][col] /= value;
    }

    for (row = 0; row < size; row++)
        for (col = 0; col < size; col++)
            m_inverse[row][col] = m_inverse[row][col + size];

/*    if ((diff=test_mult(m,m_inverse,size))>mat_tolerance) 
      { MATINV_diff= diff; send(MATINV); } OLD */
}


void eliminate(int row, int col, int row_to_sub, real **m, int n_cols) {
    int j;
    real value;

    /* Get m[row][col]=0 by subtracting the row with a 1 in col
       from m[row], where the row is scaled appropriately */

    if (m[row][col] == 0.0) return;
    if (m[row_to_sub][col] == 0.0) send(CRASH); /* used to send(SINGMAT) */
    /* but SINGMAT was undefined */
    value = m[row][col] / m[row_to_sub][col];
    for (j = 0; j < n_cols; j++) m[row_to_sub][j] *= value;
    for (j = 0; j < n_cols; j++) m[row][j] -= (m[row_to_sub][j]);
}


void
array_times_matrix(const real *a, real **b, int rows, int columns, real *c) {
    int i, j;

    for (i = 0; i < columns; i++)
        for (j = 0, c[i] = 0.0; j < rows; j++) c[i] += a[j] * b[j][i];
}
