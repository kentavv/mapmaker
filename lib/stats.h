#ifndef _STATS_H_
#define _STATS_H_

/******************************************************************************

  ####    #####    ##     #####   ####           #    #
 #          #     #  #      #    #               #    #
  ####      #    #    #     #     ####           ######
      #     #    ######     #         #   ###    #    #
 #    #     #    #    #     #    #    #   ###    #    #
  ####      #    #    #     #     ####    ###    #    #

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

typedef struct normal_test {
    real skew;
    real kurt;
    real mean;
    real sigma;
    real quart_ratio;
    real within_one_fourth;
    real within_one_half;
    real within_one;
    real within_two;
    real within_three;
} NORMAL_TEST;

NORMAL_TEST *check_normalcy(real *dist, int dist_size);

void skew(NORMAL_TEST *normal_results, real *dist, int dist_size);

void quartile(NORMAL_TEST *normal_results, real *dist, int dist_size);

void within(NORMAL_TEST *normal_results, real *dist, int dist_size);

void print_normal(NORMAL_TEST *to_be_printed, real lamda);

void print_rhisto(real *dist, int dist_size);

/* The function normal_test will do several tests on an array of
   values to see how normal of a distribution it is.  As it makes
   a test, it will then store the results in a struct of type normal
   test.  The reason for this is that another function box-cox will be 
   used on the data to determine which is the most normal representation
   of the data. */

#endif