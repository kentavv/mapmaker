#ifndef _QMAP_H_
#define _QMAP_H_

/******************************************************************************

  ####   #    #    ##    #####           #    #
 #    #  ##  ##   #  #   #    #          #    #
 #    #  # ## #  #    #  #    #          ######
 #  # #  #    #  ######  #####    ###    #    #
 #   #   #    #  #    #  #        ###    #    #
  ### #  #    #  #    #  #        ###    #    #

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

/***** qmap.h - Everything one needs to call QCTM should be in here. *****/

/****************** CALLING CONVENTIONS FOR QCTM ******************
Before entry, QCTM requires that pos_tolerance, like_tolerance, and
mat_tolerance be set, and that qctm_alloc be run (which requires that
max_intervals and max_individuals be set). Data must be processed with
prepare_data() for each set of intervals to run through qctm.  In
addition, the map should be run through initial_qctm_values(). If
certain values in the data struct do not correspond to those in the
map (data->num_intervals and map->n_intervals, for example) then god
knows what will happen.

In the map struct, qtl_weight, qtl_pos, mu, sigma_sq, null_mu, and
null_sigma_sq should be initialized, as should qtl_dominance if the
data is f2 intercross. Note that this can all be done by the canned
procedures provided. All of these variables, as well as var_explained,
chi_sq, null_log_like, abs_log_like, and log_like will be
side-effected!  Map->left and map->right are ignored, except in
print_iter (see print_iterations(), defined elsewhere).  Qtl_pos and
qtl_weight must be #intervals and #intervals+1 (IMPORTANT!) in length,
respectively.

Throughout the qctm code, qtl_pos and qtl_weight are arrays, while mu
and sigma_sq, are simply pointers to single reals. QCTM may possibly
send messages CRASH, SINGMAT, or MATMULT, which indicate failure. In
this case, the QTL_MAP results are undefined. No memory is allocated
by qctm (for efficiency, all big things are static), and thus no
cleanup after these messages is required.
*********************************************************************/

/* Pretty much everything that happens to data produced by qctm does 
   (and should) happen using the QTL_MAP struct... */

/* Constraints for the intercross weight and dominance terms:
   a * weight + b * dominance = c */

typedef struct {
    real backx_weight;  /* used for backx only, the rest are for interx only */
    int interx_type;    /* defined immediately below */
    real a, b, c;
} GENETICS;
#define FREE        0
#define DOMINANT    1
#define RECESSIVE   2
#define ADDITIVE    3
#define CONSTRAINED 4
#define FIXED       5 
#define TEST_MODELS 6  /* TEST_MODELS is only allowed inside QTL_SEQUENCEs */
#define NUM_MODELS  4  /* 0-3 as FIXED and CONSTRAINED models don't count */
#define F3DOMINANT  7
#define F3RECESSIVE 8

typedef struct { 
    int trait;
    int num_intervals;
    int num_continuous_vars;
    int max_intervals;			/* max # allocated for */
    int max_continuous_vars;
    int *left, *right;			/* [interval#] => left & right loci */
    real *interval_len; 		/* [interval#] => a rec frac */
    real *qtl_pos;			/* [interval#] also is a r.f. */
    real *fix_pos;	     	        /* [interval#] is DONT_FIX or a r.f. */
    real *qtl_weight;			/* [interval#] */
    real *qtl_dominance;                /* [interval#] only for intercross! */
    int  *cont_var;			/* [cont_var#] => a valid trait# */
    real *cont_var_weight; 		/* [cont_var#] */
    real *fix_cont_var_weight; 		/* [cont_var#] */
    GENETICS *constraint;               /* an array of #intervals structs */
    real mu, sigma_sq, var_explained, chi_sq;
    real null_mu, null_sigma_sq, null_log_like;
    real log_like, no_data_like, abs_log_like;   
} QTL_MAP;

#define DONT_FIX   OBSCURE_REAL         /* use for fix_pos or fix_weight */
#define EPISTASIS_TERM -1 		/* special cont-var */

/*** handy functions in QCTM.C ***/

/*** things in QDATA.C ***/
bool qctm_globals_avail();
void alloc_qctm_globals();
void free_qctm_globals();
real model_prediction(); /* args: QTL_MAP *map; int indiv; */

/*** things in QTOP.C ***/
QTL_MAP *alloc_qtl_map(); /* args: int n_intervals, n_continuous_vars; */
void free_qtl_map();
bool reset_map();
void really_reset_map();
int  add_interval();
void mapcpy();
void make_qtl_map(); /* args: QTL_MAP *map; sets map->trait the runs 
   prepare_data(), initial_qctm_values(), and qtl_conv_to_map() */
/* BAGGED THIS:
real qtl_map_like();  like make_qtl_map(), but it assume that the map->
   weight, dominance, and cont_var_weight values have been set, and it
   calls qtl_noconv_to_map(). The likelihood is returned AND map->log_like
   is set to it. */

void copy_genetics();
bool constrained(); /* args: GENETICS *genetics; */


/*** random things ***/
#define INF_LOCUS 	-1
#define NO_LOCUS 	-2
#define ANY_LOCUS 	-3
#define VERY_UNLIKELY	-1e30
/* #define NO_TRAIT	-1 */

#endif
