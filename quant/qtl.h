/******************************************************************************

  ####    #####  #               #    #
 #    #     #    #               #    #
 #    #     #    #               ######
 #  # #     #    #        ###    #    #
 #   #      #    #        ###    #    #
  ### #     #    ######   ###    #    #

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

/* This is the only file which understands the order in which #include
   files need to be included! */

#ifdef INC_QTOPLEVEL
#define INC_TABLE     /* NOT included by INC_LIB */
#define INC_CALLQCTM
#endif

#include "system.h"

/******************** GLOBAL DECLARATIONS SPECIFIC TO QTL ********************/
/* Generally speaking, this declares QTL commands, user accessible state 
   variables, _init procedures, and a very few globally useful functions. */

/* Extensions for MAPMAKER files... */
#ifndef DOS
#define RAW_EXT   ".raw"
#define TEMP_EXT  ".temp"
#define DATA_EXT  ".data"
#define DATA_OLD  ".xdata"
#define MAPS_EXT  ".maps"
#define MAPS_OLD  ".xmaps"
#define TRAIT_EXT ".traits"
#define TRAIT_OLD ".xtraits"
#define QTL_EXT   ".qtls"
#define QTL_OLD   ".xqtls"
#else
#define RAW_EXT   ".raw"
#define TEMP_EXT  ".tmp"
#define DATA_EXT  ".dat"
#define DATA_OLD  ".xda"
#define MAPS_EXT  ".map"
#define MAPS_OLD  ".xma"
#define TRAIT_EXT ".tra"
#define TRAIT_OLD ".xtr"
#define QTL_EXT   ".qtl"
#define QTL_OLD   ".xqt"
#endif



/* global defs */
#define SEQ_LEN MAXLINE

/* allowed data types */
#define NO_DATA     0
#define INTERCROSS  1
#define BACKCROSS   2
#define ANY_DATA  (-1)

/* These are the absolute limits for max_intervals and max_individuals. 
   The defaults are set in cmd_init() (in qtop.c) Note that max_intervals is 
   limited because 2^(2*num_intervals) must fit in a GENOTYPE. */
#define GENOTYPE 	  short
#define MAX_INTERVALS     7
#define MAX_GENOTYPE_VARS 15
#define MAX_CONTINUOUS_VARS 3
#define MAX_INDIVIDUALS   5001

#define RECFRACS 0
#define CENTIMORGANS 1

/*** functions in QTOP.C ***/
void top_init();   /* no args */

/*** in QCTM.C ***/
void qctm_init();  /* no args */
real haldane();    /* args: real rec_frac; returns centiMorgans (not Morgans)*/
real unhaldane();  /* args: centimorgans;  returns recfrac; CHANGE THESE!!!! */
real kosambi();    /* args: real rec_frac; returns centiMorgans (not Morgans)*/
real unkosambi();  /* args: centimorgans;  returns recfrac; CHANGE THESE!!!! */

/* the functions that the programmers call */

#define haldane_cm(theta) (100.0 * haldane(theta)) 
#define unhaldane_cm(cm)  (unhaldane(cm/100.0))

/*** in QDATA.C ***/
void data_init();  /* no args */
extern int map_function;
#define HALDANE 0
#define KOSAMBI 1

/*** in QSEQ.C ***/
void seq_init();   /* no args */
char *expand_named_entries(); /* takes an interval_string and returns a 
				 newly allocated one with names expanded */

/*** in QCMDS.C ***/
void cmd_init();   /* no args */
void help_init();  /* no args - includes code qtl.code generated by makehelp */

/*** in QWIGGLE.C ***/
void wiggle_init();  /* no args */

/*** QTL user state variables  ***/				/* file */
/* all (most?) of these are initialized by cmd_init() */
extern int units, print_names, print_mapm_loci, print_scans; 	/* qtop.c */
extern real pos_tolerance, like_tolerance;			/* qctm.c */
extern real mat_tolerance;					/* qctm.c */
extern int print_iter, print_rec_mat, bag_qctm, debug_qctm;	/* qctm.c */
extern bool print_brute_force, debug_newton, brute_force;		/* qctm.c */
extern int max_intervals, max_continuous_vars;			/* qdata.c */
extern int segregation_distortion; /* NOT USED YET! */
extern bool altered_chroms;                                     /* qraw.c */

/*** Messages and message variables used by the QTL program ***/
/* All are declared in qtop.c, and set up by top_init(). Some of these will 
   be obsoleted soon! */
#define SINGMAT 	USER_MESSAGE(0)
#define MATINV	 	USER_MESSAGE(1)
#define DOFREE 		USER_MESSAGE(3)
#define BADDATA		USER_MESSAGE(4)
#define QUIT_QTL	USER_MESSAGE(6)
#define NOPARSE		USER_MESSAGE(8)
#define BAD_MKINT	USER_MESSAGE(9)
#define BADSEQ		USER_MESSAGE(10)
#define BADTRAIT	USER_MESSAGE(11)

extern int   BADDATA_line_num;
extern char *BADDATA_ln; /* a string (MAXLINE+1 long) which is side-effected */
extern char *BADDATA_error; /* a char ptr which is set */
extern real  MATINV_diff;
extern int   NOPARSE_err;
extern int   BAD_MKINT_err;
extern int   BADSEQ_errpos;
extern char *BADSEQ_errmsg; /* a ptr which is set (usually to msgstr) */
extern char *BADTRAIT_errmsg; /* also a ptr which is set (usually to msgstr) */
int dum_loc;      /* Used to set proper number of loci in data set */

/* Various Kludges to go from the old to new helpers library!  */
#define INSIZE LINE*2
#define why_(str) nstrcpy(BADDATA_error,str,LINE)


/*** QTL #include files ***/

#ifdef INC_CALLQCTM
#include "qmap.h"
#endif

#ifdef INC_QTOPLEVEL
#include "qtop.h"    /* include after qmap.h */
#endif

#ifdef INC_QLOWLEVEL
#include "qlow.h"
#endif




