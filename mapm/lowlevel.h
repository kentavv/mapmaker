#ifndef _LOWLEVEL_H_
#define _LOWLEVEL_H_

/******************************************************************************

 #        ####   #    #  #       ######  #    #  ######  #               #    #
 #       #    #  #    #  #       #       #    #  #       #               #    #
 #       #    #  #    #  #       #####   #    #  #####   #               ######
 #       #    #  # ## #  #       #       #    #  #       #        ###    #    #
 #       #    #  ##  ##  #       #        #  #   #       #        ###    #    #
 ######   ####   #    #  ######  ######    ##    ######  ######   ###    #    #

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#include "map_info.h"


/********************** Raw Data Struct ************************************/
  			
/* symbols used for F2 data preparation and storage */
#define PARENTAL_TYPE_A           'A'
#define PARENTAL_TYPE_B           'B'
#define TYPE_NOT_A                'C'
#define TYPE_NOT_B                'D'
#define HYBRID_TYPE_H             'H'
#define MISSING_DATA              '-'

/* Definition for f2 raw data storage. */
typedef struct {
    int         num_indivs;
    int         cross_type;              /* F2_INTERCROSS or F2_BACKCROSS */
    real        **allelic_distribution;  /* [loci][4: AA,AB,BA,BB] */
    char        **allele;                /* [loci][indivs] */
}       F2_RAW;

typedef struct {
    char        filename[PATH_LENGTH+1]; /* raw data file */
    int         filenumber;       	 /* for concurrency checks */
    int         data_type;        	 /* defined above */
    int 	num_markers;
    char        **locus_name;     	 /* [loci][MAX_NAME_LEN] */
    union {
#ifdef HAVE_CEPH
    	CEPH_RAW      	ceph;             /* definition in ceph.h?? */
#endif
    	F2_RAW   	f2;
    } data;
} RAW_DATA;

extern RAW_DATA raw;

/****************************************************************************/


/* Guts of things used only inside converge_to_map. Much of this is for the
   old_ctm implementation of F2 etc will be obsoleted by HMM... */

#define RIGHT 		1
#define LEFT		0
//#define NEGATIVE        1
//#define POSITIVE        0

typedef real   RECVECTOR[2];
//typedef real   RECARRAY[2][2];
//typedef char   CONV_HIST;
//typedef int    (*PFI)();   /* pointer to function returning an integer */
//typedef int    (*PFB)();   /* pointer to function returning a bool(int) */

/* Locus structure for converge_to_map() down.  Used to store
   the markers that converge_to_map() is being used on. */
typedef struct {
    int         size;
    int         *Entry;
} LOCUS;

#ifdef HAVE_CEPH
/*****************************************************************************
                       COOKED CEPH-DATA
*****************************************************************************/

/* Processed data storage for phase known */
typedef struct {
    int         num_markers;
    int         num_families;
    FLDB        **fldb;      /* as in CEPH_RAW, KNOWN data needs */
}       KNOWN_COMMUNITY;             /* no preprocessing */

/* Processed data for phase unknown program */
typedef struct {
    event_vector event[2];      /* [event] */
    real prob[2];     	/* [side] */
    real ap_prob;      	/* a priori probability */
} INHER_VECTOR;
    
typedef struct {
    int eff_locus[2];		/* effective parental locus */
    int non_informative;	/* true iff non-inf for both parents */
    long index;			/* size of dist */
    INHER_VECTOR *dist;		/* probability distribution */
} INHER_STR;

typedef struct {
    INHER_STR *genes;
    int kids;		        /* number of kids in family */
} FAMILY;

typedef struct {
    FAMILY *family;     
    int num_fams;	        /* number of families */
} UNKNOWN_COMMUNITY;

/* Definitions used in the count_recs() procedures */
#define LOCUS_AT_INFINITY -1
#define TWO_TO_MAX_KIDS 32768

typedef struct {
    long size;
    int eff_locus;
    event_vector Entry[TWO_TO_MAX_KIDS];
} PREV_LOCUS;

#endif


/*****************************************************************************
			COOKED F2 TYPE DATA
*****************************************************************************/

//#define NEUTER                      0
//#define APRIORI                     2

/* Definition for processed backcross data storage. */
/* As with phase known data, backcross data needs no preprocessing. */
typedef struct {
    int         num_markers;
    int         num_indivs;
    char        **allele;         /* [loci][indivs], as in F2_RAW */
}       BACKX_GENERATION;

/* Processed intercross data storage. */
typedef struct {
    int 	index;
    bool     non_informative;
    int         event[4];        /* event: AA,AB,BA,BB */
    real      prob_dist[3][4]; /* [apriori,left,right][event: AA,AB,BA,BB] */
}		F2_INHER_STR;
    
//typedef real PROB_DIST[3][4];  /* as in F2_INHER_STR */

typedef struct {
    int             number_indivs;
    F2_INHER_STR    **data;     
}                   F2_GENERATION;


/* processed data structure - created in the various init_for_em procedures */
typedef union {
#ifdef HAVE_CEPH
    KNOWN_COMMUNITY	*known;        /* in ceph.h */
    UNKNOWN_COMMUNITY	*unknown;      /* in ceph.h */
#endif
    F2_GENERATION	*f2;           /* in f2.h */
    BACKX_GENERATION    *backcross;    /* in f2.h */
}	PROCESSED_DATA;


/* Functions */
/*in ctm.c */
int converge_to_map_mt(MAP *map, int index);  /* arg: MAP; calls map making functions and returns
                            the converged map (alters: rec_frac, log-like) */
int converge_to_map(MAP *map);  /* arg: MAP; calls map making functions and returns
                            the converged map (alters: rec_frac, log-like) */
/* in multipt.c */
//void allocate_recs();    /* These procedures are called by converge_to_map */
//void allocate_temps();   /* and assist in the map making. */
//void free_memory_temps();
//void no_history();
//bool norm_conv_rule();
//bool norm_inner_conv_rule();
//bool inner_converge_instantly();
//bool converge_instantly();
//void norm_make_new_map();

///* map making procedures (separated by data type) */
///* in interx.c */
//void f2_init_for_em();
//void f2_count_recs();
//void f2_free_memory_from_em();
///* the following are also in interx.c but are used by other functions */
//int symbol_value();
//int changes();
//real f2_prob(), f3_prob();
//real power();

///* in known.c */
//void known_init_for_em();
//void known_count_recs();
//void known_free_memory_from_em();
//
///* in backx.c */
//void backcross_init_for_em();
//void backcross_count_recs();
//void backcross_free_memory_from_em();
//
///* in unk_init.c */
//void unk_init_for_em();
//void unk_free_memory_from_em();
///*in unk_cnt.c */
//void unk_count_recs();
//
///* in new_count.c */
//void new_unk_init_for_em();
//void new_unk_free_memory_from_em();
//void new_unk_count_recs();


/* Constant definitions useful in converge_to_map() and 
   lower map making procedures. */
#define NUM_CONV_RULES 1
#define NORMAL_CONVERGENCE 0
#define NUM_INNER_CONV_RULES 1
#define NORM_INNER_CONVERGENCE 0

#define UNINFORM	-1             /* Defs for marker types */
#define INFORM		0              /*   and zygosity.  	*/
//#define HETERO		INFORM
//#define HOMO		UNINFORM
//#define EITHER		4
//#define FIRSTHETERO	2

#define REC 		1
#define NOREC 		0
//#define for_recs(r)	for(r=REC;r>=NOREC;r--)

#define SMALLRECOMBS    0.00001
#define BIGRECOMBS      0.50

#endif
