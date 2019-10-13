#ifndef _QTOP_H_
#define _QTOP_H_

/******************************************************************************

  ####    #####   ####   #####           #    #
 #    #     #    #    #  #    #          #    #
 #    #     #    #    #  #    #          ######
 #  # #     #    #    #  #####    ###    #    #
 #   #      #    #    #  #        ###    #    #
  ### #     #     ####   #        ###    #    #

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#include "table.h"

/***** qtop.h - Declarations of things which help the QTL commands
deal with data produced by QCTM - this file requires that qmap.h and
table.h be included first! *****/

/***** QTL INTERVALS LISTS (SEQUENCES) *****/

/* The compiled interval specification is a linked list... Only code in qseq.c
   should care what this really looks like. This should be a struct of unions,
   but it isn't. */

typedef struct qtl_seq_option {
    int type;            /* one of the defines below */
    union {
        struct {
            int left, right;
            real fix_pos;
        } interval;
        struct {
            int from, to;
        } range;
        struct {
            int num, dummy;    /* dummy is here for debugging purposes only! */
            real fix_weight;
        } cont_var;
    } isa;
    struct qtl_seq_option *next;
} QTL_SEQ_OPTION;

#define GARBAGE     0    /* possible types for an INTERVALS struct */
#define INTERVAL 1
#define RANGE     2
#define CONT_VAR 3

typedef struct qtl_sequence {
    QTL_SEQ_OPTION *option; /* Linked list of the possibilities for this int */
    int repeat;
    struct qtl_sequence *next;
    bool dont_free;               /* If a global struct pts here, dont free */
    int num_options;              /* total num options for this interval */
    bool isa_continuous_var;      /* If FALSE, fill in vars left...genetics */

    int *left, *right;            /* [#opts] enumerate the possibilities */
    real *fix_pos;                /* for each possibility */
    bool *contig;                 /* flag discontinuities w/FALSE */
    GENETICS genetics;

    int *cont_var;          /* iff isa_cont_var: [#opts] trait numbers */
    real *fix_cont_var_weight;      /* iff isa_cont_var: [#opts] weights */
    int *count;          /* permutation state, array <repeat> long */
    int test;                     /* permutation state for test genetics */
    real interval_cm, pos_cm;     /* permutation state for wiggle */
} QTL_SEQUENCE;


/*** things in QSEQ.C ***/
extern QTL_SEQUENCE *ints;
extern char *ints_string;

typedef struct {
    int trait;
    int seq_history_num;
    TABLE *named_sequences;
    TABLE *sequence_history;
} STATUS_CONTEXT;

extern STATUS_CONTEXT **context;
extern int num_contexts, active_context;

//void allocate_context();
void context_init(void);

//bool change_context();
//bool create_new_context();
void kill_context(STATUS_CONTEXT *con, bool save_it);

void free_context(STATUS_CONTEXT *con);

/* This the standard ways to loop through the intervals in the qtl_sequence: */
#define for_all_orders(seq, map, perm)                \
for (Omore=TRUE, perm=FIRST;                    \
     Omore && reset_map(map) && get_order(seq,FALSE,map);    \
     Omore=next_order(seq,&(perm)))

#define for_wiggle_orders(seq, map, step, perm)            \
for (Omore=TRUE, perm=FIRST;                    \
     Omore && reset_map(map) && get_order(seq,TRUE,map);    \
     Omore=next_wiggle(seq,&(perm),step))
extern int Omore;

/* These constants reflect the order in which permutations occur. It
in this case IS OK to assume the numbering, and to use >, <, etc when
testing the value of perm. In for_all_orders, perm may be set to
FIRST, CONTIG,  TEST, or ORDER.  In for_wiggle_orders, perm= FIRST,
WIGGLE, CONTIG, SKIP, TEST or ORDER */

#define FIRST  0  /* The first permutation */
#define WIGGLE 1  /* wiggle only: One step inside an interval */
#define CONTIG 2  /* A new, however contiguous, interval or set of intervals */
#define SKIP   3  /* wiggle only: A new, non-contiguous interval */
#define TEST   4  /* Permuted the genetics - considered discontiguous */
#define ORDER  5  /* A discontiguous perm. Any left seq permutation in wiggle,
		     including a test perm, is flagged as ORDER. */

bool set_qtl_sequence(char *str, char *errmsg, int *errpos); /* args: char *str, *errmsg; int *errpos;
   Side-effect globals ints and ints_string, returns TRUE if all is ok or
   FALSE on error, in which case neither ints nor ints_string are touched,
   while the args errmsg and errpos are side-effected with an error message
   and error position in the string. When ok, the old ints struct is freed. */

/* The lower level functions... */

bool reset_state(QTL_SEQUENCE *p, bool wiggle, int *pnum_intervals, int *pcont_vars, int *pnum_orders, int *pwiggle_ints); /* args: QTL_SEQUENCE *seq; bool wiggle; int *num_ints;
		       int *num_orders, *num_wiggle_ints; 
   Resets the permutation state of the seq struct and also counts the
   #intervals per order and the #orders specified (either pointer may be
   NULL). This function is usually called from qtl_ready() in qtop.c. If 
   wiggle is true, the num_orders is for the intervals to the right of the
   wiggled interval, and *num_wiggle_ints is set to point to the 
   #possibilities for the rightmost qtl. Also in this case, FALSE may 
   be returned if the sequence is not appropriate for for_wiggle_orders(),
   otherwise TRUE is returned. */

bool get_order(QTL_SEQUENCE *p, bool wiggle, QTL_MAP *map);   /* args: QTL_SEQUENCE *seq; bool wiggle; QTL_MAP *map;
   For a permutation of the seq struct, this APPENDS that sequence of
   intervals to those in a QTL_MAP - this always returns TRUE */

bool next_order(QTL_SEQUENCE *p, int *perm); /* args: QTL_SEQUENCE *seq; int *perm;
   Permutes the seq struct and returns TRUE if there is another
   permutation, or returns FALSE otherwise. perm is side-effected as 
   described above. */

bool next_wiggle(QTL_SEQUENCE *p, int *perm, real cm_step); /* args: QTL_SEQUENCE *seq; int *perm; real *cm_step;
   Like next_order(), except the rightmost qtl is wiggled through at the 
   cm_step rather than more simply permuted. */

QTL_SEQUENCE *compile_intervals(char *str); /* args: char *str; sends msg on error */

//void free_qtl_sequence(); /* args: QTL_SEQUENCE *ints; */

void getdataln(FILE *fp);

bool name_sequence(char *name, char *seq, char **why_not);

bool unname_sequence(char *name, char **why_not);

void add_to_seq_history(char *seq);

bool valid_locus_str(char *str, int *num, char *errmsg);

bool valid_locus_num(int *num);

void get_compare_nums(char *str, int *compare, int *contig);

//bool valid_locus_str();     /* which of these are obsolete? */
//bool valid_interval_num();
//bool valid_locus_num();
//
//bool name_sequence();
//bool unname_sequence();
//void add_to_seq_history (char *seq);
//void add_to_seq_history();
//void name_peaks();
#define MAX_HISTORY_SEQS 23
#define MAX_NAMED_SEQS   44
//void getdataln();
//void get_compare_nums();

/* For qctm.c */


/*** things in QTOP.C ***/
void qtl_ready(int data_type, int need_seq, bool need_trait, bool will_call_qctm); /* args: int data_type, need_seq; bool need_trait,calls_qctm;
   Checks if qtl is set up as required, and calls error (which sends
   SOFTABORT, see shell.h) if she's not. In this process, reset_state() is
   called. If she is and if min_interval_permutations>0 then
   *num_intervals and *num_orders are side-effected (if they are
   non-NULL). If need_seq==WIGGLE then *num_orders will indicate the
   #permutations of the left part of the sequence, and *num_wiggle_ints
   will indicate the #intervals specified for the rightmost qtl. */

#define SEQ         1
#define NOSEQ       0
#define WIGSEQ      (-1)
#define ONEINT      2    /* FROB */
#define TRAIT       TRUE
#define NOTRAIT     FALSE
#define QCTM        TRUE
#define NOQCTM      FALSE

/* Globals for commands - set before the command procedure is invoked */
extern int trait;
extern char *trait_string;
extern int print_maps;

/* These globals are all set by qtl_ready() */
extern int num_intervals, num_orders, num_ints_to_wiggle;
extern bool *free_genetics; /* [interval#] => TRUE if genetics are free */
extern QTL_MAP *map; /* allocated for max_intervals */
//char *trait_str();
bool valid_trait_str(char *str, int *num, char *errmsg);

bool valid_new_trait_name(char *str, char *errmsg);

bool valid_trait_num(int num);

void set_trait_spec(char *str);

//bool valid_trait_str(); /* args: char *str; int *trait_num; char *error_msg;
//   side-effect *trait_num and return TRUE if str is a valid trait number or
//   name, or return FALSE and side-effect error_msg (WHICH MAY NOT BE NULL!)
//   otherwise. str must contain only a single token, or you will get FALSE. */
//
//void set_trait_spec(); /* args: char *str; sets the global vars trait and
//   trait_string if str passes valid_trait_num, or sends BADTRAIT message */
//
//bool valid_new_trait_name(); /* args: char *str, *error_msg; return TRUE
//   if str is a valid name for a new trait. Return false and side-effect
//   error_msg (WHICH MAY NOT BE NULL!) otherwise. str must contain only a
//   single token, or you will get FALSE. */
//
//bool valid_trait_num(); /* args: int num; return TRUE if num is a valid trait
//   number FROM THE COMPUTER'S PERSPECTIVE (that is 0...#traits-1), not the
//   user's perspective (which is 1...#traits). */
//
//bool new_trait_num(); /* args: int *num; return TRUE if there is room for
//   another trait in the raw struct, and side-effect *num to be the appropriate
//   number for it. Return FALSE otherwise. */


/*** To save wiggle data... ***/

/* The WIGGLE struct holds the LOD scores, etc, for a 'wiggle' search
for qtls. It is assumed that the rightmost interval is the one
'wiggled' across the intervals specified - if various left orders are
allowed by the sequence or command, the data are stored in a 2D
matrix, indexed by the left permutation# and the wiggled interval#.
This way, we could allow a WIGGLE struct to hold the results of a
command such as: 'search the genome for a 2nd qtl for trait N,
assuming that there is one known qtl somewhere between markers 12 and
15.' Wiggles testing many possible genetic models for the rightmost
interval will be considered as multiple left perms. 

The wiggle struct is dually indexed - all data are collected by the
operation which generated them into the wiggles struct. A COMPARE
struct is used to store the same kind of data for non-wiggle
searches. Wiggles and/or compare seraches for single interval models
(allowing test_genetics) are also collected for each trait for the
entire genome into the qtls struct. All of these structs are saved. */

typedef struct {
    int trait;            /* the trait and sequence should uniquely */
    QTL_SEQUENCE *seq;        /* specify what this test did */
    char *seq_string;
    int num_intervals;    /* rightmost interval is the one "wiggled" */
    int num_orders;             /* #permutations of the left part of the seq */
    int num_wiggled_intervals;  /* #intervals for rightmost interval in seq */
    int order, wiggle_interval; /* state used for building these things */
    struct wiggle_int ***data;  /* [order#][interval#]=> ptr. If data==NULL */
} WIGGLE_OPERATION;             /* then this struct isn't filled in! */

typedef struct wiggle_int {
    bool contig;          /* FALSE marks the start of a discontinuity */
    QTL_MAP *map;         /* Map for best_point, which may be calculated by
			     max-like. It  shouldn't be NULL. In qtls struct
    			     only, if point is NULL, only *map is filled in. */
    bool max_like_map;    /* TRUE if the map is from a ML calculation */
    bool in_qtls_struct;  /* TRUE if a ptr to this struct is in qtls */
    WIGGLE_OPERATION *op; /* back-ptr to the operation which made these data */
    real cm_increment;
    int num_points;
    struct wiggle_point **point; /* [data_point#] => ptr to struct below */
    int best_point;              /* data_point# for the highest LOD */
    int point_num;               /* state for building these things */

} WIGGLE_INTERVAL;

typedef struct wiggle_point {
    real qtl_pos;
    real lod_score;
    real var_explained;
    real qtl_weight;
    real qtl_dominance;          /* set in INTERCROSS only */
} WIGGLE_POINT;

/* Wiggles is a global collection of all wiggle data: [operation#]=>ptr */
extern WIGGLE_OPERATION **wiggles;
extern int num_wiggles, max_wiggles, first_wiggle;

/* The COMPARE struct is similar...  */

typedef struct {
    int trait;
    QTL_SEQUENCE *seq;
    char *seq_string;
    int num_intervals;
    int num_orders;
    int num_contigs;
    int order;
    struct comp_data **data;
} COMPARE_OPERATION;

typedef struct comp_data {
    bool contig;    /* FALSE marks the start of a discontinuity */
    QTL_MAP *map;   /* map struct for the Max-Like qtl_pos */
    COMPARE_OPERATION *op;
} COMPARE_DATA;

extern COMPARE_OPERATION **compares;
extern int num_compares, max_compares, first_compare;

/* Qtls is a global collection of the wiggle data for single interval models,
   possibly allowing some genetics constraints. When indexed by 
   [trait#][model#][interval#] we get a pointer to a WIGGLE_INTERVAL struct 
   which is also in the wiggles collection, and has a back ptr to the 
   WIGGLE_OPERATION which created it (which is NULL if the search used 
   compare rather than wiggle. 
   NOTE: THIS HAS ONLY BEEN HALF IMPLEMENTED, MAY NOT WORK, AND SOON MAY BE
   DELETED FOR A SLIGHTLY BETTER MECHANISM. */
/* extern WIGGLE_INTERVAL ****qtls; */


void allocate_qtl_struct(int n_wiggles, int n_compares); /* args: int n_wiggles, n_compares; Allocates
   most of the global structs qtls, wiggles, and compares, using globals 
   raw.max_traits, raw.data_type, and raw.n_loci as params. Thus, this should
   happen upon data loading. Parts of these structs are filled in however. */

int allocate_wiggle_struct(int trait, QTL_SEQUENCE *seq, char *seq_str, int n_intervals, int n_orders, int n_wiggled); /* args: int trait; QTL_SEQUENCE *seq;
				 char *seq_str; int n_intervals, n_orders; 
				 int n_wiggled_intervals; 
   Returns and index into the wiggles struct, or -1 if it can't (maybe
   it's full - we need to add code to allow recycling or expansion). The
   WIGGLE_INTERVAL structs are allocated, although they are not filled in
   in any way. Note that the args n_intervals, n_left_orders, and
   n_wiggled_intervals are not maximums - they must be precisely correct
   (this is enforced). */

void store_wiggle_interval(int wiggle_num, QTL_MAP *map, bool new_left_order, bool contig, real cm_step); /* args: int wiggle_num; QTL_MAP *map;
				 bool new_left_order, contig; real cm_step;
   Initiates storing of wiggle data for the specified interval. The
   data wil be stored in both wiggles[wiggle_num] and qtls if
   appropriate. This function allocates the space for the data points,
   and thus assumes that the for_wiggle_orders() macro is being used in
   how much space it allocates. */

void store_wiggle_point(int wiggle_num, QTL_MAP *map); /* int wiggle_num; QTL_MAP *map; Stores the
   map data in the appropriate WIGGLE_DATA_POINT struct. This also 
   assumes that for_wiggle_orders() is being used! */

void bash_wiggle_struct(int n); /* args: int index; This kills the wiggle struct
   entry number index. No effort is made to recycle it or its slot in 
   the wiggles struct, and this is really intended for exception handling
   in the wiggle command. This is a KLUDGE for now. */

int allocate_compare_struct(int trait, QTL_SEQUENCE *seq, char *seq_str, int n_intervals, int n_orders); /* args: int trait; QTL_SEQUENCE *seq;
				  char *seq_str; int n_intervals, n_orders;
   Like allocate_wiggle_struct(). */

void bash_compare_struct(int n); /* args: int n; like bash_wiggle_struct(). */

void store_compare_map(int compare_num, QTL_MAP *map, bool contig); /* args: int compare_num; QTL_MAP *map;
				  bool contig;
   Kind of like store_wiggle_interval(). */
void get_seq_free_genetics(QTL_SEQUENCE *p, bool *free);

void get_wiggle_nums(char *str, int *wiggle, int *order);

void save_wiggle(FILE *fp, int n);

void load_wiggle(FILE *fp);

void save_compare(FILE *fp, int n);

void load_compare(FILE *fp);
//
//void get_wiggle_nums(); /* args: int *wiggle_num, *order_num; gets cmd args */
//bool isa_wiggle_interval(); /* args: int wiggle_num; */
//void save_wiggle();
//void save_compare();
//void load_wiggle();
//void load_compare();

/* The WIGGLE_PEAK struct is a distilation of wiggles/qtls data... */

typedef struct wiggle_peak {
    int left, right;
    real qtl_pos;
    real lod_score;
    real var_explained;
    real qtl_weight;
    real qtl_dominance;     /* dominance set for INTERCROSS only */

    int backward_left, backward_right;
    real backward_pos;
    int forward_left, forward_right;
    real forward_pos;

    QTL_MAP *map;           /* May be NULL */
    bool free_map;          /* When freeing this struct, free the map too? */

    struct wiggle_peak *next; /* a linked list */
} WIGGLE_PEAK;

void name_peaks(WIGGLE_PEAK *peak, char *prefix, bool forget);

WIGGLE_PEAK *find_wiggle_peaks(int wiggle_num, int left_order_num, real threshold, real qtl_falloff, real confidence_falloff, real min_peak_delta,
                               bool get_peak_maps);

#define OFF_END (-1)

void free_wiggle_peaks(WIGGLE_PEAK *p);

void print_qtl_map(QTL_MAP *map, bool *free_genetics);

void map_printer(QTL_MAP *map, bool print_genetics);

void print_short_title(void);

void print_short_qtl_map(QTL_MAP *map, real threshold, real scale);

void print_iteration(int iter, QTL_MAP *map, real delta_log_like);

void print_null_iteration(QTL_MAP *map);

//void do_print_E_step(real **expected_genotype, real **S_matrix, GENO_PROBS *expected_recs, int n_individuals, int n_genotype_vars, int n_continuous_vars, int n_intervals);
void print_wiggle_title(void);

void print_wiggle_interval(QTL_MAP *map);

void print_wiggle_map(QTL_MAP *map, real base_like, real scale);

void print_wiggle_genetics(GENETICS *genetics);

void print_wiggle_left_seq(QTL_MAP *map);

void print_saved_wiggles(void);

void print_saved_wiggle(int wiggle);

void print_saved_wiggle_order(int wiggle, int order, real base_like, real scale);

void print_peak(WIGGLE_PEAK *peak, int num);

void print_test_wiggle_map(WIGGLE_POINT **point, real threshold);

void print_test_wiggle_interval(QTL_MAP *map);

void print_test_wiggle_title(void);

void print_test_wiggle_order(int wiggle, int order, real threshold);

void print_trait(int for_num_maps);

void print_seq(void);

void print_old_seq(char *str);

char *trait_str(void);

char *interval_str(int left, int right, bool fill);

char *dist_str(real rec_frac, bool fill);

//char *units_str(bool fill);
char *genetics_str(GENETICS *genetics, bool verbose);

//char *left_seq_str(QTL_MAP *map);
void print_tiny_map(QTL_MAP *map, int num, real offset);

void print_saved_compares(void);

void print_best_saved_maps(int compare, int contig, real threshold, real falloff);

void print_saved_maps(int compare, int contig);

void get_fixed_qtl_weights(QTL_MAP *map);

///** The myriad of QTL_MAP printing routines: All are now in qprint.c */
void print_map_divider(void);
//void print_qtl_map();
//void print_short_qtl_map(); /* FROB */
//void print_short_title();   /* FROB */
//void print_tiny_map();  /* unused and broken */
//
//void print_wiggle_title();     /* no args */
//void print_wiggle_interval();  /* args: QTL_MAP *map; */
//void print_wiggle_map();       /* args: QTL_MAP *map; real threshold, scale; */
//void print_wiggle_genetics();  /* args: GENETICS *genetics; */
//void print_wiggle_left_seq();  /* args: QTL_MAP *map; */
//
//void print_saved_wiggle();
//void print_saved_wiggles();
//void print_saved_wiggle_order();
//void print_peak(); /* arg: WIGGLE_PEAK *peak; */
//
//void print_saved_compares();
//void print_saved_maps();
//void print_best_saved_maps();
//
//void print_test_wiggle_title();    /* no args */
////void print_test_wiggle_interval(); /* args: QTL_MAP *map; */
//void print_test_wiggle_order();    /* args: int wiggle, order; real thresh;*/

//char *genetics_str();  /* args: GENETICS *ptr; bool verbose; */
//char *interval_str();  /* args: int left, right; bool pad; */
//char *left_seq_str();  /* args: QTL_MAP *map; */
//char *dist_str();
//char *units_str();

//void print_trait();   /* arg: int for_how_many_traits; */
//void print_seq();     /* no args */
//void print_old_seq(); /* args: char *seq_string; */

//void get_fixed_qtl_weights(); /* args: QTL_MAP *map; */

#define MAP_DIVIDER \
"=============================================================\n"
#define LONG_MAP_DIVIDER \
"===========================================================================\n"
#define ITER_DIVIDER \
"-------------------------------------------------------------\n"
#define BIG_DIVIDER \
"=========================================================================\n"


/*** A struct to save lists of QTL_MAP structs... ***/
/* This code has not been used in a long time! */

typedef struct sav_map_struct {
    int max_maps, max_intervals;    /* max alloced for */
    int num_maps;            /* number of maps saved */
    real threshold;            /* min log like */
    QTL_MAP **map;            /* array of ptrs to QTL_MAPS */
    QTL_MAP *unused;        /* plus one unused for scratch */
    int unused_checked_out;        /* and we keep track of it */
} SAVE_QTL_MAPS;

SAVE_QTL_MAPS *alloc_saved_maps(int num_maps, int num_intervals, real threshold);

void free_saved_maps(SAVE_QTL_MAPS *p);

QTL_MAP *save_map(QTL_MAP *map_to_save, SAVE_QTL_MAPS *the_maps);

QTL_MAP *get_unused_map(SAVE_QTL_MAPS *the_maps);

void return_unused_map(QTL_MAP *map_to_return, SAVE_QTL_MAPS *the_maps);
//SAVE_QTL_MAPS *alloc_saved_maps();
//QTL_MAP *save_map();
//QTL_MAP *get_unused_map();
//void free_saved_maps();
//void return_unused_map();
//#define maps_saved(saved_maps) (saved_maps->num_maps)
//#define nth_map(saved_maps,i)  (saved_maps->map[i])
//#define best_map(saved_maps)   (saved_maps->map[0])

#endif
