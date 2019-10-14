#ifndef _MAP_INFO_H_
#define _MAP_INFO_H_

/******************************************************************************

 #    #    ##    #####              #    #    #  ######   ####           #    #
 ##  ##   #  #   #    #             #    ##   #  #       #    #          #    #
 # ## #  #    #  #    #             #    # #  #  #####   #    #          ######
 #    #  ######  #####              #    #  # #  #       #    #   ###    #    #
 #    #  #    #  #                  #    #   ##  #       #    #   ###    #    #
 #    #  #    #  #      #######     #    #    #  #        ####    ###    #    #

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#include "mapm.h"
//#include "toplevel.h"

/***** MAP struct - holds all information about a linkage map - maps.c *****/
typedef struct {
    char *map_name;     /* name of MAP */
    int max_loci;
    int num_loci;
    int *locus;         /* [num_loci] (markers in the map) */
    double **rec_frac;     /* [intervals][M/F], recombination fractions */
    int unlink;         /* indicates interval (if any) held unlinked */
    int *fix_interval;  /* [intervals], intervals not converging */
    int sex_specific;
    double log_like;
    bool allow_errors;   /* if TRUE, below things are valid */
    double **error_lod;    /* [num_loci][indiv] - maybe NULL */
    double *error_rate;    /* [num_loci] - NULL if error_lod==NULL */
    bool modified;
} MAP;

/* for map->unlink */
#define MANY_UNLINKED (-2)
#define NONE_UNLINKED (-1)
/* A number 1 or greater indicates the interval to be held unlinked (interval
   0 is the open interval left of the right-most locus). */

/* Functions for MAP struct in maps.c */
void allocate_mapping_data(int num_markers);

void free_mapping_data(int num_markers);

void free_map(MAP *map);

bool clean_map(MAP *map);

void mapcpy(MAP *to, const MAP *from, bool clean_it);

int insert_locus(MAP *map, int position, int locus);

MAP *allocate_map(int maxloci);
//MAP *allocate_map();   /* arg:  max_loci; alloc map for up to max_loci */
//void free_map();       /* arg:  MAP *map; frees map previously allocated */
//void mapcpy();         /* args: MAP *to, *from; bool clean_it; copy from->to */
//int insert_locus();    /* args: map, position, locus; adds locus to map */
//bool clean_map();      /* arg:  map; resets rf's etc in map, returns TRUE */
//void setup_error_matrix(); /* args: MAP *map, int n_indivs, real rate; */

/* Call init_rec_fracs only if you know that map->rec_frac[][] have been 
   initialized to either NOT_FIXED, UNLINK_ME, or a fixed value. */
void init_for_ctm(MAP *map, bool sex, bool errors, bool start);

void init_rec_fracs(MAP *map);

void init_not_fixed(MAP *map);
//void init_rec_fracs(); /* arg:  MAP *m; sets rf, fix_interval, etc for ctm */
//void init_not_fixed(); /* arg:  MAP *m; sets rf, fix_interval for ctm */
//void init_for_ctm();   /* arg:  MAP *m; bool sex, errors, hot_start; similar */




typedef struct {  /* nifty little struct used for the above */
    real like;
    real dist;
    bool zero;
    real net_error;
    real worst_error;
} PLACE;

typedef struct {  /* working data for place cmd and extend_order() */
    int locus;
    bool *excluded;  /* [position-in-order] */
    int num_places, best_pos;
    bool off_end, errors;
    int priority;    /* a random number */
    MAP *best_map;
} PLACEME;


/***** SAVED_LIST Multiple map structure - maps.c *****/

typedef struct {
    int max_maps;      /* used in allocation of structure */
    int max_loci;      /* used in allocation of structure */
    MAP **map_list;    /* list of max_maps pointers to MAPs */
    MAP *extra_map;    /* expendable map-for insertion and deletion */
    int num_maps;      /* number of maps in map_list (<= max_maps) */
    double threshold;     /* highest acceptable likelihood when sorting */
    bool sorted;        /* determines whether this is a sorted list */
} SAVED_LIST;

#define VERY_UNLIKELY -999999.0 /* useful threshold */
#define SORTED     1
#define UNSORTED   0
#define FULL_LIST  (-1)

void write_mapping_data(FILE *fp);

void read_mapping_data(FILE *fp);

SAVED_LIST *allocate_map_list(int maxmaps, int maxloci, bool sortflag, MAP **map);

MAP *get_map_to_bash(SAVED_LIST *list);

void free_map_list(SAVED_LIST *list);

void clean_list(SAVED_LIST *list);

MAP *get_best_map(SAVED_LIST *list);

int insert_map_into_list(SAVED_LIST *list, MAP **map);

void overwrite_map_num(SAVED_LIST *list, MAP **map, int chrom);

void read_map(FILE *fp, MAP *map);

void write_map(FILE *fp, MAP *map);
//SAVED_LIST *allocate_map_list();  /* args: max_maps, max_loci, sorted; */
//MAP *get_map_to_bash(); /* arg: SAVED_LIST *list; gets the extra map to bash */
//int insert_map_into_list(); /* arg: *SAVED_LIST, **MAP; inserts *MAP */
//                             /* sets *MAP to new extra_map for re-use */
//void overwrite_map_num(); /* args: SAVED_LIST *list; MAP **map; int chrom; */
//
//void free_map_list();  /* arg: *SAVED_LIST; frees SAVED_LIST structure */
//void clean_list();     /* arg: *SAVED_LIST; resets values in list for re-use */
//MAP *get_best_map();   /* arg: *SAVED_LIST; returns best map in list */
//
//void read_map(), write_map(); /* args: FILE *fp; */



/***** Global structs used to store 2pt and 3pt data - in map_23.c *****/

typedef struct {
    double lodscore[3];    /* actually, only SEXSPEC and NOSEX are used */
    double theta[3];       /* MALE, FEMALE, and NOSEX */
    bool used;        /* for garbage collection */
} TWO_PT_DATA;

#define UNLINKED_LOD   0.50  /* Because LODs less than this never get stored */
#define UNLINKED_THETA 0.40  /* we should not allow lodbounds<UNLINKED_LOD */

void allocate_two_pt(int num_loci);

void free_two_pt(int num_loci);
//void allocate_two_pt();  /* arg: num_loci; allocs and inits global 2pt data */
//void free_two_pt();      /* arg: num_loci; frees global */
//void expand_two_pt();    /* arg: num_entries; pre-allocates space in global */
#define EXPAND_DEFAULT (-1) /* can be arg to above */

void compute_two_pt(int a, int b);

bool get_two_pt(int a, int b, real *lod, real *theta);

bool get_sex_two_pt(int a, int b, real *lod, real *thetam, real *thetaf);

void allocate_three_pt(int num_total);

void bash_all_three_pt(int num_total);

void free_three_pt(int num_total);

bool restore_triple(int locus1, int locus2, int locus3, real *d1, real *d2, real *d3);

bool three_linked(int *locus, real lodbound, real thetabound, int num_links, bool sex);
//
///* these two calculate on demand WHAT ABOUT HAPS??? */
//void compute_two_pt();  /* args: int locus1, locus2; */
//bool get_two_pt();	/* args: int locus1, locus2; real *lod, *theta; */
//                        /* Returns TRUE if !unlinked */
//bool get_sex_two_pt();	/* args: int locus1, locus2; real *lod, *m, *f; */
//
//void allocate_three_pt(); /* arg: num_loci; allocs and inits global 3pt data */
//void free_three_pt();     /* arg: num_loci; frees global */
//void bash_all_three_pt(); /* arg: num_loci; frees then allocates */
//
//bool restore_triple();
///* args: int locus1, locus2, locus3; real *delta1, *delta2, *delta3;
//   returns TRUE and side-effects deltas if three-pt data exist */
//
////void compute_triple();
/////* args: int locus1, locus2, locus3; real *delta1, *delta2, *delta3;
////   log-likelihoods (0.0 or negative) are filled in. */
//
////void insert_triple();
/////* args: int locus1, locus2, locus3; real delta1, delta2, delta3; saves it */
//
//bool three_linked(); /* args: int *locus; real lod, theta; int nlinks; */
//void compute_3pt();  /* args: *seq; bool sex; real error_rate, *like; *map; */

extern bool two_pt_touched;
extern bool three_pt_touched;

void read_two_pt(FILE *fp);

void write_two_pt(FILE *fp);

void read_three_pt(FILE *fp);

void write_three_pt(FILE *fp);


/**************** Order building code in orders.c ****************/

void get_linkage_group(int *locus, int *num_loci, int *linkage_group, int *group_size, real lodbound, real thetabound);

void setup_haplo_group(int *locus, int num_loci);

bool delete_haplo_groups(int *locus, int num_loci, int *old_locus, int *num_old);

bool force_haplo_sanity(int *locus, bool verbose);

void find_haplo_group(int *locus, int *num_loci, int *haplo_group, int *num_haplo, int *old_obs, int *new_obs);
//
//void get_linkage_group();
//
//void setup_haplo_group();   /* args: loci, num_loci; */
//bool delete_haplo_groups(); /* args: loci, num_loci; TRUE if any deleted */
//void find_haplo_group();    /* args: loci *#loci group *#group obs1 obs2 */
//bool force_haplo_sanity();  /* int *locus; bool verbose; FALSE if changed */

void setup_3pt_data(int *locus, int num_loci, real threshold);
/* args: int locus[], num_loci; real three_pt_threshold;
   Sets the global matrix used by create_order. three_pt_threshold
   should be negative, or if zero, 3pt data are unused. */
//void free_3pt_data();

//bool start_3pt();
/* args: int locus[], *num_loci; real threshold; int order[], *num_ordered; 
   If TRUE returned, all args (but threshold) are side-effected. */

int three_pt_exclusions(int *order, int num_placed, int newmarker, bool *excluded);

/* args: int order[], num_placed, new_locus; bool excluded[];
   returns num_orders;
   Examines three-pt data for inserting new_locus in a order.
   excluded[i], for i=0...num_placed, and *num_orders are side-effected.
   This assumes that excluded[i] is initialized (likely to FALSE) */

bool three_pt_verify(int *locus, int num_loci, int window_size);

/* args: int loci[], num_loci, window_size;
   Examines an order to see if it is compatible with the three-point data.
   Return FALSE if order is excluded. */

void informative_subset(int *locus, int num_loci, int min_infs, real min_theta, bool codom_only, bool haplos, int *subset, int *num);

/* args: locus[], num_loci, min_infs; real min_dist; bool codom_only, haplos;
   int subset[], *num; finds a informative and not-too-closely-spaced set
   of loci. Does not understand sex_specific, nor does it really do anything
   smart with mixed typings. */

void extend_order(MAP *placed, PLACEME **unplaced, int *num_unplaced, real npt_thresh, bool print_anyway);

void find_window(int *order, int num_placed, int new_locus, int *excluded, int min_window, int *start, int *end);

/* int order[], num_placed, new_locus,  excluded[], min_window, *start, *end;
   find the appropriate region in which to try to place a marker, based on
   excluded[], and maybe (someday) informativeness */

void place_locus(MAP *map, int locus, int start, int finish, bool *excluded, PLACE **result, int *best_pos, MAP *best_map, MAP *temp);

/* args: MAP *map; int locus; int start, finish; bool excluded[];
   PLACE *result[]. int *best_pos; MAP *best_map, *temp_map;
   Tries locus into map order between locus #s i=start..finish
   (inclusive), for which !excluded[i]. result[i] is side-effected, where
   like may be 0.0 (best), NO_LIKE (untried) or ZERO_LIKE (if places on
   top of a marker, one side gets this value */

int npt_exclusions(PLACE **place, int num_loci, real like_thresh, real one_err_thresh, real net_err_thresh, bool *excluded, bool *zero_place,
                   int *off_ends, real *error_lod, bool *single_error, real *next_best, real *best_unallowed);
/* args: PLACE *placements[]; int num_loci; real like_threshold;
   real worst_error_threshold, net_error_threshold; 
   bool excluded[], *zero_placement, *placed_off_end; real *error_lod;
   bool *single_error; real *second_best_like, *best_unallowed_like; 
   Examines the likes produced by place_locus(), side-effecting excluded[i], 
   and the many flags. */

/***** Chromosome framework, assignment, and placement stuff - chroms.c *****/

extern SAVED_LIST *chromosome;  /* malloced by allocate_mapping_structs() */
#define chrom2str(x) ((x)>=0 ? chromosome->map_list[x]->map_name : "none")

bool isa_chrom(char *name, int *chrom);
//bool isa_chrom(); /* args: char *name; int *chrom; side-effected if TRUE */
#define num_chromosomes (chromosome->num_maps)
extern int current_chrom; /* set by the sequence command or reset_state() */

bool make_new_chrom(char *name, int *num);

//bool make_new_chrom();   /* args: char *name; returns # */
MAP *get_chrom_frame(int chrom, int *num_loci);

bool framework_marker(int locus);

void set_chrom_frame(int chrom, MAP *new);

void get_chrom_loci(int chrom, int *locus, int which_loci, int *num_loci, int *num_framework);
//void set_chrom_frame();  /* int chrom; char *name; MAP *map; */
///* When this is called, new MUST be equal to get_map_to_bash(chromosome), amd
//   assigned_to(*,chrom) must have be true for all loci in new map. We also
//   assume haplo_sanity is true for old framework, and force it for new one. */
//
//MAP *get_chrom_frame();  /* args: int chrom, *n_loci; n_loci may be 0 or 1 */
//bool framework_marker(); /* int locus; haplos are NOT checked: if use_haplos
//  is on, the chrom frames better only include the first marker in any haplo
//  group */
//
//void get_chrom_loci();  /* args: int chrom, *locus, *num_loci; int which_loci;
//			   copies the chrom's loci into an array. */
#define FRAMEWORK 1 /* do not change these w/o changing get_chrom_loci() */
#define NON_FRAME 2
#define ALL_LOCI  3 /* =FRAMEWORK+NON_FRAME */
#define HAPLOS    4

/* Rule on haplos: if use_haplotypes is on then only the first locus in any
   haplo group is included, unless HAPLOS is specified. if use_haplotypes is,
   off then all fw/assigned/placed markers are included. */

void count_chrom_loci(int chrom, int *n_anchor, int *n_frame, int *n_total, int *n_placed, int *n_unique, int *n_region, bool haplos_too, int *temp);
//void count_chrom_loci();
///* args: int chrom, *n_anchor, *n_frame, *n_total, *n_placed, *n_unique;
//   int *temp; bool haplos_too; temp should be aloced raw.markers long */

typedef struct {
    int chromosome;      /* framework this marker is linked to */
    int linked_to;       /* closest linked marker in that framework */
    real LODscore, theta; /* two-point values with that marker */
    bool modified;        /* needs to be saved */
    int status;
} ASSIGNMENT;

extern ASSIGNMENT **assignment;

/* assignment status values: */
/* Data which must be valid:  CHROM  LOD&LOCUS  NAME      */
#define A_PROBLEM   (-3)  /*  nyet   nyet       conflict  */
#define A_CHANGED   (-2)  /*  nyet   nyet       changed   */
#define A_UNLINKED  (-1)  /*  nyet   nyet       unlinked  */
#define A_UNKNOWN    (0)  /*  nyet   nyet       -         */
#define A_BORDERLINE (1)  /*  da     da         low-lod   */
#define A_ASSIGNED   (2)  /*  da     da         linked    */
#define A_ATTACH     (3)  /*  da     nyet       attached  */
#define A_FRAMEWORK  (4)  /*  da     nyet       framework */
#define A_ANCHOR     (5)  /*  da     nyet       anchor    */
/* note: A_FRAMEWORK means that the locus was in a framework when 
   do_assignment was called, resulting in loss of linkage info for it */

/* all four below silently force_haplo_sanity */
bool assigned(int locus);

int assignment_state(int locus);

int assignment_chrom(int locus);

bool assigned_to(int locus, int chrom);

bool anchor_locus(int locus);
//bool assigned();      /* args: int locus; */
//bool assigned_to();   /* args: int locus, chrom; TRUE if locus is on chrom */
//bool anchor_locus();  /* int locus; */
//int assignment_chrom(); /* args: int locus; */
//int assignment_state(); /* args: int locus; */

/* Macros for yucks, implicitly require force_haplo_sanity */
#define assignment_lod(locus)      (assignment[locus]->LODscore)
#define assignment_locus(locus)    (assignment[locus]->linked_to)
#define assignment_theta(locus)    (assignment[locus]->theta)

bool is_assignable(int locus, int chrom, bool fix_frames);

void assign_this(int locus, int state, int chrom, real lod, real theta, int linked_locus, char *msg);

//void unassign_this(int locus, int state);
void set_chrom_anchors(int chrom, int *locus, int num_loci);

void do_assignments(int *locus, int num_loci, real lod1, real unlinked_lod1, real theta1, real lod2, real unlinked_lod2, real theta2, bool haplo);

//bool is_assignable(); /* args: int locus, chrom; prints msg if FALSE */
///* note that is_assignable is FALSE if locus is insane, haplo-wise */
//
//void assign_this();
///* int locus, state, chrom; real lod, theta; int linked_locus; char *msg;
//   msg pre-empts any other message, only used for A_PROBLEM as yet
//   This can unassign anchor markers, but not reassign them! */
//
////void attach_this();   /* int locus, state; char *msg; just a convenience */
////void unassign_this(); /* int locus, state; char *msg; just a convenience */
//
//void set_chrom_anchors(); /* int chrom, *locus, num_loci */
///* The only right way to assign() a locus as A_ANCHOR. Will reassign if need
//   be and it can. Do not use is_assignable() beforehand, as that never lets
//   one assign an anchor. */
//
//void do_assignments(); /* args: int *locus, num_loci;
//   real lod1, unlinked_lod1, theta1, lod2, unlinked_lod2, theta2; bool haplo;
//   Does the real work - is_assignable must have been verified */

typedef struct {
    int chromosome;
    int num_intervals;
    int interval[MAX_STORED_PLACES];
    real like_ratio[MAX_STORED_PLACES];
    real distance[MAX_STORED_PLACES];
    real errors[MAX_STORED_PLACES];
    real threshold;
    real error_lod;
    bool single_error;
    bool modified;
    int status;
} PLACEMENT;

extern PLACEMENT **placement;

/* placement status values: */
/* Data which must be valid:  LIKES  ERRORS  NAME      */
#define M_PROBLEM   (-2)  /*  nyet   nyet    problem   */
#define M_FRAMEWORK (-1)  /*  nyet   nyet    framework */
#define M_UNKNOWN    (0)  /*  nyet   nyet    -         */
#define M_ERROR      (1)  /*  da     da      errors?   */
#define M_OFFEND     (2)  /*  da     da      off-end   */
#define M_REGION     (3)  /*  da     da      region    */
#define M_UNIQUE     (4)  /*  da     da      unique    */
#define M_ZERO       (5)  /*  da     da      zero      */

#define NO_CHROM  (-1)
#define ANY_CHROM (-2) /* flag for is_assignable() */
#define NEW_CHROM (-3) /* flag for set_chrom_frame */
#define NO_LOCUS  (-1)
#define NO_LOD    0.0
#define NO_THETA  0.4996  /* HMM_MAX_THETA 0.4995 */
#define NO_LIKE   (-999.98765)  /* hope that's obscure */
#define ZERO_LIKE (-999.00123)  /* hope that's obscure too */
#define NO_ERRORS NO_LIKE       /* for ->error_lod */

#define FRAMEWORK_MARKER (-1)

bool placed_locus(int locus);

bool placement_state(int locus);

bool is_placeable(int locus, int chrom);

int place_this(int locus, int chrom, PLACE **place, real like_thresh, real one_err_thresh, real net_err_thresh, bool *excluded);

void unplace_this(int locus, int chrom, int status, bool verbose);

int best_placement(int locus);

int second_best_placement(int locus, real *like);

void bash_mapping_data(int *changed, int num_changed);

void allocate_mapping_data(int num_markers);

void free_mapping_data(int num_markers);

void write_mapping_data(FILE *fp);

void read_mapping_data(FILE *fp);

void print_ps_map(FILE *fp, MAP *map);

void print_ps_chrom(FILE *fp, int chrom);

void print_all_ps_chroms(FILE *fp);
//bool placed_locus();   /* args: int locus; forces haplo_sanity w/no msg */
//int placement_state(); /* args: int locus; forces haplo_sanity w/no msg */
//
//bool is_placeable(); /* args: int locus, new_chrom; */
//  /* prints msg if FALSE, and it WILL be FALSE if !haplo_sanity */
//
//int place_this();
///* args: int locus, chrom; PLACE **placements; real like_threshold;
//   real worst_error_threshold, net_error_threshold; bool excluded[];
//   sets placement struct and returns #intervals within threshold.
//   uses npt_exclusions() to do its real work. Is chatty. assumes
//   haplo_sanity */
//void unplace_this(); /* args: int locus, chrom; bool problem, print_msg; */
//
//int  best_placement();        /* args: int locus; */
//int  second_best_placement(); /* args: int locus; real *like; */
//
//void allocate_mapping_data(), free_mapping_data(); /* args: num_markers; */
//void write_mapping_data(),    read_mapping_data(); /* args: FILE *fp; */
//void bash_mapping_data(); /*args: int changed_markers[], num_changed; */
//
//void print_ps_map(); /* args: fp, map; in ps_maps.c */
//void print_ps_chrom(); /* args: fp, chrom; in ps_maps.c */
//void print_all_ps_chroms();  /* args: fp; in ps_maps.c */


/***************** Misc Saved Mapping Info *****************/
/* support is all now in map_123.c */

extern char **note;      /* [locus]=>string, Note about a particular marker */
extern bool *modified;   /* [locus] FALSE if no, TRUE if yes */
extern real *error_rate; /* [locus] apriori, see use_error_rate, state.c */
extern int *my_group;   /* [locus] For group command, NO_GROUP if unknown */
extern int num_groups;
#define NO_GROUP (-1)

extern int *haplo_first, *haplo_next;  /* [locus] */
extern int *order_first, *unorder_first; /* [order]=> locus# or NO_LOCUS */
extern int *order_next;  /* [locus]=>locus or NO_LOCUS, for order or unorder */
extern int num_orders;
#define haplotyped(locus)  (haplo_first[locus]!=NO_LOCUS)
#define haplotype_subordinate(n) (haplo_first[n]!=NO_LOCUS &&haplo_first[n]!=n)

extern int *class;        /* [locus] class of each marker */
extern char **class_name;  /* [NUM_CLASSES][NAME_LENGTH] name of each class */
#define NO_CLASS     0
#define NUM_CLASSES  11

bool isa_class(char *name, int *num);

bool make_new_class(char *name, char **why_not);

void print_class_names(void);

void bash_order_info(int *changed, int num_changed);

void allocate_order_data(int num_markers);

void free_order_data(int num_markers);

void write_order_data(FILE *fp);

void read_order_data(FILE *fp);
//
//bool isa_class();	  /* args: char *name; int *num; */
//bool make_new_class();	  /* args: char *name; char **why_not; */
//void print_class_names(); /* lets print() auto_wrap, no nl on end */
//
//void allocate_order_data(); /* args: num_markers; */
//void free_order_data();     /* args: num_markers; */
//void write_order_data(); /* args: FILE *fp; */
//void read_order_data();  /* args: FILE *fp; */
//void bash_order_info();  /* args: int changed_marker[], num_changed; */


/*************** Haldane/Kosambi Map Functions - now in maps.c ***************/
typedef struct {
    char name[12];

    double (*add)(real, real);

    double (*apportion)(bool, real, real, real);

    double (*rec_to_dist)(real);

    double (*dist_to_rec)(real);

    double (*d_to_r_deriv)(real);
} MAP_FUNCTION;

#define TEN_CM      0.0906346 /* haldane, that is */
#define THIRTY_CM   0.2255942 /* ditto */
#define FIFTY_CM    0.3160603 /* ditto */
#define EIGHTY_CM   0.3990517 /* ditto */
#define ONE_CM      0.00900663
#define UNLINKED_CM 0.4995441 /* ditto 350 cM */
#define cm_dist(theta) ((*mapfunction->rec_to_dist)(theta)*100.0)

extern MAP_FUNCTION *mapfunction;
extern MAP_FUNCTION maps[2];
extern int num_map_functions;

void map_func(int mapnum);
//void map_func(); /* args: int num; sets map function */
#define map_func_name() (mapfunction->name)
#define map_func_num()  (mapfunction->name[0]=='H' ? HALDANE:KOSAMBI)
#define HALDANE 0
#define KOSAMBI 1

void new_print_placements(MAP *map, PLACEME **placed, int num_loci);

#endif
