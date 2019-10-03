/******************************************************************************

 #    #    ##    #####   #    #          #    #
 ##  ##   #  #   #    #  ##  ##          #    #
 # ## #  #    #  #    #  # ## #          ######
 #    #  ######  #####   #    #   ###    #    #
 #    #  #    #  #       #    #   ###    #    #
 #    #  #    #  #       #    #   ###    #    #

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

/***************** Get the common #includes first...  *****************/

#define INC_TABLE /* need TABLE to compile toplevel.h and STATUS_CONTEXT */

#ifndef command
#define command void /* needed if INC_SHELL was not defined */
#endif

#include "system.h"


/***************** MAPM Constants *****************/
/* you should be able to change these and recompile */

#ifndef _BIG_DATASETS
#define MAX_LOCI          1000  /* max #loci in file, NOT used to malloc */
#define MAX_CHROM_LOCI    100   /* is #loci per frame, used only in map_2.c */
#define MAX_MAP_LOCI      101   /* max #loci in a MAP on any call to CTM */
#define MAX_SEQ_LOCI      1001  /* in words, after expanding, #loci+1 */
#define MAX_SEQ_TOKENS    1100  /* in words, after expanding, #loci+100 */
#define MAX_SEQ_LEN       11000 /* in chars, ~10*MAX_SEQ_TOKENS */
#else
#define MAX_LOCI          5000
#define MAX_CHROM_LOCI    500   
#define MAX_MAP_LOCI      501   
#define MAX_SEQ_LOCI      5001
#define MAX_SEQ_TOKENS    5100  
#define MAX_SEQ_LEN       51000
#endif

#define F2VERSION         3     /* data file format version */
#define MAX_CHROMOSOMES   50    /* fairly obvious */
#define MAX_HISTORY_SEQS  100   /* this table is a fixed size */
#define MAX_NAMED_SEQS	  20    /* actually, this table expands as needed */
#define MAX_STORED_PLACES 20    /* for placement struct */
#define PLACEMENT_THRESHOLD (-5.0)  /* like to store placement info at all */

#define INIT_2PT_SIZE(n)  min(max(5000,((n)*(n))/10),50000)
/* at least 5000, no more than 50,000 to start */
#define INIT_3PT_SIZE(n)  ((n)*25)  /* avg num 3pt orders per raw locus */
#define INIT_3PT_GROUP 100 /* 3pt exclusion mat size, =n^3 chars - orders.c */

#define NAME_LEN 8  /* changed since V2 */
/* do not change this, as the printed name size is assumed to be 9 chars in
   many places in mapm, and use_haplotypes may add a trailing '+' to names */

#define DEFAULT_ERROR_RATE 0.01
#define ZERO_PLACE (-0.0099)   /* virtually zero likelihood for place_locus */
#define ZERO_DIST  (0.0005)    /* virtually zero distance for place_locus */
#define MAX_NOTE_LEN      80   /* in chars */

/* Extensions for MAPMAKER files... */
#ifndef _SYS_DOS
#define PS_EXT    ".ps"
#define RAW_EXT   ".raw"
#define PREP_EXT  ".prep"
#define TEMP_EXT  ".temp"
#define DATA_EXT  ".data"
#define DATA_OLD  ".xdata"
#define MAPS_EXT  ".maps"
#define MAPS_OLD  ".xmaps"
#define TWO_EXT   ".2pt"
#define TWO_OLD   ".x2pt"
#define THREE_EXT ".3pt"
#define THREE_OLD ".x3pt"
#define TRAIT_EXT ".traits"
#define TRAIT_OLD ".xtraits"
#else
#define PS_EXT    ".ps"
#define RAW_EXT   ".raw"
#define PREP_EXT  ".pre"
#define TEMP_EXT  ".tmp"
#define DATA_EXT  ".dat"
#define DATA_OLD  ".xda"
#define MAPS_EXT  ".map"
#define MAPS_OLD  ".xma"
#define TWO_EXT   ".2pt"
#define TWO_OLD   ".x2p"
#define THREE_EXT ".3pt"
#define THREE_OLD ".x3p"
#define TRAIT_EXT ".tra"
#define TRAIT_OLD ".xtr"
#endif


/****************************** Global Stuff *******************************/

/* Lowlevel constant definitions used in many places throughout mapmaker */
#define	MALE		0 
#define FEMALE		1 
#define NOSEX		2
#define SEXSPEC		0   /* for two-point lod-score */
#define for_sexes(s)	for(s=MALE;s<=FEMALE;s++)

/* Signals */
#define BADSEQ      USER_MESSAGE(2)
#define PREPERROR   USER_MESSAGE(3)
#define BADDATA	    USER_MESSAGE(4)

extern int   BADDATA_line_num;
extern char *BADDATA_text;
extern char *BADDATA_reason;

#define SELF_DELIMITING_TOKENS  "<>[]{}~|!=" /* mostly for seq parser */



/****************************** Raw Data Types *******************************/

#define CEPH 0
#define PH_KNOWN 0         /* subtypes of CEPH - for raw.use_number? */
#define PH_UNKNOWN 1
#define F2 2
#define F2_INTERCROSS 2    /* subtypes of data type F2 - for raw.cross_type */
#define F2_BACKCROSS  3
#define RI_SIB  4
#define RI_SELF 5
#define F3_SELF 6

#define NO_DATA (-1)       /* fake data types used as args to mapm_ready() */
#define ANY_DATA (-2)
#define MAYBE_DATA (-3)
#define RI_NOT (-4)
#define F2_DATA 2
#define CEPH_DATA 0

/* for raw.data_type, use only CEPH, F2, or NO_DATA */
#define NUM_DATA_TYPES 4  /* used in old_ctm select_procs */

/* Things used by the data readers in reader.c */
void getdataln();
void baddata();
void do_save_data();
void do_unload_data(); 
void do_read_data();


/***************** MAPM Status Variables  ******************/
/* Variables that the user controls - see reset_state(), state_init() */
void reset_state();  /* set all vars to sane values */
void undo_state();   /* free a little stuff */

/* NOTE! Keep the declarations of vars in the same order here, below in 
   reset_state(), in read/write_state, and in mapm.h. */

/* general */
extern bool print_names;
extern real tolerance;
extern int  units;
#define RECFRACS 0
#define CENTIMORGANS 1
extern bool auto_save;
extern int  print_maps;  /* should do something with this */

/* two-point */
extern real default_lod;
extern real default_theta;
extern bool use_haplotypes;

/* three-point */
extern bool use_three_pt;
extern real triplet_lod;
extern real triplet_theta;
extern int  triplet_num_links;
extern real three_pt_threshold;
extern int  three_pt_window;
extern real triplet_error_rate;
#define LOCUS_ERROR_RATE OBSCURE_REAL /* for above */

/* Order Maker */
extern real npt_threshold;
extern real npt_first_threshold;
extern int  npt_window;
extern int  npt_min_indivs; /* infomativeness criteria */
extern bool npt_codominant;
extern real npt_min_theta;
extern bool print_all_maps;

/* F2 only */
extern bool fake_maps;  /* just for debugging */
extern bool segregation_distortion; /* Obsolete??? */
extern bool use_hmm;
extern bool use_error_rate;
extern real error_lod_thresh;
extern real error_net_thresh;
extern real error_single_thresh;
extern bool print_all_errors;

/* CEPH only */
extern int  sex_specific; 	
extern int  compress_DNA;
extern int  print_problem_size;
extern long max_problem_size;

/* Obsolete (Old CTM?) */
extern int  time_stamping;
extern real inner_tolerance;
extern int  inner_loop;
extern int  convergence_rule;
extern int  original_markers;
extern real startrecombs;
extern int  print_dots;

/* Status Context */
typedef struct {
    int         sex_specific;
    int         compress_DNA;
    int         use_number;
    long        max_problem_size;
    int         seq_history_num;
    TABLE       *sequence_history;
    TABLE       *named_sequences;
} STATUS_CONTEXT;

extern STATUS_CONTEXT **context;
extern int num_contexts;
extern int active_context;
#define the_context          (context[active_context])
#define the_seq_history_num  (context[active_context]->seq_history_num)
#define the_sequence_history (context[active_context]->sequence_history)
#define the_named_sequences  (context[active_context]->named_sequences)

void allocate_context();
void free_context();
bool change_context();
bool create_new_context();
void read_status();  /* arg: FILE *fp; */
void write_status(); /* arg: FILE *fp; */

/* Other State Vars of note
MAP_FUNCTION *mapfunction (map_info.c, maps.c)
raw.data.ceph.use_number
SEQ_NODE *seq; char *seq_string; set by set_current_seq
*/ 


/***************** MAPMAKER Commands *****************/

#define RAWDATA_CMD  1
#define SEQUENCE_CMD 2
#define TWOPT_CMD    3
#define MAPPING_CMD  4
#define OPTION_CMD   5
#define NAMES_CMD    6
#define DATAVIEW_CMD 7
#define SYSTEM_CMD   8
#define ALIAS_CMD    9
#define WIZARD_CMD   10
#define HELPONLY_CMD 11
#define AUTO_CMD     12
#define B_TOPIC 2


/**** in state.c ****/ 
/* This is a rat's nest. Keep the order here identical to the declarations
   above, which should be identical to the order of declarations in
   state.c, which should be identical to the order of settings in
   reset_state(), which should be identical to the order of reading and
   writing in read_status() and write_status(), which should be identical
   to the order of the set_* functions in state.c. If you don't keep
   these consistent, I will kill you. */

/* general */
command set_print_names();
command set_tolerance();
command set_units();
command set_cm_func();
command set_autosave();
command set_more_mode();

/* two point */
command set_default_linkage();

/* three point */
command set_use_3pt();
command set_3pt_linkage();
command set_3pt_threshold();
command set_3pt_errors();

/* order maker */
command set_npt_threshold();
command set_inf_threshold();
command set_print_all_maps();

/* F2 and error checker stuff */
command set_fake_maps(); /* wizard */
command set_use_error_rate();
command set_error_lod_thresh();

/**** in sys_cmds.c ****/
command new_load_data();
command new_save_data();
command new_prepare();
command set_error_rate();
command make_note();		/* note */
command set_age();

/* sequence commands now in sys_cmds.c */
command sequence();
command expand_sequence();
command history();
command let();
command let_expanding();
command names();      		
command forget();     		/* forget name? */

command download();   /* downloads data from map_base for Mapmaker loading */
command save_to_database();/* sends map data (frame/assign/place) back to the database */
command import();
command export();


/**** in two_cmds.c ****/
command two_point();
command three_point();
command forget_three_point();
command haplotype();
command unhaplotype();
command group();
command order_maker();
command greedy();
command list_loci();
command list_mapping();
command list_haplotypes();

command biglods();
command lodtable();
command near_locus();
command near_chrom();
command pairwise();
command place();
command suggest_subset();

/* in cmds.c */
command make_map();		/* map */
command draw_map();             /* map */
command genotypes();		/* map */
command compare();
command try();
command likely(); 		/* likelihood */
command permsex();    		/* permute sexes */
command use();        		/* use */
command ripple();

command set_framework();  /* makes a map, stores it in chromosome list */
command set_anchors();
command print_chromosomes();  /* prints SAVED_LIST of chromosomes */
command print_short_chroms(); /* for F2 data, prints chromosomes in raw F2 form */
command attach();
command assign();
command unassign();
command save();
command meiosis_counter();
command chrom_pics();
command family_cross_count();
command delete_locus();
command append_locus();
command edit_sequence();
command new_insert();
command new_delete();
command new_append();
command show_seq_history(); /* previous sequences */
command set_class();
command make_classes();
void init_class_names(); /* could become an interactive command */

/* map_func.c */
command cm_func();

/* in auto.c */
command assign(), unassign(), attach(), chrom(), show_chrom(), list_loci();
command place_together(), subset(), place();
command make_chromosome();
command list_chroms();
command list_assignments();
command batch_place(), iter();
command draw_chromosome();
command draw_all_chromosomes();


void seg_dist();

bool mapm_save_on_exit(); /* now in main.c */

command prepdat();    		
command translate();


/* in sequence.c */
bool valid_new_name(char *str);

/* moved set_age and make_note */


/***************** Init Functions For Specific .c Files ******************/
void npt_cmds_init();
void map_init();
void sequence_init();
void state_init();
void data_init();

/***************** Other MAPM Include Files  ******************/

#include "lowlevel.h"
#include "map_info.h"
#include "toplevel.h"


/*
#ifdef INC_PREPDATA
#include "prepdata.h"
#endif
*/

