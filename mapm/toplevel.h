#ifndef _TOPLEVEL_H_
#define _TOPLEVEL_H_

/******************************************************************************

  #####   ####   #####   #       ######  #    #  ######  #               #    #
    #    #    #  #    #  #       #       #    #  #       #               #    #
    #    #    #  #    #  #       #####   #    #  #####   #               ######
    #    #    #  #####   #       #       #    #  #       #        ###    #    #
    #    #    #  #       #       #        #  #   #       #        ###    #    #
    #     ####   #       ######  ######    ##    ######  ######   ###    #    #

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

/********************** Useful things for MAPM Commands **********************/

void mapm_ready(int data_type, int min_seq_loci, bool permable_seq, int *seq_loci);
///* in toplevel.c */
//void mapm_ready(); /* args: int data_type_needed;
//		            int min_loci_in_seq; bool permutable_seq;
//			    int *loci_in_seq; */

/* Data_type_needed may be one of CEPH, F2, ANY_DATA, NO_DATA, or
   MAYBE_DATA.  Min_loci_in_seq should be >0 if a seq is required, and
   should be the minimum number of loci which MUST be in the sequence. If
   it's non-zero, then the sequence is verified and the global seq is
   set, otherwise seq may be bogus. If min_loci_in_seq==MAYBE_SEQ (-1)
   then the size of the sequence is not tested, although seq IS set!

   Permutable_seq should be TRUE if it must be, FALSE if it should not
   be (e.g. if it is used as a list-o-loci), or MAYBE if no test should
   be performed. This arg is ignored if min_loci_in_seq==0.  As a
   convenient side-effect, the number of loci in the sequence is stored
   in *loci_in_seq, if it is non-null and min_loci_in_seq is non-zero. */
#define MAYBE_SEQ (-1)

#define PERM_SEQ          1  /* TRUE  */
#define LIST_SEQ          0  /* FALSE tests assignment haplos dups */
#define MAYBE_PERM      (-1) /* MAYBE */
#define UNCRUNCHED_LIST (-2)
#define ONE_ORDER       (-3)

#define NEED_SEQ_OR_ARGS "need loci as arguments OR sequence must be set"

bool crunch_locus_list(int *locus, int *num_loci, bool verbose, bool check_assignments, bool in_sequence);
///* this is also now all in main.c */
//bool crunch_locus_list();
///* args: locus[], *num_loci; int verbosity,  check_assignments, in_sequence;
//   side-effects locus[], *num_loci; and, if verbose, may print a msg.
//   Deletes duplicates in the list, multiple loci from the same haplo
//   group, and renumbers haploed loci to their primary. If
//   verbosity==ORDER_ERROR and if duplicates existed, error() is called.
//   If verbosity==CRUCH_WARNINGS, warnings are printed, otherwise it is
//   silent. If check_assignments and the chrom is set, then chrom
//   assignments of all loci to the current_chrom are verified. */

#define CRUNCH_WARNINGS TRUE  /* for verbosity */
#define ORDER_ERRORS    MAYBE
#define SILENTLY        FALSE
#define CHECK_CHROM     TRUE  /* for check_assignments */
#define ANY_CHROMS      FALSE
#define IN_SEQ          TRUE  /* for in_sequence */
#define IN_ARGS         FALSE

int get_chrom_arg(bool allow_no_chrom);
//int get_chrom_arg(); /* args: bool allow_no_chrom_answer; */
/* tries using get_one_arg */

//bool get_markers();
//bool get_intervals();
//bool get_reals();
bool input_dist(real *dist);
//bool input_dist();

/* stuff in reader.c */
//void getdataln(FILE *fp);
//void baddata(char *reason);
int data_loaded(void);
char *data_info(bool add_nums);
FILE *start_save_to_file(char *name, char *ext, char *type, bool *exists);
void finish_save_to_file(char *name, char *oldext, bool exists);
void do_load_data(FILE *fp, char *filename, bool israw);
//void do_unload_data(void);
//void do_save_data(char *base_name, bool save_genos_too);
int read_data_file_header(FILE *fp, char *filename);
int read_raw_file_header(FILE *fp, char *filename, char *symbol);
bool read_magic_number(FILE *fp, char *type);
void write_magic_number(FILE *fp, char *type);
int new_magic_number(void);
void allocate_f2_data(int n_loci, int n_indivs);
void free_f2_data(int n_loci);
void free_traits(void);
void new_read_f2_data(FILE *fp);
void new_read_f2_locus(FILE *fp, int locus_num);
void write_f2_data(FILE *fp);
void new_read_f2_raw(FILE *fp, char *symbols);
void read_raw_f2_locus(FILE *fp, int locus_num, char *symbol);
void read_raw_trait(FILE *fp, int trait_num);
//int symbol_value(int chr, char *symb);
void write_traits(FILE *fp);
void add_to_seg_dist(char c, int locus);
void scale_seg_dist(int locus);

void mapm_data_info(FILE *fp);

void try_to_load (FILE *fp, char *name, bool prev_data, bool raw);

//void read_data();    /* args: FILE *fp; char *file_name; */
//bool data_loaded();  /* no args */
//char *data_info();   /* args: bool add_numbers_to_info; returns temp_string */
//void do_load_data();
//void do_unload_data();
//void allocate_f2_data();
//void free_f2_data();
//void mapm_data_info(); /* for mapm_update_hook() */
//void try_to_load();  /* FILE *fp; char *name; bool prev_data, is_raw;
//			in sys_cmds.c - where should I put this decl? */



/***************************** MAPM's Map Makers  ****************************/

void quick_two_pt(int locus0, int locus1, TWO_PT_DATA *two_pt, bool sex /* do both with and without sex spec */); /* ctm.c */
void f2_quick_two_pt(int loc1, int loc2, TWO_PT_DATA *two_pt, bool sexflag); /* in quick23.c, used by the above */

/* in ctm.c */
//void converge_to_map(MAP *map);
void f2_genotype(int locus, bool haplo, int *observation);
bool merge_genotypes(int locus, int *observation, int *new_observation);
int f2_count_infs(int *num_dom, int *num_het, int *observation);

//void converge_to_map(); /* args: MAP *map */
//void f2_genotype(); /* args: int locus; bool haplo; int *observation[#indiv] */
//bool merge_genotypes(); /* args: int new_locus, old_obs[], new_obs[] */
//  /* returns TRUE if no obligate recs, WARNING: will side-effect new_obs even
//     if FALSE is returned! */
//int  f2_count_infs(); /* args: *n_dominant, *n_het, obs[]; return #infs */
extern int *observations; /* [num_indivs] useful as a temp to the above */

/* this stuff is used by the new CTM and by a few other things */
#define OBS_MISSING 0
#define OBS_A 1
#define OBS_H 2
#define OBS_B 3
#define OBS_NOT_A 4  /* C */
#define OBS_NOT_B 5  /* D */
#define fully_inf(obs) ((obs)==OBS_A || (obs)==OBS_H || (obs)==OBS_B)
extern int obligate_recs[6][6];

void allocate_hmm_temps (int total_loci, int num_indivs, int cross_type);
void free_hmm_temps (int total_loci, int num_indivs, int cross_type);
//void allocate_hmm_temps(), free_hmm_temps();
/* args: int n_loci, n_indivs, cross_type; */



/******************************* MAPM Sequences ******************************/

/* The structure which holds compiled sequences. Nothing but the code in
   sequence.c should ever care what this RELLY looks like... */

typedef struct seq_struct {
    flag type;                /* one of the #defines below */
    union {
        struct {
	    int num;
	  } locus;
	struct {
	    int first, last;
	  } locus_range;
	struct {
	    struct seq_struct *contents;
	  } megalocus;
	struct {
	    struct seq_struct *contents;
	    bool flip; /* permutation state */
	  } invertable_set;
	struct {
	    struct seq_struct *contents;
	    struct seq_struct *reordered_contents; /* list using next_item */
	  } unordered_set;
    } node;
    struct seq_struct *next,*prev; /* make a bidirectionally linked list */
    real fixed_distance;           /* not always valid, usually NOT_FIXED */
    bool unlinked_now;             /* state */
    int insert_pos;                /* state for elements of an unordered_set */
    struct seq_struct *next_item;  /* also state: a differently ordered list */
} SEQ_NODE;

#define NONE		0
#define SEQLOCUS 	1
#define LOCUS_RANGE	2
#define MEGALOCUS 	3
#define INVERTABLE_SET 	4
#define UNORDERED_SET 	5

#define NOT_FIXED OBSCURE_REAL
#define UNLINK_ME OBSCURE_REAL2

void allocate_seq_stuff(int n_loci);
void free_seq_stuff(void);
//void allocate_seq_stuff(); /* args: int max_loci; call after loading data */
//void free_seq_stuff();     /* args: int max_loci; */

/* BADSEQ message */
extern char *BADSEQ_errmsg;
extern int BADSEQ_errpos;     /* index of the erroneous token in seq string */
void badseq(char *errmsg);
void print_badseq(void);
//void badseq();                /* args: char *error_message; sends BADSEQ */
//void print_badseq();	      /* no args */

/* Fun things you can do to compiled sequences... */
int count_loci (SEQ_NODE *p);
//int count_loci();	/* args: SEQ_NODE *seq; count loci in seq */
bool has_fixed_dists(SEQ_NODE *p);
//bool has_fixed_dists();	/* args: SEQ_NODE *seq; if any fixed_dists in seq */
bool unpermutable(SEQ_NODE *p);
//bool unpermutable();	/* args: SEQ_NODE *seq; TRUE if seq is unpermutable */
#define permutable(seq) !unpermutable(seq)

/* These two global arrays are ONLY for crunch_locus_list() in main.c */
extern int *seq_locus; /* [MAX_SEQ_LOCI] */
extern int *use_locus; /* [raw.num_markers] */

/* Similarly, this stuff is ONLY for use in seq_cmds.c */
extern char **seq_tokens;       /* [MAX_SEQ_TOKENS] */
extern char *new_seq, *the_seq; /* [MAX_SEQ_LEN] */

/* FIX these should move to sequence.c/toplevel.h - is now unused? */
//void set_seq_from_list();


/**** Below are the real command level interfaces to sequences ****/

/* MAPMAKER's current sequence. set_current_sequence() side-effects the global
   variables seq and seq_string. If an error occurs, the BADSEQ message is 
   sent and seq and seq_string are NOT changed. Note that set_current_sequence
   now uses delete_comments() and expand_sequence(), and it expands locus 
   names, locus range settings, etc. */

void set_current_seq();   /* args: char *str; bool expanded; str=NULL unsets */
extern SEQ_NODE *seq;     /* global: the current sequence */
extern char *seq_string;  /* ditto */
void check_current_seq(); /* for use by mapm_ready */
void parse_locus_args();  /* args: int **loci, *num_loci; may send error */
void make_compare_seq();  /* int *locus, num_loci, start_set, set_size; */
//void expand_seq_names();  /* char *str */

bool alloc_list_of_all_loci(); /* args: SEQ_NODE *seq; int **locus,*num_loci;*/
/* For convenience, it's like array(,count_loci(),int) then
   get_list_of_all_loci() followed by crunch_locus_list(), testing haplos
   and duplicates (but NOT chrom assignments, and it is set to print
   warnings). Returns TRUE if num_loci>0, FALSE (w/o allocating)
   otherwise */

void get_list_of_all_loci();
/* args: SEQ_NODE *seq; int *locus, *num_loci, max_loci;
   sets locus[...] and side-effects *num_loci to list all loci in seq
   max_loci is the malloced size of the list, we get a CRASH if it
   is exceeded. The array needs to be big enough to hold ALL loci in the seq, 
   not just haplo groups, one of each duplicate, etc */
#define sort_loci(loci,num_loci) inv_isort(loci,num_loci)

void get_one_order(); /* args: SEQ_NODE *seq; MAP *map; */
/* Like get_list_of_all_loci. Map->max_loci must be big enough. */

/* iterator macro for_all_orders(); args: SEQ_NODE *seq; MAP *map; (not **map)
   Iterates over all orders specified by seq, side-effecting map. It is
   OK for the ptr specified by map to change between iterations (this is
   a wierd thing which happens because this is a macro, not a function). 
   However, also because this is a macro, do not give it arguments which
   side-effect something (eg: for_all_orders(seq,maps[i++]) is a no-no).
   Also, do not try to put one call to for_all_orders() inside another, 
   nor put a call to unpermutable() or enumerate_loci() inside a 
   for_all_orders() loop! */

#define get_map_order(S,M) \
     get_order(S,(M)->locus,(M)->rec_frac,&(M)->num_loci,(M)->max_loci)
#define for_all_orders(S,M) \
for (Oagain=TRUE, reset_seq(S,TRUE); \
     Oagain && clean_map(M) && get_map_order(S,M); \
     Oagain=perm_seq(S,FALSE,FALSE))

extern bool Oagain;     /* internal use only */
void reset_seq();	/* internal only! */
bool perm_seq();	/* internal only! */
bool get_order();	/* internal only! */

/* iterator macro for_all_pairs()  args: int *loci, num_loci, locus1, locus2;
   Loops over all sq(num_loci)/2 unique pairs of the loci specified. 
   locus1 and locus2 are set to the aprropriate values. Never use a 
   for_all_pairs() loop inside another! */

#define for_all_locus_pairs(L,num,a,b) 	\
 for (Pi=0; Pi<num-1; Pi++) for(Pj=Pi+1,a=L[Pi],b=L[Pj]; Pj<num; Pj++,b=L[Pj])
extern int Pi, Pj;  /* internal only! */

/* iterator macro for_all_3pt_seqs() args: int *loci, num_loci; SEQ_NODE *seq;
   Sets seq equal to the three point sequences for an ambiguous order (e.g. 
   the equivalent of "{a b c}") for all unique triples of loci, a, b, and c 
   in the array of locus numbers. Never layer one of these loops inside 
   another. Calling this with num_loci<3 causes a CRASH. */

#define for_all_3pt_seqs(loci,num,seq) 			\
  for (seq=Tinit(&Tagain,&Ti,&Tj,&Tk,loci,num); Tagain; \
       Tagain=Tnext(&Ti,&Tj,&Tk,loci,num,seq))
extern int Ti, Tj, Tk;  /* internal use only! */
extern bool Tagain;     /* internal use only! */
SEQ_NODE *Tinit();      /* internal use only! */
bool Tnext();           /* internal use only! */



/************************ Names handling routines *************************/

bool is_a_locus(); /* args: char *token; int *num; char **errmsg; */
/* Gets the locus number (0...max-1, inclusive) for the locus
specified by str. Token MUST be a despace()ed and filter()ed SINGLE
TOKEN! If token is not a valid locus name or number, FALSE is
returned, *errmsg is set, and *num is left undefined. Otherwise TRUE
is returned, *num is set, and *errmsg is left undefined. This routine
checks that token matches the locus name without matching any sequence
names, etc. and thus, abbreviations work in a semantically nice way. */

bool is_a_sequence(); /* args: char *token, **str; char **errmsg; */
/* If the token is a reference to a MAPMAKER sequence (either a named
sequence or a history reference) then TRUE is returned and *str is
side-effected. Otherwise FALSE is returned and *errmsg is
side-effected. Abbreviated names are handled in the same way that
is_a_locus() handles them. Token must be despace()ed and filter()ed! */

//bool is_a_special_sequence(); /* args: char *token, **str; char **errmsg; */
///* Like above, and called by is_a_sequence(). The only difference is that
//   TRUE can be returned with errmsg set (!nullstr()), meaning that the name
//   WOULD be valid but doesn't happen to be, given MAPM's state. */

void print_special_sequences(); /* no args */
void print_user_sequences();    /* no args - print user names */
void print_history_seqs();      /* args: int max_num_to_print; */

bool name_sequence();      /* args: char *name, *seq, **errmsg; */
bool unname_sequence();    /* args: char *name, **errmsg; */
/* Thes both return TRUE if successful, and otherwise return FALSE and 
   side-effect *errmsg. */

void add_to_seq_history (char *seq, bool is_next_entry);
//void add_to_seq_history(); /* args: char *seq; */

/* bool valid_name(); */
/* args: char *str; TRUE if str is a valid name token*/
/* This does NOT necessarily mean that str is defined*/

void tokenize_seq(); /* args: char *seq, **token; int *num_tokens; */
  /* parses seq into into an array of words (using SELF_DELIMITING) */

void untokenize_seq(); /* args: char *str, **token; int *num_tokens; */
  /* unparses tokens side-effecting str, but does not set_current_seq() */

void print_sequence(); /* no args */


/**************************** Stuff in print.c ****************************/

char *rf2str();   /* arg: rec_frac; returns cM if print_cm, else itself */
char *loc2str();  /* arg: locus; returns name if print_names, or number */
char *locs2str(); /* args: locus1, locus2; */
char *rag();      /* args str; wrap around the above to get ragged result */
char *locname();  /* int locus; bool print_haplo_mark; RAGGED output */

void print_tiny_map();  /* arg: map; one liner (loci, log-like) */
//void print_short_map(); /* arg: map; three liner (loci, rec_fracs, log-like) */
void print_long_map();  /* arg: map; expanded map output */
void print_special_map();  /* see print.c */
void print_list();      /* arg: list; prints list */
void print_trys();
//void print_permsex();

void print_f2_map_genotypes(); /* prints genotypes w/ X-overs and errors
  args: MAP *map; bool explode_haplos; int num_old; int *old_locus; 
  can omit old_locus if num_old==0 old_locus are printed as is, non-old 
  get parentheses */

void print_haplo_summary();   /* int *loci, num_loci; */
void print_locus_summary();   /* int *loci, num_loci; bool haplos_too; */
void print_mapping_summary(); /* int *loci, num_loci; bool haplos_too; */

void new_print_placements();
/* args: int *order, num_order, *unplaced_loci, num_unplaced; bool **excluded;
   index into unplaced_loci[] is the same as the 1st index into excluded[][] */

#define MAP_DIVIDER "===============================================================================\n"
#define SUB_DIVIDER "-------------------------------------------------------------------------------\n"

#endif
