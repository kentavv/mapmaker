#ifndef _QLOW_H_
#define _QLOW_H_

/******************************************************************************

  ####   #        ####   #    #          #    #
 #    #  #       #    #  #    #          #    #
 #    #  #       #    #  #    #          ######
 #  # #  #       #    #  # ## #   ###    #    #
 #   #   #       #    #  ##  ##   ###    #    #
  ### #  ######   ####   #    #   ###    #    #

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

/* qlow.h - declarations needed only deep within the bowls of QCTM */

/* These locus and interval genotype codes are used in lots of places.
   Do NOT change these definitions: if you do, all hell may break loose! 
   Here A means an A/A individual, B means B/B, H means A/B, and 
   I means B/A. For BACKCROSS, use A and H, not A and B for consistency 
   (in particular, so that count_recs and a few other things work). 
   The possible genotype encoding software in qraw.c presumes the 
   codes defined here! */

#define MAX_LOCUS_GENOTYPES 3
#define A 0
#define H 1
#define B 2
#define I 3 /* not used */

#define for_locus_genotypes(data_type,geno) \
 for (geno=0; geno<(data_type==BACKCROSS ? 2:3); geno++)

#define MAX_INTERVAL_GENOTYPES 9
#define AA 0
#define AH 1
#define HA 2
#define HH 3
#define AB 4
#define HB 5
#define BA 6
#define BB 7
#define BH 8

#define for_interval_genotypes(data_type,geno) \
 for (geno=0; geno<(data_type==BACKCROSS ? 4:9); geno++)
/**** used to be 3 for BACKCROSS - bug fix 6/7/90 ****/

/* these are declared and set in qraw.c */
extern int *left_genotype, *right_genotype; 
    /* index with an interval genotype code => locus genotype code */
extern int **make_interval_genotype;
    /* index with [left_genotype_code][right_genotype_code]=> interval code */

typedef real INTERVAL_GENOTYPE_PROBS[16], LOCUS_GENOTYPE_PROBS[4];
typedef real GENO_PROBS[9][4];    /* 16 and 4 for indexing efficiency */
typedef int  INT_POSS_GENOTYPES[16];

#define MAX_RECS        2
#define MAX_REC_FRAC	0.49995
#define MIN_REC_FRAC	0.00004
#define MAX_FRAC_OF_RF	0.99999
#define MAX_CM		998.5
#define MIN_CM		0.04
#define MAX_FRAC_OF_CM  0.99999
#define MAX_CHROM_LOC   100

/*** This struct holds the data needed by QCTM ***/

typedef struct {
    int  num_intervals;        
    int  num_individuals;	/* < raw.n_indivs if some are missing pheno */
    int  num_continuous_vars;	
    int  max_intervals;		/* max # this struct is allocated for */
    int  max_individuals;	
    int  max_continuous_vars;   
    real *interval_len; 	/* [interval#] => length (a rec-frac really) */
    INTERVAL_GENOTYPE_PROBS **genotype_prob;    /* [indiv#][int#][geno-code] */
    real *phenotype; 		/* [individual#] */
    real **cont_var;		/* [var#][individual#] */
    int  **num_genotypes;	/* [individual#][interval#]=> #of poss genos */
    INT_POSS_GENOTYPES **genotype; /* [individual#][interval#][genotype_num] */
} DATA;

extern DATA *data; /* a global used by make_qtl_map() */

/*** Functions in QDATA.C ***/
DATA *alloc_data();  /* args: int max_intervals, max_continuous_vars; */
void free_data();    /* args: DATA *ptr; */

void prepare_data(); /* args: QTL_MAP *map; DATA *data; in qdata.c */
/* Prepare_data() uses the trait, left, right and num_intervals elements
   of the map struct, along with the global raw struct, to side-effect an 
   existing data struct. */	   

void initial_qctm_values(); /* args: DATA *data; QTL_MAP *map; in qdata.c */
/* map->fix_pos, ->fix_weight, ->fix_dominance, ->trait, ->left, and
   ->right must be set, and prepare_data() must have been run on data;
   map is side-effected. */

/*** In QCTM.C ***/
void qtl_conv_to_map(); /* args: DATA *data; QTL_MAP *map; in qctm.c */
/* side-effects map */

/* Those three routines do all the work needed to make QTL maps. */

/* These are in qprint.c, and are called in the guts of QCTM */
void print_iteration();        
void print_null_iteration();
void do_print_E_step();


/*** QCTM's state and global variables - declared and used in qctm.c, but
     alloced in qdata.c ***/

extern int max_intervals; /* this var is used all over the code */
extern int max_genotype_vars; 
extern int max_interx_genotypes, max_backx_genotypes;
extern int **lookup_genotype;
extern real **lookup_coded_genotype;

extern real **expected_genotype;		 
extern real **S_matrix, **S_inverse;         	 
extern real **indiv_S_matrix;                	
extern real *int_like;	                         
extern real *qctm_qtl_pos;  	               
extern real *qctm_qtl_weight, *null_qtl_weight, *fix_qtl_weight;
extern real *temp_row; 	
extern GENO_PROBS *indiv_rec_like;             
extern GENO_PROBS *expected_recs, *trans_prob;
extern INTERVAL_GENOTYPE_PROBS *rec_like;    
extern char geno_chars[10],default_backcross_chars[10];
extern char default_intercross_chars[10];


/***** Stuff one needs to muck around with RAW data *****/

typedef struct {
        int data_type;  	/* valid data types are defined in qg.h */
	bool f3;		
	char file[PATH_LENGTH+1];
	int filenumber;
	int n_loci, n_traits, max_traits;
	int n_indivs, name_len;
	int n_chroms;
	int *original_locus;   /* [qtl_locus_number] */
	int *chrom_start, *chrom_n_loci;
	char **locus;		/* [indvidual][locus#] => one char genotype */
	char **locus_name;	/* [locus#] => string */
	LOCUS_GENOTYPE_PROBS **left_cond_prob, **right_cond_prob; 
	                        /* [indiv#][locus#][locus-geno-code] */
	real **trait;		/* [individual#][trait#] */
	char **trait_name;	/* [trait#] => string */
	char **trait_eqn;	/* [trait#] => string */
	real *map_dist;		/* [locus#] for dist to locus[locus#+1] */
} RAW;
/* FROB */
#define MAX_TRAITS(traits_in_file)(traits_in_file>10 ? 6*(traits_in_file) : 20)
#define MISSING_PHENO	-100000.0

extern RAW raw;

/*** Functions defined in QRAW.C ***/
void raw_init();             /* no args */
void read_data();            /* args FILE *fp; char *path; */
void read_olddata();         /* same args as read_data     */
void crunch_data();          /* no args; makes L&R cond_probs in raw struct */
bool data_loaded();          /* no args; return TRUE if data is loaded */
real map_length();           /* args: int left_locus, right_locus; */
real transition_prob();      /* args: int data_type, left, right; real theta;*/
real apriori_prob();         /* args: int data_type, geno; */
void indiv_interval_probs(); /* args: INTERVAL_GENOTYPE_PROBS *prob; 
			        int index, indiv, left, right; real theta; */
void save_traitfile();       /* command call "store status" */
void free_raw();             /* no args; never been used! */



/* Definitions from MAPMAKER needs for prep_dat */
typedef struct {
    int         max_loci;
    int         num_loci;
    int         *locus;         /* [num_loci] (markers in the map) */
    double      **rec_frac;     /* [intervals][M/F], recombination fractions */
    int         unlink;         /* indicates interval (if any) held unlinked */
    int         *fix_interval;  /* [intervals], intervals not converging */ 
    int         sex_specific;
    double      log_like;
}       MAP;

extern MAP *raw_map;

typedef struct {
    int         max_maps;      /* used in allocation of structure */
    int         max_loci;      /* used in allocation of structure */
    MAP         **map_list;    /* list of max_maps pointers to MAPs */
    MAP         *extra_map;    /* expendable map-for insertion and deletion */
    int         num_maps;      /* number of maps in map_list (<= max_maps) */
    double      threshold;     /* highest acceptable likelihood when sorting */
    bool        sorted;        /* determines whether this is a sorted list */
} SAVED_LIST;

extern SAVED_LIST *list;

#define INTERCROSS_ALLELE_CHARS "ABCDH-"
#define NAME_LEN 10
#define TRAIT_EQN_LEN 250
#define UNSORTED 0
#define NONE_UNLINKED -1

void create_status();
SAVED_LIST *allocate_map_list();
MAP *allocate_map();
void sort_last();
void free_map();
void free_map_list();
int symbol_value();

#define PARENTAL_TYPE_A           'A'
#define PARENTAL_TYPE_B           'B'
#define TYPE_NOT_A                'C'
#define TYPE_NOT_B                'D'
#define HYBRID_TYPE_H             'H'
#define MISSING_DATA              '-'

extern bool fix_weight_kludge; /* BIG KLUDGE for tweak-weight cmd */

#endif