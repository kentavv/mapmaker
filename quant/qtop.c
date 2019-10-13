/******************************************************************************

  ####    #####   ####   #####            ####
 #    #     #    #    #  #    #          #    #
 #    #     #    #    #  #    #          #
 #  # #     #    #    #  #####    ###    #
 #   #      #    #    #  #        ###    #    #
  ### #     #     ####   #        ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

/* QTOP.C - The toplevel for MAPMAKER/QTL, including lots of utility functions
   for the commands, and the main() procedure. */

//#define INC_LIB
//#define INC_SHELL
//#define INC_CALLQCTM
//#define INC_QTOPLEVEL
//#define INC_QLOWLEVEL
#include "qtl.h"

/***** User State Variables - global *****/
int units, print_names, print_mapm_loci, print_scans, print_maps; 
int trait;
char *trait_string;
QTL_MAP *map;  /* used by many commands */
int num_intervals, num_orders, num_ints_to_wiggle, num_continuous_vars;
bool *free_genetics;

/***** declarations to keep shared code happy *****/
/*********** QTL should never use these ***********/
bool sex_specific = 0;
real startrecombs = .05;
/**************************************************/
 
/***** for internal use only *****/
//int main();
//int get_pair_entry();
//void qtl_top();
//char *genetics_str();
DATA *map_data; /* allocated max_intervals long - used for scratch space */

bool isa_locus_name(char *str, int *num /* the locus# if TRUE, else the #matched if FALSE */, bool *exact);
bool isa_seq_name(char *str, int *num /* the #matched if FALSE - undefined if true */, bool *exact);
bool isa_trait_name(char *str, int *num /* the trait# if TRUE, else the #matched if FALSE */, bool *exact);

void make_help_entries(void);

/********** Message definitions for QTL **********/
char *BADDATA_ln;
char *BADDATA_error;
int   BADDATA_line_num;
real  MATINV_diff;
int   NOPARSE_err;
int   BAD_MKINT_err;
int   BADSEQ_errpos;
char *BADSEQ_errmsg;
char *BADTRAIT_errmsg;

//void ps_MATINV();
//void ps_BADDATA();

void
ps_MATINV (char *a)  { sprintf(ps, "diff= %lf", MATINV_diff); }
void 
ps_BADDATA (char *a) { sprintf(ps, "line=\"%s\"\nerror=%s",
                            truncstr(BADDATA_ln,60), truncstr(BADDATA_error,70)); }

void 
top_init (void)
{ 
  BADDATA_ln= BADDATA_error= NULL; 
  MATINV_diff= 0.0; NOPARSE_err= 0;
  BADSEQ_errpos=0; BADSEQ_errmsg= BADTRAIT_errmsg= NULL; 

  array(BADDATA_error,MAXLINE+1,char);
  map=NULL; map_data= NULL; 
  trait= -1; array(trait_string,MAXLINE,char);
  num_intervals= num_orders= num_ints_to_wiggle= num_continuous_vars= 0;
  free_genetics= NULL;

  setmsg(SINGMAT,	"attempt to invert singular matrix",sender,NULL); 
  setmsg(MATINV,	"matrix invert failed",             sender,ps_MATINV);
  setmsg(BADDATA,	"bad quant data file",              sender,ps_BADDATA);
  setmsg(QUIT_QTL,	"quit program",                     sender,NULL);
  setmsg(NOPARSE,	"parse failed",                     sender,NULL);
  setmsg(BAD_MKINT,	"parse failed",                     sender,NULL);
  setmsg(BADSEQ, 	"illegal QTL interval list",        sender,NULL);
  setmsg(BADTRAIT,	"illegal QTL trait spec",           sender,NULL);
}


/********** Welcome to the Top of the World **********/

int main(int argc, char *argv[])
{
    char help_filename[PATH_LENGTH+1];

    custom_lib_init();
    get_cmd_line_args(&argc,argv);
    tty_hello();
    seedrand(RANDOM);

    strcpy(help_filename,"qtl");
    shell_init("MAPMAKER/QTL","1.1b","1988-1992",help_filename);
    banner();
    photo_banner_hook = NULL;

    /* Initialize the various pieces of QTL */
    top_init();		/* qtop.c (you are there) */
    seq_init();		/* qseq.c */
    data_init();	/* qdata.c */
    raw_init();		/* qraw.c */
    qctm_init();	/* qctm.c */
    context_init();     /* qcontext.c */
    wiggle_init();      /* qwiggle.c */
    
    self_delimiting= mkstrcpy("[]<>");

    cmd_init(); 	/* qcmds.c - among other things, this initializes */
                        /*	     the user accessible variables */
    make_help_entries(); 

    run { help_file=open_file("qtl.help",READ); } except_when(CANTOPEN) {}

    if (!nullstr(file_arg[PHOTO_FILE_ARG])) {
        run {
                nl();
                if ((append_it && !photo_to_file(file_arg[PHOTO_FILE_ARG], APPEND)) ||
                    (!append_it && !photo_to_file(file_arg[PHOTO_FILE_ARG], WRITE)))
                    send(CANTOPEN);
                sprintf(ps, "photo is on: file is '%s'\n", photo_file);
                pr();
            } on_error { print("error opening photo file\n"); }
    }

    if (!nullstr(file_arg[RUN_FILE_ARG])) {
        run {
                redirect_input(file_arg[RUN_FILE_ARG], TRUE);
            } on_error { print("error opening run file\n"); }
    }

    command_loop();
    /* screen_end(); */
    exit_main();
}


/* The screen header:
         1         2         3         4         5         6         7
1234567890123456789012345678901234567890123456789012345678901234567890123456789
-------------------------------------------------------------------------------
MAPMAKER/QTL v1.1b     genetics: constrained     data: 123456789012.dat
trait: 1 f(1234567890)    units: kossambi CM    photo: 123456789012.out
sequence #123: blah blah blah
-------------------------------------------------------------------------------
*/

#define LINE0 "MAPMAKER/QTL V1.1b     genetics: %-11s  data: %s"
#define LINE1 "trait: %-15s units: %-11s     photo: %s"
#define LINE2 "sequence #%d: %-64s"

void 
qtl_top (char **line, int lines, int cols)
{
    char tr[16], file[PATH_LENGTH+1];

    if (!data_loaded()) { sprintf(line[0], LINE0, "", "<none>"); } else { 
	nstrcpy(file,raw.file,cols-55); 
	sprintf(line[0], LINE0, "obsolete", file); 
    }

    if (!data_loaded() || trait<0 || trait>raw.n_traits) strcpy(tr,"");
      else sprintf(tr, "%d (%s)", trait + 1, raw.trait_name[trait]);
    if (!log_open) strcpy(file,"<off>"); 
      else nstrcpy(file,photo_file,cols-55);
    sprintf(line[1], LINE1, tr, (units ? "haldane cM" : "rec-fracs"), file);

    if (!data_loaded() || ints==NULL) strcpy(line[2],"sequence: ");
      else sprintf(line[2], LINE2, context[active_context]->seq_history_num,
                   ints_string);
}    



/********** QTL_MAP handling routines **********/

QTL_MAP *
alloc_qtl_map (int num_intervals, int num_cont_vars)
{
    QTL_MAP *map=NULL;
        
    run {
	if(num_intervals<1) send(CRASH);
	if(num_cont_vars==0) num_cont_vars=1; /* alloc a space, don't use it */
	single(map, QTL_MAP);

	map->max_intervals=0; 
	map->left= map->right= NULL; map->interval_len= NULL;
	map->qtl_pos= map->qtl_weight= NULL;
	map->fix_pos= NULL; map->constraint=(GENETICS*)NULL;

	map->max_continuous_vars=0; map->cont_var=NULL; 
	map->cont_var_weight=NULL; map->fix_cont_var_weight=NULL;
	
	array(map->left, num_intervals, int);
	array(map->right, num_intervals, int);
	array(map->interval_len, num_intervals, real);
	array(map->qtl_pos, num_intervals, real);
	array(map->qtl_weight, num_intervals, real); 
	array(map->qtl_dominance, num_intervals, real); 
	array(map->fix_pos, num_intervals, real);
	array(map->constraint, num_intervals, GENETICS);
	array(map->cont_var, num_cont_vars, int);
	array(map->cont_var_weight, num_cont_vars, real);
	array(map->fix_cont_var_weight, num_cont_vars, real);
	
	map->max_intervals= num_intervals; 
	map->max_continuous_vars= num_cont_vars;
	reset_map(map);

    } when_aborting { free_qtl_map(map); relay; }
    return(map);
}


void 
free_qtl_map (QTL_MAP *map)
{
    if (map==NULL) return;

    unarray(map->left, int);
    unarray(map->right, int);
    unarray(map->interval_len, real);
    unarray(map->qtl_pos, real);
    unarray(map->qtl_weight, real);
    unarray(map->qtl_dominance, real);
    unarray(map->fix_pos, real);
    unarray(map->constraint, GENETICS);
    unarray(map->cont_var, int);
    unarray(map->cont_var_weight, real);
    unarray(map->fix_cont_var_weight, real);

    map->left= map->right= NULL; map->interval_len= NULL;
    map->qtl_pos= map->qtl_weight= map->qtl_dominance= NULL;
    map->fix_pos= NULL; map->constraint= (GENETICS*)NULL;
    map->num_intervals= map->max_intervals= 0;
    map->max_continuous_vars=0; map->cont_var=NULL; 
    map->cont_var_weight=NULL; map->fix_cont_var_weight=NULL;
    unsingle(map, QTL_MAP);
}


/* Reset_map() does NOT change ->trait, ->qtl_weight, ->qtl_dominance,
   ->fix_weight, ->fix_dominance, ->interval_len, ->qtl_pos, ->fix_pos,
   ->left, ->right, cont_var, ->fix_cont_var_weight, ->cont_var_weight, etc,
   except to the extent that ->num_intervals=0 and
   ->num_continuous_vars=0. */

bool 
reset_map (QTL_MAP *map)    
{
    if (map==NULL) send(CRASH);

    map->num_intervals= 0; map->num_continuous_vars= 0; 
    map->mu= map->sigma_sq= map->abs_log_like= 0.0;
    map->null_mu= map->null_sigma_sq= map->null_log_like= 0.0;
    map->var_explained= map->chi_sq= 0.0;
    map->log_like= map->no_data_like= 0.0;

    return(TRUE);
}


int 
add_interval (
    QTL_MAP *map,
    int left,
    int right, /* loci */
    real fix_pos,
    GENETICS *genetics
)
{
    int i;
    if (map==NULL || map->num_intervals>=map->max_intervals) send(CRASH); 
/*  if (!check_interval(&left,&right,&fix_pos)) send(CRASH); KLUDGE */
    
    i= (map->num_intervals)++;
    map->left[i]= left; 
    map->right[i]= right;
    map->fix_pos[i]= fix_pos;
    map->constraint[i].backx_weight= genetics->backx_weight;
    map->constraint[i].interx_type=  genetics->interx_type;
    map->constraint[i].a= genetics->a;
    map->constraint[i].b= genetics->b;
    map->constraint[i].c= genetics->c;
    return(i);
}


int add_continuous_var(QTL_MAP *map, int trait, real fix_weight)
{
    int i;
    if (map==NULL || map->num_continuous_vars>=map->max_continuous_vars) 
      send(CRASH); 
    if (trait!=EPISTASIS_TERM && !valid_trait_num(trait)) send(CRASH);
    
    i= (map->num_continuous_vars)++;
    map->cont_var[i]= trait; 
    map->fix_cont_var_weight[i]= fix_weight;
    return(i);
}


/* This is (I think) as yet unused, and hence, really untested */
void 
mapcpy (   /* strcpy()-ish backasswards notation */
    QTL_MAP *to,
    QTL_MAP *from
)
{
    int i;
	
    if (to->max_intervals < from->num_intervals) send(CRASH);
    if (to->max_continuous_vars < from->num_continuous_vars) send(CRASH);

	to->trait=		from->trait;
	to->num_intervals=	from->num_intervals;
	to->num_continuous_vars=from->num_continuous_vars;

    for (i=0; i<from->num_intervals; i++) {
	to->left[i]= 		from->left[i];
	to->right[i]= 		from->right[i];
	to->interval_len[i]= 	from->interval_len[i];
	to->qtl_pos[i]= 	from->qtl_pos[i];
	to->fix_pos[i]= 	from->fix_pos[i];
	to->qtl_weight[i]= 	from->qtl_weight[i];
	to->qtl_dominance[i]= 	from->qtl_dominance[i];
	copy_genetics(&to->constraint[i],&from->constraint[i]);
    }

    for (i=0; i<from->num_continuous_vars; i++) {
	to->cont_var[i]=		from->cont_var[i];
	to->cont_var_weight[i]=		from->cont_var_weight[i];
	to->fix_cont_var_weight[i]=	from->fix_cont_var_weight[i];
    }		

    to->mu= 		from->mu;
    to->sigma_sq= 	from->sigma_sq;
    to->var_explained= 	from->var_explained;
    to->chi_sq= 	from->chi_sq;
    to->null_mu= 	from->null_mu;
    to->null_sigma_sq= 	from->null_sigma_sq;
    to->null_log_like= 	from->null_log_like;
    to->log_like= 	from->log_like;
    to->no_data_like=	from->no_data_like;
    to->abs_log_like= 	from->abs_log_like;
}


void 
copy_genetics (GENETICS *to, GENETICS *from)
{
    to->backx_weight= from->backx_weight;
    to->interx_type=  from->interx_type;
    to->a=from->a; to->b=from->b; to->c=from->c;
}


bool 
constrained (GENETICS *genetics)
{
    if (raw.data_type==BACKCROSS) return(genetics->backx_weight!=DONT_FIX);
    else return(genetics->interx_type!=FREE);
}


#ifdef OLD_CODE /* This was never updated to deal with continuous_vars, and
		   maybe some other things. */

/******************** QTL_MAP SAVING SOFTWARE ********************/

/* DIRECTIONS: Each MAP_LIST contains ONE QTL_MAP that is "free": that
is, not part of the sorted list and free to be side-effected and used.
Use get_unused_map() ONCE to get a ptr to this QTL_MAP. You may return
this map to the struct using either with return_unused_map() or save_map().
In the latter case, a ptr to a new unused QTL_MAP is returned to you.
Only attempt to use save_map() on QTL_MAPS gotten in one of these ways! Also,
only side-effect QTL_MAPs gotten in this way: the others are in the
MAP_LIST! If you free the MAP_LIST without returning the unused QTL_MAP,
it will not be freed and should be freed explicitly with free_qtl_map(). 
Otherwise, if you do return it, it will be freed, and thus any ptrs to it 
should be thrown away. */


SAVE_QTL_MAPS *
alloc_saved_maps (int num_maps, int num_intervals, real threshold)
{
	SAVE_QTL_MAPS *p;
	int i;
	
	if (num_intervals<0 || num_maps<0 ) send(CRASH);
	single(p, SAVE_QTL_MAPS);
	p->unused_checked_out= FALSE;
	p->num_maps= 0;
	p->max_maps= num_maps;
	p->max_intervals= num_intervals;
	p->threshold= threshold;
	array(p->map, num_maps, QTL_MAP*);
	for (i=0; i<num_maps; i++) 
		p->map[i]= alloc_qtl_map(num_intervals);
	p->unused= alloc_qtl_map(num_intervals);
	return(p);
}
	

void 
free_saved_maps (SAVE_QTL_MAPS *p)
{
	int i;
	
	for (i=0; i<p->max_maps; i++)
		free_qtl_map(p->map[i]);
	if (!p->unused_checked_out) 
	    { free_qtl_map(p->unused); p->unused= NULL; }
	p->map= NULL;
	p->max_maps=0;
	p->max_intervals=0;
	unsingle(p, SAVE_QTL_MAPS);
}


QTL_MAP *
save_map (QTL_MAP *map_to_save, SAVE_QTL_MAPS *the_maps)
{
	int i;
	QTL_MAP *prev, *temp;
	
	if (!the_maps->unused_checked_out || map_to_save!=the_maps->unused) 
		send(CRASH);
	if (map_to_save->log_like < the_maps->threshold) return(FALSE);
	for (i=0, prev=NULL; i<the_maps->num_maps; i++) {
	    if (prev!=NULL) { 
		    temp=the_maps->map[i];
		    the_maps->map[i]= prev;
		    prev= temp;
	    } else if (map_to_save->log_like > the_maps->map[i]->log_like) {
		    prev= the_maps->map[i];
		    the_maps->map[i]= map_to_save;
	    }
        }
	if (the_maps->num_maps < the_maps->max_maps && prev==NULL) {
		prev= the_maps->map[the_maps->num_maps];
		the_maps->map[the_maps->num_maps]= map_to_save;
		(the_maps->num_maps)++;
	} else if (the_maps->num_maps< the_maps->max_maps && prev!=NULL) {
		temp= the_maps->map[the_maps->num_maps];
		the_maps->map[the_maps->num_maps]= prev;
		prev= temp;
		(the_maps->num_maps)++;
	} else if (prev==NULL) return(map_to_save); /* don't save it */

	the_maps->unused= prev; 
	return(prev); 
}


QTL_MAP *
get_unused_map (SAVE_QTL_MAPS *the_maps)
{
	if (the_maps->unused_checked_out) send(CRASH);
	the_maps->unused_checked_out= TRUE;
	return(the_maps->unused);
}


void 
return_unused_map (QTL_MAP *map_to_return, SAVE_QTL_MAPS *the_maps)
{
  	if (!the_maps->unused_checked_out || 
	    map_to_return != the_maps->unused) send(CRASH);
	the_maps->unused_checked_out= FALSE;
}

#endif

	
/********** OTHER USEFUL THINGS FOR QTL COMMANDS **********/

#define eNEED_DATA \
"No data have yet been loaded.\nUse the 'load data' command first."
#define eNEED_INTERX \
"Backcross data have been loaded.\nThe '%s' command requires intercross data."
#define eNEED_BACKX \
"Intercross data have been loaded.\nThe '%s' command requires backcross data."
#define eNEED_INTS \
"The interval sequence has not been set.\nUse the 'sequence' command first."
#define eNEED_TRAIT \
"A trait has not yet been selected.\nUse the 'trait' command first."
#define eWRONG_SEQ \
"The current sequence is not apropriate for the '%s' command.\nThe rightmost interval may not have a repeat count nor any '@' symbols."
#define eONE_INT \
"The current sequence is not apropriate for the '%s' command.\nOnly one interval may be specified."  /* FROB */


void 
qtl_ready (int data_type, int need_seq, bool need_trait, bool will_call_qctm)
/* Side-effects globals num_intervals, num_orders, num_ints_to_wiggle, 
   free_genetics, and allocate map and map_data, if needed. */
{

    if (data_type!=NO_DATA && !data_loaded()) error(eNEED_DATA);
    if (data_type==INTERCROSS && raw.data_type!=INTERCROSS) 
      { sprintf(ps, eNEED_INTERX, com); error(ps); }
    if (data_type==BACKCROSS && raw.data_type!=BACKCROSS) 
      { sprintf(ps, eNEED_BACKX, com); error(ps); }

    if (need_seq!=NOSEQ) {
	if (ints==NULL) error(eNEED_INTS);
	if (!reset_state(ints,(need_seq==WIGSEQ),&num_intervals,
		&num_continuous_vars,&num_orders,&num_ints_to_wiggle)) 
	  { sprintf(ps, eWRONG_SEQ, com); error(ps); }
	if (need_seq==ONEINT && num_intervals!=1) 
	  { sprintf(ps, eONE_INT, com); error(ps); } /* FROB */
	if (free_genetics==NULL) array(free_genetics,abs(max_intervals),bool);
	get_seq_free_genetics(ints,free_genetics);
    }

    if (need_trait) {
	if (trait<0 || trait>raw.n_traits-1) error(eNEED_TRAIT);
	if (map!=NULL) map->trait= trait;
    }

    if (will_call_qctm) {
	if (data_type==NO_DATA) error(eNEED_DATA);
	if (max_intervals<0) max_intervals= -max_intervals;
	if (max_continuous_vars<0) max_continuous_vars= -max_continuous_vars;
	if (!qctm_globals_avail()) alloc_qctm_globals();
	if (map==NULL) map=alloc_qtl_map(max_intervals,max_continuous_vars);
	if (map_data==NULL) map_data=alloc_data(max_intervals,
						max_continuous_vars);
    }
}


#define eNOTRAIT    "Missing trait name or number following '*'."
#define eNOLOCUS    "Missing locus name or number following '*'."
#define eTRAITKILLED "Trait %d has been deleted from the data set."
#define eNTRAITTOKS "A trait must be specified by a single name or number."
#define eNLOCUSTOKS "A locus must be specified by a single name or number."
#define eONETRAIT \
 "Trait number is out of range - There is only one trait in the data set."
#define eNTRAITS    "Trait number is out of range - Use a number from 1 to %d."
#define eNLOCI      "Locus number is out of range - Use a number from 1 to %d."
#define eTRAITNAME  "There is no trait named '%s'."
#define eLOCUSNAME  "There is no locus named '%s'."
#define eNAMEAMBIG  "The name '%s' is ambiguous - Supply more characters."
#define eNAMEINUSE  "The name '%s' is already in use - Try a different name."
#define eNAMETOOLONG \
  "The name '%s' is too long - limit names to %d characters."
#define eTRAITNEEDNAME "You must specify a name (not a number) for the trait."

bool 
valid_locus_str (  /* return TRUE or FALSE */
    char *str,     /* assume this is a despaced/lowercased/filtered token */
    int *num,      /* side-effect with the number if TRUE */
    char *errmsg  /* side-effect with a message if FALSE */
)
{
    int match, exact;
    char *name= get_temp_string();
    if (errmsg==NULL) send(CRASH);
    
    if (nullstr(str)) {
	strcpy(errmsg,eNOTRAIT);
	return(FALSE);
	
    } else if (itoken(&str,iREQUIRED,num)) {
	/* --*num; valid_locus_num now handles this */
	if (!valid_locus_num(num)) {
	    sprintf(errmsg, eNLOCI, raw.n_loci); 
	    return(FALSE);
	} else if (!nullstr(str)) { /* check for only one token */
	    strcpy(errmsg,eNLOCUSTOKS); return(FALSE);
	} else return(TRUE);

    } else if(stoken(&str,sREQUIRED,name)) {
	if (name[0]=='*') { 
	    name++; 
	    if (nullstr(name)) { strcpy(errmsg,eNOLOCUS); return(FALSE); }
	}
	if (!valid_name(name)) 
	  { strcpy(errmsg,"illegal name"); return(FALSE); }
	if (len(name)>NAME_LEN) 
	  { sprintf(errmsg, eNAMETOOLONG, name, NAME_LEN); return(FALSE); }
	if (!nullstr(str)) { /* check for only one token */
	    strcpy(errmsg,eNLOCUSTOKS); return(FALSE);
	} else {
	    if (!isa_locus_name(name,num,&exact)) {
		if (*num==0) sprintf(errmsg, eLOCUSNAME, name);
		else /* *num>=2 */ sprintf(errmsg, eNAMEAMBIG, name); 
		return(FALSE);
	    }
	    if (!exact && (isa_trait_name(name,&match,&exact) ||
			   isa_seq_name(name,&match,&exact))) {
		    if (match>0) sprintf(errmsg, eNAMEAMBIG, name);
		    else sprintf(errmsg, eLOCUSNAME, name);
		    return(FALSE);
	    }
	    return(TRUE); 
	}
	
    } else { /* !itoken && !stoken */
	strcpy(msgstr,"Illegal token???"); 
	return(FALSE);
    }
}


bool 
valid_trait_str (  /* return TRUE or FALSE */
    char *str,    /* assume this is a despaced/lowercased/filtered token */
    int *num,     /* side-effect with the number if TRUE */
    char *errmsg /* side-effect with a message if FALSE */
)
{
    int match, exact;
    char *name= get_temp_string();
    
    if (errmsg==NULL) send(CRASH);
    
    if (nullstr(str)) {
	strcpy(errmsg,eNOTRAIT);
	return(FALSE);
	
    } else if (itoken(&str,iREQUIRED,num)) {
        --*num; 
	if (!valid_trait_num(*num)) {
	    if (raw.n_traits==1) strcpy(errmsg,eONETRAIT); 
	    else sprintf(errmsg, eNTRAITS, raw.n_traits); 
	    return(FALSE);
	} else if (nullstr(raw.trait_name[*num])) { 
	    sprintf(errmsg, eTRAITKILLED, *num + 1); return(FALSE);
	} else if (!nullstr(str)) { /* check for only one token */
	    strcpy(errmsg,eNTRAITTOKS); return(FALSE);
	} else  return(TRUE);
	
    } else if(stoken(&str,sREQUIRED,name)) {
	if (name[0]=='*') { 
	    name++; 
	    if (nullstr(name)) { strcpy(errmsg,eNOTRAIT); return(FALSE); }
	}
	if (!valid_name(name)) 
	  { strcpy(errmsg,"illegal name"); return(FALSE); }	  
	if (len(name)>NAME_LEN) 
	  { sprintf(errmsg, eNAMETOOLONG, name, NAME_LEN); return(FALSE); }
	if (!nullstr(str)) { /* check for only one token */
	    strcpy(errmsg,eNTRAITTOKS); return(FALSE);
	} else {
	    if (!isa_trait_name(name,num,&exact)) {
		if (*num==0) sprintf(errmsg, eTRAITNAME, name);
		else /* *num>=2 */ sprintf(errmsg, eNAMEAMBIG, name); 
		return(FALSE);
	    }
	    if (!exact && (isa_locus_name(str,&match,&exact) ||
			   isa_seq_name(str,&match,&exact))) {
		    if (match>0) sprintf(errmsg, eNAMEAMBIG, name);
		    else sprintf(errmsg, eTRAITNAME, name);
		    return(FALSE);
	    }
	    return(TRUE); 
	}
	
    } else { /* !itoken && !stoken */
	strcpy(msgstr,"Illegal token???"); 
	return(FALSE);
    }
}


bool 
valid_new_trait_name (  /* return TRUE or FALSE */
    char *str,     /* assume this is a despaced/lowercased/filtered token */
    char *errmsg  /* side-effect with a message if FALSE */
)
{
    int match, exact;
    char *name= get_temp_string();
    
    if (errmsg==NULL) send(CRASH);
    
    if (nullstr(str)) {
	strcpy(errmsg,eNOTRAIT);
	return(FALSE);
	
    } else if (itoken(&str,iREQUIRED,&match)) {
	strcpy(errmsg,eTRAITNEEDNAME);
	return(FALSE);
	
    } else if(stoken(&str,sREQUIRED,name)) {
	if (name[0]=='*') { 
	    name++; 
	    if (nullstr(name)) { strcpy(errmsg,eNOTRAIT); return(FALSE); }
	}
	if (!valid_name(name)) 
	  { strcpy(errmsg,"illegal name"); return(FALSE); }
	if (len(name)>NAME_LEN) 
	  { sprintf(errmsg, eNAMETOOLONG, name, NAME_LEN); return(FALSE); }
	else if (!nullstr(str)) /* check for only one token */
	  { strcpy(errmsg,eNTRAITTOKS); return(FALSE); }
	else if ((isa_trait_name(name,&match,&exact) && exact) ||
		 (isa_locus_name(name,&match,&exact) && exact) || 
		 (isa_seq_name(name,&match,&exact) && exact)) 
	  { sprintf(errmsg, eNAMEINUSE, name); return(FALSE); }
	return(TRUE);
	
    } else { /* !itoken && !stoken */
	strcpy(msgstr,"Illegal token???"); 
	return(FALSE);
    }
}


bool 
valid_locus_num (int *num)
{ 
    int i, x, new_num = -1;

    for(i=0; i<raw.n_loci; i++) {
	x=raw.original_locus[i];
	if (x == *num) 
	  new_num = i;
    }
    *num = new_num;
    return(*num>=0);
}

bool 
valid_trait_num (int num)
{ return(irange(&num,0,raw.n_traits-1) && !nullstr(raw.trait_name[num])); }

/* bool new_trait_num(num) UNUSED?
int *num;
{ if (raw.n_traits<raw.max_traits)
    { *num=raw.n_traits++; return(TRUE); }
  else return(FALSE);
} */





/* Lower-level */

bool isa_trait_name(char *str, int *num /* the trait# if TRUE, else the #matched if FALSE */, bool *exact)
{
    int i, n_matched;

    *num=0; *exact=FALSE;
    /* see if we uniquely match a trait name */
    for (i=0, n_matched=0; i<raw.n_traits; i++)
      if (!nullstr(raw.trait_name[i]) && matches(str,raw.trait_name[i])) { 
	  n_matched++; *num=i; 
	  if (len(str)==len(raw.trait_name[i])) 
	    { *exact=TRUE; n_matched=1; break; }
      }
    return(n_matched==1 ? TRUE:FALSE);
}


bool isa_locus_name(char *str, int *num /* the locus# if TRUE, else the #matched if FALSE */, bool *exact)
{
    int i, n_matched;

    *num=0; *exact=FALSE;
    /* see if we uniquely match a locus name */
    for (i=0, n_matched=0; i<raw.n_loci; i++)
      if (matches(str,raw.locus_name[i])) { 
	  n_matched++; *num=i; 
	  if (len(str)==len(raw.locus_name[i])) 
	    { *exact=TRUE; n_matched=1; break; }
      }
/*  if (matches(str,"inf"))  { *num= INF_LOCUS; n_matched++; }
    if (matches(str,"ter"))  { *num= INF_LOCUS; n_matched++; }
    if (matches(str,"pter")) { *num= INF_LOCUS; n_matched++; }
    if (matches(str,"qter")) { *num= INF_LOCUS; n_matched++; } */
    return(n_matched==1 ? TRUE:FALSE);
}


bool isa_seq_name(char *str, int *num /* the #matched if FALSE - undefined if true */, bool *exact)
{
    char *seq, *full_name;
    int err;

    if (!get_named_entry(str,&seq,&full_name,
	  context[active_context]->named_sequences,&err)) {
	if (err==NAME_DOESNT_MATCH) *num=0;
	else /* err==NAME_IS_AMBIGUOUS */ *num=2;
	return(FALSE);
    } 
    *exact= (len(full_name)==len(str) ? TRUE:FALSE);
    return(TRUE);
}


/* Upper level */

void 
set_trait_spec (char *str)
{
    int trait_num;

    if (!valid_trait_str(str,&trait_num,msgstr)) {
	BADTRAIT_errmsg= msgstr;
	send(BADTRAIT);
    } 
    trait= trait_num;
    sprintf(trait_string, "%d (%s)", trait + 1, raw.trait_name[trait]);
}


void 
make_qtl_map (QTL_MAP *map)
{
    map->trait = trait;
    /* map_data is a global we keep handy just for this! */
    prepare_data(map,map_data);
    initial_qctm_values(map_data,map);
    qtl_conv_to_map(map_data,map);
}
