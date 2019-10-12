/******************************************************************************

  ####    #####    ##     #####  ######           ####
 #          #     #  #      #    #               #    #
  ####      #    #    #     #    #####           #
      #     #    ######     #    #        ###    #
 #    #     #    #    #     #    #        ###    #    #
  ####      #    #    #     #    ######   ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#define INC_LIB
#define INC_SHELL
#include "mapm.h"

void state_init_values();
#define MAX_CONTEXTS 1     

/* NOTE: keep the declarations of vars in the same order here, below in 
   reset_state(), in read/write_state, and in mapm.h. */

/* global state variables THESE SHOULD MOVE */

/**** general ****/
bool print_names;
real tolerance;
int  units;
bool auto_save;
/* also more_mode */

/**** two point ****/
real default_lod;
real default_theta;
bool use_haplotypes;

/**** three point ****/
bool use_three_pt;
real triplet_lod;
real triplet_theta;
int  triplet_num_links;
real three_pt_threshold;
int  three_pt_window;
real triplet_error_rate;

/**** order maker ****/
real npt_threshold;
real npt_first_threshold;
int  npt_window;
int  npt_min_indivs; /* infomativeness criteria */
bool npt_codominant;
real npt_min_theta;
bool print_all_maps;

/* F2 only */
bool fake_maps;  /* just for debugging */
bool use_error_rate;
real error_lod_thresh;
real error_single_thresh;
real error_net_thresh;

/* Contexts */
STATUS_CONTEXT **context;
int active_context;
int num_contexts;

/* CEPH only or Obsolete */
int  sex_specific; 	
int  triplet_sex;
int  compress_DNA;
int  print_problem_size;
long max_problem_size;
bool segregation_distortion;
bool use_hmm;
int  print_maps;
bool print_all_errors;
int  time_stamping;
real inner_tolerance;
int  inner_loop;
int  convergence_rule;
int  original_markers;
real startrecombs;
int  print_dots;


void 
state_init (void)
{
    parray(context,MAX_CONTEXTS,STATUS_CONTEXT);
}


void 
reset_state (void)
{
    allocate_context(context[0]);
    state_init_values();
    set_current_seq(NULL,FALSE);
}


void 
undo_state (void) 
{ free_context(context[0]); state_init_values(); }
    

void 
state_init_values (void)
{

/**** General ****/

    print_names=FALSE;
    /* if TRUE, proper locus names will be displayed instead of locus numbers 
       in maps */

    tolerance=0.001; 
    /* convergence tolerance: a map is considered converged when two
       successive maps have a log-likelihood difference less than the
       tolerance */

    units=CENTIMORGANS; map_func(HALDANE);
    /* if RECFRACS, recombination fractions will be displayed 
       rather than cM distances in short maps and other functions */

    auto_save=TRUE;
    /* everything is saved upon quitting if TRUE */

/**** two point ****/

    default_lod=3.0; default_theta=FIFTY_CM;
    /* default linkage criteria used in two point */

    use_haplotypes=TRUE;
    /* whether we should use the haplotype info or not - should be set only by
       the haplotype command, which needs to do many other things when this 
       changes. NO LONGER REALLY A STATE VAR - IS FIXED AT TRUE. */

/**** three point ****/

    use_three_pt=TRUE;
    /* use fast 3pt extensions to place, etc by default */

    triplet_lod=3.0; triplet_theta=FIFTY_CM; triplet_num_links=2;
    /* default linkage criteria used for three point triples */

    three_pt_threshold=4.0; three_pt_window=9;
    /* likelihood diff for which a three point order is considered excluded.
       NOTE THRESH IS >0, ALTHOUGH MOST OF THE TIME IT SHOULD BE NEGATED! */

    triplet_error_rate=0.0;
    /* error rate to use in 3pt calculations, 0.0 means use_error_rate=FALSE,
       LOCUS_ERROR_RATE means to use the per locus rates (as Npt does), or
       a real# is a fixed error rate (must be 0.0 or fixed to use fast3pt) */
    
/**** order maker ****/

    npt_threshold=2.0; npt_first_threshold=3.0; npt_window=7;
    /* likelihood diffs for which a npt order is considered excluded.
       NOTE THRESHS ARE >0, ALTHOUGH MOST OF THE TIME IT SHOULD BE NEGATED! */

    npt_min_indivs=1; npt_codominant=FALSE; npt_min_theta=ONE_CM;
    /* default informativeness criteria */

    print_all_maps=FALSE;
    /* for place, order, build, together, etc */

/**** F2 Stuff, including error checker ****/
    
    fake_maps=FALSE;
    /* bypass hmm ctm - for testing */

    use_error_rate=FALSE;
    /* whether ctm should use the error kludge or not */
    
    error_lod_thresh=1.0;
    /* error lod at which print_long_map barks, and which get summed into net*/

    error_single_thresh=2.0; error_net_thresh=3.0; 
    /* error lod at which order making barks */

/**** Contexts, Partially Obsolete ****/

    num_contexts=1; active_context=0;
    /* these never change as yet */

/**** CEPH and Other Obsolete Stuff ****/

    sex_specific= context[active_context]->sex_specific= FALSE;
    /* set to TRUE if sex specific maps to be made */

    triplet_sex=FALSE;
    /* sex-specific flag for 3pt calculations - always FALSE for now */
    
    compress_DNA= context[active_context]->compress_DNA= TRUE;
    /* always on - might be used in future */

    max_problem_size= context[active_context]->max_problem_size= 10000000;
    /* (default = 10000000) above this threshold, maps will not 
       be computed, as they will take too long */

    inner_tolerance= 0.01; 
    /* inner convergence tolerance */

    inner_loop= TRUE;  
    /* inner looping switch */

    print_problem_size= FALSE;
    /* problem sizes appear as each ceph map is computed */

    startrecombs= 0.05;  
    /* starting value for all recombination fractions */

    convergence_rule= 0;  
    /* unchangeable - perhaps useful in the future */

    time_stamping= FALSE;  
    /* if TRUE, each ctm map is timed and the time is printed */

    print_maps= FALSE; /* obsolete? */
    /* if TRUE, full maps will be displayed as each map is computed in 
       some multi map functions (compare,try) */

    segregation_distortion=FALSE;  /* presently obsolete */
    /* if TRUE, segregation distortion is accounted for in F2 map making */

    use_hmm=TRUE; /* also obsolete */
    /* if TRUE use hmm ctm by default */

    print_all_errors= FALSE; /* unused? */
    /* if TRUE, all candidate errors will be printed by commands like map,
       otherwise only the worst N may be printed...  */
}



/**** General ****/

command 
set_print_names (void)
{
    mapm_ready(ANY_DATA,0,0,NULL);
    maybe_set_bool(&print_names);
}

command 
set_tolerance (void)
{
    mapm_ready(ANY_DATA,0,0,NULL);
    maybe_set_real(&tolerance,0.0000001,1.0,8.6);
}

command 
set_units (void)
{
    char type[TOKLEN+1];
    mapm_ready(ANY_DATA,0,0,NULL);
    get_one_arg(stoken,"",type);

    if (nullstr(type)) {
	if (units==RECFRACS)
	  print("the 'units' are currently recombination-fractions\n");
	else if (map_func_num()==HALDANE)
	  print("the 'units' are currently (Haldane) centimorgans\n");
	else if (map_func_num()==KOSAMBI)
	  print("the 'units' are currently (Kosambi) centimorgans\n");
	else send(CRASH);

    } else { /* have args */
	if (matches(type,"recombination fractions") || 
	    matches(type,"rec-fracs") ||
	    matches(type,"recombination-fractions") ||
	    matches(type,"rec fracs") ||
	    matches(type,"rf")) units=RECFRACS;
	else if (matches(type,"centimorgans") || 
		 matches(type,"cm")) units=CENTIMORGANS;
	else usage_error(1);

	if (units==RECFRACS)
	  print("the 'units' are now recombination-fractions.\n");
	else if (map_func_num()==HALDANE)
	  print("the 'units' are now set to (Haldane) centimorgans.\n");
	else if (map_func_num()==KOSAMBI)
	  print("the 'units' are now set to (Kosambi) centimorgans.\n");
	else send(CRASH);
    }
}

command 
set_cm_func (void)
{
    int i;

    char type[TOKLEN+1];
    mapm_ready(ANY_DATA,0,FALSE,NULL);
    get_one_arg(stoken,"",type);

    if (nullstr(type)) {
	sprintf(ps, "centimorgan function: %s\n", map_func_name()); pr();
	return;
    }
    /* assume here that the mapfunc names have different first letters */
    for (i=0; i<num_map_functions; i++)
      if (matches(type,maps[i].name)) {
	  map_func(i);
	  sprintf(ps, "centimorgan function: %s\n", map_func_name()); pr();
	  return;
      }
    error("unknown function name\nexpected either 'Haldane' or 'Kosambi'");
}


command 
set_autosave (void)
{
    mapm_ready(ANY_DATA,0,0,NULL);
    maybe_set_bool(&auto_save);
}

command 
set_more_mode (void)
{
    mapm_ready(MAYBE_DATA,0,0,NULL);
    maybe_set_bool(&more_mode);
}


/**** two-point ****/

#define DEF_LINK_LOD    "default LOD score threshold is %.2lf\n"
#define DEF_LINK_THETA  "default %s threshold is %.2lf\n"
command 
set_default_linkage (void)
{
    real lod, theta;
    bool have_theta=FALSE;

    if (!nullstr(args)) {
	if (!rtoken(&args,rREQUIRED,&lod) || !rrange(&lod,0.0,100.0))
	  set_usage_error("a real number from 0 to 100");
	if (!nullstr(args)) {
	    if (!rtoken(&args,rREQUIRED,&theta) || !input_dist(&theta))
	      set_usage_error("a real number, greater than zero (cm or rf)");
	    else have_theta= TRUE;
	}
	default_lod= lod; 
	if (have_theta) default_theta=theta;
    }

    sprintf(ps, DEF_LINK_LOD, default_lod); pr();
    if (units!=RECFRACS) sprintf(ps, DEF_LINK_THETA, "centimorgan distance",
		      ((*mapfunction->rec_to_dist)(default_theta))*100.0);
    else sprintf(ps, DEF_LINK_THETA, "recombination-fraction", default_theta);
    pr();
}


/**** three-point ****/

command 
set_use_3pt (void)
{
    mapm_ready(ANY_DATA,0,0,NULL);
    maybe_set_bool(&use_three_pt);
}

#define DEF_3PT_LOD   "triplet LOD score threshold is %.2lf\n"
#define DEF_3PT_THETA "triplet %s threshold is %.2lf\n"
#define DEF_3PT_NUM   "number of linkages required is %2d\n"
command 
set_3pt_linkage (void)
{
    real lod, theta;
    bool have_theta=FALSE, have_num=FALSE; 
    int num;

    mapm_ready(ANY_DATA,0,0,NULL);
    if (!nullstr(args)) {
	if (!rtoken(&args,rREQUIRED,&lod) || !rrange(&lod,0.0,100.0))
	  set_usage_error("a real number from 0 to 100");
	if (!nullstr(args)) {
	    if (!rtoken(&args,rREQUIRED,&theta) || !input_dist(&theta))
	      set_usage_error("a real number, greater than zero (cm or rf)");
	    else have_theta= TRUE;
	}
	if (!nullstr(args)) { /* still */
	    if (!itoken(&args,0,&num) || !irange(&num,2,3))
	      set_usage_error("an integer number, either 2 or 3");
	    else have_num= TRUE;
	}
	triplet_lod= lod;
	if (have_theta) triplet_theta=theta;
	if (have_num)   triplet_num_links=num;
    }
    sprintf(ps, DEF_3PT_LOD, triplet_lod); pr();
    if (units!=RECFRACS) sprintf(ps, DEF_3PT_THETA, "centimorgan distance",
		      ((*mapfunction->rec_to_dist)(triplet_theta))*100.0);
    else sprintf(ps, DEF_3PT_THETA, "recombination-fraction", triplet_theta);
    pr();
    sprintf(ps, DEF_3PT_NUM, triplet_num_links); pr();
}

#define TRIP_THRESH \
  "triplet log-likelihood exclusion threshold is %.2lf.\nwindow size is %d.\n"
command 
set_3pt_threshold (void)
{
    real like;
    int wind, have_num=FALSE;

    mapm_ready(ANY_DATA,0,0,NULL);
    if (!nullstr(args)) {
	if (!rtoken(&args,rREQUIRED,&like) || !rrange(&like,0.0,100.0))
	  set_usage_error("a real number from 0 to 100");
	if (!nullstr(args)) { /* still */
	    if (!itoken(&args,0,&wind) || !irange(&wind,3,99) || (wind%2)!=1)
	      set_usage_error("an odd integer number from 3 to 99");
	    else have_num=TRUE;
	}
	three_pt_threshold=like;
	if (have_num) three_pt_window=wind;
    }
    sprintf(ps, TRIP_THRESH, three_pt_threshold, three_pt_window); pr();
}

#define SET_3ERR_USAGE \
  "'on', 'off', or a real number from 0 to 10 (percent chance of error)"
command 
set_3pt_errors (void)
{
    real rate;

    mapm_ready(ANY_DATA,0,0,NULL);
    if (!nullstr(args)) {
	if (rtoken(&args,rREQUIRED,&rate)) {
	    if (!rrange(&rate,0.0,10.0)) set_usage_error(SET_3ERR_USAGE);
	    triplet_error_rate=rate/100.0;
	} else if (nmatches(args,"on",2))  triplet_error_rate=LOCUS_ERROR_RATE;
	else if (nmatches(args,"off",2)) triplet_error_rate=0.0;
	else set_usage_error(SET_3ERR_USAGE);
    }
    if (triplet_error_rate==0.0) print("'triple error detection' is off.\n");
    else if (triplet_error_rate==LOCUS_ERROR_RATE)
      print("'triple error detection' is on.\n");
    else {
	print("'triple error detection' in three-point analysis is on.\n");
	sprintf(ps, "'error probability' for all loci is fixed at %.2lf%%.\n",
	   triplet_error_rate*100.0);
    }
}


/**** Order Maker ****/

#define NPT_THRESH1 \
"multipoint exclusion threshold is %.2lf.\nwindow size is %d.\n"
#define NPT_THRESH2 \
"multipoint exclusion threshold is %.2lf.\nwindow size is %d.\n\
strict multipoint exclusion threshold is %.2lf.\n"

command 
set_npt_threshold (void)
{
    real like, like2;
    int wind, have_num=FALSE, have_2nd=FALSE;

    mapm_ready(ANY_DATA,0,0,NULL);
    if (!nullstr(args)) {
	if (!rtoken(&args,rREQUIRED,&like) || !rrange(&like,0.0,100.0))
	  set_usage_error("a real number from 0 to 100");
	if (!nullstr(args)) { /* still */
	    if (!itoken(&args,0,&wind) || !irange(&wind,3,99) || (wind%2)!=1)
	      set_usage_error("an odd integer number from 3 to 99");
	    else have_num=TRUE;
	    if (!nullstr(args)) { /* still */
		if (!rtoken(&args,rREQUIRED,&like2)||!rrange(&like2,0.0,100.0))
		  set_usage_error("a real number from 0 to 100");  
		else have_2nd=TRUE;
	    }
	}
	npt_threshold=like;
	if (have_num) npt_window=wind;
	if (have_2nd) npt_first_threshold=like2;
	  else npt_first_threshold=max(npt_first_threshold,npt_threshold);
    }
    if (npt_threshold==npt_first_threshold)
      { sprintf(ps, NPT_THRESH1, npt_threshold, npt_window); pr(); }
    else 
      { sprintf(ps, NPT_THRESH2, npt_threshold, npt_window, npt_first_threshold);pr(); }
}


#define INF_THRESHS \
"Informativeness Criteria: min Distance %s, min #Individuals %d%s\n"
#define CD_NOT \
"slecting codominant markers (or not) does not apply to the loaded data"
command 
set_inf_threshold (void)
{
    real theta;
    int indivs, codom, have_indivs=FALSE, have_cd=FALSE;
    
    mapm_ready(F2_DATA,0,0,NULL);
    if (!nullstr(args)) {
	if (!rtoken(&args,rREQUIRED,&theta) || !input_dist(&theta))
	  set_usage_error("a map distance, cM or RF");
	if (!nullstr(args)) { /* still */
	    if (!itoken(&args,0,&indivs) || 
		!irange(&indivs,3,raw.data.f2.num_indivs))
	      set_usage_error("an integer from 1 to #individuals in data set");
	    else have_indivs=TRUE;
	    if (!nullstr(args)) { /* still */
		if (matches(args,"codominant") || streq(args,"cd")) 
		  { codom=TRUE; have_cd=TRUE; }
		else if (matches("any",args) || matches("all",args))
		  { codom=FALSE; have_cd=TRUE; }
		else set_usage_error("either 'codominant' or 'any'");
	    }
	}
	if (have_cd &&
	    raw.data.f2.cross_type!=F2_INTERCROSS &&
	    raw.data.f2.cross_type!=F3_SELF) error(CD_NOT);
	npt_min_theta=theta;
	if (have_indivs) npt_min_indivs=indivs;
	if (have_cd)     npt_codominant=codom;
    }
    sprintf(ps, INF_THRESHS, rag(rf2str(npt_min_theta)), npt_min_indivs,
            (npt_codominant ? "\n(codominant markers only)":"")); pr();
}


command 
set_print_all_maps (void)
{
    mapm_ready(ANY_DATA,0,0,NULL);
    maybe_set_bool(&print_all_maps);
}



/**** F2 and error checker ****/

command 
set_fake_maps (void)
{
    mapm_ready(F2,0,0,NULL);
    maybe_set_bool(&fake_maps);
}

command 
set_use_error_rate (void)
{
    mapm_ready(F2,0,0,NULL);
    maybe_set_bool(&use_error_rate);
}

#define ERROR_THRESHES "error LOD threshold to print: %.2lf\n\
error LOD thresholds for mapping: single: %.2lf, net: %.2lf\n"
command 
set_error_lod_thresh (void)
{
    real t1, t2, t3;

    mapm_ready(F2,0,0,NULL);
    if (!nullstr(args)) {
	if (!rtoken(&args,rREQUIRED,&t1) || !rrange(&t1,0.0,99.0))
	  set_usage_error("a real number from 0 to 99");
	if (!nullstr(args)) {
	    if (!rtoken(&args,rREQUIRED,&t2) || !rrange(&t2,0.0,99.0))
	      set_usage_error("a real number from 0 to 99");
	    if (!nullstr(args)) {
		if (!rtoken(&args,rREQUIRED,&t3) || !rrange(&t2,0.0,99.0))
		  set_usage_error("a real number from 0 to 99");
	    }
	}
	error_lod_thresh=t1; error_single_thresh=t2; error_net_thresh=t3;
    }
    sprintf(ps, ERROR_THRESHES, error_lod_thresh, error_single_thresh,
            error_net_thresh); pr();
}


/**************** State Contexts ****************/

void 
allocate_context (STATUS_CONTEXT *con)
{
    run {
	con->sex_specific= FALSE;
	con->compress_DNA= TRUE;
	con->max_problem_size= 10000000;
#ifdef HAVE_CEPH
	con->use_number= USE_NO_DATA;
#endif
	con->named_sequences= allocate_table(MAX_NAMED_SEQS,MAX_SEQ_LEN,
					     EXPANDS_BY(MAX_NAMED_SEQS),
					     INDEX_BY_NAME);
	con->sequence_history= allocate_table(MAX_HISTORY_SEQS,MAX_SEQ_LEN,
					      CANT_EXPAND,INDEX_BY_NUMBER);
	con->seq_history_num= next_entry_number(con->sequence_history);

    } on_exit {
	free_context(con);
	relay_messages;
    }
}

void 
free_context (STATUS_CONTEXT *con)
{
    free_table(con->named_sequences);
    free_table(con->sequence_history);
}


bool change_context(new_context)
int new_context;
{
    if(context[new_context] != NULL) {
	active_context = new_context;
	/* set the globals */
	sex_specific = context[active_context]->sex_specific;
	compress_DNA = context[active_context]->compress_DNA;
	max_problem_size = context[active_context]->max_problem_size;
#ifdef HAVE_CEPH
	if(raw.data_type == CEPH)
	  raw.data.ceph.use_number = context[active_context]->use_number;
#endif
	return(TRUE);
    } else {
	return(FALSE);
    }
}


bool create_new_context(new_context)
int new_context;
{
    if(context[new_context] != NULL || new_context > MAX_CONTEXTS)
      return(FALSE);
    
    allocate_context(context[new_context]);
    
    /* take current values as defaults */
    context[new_context]->sex_specific = context[active_context]->sex_specific;
    context[new_context]->compress_DNA = context[active_context]->compress_DNA;
    context[new_context]->use_number = context[active_context]->use_number;
    context[new_context]->max_problem_size = 
      context[active_context]->max_problem_size;

    return(TRUE);
}


#ifdef OBSOLETE
void kill_context(con,save_it)
STATUS_CONTEXT *con;
bool save_it;
{
    char *name,*seqnce,*err;

    if(save_it) {
	for(Te=con->named_sequences->list; Te!=NULL; Te=Te->next) {
	    if(!name_sequence(Te->id.name,Te->string,&err))
	      error(err);
	}
    }
    free_context(con);
}
#endif


/**** CEPH/OBSOLETE? ****/

#ifdef OBSOLETE
command 
set_sex_specific (void)
{
    maybe_set_bool(&sex_specific);
    context[active_context]->sex_specific = sex_specific;
}

command 
set_3pt_sex (void) /* unused */
{
    mapm_ready(ANY_DATA,0,0,NULL);
    maybe_set_bool(&triplet_sex);
}

command 
set_segregation_distortion (void)
{
    maybe_set_bool(&segregation_distortion);
}

command 
set_print_maps (void)
{
    maybe_set_bool(&print_maps);
}

command 
set_inner_tolerance (void)
{
    maybe_set_real(&inner_tolerance,0.0000001,1.0,8.6);
}

command 
set_startrecombs (void)
{
    maybe_set_real(&startrecombs,0.0,0.5,4.2);
}

command 
set_inner_loop (void)
{
    maybe_set_bool(&inner_loop);
}

command 
set_print_problem_size (void)
{
    maybe_set_bool(&print_problem_size);
}

command 
set_max_problem_size (void)
{
    maybe_set_long(&max_problem_size,100,100000000);
    context[active_context]->max_problem_size = max_problem_size;
}

command 
set_time_stamping (void)
{
    maybe_set_bool(&time_stamping);
}

command 
set_print_dots (void)
{
    maybe_set_bool(&print_dots);
}

command 
set_use_hmm (void)
{
    maybe_set_bool(&use_hmm);
}

command 
set_use_haplotypes (void)
{
    maybe_set_bool(&use_haplotypes);
}
#endif /* OBSOLETE */



/**************** Save and Load Status info ****************/


void write_status(fp)
FILE *fp;
{
    int i, usenum;

    fprint(fp,"*MapmakerStatusInfo:\n");

    sprintf(ps, "*PrintNames: %d\n", print_names); fpr(fp);
    sprintf(ps, "*Tolerance: %lf\n", tolerance); fpr(fp);
    sprintf(ps, "*Units: %d\n", units); fpr(fp);
    sprintf(ps, "*MapFunc: %d\n", map_func_num()); fpr(fp);

    sprintf(ps, "*DefaultLOD: %lf\n", default_lod); fpr(fp);
    sprintf(ps, "*DefaultTheta: %lf\n", default_theta); fpr(fp);

    sprintf(ps, "*UseThreePoint: %d\n", use_three_pt); fpr(fp);
    sprintf(ps, "*TripletLOD: %lf\n", triplet_lod); fpr(fp);
    sprintf(ps, "*TripletTheta: %lf\n", triplet_theta); fpr(fp);
    sprintf(ps, "*TripletNumLinks: %d\n", triplet_num_links); fpr(fp);
    sprintf(ps, "*ThreePointThreshold: %lf\n", three_pt_threshold); fpr(fp);
    sprintf(ps, "*ThreePointWindow: %d\n", three_pt_window); fpr(fp);
    sprintf(ps, "*TripletErrorRate: %lf\n", triplet_error_rate); fpr(fp);

    sprintf(ps, "*NptThreshold: %lf\n", npt_threshold); fpr(fp);
    sprintf(ps, "*NptFirstThreshold: %lf\n", npt_first_threshold); fpr(fp);
    sprintf(ps, "*NptWindow: %d\n", npt_window); fpr(fp);
    sprintf(ps, "*NptMinIndivs: %d\n", npt_min_indivs); fpr(fp);
    sprintf(ps, "*NptCodominant: %d\n", npt_codominant); fpr(fp);
    sprintf(ps, "*NptMinTheta: %lf\n", npt_min_theta); fpr(fp);
    sprintf(ps, "*PrintAllMaps: %d\n", print_all_maps); fpr(fp);

    sprintf(ps, "*UseErrorRate: %d\n", use_error_rate); fpr(fp);
    sprintf(ps, "*ErrorLodThreshold: %lf\n", error_lod_thresh); fpr(fp);
    sprintf(ps, "*ErrorSingleThreshold: %lf\n", error_single_thresh); fpr(fp);
    sprintf(ps, "*ErrorNetThreshold: %lf\n", error_net_thresh); fpr(fp);

    sprintf(ps, "*Contexts: %d\n", num_contexts); fpr(fp);
    sprintf(ps, "*ActiveContext: %d\n", active_context); fpr(fp);
    fnl(fp);

    for (i=0; i<num_contexts; i++) {
      sprintf(ps, "*Context %d\n", i + 1); fpr(fp);
      sprintf(ps, "*SexSpecific: %d\n", context[i]->sex_specific); fpr(fp);
      sprintf(ps, "*CompressDNA: %d\n", context[i]->compress_DNA); fpr(fp);
      if (raw.data_type==CEPH) usenum=context[i]->use_number; else usenum=0;
      sprintf(ps, "*UseNum: %d\n", usenum); fpr(fp);
      sprintf(ps, "*MaxProblemSize: %ld\n", context[i]->max_problem_size); fpr(fp);
      
      sprintf(ps, "*SavedNames: %d\n",
              count_table_entries(context[i]->named_sequences)); fpr(fp);
      save_table(context[i]->named_sequences,fp,INDEX_BY_NAME);
	    
      sprintf(ps, "*SequenceHistory: %d\n",
              count_table_entries(context[i]->sequence_history)); fpr(fp);
      save_table(context[i]->sequence_history,fp,INDEX_BY_NUMBER);
    }
}


void read_status(fp)
FILE *fp;
{
    int i, mapnum, usenum, num;
    char word[TOKLEN+1];

    fscanf(fp,"%30s\n",word); 
    if (!streq(word,"*MapmakerStatusInfo:")) send(CRASH);

    fscanf(fp,"%*s %d\n", &print_names);
    fscanf(fp,"%*s %lf\n",&tolerance);
    fscanf(fp,"%*s %d\n", &units);
    fscanf(fp,"%*s %d\n", &mapnum); map_func(mapnum);

    fscanf(fp,"%*s %lf\n",&default_lod);
    fscanf(fp,"%*s %lf\n",&default_theta);

    fscanf(fp,"%*s %d\n", &use_three_pt);
    fscanf(fp,"%*s %lf\n",&triplet_lod);
    fscanf(fp,"%*s %lf\n",&triplet_theta);
    fscanf(fp,"%*s %d\n", &triplet_num_links);
    fscanf(fp,"%*s %lf\n",&three_pt_threshold);
    fscanf(fp,"%*s %d\n", &three_pt_window);
    fscanf(fp,"%*s %lf\n",&triplet_error_rate);

    fscanf(fp,"%*s %lf\n",&npt_threshold);
    fscanf(fp,"%*s %lf\n",&npt_first_threshold);
    fscanf(fp,"%*s %d\n", &npt_window);
    fscanf(fp,"%*s %d\n", &npt_min_indivs);
    fscanf(fp,"%*s %d\n", &npt_codominant);
    fscanf(fp,"%*s %lf\n",&npt_min_theta);
    fscanf(fp,"%*s %d\n", &print_all_maps);

    fscanf(fp,"%*s %d\n", &use_error_rate);
    fscanf(fp,"%*s %lf\n",&error_lod_thresh);
    fscanf(fp,"%*s %lf\n",&error_single_thresh);
    fscanf(fp,"%*s %lf\n",&error_net_thresh);

    fscanf(fp,"%20s %d\n",word,&num_contexts);
    if (!streq(word,"*Contexts:")) send(CRASH);
    fscanf(fp,"%*s %d",&active_context);

    for (i=0; i<num_contexts; i++) { /* Fix this... */
	fscanf(fp,"%*s %*d\n");
	fscanf(fp,"%*s %d\n",&sex_specific);
	context[i]->sex_specific= sex_specific;
	fscanf(fp,"%*s %d\n",&compress_DNA);
	context[i]->compress_DNA = compress_DNA;
	fscanf(fp,"%*s %d\n",&usenum);

#ifdef HAVE_CEPH
	if (raw.data_type==CEPH) {
	    raw.data.ceph.use_number=usenum;
	    context[i]->use_number=usenum;
	}
#endif
	fscanf(fp,"%*s %ld\n",&max_problem_size);
	context[i]->max_problem_size=max_problem_size;
	fscanf(fp,"%*s %d\n",&num);  /* "Saved names: %d" header */
	load_table(context[i]->named_sequences,fp,INDEX_BY_NAME,num);
	fscanf(fp,"%*s %d\n",&num); context[i]->seq_history_num=num; 
	load_table(context[i]->sequence_history,fp,INDEX_BY_NUMBER,num);
    }

    if (the_seq_history_num>0) set_current_seq("none",FALSE);
}
