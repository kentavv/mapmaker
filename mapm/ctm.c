/******************************************************************************

  ####    #####  #    #           ####
 #    #     #    ##  ##          #    #
 #          #    # ## #          #
 #          #    #    #   ###    #
 #    #     #    #    #   ###    #    #
  ####      #    #    #   ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

/************** Hiden Markov Model based Converge to Map for F2 **************/

#define INC_LIB
#define INC_SHELL
#include "mapm.h"

/**** statics all setup by setup_hmm() - arrays are freed by free_hmm() ****/
int n_indivs, n_loci, n_intervals;
int max_observations, max_states;
int  **n_poss_states, ***poss_state;     /* [indiv][locus],[indiv][locus][#] */
real ***trans_prob, ***exp_transitions;  /* [interval][from_state][to_state] */
real ***implied_recs, ***implied_norecs; /* [interval][from_state][to_state] */
real ***obs_prob, ***exp_observations;   /* [locus][observation/MD][state] */
bool *obs_prob_fixed;  /* [locus] - if all are fixed then exp_obs==NULL */
int  **observation;    /* [indiv][locus] => OBS_ code */
int  *the_observation; /* [indiv] a temp used only in f2_setup_hmm */
real **apriori_prob;   /* [indiv][state] - for leftmost locus */
int  no_data;          /* index [indiv] into apriori_prob for right locus */
real *e_step_temp;     /* [state]- used in LCP and RCP calculations */
real ***exp_genotype;  /* [indiv][locus][state] for ERROR-KLUDGE */

int *observations; /* [indiv] externed for use as a temp */
void hmm_converge_to_map();

void (*hmm_set_probs)();
bool error_kludge; /* TRUE for ERROR-KLUDGE, e.g. error_rate>0 */
bool separate_recs; 
/* TRUE if distinct implied_rec/norec entries for each interval */

typedef struct {
    real recs;   /* for now */
    real norecs;
} SUFF_STATS;
SUFF_STATS **sufficient_stats; /* [locus] => ptr to struct */

/* Internals for HMM E-Step (obs stuff moved to toplevel,h) */
real **left_prob, **right_prob;  /* [locus#][state#] */
real ***indiv_transitions;       /* [interval][from_state][to_state] */
real ***indiv_observations;      /* [locus][observation][state], maybe NULL */

int obligate_recs[6][6]= { /* externed in toplevel.h */ 
    /* -  A  H  B  C  D */
    {  0, 0, 0, 0, 0, 0 }, /* - */
    {  0, 0, 1, 2, 1, 0 }, /* A */
    {  0, 1, 0, 1, 0, 0 }, /* H */
    {  0, 2, 1, 0, 0, 1 }, /* B */
    {  0, 1, 0, 0, 0, 0 }, /* C */
    {  0, 0, 1, 0, 0, 0 }  /* D */
};

#define HMM_MAX_THETA 0.4995
#define HMM_MIN_THETA 0.0004

/* F2 backcross/intercross, RI states */
#define STATE_A 0
#define STATE_B 1 
#define STATE_H 2   /* intercross only */

int ri_type;  /* RI_SELF, RI_SIB, or RI_NOT */


/**** F3 specific stuff Follows: ****/

#define A 0  /* F2 States... */
#define H 1
#define J 2
#define B 3

/* States for one F3 chromosome (grand- maternal or paternal origin) */
#define M 0
#define P 1

/* There are 4 F2 states and 16 hidden F2/F3 states: */
int A_AMAM, A_AMAP, A_APAM, A_APAP, H_AMAM, H_AMBP, H_BPAM, H_BPBP;
int J_BMBM, J_BMAP, J_APBM, J_APAP, B_BMBM, B_BMBP, B_BPBM, B_BPBP;
/* These are temps for setup_f3_self_hmm(), although f2_recs is also 
   used in hmm_f3_self_set_probs() */
int f3_states, mat_chrom[16], pat_chrom[16], f2_state[16];
real f2_recs[4][4];
int  f3_state();

/**** The internal routines ****/

real hmm_e_step();
void hmm_count_stats(), hmm_make_new_map();
void setup_hmm();

void setup_bc_like_hmm(), hmm_bc_like_set_probs();
void setup_f2_hmm(), hmm_f2_set_probs();
void setup_f3_self_hmm(), hmm_f3_self_set_probs();

void test_dump(), test_dump0();
real null_like;
bool twopt_kludge;
MAP *map2;
void count_error_lods();
void hmm_fake_converge_to_map();


/**** External Routines, mostly ****/

void 
allocate_hmm_temps (int total_loci, int num_indivs, int cross_type)
{
    int i, num_states=0, num_observations=0, num_loci, num_intervals;
    bool special_recs_hack;
    
    switch (cross_type) {
	case F2_BACKCROSS:  num_states=2;  num_observations=3; break;
	case F2_INTERCROSS: num_states=3;  num_observations=6; break;
	case RI_SELF:       num_states=2;  num_observations=3; break;
	case RI_SIB:        num_states=2;  num_observations=3; break;
	case F3_SELF:       num_states=16; num_observations=6; break;
	default: send(CRASH);
    }
    num_loci=MAX_MAP_LOCI; num_intervals=num_loci-1;
    special_recs_hack= (cross_type=F2_INTERCROSS ? TRUE:FALSE);

    parray(sufficient_stats,num_loci,SUFF_STATS);
    matrix(left_prob,num_loci,num_states,real);
    matrix(right_prob,num_loci,num_states,real);
    indiv_observations=NULL;

    array(trans_prob,num_intervals,real**);
    array(exp_transitions,num_intervals,real**);
    array(indiv_transitions,num_intervals,real**);
    array(implied_recs,num_intervals,real**);
    array(implied_norecs,num_intervals,real**);
    for (i=0; i<num_intervals; i++) {
	matrix(trans_prob[i],num_states,num_states,real);
	matrix(exp_transitions[i],num_states,num_states,real);
	matrix(indiv_transitions[i],num_states,num_states,real);
	if (special_recs_hack || i==0) { 
	    matrix(implied_recs[i],num_states,num_states,real);
	    matrix(implied_norecs[i],num_states,num_states,real);
	}
    }
    if (!special_recs_hack) for (i=1; i<num_intervals; i++) {
	/*  Do the funky pointer KLUDGE */
	implied_recs[i]= implied_recs[0];
	implied_norecs[i]= implied_norecs[0];
    }

    array(obs_prob_fixed,num_loci,bool);
    array(obs_prob,num_loci,real**);
    for (i=0; i<num_loci; i++) {
	matrix(obs_prob[i],num_observations,num_states,real);
    }

    matrix(n_poss_states,num_indivs,num_loci,int);
    array(poss_state,num_indivs,int**);
    for (i=0; i<num_indivs; i++) {
	matrix(poss_state[i],num_loci,num_states,int);
    }

    exp_observations=NULL;  /* unused at present */

    matrix(observation,num_indivs,num_loci,int);
    array(the_observation,num_indivs,int);
    matrix(apriori_prob,num_indivs+1,num_states,real);
    array(e_step_temp,num_states,real);

    array(exp_genotype,num_indivs,real**);  
    for (i=0; i<num_indivs; i++) {
	matrix(exp_genotype[i],num_loci,num_states,real);
    }
    map2= allocate_map(2); /* should have NULL error_lod stuff */
    array(observations,num_indivs,int); /* extern */
}


void 
free_hmm_temps (int total_loci, int num_indivs, int cross_type)
{
    int i, num_states=0, num_observations=0, num_loci, num_intervals;
    bool special_recs_hack;
    
    switch (cross_type) {
	case F2_BACKCROSS:  num_states=2;  num_observations=3; break;
	case F2_INTERCROSS: num_states=3;  num_observations=6; break;
	case RI_SELF:       num_states=2;  num_observations=3; break;
	case RI_SIB:        num_states=2;  num_observations=3; break;
	case F3_SELF:       num_states=16; num_observations=6; break;
	default: send(CRASH);
    }
    num_loci=MAX_MAP_LOCI; num_intervals=num_loci-1;
    special_recs_hack= (cross_type=F2_INTERCROSS ? TRUE:FALSE);

    unparray(sufficient_stats,num_loci,SUFF_STATS);
    unmatrix(left_prob,num_loci,real);
    unmatrix(right_prob,num_loci,real);
    
    for (i=0; i<num_intervals; i++) {
	unmatrix(trans_prob[i],num_states,real);
	unmatrix(exp_transitions[i],num_states,real);
	unmatrix(indiv_transitions[i],num_states,real);
	if (special_recs_hack || i==0) {
	    unmatrix(implied_recs[i],num_states,real);
	    unmatrix(implied_norecs[i],num_states,real);
	}
    }
    unarray(trans_prob,real**);
    unarray(exp_transitions,real**);
    unarray(indiv_transitions,real);
    unarray(implied_recs,real**);
    unarray(implied_norecs,real**);
    /* indiv_observations was never allocated */

    unarray(obs_prob_fixed,bool);
    for (i=0; i<num_loci; i++) {
	unmatrix(obs_prob[i],num_observations,real); 
    }
    unarray(obs_prob,real**);

    unmatrix(n_poss_states,num_indivs,int);
    for (i=0; i<num_indivs; i++) {
	unmatrix(poss_state[i],num_loci,int);
    }
    unarray(poss_state,int);

    /* exp_observations was never allocated */

    unmatrix(observation,num_indivs,int);
    unarray(the_observation,int);
    unmatrix(apriori_prob,num_indivs+1,real);
    unarray(e_step_temp,real);

    for (i=0; i<num_indivs; i++) {  
	unmatrix(exp_genotype[i],num_loci,real);
    }
    unarray(exp_genotype,real**);  
    free_map(map2); /* should have NULL error_lod stuff */
    unarray(observations,int); /* extern */
}


void 
converge_to_map (MAP *map)
{
    if (map==NULL || map->num_loci<2 || raw.data_type!=F2) send(CRASH);
    if (map->num_loci>MAX_MAP_LOCI) 
      { strcpy(ps,"too many loci in one map"); error(ps); }

    setup_hmm(map); /* set statics above, and set the initial parameters */
    if (fake_maps) hmm_fake_converge_to_map(map);
    else hmm_converge_to_map(map);
    /* else old_converge_to_map(map); */
}


void f2_genotype(locus,haplo,observation)
int locus;
bool haplo; /* merge haplotypes */
int *observation; /* array of raw.num_indivs observations */
{
    int j;
    
    if (raw.data_type!=F2) send(CRASH);

    if (haplo && haplotyped(locus))
      locus=haplo_first[locus];  /* go to the head of the list */

    for (j=0; j<raw.data.f2.num_indivs; j++)
      switch (raw.data.f2.allele[locus][j]) {
	  case PARENTAL_TYPE_A:	observation[j]=OBS_A; break;
	  case PARENTAL_TYPE_B:	observation[j]=OBS_B; break;
	  case HYBRID_TYPE_H:	observation[j]=OBS_H; break;
	  case TYPE_NOT_A:	observation[j]=OBS_NOT_A; break;
	  case TYPE_NOT_B:	observation[j]=OBS_NOT_B; break;
	  default: 		observation[j]=OBS_MISSING; break;
      }
    
    if (haplo && haplo_first[locus]!=NO_LOCUS) 
      while ((locus=haplo_next[locus])!=NO_LOCUS)
	merge_genotypes(locus,observation,observation);
}


bool merge_genotypes(locus,observation,new_observation)
int locus;
int *observation, *new_observation; /* array of raw.num_indivs observations */
{
    int j, the_obs;
    
    for (j=0; j<raw.data.f2.num_indivs; j++) {
	switch (raw.data.f2.allele[locus][j]) {
	    case PARENTAL_TYPE_A:	the_obs=OBS_A; break;
	    case PARENTAL_TYPE_B:	the_obs=OBS_B; break;
	    case HYBRID_TYPE_H:		the_obs=OBS_H; break;
	    case TYPE_NOT_A:		the_obs=OBS_NOT_A; break;
	    case TYPE_NOT_B:		the_obs=OBS_NOT_B; break;
	    default: 			the_obs=OBS_MISSING; break;
	}
	if (obligate_recs[observation[j]][the_obs]) 
	  return(FALSE);
	if ((observation[j]==OBS_MISSING && the_obs!=OBS_MISSING) ||
	    (!fully_inf(observation[j]) && fully_inf(the_obs)))
	  new_observation[j]=the_obs;
	else new_observation[j]= observation[j];
    }
    return(TRUE);
}


int 
f2_count_infs ( /* return #inf */
    int *num_dom,
    int *num_het, /* side-effected */
    int *observation /* array of raw.num_indivs observations */
)
{
    int homo, het, dom, missing, j;

    homo= het= dom= missing= 0;
    for (j=0; j<raw.data.f2.num_indivs; j++)
      switch (observation[j]) {
	  case OBS_A: homo++; break;
	  case OBS_B: homo++; break;
	  case OBS_H: het++; break;
	  case OBS_NOT_A: dom++; break;
	  case OBS_NOT_B: dom++; break;
	  default: missing++; break;
      }
    *num_dom= dom; *num_het= het;
    return(raw.data.f2.num_indivs - missing);
}


/**** The real Mapmaker itself ****/

void 
hmm_converge_to_map (MAP *map)
{
    int iter;
    int testing, one;
    real old_like, new_like;

    old_like= new_like= 0.0; iter=0;
    testing=0; one=(n_loci>=3 ? 1:0); /* just debugging vars */
    for (iter=0; TRUE; ++iter) {

	(*hmm_set_probs)(map,trans_prob,obs_prob,implied_recs);
	old_like= new_like;
	new_like= hmm_e_step(map,trans_prob,obs_prob,exp_transitions,
			     exp_observations);
	if (n_loci==2 && iter==0) null_like=new_like;

	/* sprintf(ps,"iter: %2d  like=%lf  thetas: %lf %lf\n",iter,new_like,
	   map->rec_frac[0][0],map->rec_frac[one][0]); pr(); */
if (testing) { /********** DEBUGGING CODE **********/
	test_dump(obs_prob,max_observations,"obs_prob");
	test_dump(implied_recs,max_states,"implied_recs");
	test_dump(trans_prob,max_states,"trans_prob");
	test_dump(exp_transitions,max_states,"from->to exp_trans: "); 
} /********** END OF DEBUGGING CODE **********/

	if (fabs(new_like - old_like) < tolerance) break; /* converged! */
	hmm_count_stats(exp_transitions,exp_observations,implied_recs,
			   sufficient_stats);

if (testing) { /********** MORE DEBUGGING CODE **********/
#define ss sufficient_stats
	sprintf(ps, "ss: R/NR %lf N %lf | %lf %lf\n", ss[0]->recs, ss[0]->norecs,
            ss[one]->recs, ss[one]->norecs); pr();
} /********** END OF DEBUGGING CODE **********/

	hmm_make_new_map(sufficient_stats,new_like,map);

    } /* iterate */

    if (error_kludge)
      count_error_lods(map->error_rate,obs_prob,exp_genotype,
		       map->error_lod);
/*    if (n_loci==2) {
	map->rec_frac[0][0]= 0.499999; 
	(*hmm_set_probs)(map,trans_prob,obs_prob,implied_recs);
	null_like= hmm_e_step(map,trans_prob,obs_prob,exp_transitions,
				  exp_observations);
    } */    
}


void 
setup_hmm (MAP *map)
{
    int type;

    n_loci=map->num_loci; 
    n_intervals=n_loci-1;
    n_indivs=raw.data.f2.num_indivs;
    no_data=raw.data.f2.num_indivs; /* extra index into apriori_prob */
    error_kludge= twopt_kludge= FALSE;
    ri_type=RI_NOT;
    type=raw.data.f2.cross_type;

    if (type==F2_INTERCROSS) setup_f2_hmm(map->locus,map->error_rate);
    else if (type==F3_SELF)  setup_f3_self_hmm(map->locus);
    else if (type==F2_BACKCROSS || type==RI_SELF || type==RI_SIB)
      setup_bc_like_hmm(map->locus,map->error_rate,type);
    else send(CRASH);
}
    

void 
setup_bc_like_hmm (
    int *locus, /* The locus numbers in order */
    double *error_rate,
    int cross_type
)
{
    int i, j, k, obs;
    real **recs, **norecs, error_prob;

    max_states=2; 	/* STATE_A or _B */
    max_observations=3; /* OBS_A, _H, or _MISSING (yup, this is odd for BC) */
    separate_recs=FALSE;
    hmm_set_probs=hmm_bc_like_set_probs;
    if (error_rate!=NULL) error_kludge=TRUE;

    if (cross_type==F2_BACKCROSS) ri_type=RI_NOT;
    else if (cross_type==RI_SELF) ri_type=RI_SELF;
    else if (cross_type==RI_SIB)  ri_type=RI_SIB;
    else send(CRASH);

    error_prob=0.0;
    for (i=0; i<n_loci; i++) {
	if (error_kludge) error_prob= error_rate[i]; /* else =0.0 */
	obs_prob[i][OBS_A][STATE_A]= 1.0 - error_prob;
	obs_prob[i][OBS_H][STATE_A]= 0.0 + error_prob;
	obs_prob[i][OBS_A][STATE_B]= 0.0 + error_prob;
	obs_prob[i][OBS_H][STATE_B]= 1.0 - error_prob;
	obs_prob[i][OBS_MISSING][STATE_A]= 1.0; 
	obs_prob[i][OBS_MISSING][STATE_B]= 1.0;
	obs_prob_fixed[i]=TRUE;

	f2_genotype(locus[i],use_haplotypes,the_observation);
	for (j=0; j<n_indivs; j++) {
	    obs=observation[j][i]=the_observation[j];
	    k=0;
	    if (obs==OBS_MISSING || error_prob>0.0) {
		poss_state[j][i][k++]=STATE_A;
		poss_state[j][i][k++]=STATE_B;
	    } else if (obs==OBS_A) {
		poss_state[j][i][k++]=STATE_A;
	    } else if (obs==OBS_H) {
		poss_state[j][i][k++]=STATE_B;
	    } else send(CRASH);
	    n_poss_states[j][i]=k;
	}
    }

    for (i=0; i<n_intervals; i++) {
	/* set the implied_recs */
	recs= implied_recs[i]; norecs= implied_norecs[i];
	recs[STATE_A][STATE_A]=0.0; norecs[STATE_A][STATE_A]=1.0;
	recs[STATE_A][STATE_B]=1.0; norecs[STATE_A][STATE_B]=0.0;
	recs[STATE_B][STATE_A]=1.0; norecs[STATE_B][STATE_A]=0.0;
	recs[STATE_B][STATE_B]=0.0; norecs[STATE_B][STATE_B]=1.0; 
    } /* end loop over intervals (i) */

    for (j=0; j<n_indivs; j++) {
	obs= observation[j][0];
	if (obs!=OBS_MISSING) {
	    apriori_prob[j][STATE_A]=0.50 * obs_prob[0][obs][STATE_A];
	    apriori_prob[j][STATE_B]=0.50 * obs_prob[0][obs][STATE_B];
	} else {
	    apriori_prob[j][STATE_A]=0.50;
	    apriori_prob[j][STATE_B]=0.50;
	}
    }
    apriori_prob[no_data][STATE_A]=0.50;
    apriori_prob[no_data][STATE_B]=0.50;
}


void 
setup_f2_hmm (
    int *locus, /* The locus numbers in order */
    double *error_rate
) 
{
    int i, j, k, obs;
    real total, **recs, **norecs, error_prob;

    max_states=3; 	/* STATE_A,_H, or _B */
    max_observations=6; /* OBS_A, _B, _H, _NOT_A, _NOT_B, or _MISSING */
    separate_recs=TRUE;
    hmm_set_probs=hmm_f2_set_probs;
    if (error_rate!=NULL) error_kludge=TRUE;

    error_prob=0.0;
    for (i=0; i<n_loci; i++) {
	if (error_kludge) error_prob= error_rate[i]; /* else =0.0 */
	obs_prob[i][OBS_A][STATE_A]= 1.0 - error_prob;
	obs_prob[i][OBS_B][STATE_A]= 0.0 + error_prob/2.0;
	obs_prob[i][OBS_H][STATE_A]= 0.0 + error_prob/2.0;

	/* Fix the NOT_A/B cases */
	obs_prob[i][OBS_NOT_A][STATE_A]= 0.0 + (2.0 * error_prob); /* KLUDGE */
	obs_prob[i][OBS_NOT_B][STATE_A]= 1.0 - error_prob;

	obs_prob[i][OBS_A][STATE_H]= 0.0 + error_prob/2.0;
	obs_prob[i][OBS_B][STATE_H]= 0.0 + error_prob/2.0;
	obs_prob[i][OBS_H][STATE_H]= 1.0 - error_prob;
	obs_prob[i][OBS_NOT_A][STATE_H]= 1.0 - error_prob;;
	obs_prob[i][OBS_NOT_B][STATE_H]= 1.0 - error_prob;;

	obs_prob[i][OBS_A][STATE_B]= 0.0 + error_prob/2.0;
	obs_prob[i][OBS_B][STATE_B]= 1.0 - error_prob;
	obs_prob[i][OBS_H][STATE_B]= 0.0 + error_prob/2.0;
	obs_prob[i][OBS_NOT_A][STATE_B]= 1.0 - (2.0 * error_prob);;
	obs_prob[i][OBS_NOT_B][STATE_B]= 0.0 + error_prob;

	obs_prob[i][OBS_MISSING][STATE_A]= 1.0; 
	obs_prob[i][OBS_MISSING][STATE_B]= 1.0;
	obs_prob[i][OBS_MISSING][STATE_H]= 1.0;

	obs_prob_fixed[i]=TRUE;

	/* make a copy of the raw data, translated for speed. Now also
	   adds the haplotype info, when available - note order
	   of indecies for the observation matrix is opposite that 
	   for the raw struct! */

	f2_genotype(locus[i],use_haplotypes,the_observation);
	for (j=0; j<n_indivs; j++) {
	    obs=observation[j][i]=the_observation[j];
	    k=0;
	    if (obs==OBS_MISSING || error_prob>0.0) {
		poss_state[j][i][k++]= STATE_A;
		poss_state[j][i][k++]= STATE_B;
		poss_state[j][i][k++]= STATE_H;
	    } else if (obs==OBS_NOT_A) {
		poss_state[j][i][k++]= STATE_B;
		poss_state[j][i][k++]= STATE_H;
	    } else if (obs==OBS_NOT_B) {
		poss_state[j][i][k++]= STATE_A;
		poss_state[j][i][k++]= STATE_H;
	    } else if (obs==OBS_A) {
		poss_state[j][i][k++]= STATE_A;
	    } else if (obs==OBS_B) {
		poss_state[j][i][k++]= STATE_B;
	    } else if (obs==OBS_H) {
		poss_state[j][i][k++]= STATE_H;
	    } else send(CRASH);
	    n_poss_states[j][i]=k;
	}
    }

    for (i=0; i<n_intervals; i++) {
	/* set most of the implied_recs, only the [H][H] entries vary */
	recs= implied_recs[i]; norecs= implied_norecs[i];
	recs[STATE_A][STATE_A]= 0.0; norecs[STATE_A][STATE_A]= 2.0; 
	recs[STATE_A][STATE_B]= 2.0; norecs[STATE_A][STATE_B]= 0.0;
	recs[STATE_A][STATE_H]= 1.0; norecs[STATE_A][STATE_H]= 1.0;
	recs[STATE_B][STATE_A]= 2.0; norecs[STATE_B][STATE_A]= 0.0;
	recs[STATE_B][STATE_B]= 0.0; norecs[STATE_B][STATE_B]= 2.0;
	recs[STATE_B][STATE_H]= 1.0; norecs[STATE_B][STATE_H]= 1.0;
	recs[STATE_H][STATE_A]= 1.0; norecs[STATE_H][STATE_A]= 1.0;
	recs[STATE_H][STATE_B]= 1.0; norecs[STATE_H][STATE_B]= 1.0;
	recs[STATE_H][STATE_H]= 0.0; norecs[STATE_H][STATE_H]= 0.0; 
	/* _H _H entries are a function of theta[i] and are changed 
	   each iteration by set_probs() */
    } /* end loop over intervals (i) */

    /* Lastly, we fill in apriori state probs for leftmost locus only:
       THIS HAS BEEN CHANGED TO CALCULATE THE UN-NORMALIZED PRIORS! */

    for (j=0; j<n_indivs; j++) {
	obs= observation[j][0]; total=0.0;
	if (obs != OBS_MISSING) { /* NOT MISSING  - general version??? */

#ifdef HMMMM
	    total+= apriori_prob[j][STATE_A]= obs_prob[0][obs][STATE_A];
	    total+= apriori_prob[j][STATE_B]= obs_prob[0][obs][STATE_B];
	    total+= apriori_prob[j][STATE_H]= obs_prob[0][obs][STATE_H];
#else
	    total+= apriori_prob[j][STATE_A]= 0.25 * obs_prob[0][obs][STATE_A];
	    total+= apriori_prob[j][STATE_B]= 0.25 * obs_prob[0][obs][STATE_B];
	    total+= apriori_prob[j][STATE_H]= 0.50 * obs_prob[0][obs][STATE_H];
#endif

/*	    apriori_prob[j][STATE_A]/= total;
	    apriori_prob[j][STATE_B]/= total;
	    apriori_prob[j][STATE_H]/= total; */

	} else { /* OBSERVATION IS MISSING */
	    apriori_prob[j][STATE_A]=0.25;
	    apriori_prob[j][STATE_B]=0.25;
	    apriori_prob[j][STATE_H]=0.50;
	}
    }
    apriori_prob[no_data][STATE_A]=0.25;
    apriori_prob[no_data][STATE_B]=0.25;
    apriori_prob[no_data][STATE_H]=0.50;
}


void 
setup_f3_self_hmm (
    int *locus /* The locus numbers in the sequence */
)
{
    int i, j, k, the_obs;
    real recs;
    int a[3], b[3], h[3], n[3], aa, ab, ah, ba, bb, bh, ha, hb, hh;
    
    max_states=16; f3_states=0;	 
    max_observations=6; /* OBS_A, _B, _H, _NOT_A, _NOT_B, or _MISSING */
    separate_recs=FALSE;
    hmm_set_probs=hmm_f3_self_set_probs;

    A_AMAM= f3_state("A-AmAm",A,M,M,OBS_A); /* 0  */
    A_AMAP= f3_state("A-AmAp",A,M,P,OBS_A); /* 1  */
    A_APAM= f3_state("A-ApAm",A,P,M,OBS_A); /* 2  */
    A_APAP= f3_state("A-ApAp",A,P,P,OBS_A); /* 3  */

    H_AMAM= f3_state("H-AmAm",H,M,M,OBS_A); /* 4  */
    H_AMBP= f3_state("H-AmBp",H,M,P,OBS_H); /* 5  */
    H_BPAM= f3_state("H-BpAm",H,P,M,OBS_H); /* 6  */
    H_BPBP= f3_state("H-BpBp",H,P,P,OBS_B); /* 7  */

    J_BMBM= f3_state("J-BmBm",J,M,M,OBS_B); /* 8  */
    J_BMAP= f3_state("J-BmAp",J,M,P,OBS_H); /* 9  */
    J_APBM= f3_state("J-ApBm",J,P,M,OBS_H); /* 10 */
    J_APAP= f3_state("J-ApAp",J,P,P,OBS_A); /* 11 */

    B_BMBM= f3_state("B-BmBm",B,M,M,OBS_B); /* 12 */
    B_BMBP= f3_state("B-BmBp",B,M,P,OBS_B); /* 13 */
    B_BPBM= f3_state("B-BpBm",B,P,M,OBS_B); /* 14 */
    B_BPBP= f3_state("B-BpBp",B,P,P,OBS_B); /* 15 */
    
    for (i=0; i<2; i++) a[i]=b[i]=h[i]=n[i]=0;
    aa=ab=ah=ba=bb=bh=ha=hb=hh=0;
    for (i=0; i<n_loci; i++) {
	/* Here we have a KLUDGE, in that the poss_state and n_poss_state 
	   vars could be used to reduce the computation. However, as yet
	   we don't */

	obs_prob_fixed[i]= TRUE;

	/* Make a copy of the raw data, translated for speed
	   KLUDGE update to use f2_genotype() */

	for (j=0; j<n_indivs; j++) {
	    n_poss_states[j][i]= f3_states; /* the slow way */
	    for (k=0; k<f3_states; k++) poss_state[j][i][k]=k;
	    switch (raw.data.f2.allele[locus[i]][j]) {
		case PARENTAL_TYPE_A:	observation[j][i]=OBS_A; break;
		case PARENTAL_TYPE_B:	observation[j][i]=OBS_B; break;
		case HYBRID_TYPE_H:	observation[j][i]=OBS_H; break;
		case TYPE_NOT_A:	observation[j][i]=OBS_NOT_A; break;
		case TYPE_NOT_B:	observation[j][i]=OBS_NOT_B; break;
		case MISSING_DATA:	observation[j][i]=OBS_MISSING; break;
		default:		send(CRASH);
	    }
	    if (i<2) { 
		if (observation[j][i]==OBS_A) { a[i]+=1; n[i]+=1; }
		if (observation[j][i]==OBS_B) { b[i]+=1; n[i]+=1; }
		if (observation[j][i]==OBS_H) { h[i]+=1; n[i]+=1; }
	    }
	    if (i==0) {
if (observation[j][i]==OBS_A && observation[j][i+1]==OBS_A) aa++;
if (observation[j][i]==OBS_A && observation[j][i+1]==OBS_B) ab++;
if (observation[j][i]==OBS_A && observation[j][i+1]==OBS_H) ah++;
if (observation[j][i]==OBS_B && observation[j][i+1]==OBS_A) ba++;
if (observation[j][i]==OBS_B && observation[j][i+1]==OBS_B) bb++;
if (observation[j][i]==OBS_B && observation[j][i+1]==OBS_H) bh++;
if (observation[j][i]==OBS_H && observation[j][i+1]==OBS_A) ha++;
if (observation[j][i]==OBS_H && observation[j][i+1]==OBS_B) hb++;
if (observation[j][i]==OBS_H && observation[j][i+1]==OBS_H) hh++;
	    }
	} /* end loop over indivs (j) */
    } /* end loop over loci (i) */

/*    printf("\n\nL: A=%4d  H=%4d  B=%4d\nR: A=%4d  H=%4d  B=%4d\n\n",
	   a[0],h[0],b[0],a[1],h[1],b[1]); printf(
"AA=%4d  AH=%4d  AB=%4d\nHA=%4d  HH=%4d  HB=%4d\nBA=%4d  BH=%4d  BB=%4d\n\n",
	   aa,ah,ab,ha,hh,hb,ba,bh,bb); */

    /* For F3, the implied recs entries are all constant. Thus, we do
       the same pointer KLUDGE we do for the obs_probs to save space. */

    /* Number of recs in the F2 generation: norecs= 2-recs */
    f2_recs[A][A]=0.0; f2_recs[H][A]=1.0; f2_recs[J][A]=1.0; f2_recs[B][A]=2.0;
    f2_recs[A][H]=1.0; f2_recs[H][H]=0.0; f2_recs[J][H]=2.0; f2_recs[B][H]=1.0;
    f2_recs[A][J]=1.0; f2_recs[H][J]=2.0; f2_recs[J][J]=0.0; f2_recs[B][J]=1.0;
    f2_recs[A][B]=2.0; f2_recs[H][B]=1.0; f2_recs[J][B]=1.0; f2_recs[B][B]=0.0;

    for (i=0; i<f3_states; i++) for (j=0; j<f3_states; j++) {
	implied_recs[0][i][j]= recs=
	    f2_recs[f2_state[i]][f2_state[j]] +
	    (mat_chrom[i]==mat_chrom[j] ? 0.0:1.0) +
	    (pat_chrom[i]==pat_chrom[j] ? 0.0:1.0);
	implied_norecs[0][i][j]= 4.0 - recs;
    }

    /* Lastly, we fill in apriori state probs for leftmost locus only */
    /* In F3 HMM, we enumerate every possible situtation: Thus, 4 F2 genotypes
       each can give rise to 4 F3-Self genotypes. Therefore, the apriori for
       any state is 1/4 * 1/4. These are not normalized. */

    for (j=0; j<n_indivs; j++) {
	the_obs= observation[j][0];
	if (the_obs==OBS_MISSING) {
	    for (k=0; k<f3_states; k++) 
	      apriori_prob[j][k]=0.0625;
	} else { /* not OBS_MISSING */
	    for (k=0; k<f3_states; k++) 
	      apriori_prob[j][k]=0.0625 * obs_prob[0][the_obs][k];
	}
    }
    for (k=0; k<f3_states; k++) apriori_prob[no_data][k]=0.0625;
}


int 
f3_state (
    char *name, /* for debugging */
    int f2,
    int mat,
    int pat,
    int obs
)
{
    int i,j;
    i= f3_states++; /* assign f3 state number, f3_states is a global */
    f2_state[i]= f2;
    mat_chrom[i]= mat;
    pat_chrom[i]= pat;

    /* Assume the KLUDGE allowing us to only have to set obs_prob[0] */
    /* well, actually, let's not and see what happens */

    for (j=0; j<n_loci; j++) {
      obs_prob[j][OBS_A][i]= (obs==OBS_A ? 1.0:0.0);
      obs_prob[j][OBS_B][i]= (obs==OBS_B ? 1.0:0.0);
      obs_prob[j][OBS_H][i]= (obs==OBS_H ? 1.0:0.0);
      obs_prob[j][OBS_NOT_A][i]= (obs==OBS_H || obs==OBS_B ? 1.0:0.0);
      obs_prob[j][OBS_NOT_B][i]= (obs==OBS_H || obs==OBS_H ? 1.0:0.0);
      obs_prob[j][OBS_MISSING][i]= 1.0;
    }

    return(i);
}


void hmm_bc_like_set_probs(map,trans_prob,obs_prob,implied_recs)
MAP *map;
real ***trans_prob;    /* [interval][from_state][to_state] */
real ***obs_prob;      /* [locus][observation/MD][state] side-effected */
real ***implied_recs;  /* [interval][from_state][to_state] side-effected */
/* for now, assume obs_probs are fixed and initialized by setup_hmm() */
{
    int i;
    real rec=0., norec, theta;
    real **prob;

    /* set the trans_probs (right geno, conditional on left) given theta */
    for (i=0; i<n_intervals; i++) {
	theta=map->rec_frac[i][0];
	if (ri_type==RI_NOT) rec=theta;
	  else if (ri_type==RI_SELF) rec=(2.0*theta)/(1.0+2.0*theta);
	  else if (ri_type==RI_SIB)  rec=(4.0*theta)/(1.0+6.0*theta);
	  else send(CRASH); /* From Haldane and Wadington (inverse of y(x)) */
	prob=trans_prob[i]; norec=1.0-rec;
	prob[STATE_A][STATE_A]=norec;
	prob[STATE_A][STATE_B]=rec;
	prob[STATE_B][STATE_A]=rec;
	prob[STATE_B][STATE_B]=norec;
    }

    /* For BC/RI we do not need to set any implied_recs/norecs - all entries
       are constant and are initialized by setup_bc1_hmm(). For RI this
       means that we are counting observed genotypic classes (as rec/norec)
       not really the implied_recs. See make_new_map. */
    /* Set obs_probs based on penetrance, etc. - not implemented yet
       exp_observations!=NULL means that obs_probs are not all fixed */
    if (exp_observations!=NULL) send(CRASH); 
}


void hmm_f2_set_probs(map,trans_prob,obs_prob,implied_recs)
MAP *map;
real ***trans_prob;    /* [interval][from_state][to_state] */
real ***obs_prob;      /* [locus][observation/MD][state] side-effected */
real ***implied_recs;  /* [interval][from_state][to_state] side-effected */
/* for now, assume obs_probs are fixed and initialized by setup_hmm() */
{
    int i;
    real theta, one_minus_theta, no_recs, one_rec, two_recs;
    real **prob, recs;

    /* set the trans_probs (right geno, conditional on left) given theta */
    for (i=0; i<n_intervals; i++) {
	theta= map->rec_frac[i][0]; one_minus_theta= 1.0-theta;
	no_recs= one_minus_theta * one_minus_theta;
	one_rec= theta * one_minus_theta;
	two_recs= theta * theta;
	prob= trans_prob[i]; 
	prob[STATE_A][STATE_A]= no_recs;	
	prob[STATE_A][STATE_B]= two_recs;
	prob[STATE_A][STATE_H]= 2.0 * one_rec;
	prob[STATE_B][STATE_A]= two_recs;
	prob[STATE_B][STATE_B]= no_recs;
	prob[STATE_B][STATE_H]= 2.0 * one_rec;
	prob[STATE_H][STATE_A]= one_rec;
	prob[STATE_H][STATE_B]= one_rec;
	prob[STATE_H][STATE_H]= no_recs + two_recs;

	/* as an optimization we only set implied_recs[H][H] - the 
	   other entries are constant and are initialized by setup_f2_hmm() */
	recs= (2.0 * two_recs)/(two_recs + no_recs);
	implied_recs[i][STATE_H][STATE_H]= recs;
	implied_norecs[i][STATE_H][STATE_H]= 2.0 - recs;
    }

    /* Set obs_probs based on penetrance, etc. - not implemented yet
       exp_observations!=NULL means that obs_probs are not all fixed */
    if (exp_observations!=NULL) send(CRASH); 
}


void hmm_f3_self_set_probs(map,trans_prob,obs_prob,implied_recs)
MAP *map;
real ***trans_prob;    /* [interval][from_state][to_state] */
real ***obs_prob;      /* [locus][observation/MD][state] side-effected */
real ***implied_recs;  /* [interval][from_state][to_state] side-effected */
/* for now, assume obs_probs are fixed and initialized by setup_hmm() */
{
    int i, k, j;
    real rec, norec, r0, r1, r2;
    real **prob;

    /* set the trans_probs (right geno, conditional on left) given theta */
    for (i=0; i<n_intervals; i++) {
	rec= map->rec_frac[i][0]; norec= 1.0 - rec;
	r0= norec * norec;
	r1= rec * norec;
	r2= rec * rec;

	/* f2_recs[][] is really the trans_prob for the F2 generation. */
	f2_recs[A][A]=r0; f2_recs[H][A]=r1; f2_recs[J][A]=r1; f2_recs[B][A]=r2;
	f2_recs[A][H]=r1; f2_recs[H][H]=r0; f2_recs[J][H]=r2; f2_recs[B][H]=r1;
	f2_recs[A][J]=r1; f2_recs[H][J]=r2; f2_recs[J][J]=r0; f2_recs[B][J]=r1;
	f2_recs[A][B]=r2; f2_recs[H][B]=r1; f2_recs[J][B]=r1; f2_recs[B][B]=r0;

	prob= trans_prob[i]; 
	for (j=0; j<f3_states; j++) for (k=0; k<f3_states; k++) 
	  prob[j][k]= f2_recs[f2_state[j]][f2_state[k]] *
	      (mat_chrom[j]==mat_chrom[k] ? norec:rec) *
	      (pat_chrom[j]==pat_chrom[k] ? norec:rec);

	/* For F3 we do not need to set any implied_recs/norecs - all entries
	   are constant and are initialized by setup_f3_self_hmm() */
    }

    /* Set obs_probs based on penetrance, etc. - not implemented yet
       exp_observations!=NULL means that obs_probs are not all fixed */
    if (exp_observations!=NULL) send(CRASH); 
}


void hmm_count_stats(exp_transitions,exp_observations,implied_recs,
			sufficient_stats)
real ***exp_transitions;   /* [interval][from_state][to_state] */
real ***exp_observations;  /* [locus][observation/MD][state] */
real ***implied_recs;	   /* [interval][from_state][to_state] */
/* for now, assume exp_observations==NULL, and is not used */
SUFF_STATS **sufficient_stats; /* [locus] => ptr to struct (recs,norecs)... */
{
    int i, k, j;
    real **transitions_exp, **recs_implied=NULL, **norecs_implied=NULL;
    real recs, norecs;

    if (!separate_recs) {
	recs_implied= implied_recs[0];
	norecs_implied= implied_norecs[0];
    }
    for (i=0; i<n_intervals; i++) {
	if (separate_recs) {
	    recs_implied= implied_recs[i];
	    norecs_implied= implied_norecs[i];
	} 
	transitions_exp= exp_transitions[i];
	recs= norecs= 0.0;
	for (k=0; k<max_states; k++)
	  for (j=0; j<max_states; j++) {
	      recs+=   transitions_exp[j][k] * recs_implied[j][k];
	      norecs+= transitions_exp[j][k] * norecs_implied[j][k];
	  }
	sufficient_stats[i]->recs= recs;
	sufficient_stats[i]->norecs= norecs;
    }
}


void hmm_make_new_map(sufficient_stats,new_like,map)
SUFF_STATS **sufficient_stats; /* [locus] => ptr to struct (recs,norecs)... */
real new_like;
MAP *map;
/* for now, assume no penetrance stuff */
{
    int i;
    real recs, norecs, sum, theta=0., ratio;
    
    for (i=0; i<n_intervals; i++) {
	if (map->fix_interval[i]) continue;
	sum=  recs=   sufficient_stats[i]->recs;
	sum+= norecs= sufficient_stats[i]->norecs;
	if (sum>=0.00001) ratio=recs/sum; else ratio=0.0; /* WHAT TO DO? */
	if (ri_type==RI_NOT) theta=ratio;
	  else if (ri_type==RI_SELF) theta=ratio/(2.0*(1.0-ratio));
	  else if (ri_type==RI_SIB)  theta=ratio/(2.0*(2.0-3.0*ratio));
	  else send(CRASH); /* From Haldane and Wadington (inverse of y(x)) */
	
	if (theta>HMM_MAX_THETA) theta=HMM_MAX_THETA;
	  else if (theta<HMM_MIN_THETA) theta=HMM_MIN_THETA;
	map->rec_frac[i][0]=theta;
    }
    map->log_like= new_like;
}



real hmm_e_step(map,trans_prob,obs_prob,exp_transitions,exp_observations)
MAP *map;
real ***trans_prob, ***exp_transitions; /* [interval][from_state][to_state] */
real ***obs_prob, ***exp_observations;  /* [locus][observation][state] */
/* for now, assume exp_observations==NULL, meaning "don't bother" */
/* return likelihood */
{
    int indiv, i, the_observation, lstate, rstate, state, leftk, rightk;
    int locus, prev, num_left_states, num_right_states, num_inf;
    int *the_left_state, *the_right_state;
    real sum, total, val, *pheno_given_geno, *the_prob, *prev_prob;
    real indiv_like, total_like;
    real **the_transitions, **my_transitions, *temp, **the_trans_prob;

    /* Sometime unroll array refs in here for a little speed */

    /* For the moment, we don't ned all of the indiv_transitions matrix - 
       so we KLUDGE and use entry [0] as BOTH my_transitions, and as temp.
       I will fix this someday  - I THINK I JUST DID IT */

    temp= e_step_temp;

    total_like= 0.0;
    for (i=0; i<n_intervals; i++)
      for (lstate=0; lstate<max_states; lstate++) {
	  for (rstate=0; rstate<max_states; rstate++)
	    exp_transitions[i][lstate][rstate]=0.0;
	  if (error_kludge)
	    for (indiv=0; indiv<n_indivs; indiv++)
	      exp_genotype[indiv][i][lstate]=0.0;
      }

    for (indiv=0; indiv<n_indivs; indiv++) {
	for (locus=0, num_inf=0; locus<n_loci; locus++)
	  if (observation[indiv][locus]!=OBS_MISSING) num_inf++;
	if (num_inf<2) continue;

	/**** first, compute the left conditioned probs & loglike ****/
	locus=0; /* leftmost locus */
	total=0.0;
	/* change this to use n_poss_states */
	for (state=0; state<max_states; state++) 
	  total+= left_prob[locus][state]= apriori_prob[indiv][state];
	for (state=0; state<max_states; state++) 
	  left_prob[locus][state]/=total;
	indiv_like= log10(total);

	/* loop over intervals - prev is marker on left and locus is 
	   that on right */
	for (prev=0, locus=1; locus<n_loci; prev=locus, locus++) {
	    num_left_states=  n_poss_states[indiv][prev];
	    num_right_states= n_poss_states[indiv][locus];
	    the_left_state=   poss_state[indiv][prev];
	    the_right_state=  poss_state[indiv][locus];
	    the_prob=  left_prob[locus];
	    prev_prob= left_prob[prev];
	    the_trans_prob= trans_prob[prev];
	    the_observation= observation[indiv][locus];
	    pheno_given_geno= obs_prob[locus][the_observation];
	    for (state=0; state<max_states; state++) the_prob[state]=0.0;
	    total= 0.0; 
	    for (rightk=0; rightk<num_right_states; rightk++) {
		rstate=the_right_state[rightk];
		sum=0.0;
		for (leftk=0; leftk<num_left_states; leftk++) {
		    lstate=the_left_state[leftk];
		    sum+= prev_prob[lstate] * the_trans_prob[lstate][rstate];
		}
		total+= temp[rstate]= sum * pheno_given_geno[rstate];
	    }
	    if (fabs(total)<0.00001) send(CRASH);
	    for (rightk=0; rightk<num_right_states; rightk++) {
		rstate= the_right_state[rightk];
		the_prob[rstate]= temp[rstate]/total;
	    }
	    indiv_like+= log10(total);
	}


	/**** now, the right conditioned probs ****/
	locus= n_loci-1;
	val= 1.0/max_states;
	for (state=0; state<max_states; state++) right_prob[locus][state]=val;

	/* Not the right way: for an F2: 0.25 0.50 0.25, etc
	   for (state=0; state<max_states; state++) 
	   right_prob[locus][state]= apriori_prob[no_data][state]; */
	   
	/* this time, prev is on the right and locus is on the left */
	for (prev=locus, locus=n_loci-2; locus>=0; prev=locus, locus--) {
	    num_left_states=  n_poss_states[indiv][locus];
	    num_right_states= n_poss_states[indiv][prev];
	    the_left_state=  poss_state[indiv][locus];
	    the_right_state= poss_state[indiv][prev];
	    the_prob=  right_prob[locus];
	    prev_prob= right_prob[prev];
	    the_trans_prob= trans_prob[locus];
	    the_observation= observation[indiv][prev]; /* was [locus] */
	    pheno_given_geno= obs_prob[prev][the_observation]; /* ditto */
	    for (state=0; state<max_states; state++) the_prob[state]=0.0;
	    total= 0.0;
	    for (leftk=0; leftk<num_left_states; leftk++) {
		lstate=the_left_state[leftk];
		sum=0.0;
		for (rightk=0; rightk<num_right_states; rightk++) {
#ifdef NOT_THIS /* CHECK: what is correct PGG below? left or right? */
		    sum+= right_prob[prev][right]
		      * trans_prob[locus][left][right]
		      * pheno_given_geno[left]; 
#else 
		    rstate=the_right_state[rightk];
		    sum+= prev_prob[rstate] * the_trans_prob[lstate][rstate]
		        * pheno_given_geno[rstate]; 
#endif
		}
		total+=sum; temp[lstate]=sum;
	    }
	    if (fabs(total)<0.00001) send(CRASH);
	    for (leftk=0; leftk<num_left_states; leftk++) {
		lstate=the_left_state[leftk];
		the_prob[lstate]= temp[lstate]/total;
	    }
	}
          

	/**** and finally, the expected transition matrix for this indiv ****/
	for (i=0; i<n_intervals; i++) {
	    my_transitions=   indiv_transitions[i];
	    num_left_states=  n_poss_states[indiv][i];
	    num_right_states= n_poss_states[indiv][i+1];
	    the_left_state=   poss_state[indiv][i];
	    the_right_state=  poss_state[indiv][i+1];
	    prev_prob= left_prob[i];
	    the_prob=  right_prob[i+1];
	    the_trans_prob= trans_prob[i];
	    the_observation= observation[indiv][i+1];
	    pheno_given_geno= obs_prob[i+1][the_observation];
	    total= 0.0;
	    for (leftk=0; leftk<num_left_states; leftk++) {
		lstate=the_left_state[leftk];
		for (rightk=0; rightk<num_right_states; rightk++) {
		    rstate=the_right_state[rightk];
		    my_transitions[lstate][rstate]= 
		      prev_prob[lstate] * the_trans_prob[lstate][rstate]
			* pheno_given_geno[rstate] * the_prob[rstate];
		    total+= my_transitions[lstate][rstate];
		}
	    }
	
	    /* normalize and massage these data in with that for all indivs
	       also massage into the exp_genotypes global for ERROR KLUDGE 
	       NB: THE EXP-GENOTYPE FOR THE RIGHTMOST INTERVAL IS NOT SET! */
	    the_transitions= exp_transitions[i];
	    for (leftk=0; leftk<num_left_states; leftk++) {
		lstate=the_left_state[leftk];
		for (rightk=0; rightk<num_right_states; rightk++) {
		    rstate=the_right_state[rightk];
		    my_transitions[lstate][rstate]/=total; /* clean this? */
		    the_transitions[lstate][rstate]+=
		      my_transitions[lstate][rstate];
		    if (error_kludge) exp_genotype[indiv][i][lstate]+=
		      my_transitions[lstate][rstate];
		}
	    }
	}
    
#ifdef DEBUGGING
	test_dump0(left_prob,"left_prob");
	test_dump0(right_prob,"right_prob");
	sprintf(ps,"my_trans-indiv %d ",indiv);
	test_dump(indiv_transitions,max_states,ps); 
	r= n= 0.0;
	for (lstate=0; lstate<max_states; lstate++) 
	  for (rstate=0; rstate<max_states; rstate++) {
	      r+= my_transitions[lstate][rstate] 
		* implied_recs[0][lstate][rstate];
	      n+= my_transitions[lstate][rstate] 
		* implied_norecs[0][lstate][rstate];
	  }
	printf("indiv recs=%lf  norecs=%lf\n",r,n); 
#endif
	
	total_like+= indiv_like;
    } /* loop over indivs */

    return(total_like);
}


void test_dump(p,rows,name) 
real ***p;
int rows;
char *name;
{
    int i, j, k, n;

    n= n_intervals > 1 ? 2:1;
    for (k=0; k<n; k++) {

	printf("=== %s[%d] ===\n    ",name,k);
	for (j=0; j<max_states; j++) { printf("    %2d     ",j); }
	printf("\n");

	for (i=0; i<rows; i++) {
	    printf("%2d: ",i); 
	    for (j=0; j<max_states; j++) { printf("%10.8lf ",p[k][i][j]);  }
	    printf("\n");
	}
    }
}


void test_dump0(p,name) 
real **p;
char *name;
{
    int j, k, n;

    n= n_intervals > 1 ? 3:2;
    printf("=== %s ===\n    ",name);
    for (j=0; j<max_states; j++) { printf("    %2d     ",j); }
    printf("\n");

    for (k=0; k<n; k++) {
	printf("%2d: ",k); 
	for (j=0; j<max_states; j++) { printf("%10.8lf ",p[k][j]);  }
	printf("\n");
    }
}



#define state_not_allowed(locus,obs,state) \
  (obs!=OBS_MISSING && obs_prob[locus][obs][state]<0.24)

void count_error_lods(apriori_rate,obs_prob,exp_genotype,
		      error_lod)
real *apriori_rate;       /* [locus#] */
real ***exp_genotype;     /* [indiv][locus][state] */
real ***obs_prob;         /* [locus#][observation][state] */
real **error_lod;         /* [locus#][indiv] side-effected */
{
    real apostiori_rate, e1, e2, like;
    int i, j, state, the_obs;

    for (j=0; j<n_indivs; j++) {
	/* don't bother with first or last locus */
	error_lod[0][j]= error_lod[n_loci-1][j]= 0.0;
	for (i=1; i<n_intervals; i++) {
	    error_lod[i][j]= 0.0;
	    if (apriori_rate[i]<VERY_SMALL) continue;
	    apostiori_rate= 0.0;
	    the_obs= observation[j][i];
	    for (state=0; state<max_states; state++) {
		if (state_not_allowed(i,the_obs,state)) {
		    apostiori_rate+= exp_genotype[j][i][state];
		}
	    }
	    e1= apriori_rate[i]; e2= apostiori_rate;
	    like= (e2/e1) * ((1.0-e1)/(1.0-e2));
	    if (like>VERY_SMALL) error_lod[i][j]= log10(like);
	      else error_lod[i][j]= 0.0;
	}
    }
}


void quick_two_pt(int locus0, int locus1, TWO_PT_DATA *two_pt, bool sex /* do both with and without sex spec */)
{
    if (raw.data_type==F2 &&
	(raw.data.f2.cross_type==F2_INTERCROSS || 
	 raw.data.f2.cross_type==F2_BACKCROSS)) {
      f2_quick_two_pt(locus0,locus1,two_pt,FALSE);

#ifdef HAVE_CEPH
    } else if (raw.data_type==CEPH) {
	ceph_quick_two_pt(locus0,locus1,two_pt,FALSE);
#endif

    } else { /* bail out - use HMM */
	map2->num_loci= 2;
	map2->locus[0]= locus0;
	map2->locus[1]= locus1;
	map2->rec_frac[0][0]= 0.499;
	map2->fix_interval[0]= FALSE;    /* map2->error_rate=NULL */
	setup_hmm(map2);
	twopt_kludge=TRUE;
	hmm_converge_to_map(map2);
	two_pt->theta[NOSEX]= map2->rec_frac[0][0];
	two_pt->lodscore[NOSEX]= map2->log_like - null_like;
	if (two_pt->lodscore[NOSEX]<0.0) two_pt->lodscore[NOSEX]=0.0;
    }
}



void 
hmm_fake_converge_to_map (MAP *map)
{
    int i, j;

    for (i=0; i<map->num_loci-1; i++) map->rec_frac[i][MALE]= randnum()/2.0;
    map->log_like= -100.0 - 20.0*randnum();
    
    if (map->error_lod!=NULL) for (i=1; i<map->num_loci-1; i++) {
	if (map->error_rate[i]>0.0) for (j=0; j<raw.data.f2.num_indivs; j++) {
	    if (randnum()<0.01) map->error_lod[i][j]= 1.0+ 4.1*randnum();
	    else map->error_lod[i][j]= 10.0*randnum();
	}
    }	
}
