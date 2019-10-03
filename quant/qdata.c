/******************************************************************************

  ####   #####     ##     #####    ##             ####
 #    #  #    #   #  #      #     #  #           #    #
 #    #  #    #  #    #     #    #    #          #
 #  # #  #    #  ######     #    ######   ###    #
 #   #   #    #  #    #     #    #    #   ###    #    #
  ### #  #####   #    #     #    #    #   ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#define INC_LIB
#define INC_SHELL
#define INC_CALLQCTM
#define INC_QLOWLEVEL
#include "qtl.h"


/* internal stuff */
void assign_probs();
void make_genotype_arrays();
int map_function;


/********* FUNCTIONS TO DEAL WITH THE DATA AND MAP STRUCTS FOR QCTM *********/
	
void data_init() {}

DATA *alloc_data(num_intervals,num_cont_vars)
int num_intervals, num_cont_vars;
{
    DATA *data;
    
    if (num_intervals<1 || num_cont_vars<0 || !data_loaded()) send(CRASH);
    run {
	single(data, DATA);
	/* Here, the +1s are for paranoia ONLY */
	data->genotype_prob= NULL; data->cont_var=NULL;
	data->phenotype= data->interval_len= NULL;
	data->max_individuals= raw.n_indivs;
	data->max_intervals= num_intervals;
	data->max_continuous_vars= num_cont_vars;

	matrix(data->genotype_prob,raw.n_indivs,num_intervals+1,
	       INTERVAL_GENOTYPE_PROBS);
	matrix(data->cont_var,num_cont_vars,raw.n_indivs,real);
	array(data->phenotype, raw.n_indivs, real);
	array(data->interval_len, num_intervals+1, real);
	matrix(data->num_genotypes, raw.n_indivs, num_intervals, int);
	matrix(data->genotype,raw.n_indivs,
	       MAX_INTERVAL_GENOTYPES,INT_POSS_GENOTYPES);

    } when_aborting { free_data(data); relay; }
    return(data);
}


void free_data(data)
DATA *data;
{
    if (data==NULL) return;
    unmatrix(data->genotype_prob, data->max_individuals,
	     INTERVAL_GENOTYPE_PROBS);
    unmatrix(data->cont_var,data->num_continuous_vars,real);
    unarray(data->interval_len, real);
    unarray(data->phenotype, real);
    unmatrix(data->num_genotypes, raw.n_indivs, int);
    unmatrix(data->genotype, raw.n_indivs, INT_POSS_GENOTYPES);

    /* Just in case a bad pointer to a DATA struct gets dereferenced */
    data->max_intervals= data->max_individuals= data->num_continuous_vars= 0; 
    data->num_intervals= data->num_individuals= data->max_continuous_vars= 0; 
    data->genotype_prob= NULL; data->cont_var=NULL;
    data->phenotype= data->interval_len= NULL;

    unsingle(data, DATA);
}

 
void prepare_data(map,data)
QTL_MAP *map;
DATA *data;		/* side-effected */
{
    int i, j, k, indiv;
    bool missing;

    if (map->num_intervals>data->max_intervals || 
	raw.n_indivs>data->max_individuals ||
	map->num_continuous_vars>data->max_continuous_vars || 
	!valid_trait_num(map->trait)) send(CRASH); 

    data->num_intervals= map->num_intervals;
    data->num_continuous_vars= map->num_continuous_vars;

    for (i=0; i<map->num_intervals; i++) {
	/* if (!check_interval(&map->left[i],&map->right[i],&map->fix_pos[i]))
	   send(CRASH);  A KLUDGE - has this been obsoleted forever? */
	data->interval_len[i]= map_length(map->left[i],map->right[i]);
    }

    for (k=0; k<map->num_continuous_vars; k++) /* CONT_VAR hooks */
      if (map->cont_var[k]!=EPISTASIS_TERM && !valid_trait_num(map->cont_var[k])) 
	send(CRASH);

    for (i=0, indiv=0; i<raw.n_indivs; i++) {
	if (raw.trait[i][map->trait]==MISSING_PHENO) continue;
	for (missing=FALSE, k=0; k<map->num_continuous_vars; k++)
	  if (raw.trait[i][map->cont_var[k]]==MISSING_PHENO) 
	    { missing=TRUE; break; }
	if (missing) continue;

	data->phenotype[indiv]= raw.trait[i][map->trait];
	for (j=0; j<map->num_intervals; j++) 
	  assign_probs(data,i,indiv,j,map->left[j],map->right[j]);
	for (k=0; k<map->num_continuous_vars; k++)
	  if (map->cont_var[k]==EPISTASIS_TERM) data->cont_var[k][indiv]= 0.0;
	  else data->cont_var[k][indiv]= raw.trait[i][map->cont_var[k]];
	indiv++;
    }
    data->num_individuals= indiv;
    if (indiv==0) 
      { error("trait(s) are all missing data for all individuals."); }

   /* for (i=0; i<raw.n_indivs; i++) {   DEBUGGING CODE
	sf(ps,"indiv %d\n",i); pr();
	for (j=0; j<9; j++) {
	    sf(ps,"%6.4lf ",data->genotype_prob[i][0][j]); pr();
	}
	nl();
    } */
}


void assign_probs(data,raw_i,data_i,interval,left,right)
DATA *data;        /* Needs data->interval_len[interval] */
int raw_i, data_i; /* The individual's numbers in the data and raw structs */
int interval;      /* The interval number */
int left, right;   /* The left and right locus numbers */
{
    int geno, num;

    if (right==INF_LOCUS) send(CRASH); /* KLUDGE for now */

    indiv_interval_probs(data->genotype_prob[data_i],interval,raw_i,
			 left,right,data->interval_len[interval]);
    num=0;
    for_interval_genotypes(raw.data_type,geno) 
      if (data->genotype_prob[data_i][interval][geno] > VERY_SMALL) 
	data->genotype[data_i][interval][num++]=geno;
    data->num_genotypes[data_i][interval]=num;
}

	
void initial_qctm_values(data, map) 
DATA *data;
QTL_MAP *map; 	/* side-effected */
/* map->fix_pos, fix_weight, fix_dominance, trait, left, and right must be 
   set, and prepare_data() must have been run on data */
{
    int i, k;
    real a, b, c, mu, sigma_sq;

    if (map==NULL || data==NULL) send(CRASH);

    for (i=0; i<map->num_intervals; i++) {
	map->interval_len[i]= data->interval_len[i]; 
	
	/* set initial qtl_pos, weight, and dominance */
	if (map->fix_pos[i]!=DONT_FIX)
	  map->qtl_pos[i]= map->fix_pos[i];
	else if (data->interval_len[i]>=MAX_REC_FRAC) /* is unlinked */
	  map->qtl_pos[i]= 0.2;
	else
	  map->qtl_pos[i]=unhaldane(haldane(data->interval_len[i])/2.0);
	if (map->qtl_pos[i] > data->interval_len[i]) 
	  map->qtl_pos[i] = data->interval_len[i];
	
	if (!fix_weight_kludge)	map->qtl_weight[i]= 0.0; 
	map->qtl_dominance[i]= 0.0;
	
	if (raw.data_type==BACKCROSS) {
	    if (map->constraint[i].backx_weight!=DONT_FIX) 
	      map->qtl_weight[i]= map->constraint[i].backx_weight;
	} else if (raw.data_type==INTERCROSS) {
	  switch (map->constraint[i].interx_type) { 
	      case FREE:	a=0.0; b=  0.0; c=0.0; break;
	      case DOMINANT:    a=1.0; b= -1.0; c=0.0; break;
	      case RECESSIVE:   a=1.0; b=  1.0; c=0.0; break;
	      case F3DOMINANT:  a=1.0; b= -2.0; c=0.0; break;
	      case F3RECESSIVE: a=1.0; b=  2.0; c=0.0; break;
	      case ADDITIVE:    a=0.0; b=  1.0; c=0.0; break;
	      case TEST_MODELS: send(CRASH); /* should never make it to here */
	      case FIXED:	/* a and b should be set, and c is ignored */
	      			map->qtl_weight[i]= 0.0; 
	      			map->qtl_dominance[i]= 0.0; break;
	      case CONSTRAINED: break; /* a,b, and c should be set! */
	      default: send(CRASH);
	  }
	  if (map->constraint[i].interx_type!=FIXED && 
	      map->constraint[i].interx_type!=CONSTRAINED) { 
	    	map->constraint[i].a= a;
	    	map->constraint[i].b= b;
	    	map->constraint[i].c= c;
	    }
      } else send(CRASH);
    }
    
    for (k=0; k<map->num_continuous_vars; k++) {
	if (map->fix_cont_var_weight[k]==DONT_FIX) map->cont_var_weight[k]=0.0;
	else map->cont_var_weight[k]= map->fix_cont_var_weight[k];
    }	
    
    for (i=0, mu= 0.0; i<data->num_individuals; i++) 
      mu+= data->phenotype[i];
    mu/= (real)(data->num_individuals);
    
    for (i=0, sigma_sq=0.0; i<data->num_individuals; i++) 
      sigma_sq+= sq(mu - data->phenotype[i]);
    sigma_sq/= (real)(data->num_individuals);
    
    map->mu= map->null_mu= mu; 
    map->sigma_sq= map->null_sigma_sq= sigma_sq;
    map->abs_log_like= map->null_log_like= map->log_like= 0.0;
    map->var_explained= map->chi_sq= 0.0;
}




/*************************** Setup stuff for QCTM ***************************/

void alloc_qctm_globals()
{
    int j;
    
    /* defaults set in qcmds.c ??? */
    if (max_intervals<0 || max_continuous_vars<0) send(CRASH);
    if (max_intervals>MAX_INTERVALS || 
	max_continuous_vars>MAX_CONTINUOUS_VARS) send(CRASH);

    /* Vars: 1/2 for each int + 1 for each CV + 1 for epistasis + mu */
    if (raw.data_type==BACKCROSS) 
      max_genotype_vars= max_intervals + max_continuous_vars + 2; 
    else if (raw.data_type== INTERCROSS) 
      max_genotype_vars= 2*max_intervals + max_continuous_vars + 2; 
    else send(CRASH);
    
    run { /* NOTE: CHECK SIZES OF ARRAYS many need max_cont_vars */
	make_genotype_arrays(raw.data_type,max_intervals);

matrix(expected_genotype,raw.n_indivs,		max_genotype_vars+1,	real);
matrix(S_matrix,	 max_genotype_vars,	max_genotype_vars+1,	real);
matrix(S_inverse,	 max_genotype_vars,	2*(max_genotype_vars+1),real);
matrix(indiv_S_matrix,	 max_genotype_vars,	max_genotype_vars,	real);
array(qctm_qtl_weight,	 max_genotype_vars,	real); 
array(fix_qtl_weight,	 max_genotype_vars, 	real);
array(null_qtl_weight,	 max_genotype_vars,	real);
array(temp_row,	 	 max_genotype_vars, 	real);
array(qctm_qtl_pos,	 max_intervals, 	real);
array(int_like,	 	 max_intervals,		real);
array(expected_recs,	 max_intervals,	 	GENO_PROBS);
array(indiv_rec_like,	 max_intervals, 	GENO_PROBS);
array(trans_prob,	 max_intervals,		GENO_PROBS);
array(rec_like,		 max_intervals, 	INTERVAL_GENOTYPE_PROBS);

for (j=0; j<max_genotype_vars+1; j++) null_qtl_weight[j]= 0.0;

    } when_aborting {
	free_qctm_globals();
	relay;
    }
}


#define RIGHT_BIT_SET(geno_vector) ((int)(geno_vector & (GENOTYPE)1))

void make_genotype_arrays(data_type,num_intervals)
int data_type,num_intervals;
{
    int i, j, this_genotype, N;
    GENOTYPE geno;

    /* The array and matrix macros don't handle things bigger with a dimension
       bigger than 32K. For now we assume that GENOTYPE is a short, and that 
       num_intervals <= 7. If this must be changed later, it can be. */

    N=num_intervals; 
 
    if (data_type==BACKCROSS) {
	max_backx_genotypes=  (GENOTYPE)lpow2(N);
	matrix(lookup_genotype,       (int)max_backx_genotypes, N,  int);
	matrix(lookup_coded_genotype, (int)max_backx_genotypes, N,  real);

	for (i=0; i<max_backx_genotypes; i++) /* all bit vectors */
	  for (geno=(GENOTYPE)i, j=0; j<num_intervals; j++, geno>>=1) { 
	      lookup_genotype[i][j]=(RIGHT_BIT_SET(geno) ? H:A);
	      lookup_coded_genotype[i][j]= (RIGHT_BIT_SET(geno) ? 1.0:0.0);
      }

    } else if (data_type==INTERCROSS) {
	/* max_interx_genotypes= 3^N - NEW [was lpow2(2*N)] */
	for (i=0, max_interx_genotypes=1; i<N; i++) max_interx_genotypes*=3;
	/* 2*N+1 is max genotype vars */
	matrix(lookup_genotype,       (int)max_interx_genotypes,N,int);
	matrix(lookup_coded_genotype, (int)max_interx_genotypes,2*N+1,real);

	for (i=0; i<max_interx_genotypes; i++) { /* all genotype vectors */
	  /* Loop over intervals */
	  for (geno=i, j=0; j<num_intervals; j++, geno/=3) {
	      /* j is interval# and k is a genotype_var# */
	      this_genotype= ((int)geno)%3;
	      if (this_genotype==0)      lookup_genotype[i][j]= A;
	      else if (this_genotype==1) lookup_genotype[i][j]= H;
	      else if (this_genotype==2) lookup_genotype[i][j]= B;
	      else send(CRASH);
	  }
      }	

    } /* end for INTERCROSS */
}


void free_qctm_globals() {}  /* KLUDGE: CURRENTLY A NOP */
bool qctm_globals_avail()  { return(null_qtl_weight!=((real*)NULL)); }



/********************* RANDOM FUNCTIONS *********************/

real haldane(theta)
real theta;
{
	if (theta==0.0) return(0.0);
	else if (theta>=MAX_REC_FRAC) return(MAX_CM);
	else return(-0.50 * log(1-2*theta)); 
}


real unhaldane(morgans)
real morgans;
{
    if (morgans==0.0) return(0.0);
    else if (morgans>=(MAX_CM/100.0)) return(MAX_REC_FRAC);
    else return((1.0-exp(-2.0*morgans))/2.0);    
}

real kosambi(theta)
real theta;
{
	if (theta==0.0) return(0.0);
	else if (theta>=MAX_REC_FRAC) return(MAX_CM);
	else return(0.25 * log((1.0 + 2.0 * theta)/(1.0 - 2.0 * theta)));
}


real unkosambi(morgans)
real morgans;
{
  	real exp_4_morgans;

	if (morgans==0.0) return(0.0);
	else if (morgans>=(MAX_CM/100.0)) return(MAX_REC_FRAC);
	else {
	  exp_4_morgans= exp(4.0 * morgans);
	  return(0.5 * (exp_4_morgans - 1.0) / (exp_4_morgans + 1.0));
	}
}


#ifdef UNUSED_CODE
/* name conflicts with map_func in def.h, included by f2_prep and map -
   will need to be changed before we can use this! */
real map_func(theta)
real theta;
{ switch(map_function) {
    case HALDANE: return(haldane(theta));
    case KOSAMBI: return(kosambi(theta));  
    default: send(CRASH); 
  }
}


real unmap_func(morgans)
real morgans;
{ switch(map_function) {
    case HALDANE: return(unhaldane(morgans));
    case KOSAMBI: return(unkosambi(morgans));  
    default: send(CRASH); 
  }
}
#endif


real model_prediction(map,indiv)
QTL_MAP *map;
int indiv;
{
    int j, k, geno, qtl, left, right;
    real pheno, interval_rf, left_rf, right_rf, sum, num, dom;
    INTERVAL_GENOTYPE_PROBS *interval_prob; /* p[0][genotype-code] => real*/
    LOCUS_GENOTYPE_PROBS qtl_geno_prob, trans_prob; /* [qtl-geno] => real */

    if (raw.data_type!=INTERCROSS) send(CRASH);
    if (raw.trait[indiv][map->trait]==MISSING_PHENO) 
      return(MISSING_PHENO);
    for (k=0; k<map->num_continuous_vars; k++)
      if (raw.trait[indiv][map->cont_var[k]]==MISSING_PHENO) 
	return(MISSING_PHENO);

    array(interval_prob,1,INTERVAL_GENOTYPE_PROBS);
    pheno=0.0;
    for (k=0; k<map->num_continuous_vars; k++) {
	if (raw.trait[indiv][map->cont_var[k]]==MISSING_PHENO)
	  return(MISSING_PHENO);
	pheno+= map->cont_var_weight[k] * raw.trait[indiv][map->cont_var[k]];
    }

    /* this is a bit of a KLUDGE for now, in term of how we find the 
       a priori genotype probs. */
    for (j=0; j<map->num_intervals; j++) { /* fill in prob[] */
	interval_rf= map->interval_len[j];
	left_rf= rmaxf(MIN_REC_FRAC,min(map->qtl_pos[j],
					MAX_FRAC_OF_RF*interval_rf));
	right_rf= (interval_rf - left_rf)/(1 - 2*left_rf);
	indiv_interval_probs(interval_prob,0,indiv,map->left[j],map->right[j],
			     interval_rf); /* see qraw.c */

	for_locus_genotypes(raw.data_type,qtl) 
	    qtl_geno_prob[qtl]=0.0;
	for_interval_genotypes(raw.data_type,geno) {
	    left= left_genotype[geno];
	    right= right_genotype[geno];
	    sum= 0.0;
	    for_locus_genotypes(raw.data_type,qtl) {
		sum+= trans_prob[qtl]=
		  (transition_prob(raw.data_type,left,qtl,left_rf) *
		   transition_prob(raw.data_type,qtl,right,right_rf)) /
		     transition_prob(raw.data_type,left,right,interval_rf);
		qtl_geno_prob[qtl]+= trans_prob[qtl]*interval_prob[0][geno];
	    }
	    if (fabs(sum-1.0)>0.000001)
	      print("*** warning: the sum of the geno_probs is not 1.0\n");
	}

	sum=0.0;
	for_locus_genotypes(raw.data_type,qtl) {
	    switch(qtl) {
		case A: num=0.0; dom=0.0; break;
		case H: num=1.0; dom=1.0; break;
		case B: num=2.0; dom=0.0; break;
	    }
	    pheno+= (num * map->qtl_weight[j] + dom * map->qtl_dominance[j])
	            * qtl_geno_prob[qtl];
	    sum+= qtl_geno_prob[qtl];
	}
	if (fabs(sum-1.0)>0.000001)
	  print("*** warning: the sum of the qtl_probs is not 1.0\n");
    } /* do next interval */
    return(pheno);
}


