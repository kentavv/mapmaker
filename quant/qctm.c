/******************************************************************************

  ####    ####    #####  #    #           ####
 #    #  #    #     #    ##  ##          #    #
 #    #  #          #    # ## #          #
 #  # #  #          #    #    #   ###    #
 #   #   #    #     #    #    #   ###    #    #
  ### #   ####      #    #    #   ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#define INC_LIB
#define INC_MISC
#define INC_CALLQCTM
#define INC_QLOWLEVEL
#include "qtl.h"

/* Globals available to the outside world: */
real like_tolerance, pos_tolerance, mat_tolerance;
bool print_iter, print_rec_mat, bag_qctm, brute_force, print_brute_force;
bool debug_newton;
int  max_intervals, max_genotype_vars, max_continuous_vars;
/* functions alloc_qctm_globals(), free_qctm_globals(), 
   qctm_globals_avail(), qctm_init(), qtl_conv_to_map(); */

/* Local to this file: */
void qctm();
void set_qctm_globals();
void fill_in_qtl_map();

real E_step();
void likelihood();	
void make_rec_probs();

void ML_qtl_weight();
void kill_entry();

void ML_qtl_pos();
real do_brute_force();
real pos_like();

real variance();
real normal_density();

bool debug_qctm;
int n_individuals, n_intervals, n_qtl_genotypes;
int n_genotype_vars, n_continuous_vars;

/* WARNING: Make sure that the (wrong) globals from qtop.c, 
   num_continuous_vars and num_intervals, aren't used in this file! */

/* All of these variables are set & alloced by routines in qdata.c */
int max_interx_genotypes, max_backx_genotypes;
int **lookup_genotype;                   /* [bit-vect][interval#]=> A,B,H */
real **lookup_coded_genotype;            /* [bit-vect][geno-var#]=> real# */

/* State variables for the QCTM E-M loop */
real **expected_genotype;		 /* [indivs][geno_vars+1] */
real *int_like;	                         /* [intervals] */
real **S_matrix, **S_inverse;         	 /* [geno_vars+1][geno_vars+1] */
real **indiv_S_matrix;                	 /* ditto */
real *qctm_qtl_pos;  	                 /* [intervals] */ 
real *null_qtl_weight, *qctm_qtl_weight; /* [geno_vars+1] */
real *fix_qtl_weight, *temp_row;	 /* [geno_vars+1] see ML_qtl_weight */
GENO_PROBS *indiv_rec_like;              /* [interval][int-geno][qtl-geno] */
GENO_PROBS *expected_recs, *trans_prob;	 /* ditto */
INTERVAL_GENOTYPE_PROBS *rec_like;       /* [interval][int-geno] */
bool epistasis_kludge;

#define SMALL_VARIANCE 	 1e-5
#define TINY_VARIANCE 	 0.001
#define TINY_LIKELIHOOD  1e-40
#define SQRT_2_PI 	 2.506628
#define BIGGEST_EXPONENT 1e20
#define NOREC 0
#define REC 1
#define RECS 2




/****************************** PRELIMINARIES ******************************/

void qctm_init() 
{
    max_interx_genotypes= max_backx_genotypes= 0;
    max_intervals= max_genotype_vars= max_continuous_vars= 0;
    pos_tolerance= like_tolerance=  mat_tolerance= -1.0;
    debug_qctm= FALSE;

    lookup_genotype= NULL; lookup_coded_genotype= NULL;
    int_like= NULL;    temp_row= NULL;
    expected_genotype= NULL; 
    S_matrix= S_inverse= indiv_S_matrix= NULL;
    expected_recs= indiv_rec_like= trans_prob= NULL;
    rec_like= NULL;
    qctm_qtl_weight= qctm_qtl_pos= NULL;
    null_qtl_weight= NULL;
}


void set_qctm_globals(data,map)	/* called only by qctm itself */
DATA *data;
QTL_MAP *map;
{
    int i, j, n, k, last;
    int n1, n2, g1, g2;
    real het_additive_coeff, a, b, c;

    /* debug_qctm=TRUE; */

    n_individuals= data->num_individuals;
    n_intervals=   data->num_intervals;
    n_genotype_vars= (raw.data_type==INTERCROSS ? 2:1) * n_intervals;
    if (n_intervals==0) n_qtl_genotypes=0; else 
      n_qtl_genotypes= ipow((raw.data_type==INTERCROSS ? 3:2),n_intervals);
    n_continuous_vars= data->num_continuous_vars;
    epistasis_kludge= FALSE;
    
    /* Paranoia: Check if call to QCTM is OK? */
    /* The n_inter<n_indiv test is to keep variance() from crashing. */
    
    if (like_tolerance<=0.0 || pos_tolerance<=0.0 || mat_tolerance<=0.0 || 
	n_intervals>=n_individuals || !qctm_globals_avail() ||
	n_intervals>max_intervals || n_continuous_vars>max_continuous_vars)
      send(CRASH); 

    if (raw.data_type==BACKCROSS) { 
	/* send(CRASH); Needs to be modified for epistasis and lookups */
	for (i=0; i<n_intervals; i++) {
	    qctm_qtl_weight[i]= map->qtl_weight[i];
	    fix_qtl_weight[i]=  map->constraint[i].backx_weight;
	}
	for (k=0; k<n_continuous_vars; k++) {
	    qctm_qtl_weight[n_intervals+k]= map->cont_var_weight[k];
	    fix_qtl_weight[n_intervals+k]=  map->fix_cont_var_weight[k];
	}
	fix_qtl_weight[n_intervals+n_continuous_vars]= DONT_FIX;

    } else if (raw.data_type==INTERCROSS) {
	for (j=0, n=0; j<n_intervals; j++, n+=2) {
	    qctm_qtl_weight[n]= map->qtl_weight[j];
	    qctm_qtl_weight[n+1]= map->qtl_dominance[j];

	    a=map->constraint[j].a; 
	    b=map->constraint[j].b; 
	    c=map->constraint[j].c; 
	    
	    if (fix_weight_kludge) { /* for scan weight KLUDGE */
	        fix_qtl_weight[n]= map->qtl_weight[j]; fix_qtl_weight[n+1]=0.0;
	    } else if (a==0.0 && b==0.0) { /* nothing fixed */
	        fix_qtl_weight[n]= DONT_FIX; fix_qtl_weight[n+1]= DONT_FIX;
	    } else if (b==0.0) { /* fix additive term */ 
		fix_qtl_weight[n]= c/a;      fix_qtl_weight[n+1]= DONT_FIX;
	    } else { /* b!=0.0 -> fix dominance term */ 
		fix_qtl_weight[n]= DONT_FIX; fix_qtl_weight[n+1]= c/b; 
	    }

	    het_additive_coeff= (b!=0.0 ? (b-a)/b : 1.0);
	    for (i=0; i<n_qtl_genotypes; i++) { 
		if (lookup_genotype[i][j]==A) {
		    lookup_coded_genotype[i][n]=   0.0; 
		    lookup_coded_genotype[i][n+1]= 0.0; 
		} else if (lookup_genotype[i][j]==B) {
		    lookup_coded_genotype[i][n]=   2.0; 
		    lookup_coded_genotype[i][n+1]= 0.0; 
		} else if (lookup_genotype[i][j]==H) {
		    lookup_coded_genotype[i][n]= het_additive_coeff;
		    lookup_coded_genotype[i][n+1]= 1.0; 
		}
	    }
	}

	last= n_genotype_vars; /* actually it's the 1st unused genotype var */

	for (k=0; k<n_continuous_vars; k++) { /* Includes epistasis term */
	    qctm_qtl_weight[last+k]= map->cont_var_weight[k];
	    fix_qtl_weight[last+k]=  map->fix_cont_var_weight[k];
	}

	if (n_continuous_vars>=1 && map->cont_var[0]==EPISTASIS_TERM) {
	    if (n_intervals<2) send(CRASH);
	    n1= n_intervals-2; n2= n_intervals-1; /* the 2 epistatic loci */
	    for (i=0; i<n_qtl_genotypes; i++) { 
		g1=lookup_genotype[i][n1]; g2=lookup_genotype[i][n2];
		if      (g1==A || g2==A)    lookup_coded_genotype[i][last]=0.0;
		else if (g1==H && g2==H)    lookup_coded_genotype[i][last]=1.0;
		else if (g1==B && g2==H)    lookup_coded_genotype[i][last]=2.0;
		else if (g1==H && g2==B)    lookup_coded_genotype[i][last]=2.0;
		else /* (g1==B && g2==B) */ lookup_coded_genotype[i][last]=4.0;
		
	    }
	    epistasis_kludge= TRUE; n_continuous_vars-= 1;
	    last+=1; n_genotype_vars= last;
	}

	/* Mu term - with or without epistasis */
	qctm_qtl_weight[last+n_continuous_vars]=map->mu; 
	fix_qtl_weight[last+n_continuous_vars]= DONT_FIX; 

    } else send(CRASH);
}


void fill_in_qtl_map(map,new_like,qtl_weight)
QTL_MAP *map;
real new_like, *qtl_weight;
{
    int j, k, n, num, last;
    real a, b, c;

    map->abs_log_like= new_like; 
    map->log_like= new_like - map->null_log_like;
    map->var_explained= 1.0 - (map->sigma_sq/map->null_sigma_sq);
    map->chi_sq= 2.0 * map->log_like * log(10.0);
    
    if (raw.data_type==BACKCROSS) {
	for (j=0; j<n_intervals; j++) {
	    map->qtl_weight[j]= qtl_weight[j];
	    map->qtl_dominance[j]= 0.0;
	}
	for (k=0; k<n_continuous_vars; k++) {
	    map->cont_var_weight[k]= qctm_qtl_weight[n_intervals+k];
	}

    } else if (raw.data_type== INTERCROSS) {
	for (j=0, n=0; j<n_intervals; j++, n+=2) {
	    a=map->constraint[j].a; 
	    b=map->constraint[j].b; 
	    c=map->constraint[j].c; 
	    
	    map->qtl_weight[j]= qtl_weight[n];
	    if (b==0.0 || fix_weight_kludge) 
	      map->qtl_dominance[j]= qtl_weight[n+1];
	    else 
	      map->qtl_dominance[j]= (c - a*qtl_weight[n])/b;
	}

	last= (!epistasis_kludge ? n_genotype_vars:   n_genotype_vars-1);
	num=  (!epistasis_kludge ? n_continuous_vars: n_continuous_vars+1);
	for (k=0; k<num; k++) { /* with or without epistasis */
	    map->cont_var_weight[k]= qctm_qtl_weight[last+k];
	}

    } else send(CRASH);
}


/************************** QTL_CONV_TO_MAP (QCTM) **************************/

void qtl_conv_to_map(data,map)
DATA *data;
QTL_MAP *map; 	/* many parts of this are side-effected */
{
    real old_like, new_like;
    int i;
    bool done;
    /* use externs null_qtl_weight, qctm_qtl_weight, qctm_qtl_pos */

    set_qctm_globals(data,map);
    map->chi_sq= map->var_explained= map->log_like= 0.0;

    /*** BAG_QCTM - SPEED THINGS UP FOR DEBUGGING OTHER PARTS OF THE CODE ***/
    if (bag_qctm) {
	map->log_like= map->abs_log_like= map->chi_sq= 1.0;
	map->null_log_like= map->no_data_like= 0.0; map->var_explained= 0.001;
	return;
    }
	
    /*** COMPUTE THE NULL QTL MAPS ***/
    map->null_log_like= 
      E_step(map->qtl_pos,null_qtl_weight,&map->null_mu,&map->null_sigma_sq,
	     data,expected_genotype,S_matrix,expected_recs);
    map->no_data_like= 0.0;
    /* was no_data_like(data,qctm_qtl_weight,mu,sigma_sq) - see Obsolete.c */
    if (print_iter) print_null_iteration(map);
    
    /*** COMPUTE THE ACTUAL QTL MAP ***/
    new_like= VERY_UNLIKELY; i=0; done=FALSE;
    do {
	old_like=new_like; 
	new_like= 
	  E_step(map->qtl_pos,qctm_qtl_weight,&map->mu,&map->sigma_sq,data,
		 expected_genotype,S_matrix,expected_recs); 

	if ((fabs(new_like-old_like) < like_tolerance) || 
	    (map->sigma_sq < SMALL_VARIANCE) ||
	    (fix_weight_kludge && i==200)) done=TRUE; /* converged! */
	else if (new_like < old_like) {
	    print("*** warning: likelihood decreased, quitting...\n");
	    done=TRUE;
	} 

	if (print_iter && (!done || print_rec_mat)) {   
	    fill_in_qtl_map(map,new_like,qctm_qtl_weight);
	    print_iteration(i,map,(i==0 ? 0.0 :	new_like-old_like));
	    if (print_rec_mat) 
	      do_print_E_step(expected_genotype,S_matrix,expected_recs,
		n_individuals,n_genotype_vars,n_continuous_vars,n_intervals);
	}

	if (done) break;

	ML_qtl_weight(S_matrix,expected_genotype,data->phenotype,
	    fix_qtl_weight,&map->mu,qctm_qtl_weight,&map->sigma_sq); 
	ML_qtl_pos(expected_recs,data->interval_len,map->fix_pos,map->qtl_pos);
	i++; /* #iterations counter */

    } while(TRUE); /* break is in the middle of the loop */
    fill_in_qtl_map(map,new_like,qctm_qtl_weight);
}



/********************************** E-STEP **********************************/

real E_step(qtl_pos,qtl_weight,mu,sigma_sq,data,expected_genotype,
	    S_matrix,expected_recs) /* return the log-likelihood */
real *qtl_pos, *qtl_weight;
real *mu, *sigma_sq; /* just ptrs to single reals */
DATA *data;
real **expected_genotype, **S_matrix;  /* both side-effected */
GENO_PROBS *expected_recs;             /* side-effected */
{ 
    int *poss_genotype; 
    real *genotype_contribution;
    real sum_exp_genotype, prediction;
    real like, indiv_like, int_like, sum_int_like, total_log_like;
    int i, j, n, k;   /* use i=individuals; j=intervals; n=genotype_vars */
    int x, y, z, geno, qtl, last, c, d, entry;

    /* externs side-effected: indiv_S_matrix, indiv_rec_like, trans_prob, 
       and rec_like */

if (debug_qctm) { /* DEBUGGING CODE */
    sprintf(ps, "qtl_pos: %lf\n", qtl_pos[0]); pr();
    sprintf(ps, "qtl_weight: %lf\n", qtl_weight[0]); pr();
    sprintf(ps, "mu: %lf\n", *mu); pr();
    sprintf(ps, "sigma_sq: %lf\n", *sigma_sq); pr(); 
} 

    /*** Init expected_genotype, expected_recs, S_matrix and rec_probs ***/
    for (j=0; j<n_intervals; j++) 
      for_interval_genotypes(raw.data_type,geno) 
	for_locus_genotypes(raw.data_type,qtl)
	  expected_recs[j][geno][qtl]= 0.0;
    for (n=0; n<n_genotype_vars; n++) {
	for (k=0; k<n_genotype_vars; k++) S_matrix[n][k]= 0.0;
	for (i=0; i<n_individuals; i++)	  expected_genotype[i][n]= 0.0;
    }
    make_rec_probs(qtl_pos,data->interval_len,trans_prob);
		
    /*** Now compute expected genotypes, recs, S_matrix and likelihood ***/
    total_log_like = 0.0;
    for (i=0; i<n_individuals; i++) {

if (debug_qctm) {
    sprintf(ps, "===== indiv %d \n", i); print(ps);
    for (j=0; j<n_intervals; j++) {
	sprintf(ps, "interval %d probs:\n", j); pr(); z=0;
	for_interval_genotypes(raw.data_type,y) {
	    sprintf(ps, "%2d %lf   ", y, data->genotype_prob[i][j][y]); pr(); z++;
	    if (z%4==0 && y!=0) nl();
	}
	if (z%4!=0) nl();
    }
}

	/*** For each indiv we calculate indiv_like, indiv_rec_like and 
	     indiv_S_matrix. These are then used in the running computation of
	     expected_genotype, S_matrix, expected_recs and the log like.
	     First, we initialize the indiv... variables. ***/

	indiv_like=0.0; sum_int_like=0.0;
	for (j=0; j<n_intervals; j++) 
	  for_interval_genotypes(raw.data_type,geno) 
	    for_locus_genotypes(raw.data_type,qtl)
	      indiv_rec_like[j][geno][qtl]= 0.0;
	for (n=0; n<n_genotype_vars; n++) 
	    for (k=0; k<=n; k++) indiv_S_matrix[n][k]= 0.0;

	/*** Now we loop over all possible qtl-genotypes, getting a like and 
	     the rec_like[][][] for each. These are mashed into indiv_like,
	     indiv_rec_like and indiv_S_matrix. ***/
		
	for (x=0; x<n_qtl_genotypes; x++) {
	    poss_genotype= lookup_genotype[x];
	    genotype_contribution= lookup_coded_genotype[x]; 
	      
	    /* this sets like, int_like, and rec_like */
	    likelihood(poss_genotype,genotype_contribution,trans_prob,
		      qtl_weight,mu,sigma_sq,data,i,&like,&int_like,rec_like); 

	    indiv_like+=like; 
	    sum_int_like+=int_like;
	    for (j=0; j<n_intervals; j++) 
	      for (y=0; y<data->num_genotypes[i][j]; y++) { 
		  geno= data->genotype[i][j][y]; 
		  indiv_rec_like[j][geno][poss_genotype[j]]+=rec_like[j][geno];
	      }
	    for (n=0; n<n_genotype_vars; n++) {
		expected_genotype[i][n]+= genotype_contribution[n] * like;
		for (k=0; k<=n; k++) indiv_S_matrix[n][k]+= 
		  genotype_contribution[n] * genotype_contribution[k] * like;
	    } 
	} /* end of for possible_genotypes (x<n_qtl_genotypes) */

	/* For regression on cont-vars alone, n_qtl_genotypes==0 and thus
	   likelihood() never gets called. Thus, we need to calculate the 
	   indiv_like here... */

	/*** For this indiv, we now normalize the indiv_rec_like and 
	     expected_genotype, then we use this to calculate expected_recs, 
	     S_matrix, and total_log_like. ***/

	if (n_genotype_vars>0) {
	    if (fabs(sum_int_like-1.0)>0.01) 
	      print("*** warning: sum_int_like in E_step is not 1.0\n"); 
	    if (indiv_like > VERY_SMALL) {
		total_log_like+= log10(indiv_like);

		for (j=0; j<n_intervals; j++) 
		  for_interval_genotypes(raw.data_type,geno) 
		    for_locus_genotypes(raw.data_type,qtl)
		      indiv_rec_like[j][geno][qtl]/= indiv_like;
		for (n=0; n<n_genotype_vars; n++) {
		    expected_genotype[i][n]/= indiv_like;
		    for (k=0; k<=n; k++) 
		      S_matrix[n][k]+= indiv_S_matrix[n][k]/indiv_like;
		}
		for (j=0; j<n_intervals; j++) 
		  for_interval_genotypes(raw.data_type,geno) 
		    for_locus_genotypes(raw.data_type,qtl)
		      expected_recs[j][geno][qtl]+= 
			indiv_rec_like[j][geno][qtl];
	    } /* end of if (indiv_like > VERY_SMALL) */

	} else {  /* n_genotype_vars==0 */
	    for (prediction= *mu, c=0; c<n_continuous_vars; c++) 
		prediction+= data->cont_var[c][i] * qtl_weight[c];
	    indiv_like= 
	      normal_density((data->phenotype[i]-prediction),sigma_sq);
	    if (indiv_like > VERY_SMALL) total_log_like+=log10(indiv_like);
	}

if (debug_qctm) { /* DEBUGGING CODE */
    int g, q;
    print("expected_geno: "); for (n=0; n<n_genotype_vars; n++) 
      { sprintf(ps, "%lf  ", expected_genotype[i][n]); print(ps); } nl();
    sprintf(ps, "indiv_like:    %lf   indiv_rec_like:\n", indiv_like); print(ps); 
    for (j=0; j<n_intervals; j++) {
	for_locus_genotypes(raw.data_type,q) {
	    for_interval_genotypes(raw.data_type,g)
	      { sprintf(ps, "[%d][%d]  ", g, q); pr(); } nl();
	    for_interval_genotypes(raw.data_type,g)
	      { sprintf(ps, "%7.5lf ", indiv_rec_like[j][g][q]); pr(); } nl();
	}
    }
}

    } /* end of for i (individuals) */


    /*** Finally, we fill in the upper right triangle of the S_matrix and 
         the constant entries of the S_matrix and the expected_genotypes.
	 These entries correspond to the continuous variables and the base 
	 trait value (mu) in the regression equation, and are used by 
	 ML_qtl_weight(), but not ML_qtl_pos(). ***/
    
    for (n=0; n<n_genotype_vars; n++) /* fill in upper right */
      for (k=0; k<n; k++) S_matrix[k][n]= S_matrix[n][k];   
    last= n_genotype_vars + n_continuous_vars;

    for (i=0; i<n_individuals; i++) expected_genotype[i][last]= 1.0;

    for (c=0; c<n_continuous_vars; c++) {
	entry=n_genotype_vars+c;
	for (i=0; i<n_individuals; i++) 
	  expected_genotype[i][entry]= data->cont_var[c][i];

	for (n=0; n<n_genotype_vars; n++) {
	    for (sum_exp_genotype=0.0, i=0; i<n_individuals; i++) 
	      sum_exp_genotype+=data->cont_var[c][i] * expected_genotype[i][n];
	    S_matrix[n][entry]=S_matrix[entry][n]= sum_exp_genotype;
	}

	for (sum_exp_genotype=0.0, i=0; i<n_individuals; i++) 
	  sum_exp_genotype+=data->cont_var[c][i];
	S_matrix[entry][last]=S_matrix[last][entry]= sum_exp_genotype;

	x=n_genotype_vars+c;
	for (d=0; d<=c; d++) {
	    y=n_genotype_vars+d;
	    for (sum_exp_genotype=0.0, i=0; i<n_individuals; i++) 
	      sum_exp_genotype+=data->cont_var[c][i] * data->cont_var[d][i];
	    S_matrix[x][y]= S_matrix[y][x]= sum_exp_genotype;
	}
    }

    for (n=0; n<last; n++) { /* for all genotype_vars and continuous_vars */
	for (sum_exp_genotype=0.0, i=0; i<n_individuals; i++) 
	  sum_exp_genotype+= expected_genotype[i][n];  /* sum for all indivs */
	S_matrix[last][n]= S_matrix[n][last]= sum_exp_genotype;  /* mu entry */
    }
    S_matrix[last][last]= (real)(n_individuals);

    return(total_log_like);
}


void likelihood(qtl_genotype,contribution,trans_prob,qtl_weight,mu,sigma_sq,
    data,indiv,total_like,total_int_like,rec_like)
int *qtl_genotype;
real *contribution;
GENO_PROBS *trans_prob;
real *qtl_weight;
real *mu, *sigma_sq;  	  /* mu and sigma_sq are pointers to single numbers */
DATA *data;
int indiv;
real *total_like;	  /* side-effected - is a  ptr to a single num */
real *total_int_like;	  /* side-effected - is also a  ptr to a single num */
INTERVAL_GENOTYPE_PROBS *rec_like; /* [interval][int-geno] - side-effected */
{
    real normal_like, prediction, ratio;
    int j, y, n, c, geno;
	
    /* Some elements of the data struct we use in here... */
    INTERVAL_GENOTYPE_PROBS *indiv_prob= data->genotype_prob[indiv];
    int *num_int_genos= data->num_genotypes[indiv];
    INT_POSS_GENOTYPES *poss_geno= data->genotype[indiv];
        /* poss_geno[interval#][poss-geno#]=> genotype code */

    /* int_like is used as a temp in here */

    *total_int_like= 1.0; 
    for (j=0; j<n_intervals; j++) {
	int_like[j]= 0.0;
     /*	for_interval_genotypes(data_type,geno) { EFFECTIVELY */
 	for (y=0; y<num_int_genos[j]; y++) { 
	    geno= poss_geno[j][y]; 
	    int_like[j]+= rec_like[j][geno]= 
	      trans_prob[j][geno][qtl_genotype[j]] * indiv_prob[j][geno];
	}
	*total_int_like *= int_like[j];
    }
    
    for (n=0, prediction= *mu; n<n_genotype_vars; n++) 
      prediction+= contribution[n] * qtl_weight[n];
    for (c=0; c<n_continuous_vars; c++) /* ADDED for cont vars */
      prediction+= data->cont_var[c][indiv] * qtl_weight[n_genotype_vars+c];
    normal_like= normal_density((data->phenotype[indiv]-prediction),sigma_sq);
    *total_like= *total_int_like * normal_like;

    for (j=0; j<n_intervals; j++) {
	ratio= (*total_like>VERY_SMALL ? *total_like/int_like[j] : 0.0);
 	for (y=0; y<num_int_genos[j]; y++) { /* for_interval_genotypes y */
	    geno= poss_geno[j][y]; 
	    rec_like[j][geno] *= ratio;
	}
    }

if (debug_qctm) { /* DEBUGGING CODE */
    int x;
    x=(raw.data_type==INTERCROSS ? 2:1);
    for (j=0; j<n_intervals; j++) {
	sprintf(ps, "poss_geno=%d cont=(%3.1lf %3.1lf) int_like=%lf rec_like:\n",
            qtl_genotype[j], contribution[x*j], contribution[x*j+1], int_like[j]);
	pr(); 
    }
    sprintf(ps, "total_int_like=%lf norml=%lf total=%lf\n=====\n",
            *total_int_like, normal_like, *total_like); pr();
} 

}


void make_rec_probs(qtl_pos,interval_rf,trans_prob)
real *qtl_pos, *interval_rf;
GENO_PROBS *trans_prob;
{
    int j, geno, qtl, left, right;
    real left_rf, right_rf, sum;

    for (j=0; j<n_intervals; j++) {
	left_rf= 
	  rmaxf(MIN_REC_FRAC,min(qtl_pos[j],MAX_FRAC_OF_RF*interval_rf[j]));
	right_rf= (interval_rf[j] - left_rf)/(1 - 2*left_rf);
	for_interval_genotypes(raw.data_type,geno) {
	    sum= 0.0;
	    left= left_genotype[geno];
	    right= right_genotype[geno];
	    for_locus_genotypes(raw.data_type,qtl) {
		sum+= trans_prob[j][geno][qtl]= 
		  (transition_prob(raw.data_type,left,qtl,left_rf) *
		   transition_prob(raw.data_type,qtl,right,right_rf)) /
		 transition_prob(raw.data_type,left,right,interval_rf[j]);
	    }		
	    if (fabs(sum-1.0)>0.01)
	      print("*** warning: the sum of the trans_probs is not 1.0\n");
	}
    }

if (debug_qctm) {
    for (j=0; j<n_intervals; j++) {
	sprintf(ps, "rec_probs for interval %d rf=%lf\n", j, interval_rf[j]); pr();
	for_interval_genotypes(raw.data_type,geno) { 
	    sprintf(ps, "%2d ", geno); pr();
	    for_locus_genotypes(raw.data_type,qtl) 
	      { sprintf(ps, "[%d] %lf ", qtl, trans_prob[j][geno][qtl]); pr(); } nl();
	}
    }
}

}


real normal_density(x,sigma_sq)
real x,*sigma_sq;
{
    real exponent;
    
    if (*sigma_sq < TINY_VARIANCE) {
	print("warning: sigma_sq exceeds tiny_variance...\n");
	*sigma_sq= TINY_VARIANCE; 
    }
    exponent = ((x*x)/(2.0 * *sigma_sq))
             + (0.5 * (log(*sigma_sq))) + log(SQRT_2_PI); 
    
    if (exponent < BIGGEST_EXPONENT) return(exp(-exponent)); 
    else return(TINY_LIKELIHOOD);
}



/******************** ML_QTL_WEIGHT ********************/
    
void ML_qtl_weight(S_matrix,expected_genotype,phenotype,fix_weight,
	mu,qtl_weight,sigma_sq)
real **S_matrix, **expected_genotype, *phenotype, *fix_weight;
real *mu, *qtl_weight, *sigma_sq;  /* side effect these three */
{
    int n, size, i;
    real total, prediction;
    /* S_inverse and row are side-effected */

    size= n_genotype_vars+n_continuous_vars+1;


    if (!fix_weight_kludge) {
	for (n=0; n<size; n++) 
	  if (fix_weight[n]!=DONT_FIX) kill_entry(S_matrix,size,n);

	mat_invert(S_matrix,size,S_inverse);
	array_times_matrix(phenotype,expected_genotype,n_individuals,size,
			   temp_row);
	array_times_matrix(temp_row,S_inverse,size,size,qtl_weight);
	
	for (n=0; n<size; n++) 
	  if (fix_weight[n]!=DONT_FIX) qtl_weight[n]=fix_weight[n];
	*mu= qtl_weight[size-1]; /* last entry is mu */

    } else { /* if fix_weight_kludge */
	for (n=0; n<size-1; n++) qtl_weight[n]= fix_weight[n];
	for (i=0, total=0.0; i<n_individuals; i++) {
	    for (n=0, prediction= 0.0; n<size-1; n++)
	      prediction+= qtl_weight[n] * expected_genotype[i][n];
	    total+= phenotype[i] - prediction;
	}
	*mu= total/((real)n_individuals);
    }

    *sigma_sq= 
      variance(phenotype,qtl_weight,expected_genotype,S_matrix,mu,size);
}


void kill_entry(S_matrix,size,i)
real **S_matrix;
int size, i; /* index of entry to kill */
{
    int n;

    for (n=0; n<size; n++) S_matrix[i][n]= S_matrix[n][i]= 0.0;
    S_matrix[i][i]= 1.0;
}


real variance(phenotype,qtl_weight,expected_genotype,S_matrix,mu,size)
real *phenotype, *qtl_weight, **expected_genotype, **S_matrix, *mu;
int size;
{
    int i, n, m;
    real var, prediction, term1, term2, pheno;
    
    term1= term2= 0.0; 

    for (i=0; i<n_individuals; i++) {
	for (n=0, prediction= 0.0; n<size-1; n++)
	  prediction+= qtl_weight[n] * expected_genotype[i][n];
	pheno = phenotype[i] - *mu;
	term1+= pheno*pheno - 2.0*pheno*prediction;
    }
    
    for (n=0; n<size-1; n++) for (m=0; m<size-1; m++) 
      term2+= qtl_weight[n] * qtl_weight[m] * S_matrix[n][m];
    
    var = (term1 + term2) / ((real)n_individuals);
    return(var);
}


/******************** ML_QTL_POS ********************/

#define DIVS 100
#define STEP 0.02

void ML_qtl_pos(expected_recs,interval_rf,fix_pos,qtl_pos)
GENO_PROBS *expected_recs;
real *interval_rf;
real *fix_pos; /* []==DONT_FIX if we don't want to fix it */
real *qtl_pos; /* side-effected */
{
    int i;
    real start, end, interval_cm, pos, max_cm, min_cm;

    for (i=0; i<n_intervals; i++) {
	
	if (fix_pos[i]!=DONT_FIX) qtl_pos[i]=fix_pos[i];
	
	else if (brute_force) {
	/***** Brute force computation to find a maximum ******/
	    
	    /* First pass */
	    interval_cm= haldane_cm(interval_rf[i]);
	    min_cm= start= MIN_CM;
	    max_cm= end= haldane_cm(MAX_FRAC_OF_RF*interval_rf[i]);
	    pos= do_brute_force(expected_recs,i,interval_rf[i],start,end,DIVS);

	    /* Second pass */
	    start= rmaxf(pos-STEP*interval_cm,min_cm);
	    end= rminf(pos+STEP*interval_cm,max_cm);
	    pos= do_brute_force(expected_recs,i,interval_rf[i],start,end,DIVS);

	    qtl_pos[i]= unhaldane_cm(pos);
	} else /* !brute_force */ send(CRASH);
    } /* for i (all intervals) */
}


real do_brute_force(expected_recs,i,interval_rf,start,end,steps)
  /* return the best pos in cM */ 
GENO_PROBS *expected_recs;
int i; /* interval# */
real interval_rf, start, end; /* start and end are positions in cM */
int steps;
{
    real inc, pos, like, max_like, best_pos=0.;
    int j;

    inc=(end-start)/((real)(steps-1)); max_like= -VERY_BIG;
    for (j=0, pos=start; j<steps; j++, pos+=inc) {
	like= pos_like(unhaldane_cm(pos),interval_rf,expected_recs,i);
	if (like>max_like) { max_like=like; best_pos=pos; }
    }
    return(best_pos);
}


real pos_like(left_theta,interval_theta,expected_recs,i)
GENO_PROBS *expected_recs;
real left_theta, interval_theta;
int i; /* the interval# */
{
    real event_like, pos_like, right_theta;
    int interval_geno, left_geno, right_geno, qtl_geno;

    right_theta= (interval_theta-left_theta)/(1.0-2.0*left_theta);
    pos_like=0.0;

    for_interval_genotypes(raw.data_type,interval_geno) {
	left_geno= left_genotype[interval_geno];
	right_geno= right_genotype[interval_geno];
	for_locus_genotypes(raw.data_type,qtl_geno) {
	    event_like= apriori_prob(raw.data_type,left_geno) *
	      transition_prob(raw.data_type,left_geno,qtl_geno,left_theta) *
	      transition_prob(raw.data_type,qtl_geno,right_geno,right_theta);
	    pos_like+=
	      log(event_like) * expected_recs[i][interval_geno][qtl_geno];
	}
    }
    return(pos_like);
}


/*    real left_recs, right_recs, left_norecs, right_norecs;
    left_recs= left_norecs= right_recs= right_norecs= 0.0;
    if (raw.data_type!=BACKCROSS) send(CRASH);
	    if (left_geno==qtl_geno) 
	       left_norecs+= expected_recs[i][interval_geno][qtl_geno];
	    else left_recs+= expected_recs[i][interval_geno][qtl_geno];
	    if (right_geno==qtl_geno) 
	       right_norecs+= expected_recs[i][interval_geno][qtl_geno];
	    else right_recs+= expected_recs[i][interval_geno][qtl_geno];
    like= left_recs*log(left_theta)   + left_norecs*log(1.0-left_theta) + 
        + right_recs*log(right_theta) + right_norecs*log(1.0-right_theta); */

#ifdef VERY_OLD_CODE

	    start= 0.0;
	    inc= (haldane_cm(interval_rf[i])/(DIVS-1)); 
	    if (print_brute_force) {
		sprintf(ps,PrBF1,i,interval_rf[i],qtl_pos[i]); 
		print(ps);
	    }
	    max_pos=do_brute_force(start,inc,interval_rf[i],
				   expected_recs,i,print_brute_force,&unimodal);
	    
	    if (!unimodal) {
		print("*** warning: pos_likes is not unimodal...\n");
		qtl_pos[i]= max_pos;
	    } else { /* unimodal */
		/* Do it again, with interval= max_pos+/-inc */
		start= haldane_cm(max_pos)-inc;  /* rrange takes care */
		inc= (2*inc)/(DIVS-1);        /* of fencepost error */
		/* if (print_brute_force) { print("second pass:\n"); } */
		qtl_pos[i]=do_brute_force(start,inc,interval_rf[i],
					  expected_recs,i,FALSE,&unimodal);
	    }
	    
	    if (print_iter && !print_brute_force) {
		sprintf(ps,PrI,i,interval_rf[i],
		   guess_pos(expected_recs,i,interval_rf[i]),
		   qtl_pos[i]); print(ps);
	    }
	    
real do_brute_force(start,inc,theta,expected_recs,interval,do_print,unimodal)
real start, inc, theta;
GENO_PROBS *expected_recs;
int interval, do_print, *unimodal;
{
	real pos, cm, max_pos, max_like, d, d2, f, prev_f;
	int j, max_j, dip, got_f;	
	
	max_like= VERY_UNLIKELY; dip= FALSE; *unimodal=TRUE; 
	for (j=0, cm=start; j<DIVS; j++,cm+=inc) {
	    pos=unhaldane_cm(cm);
	    if (!rrange(&pos, MIN_REC_FRAC, MAX_FRAC_OF_RF*theta)) continue;
	    pos_likes(pos,expected_recs,interval,theta,&d,&d2,&f);
	    if (do_print) {
		sprintf(ps,PrBF2,pos,f,dip);print(ps); 
	    }
	    if (f>max_like) { 
		max_like=f; max_pos=pos; max_j=j; 
		if (dip) *unimodal= FALSE;
	    } else if (f<max_like) dip= TRUE;
	}
	if (do_print) {
	    sprintf(ps,PrBF3,max_pos,max_like); 
	    print(ps); 
        }
	if (max_like==VERY_UNLIKELY) { send(CRASH); }
	return(max_pos);
}


void pos_likes(lambda,expected_recs,i,theta,deriv,deriv2,f)
real lambda;
GENO_PROBS *expected_recs;
int i;
real theta;
real *deriv, *deriv2, *f;  
{
	real w1,x1,y1,z1,w2,x2,y2,z2;
	real lambda_star, one_minus_lambda, one_minus_2_lambda;
	real sq_one_minus_2_lambda, theta_minus_lambda;
	real one_minus_theta_lambda;
	real nr1, nr2, r_1, r_2;
	
	rrange(&lambda,0.001,0.499);
	
	r_1= expected_recs[i][REC][REC]   + expected_recs[i][REC][NOREC];
	r_2= expected_recs[i][REC][REC]   + expected_recs[i][NOREC][REC];
	nr1= expected_recs[i][NOREC][REC] + expected_recs[i][NOREC][NOREC];
	nr2= expected_recs[i][REC][NOREC] + expected_recs[i][NOREC][NOREC];
	
	lambda_star= (theta-lambda)/(1.0-2.0*lambda);
	*f= r_1*log(lambda) + r_2*log(lambda_star) +
	    nr2*log(1.0-lambda_star) + nr1*log(1.0-lambda);
/*
	one_minus_lambda=	1.0 - lambda;
	one_minus_2_lambda= 	1.0 - (2.0*lambda);
	sq_one_minus_2_lambda= 	sq(one_minus_2_lambda);
	theta_minus_lambda= 	theta - lambda;
	one_minus_theta_lambda= 1.0 - theta - lambda;

	w1= r_1 / lambda;
	x1= r_2 * (2.0/one_minus_2_lambda - 1.0/theta_minus_lambda);
	y1= nr2 * (-(2.0*lambda_star-1.0) / one_minus_theta_lambda);
	z1= nr1 * (-1.0 / one_minus_lambda);
	
	w2= -r_1 / sq(lambda);
	x2=  r_2 * (4.0/sq_one_minus_2_lambda - 1.0/sq(theta_minus_lambda));
	y2= -nr2 * ((-2.0/sq_one_minus_2_lambda) + 
		    (4.0*theta - 2.0)/
			(sq(one_minus_theta_lambda) * one_minus_2_lambda));
	z2=  nr1 * (-1.0/sq(one_minus_lambda));

	*deriv= w1 + x1 + y1 + z1;
	*deriv2= w2 + x2 + y2 + z2;
*/
}

real guess_pos(expected_recs,i,theta_interval)
GENO_PROBS *expected_recs;
int i;
real theta_interval;
{
	real total_recs, left_recs, right_recs, L, R, X;
	real left_cm, right_cm, cm_pos, pos;
	int left, right;

	total_recs= left_recs= right_recs= 0.0;
	{
	    L= ((real)left); R= ((real)right);
	    X= expected_recs[i][left][right];
	    left_recs+=  L * X;
	    right_recs+= R * X;
	    total_recs+= (L+R) * X;
	}
	
	left_cm=  haldane_cm(rminf(left_recs/total_recs, 0.49));
	right_cm= haldane_cm(rminf(right_recs/total_recs,0.49));
        cm_pos =  haldane_cm(theta_interval) * (left_cm /(left_cm +right_cm));
	pos = unhaldane_cm(cm_pos);
	return(pos);
}   

#endif

void qtl_noconv_to_map(data,map) /* OBSOLETE, I THINK */
DATA *data;
QTL_MAP *map; 	/* many parts of this are side-effected */
{
    real new_like;
    /* use externs null_qtl_weight, qctm_qtl_weight, qctm_qtl_pos */

    set_qctm_globals(data,map);
    map->chi_sq= map->var_explained= map->log_like= 0.0;

    /*** COMPUTE THE NULL QTL MAPS ***/
    map->null_log_like= 
      E_step(map->qtl_pos,null_qtl_weight,&map->null_mu,&map->null_sigma_sq,
	     data,expected_genotype,S_matrix,expected_recs);
    map->no_data_like= 0.0;
    
    /*** COMPUTE THE ACTUAL QTL MAP ***/
    /* qctm_qtl_weight[] was initialized from the map->qtl_weight[] and
       map->qtl_dominance[] entries by set_qctm_globals() */
    new_like= 
      E_step(map->qtl_pos,qctm_qtl_weight,&map->mu,&map->sigma_sq,data,
	     expected_genotype,S_matrix,expected_recs); 
    fill_in_qtl_map(map,new_like,qctm_qtl_weight);
}



#ifdef DEBUGGING_CODE
    printf("Expected Genotype:\n");
    for (i=0; i<n_individuals; i++) {
	printf("%d: ",i);
	for (j=0; j<size; j++) printf("%8.5lf ",expected_genotype[i][j]);
	printf("\n");
    }
    printf("\nS Matrix:\n");
    for (i=0; i<size; i++) {
	printf("%d: ",i);
	for (j=0; j<size; j++) printf("%8.5lf ",S_matrix[i][j]);
	printf("\n");
    }
#endif
      
