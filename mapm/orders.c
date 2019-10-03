/******************************************************************************

  ####   #####   #####   ######  #####    ####            ####
 #    #  #    #  #    #  #       #    #  #               #    #
 #    #  #    #  #    #  #####   #    #   ####           #
 #    #  #####   #    #  #       #####        #   ###    #
 #    #  #   #   #    #  #       #   #   #    #   ###    #    #
  ####   #    #  #####   ######  #    #   ####    ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#define INC_LIB
#define INC_SHELL
#define INC_MISC
#include "mapm.h"

#define set_excluded(a,b,c) \
  three_pt_excluded[three_pt_index[a]][three_pt_index[b]][three_pt_index[c]] \
    =EXCLUDED /* OBSOLETE? */

void left_exclusion(), right_exclusion(), middle_exclusion();
bool is_excluded(), is_computed();

void add_to_order();
int  compare_markers_to_try();
void randomize_markers_to_try();
void rank_markers_to_try();
void add_order();
bool extend_order_by_1();
int  **add_orders();
bool middles_excluded();

char ***three_pt_excluded=NULL;
int  three_pt_size=0;
int *three_pt_index=NULL;

#define EXCLUDED   0
#define ALLOWED    1
#define UNCOMPUTED 2

void alloc_3pt_matrix(); /* args: int num_loci; sets a global, will free 
   and realloc if it's too small, use to prealloc for efficiency */
void free_3pt_matrix(); 



/**************** Two-point, Haplotypes ****************/

void get_linkage_group(locus,num_loci,linkage_group,group_size,lodbound,
		       thetabound)
int *locus, *num_loci, *linkage_group, *group_size;
real lodbound,thetabound;
{
    int i, j, linked, *unlinked_list, unlinked_size, oldsize, a, b;
    real lod, theta, thetam, thetaf;

    /* start linkage_group with first marker in the list of loci */
    linkage_group[0]=locus[0];
    *group_size= 1;
    unlinked_list= locus;
    *num_loci-= 1;  
    unlinked_size= *num_loci;
    locus= &unlinked_list[1];
    
    do {
	oldsize= unlinked_size;
	unlinked_size= 0;
	for (i=0; i<oldsize; i++) { /* for all loci a NOT in the group */
	    a= locus[i];
	    for (j=0, linked=FALSE; j<*group_size; j++) { /* for all in group*/
		b= linkage_group[j];
		if (!sex_specific) {
		    get_two_pt(a,b,&lod,&theta);
		    if (lod>=lodbound && theta<=thetabound)
		      { linked=TRUE; break; }
		} else {
		    get_sex_two_pt(a,b,&lod,&thetam,&thetaf);
		    if (lod>lodbound && 
			(thetam<=thetabound || thetaf<=thetabound))
		      { linked=TRUE; break; }
		}
	    }
	    if (linked) {
		linkage_group[*group_size]= a;
		++*group_size;
	    } else {
		/* Effectively, write back into the locus array, relying on the
		   fact that unlinked_size<=*num_loci */
		unlinked_list[unlinked_size]= a;
		++unlinked_size;
	    }
	} /* for i */
	locus= unlinked_list;
	*num_loci= unlinked_size;
    } while (unlinked_size!=oldsize && unlinked_size!=0);
}


void find_haplo_group(locus,num_loci,haplo_group,num_haplo,old_obs,new_obs)
int *locus, *num_loci, *haplo_group, *num_haplo;
int *old_obs, *new_obs;
{
    int i, num_inf, num_dom, best_i, foo;
    int best_score, old_score, score;
    
    /* get the most informative marker */
    best_score= 0;
    for (i=0; i<*num_loci; i++) {
	f2_genotype(locus[i],FALSE,old_obs);
	num_inf= f2_count_infs(&num_dom,&foo,old_obs);
	score= 2*num_inf - num_dom;
	if (score>best_score) { best_score=score; best_i=i; }
    }
    
    /* delete it, and make it the haplo group */
    haplo_group[0]=locus[best_i]; *num_haplo=1;
    f2_genotype(locus[best_i],FALSE,old_obs);
    locus[best_i]= locus[--*num_loci]; 
    if (*num_loci==0) return;
    
    do { /* find one to add to the group */
	old_score= best_score;
	best_score= 0;
	for (i=0; i<*num_loci; i++)
	  if (merge_genotypes(locus[i],old_obs,new_obs)) { /* if no recs */
	      num_inf= f2_count_infs(&num_dom,&foo,new_obs);
	      score= 2*num_inf - num_dom;
	      if (score>best_score) { best_score=score; best_i=i; }
	  }
	if (best_score !=0 ) { /* if we found any with no recs */
	    merge_genotypes(locus[best_i],old_obs,old_obs); /*updates old_obs*/
	    haplo_group[*num_haplo]=locus[best_i]; ++*num_haplo; /* add it */
	    locus[best_i]= locus[--*num_loci]; /* delete it */
	}
    } while (*num_loci>0 && best_score!=0); /* got one and more to try */
}



/**************** Three Point ****************/


void alloc_3pt_matrix(num_group)
int num_group;
{ 
    if (three_pt_excluded!=NULL && num_group<=three_pt_size) return;

    run {
	if (three_pt_excluded!=NULL) free_3pt_matrix();
	three_pt_size=num_group;
	three_pt_excluded= alloc_char_3d_matrix(num_group,num_group,num_group);
	array(three_pt_index,raw.num_markers,int);
    } when_aborting {
	free_3pt_matrix();
	relay_messages;
    }
}
      


void free_3pt_matrix()
{
    free_char_3d_matrix(three_pt_excluded,three_pt_size,three_pt_size);
    unarray(three_pt_index,int);
    three_pt_excluded=NULL; three_pt_size=0;
}

    

void setup_3pt_data(locus,num_loci,threshold)
int *locus, num_loci;
real threshold;
{
    int i, j, k;
    real d1, d2, d3;

    alloc_3pt_matrix(num_loci);

    /* three_pt_index[] is lookup by raw locus# to three_pt_exclusion index */
    for (i=0; i<raw.num_markers; i++) three_pt_index[i]= -1;
    for (i=0; i<num_loci; i++) three_pt_index[locus[i]]=i;
    
    for (i=0; i<num_loci; i++) 
      for (j=0; j<num_loci; j++) 
	for (k=0; k<num_loci; k++) 
	  three_pt_excluded[i][j][k]= UNCOMPUTED;
    
    for (i=0; i<num_loci-2; i++) {
	for (j=i+1; j<num_loci-1; j++) {
	    for (k=j+1; k<num_loci; k++) {

		if (!restore_triple(locus[i],locus[j],locus[k],&d1,&d2,&d3)) 
		  continue;

		/* maybe later compute on the fly ADD keep_user_amused 
		   compute_triple(locus[i],locus[j],locus[k],&d1,&d2,&d3); */

		if (d1<threshold) {
		    three_pt_excluded[i][j][k]=EXCLUDED;
		    three_pt_excluded[k][j][i]=EXCLUDED;
		} else {
		    three_pt_excluded[i][j][k]=ALLOWED;
		    three_pt_excluded[k][j][i]=ALLOWED;
		}
		
		if (d2<threshold) {
		    three_pt_excluded[i][k][j]=EXCLUDED;
		    three_pt_excluded[j][k][i]=EXCLUDED;
		} else {
		    three_pt_excluded[i][k][j]=ALLOWED;
		    three_pt_excluded[j][k][i]=ALLOWED;
		}
		
		if (d3<threshold) {
		    three_pt_excluded[j][i][k]=EXCLUDED;
		    three_pt_excluded[k][i][j]=EXCLUDED;
		} else {
		    three_pt_excluded[j][i][k]=ALLOWED;
		    three_pt_excluded[k][i][j]=ALLOWED;
		}
	    }
	}
    }
}

void free_3pt_data() { free_3pt_matrix(); }


bool start_3pt(locus,num_loci,threshold,order,num_ordered) /* improve this! */
int *locus, *num_loci;     /* starters are deleted */
real threshold;
int *order, *num_ordered; /* the starting order */
{
    int i, j, k, num_left, starter[3], temp;
    real d1, d2, d3, best;

    starter[0]= starter[1]= starter[2]= NO_LOCUS;
    best= 0.0;

    for (i=0; i<*num_loci-2; i++) {
	for (j=i+1; j<*num_loci-1; j++) {
	    for (k=j+1; k<*num_loci; k++) {
		if (!restore_triple(locus[i],locus[j],locus[k],&d1,&d2,&d3))
		  continue;
		
		if (d2<best && d3<best) {
		    best=min(d2,d3); if (best>threshold) { best=0.0; continue;}
		    starter[0]=i; starter[1]=j; starter[2]=k; 
		} else if (d1<best && d3<best) {
		    best=min(d1,d3); if (best>threshold) { best=0.0; continue;}
		    starter[0]=i; starter[1]=k; starter[2]=j;
		} else if (d1<best && d2<best) {
		    best=min(d1,d2); if (best>threshold) { best=0.0; continue;}
		    starter[0]=j; starter[1]=i; starter[2]=k;
		}
	    }
	}
    }

    if (best==0.0) return(FALSE);

    /* else, setup order list and delete starters from locus list */
    order[0]=locus[starter[0]]; locus[starter[0]]=NO_LOCUS;
    order[1]=locus[starter[1]]; locus[starter[1]]=NO_LOCUS;
    order[2]=locus[starter[2]]; locus[starter[2]]=NO_LOCUS;
    *num_ordered= 3;

    num_left=0;
    for (i=0; i<*num_loci; i++) {
	temp=locus[i];
	if (temp!=NO_LOCUS) locus[num_left++]=temp;
    }
    *num_loci= num_left;
    return(TRUE);
}



bool is_excluded(a,b,c)
int a, b, c;
{ return(three_pt_excluded[three_pt_index[a]][three_pt_index[b]]
	                  [three_pt_index[c]]==EXCLUDED); }
 

bool is_computed(a,b,c)
int a, b, c;
{ return(three_pt_excluded[three_pt_index[a]][three_pt_index[b]]
                          [three_pt_index[c]]!=UNCOMPUTED); }
	 

void left_exclusion(excluded,locus,num_loci,newmarker)
bool *excluded;
int *locus, num_loci, newmarker;
{
    int middle, end, i;

    for (middle=num_loci-2; middle>=0; middle--) {
	for (end=num_loci-1; end>middle; end--) {
	    if (is_excluded(newmarker,locus[middle],locus[end])) {
		for (i=0; i<=middle; i++)
		  excluded[i]=TRUE;
		return;
	    }
	}
    }
}


void right_exclusion(excluded,locus,num_loci,newmarker)
bool *excluded;
int *locus, num_loci, newmarker;
{
    int middle,end,i;

    for (middle=1; middle<=num_loci-1; middle++) {
	for (end=0; end<middle; end++) {
	    if (is_excluded(locus[end],locus[middle],newmarker)) {
		for (i=num_loci; i>middle; i--)
		  excluded[i]=TRUE;
		return;
	    }
	}
    }
}


void middle_exclusion(excluded,locus,num_loci,leftmost,rightmost,
		      init_window_size,newmarker)
bool *excluded;
int *locus, num_loci;
int leftmost, rightmost;  /* presumably, the left and rightmost intervals? */
int init_window_size, newmarker;
{
    int first, last, start, end, left, width, i, j;
    
    for (width=init_window_size; width>=1; width--) {
	first=max(1,leftmost-width+1);
	last= min(rightmost,num_loci-width);
	for (i=first; i<=last; i++) {
	    start=i;
	    end=i+width-1;
	    if (is_excluded(locus[start-1],newmarker,locus[end])) {
		for (j=start; j<=end; j++)
		  excluded[j]=TRUE;
		if (start>leftmost && width>1)
		  middle_exclusion(excluded,locus,num_loci,leftmost,start-1,
				   width-1,newmarker);
		if (end==rightmost) return;
		leftmost=end+1;
	    }
	}
    }
}


/* orig frame:         a - b - c - d - e - f - g
   orig frame index:   0 - 1 - 2 - 3 - 4 - 5 - 6    num_placed=7
   exclusion index:  0 - 1 - 2 - 3 - 4 - 5 - 6 - 7                */

int  three_pt_exclusions(order,num_placed,newmarker,excluded)
int  *order, num_placed, newmarker;
bool *excluded;  /* side-effected */
{
    int i, j, left, right, num_places;
    for (i=0; i<num_placed+1; i++) excluded[i]=FALSE;

    left_exclusion(excluded,order,num_placed,newmarker);
    right_exclusion(excluded,order,num_placed,newmarker);
	
    left=1; right=num_placed;
    j=1; while (j<num_placed && excluded[j]) j++;
    if (j!=num_placed)
      middle_exclusion(excluded,order,num_placed,left,right,right-left,
		       newmarker);
    for (num_places=0, i=0; i<num_placed+1; i++) 
      if (!excluded[i]) ++num_places;
    return(num_places);
}


bool three_pt_verify(locus,num_loci,window_size)
int *locus, num_loci, window_size;
{
    int i, j, k, left, right;
    if (num_loci<3 || window_size<3) return(TRUE); /* crash? */

    for (i=0; i<num_loci-2; i++)
      for (k=i+2; k<num_loci && k<i+window_size; k++)
	for (j=i+1; j<k; j++)
	  if (is_excluded(locus[i],locus[j],locus[k])) return(FALSE);

    return(TRUE);
}


/**************** N-Pt Maker Placer ****************/

void place_locus(map,locus,start,finish,excluded,result,best_pos,best_map,temp)
MAP *map;          /* map to place locus into, should be ctm-ed! */
int locus;
int start, finish; /* The leftmost and rightmost frame markers (i) to use */
bool *excluded; /* [interval#], one must be FALSE, NOT set */
PLACE **result; /* [interval#]->foo side-effected like=NO_LIKE if untried */
int *best_pos;   /* side effected */
MAP *best_map, *temp;  /* max_loci must be >=finish-start+2 */
{
    int i, j, k, num_allowed, last, *zero=NULL;
    real best_like, lod, theta;
    real theta1, theta2, dist1, dist2, distframe, distscale;
    
    /* Generate N-point maps from left to right markers with locus 
       inserted into each of the allowed positions, even one */

    for (j=0, num_allowed=0; j<=map->num_loci; j++) {
	if (!excluded[j]) num_allowed++;
	result[j]->like=NO_LIKE;   result[j]->dist=NO_THETA; 
	result[j]->net_error= 0.0; result[j]->worst_error= 0.0;
	result[j]->zero=FALSE;
    }
	
    if (sex_specific || num_allowed==0 || best_map->max_loci<finish-start+2 ||
	temp->max_loci<finish-start+2) send(CRASH);
      
    last=finish+1;
    best_like=VERY_UNLIKELY;
    for (i=start; i<=last; i++) /* i=the interval # to place in */
      if (!excluded[i]) {
	  clean_map(temp); k=0;
	  if (i==0) temp->locus[k++]=locus;
	  for (j=start; j<=finish; j++) { /* locus #s */
	      temp->locus[k++]=map->locus[j];
	      if (j+1==i) temp->locus[k++]=locus;
	  }
	  temp->num_loci= k;
	  init_rec_fracs(temp);
	  converge_to_map(temp);
	  result[i]->like=temp->log_like;

	  if (temp->log_like>best_like) { /* keep the best */
	      best_like=temp->log_like; *best_pos=i; 
	      mapcpy(best_map,temp,FALSE); 
	  }
	  
	  /* normalize placement dists by any expansion of this interval
	     if i=first or last interval in the FRAMEWORK, then 
	     dist= the correct rec frac, otherwise the placement dist = 
	     DL*(Dframe/(DL+DR)) */
	  if (i==start) {
	      result[i]->dist= theta= temp->rec_frac[0][0];
	      if (theta<ZERO_DIST) result[i]->zero=TRUE;
	  } else if (i==last) {
	      result[i]->dist= theta= temp->rec_frac[i-start-1][MALE];
	      if (theta<ZERO_DIST) result[i]->zero=TRUE;
	  } else {
	      theta1= temp->rec_frac[i-start-1][MALE];
	      theta2= temp->rec_frac[i-start][MALE];
	      if (theta1<ZERO_DIST || theta2<ZERO_DIST) result[i]->zero=TRUE;
	      dist1= (*mapfunction->rec_to_dist)(theta1);
	      dist2= (*mapfunction->rec_to_dist)(theta2);
	      distframe= (*mapfunction->rec_to_dist)(map->rec_frac[i-1][0]);
	      distscale= dist1*(distframe/(dist1+dist2));
	      result[i]->dist= (*mapfunction->dist_to_rec)(distscale);
	  }

	  /* determine the error lods for any particular placement. 
	     Here is a stupid algorithm which looks only at this locus's
	     lods. we realy need to look at flankers too! */
	  result[i]->net_error=0.0; result[i]->worst_error=0.0;
	  if (temp->allow_errors && (i>start || i<last) &&
	      temp->error_rate[i-start]>0.0) {
	      for (j=0; j<raw.data.f2.num_indivs; j++) { /* must be F2 */
		  lod= temp->error_lod[i-start][j];
		  if (lod>result[i]->worst_error) result[i]->worst_error=lod;
		  if (lod>=error_lod_thresh)      result[i]->net_error+=lod;
	      }
	  }
      } /* next interval i */
    if (best_like==VERY_UNLIKELY) send(CRASH);

    /* normalize the likes */
    for (i=start; i<=last; i++) if (!excluded[i]) {
	if (i==*best_pos) result[i]->like=0.0; 
	else result[i]->like-=best_like;
    }
	
    /* if all placements with nearly 0 like are zero dist from a framework 
       locus, remove all but the best (thus place_this will place the marker 
       uniquely) */
    /* here we have a simple (and likely not complete enough) test... 
       are the likelihoods in consecutive intervals really zero, 
       and does each placement have SOME zero dist. */
    if (result[*best_pos]->zero) {
	for (i=*best_pos+1; !excluded[i] && i<=last &&
	     result[i]->like>ZERO_PLACE && result[i]->zero; i++) 
	  result[i]->like=ZERO_LIKE;
	for (i=*best_pos-1; !excluded[i] && i>=start && 
	     result[i]->like>ZERO_PLACE && result[i]->zero; i--)
	  result[i]->like=ZERO_LIKE;
    }
}


int npt_exclusions(place,num_loci,like_thresh,one_err_thresh,net_err_thresh,
		   excluded,zero_place,off_ends,error_lod,single_error,
		   next_best,best_unallowed)
/* return #ints ok */
PLACE **place;  /* results from above */
int  num_loci;  /* in order */
real like_thresh, one_err_thresh, net_err_thresh;
bool *excluded;     /* contents side-effected */
bool *zero_place;   /* side-effected, T if all acceptable places are zeros */
int  *off_ends;     /* side-effected, >0 if all acceptable places are end(s) */
real *error_lod;    /* side-effected, is greatest for all acceptable places */
bool *single_error; /* side-effected, T if error_place && err due to summing */
real *next_best;    /* side-effected, best acceptable but non_zero like (<0) */
real *best_unallowed; /* side-effected, best like >threshold */
{
    int i, best_i, num, zero, one, ends;
    real second, unallowed, this_like, err, x;
    if (like_thresh>=0.0) send(CRASH);

    num=0; zero=TRUE; err=NO_ERRORS; one=FALSE, ends=2; best_i= -1;
    second=NO_LIKE; unallowed=NO_LIKE;  /* may stay this way */

    for (i=0; i<=num_loci; i++) {
	this_like=place[i]->like;
	if (this_like!=NO_LIKE && this_like!=ZERO_LIKE && 
	      this_like>like_thresh) { /* it's acceptable */
	    num++;
	    if (excluded!=NULL) excluded[i]=FALSE;
	    if (this_like==0.0 && best_i==-1) best_i=i;
	      else if (this_like>second) second=this_like;
	    if ((x=place[i]->worst_error)>=one_err_thresh) { err=x; one=TRUE; }
	      else if ((x=place[i]->net_error)>=net_err_thresh) err=x;

	} else { /* it's unacceptable */
	    if (excluded!=NULL) excluded[i]=TRUE;
	    if (i==0 || i==num_loci) ends--;
	    if (this_like!=NO_LIKE) {
		if (this_like!=ZERO_LIKE) zero=FALSE;
		if (this_like>unallowed)  unallowed=this_like;
	    }
	}
    }
    if (best_i==-1) send(CRASH);
    if (place[best_i]->worst_error>=one_err_thresh) one=TRUE; /* override */
    
    if (zero_place!=NULL)     *zero_place= (num==1 && zero);
    if (error_lod!=NULL)      *error_lod=err;
    if (single_error!=NULL)   *single_error=one;
    if (off_ends!=NULL)       *off_ends=(ends>0 ? TRUE:FALSE);
    if (next_best!=NULL)      *next_best=second;
    if (best_unallowed!=NULL) *best_unallowed= unallowed;
    return(num);
}


/**************** Automatic Sub-Order Finder and Stuff ****************/


void informative_subset(locus,num_loci,min_infs,min_theta,codom_only,haplos,
			subset,num)
int *locus, num_loci;
int min_infs;
real min_theta;
bool codom_only, haplos;
int *subset, *num;
{
    int n_infs, n_dom_obs, n_het_obs, type, i, j, a, b, delete, deleted;
    int a_infs, b_infs;
    real lod, theta;
    *num=0;

    /* filter for informativeness */
    for (i=0; i<num_loci; i++) {
	f2_genotype(locus[i],haplos,observations);
	n_infs=f2_count_infs(&n_dom_obs,&n_het_obs,observations);
	type=raw.data.f2.cross_type;

	if (n_infs<min_infs) continue;
	if (type==F3_SELF || type==F2_INTERCROSS) {
	    if (codom_only && n_dom_obs!=0) continue;
	}
	subset[(*num)++]=locus[i];
    }
    if (*num<2) return;

    /* find zeros - a bit kludgey */
    do {
	deleted=FALSE;
	for (i=0; i<*num-1; i++) {
	    for (j=i+1; j<*num; j++) {
		a=subset[i]; b=subset[j];
		get_two_pt(a,b,&lod,&theta);
		if (theta<min_theta) {
		    /* delete least informative marker, to bad we don't save
		       the infomativeness values */
		    f2_genotype(a,haplos,observations);
		    a_infs=f2_count_infs(&n_dom_obs,&n_het_obs,observations);
		    f2_genotype(b,haplos,observations);
		    b_infs=f2_count_infs(&n_dom_obs,&n_het_obs,observations);
		    /* do something smarter with mix/dr markers */
		    if (a_infs>b_infs) delete=(coin_flip() ? i:j);
		      else if (a_infs>b_infs) delete=j;
		      else delete=i;
		    if (delete<*num-1) subset[delete]=subset[*num-1];
		    *num-=1; deleted=TRUE; break;
		}
	    }
	    if (deleted) break;
	}
    } while (deleted);
}



#define irand(n) (int)(randnum()*(real)n)
#define SEED_FAILED \
  "Failed to find a starting order at log-likelihood threshold %.2lf.\n"
#define SEARCH_INF \
"Searching for a unique starting order containing %d of %d informative loci...\n"
#define SEARCH_ALL \
  "Searching for a starting order containing %d of all %d loci...\n"
#define MAX_SAFETY 100

bool find_seed_order(is_subset,locus,num_loci,size,max_tries,thresh,map,
		     temp_map,temp)
bool is_subset;
int *locus, num_loci, size, max_tries;
real thresh;
MAP *map, *temp_map;
bool **temp; /* [num_loci][num_loci] */
{
    int try, pick, match, safety, i, j, k;
    /* size must be >=num_loci; if size==num_loci, there is only one choice.
       if size==num_loci-1, there are only num_loci choices, else there are
       many orders. if size>num_loci it's an error */

    if (num_loci<3) send(CRASH);
    if (size>num_loci) size=num_loci;
    try=0;
    temp_map->num_loci=size; map->num_loci=0;

    if (is_subset) sf(ps,SEARCH_INF,size,num_loci); 
      else sf(ps,SEARCH_ALL,size,num_loci);
    pr();

    if (size==num_loci) {
	for (i=0; i<num_loci; i++) temp_map->locus[i]=locus[i];
	keep_user_amused("subset",try+1,0);
	if (good_seed(map,temp_map,thresh)) return(TRUE);
	sf(ps,SEED_FAILED,-thresh); pr();
	return(FALSE);

    } else if (size==num_loci-1) {
	for (j=0; j<num_loci; j++) { /* locus to punt */
	    k=0;
	    for (i=0; i<num_loci; i++)
	      if (i!=j) temp_map->locus[k++]=locus[i];
	    keep_user_amused("subset",++try,0);
	    if (good_seed(map,temp_map,thresh)) return(TRUE);
	    if (try>max_tries) break;
	}
	sf(ps,SEED_FAILED,-thresh); pr();
	return(FALSE);
	
    } else /* size << num_loci */ {
	do { /* try a subset */
	    safety=0;
	    do { /* find an untried subset */
		/* first use temp as a flag to indicate loci in order */
		for (i=0; i<num_loci; i++) temp[try][i]=FALSE;
		for (i=0; i<size; i++) {
		    do pick=irand(num_loci); while (temp[try][pick]);
		    temp_map->locus[i]=locus[pick]; temp[try][pick]=TRUE;
		}
		/* have we tried this subset before? */
		sort_loci(temp_map->locus,size);
		for (j=0, match=FALSE; j<try; j++) {
		    for (match=TRUE, k=0; k<size; k++) 
		      if (temp_map->locus[k]!=temp[j][k]) { match=FALSE;break;}
		    if (match) break;
		}
		/* safety is a very lame way to see when we should give up, 
		   because we have tried all possible subsets */
		safety++;
	    } while (match && safety<=MAX_SAFETY);
	    if (match) break;
	    for (k=0; k<size; k++) temp[try][k]=temp_map->locus[k];
	    keep_user_amused("subset",++try,0);
	    /* for (i=0; i<size; i++)
	       { sf(ps,"%d ",temp_map->locus[i]+1); pr(); } nl(); */
	    if (good_seed(map,temp_map,thresh)) return(TRUE);
	} while(try<=max_tries);
	sf(ps,SEED_FAILED,-thresh); pr();
	return(FALSE);
    }
}


bool good_seed(map,temp_map,thresh)
MAP *map, *temp_map;
real thresh;
{
    real best2=VERY_UNLIKELY, best=VERY_UNLIKELY;
    
    make_compare_seq(temp_map->locus,temp_map->num_loci,0,temp_map->num_loci);
    for_all_orders(seq,temp_map) {
	if (use_three_pt &&
	    !three_pt_verify(temp_map->locus,temp_map->num_loci,
			     three_pt_window)) continue;
	init_rec_fracs(temp_map);
	converge_to_map(temp_map);
	if (temp_map->log_like>best) 
	  { best2=best; best=temp_map->log_like; mapcpy(map,temp_map,TRUE); }
	else if (temp_map->log_like>best2) { best2=temp_map->log_like; }
    }
    if (best==VERY_UNLIKELY) return(FALSE);
    else if (best2==VERY_UNLIKELY || best2-best<thresh) {
	sf(ps,"Got one at log-likelihood %.2lf\n",best-best2); pr();
	return(TRUE); 
    } else return(FALSE);
}
    


/**************** Automatic Order Maker ****************/

#define BY_NPT 1
#define BY_3PT 0
#define BY_NPT_ENDS (-1)
#define BY_NPT_ERRORS (-2)
#define BY_3PT_ENDS (-3)
#define NO_PLACES 9999

#define ALL_EXCLUDED \
"Three-point analysis excludes %d %smarker%s from all intervals\n"

#define sf_locnum(locus,paren) \
  sf(ps,"%s%d%s%s",(paren ? "(":""),(locus)+1, \
     (use_haplotypes && haplotyped(locus) ? "+":""),(paren ? ")":""))
     
void extend_order(placed,unplaced,num_unplaced,npt_thresh,print_anyway)
MAP *placed;  		/* starting order, must be total long */
PLACEME **unplaced;
int *num_unplaced;
real npt_thresh;        /* both are side-effected, and may have NO_LOCUS */
bool print_anyway;
{
    int i, j, k, num, total;
    bool placed_any, contradiction;
    PLACE   **placements=NULL;
    MAP     *temp_map=NULL;
    
    run {
	total=placed->num_loci + *num_unplaced;
	/* if (*num_unplaced==0) return; let it fall through */
	if (placed->num_loci<2) send(CRASH);
	
	temp_map=allocate_map(total);
	parray(placements,total,PLACE);
	for (i=0; i<*num_unplaced; i++) {
	    unplaced[i]->best_pos=(-1);
	    unplaced[i]->num_places=0;
	    unplaced[i]->best_map->num_loci=0;
	}

	print("Start:   ");
	for (j=0; j<placed->num_loci; j++) {
	    sf_locnum(placed->locus[j],FALSE); pr();
	    if (j!=placed->num_loci-1) print(" ");
	} nl();
	
	placed_any=FALSE;
	while (*num_unplaced>0)
	  if (!extend_order_by_1(placed,unplaced,num_unplaced,npt_thresh,
				 &contradiction,placements,temp_map)) {
	      sf(ps,"No unique placements for %d %smarker%s\n",
		 *num_unplaced,(placed_any ? "remaining ":""),
		 maybe_s(*num_unplaced)); pr();
	      if (contradiction==TRUE) {
		  sf(ps,"Three-point analysis allows no %s.\n",
		     (placed_any ? "further ordering":"orders at all")); pr();
	      } else if (contradiction==MAYBE) {
		  print("Maximum number of loci in one map reached.\n");
	      }
	      break;
	  } else placed_any=TRUE;
	if (*num_unplaced==0)
	  { sf(ps,"Uniquely ordered all %d markers\n\n",total); pr(); }
	
	if (print_anyway || placed_any) {
	    init_not_fixed(placed); 
	    converge_to_map(placed);
	    nl();
	    print_long_map(placed,"Map:");
	    if (*num_unplaced>0) {
		nl();
		print("Markers placed relative to above map:\n");
		new_print_placements(placed,unplaced,*num_unplaced);
	    }
	}
    } on_exit {
	unparray(placements,total,PLACE);
	free_map(temp_map);
	relay_messages;
    }
}


bool extend_order_by_1(placed,unplaced,num_unplaced,npt_thresh,contradiction,
		       placements,temp_map)
MAP *placed;
PLACEME **unplaced;
int  *num_unplaced;
real npt_thresh;
bool *contradiction; /* iff return FALSE, is TRUE if 3pt excludes all */
PLACE **placements; /* these are all temps */
MAP *temp_map;
{
    int places, best_places, best_pos, best_i, i, j;
    int left, right, pos, how;
    real next_best, best_unallowed, error_lod, e1, e2; /* e1, e2=threshs */
    bool zero_place, off_ends, single_error;

    if (placed->num_loci>=MAX_MAP_LOCI) 
      { *contradiction=MAYBE; return(FALSE); }
    randomize_markers_to_try(unplaced,*num_unplaced);
    for (i=0; i<*num_unplaced; i++) unplaced[i]->best_map->num_loci=0;

    /* Try to place remaining loci, into order, in their order */
    if (use_three_pt) {
	best_places=NO_PLACES; best_i=-1;
	for (i=0; i<*num_unplaced; i++) { /* sorted */
	    unplaced[i]->num_places= places=
	      three_pt_exclusions(placed->locus,placed->num_loci,
				  unplaced[i]->locus,unplaced[i]->excluded);
	    unplaced[i]->best_pos=(-1); /* NO_POS? */
	    if (places>0 && places<best_places) {
		best_places=places; best_i=i;
		if (places==1 && !middles_excluded(unplaced[i]->excluded,
						   placed->num_loci)) {
		    /* It doesn't get any better than this! */
		    add_to_order(BY_3PT,placed->locus,&placed->num_loci,
				 unplaced,best_i,num_unplaced);
		    return(TRUE);
		}
	    }
	} /* for i: loop over num_to_place trying 3pt */
	if (best_places==NO_PLACES) { *contradiction=TRUE; return(FALSE); }

	/**** We get to here if all allowed positions are off-end, or if no 
	      unique placements were found ****/
#ifdef DONT_DO_THIS
	for (i=0; i<*num_unplaced; i++) if (unplaced[i]->num_places==1) {
	    /* a troublesome however unique placement */
	    add_to_order(BY_3PT_ENDS,placed->locus,&placed->num_loci,unplaced,
			 i,num_unplaced);
	    return(TRUE);
	}   
#endif
	/* otherwise fall-through to Npt tests below */
	
    } else { /* !use_three-pt */
	for (i=0; i<*num_unplaced; i++) {
	    /* unplaced[i]->num_places= placed->num_loci+1; */
	    for (j=0; j<placed->num_loci; j++) unplaced[i]->excluded[j]=FALSE;
	}
    }

    /**** Have No UNIQUE Placements, so try N-pt analysis ****/
    for (i=0; i<placed->num_loci-1; i++) 
      placed->rec_frac[i][MALE]=placed->rec_frac[i][FEMALE]=NOT_FIXED;
    init_rec_fracs(placed); converge_to_map(placed);
    rank_markers_to_try(unplaced,*num_unplaced);

    best_places=NO_PLACES; best_i=-1;
    for (i=0; i<*num_unplaced; i++) { /* sorted */
	keep_user_amused("marker",i+1,*num_unplaced);
	if (unplaced[i]->num_places==0) continue;
	find_window(placed->locus,placed->num_loci,
		    unplaced[i]->locus,unplaced[i]->excluded,
		    npt_window,&left,&right);
	place_locus(placed,unplaced[i]->locus,left,right,
		    unplaced[i]->excluded,placements,&unplaced[i]->best_pos,
		    unplaced[i]->best_map,temp_map);
	unplaced[i]->num_places= places=
	  npt_exclusions(placements,placed->num_loci,npt_thresh,
			 error_single_thresh,error_net_thresh, /* GLOBALS */
			 unplaced[i]->excluded,&zero_place,&off_ends,
			 &error_lod,&single_error,&next_best,&best_unallowed);
	unplaced[i]->off_end= off_ends;
	unplaced[i]->errors= (error_lod!=NO_ERRORS);

	/* preferentially add loci which place well */
	if (places<best_places && !unplaced[i]->off_end && 
	        !unplaced[i]->errors) {
	      best_places=places; best_i=i;
	    if (places==1) { /* Greedy: It doesn't get any better than this! */
		add_to_order(i+1,placed->locus,&placed->num_loci,unplaced,
			     best_i,num_unplaced);
		return(TRUE);
	    }
	}
    } /* for i: loop over num_to_place trying 3pt */
    
    /**** We get to here if all placements are troublesome, or if no unique
          placements were found ****/
    for (i=0; i<*num_unplaced; i++) if (unplaced[i]->num_places==1) {
	/* a troublesome however unique placement */
	how= (unplaced[i]->off_end ? BY_NPT_ENDS:BY_NPT_ERRORS);
	add_to_order(how,placed->locus,&placed->num_loci,unplaced,i,
		     num_unplaced);
	return(TRUE);
    }

    /**** Else we have no unique placements, we fail ****/
    *contradiction=FALSE; return(FALSE);
}


void randomize_markers_to_try(unplaced,num_unplaced)
PLACEME **unplaced;
int num_unplaced;
{
    int i, n_dom, n_het, num;

    for (i=0; i<num_unplaced; i++) {
	unplaced[i]->num_places=0;
	f2_genotype(unplaced[i]->locus,use_haplotypes,observations);
	num=f2_count_infs(&n_dom,&n_het,observations);
	unplaced[i]->priority= 
	  imaxf(num-2*n_dom,0) +
	  irand(raw.data.f2.num_indivs/10);
    }
    qsort((QSORT_DATA_PTR_TO*)unplaced,(QSORT_LENGTH)num_unplaced,
	  sizeof(PLACEME*),compare_markers_to_try);
}


void rank_markers_to_try(unplaced,num_unplaced)
PLACEME **unplaced;
int num_unplaced;
{
    qsort((QSORT_DATA_PTR_TO*)unplaced,(QSORT_LENGTH)num_unplaced,
	  sizeof(PLACEME*),compare_markers_to_try);
}


int compare_markers_to_try(a,b)
QSORT_COMPARE_PTR_TO(PLACEME*) *a;
QSORT_COMPARE_PTR_TO(PLACEME*) *b;
{
    PLACEME *x, *y;
    
    x=*a; y=*b;
    /* num_places is three_pt criteria */
    if      (x->num_places < y->num_places) return(-1);
    else if (x->num_places > y->num_places) return(1);
    else if (x->priority > y->priority) return(-1);
    else if (x->priority < y->priority) return(1);
    else return(0);
}


void add_to_order(how,order,num_ordered,unplaced,index,num_unplaced)
int how, *order, *num_ordered;
PLACEME **unplaced;
int index, *num_unplaced;
{
    int i, pos, last, new_marker;
    PLACEME *temp;
    new_marker=unplaced[index]->locus;

    for (pos=0; pos<=*num_ordered; pos++)
      if (!unplaced[index]->excluded[pos]) break;
    if (pos==*num_ordered+1) send(CRASH);

    new_marker=unplaced[index]->locus;
    if (pos<*num_ordered)
      for (i=*num_ordered-1; i>=pos; i--) order[i+1]=order[i];
    order[pos]=new_marker;
    *num_ordered+=1;

    last=*num_unplaced-1;
    if (*num_unplaced==1) *num_unplaced=0;
    else if (index==last) *num_unplaced-=1;
    else {
	temp=unplaced[index];
	unplaced[index]=unplaced[last];
	unplaced[last]=temp;
	*num_unplaced-=1;
    }
    
#ifndef DEBUGGING
    if (how>0)                        { sf(ps,"Npt-%d:\t ",how); pr(); }
#else
    if (how>0)                        { sf(ps,"Npt:\t "); pr(); }
#endif
      else if (how==BY_3PT) 		print("3pt:     ");
      else if (how==BY_NPT_ENDS) 	print("Npt-End: ");
      else if (how==BY_NPT_ERRORS) 	print("Npt-Err: ");
      else if (how==BY_3PT_ENDS) 	print("3pt-End: ");
      else send(CRASH);
    to_column(9);

    for (i=0; i<*num_ordered; i++) { /* print new order */
	sf_locnum(order[i],(order[i]==new_marker)); pr();
	if (i!=*num_ordered-1) print(" ");
    } nl();
}


/* to keep this all straight:
   orig frame:         a - b - c - d - e - f - g
   orig frame index:   0 - 1 - 2 - 3 - 4 - 5 - 6    num_placed=7
   exclusion index:  0 - 1 - 2 - 3 - 4 - 5 - 6 - 7 
   new frame:              b - c - d - e - f    
   new frame index:        0 - 1 - 2 - 3 - 4        */

void find_window(order,num_placed,new_locus,excluded,min_window,start,end)
int *order; /*NOTUSED*/
int num_placed, new_locus, *excluded;
int min_window;    /* #of loci to left and right, should be odd, e.g. 5 */
int *start, *end;  /* side-effected with interval indecies */
{
    int j, step;
    /* simple implementation for now - this should later use informativeness */

    step= max((min_window-1)/2,1);
    for (j=0; j<=num_placed && excluded[j]; j++) {}
    *start= max(j-step,0);
    for (j=num_placed; j>0 && excluded[j]; j--) {}
    *end= min(j+step-1,num_placed-1);   /* Is this right? */
}


int count_positions(excluded,num_order,first_place)
int *excluded, num_order, *first_place;
{
    int j, num_places;

    for (j=0, num_places=0; j<=num_order; j++) if (!excluded[j]) {
	if (num_places==0 && first_place!=NULL) *first_place=j;
	num_places++;
    }
    return(num_places);
}


bool middles_excluded(excluded,num_order)
int *excluded, num_order;
{
    int j;

    for (j=1; j<num_order; j++) if (!excluded[j]) return(FALSE);
    return(TRUE);
}


#ifdef OBSOLETE /************************************************************/
void make_npt_maps(order,start,finish,newmarker,threshold,excluded,map,like)
int *order, start, finish, newmarker;
real threshold;    /* threshold is negative */
bool *excluded;
MAP *map;
real *like;
{
    int i, j, k, best_i, maybe_zero, first;
    real best_like;

    best_like= VERY_UNLIKELY;
    for (i=start; i<=finish+1; i++) /* the interval #s */
      if (!excluded[i]) {
	  clean_map(map); k=0;
	  if (i==0) map->locus[k++]=newmarker;
	  for (j=start; j<=finish; j++) { /* locus #s */
	      map->locus[k++]=order[j];
	      if (j+1==i) map->locus[k++]=newmarker;
	  }
	  map->num_loci= k;
	  
	  init_rec_fracs(map);
	  converge_to_map(map);
	  like[i]=map->log_like;
	  if (like[i]>best_like) { best_like=like[i]; best_i=i; }
      }
	
    for (i=start; i<=finish+1; i++) if (!excluded[i]) {
	if (like[i] < best_like+threshold) excluded[i]=TRUE;
	like[i]-= best_like; /* check this! */
    }

    /* if all likes are zero just pick one - NEED ALSO TEST FOR DIST */
    maybe_zero=TRUE;
    for (i=start; i<=finish+1; i++) if (!excluded[i]) {
	if (like[i]<=ZERO_PLACE) maybe_zero=FALSE;
    } 
    first=TRUE;
    if (maybe_zero) for (i=start; i<=finish+1; i++) if (!excluded[i]) {
	if (first) first=FALSE; else excluded[i]=TRUE;
    }
}
#endif /* OBSOLETE */


#ifdef EXPERIMENTAL /*******************************************************/

/* to keep this all straight:
   orig frame:         a - b - c - d - e - f - g
   orig frame index:   0 - 1 - 2 - 3 - 4 - 5 - 6
   exclusion index:  0 - 1 - 2 - 3 - 4 - 5 - 6 - 7
   new frame:              y - c - d - e - z    
   new frame index:        0 - 1 - 2 - 3 - 4   
   place_in:                 0 - 1 - 2 - 3               */

void make_pair_maps(order,start,finish,locus,i1,i2,threshold,excluded,
		   map,temp_order)
int *order, start, finish;
int *locus, i1, i2;  /* the locus array and indecies of the two to place */
real threshold;      /* threshold is negative */
bool **excluded;
MAP *map;
int **temp_order;     /* used as a temp 1,2=frame, 3..., n... w/places */
{
    int j, k, n, n_in_order, n_start, n_orders;
    real best;

    temp_order[1][0]=NO_LOCUS; /* add two dummy loci */
    for (k=1, j=start; j<=finish; j++) temp_order[1][k++]=order[j];
    new_order[1][k]=NO_LOCUS; n_start=k+1;

    for (k=0, j=start; j<=finish+1; j++) place_in[k++]=!excluded[i1][j];
    for (k=0, j=start; j<=finish; j++)   temp_order[2][k++]=order[j];
    n_in_order=k;

    add_orders(temp_order[1],n_start,locus[i1],place_in,&temp_order[2],1,
	       n_in_order,&temp_order[3],&n_orders);

    ++n_in_order;
    for (k=0, j=start; j<=finish+1; j++) place_in[k++]=!excluded[i2][j];
    n= 3+n_orders;

    add_orders(temp_order[1],n_start,locus[i2],place_in,&temp_order[3],
	       n_orders,n_in_order,&temp_order[n],&n_orders);

    for (i=0; i<n_orders; i++) { 
	sf(ps,"   %d->\t",i+1); pr();
	for (j=0; j<n_in_order+1; j++) 
	  { sf(ps,"%d ",temp_order[n+i][j]); pr(); }
	nl();
    }

    n_in_order-=1; /* ignore dummies on end */
    best= VERY_UNLIKELY;
    for (i=0; i<n_orders; i++) {
	clean_map(map);
	for (j=0; j<n_in_order; j++) 
	  map->locus[j]= temp_order[n+i][j+1]; /* watch dummies on end */
	map->num_loci= n_in_order;
	init_rec_fracs(map);
	converge_to_map(map);
    }




void add_orders(frame,n_frame,new_marker,place_in,old_order,n_orders,n_loci,
		new_order,num_orders)
int *frame;      /* frame must have LEFT_DUMMY@0 and RIGHT_DUMMY@n-1 */
int n_frame;     /* thus n_frame is real frame len+2 */
int new_marker;  
bool *place_in;  /* [n]=>T/F, n is the interval right of frame[n], 0..n_f-1 */
int **old_order; /* [order#][position#] */
int n_orders;    /* #old orders */
int n_loci;      /* in each old order */
int **new_order; /* [order#][position#] */  
int *new_orders; /* ptr to int - side-effected */
{
    int n, i, j, k, left, right;
    
    k=0;
    for (n=0; n<n_orders; n++) {
	for (i=0; i<n_frame-1; i++) if (place_in[i]) { /* never off end */
	    left=frame[i]; right=frame[i+1];
	    /* j= position of left marker in old_order[n] */
	    for (j=0; old_order[n][j]!=left; j++) {}

	    do { /* until j hits right marker in old_order[n] */
		add_order(old_order[n],n_loci,j,new_marker,new_order[k++]);
	    } while (old_order[n][++j]!=right);
	}
    }
    *new_orders= k;
}


void add_order(old_order,n_loci,pos,new_marker,new_order)
int *old_order, n_loci, pos, new_marker, *new_order;
{
    int i, d;

    for (i=0, d=0; i<n_loci; i++) {
	new_order[i+d]= old_order[i];
	if (i==pos) { new_order[i+1]= new_marker; d=1; }
    }
}


/* needs excluded, a little tweaking */
int **count_orders(frame,n_frame,new_marker,place_in,old_order,n_orders,n_loci,
		   num_orders)
int *frame;      /* frame must have LEFT_DUMMY@0 and RIGHT_DUMMY@n-1 */
int n_frame;     /* thus n_frame is real frame len+2 */
int new_marker;  
bool *place_in;  /* [n]=>T/F, n is the interval right of frame[n], 0..n_f-1 */
int **old_order; /* [order#][position#] */
int n_orders;    /* #old orders */
int n_loci;      /* in each old order */
int **new_order; /* [order#][position#] */  
int *new_orders; /* ptr to int - side-effected */
{
    int n, i, j, k, left, right;
    
    k=0;
    for (n=0; n<n_orders; n++) {
	for (i=0; i<n_frame-1; i++) if (place_in[i]) { /* never off end */
	    left=frame[i]; right=frame[i+1];
	    /* j= position of left in old_order[n] */
	    for (j=0; old_order[n][j]!=left; j++) {}
	    
	    do { /* until we hit right in old_order[n] */
		add_order(old_order[n],n_loci,j,new_marker,new_order[k++]);
	    } while (old_order[n][++j]!=right);
	}
    }
    *new_orders= k;
}
#endif  /* EXPERIMENTAL */

