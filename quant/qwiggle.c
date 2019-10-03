/******************************************************************************

  ####   #    #     #     ####    ####   #       ######           ####
 #    #  #    #     #    #    #  #    #  #       #               #    #
 #    #  #    #     #    #       #       #       #####           #
 #  # #  # ## #     #    #  ###  #  ###  #       #        ###    #
 #   #   ##  ##     #    #    #  #    #  #       #        ###    #    #
  ### #  #    #     #     ####    ####   ######  ######   ###     ####

******************************************************************************/

/* qwiggle.c - code for global wiggle/compare data storage, save & load */

#define INC_LIB
#define INC_SHELL
#define INC_CALLQCTM
#define INC_QTOPLEVEL
#define INC_QLOWLEVEL
#include "qtl.h"

/* globals */
/* WIGGLE_INTERVAL ****qtls; */
WIGGLE_OPERATION **wiggles;   
int num_wiggles, max_wiggles, first_wiggle;
COMPARE_OPERATION **compares;   
int num_compares, max_compares, first_compare;

/* local to this file */
bool step();
WIGGLE_PEAK *peak_finder();
void save_wiggle(),save_interval(),save_qtl_map(),save_compare();
void load_wiggle(),load_interval(),load_compare();
QTL_MAP *load_qtl_map();

/* functions */

void wiggle_init()
{
    /* qtls=NULL; */
    wiggles=NULL;
    compares=NULL;
    num_wiggles= num_compares= 0;
    first_wiggle= first_compare= 0;
    max_wiggles= max_compares= 0;  /* WHAT DEFAULT TO USE? */
}



void allocate_qtl_struct(n_wiggles,n_compares) /* allocate globals */
int n_wiggles, n_compares;
/* use globals raw.max_traits, raw.data_type, and raw.n_loci as params */
{
    int i, j, k, n_models, n_intervals;

    if (raw.data_type==NO_DATA) send(CRASH);
    wiggles=NULL; compares=NULL; /* qtls=NULL; */

    run {
	parray(wiggles,n_wiggles,WIGGLE_OPERATION); 
	max_wiggles=n_wiggles; num_wiggles=0;
	for (i=0; i<n_wiggles; i++) wiggles[i]->data=NULL;

	parray(compares,n_compares,COMPARE_OPERATION); 
	max_compares=n_compares; num_compares=0;
	for (i=0; i<n_compares; i++) compares[i]->data=NULL;

	n_models= (raw.data_type==INTERCROSS ? NUM_MODELS:1);
	n_intervals= raw.n_loci-1;
	/*  matrix(qtls,raw.n_traits,n_models,WIGGLE_INTERVAL**);
	    for  (i=0; i<raw.n_traits; i++) for (j=0; j<n_models; j++) {
	    array(qtls[i][j],n_intervals,WIGGLE_INTERVAL*);
	    for (k=0; k<n_intervals; k++) qtls[i][j][k]=NULL;
	    }  */

    } except_when(NOMEMORY) {
	unparray(wiggles,n_wiggles,WIGGLE_OPERATION);
	unparray(compares,n_compares,compare_operation);
	/* if (qtls!=NULL) 
	   for (i=0; i<raw.n_traits; i++) for (j=0; j<n_models; j++)
	   unparray(qtls[i][j],n_intervals,WIGGLE_INTERVAL);
	   unmatrix(qtls,raw.n_traits,WIGGLE_INTERVAL**); */
	wiggles=NULL; compares=NULL; /* qtls=NULL; */
	num_wiggles=max_wiggles=num_compares=max_compares=0;
	relay;
    }
}


/***** Stuff for the wiggle struct to save wiggle output... *****/

int allocate_wiggle_struct(trait,seq,seq_str,n_intervals,n_orders,n_wiggled)
int trait;
QTL_SEQUENCE *seq;
char *seq_str;
int n_intervals, n_orders, n_wiggled;
/* return wiggles struct entry#, or -1 if error */
{
    int n, i, j;

    /* KLUDGE - How to recycle these things? */
    if (num_wiggles==max_wiggles) return(-1);

    n= num_wiggles++;
    wiggles[n]->data=NULL; 
    wiggles[n]->seq_string=NULL; 
    run {
	if (trait == -1) {
	    wiggles[n]->trait= -1;
	    wiggles[n]->seq = NULL;
	    wiggles[n]->num_intervals = -1;
	    wiggles[n]->num_orders = -1;
	    wiggles[n]->num_wiggled_intervals = -1;
	    wiggles[n]->order = 0;
	    wiggles[n]->wiggle_interval = -1;
	    wiggles[n]->seq_string = NULL;
	}
	else {
	    wiggles[n]->trait=trait;
	    wiggles[n]->seq= seq; seq->dont_free=TRUE;
	    wiggles[n]->num_intervals=n_intervals;
	    wiggles[n]->num_orders= n_orders;
	    wiggles[n]->num_wiggled_intervals= n_wiggled;
	    wiggles[n]->order= -1;  /* For store_wiggle_interval */
	    wiggles[n]->wiggle_interval= 0;
            /********** 
              would like to expand the sequence so that changes to user-defined
              names do not ruin any saved scan information, but this can't be 
	      done here - too many problems - just DON'T edit user-def names!
              wiggles[n]->seq_string= expand_named_entries(seq_str);
	    **********/
	    wiggles[n]->seq_string= mkstrcpy(seq_str);
	
	    matrix(wiggles[n]->data,n_orders,n_wiggled,WIGGLE_INTERVAL*);
	    for (i=0; i<n_orders; i++) for (j=0; j<n_wiggled; j++) {
		single(wiggles[n]->data[i][j],WIGGLE_INTERVAL);
		wiggles[n]->data[i][j]->point=NULL; 
		wiggles[n]->data[i][j]->map=NULL;
	    }
	}
    } when_aborting { bash_wiggle_struct(n); relay; }
    return(n);
}


void bash_wiggle_struct(n) 
/* this is a KLUDGE for now - we need to recycle these better */
int n;
{
    unmatrix(wiggles[n]->data,wiggles[n]->num_orders,WIGGLE_INTERVAL);
    unarray(wiggles[n]->seq_string,char);

    /* change this if any functions besides allocate_wiggle_struct()
       or allocate_compare_struct() could set seq->dont_fix! */
    wiggles[n]->seq->dont_free=FALSE; 
    wiggles[n]->data=NULL; /* This is the flag that says 'not filled in' */
    wiggles[n]->seq_string=NULL;
    if (n==num_wiggles-1) num_wiggles--;
} 


void store_wiggle_interval(wiggle_num,map,new_left_order,contig,cm_step)
int wiggle_num;
QTL_MAP *map;
bool new_left_order;
bool contig;
real cm_step;
{
    WIGGLE_OPERATION *op;
    WIGGLE_INTERVAL *data, *prev;
    int i, j, k, n, t;
    real cm, interval_cm;
    bool new;

    if (wiggle_num<0) return;

    if (wiggle_num>=max_wiggles || wiggles[wiggle_num]==NULL ||
	(op=wiggles[wiggle_num])->data==NULL || map==NULL) send(CRASH);

    if (op->order>=0) { /* prev order exists */
	prev=op->data[op->order][op->wiggle_interval];
	if (prev==NULL || prev->point_num!=prev->num_points-1) send(CRASH);
    }

    if (new_left_order) {
	if (op->order>=0) {
	    if (op->order>=op->num_orders-1) send(CRASH);
	    if (op->wiggle_interval!=op->num_wiggled_intervals-1) send(CRASH);
	}
	j=(op->wiggle_interval=0); i=(op->order+=1); 	
    } else {
	if (op->wiggle_interval>=op->num_wiggled_intervals-1) send(CRASH);
	i=op->order; j=(op->wiggle_interval+=1);
    }
    
    k= map->num_intervals-1; /* the rightmost interval */
    data=op->data[i][j]; if (data==NULL) send(CRASH);
    data->point=NULL; data->map=NULL;
    run { 
	data->contig= contig;
	data->map= alloc_qtl_map(map->num_intervals,map->num_continuous_vars);
	mapcpy(data->map,map);
	data->max_like_map= FALSE;
	data->in_qtls_struct= FALSE;
	data->op= op; 
	data->cm_increment= cm_step;
	/* don't change this without changing the code in wiggle_perm() */
	interval_cm= haldane_cm(map_length(map->left[k],map->right[k]));
	for (cm=0.0, n=0; cm<interval_cm; cm+=cm_step) n++;
	data->num_points=n;
	parray(data->point,data->num_points,WIGGLE_POINT);
	data->point_num= -1; /* so that store_wiggle_point can increment */
	data->best_point=0;

    } when_aborting {
	free_qtl_map(data->map);
	unparray(data->point,data->num_points,WIGGLE_POINT);
	data->point= NULL; data->map=NULL; 
	relay;
    }

    /* This wiggle interval MAY want to go into the qtls struct also... */
    do { 
	if (map->num_intervals!=1) break;
	t= map->trait;
	if (raw.data_type==INTERCROSS) {
	    if ((i=map->constraint[0].interx_type)>=NUM_MODELS) break;
	} else { i=0; if (map->constraint[0].backx_weight!=DONT_FIX) break; }
	if ((j=map->left[0])!=map->right[0]-1) break;
#ifdef DELETED
	/* It does want to go in... */
	/* if (qtls[t][i][j]!=NULL) {  but these wiggle data already exist! */
	    /* prev= qtls[t][i][j]; new=FALSE;
	    /* If the old one has points, we keep the higher resolution one */
	    if (prev->point==NULL || data->cm_increment<=prev->cm_increment) {
		prev->in_qtls_struct=FALSE; 
		qtls[t][i][j]=data;
		data->in_qtls_struct=TRUE;
		new=TRUE;
	    }
	    /* If the max_like map, but not the points, was filled in 
	       in the old struct we merge the ML map into the new struct. */
	    if (new && prev->max_like_map) {
		mapcpy(data->map,prev->map);
		data->max_like_map= TRUE;
	    }
	    /* If ONLY the max_like_map exists in the old wiggle_interval 
	       struct then it should be freed. */
	   if (new   && prev->point==NULL) { 
		free_qtl_map(prev->map); prev->map=NULL;
		unsingle(prev,WIGGLE_INTERVAL);
	    }

	    else { /* No previous data exist. */
		qtls[t][i][j]=data;
		data->in_qtls_struct=TRUE;
	}
#endif
    } while(FALSE); /* not really a loop, just used do {} for convenience */
}


void store_wiggle_point(wiggle_num,map)
int wiggle_num;
QTL_MAP *map;
{
    WIGGLE_OPERATION *op;
    WIGGLE_INTERVAL *data;
    int i, j;

    if (wiggle_num<0) return;

    if (wiggle_num>=max_wiggles || wiggles[wiggle_num]==NULL ||
	(op=wiggles[wiggle_num])->data==NULL || map==NULL) send(CRASH);

    data=op->data[op->order][op->wiggle_interval]; 
    if (data==NULL || data->map==NULL) send(CRASH);
    if (data->point_num>=data->num_points-1) send(CRASH);
    i= (data->point_num+=1);
    j= map->num_intervals-1; /* the rightmost interval */

    data->point[i]->qtl_pos= map->qtl_pos[j];
    data->point[i]->lod_score= map->log_like;
    data->point[i]->var_explained= map->var_explained;
    data->point[i]->qtl_weight= map->qtl_weight[j];
    if (raw.data_type!=INTERCROSS) data->point[i]->qtl_dominance= 0.0;
      else data->point[i]->qtl_dominance= map->qtl_dominance[j];

    /* the "|| i == 0" was added to copy the first map into the best map
       as an initial condition - mjd 4/3/90 */

    if (data->point[data->best_point]->lod_score < map->log_like || i == 0)
      { data->best_point=i; mapcpy(data->map,map); }
}




/***** Stuff for the COMPARE struct to save compare output... *****/

int allocate_compare_struct(trait,seq,seq_str,n_intervals,n_orders)
int trait;
QTL_SEQUENCE *seq;
char *seq_str;
int n_intervals, n_orders;
/* return compare struct entry#, or -1 if error */
{
    int n, i, j;

    /* KLUDGE - How to recycle these things? */
    if (num_compares==max_compares) return(-1);

    n= num_compares++;
    compares[n]->data=NULL; compares[n]->seq_string=NULL;
    run {
	if (trait== -1){
	    compares[n]->trait= -1;
	    compares[n]->seq = NULL;
	    compares[n]->num_intervals= -1;
	    compares[n]->num_orders = -1;
	    compares[n]->num_contigs = -1;
	    compares[n]->order = 0;
	    compares[n]->seq_string = NULL;
	    compares[n]->data = NULL;
	}
	else {
	    compares[n]->trait=trait;
	    compares[n]->seq= seq; seq->dont_free=TRUE;
	    compares[n]->num_intervals=n_intervals;
	    compares[n]->num_orders= n_orders;
	    compares[n]->num_contigs= 0;
	    compares[n]->order= -1;  /* For store_compare_interval */
	    compares[n]->seq_string= mkstrcpy(seq_str);
	    parray(compares[n]->data,n_orders,COMPARE_DATA);
	}
    } when_aborting { bash_compare_struct(n); relay; }
    return(n);
}	


void bash_compare_struct(n) 
/* this is a KLUDGE for now - we need to recycle these better */
int n;
{
    free_qtl_map(compares[n]->data[0]->map);
    compares[n]->data[0]->map=NULL; 
    unparray(compares[n]->data,compares[n]->num_orders,COMPARE_DATA);
    unarray(compares[n]->seq_string,char);
    /* change this if any functions besides allocate_wiggle_struct()
       or allocate_compare_struct() could set seq->dont_fix! */
    compares[n]->seq->dont_free=FALSE; 
    compares[n]->data=NULL; /* This is the flag that says 'not filled in' */
    compares[n]->seq_string=NULL;
    if (n==num_compares-1) num_compares--;
} 


void store_compare_map(compare_num,map,contig)
int compare_num;
QTL_MAP *map;
bool contig;
{
    COMPARE_OPERATION *op;
    COMPARE_DATA *data;
    WIGGLE_INTERVAL *wig, *prev;
    int i, j, k, n, t;

    if (compare_num<0) return;
    if (compare_num>=max_compares || compares[compare_num]==NULL ||
	(op=compares[compare_num])->data==NULL || map==NULL) send(CRASH);

    if (op->order<0 || !contig) op->num_contigs++;
    i=(op->order+=1); if (i>=op->num_orders) send(CRASH);
    data=op->data[i]; if (data==NULL) send(CRASH);
    data->map=NULL;

    run {
	data->contig= contig;
	data->map= alloc_qtl_map(map->num_intervals,map->num_continuous_vars);
	mapcpy(data->map,map);
	data->op= op;
    } when_aborting {
	free_qtl_map(data->map);
	data->map=NULL;
	relay;
    }

    /* This map MAY want to go into the qtls struct also... */
    wig=NULL;
    do { run { 
	if (map->qtl_pos[0]!=DONT_FIX) break;
	if (map->num_intervals!=1) break;
	t= map->trait;
	if (raw.data_type==INTERCROSS) {
	    if ((i=map->constraint[0].interx_type)>=NUM_MODELS) break;
	} else { i=0; if (map->constraint[0].backx_weight!=DONT_FIX) break; }
	if ((j=map->left[0])!=map->right[0]-1) break;
#ifdef DELETED
	/* It does want to go in... */
	if (qtls[t][i][j]!=NULL) { /* but data for this qtl already exists! */
	    prev= qtls[t][i][j];
	    /* Regardless of whether the old one has a max-like qtl_map, 
	       we give it this one. */
	    prev->max_like_map= TRUE;
	    mapcpy(prev->map,data->map);
	} else { /* No previous data exist, so we alloc a dummy wiggle_int. */
	    wig=NULL;
	    single(wig,WIGGLE_INTERVAL); 
	    wig->map=NULL; wig->point=NULL;
	    wig->map= 
	      alloc_qtl_map(map->num_intervals,map->num_continuous_vars);
	    mapcpy(wig->map,data->map);
	    wig->in_qtls_struct=TRUE;
	    wig->max_like_map= TRUE;
	    qtls[t][i][j]=wig;
	}
#endif
    } when_aborting {
	if (wig!=NULL) free_qtl_map(wig->map);
	unsingle(wig,WIGGLE_INTERVAL);
	relay;
    } } while(FALSE); /* not really a loop, just used do {} for convenience */
}



bool step(interval,point,contig,wig,n_intervals,forwards)
int *interval, *point;
bool *contig;
WIGGLE_INTERVAL **wig;
int n_intervals;
bool forwards;
{
    if (*interval >= n_intervals) return(FALSE);
    if (forwards) {
	if (++*point<wig[*interval]->num_points) { *contig=TRUE; return(TRUE);}
	else if (++*interval<n_intervals) { 
	    *contig= wig[*interval]->contig;
	    *point=0; return(TRUE); 
	} else {
	    *contig = FALSE;
	    return(FALSE);
	}
    } else { /* backwards */
	if (--*point>=0) { *contig=TRUE; return(TRUE); }
	else if (--*interval>=0) {
	    if(*interval+1 == n_intervals) *contig = TRUE;
	    else *contig= wig[*interval + 1]->contig;
	    *point=wig[*interval]->num_points-1; return(TRUE); 
	} else {
	    *contig = FALSE;
	    return(FALSE);
	}
    }
}


WIGGLE_PEAK *find_wiggle_peaks(wiggle_num,left_order_num,threshold,qtl_falloff,
			       confidence_falloff,min_peak_delta,get_peak_maps)
int wiggle_num, left_order_num;
real threshold, qtl_falloff, confidence_falloff, min_peak_delta;
bool get_peak_maps;
{
    int i, j, peak_i, peak_j, n_intervals;
    bool still_falling, off_end, contig;
    real lod, local_maxima, pos, cm;
    real local_minima, prev_lod, starting_value;
    WIGGLE_INTERVAL **interval; 
    WIGGLE_PEAK *first, *last, *peak;

    if (wiggle_num<0) return(NULL);

    if (wiggles[wiggle_num]==NULL || wiggles[wiggle_num]->data==NULL || 
	(interval=wiggles[wiggle_num]->data[left_order_num])==NULL)
      send(CRASH);
    
    /* Look for a local maxima over the LOD threshold... */
    i=0; j=0; 
    contig=FALSE; first=last=NULL; 
    n_intervals= wiggles[wiggle_num]->num_wiggled_intervals; 

    do {
	if (!contig) { 
	    local_maxima= -10.0; local_minima= 1000.0; lod= 1000.0; 
	    peak_i=i; peak_j=j; still_falling=FALSE; off_end=FALSE;
	    starting_value = interval[i]->point[j]->lod_score;
	}
	prev_lod=lod;
	lod=interval[i]->point[j]->lod_score;
	if (!still_falling && lod>local_maxima && local_maxima>threshold &&
	    (local_maxima>=local_minima+min_peak_delta || 
	     local_minima >= starting_value + confidence_falloff)) {
	    /* Found a local maxima - now get peak and confidence interval. */
	    peak=peak_finder(&peak_i,&peak_j,&off_end,get_peak_maps,
			     interval,n_intervals,threshold,
			     confidence_falloff);
	    /* Store the data for this peak */
	    if (first==NULL) first=last=peak;
	      else { last->next=peak; last=peak; }
	    /* peak_finder sets the counters to be one step beyond peak. */
	    i=peak_i; j=peak_j; local_maxima= -10.0; 
	    still_falling=TRUE; lod= 1000.0; starting_value = 2000.0;
	    if(peak_i >= n_intervals) local_minima = 1000.0;
	    else local_minima = interval[i]->point[j]->lod_score;
	} else { /* Not a local maxima */
	    if (still_falling && lod>prev_lod) still_falling=FALSE;
	    if (!still_falling && lod>local_maxima) 
	      { local_maxima=lod; peak_i=i; peak_j=j; }
	    if (lod<local_minima) { local_minima=lod; }
	}
    } while (step(&i,&j,&contig,interval,n_intervals,TRUE));

    /* Deal with with fencepost error - have a max if lod is rising at end! */
    if (!off_end && lod>threshold && lod==local_maxima) {
	peak=peak_finder(&peak_i,&peak_j,&off_end,get_peak_maps,
			 interval,n_intervals,threshold,confidence_falloff);
	if (first==NULL) first=last=peak;
	  else { last->next=peak; last=peak; }
    }
      
    return(first);
}			


WIGGLE_PEAK *peak_finder(interval,point,off_end,get_peak_maps,wig,n_intervals,
			 threshold,falloff)
int *interval, *point;
bool *off_end;
bool get_peak_maps; /* if TRUE, fill in peak->map MAY REQUIRE CALCULATION! */
WIGGLE_INTERVAL **wig;
int n_intervals;
real threshold, falloff;
{
    /* i's are intervals and j's are points */
    int peak_i, peak_j, left_i, left_j, right_i, right_j, i, j, k;
    real lod, maxima; 
    bool off_right, off_left, contig;
    WIGGLE_PEAK *peak=NULL;

    run {
	*off_end=TRUE;

	/* First go to right boundary, finding the peak in the process. The args
	   will be side-effected to be one step beyond the right boundary. */

	right_i=peak_i= *interval; right_j=peak_j= *point;    
	maxima= wig[*interval]->point[*point]->lod_score; 
	off_right=TRUE; *off_end=TRUE; 

	while (step(interval,point,&contig,wig,n_intervals,TRUE)) {
	    if (!contig) { *off_end=FALSE; break; }
	    lod= wig[*interval]->point[*point]->lod_score;
	    if (lod>maxima) { maxima=lod; peak_i= *interval; peak_j= *point; }
	    else if ((falloff<0.0 && lod<maxima+falloff) ||
		     (falloff>0.0 && lod<falloff)) 
	      { off_right=FALSE; *off_end=FALSE; break; }
	    else { right_i= *interval; right_j= *point; }
	}
	    
	/* Now go to left boundary... */

	i=left_i=peak_i; j=left_j=peak_j; off_left=TRUE;    
	while (step(&i,&j,&contig,wig,n_intervals,FALSE)) {
	    if (!contig) { break; }
	    lod= wig[i]->point[j]->lod_score;
	    if ((falloff<0.0 && lod<maxima+falloff) ||
		(falloff>0.0 && lod<falloff)) { off_left=FALSE; break; }
	    else { left_i=i; left_j=j; }
	}

	/* And finally pack these data away... */

	single(peak,WIGGLE_PEAK);
	peak->next=NULL; peak->map=NULL; 
	k= wig[peak_i]->map->num_intervals-1;

	peak->left= wig[peak_i]->map->left[k];
	peak->right= wig[peak_i]->map->right[k];
	peak->qtl_pos= wig[peak_i]->point[peak_j]->qtl_pos;
	peak->lod_score= wig[peak_i]->point[peak_j]->lod_score;
	peak->var_explained= wig[peak_i]->point[peak_j]->var_explained;
	peak->qtl_weight= wig[peak_i]->point[peak_j]->qtl_weight;
	peak->qtl_dominance= wig[peak_i]->point[peak_j]->qtl_dominance;
	
	peak->forward_left= wig[right_i]->map->left[k];
	peak->forward_right= wig[right_i]->map->right[k];
	if (off_right) peak->forward_pos=OFF_END;
	else peak->forward_pos= wig[right_i]->point[right_j]->qtl_pos;

	peak->backward_left= wig[left_i]->map->left[k];
	peak->backward_right= wig[left_i]->map->right[k];
	if (off_left) peak->backward_pos=OFF_END;
	else peak->backward_pos= wig[left_i]->point[left_j]->qtl_pos;

	if (get_peak_maps) {
	    /* We assume here that best_point is pretty much near the 
	       ML_qtl_pos, in that we ignore max_like_map. */
	    if (wig[peak_i]->best_point==peak_j) {
		peak->map=wig[peak_i]->map;
		peak->free_map=FALSE;
	    } else { /* Shit! no map available for this pos, so calculate it */
		peak->map= alloc_qtl_map(wig[peak_i]->map->num_intervals,
				 wig[peak_i]->map->num_continuous_vars);
		peak->free_map=TRUE;
		mapcpy(peak->map,wig[peak_i]->map);
		peak->map->fix_pos[k]=peak->qtl_pos;
		temp_print("Calculating QTL map...");
		make_qtl_map(peak->map);
		temp_print(NULL);
	    }
	}

    } when_aborting { free_wiggle_peaks(peak); relay; }
    return(peak);
}
			     

void free_wiggle_peaks(p)
WIGGLE_PEAK *p;
{
    if (p==NULL) return;
    free_wiggle_peaks(p->next);
    if (p->free_map && p->map!=NULL) free_qtl_map(map);
    unsingle(p,WIGGLE_PEAK);
}


bool isa_test_wiggle(wiggle_num)
int wiggle_num;
{
    QTL_SEQUENCE *p;
    
    if (raw.data_type==BACKCROSS) return(FALSE);
    p=wiggles[wiggle_num]->seq; while (p->next!=NULL) p=p->next;
    if (p->genetics.interx_type!=TEST_MODELS) return(FALSE);
    return(TRUE);
}
	
    
#define NO_WIGGLES \
"No scan results have yet been saved.\nUse the 'scan' command first."
#define WIGGLE_NUM \
"%d is not a valid number for a saved scan result.\nUse a number from %d-%d.\n"
#define ORDER_ONE \
"%d.%d is not a valid number for a saved scan result.\nScan %d only has one entry, numbered %d.1.\n"
#define ORDER_NUM \
"%d.%d is not a valid number for a saved scan result.\nUse a number from %d.1-%d.%d.\n"

void get_wiggle_nums(str,wiggle,order)
char *str;
int *wiggle, *order;
{
    bool order_set;
    int i;

    if (num_wiggles==0) error(NO_WIGGLES);
    if (nullstr(str)) { *wiggle= *order= -1; return; }

    if ((i=strfinder(str,'.'))>0) { 
	str[i++]='\0'; if (sscanf(str,"%d",wiggle)<1) usage_error(num_args);
	order_set=TRUE; if (sscanf(&str[i],"%d",order)<1) 
	  usage_error(num_args);
	if ((--*wiggle)<0 || (--*order)<0) usage_error(num_args);
    } else { 
	*order= -1; if (sscanf(str,"%d",wiggle)<1) usage_error(num_args); 
	order_set=FALSE; if ((--*wiggle)<0) usage_error(num_args);
    }
	
    if (*wiggle<first_wiggle || *wiggle>=num_wiggles) {
	sf(ps,WIGGLE_NUM,*wiggle+1,first_wiggle+1,num_wiggles); 
	error(ps); 
    }
    if (wiggles[*wiggle]==NULL) send(CRASH);
	
    if (order_set) 
	if (*order>=wiggles[*wiggle]->num_orders) {
	    if (wiggles[*wiggle]->num_orders==1) 
	      sf(ps,ORDER_ONE,*wiggle+1,*order+1,*wiggle+1,*wiggle+1);
	    else sf(ps,ORDER_NUM,*wiggle+1,*order+1,*wiggle+1,*wiggle+1,
	       wiggles[*wiggle]->num_orders); 
	    error(ps); 
	}
}


#define NO_COMPARES \
"No compare results have yet been saved.\nUse the 'compare' command first."
#define COMPARE_NUM \
"%d is not a valid number for a saved compare result.\nUse a number from %d-%d.\n"
#define ORDER_ONE_COMP \
"%d.%d is not a valid number for a saved compare result.\nCompare %d only has one contiguous entry, numbered %d.1.\n"
#define ORDER_NUM_COMP \
"%d.%d is not a valid number for a saved compare result.\nUse a number from %d.1-%d.%d.\n"


void get_compare_nums(str,compare,contig)
char *str;
int *compare,*contig;
{
    int i;
    bool contig_set;

    if(num_compares == 0) error(NO_COMPARES);
    if (nullstr(str)) { *compare= *contig= -1; return; }
    
    if((i=strfinder(str,'.')) > 0) {
	str[i++] = '\0'; if(sscanf(str,"%d",compare) < 1) usage_error(num_args);
	contig_set = TRUE; 
	if(sscanf(&str[i],"%d",contig)<1) usage_error(num_args);
	if ((--*compare)<0 || (--*contig)<0) usage_error(num_args);
    } else {
	*contig = -1; if(sscanf(str,"%d",compare)<1) usage_error(num_args);
	contig_set = FALSE; if ((--*compare)<0) usage_error(num_args);
    }

    if(*compare == first_compare || *compare >= num_compares) {
	sf(ps,COMPARE_NUM,*compare+1,first_compare+1,num_compares);
	error(ps);
    }
    if(compares[*compare]==NULL) send(CRASH);

    if(contig_set)
      if(*contig>=compares[*compare]->num_contigs) {
	  if(compares[*compare]->num_contigs==1)
	    sf(ps,ORDER_ONE_COMP,*compare+1,*contig+1,*compare+1,*compare+1);
	  else
	    sf(ps,ORDER_NUM_COMP,*compare+1,*contig+1,*compare+1,*compare+1,
	       compares[*compare]->num_contigs);
	  error(ps);
      }
}


void save_wiggle(fp,n)
FILE *fp;
int n;
{
    int i,j;
    if (wiggles[n]->data == NULL)
      fprintf(fp,"%d\n",-1);
    else {
	fprintf(fp,"%d %s\n",wiggles[n]->trait,wiggles[n]->seq_string);
	
	for(i = 0; i < wiggles[n]->num_orders; i++) {
	    for(j = 0; j < wiggles[n]->num_wiggled_intervals; j++) {
		if(wiggles[n]->data[i][j] != NULL) {
		    save_interval(fp, wiggles[n]->data[i][j]);
		}
	    }
	}
    }
}


void save_interval(fp,interval)
FILE *fp;
WIGGLE_INTERVAL *interval;
{
    int i;

    
    sprintf(ps,"%d %.3lf %d %d %d %d\n",interval->num_points,
	    interval->cm_increment, interval->contig, interval->max_like_map,
	    interval->in_qtls_struct, interval->best_point);
    fpr(fp);

    /* save the peak map */
    save_qtl_map(fp,interval->map);

    /* save other data points (if any) */
    for(i = 0; i < interval->num_points; i++) {
	sf(ps,"%.4lf %.3lf %.3lf %.3lf %.4lf\n",interval->point[i]->qtl_pos,
	   interval->point[i]->qtl_weight, interval->point[i]->qtl_dominance,
	   interval->point[i]->lod_score, interval->point[i]->var_explained);
	fpr(fp);
    }
}


void save_qtl_map(fp,map)
FILE *fp;
QTL_MAP *map;
{
    int i;

    sf(ps,"%d %d %d %d %d\n",map->trait,map->num_intervals,map->max_intervals,
       map->num_continuous_vars,map->max_continuous_vars);
    fpr(fp);

    for(i = 0; i < map->num_intervals; i++) {
	sf(ps,"%d %d %.3lf %.3lf %.3lf %.3lf %.3lf\n", map->left[i],
	   map->right[i], map->interval_len[i], map->qtl_pos[i],
	   map->qtl_weight[i], map->qtl_dominance[i], map->fix_pos[i]);
	fpr(fp);
	sf(ps,"%.4lf %d %.4lf %.4lf %.4lf\n", map->constraint[i].backx_weight,
	   map->constraint[i].interx_type, map->constraint[i].a,
	   map->constraint[i].b, map->constraint[i].c);
	fpr(fp);
    }

    for(i = 0; i < map->num_continuous_vars; i++) {
	sf(ps,"%d %lf %lf\n",map->cont_var[i],map->cont_var_weight[i],
	   map->fix_cont_var_weight[i]);
	fpr(fp);
    }


    sf(ps,"%.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf\n",
       map->mu, map->sigma_sq, map->var_explained, map->chi_sq,
       map->null_mu, map->null_sigma_sq, map->null_log_like,
       map->log_like, map->no_data_like, map->abs_log_like);
    fpr(fp);
}


void load_wiggle(fp)
FILE *fp;
{
    int trait,i,j,n;
    int n_orders,n_ints,n_wiggled,foo;
    QTL_SEQUENCE *seqnce;

    fgetdataln(fp,NULL);
    if(!itoken(&ln,iREQUIRED,&trait))  return;
    if (trait != -1) {
	_filter(ln);
	seqnce = compile_intervals(ln);
	reset_state(seqnce,TRUE,&n_ints,&foo,&n_orders,&n_wiggled);
    }
    n = allocate_wiggle_struct(trait,seqnce,ln,n_ints,n_orders,n_wiggled);

    if(n == -1) return;
    if (trait != -1) {
	for(i = 0; i < wiggles[n]->num_orders; i++) {
	    for(j = 0; j < wiggles[n]->num_wiggled_intervals; j++) {
		if(wiggles[n]->data[i][j] != NULL) {
		    load_interval(fp, wiggles[n]->data[i][j]);
		}
	    }
	}
    }
}


void load_interval(fp,interval)
FILE *fp;
WIGGLE_INTERVAL *interval;
{
    int num_points,i,matched,i1,i2,i3,i4;
    real pos,weight,dom,lod,var_exp,inc;

    run {
	if(fscanf(fp,"%d %lf %d %d %d %d\n",&num_points,&inc,&i1,&i2,&i3,&i4) 
	   != 6)
	  error("error in loading saved qtls file");

	interval->num_points = num_points;  interval->cm_increment = inc;
	interval->contig = i1;  interval->max_like_map = i2;
	interval->in_qtls_struct = i3;  interval->best_point = i4;

	/* load map */
	interval->map = load_qtl_map(fp);

	/* load other data points */
	if(interval->num_points > 0)
	    parray(interval->point, interval->num_points, WIGGLE_POINT);

	for(i = 0; i < interval->num_points; i++) {
	    matched = fscanf(fp,"%lf %lf %lf %lf %lf\n",
	                     &pos,&weight,&dom,&lod,&var_exp);
			     
	    if(matched != 5)
	      error("error in loading saved qtls file");

	    interval->point[i]->qtl_pos = pos;  
	    interval->point[i]->qtl_weight = weight;
	    interval->point[i]->qtl_dominance = dom;
	    interval->point[i]->lod_score = lod;  
	    interval->point[i]->var_explained = var_exp;
	}
    } when_aborting {
	if(msg == ENDOFILE) { 
	    print("Error while loading saved qtls file...data not loaded.\n");
	}
	else relay;
    }	
}


QTL_MAP *load_qtl_map(fp)
FILE *fp;
{
    int i1,i2,i3,i4,i5, i;
    real r1,r2,r3,r4,r5,r6,r7,r8,r9,r10;
    QTL_MAP *map;

    if(fscanf(fp,"%d %d %d %d %d\n",&i1,&i2,&i3,&i4,&i5) != 5)
      error("error in loading qtls file");

    map = alloc_qtl_map(i3,i5);
    
    map->trait = i1;  map->num_intervals = i2;
    map->max_intervals = i3;
    map->num_continuous_vars = i4;  map->max_continuous_vars = i5;

    for(i = 0; i < map->num_intervals; i++) {
        if(fscanf(fp,"%d %d %lf %lf %lf %lf %lf\n",&i1,&i2,&r1,&r2,&r3,&r4,&r5)
	   != 7)
	  error("error in loading qtls file");

	map->left[i] = i1;  map->right[i] = i2;
	map->interval_len[i] = r1;  map->qtl_pos[i] = r2;
	map->qtl_weight[i] = r3;  map->qtl_dominance[i] = r4;
	map->fix_pos[i] = r5;

	if(fscanf(fp,"%lf %d %lf %lf %lf\n",&r1,&i1,&r2,&r3,&r4) != 5)
	  error("error in loading qtls file");
    
	map->constraint[i].backx_weight = r1;
	map->constraint[i].interx_type = i1;  map->constraint[i].a = r2;
	map->constraint[i].b = r3;  map->constraint[i].c = r4;
    }

    for(i = 0; i < map->num_continuous_vars; i++) {
	if(fscanf(fp,"%d %lf %lf\n",&i1,&r1,&r2) != 3)
	  error("error loading qtls file");
	
	map->cont_var[i] = i1;  map->cont_var_weight[i] = r1;
	map->fix_cont_var_weight[i] = r2;
    }

    if(fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	      &r1,&r2,&r3,&r4,&r5,&r6,&r7,&r8,&r9,&r10) != 10)
      error("error in loading qtls file");

    map->mu = r1;  map->sigma_sq = r2;
    map->var_explained = r3;  map->chi_sq = r4;
    map->null_mu = r5;  map->null_sigma_sq = r6; 
    map->null_log_like = r7;  map->log_like = r8;
    map->no_data_like = r9;  map->abs_log_like = r10;

    return(map);
}


void save_compare(fp,n)
FILE *fp;
int n;
{
    int i;

    if(compares[n] == NULL) return;
    if (compares[n]->data == NULL) {
	fprintf(fp,"%d\n", -1);
	return;
      }
    fprintf(fp,"%d %s\n",compares[n]->trait,compares[n]->seq_string);

    for(i = 0; i < compares[n]->num_orders; i++) {
	if(compares[n]->data[i] != NULL) {
	    fprintf(fp,"%d\n",compares[n]->data[i]->contig);
	    save_qtl_map(fp,compares[n]->data[i]->map);
	}
    }
}


void load_compare(fp)
FILE *fp;
{
    int trait,i,contig,n;
    int n_orders,n_ints,n_wiggled,foo;
    QTL_SEQUENCE *seqnce;

    fgetdataln(fp,NULL);
    if(!itoken(&ln,iREQUIRED,&trait))  return;
    if(trait != -1) {
    	_filter(ln);
	seqnce = compile_intervals(ln);
	reset_state(seqnce,FALSE,&n_ints,&foo,&n_orders,&n_wiggled);
    }
    n = allocate_compare_struct(trait,seqnce,ln,n_ints,n_orders);
    if(n == -1) return;
    if(compares[n]->num_contigs != -1) {
	compares[n]->num_contigs = 1;
	for(i = 0; i < compares[n]->num_orders; i++) {
	    fscanf(fp,"%d\n",&contig);
	    compares[n]->data[i]->contig = contig;
	    if(!contig) compares[n]->num_contigs += 1;
	    compares[n]->data[i]->map = load_qtl_map(fp);
	}
    }
}




