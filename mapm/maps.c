/******************************************************************************

 #    #    ##    #####    ####            ####
 ##  ##   #  #   #    #  #               #    #
 # ## #  #    #  #    #   ####           #
 #    #  ######  #####        #   ###    #
 #    #  #    #  #       #    #   ###    #    #
 #    #  #    #  #        ####    ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

//#define INC_LIB
//#define INC_SHELL
#include "mapm.h"
//#include "toplevel.h"
//#include "lowlevel.h"

/****************** Support for the MAP and MAP_LIST structs *****************/

/* internal only */
static void allocate_error_matrix(MAP *map);
static int insert_sorted_map(SAVED_LIST *list);
static void sort_last(SAVED_LIST *list);
static int insert_unsorted_map(SAVED_LIST *list);

/******************************* Map Functions *******************************/
static real apportion(bool rec_flag, real both, real first, real second);
static real poisson_add(real first, real second);
static real poisson_r_d(real rec_frac);
static real poisson_d_r(real dist);
static real kosa_add(real first, real second);
static real kosa_r_d(real rec_frac);
static real kosa_d_r(real dist);
static real poisson_d_r_deriv(real dist);
static real kosa_d_r_deriv(real dist);

/* external */
MAP_FUNCTION *mapfunction;
MAP_FUNCTION maps[2];
int num_map_functions;


MAP *
allocate_map (int maxloci)
{
    MAP *map;
    run {
        map= NULL;
	single(map, MAP);
	map->locus= NULL;
	map->rec_frac= NULL;
	map->fix_interval= NULL;
        map->num_loci= 0;
        map->max_loci= maxloci;
        map->sex_specific= FALSE;
	map->unlink= NONE_UNLINKED;
        map->log_like= 0.0;
	array(map->locus, map->max_loci, int);
	matrix(map->rec_frac, map->max_loci - 1, 2, real);
	array(map->fix_interval, map->max_loci - 1, int);
	array(map->map_name, NAME_LEN+1, char);
	map->allow_errors= FALSE;
        map->error_rate= NULL; 
	map->error_lod= NULL;
	clean_map(map);
    } when_aborting {
	free_map(map);
	relay;
    }
    return(map);
}


void 
allocate_error_matrix (MAP *map)
{
    if (map->error_rate!=NULL) return;
    run {
	array(map->error_rate,map->max_loci,real);
	matrix(map->error_lod,map->max_loci,raw.data.f2.num_indivs,real);
    } when_aborting {
	unarray(map->error_rate,real);
	unmatrix(map->error_lod,map->max_loci,real);
	relay;
    }
}


void 
free_map (MAP *map)
{
    if(map == NULL) return;
    unarray(map->fix_interval, int);
    unmatrix(map->rec_frac, map->max_loci - 1, real);
    unarray(map->locus, int);
    if (map->error_rate != NULL) {
        unarray(map->error_rate, real);
        unmatrix(map->error_lod, map->max_loci, real);
    }
    map->num_loci= 0;
    map->max_loci= 0;
    map->sex_specific= FALSE;
    map->allow_errors= FALSE;
    map->log_like= 0.0;
    unsingle(map, MAP);
}


bool 
clean_map (MAP *map)
{
    int i, j;

    map->num_loci= 0;
    map->unlink= NONE_UNLINKED; /* now unused */
    map->sex_specific= FALSE;
    map->allow_errors= FALSE;
    map->log_like= 0.0;
    for (i=0; i<map->max_loci-1; i++) {
	map->fix_interval[i]= FALSE;
	map->rec_frac[i][MALE]= map->rec_frac[i][FEMALE]= NOT_FIXED;
    }
    if (map->error_rate!=NULL) {
        for (i=0; i<map->max_loci; i++) {
	    map->error_rate[i]= 0.0;
	    for (j=0; j<raw.data.f2.num_indivs; j++) map->error_lod[i][j]=0.0;
	}
    }
    return(TRUE);
}


void 
init_for_ctm ( /* what else */
    MAP *map,
    bool sex,
    bool errors,
    bool start /* deal with start */
)
{ 
    int i, j;

    map->sex_specific= sex;    if (sex && raw.data_type!=CEPH)  send(CRASH);
    map->allow_errors= errors; if (errors && raw.data_type!=F2) send(CRASH);
    map->log_like= 0.0;

    for (i=0; i<map->num_loci-1; i++) {
	if (map->unlink==i+1) { /* now unused */
	    map->fix_interval[i]= TRUE;
	    map->rec_frac[i][MALE]= 0.50;
	    map->rec_frac[i][FEMALE]= 0.50;
	} 
	if (map->rec_frac[i][MALE]==UNLINK_ME) { 
	    map->fix_interval[i]= TRUE;
	    map->rec_frac[i][MALE]= 0.50;
	    map->rec_frac[i][FEMALE]= 0.50;
	} else if (map->rec_frac[i][MALE]!=NOT_FIXED) {
	    map->fix_interval[i]= TRUE;
	} else {
	    map->fix_interval[i]= FALSE;
	    map->rec_frac[i][MALE]= map->rec_frac[i][FEMALE]= startrecombs;
	}
    }
    if (map->allow_errors) {
	if (map->error_rate==NULL) allocate_error_matrix(map);
	for (i=1; i<map->num_loci-1; i++) {
	    map->error_rate[i]= error_rate[map->locus[i]];
	    for (j=0; j<raw.data.f2.num_indivs; j++) map->error_lod[i][j]=0.0;
	}
	map->error_rate[0]= map->error_rate[map->num_loci-1]= 0.0;
    }
}


void 
init_rec_fracs (MAP *map)
{ init_for_ctm(map,sex_specific,use_error_rate,TRUE); }


void 
init_not_fixed (MAP *map)
{ 
    int i;
    for (i=0; i<map->num_loci-1; i++)
      map->rec_frac[i][FEMALE]= map->rec_frac[i][MALE]= NOT_FIXED;
    init_for_ctm(map,sex_specific,use_error_rate,TRUE);
}


void 
mapcpy (MAP *to, MAP *from, bool clean_it)
{
    int i, j;
    if (to->max_loci<from->num_loci) send(CRASH);

    to->num_loci= from->num_loci;
    to->sex_specific= from->sex_specific;
    to->unlink= from->unlink; /* now unused */
    to->allow_errors= from->allow_errors;
    to->log_like= from->log_like;
    for (i=0; i<from->num_loci; i++) {
	to->locus[i]= from->locus[i];
    }
    for (i=0; i<from->num_loci-1; i++) {
	if (!clean_it) {
	    to->fix_interval[i]= from->fix_interval[i];
	    to->rec_frac[i][MALE]=   from->rec_frac[i][MALE];
	    to->rec_frac[i][FEMALE]= from->rec_frac[i][FEMALE];
	} else if (from->fix_interval[i]) {
	    to->fix_interval[i]= FALSE;
	    to->rec_frac[i][MALE]=   from->rec_frac[i][MALE]; /* fixes rfs */
	    to->rec_frac[i][FEMALE]= from->rec_frac[i][FEMALE];
	} else {
	    to->fix_interval[i]= FALSE;
	    to->rec_frac[i][MALE]=   NOT_FIXED;
	    to->rec_frac[i][FEMALE]= NOT_FIXED;
	}
    }
    if (from->allow_errors) {
	if (to->error_rate==NULL) allocate_error_matrix(to);
        for (i=0; i<from->num_loci; i++) {
	    to->error_rate[i]= from->error_rate[i];
	    for (j=0; j<raw.data.f2.num_indivs; j++) 
	      to->error_lod[i][j]= from->error_lod[i][j];
	}
    }
    if (!nullstr(from->map_name)) strcpy(to->map_name,from->map_name);
      else strcpy(to->map_name,"");
}


int 
insert_locus (    /* Returns TRUE if successful */
    MAP *map,
    int position,
    int locus
)
{
    int i;

    if (map->num_loci==map->max_loci) return(FALSE);

    /* this needs testing */
    if ((position>0 && position<map->num_loci) &&
	(map->rec_frac[position-1][MALE]!=NOT_FIXED ||
	 map->rec_frac[position-1][FEMALE]!=NOT_FIXED))
      error("attempt to insert a locus into a interval with fixed distance");

    if (map->unlink==position) /* not used */
      error("attempt to insert a locus into an unlinked interval");

    map->num_loci+= 1;
    for (i=map->num_loci-1; i>position; i--) {
	map->locus[i]= map->locus[i-1];
	if (i==map->num_loci-1) continue;
	map->rec_frac[i][MALE]=   map->rec_frac[i-1][MALE];
	map->rec_frac[i][FEMALE]= map->rec_frac[i-1][FEMALE];
    }
    map->locus[position]= locus;
    return(TRUE);
}


SAVED_LIST *
allocate_map_list (int maxmaps, int maxloci, bool sortflag, MAP **map)
{
    SAVED_LIST *list;
    int i;

    run {
        list= NULL;
        single(list, SAVED_LIST);
        list->max_maps= maxmaps;
        list->max_loci= maxloci;
	list->map_list= NULL;
	list->extra_map= NULL;
        array(list->map_list, list->max_maps, MAP*);
        for(i= 0; i < list->max_maps; i++)
            list->map_list[i]= allocate_map(list->max_loci);
        list->extra_map= allocate_map(list->max_loci);
        list->num_maps= 0;
        list->threshold= VERY_UNLIKELY;
        list->sorted= sortflag;
	if (map!=NULL) *map= list->extra_map;
    } except_when(ABORT)
	    free_map_list(list);
    return(list);
}


MAP *
get_map_to_bash (SAVED_LIST *list)
{ 
    clean_map(list->extra_map); 
    return(list->extra_map); 
}


void 
free_map_list (SAVED_LIST *list)
{
    int i;
    if(list == NULL) return;
    free_map(list->extra_map);
    if(list->map_list != NULL) {
        for(i= list->max_maps-1; i >= 0; i--)
	    free_map(list->map_list[i]);
	unarray(list->map_list, MAP*);
	unsingle(list, SAVED_LIST);
    }
}


void 
clean_list (SAVED_LIST *list)
{
    int i;

    list->num_maps= 0;
    list->threshold= VERY_UNLIKELY;
    for (i=0; i<list->max_maps-1; i++)
      clean_map(list->map_list[i]);
    clean_map(list->extra_map);
}


MAP *
get_best_map (SAVED_LIST *list)
{
    int i;
    real best;
    MAP *map=NULL;

    best=VERY_UNLIKELY; 
    if (list->num_maps==0) send(CRASH);
    for (i=0; i<list->num_maps; i++) {
        if (list->map_list[i]->log_like>best) 
	  { map=list->map_list[i]; best=map->log_like; }
    }
    return(map);
}
       
int 
insert_map_into_list (SAVED_LIST *list, MAP **map)
{
    int val; 

    list->extra_map= *map;
    if(list->sorted)
        val= insert_sorted_map(list);
    else 
        val= insert_unsorted_map(list);
    *map= list->extra_map;
    return(val);
}

void 
overwrite_map_num (SAVED_LIST *list, MAP **map, int chrom)
{
    MAP *tempmap;
    
    tempmap= list->map_list[chrom];
    list->extra_map= *map;
    list->map_list[chrom]= list->extra_map;
    list->extra_map= tempmap;
    *map= list->extra_map;
    return;
}


int 
insert_sorted_map (SAVED_LIST *list)
{
    int nmaps;
    MAP *tempmap;

    nmaps= list->num_maps;
    if(nmaps == 0) {
        tempmap= list->map_list[0];
        list->map_list[0]= list->extra_map;
	list->extra_map= tempmap;
	list->num_maps= 1;
	return(TRUE);
    }
    else if(list->extra_map->log_like < list->threshold)
        return(FALSE); /* extra_map->log_like is unacceptable */
    else {
        if(nmaps < list->max_maps) {
	    tempmap= list->map_list[nmaps];
	    list->map_list[nmaps]= list->extra_map;
	    list->extra_map= tempmap;
	    list->num_maps += 1;
	    sort_last(list);
	    return(TRUE);
	}
	else {
	    tempmap= list->map_list[nmaps-1];
	    list->map_list[nmaps-1]= list->extra_map;
	    list->extra_map= tempmap;
	    sort_last(list);
	    list->threshold= list->map_list[nmaps-1]->log_like;
    /* when map_list is full, threshold is set to worst log_like in list */
	    return(TRUE);
	}
    }
}


void 
sort_last (   /* Sorts last entry in the list */
    SAVED_LIST *list
)
{
    int i;
    MAP *movable,*tempmap;

    movable= list->map_list[list->num_maps-1]; /* last in list */
    for(i= list->num_maps-2; i >= 0; i--) {
        if(movable->log_like < list->map_list[i]->log_like)
	    return; /* it's in the proper place */
	else {
	    tempmap= list->map_list[i];
	    list->map_list[i]= movable;
	    list->map_list[i+1]= tempmap; 
	}
    }
    return;
}


int 
insert_unsorted_map (SAVED_LIST *list)
{
    MAP *tempmap;

    if(list->num_maps == list->max_maps)
        return(FALSE); /* no vacancy */
    else {
        tempmap= list->map_list[list->num_maps];
        list->map_list[list->num_maps]= list->extra_map;
	list->extra_map= tempmap;
	list->num_maps += 1;
	return(TRUE);
    }
}
	    

/******************************* Map Functions *******************************/

//real apportion();
//real poisson_add(), poisson_d_r();
//real poisson_r_d(), kosa_r_d();
//real kosa_add(), kosa_d_r();
//real poisson_d_r_deriv(), kosa_d_r_deriv();


real 
apportion (
    bool rec_flag, /* assume REC==1 */
    real both,
    real first,
    real second
)
{
    if (rec_flag) return ((first - second + both) / (2.0 * both));
    else return ((first + second - both) / (2.0 * (1.0 - both)));
}

real 
poisson_add (real first, real second)
{
    return first * (1.0 - second) + (1.0 - first) * second;
}

real 
poisson_r_d (real rec_frac)
{
#ifdef DONT_DO_THIS    
    if (raw.data_type == F2) {
	if (raw.data.f2.cross_type == RI_SIB)
	    rec_frac= rec_frac/(4 - (6*rec_frac));
	else if (raw.data.f2.cross_type == RI_SELF)
	    rec_frac= rec_frac/(2 - (2*rec_frac));
	else if (raw.data.f2.cross_type == F3_SELF)
	    rec_frac= (1.0 - sqrt(1.0-(2*rec_frac))) / 2.0;
    }
#endif

    return ((real) (rec_frac < 0.49999 ?
     -0.5 * log ((real)(1.0 - 2.0 * rec_frac)) :
     9.999));
}

real 
poisson_d_r (real dist)
{
    real rec_frac;

    rec_frac= 0.5 * (1.0 - exp((real) (-2.0 * dist)));
    return(rec_frac);

#ifdef DONT_DO_THIS    
    if(raw.data_type == F2) {
	if(raw.data.f2.cross_type == RI_SIB)
	    rec_frac= (4*rec_frac)/(1 + (6*rec_frac));
	else if(raw.data.f2.cross_type == RI_SELF)
	    rec_frac= (2*rec_frac)/(1 + (2*rec_frac));
	else if(raw.data.f2.cross_type == F3_SELF)
	    rec_frac= 2*rec_frac*(1.0-rec_frac);
    }
#endif
}

real 
kosa_add (real first, real second)
{
    return (first + second) / (1.0 + 4.0 * first * second);
}

real 
kosa_r_d (real rec_frac)
{
#ifdef DONT_DO_THIS    
    if(raw.data_type == F2) {
	if(raw.data.f2.cross_type == RI_SIB)
	    rec_frac= rec_frac/(4 - (6*rec_frac));
	else if(raw.data.f2.cross_type == RI_SELF)
	    rec_frac= rec_frac/(2 - (2*rec_frac));
	else if(raw.data.f2.cross_type == F3_SELF) 
	    rec_frac= (1.0 - sqrt(1.0 - (2*rec_frac))) / 2.0;
    }
#endif

    return ((real) (rec_frac < .49999 ?
     0.25 * log ((real)((1.0 + 2.0 * rec_frac) / (1.0 - 2.0 * rec_frac))) :
     9.999));
}

real 
kosa_d_r (real dist)
{
    real rec_frac;

    rec_frac= 0.5 * tanh ((real)(2.0 * dist));
    return(rec_frac);

#ifdef DONT_DO_THIS
    if(raw.data_type == F2) {
	if(raw.data.f2.cross_type == RI_SIB)
	    rec_frac= (4*rec_frac)/(1 + (6*rec_frac));
	else if(raw.data.f2.cross_type == RI_SELF)
	    rec_frac= (2*rec_frac)/(1 + (2*rec_frac));
	else if(raw.data.f2.cross_type == F3_SELF)
	    rec_frac= 2*rec_frac*(1.0-rec_frac);
    }
#endif
}

real 
poisson_d_r_deriv (real dist)
{
    return ((real) (exp(-2.0*dist)));
}

real 
kosa_d_r_deriv (real dist)
{
    return ((real) 4.0/( (exp(2.0*dist) + exp(-2.0*dist)) * 
			    (exp(2.0*dist) + exp(-2.0*dist)) ) );
}


void 
map_func (int mapnum)
{ mapfunction= &maps[mapnum]; }


void 
map_init (void)
{
    strcpy(maps[HALDANE].name,"Haldane");
    maps[HALDANE].add= poisson_add;
    maps[HALDANE].apportion= apportion;
    maps[HALDANE].rec_to_dist= poisson_r_d;
    maps[HALDANE].dist_to_rec= poisson_d_r;
    maps[HALDANE].d_to_r_deriv= poisson_d_r_deriv;

    strcpy(maps[KOSAMBI].name,"Kosambi");
    maps[KOSAMBI].add= kosa_add;
    maps[KOSAMBI].apportion= apportion;
    maps[KOSAMBI].rec_to_dist= kosa_r_d;
    maps[KOSAMBI].dist_to_rec= kosa_d_r;
    maps[KOSAMBI].d_to_r_deriv= kosa_d_r_deriv;

    num_map_functions= 2;
    map_func(HALDANE);
}



/******************************* Save/Load *******************************/


void 
read_map (FILE *fp, MAP *map)
{
    int i, j, num_loci, num, unlink, sex, errors;
    real rnum, like;
    char name[NAME_LEN+2], str[TOKLEN+1];
 
    clean_map(map);
    
    fgetln(fp);
    if (!nstoken(&ln,sREQUIRED,name,NAME_LEN+1) || name[0]!='*' ||
	sscanf(ln,"%d %d %d %d %lf",&num_loci,&unlink,&sex,&errors,&like)!=5)
      baddata("expected *chrom-name # # # #");
    else if (num_loci>map->max_loci) /* num_loci may be 0 */
      baddata("num_loci too large for map");
    
    map->unlink=unlink;           map->sex_specific=sex;
    strcpy(map->map_name,name+1); map->allow_errors= errors;
    map->num_loci=num_loci;       map->log_like=like;

    fgetln(fp);
    for (i=0; i<num_loci; i++) {
	if (nullstr(ln)) fgetln(fp);
	if (!itoken(&ln,iREQUIRED,&num)) send(CRASH);
	map->locus[i]=num;
    }
    fgetln(fp);
    for (i=0; i<num_loci-1; i++) {
	if (nullstr(ln)) fgetln(fp);
	if (!rtoken(&ln,rREQUIRED,&rnum)) send(CRASH);
	map->rec_frac[i][MALE]=rnum;
    }

    if (raw.data_type==CEPH && sex) {
	fgetln(fp);
	for(i=0; i<num_loci-1; i++) {
	    if (nullstr(ln)) fgetln(fp);
	    if (!rtoken(&ln,rREQUIRED,&rnum)) send(CRASH);
	    map->rec_frac[i][FEMALE]=rnum;
	}
    } else for (i=0; i<num_loci-1; i++) map->rec_frac[i][FEMALE]=0.0;

    fgetln(fp);
    for(i=0; i<num_loci-1; i++) {
	if (nullstr(ln)) fgetln(fp);
	if (!itoken(&ln,iREQUIRED,&num)) send(CRASH);
	map->fix_interval[i]=num;
    }   
    if (raw.data_type==F2 && errors) {
	if (map->error_rate==NULL) {
	    array(map->error_rate,map->max_loci,real);
	    matrix(map->error_lod,map->max_loci,raw.data.f2.num_indivs,real);
	}
	for (i=0; i<map->num_loci; i++) {
	    fgetln(fp);
	    if (!stoken(&ln,sREQUIRED,str) || !streq(str,"errors")) 
	      send(CRASH);
	    if (!rtoken(&ln,rREQUIRED,&rnum)) send(CRASH);
	    map->error_rate[i]=rnum;
	    if (i==0 || i==map->num_loci-1) {  /* not ends */
		for(j=0; j<raw.data.f2.num_indivs; j++) 
		  map->error_lod[i][j]= 0.0;
	    } else {
		fgetln(fp);
		for (j=0; j<raw.data.f2.num_indivs; j++) {
		    if (nullstr(ln)) fgetln(fp);
		    if (!rtoken(&ln,rREQUIRED,&rnum)) send(CRASH);
		    map->error_lod[i][j]=rnum;
		}
	    }
	}
    } /* else do not touch errors */
}


void 
write_map (FILE *fp, MAP *map)
{
    int i, j;

    sprintf(ps, "*%s %d %d %d %d %.3lf\n",
            map->map_name, map->num_loci, map->unlink, map->sex_specific,
            map->allow_errors, map->log_like); fpr(fp);

    for (i=0; i<map->num_loci; i++) {
        if (i%10==0 && i!=0 && i!=map->num_loci-1) fnl(fp);
	sprintf(ps, "%d ", map->locus[i]); fpr(fp);
    }
    fnl(fp);

    for (i=0; i<map->num_loci-1; i++) {
        if (i%10==0 && i!=0) fnl(fp);
	sprintf(ps, "%6.4lf ", map->rec_frac[i][MALE]); fpr(fp);
    }
    fnl(fp);

    if (raw.data_type==CEPH && map->sex_specific) {
        for(i=0; i<map->num_loci-1; i++) {
	    if (i%10==0 && i!=0) fnl(fp);
	    sprintf(ps, "%6.4lf ", map->rec_frac[i][FEMALE]); fpr(fp);
	}
	fnl(fp);
    }

    for(i=0; i<map->num_loci-1; i++) {
        if (i%20==0 && i!=0 && i!=map->num_loci-1) fnl(fp);
	sprintf(ps, "%d ", map->fix_interval[i]); fpr(fp);
    }
    fnl(fp);

    if (raw.data_type==F2 && map->allow_errors) {
        for (i=0; i<map->num_loci; i++) {
	    sprintf(ps, "errors %.7lf\n", map->error_rate[i]); fpr(fp);
	    if (i==0 || i==map->num_loci-1) continue;  /* not ends */
	    for (j=0; j<raw.data.f2.num_indivs; j++) {
	        if (j%10==0 && j!=0) fnl(fp);
		sprintf(ps, "%.3lf ", map->error_lod[i][j]);
		fpr(fp);
	    }
	    fnl(fp);
	}
    }
}
