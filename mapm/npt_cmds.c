/******************************************************************************

 #    #  #####    #####           ####   #    #  #####    ####            ####
 ##   #  #    #     #            #    #  ##  ##  #    #  #               #    #
 # #  #  #    #     #            #       # ## #  #    #   ####           #
 #  # #  #####      #            #       #    #  #    #       #   ###    #
 #   ##  #          #            #    #  #    #  #    #  #    #   ###    #    #
 #    #  #          #   #######   ####   #    #  #####    ####    ###     ####

******************************************************************************
   This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

//#define INC_LIB
//#define INC_SHELL
//#define INC_MISC
#include "mapm.h"
//#include "toplevel.h"
//#include "lowlevel.h"

/* internal procedures and variables */
static bool try_marker (
        int *marker,
        MAP *original_map, /* contains the list of loci in the order */
        SAVED_LIST *list,
        bool *excluded,
        bool *zero,
        int *count
);

bool find_seed_order(bool is_subset, int *locus, int num_loci, int size, int max_tries, real thresh, MAP *map, MAP *temp_map, bool **temp /* [num_loci][num_loci] */);

/* ERROR KLUDGE Stuff */
#define ESTEPS  200
#define ESTART  (-9.9)
#define ESTEP   0.1
int errs[ESTEPS], noterrs[ESTEPS], lod[ESTEPS];
int cum_not[ESTEPS], cum_err[ESTEPS];
real ratio[ESTEPS];


void 
npt_cmds_init (void)
{ 
    int i;
    for (i=0; i<ESTEPS; i++) errs[i]=noterrs[i]=0;
}


command 
make_map (void)
{
    MAP *map;
    int num_loci;

    map= NULL;
    mapm_ready(ANY_DATA,2,MAYBE_PERM,&num_loci);
    nomore_args(0);

    run {
        map=allocate_map(num_loci);
	print(MAP_DIVIDER);
	for_all_orders(seq, map) {
	    init_rec_fracs(map); 	
	    converge_to_map(map);
	    print_long_map(map,"Map:");
	    print(MAP_DIVIDER);
	}
    } on_exit { 
	free_map(map); 
	relay_messages; 
    }
}


command 
draw_map (void)
{
    MAP *map;
    int num_loci;
    char name[PATH_LENGTH+1];
    FILE *fp=NULL;

    /* ONE_ORDER only - really don't want this to try to create a 
       60 or more page PostScript file if someone inadvertantly executes
       it with a compare type sequence of many permutations */
 
    map= NULL;
    mapm_ready(ANY_DATA,2,ONE_ORDER,&num_loci);
    use_uncrunched_args();
    get_arg(stoken,"map",name);
    nomore_args(num_args);

    run {
	if (!make_filename(name,FORCE_EXTENSION,PS_EXT))
	  { sprintf(ps, "illegal filename '%s'", name); error(ps); }
	fp= open_file(name,WRITE);
	sprintf(ps, "Drawing map in PostScript file '%s'... \n", name); pr();
        map=allocate_map(num_loci);
	get_one_order(seq,map);
	init_rec_fracs(map); 	
	converge_to_map(map);
	print_ps_map(fp,map);
	close_file(fp);
	print("ok\n");

    } on_exit { 
        if (msg==CANTOPEN) { 
	    sprintf(ps, "Can't create output file '%s'", name); error(ps);
	} else if (msg==CANTCLOSE) {
	    error("\nCan't close file - disk is full?"); 
	} else {
	    free_map(map); 
	    relay_messages;
	}
    }
}


command 
genotypes (void)
{
    MAP *map;
    int num_loci;

    map= NULL;
    mapm_ready(F2_DATA,2,MAYBE_PERM,&num_loci);
    nomore_args(0);

    run {
        map=allocate_map(num_loci);
	print(MAP_DIVIDER); print("Genotypes:\n");
	for_all_orders(seq, map) {
	    init_rec_fracs(map); 	
	    converge_to_map(map);
	    print_f2_map_genotypes(map,use_haplotypes,TRUE,0, NULL);
        print(MAP_DIVIDER);
	}
    } on_exit {
	free_map(map); 
	relay_messages; 
    }
}


#define COMPARE_NONE \
  "no orders allowed by three-point analysis at threshold %.2lf\n"
#define BEST_THRESH "\nBest order%s at threshold %.2lf"
#define BEST_ORDERS "\nBest %d orders"
#define BEST_ORDER  "\nBest order"
#define N_EXCLUDED  " (%d excluded by three-point analysis)"

command 
compare (void)
{
    SAVED_LIST *list=NULL;
    MAP *map;
    int maps_to_save, num_to_print, loci_per_map, n, total, excluded, tried, i;
    int *locus, num_loci, prev;
    real threshold, best;

    mapm_ready(ANY_DATA,2,PERM_SEQ,&num_loci);
    get_arg(itoken,20,&maps_to_save);
    get_arg(rtoken,0.0,&threshold);
    if (maps_to_save<1 || maps_to_save>999) usage_error(1);
    if (threshold<0.0 || threshold>100.0) usage_error(2);
    if (threshold==0.0) threshold=VERY_UNLIKELY; else threshold=-threshold;
    nomore_args(num_args);

    run {
       list=allocate_map_list(maps_to_save,num_loci,SORTED,&map);
       array(locus,num_loci,int); /* NO CRUNCH?, #loci shouldn't change */
       get_list_of_all_loci(seq,locus,&loci_per_map,num_loci);
       num_orders=0; /* global - deletes all pre-existing orders */
       for (i=0; i<raw.num_markers; i++) order_next[i]=NO_LOCUS;
       if (use_three_pt)
	 setup_3pt_data(locus,loci_per_map,-three_pt_threshold);

       total=excluded=tried=0; n=0;
       for_all_orders(seq,map) total++;
       for_all_orders(seq,map) {
	   keep_user_amused("map",++n,total);
	   if (use_three_pt &&
	       !three_pt_verify(map->locus,map->num_loci,three_pt_window))
	     { excluded++; continue; }
	   init_rec_fracs(map);
	   converge_to_map(map);
	   insert_map_into_list(list,&map);
	   tried++;
       }
       if (tried==0)
	 { sprintf(ps, COMPARE_NONE, three_pt_threshold); pr(); abort_command(); }

       map=get_best_map(list); best=map->log_like; i=1;
       if (threshold!=VERY_UNLIKELY) {/* can you say "abstraction violation" */
	   while (i<list->num_maps && 
		  list->map_list[i]->log_like>(best+threshold)) i++;
	   num_to_print=i;
	   sprintf(ps, BEST_THRESH, maybe_s(num_to_print), -threshold); pr();
       } else {
	   num_to_print=min(maps_to_save,tried);
	   if (num_to_print>1) sprintf(ps, BEST_ORDERS, num_to_print);
	     else strcpy(ps,BEST_ORDER);
	   pr();
       }
       if (excluded>0) { sprintf(ps, N_EXCLUDED, excluded); pr(); }
       print(":\n");
       print_list(list,num_to_print);

       print("order1 is set");
       for (i=0; i<map->num_loci; i++) {
	   /* sprintf(ps,"%s ",rag(loc2str(map->locus[i]))); pr(); */
	   if (i==0) order_first[0]=map->locus[i]; 
	     else order_next[prev]= map->locus[i]; /* lint warning is OK */
	   prev=map->locus[i];
       }
       nl(); 
       unorder_first[0]=NO_LOCUS; order_next[prev]=NO_LOCUS; num_orders=1;
       
   } on_exit { 
       free_map_list(list); 
       relay_messages; 
   }
}


command 
ripple (void)
{
    MAP *map0=NULL, *map;
    SAVED_LIST *list;
    int n_loci, window, excluded, tried, i, j, k, n;
    real thresh, best;
    bool same;

    mapm_ready(ANY_DATA,3,ONE_ORDER,&n_loci);
    get_arg(itoken,5,&window);   if (!irange(&window,3,99)) usage_error(1);
    get_arg(rtoken,2.0,&thresh);if (!rrange(&thresh,0.0,100.0)) usage_error(1);
    if (window>n_loci) 
      error("window size is greater than number of loci in seq");
    if (has_fixed_dists(seq)) 
      error("can not use sequence with fixed distances");
    nomore_args(num_args);
    for (i=window, n=1; i>1; i--) n*=i; /* window factorial */

    run {
	map0=allocate_map(n_loci);
	list=allocate_map_list(n,n_loci,SORTED,&map);
	get_one_order(seq,map0);
	init_rec_fracs(map0);
	print(MAP_DIVIDER);
	converge_to_map(map0);
	print_long_map(map0,"Map To Test:");
	if (use_three_pt)
	  setup_3pt_data(map0->locus,map0->num_loci,-three_pt_threshold);

	print(MAP_DIVIDER);
	sprintf(ps, "Window-size: %d  Log-likelihood Threshold: %.2lf\n",
            window, thresh); pr();
	print("Comparing maps with ALL flanking markers...\n\n");
	for (i=0; i<=n_loci-window; i++) {
	    clean_list(list); map=get_map_to_bash(list);
	    make_compare_seq(map0->locus,map0->num_loci,i,window);

	    print("compare ");
	    if (i>0) print("...{"); else print("   {");
	    for (j=i; j<i+window; j++) {
		print(rag(loc2str(map0->locus[j]))); 
		if (j!=i+window-1) print(" ");
	    }
	    if (i+window-1<map0->num_loci-1) print("}... "); 
	      else print("}    ");
	    /* to_column(17+(print_names ? 10:5)*window); */
	    flush();

	    excluded=tried=0; 
	    for_all_orders(seq,map) {
		if (use_three_pt && /* only threept the compared part */
		    !three_pt_verify(map->locus+i,window,three_pt_window))
		  { excluded++; continue; }
		init_rec_fracs(map);
		converge_to_map(map);
		insert_map_into_list(list,&map); tried++;
	    }
	    if (tried==0) { print("(all excluded)\n"); continue; }

	    same=TRUE; map=get_best_map(list); 
	    for (j=0; j<map0->num_loci; j++)
	      if (map0->locus[j]!=map->locus[j]) { same=FALSE; break; }
	    best=map->log_like; k=1;
	    while (k<list->num_maps && /* thresh>0 here */
		   list->map_list[k]->log_like>(best-thresh)) k++;

	    if (!same || k>1) {
		sprintf(ps, "\nbest order%s:", maybe_s(k)); pr();
		if (excluded>0) { sprintf(ps, " (%d excluded)", excluded); pr(); }
		nl(); print_list(list,k); nl();	
	    } else print("ok\n");
	}
	print(MAP_DIVIDER);

    } on_exit {
	free_map_list(list); 
	free_map(map0);
	relay_messages;
    }
}


#define TRY_NONE "%s - all orders excluded by three-point analysis\n"
#define TRY_TOO_MANY "too many loci to try"
#define MAX_PAIRED 20

command 
try (void)
{
    SAVED_LIST **list=NULL;
    MAP *map=NULL;
    int *marker_to_try=NULL, **new_marker=NULL, num_to_try_at_once, num_tries;
    int num_seq_loci, num_intervals, first_marker, i, j, n, m, next;
    int num_total, max_paired;
    bool **exclude_interval=NULL, **zero_placement=NULL;
    char *err, token[TOKLEN+1];

    void expand_seq_names();      /* KLUDGE from sequence.c, gotta do better */
    char *markers, str[MAX_SEQ_LEN+99]; /* KLUDGE KLUDGE KLUDGE */

    mapm_ready(ANY_DATA,2,MAYBE_PERM,&num_seq_loci);
    strcpy(str,args); markers=str;
    if (nullstr(args)) usage_error(0);
    /* NEED SOMETHING LIKE PARSE_LOCUS_ARGS AND CRUNCH */

    run {
	array(marker_to_try,raw.num_markers,int); /* w/paired, num_total */
	matrix(new_marker,raw.num_markers,MAX_PAIRED,int);
	
	num_tries=0; num_total=0; max_paired=1;
	for (i=0; i<raw.num_markers; i++) for (j=0; j<MAX_PAIRED; j++)
	  new_marker[i][j]=NO_LOCUS;

	/* preprocess out the paired loci */
	expand_seq_names(markers);
	while (stoken(&markers,sREQUIRED,token)) {
	    if (streq(token,"[")) {
		if (!stoken(&markers,sREQUIRED,token))
		  error("expected a locus after '['");
		if (!is_a_locus(token,&n,&err) || !nullstr(err))
		  { sprintf(ps, "'%s' is not a locus", token); error(ps); }
		if (num_total==raw.num_markers) error(TRY_TOO_MANY);
		next=0;
		new_marker[num_tries][next++]=n;
		marker_to_try[num_total++]=n;
		while (stoken(&markers,sREQUIRED,token)) {
		    if (streq(token,"]")) break;
		    if (!is_a_locus(token,&m,&err) || !nullstr(err))
		      { sprintf(ps, "'%s' is not a locus", token); error(ps); }
		    if (num_total==raw.num_markers) error(TRY_TOO_MANY);
		    marker_to_try[num_total++]=m;
		    if (next==MAX_PAIRED) error(TRY_TOO_MANY);
		    new_marker[num_tries][next++]=m;
		}
		if (!streq(token,"]")) error("missing ']'");
		max_paired=max(next,max_paired);
		num_tries++;
		
	    } else if (streq(token,"<")) {
		if (!stoken(&markers,sREQUIRED,token))
		  error("expected a locus after '<'");
		if (!is_a_locus(token,&n,&err) || !nullstr(err))
		  { sprintf(ps, "'%s' is not a locus", token); error(ps); }
		if (num_total==raw.num_markers) error(TRY_TOO_MANY);
		next=0;	
		new_marker[num_tries][next++]=n;
		marker_to_try[num_total++]=n;
		while (stoken(&markers,sREQUIRED,token)) {
		    if (streq(token,">")) break;
		    if (!is_a_locus(token,&m,&err) || !nullstr(err))
		      { sprintf(ps, "'%s' is not a locus", token); error(ps); }
		    if (num_total==raw.num_markers) error(TRY_TOO_MANY);
		    marker_to_try[num_total++]=m;
		    if (next==MAX_PAIRED) error(TRY_TOO_MANY);
		    new_marker[num_tries][next++]=m;
		}
		if (!streq(token,">")) error("missing '>'");
		if (next>1) { /* then it's paired */
		    for (i=0; i<next; i++)
		      new_marker[num_tries+1][i]=
			new_marker[num_tries][next-i-1];
		    max_paired=max(next,max_paired);
		    num_tries+=2;
		} else num_tries++;
		
	    } else {
		if (!is_a_locus(token,&n,&err) || !nullstr(err))
		  { sprintf(ps, "'%s' is not a locus", token); error(ps); }
		if (num_total==raw.num_markers) error(TRY_TOO_MANY);
		new_marker[num_tries][0]=n;
		marker_to_try[num_total++]=n;
		num_tries++;
	    }
	}

	/* test some limits here?? */
	num_intervals= num_seq_loci+1;
	num_to_try_at_once= min((print_names ? 6:8),num_tries);
	matrix(exclude_interval,num_to_try_at_once,num_intervals,bool);
	matrix(zero_placement,num_to_try_at_once,num_intervals,bool);
	map= allocate_map(num_seq_loci);
	array(list,num_to_try_at_once,SAVED_LIST*); 
	for (i=0; i<num_to_try_at_once; i++) 
	  list[i]= allocate_map_list(num_intervals+1,num_seq_loci+max_paired,
				     UNSORTED,NULL);

	/* list all loci in seq and args in map->locus, then setup 3pt */
	get_one_order(seq,map);
	init_rec_fracs(map); /* no ctm */
	if (use_three_pt) {
	    for (j=0; j<map->num_loci; j++) {
		if (num_total==raw.num_markers) error(TRY_TOO_MANY);
		marker_to_try[num_total++]=map->locus[j];
	    }
	    setup_3pt_data(marker_to_try,num_total,-three_pt_threshold);
	}

	/* Here we loop over 8 markers to try at a time, calculating the
	   8 map_lists and printing the results. Each map_list has the 
	   results of sticking one locus into all selected intervals of
	   one order. Note that i counts 0..7 and j counts from n...n+7, 
	   where n is the 1st-marker of this set of 8 to try. */
	
	n=0;
	for_all_orders(seq,map) {
	    for (first_marker=0; first_marker<num_tries; 
		 first_marker+=num_to_try_at_once) {
		/* For each 8 markers to try... */ 
		for (i=0, j=first_marker; j<num_tries && i<num_to_try_at_once;
		     j++, i++) {
		    if (!try_marker(new_marker[j],map,list[i],
				    exclude_interval[i],zero_placement[i],&n))
		      { sprintf(ps, TRY_NONE, rag(loc2str(marker_to_try[j]))); pr(); }
		}
		nl();
		print_trys(list,map,exclude_interval,new_marker,i,
			   first_marker);
		nl(); n=0;
		for (i=0; i<num_to_try_at_once; i++) clean_list(list[i]);
	    }
        }

    } on_exit {
	unarray(marker_to_try,int);
	unarray(new_marker,int);
	unmatrix(exclude_interval,num_to_try_at_once,bool);
	unmatrix(zero_placement,num_to_try_at_once,bool);
	free_map(map);
	if (list!=NULL) {
	    for (i=0; i<num_to_try_at_once; i++) free_map_list(list[i]);
	    unarray(list,SAVED_LIST*);
	}
	relay_messages;
    }
}


/* really should replace this all with npt_exclusions, or some such thing */
#define is_zero(rec_frac,sex) \
  (rec_frac[MALE]<ZERO_DIST && (!sex || rec_frac[FEMALE]<ZERO_DIST))

bool 
try_marker (
    int *marker,
    MAP *original_map, /* contains the list of loci in the order */
    SAVED_LIST *list,
    bool *excluded,
    bool *zero,
    int *count
)
{
    int i, j, num_ok, last, sex, num_paired, best_i=0;
    MAP *map;
    real best;
    
    clean_list(list);
    map=get_map_to_bash(list);
    original_map->unlink= NONE_UNLINKED;

    for (i=0; i<original_map->num_loci; i++) excluded[i]=FALSE;
    num_ok= original_map->num_loci+1;
    for (num_paired=0; marker[num_paired]!=NO_LOCUS; num_paired++)
      if (use_three_pt)
	num_ok=three_pt_exclusions(original_map->locus,original_map->num_loci,
				   marker[num_paired],excluded);
    if (num_paired==0) send(CRASH);
    
    last=original_map->num_loci; best=VERY_UNLIKELY;
    for (i=0; i<=last; i++) if (!excluded[i]) {
	mapcpy(map,original_map,TRUE);
	for (j=0; j<num_paired; j++) 
	  if (!insert_locus(map,i+j,marker[j])) send(CRASH);
	keep_user_amused("map",++*count,0);
	init_rec_fracs(map);
	converge_to_map(map);
	if (map->log_like>best) { best=map->log_like; best_i=i; }
	
	/* find zero placements. indexing is:
	   orig-loci:     0   1     |     0   1     |     0   1
	   intervals:   0   1   2   |   0   1   2   |   0   1   2
	   new-seq:   x   0   1     |     0 x 1     |     0   1   x
	   recfracs:    0   1       |      0 1      |       0   1   
	   index:        i=0        |      i=1      |        i=2       */
	sex=map->sex_specific;
	if (i>0    && is_zero(map->rec_frac[i-1],sex)) zero[i]=TRUE;
	if (i<last && is_zero(map->rec_frac[i+num_paired-1],sex)) zero[i]=TRUE;
	if (zero[i]) for (j=1; j<num_paired; j++)
	  if (i+j<last && !is_zero(map->rec_frac[i+j],sex)) zero[i]=FALSE;
	insert_map_into_list(list,&map);
    } else keep_user_amused("map",++*count,0);

    if (num_ok>0 && zero[best_i]) {
	for (i=best_i+1; !excluded[i] && i<=last &&
	     list->map_list[i]->log_like-best>ZERO_PLACE && zero[i]; i++)
	  list->map_list[i]->log_like=ZERO_LIKE;
	for (i=best_i-1; !excluded[i] && i>=0 && 
	     list->map_list[i]->log_like-best>ZERO_PLACE && zero[i]; i--)
	  list->map_list[i]->log_like=ZERO_LIKE;
    }
    /* INF dist */
    mapcpy(map,original_map,TRUE); /* cleans map */
    i=map->num_loci-1; /* interval to unlink */
    for (j=0; j<num_paired; j++) 
      if (!insert_locus(map,map->num_loci+j,marker[j])) send(CRASH);
    keep_user_amused("map",++*count,0);
    map->rec_frac[i][MALE]=  UNLINK_ME;
    map->rec_frac[i][FEMALE]=UNLINK_ME;
    init_rec_fracs(map);
    converge_to_map(map);
    insert_map_into_list(list,&map);

    return(num_ok==0 ? FALSE:TRUE);
}



#define ORDER_NO_GROUPS \
 "\nNo linkage groups found at specified threshold.\n"
#define ORDER_AT "Placing at log-likelihood threshold %.2lf...\n"
#define START_CRITERIA \
 "Starting Orders: Size %d, Log-Likelihood %.2lf, Searching up to %d subsets\n"
#define INF_CRITERIA \
 "Informativeness: min #Individuals %d%s, min Distance %s\n"
#define LINK_CRITERIA \
 "Linkage Groups at min LOD %.2lf, max Distance %s\n"

command 
order_maker (void)
{
    int *loci=NULL, num_loci, *linkage_group=NULL, num_unlinked, group_size;
    int *subset=NULL, **seed_temp=NULL, subset_size, groups_done;
    int i, j, num, prev, num_unplaced, seed_size, seed_tries;
    bool seed_ok, found;
    MAP  *order=NULL, *seed_map=NULL;
    PLACEME **unplaced;
    real lodbound, thetabound, seed_like;

    /* args: lod theta start-window start-like max-to-try */
    mapm_ready(F2_DATA,3,LIST_SEQ,&num_loci); /* F2 because #indivs */
    if (!rtoken(&args,default_lod,&lodbound) || !rrange(&lodbound,0.0,1000.0)) 
      usage_error(1);
    if (!rtoken(&args,default_theta,&thetabound) || !input_dist(&thetabound))
      usage_error(2);
    get_arg(itoken,5,&seed_size);   if (seed_size<3) usage_error(0);
    get_arg(rtoken,3.0,&seed_like); if (seed_like<=0.0) usage_error(0);
    get_arg(itoken,50,&seed_tries); if (seed_tries<1) usage_error(0);
    nomore_args(0);
    seed_like= -seed_like;

    sprintf(ps, LINK_CRITERIA, lodbound, rag(rf2str(thetabound))); pr();
    sprintf(ps, START_CRITERIA, seed_size, -seed_like, seed_tries); pr();

    sprintf(ps, INF_CRITERIA, npt_min_indivs,
            (npt_codominant ? " (codominant only)":""),
            rag(rf2str(npt_min_theta))); pr();
    if (npt_threshold==npt_first_threshold) /* globals */
      sprintf(ps, "Placement Threshold %.2lf, Npt-Window %d\n", npt_threshold,
              npt_window);
    else 
      sprintf(ps, "Placement Threshold-1 %.2lf, Threshold-2 %.2lf, Npt-Window %d\n",
              npt_first_threshold, npt_threshold, npt_window);
    pr();
    
    run {
	alloc_list_of_all_loci(seq,&loci,&num_loci);
	array(linkage_group,num_loci,int);
	order=allocate_map(MAX_MAP_LOCI);
	seed_map=allocate_map(seed_size);
	array(subset,num_loci,int);
	matrix(seed_temp,seed_tries+1,num_loci,int);
	parray(unplaced,num_loci,PLACEME);
	for (i=0; i<num_loci; i++) {
	    array(unplaced[i]->excluded,MAX_MAP_LOCI+1,bool);
	    unplaced[i]->best_map=allocate_map(MAX_MAP_LOCI);
	}

	num_unlinked=num_loci;
	num_orders=0; /* global - deletes all pre-existing orders */
	for (i=0; i<raw.num_markers; i++) order_next[i]=NO_LOCUS;

	groups_done=0;
	do {
	    get_linkage_group(loci,&num_unlinked,linkage_group,&group_size,
			      lodbound,thetabound);
	    if (group_size>=2) {
		print(MAP_DIVIDER);
		inv_isort(linkage_group,group_size);
		groups_done++;
		sprintf(ps, "Linkage group %d, %d Markers:\n", groups_done,
                group_size); pr();
		for (i=0; i<group_size; i++) {
		    num=linkage_group[i];
		    sprintf(ps, "%4d %s%s", num + 1, raw.locus_name[num],
                    (use_haplotypes && haplotyped(num) ? "+":""));
		    pad_to_len(ps,14); pr();
		    if (i==group_size-1 || i%5==4) nl(); else print("  ");
		}
		nl();

		if (group_size<3) { 
		    print("Group is too small to map.\n"); 
		    sprintf(ps, "order%d= ", groups_done); pr();
		    for (i=0; i<group_size; i++) {
			sprintf(ps, "%s ", rag(loc2str(linkage_group[i]))); pr();
			if (i==0) order_first[groups_done-1]=linkage_group[i];
			else order_next[prev]=linkage_group[i]; 
			prev=linkage_group[i]; /* lint warning is OK */
		    }
		    nl();
		    num_orders++;
		    continue;
		}

		if (use_three_pt) 
		  setup_3pt_data(linkage_group,group_size,-three_pt_threshold);
		seed_ok=FALSE;
		informative_subset(linkage_group,group_size,
				   npt_min_indivs,npt_min_theta,npt_codominant,
				   use_haplotypes,subset,&subset_size);

		if (subset_size<3 || subset_size<seed_size) 
		  print("Most informative subset is too small... \n");
		else if (subset_size==group_size) 
		  print("All markers are informative... \n");
		else {
		    print("Most informative subset: ");
		    for (i=0; i<subset_size; i++) {
			sprintf(ps, "%d%s", (subset[i]) + 1,
                    (use_haplotypes && haplotyped(subset[i]) ? "+":""));
			pr(); if (i<subset_size-1) print(" ");
		    }
		    nl();
		    if (find_seed_order(TRUE,subset,subset_size,seed_size,
					seed_tries,seed_like,order,seed_map,
					seed_temp))
		      seed_ok=TRUE;
		      else if (group_size==subset_size) continue; /* punt LG */
		}

		if (!seed_ok) {
		    if (find_seed_order(FALSE,linkage_group,group_size,
					seed_size,seed_tries,seed_like,order,
					seed_map,seed_temp))
		      seed_ok=TRUE;
		      else continue; /* punt this LG */
		}

		/* delete seed markers from linkage group */
		num_unplaced=0;
		for (j=0; j<group_size; j++) {
		    for (i=0, found=FALSE; i<order->num_loci; i++)
		      if (order->locus[i]==linkage_group[j])
			{ found=TRUE; break; }
		    if (!found) 
		      unplaced[num_unplaced++]->locus=linkage_group[j];
		}

		if (num_unplaced>0) {
		    /* extend_order is chatty, prints map and places at end */
		    nl(); sprintf(ps, ORDER_AT, npt_first_threshold); pr();
		    extend_order(order,unplaced,&num_unplaced,
				 -npt_first_threshold,TRUE);

		    if (npt_first_threshold!=npt_threshold && num_unplaced>0) {
			print(SUB_DIVIDER);
 			sprintf(ps, ORDER_AT, npt_threshold); pr();
			extend_order(order,unplaced,&num_unplaced,
				     -npt_threshold,FALSE);
		    }
		}

		nl();
		sprintf(ps, "order%d= ", groups_done); pr();
		for (i=0; i<order->num_loci; i++) {
		    sprintf(ps, "%s ", rag(loc2str(order->locus[i]))); pr();
		    if (i==0) order_first[groups_done-1]=order->locus[i];
		      else order_next[prev]=order->locus[i]; 
		    prev=order->locus[i]; /* lint warning is OK */
		}
		nl();
		sprintf(ps, "other%d= ", groups_done); pr();
		for (i=0; i<num_unplaced; i++) {
		    sprintf(ps, "%s ", rag(loc2str(unplaced[i]->locus))); pr();
		    if (i==0) unorder_first[groups_done-1]=unplaced[i]->locus;
		      else order_next[prev]=unplaced[i]->locus;
		    prev=unplaced[i]->locus; /* lint warning is OK */
		}
		nl();
		num_orders++;

		if (print_all_maps) for (i=0; i<num_unplaced; i++)
		  if (unplaced[i]->best_map->num_loci>0) {
		      print(SUB_DIVIDER);
		      sprintf(ps, "Best placement of %s:",
                      rag(loc2str(unplaced[i]->locus)));
		      print_special_map(unplaced[i]->best_map,ps,
					order->num_loci,order->locus);
		  }

	    } /* if group size */
	} while (num_unlinked>0);
	
        if (groups_done==0) print(ORDER_NO_GROUPS);
	else print(MAP_DIVIDER);

    } on_exit {
	unarray(loci,int);
	unarray(linkage_group,int);
	free_map(order);
	free_map(seed_map);
	unarray(subset,int);
	unmatrix(seed_temp,seed_tries,int);
	for (i=0; i<num_loci; i++) {
	    unarray(unplaced[i]->excluded,bool);
	    free_map(unplaced[i]->best_map);
	}
	unparray(unplaced,num_loci,PLACEME);
	relay_messages;
    }
}


command 
greedy (void)
{
    int *locus=NULL, num_loci=0, num_order=0, *all_loci=NULL, total;
    int i=0, prev=0, num=0, num_unplaced=0;
    MAP *order=NULL;
    PLACEME **unplaced=NULL;

    /* args: lod theta start-window start-like max-to-try */
    mapm_ready(F2_DATA,2,LIST_SEQ,&num_order); /* F2 because place_locus */
    if (nullstr(args)) usage_error(0);
    parse_locus_args(&locus,&num_loci); /* error if fails, >=1 */
    crunch_locus_list(locus,&num_loci,CRUNCH_WARNINGS,ANY_CHROMS,IN_ARGS);

    if (npt_threshold==npt_first_threshold) /* globals */
      sprintf(ps, "Placement Threshold %.2lf, Npt-Window %d\n",
              npt_threshold, npt_window);
    else 
      sprintf(ps, "Placement Threshold-1 %.2lf, Threshold-2 %.2lf, Npt-Window %d\n",
              npt_first_threshold, npt_threshold, npt_window);
    pr();

    if (num_order>=MAX_MAP_LOCI) error("starting order is too big");
    run {
            if (num_loci == 0) {
                printf("\n\nWARNING: greedy(): num_loci == 0, returning early.\n\n"); fflush(NULL);
                // Not returning early will lead to an exception in: parray(unplaced,num_loci,PLACEME);
                send(CRASH);
            }

    total= num_order+num_loci;
	order=allocate_map(total); /* maybe MAX_MAP */
	get_list_of_all_loci(seq,order->locus,&order->num_loci,total);
	array(all_loci,total,int);
	parray(unplaced,num_loci,PLACEME);
	for (i=0; i<num_loci; i++) {
	    array(unplaced[i]->excluded,total+1,bool);
	    unplaced[i]->best_map=allocate_map(total);
	}

	num_orders=0; /* global - deletes all pre-existing orders */
	for (i=0; i<raw.num_markers; i++) order_next[i]=NO_LOCUS;
	print(MAP_DIVIDER);

	total=0;
	for (i=0; i<order->num_loci; i++) all_loci[total++]=order->locus[i];
	for (i=0; i<num_loci; i++) all_loci[total++]=locus[i];
	sort_loci(all_loci,total);
	if (use_three_pt)
	  setup_3pt_data(all_loci,total,-three_pt_threshold);
	
	for (i=0; i<num_loci; i++) unplaced[i]->locus=locus[i];
	num_unplaced=num_loci;

	sprintf(ps, "%d Markers to order:\n", total); pr();
	for (i=0; i<total; i++) {
	    num=all_loci[i];
	    sprintf(ps, "%4d %s%s", num + 1, raw.locus_name[num],
                (use_haplotypes && haplotyped(num) ? "+":""));
	    pad_to_len(ps,14); pr();
	    if (i==total-1 || i%5==4) nl(); else print("  ");
	}

	nl(); sprintf(ps, ORDER_AT, npt_first_threshold); pr();
	extend_order(order,unplaced,&num_unplaced,-npt_first_threshold,TRUE);

	if (npt_first_threshold!=npt_threshold && num_unplaced>0) {
	    print(SUB_DIVIDER);
	    sprintf(ps, ORDER_AT, npt_threshold); pr();
	    extend_order(order,unplaced,&num_unplaced,-npt_threshold,FALSE);
	}
	nl();
	print("order1= ");
	for (i=0; i<order->num_loci; i++) {
	    sprintf(ps, "%s ", rag(loc2str(order->locus[i]))); pr();
	    if (i==0) order_first[0]=order->locus[i];
	      else order_next[prev]=order->locus[i]; 
	    prev=order->locus[i]; /* lint warning is OK */
	}
	nl();
	print("other1= ");
	for (i=0; i<num_unplaced; i++) {
	    sprintf(ps, "%s ", rag(loc2str(unplaced[i]->locus))); pr();
	    if (i==0) unorder_first[0]=unplaced[i]->locus;
	      else order_next[prev]=unplaced[i]->locus;
	    prev=unplaced[i]->locus;
	}
	nl();
	num_orders=1;

	if (print_all_maps) for (i=0; i<num_unplaced; i++)
	  if (unplaced[i]->best_map->num_loci>0) {
	      print(SUB_DIVIDER);
	      sprintf(ps, "Best placement of %s:", rag(loc2str(unplaced[i]->locus)));
	      print_special_map(unplaced[i]->best_map,ps,
				order->num_loci,order->locus);
	  }
	print(MAP_DIVIDER);	

    } on_exit {
	unarray(locus,int);
	unarray(all_loci,int);
	free_map(order);
	for (i=0; i<num_loci; i++) {
	    unarray(unplaced[i]->excluded,bool);
	    free_map(unplaced[i]->best_map);
	}
	unparray(unplaced,num_loci,PLACEME);
	relay_messages;
    }
}
