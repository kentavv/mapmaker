/******************************************************************************

   ##    #    #   #####   ####            ####   #    #  #####            ####
  #  #   #    #     #    #    #          #    #  ##  ##  #    #          #    #
 #    #  #    #     #    #    #          #       # ## #  #    #          #
 ######  #    #     #    #    #          #       #    #  #    #   ###    #
 #    #  #    #     #    #    #          #    #  #    #  #    #   ###    #    #
 #    #   ####      #     ####  #######   ####   #    #  #####    ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

//#define INC_LIB
//#define INC_SHELL
//#define INC_MISC
#include "mapm.h"
//#include "toplevel.h"
//#include "lowlevel.h"

/* external proccedures WHY WHY WHY FIX FIX FIX */
//void try_marker();

/* local */
//void maybe_print_together();

#define ASS_ALREADY "%s- already assigned to %s at LOD %4.1lf...not changing\n"


/**************** COMMANDS ****************/

#define CANT_MAKE "warning: can't make chromosome '%s'... %s\n"

command 
make_chromosome (void)
{
    char name[TOKLEN+1];
    int i, num;

    mapm_ready(ANY_DATA,0,0,NULL);

    if (!nullstr(args)) {
	while (stoken(&args,sREQUIRED,name)) {
	    if (!valid_name(name)) 
	      { sprintf(ps, CANT_MAKE, name, "illegal name"); pr(); }
	    else if (!valid_new_name(name))
	      { sprintf(ps, CANT_MAKE, name, "name is already in use"); pr(); }
	    else if (!make_new_chrom(name,&num))
	      { sprintf(ps, CANT_MAKE, name, "chromosome already exists"); pr(); }
	}
    }
    print("chromosomes defined: ");
    if (num_chromosomes==0) print("none");
    else for (i=0; i<num_chromosomes; i++)
      { print(chrom2str(i)); print(" "); }
    nl();
}


command 
set_anchors (void)
{
    int i, *locus=NULL, num_loci, chrom;

    mapm_ready(ANY_DATA,MAYBE_SEQ,LIST_SEQ,&num_loci);
    chrom=get_chrom_arg(FALSE);
    nomore_args(num_args);

    run {
	if (num_loci>0)	{
	    alloc_list_of_all_loci(seq,&locus,&num_loci);
	    for (i=0; i<num_loci; i++)
	      if (!is_assignable(locus[i],chrom,FALSE))
		abort_command(); /* prints msg - is this a sufficient test? */
	    set_chrom_anchors(chrom,locus,num_loci);
	    sprintf(ps, "chromosome %s anchor(s): ", chrom2str(chrom)); pr();
	    for (i=0; i<raw.num_markers; i++) /* requires haplo sanity */
	      if (assigned_to(i,chrom) && anchor_locus(i))
		{ print(locname(i,TRUE)); print(" "); }
	    nl();
	} else {
	    set_chrom_anchors(chrom,locus,0);
	    sprintf(ps, "chromosome %s now has no anchors\n", chrom2str(chrom)); pr();
	}
    } on_exit { 
	unarray(locus,int);
	relay_messages; 
    }
}


#define CHROM_SETTING "setting framework for chromosome %s..."
#define CHROM_FRAME_EXISTS \
  "chromosome %s framework already exists, changing it..."
#define CHROM_NOFRAME "chromosome %s framework is %sempty..."

#define CHROM_FRAME "%s framework:\n"
#define CHROM_MARKS "%s Markers:\n"

command 
set_framework (void)
{
    MAP *map, *old;
    char title[TOKLEN+1];
    int chrom, num_loci;

    mapm_ready(ANY_DATA,MAYBE_SEQ,ONE_ORDER,&num_loci);
    chrom=get_chrom_arg(FALSE);
    nomore_args(num_args);
    if (num_loci>MAX_CHROM_LOCI)
      error("too many loci in the sequence for a chromosome framework");

    old=get_chrom_frame(chrom,NULL);
    map=get_map_to_bash(chromosome);
    if (num_loci==0) {
	sprintf(ps, CHROM_NOFRAME, chrom2str(chrom), (old->num_loci > 0 ? "now " : ""));
	pr(); nl();
	clean_map(map);
	set_chrom_frame(chrom,map);
	return;
    }
    run {
	get_one_order(seq,map);
	if (old->num_loci>0)
	  { sprintf(ps, CHROM_FRAME_EXISTS, chrom2str(chrom)); pr(); nl(); }
	else { sprintf(ps, CHROM_SETTING, chrom2str(chrom)); pr(); nl(); }
    
	init_rec_fracs(map);
	converge_to_map(map);
	set_chrom_frame(chrom,map);
	
	print(MAP_DIVIDER);
	sprintf(title, CHROM_FRAME, chrom2str(chrom));
	print_long_map(map,title);
	print(MAP_DIVIDER);
	
    } on_exit { relay_messages; }
}


#define CHROMS_TITLE \
" Chromosome:  #Total  #Frame  #Anchors  #Placed  #Unique  #Region\n"
/*2345678.......234.....234......234.......234......234......234*/
#define CHROMS_FORMAT \
"   %-8s    %4d    %4d     %4d      %4d     %4d     %4d\n"

command 
list_chroms (void)
{
    int frame, total, anchor, placed, unique, region, i, *loci=NULL;
    int sum_frame, sum_total, sum_anchor, sum_placed, sum_unique, sum_region;
    
    mapm_ready(ANY_DATA,0,0,NULL);
    nomore_args(0);
    if (num_chromosomes==0) error("no chromosomes presently defined\n");

    array(loci,raw.num_markers,int);
    sum_frame=sum_total=sum_anchor=sum_placed=sum_unique=sum_region= 0;
    nl(); print(CHROMS_TITLE);
    
    for (i=0; i<chromosome->num_maps; i++) {
	count_chrom_loci(i,&anchor,&frame,&total,&placed,&unique,&region,
			 TRUE,loci);
	sprintf(ps, CHROMS_FORMAT, chrom2str(i), total, frame, anchor, placed,
            unique, region); pr();
	sum_frame+=frame; sum_total+=total; sum_anchor+=anchor;
	sum_placed+=placed; sum_unique+=unique; sum_region+=region;
    }
    sprintf(ps, CHROMS_FORMAT, "Total:", sum_total, sum_frame, sum_anchor,
            sum_placed, sum_unique, sum_region); pr();
    unarray(loci,int);
}


command 
list_assignments (void)
{
    int chrom, i;
    
    mapm_ready(ANY_DATA,0,0,NULL);
    nomore_args(0);
    if (num_chromosomes==0) error("no chromosomes presently defined\n");
    print("Chromosome assignments of all loci:\n"); nl();

    for (chrom=0; chrom<chromosome->num_maps; chrom++) {
	sprintf(ps, "%s= ", chrom2str(chrom)); pr();
	for (i=0; i<raw.num_markers; i++) if (assigned_to(i,chrom))
	  { print(rag(loc2str(i))); print(" "); }
	nl();
	print("-------\n");
    }
    print("unassign= ");
    for (i=0; i<raw.num_markers; i++) if (!assigned(i))
      { print(rag(loc2str(i))); print(" "); }
    nl();
}


command 
list_mapping (void)
{
    int source, num_loci, *locus=NULL;

    mapm_ready(ANY_DATA,MAYBE_SEQ,UNCRUNCHED_LIST,&num_loci);
    run {
	if (!nullstr(args)) {
	    parse_locus_args(&locus,&num_loci); /* error if fails */
	    if (num_loci==0) error("no loci in arguments\n");
	    source=IN_ARGS;
	} else {
	    if (!alloc_list_of_all_loci(seq,&locus,&num_loci))
	      error(NEED_SEQ_OR_ARGS);
	    source=IN_SEQ;
	}
	crunch_locus_list(locus,&num_loci,CRUNCH_WARNINGS,ANY_CHROMS,source);
	nl(); 
	hold(more_mode) print_mapping_summary(locus,num_loci,use_haplotypes);
    } on_exit {
	unarray(locus,int);
	relay_messages;
    }
}


command 
assign (void)
{
    int i, n_loci, *locus=NULL;
    real theta, lod, min_lod, unlinked_lod;

    mapm_ready(ANY_DATA,1,UNCRUNCHED_LIST,&n_loci);

    if (!rtoken(&args,default_lod,&lod) || !rrange(&lod,0.0,1000.0))
      error("Bad value for LOD score bound");
    if (!rtoken(&args,default_theta,&theta) || !input_dist(&theta))
      error("Bad value for theta bound");
    if (!rtoken(&args,lod,&unlinked_lod) || !rrange(&unlinked_lod,0.0,lod))
      error("Bad value for maximum unlinked lod bound");
    if (!rtoken(&args,lod,&min_lod) || !rrange(&min_lod,0.0,lod))
      error("Bad value for marginal LOD bound");
    nomore_args(0);

    run {
	alloc_list_of_all_loci(seq,&locus,&n_loci);
	crunch_locus_list(locus,&n_loci,CRUNCH_WARNINGS,ANY_CHROMS,IN_SEQ);
	for (i=0; i<n_loci; i++)
	  if (!is_assignable(locus[i],ANY_CHROM,TRUE)) locus[i]=NO_LOCUS;
	do_assignments(locus,n_loci,lod,unlinked_lod,theta,
		       min_lod,min_lod,theta,FALSE);  /* FIX */
	
    } on_exit {
	unarray(locus, int);
	relay_messages;
    }
}


command 
attach (void)
{
    int i, chrom, n_loci, *locus=NULL;

    mapm_ready(ANY_DATA,1,UNCRUNCHED_LIST,&n_loci);
    chrom=get_chrom_arg(FALSE);
    nomore_args(num_args);

    run {
	alloc_list_of_all_loci(seq,&locus,&n_loci);
	crunch_locus_list(locus,&n_loci,CRUNCH_WARNINGS,ANY_CHROMS,IN_SEQ);
	for (i=0; i<n_loci; i++) 
	  if (is_assignable(locus[i],ANY_CHROM,TRUE)) {
	      if (assigned_to(locus[i],chrom) && 
		  (assignment_state(locus[i])==A_ASSIGNED || 
		   assignment_state(locus[i])==A_BORDERLINE))
		{ sprintf(ps, ASS_ALREADY, loc2str(locus[i]),
                  chrom2str(assignment_chrom(locus[i])),
                  assignment_lod(locus[i])); pr(); continue; }
	      assign_this(locus[i],A_ATTACH,chrom,NO_LOD,NO_THETA,NO_LOCUS,"");
	  }
	print("ok\n");
	
    } on_exit {
	unarray(locus, int);
	relay_messages;
    }
}


command 
unassign (void)
{
    int n_loci, i, *locus=NULL;

    mapm_ready(ANY_DATA,1,UNCRUNCHED_LIST,&n_loci);
    nomore_args(0);

    run {
	alloc_list_of_all_loci(seq,&locus,&n_loci);
	crunch_locus_list(locus,&n_loci,CRUNCH_WARNINGS,ANY_CHROMS,IN_SEQ);
	for (i=0; i<n_loci; i++) {
	    if (is_assignable(locus[i],ANY_CHROM,TRUE))
	      assign_this(locus[i],A_UNKNOWN,NO_CHROM,NO_LOD,
			  NO_THETA,NO_LOCUS,"");
	}
    } on_exit {
	unarray(locus,int);
	relay_messages;
    }
}


#define PLACE_HEADER  "Placing markers at log-likelihood threshold %.2lf...\n"
#define PLACE_SOME    "%d marker%s to place on chromosome %s...\n"
#define PLACE_NOFRAME "chromosome has no framework order - no markers placed\n"
#define PLACE_NONE    "unable to place any markers in the sequence"
#define NONE_OK "no orders allowed by three-point analysis, not placed\n"
#define NPT_WINDOW 5
command 
place (void)
{
    int  chrom, i, left, right;
    int  *locus=NULL, n_loci, *chrom_locus=NULL, n_to_place, n_allowed;
    int  *chrom_all=NULL, n_chrom, n_frame, n_to_punt;
    real npt_thresh;
    bool got_any;
    PLACE **placements=NULL;
    PLACEME **unplaced=NULL;
    MAP *temp=NULL, *frame; /* frame is NOT malloced */
  
    mapm_ready(ANY_DATA,1,LIST_SEQ,&n_loci);
    if (!rtoken(&args,2.0,&npt_thresh) || 
	!rrange(&npt_thresh,0.0,-1.0 * PLACEMENT_THRESHOLD))
      error("Bad value for multipoint likelihood threshold");
    
    run {
	alloc_list_of_all_loci(seq,&locus,&n_loci);
	array(chrom_locus,raw.num_markers,int);
	array(chrom_all,raw.num_markers,int);
	parray(placements,MAX_CHROM_LOCI+1,PLACE);
	parray(unplaced,n_loci,PLACEME);
	for (i=0; i<n_loci; i++) {
	    array(unplaced[i]->excluded,MAX_CHROM_LOCI+1,bool);
	    unplaced[i]->best_map=allocate_map(MAX_CHROM_LOCI+1);
	}
	temp=allocate_map(MAX_CHROM_LOCI+1);

	sprintf(ps, PLACE_HEADER, npt_thresh); pr();
	for (i=0; i<n_loci; i++) /* is_placable may print a msg & unplace */
	  if (!is_placeable(locus[i],ANY_CHROM)) locus[i]=NO_LOCUS;

	got_any= FALSE;
	for (chrom=0; chrom<chromosome->num_maps; chrom++) {
	    /* chrom_locus[], n_to_place is all loci in seq on this chrom */
	    n_to_place=0;
	    for (i=0; i<n_loci; i++) if (locus[i]!=NO_LOCUS)
	      if (assigned_to(locus[i],chrom) && is_placeable(locus[i],chrom))
		{ chrom_locus[n_to_place++]=locus[i]; locus[i]=NO_LOCUS; }
	    if (n_to_place==0) continue; /* next chrom */
	    
	    frame= get_chrom_frame(chrom,&n_frame);
	    /* chrom_all has the framework first, then chrom_locus loci */
	    get_chrom_loci(chrom,chrom_all,ALL_LOCI,&n_chrom,NULL);
	    sort_loci(chrom_all,n_chrom);

	    /* if (!got_any) nl(); */
	    print(MAP_DIVIDER);
	    sprintf(ps, PLACE_SOME, n_to_place, maybe_s(n_to_place), chrom2str(chrom));
	    pr(); nl();
	    if (n_frame<2) { sprintf(ps, PLACE_NOFRAME); pr(); nl(); continue; }
	    got_any=TRUE; /* yes, do this here */
	    /* if (use_haplotypes) 
	       { print_haplo_summary(chrom_all,n_chrom); nl(); } */
	    print_long_map(frame,"Framework Map:"); nl();

	    if (use_three_pt)
	      setup_3pt_data(chrom_all,n_chrom,-three_pt_threshold);

	    n_to_punt=0;
	    print("Placing Markers:\n");
	    for (i=0; i<n_to_place; i++) {
		unplaced[i]->locus=chrom_locus[i];

		if (use_three_pt) {
		    /* three_pt_exclusions sets excluded[x] to FALSE first */
		    n_allowed=three_pt_exclusions(frame->locus,frame->num_loci,
						  unplaced[i]->locus,
						  unplaced[i]->excluded);
		    find_window(frame->locus,frame->num_loci,
				unplaced[i]->locus,unplaced[i]->excluded,
				NPT_WINDOW,&left,&right);
		} else { /* no three-pt -> all orders are OK */
		    n_allowed=n_frame+1; 
		    left=0; right=n_frame; /* window= all intervals */
		}
		
		if (n_allowed==0) { /* (un)place_this() does all of the i/o */
		    unplace_this(unplaced[i]->locus,chrom,M_PROBLEM,TRUE);
		    unplaced[i]->locus=NO_LOCUS; n_to_punt++;

		} else {
		    place_locus(frame,unplaced[i]->locus,left,right,
				unplaced[i]->excluded,placements,
				&unplaced[i]->best_pos,unplaced[i]->best_map,
				temp);
		    place_this(unplaced[i]->locus,chrom,placements,-npt_thresh,
			       10000.0,error_net_thresh,unplaced[i]->excluded);
		    /* place_this() calls npt_exclusions() on placements */
		}
	    }

	    nl();
	    if (n_to_place-n_to_punt==0) print("No Valid Placements\n");
	    else {
		print("Placements:\n");
		new_print_placements(frame,unplaced,n_to_place);
		for (i=0; i<n_to_place; i++)
		  if (print_all_maps && unplaced[i]->locus!=NO_LOCUS && 
		        placed_locus(unplaced[i]->locus)) {
		      print(SUB_DIVIDER);
		      sprintf(ps, "Best placement of %s:",
                      rag(loc2str(unplaced[i]->locus)));
		      print_special_map(unplaced[i]->best_map,ps,
					n_frame,frame->locus);
		  }
	    }

	} /* next chrom */
	if (got_any) print(MAP_DIVIDER);
	else error(PLACE_NONE);

    } on_exit {
	unarray(locus,int);
	unparray(placements,MAX_CHROM_LOCI+1,real);
	unarray(chrom_locus,int);
	unarray(chrom_all,int);
	free_map(temp);
	for (i=0; i<n_loci; i++) {
	    free_map(unplaced[i]->best_map);
	    unarray(unplaced[i]->excluded,bool);
	}
	unparray(unplaced,n_loci,PLACEME);
	relay_messages;
    }
}


command 
show_chrom (void)
{
    int chrom;
    int *marker=NULL, num_markers, **placement_state=NULL;
    char title[TOKLEN+1];

    mapm_ready(ANY_DATA,0,0,NULL);
    chrom=get_chrom_arg(FALSE);
    nomore_args(num_args);

    run {
	array(marker,raw.num_markers,int); /* non-frame markers */
	get_chrom_loci(chrom,marker,NON_FRAME,&num_markers,NULL);
	matrix(placement_state,num_markers,
	       chromosome->map_list[chrom]->num_loci+1,int);

	nl();
	print(MAP_DIVIDER);
	sprintf(title, CHROM_FRAME, chrom2str(chrom));
	print_long_map(chromosome->map_list[chrom],title);
	print(MAP_DIVIDER);
	
	sprintf(ps, CHROM_MARKS, chrom2str(chrom)); pr(); nl();
	get_chrom_loci(chrom,marker,ALL_LOCI,&num_markers,NULL);
	print_locus_summary(marker,num_markers,TRUE); /* change to FALSE? */
	if (use_haplotypes) { nl(); print_haplo_summary(marker,num_markers); }
	print(MAP_DIVIDER);

/*	get_chrom_loci(chrom,marker,&num_markers,...);
	print_3pt_placements(chromosome->map_list[chrom]->locus,
			     chromosome->map_list[chrom]->num_loci,
			     marker,num_markers,placement_state);    
    	print(MAP_DIVIDER);  */
	
    } on_exit {
	unarray(marker,int);
	relay_messages;
    }
}


#define PLACE_AT \
"Placing markers at log-likelihood threshold %.1lf:\n"

command 
place_together (void)
{
    int  chrom, i;
    int  *locus=NULL, n_loci, *chrom_locus=NULL, n_to_place;
    int  *chrom_all=NULL, n_chrom, n_frame, n_unplaced, num, prev;
    real npt_thresh = 0.;
    bool got_any;
    PLACEME **unplaced=NULL;
    MAP *order=NULL, *frame; /* frame is NOT malloced */
  
    mapm_ready(ANY_DATA,1,LIST_SEQ,&n_loci);
    nomore_args(0);
    if (npt_threshold==npt_first_threshold) /* globals */
      sprintf(ps, "Placement Threshold %.2lf, Npt-Window %d\n",
              npt_threshold, npt_window);
    else 
      sprintf(ps, "Placement Threshold-1 %.2lf, Threshold-2 %.2lf, Npt-Window %d\n",
              npt_first_threshold, npt_threshold, npt_window);
    pr();

    run {
	alloc_list_of_all_loci(seq,&locus,&n_loci);
	array(chrom_locus,raw.num_markers,int);
	array(chrom_all,raw.num_markers,int);
	parray(unplaced,n_loci,PLACEME);
	for (i=0; i<n_loci; i++) {
	    array(unplaced[i]->excluded,MAX_CHROM_LOCI+1,bool);
	    unplaced[i]->best_map=allocate_map(MAX_CHROM_LOCI+1);
	}
	order=allocate_map(MAX_CHROM_LOCI+1);

	sprintf(ps, PLACE_HEADER, npt_thresh); pr();
	for (i=0; i<n_loci; i++) /* is_placable may print a msg & unplace */
	  if (!is_placeable(locus[i],ANY_CHROM)) locus[i]=NO_LOCUS;

	got_any= FALSE;
	for (chrom=0; chrom<chromosome->num_maps; chrom++) {
	    /* chrom_locus[], n_to_place is all loci in seq on this chrom */
	    n_to_place=0;
	    for (i=0; i<n_loci; i++) if (locus[i]!=NO_LOCUS)
	      if (assigned_to(locus[i],chrom) && is_placeable(locus[i],chrom))
		{ chrom_locus[n_to_place++]=locus[i]; locus[i]=NO_LOCUS; }
	    if (n_to_place==0) continue; /* next chrom */
	    
	    frame= get_chrom_frame(chrom,&n_frame);
	    /* chrom_all has the framework first, then chrom_locus loci */
	    get_chrom_loci(chrom,chrom_all,ALL_LOCI,&n_chrom,NULL);
	    sort_loci(chrom_all,n_chrom);

	    /* if (!got_any) nl(); */
	    print(MAP_DIVIDER);
	    sprintf(ps, PLACE_SOME, n_to_place, maybe_s(n_to_place), chrom2str(chrom));
	    pr(); nl();
	    if (n_frame<2) { sprintf(ps, PLACE_NOFRAME); pr(); nl(); continue; }
	    got_any=TRUE; /* yes, do this here */

	    /* if (use_haplotypes) 
	       { print_haplo_summary(chrom_all,n_chrom); nl(); } */
	    if (use_three_pt)
	      setup_3pt_data(chrom_all,n_chrom,-three_pt_threshold);

	    for (i=0; i<n_chrom; i++) {
		num=chrom_all[i];
		sprintf(ps, "%4d %s%s", num + 1, raw.locus_name[num],
                (use_haplotypes && haplotyped(num) ? "+":""));
		pad_to_len(ps,14); pr();
		if (i==n_chrom-1 || i%5==4) nl(); else print("  ");
	    }

	    mapcpy(order,frame,TRUE);
	    for (i=0; i<n_to_place; i++) unplaced[i]->locus=chrom_locus[i];
	    n_unplaced= n_to_place;

	    nl(); sprintf(ps, PLACE_AT, npt_first_threshold); pr();
	    extend_order(order,unplaced,&n_unplaced,-npt_first_threshold,TRUE);

	    if (npt_first_threshold!=npt_threshold && n_unplaced>0) {
		print(SUB_DIVIDER);
		sprintf(ps, PLACE_AT, npt_threshold); pr();
		extend_order(order,unplaced,&n_unplaced,-npt_threshold,FALSE);
	    }

	    nl();
	    sprintf(ps, "order%d= ", num_orders + 1); pr();
	    for (i=0; i<order->num_loci; i++) {
		sprintf(ps, "%s ", rag(loc2str(order->locus[i]))); pr();
		if (i==0) order_first[0]=order->locus[i];
		else order_next[prev]=order->locus[i]; 
		prev=order->locus[i]; /* lint warning is OK */
	    }
	    nl();
	    sprintf(ps, "other%d= ", num_orders + 1); pr();
	    for (i=0; i<n_unplaced; i++) {
		sprintf(ps, "%s ", rag(loc2str(unplaced[i]->locus))); pr();
		if (i==0) unorder_first[0]=unplaced[i]->locus;
		else order_next[prev]=unplaced[i]->locus;
		prev=unplaced[i]->locus;
	    }
	    nl();
	    num_orders++;

	    if (print_all_maps) for (i=0; i<n_unplaced; i++)
	      if (unplaced[i]->best_map->num_loci>0) {
		  print(SUB_DIVIDER);
		  sprintf(ps, "Best placement of %s:",
                  rag(loc2str(unplaced[i]->locus)));
		  print_special_map(unplaced[i]->best_map,ps,
				    order->num_loci,order->locus);
	      }

	} /* next chrom */
	if (got_any) print(MAP_DIVIDER); else error(PLACE_NONE);

    } on_exit {
	/* add FREEs */
	relay_messages;
    }
}


command 
draw_chromosome (void)
{
    int chrom;
    char name[PATH_LENGTH+1];
    FILE *fp=NULL;

    mapm_ready(ANY_DATA,0,0,NULL);
    use_uncrunched_args();
    chrom=get_chrom_arg(FALSE);
    get_arg(stoken,chrom2str(chrom),name);
    nomore_args(num_args);

    run {
	if (!make_filename(name,FORCE_EXTENSION,PS_EXT))
	  { sprintf(ps, "illegal filename '%s'", name); error(ps); }
	fp= open_file(name,WRITE);
	sprintf(ps, "Drawing chromosome %s in PostScript file '%s'... \n",
            chrom2str(chrom), name); pr();
	print_ps_chrom(fp,chrom);
	close_file(fp);
	print("ok\n");

    } except {
	when CANTOPEN: sprintf(ps, "Can't create output file '%s'", name); error(ps); break;
	when CANTCLOSE: error("\nCan't close file - disk is full?"); break;
    default: close_file(fp); relay_messages; break;
    }
}


command 
draw_all_chromosomes (void)
{
    char name[PATH_LENGTH+1];
    FILE *fp=NULL;

    mapm_ready(ANY_DATA,0,0,NULL);
    use_uncrunched_args();
    get_arg(stoken,raw.filename,name);
    nomore_args(num_args);
    
    run {
	/* print all chroms in same file to same scale */
	if (!make_filename(name,FORCE_EXTENSION,PS_EXT)) 
	  { sprintf(ps, "illegal filename '%s'", name); error(ps); }
	fp=open_file(name,WRITE);
	sprintf(ps, "Drawing all chromosomes in PostScript file '%s'... \n", name); 
	pr();
	print_all_ps_chroms(fp);
	close_file(fp);
	print("ok\n");
	
    } except {
	when CANTOPEN: sprintf(ps, "Can't create output file '%s'", name); error(ps); break;
    when CANTCLOSE: error("\nCan't close file - disk is full?"); break;
	default: close_file(fp); relay_messages; break;
    }
}
