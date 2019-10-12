/******************************************************************************

  ####   #    #  #####    ####   #    #   ####            ####
 #    #  #    #  #    #  #    #  ##  ##  #               #    #
 #       ######  #    #  #    #  # ## #   ####           #
 #       #    #  #####   #    #  #    #       #   ###    #
 #    #  #    #  #   #   #    #  #    #  #    #   ###    #    #
  ####   #    #  #    #   ####   #    #   ####    ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#define INC_LIB
#define INC_MISC
#define INC_SHELL
#include "mapm.h"

SAVED_LIST *chromosome;
ASSIGNMENT **assignment;
PLACEMENT  **placement;
int current_chrom;
bool do_assignment();
int  get_anchors();


/* do_assignments anchor_locus framework_marker set_anchor_loci */

/**************** Chromosomes ****************/

bool make_new_chrom(name,num)
char *name;
int *num;
{
    int n, m;

    if ((n=chromosome->num_maps)==MAX_CHROMOSOMES-1)
      error("too many chromosomes - can't make any more");
    if (isa_chrom(name,&m)) return(FALSE);

    strcpy(chromosome->map_list[n]->map_name,name);
    chromosome->map_list[n]->num_loci=0;
    chromosome->num_maps+=1;
    *num=n;
    return(TRUE);
}


bool isa_chrom(name,chrom)
char *name;
int *chrom; /* side-effected */
{
    int i;

    *chrom= NO_CHROM;
    if (xstreq(name,"any")) { *chrom=NO_CHROM; return(TRUE); }
    for (i=0; i<chromosome->num_maps; i++)
      if (xstreq(name,chromosome->map_list[i]->map_name)) { *chrom=i; break; }
    return (*chrom==NO_CHROM ? FALSE:TRUE);
}


MAP *
get_chrom_frame (
    int chrom,
    int *num_loci /* side-effected,  can be 0, or 1, or greater */
)
{
    if (num_loci!=NULL) *num_loci= chromosome->map_list[chrom]->num_loci;
    return(chromosome->map_list[chrom]);
}	


bool framework_marker(locus)
int locus;
{
    int j;
    MAP *frame;

    /* force_haplo_sanity(&locus,FALSE); */
    if (!assigned(locus)) return(FALSE);
    
    frame=chromosome->map_list[assignment_chrom(locus)];
    for (j=0; j<frame->num_loci; j++) {
	if (locus==frame->locus[j]) return(TRUE);
	/* if (!check_haplos) continue;
	   for (k=haplo_first[frame->locus[j]]; k!=NO_LOCUS; k=haplo_next[k]) 
	   if (locus==k) return(TRUE); */
    }
    return(FALSE);
}

   
#define FRAME_OTHER_CHROM \
  "bad framework for %s - %s is assigned to %s"
#define FRAME_NO_CHROM \
  "bad framework - %s is not assigned to a chromosome"
#define FRAME_DUPS \
  "bad framework - one or more loci are repeated in sequence"

void 
set_chrom_frame (
    int chrom,
    MAP *new   /* warning: side-effected! */
)
/* When this is called, new MUST be equal to get_map_to_bash(chromosome), amd
   assigned_to(*,chrom) must have be true for all loci in new map. We also 
   assume haplo_sanity is true for old framework, and force it for new one. */
{
    int i, j, found, other;
    MAP *old;

    old=get_chrom_frame(chrom,NULL);

    /* check for problems: markers on other chroms, haplo_sanity, duplicates */
    for (j=0; j<new->num_loci; j++) {
	/* force_haplo_sanity(new->locus+j,TRUE); */
	if (!assigned(new->locus[j])) {
	    sprintf(ps, FRAME_NO_CHROM, loc2str(new->locus[j])); error(ps);
	} else if ((other=assignment_chrom(new->locus[j]))!=chrom) {
	    sprintf(ps, FRAME_OTHER_CHROM, chrom2str(chrom), loc2str(new->locus[j]),
                chrom2str(other)); error(ps);
	}
    }

    /* unplace markers which WERE on the framework */
    for (i=0; i<old->num_loci; i++) {
	found=FALSE; 
	for (j=0; j<new->num_loci; j++)
	  if (old->locus[i]==new->locus[j]) { found=TRUE; break; }
	if (!found) unplace_this(old->locus[i],NO_CHROM,M_UNKNOWN,FALSE);
    }

    /* unplace ALL markers which were placed relative to the old framework */
    for (j=0; j<raw.num_markers; j++)
      if (assigned_to(j,chrom) && placed_locus(j))
	unplace_this(j,NO_CHROM,M_UNKNOWN,FALSE);

    /* quasi-place new markers in the new framework on this chrom */
    for (i=0; i<new->num_loci; i++) {
	found=FALSE; 
	for (j=0; j<old->num_loci; j++)
	  if (old->locus[i]==new->locus[j]) { found=TRUE; break; }
	if (!found) unplace_this(new->locus[i],chrom,M_FRAMEWORK,FALSE);
    }

    /* save the framework */
    strcpy(new->map_name,chromosome->map_list[chrom]->map_name);
    overwrite_map_num(chromosome,&new,chrom);
    chromosome->map_list[chrom]->modified= TRUE;
}


void 
get_chrom_loci (int chrom, int *locus, int which_loci, int *num_loci, int *num_framework)
{
    int i, j, k, n, primary;
    bool framework=FALSE, non_frame=FALSE, haplos_anyway=FALSE;

    if (which_loci>=HAPLOS)    { haplos_anyway=TRUE; which_loci-=HAPLOS; }
    if (which_loci>=NON_FRAME) { non_frame=TRUE; which_loci-=NON_FRAME;  }
    if (which_loci>=FRAMEWORK) { framework=TRUE; which_loci-=FRAMEWORK;  }

    n=0;
    if (framework) 
      for (j=0; j<chromosome->map_list[chrom]->num_loci; j++) {
	  i=chromosome->map_list[chrom]->locus[j];
	  if (use_haplotypes && haplos_anyway && haplotyped(i)) {
	      for (k=haplo_first[j]; k!=NO_LOCUS; k=haplo_next[k])
		locus[n++]=k;
	  } else { /* ignore haplos, and assume chrom frame is sane */
	      locus[n++]= i;
	  }
      }
    if (num_framework!=NULL) *num_framework=n;

    if (non_frame)
      for (i=0; i<raw.num_markers; i++) {
	  if (use_haplotypes && !haplos_anyway && haplotyped(i) &&
	      haplo_first[i]!=i) continue; /* then skip non-primaries */
	  else if (use_haplotypes && haplotyped(i)) primary=haplo_first[i];
	  else primary=i;
	  if (assigned_to(primary,chrom) && !framework_marker(primary))
	    locus[n++]=i;
      }
    if (num_loci!=NULL) *num_loci=n; /* the total! */
}	


void count_chrom_loci(chrom,n_anchor,n_frame,n_total,n_placed,n_unique,
		      n_region,haplos_too,temp)
int chrom;
int *n_anchor, *n_frame, *n_total, *n_placed, *n_unique, *n_region;
bool haplos_too;
int *temp; /* at most raw.num_markers long, alloced outside */
{
    int which, i, locus, primary;
    int frame, total, anchor, placed, unique, region;
    
    if (use_haplotypes && haplos_too) which=ALL_LOCI+HAPLOS; 
      else which=ALL_LOCI;
    get_chrom_loci(chrom,temp,which,&total,&frame);
    
    anchor=placed=unique=region=0;
    for (i=frame; i<total; i++) { /* Loop over non-frame markers */
	locus=temp[i];
	/* if !haplos_too then only primaries will be in the list */
	if (use_haplotypes && haplotyped(locus)) primary=haplo_first[locus];
	  else primary=locus; 
	if (anchor_locus(primary)) anchor++;
	if (placed_locus(primary)) {
	    placed++;
	    if (placement_state(primary)==M_UNIQUE) unique++;
	      else if (placement_state(primary)==M_REGION) region++;
	}
    }
    if (n_anchor!=NULL) *n_anchor=anchor;
    if (n_frame!=NULL)  *n_frame= frame;
    if (n_total!=NULL)  *n_total= total;
    if (n_placed!=NULL) *n_placed=placed;
    if (n_unique!=NULL) *n_unique=unique;
    if (n_region!=NULL) *n_region=region;
}


/**************** Assignments ****************/

bool assigned(locus)
int locus; 
{ force_haplo_sanity(&locus,FALSE); return(assignment[locus]->status>0); }

int 
assignment_state (int locus)
{ force_haplo_sanity(&locus,FALSE); return(assignment[locus]->status); }

int 
assignment_chrom (int locus)
{ 
  force_haplo_sanity(&locus,FALSE);
  return(assignment[locus]->status>0 ? assignment[locus]->chromosome:NO_CHROM);
}

bool assigned_to(locus,chrom)
int locus, chrom; 
{
  force_haplo_sanity(&locus,FALSE);
  return(assignment[locus]->status>0 && assignment[locus]->chromosome==chrom);
}

bool anchor_locus(locus)
int locus;
{ return(assignment_state(locus)==A_ANCHOR); } /* does force_haplo_sanity */

    
#define ASS_ANCHOR  "%s- anchor locus on %s...cannot re-assign\n"
#define ASS_REFRAME "%s- framework marker on %s...cannot re-assign to %s\n"
#define ASS_UNFRAME "%s- framework marker on %s...cannot un-assign\n"
#define ASS_ISFRAME "%s- framework marker on %s...remains attached\n"
#define ASS_INSANE  \
  "%s- not the name of its haplotype-group, can't assign or unassign\n"

bool is_assignable(locus,chrom,fix_frames)
int locus, chrom; /* maybe NO_CHROM */
bool fix_frames;
{
    int old;

    /* if (!force_haplo_sanity(&locus,FALSE)) { 
	sprintf(ps,ASS_INSANE,loc2str(locus)); pr();
	return(FALSE);
    } */
      
    if (!assigned(locus)) return(TRUE);
    old= assignment_chrom(locus);

    if (anchor_locus(locus) && old!=chrom) {
	sprintf(ps, ASS_ANCHOR, loc2str(locus), chrom2str(old)); pr();
	return(FALSE); /* don't touch assignment */
    }
    if (framework_marker(locus) && old!=chrom) {
	/* includes chrom==NO_CHROM case, eg unassigning */
	if (!fix_frames) {
	    if (chrom==NO_CHROM) 
	      sprintf(ps, ASS_UNFRAME, loc2str(locus), chrom2str(old));
	    else 
	      sprintf(ps, ASS_REFRAME, loc2str(locus), chrom2str(old),
                  chrom2str(chrom));
	    pr();
	    return(FALSE);
	} else {
	    /* need to change assignment - lod has disappeared */
	    /* keep old assignment[locus]->chromosome */
	    assignment[locus]->status= A_FRAMEWORK;
	    assignment[locus]->LODscore= NO_LOD;
	    assignment[locus]->theta= NO_THETA;
	    assignment[locus]->linked_to= NO_LOCUS;
	    assignment[locus]->modified= TRUE;
	    /* no need to unplace - it's (still) a framework marker */
	    sprintf(ps, ASS_ISFRAME, loc2str(locus), chrom2str(old)); pr();
	    return(FALSE); /* ??? */
	}
    }
    return(TRUE);
}


#define ASS_NONE    "%s- unassigned\n"
#define ASS_UNASS   "%s- unassigned (previously assigned to chrom %s)\n"
#define ASS_UNANCH  "%s- unassigned (previously anchor locus on %s)\n"
#define ASS_NEWANCH "%s- unassigned (anchor loci on %s have been changed)\n"

#define ASS_FRAME   "%s- framework locus...remains attached to %s\n"
#define ASS_ASSIGN  "%s- assigned to %s at LOD %4.1f\n"
#define ASS_BORDER  "%s- assigned to %s at LOD *%4.1f*\n"
#define ASS_ATTACH  "%s- attached to %s\n"

#define CHG_ASSIGN "%s- assigned to %s at LOD %4.1f (previously assigned to %s)\n"
#define CHG_ATTACH "%s- attached to %s (previously assigned to %s)\n"
#define CHG_BORDER "%s- assigned to %s at LOD *%.1f* (previously assigned to %s)\n"

/* this will unassign anchors, but not reassign them correctly, and it will 
   certainly not make loci anchors - this is OK for now */
void assign_this(locus,state,chrom,lod,theta,linked_locus,msg)
int locus, state, chrom;
real lod, theta;
int linked_locus;
char *msg; /* pre-empts other message, only for A_PROBLEM as yet */
{
    if (state<=0) { /* set to an unassigned state */
	if (!nullstr(msg)) /* problem state: use this msg */
	  strcpy(ps,msg);
	else if (anchor_locus(locus))
	  sprintf(ps, ASS_UNANCH, loc2str(locus), chrom2str(assignment_chrom(locus)));
	else if (assigned(locus)) 
	  sprintf(ps, ASS_UNASS, loc2str(locus), chrom2str(assignment_chrom(locus)));
	else 
	  sprintf(ps, ASS_NONE, loc2str(locus));
	pr();
	unplace_this(locus,NO_CHROM,M_UNKNOWN,FALSE);

    } else { /* assign to a chrom */
	if (assigned(locus) && assignment_chrom(locus)!=chrom) { /* changed */
	    if (framework_marker(locus)) send(CRASH);
	    if (state==A_ATTACH) /* no lod state */
	      sprintf(ps, CHG_ATTACH, loc2str(locus), chrom2str(chrom),
                  chrom2str(assignment_chrom(locus)));
	    else if (state==A_FRAMEWORK) send(CRASH);
	    else if (state==A_BORDERLINE)
	      sprintf(ps, CHG_BORDER, loc2str(locus), chrom2str(chrom), lod,
                  chrom2str(assignment_chrom(locus)));
	    else /* state==A_ASSIGNED */
	      sprintf(ps, CHG_ASSIGN, loc2str(locus), chrom2str(chrom), lod,
                  chrom2str(assignment_chrom(locus)));
	    unplace_this(locus,NO_CHROM,M_UNKNOWN,FALSE);

	} else { /* assignment has not been changed */
	    if (state==A_ATTACH) /* no lod states */
	      sprintf(ps, ASS_ATTACH, loc2str(locus), chrom2str(chrom)); 
	    else if (state==A_FRAMEWORK) /* no lod states */
	      sprintf(ps, ASS_FRAME, loc2str(locus), chrom2str(chrom)); 
	    else if (state==A_BORDERLINE) 
	      sprintf(ps, ASS_BORDER, loc2str(locus), chrom2str(chrom), lod);
	    else
	      sprintf(ps, ASS_ASSIGN, loc2str(locus), chrom2str(chrom), lod);
	}
	pr();
    }

    assignment[locus]->chromosome= chrom;
    assignment[locus]->status= state;
    assignment[locus]->LODscore= lod;
    assignment[locus]->theta= theta;
    assignment[locus]->linked_to= linked_locus;
    assignment[locus]->modified= TRUE;
}


void 
unassign_this (int locus, int state)
{ assign_this(locus,state,NO_CHROM,NO_LOD,NO_THETA,NO_LOCUS,""); }


void 
attach_this (int locus, int state, int chrom)
{ assign_this(locus,state,chrom,NO_LOD,NO_THETA,NO_LOCUS,""); }


#define ANCHOR_ISFRAME \
  "%s is a framework marker on %s...cannot re-assign to %s\n"
#define ANCHOR_ISANCH \
  "%s is a framework marker on %s...cannot re-assign to %s\n"
#define ANCHOR_ISNOW "%s- anchor locus on %s\n"
#define ANCHOR_DUPS  "bad framework\none or more loci are repeated in sequence"

void 
set_chrom_anchors (
    int chrom,
    int *locus,
    int num_loci /* maybe 0 */
)
{ 
    int i, j, old_chrom;
    bool still_got_it;

    /* check for problems: FW/anchor markers on other chroms, 
       need to assume haplo_sanity */
    for (j=0; j<num_loci; j++) {
	old_chrom=(assigned(locus[j]) ? assignment_chrom(locus[j]):NO_CHROM);
	if (old_chrom!=chrom && framework_marker(locus[j])) {
	    sprintf(ps, ANCHOR_ISFRAME, loc2str(locus[j]), chrom2str(old_chrom), chrom2str(chrom));
	    error(ps);
	}
	if (old_chrom!=chrom && anchor_locus(locus[j])) {
 	    sprintf(ps, ANCHOR_ISANCH, loc2str(locus[j]), chrom2str(old_chrom), chrom2str(chrom));
	    error(ps);
	}
    }

    /* unassign old anchors */
    for (i=0; i<raw.num_markers; i++)
      if (assigned_to(i,chrom) && anchor_locus(i)) {
	  still_got_it=FALSE;
	  for (j=0; j<num_loci; j++) 
	    if (locus[j]==i) { still_got_it=TRUE; break; }
 	  if (!still_got_it) {
	      if (framework_marker(i)) attach_this(i,A_FRAMEWORK,chrom);
	        else unassign_this(i,A_CHANGED);
	  }
      }

    /* assign new anchors */
    for (j=0; j<num_loci; j++) {
	if (assigned(locus[j]) && assignment_chrom(locus[j])!=chrom)
	  unplace_this(locus[j],NO_CHROM,M_UNKNOWN,FALSE);
	sprintf(ps, ANCHOR_ISNOW, loc2str(locus[j]), chrom2str(chrom)); pr();
	assignment[locus[j]]->chromosome= chrom;
	assignment[locus[j]]->status=     A_ANCHOR;
	assignment[locus[j]]->LODscore=   NO_LOD;
	assignment[locus[j]]->theta=      NO_THETA;
	assignment[locus[j]]->linked_to=  NO_LOCUS;
	assignment[locus[j]]->modified=   TRUE;
    }
}


#define DOASS_NOANCHS "no anchor loci have yet been assigned to chromosomes\n" 
#define DOASS_NOCHROMS "no chromosomes have yet been defined\n"

void do_assignments(locus,num_loci,lod1,unlinked_lod1,theta1,
		    lod2,unlinked_lod2,theta2,haplo)
int *locus, num_loci;  /* is_assignable must have been verified */
real lod1, unlinked_lod1, theta1, lod2, unlinked_lod2, theta2;
bool haplo;
{
    int i, state, num_chroms;
    bool assigned_any;
    int **anchor=NULL, *count=NULL;
    real lod, unlinked_lod, theta;

    run {
	if ((num_chroms=chromosome->num_maps)==0) error(DOASS_NOCHROMS);
	matrix(anchor,num_chroms,raw.num_markers,int);
	array(count,num_chroms,int); for (i=0; i<num_chroms; i++) count[i]=0;
	if (!get_anchors(anchor,count,locus,num_loci,haplo))
	  error(DOASS_NOANCHS);

	lod=lod1; unlinked_lod=unlinked_lod1; theta=theta1;
	assigned_any=TRUE; 
	/* when we start loop with assigned_any=FALSE, we switch lods */
	do {
	    if (!assigned_any) {
		lod=lod2; unlinked_lod=unlinked_lod2; theta=theta2;
	    }
	    assigned_any=FALSE;
	    for (i=0; i<num_loci; i++) {
		if (locus[i]!=NO_LOCUS && 
		    do_assignment(locus[i],lod,unlinked_lod,theta,
				  anchor,count,num_chroms)) {
		    /* msg is printed */
		    locus[i]=NO_LOCUS;
		    assigned_any=TRUE;
		}
	    }
	} while (assigned_any || lod!=lod2);
	
	/* if a locus made it to here, its unlinked */
	for (i=0; i<num_loci; i++) if (locus[i]!=NO_LOCUS) {
	    state=assignment_state(i);
	    if (state!=A_FRAMEWORK && state!=A_ATTACH)
	      unassign_this(locus[i],A_UNLINKED);
	    else attach_this(locus[i],state,assignment_chrom(locus[i]));
	}
    } on_exit {
	unmatrix(anchor,num_chroms,int);
	unarray(count,int);
    }
}



#define ASS_MULTI  "%s- conflicting data (%s LOD %.2lf, %s LOD %.2lf)\n"

bool do_assignment(locus,lodbound,minlodbound,thetabound,
		   anchor,count,num_groups)
int locus;
real lodbound, minlodbound, thetabound;
int **anchor;   /* [num_groups][0..count[this_group]-1] */
int *count;     /* [num_groups] */
int num_groups;
{
    int assigned_locus, j, k, state;
    int pos_locus=0, this_chrom=0, pos_chrom, alt_chrom, pos_group=0;
    real lod, theta, this_lod, pos_lod, alt_lod, this_theta=0., pos_theta=0.;
    bool this_group;
    char *temp= get_temp_string();

    /* this actually does a bit of i/o (well, one line per locus) by
       virtue of calling assign_this(), which is gabby */ 

    /* this has slightly special behavior for loci with A_FRAMEWORK or
       A_ATTACH, in that it will not unlink them if there is no evidence for
       assignment (it will try to RELINK them however) THIS NEEDS WORK */

    if (sex_specific) send(CRASH);

    pos_chrom= alt_chrom= NO_CHROM;
    pos_lod= alt_lod= 0.0;
	    
    for (j=0; j<num_groups; j++) {
	/* if (do_changed_only && !changed[j]) continue; */
	this_group=FALSE; this_lod= 0.0;
	for (k=0; k<count[j]; k++) {
	    assigned_locus= anchor[j][k];
	    get_two_pt(locus,assigned_locus,&lod,&theta);
	    if (lod>=minlodbound && theta<=thetabound) {
		this_group= TRUE;
		this_chrom= j;
		if (lod>this_lod) {
		    this_lod= lod;
		    this_theta= theta;
		}
	    }
	} 
	if (this_group) { /* matched one in this group */
	    if (pos_chrom==NO_CHROM || 
		   (pos_chrom==this_chrom && this_lod>pos_lod)) {
		/* we've matched no other chroms in this or other groups */
		pos_chrom= this_chrom; pos_locus= assigned_locus; pos_group= j;
		pos_lod= this_lod;     pos_theta= this_theta;

	    } else if (this_lod>pos_lod) {
		/* we've matched another chrom, but this one is best */
		if (pos_lod>alt_lod) /* if ex-best chr is now second best */
		  { alt_chrom= pos_chrom; alt_lod= pos_lod; }
		pos_chrom= this_chrom; pos_locus= assigned_locus; pos_group= j;
		pos_lod= this_lod;     pos_theta= this_theta;
		
	    } else if (alt_chrom!=NO_CHROM && this_lod>alt_lod) {
		/* we've matched another, better, chrom, this is 2nd best */
		alt_chrom= this_chrom; alt_lod= this_lod;
	    }
	}
    }
    
    if (pos_lod>=lodbound && alt_chrom!=NO_CHROM && 
	((state=assignment_state(locus))!=A_FRAMEWORK && state!=A_ATTACH)) {
	/* really linked, but there is evidence of linkage to other chroms */
	sprintf(temp, ASS_MULTI, loc2str(locus),
            chrom2str(pos_chrom), pos_lod,
            chrom2str(alt_chrom), alt_lod);
	assign_this(locus,A_PROBLEM,NO_CHROM,NO_LOD,NO_THETA,NO_LOCUS,temp);
	return(TRUE);

    } else if (pos_lod>=lodbound) {
	/* really linked, only to this chrom */
	assign_this(locus,A_ASSIGNED,pos_chrom,pos_lod,pos_theta,pos_locus,"");
	anchor[pos_group][(count[pos_group])++]=locus;
	return(TRUE);

    } else return(FALSE);
}
    


bool get_anchors(anchor,count,locus,num_loci,haplo) /* internal use only */
int **anchor, *count; /* side-effected if non-null */
int *locus, num_loci; /* markers to assign now */
bool haplo;   /* take only haplo_firsts */
{
    int i, j, state, chrom, found, n=0;

    for (i=0; i<raw.num_markers; i++) if (assigned(i)) {
	state=assignment_state(i); chrom=assignment_chrom(i);
	if ((state==A_ANCHOR || state==A_BORDERLINE || state==A_ASSIGNED) &&
	      (!haplo || haplo_first[i]==NO_CHROM || haplo_first[i]==i)) {
	    for (found=FALSE, j=0; j<num_loci; j++) 
	      if (locus[j]==i) { found=TRUE; break; }
	    if (!found)	{ anchor[chrom][(count[chrom])++]=i; n++; }
	}
    }
    return(n);
}



/**************** Placements ****************/

bool placed_locus(locus)
int locus;
{ force_haplo_sanity(&locus,FALSE); return(placement[locus]->status>0); }

bool placement_state(locus)
int locus;
{ force_haplo_sanity(&locus,FALSE); return(placement[locus]->status); }


#define PLACE_ISFRAME \
  " %4d  %-9s - framework marker on chromosome %s, can't place\n"
#define PLACE_NOTASS  \
  " %4d  %-9s - not assigned to a chromosome, can't place\n"
#define PLACE_NOTCHR  \
  " %4d  %-9s - not assigned to chromosome %s, can't place\n"
#define PLACE_INSANE  \
  " %4d  %-9s - not the name of its haplotype-group, can't place\n"

#define PLACE_PROBLEM \
  " %4d  %-9s - no locations allowed, can't place\n"
#define PLACE_UNKNOWN \
  " %4d  %-9s - unplaced"
#define PLACE_FRAMEWORK \
  " %4d  %-9s - framework marker"
#define PLACE_TOOMANY \
  " %4d  %-9s - too many possible locations, can't place\n"

#define PLACE_OFFEND  \
  " %4d  %-9s - placed %soff %send of framework, like%s%.2f\n"
#define PLACE_UNIQUE  \
  " %4d  %-9s - placed %sin one interval, like%s%.2f\n"
#define PLACE_ZERO    \
  " %4d  %-9s - placed %sat zero distance, like%s%.2f\n"
#define PLACE_REGION  \
  " %4d  %-9s - placed %sin %d intervals, like%s%.2f 2nd%s%.2f\n"
/*234  123456789 - placed *with errors* in 24 intervals, like=-9.24, 2nd=-5.*/

bool is_placeable(locus,chrom)
int locus, chrom;
{
    /* unplace_this() helps to clean up the placement structure, when this
       is called by place or place_together */

    if (!force_haplo_sanity(&locus,FALSE)) {
	sprintf(ps, PLACE_INSANE, locus + 1, locname(locus, TRUE));
	pr(); return(FALSE);

    } else if (framework_marker(locus)) {
	sprintf(ps, PLACE_ISFRAME, locus + 1, locname(locus, TRUE),
            chrom2str(assignment_chrom(locus)));	pr(); 
	/* unplace_this(locus,NO_CHROM,M_UNKNOWN,FALSE); WHY? */
	return(FALSE);

    } else if (!assigned(locus)) {
	sprintf(ps, PLACE_NOTASS, locus + 1, locname(locus, TRUE)); pr();
	/* unplace_this(locus,NO_CHROM,M_UNKNOWN,FALSE); */
	return(FALSE);

    } else if (chrom!=ANY_CHROM && !assigned_to(locus,chrom)) {
	/* not used yet */
	sprintf(ps, PLACE_NOTCHR, locus + 1, locname(locus, TRUE), chrom2str(chrom)); pr();
	/* unplace_this(locus,NO_CHROM,M_UNKNOWN,FALSE); */
	return(FALSE);

    } else return(TRUE);
}


int place_this(locus,chrom,place,like_thresh,one_err_thresh,net_err_thresh,
	       excluded)
int locus, chrom;
PLACE **place;   /* better be ralative to chrom's fw order */
real like_thresh, one_err_thresh, net_err_thresh;
bool *excluded;  /* used as a temp for npt_exclusions() */
{
    int j, n, num, good;
    real second_best, support, error_lod;
    bool zero, off_end, is_single_error;
    char *sgn, *sgn2, *witherr;

    num= chromosome->map_list[chrom]->num_loci;

    good= npt_exclusions(place,num,like_thresh,one_err_thresh,net_err_thresh,
			 excluded,&zero,&off_end,&error_lod,&is_single_error,
			 &second_best,&support);

    placement[locus]->modified= TRUE;
    placement[locus]->chromosome=chrom;
    placement[locus]->threshold= like_thresh;
    placement[locus]->error_lod= error_lod; /* maybe NO_ERRORS */
    placement[locus]->single_error= is_single_error;

    n=0;
    for (j=0; j<=num; j++) if (!excluded[j]) {
	if (place[j]->like>PLACEMENT_THRESHOLD) {
	    if (n<MAX_STORED_PLACES) {
		placement[locus]->interval[n]= j;
		placement[locus]->like_ratio[n]= place[j]->like;
		placement[locus]->distance[n]=   place[j]->dist;
		placement[locus]->num_intervals= ++n;
	    } else { /* we are full */
		placement[locus]->num_intervals= 0;
		placement[locus]->threshold= 0.0;
		placement[locus]->status= M_PROBLEM;
		sprintf(ps, PLACE_TOOMANY, locus + 1, locname(locus, TRUE)); pr();
		return(0);
	    }
	}
    }

    sgn=ptr_to("=");
    if (support==NO_LIKE) { sgn=ptr_to("<"); support=-5.0; }
    witherr=ptr_to("");
    if (error_lod!=NO_ERRORS) witherr=ptr_to("*with errors* "); 

    if (off_end) {
	placement[locus]->status= (error_lod!=NO_ERRORS ? M_ERROR:M_OFFEND);
	sprintf(ps, PLACE_OFFEND, locus + 1, locname(locus, TRUE), witherr,
            (off_end==2 ? "*either* ":""), sgn, support);
	pr();
	/* print the errors? */

    } else if (good==1) {
	if (zero) {
	    placement[locus]->status=(error_lod!=NO_ERRORS ? M_ERROR:M_ZERO);
	    sprintf(ps, PLACE_ZERO, locus + 1, locname(locus, TRUE), witherr,
                sgn, support); pr();

	} else {
	    placement[locus]->status=(error_lod!=NO_ERRORS ? M_ERROR:M_UNIQUE);
	    sprintf(ps, PLACE_UNIQUE, locus + 1, locname(locus, TRUE), witherr,
                sgn, support); pr();
	}

    } else { /* good>1 */
	sgn2=ptr_to("=");
	if (second_best==NO_LIKE) { sgn2=ptr_to("<"); second_best=-5.0; }
	placement[locus]->status= (error_lod!=NO_ERRORS ? M_ERROR:M_REGION);
	sprintf(ps, PLACE_REGION, locus + 1, locname(locus, TRUE), witherr, good,
            sgn, support, sgn2, second_best); 
	pr();
	/* print possible errors? */
    }

    return(good);
}


void unplace_this(locus,chrom,status,verbose)
int locus, chrom, status;
bool verbose;
{
    placement[locus]->modified= TRUE;
    placement[locus]->num_intervals= 0;
    placement[locus]->threshold= 0.0;
    placement[locus]->error_lod=  NO_ERRORS;
    placement[locus]->single_error= TRUE;
    placement[locus]->chromosome= chrom;
    placement[locus]->status= status;
    if (verbose) {
	if (status==M_PROBLEM)
	  { sprintf(ps, PLACE_PROBLEM, locus + 1, locname(locus, TRUE)); pr(); }
	else if (status==M_FRAMEWORK)
	  { sprintf(ps, PLACE_FRAMEWORK, locus + 1, locname(locus, TRUE)); pr(); }
	else
	  { sprintf(ps, PLACE_UNKNOWN, locus + 1, locname(locus, TRUE)); pr(); }
    }
}


/* These have no error traps - only use them if valid data are stored! */
int 
best_placement (int locus)
{ 
    int i; 
    
    if (placement[locus]->num_intervals==0) send(CRASH);
    for (i=0; i<placement[locus]->num_intervals; i++) 
      if (placement[locus]->like_ratio[i]==0.0) break;
    if (i==placement[locus]->num_intervals) send(CRASH);
    return(i);
}


int second_best_placement(locus,like) /* only ok for M_REGION */
int locus;
real *like;
{ 
    int i, j; /* need to add many error traps */
    real best;
    
    for (i=0; i<placement[locus]->num_intervals; i++)
      if (placement[locus]->like_ratio[i]==0.0) break;
    best= -VERY_BIG;
    for (j=0; j<placement[locus]->num_intervals; j++)
      if (j!=i && placement[locus]->like_ratio[j]>best) 
	{ best=placement[locus]->like_ratio[j]; ; }
    if (like!=NULL) *like=best;
    if (best== -VERY_BIG) return(-1);
    return(j);
}



/************ Load/Save/Reset ************/


void 
bash_mapping_data (
    int *changed,
    int num_changed /* changed markers */
)
{
    int i, j, n;
    MAP *map;
    bool bash[MAX_CHROMOSOMES];
    for (j=0; j<chromosome->num_maps; j++) bash[j]=FALSE;

    for (i=0; i<num_changed; i++) {
	n=changed[i];
	sprintf(ps, "%-8s - unassigned and unplaced (changed)\n",
            locname(n,FALSE)); pr();
    }
    for (i=0; i<num_changed; i++) {
	n=changed[i];
	if (assignment[n]->status==A_ANCHOR) {
	    sprintf(ps, "warning: %s WAS an anchor locus for chromosome %s\n",
                locname(n,FALSE), chrom2str(assignment[n]->chromosome)); pr();
	}
    }
    for (i=0; i<num_changed; i++) {
	n=changed[i];
	if (assignment[n]->status!=A_UNKNOWN)
	  assignment[n]->status=   A_CHANGED;
	assignment[n]->chromosome= NO_CHROM; 
	assignment[n]->linked_to=  NO_LOCUS;
	assignment[n]->LODscore=   NO_LOD;
	assignment[n]->theta=      NO_THETA;
	assignment[n]->modified=   FALSE;

	if (placement[n]->status==M_FRAMEWORK)
	  bash[assignment[n]->chromosome]=TRUE;
	placement[n]->status= M_UNKNOWN; 
	placement[n]->num_intervals= 0;
	placement[n]->threshold= 0.0;
	placement[n]->error_lod= NO_ERRORS;
	placement[n]->single_error= TRUE;
	placement[n]->chromosome= NO_CHROM;
	placement[n]->modified= FALSE;
    }

    for (j=0; j<chromosome->num_maps; j++) 
      if (bash[j] && chromosome->map_list[j]->num_loci>0) {
	  sprintf("chromosome %s framework cleared, all loci unplaced\n",
              chrom2str(j)); pr();
	map=get_map_to_bash(chromosome);
	clean_map(map); /* num_loci -> 0 */
	set_chrom_frame(j,map);
    }
}


void 
allocate_mapping_data (int num_markers)
{
    int locus;

    chromosome= 
      allocate_map_list(MAX_CHROMOSOMES,MAX_CHROM_LOCI,UNSORTED,NULL);
    parray(assignment,num_markers,ASSIGNMENT);
    parray(placement, num_markers,PLACEMENT);

    for (locus=0; locus<num_markers; locus++) {
	assignment[locus]->status=    A_UNKNOWN;
	assignment[locus]->chromosome=NO_CHROM;
	assignment[locus]->linked_to= NO_LOCUS; 
	assignment[locus]->LODscore=  NO_LOD;
	assignment[locus]->theta=     NO_THETA;
	assignment[locus]->modified= FALSE;

	placement[locus]->status= M_UNKNOWN;
	placement[locus]->num_intervals= 0;
	placement[locus]->threshold= 0.0;
	placement[locus]->error_lod= NO_ERRORS;
	placement[locus]->single_error= TRUE;
	placement[locus]->chromosome= NO_CHROM;
	placement[locus]->modified= FALSE;
    }
}


void 
free_mapping_data (int num_markers)
{
    free_map_list(chromosome);
    unparray(assignment,num_markers,ASSIGNMENT);
    unparray(placement, num_markers,PLACEMENT);
}


void write_mapping_data(fp)
FILE *fp;
{
    int locus, i, j;

    sprintf(ps, "*Chromosomes: %d\n", chromosome->num_maps); fpr(fp);
    for (i=0; i < chromosome->num_maps; i++) {
        write_map(fp, chromosome->map_list[i]);
    }
   
    fprint(fp,"*Assignments and Placements:\n");
    for (locus=0; locus< raw.num_markers; locus++) {
        sprintf(ps, "*%-8s %2d", raw.locus_name[locus], assignment[locus]->status); 
	fpr(fp);
	if (assignment[locus]->status==M_UNKNOWN) { fnl(fp); continue; }

	sprintf(ps, " %d %d %f %f | %d %f %f %d %d \n",
            assignment[locus]->chromosome, assignment[locus]->linked_to,
            assignment[locus]->LODscore, assignment[locus]->theta,
            placement[locus]->status, placement[locus]->threshold,
            placement[locus]->error_lod, placement[locus]->single_error,
            placement[locus]->num_intervals); fpr(fp);

	if (placement[locus]->num_intervals==0) continue;

	for (j=0; j < placement[locus]->num_intervals; j++) {
	    sprintf(ps, " %d %f %f",
                placement[locus]->interval[j], placement[locus]->like_ratio[j],
                placement[locus]->distance[j]); fpr(fp);
	}
	fnl(fp);
    }
}


void read_mapping_data(fp)
FILE *fp;
{
    int locus, num, i, j, num_chroms;
    real rnum;
    char word[TOKLEN+1], temp_locus_name[NAME_LEN+2];
    MAP  *map;

    /* chromosomes */
    fgetln(fp);
    stoken(&ln,sREQUIRED,word);
    if (!streq(word,"*Chromosomes:") || !itoken(&ln,iREQUIRED,&num_chroms)) {
        baddata("error finding *Chromosomes:");
    }
    for (i=0; i < num_chroms; i++) {
        map= get_map_to_bash(chromosome);
	read_map(fp,map);
	insert_map_into_list(chromosome,&map);
    }
    if (chromosome->num_maps != num_chroms) {
        baddata("listed number of chromosomes and actual number do not agree");
    }

    fgetln(fp);
    if (!streq(ln,"*Assignments and Placements:"))
        baddata("error finding *Assignments and Placements:");
    for (locus=0; locus < raw.num_markers; locus++) {
        getdataln(fp);
    
	if (!nstoken(&ln,sREQUIRED,temp_locus_name,NAME_LEN+1) || 
	    temp_locus_name[0]!='*' || len(temp_locus_name)<2)
	  baddata("expected *name");
	else if (!streq(raw.locus_name[locus],&temp_locus_name[1]))
	  baddata("locus names don't match");

	/* assignment info */
	itoken(&ln,iREQUIRED,&num);
	assignment[locus]->status=num;
	if (num==A_UNKNOWN) continue; /* next locus */

	if (sscanf(ln,"%d %d %lf %lf %s %d %lf %lf %d %d",
		   &assignment[locus]->chromosome,
		   &assignment[locus]->linked_to,
		   &assignment[locus]->LODscore,
		   &assignment[locus]->theta,    
		   word,
		   &placement[locus]->status,
		   &placement[locus]->threshold,
		   &placement[locus]->error_lod,
		   &placement[locus]->single_error,
		   &placement[locus]->num_intervals)!=10 || !streq(word,"|"))
	  baddata("unable to parse assignment/placement line");

	placement[locus]->chromosome= assignment[locus]->chromosome;
	if (placement[locus]->num_intervals==0) continue; /* next locus */
        
	getdataln(fp); /* get next line */
	for (j=0; j < placement[locus]->num_intervals; j++) {
	    itoken(&ln,iREQUIRED,&num); 
	    placement[locus]->interval[j]=num;
	    rtoken(&ln,rREQUIRED,&rnum);
	    placement[locus]->like_ratio[j]=rnum;
	    rtoken(&ln,rREQUIRED,&rnum);
	    placement[locus]->distance[j]=rnum;
	}
    }
}
