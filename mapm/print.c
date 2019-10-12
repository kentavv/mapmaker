/******************************************************************************

 #####   #####      #    #    #   #####           ####
 #    #  #    #     #    ##   #     #            #    #
 #    #  #    #     #    # #  #     #            #
 #####   #####      #    #  # #     #     ###    #
 #       #   #      #    #   ##     #     ###    #    #
 #       #    #     #    #    #     #     ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#define INC_MEM
#define INC_IO
#define INC_MATH
#define INC_STR
#define INC_MSG
#include "mapm.h"

#define empty_loc_str loc2str(-1)
#define empty_locs_str locs2str(-1,-1)
void print_err();

#define MAX_LOG_LIKE_COLUMN 66

void print_tiny_map(MAP *map, char *text  /* make sure its OK to use ps for this */, real like_base)
{
    int i, stop_column; 
    bool wrapped, had_unlinked=FALSE;

    stop_column= MAX_LOG_LIKE_COLUMN-len(empty_loc_str)-3;
    wrapped= FALSE;
    
    print(text); print("  ");
    for (i=0; i<map->num_loci; i++) {
	if (at_column()>stop_column) { 
	    if (!wrapped) {
		to_column(MAX_LOG_LIKE_COLUMN);
		if (like_base!=0.0)
		  sprintf(ps,"Like: %s",rsn(6.2,map->log_like-like_base));
		  else sprintf(ps, "Like: %s", rsn(8.2, map->log_like));
		pr();
	    }
	    nl(); space(len(text)+2); 
	    wrapped=TRUE;
	}
	if (!print_names) sprintf(ps, "%d", map->locus[i] + 1);
	    else strcpy(ps,raw.locus_name[map->locus[i]]);
	pr();
	if (i<map->num_loci-1 && map->rec_frac[i][0]==UNLINKED_CM) 
	  { print(" *"); had_unlinked=TRUE; }
	print(" ");
    }
    if (!wrapped) {
	if (!had_unlinked) print("  ");
	if (like_base!=0.0) 
	  sprintf(ps,"Like: %s",rsn(6.2,map->log_like-like_base)); 
	  else sprintf(ps,"Like: %s",rsn(6.2,map->log_like));
	pr();
    }
    nl();
}



void 
print_err (int indiv, int locus, MAP *map)
{
    sprintf(ps, "[#%d %c-%c-%c %.2lf]", indiv + 1,
            raw.data.f2.allele[map->locus[locus-1]][indiv],
            raw.data.f2.allele[map->locus[locus]]  [indiv],
            raw.data.f2.allele[map->locus[locus+1]][indiv],
            map->error_lod[locus][indiv]); 
    pad_to_len(ps,18); pr();
}


void 
print_long_map (MAP *map, char *title)
{ print_special_map(map,title,0,NULL); }


void 
print_special_map (
    MAP *map,
    char *title,  /* ps MAY be used for this! */
    int num_old,
    int *old_locus  /* can omit if num_args==0 - should do this with VARARGS */
)
{
    int  i, j, n, m;
    real length, male_cm, female_cm, male_length, female_length;
    real best, next, lod;
    int  bestj, nextj;
    char p1, p2, num[TOKLEN+1], my_title[TOKLEN+1];

    if (!map->sex_specific) { /**************** CEPH or F2 ****************/
	length= 0.0;
	/* don't use ps until the title is printed */
	sprintf(my_title, "%-26s  ", (!nullstr(title) ? title : "Map:"));
	print(my_title);
	if (map->error_lod!=NULL) { sprintf(ps, "Apriori"); pr(); }
	nl();

	sprintf(ps, "  Markers          Distance "); pr(); /* that's 28 chars */
	if (map->error_lod!=NULL) { sprintf(ps, "  Prob  Candidate Errors"); pr(); }
	nl();

	for (i=0; i<map->num_loci-1; ++i) {
	    p1=p2=' ';
	    if (num_old>0) for (p1='(', p2=')', j=0; j<num_old; j++) 
	      if (old_locus[j]==map->locus[i]) { p1=p2=' '; break; }
	    sprintf(num, "%c%d", p1, map->locus[i] + 1);
	    sprintf(ps, "%5s%c %-9s  ", num, p2, locname(map->locus[i], use_haplotypes));
	    pr();

	    sprintf(ps, " %s %2s", rf2str(map->rec_frac[i][0]),
                (units==CENTIMORGANS ? "cM":"rf")); pr();
	    length+= cm_dist(map->rec_frac[i][0]);
	    
	    if (map->error_lod!=NULL && i>0) { /* must be F2 data */
		sprintf(ps, "   %s%%  ", rsd(3.1, map->error_rate[i] * 100.0)); pr();

		if (!print_all_errors) { /* print the two worst and a count */
		    best= next= 0.0; bestj= nextj=0; n=0;
		    for (j=0; j<raw.data.f2.num_indivs; j++)
		      if ((lod=map->error_lod[i][j])>=error_lod_thresh) {
			  n++;
			  if (lod>best) 
			    { next=best; nextj=bestj; best=lod; bestj=j; }
			  else if (lod>next) { next=lod;  nextj=j; }
		      }
		    if (n==0) print("-");
		    if (n>=1) print_err(bestj,i,map);
		    if (n>=2) print_err(nextj,i,map);
		    if (n>=3) { sprintf(ps, "%d more", n - 2); pr(); }

		} else { /* print all errors, unsorted */
		    for (n=0, j=0; j<raw.data.f2.num_indivs; j++)
		      if ((lod=map->error_lod[i][j])>=error_lod_thresh) n++;
		    if (n==0) print("-");
		    else for (m=0, j=0; j<raw.data.f2.num_indivs; j++)
		      if ((lod=map->error_lod[i][j])>=error_lod_thresh) {
			  if (m>0 && m%2==0) /* nl + 28 (to end of cM) + 9 */
			    print("\n                                     "); 
			  print_err(j,i,map);
		    }
		}
	    }
	    nl();
	}
	i=map->num_loci-1;
	p1=p2=' ';
	if (num_old>0) for (p1='(', p2=')', j=0; j<num_old; j++) 
	  if (old_locus[j]==map->locus[i]) { p1=p2=' '; break; }
	sprintf(num, "%c%d", p1, map->locus[i] + 1);
	sprintf(ps, "%5s%c %-9s  ", num, p2, locname(map->locus[i], use_haplotypes));
	pr();

	print("----------\n");
	sprintf(ps,
            "                   %5s cM   %d markers   log-likelihood= %.2lf\n",
            rsd(5.1,length), map->num_loci, map->log_like); pr();

    } else { /**************** Sex-Specific CEPH ****************/
	/* needs fixing! */
	male_length= 0.0;
	female_length= 0.0;
	print(empty_locs_str);
	print("    MALE-MAP:             FEMALE-MAP:\n");
	for (i=0; i<map->num_loci-1; ++i) {
	    
	    male_cm= cm_dist(map->rec_frac[i][MALE]);
	    female_cm= cm_dist(map->rec_frac[i][FEMALE]);
	    sprintf(ps, "  %s  %s cM   %s %%     %s cM   %s %%\n",
                locs2str(map->locus[i],map->locus[i+1]),
                rsd(5.1,male_cm), rsd(5.1,map->rec_frac[i][MALE]*100.0),
                rsd(5.1,female_cm), rsd(5.1,map->rec_frac[i][FEMALE]*100.0));
            print(ps);
	    male_length += male_cm;
	    female_length += female_cm;
        }
	print(empty_locs_str);
	print("    --------               --------\n"); 
	print(empty_locs_str);
	sprintf(ps, "    %s cM               %s cM\n\n",
            rsd(5.1,male_length), rsd(5.1,female_length)); pr();
	sprintf(ps, "  log-likelihood = %lf\n", map->log_like);
	print(ps);
    }
}


void 
print_short_map (MAP *map, char *text)
{
    int i, col;

    if (len(text)<5) col=6; else col=len(text)+1;
    print(text); print(" "); to_column(col+2); 
    for (i=0; i<map->num_loci; i++) {
	print(loc2str(map->locus[i]));
	if (map->unlink==i+1) print(" * ");
	else if (!print_names) print("   "); else print(" ");
    }
    if (map->unlink==NONE_UNLINKED && print_names) print("  ");
    sprintf(ps, " log-likelihood: %s\n", rsn(7.2, map->log_like)); pr();

    if (!map->sex_specific) {
	print("map:"); to_column(col+3); 
	if (print_names) space((len(empty_loc_str)/2) + 2); else print(" ");
	for (i=0; i<map->num_loci-1; i++) {
	    print(rf2str(map->rec_frac[i][MALE]));
	    if (map->unlink==i+1 && print_names) print("  ");
	    if (print_names) space(len(empty_loc_str)-5+1); else print(" ");
	}
	nl();

    } else { /* sex_spec */
	print("male:"); to_column(col+3); 
	if (print_names) space(len(empty_loc_str)/2); else print(" ");
	for (i=0; i<map->num_loci-1; i++) {
	    print(rf2str(map->rec_frac[i][MALE]));
	    if (map->unlink==i+1 && print_names) print("  ");
	    if (print_names) space(len(empty_loc_str)-5+1); else print(" ");
	}
	nl();
	print("female:"); to_column(col+3); 
	if (print_names) space(len(empty_loc_str)/2); else print(" ");
	for (i=0; i<map->num_loci-1; i++) {
	    print(rf2str(map->rec_frac[i][FEMALE]));
	    if (map->unlink==i+1 && print_names) print("  ");
	    if (print_names) space(len(empty_loc_str)-5+1); else print(" ");
	}
	nl();
    }
}


void 
print_list (SAVED_LIST *list, int how_many)
{
    int i,num_to_print;
    char str[TOKLEN+1];
    real best;

    if (how_many==FULL_LIST || how_many>list->num_maps || how_many<0)
      num_to_print= list->num_maps;
    else num_to_print= how_many;

    best= -VERY_BIG;
    for (i=0; i<list->num_maps; i++) 
      best= rmaxf(list->map_list[i]->log_like,best);

    for (i=0; i<num_to_print; i++) {
	if (i!=0 || num_to_print!=1) sprintf(str, "%d:", i + 1);
	    else strcpy(str,"1:");
	pad_to_len(str,4);
        print_tiny_map(list->map_list[i],str,best);
    }
}


char *
rf2str (real rec_frac)
{
    real d;

    rrange(&rec_frac,0.0,0.499999);
    if (units==CENTIMORGANS) {
	d=(*mapfunction->rec_to_dist)(rec_frac); 
	if (d<0.0) d=0.0; /* maybe some rounding error */
	return(rsd(5.1,100.0*d));
    } else return(rsd(5.3,rec_frac));
}


char *
loc2str (int locus)
{
    char *str= get_temp_string(), haplo=' ';
    
    if (locus<0) {
	if (print_names) strcpy(str,"         "); else strcpy(str,"     ");
    } else {
	if (use_haplotypes && haplotyped(locus) && 
	    !haplotype_subordinate(locus)) haplo='+';
	if (print_names)
	    { sprintf(str, "%s%c", raw.locus_name[locus], haplo); pad_to_len(str, 9); }
	  else { sprintf(str, "%d%c", locus + 1, haplo); pad_to_len(str, 5); }
    }
    return(str);
}


char *
locname ( /* ragged */
    int locus,
    bool haplo_mark
)
{
    char *str= get_temp_string();
    bool haplo=FALSE;
    
    if (haplo_mark && haplotyped(locus)) haplo=TRUE;
    sprintf(str, "%s%s", raw.locus_name[locus], (haplo ? "+" : ""));
    return(str);
}


char *
locs2str (int locus1, int locus2)
{
    char *str, *haplo1=get_temp_string(), *haplo2=get_temp_string();
    int i;

    str= get_temp_string();
    if (locus1<0 || locus2<0) {
	str[0]='\0';
    } else {
	strcpy(haplo1,(use_haplotypes && haplotyped(locus1)) ? "+":"");
	strcpy(haplo2,(use_haplotypes && haplotyped(locus2)) ? "+":"");
	if (print_names) 
	  sprintf(str,"%s%s - %s%s",
		  raw.locus_name[locus1],haplo1,raw.locus_name[locus2],haplo2);
	  else sprintf(str,"%d%s - %d%s",locus1+1,haplo1,locus2+1,haplo2);
    }
    for (i=len(str); i< (print_names ? 25:13 ); i++) str[i]=' ';
    str[i]='\0';
    return(str);
}


char *
rag (char *str)
{
    int i=0;
    while (*str==' ') str++;
    while (str[i]!='\0') i++;
    i--;
    while (str[i]==' ') i--;
    str[++i]='\0';
    return(str);
}


void 
print_trys (
    SAVED_LIST **list,
    MAP *base_map,	   /* The map the tried markers were added to */
    bool **excluded,   /* [i][n]=FALSE if we tried try_marker[i] in excluded */
    int **new_marker,
    int num_tried,
    int first
)
{
    int i, j, k, q, width_ea, last;
    bool any_paired=FALSE;
    real best[8];
    char *line1= get_temp_string(), *line2= get_temp_string();
    MAP *map;

    if (num_tried>8) send(CRASH);
    if (print_names) width_ea=11; else width_ea=8;

    sprintf(line2, empty_loc_str); sprintf(line1, empty_loc_str); 
    strcat(line2,"   "); strcat(line1,"   ");
    
    for (i=0; i<num_tried; i++) {
	if (list[i]->num_maps==0) send(CRASH);
        else { map=get_best_map(list[i]); best[i]=map->log_like; }
	if (new_marker[i+first][1]!=NO_LOCUS) { /* then paired */
	    last=1; while(new_marker[i+first][last+1]!=NO_LOCUS) last++;
	    sprintf(ps, "[%s", rag(loc2str(new_marker[i + first][0])));
	    strcat(line1,pad_to_len(ps,width_ea));
	    sprintf(ps, "%s%s]", (last > 1 ? "-" : " "),
                rag(loc2str(new_marker[i+first][last])));
	    strcat(line2,pad_to_len(ps,width_ea));
	    any_paired=TRUE;
	} else {
	    strcat(line1,pad_to_len(empty_loc_str,width_ea));
	    strcat(line2,pad_to_len(loc2str(new_marker[i+first][0]),width_ea));
	}
    }
    if (any_paired) { print(line1); nl(); }
    print(line2); nl();

    print(empty_loc_str); print("  "); 
    for (i=0; i<num_tried*width_ea-1; i++) print("-");
    print("\n");

    for (i=0; i<=base_map->num_loci; i++) { /* i=interval */
	print(empty_loc_str); print(" |");
	for (j=0; j<num_tried; j++)  /* j=try_marker */
	  if (!excluded[j][i]) {
	      for (k=q=0; k<i; k++) if (!excluded[j][k]) q++;
	      if (list[j]->map_list[q]->log_like==ZERO_LIKE) sprintf(ps, "  zero");
	      else sprintf(ps, "%s", rsd(6.2, list[j]->map_list[q]->log_like - best[j]));
	      print(pad_to_len(ps,width_ea-(j==num_tried-1 ? 1:0)));
	  } else {
	      strcpy(ps,"   x");
	      print(pad_to_len(ps,width_ea-(j==num_tried-1 ? 1:0)));
	  }
	print("|\n");
	if (i<base_map->num_loci) {
	    print(loc2str(base_map->locus[i])); print(" |");
	    for (j=0; j<num_tried*width_ea-1; j++) print(" ");
	    print("|\n");
	}
    }

    space(len(empty_loc_str)); print(" |");
    for (i=0; i<num_tried*width_ea-1; i++) print("-");
    print("|\n");

    print("INF"); space(len(empty_loc_str)-3); print(" |");
    i=base_map->num_loci+1;
    for (j=0; j<num_tried; j++) {
	for (k=q=0; k<i; k++) if (!excluded[j][k]) q++;
        sprintf(ps, "%s", rsd(6.2, list[j]->map_list[q]->log_like - best[j]));
	print(pad_to_len(ps,width_ea-(j==num_tried-1 ? 1:0)));
    } 
    print("|\n");

    print(empty_loc_str); print("  "); 
    for (i=0; i<num_tried*width_ea-1; i++) print("-");
    print("\n");

    print("BEST"); space(len(empty_loc_str)+1-4); 
    for (j=0; j<num_tried; j++) {
        sprintf(ps, "%s", rsd(7.2, best[j]));
	print(pad_to_len(ps,width_ea));
    }
    nl();
}


void 
print_haplo_summary (int *locus, int num_loci)
{
    int i, j, any=FALSE;
    
    for (i=0; i<num_loci; i++) {
	if (haplotyped(locus[i])) {
	    if (!any) { print("Haplotype groups:\n"); any=TRUE; }
	    sprintf(ps, " %4d  %-9s = [", haplo_first[locus[i]] + 1,
                locname(haplo_first[locus[i]],TRUE)); pr();
	    for (j=haplo_first[locus[i]]; j!=NO_CHROM; j=haplo_next[j]) {
		if (j!=haplo_first[locus[i]]) print(" ");
		if (print_names) sprintf(ps, "%s", locname(j, FALSE));
		else sprintf(ps, "%d", j + 1);
		pr();
	    }
	    print("]\n");
	}
    }
    if (!any) print("No Haplotype Groups\n");
}



char *
region2str (int locus, char **errs)
{
    *errs= ptr_to("0.0");
    return("abc1-xyz999");
}


char *
genetics2str (int locus, bool haplo)
{
    int type, n_infs, n_dom_obs, n_het_obs;
    char *retoin= get_temp_string();

    if (raw.data_type!=F2) { strcpy(retoin,"ceph"); return(retoin); }

    f2_genotype(locus,haplo,observations);
    n_infs= f2_count_infs(&n_dom_obs,&n_het_obs,observations);
    type=raw.data.f2.cross_type;

    if (type!=F2_INTERCROSS && type!=F3_SELF) 
      { sprintf(retoin, "%4d", n_infs); return(retoin); }
    
    if (n_dom_obs==0)      sprintf(retoin, "%4d codom", n_infs);
    else if (n_het_obs==0) sprintf(retoin, "%4d +/-  ", n_infs);
    else                   sprintf(retoin, "%4d mixed", n_infs);
    return(retoin);
}
    

#define LIST_HEAD1 \
"                            Error           Linkage   Haplotype"
#define LIST_HEAD2 \
" Num  Name       Genotypes  Prob   Chrom    Group     Group     Class    New?"
/*234  12345678+ 1234567890 12.45%  12345678 group123  12345678+ 12345678 new */
#define LIST_LINE \
   "%4d  %-9s %10s %5.2lf%%  %-8s %-8s  %-9s"
/*  num  name typ  err       chr  lg    hap */


void 
print_locus_summary (int *locus, int n_loci, bool haplo)
{
    int i, g;
    char *chrom, *hap, lg[TOKLEN+1];

    for (i=0; i<n_loci; i++)
      if (modified[locus[i]]) { break; }
    
    print(LIST_HEAD1); nl();
    print(LIST_HEAD2); nl();
    
    for (i=0; i<n_loci; i++) {
	if (!assigned(locus[i]))chrom= ptr_to("-");
	else chrom= chrom2str(assignment_chrom(locus[i]));
	
	if (!haplotyped(locus[i])) hap= ptr_to("-");
	else hap= locname(haplo_first[locus[i]],TRUE);
	
	g= my_group[locus[i]];
	if (g==NO_GROUP) strcpy(lg,"-");
	else if (g==num_groups-1) strcpy(lg,"unlinked");
	else sprintf(lg,"group%d",g+1);
	
	sprintf(ps, LIST_LINE, locus[i] + 1, locname(locus[i], haplo),
            genetics2str(locus[i],haplo),
	   error_rate[locus[i]]*100.0, chrom, lg, hap); pr(); 
	if (class[locus[i]]==NO_CLASS) sprintf(ps, " -       "); 
	else sprintf(ps, "%-8s ", class_name[class[locus[i]]]);
	pr();
	if (modified[locus[i]]) print("new"); else print(" - ");
	nl();
    }    
}

#define MAPPING_HEAD1 \
"                                                      2nd    Left "
#define MAPPING_HEAD2 \
" Num  Name      Assignment Chrom     LOD   Mapping    Like   Locus     Errors"
/*234  12345678+ 123456789  12345678 12.45  123456789 -23.45  12345678+  1.34* */
#define MAPPING_LINE \
   "%4d  %-9s %-9s  %-8s %-5s  %-9s %-6s  %-9s %s%s\n"
/*  num  name ass   chr  lod   map  like loc  err */


void 
print_mapping_summary (int *locus, int n_loci, bool haplo)
{
    int i, j, k, pos, state;
    real val;
    char *chrom, *lod, *assign;
    char *mapping, *left, *like, *errs, *star;

    print(MAPPING_HEAD1); nl();
    print(MAPPING_HEAD2); nl();
    
    for (i=0; i<n_loci; i++) {
	chrom=  ptr_to("-");
	lod=    ptr_to("  -  ");
//	char *theta=  ptr_to("   - ");  /* theta and linked are not printed */
//	char *linked= ptr_to("-");
	state= assignment_state(locus[i]);
	
	if (state==A_PROBLEM) {
	    assign= ptr_to("conflict!");
	} else if (state==A_CHANGED) {
	    assign= ptr_to("changed");
	} else if (state==A_UNLINKED) {
	    assign= ptr_to("unlinked");
	} else if (state==A_UNKNOWN) { 
	    assign= ptr_to("-");
	} else if (state==A_BORDERLINE) {
	    assign= ptr_to("low-lod");
	    chrom=  chrom2str(assignment_chrom(locus[i])); 
	    lod=    rsd(5.2,assignment_lod(locus[i]));
//	    theta=  rsd(5.1,assignment_theta(locus[i]));
//	    linked= loc2str(assignment_locus(locus[i]));
	} else if (state==A_ASSIGNED) {
	    assign= ptr_to("linked");
	    chrom=  chrom2str(assignment_chrom(locus[i])); 
	    lod=    rsd(5.2,assignment_lod(locus[i]));
//	    theta=  rsd(5.1,assignment_theta(locus[i]));
//	    linked= loc2str(assignment_locus(locus[i]));
	} else if (state==A_ATTACH) {
	    assign= ptr_to("attached");
	    chrom=  chrom2str(assignment_chrom(locus[i]));
	} else if (state==A_FRAMEWORK) {
	    assign= ptr_to("framework");
	    chrom=  chrom2str(assignment_chrom(locus[i]));
	} else if (state==A_ANCHOR) {
	    assign= ptr_to("anchor");
	    chrom=  chrom2str(assignment_chrom(locus[i]));
	} else {
	    assign= ptr_to("???");
	}
	
	like= ptr_to("  -");
	left= ptr_to("-");
	errs= ptr_to(" - ");
	star= ptr_to("");
	state=placement_state(locus[i]);
	
	if (state==M_PROBLEM) {
	    mapping= ptr_to("problem!");
	} else if (state==M_FRAMEWORK) {
	    mapping= ptr_to("framework");
	} else if (state==M_UNKNOWN) {
	    mapping= ptr_to("-");
	} else if (state==M_ERROR) {
	    mapping= ptr_to("errors?");
	} else if (state==M_OFFEND) {
	    mapping= ptr_to("off-end");
	} else if (state==M_REGION) {
	    mapping= ptr_to("region");
	} else if (state==M_UNIQUE) {
	    mapping= ptr_to("unique");
	} else if (state==M_ZERO) {
	    mapping= ptr_to("zero");
	} else {
	    mapping= ptr_to("???");
	}

	if (placed_locus(locus[i])) {
	    j=best_placement(locus[i]);
	    k= assignment_chrom(locus[i]);
	    pos= placement[locus[i]]->interval[j];
	    if (pos==0) left=ptr_to("leftmost");
	    else left=loc2str(chromosome->map_list[k]->locus[pos-1]);
	    if (placement[locus[i]]->error_lod!=NO_ERRORS) {
		errs= rsd(5.2,placement[locus[i]]->error_lod);
		if (!placement[locus[i]]->single_error) star=ptr_to("*");
	    }
	    if (second_best_placement(locus[i],&val)>=0) like=rsd(6.2,val);
	}
	
	sprintf(ps, MAPPING_LINE, locus[i] + 1, locname(locus[i], haplo),
            assign, chrom, lod, mapping, like, left, errs, star);
	pr();
    }
}



#ifdef OBSOLETE
void 
print_placements (
    int *order,
    int num_order,
    int *locus,
    int num_loci,  /* unplaced, maybe NO_LOCUS */
    bool **excluded       /* first index is that into locus[] */
)
{
    int num_left, num_across, num_done, i, j;

    num_left= num_order+1;
    num_done=0;
    do {
	num_across=(num_left > 17) ? 18 : num_left;
	print("    ");
	for (j=num_done; j<num_across+num_done-1; j++) {
	    sprintf(ps," %3d",order[j]+1);  pr();
	}
	nl();
	
	for (i=0; i<num_loci; i++) if (locus[i]!=NO_LOCUS) {
	    sprintf(ps,"%3d ",locus[i]+1);  pr();
	    for (j=num_done; j<num_across+num_done; j++) {
		switch(excluded[i][j]) {
		    case TRUE:  print("."); break; /* kind of rely on TRUE=1 */
		    case FALSE: print("*"); break;
		    case MAYBE: print("^"); break;
		}
		if (j!=num_across+num_done-1) print("...");
	    }
	    nl();
	}
	nl();
	num_left-= 18;
	num_done+= 18;
    } while(num_left > 0);
}	
#endif


/*345678901234567890123456789012345678901234567890123456789012345678901234567*/
/* 15 intervals, or 14 middle intervals + both ends =16

12345678 num_done=0  locus fmts: xx1x_  x12x_  x123_  1234_
   xxxx 1234 1234 1234   1  1234  34  1234  234 _234 _234 _234 _234 _234 _234
      xxxx:__3 :123 :_23 :_23 :_23 :_23 :_23 :_23 :_23 :_23 :_23 :_23 :_23 :...
 234 1 .*.:..*.:..*.:.*..:....:....:....:....:....:....:....:....:....:....:...
1234 12...:....:.**.:....:....:....:....:....:....:....:....:....:....:....:...

12345 num_done>0
     1234 1234 _234 _234 _234 _234 _234 _234 _234 _234 _234 _234 _234 _234 1234
   xxxx:__3 :123 :_23 :_23 :_23 :_23 :_23 :_23 :_23 :_23 :_23 :_23 :_23 :.23.:
 234 1 :..*.:..*.:.*..:....:....:....:....:....:....:....:....:....:....:....:
1234 12:....:.**.:....:....:....:....:....:....:....:....:....:....:....:....:
*/

void 
pr_placement (bool *excluded, int interval, int best, int num_intervals, int rightmost)
{
    if (interval==best) print("**.");
      else if (excluded[interval]) print("..."); else print(".*.");
    if (interval<num_intervals-1) print(":");
    if (interval<num_intervals-2 && interval<rightmost) print(".");
}


void 
pr_locus (int n)
{ 
    if (n<99) sprintf(ps, " %2d  ", n + 1);
    else sprintf(ps, "%4d ", n + 1);
    pr();
}


void 
pr_dist (real rf)
{ 
    sprintf(ps, "%3.0lf-:", cm_dist(rf));
    if (ps[0]==' ') ps[0]='-';
    if (ps[1]==' ') ps[1]='-';
    pr();
}


void 
new_print_placements (
    MAP *map,		/* framework */
    PLACEME **placed,
    int num_loci		/* placed loci, maybe NO_LOCUS */
)
{
    int num_remaining, num_across, num_done, num_intervals;
    int i, j, num_places, rightmost;

    if (map->sex_specific) send(CRASH);

    num_done=0; 
    num_intervals=map->num_loci+1;              /* interval #s */
    do {
	num_remaining= num_intervals-num_done;
	num_across= min(14,num_remaining);      /* interval #s */
	if (num_across==num_remaining-1) num_across++;
	rightmost=  num_done+num_across-1;      /* num_done is start */
	
	/**** draw loci, i=interval#, i-1, i= flanking-locus#'s ****/
	if (num_done==0) print("   "); /*3*/ else print("     "); /*5*/
	for (i=num_done; i<=rightmost; i++) { /* draw interval left locus */
	    if (i==0) print("     ");
	    else pr_locus(map->locus[i-1]);
	} /* now draw interval right locus */
	if (rightmost<num_intervals-1) pr_locus(map->locus[rightmost]);
	nl();

	/**** draw dists, i=interval# ****/
	if (num_done>0) print("       :"); /*3+4+:*/ else print("      ");/*5*/
	for (i=num_done; i<=rightmost; i++) { /* draw interval dist */
	    if (i==num_intervals-1) continue;
	    else if (i==0 || map->rec_frac[i-1][0]==NO_THETA) print("    :");
	    else pr_dist(map->rec_frac[i-1][0]);
	}
	nl();

	for (j=0; j<num_loci; j++) if (placed[j]->locus!=NO_LOCUS) {
	    /**** draw placements ****/
	    for (num_places=0, i=0; i<map->num_loci+1; i++) 
	      if (!placed[j]->excluded[i]) num_places++;
	    sprintf(ps, "%4d %-2d", placed[j]->locus + 1, num_places); pr();
	    if (num_done>0) print(":.");
	    
	    for (i=num_done; i<=rightmost; i++) /* draw interval dist */
	      pr_placement(placed[j]->excluded,i,placed[j]->best_pos,
			   num_intervals,rightmost);
	    nl();
	}
	num_done+= num_across;
	if (num_done<num_intervals) nl();
    } while (num_done<num_intervals);
}


#define NOER OBSCURE_REAL

void 
print_geno_line (int locus, real error_rate, int firsti, int lasti, int *obs, bool is_old, bool isa_haplo)
{
    char p1, p2, name[TOKLEN+1];
    int indiv;

    if (is_old) p1=p2=' '; else { p1='('; p2=')'; }

    sprintf(name, "%s%c", locname(locus, isa_haplo), p2);
    if (error_rate==NOER) sprintf(ps, "%c%-11s      ", p1, name);
      else sprintf(ps, "%c%-11s %s%% ", p1, name, rsd(3.1, error_rate * 100.0));
    pr();

    for (indiv=firsti; indiv<=lasti; indiv++) {
	switch (obs[indiv]) {
	    case OBS_A:       print("A"); break;
	    case OBS_B:       print("B"); break;
	    case OBS_H:       print("H"); break;
	    case OBS_NOT_A:   print("C"); break;
	    case OBS_NOT_B:   print("D"); break;
	    case OBS_MISSING: print("-"); break;
	    default:          print("?"); break;
	}
    }
    nl();
}


#define OBLIG_REC "%5.1lf cM %3d loci in indiv %-3d  %s(%c) - %s(H) - %s(%c)\n"

void 
print_f2_map_genotypes (
    MAP *map,
    bool use_haplos,
    bool explode_haplos,
    int num_old,
    int *old_locus  /* can omit if num_old==0 - should do this with VARARGS */
)
{
    int indiv, n_indivs, firsti, lasti;
    int locus, i, j, *obs=NULL, *prev_obs=NULL, *num_recs=NULL, old=0;
    int *last_homo=NULL, *last_het=NULL, *homo_was=NULL;
    real error_rate;

    if (map->sex_specific || raw.data_type!=F2) send(CRASH);

    firsti=0; n_indivs=raw.data.f2.num_indivs;
    array(obs,n_indivs,int);
    array(prev_obs,n_indivs,int);
    array(num_recs,n_indivs,int);
    array(last_homo,n_indivs,int);
    array(last_het,n_indivs,int);
    array(homo_was,n_indivs,int);

    do { /* loop over sets of indivs */
	if (n_indivs-firsti>=60) lasti=firsti+49; else lasti=n_indivs-1;
	
	if (n_indivs>=1000) {
	          /* 123456789012 */
	    print("                  ");
	    for (indiv=firsti; indiv<=lasti; indiv++)
	      { sprintf(ps, "%d", ((indiv + 1) / 1000) % 10); pr(); }
	    nl();
	}
	if (n_indivs>=100) {
	          /* 123456789012 */
	    print("                  ");
	    for (indiv=firsti; indiv<=lasti; indiv++)
	      { sprintf(ps, "%d", ((indiv + 1) / 100) % 10); pr(); }
	    nl();
	}
	print("                  ");
	for (indiv=firsti; indiv<=lasti; indiv++)
	  { sprintf(ps, "%d", ((indiv + 1) / 10) % 10); pr(); }
	print("\n                  ");
	for (indiv=firsti; indiv<=lasti; indiv++)
	  { sprintf(ps, "%d", (indiv + 1) % 10); pr(); }
	print("\n                  ");
	for (indiv=firsti; indiv<=lasti; indiv++) print("-");
	nl();

	for (indiv=firsti; indiv<=lasti; indiv++) {
	    prev_obs[indiv]=OBS_MISSING; num_recs[indiv]=0; 	    
	}

	for (i=0; i<map->num_loci; i++) {
	    f2_genotype(map->locus[i],use_haplos,obs);

	    if (i>0) {
		sprintf(ps, "         %s %s ", rf2str(map->rec_frac[i - 1][0]),
                (units==CENTIMORGANS ? "cM":"rf")); pr();
		for (indiv=firsti; indiv<=lasti; indiv++) {
		    if (map->error_lod!=NULL &&
			(map->error_lod[i][indiv]>=error_lod_thresh ||
			 map->error_lod[i-1][indiv]>=error_lod_thresh)) {
			print("|");
			num_recs[indiv]+= 
			  obligate_recs[obs[indiv]][prev_obs[indiv]];
		    } else if (obligate_recs[obs[indiv]][prev_obs[indiv]]==2) 
		      { print("O"); num_recs[indiv]+=2; }
		    else if (obligate_recs[obs[indiv]][prev_obs[indiv]]==1) 
		      { print("X"); num_recs[indiv]+=1; }
		    else print(" ");
		}
		nl();
	    }

	    for (indiv=firsti; indiv<=lasti; indiv++)
	      if (obs[indiv]!=OBS_NOT_A && obs[indiv]!=OBS_NOT_B)
		prev_obs[indiv]=obs[indiv];

	    if (num_old>0) for (old=FALSE, j=0; j<num_old; j++) 
	      if (old_locus[j]==map->locus[i]) { old=TRUE; break; }
	    if (map->error_rate!=NULL && i!=0 && i<map->num_loci-1)
	      error_rate=map->error_rate[i]; else error_rate=NOER;
	    print_geno_line(map->locus[i],error_rate,firsti,lasti,obs,
			    old,use_haplos && haplotyped(map->locus[i]));
	    
	    if (use_haplos && explode_haplos && haplotyped(map->locus[i]))
	      for (locus=haplo_first[map->locus[i]]; locus!=NO_LOCUS;
		   locus=haplo_next[locus]) {
		  f2_genotype(locus,FALSE,obs);
		  print_geno_line(locus,NOER,firsti,lasti,obs,old,FALSE);
	      }
	}

	print("                  ");
	for (indiv=firsti; indiv<=lasti; indiv++) print("-");
	print("\n           #Recs: ");
	for (indiv=firsti; indiv<=lasti; indiv++) 
	  if (num_recs[indiv]>9) print("*");
	  else { sprintf(ps, "%d", num_recs[indiv] % 10); pr(); }
	nl(); nl();

	firsti=lasti+1;
    } while (firsti<n_indivs);

#ifdef FIX_THIS_SOMEDAY
    nl();
    for (indiv=0; indiv<n_indivs; indiv++) {
	last_homo[indiv]= last_het[indiv]= NO_LOCUS;
    }
    got_one=FALSE;
    for (indiv=0; indiv<n_indivs; indiv++) {
	for (i=0; i<map->num_loci; i++) {
	    f2_genotype(map->locus[i],use_haplos,obs);
	    
	    /* find obligate doubles - ignoring all C's D's for now */
	    if (obs[indiv]==OBS_A || obs[indiv]==OBS_B) {
		if (last_het[indiv]==NO_LOCUS || obs[indiv]!=homo_was[indiv]) {
		    last_homo[indiv]=i; last_het[indiv]=NO_LOCUS;
		    homo_was[indiv]=obs[indiv];
		} else {
		    /* got one a-h-a, or b-h-b */
		    for (cm=0.0, j=last_homo[indiv]; j<i; j++) 
		      cm+=(*mapfunction->rec_to_dist)(map->rec_frac[j][0]);
		    cm*=100.0;
		    ch= (obs[indiv]==OBS_A ? 'A':'B');
		    if (!got_one) 
		      { print("Obligate Doubles:\n"); got_one=TRUE; }
		    num= i-last_homo[indiv]+1;
 		    sprintf(ps,OBLIG_REC,cm,num,indiv+1,
		       rag(loc2str(map->locus[last_homo[indiv]])),ch,
		       rag(loc2str(map->locus[last_het[indiv]])),
		       rag(loc2str(map->locus[i])),ch);
		    pr();
		}
	    }
	    if (obs[indiv]==OBS_H && last_homo[indiv]!=NO_LOCUS &&
		last_het[indiv]==NO_LOCUS) {
		  last_het[indiv]=i;
	    }
	}
    }
    if (!got_one) print("No obligate doubles anywhere.\n");
#endif
    unarray(obs,int); 
    unarray(prev_obs,int);
    unarray(num_recs,int);
    unarray(last_homo,int);
    unarray(last_het,int);
    unarray(homo_was,int);
}

