/******************************************************************************

  ####   #####   #####      #    #    #   #####           ####
 #    #  #    #  #    #     #    ##   #     #            #    #
 #    #  #    #  #    #     #    # #  #     #            #
 #  # #  #####   #####      #    #  # #     #     ###    #
 #   #   #       #   #      #    #   ##     #     ###    #    #
  ### #  #       #    #     #    #    #     #     ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

/********** QPRINT.C - QTL's PRINTING ROUTINES **********/

#define INC_LIB
#define INC_SHELL
#define INC_CALLQCTM
#define INC_QTOPLEVEL
#define INC_QLOWLEVEL /* why? */
#include "qtl.h"

char *interval_str();
char *trait_str();

void print_test_wiggle_map();      
void map_printer();


/* QTL MAP OUTPUT FORMAT:   STEVE changed! Again!

==============================================================
QTL map for trait 1 (no_error):

INTERVALS  LENGTH  QTL-POS  WEIGHT  DOMINANCE                  
234-678    100.0   100.0   -0.000  -0.000
1-2        100.0   100.0    0.000   0.000

CONTINUOUS-VARIABLE         WEIGHT
123 (1234567890)           -0.000

INTERVALS  LENGTH  POSITION  GENETICS    WEIGHT  DOMINANCE     <- !print_names
234-678    100.0   100.0     recessive  -0.000  -0.000
1-2        100.0   100.0     free        0.000   0.000

INTERVALS              LENGTH  POSITION  GENETICS    WEIGHT  DOMINANCE
1234567890 1234567890  100.0   100.0     recessive   0.000   0.000

chi^2= 12.456 (12 D.F.)        log-likelihood= 1234.67          
mean= 12.456  sigma^2= 12.456  variance-explained= 0.1 %
==============================================================  
*/

#define INT_TITLE_LINE	"INTERVAL%s%s  LENGTH  QTL-POS  "
#define INT_LINE 	"%s  %-5s   %-5s   "
#define GENETICS_TITLE  "GENETICS    "
#define GENETICS_FORMAT " %-9s  "
#define WEIGHT_TITLE    "WEIGHT  "
#define WEIGHT_FORMAT   "%s "  
#define DOM_TITLE       "DOMINANCE"
#define DOM_FORMAT      "%s"
#define LOD_TITLE       "  %VAR     LOD" /* FROB */
#define LOD_FORMAT      "  %4.2lf%%  %-5.2lf"
#define CV_TITLE        "CONTINUOUS-VARIABLE         WEIGHT\n"
#define CV_LINE         "%-20s       %s\n"

#define LIKE_LINE 	\
"chi^2= %-6.3lf (%-2d D.F.)        log-likelihood= %-7.2lf\n"
#define VAR_LINE 	\
"mean= %-6.3lf  sigma^2= %-6.3lf  variance-explained= %-4.1lf%%\n"

#define ITER_LINE 	\
"log-likelihood= %-10.5lf     delta-log-like= %-10.5lf\n"
#define EXTRA_LINE1  	\
"abs-log-like= %-10.5lf       no-data-like= %-10.5lf\n"
#define EXTRA_LINE2	\
"null-log-like= %-10.5lf      like-diff= %-10.5lf\n"

#define SHORT_MAP_MAX_STARS 10

void print_map_divider()
{ 
    if (print_names && raw.data_type==INTERCROSS) print(LONG_MAP_DIVIDER);
    else print(MAP_DIVIDER);
}


void print_qtl_map(map,free_genetics) /* has CONT_VAR hooks */
QTL_MAP *map;
bool *free_genetics;
{
    int i, print_genetics, df_per;

    if (free_genetics==NULL) print_genetics=TRUE; 
    else for (i=0, print_genetics=FALSE; i<map->num_intervals; i++) 
      if (!free_genetics[i]) { print_genetics=TRUE; break; }

    /* hold(6+map->num_intervals); */
    map_printer(map,print_genetics); nl();
    df_per= (raw.data_type==INTERCROSS ? 2:1);
    sprintf(ps, LIKE_LINE, map->chi_sq, df_per * map->num_intervals, map->log_like); pr();
    sprintf(ps, VAR_LINE, map->mu, map->sigma_sq, map->var_explained * 100.0); pr();
    
    /*	sprintf(ps,EXTRA_LINE1,map->abs_log_like,map->no_data_like); pr();
	sprintf(ps,EXTRA_LINE2,map->null_log_like,
	map->no_data_like-map->null_log_like); pr(); */
}


void map_printer(map,print_genetics)
QTL_MAP *map;
bool print_genetics;
{
    int i;
    char *trait_string= get_temp_string();
	
    if (map->num_intervals>0) {
	sprintf(ps, INT_TITLE_LINE, map->num_intervals > 1 ? "S" : " ",
	   print_names ? "            ":""); pr();
	if (print_genetics) print(GENETICS_TITLE);
	print(WEIGHT_TITLE);
	if (raw.data_type==INTERCROSS) print(DOM_TITLE);
	nl();
    }
     for (i=0; i<map->num_intervals; i++) {
	sprintf(ps, INT_LINE, interval_str(map->left[i], map->right[i], TRUE),
            dist_str(map->interval_len[i],FALSE),
            dist_str(map->qtl_pos[i],FALSE)); pr();
	if (print_genetics) { 
	    sprintf(ps, GENETICS_FORMAT, genetics_str(&map->constraint[i], FALSE)); 
	    pr();
	}
	sprintf(ps, WEIGHT_FORMAT, rsn(7.5, map->qtl_weight[i])); pr();
	if (raw.data_type==INTERCROSS) 
	  { sprintf(ps, DOM_FORMAT, rsn(7.5, map->qtl_dominance[i])); pr(); }
	nl();
    }

    if (map->num_continuous_vars>0) {
	if (map->num_intervals>0) nl(); /* divider line */
	print(CV_TITLE);
    }
    for (i=0; i<map->num_continuous_vars; i++) {
	if (map->cont_var[i]==EPISTASIS_TERM) strcpy(trait_string,"epistasis");
	else sprintf(trait_string, "%d (%s)", map->cont_var[i] + 1,
                 raw.trait_name[map->cont_var[i]]);
	   
	sprintf(ps, CV_LINE, trait_string, rsn(7.5, map->cont_var_weight[i])); pr();
    }
}


void print_short_title() /* FROB */
{
    sprintf(ps, INT_TITLE_LINE, " ", print_names ? "            " : ""); pr();
    print(WEIGHT_TITLE);
    if (raw.data_type==INTERCROSS) print(DOM_TITLE);
    print(LOD_TITLE);
    nl();
}


void print_short_qtl_map(map,threshold,scale) /* FROB */
QTL_MAP *map;
real threshold, scale;
{
    int i, n, last;
    last= map->num_intervals-1;

    sprintf(ps, INT_LINE, interval_str(map->left[last], map->right[last], TRUE),
            dist_str(map->interval_len[last],FALSE),
            dist_str(map->qtl_pos[last],FALSE)); pr();
    sprintf(ps, WEIGHT_FORMAT, rsn(7.5, map->qtl_weight[last])); pr();
    if (raw.data_type==INTERCROSS) 
      { sprintf(ps, DOM_FORMAT, rsn(7.5, map->qtl_dominance[last])); pr(); }
    sprintf(ps, LOD_FORMAT, map->var_explained * 100.0, map->log_like); pr();
    if (map->log_like>=threshold) {
	n= ((int)((map->log_like-threshold)/scale))+1;
	if (n>SHORT_MAP_MAX_STARS) n=SHORT_MAP_MAX_STARS;
	for (i=0; i<n; i++) print("*");
    }
    nl();
}


void print_iteration(iter,map,delta_log_like)
int iter;
QTL_MAP *map;
real delta_log_like;
{
	print(ITER_DIVIDER);
	sprintf(ps, "iteration #%d:\n", iter + 1); print(ps);
	map_printer(map,FALSE);
	sprintf(ps, ITER_LINE, map->log_like, delta_log_like);
	print(ps);
	sprintf(ps, VAR_LINE, map->mu, map->sigma_sq, map->var_explained * 100.0); 
	print(ps);
}


#define NULL_MAP_LINE \
"null-like= %-8.3lf  null_sigma^2= %-6.3lf \nno-data-like=%-8.3lf\n"

void print_null_iteration(map)
QTL_MAP *map;
{
	print(ITER_DIVIDER);
	sprintf(ps, NULL_MAP_LINE, map->null_log_like, map->null_sigma_sq,
            map->no_data_like);
	print(ps);
}


/* has CONT_VAR hooks */
void do_print_E_step(expected_genotype,S_matrix,expected_recs,n_individuals,
		     n_genotype_vars,n_continuous_vars,n_intervals)
real **expected_genotype, **S_matrix;
GENO_PROBS *expected_recs;
int n_individuals, n_genotype_vars, n_continuous_vars, n_intervals;
{
    int i, j, n, g, q, m, num;

    num= n_genotype_vars + n_continuous_vars;

    nl();
    print("expected_genotype:\nindiv");
    for (n=0; n<num+1; n++) 
      { sprintf(ps, "  %-2d   ", n); pr(); } nl();
    for (i=0; i<n_individuals; i++) {
	sprintf(ps, "%-4d ", i); pr();
	for (n=0; n<num+1; n++)
	    { sprintf(ps, "%-6.3lf ", expected_genotype[i][n]); pr(); }
	nl();
    }

    print("\nS_matrix:\n   ");
    for (n=0; n<num+1; n++) { sprintf(ps, "  %-2d   ", n); pr(); } nl();
    for (n=0; n<num+1; n++) {
	sprintf(ps, "%-2d ", n); pr();
	for (m=0; m<num+1; m++) 
	  { sprintf(ps, "%-6.2lf ", S_matrix[n][m]); pr(); }
	nl();
    }

    print("\nexpected recs: ");
    for (j=0; j<n_intervals; j++) {
	sprintf(ps, "interval %d:\n", j); pr();
	for_locus_genotypes(raw.data_type,q) {
	    for_interval_genotypes(raw.data_type,g)      
	      { sprintf(ps, "[%d][%d]  ", g, q); pr(); } nl();
	    for_interval_genotypes(raw.data_type,g)
	      { sprintf(ps, "%7.5lf ", expected_recs[j][g][q]); pr(); } nl();
	}
    }
    nl();
}


/* WIGGLES OUTPUT FORMAT.... (B1 AND F2, RESPECTIVELY)

POS     WEIGHT  %VAR  LOG-LIKE |
-------------------------------| (1 2) 10.0 cm
123.5  -2.456  12.4%  123.567  | 1234567890123456789012345678901234567890123456
-------------------------------|

POS     WEIGHT  DOM     %VAR  LOG-LIKE |
---------------------------------------| (1234567890 1234567890) 100.0 cm
123.5  -2.456  -2.456  12.4%  123.567  | 12345678901234567890123456789012345678
---------------------------------------|
*/

#define B1_WIGGLE_TITLE "POS     WEIGHT  %VAR  LOG-LIKE | "
#define B1_WIGGLE_MAP   "%5s  %s  %4.1lf%%  %s  | "
#define B1_NULL_WIGGLE  "-------------------------------|\n"
#define B1_WIGGLE_INTERVAL	\
"-------------------------------| %s %s %s\n"
#define B1_SPACES_LEFT 43

#define F2_WIGGLE_TITLE "POS     WEIGHT  DOM     %VAR  LOG-LIKE | "
#define F2_WIGGLE_MAP   "%5s  %s  %s  %4.1lf%%  %s  | "
#define F2_NULL_WIGGLE  "---------------------------------------|\n"
#define F2_WIGGLE_INTERVAL	\
"---------------------------------------| %s %s %s\n"
#define F2_SPACES_LEFT 35


void print_wiggle_title() 
{ 
    if (raw.data_type==BACKCROSS) { print(B1_WIGGLE_TITLE); }
    else { print(F2_WIGGLE_TITLE); }
    nl();
}


void print_wiggle_interval(map)
QTL_MAP *map;
{ 
    int i;

    if (map==NULL) { 
	if (raw.data_type==BACKCROSS) strcpy(ps,B1_NULL_WIGGLE); 
	else strcpy(ps,F2_NULL_WIGGLE); 
    } else {
	i= map->num_intervals-1;
	sprintf(ps,
            (raw.data_type==BACKCROSS ? B1_WIGGLE_INTERVAL:F2_WIGGLE_INTERVAL),
            interval_str(map->left[i],map->right[i],FALSE),
            dist_str(map_length(map->left[i],map->right[i]),FALSE),
            units_str(FALSE)); 
    }
    pr();
}
  

void print_wiggle_map(map,base_like,scale)
QTL_MAP *map;
real base_like, scale;
{
    int i, j, stars;

    if (map==NULL || scale<=0.0) send(CRASH);
    i= map->num_intervals-1;

    if (raw.data_type==BACKCROSS)
      sprintf(ps, B1_WIGGLE_MAP, dist_str(map->qtl_pos[i], TRUE),
              rsn(6.3,map->qtl_weight[i]),
	 map->var_explained*100.0, rsd(7.3,map->log_like));
    else /* data_type==INTERCROSS */
      sprintf(ps, F2_WIGGLE_MAP, dist_str(map->qtl_pos[i], TRUE),
              rsn(6.3,map->qtl_weight[i]), rsn(6.3,map->qtl_dominance[i]),
	    map->var_explained*100.0, rsd(7.3,map->log_like));
    pr();

    if (map->log_like>base_like) {
	stars=iminf((int)((map->log_like-base_like+scale+0.001)/scale),
	      (raw.data_type==BACKCROSS ? B1_SPACES_LEFT:F2_SPACES_LEFT)+1);
	for (j=0; j<stars; j++) print("*");
    }
    nl();
}


void print_wiggle_genetics(genetics)
GENETICS *genetics;
{
    if (constrained(genetics)) {
	print("Scanned QTL genetics are constrained to be: ");
	print(genetics_str(genetics,TRUE)); nl();
    } else print("Scanned QTL genetics are free.\n");
}


/* left_seq_str() in this file STILL needs CONT_VAR hooks! */

void print_wiggle_left_seq(map) 
QTL_MAP *map;
{
    char *left_seq;
    
    if (nullstr(left_seq=left_seq_str(map))) print("No fixed-QTLs.\n");
    else { print("Fixed-QTLs: "); print(left_seq); nl(); }
}


#define WIG_LIST_TITLE "  NUM  TRAIT%s  SEQUENCE\n"
#define WIG_LIST_NAMES "%3d.%s  %-10s  %s\n"
#define WIG_LIST_NUMS  "%3d.%s  %3d    %s\n"

void print_saved_wiggles()
{
    int i;
    WIGGLE_OPERATION *op;

    if (first_wiggle==num_wiggles) return;
    sprintf(ps, WIG_LIST_TITLE, (print_names ? "     " : "")); pr();
    for (i=first_wiggle; i<num_wiggles; i++) {
	if ((op=wiggles[i])==NULL) send(CRASH);
	if (op->data == NULL)
	  sprintf(ps, "  <deleted>\n");
	else if (print_names) 
	  sprintf(ps, WIG_LIST_NAMES, i + 1, (op->num_orders == 1 ? "1" : "x"),
              raw.trait_name[op->trait], op->seq_string);
	else sprintf(ps, WIG_LIST_NUMS, i + 1, (op->num_orders == 1 ? "1" : "x"),
		op->trait+1, op->seq_string);
	pr(); 
    }
} 


#define ORD_LIST_TITLE "  NUM    GENETICS   FIXED-QTLS\n"
#define ORD_LIST_FORM  "%3d.%-3d  %-9s  %s\n"

void print_saved_wiggle(wiggle)
int wiggle;
{
    int i, k;
    WIGGLE_OPERATION *op;

    if (wiggle<0) return;
    if (wiggle>=num_wiggles || (op=wiggles[wiggle])==NULL) send(CRASH);

    print(ORD_LIST_TITLE);
    k=op->num_intervals-1;
    for (i=0; i<op->num_orders; i++) {
	sprintf(ps, ORD_LIST_FORM, wiggle + 1, i + 1,
            genetics_str(&op->data[i][0]->map->constraint[k],FALSE),
            left_seq_str(op->data[i][0]->map));
	pr();
    }
}
	

void print_saved_wiggle_order(wiggle,order,base_like,scale)
int wiggle, order;
real base_like, scale;
{
    int i, j, k, stars;
    WIGGLE_OPERATION *op;
    WIGGLE_INTERVAL **data;
    WIGGLE_POINT *point;

    if (wiggle<0) return;
    if (wiggle>=num_wiggles || (op=wiggles[wiggle])==NULL) send(CRASH);
    if (order<0 || order>=op->num_orders || base_like<0.0 || scale<=0.0) 
      send(CRASH);
    data=op->data[order];

    print_wiggle_title();
    for (i=0; i<op->num_wiggled_intervals; i++) {
	if (i>0 && !data[i]->contig) { print_wiggle_interval(NULL); nl(); }
	print_wiggle_interval(data[i]->map);
	for (j=0; j<data[i]->num_points; j++) {
	    point=data[i]->point[j];
	    if (raw.data_type==BACKCROSS)
	      sprintf(ps, B1_WIGGLE_MAP, dist_str(point->qtl_pos, TRUE),
                  rsn(6.3,point->qtl_weight),point->var_explained*100.0,
                  rsd(7.3,point->lod_score));
	    else /* data_type==INTERCROSS */
	      sprintf(ps, F2_WIGGLE_MAP, dist_str(point->qtl_pos, TRUE),
                  rsn(6.3,point->qtl_weight), rsn(6.3,point->qtl_dominance),
		 point->var_explained*100.0, rsd(7.3,point->lod_score));
	    pr();
	    if (point->lod_score>base_like) {
		stars=
		 iminf((int)((point->lod_score-base_like+scale+0.001)/scale),
		 (raw.data_type==BACKCROSS ?B1_SPACES_LEFT:F2_SPACES_LEFT)+1);
		for (k=0; k<stars; k++) print("*");
	    }
	    nl();
	} 
    }
    print_wiggle_interval(NULL); 
}


#define PEAK_TITLE   "QTL-Map for peak %d:\n"
#define PEAK_L_CONF  "Confidence Interval:  Left Boundary= "
#define PEAK_R_CONF  "                     Right Boundary= "
#define PEAK_BOUND   "%s + %s\n"
#define PEAK_NOBOUND "%s (off end)\n"

void print_peak(peak,num)
WIGGLE_PEAK *peak;
int num;
{
    char *i;
    
    print_map_divider();
    sprintf(ps, PEAK_TITLE, num + 1); pr();
    print(PEAK_L_CONF);
    if (peak->backward_pos==OFF_END) {
	i=interval_str(peak->backward_left,peak->backward_right,FALSE);
	sprintf(ps, PEAK_NOBOUND, i); pr();
    } else { 
	i=interval_str(peak->backward_left,peak->backward_right,FALSE);
	sprintf(ps, PEAK_BOUND, i, dist_str(peak->backward_pos, FALSE)); pr();
    }
    print(PEAK_R_CONF);
    if (peak->forward_pos==OFF_END) {
	i=interval_str(peak->forward_left,peak->forward_right,FALSE);
	sprintf(ps, PEAK_NOBOUND, i); pr();
    } else { 
	i=interval_str(peak->forward_left,peak->forward_right,FALSE);
	sprintf(ps, PEAK_BOUND, i, dist_str(peak->forward_pos, FALSE)); pr();
    }
    nl();
    if (peak->map==NULL) send(CRASH); /* KLUDGE for now */
    print_qtl_map(peak->map,NULL);
}




/* SAVED TEST WIGGLES OUTPUT FORMAT... (F2 DATA ONLY)
-----------------------------------------------------------------------------|
GENETICS:             FREE         | DOMINANT    | RECESSIVE   | ADDITIVE    |
-----------------------------------------------------------------------------|
                              LOG  |        LIKE |        LIKE |        LIKE |
POS    WEIGHT  DOM    %VAR    LIKE | %VAR   DIFF | %VAR   DIFF | %VAR   DIFF |
-----------------------------------------------------------------------------|
interval= xxxxx length= 12.45 cm                                             |
-----------------------------------------------------------------------------|
123.4 -2.345  -2.345  12.4  123.56 | 12.4 -12.67 | 12.4 -12.56 | 12.4 -12.56 |
-----------------------------------------------------------------------------|
*/

#define TEST_WIGGLE_DIVIDER \
"-----------------------------------------------------------------------------|"
#define TEST_WIGGLE_TITLE \
"GENETICS:             FREE         | DOMINANT    | RECESSIVE   | ADDITIVE    |"
#define TEST_WIGGLE_HEADER1 \
"                              LOG  |        LIKE |        LIKE |        LIKE |"
#define TEST_WIGGLE_HEADER2 \
"POS    WEIGHT  DOM    %VAR    LIKE | %VAR   DIFF | %VAR   DIFF | %VAR   DIFF |"

#define TEST_WIGGLE_INTERVAL1 "interval= %s length= %s %s"
#define TEST_WIGGLE_MAP1 "%5s %s  %s  %4.1lf  %s | "
#define TEST_WIGGLE_MAP2 "%4.1lf %s | %4.1lf %s | %4.1lf %s |%s"

/* We had to split up the TEST_WIGGLE_MAP line because of a (compiler?) bug 
   on the HP800. Leave it as is. */


void print_test_wiggle_map(point,threshold)
WIGGLE_POINT **point; 
real threshold;
{
    real lod; lod= point[FREE]->lod_score;

    sprintf(ps, TEST_WIGGLE_MAP1,
            dist_str(point[FREE]->qtl_pos,TRUE),
            rsn(6.3,point[FREE]->qtl_weight),
            rsn(6.3,point[FREE]->qtl_dominance),
       point[FREE]->var_explained*100.0,
            rsd(6.2,point[FREE]->lod_score)); pr();

    sprintf(ps, TEST_WIGGLE_MAP2,
       point[DOMINANT]->var_explained*100.0,
            rsd(6.2,point[DOMINANT]->lod_score-lod),
       point[RECESSIVE]->var_explained*100.0,
            rsd(6.2,point[RECESSIVE]->lod_score-lod),
       point[ADDITIVE]->var_explained*100.0,
            rsd(6.2,point[ADDITIVE]->lod_score-lod),
            (lod+0.00499>=threshold ? "*\n":"\n")); pr(); 
}


void print_test_wiggle_interval(map)
QTL_MAP *map;
{ 
    int i, last; 
    print(TEST_WIGGLE_DIVIDER); nl(); 
    if (map==NULL) return;

    last= map->num_intervals-1;
    sprintf(ps, TEST_WIGGLE_INTERVAL1,
            interval_str(map->left[last],map->right[last],FALSE),
            dist_str(map->interval_len[last],FALSE), units_str(FALSE)); pr();
    for (i=len(ps);i<77; i++) print(" ");
    print("|\n"); print(TEST_WIGGLE_DIVIDER); nl();
}


void print_test_wiggle_title()
{ print(TEST_WIGGLE_DIVIDER); nl(); print(TEST_WIGGLE_TITLE); nl(); 
  print(TEST_WIGGLE_DIVIDER); nl(); print(TEST_WIGGLE_HEADER1); nl();  
  print(TEST_WIGGLE_HEADER2); nl(); } 


void print_test_wiggle_order(wiggle,order,threshold)
int wiggle, order;
real threshold;
{
    int i, j, test;
    WIGGLE_OPERATION *op;
    WIGGLE_INTERVAL **data[NUM_MODELS];
    WIGGLE_POINT *point[NUM_MODELS];

    if (wiggle<0) return;
    if (wiggle>=num_wiggles || (op=wiggles[wiggle])==NULL) send(CRASH);
    if (order<0 || order>=op->num_orders || threshold<0.0)
      send(CRASH);
    for (test=0; test<NUM_MODELS; test++) data[test]=op->data[order+test];
      
    print_test_wiggle_title();
    for (i=0; i<op->num_wiggled_intervals; i++) {
	if (i>0 && !data[FREE][i]->contig) 
	  { print_test_wiggle_interval(NULL); nl(); }
	print_test_wiggle_interval(data[FREE][i]->map);

	for (j=0; j<data[FREE][i]->num_points; j++) {
	    for (test=0; test<NUM_MODELS; test++) 
	      point[test]=data[test][i]->point[j];
	    print_test_wiggle_map(point,threshold);
	} 
    }
    print_test_wiggle_interval(NULL); 
}



/*****************************************************************************/

void print_trait(for_num_maps) 
int for_num_maps;
{ sprintf(ps, "QTL map%s for trait %s:\n", maybe_s(for_num_maps), trait_str()); pr(); }


void print_seq()
{ print("Sequence: "); print(ints_string); nl(); }


void print_old_seq(str)
char *str;
{ print("Sequence: "); print(str); nl(); }


char *trait_str()
{
    char *str= get_temp_string(); 
    sprintf(str, "%d (%s)", trait + 1, raw.trait_name[trait]); 
    return(str);
}
    

char *interval_str(left,right,fill)
int left, right;
bool fill;
{
    char *str= get_temp_string();
    if (print_mapm_loci) {
	if (!print_names) sprintf(str, "%d-%d", raw.original_locus[left],
                              raw.original_locus[right]);
	else sprintf(str, "%s-%s", raw.locus_name[left],
                 raw.locus_name[right]);
	if (fill) pad_to_len(str,print_names ? 21:9);
	return(str);
    }
    else {
	if (!print_names) sprintf(str, "%d-%d", left + 1, right + 1);
	else sprintf(str, "%s-%s", raw.locus_name[left], raw.locus_name[right]);
	if (fill) pad_to_len(str,print_names ? 21:9);
	return(str);
    }
}

char *dist_str(rec_frac,fill)  /* Always takes 5 spaces */
real rec_frac;
bool fill; /* unused? */
{
    char *str; str=get_temp_string();

    if (units == CENTIMORGANS) {
	if (fill) sprintf(str, "%-5.1lf", haldane_cm(rec_frac));
	  else sprintf(str, "%-3.1lf", haldane_cm(rec_frac));
	str[5]='\0'; /* just in case */
    } else sprintf(str, "%-5.3lf", rec_frac);

    return(str);
}


char *units_str(fill)  /* Takes 5 spaces, if fill, 2 otherwise */
bool fill;
{
    char *str; str=get_temp_string();

    if (units == CENTIMORGANS) strcpy(str,"cM   ");
    else strcpy(str,"RF   ");
    if (!fill) str[2]='\0';

    return(str);
}


char *genetics_str(genetics,verbose)
GENETICS *genetics;
bool verbose;
{ 
    char *str= get_temp_string();

    if (!data_loaded()) strcpy(str,"");

    else if (raw.data_type==BACKCROSS) {
	if (genetics->backx_weight==DONT_FIX) strcpy(str,"free");
	else if (verbose) sprintf(str, "W=%-5.3lf", genetics->backx_weight);
	else strcpy(str,"fixed");

    } else if (raw.data_type==INTERCROSS)
      switch (genetics->interx_type) {
	  case FREE: 		strcpy(str,"free"); 	   break;
	  case DOMINANT: 	strcpy(str,"dominant");    break;
	  case ADDITIVE: 	strcpy(str,"additive");    break;
	  case RECESSIVE: 	strcpy(str,"recessive");   break;
	  case F3DOMINANT: 	strcpy(str,"f3dom");       break;
	  case F3RECESSIVE: 	strcpy(str,"f3rec");       break;
	  case TEST_MODELS:     strcpy(str,"try");         break;
	  case CONSTRAINED: 	
	      if (!verbose) strcpy(str,"constrain");
	      else sprintf(str, "constraints: A=%-4.2lf B=%-4.2lf C=%-4.2lf",
                       genetics->a, genetics->b, genetics->c);
	      break;
	  case FIXED:
	      if (!verbose) strcpy(str,"fixed"); 
	      else sprintf(str, "fixed: A=%-5.3lf D=%-5.3lf", genetics->a,
                       genetics->b);
	      break;
	  default: send(CRASH);
      }

    else send(CRASH);
    return(str);
}


char *left_seq_str(map)
QTL_MAP *map;
{
    char *str=get_temp_string(), *qtl=get_temp_string();
    char *interval, *genetics, pos[11];
    int j;

    for (j=0; j<map->num_intervals-1; j++) { /* all but rightmost interval */
	interval= interval_str(map->left[j],map->right[j],FALSE);
	genetics= genetics_str(&map->constraint[j],FALSE);
	if (streq(genetics,"free")) genetics[0]='\0';
	  else strins(genetics," :");
	if (map->fix_pos[j]==DONT_FIX) pos[0]='\0';
	  else sprintf(pos,"+%s",dist_str(map->fix_pos[j],FALSE));
	sprintf(qtl,"%s[%s%s%s]",(str[0]!='\0' ?" ":""),interval,pos,genetics);
	maxstrcat(str,qtl,TEMP_STRING_LEN);
    }
    return(str);
}
		


/********************* TINY QTL-MAP FORMAT: ****************************

     INTERVALS                                                %VAR  LOG-LIKE
#123 13-56 23-56 33-56 43-56 53-56 63-56 73-56 83-56 93-56    12.4  123.56
 POS: 00.0 
   W: 0.000
   D: 0.000

#123 123-567 223-567 323-567 423-567 523-567 623-567 723-567  12.4  123.56
#123 1234567890-1234567890 1234567890-1234567890  12.4  1234.67 */

#define TINY_TITLE1 "     INTERVAL%s:"
#define TINY_TITLE2 " %%VAR  %s\n"
#define TINY_LINE1  "#%-3d "   
#define TINY_LINE2  " %4.1lf  %-7.2lf\n"

void print_tiny_map(map,num,offset) /* unused, and broken */
QTL_MAP *map;
int num; /* number for map (starts at 0) or num_intervals if map is null */
real offset;
{
    int i, spaces_per_int, ints_per_line, columns, right, n_ints=0;
    char *interval;
    
    if (map==NULL || map->num_continuous_vars>0) send(CRASH);

    if (!print_names) {
	spaces_per_int=raw.n_loci<=99 ? 6:8; 
	ints_per_line= raw.n_loci<=99 ? 9:7; 
    } else { spaces_per_int=22; ints_per_line=2; }
    
    columns= (spaces_per_int * map->num_intervals);
    if (columns>60) columns=60;
    right=imaxf(5+columns,15);

    if (num<0) { /* header only */
	sprintf(ps, TINY_TITLE1, n_ints > 1 ? "S" : " "); pr(); to_column(right); 
	sprintf(ps, TINY_TITLE2, (offset > 0.0 ? "LIKE-DIFF" : "LOG-LIKE")); pr();

    } else { /* print the map */
	sprintf(ps, TINY_LINE1, num + 1); pr();
	for (i=0; i<map->num_intervals; i++) {
	    interval= interval_str(map->left[i],map->right[i],FALSE);
	    print(interval); print(" ");
	    if ((i+1)%ints_per_line==0 && i!=map->num_intervals-1) 
	      print("\n    ");
	}
	to_column(right);	    
	sprintf(ps, TINY_LINE2, map->var_explained * 100.0, map->log_like - offset); pr();
    }
}



/***** Show-compares printing functions *****/

#define COMP_LIST_TITLE "  NUM  TRAIT%s  SEQUENCE\n"
#define COMP_LIST_NAMES "%3d.%s  %-10s  %s\n"
#define COMP_LIST_NUMS  "%3d.%s  %3d    %s\n"

void print_saved_compares()
{
    int i;
    COMPARE_OPERATION *op;

    if(first_compare == num_compares)  return;
    sprintf(ps, COMP_LIST_TITLE, (print_names ? "     " : "")); pr();
    for(i=first_compare; i<num_compares; i++) {
	if ((op=compares[i]) == NULL) send(CRASH);
	if (op->data==NULL) 
	  sprintf(ps, "  This compare (number%d) has been deleted!\n", i + 1);
	else if(print_names)
	  sprintf(ps, COMP_LIST_NAMES, i + 1, (op->num_contigs == 1 ? "1" : "x"),
              raw.trait_name[op->trait], op->seq_string);
	else
	  sprintf(ps, COMP_LIST_NUMS, i + 1, (op->num_contigs == 1 ? "1" : "x"),
	     op->trait+1, op->seq_string);
	pr();
    }
}


#define GOOD_MAPS_TITLE "QTL-maps above LOD %sfalloff:\n"
#define BAD_MAPS_TITLE  "QTL-maps below LOD %sfalloff:\n"
#define COMP_MAP_TOP    "QTL-map #%d:  LOD score difference= %-5.2lf\n"

void print_best_saved_maps(compare,contig,threshold,falloff)
int compare, contig;
real threshold, falloff;
{
    COMPARE_OPERATION *op;
    int start, i, *like_index, num_in_list, j;
    bool bad_maps_yet, good_maps_printed;
    real best, *like_list;

    like_index= NULL;  like_list= NULL;
    run {
	best= 0.0;
	array(like_list,raw.n_loci,real);
	array(like_index,raw.n_loci,int);

	if ((op=compares[compare])==NULL) send(CRASH);

	for (i=0, start=0; i<contig; i++) 
	  { while (op->data[start]->contig) start++; }

	like_list[0]=op->data[start]->map->log_like;
	like_index[0]=start;
	num_in_list=1;
	start++;
	
/* while (op->data[start]->contig) {  */

	for(j = start; j<op->num_orders; j++) {
	    i=num_in_list;
	    while (op->data[start]->map->log_like > like_list[i-1]) {
		like_list[i]=like_list[i-1];
		like_index[i]=like_index[i-1];
		i--;
		if (i==0) break;
	    }
	    like_list[i]=op->data[start]->map->log_like;
	    like_index[i]=start;
	    num_in_list++; start++;
	}
	best=like_list[0];

	bad_maps_yet=good_maps_printed=FALSE;
	if (best>=threshold) { 
	    sprintf(ps, GOOD_MAPS_TITLE, (falloff < 0.0 ? "threshold minus " : "")); pr();
	    print_map_divider();
	    good_maps_printed=TRUE;
	}

	for (i=0; i<num_in_list; i++) {
	    map=op->data[like_index[i]]->map;
	    if (good_maps_printed && ((falloff>=0.0 && map->log_like>falloff)||
		    (falloff<0.0 && map->log_like > best+falloff))) {
		sprintf(ps, COMP_MAP_TOP, i + 1, map->log_like - best); pr(); nl();
		print_qtl_map(map,NULL);
		print_map_divider();
	    } else {
		if (!bad_maps_yet) {
		    if (good_maps_printed) nl(); 
		    sprintf(ps, BAD_MAPS_TITLE,
                    (falloff<0.0 ? "threshold minus ":"")); pr();
		    print_tiny_map(map,-1,0.0); /* title */
		    bad_maps_yet=TRUE;
		}
		print_tiny_map(map,i,0.0);
	    }
	}

    } on_exit {
	unarray(like_index,int);
	unarray(like_list,real);
	relay_messages;
    }
}


#define SAVED_BEST_LIKE \
"LOD score maximum for these maps: %-4.2lf\n"

void print_saved_maps(compare,contig)
int compare, contig;
{
    COMPARE_OPERATION *op;
    int i, start;
    real best, like;

    if ((op=compares[compare])==NULL) send(CRASH);
    best= 0.0;

    for (i=0, start=0; i<contig; i++) 
      { while (op->data[start]->contig) start++; }
    for (i=start; i<op->num_orders; i++) {
	if (i>start && !op->data[i]->contig) break;
	if ((like=op->data[i]->map->log_like)>best) best=like;
    }
    
    sprintf(ps, SAVED_BEST_LIKE, best); pr(); nl();
    print_tiny_map(op->data[0]->map,-1,1.0); /* title */
    
    for (i=start; i<op->num_orders; i++) {
	if (i>start && !op->data[i]->contig) break;
	print_tiny_map(op->data[i]->map,i,best);
    }
}

void get_fixed_qtl_weights(map)
QTL_MAP *map;
{
    int i;

    if (raw.data_type==BACKCROSS) print("Enter weight for QTL(s):\n");
    else print("Enter weight and dominance for QTL(s):\n");
    for (i=0; i<map->num_intervals; i++) {
	print(interval_str(map->left[i],map->right[i],FALSE));
	getln(": ");
	if (!rtoken(&ln,rREQUIRED,&map->qtl_weight[i])) error("bad weight!");
	if (raw.data_type==INTERCROSS &&
	    !rtoken(&ln,rREQUIRED,&map->qtl_dominance[i]))error("bad weight!");
    }
    print("Enter weight for cont-vars:\n");
    for (i=0; i<map->num_continuous_vars; i++) {
	if (map->cont_var[i]==EPISTASIS_TERM) strcpy(ps,"epistasis: ");
	else sprintf(ps, "%d (%s): ", map->cont_var[i] + 1, raw.trait_name[map->cont_var[i]]);
	getln(ps);
	if (!rtoken(&ln,rREQUIRED,&map->cont_var_weight[i]))
	  error("bad weight!");
    }
/*    getln("Enter mu and sigma-squared: ");
    if (!rtoken(&ln,rREQUIRED,&map->mu) ||
	!rtoken(&ln,rREQUIRED,&map->sigma_sq)) 	  error("bad weight!"); */
}




