/******************************************************************************

  ####   #####     ##    #    #           ####
 #    #  #    #   #  #   #    #          #    #
 #    #  #    #  #    #  #    #          #
 #  # #  #####   ######  # ## #   ###    #
 #   #   #   #   #    #  ##  ##   ###    #    #
  ### #  #    #  #    #  #    #   ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

/* QTL raw data */

//#define INC_LIB
//#define INC_SHELL
//#define INC_QLOWLEVEL
//#define INC_CALLQCTM
//#define INC_QTOPLEVEL
#include <stdlib.h>
//#undef CRASH
#include "qtl.h"

int field(char **p_str, int length, char *val);   /* OBSOLETE */

/* External statics */
RAW raw;
int **make_interval_genotype; 
int *left_genotype, *right_genotype;
real *map_space;
int *order;
int segregation_distortion;
bool altered_chroms;
bool already_printed;

/* Internal interfaces */
char name_chars[70];
int nam_len;
#define LEFT 0
#define RIGHT 1
#define F2VERSION 3

void make_count_recs_etc(void);
void read_map_locus(FILE *fp, int indivs, int t_loc, int *order, int n_loci);
real read_map_distance(FILE *fp);
void read_quant_trait(FILE *fp, int num, int indivs);
void read_olddata(FILE *fp, char *path);
void read_oldmap_locus(FILE *fp, int num, int indivs);
void read_oldquant_trait(FILE *fp, int num, int indivs);
void read_oldmap_distance(FILE *fp, int num);
void initial_prob(LOCUS_GENOTYPE_PROBS *prob, int locus, int this_genotype);
real pheno_given_geno(int observation, int genotype);
void condition(LOCUS_GENOTYPE_PROBS *prob, int locus, int prev_locus, int observation, real rec_frac, int side);
void altered_chroms_message(void);
//void getdataln();
//void read_map_locus();
//void read_quant_trait();
//real interval_likelihood();
//real pheno_given_geno();
//
//void make_count_recs_etc();
//
//bool missing_data(); /* OBSOLETE, NOT USED */
//void initial_prob();
//void condition();
//void altered_chroms_message();

#define MAXEXCEED "max # indiv or inter exceeded"

#define FILE_MISMATCH "\
The '.data' and '.trait' files were not prepared from the same raw\n\
data file. You will need to re-prepare your data in MAPMAKER."

#define OLD_FILENUM "\
This file was prepared using an older version of MAPMAKER. You will\n\
need to re-prepare it using the 'prepare data' command in MAPMAKER."

/********** PRELIMINARIES **********/

void 
raw_init (void)
{
    raw.n_loci= 0; 		raw.n_traits= 0;
    raw.trait= NULL; 		raw.trait_name= NULL;
    raw.locus= NULL; 		raw.locus_name= NULL;
    raw.map_dist= NULL;		raw.data_type= NO_DATA; raw.f3=FALSE;
    raw.left_cond_prob= NULL;	raw.right_cond_prob= NULL;
    raw.file[0]='\0';	        raw.trait_eqn= NULL; 
    strcpy(default_backcross_chars, "ABH-");
    strcpy(default_intercross_chars,"ABCDH-");
    altered_chroms = FALSE;
    already_printed = FALSE;
    make_count_recs_etc();
}

void 
make_count_recs_etc (void) 
{
    int **x, *y;

   matrix(make_interval_genotype,MAX_LOCUS_GENOTYPES,MAX_LOCUS_GENOTYPES,int);
   array(left_genotype,MAX_INTERVAL_GENOTYPES,int);
   array(right_genotype,MAX_INTERVAL_GENOTYPES,int);

    x=make_interval_genotype;
    x[A][A]=AA; x[A][H]=AH; x[A][B]=AB;
    x[H][A]=HA; x[H][H]=HH; x[H][B]=HB;
    x[B][A]=BA; x[B][H]=BH; x[B][B]=BB;

    y=left_genotype;
    y[AA]= y[AB]= y[AH]= A;   
    y[HA]= y[HB]= y[HH]= H;   
    y[BA]= y[BB]= y[BH]= B; 

    y=right_genotype;
    y[AA]= y[BA]= y[HA]= A;   
    y[AH]= y[BH]= y[HH]= H;   
    y[AB]= y[BB]= y[HB]= B; 
}


void 
free_raw (void)  /* alloced by read_data - NEVER BEEN TESTED!! */
{
	unmatrix(raw.locus, raw.n_indivs, char);
	unmatrix(raw.locus_name, raw.n_loci, char);
	unmatrix(raw.trait, raw.n_indivs, real);
	unmatrix(raw.trait_name, raw.max_traits, char);
	unmatrix(raw.trait_eqn, raw.max_traits, char);	
	unmatrix(raw.left_cond_prob, raw.n_indivs, LOCUS_GENOTYPE_PROBS);
	unmatrix(raw.right_cond_prob, raw.n_indivs, LOCUS_GENOTYPE_PROBS);
	unarray(raw.map_dist, real);
	raw.n_loci= raw.n_traits= raw.n_indivs= raw.name_len= 0;
	raw.file[0]='\0';
}

int 
data_loaded (void) { return(raw.file[0]!='\0'); }



/********** THE DATA READER **********/

void 
getdataln ( /* get next nonblank/noncomment data file line */
    FILE *fp
)
{ do { fgetln(fp); BADDATA_line_num++; } while(nullstr(ln) || ln[0]=='#'); 
  BADDATA_ln= ln; }

#define DUMMY_LOCI 1

void 
read_data (FILE *fpa, FILE *fpb, FILE *fpm, char *temp)
{
 
    int k=0,v,l,i,num_chrom,num_loc=0,t_loc=0,j=0,loc;
    int name_len=80, checker, mapm_loci;
    int n_indivs, n_loci, n_traits, random_check1, random_check2, num_entries;
    real rf;
    char *name,*err;
    char num_of_chroms[TOKLEN+1];
    char all_str[5*SEQ_LEN], current_chrom[SEQ_LEN];
    char default_intercross_chars[10], default_backcross_chars[10];

    name = NULL;
    run {
	getdataln(fpa);
	sscanf(ln,"%*d %*d %d\n",&mapm_loci);
	mapm_loci *= 2;
	if(mapm_loci < 1) error("insufficient genetic data for analysis\n");

	array(map_space,mapm_loci,real);
	array(order,mapm_loci,int);
	array(name, NAME_LEN+1, char);

	strcpy(default_backcross_chars, "ABH-");
	strcpy(default_intercross_chars,"ABCDH-");

	getdataln(fpm);
	while(nstrcmp(ln, "*Chromosomes:", 13) != 0) getdataln(fpm);
	

	/***** DON'T BOTHER WITH THESE MAPM VARIABLES - JUST USE QTL DEFAULTS *****/
	/**********
	while(strcmp(ln, "*Print names: 0") != 0 && strcmp(ln,"*Print names: 1") != 0) 
	  getdataln(fpa);
	sscanf(ln,"%*s %*s %d",&print_names);
	fscanf(fpa,"%*s %*s %*d\n");
	fscanf(fpa,"%*s %*s %d\n",&print_maps);
	fscanf(fpa,"%*s %*s %d\n",&segregation_distortion);
	fscanf(fpa,"%*s %d\n",&mapnum);
	fscanf(fpa,"%*s %*s %d\n",&wizard_mode);

	getdataln(fpa); 
	getdataln(fpa);
	getdataln(fpa); 
	**********/

	stoken(&ln, sREQUIRED, num_of_chroms); /* get rid of "*Chromosomes:" */
	itoken(&ln, iREQUIRED, &num_chrom);
	if(num_chrom == 0)
	  error("no linkage maps saved...QTL cannot use this data yet");
	raw.n_chroms = num_chrom;
	array(raw.chrom_start,num_chrom,int);
	array(raw.chrom_n_loci,num_chrom,int);
	array(raw.original_locus,mapm_loci,int);

	/* Now I loop through the chromosomes in order to get all the
	   important loci, and the map distances between them */
	strcpy(all_str,"");
	dum_loc=0;
	for(i=0;i<num_chrom;i++) {
	    current_chrom[0]='\0';
	    getdataln(fpm); /* should be name and sequence line */
	    while(ln[0] != '*') getdataln(fpm); /* just to make sure */
	    stoken(&ln,sREQUIRED,name);  /* keep for named sequence */
	    itoken(&ln, iREQUIRED, &num_loc);
	    raw.chrom_n_loci[i] = num_loc;
	    if (!print_mapm_loci) {
		if (DUMMY_LOCI) sprintf(ps, "%d-%d", t_loc + 2, t_loc + num_loc + 1);
		else sprintf(ps, "%d-%d", t_loc + 1, t_loc + num_loc);
		strcat(all_str,ps);
		strcat(all_str," ");
	    }
	    if (DUMMY_LOCI) {
		dum_loc++;
		num_loc++;
	    }
	    if (!print_mapm_loci)
	      name_sequence(name,ps,&err); /* set chrom names */
	    
	    t_loc += num_loc;
	    getdataln(fpm); 
	    if (DUMMY_LOCI) {
		checker = -1;
		order[j] = -1;
		j++;
		for(v=0;v<num_loc-1;v++) {
		    if (nullstr(ln)) getdataln(fpm);
		    itoken(&ln, iREQUIRED, &loc); 
		    if (checker == -1) {
			raw.chrom_start[i] = loc;
			checker = 1;
		    }
		    order[j] = loc;
		    if (print_mapm_loci) {
			sprintf(ps, "%d ", loc + 1);
			strcat(current_chrom, ps);
		    } 
		    j++;
		}
		if (print_mapm_loci) {
		    name_sequence(name,current_chrom,&err);
		    strcat(all_str,&name[1]);
		    strcat(all_str," ");
		} 
		map_space[k] = .499; k++;
		for (l=0; l<num_loc-2; l++,k++) 
		  map_space[k] = read_map_distance(fpm/*,l*/);
		map_space[k] = 0.499; k++;
		getdataln(fpm);
	    }
	    else {
		checker = -1;
		for(v=0;v<num_loc;v++) {
		    if (nullstr(ln)) getdataln(fpm);
		    itoken(&ln, iREQUIRED, &loc); 
		    if (checker == -1) {
			raw.chrom_start[i] = loc;
			checker = 1;
		    }
		    order[j] = loc;
		    if (print_mapm_loci) {
		       sprintf(ps, "%d ", loc);
		       strcat(current_chrom, ps);
		       } 
		    j++;
		}
		if (print_mapm_loci) {
		    name_sequence(name,current_chrom,&err);
		    strcat(all_str,current_chrom);
		    strcat(all_str,current_chrom);
		} 
		if (DUMMY_LOCI) {
		    map_space[k] = .499; k++;
		}
		for (l=0; l<num_loc-1; l++,k++) 
		  map_space[k] = read_map_distance(fpm/*,l*/);
		map_space[k] = 0.499; k++;
		getdataln(fpm);
	    }
	}
	if (DUMMY_LOCI) {
	    order[j]= -1; /* Add dummy locus to end */
	    map_space[k] = 0.499; k++; t_loc++;
	}
	name_sequence("*all",all_str,&err);
	frewind(fpa);
	getdataln(fpa);getdataln(fpa);
	strcpy(raw.file,temp);
	raw.n_loci = t_loc;
	if (!(itoken(&ln,iREQUIRED,&random_check1) && random_check1>0 &&
	      itoken(&ln,iREQUIRED,&n_indivs) && n_indivs>0 && 
	      itoken(&ln,iREQUIRED,&n_loci) && n_loci>0))
	  send(BADDATA);
	if (random_check1%10 != F2VERSION)
	  error(OLD_FILENUM);
	raw.n_indivs = n_indivs;
	
	
	/*if(max_individuals<0) {
	  if (max_individuals>MAX_INDIVIDUALS) 
	  { why_(MAXEXCEED); send(BADDATA); }
	  else max_individuals= imaxf(n_indivs+1,0-max_individuals);
	  } else if (max_individuals>0 && n_indivs>max_individuals)
	  { why_(MAXEXCEED); send(BADDATA); }*/
	
	matrix(raw.locus_name,t_loc,name_len,char);
	matrix(raw.locus,n_indivs,t_loc,char);
	array(raw.map_dist,t_loc,real);
	matrix(raw.left_cond_prob, n_indivs, t_loc, 
	       LOCUS_GENOTYPE_PROBS); 
	matrix(raw.right_cond_prob, n_indivs, t_loc,
	       LOCUS_GENOTYPE_PROBS); 
	for(i=0;i<t_loc;i++) raw.map_dist[i] = map_space[i];
	read_map_locus(fpa,n_indivs,t_loc,order,n_loci);
	
	/* Now read the traits and QTL state variables from .trait file */

	getdataln(fpb);
	itoken(&ln, iREQUIRED, &random_check2);
	if (random_check2 != random_check1)
	  error(FILE_MISMATCH);
	getdataln(fpb);
	itoken(&ln, iREQUIRED, &n_traits);
	raw.filenumber = random_check2;
	raw.n_traits = n_traits;
	raw.max_traits= MAX_TRAITS(n_traits);
	matrix(raw.trait,n_indivs,raw.max_traits,real);
	matrix(raw.trait_name,raw.max_traits,name_len,char);
	matrix(raw.trait_eqn,raw.max_traits,TRAIT_EQN_LEN+1,char);
	for (j=0; j<n_traits; j++) 
	  read_quant_trait(fpb,j,n_indivs); 
	getdataln(fpb); /* Skips over line "#QTL only variables" */
	sscanf(ln,"%*s %*s %*s %d\n",&print_mapm_loci);
	fscanf(fpb,"%*s %*s %lf\n",&like_tolerance);
	fscanf(fpb,"%*s %*s %d\n", &brute_force);
	fscanf(fpb,"%*s %*s %d\n", &max_intervals);
	fscanf(fpb,"%*s %*s %*s %d\n", &max_continuous_vars);
	fscanf(fpb,"%*s %*s %d\n",&max_wiggles);
	fscanf(fpb,"%*s %*s %d\n",&max_compares);
	fscanf(fpb,"%*s %*s %d\n",&units);
	fscanf(fpb,"%*s %d\n",&num_chrom);

	/* Loop through the chromosomes, checking to ensure
	   no changes since QTL was last used */


	if (num_chrom == 0) {}
	else {  
	    if(num_chrom != raw.n_chroms)  {
	        altered_chroms_message();
	    }
	    if (DUMMY_LOCI) k = 1;
	    else k=0;
	    for(i=0;i<num_chrom;i++) {
		j = 0;
		getdataln(fpb);  
		sscanf(ln, "%*s %d", &num_loc);
		if(num_loc != raw.chrom_n_loci[i])
		  altered_chroms_message();
		getdataln(fpb);  
		checker = -1;
		for(v=0;v<num_loc;v++) {
		    if (nullstr(ln)) getdataln(fpb);
		    itoken(&ln, iREQUIRED, &loc); 
		    if (checker == -1) {
			if(loc != raw.chrom_start[i])
			  altered_chroms_message();
			checker = 1;
		    }
		}
		for(v=0;v<num_loc-1;v++,k++) {
		    if(nullstr(ln)) getdataln(fpb);
		    rtoken(&ln, rREQUIRED, &rf);
		    if(rf != map_space[k]) {
		        if(!(map_space[k] == .499 && rf >= .49)) 
			  altered_chroms_message();
		    }
		}
		if (DUMMY_LOCI) k++;
		k++;
		getdataln(fpb);getdataln(fpb);
	    }
	}


	
	fscanf(fpb,"%*s %*s %*s %d\n",&num_contexts);
	fscanf(fpb,"%*s %*s %d\n",&active_context);
	for(i = 0; i < num_contexts; i++) {
	    fscanf(fpb,"%*s %*d\n");
	    fscanf(fpb,"%*s %d\n", &trait);
	    if(!altered_chroms) {
		fscanf(fpb,"%*s %*s %d\n",&num_entries);
		load_table(context[i]->named_sequences,
			   fpb,INDEX_BY_NAME,num_entries);
		fscanf(fpb,"%*s %*s %d\n",&num_entries);
		load_table(context[i]->sequence_history,
			   fpb,INDEX_BY_NUMBER,num_entries);
		context[i]->seq_history_num= 
		  context[i]->sequence_history->next_entry_num;
	    } else {
		context[i]->seq_history_num = 0;
	    }
	}
    } on_exit {
	unarray(name, char);
	relay_messages;
    } return;
}


void 
save_traitfile (FILE *fp)
{
    int i,j,loci_tot,map_tot;
    
    run {
        sprintf(ps, "%d mapmaker trait data\n", raw.filenumber);
	fpr(fp);
	sprintf(ps, "%d\n", raw.n_traits);
	fpr(fp);

	for(i=0;i<raw.n_traits;i++) {
	    sprintf(ps, "*%-9s", raw.trait_name[i]);fpr(fp);

	    if (!nullstr(raw.trait_eqn[i])) 
	      { sprintf(ps, " =%s\n", raw.trait_eqn[i]);fpr(fp); }

	    for(j=0;j<raw.n_indivs;j++) {
		if(j % 5 == 0 && j != 0)
		  fprint(fp,"\n          ");
		sprintf(ps, "%12.3lf ", raw.trait[j][i]);fpr(fp);
	    }
	    fnl(fp);
	}

	fprint(fp,"#QTL only variables:\n");
	sprintf(ps, "*Print mapm loci: %d\n", print_mapm_loci);fpr(fp);
	sprintf(ps, "*Like tolerance: %lf\n", like_tolerance);fpr(fp);
	sprintf(ps, "*Brute force: %d\n", brute_force);fpr(fp);
	sprintf(ps, "*Max intervals: %d\n", max_intervals);fpr(fp);
	sprintf(ps, "*Max continuous vars: %d\n", max_continuous_vars);fpr(fp);
	sprintf(ps, "*Max wiggles: %d\n", max_wiggles);fpr(fp);
	sprintf(ps, "*Max compares: %d\n", max_compares);fpr(fp);
	sprintf(ps, "*Default units: %d\n", units); fpr(fp);
	sprintf(ps, "*Chromosomes: %d\n", raw.n_chroms);
	fpr(fp);
	loci_tot = 0;map_tot = 0;
	for (i=0;i<raw.n_chroms;i++) {
	    if (DUMMY_LOCI) {
		loci_tot++;
		map_tot++;
	    }
	    sprintf(ps, "chr%d %d\n", i + 1, raw.chrom_n_loci[i]); fpr(fp);
	    for(j=0;j<raw.chrom_n_loci[i];j++) {
		if(j % 18 == 0 && j != 0)
		  fprint(fp,"\n");
		sprintf(ps, "%d ", order[loci_tot]);fpr(fp);
		loci_tot++;
	    }
	    fnl(fp);
	    for(j=0;j<raw.chrom_n_loci[i]-1;j++) {
		if(j % 12 == 0 && j != 0)
		  fprint(fp,"\n");
		sprintf(ps, "%.4lf ", map_space[map_tot]);fpr(fp);
		map_tot++;
	    }
	    fnl(fp);
	    map_tot++;
	    fprint(fp,"0\n");
	    fprint(fp,"0\n");
	}
	sprintf(ps, "*Number of contexts: %d\n", num_contexts);
	fpr(fp);
	sprintf(ps, "*Active context: %d\n", active_context);
	fpr(fp);
	fnl(fp);
	for(i = 0; i < num_contexts; i++) {
	    sprintf(ps, "*Context %d\n", i + 1);  fpr(fp);
	    sprintf(ps, "*Trait: %d\n", trait); fpr(fp);
	    sprintf(ps, "*Named Sequences: %d\n",
                count_table_entries(context[i]->named_sequences));
	    fpr(fp);
	    
	    save_table(context[i]->named_sequences,fp,INDEX_BY_NAME);
	    sprintf(ps, "*Sequence history: %d\n",
                count_table_entries(context[i]->sequence_history));
	    fpr(fp);
	    save_table(context[i]->sequence_history,fp,INDEX_BY_NUMBER);
	}
    } on_exit {
	relay_messages;
    }
}


void 
read_map_locus (FILE *fp, int indivs, int t_loc, int *order, int n_loci)
/* could send BADDATA or IOERROR */
{
    
    char **temp_set,c, nam[TOKLEN+1];
    int j,k,i,t,num_of_terms = 0;
    matrix(temp_set,n_loci,80,char);
    for(i=0;i<t_loc;i++) {
	if (order[i] == -1) {
	    sprintf(ps, "ter%d", num_of_terms);
	    num_of_terms++;
	    strcpy(raw.locus_name[i],ps);
	    for(k=0;k<indivs;k++)
	      raw.locus[k][i] = '-';
	}
    }
    for (j=0;j<n_loci;j++) {
	getdataln(fp);
	while(ln[0] != '*') { getdataln(fp); }
	t=j;
	for(i=0;i<t_loc;i++) {
	    if (t == order[i]) {
		t = i;
		raw.original_locus[t] = order[i]+1; 
		self_delimiting = ptr_to("");
		if (!stoken(&ln, sREQUIRED, nam))
		  { send(BADDATA); }
		strcpy(raw.locus_name[t],&nam[1]);
		/* read the locus' data */
		for (k=0; k<indivs; k++) {
		    if (nullstr(ln)) getdataln(fp);
		    if (!parse_char(&ln,geno_chars,TRUE,&c)) 
		      { send(BADDATA); }
		    else if (c == 'B' && raw.data_type == BACKCROSS)
		      { raw.locus[k][t]= 'H'; }
		    else 
		      { raw.locus[k][t]= c; }
		}
		if (!nullstr(ln)) {	send(BADDATA); }
		break;
	    }
	}
    }
}
		

real read_map_distance(FILE *fp/*int num*/)
{
	real map_dist;
	if (nullstr(ln)) getdataln(fp);
	if (!rtoken(&ln, rREQUIRED, &map_dist))

/* REMOVED	    !rrange(&raw.map_dist[num],0.0,0.5))	*/
		{ why_("bad map distance"); send(BADDATA); }

/* THIS WAS A MAJOR KLUDGE - BUT IT IS TEMPORARY! */
        if (map_dist > 0.5) {
	    map_dist= unhaldane_cm(map_dist);
      }
/*	sprintf(ps,"%-3d %lf\n",num,raw.map_dist[num]); print(ps);   */
	rrange(&map_dist,MIN_REC_FRAC,MAX_REC_FRAC);	
    return(map_dist);
    }


void 
read_quant_trait (FILE *fp, int num, int indivs)
/* could send IOERROR or BADDATA */
{
    
    char tok[TOKLEN+1], c;
    int i;

    getdataln(fp); 
    /* read the locus name */
    if (!nstoken(&ln, sREQUIRED, tok, NAME_LEN)) { send(BADDATA); }
    strcpy(raw.trait_name[num],&tok[1]); raw.trait_eqn[num][0]='\0';
    if (parse_char(&ln,"=",TRUE,&c)) {
	if (nullstr(ln)) send(BADDATA); 
	nstrcpy(raw.trait_eqn[num],ln,TRAIT_EQN_LEN);
	ln = NULL;
    } 
    /* read its data */
    for (i=0; i<indivs; i++) {
	if (nullstr(ln)) getdataln(fp);
	if (!rtoken(&ln,rREQUIRED,&raw.trait[i][num])) {
	    if (!stoken(&ln,sREQUIRED,tok) || strcmp(tok,"-")) 
	      { send(BADDATA); }
	    else raw.trait[i][num]= MISSING_PHENO;
	}
    }
    if (!nullstr(ln)) { send(BADDATA); }
}
  

void
read_olddata (FILE *fp, char *path)
{
    int i, j, n_indivs, n_loci, n_traits;
    char tok[TOKLEN + 1];
    extern int nam_len;
	
    BADDATA_line_num= 0;
    if (!nullstr(raw.file)) free_raw();
    
    run {
	getdataln(fp); 
	if (!(itoken(&ln,iREQUIRED,&n_indivs) && n_indivs>0 && 
	      itoken(&ln,iREQUIRED,&n_loci) && n_loci>0 &&
	      itoken(&ln,iREQUIRED,&n_traits) && n_traits>0 &&
	      itoken(&ln,iREQUIRED,&nam_len) && 
	      nam_len>0 && nam_len<=TOKLEN)) 
	  { why_("bad header or header entries");
	    send(BADDATA); }
	
	
	raw.n_loci= n_loci;
	raw.n_traits= n_traits;
	raw.n_indivs= n_indivs;
	raw.name_len= nam_len;
	raw.left_cond_prob= NULL;
	raw.right_cond_prob=NULL;
	
	matrix(raw.locus, n_indivs, n_loci, char); 
	matrix(raw.locus_name, n_loci, nam_len, char);
	array(raw.map_dist, n_loci, real);
	raw.max_traits=MAX_TRAITS(raw.n_traits);
	matrix(raw.trait,n_indivs,raw.max_traits,real);
	matrix(raw.trait_name,raw.max_traits,nam_len,char);
	matrix(raw.left_cond_prob, n_indivs, n_loci, LOCUS_GENOTYPE_PROBS); 
	matrix(raw.right_cond_prob, n_indivs, n_loci,LOCUS_GENOTYPE_PROBS); 
	
	for (j=0; j<n_loci; j++) 
	  read_oldmap_locus(fp,j,n_indivs);
	for (j=0; j<n_traits; j++) 
	  read_oldquant_trait(fp,j,n_indivs);
	
	getdataln(fp); stoken(&ln,sREQUIRED,tok);
	if (!streq(tok,"map:")) 
	  { why_("expected 'MAP:'"); send(BADDATA); }
	
	for (i=0; i<n_loci-1; i++) read_oldmap_distance(fp,i);
	
	if (!nullstr(ln) || !end_of_text(fp)) 
	  { why_("data after logical end of file"); send(BADDATA); }
	
	close_file(fp);
	nstrcpy(raw.file,path,PATH_LENGTH);
	
    } when_aborting {
	free_raw();
	raw.file[0]='\0';
	relay;
    }
}


void 
read_oldmap_locus (FILE *fp, int num, int indivs)
/* could send BADDATA or IOERROR */
{
    char c;
    int i;

    /* read the locus name */
    getdataln(fp); 
    if (!field(&ln, nam_len, raw.locus_name[num])) 
        { why_("missing locus name"); send(BADDATA); }

    /* read the locus' data */
    for (i=0; i<indivs; i++) {
	if (nullstr(ln)) getdataln(fp);
	if (!parse_char(&ln,geno_chars,TRUE,&c)) 
	    { why_("expected a genotype character"); send(BADDATA); }
	else raw.locus[i][num]= c;
    }
    if (!nullstr(ln)) {	why_("extra data"); send(BADDATA); }
}
		

void 
read_oldquant_trait (FILE *fp, int num, int indivs)
/* could send IOERROR or BADDATA */
{
    char tok[TOKLEN+1];
    int i;
    
    getdataln(fp); 
    /* read the locus name */
    if (!field(&ln, nam_len, raw.trait_name[num])) 
	{ why_("bad locus name field"); send(BADDATA); }
    
    /* read its data */
    for (i=0; i<indivs; i++) {
	if (nullstr(ln)) getdataln(fp);
	if (!rtoken(&ln,rREQUIRED,&raw.trait[i][num])) {
	    if (!stoken(&ln,sREQUIRED,tok) || strcmp(tok,"-")) 
	      { why_("bad trait data"); send(BADDATA); }
	    else raw.trait[i][num]= MISSING_PHENO;
	}
    }
    if (!nullstr(ln)) { why_("extra data"); send(BADDATA); }
}

void 
read_oldmap_distance (FILE *fp, int num)
{
	int bar;

	if (nullstr(ln)) getdataln(fp);
	if (!rtoken(&ln, rREQUIRED, &raw.map_dist[num]))
/* REMOVED	    !rrange(&raw.map_dist[num],0.0,0.5))	*/
		{ why_("bad map distance"); send(BADDATA); }

/* THIS WAS A MAJOR KLUDGE - BUT IT IS TEMPORARY! */
        if (raw.map_dist[num] > 0.5) {
	  bar= (int) (raw.map_dist[num] + 0.499);
	  raw.map_dist[num]= unkosambi((real) bar);
	}
/*	sprintf(ps,"%-3d %lf\n",num,raw.map_dist[num]); print(ps);   */
	rrange(&raw.map_dist[num],MIN_REC_FRAC,MAX_REC_FRAC);	
}








/*** THE NEW AND IMPROVED RAW DATA CRUNCHER - USED AFTER LOADING DATA ****/

#define PROBS_NOT_EQ \
"*** WARNING: left and right probs unequal: indiv %d locus %d genotype %d\n"
#define PROBS_NOT_1 \
"*** WARNING: interval genotype probs not 1: indiv %d interval %d\n"

void 
crunch_data (void) /* side effects the raw data struct */
{
    int i, k, first, last;

    real r[3], l[3], right, left;
    INTERVAL_GENOTYPE_PROBS *indiv_probs;
    int geno;

    first=0; last=raw.n_loci-1;
    for (i=0; i<raw.n_indivs; i++) {

	/* get left conditioned probs */
	initial_prob(raw.left_cond_prob[i],first,raw.locus[i][first]);
	for (k=1; k<=last; k++) 
	  condition(raw.left_cond_prob[i],k,k-1,raw.locus[i][k],
		    raw.map_dist[k-1],LEFT);

	/* and now the right conditioned probs */
	for_locus_genotypes(raw.data_type,geno)
	  raw.right_cond_prob[i][last][geno] = 
	    pheno_given_geno(raw.locus[i][last],geno);
	for (k=last-1; k>=first; k--) 
	  condition(raw.right_cond_prob[i],k,k+1,raw.locus[i][k],
		    raw.map_dist[k],RIGHT);
    }

    /* DEBUGGING STUFF: We now check the L-R- conditioning */
    
    array(indiv_probs,2,INTERVAL_GENOTYPE_PROBS);
    for (i=0; i<raw.n_indivs; i++) 
      for (k=first; k<last-1; k++) {
	  r[A]=r[B]=r[H]=l[A]=l[B]=l[H]=left=right=0.0;

	  indiv_interval_probs(indiv_probs,0,i,k,k+1,raw.map_dist[k]);
	  for_interval_genotypes(raw.data_type,geno) 
	    l[right_genotype[geno]]+= indiv_probs[0][geno];

	  indiv_interval_probs(indiv_probs,1,i,k+1,k+2,raw.map_dist[k+1]);
	  for_interval_genotypes(raw.data_type,geno) 
	    r[left_genotype[geno]]+= indiv_probs[1][geno];

	  for_locus_genotypes(raw.data_type,geno) {
	      if (fabs(l[geno]-r[geno])>=0.01) 
		{ sprintf(ps, PROBS_NOT_EQ, i + 1, k + 2, geno); pr(); }
	      left+=l[geno]; right+=r[geno];
	  }
	  if (fabs(left-1.0)>=0.01) { sprintf(ps, PROBS_NOT_1, i + 1, k + 1); pr(); }
	  if (fabs(right-1.0)>=0.01) { sprintf(ps, PROBS_NOT_1, i + 1, k + 2); pr(); }

	 /* sprintf(ps,"Indiv = %d Locus = %d left=%lf right=%lf\n",i+1,k+2,
	     left,right); pr();*/
	  /* for_interval_genotypes(raw.data_type,foo) 
	    { sprintf(ps,"geno %d: left_int_prob=%lf right_int_prob=%lf\n",
		 foo,indiv_probs[0][foo],indiv_probs[1][foo]); pr(); }
		 */
	  /* for_locus_genotypes(raw.data_type,foo) 
	    { sprintf(ps,"geno %d lcp=%lf rcp=%lf\n",foo,
		 raw.left_cond_prob[i][k+1][foo],
		 raw.right_cond_prob[i][k+1][foo]); pr(); }
		 */
      }

}


void initial_prob(prob,locus,this_genotype)
LOCUS_GENOTYPE_PROBS *prob;
int locus;
char this_genotype;
{
    prob[locus][A]= prob[locus][B]= prob[locus][H]= 0.0;

    if (raw.data_type==BACKCROSS) 
      switch (this_genotype) {
	  case 'A': prob[locus][A]=1.0; break;
	  case 'B':                             
	  case 'H': prob[locus][H]=1.0; break;
	  case '-': prob[locus][A]=0.5; prob[locus][H]=0.5; break;
	  default:  send(CRASH);
      }
    else if (raw.data_type==INTERCROSS && !raw.f3) 
      switch (this_genotype) {
	  case 'A': prob[locus][A]=1.0; break;
	  case 'B': prob[locus][B]=1.0; break;
	  case 'H': prob[locus][H]=1.0; break;
	  case 'C': prob[locus][H]=2.0/3.0; prob[locus][B]=1.0/3.0; break;
	  case 'D': prob[locus][H]=2.0/3.0; prob[locus][A]=1.0/3.0; break;
	  case '-': prob[locus][A]=0.25; prob[locus][B]=0.25;
	            prob[locus][H]=0.5;  break;
	  default:  send(CRASH);
      }
    else if (raw.data_type==INTERCROSS && raw.f3) 
      switch (this_genotype) {
	  case 'A': prob[locus][A]=1.0; break;
	  case 'B': prob[locus][B]=1.0; break;
	  case 'H': prob[locus][H]=1.0; break;
	  case 'C': prob[locus][H]=0.4; prob[locus][B]=0.6; break;
	  case 'D': prob[locus][H]=0.4; prob[locus][A]=0.6; break;
	  case '-': prob[locus][A]=3.0/8.0; prob[locus][B]=3.0/8.0;
	            prob[locus][H]=0.25;    break;
	  default:  send(CRASH);
      }
    else send(CRASH);
}

real pheno_given_geno(observation,genotype)
char observation; /* eg the 'phenotype' */
int genotype;
{
    if (raw.data_type==BACKCROSS)
      switch (observation) {
	  case 'A': return((genotype==A) ? 1.0 : 0.0);
	  case 'B':
	  case 'H': return((genotype==H) ? 1.0 : 0.0);
	  case '-': return(1.0);
	  default:  send(CRASH);
      }
    else if (raw.data_type==INTERCROSS && !raw.f3)
      switch (observation) {
	  case 'A': return((genotype==A) ? 1.0 : 0.0);
	  case 'B': return((genotype==B) ? 1.0 : 0.0);
	  case 'H': return((genotype==H) ? 1.0 : 0.0);
	  case 'C': return((genotype!=A) ? 1.0 : 0.0);
	  case 'D': return((genotype!=B) ? 1.0 : 0.0);
	  case '-': return(1.0);
	  default:  send(CRASH);
      }
    else if (raw.data_type==INTERCROSS && raw.f3) /* no different??? */
      switch (observation) {
	  case 'A': return((genotype!=B) ? 1.0 : 0.0);
	  case 'B': return((genotype!=A) ? 1.0 : 0.0);
	  case 'H': return((genotype==H) ? 1.0 : 0.0);
	  case 'C': return((genotype!=A) ? 1.0 : 0.0);
	  case 'D': return((genotype!=B) ? 1.0 : 0.0);
	  case '-': return(1.0);
	  default:  send(CRASH);
      }
    send(CRASH);
    abort();
}


void condition(prob,locus,prev_locus,observation,rec_frac,side)
LOCUS_GENOTYPE_PROBS *prob;
int locus, prev_locus;
char observation;
real rec_frac;
int side;
{
    int i, geno_was, geno_is, geno1,geno2;
    real total;

    for (i=0; i<MAX_LOCUS_GENOTYPES; i++) prob[locus][i]=0.0;
    total= 0.0;
    
    for_locus_genotypes(raw.data_type,geno_was)
      for_locus_genotypes(raw.data_type,geno_is) {
	  geno1 = (side==LEFT) ? geno_was:geno_is;
	  geno2 = (side==LEFT) ? geno_is :geno_was;
	  prob[locus][geno_is]+= 
	    prob[prev_locus][geno_was] *
	    transition_prob(raw.data_type,geno1,geno2,rec_frac);
      }
    for_locus_genotypes(raw.data_type,geno_is) 
      total += prob[locus][geno_is] *= pheno_given_geno(observation,geno_is);
    for_locus_genotypes(raw.data_type,geno_is) prob[locus][geno_is]/=total;
}		  
		
	
real transition_prob(data_type,geno_was,geno_is,theta)
int data_type, geno_was, geno_is;
real theta;
{
    real sum;
    real C2, D2, E2, F2, G2, C3, D3, E3, F3, G3;

    if (data_type==BACKCROSS) {
	if (geno_was==geno_is) return(1.0-theta);
	else return(theta);

    } else if (data_type==INTERCROSS && !raw.f3) {
	switch(geno_was) {
	    case A: switch(geno_is) {
		case A: return(sq(1.0-theta));
		case B: return(sq(theta)); 
		case H: return(2.0*theta*(1.0-theta));
	    }
	    break;
	    case B: switch(geno_is) {
		case A: return(sq(theta)); 
		case B: return(sq(1.0-theta));
		case H: return(2.0*theta*(1.0-theta)); 
	    }
	    break;
	    case H: switch(geno_is) {
		case A: return(theta*(1.0-theta)); 
		case B: return(theta*(1.0-theta));
		case H: return(sq(theta)+sq(1.0-theta)); 
	    }
	    break;
	}

    } else if (data_type==INTERCROSS && raw.f3) { /* haldane & waddington */
	/* OK BONEHEAD - Why did I put I's in here??? */
	C2= 0.5*sq(1.0-theta); 
	D2= 0.5*sq(theta);
	E2= theta*(1.0-theta);
	F2= sq(1.0-theta);
	G2= sq(theta);
	
	C3= C2 + 0.5*E2 + 0.25*sq(1.0-theta)*F2 + 0.25*sq(theta)*G2;
	D3= D2 + 0.5*E2 + 0.25*sq(theta)*F2 + 0.25*sq(1.0-theta)*G2;
	E3= 0.5*E2 + 0.5*theta*(1.0-theta)*(F2 + G2);
	F3= 0.5*sq(1.0-theta)*F2 + 0.5*sq(theta)*G2;
	G3= 0.5*sq(theta)*F2 + 0.5*sq(1.0-theta)*G2;

	switch(geno_was) {
	    case A: 
	        sum= C3 + D3 + E3;
		switch(geno_is) {
		    case A:  return(C3/sum);
		    case B:  return(D3/sum);
		    case H:  case I: return(0.5*(E3/sum));
		    default: send(CRASH); break;
		}
		break;
	    case B: 
	        sum= C3 + D3 + E3;
		switch(geno_is) {
		    case B:  return(C3/sum);
		    case A:  return(D3/sum);
		    case H:  case I: return(0.5*(E3/sum));
		    default: send(CRASH); break;
		}
            break;
		case H: case I:
		sum= 2.0*E3 + F3 + G3;
		switch(geno_is) {
		    case A:  return(E3/sum);
		    case B:  return(E3/sum);
		    case H: case I: 
		       return( (geno_was==geno_is ? F3:G3)/sum);
		    default: send(CRASH); break;
		}
		break;
	    default: send(CRASH); break;
	}

#ifdef ANOTHER_SHOT
    } else if (data_type==INTERCROSS && raw.f3) {  /* new look */
	x= sq(1.0-theta); y=theta*(1.0-theta); z=sq(theta); 
	switch(geno_was) {
	    case A: 
		switch(geno_is) {
		    case A:  return((1.0/3.0)*(sq(x) + sq(z) + 2.0*y + 2.0*x));
		    case B:  return((2.0/3.0)*(x*z + y + z));
		    case H:  return((2.0/3.0)*(y*x + y*z + y));
		    default: send(CRASH);
		}
	    case B:
		switch (geno_is) {
		    case B:  return((1.0/3.0)*(sq(x) + sq(z) + 2.0*y + 2.0*x));
		    case A:  return((2.0/3.0)*(x*z + y + z));
		    case H:  return((2.0/3.0)*(y*x + y*z + y));
		    default: send(CRASH);
		}
	    case H:
		switch(geno_is) {
		    case A:  return(y*(1.0+x+z));
		    case B:  return(y*(1.0+x+z));
		    case H:  return(x*(x+z) + z*(x+z));
		    default: send(CRASH);
		}
	    default: send(CRASH);
	}
#endif
#ifdef YET_ANOTHER
    } else if (data_type==INTERCROSS && raw.f3) {  /* From Transition Matrix */
    INTERVAL_GENOTYPE_PROBS prob;
	x= sq(1.0-theta); y=theta*(1.0-theta); z=sq(theta);
	sum=0.0;  
	switch(geno_was) {
	    case A: 
	    	sum+= prob[AA]= 0.25*x + 0.25*y + 0.125*x*x +0.125*z*z;
		sum+= prob[AH]= 0.25*y + 0.25*x*y + 0.25*z*y;
		sum+= prob[AB]= 0.25*y + 0.25*z + 0.25*x*z;
		switch(geno_is) {
		    case A:  return(prob[AA]/sum);
		    case H:  return(prob[AH]/sum);
		    case B:  return(prob[AB]/sum);
		    default: send(CRASH);
		}
	    case H:
	    	sum+= prob[HA]= 0.25*y + 0.25*x*y + 0.25*z*y;
		sum+= prob[HH]= 0.25*x*x + 0.25*z*z + 0.25*x*z + 0.25*z*x;
		sum+= prob[HB]= 0.25*y + 0.25*x*y + 0.25*z*y;
		switch(geno_is) {
		    case A:  return(prob[HA]/sum);
		    case H:  return(prob[HH]/sum);
		    case B:  return(prob[HB]/sum);
		    default: send(CRASH);
		}
	    case B:
	    	sum+= prob[BB]= 0.25*x + 0.25*y + 0.125*x*x +0.125*z*z;
		sum+= prob[BH]= 0.25*y + 0.25*x*y + 0.25*z*y;
		sum+= prob[BA]= 0.25*y + 0.25*z + 0.25*x*z;
		switch(geno_is) {
		    case A:  return(prob[BA]/sum);
		    case H:  return(prob[BH]/sum);
		    case B:  return(prob[BB]/sum);
		    default: send(CRASH);
		}
	    default: send(CRASH);
	}
#endif 
    }

    send(CRASH);
    abort();
}


real map_length(left,right) 	/* return rf dist from left to right */
int left, right;		/* Assumes check_interval() has succeeded */
{
    real total_cm, total_rf;
    int i;
    
	
/*  if (zero_interval(left,right)) return(MIN_REC_FRAC);  */
    if (right==INF_LOCUS) return(MAX_REC_FRAC); 
    if (left+1==right) return(raw.map_dist[left]);
    
    for (i=left, total_cm=0.0; i<right; i++) {
	if (raw.map_dist[i]>MAX_REC_FRAC) return(MAX_REC_FRAC);
	else total_cm+= haldane_cm(raw.map_dist[i]);
	if (total_cm>MAX_CM) return(MAX_REC_FRAC);
    }
    total_rf= unhaldane_cm(total_cm); 
    rrange(&total_rf,MIN_REC_FRAC,MAX_REC_FRAC);
    return(total_rf);
}


void indiv_interval_probs(prob,data_indiv_num,raw_indiv_num,left_locus_num,
  right_locus_num,theta)
INTERVAL_GENOTYPE_PROBS *prob; /* [num][geno-code] => real */
int data_indiv_num, raw_indiv_num;
int left_locus_num, right_locus_num;
real theta; /* provided as an argument for efficiency */
{
    real sum;
    int interval_geno, left, right;

    sum= 0.0;
    for_interval_genotypes(raw.data_type,interval_geno) {
	left= left_genotype[interval_geno];
	right=right_genotype[interval_geno];
	sum+= prob[data_indiv_num][interval_geno]= 
	  raw.left_cond_prob[raw_indiv_num][left_locus_num][left] *
	    transition_prob(raw.data_type,left,right,theta) *
	      raw.right_cond_prob[raw_indiv_num][right_locus_num][right];
    }

    for_interval_genotypes(raw.data_type,interval_geno) {
	if(sum == 0.0) print("sum = 0\n");
	prob[data_indiv_num][interval_geno]/= sum;
    }
}


real apriori_prob(data_type,geno)
int data_type, geno;
{
    if (data_type==BACKCROSS) return(0.5);

    else switch(geno) {
	case A: return(0.25);
	case B: return(0.25);
	case H: return(0.50);
    }

    abort();
}


void 
altered_chroms_message (void)
{
if(!already_printed) {
    print("The linkage map data in the '.data' file has been altered since its\n");
    print("last usage. The saved scan results and saved names are obsolete and\n");
    print("will not be loaded.\n");
    already_printed = TRUE;
}
altered_chroms = TRUE;
}

