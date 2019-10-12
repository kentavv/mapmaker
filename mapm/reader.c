/******************************************************************************

 #####   ######    ##    #####   ######  #####            ####
 #    #  #        #  #   #    #  #       #    #          #    #
 #    #  #####   #    #  #    #  #####   #    #          #
 #####   #       ######  #    #  #       #####    ###    #
 #   #   #       #    #  #    #  #       #   #    ###    #    #
 #    #  ######  #    #  #####   ######  #    #   ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#define INC_LIB
#define INC_SHELL
#define INC_EQN
#include "mapm.h"

int symbol_value (int chr, char *symb);

RAW_DATA raw;
int original_markers;

#define end_of_text(fp) feof(fp)
int  data_line; /* Remember to reset when you open the data file! */

/* for the BADDATA message */
int BADDATA_line_num;
char *BADDATA_text;
char *BADDATA_reason;
char *badstr;

int read_data_file_header();
int read_raw_file_header();
FILE *start_save_to_file();
void finish_save_to_file();
bool read_magic_number();
void write_magic_number();
int  new_magic_number();

void new_read_f2_data();
void allocate_f2_data();
void free_f2_data();
void free_traits();
void new_read_f2_locus();
void write_f2_data();

void new_read_f2_raw();
void read_raw_f2_locus();
void read_raw_trait();
void write_traits();

void add_to_seg_dist(char c, int locus);
void scale_seg_dist(int locus);

bool uppercase_genotypes; /* set by read_raw_header() for read_raw_f2_locus()*/

real **trait=NULL;
char **trait_name=NULL;
int num_traits=0;

#define MISSING_PHENO -100000.0

#define F2_ONLY \
  "only experimental data may be loaded into this version of MAPMAKER"

#define LOADING_FROM     "loading %s data from file '%s'... "
#define SAVING_TO        "saving %s data in file '%s'... "
#define NOT_SAVING_TO    "can't write file '%s'... "
#define NOT_LOADING_FROM "%s data in file '%s' is old... not loading\n"

#define OLD_FILENUM "\
This file was prepared using an older version of MAPMAKER.\n\
You will need to re-prepare it using the 'prepare data' command."


void 
data_init (void)
{
    strcpy(raw.filename,"");
    raw.data_type= NO_DATA;
    setmsg(BADDATA,"error in data file",sender,strmsg_default);
    array(BADDATA_text,MAXLINE+1,char); array(BADDATA_reason,MAXLINE+1,char);
    array(badstr,LINE+1,char);
}


void 
getdataln ( /* get next nonblank,noncomment data file line */
    FILE *fp
)
{ do { fgetln(fp); data_line++; } while (nullstr(ln)||ln[0]=='#'); }

  
void 
baddata ( /* send data reading error message */
    char *reason /* NOTE: should be a constant or a global */
)
{ nstrcpy(BADDATA_reason,reason,MAXLINE); BADDATA_line_num=data_line;
  nstrcpy(BADDATA_text,ln,MAXLINE); send(BADDATA); }


int 
data_loaded (void)
{ return (raw.data_type!=NO_DATA); }


char *
data_info (bool add_nums)
{
    char *str1= get_temp_string(), *str2= get_temp_string();

#ifdef HAVE_CEPH
    if (raw.data_type==CEPH) {
	if (add_nums) 
	 sprintf(str,"CEPH data (%d families, %d loci)",
	     raw.data.ceph.num_families,raw.num_markers);
       else 
	  strcpy(str,"CEPH data");
       return(str);
   } 
#endif /* else */

    switch (raw.data.f2.cross_type) {
	case F2_INTERCROSS: strcpy(str1,"F2 intercross data "); break;
	case F2_BACKCROSS:  strcpy(str1,"F2 backcross data "); break;
	case RI_SIB:        strcpy(str1,"RI (sib-mating) data "); break;
	case RI_SELF:       strcpy(str1,"RI (selfing) data "); break;
	case F3_SELF:       strcpy(str1,"F3 intercross (selfing) data"); break;
	default:            send(CRASH);
    }
    
    if (!add_nums) return(str1);
    else { sprintf(str2,"%s (%d individuals, %d loci)",str1,
		   raw.data.f2.num_indivs,raw.num_markers);
	   return(str2); 
    }
}


FILE *
start_save_to_file (char *name, char *ext, char *type, bool *exists)
{
    FILE *fp;
    char tmpname[PATH_LENGTH+1];

    make_filename(name,FORCE_EXTENSION,ext);
    strcpy(tmpname,name); make_filename(tmpname,FORCE_EXTENSION,TEMP_EXT);
   
    fp=NULL; *exists=FALSE;
    run {
	fp=open_file(name,READ);
	*exists=TRUE;
	close_file(fp);
    } trapping(CANTOPEN) {} 

    fp=NULL;
    run {
	fp=open_file(tmpname,WRITE);
	sprintf(ps, SAVING_TO, type, name); pr(); flush();
    } except_when(CANTOPEN) {
	sprintf(ps, NOT_SAVING_TO, name); pr(); nl();
	abort_command();
    }

    return(fp);
}


void 
finish_save_to_file (char *name, char *oldext, bool exists)
{
    char tmpname[PATH_LENGTH+1], oldname[PATH_LENGTH+1];
    strcpy(tmpname,name); make_filename(tmpname,FORCE_EXTENSION,TEMP_EXT);
    strcpy(oldname,name); make_filename(oldname,FORCE_EXTENSION,oldext);

    if (exists) { rename_file(name,oldname); }
    rename_file(tmpname,name);
    print("ok\n");
}
    

/**** Do Real Work Now ****/

void 
do_load_data (FILE *fp, char *filename, bool israw)
{
    char name[PATH_LENGTH+1], symbols[7], type[TOKLEN+1];
    FILE *fp2;
    int ntype;
    if (data_loaded()) send(CRASH);
    num_traits=0; trait=NULL; trait_name=NULL;

    sprintf(ps, "%sing data from file '%s'... ", (israw ? "prepar" : "load"), filename);
    pr(); flush();

    if (!israw) ntype=read_data_file_header(fp,filename); /* sets raw.stuff */
      else ntype=read_raw_file_header(fp,filename,symbols); /* ditto */

    if (ntype==F2) { /* read_f2_frob allocates the matrix first */
	if (!israw) new_read_f2_data(fp); else new_read_f2_raw(fp,symbols);
    }
#ifndef HAVE_CEPH
    else baddata(F2_ONLY);
#else
    else read_ceph_data(fp);
    allocate_sex_choose();
#endif

    print("ok\n  "); print(data_info(TRUE)); /*print(".");*/ flush();
    strcpy(raw.filename,filename); raw.data_type=ntype; /* got it */
    strcpy(name,filename);

    reset_state();
    allocate_order_data(raw.num_markers);   /* print("."); flush(); */
    allocate_seq_stuff(raw.num_markers);       print("."); flush();
    allocate_mapping_data(raw.num_markers);    print("."); flush();
    allocate_two_pt(raw.num_markers);       /* print("."); flush(); */
    allocate_three_pt(raw.num_markers);        print("."); flush();
    allocate_hmm_temps(raw.num_markers,raw.data.f2.num_indivs,
		raw.data.f2.cross_type);
    print(" ok\n");

    make_filename(name,FORCE_EXTENSION,MAPS_EXT); strcpy(type,"map");
    run {
	fp2=NULL;
        fp2=open_file(name,READ);
	if (read_magic_number(fp2,type)) {
	    sprintf(ps, LOADING_FROM, type, name); pr(); flush();
	    read_order_data(fp2);
	    read_mapping_data(fp2);
	    read_status(fp2);
	    print("ok\n");
	} else { sprintf(ps, NOT_LOADING_FROM, type, name); pr(); }
	close_file(fp2);
    } except_when(CANTOPEN) { } /* need a better handler */

    make_filename(name,FORCE_EXTENSION,TWO_EXT); strcpy(type,"two-point");
    run {
	fp2=NULL;
	fp2=open_file(name,READ);
	if (read_magic_number(fp2,type)) {
	    sprintf(ps, LOADING_FROM, type, name); pr(); flush();
	    read_two_pt(fp2);
	    print("ok\n");
	} else { sprintf(ps, NOT_LOADING_FROM, type, name); pr(); }
	close_file(fp2);
    } except_when(CANTOPEN) { } /* need a better handler */

    make_filename(name,FORCE_EXTENSION,THREE_EXT); strcpy(type,"three-point");
    run {
	fp2=NULL;
        fp2=open_file(name,READ);
        if (read_magic_number(fp2,type)) {
	    sprintf(ps, LOADING_FROM, type, name); pr(); flush();
	    read_three_pt(fp2);
	    print("ok\n");
	} else { sprintf(ps, NOT_LOADING_FROM, type, name); pr(); }
	close_file(fp2);
    } except_when(CANTOPEN) { } /* need a better handler */
}


void 
do_unload_data (void)
{
    if (raw.data_type==F2)
      free_f2_data(raw.num_markers/*,raw.data.f2.num_indivs*/);
#ifdef HAVE_CEPH
      else free_ceph_data(raw.num_markers,raw.data.ceph,num_families);
    free_sex_choose();
#endif
    strcpy(raw.filename,""); raw.data_type=NO_DATA; /* use cross_type below */
    undo_state();
    free_order_data(raw.num_markers);
    free_mapping_data(raw.num_markers);
    free_seq_stuff(/*raw.num_markers*/);
    free_two_pt(raw.num_markers);
    free_three_pt(raw.num_markers);
    free_hmm_temps(raw.num_markers,raw.data.f2.num_indivs,
		   raw.data.f2.cross_type);
}


void 
do_save_data (char *base_name, bool save_genos_too)
{
    FILE *fp=NULL;
    char name[PATH_LENGTH+1];
    bool exists;

    run {
	strcpy(name,base_name);
	if (save_genos_too) {
	    fp= start_save_to_file(name,DATA_EXT,"genotype",&exists);
	    write_f2_data(fp); /* deals with header and magic number */
	    close_file(fp); fp=NULL;
	    finish_save_to_file(name,DATA_OLD,exists);
	}
	fp= start_save_to_file(name,MAPS_EXT,"map",&exists);
	write_magic_number(fp,"map");
	write_order_data(fp);
	write_mapping_data(fp);
	write_status(fp);
	close_file(fp); fp=NULL;
	finish_save_to_file(name,MAPS_OLD,exists);

	if (two_pt_touched) {
	    fp= start_save_to_file(name,TWO_EXT,"two-point",&exists);
	    write_magic_number(fp,"two-point");
	    write_two_pt(fp);
	    close_file(fp); fp=NULL;
	    finish_save_to_file(name,TWO_OLD,exists);
	}
	if (three_pt_touched) {
	    fp= start_save_to_file(name,THREE_EXT,"three-point",&exists);
	    write_magic_number(fp,"three-point");
	    write_three_pt(fp);
	    close_file(fp); fp=NULL;
	    finish_save_to_file(name,THREE_OLD,exists);
	}
	if (num_traits>0) {
	    fp= start_save_to_file(name,TRAIT_EXT,"traits",&exists);
	    write_magic_number(fp,"trait");
	    write_traits(fp);
	    close_file(fp); fp=NULL;
	    finish_save_to_file(name,TRAIT_OLD,exists);
	    free_traits();
	}
	two_pt_touched= FALSE;
	three_pt_touched= FALSE;

    } on_error {
	if (fp!=NULL && msg!=CANTCLOSE) close_file(fp);
	relay;
    }
}


/**** Lower level stuff ****/

int 
read_data_file_header (FILE *fp, char *filename)
{
    int type=NO_DATA;
    if (data_loaded() || fp==NULL) send(CRASH);
    
    getdataln(fp);

    if (streq(ln,"prepared data f2 intercross")) {
        type=F2; raw.data.f2.cross_type=F2_INTERCROSS;
    } else if (streq(ln,"prepared data f2 backcross")) {
        type=F2; raw.data.f2.cross_type=F2_BACKCROSS;
    } else if (streq(ln,"prepared data f2 ri-sib")) {
        type=F2; raw.data.f2.cross_type=RI_SIB;
    } else if (streq(ln,"prepared data f2 ri-self")) {
        type=F2; raw.data.f2.cross_type=RI_SELF;
    } else if (streq(ln,"prepared data f3")) {
        type=F2; raw.data.f2.cross_type=F3_SELF;
#ifdef HAVE_CEPH
    } else if (streq(ln,"data type ceph")) {
        type=CEPH;
#endif
    } else baddata("can't recognize the data type");
    if (type!=F2) send(CRASH);

    fscanf(fp,"%d %d %d",&raw.filenumber,&raw.data.f2.num_indivs,
	   &raw.num_markers);
    if (raw.filenumber%10 != F2VERSION) baddata(OLD_FILENUM);
    return(type);
}


#define OLD_CHROM_WARN \
"warning: chromosome definitions in the raw file are ignored by this version\n\
of MAPMAKER. Type 'help prepare data' for instructions on how to accomplish\n\
the same thing using a '.prep' file."
#define SYM_SYMBOLS "found 'symbols' without genotype symbols"
#define SYM_FORMAT \
"the format of the genotype symbol definitions is wrong\n\
Type 'help data formats' for instructions."
#define SYM_REPEAT "a genotype symbol is defined more than once" 
#define OK_SYMBOLS \
  "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890+-"

int 
read_raw_file_header (FILE *fp, char *filename, char *symbol)
{
    int type, n_indivs, n_loci, n_chroms, n_traits, i, n=0, first;
    char str[TOKLEN+1], *dflt=NULL, user[2];
    if (data_loaded() || fp==NULL) send(CRASH);

    getdataln(fp); while(nullstr(ln)); crunch(ln);

    if (streq(ln,"data type f2 intercross")) {
        type=F2; raw.data.f2.cross_type=F2_INTERCROSS; dflt=ptr_to("-AHBCD");
    } else if (streq(ln,"data type f3 intercross")) {
        type=F2; raw.data.f2.cross_type=F3_SELF; dflt=ptr_to("-AHBCD");
    } else if (streq(ln,"data type f2 backcross")) {
        type=F2; raw.data.f2.cross_type=F2_BACKCROSS; dflt=ptr_to("-AH");
    } else if (streq(ln,"data type ri sib")) {
        type=F2; raw.data.f2.cross_type=RI_SIB; dflt=ptr_to("-AB");
    } else if (streq(ln,"data type ri self")) {
	type=F2; raw.data.f2.cross_type=RI_SELF; dflt=ptr_to("-AB");
#ifdef HAVE_CEPH
    } else if (streq(ln,"data type ceph")) {
        raw.data_type=CEPH;
#endif
    } else baddata("bad header first line - can't recognize data type");
    if (type!=F2) send(CRASH);

    do getdataln(fp); while(nullstr(ln)); /* NOT crunched */
    if (!(itoken(&ln,iREQUIRED,&n_indivs) && n_indivs > 0 &&
	  itoken(&ln,iREQUIRED,&n_loci) && n_loci > 0 &&
	  itoken(&ln,iREQUIRED,&n_traits) && n_traits >= 0))
      error("bad #individuals, #loci or #traits");
    
    uppercase_genotypes=TRUE; first=TRUE; user[1]='\0';
    if (len(dflt)==3) strcpy(symbol,"   "); else strcpy(symbol,"      ");
    if (!nullstr(ln)) { /* symbol definitions, old chrom#, case thing */
	if (itoken(&ln,iREQUIRED,&n_chroms))
	  if (n_chroms>0) { print(OLD_CHROM_WARN); } 
	/* BIG KLUDGE BECAUSE SELF_DELIMITING INCLUDES '=' */
	while (!nullstr(ln) && sscanf(ln,"%s",str)==1) {
	    while (!white(ln[0]) && ln[0]!='\0') ln++;
	    while (white(ln[0])  && ln[0]!='\0') ln++;
	    if (first && xstreq(str,"case")) 
	      { uppercase_genotypes=FALSE; continue; }
	    if ((first && xstreq(str,"symbols")) || xstreq(str,"symbols:")) {
		if (nullstr(ln)) baddata(SYM_SYMBOLS);
		  else continue;
	    }
	    first=FALSE;
	    if (len(str)!=3 || str[1]!='=' || !strin(OK_SYMBOLS,str[0]) ||
		!strin(dflt,str[2])) baddata(SYM_FORMAT);
	    switch (str[2]) {
		case '-':           n=0; break;
		case 'A': case 'a': n=1; break;
		case 'H': case 'h': n=2; break;
		case 'B': case 'b': n=3; break;
		case 'C': case 'c': n=4; break;
		case 'D': case 'd': n=5; break;
		default: baddata(SYM_FORMAT);
	    }
	    user[0]=str[0]; if (uppercase_genotypes) uppercase(user);

	    /* trick for RIs, want to use -AB as default but convert it 
	       to -AH surreptitiously */
	    if (n==3 && (raw.data.f2.cross_type==RI_SIB ||
			 raw.data.f2.cross_type==RI_SELF)) n=2;

	    if (strin(symbol,*user) || symbol[n]!=' ') baddata(SYM_REPEAT);
	    symbol[n]=user[0];
	  }
      }
    for (i=0; symbol[i]!='\0'; i++)
      if (symbol[i]!=' ') continue;
      else if (strin(symbol,dflt[i])) baddata(SYM_REPEAT);
      else symbol[i]=dflt[i];

    raw.num_markers= n_loci;
    raw.data.f2.num_indivs= n_indivs;
    raw.filenumber= new_magic_number();
    num_traits= n_traits;
    return(type);
}


bool 
read_magic_number (
    FILE *fp,
    char *type /* no spaces allowed, must parse with %s */
)
{
    int magic_number;
    char id[TOKLEN+1];

    getdataln(fp);
    if (sscanf(ln,"%d %*s %s %*s",&magic_number,id)!=2) return(FALSE); 
    else if (magic_number!=raw.filenumber) return(FALSE);
    else if (!streq(id,type)) return(FALSE);
    else return(TRUE);
}


void 
write_magic_number (FILE *fp, char *type)
{ fprintf(fp,"%d mapmaker %s data\n",raw.filenumber,type); }


int 
new_magic_number (void)
{ return(((int)(randnum()*3275.0))*10 + F2VERSION); }



/**** Actually Get the Genotypes ****/

void 
allocate_f2_data (int n_loci, int n_indivs)
{
    matrix(raw.locus_name,n_loci,NAME_LEN+1,char);
    matrix(raw.data.f2.allele,n_loci,n_indivs,char);
    matrix(raw.data.f2.allelic_distribution,n_loci,4,real);
    if (num_traits>0) {
	matrix(trait,num_traits, raw.data.f2.num_indivs,real);
	matrix(trait_name,num_traits,NAME_LEN+1,char);
    }
}


void free_f2_data(int n_loci/*,int n_indivs*/)
{
    unmatrix(raw.locus_name, n_loci, char);
    unmatrix(raw.data.f2.allele, n_loci, char);
    unmatrix(raw.data.f2.allelic_distribution, n_loci, real);
    raw.data_type=NO_DATA; raw.filename[0]='\0';

    /* Just for sanity, although if this happens we get a big memory leak */
    num_traits=0; 
}


void 
free_traits (void)
{
    if (num_traits>0) {
	unmatrix(trait,num_traits,real);
	unmatrix(trait_name,num_traits,char);
	num_traits=0;
    }
}

void 
new_read_f2_data (
    FILE *fp /* sends BADDATA, or MALLOC, INTERRUPT etc. if an error */
)
{
    int j;
    run {
	allocate_f2_data(raw.num_markers,raw.data.f2.num_indivs);
	for (j=0; j<raw.num_markers; j++) new_read_f2_locus(fp,j);
    } on_error {
	free_f2_data(raw.num_markers/*,raw.data.f2.num_indivs*/);
	if (msg==ENDOFILE) baddata("unexpected end of file");
	relay;
    }
}	


void 
new_read_f2_locus (FILE *fp, int locus_num)
{
    int i;
    char c, name[NAME_LEN+2];

    getdataln(fp); /* Must be at the start of a line */
  
    /* This should never fail, as ln should never be empty! */
    nstoken(&ln,sREQUIRED,name,NAME_LEN+1);/* truncates */
    if (name[0]!='*' || len(name)<2) baddata("expected *name");
    strcpy(raw.locus_name[locus_num],name+1);
	
    raw.data.f2.allelic_distribution[locus_num][0]= 0.0;
    raw.data.f2.allelic_distribution[locus_num][1]= 0.0;
    raw.data.f2.allelic_distribution[locus_num][2]= 0.0;
    raw.data.f2.allelic_distribution[locus_num][3]= 0.0;

    for (i=0; i<raw.data.f2.num_indivs; i++) {
        if (nullstr(ln)) getdataln(fp);
	if (!parse_char(&ln,"-ABCDH",SKIPWHITE,&c))
	  baddata("illegal genotype character");
	raw.data.f2.allele[locus_num][i]= c;
	/* add_to_seg_dist(c,locus_num); */
    }
    /* scale_seg_dist(locus_num); */
}


void 
write_f2_data (FILE *fp)
{
    int i, j;
    char header[MAXLINE+1];

    if (raw.data.f2.cross_type==F2_INTERCROSS) 
      fprint(fp,"prepared data f2 intercross\n");
    else if (raw.data.f2.cross_type==F2_BACKCROSS) 
      fprint(fp,"prepared data f2 backcross\n");
    else if (raw.data.f2.cross_type==RI_SIB) 
      fprint(fp,"prepared data f2 ri-sib\n");
    else if (raw.data.f2.cross_type==RI_SELF) 
      fprint(fp,"prepared data f2 ri-self\n");
    else if (raw.data.f2.cross_type==F3_SELF)
      fprint(fp,"prepared data f3\n");
    
    sprintf(header, "%d %d %d\n", raw.filenumber, raw.data.f2.num_indivs,
            raw.num_markers); fprint(fp,header);
    
    for (i=0; i<raw.num_markers; i++) {
	sprintf(ps, "*%-10s   ", raw.locus_name[i]);
	fpr(fp);
	for (j=0; j<raw.data.f2.num_indivs; j++) {
	    if (j%50==0 && j!=0) fprint(fp,"\n        ");
	    sprintf(ps, "%c", raw.data.f2.allele[i][j]);
	    fpr(fp);
	    }
	fnl(fp);
    }
}


void 
new_read_f2_raw (
    FILE *fp, /* sends BADDATA, or MALLOC, INTERRUPT etc. if an error */
    char *symbols
)
{
    int i;

    run {
	allocate_f2_data(raw.num_markers,raw.data.f2.num_indivs);
	for (i=0; i<raw.num_markers; i++) read_raw_f2_locus(fp,i,symbols);
	for (i=0; i<num_traits; i++) read_raw_trait(fp,i);
    } on_error {
	free_f2_data(raw.num_markers/*,raw.data.f2.num_indivs*/);
	if (msg==ENDOFILE) baddata("unexpected end of file");
	relay;
    }
}



void 
read_raw_f2_locus (FILE *fp, int locus_num, char *symbol)
{
    int i, j;
    char c, *name, converted;

    name= get_temp_string();

    getdataln(fp);
    nstoken(&ln,sREQUIRED,name,NAME_LEN+1);
    while (name[0]=='#') {
        /* if there's a comment after the locus data */
        getdataln(fp); nstoken(&ln,sREQUIRED,name,NAME_LEN+1);
    }

    /* make sure name starts with '*' */
    if (name[0]!='*') {
        sprintf(badstr, "name does not have '*' or too many indivs in locus %s\n",
                raw.locus_name[locus_num > 0 ? locus_num-1 :0]);
	baddata(badstr);
    }
    if (len(name)<2) {
        sprintf(badstr, "expected locus name after '*', following locus %s\n",
                raw.locus_name[locus_num > 0 ? locus_num-1 : 0]);
	baddata(badstr);
    }
    if (!strin(NAME_FIRST_CHARS,name[1])) {
        sprintf(badstr, "illegal first character in locus name %s\n", name + 1);
	baddata(badstr);
    }
    for (i=2; name[i]!='\0'; i++) {
	if (name[i]=='-') name[i]='_';
        else if (!strin(NAME_CHARS,name[i])) {
	    sprintf(badstr, "illegal character '%c' in locus name %s",
                name[i],name+1);
	    baddata(badstr);
	}
    }

    /* make sure marker name is not a duplicate */
    for (j=0; j<locus_num; j++) {
	if (xstreq(name+1,raw.locus_name[locus_num])) {
	    sprintf(badstr, "two loci with the same name: %s", name + 1);
	    baddata(badstr);
	}
    }
    strcpy(raw.locus_name[locus_num],name+1);
    if (uppercase_genotypes) uppercase(ln);

    raw.data.f2.allelic_distribution[locus_num][0]= 0.0;
    raw.data.f2.allelic_distribution[locus_num][1]= 0.0;
    raw.data.f2.allelic_distribution[locus_num][2]= 0.0;
    raw.data.f2.allelic_distribution[locus_num][3]= 0.0;

    i=0;
    while (i<raw.data.f2.num_indivs) {
        while (nullstr(ln)) {
	    getdataln(fp);
	    if (uppercase_genotypes) uppercase(ln);
	}
	if (!parse_char(&ln,symbol,SKIPWHITE,&c)) {
	    if (c=='#') { strcpy(ln,""); continue; }
	    else if (c=='*') sprintf(badstr, "not enough data points for locus %s",
                                 raw.locus_name[locus_num]);
	    else sprintf(badstr, "illegal genotype character '%c' in locus %s",
                     c, raw.locus_name[locus_num]);
	    baddata(badstr);
	} else {
	    converted= symbol_value(c,symbol);
	    raw.data.f2.allele[locus_num][i++]= converted;
	    add_to_seg_dist(converted,locus_num);
	}
    }
    if (!nullstr(ln) && !parse_char(&ln,"#",SKIPWHITE,&c)) {
	sprintf(badstr, "too many data points for locus %s",
            raw.locus_name[locus_num]); baddata(badstr);
    }
    scale_seg_dist(locus_num);
}

#define ACTIVATE_EQUATION

void 
read_raw_trait (FILE *fp, int trait_num)
{
    int i, j, k, l, did_we_do_eqn, templen;
    char *tempstr, *temp2, *name, *tok;
    real neat, coc;
    EQUATION **postfixed;

    tempstr = get_temp_string();
    temp2 = get_temp_string();
    name = get_temp_string();
    tok = get_temp_string();

    while(nullstr(ln)) getdataln(fp);
    nstoken(&ln,sREQUIRED,name,NAME_LEN+1);
    while(name[0] == '#') {
        /* if there's a comment after the locus data */
        getdataln(fp); nstoken(&ln,sREQUIRED,name,NAME_LEN+1);
    }
    strcpy(tempstr,name);
    /* make sure name starts with '*' */
    if (tempstr[0] != '*') {
        sprintf(badstr, "trait does not have '*' or too many indivs in trait %s.\n",
                trait_name[trait_num > 0 ? trait_num-1 : 0]);
	baddata(name);
    }
    if ((templen=len(tempstr)) < 2) {
        sprintf(badstr, "expected trait name after '*', following trait %s\n",
                trait_name[trait_num > 0 ? trait_num-1 : 0]);
	baddata(badstr);
    }
    if (!strin(NAME_FIRST_CHARS,tempstr[1])) {
        sprintf(badstr, "illegal first character in trait name %s\n", tempstr);
	baddata(badstr);
    }

    for (i=2; i<templen; i++) {
        if (!strin(NAME_CHARS,tempstr[i])) {
	    sprintf(badstr, "illegal character '%c' in trait name %s",
                tempstr[i], tempstr);
	    baddata(badstr);
	}
    }
    tempstr=name+1;

    /* make sure trait name is not a duplicate */
    for (j=0; j<trait_num; j++) {
        strcpy(temp2,trait_name[j]);
	crunch(tempstr); crunch(temp2);
	if (!strcmp(tempstr,temp2)) {
	    sprintf(badstr, "two traits with the same name: %s", tempstr);
	    baddata(badstr);
	}
    }
	    
    strcpy(trait_name[trait_num], name+1);

    self_delimiting = ptr_to("=");
    did_we_do_eqn = 0;
    for(j = 0; j < raw.data.f2.num_indivs; j++) {
        if (nullstr(ln)) getdataln(fp);
	if (!rtoken(&ln,rREQUIRED,&neat)) {
	    if (stoken(&ln,sREQUIRED,tok) && strcmp(tok,"-")) {
	        
	        if (tok[0] == '*') {
		    sprintf(badstr, "not enough indivs in trait %s\n", name + 1);
		    baddata(badstr);
		}
		if (tok[0] == '#') {
		    /* rest of line is a comment, ignore it - decrement j, 
		       will be incremented by the for loop in the next step 
		       but this individual has not been filled yet */
		    strcpy(ln,"");
		    j--; 
		    continue;
		}
#ifdef ACTIVATE_EQUATION
		if (strcmp(tok,"=")==0) {
		    table_size = trait_num;
		    for(k=0;k<trait_num;k++) {
		        strcpy(variable_table[k],trait_name[k]);
		    }
		    /* get rid of comment after equation before sending it
		       to the equation parser (which is easily confused) */
		    l = 0;
		    while(ln[l] != '\0') {
		      if (ln[l]=='#') { ln[l]='\0'; break; }
		      l++;
		    }
		    postfixed = make_equation(ln,variable_lookup);
		    for(l=0;l<raw.data.f2.num_indivs;l++) {
		        if (did_we_do_eqn == 1) break;
			for(k=0;k<trait_num;k++) value_table[k] = trait[k][l];
			k=trait_num;
			coc=evaluate_equation(postfixed,value_lookup);
			trait[trait_num][l] =coc;
		    }
		    did_we_do_eqn = 1;
		}
#endif
		if (did_we_do_eqn == 1) { strcpy(ln,""); break; }
		else {
		  sprintf(badstr, "error in reading trait information for trait %s",
                  name);
		  baddata(badstr);
		}
	    }
	    else trait[trait_num][j] = MISSING_PHENO;
	}
	else
	    trait[trait_num][j] = neat;
    }
    self_delimiting = ptr_to(SELF_DELIMITING_TOKENS);
}


int 
symbol_value (int chr, char *symb)
{
    /* CHANGED FOR THIS VERION - NOW READS "-AHBCD" */

    if      (chr==symb[0]) return(MISSING_DATA);
    else if (chr==symb[1]) return(PARENTAL_TYPE_A);
    else if (chr==symb[2]) return(HYBRID_TYPE_H);
    else if (chr==symb[3]) return(PARENTAL_TYPE_B);
    else if (chr==symb[4]) return(TYPE_NOT_A);
    else if (chr==symb[5]) return(TYPE_NOT_B);
    
    else return(-10);
}


void 
write_traits (FILE *fp)
{
    int i, j;

    fprintf(fp,"%d\n",num_traits);

    for(i = 0; i < num_traits; i++) {
        sprintf(ps, "*%-10s", trait_name[i]); fpr(fp);
	for(j = 0; j < raw.data.f2.num_indivs; j++) {
	    if(j % 5 == 0 && j != 0) fprintf(fp,"\n          ");
	    sprintf(ps, "%12.3lf ", trait[i][j]); fpr(fp);
	}
	fnl(fp);
    }
    fnl(fp);
    fprint(fp,"#QTL only variables:\n");
    fprint(fp,"*Print mapm loci: 1\n");
    fprint(fp,"*Like tolerance: 0.001\n");
    fprint(fp,"*Brute force: 1\n");
    fprint(fp,"*Max intervals: -7\n");
    fprint(fp,"*Max continuous vars: -3\n");
    fprint(fp,"*Max wiggles: 0\n");
    fprint(fp,"*Max compares: 0\n");
    fprint(fp,"*Default units: 1\n");
    fprint(fp,"*Chromosomes: 0\n");
    fprint(fp,"*Number of contexts: 1\n");
    fprint(fp,"*Active context: 0\n");
    fnl(fp);
    fprint(fp,"*Context 1\n");
    if(num_traits == 1)
      fprint(fp,"*Trait: 0\n");
    else 
      fprint(fp,"*Trait: -1\n");
    fprint(fp,"*Named sequences: 0\n");
    fprint(fp,"*Sequence history: 0\n");
}




void add_to_seg_dist(char c, int locus)
{
    switch(c) {
	case PARENTAL_TYPE_A: {
	    raw.data.f2.allelic_distribution[locus][0] += 1.0;
	    return;
	}
	case PARENTAL_TYPE_B: {
	    raw.data.f2.allelic_distribution[locus][3] += 1.0;
	    return;
	}
	case TYPE_NOT_A: {
	    if(raw.data.f2.cross_type == F3_SELF) {
		raw.data.f2.allelic_distribution[locus][1] += (1.0/5.0);
		raw.data.f2.allelic_distribution[locus][2] += (1.0/5.0);
		raw.data.f2.allelic_distribution[locus][3] += (3.0/5.0);
	    } else {
		raw.data.f2.allelic_distribution[locus][1] += (1.0/3.0);
		raw.data.f2.allelic_distribution[locus][2] += (1.0/3.0);
		raw.data.f2.allelic_distribution[locus][3] += (1.0/3.0);
	    }
	    return;
	}
	case TYPE_NOT_B: {
	    if(raw.data.f2.cross_type == F3_SELF) {
		raw.data.f2.allelic_distribution[locus][0] += (3.0/5.0);
		raw.data.f2.allelic_distribution[locus][1] += (1.0/5.0);
		raw.data.f2.allelic_distribution[locus][2] += (1.0/5.0);
	    } else {
		raw.data.f2.allelic_distribution[locus][0] += (1.0/3.0);
		raw.data.f2.allelic_distribution[locus][1] += (1.0/3.0);
		raw.data.f2.allelic_distribution[locus][2] += (1.0/3.0);
	    }
	    return;
	}
	case HYBRID_TYPE_H: {
	    raw.data.f2.allelic_distribution[locus][1] += 0.5;
	    raw.data.f2.allelic_distribution[locus][2] += 0.5;
	    return;
	}
    }
}

void 
scale_seg_dist (int locus)
{
    real total;
    int i;
    
    total = 0.0;
    for(i=0; i<4; i++) {
	total += raw.data.f2.allelic_distribution[locus][i];
    }
    for(i=0; i<4; i++) {
	raw.data.f2.allelic_distribution[locus][i] =
	  raw.data.f2.allelic_distribution[locus][i] / total;
    }
}
