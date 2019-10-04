/******************************************************************************

  ####   ######   ####   #    #  ######  #    #   ####   ######           ####
 #       #       #    #  #    #  #       ##   #  #    #  #               #    #
  ####   #####   #    #  #    #  #####   # #  #  #       #####           #
      #  #       #  # #  #    #  #       #  # #  #       #        ###    #
 #    #  #       #   #   #    #  #       #   ##  #    #  #        ###    #    #
  ####   ######   ### #   ####   ######  #    #   ####   ######   ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#define INC_LIB
#define INC_SHELL
#define INC_MISC
#define INC_TABLE
#include "mapm.h"

#define SEQ_EXPAND_MAX	    10
#define NO_TERM_CHAR	    '\0'
#define LEFT_PAREN 	    "("
#define LEFT_ANGLE	    "<"
#define LEFT_SET   	    "{"
#define LEFT_BRACKET 	    "["
#define LEFT_PAREN_char	    '('
#define LEFT_ANGLE_char	    '<'
#define LEFT_SET_char	    '{'
#define LEFT_BRACKET_char   '['
#define RIGHT_PAREN_char    ')'
#define RIGHT_ANGLE_char    '>'
#define RIGHT_SET_char      '}'
#define RIGHT_BRACKET_char  ']'
#define LEFT_ANYTHING       "[{(<"
#define RIGHT_ANYTHING      "]})>"
#define OUT_OF_SYNTAX_CHARS "=-]})>"
#define EQUALS_SIGN         "="
#define QUESTION_MARK       "?"
#define SURROGATE_DASH      "~"
#define SURROGATE_DASH_char '~'

/* BADSEQ message */
int BADSEQ_errpos;
char *BADSEQ_errmsg; 

/* MAPMAKER's current sequence, etc. */
SEQ_NODE *seq=NULL, **seq_nodes=NULL;
char *seq_string=NULL;
int used_nodes;

/* Stuff for use in other files, but allocated here (BE CAREFUL STEVE!) */
int *seq_locus=NULL; /* for is_a_special_seq() (eg compile_ or check_seq())*/
int *use_locus=NULL; /* and for crunch_locus_list() in main.c */ 
char **seq_tokens=NULL;               /* for insert and delete in sys_cmds.c */
char *new_seq=NULL, *the_seq=NULL;   /* for commands in sys_cmds.c */

/* DECLARATIONS FOR INTERNAL USE ONLY! */
void compile_sequence();   /* args: char *str; sends BADSEQ on error */
int  delete_comments();
void expand_seq_hist();
void expand_seq_names();
void next_token();
void seqcomp();
void add_back_pointers();
void make_reordered_list();
int matching();
SEQ_NODE *mk_blank_node();
SEQ_NODE *mk_locus();
SEQ_NODE *mk_locus_range();
SEQ_NODE *mk_megalocus();
SEQ_NODE *mk_invertable_set();
SEQ_NODE *mk_unordered_set();
SEQ_NODE *get_seq_node();
SEQ_NODE *alloc_three_pt_seq();
void parse_fixed_distance();
void swap_for_dash();
void unswap_for_dash();
bool parse_seq_chrom();
void seq_expand();
void add_chr();
bool is_a_named_sequence();
bool is_an_old_sequence();
bool is_a_special_sequence();
void print_special_sequence();
bool do_get_order();	/* internal only! */
#define FORWARDS  1	/* for direction parameter to above */
#define BACKWARDS 2
#define REORDERED 3

bool is_a_named_locus(char *str /* must be a single token, downcased */, int *n);  /* internal use only */
bool valid_new_name(char *str); /* now does not check valid_name() */

char *maybe_seq=NULL;
char *seq_temp=NULL;
char *err_temp=NULL;

char token[TOKLEN+1]; /* global for compile_sequence(), expand_seq_names() */
bool Oagain;     /* used in for_all_orders() macro in sequence.h */
int Ti, Tj, Tk;  /* used in for_all_3pt_seqs macro in sequence.h */
SEQ_NODE *Tseq;  /* ... as is this */
int Pi, Pj;      /* used in for_all_locus_pairs macro in sequence.h */
bool Tagain;     /* ... as is this */


void sequence_init()
{ 
    self_delimiting=mkstrcpy(SELF_DELIMITING_TOKENS);
    setmsg(BADSEQ,"attempt to parse bad sequence",sender,NULL);
    array(seq_string,MAX_SEQ_LEN+1,char); seq_string[0]='\0';
    Tseq=alloc_three_pt_seq();
}


void allocate_seq_stuff(n_loci)
int n_loci;
{
    array(seq_temp,  MAX_SEQ_LEN+99,char); seq_temp[0]='\0';
    array(err_temp,  MAXLINE,char);        err_temp[0]='\0';
    array(use_locus, n_loci,int);
    array(seq_locus, MAX_SEQ_LOCI,int); /* MAX_SEQ_LOCI > n_loci */
    array(maybe_seq, MAX_SEQ_LEN+99,char); maybe_seq[0]='\0';
    array(new_seq,   MAX_SEQ_LEN+99,char); new_seq[0]='\0'; /* +99 for chrom */
    array(the_seq,   MAX_SEQ_LEN+99,char); the_seq[0]='\0';
    parray(seq_nodes,MAX_SEQ_TOKENS,SEQ_NODE); /* overestimate, but checked */
    matrix(seq_tokens,MAX_SEQ_TOKENS,TOKLEN+1,char);
}


void free_seq_stuff(/*int n_loci*/)
{
    unarray(seq_temp,char);
    unarray(err_temp,char);    
    unparray(seq_nodes,MAX_SEQ_TOKENS,SEQ_NODE);

    unarray(seq_locus,int);
    unarray(use_locus,int);
    unarray(maybe_seq,char);
    unarray(new_seq,char);
    unarray(the_seq,char);
    unmatrix(seq_tokens,MAX_SEQ_TOKENS+1,char);
}


void badseq(errmsg) char *errmsg; { BADSEQ_errmsg=errmsg; send(BADSEQ); }

void print_badseq()
{
    print("error in sequence: "); 
    /* print("Error in sequence '"); print(new_seq); print("'\n");
       print("                   "); space(BADSEQ_errpos); print("^\n"); */
    print(BADSEQ_errmsg); nl();
}

void print_sequence()
{
    if (the_seq_history_num==0) { print("sequence not set\n"); return; }
    sprintf(ps, "sequence #%d= ", the_seq_history_num); pr();
    if (current_chrom!=NO_CHROM) 
      { sprintf(ps, "%s: ", chrom2str(current_chrom)); pr(); }
    if (!nullstr(seq_string)) print(seq_string); else print("none");
    nl();
}


/******************** TO MAKE SEQ_NODE TREES ********************/

SEQ_NODE *get_seq_node()
{
    if (used_nodes>=MAX_SEQ_TOKENS) badseq("sequence too long");
    return(seq_nodes[used_nodes++]);
}

SEQ_NODE *mk_blank_node() 
{
    SEQ_NODE *node=get_seq_node();
    node->type= NONE;
    node->next= NULL;
    node->prev= NULL;
    node->fixed_distance= NOT_FIXED;
    return(node);
}


SEQ_NODE *mk_locus(num)
int num;
{
    SEQ_NODE *node=get_seq_node();
    node->type= SEQLOCUS;
    node->next= NULL;
    node->node.locus.num= num;
    node->prev= NULL;
    node->fixed_distance= NOT_FIXED;
    return(node);
}


SEQ_NODE *mk_locus_range(first,last)
int first, last;
{
    SEQ_NODE *node=get_seq_node();

    node=mk_locus(first);
    if (first==last) {
	node->next= NULL;
	return(node);
    }
    else if (first<last)
      node->next= mk_locus_range(first+1,last);
    else if (first>last)
      node->next= mk_locus_range(first-1,last);
    
    node->prev= NULL;
    node->fixed_distance= NOT_FIXED;
    return(node);
}


SEQ_NODE *mk_megalocus(contents)
SEQ_NODE *contents;
{
    SEQ_NODE *node=get_seq_node();
    node->type= MEGALOCUS;
    node->next= NULL;
    node->node.megalocus.contents= contents;
    node->prev= NULL;
    node->fixed_distance= NOT_FIXED;
    return(node);
}


SEQ_NODE *mk_invertable_set(contents)
SEQ_NODE *contents;
{
    SEQ_NODE *node=get_seq_node();
    node->type= INVERTABLE_SET;
    node->next= NULL;
    node->node.invertable_set.contents= contents;
    node->prev= NULL;
    node->fixed_distance= NOT_FIXED;
    return(node);
}


SEQ_NODE *mk_unordered_set(contents)
SEQ_NODE *contents;
{
    SEQ_NODE *node=get_seq_node();
    node->type= UNORDERED_SET;
    node->next= NULL;
    node->node.unordered_set.contents= contents;
    node->prev= NULL;
    node->fixed_distance= NOT_FIXED;
    return(node);
}


SEQ_NODE *alloc_three_pt_seq()
{
    int i;
    SEQ_NODE **node=NULL;
    parray(node,4,SEQ_NODE);
    
    node[0]->type= UNORDERED_SET;
    node[0]->next= NULL;
    node[0]->prev= NULL;
    for (i=1; i<=3; i++) {
	node[i]->type= SEQLOCUS;
	node[i]->node.locus.num= 0;
	node[i]->prev= NULL;
	node[i]->fixed_distance= NOT_FIXED;
    }
    node[0]->node.unordered_set.contents= node[1];
    node[0]->fixed_distance= NOT_FIXED;
    node[1]->next=node[2]; node[2]->next=node[3]; node[3]->next=NULL;
    return(node[0]);
}


/* compile_sequence() makes the big SEQ_NODE trees, sets global seq */

void compile_sequence(str)
char *str; /* may side-effect str? a little */
{
  	SEQ_NODE *result;
	char *copy; /* copy of pointer, not the string itself */
	
	BADSEQ_errmsg= null_string; BADSEQ_errpos= 0;
	run {
	    seq=NULL; used_nodes=0; result=NULL;
	    swap_for_dash(str);
	    copy=str;
	    seqcomp(&copy,&result,NO_TERM_CHAR); /* BADSEQ sent if err */
	    if (result!=NULL) add_back_pointers(result,NULL);
	    unswap_for_dash(str);

	} except_when(BADSEQ) { 
	    seq=NULL; used_nodes=0;
	    unswap_for_dash(str);
	    BADSEQ_errpos= len(str) - len(copy);
	    if (BADSEQ_errpos>0) BADSEQ_errpos--;
	    relay;
	}

	seq=result;
	return;
}	    
	  

void next_token(p_str) /* side-effects global variable token - internal use */
/* IMPORTANT: Keep in mind that token is static: therefore recursive calls to 
   seqcomp() will bash the contents of token before they return! */
char **p_str;
{ 
    if (!nstoken(p_str,sREQUIRED,token,TOKLEN)) 
        badseq("unexpected end of sequence"); 
}


void seqcomp(p_str,first_node,term_char)
char **p_str; /* str may not be empty */
SEQ_NODE **first_node;
char term_char;
{
    char c, *not_a_locus_errmsg;
    int n, m;
    SEQ_NODE **p, *q;

    p=NULL; 

    if (nullstr(*p_str)) { 
        if (term_char==NO_TERM_CHAR) return;
	  else badseq("empty megalocus or set");

    } else do {

        if (p==NULL) p= first_node; 
	  else do { p= &((*p)->next); } while(*p != NULL);

	next_token(p_str); /* Remember, self-delimiting= "<>[]{}()" */

	if (len(token)==1 && strin(OUT_OF_SYNTAX_CHARS,token[0])) {
	    if (strin(RIGHT_ANYTHING,token[0])) {
		if (term_char!=NO_TERM_CHAR && token[0]==term_char) {
		    *p=NULL; return; 
		} else { 
		    sprintf(msgstr, "mismatched '%c'", token[0]);
		    badseq(msgstr); 
		}
	    } else { /* not a RIGHT_ANYTHING, must be something else wrong */
		sprintf(msgstr, "did not expect '%c'", token[0]);
		badseq(msgstr); 
	    }

	} else if (streq(token,LEFT_BRACKET)) {
	    seqcomp(p_str,&q,RIGHT_BRACKET_char);
	    *p= mk_megalocus(q);
	    parse_fixed_distance(p_str,(term_char!=NO_TERM_CHAR),*p);
	    continue;

	} else if (streq(token,LEFT_ANGLE)) {
	    seqcomp(p_str,&q,RIGHT_ANGLE_char);
	    *p= mk_invertable_set(q);
	    parse_fixed_distance(p_str,(term_char!=NO_TERM_CHAR),*p);
	    continue;

	} else if (streq(token,LEFT_SET)) {
	    seqcomp(p_str,&q,RIGHT_SET_char);
	    *p= mk_unordered_set(q);
	    parse_fixed_distance(p_str,(term_char!=NO_TERM_CHAR),*p);
	    continue;

	} else if (is_a_locus(token,&n,&not_a_locus_errmsg)) {
	    if (parse_char(p_str,SURROGATE_DASH,SKIPWHITE,&c)) { /* "a - b" */
		if (nullstr(*p_str)) 
		  badseq("expected locus number or name after '-'");
		next_token(p_str);
		if (!is_a_locus(token,&m,&not_a_locus_errmsg))
		  badseq(not_a_locus_errmsg);
		else *p=mk_locus_range(n,m);
	    } else *p=mk_locus(n); /* just a "locus" */
	    parse_fixed_distance(p_str,(term_char!=NO_TERM_CHAR),*p);
	    continue;

	} else if (streq(token,"none")) {
	    continue;

	} else badseq(not_a_locus_errmsg); /* assumed it was... */

    } while(!nullstr(*p_str));
    
    return;
}


void parse_fixed_distance(str,inside,p)
char **str;
bool inside;
SEQ_NODE *p;
{
    char c;

    /* would like to make these self delimiting, but... */
    if (parse_char(str,"=*?",SKIPWHITE,&c)) { /* "=?*" */
	if (inside) badseq("can't fix distances inside {..}, <..>, or [..]");
	if (c=='*')      p->fixed_distance=UNLINKED_CM;
	else if (c=='?') p->fixed_distance=UNLINK_ME;
	else /* c=='=' */ {
	    if (!rtoken(str,rREQUIRED,&p->fixed_distance) ||
		!input_dist(&p->fixed_distance))
	      badseq("bad or missing map distance after '='");
	}
	if (nullstr(*str)) badseq("can't end sequence with a distance");
    }
}


void add_back_pointers(this,prev)
SEQ_NODE *this, *prev;
{
  this->prev= prev;
  switch (this->type) {
    case MEGALOCUS:      
    	add_back_pointers(this->node.megalocus.contents,(SEQ_NODE*) NULL);
    	break;
    case INVERTABLE_SET: 
    	add_back_pointers(this->node.invertable_set.contents,(SEQ_NODE*) NULL);
    	break;
    case UNORDERED_SET:
        add_back_pointers(this->node.unordered_set.contents,(SEQ_NODE*) NULL);
	    break;
  }
  if (this->next!=NULL) add_back_pointers(this->next,this);
}


#ifdef OBSOLETE
void free_seq_nodes(p)
SEQ_NODE *p;
{
  if (p==NULL) return;
  switch (p->type) {
    case MEGALOCUS:      free_seq_nodes(p->node.megalocus.contents); break;
    case INVERTABLE_SET: free_seq_nodes(p->node.invertable_set.contents);break;
    case UNORDERED_SET:  free_seq_nodes(p->node.unordered_set.contents); break;
  }
}
#endif



/********************* THINGS ONE CAN TO TO SEQUENCES *********************/

/* get_order() gets the current order of loci into an array, and perm_seq()
   "increments" the sequence to the next order. These functions could have 
   been wrapped up into one routine (eg. one decent of the seq tree), however
   this is easier to deal with. */

bool get_order(p,locus,theta,num_loci,max_loci) /* TRUE or sends msg */
SEQ_NODE *p;
int *locus;
real **theta;
int *num_loci, max_loci; 
/* locus is an array of locus nums, i is a ptr to its index, and max_loci is 
   the maximum number of loci the array holds (CRASH if it's exceeded) */
{
    int i, val;
    
    if (p==NULL) send(CRASH);
    *num_loci=0;
    if (theta!=NULL) for (i=0; i<max_loci-1; i++)
      theta[i][MALE]= theta[i][FEMALE]= NOT_FIXED;
    val= do_get_order(p,locus,theta,num_loci,max_loci,FORWARDS);
    for (i=0; i<*num_loci; i++) force_haplo_sanity(&locus[i],FALSE);
    return(val);
}


void get_one_order(p,map)
SEQ_NODE *p;
MAP *map;
{ 
    reset_seq(p,TRUE);
    clean_map(map); /* resets map->num_loci, among other things */
    get_order(p,map->locus,map->rec_frac,&map->num_loci,map->max_loci);
}


void get_list_of_all_loci(p,loci,num_loci,max_loci)
SEQ_NODE *p;
int *loci, *num_loci; /* side-effect *num and set loci[i] for 0<=i<*num */
int max_loci;         /* the max the array holds */
/* BUGS: get CRASH if array is exceeded, needs to be big enough to hold all 
   loci, not just haplo groups, even if use_haplotypes==TRUE. The size 
   returned by mapm_ready should be right for this. */
{ 
    *num_loci=0;
    reset_seq(p,TRUE); 
    get_order(p,loci,NULL,num_loci,max_loci);
}


bool alloc_list_of_all_loci(p,loci,num_loci)
SEQ_NODE *p;
int **loci, *num_loci; /* side-effect *loci and *num_loci setting (*loci)[i] */
{ 
    int total, *retoin=NULL; 

    *loci=NULL; 
    total=count_loci(p); 
    if (total==0) return(FALSE);
    array(retoin,total,int);

    *num_loci=0;
    reset_seq(p,TRUE);
    get_order(p,retoin,NULL,num_loci,total); 

    *loci=retoin; 
    return(TRUE);
}


bool do_get_order(p,locus,theta,i,max_loci,direction) /* TRUE or sends msg */
SEQ_NODE *p;
int *locus;
real **theta;
int *i, max_loci; 
int direction;
/* locus is an array of locus nums, i is a ptr to its index, and max_loci is 
   the maximum number of loci the array holds (CRASH if it's exceeded) */
{
    int j; 

    if (p==NULL) return(TRUE); /* should never happen? */
    if (direction==BACKWARDS) while(p->next!=NULL) p=p->next; /* to end */
    
    do {                   /* recursively loop over all nodes */
	switch (p->type) {
	  case SEQLOCUS:		
	    if (*i>=max_loci) send(CRASH);
	    locus[(*i)++]=p->node.locus.num; /* add this locus to the order */
	    break;
	  case LOCUS_RANGE:	
	    if ((*i+ p->node.locus_range.last -p->node.locus_range.first)
		>= max_loci) send(CRASH);
	    for (j=p->node.locus_range.first; j<=p->node.locus_range.last;j++) 
	      locus[(*i)++]= j; /* add all loci in the range */
	    break;
	  case MEGALOCUS: 	
	    do_get_order(p->node.megalocus.contents,locus,theta,i,max_loci,
		      FORWARDS); 
	    break;
	  case INVERTABLE_SET:	
	    /* add contents, maybe backwards*/
	    if ((!p->node.invertable_set.flip && direction!=BACKWARDS)  ||
		(p->node.invertable_set.flip  && direction==BACKWARDS)) 
	      do_get_order(p->node.invertable_set.contents,locus,theta,i,
			   max_loci,FORWARDS);
	    else 
	      do_get_order(p->node.invertable_set.contents,locus,theta,i,
			   max_loci,BACKWARDS);
	    break;
	  case UNORDERED_SET: 
	    make_reordered_list(p->node.unordered_set.contents,
				&p->node.unordered_set.reordered_contents);
	    do_get_order(p->node.unordered_set.reordered_contents,
		      locus,theta,i,max_loci,REORDERED);
	    break;
	  default: send(CRASH);
	}
	
	/* if not the last locus and if this iter added a locus to the list */
	if (theta!=NULL && p!=NULL && *i<max_loci) {
	    if (p->fixed_distance==UNLINK_ME && p->unlinked_now) {
		theta[(*i)-1][MALE]=   UNLINKED_CM;
		theta[(*i)-1][FEMALE]= UNLINKED_CM;
	    } else if (p->fixed_distance==UNLINK_ME) {
		theta[(*i)-1][MALE]=   NOT_FIXED;
		theta[(*i)-1][FEMALE]= NOT_FIXED;
	    } else if (p->fixed_distance!=NOT_FIXED) {
		theta[(*i)-1][MALE]=   p->fixed_distance;
		theta[(*i)-1][FEMALE]= p->fixed_distance;
	    } else {
		theta[(*i)-1][MALE]=   NOT_FIXED;
		theta[(*i)-1][FEMALE]= NOT_FIXED;
	    }
	}
	switch (direction) {
	    case FORWARDS:  p=p->next; break;
	    case BACKWARDS: p=p->prev; break;
	    case REORDERED: p=p->next_item; break;
	    default: send(CRASH);
	}
    } while (p!=NULL);
    return(TRUE); 
}


void make_reordered_list(p,list) /* intenal use only */
SEQ_NODE *p;
SEQ_NODE **list;
{
    int j;
    SEQ_NODE *ins_here, *q;
    
    /* We use the items list to hold the SEQ_NODE pointers of the contents
       of the set in the prescribed order. The insert_pos variables 
       dictate where each items get inserted: they must be set by
       reset_seq() and perm_seq(). */
    
    if (p->next!=NULL) make_reordered_list(p->next,list);
      else { *list=p; (*list)->next_item= (SEQ_NODE*) NULL; return; }
    
    if (p->insert_pos==0) { q= *list; *list=p; (*list)->next_item=q; }
    else {
	for (j=0, ins_here= *list; j<p->insert_pos-1; j++) 
	  ins_here= ins_here->next_item;
	q=ins_here->next_item; 
	ins_here->next_item=p; 
	ins_here->next_item->next_item=q;
    }
    return;
}


void reset_seq(p,reset_insert_pos)
SEQ_NODE *p;
bool reset_insert_pos;
{
  if (p==NULL) return;
  if (reset_insert_pos) p->insert_pos= 0;
  switch (p->type) {
    case MEGALOCUS: 	 reset_seq(p->node.megalocus.contents,TRUE); 
                         break;
    case INVERTABLE_SET: reset_seq(p->node.invertable_set.contents,TRUE); 
                         p->node.invertable_set.flip= FALSE; 
                         break;
    case UNORDERED_SET:	 reset_seq(p->node.unordered_set.contents,TRUE); 
                         break;
  }
  if (p->fixed_distance==UNLINK_ME) p->unlinked_now=FALSE;
  reset_seq(p->next,reset_insert_pos);
}


bool perm_seq(p,order_matters,in_unordered_set) 
/* return TRUE if there is another order, side-effecting the SEQ_NODE tree */
SEQ_NODE *p;
bool order_matters;
bool in_unordered_set;
{
  int max_position;
  SEQ_NODE *q;

  if (p==NULL) send(CRASH);
  if (p->next!=NULL) order_matters= TRUE; /* anytime >= 2 neighboring nodes */
  
  if (p->fixed_distance==UNLINK_ME) {
      if (p->unlinked_now==FALSE) { p->unlinked_now=TRUE; return(TRUE); }
      if (p->next!=NULL && perm_seq(p->next,order_matters,in_unordered_set))
	{ p->unlinked_now=FALSE; return(TRUE); }
  } else if (p->next!=NULL && perm_seq(p->next,order_matters,in_unordered_set))
      return(TRUE);

  switch (p->type) {
      case SEQLOCUS:	 return(FALSE);
      case LOCUS_RANGE:	 return(FALSE);
      
      case MEGALOCUS:   
      if (perm_seq(p->node.megalocus.contents,order_matters,FALSE)) {
	  reset_seq(p->next,!in_unordered_set);
	  p->unlinked_now=FALSE;
	  return(TRUE); 
      } else return(FALSE);

      case INVERTABLE_SET: 
      if (perm_seq(p->node.invertable_set.contents,order_matters,FALSE)) {
	  reset_seq(p->next,!in_unordered_set); 
	  p->unlinked_now=FALSE;
	  return(TRUE); 
      } else if (order_matters && !p->node.invertable_set.flip) {
	  p->node.invertable_set.flip= TRUE; 
	  reset_seq(p->node.invertable_set.contents,FALSE); 
	  reset_seq(p->next,!in_unordered_set); 
	  p->unlinked_now=FALSE;
	  return(TRUE);
      } else return(FALSE);
      
    case UNORDERED_SET:
      if (perm_seq(p->node.unordered_set.contents,order_matters,TRUE)) { 
	  reset_seq(p->next,!in_unordered_set); 
	  p->unlinked_now=FALSE;
	  return(TRUE);
      } else {
	  q=p->node.unordered_set.contents; 
	  max_position=1; 
	  while (q->next!=NULL) q=q->next;  /* go to the end */
	  while ((q=q->prev)!=NULL)         /* start looping with 2nd to end */
	      if (q->insert_pos < max_position && 
		  (max_position!=1 || order_matters)) {
		    ++(q->insert_pos);
		    while ((q=q->next)!=NULL) q->insert_pos= 0;
		    reset_seq(p->node.unordered_set.contents,FALSE);
		    reset_seq(p->next, !in_unordered_set);
		    p->unlinked_now=FALSE;
		    return(TRUE);
	      } else ++max_position;
	  return(FALSE);
      }
    }
    return(FALSE);
}


int count_loci(p)
SEQ_NODE *p;
{ 
if (p==NULL) return(0);
else switch (p->type) {
  case SEQLOCUS: 		
    return(count_loci(p->next) + 1);
  case LOCUS_RANGE: 	
    return(count_loci(p->next) + p->node.locus_range.last 
	   - p->node.locus_range.first + 1);
  case MEGALOCUS: 
    return(count_loci(p->next) + count_loci(p->node.megalocus.contents));
  case INVERTABLE_SET:  
    return(count_loci(p->next) + count_loci(p->node.invertable_set.contents));
  case UNORDERED_SET:
    return(count_loci(p->next) + count_loci(p->node.unordered_set.contents));
  default: send(CRASH);
  }
  send(CRASH); return(0); /* never reached */
}


bool has_fixed_dists(p)
SEQ_NODE *p;
{ 
if (p==NULL) return(FALSE);
else switch (p->type) {
  case SEQLOCUS: 		
  case LOCUS_RANGE:
    return(p->fixed_distance!=NOT_FIXED || has_fixed_dists(p->next));
  case MEGALOCUS:
    return(p->fixed_distance!=NOT_FIXED || has_fixed_dists(p->next) ||
	   has_fixed_dists(p->node.megalocus.contents));
  case INVERTABLE_SET:  
    return(p->fixed_distance!=NOT_FIXED || has_fixed_dists(p->next) ||
	   has_fixed_dists(p->node.invertable_set.contents));
  case UNORDERED_SET:
    return(p->fixed_distance!=NOT_FIXED || has_fixed_dists(p->next) ||
	   has_fixed_dists(p->node.unordered_set.contents));
  default: send(CRASH);
  }
  send(CRASH); return(0); /* never reached */
}


bool unpermutable(p)
SEQ_NODE *p;
{ reset_seq(p,TRUE); return(!perm_seq(p,FALSE,FALSE)); }



/***************** rather crufty stuff for 3pt *****************/

SEQ_NODE *Tinit(ok,i_cnt,j_cnt,k_cnt,locus,num_loci)
bool *ok;
int *i_cnt, *j_cnt, *k_cnt;
int *locus, num_loci;
{
    if (num_loci<3) send(CRASH);
    *ok= TRUE;
    *i_cnt= 0; 
    *j_cnt= 1;
    *k_cnt= 2;
    Tseq->node.unordered_set.contents->node.locus.num= locus[0];
    Tseq->node.unordered_set.contents->next->node.locus.num= locus[1];
    Tseq->node.unordered_set.contents->next->next->node.locus.num= locus[2];
    add_back_pointers(Tseq,(SEQ_NODE*)NULL);
    return(Tseq);
}
    

bool Tnext(i_cnt,j_cnt,k_cnt,locus,num_loci,p)
int *i_cnt, *j_cnt, *k_cnt;
int *locus, num_loci;
SEQ_NODE *p;
{
    if (++*k_cnt > num_loci-1) {		/* k is innermost loop */
	if (++*j_cnt > num_loci-2) {
	    if (++*i_cnt > num_loci-3) {
		return(FALSE);
	    } 
	    *j_cnt= *i_cnt+1;
	}
	*k_cnt= *j_cnt+1;
    }
    p->node.unordered_set.contents->node.locus.num= locus[*i_cnt];
    p->node.unordered_set.contents->next->node.locus.num= locus[*j_cnt];
    p->node.unordered_set.contents->next->next->node.locus.num= locus[*k_cnt];
    return(TRUE);
}
    

/***************** MAPMAKER's Current Sequence *****************/

void set_current_seq(str,expanded) /* use nullstr to unset current sequence */
char *str; /* side-effected: assumed to be MAX_SEQ_LEN long (use new_seq) */
bool expanded; /* if seq_str is set to the name-expanded seq or not */
/* if error relay BADSEQ message without changing anything */
{
    int save_chrom, num_loci;
    char *rest, chrom_name[99];

    if (nullstr(str)) { 
	seq_string[0]='\0'; seq=NULL; 
	current_chrom=NO_CHROM;
	return; 
    }
    run {
	save_chrom=current_chrom;
	if (!parse_seq_chrom(str,&current_chrom,&rest)) 
	  current_chrom=save_chrom;
	expand_seq_hist(str);
	strcpy(maybe_seq,rest);
	expand_seq_names(maybe_seq); /* needs current_chrom set right */
	compile_sequence(maybe_seq); /* does the right thing w/"none" or "" */
	/* compile_sequence() sends BADSEQ message or sets global seq */
	num_loci=count_loci(seq);
	if (num_loci>MAX_SEQ_LOCI) badseq("too many loci in sequence");
	/* if we fall to here, we got it */
	if (expanded) strcpy(seq_string,maybe_seq);
	  else strcpy(seq_string,rest);
	despace(seq_string);

	if (current_chrom!=NO_CHROM) {
	    sprintf(chrom_name, "%s: ", chrom2str(current_chrom));
	    strins(maybe_seq,chrom_name);
	}
	add_to_seq_history(maybe_seq,TRUE);

	/* Any changes in here should be reflected also in check_current_seq()
	   and name_sequence() */
    } when_aborting { current_chrom=save_chrom; relay; }
}


void check_current_seq(num_loci)
int *num_loci;
{
    int foo;
    char *rest, chrom_name[99];

    if (nullstr(seq_string)) { seq=NULL; return; }
    strcpy(maybe_seq,seq_string);
    parse_seq_chrom(maybe_seq,&foo,&rest); /* current_chrom shouldn't change */
    expand_seq_names(rest); /* needs current_chrom set right */
    compile_sequence(rest); /* sends BADSEQ message or sets global seq */
    *num_loci=count_loci(seq);
    if (*num_loci>MAX_SEQ_LOCI) badseq("too many loci in sequence");

    if (current_chrom!=NO_CHROM) {
	sprintf(chrom_name, "%s: ", chrom2str(current_chrom));
	strins(rest,chrom_name);
    }
    add_to_seq_history(rest,FALSE);
}


#define NO_ORDERED_ARGS \
  "can't use ordered/special sequences ([..], {..}, <..>, etc) as argument"

void parse_locus_args(loci,num_loci) /* send error if fail */
int **loci, *num_loci; /* both side-effected */
{
    char *copy, *errmsg, *rest, name[TOKLEN+1];
    int locus=0, *result, n, next, i, save_chrom;
    bool was_range;

    if (nullstr(args)) { num_loci=0; return; }
    nstrcpy(maybe_seq,args,MAX_SEQ_LEN);

    run {
	save_chrom= current_chrom;
	parse_seq_chrom(maybe_seq,&current_chrom,&rest);
	expand_seq_hist(rest);
	expand_seq_names(rest); /* needs current_chrom set right */
	swap_for_dash(rest);
	current_chrom= save_chrom;
	if (nullstr(rest)) { *num_loci=0; return; }
    } except_when(BADSEQ) {
	current_chrom= save_chrom;
	sprintf(ps, "bad locus (sequence) arguments\n%s", BADSEQ_errmsg);
	error(ps);
    }

    copy=rest; n=0; was_range=FALSE; /* count and check pass */
    while (stoken(&copy,sREQUIRED,name)) {
	/* remember the self delimiting chars are in effect */
	if (strin(self_delimiting,name[0]) && !streq(name,SURROGATE_DASH))
	  error(NO_ORDERED_ARGS);
	else if (streq(name,SURROGATE_DASH)) {
	    if (n==0 || !stoken(&copy,sREQUIRED,name) || was_range)
	      error ("illegal range:\ncorrect usage is locus1-locus2");
	    if (!is_a_locus(name,&next,&errmsg) || !nullstr(errmsg)) 
	      error(errmsg);
	    was_range=TRUE; n+=abs(next-locus); /* +1-1 */
	} else if (!is_a_locus(name,&locus,&errmsg) || !nullstr(errmsg))
	  error(errmsg);
	else { n++; was_range=FALSE; }
    }
    *num_loci=n;

    if (loci==NULL) return;
    array(result,n,int);
    copy=rest; n=0; /* get'em */
    while (stoken(&copy,sREQUIRED,name))
      if (streq(name,SURROGATE_DASH)) {
	  stoken(&copy,sREQUIRED,name); is_a_locus(name,&next,&errmsg);
	  if (next>=locus) for (i=locus+1; i<=next; i++) result[n++]=i;
	    else for (i=locus-1; i>=next; i--) result[n++]=i;
      } else {
	  is_a_locus(name,&locus,&errmsg);
	  result[n++]=locus;
      }
    if (n!=*num_loci) send(CRASH);
    *loci=result;
}


#define CHROM_NOT_EXISTS "there is no chromosome named '%s'"
#define SHORTHAND \
  "can't list other loci with chromosome-name, 'any', 'all', or 'none'"

bool parse_seq_chrom(str,chrom,rest) /* sets current_chrom */
char *str; /* bashed */
int *chrom;
char **rest;
{
    char *copy, name[TOKLEN+1];

    copy=str;
    if (split_string(str,rest,':')) {
	if (nullstr(str)) {
	    *chrom=NO_CHROM;
	    if (nullstr(*rest)) strcpy(*rest,seq_string);
	} else if (!stoken(&copy,sREQUIRED,name)) {
	    badseq("bad chromosome name name before ':'");
	} else if (streq(name,"any")) {
	    *chrom=NO_CHROM;
	    if (nullstr(*rest)) strcpy(*rest,seq_string);
	} else if (isa_chrom(name,chrom)) {
	    if (nullstr(*rest)) strcpy(*rest,seq_string);
	} else { 
	    sprintf(ps, CHROM_NOT_EXISTS, name); badseq(ps);
	}
	if (strin(*rest,':')) badseq("':' is only used to select chromosome");

    } else if (stoken(&copy,sREQUIRED,name)) { /* can't fail, can it? */
	/* maybe shorthand c-name w/no colon - these we substitute for */
	*rest=str;
	if (streq(name,"all") || streq(name,"any")) {
	    if (!nullstr(copy)) badseq(SHORTHAND);
	    *chrom=NO_CHROM;
	    strcpy(str,"all");
	} else if (streq(name,"none")) {
	    if (!nullstr(copy)) badseq(SHORTHAND);
	    *chrom=NO_CHROM;
	    strcpy(str,"none");
	} else if (isa_chrom(name,chrom)) {
	    if (!nullstr(copy)) badseq(SHORTHAND);
	    strcpy(str,"assigned");
	} else {
	    return(FALSE);
	}

    } else /* nullstr */ {
	*rest=str;
	return(FALSE);
    }
    return(TRUE);
}
	    

void make_compare_seq(locus,num_loci,start_set,set_size)
int *locus, num_loci, start_set, set_size;
{
    int i;
    char *last;
    
    last=seq_temp; last[0]='\0';
    for (i=0; i<num_loci; i++) {
	if (i==start_set) { strcpy(last,"{"); last++; }
	sprintf(last,"%d ",locus[i]+1); 
	while (*last!='\0') ++last;
	if (i==start_set+set_size-1) { strcpy(last,"}"); last++; }
    }
    compile_sequence(seq_temp);
}



/******************* SEQUENCE TEXT HACKING ********************/

#define mismatch_(pos,ch) \
 { BADSEQ_errpos=pos; sprintf(err_temp,"mismatched '%c'",ch); badseq(err_temp); }

int delete_comments(match_chr,str)
/* Delete all parenthesized comments, dealing with levels of delimiters 
   correctly. Side-effect both str and the pointer to it (to indicate where 
   an error occured). The origonal pointer to str should be saved by the 
   caller, as it will point to the edited sequence: *p_str may not. 
   Also now checks for all balanced punctuation. Signal BADSEQ if fails. */
char match_chr, *str;
{
    int i,j; 
    
    for (i=0; str[i]!='\0'; i++)
      if (match_chr!='\0' && str[i]==match_chr) return(i);
      else if (str[i]==LEFT_BRACKET_char) {
	  if ((j=delete_comments(RIGHT_BRACKET_char,str+i+1))<0)
	    mismatch_(i,LEFT_BRACKET_char) else i=j;
      } else if (str[i]==LEFT_ANGLE_char) {
	  if ((j=delete_comments(RIGHT_ANGLE_char,str+i+1))<0)
	    mismatch_(i,LEFT_ANGLE_char) else i=j;
      } else if (str[i]==LEFT_SET_char) {
	  if ((j=delete_comments(RIGHT_SET_char,str+i+1))<0)
	    mismatch_(i,LEFT_SET_char) else i=j;
      } else if (str[i]==LEFT_PAREN_char) {
	  if ((j=matching(RIGHT_PAREN_char,str+i+1))<0) 
	    mismatch_(i,LEFT_PAREN_char) 
	  else {
	      strdel(str+i,j+2);
	      despace(str);
	      continue; /* with i+1 */
	  }
      }
    return(-1);
}


int matching(match_chr,str) /* internal use only */
/* return index (or -1 if none) of char matching the match_chr */
char match_chr, *str;
{
    int i;

    for (i=0; str[i]!=match_chr; i++)
        if (str[i]=='\0') return(-1);
	else if (str[i]==LEFT_PAREN_char) {
	    if ((i=matching(RIGHT_PAREN_char,str+i+1))==-1) return(-1);
	} else if (str[i]==LEFT_BRACKET_char) {
	    if ((i=matching(RIGHT_BRACKET_char,str+i+1))==-1) return(-1);
	} else if (str[i]==LEFT_ANGLE_char) {
	    if ((i=matching(RIGHT_ANGLE_char,str+i+1))==-1) return(-1);
	} else if (str[i]==LEFT_SET_char) {
	    if ((i=matching(RIGHT_SET_char,str+i+1))==-1) return(-1);
	}
    return(i);
}


/* We can't have the '-' symbol be SELF_DELIMITING, because then negative 
   numbers parse wrong. So, we swap '-' for SURROGATE_DASH which should be
   SELF_DELIMITING */

void swap_for_dash(str) /* internal use only */
char *str;
{ int i; 
  for (i=0; str[i]!='\0'; i++) if (str[i]=='-') str[i]=SURROGATE_DASH_char; }


void unswap_for_dash(str) /* internal use only */
char *str;
{ int i; 
  for (i=0; str[i]!='\0'; i++) if (str[i]==SURROGATE_DASH_char) str[i]='-'; }


#define BADSEQ_TOOLONG \
"expansion of names or history results in a sequence which is too long"
#define BADSEQ_LOOPING \
"named sequences are defined in terms of each other"
#define BADSEQ_BADCHROM \
"bad chromosome name with ':' in history or name expansion"
#define BADSEQ_NOTCHROM \
"chromosome name in history or name expansion is different than sequence's"

void expand_seq_names(str) /* names, and any hist refs they contain BUT... */
char *str; /* side-effected! */
/* send BADSEQ if expansion >MAX_SEQ_LEN chars or if it contains circularly 
   defined names- this should never happen */
{ 
    char *token_start, *expansion=NULL, *errmsg=NULL, *all, *rest;
    int loop_count, chrom;

    loop_count=0; all=str; BADSEQ_errpos= -1; 
    while (!nullstr(str)) {
	while (white(*str)) str++;    /* gobble whitespace */
	token_start=str; next_token(&str);
	if (is_a_sequence(token,&expansion,&errmsg)) {
	    if (!nullstr(errmsg)) badseq(errmsg);
	    if(seq_temp != expansion) strcpy(seq_temp,expansion); /* if-statement to avoid valgrind warning of overlapping strings when pointers seq_temp = expansion, but why are they equal?) */
	    if (parse_seq_chrom(seq_temp,&chrom,&rest)) {
		if (chrom!=current_chrom && current_chrom!=NO_CHROM)
		  badseq(BADSEQ_NOTCHROM);
		current_chrom=chrom;
	    }
	    strdel(token_start,len(token));              /* delete the token */
	    if (len(all)+len(rest)>=MAX_SEQ_LEN) badseq(BADSEQ_TOOLONG);
	    strins(token_start,rest);  /* replace it with its expansion */
	    str=token_start;           /* step back to expand the expansion */
	    if (loop_count++==SEQ_EXPAND_MAX) badseq(BADSEQ_LOOPING);
	}
    }
}


void expand_seq_hist(str) 
char *str; /* side-effected! */
/* send BADSEQ if expansion >MAX_SEQ_LEN chars or if it contains circularly 
   defined names- this should never happen */
{ 
    char *token_start, *rest, *value, *errmsg=NULL, *all;
    int loop_count, chrom;

    loop_count=0; all=str; BADSEQ_errpos= -1; 
    while (!nullstr(str)) {
	while (white(*str)) str++;    /* gobble whitespace */
	token_start=str; next_token(&str);
	if (is_an_old_sequence(token,&value,&errmsg)) {
	    if (!nullstr(errmsg)) badseq(errmsg);
	    strcpy(seq_temp,value);
	    if (parse_seq_chrom(seq_temp,&chrom,&rest)) {
		if (chrom!=current_chrom && current_chrom!=NO_CHROM)
		  badseq(BADSEQ_NOTCHROM);
		current_chrom=chrom;
	    }
	    strdel(token_start,len(token));              /* delete the token */
	    if (len(all)+len(rest)>=MAX_SEQ_LEN) badseq(BADSEQ_TOOLONG);
	    strins(token_start,rest); 
	    str=token_start;            /* step back to expand the expansion */
	    if (loop_count++==SEQ_EXPAND_MAX) badseq(BADSEQ_LOOPING);
	}
    }
}


void tokenize_seq(seq,token,num_tokens)
char *seq, **token;
int *num_tokens;
{
    swap_for_dash(seq);  /* makes '-' a delimiter ('~') */
    expand_seq_names(seq);

    *num_tokens= 0; /* remember self_delimiting! */
    while (stoken(&seq,sREQUIRED,token[++*num_tokens])>0) {
	if (*num_tokens>=MAX_SEQ_TOKENS) 
	  error("sequence too long (too many words)");
    }
}


void untokenize_seq(seq,token,num_tokens)
char *seq, **token;
int num_tokens;
{
    int i, n, length;
    char *s;

    s=seq; *s='\0'; n=0;
    for (i=0; i<num_tokens; i++) {
	length= len(token[i]);
	if (n+length+1>MAX_SEQ_LEN) error("sequence too long");
	strins(s,token[i]); n+=length; while (*s!='\0') s++;
	if (i!=num_tokens-1) { *s= ' '; s++; *s='\0'; }
    }
    despace(seq);
    unswap_for_dash(seq);
}


/************** MAPMAKER sequence names and history lookups **************/
/* As of V3, the abbreviation mechanism is disabled, as it has been causing 
all sorts of problems, and as the special_names will make these problems worse.

To keep repeditive tests from taking forever, we assume that single
tokens are parsed out before they are tested. We also assume that all
names inserted into the tables and so forth have been checked (that
is, using valid_name()) so that any token matching one must also be
valid. In addition we assume that stored names must start with an
alphabetic character (eg, they alone may not start with a symbol and
may not parse into an int). These assumptions are used implicitly in
the following two routines: it makes the code somewhat intertwined in
an ugly fashion, but we don't want to waste too much time here. 

Wierd error behavior of these two routines (and others below?) is that
they may return TRUE *EVEN IF THERE IS AN ERROR*, in which case
!nullstr(*why_not).  This is to say when a name *should* be ok but it
isn't (e.g. "sequence new" with no new markers) */

bool is_a_locus(str,n,why_not)
char *str;      /* must be a non-null token and a valid_name */
int *n;         /* side-effected with locus# iff TRUE is returned */
char **why_not; /* side-effected iff FALSE is returned */
{ 
    int num;

    if (!data_loaded() || !is_a_token(str)) send(CRASH);
    *why_not=ptr_to("");
    
    if (itoken(&str,iREQUIRED,&num)) {
	*n=num-1;
	if (!irange(n,0,raw.num_markers-1)) {
	    sprintf(err_temp, "locus number '%d' out of valid range", num);
	    *why_not=ptr_to(err_temp);
	    return(FALSE);
	} else return(TRUE);

    } else if (!valid_name(str)) {
	sprintf(err_temp, "invalid name '%s'", str);
	*why_not=ptr_to(err_temp);
	return(FALSE);

    } else if (!is_a_named_locus(str,n)) {
	sprintf(err_temp, "name '%s' is not defined", str);
	*why_not=ptr_to(err_temp);
	return(FALSE);
    }
    return(TRUE);
}


bool is_a_sequence(str,seq,why_not)
char *str; 	 /* must be a non-null token and a valid_name */
char **seq;      /* side-effected */
char **why_not;  /* side-effected iff FALSE is returned */
{
    char name[TOKLEN+1];

    if (!data_loaded() || !is_a_token(str)) send(CRASH);
    *why_not=ptr_to("");
    strcpy(name,str); lowercase(name);

    if (is_an_old_sequence(str,seq,why_not)) {
	return(TRUE);    /* but why_not may be set */
	
    } else if (!valid_name(str)) {
	*why_not= ptr_to("illegal name");
	return(FALSE);

    } else if (is_a_special_sequence(name,seq,why_not)) {
	return(TRUE);    /* why_not may be set */
	
    } else if (!is_a_named_sequence(str,seq)) {
	*why_not=ptr_to("name is not defined");
	return(FALSE);

    } else return(TRUE);
}


bool name_sequence(name,str,why_not,expanded)
char *name;
char *str; /* assume its MAX_SEQ_LEN long, for expansions */
char **why_not;
bool expanded;
{ 
    char *foo, *rest;
    int save_chrom;
    /* if (name[0]=='*') name++; */

    despace(str);
    if (expanded) {
	run {
	    save_chrom=current_chrom;
	    if (parse_seq_chrom(str,&current_chrom,&rest)) {
		expand_seq_hist(rest);
		expand_seq_names(rest); /* needs current_chrom set */
		if (current_chrom!=NO_CHROM) 
		  sprintf(maybe_seq, "%s: %s", chrom2str(current_chrom), rest);
		else strcpy(maybe_seq,rest);
	    } else {
		expand_seq_hist(str);
		expand_seq_names(str); /* needs current_chrom set */
		if (current_chrom!=NO_CHROM)
		  sprintf(maybe_seq, "%s: %s", chrom2str(current_chrom), str);
		else strcpy(maybe_seq,str);
	    }
	} on_exit { current_chrom=save_chrom; relay_messages; }
    } else strcpy(maybe_seq,str);

    if (!valid_name(name) || len(name)>NAME_LEN) {
	*why_not= ptr_to("illegal name"); 
	return(FALSE);
    } else if (!valid_new_name(name) && !is_a_named_sequence(name,&foo)) {
	*why_not= ptr_to("can not re-define special named sequences");
	return(FALSE);
    }
    put_named_entry(name,maybe_seq,context[active_context]->named_sequences);
    return(TRUE);
}
  

bool unname_sequence(name,why_not)
char *name;
char **why_not;
{ 
    int fail;
    char *seq;

    if (!valid_name(name)) {	
	*why_not= ptr_to("illegal name");
	return(FALSE);
    }
    /* if (name[0]=='*') name++; */
    if (is_a_named_sequence(name,&seq)) {
	*why_not= ptr_to("name is not defined");
	return(FALSE); 
    }
    if (!delete_named_entry(name,
			    context[active_context]->named_sequences,&fail)){
	if (fail==NAME_DOESNT_MATCH)
	    *why_not= ptr_to("name is not defined");
	else 
	    *why_not= ptr_to("name is ambiguous");
	return(FALSE); 
    }
    return(TRUE);
}
  

void add_to_seq_history(seq,is_next_entry)
char *seq;
bool is_next_entry;
{ 
    int num;
    char *the_sequence;

    /* note: the_seq_history_num is the next number, not this one! */
    if (is_next_entry) {
	num=the_seq_history_num;
	if (!nullstr(seq)) put_numbered_entry(seq,the_sequence_history,&num);
	  else put_numbered_entry("none",the_sequence_history,&num);
	num++; 
	the_seq_history_num=num;
    } else {	
	if (the_seq_history_num==0) return;
	num=the_seq_history_num-1;
	if (!get_numbered_entry(num,&the_sequence,the_sequence_history)) 
	  send(CRASH);
	if (!nullstr(seq)) strcpy(the_sequence,seq);
	  else strcpy(the_sequence,"none");
    }
}


bool valid_new_name(char *str) /* now does not check valid_name() */
{
    int foo; char *x, *y;
    return (!is_a_named_locus(str,&foo) && 
	    !is_a_special_sequence(str,&x,&y) &&
	    !isa_chrom(str,&foo) &&
	    !isa_class(str,&foo));
}


/********** Internal procedures to implement these... **********/

bool is_a_named_locus(char *str /* must be a single token, downcased */, int *n)  /* internal use only */
{
    int i;

    for (i=0; i<raw.num_markers; i++)
      if (xstreq(str,raw.locus_name[i])) { *n=i; return(TRUE); }
    return(FALSE);
}


bool is_a_named_sequence(char *str, char **seq)  /* internal use only */
{
    char *foo; int fail;
    return (get_named_entry(str,&foo,seq,the_context->named_sequences,&fail));
}


bool is_an_old_sequence(str,seq,why_not)  /* internal use only */
char *str;
char **seq;      /* side-effected iff TRUE returned */
char **why_not;  /* side-effected iff TRUE returned and *str=="" */
{
    int n;
    char *rest;

    *why_not=ptr_to("");
    if (str[0]=='#')  {
	if (the_seq_history_num==0) {
	    *why_not= ptr_to("no sequences - can't use '#'");
	    return(TRUE);
	}
	rest=str+1;
	if (!itoken(&rest,iREQUIRED,&n)) {
	    *why_not= ptr_to("invalid sequence number following '#'");
	    return(TRUE);
	}
	if (n<=0) n=the_seq_history_num-n-1;
	  else --n;
	if (!get_numbered_entry(n,seq,the_sequence_history)) {
	    *why_not= ptr_to("sequence number following '#' is out of range");
	    return(TRUE);
	} else 
	  return(TRUE);
    } else return(FALSE);
}


bool is_a_special_sequence(str,seq,why_not) /* internal use only? */
char *str; 
char **seq;      /* side-effected iff TRUE returned */
char **why_not;  /* side-effected iff TRUE returned and *str=="" */
{
    int chrom, one_chrom, i, j, k, n, contig, num, ordered_seq;
    char *last;

    one_chrom=FALSE; chrom=0;
    if (current_chrom!=NO_CHROM) { one_chrom=TRUE; chrom=current_chrom; }
    n=0; ordered_seq=FALSE;
    *seq= *why_not= ptr_to("");

    if (streq(str,"new")) {
	for (i=0; i<raw.num_markers; i++) {
	    if (one_chrom && !assigned_to(i,chrom)) continue;
	    if (!modified[i]) continue;
	    if (use_haplotypes && haplotype_subordinate(i)) continue;
	    seq_locus[n++]= i;
	}

    } else if (streq(str,"assigned") || streq(str,"ass")) {
	for (i=0; i<raw.num_markers; i++) {
	    if (!assigned(i)) continue;
	    if (one_chrom && assignment_chrom(i)!=chrom) continue;
	    if (use_haplotypes && haplotype_subordinate(i)) continue;
	    seq_locus[n++]= i;
	}

    } else if (streq(str,"anchors") || streq(str,"anchor")) {
	for (i=0; i<raw.num_markers; i++) {
	    if (!anchor_locus(i)) continue;
	    if (one_chrom && assignment_chrom(i)!=chrom) continue;
	    if (use_haplotypes && haplotype_subordinate(i)) continue;
	    seq_locus[n++]= i;
	}

    } else if (streq(str,"changed")) {
	for (i=0; i<raw.num_markers; i++) {
	    if (!assigned(i) || assignment_state(i)!=A_CHANGED) continue;
	    if (one_chrom && assignment_chrom(i)!=chrom) continue;
	    if (use_haplotypes && haplotype_subordinate(i)) continue;
	    seq_locus[n++]= i;
	}

    } else if (isa_class(str,&num)) {
	for (i=0; i<raw.num_markers; i++) {
	    if (one_chrom && !assigned_to(i,chrom)) continue;
	    if (class[i]!=num) continue;
	    if (use_haplotypes && haplotype_subordinate(i)) continue;
	    seq_locus[n++]= i;
	}

    } else if (streq(str,"framework") || streq(str,"frame")) {
	/* frameworks have forced haplo_sanity */
	for (j=0; j<chromosome->num_maps; j++) {
	    if (one_chrom && j!=chrom) continue;
	    for (k=0; k<chromosome->map_list[j]->num_loci; k++) {
		i=chromosome->map_list[j]->locus[k];
		seq_locus[n++]= i;
	    }
	}
	if (one_chrom) ordered_seq=TRUE;
	  else { ordered_seq=FALSE; sort_loci(seq_locus,n); }


    /***** not chrom specific names *****/

    } else if (streq(str,"all")) {
	for (i=0; i<raw.num_markers; i++) {
	    if (use_haplotypes && haplotype_subordinate(i)) continue;
	    seq_locus[n++]= i;
	}

    } else if (streq(str,"unassign") || streq(str,"unassigned")) {
	for (i=0; i<raw.num_markers; i++) {
	    if (assigned(i)) continue;
	    if (use_haplotypes && haplotype_subordinate(i)) continue;
	    seq_locus[n++]= i;
	}

    } else if (nstreq(str,"group",5)) { /* insane */
	if (num_groups==0) {
	    *why_not=ptr_to("no grouped markers in the data set");
	    return(TRUE);
	} else if (sscanf(str+5,"%d",&num)!=1 || num<1 || num>num_groups-1) {
	    *why_not=ptr_to("that group is not defined"); /* fix this */
	    return(TRUE);
	}
	num--;
	for (i=0; i<raw.num_markers; i++) {
	    if (my_group[i]!=num) continue;
	    seq_locus[n++]= i;
	}
 
    } else if (streq(str,"unlinked") || streq(str,"unlink")) { /* insane */
	if (num_groups==0) {
	    *why_not=ptr_to("no grouped markers in the data set");
	    return(TRUE);
	}
	num--;
	for (i=0; i<raw.num_markers; i++) {
	    if (my_group[i]!=num_groups-1) continue;
	    seq_locus[n++]= i;
	}
 
    } else if (nstreq(str,"order",5)) { /* insane */
	ordered_seq=TRUE;
	if (num_orders==0) {
	    *why_not=ptr_to("no ordered markers in the data set");
	    return(TRUE);
	} else if (sscanf(str+5,"%d",&num)!=1 || !irange(&num,1,num_orders)) {
	    *why_not=ptr_to("that order is not defined"); /* fix this */
	    return(TRUE);
	}
	num--;
	for (i=order_first[num]; i!=NO_LOCUS; i=order_next[i])
	  seq_locus[n++]= i;

    } else if (nstreq(str,"other",5)) { /* insane */
	ordered_seq=TRUE;
	if (num_orders==0) {
	    *why_not=ptr_to("no ordered markers in the data set");
	    return(TRUE);
	} else if (sscanf(str+5,"%d",&num)!=1 || !irange(&num,1,num_orders)) {
	    *why_not=ptr_to("that order is not defined"); /* fix this */
	    return(TRUE);
	}
	num--;
	for (i=unorder_first[num]; i!=NO_LOCUS; i=order_next[i])
	  seq_locus[n++]= i;

    } else if (streq(str,"none")) {
	n=0;
	return(TRUE);

    } else return(FALSE);
    

    /***** delete haplos, find contigs, and make a string from it *****/

#ifdef OBSOLETE
    for (i=0; i<n; i++) use_locus[i]=TRUE;
    /* ordered_seq means seq loci must be in ascending order, eg for haplo */
    if (use_haplotypes && !ordered_seq) {
	for (i=0; i<raw.num_markers; i++) use_locus[i]=FALSE;
	for (i=0; i<n; i++) {
	    locus= seq_locus[i];
	    if (!haplotyped(locus)) use_locus[locus]=TRUE;
	    else if (locus==haplo_first[locus]) use_locus[locus]=TRUE;
	    else use_locus[haplo_first[locus]]=TRUE;
	}
	n=0;
	for (i=0; i<raw.num_markers; i++) if (use_locus[i]) seq_locus[n++]=i;
    }
#endif

    last=seq_temp; last[0]='\0'; ordered_seq=TRUE; /* KLUDGE */
    for (i=0; i<n; i++) {
	if (ordered_seq) {
	    if (!print_names) sprintf(last, "%d", seq_locus[i] + 1);
	      else sprintf(last, "%s", raw.locus_name[seq_locus[i]]);
	} else {
	    contig=0; j=i;
	    while (j!=n-1 && seq_locus[j+1]==seq_locus[j]+1) { j++; contig++; }
	    if (contig<2) {
		if (!print_names) sprintf(last, "%d", seq_locus[i] + 1);
		  else sprintf(last, "%s", raw.locus_name[seq_locus[i]]);
	    } else {
		if (!print_names) 
		  sprintf(last, "%d-%d", seq_locus[i] + 1, seq_locus[j] + 1);
		  else sprintf(last, "%s-%s", raw.locus_name[seq_locus[i]],
                       raw.locus_name[seq_locus[j]]);
		i=j;
	    }
	}
	while (*last!='\0') ++last;
	*last=' '; ++last; *last='\0';
    }

    *seq=seq_temp; 
    return(TRUE);
}


void print_special_sequences()
{
    int i;
    char name[10];
    
    if (current_chrom!=NO_CHROM) {
	sprintf(ps, "Special Names (chromosome=%s):\n", chrom2str(current_chrom));
    } else sprintf(ps, "Special Names:\n");
    pr();

    print_special_sequence("changed",FALSE);
    print_special_sequence("new",TRUE);
    print_special_sequence("anchors",TRUE);
    print_special_sequence("assigned",TRUE);
    print_special_sequence("framework",TRUE);
    print_special_sequence("unassign",TRUE);

    for (i=0; i<num_groups-1; i++) /* groupN */
      { sprintf(name, "group%d", i + 1); print_special_sequence(name, FALSE); }
    print_special_sequence("unlinked",FALSE);

    for (i=0; i<num_orders; i++) { /* orderN */
	sprintf(name, "order%d", i + 1); print_special_sequence(name, TRUE);
	sprintf(name, "other%d", i + 1); print_special_sequence(name, FALSE);
    }

    print(" others:    all, none, <chromosome-name>, etc.\n");
}


void print_user_sequences()
{
    bool i, locus, any=FALSE, in_class=FALSE;
    char *seq, line[MAX_SEQ_LEN+99];

    for (Te=context[active_context]->named_sequences->list; Te!=NULL; 
	 Te=Te->next) {
	if (!any) { print("User Defined Names:\n"); any=TRUE; }
	sprintf(ps, " %s= ", Te->id.name); pad_to_len(ps, 12); pr(); 
	seq=Te->string;
	if (nullstr(seq)) strcpy(line,"none"); else strcpy(line,seq);
	if (len(line)>64) { truncstr(line,64); strcat(line,"..."); }
	print(line); nl(); 
    }
    for (i=0; i<NUM_CLASSES; i++) if (!nullstr(class_name[i])) {
	for (in_class=FALSE, locus=0; locus<raw.num_markers; locus++) {
	    if (use_haplotypes && haplotype_subordinate(locus)) continue;
	    if (class[locus]==i) { in_class=TRUE; break; }
	}
	if (!in_class) continue;
	if (!any) { print("User Defined Names:\n"); any=TRUE; }
	print_special_sequence(class_name[i],TRUE);
    }
    if (!any) print("No User Defined Names.\n");
}


void print_special_sequence(name,print_if_none)
char *name;
bool print_if_none;
{
    char *seq, *errmsg, line[MAX_SEQ_LEN+99];

    if (!is_a_special_sequence(name,&seq,&errmsg)) send(CRASH);
    /* or ? if (nullstr(errmsg)) */

    if (print_if_none || !nullstr(seq)) {
	sprintf(line, " %s= ", name); pad_to_len(line, 12); print(line);
	if (nullstr(seq)) strcpy(line,"none"); else strcpy(line,seq);
	if (len(seq)>64) { truncstr(line,64); strcat(line,"..."); }
	print(line); nl();
    }
}


void print_history_seqs(num_to_do)
int num_to_do;
{
    int i, start, any=FALSE;
    char *str, line[MAX_SEQ_LEN+99];

    if (context[active_context]->seq_history_num<=num_to_do) start=0;
      else start=context[active_context]->seq_history_num-num_to_do;

    for (i=start; i<context[active_context]->seq_history_num-1; i++) {
	get_numbered_entry(i,&str,context[active_context]->sequence_history);
	sprintf(ps, "#%d= ", i + 1); pad_to_len(ps, 6); pr();
	strcpy(line,str);
	if (len(str)>70) { truncstr(line,70); strcat(line,"..."); }
	print(line); nl(); any=TRUE;
    }
    if (!any) print("none\n");
    /* print("##=   "); print(seq_string); nl(); */
}
