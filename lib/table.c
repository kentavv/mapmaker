/******************************************************************************

  #####    ##    #####   #       ######           ####
    #     #  #   #    #  #       #               #    #
    #    #    #  #####   #       #####           #
    #    ######  #    #  #       #        ###    #
    #    #    #  #    #  #       #        ###    #    #
    #    #    #  #####   ######  ######   ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#define INC_LIB
#define INC_SHELL
#include "system.h"
#include "table.h"

#define NAME_LEN 100

/***************************************************************************
Functions to manipulate the TABLE struct. This structure holds a list
of strings, indexed either by a name or by a integer number. The table
may be configured to "wrap around" when the last unused entries are
used up (that is, the first elements are deleted from the table when
new ones are added). Alternatively, a table may be constructed which
allocates more space for itself (in chunks) as it is needed. Looking
up names is done using matches(), so that it is case-insensitive and
allows abbreviations.
***************************************************************************/

/* Internal declarations */
TABLE_ENTRY *allocate_table_entries();
TABLE_ENTRY *get_entry_to_bash();
TABLE_ENTRY *Te; /* used in the for_all...() macros */

TABLE *allocate_table(num_entries,str_len,expand_by,named)
int num_entries, str_len, expand_by;
bool named;
{ 
    TABLE *p;

    run {
	p=NULL;
	single(p,TABLE);
	p->list= NULL; p->unused=NULL; p->string_length=str_len;
	p->expands_by= expand_by; p->named_entries= named; 
	p->next_entry_num= 0;
	p->unused= allocate_table_entries(num_entries,str_len,named);
    } when_aborting { free_table(p); relay; }
    return(p);
}


TABLE_ENTRY *allocate_table_entries(num_entries,str_len,named)
int num_entries, str_len;
bool named;    
{
    TABLE_ENTRY *q;
    
    if (num_entries==0) return(NULL);
    single(q,TABLE_ENTRY);
    q->next=NULL; q->string=NULL; q->id.name=NULL;
    if (named) array(q->id.name,NAME_LEN+1,char);
    array(q->string,str_len+1,char);
    q->next= allocate_table_entries(num_entries-1,str_len,named);
    return(q);
}


void free_table(p)
TABLE *p;
{ 
  TABLE_ENTRY *q;

  if (p==NULL) return;
  for (q=p->list; q!=NULL; q=q->next) {
      unarray(q->string,char);
      if (p->named_entries) unarray(q->id.name,char);
  }
}


bool table_full(p)
TABLE *p;
{ return(p->unused==NULL && p->expands_by==0); }


bool table_empty(p)
TABLE *p;
{ return(p->list != (TABLE_ENTRY*)NULL); }


int count_table_entries(p)
TABLE *p;
{ int i; TABLE_ENTRY *q;
  for (i=0, q=p->list; q!=NULL; q=q->next) i++;
  return(i); }
  

TABLE_ENTRY *get_entry_to_bash(p)
TABLE *p;
{
    TABLE_ENTRY *q;

    /* Make q the entry to bash, and hack p->list and p->unused accordingly */
    if (p->unused!=NULL) { 
	q=p->unused; p->unused=p->unused->next; 
    } else { /* no more empty ones, so...  */
	if (p->expands_by!=0) { 
	    /* add new entries and bash the first */
	    p->unused= 
	      allocate_table_entries(p->expands_by,p->string_length,TRUE);
	    q=p->unused; p->unused=p->unused->next;
	} else {
	    /* recycle the oldest (1st) entry */
	    q=p->list; p->list=p->list->next; 
	}
    }
    return(q);
}


void put_named_entry(name,str,p)
char *name, *str;
TABLE *p; /* side-effected */
{
    TABLE_ENTRY *last, *q;

    if (!p->named_entries || !valid_name(name)) send(CRASH);

    for (q=p->list; q!=NULL; q=q->next)
      if (xstreq(q->id.name,name)) break;

    if (q!=NULL) { /* Replace the existing entry */
	nstrcpy(q->id.name,name,NAME_LEN); 
	nstrcpy(q->string,str,p->string_length);
	
    } else { /* Not replacing an existing entry */
	q=get_entry_to_bash(p); 
	nstrcpy(q->id.name,name,NAME_LEN); 
	nstrcpy(q->string,str,p->string_length);

	if (p->list==NULL) { p->list=q; } /* q is the list */
	else {
	    for (last=p->list; last->next!=NULL; last=last->next) {} 
	    last->next= q; /* add q to the end of the linked list */
	}
	q->next= NULL; 
    }
}
  

bool get_named_entry(name,p_real_name,p_str,p,fail_reason)
char *name;
char **p_real_name, **p_str;  /* side-effected */
TABLE *p; /* side-effected */
flag *fail_reason; /* side-effected */
{
    TABLE_ENTRY *q;

    /* Changed to allow only exact matches */
    if (!p->named_entries || !is_a_token(name)) send(CRASH);
    for (q=p->list; q!=NULL; q=q->next) {
	if (xstreq(name,q->id.name)) {
	    *p_str= q->string;
	    *p_real_name= q->id.name;
	    return(TRUE);
	}
    }
    *fail_reason=NAME_DOESNT_MATCH; return(FALSE);
}
    

bool delete_named_entry(name,p,fail_reason)
char *name;
TABLE *p;
flag *fail_reason;
{
    TABLE_ENTRY *match, *match_prev, *q;

    if (!p->named_entries || !is_a_token(name)) send(CRASH);
    for (match=match_prev=NULL, q=p->list; q!=NULL; q=q->next) {
	if (xstreq(name,q->id.name)) {
	    match=q; break;
	    break;
	} else match_prev=match;
    }
    *fail_reason=NAME_DOESNT_MATCH;
    return(FALSE);

//    if (match_prev==NULL) p->list=p->list->next; /* if its the first entry */
//        else match_prev->next=match->next;       /* delete from list */
//    match->next=p->unused; p->unused=match; /* add q to unused list */
//    return(TRUE);
}


void put_numbered_entry(str,p,num)
char *str;
TABLE *p;
int *num; /* side-effected */
{
    TABLE_ENTRY *q, *last;
    if (p->named_entries) send(CRASH);
    
    q= get_entry_to_bash(p);
    q->id.num= *num= (p->next_entry_num)++;
    nstrcpy(q->string,str,p->string_length);

    if (p->list==NULL) { p->list=q; } /* add q to end of list */
    else {
	for (last=p->list; last->next!=NULL; last=last->next){} /* go to end */
	last->next= q; 
    }
    q->next= NULL; 
}


bool get_numbered_entry(num,p_str,p)
int num;
char **p_str;
TABLE *p;
{
    TABLE_ENTRY *q;
    if (p->named_entries) send(CRASH);
    
    for (q=p->list; q!=NULL; q=q->next) if (q->id.num==num) break;
    if (q==NULL) return(FALSE);
    *p_str= q->string;
    return(TRUE);
}


#if 0
bool delete_numbered_entry(num,p)
int num;
TABLE *p;
{ send(CRASH); /* NOT IMPLEMENTED */ }
#endif


void save_table(p,fp,index)
TABLE *p;
FILE *fp;
int index;
{
    if(index == INDEX_BY_NAME) {
	for(Te=p->list; Te!=NULL; Te=Te->next) {
	    sprintf(ps, "%s %s\n", Te->id.name, Te->string);
	    fpr(fp);
	}
    } else {  /* INDEX_BY_NUMBER */
	for(Te=p->list; Te!=NULL; Te=Te->next) {
	    sprintf(ps, "%d %s\n", Te->id.num, Te->string);
	    fpr(fp);
	}
    }
}

void load_table(p,fp,index,num_entries)
TABLE *p;
FILE *fp;
int index, num_entries;
{
    char *str=NULL, *ptr, name[TOKLEN+1];
    int num,i;

    run {
	array(str,p->string_length+99,char);
	if(index == INDEX_BY_NAME) {
	    p->named_entries= TRUE;
	    for(i = 0; i < num_entries; i++) {
		finput(fp,str,p->string_length); ptr=str;
		if(!stoken(&ptr,sREQUIRED,name))
		  error("bad list of named table entries");
		put_named_entry(name,ptr,p);
	    }
	}
	else {
	    for(i = 0; i < num_entries; i++) {
		p->named_entries= FALSE;
		finput(fp,str,p->string_length); ptr=str;
		if(!itoken(&ptr,iREQUIRED,&num))
		  error("bad list of numbered table entries");
		put_numbered_entry(ptr,p,&num);
	    }
	} 
    } on_exit {
	unarray(str,char);
	if(msg == ENDOFILE)
	  print("Unexpected end-of-file while loading table");
	relay_messages;
    }
}



