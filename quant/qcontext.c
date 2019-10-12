/******************************************************************************

  ####    ####    ####   #    #   #####  ######  #    #   #####           ####
 #    #  #    #  #    #  ##   #     #    #        #  #      #            #    #
 #    #  #       #    #  # #  #     #    #####     ##       #            #
 #  # #  #       #    #  #  # #     #    #         ##       #     ###    #
 #   #   #    #  #    #  #   ##     #    #        #  #      #     ###    #    #
  ### #   ####    ####   #    #     #    ######  #    #     #     ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#define INC_LIB
#define INC_SHELL
#define INC_TABLE
#define INC_CALLQCTM
#define INC_QTOPLEVEL
#include "qtl.h"

void allocate_context (STATUS_CONTEXT *con);

#define MAX_CONTEXTS 10
STATUS_CONTEXT **context;
int active_context, num_contexts;

void 
context_init (void)
{
    context = NULL;
    run {
	parray(context,MAX_CONTEXTS,STATUS_CONTEXT);

	num_contexts = 1;
	active_context = 0;

	allocate_context(context[0]);

    } on_exit {
	relay_messages;
    }
}


void 
allocate_context (STATUS_CONTEXT *con)
{
    TABLE *seqhist,*names;

    seqhist = names = NULL;
    run {
	con->trait = 0;
	seqhist = 
	  allocate_table(MAX_HISTORY_SEQS,SEQ_LEN,CANT_EXPAND,INDEX_BY_NUMBER);
	con->sequence_history = seqhist;
	con->seq_history_num= next_entry_number(con->sequence_history);
	names=allocate_table(MAX_NAMED_SEQS,SEQ_LEN,EXPANDS_BY(MAX_NAMED_SEQS),
			     INDEX_BY_NAME);
	con->named_sequences = names;
    } on_exit {
	relay_messages;
    }
}

void 
free_context (STATUS_CONTEXT *con)
{
    free_table(con->named_sequences);
    free_table(con->sequence_history);
}


bool 
change_context (int new_context)
{
    if(context[new_context] != NULL) {
	active_context = new_context;
	/* set the globals */
	trait = context[active_context]->trait;
	return(TRUE);
    } else {
	return(FALSE);
    }
}


bool 
create_new_context (int new_context)
{
    if(context[new_context] != NULL || new_context > MAX_CONTEXTS)
      return(FALSE);
    
    allocate_context(context[new_context]);
    
    /* take current values as defaults */
    context[new_context]->trait = context[active_context]->trait;

    return(TRUE);
}


void 
kill_context (STATUS_CONTEXT *con, bool save_it)
{
    char *err;

    if(save_it) {
	for(Te=con->named_sequences->list; Te!=NULL; Te=Te->next) {
	    if(!name_sequence(Te->id.name,Te->string,&err))
	      error(err);
	}
    }
    free_context(con);
}



