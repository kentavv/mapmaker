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

#include "qtl.h"

void allocate_context(STATUS_CONTEXT *con);

#define MAX_CONTEXTS 10
STATUS_CONTEXT **context;
int active_context, num_contexts;

void context_init(void) {
    context = NULL;
    run {
            parray(context, MAX_CONTEXTS, STATUS_CONTEXT);

            num_contexts = 1;
            active_context = 0;

            allocate_context(context[0]);

        } on_exit {
        relay_messages;
    }
}


void
allocate_context(STATUS_CONTEXT *con) {
    TABLE *seqhist, *names;

    seqhist = names = NULL;
    run {
            con->trait = 0;
            seqhist =
                    allocate_table(MAX_HISTORY_SEQS, SEQ_LEN, CANT_EXPAND, INDEX_BY_NUMBER);
            con->sequence_history = seqhist;
            con->seq_history_num = next_entry_number(con->sequence_history);
            names = allocate_table(MAX_NAMED_SEQS, SEQ_LEN, EXPANDS_BY(MAX_NAMED_SEQS),
                                   INDEX_BY_NAME);
            con->named_sequences = names;
        } on_exit {
        relay_messages;
    }
}
