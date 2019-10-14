/******************************************************************************

  ####    #   #   ####            ####   #    #  #####    ####            ####
 #         # #   #               #    #  ##  ##  #    #  #               #    #
  ####      #     ####           #       # ## #  #    #   ####           #
      #     #         #          #       #    #  #    #       #   ###    #
 #    #     #    #    #          #    #  #    #  #    #  #    #   ###    #    #
  ####      #     ####  #######   ####   #    #  #####    ####    ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#include "mapm.h"

/* Auxilliary stuff for load/save/prep */
FILE *try_to_open(char *name, int mode, char *ext, bool prev_data) /* return fp or send error */;

void
try_to_unload(FILE *fp, /* file to READ just so we can close it if an error occurs */ bool ask_first, bool do_save, bool do_unload, bool genos_too);

bool is_an_old_sequence(char *str, char **seq, char **why_not);

bool is_a_named_sequence(char *str, char **seq);  /* internal use only */

#define SAVE_ERR1 \
 "\nAn error occured while saving data: previous data remain loaded\n"
#define SAVE_ERR2 \
"If any files became corrupted, the old data should be available in\nbackup file(s).\n"
#define LOAD_CANTOPEN "unable to read from file '%s'%s"
#define SAVE_CANTOPEN "unable to write to file '%s'%s"
#define LOAD_ERROR "\nan error occured while loading data%s\n"
#define PREV_LOST  ": previous data were lost"
#define PREV_KEPT  "\nprevious data remain loaded"
#define BAD_FILENAME "illegal filename '%s'"
#define DATA_LOADED "data from '%s' are loaded\n  %s\n"
#define LOAD_MAYBE  "data MAY be correctly partially loaded\n"
#define NO_AUTOSAVE "warning: 'auto save' is off - not saving data\n"
#define RAW_EXTENSION "raw file has a name which will be overwritten"

bool just_prepared = FALSE;

FILE *try_to_open(char *name, int mode, char *ext, bool prev_data) {
    /* return fp or send error */
    FILE *fp;
    run {
            if (!make_filename(name, mode, ext)) {
                sprintf(ps, BAD_FILENAME, name);
                error(ps);
            }
            fp = open_file(name, READ);
        } except_when(CANTOPEN) {
        sprintf(ps, LOAD_CANTOPEN, name, (prev_data ? PREV_KEPT : ""));
        error(ps);
    }
    return (fp);
}


void try_to_unload(FILE *fp, /* file to READ just so we can close it if an error occurs */ bool ask_first, bool do_save, bool do_unload,
                   bool genos_too) {
    bool doit;
    char token[TOKLEN + 1];

    run {
            if (do_save && ask_first && interactive && !redirecting_input) {
                run getln("save current data set first? [yes] ");
                    except {
                        when ENDOINPUT:
                            ln[0] = 'y';
                            break;
                        default:
                            relay;
                    }
                if (stoken(&ln, sREQUIRED, token) && matches(token, "no"))
                    doit = FALSE;
                else doit = TRUE;
            } else doit = TRUE;
            if (do_save && doit) do_save_data(raw.filename, genos_too);
            else if (do_save) print("ok - not saving data\n");
            else /* !do_save */ print(NO_AUTOSAVE);
            if (do_unload) do_unload_data();

        } on_error {
        if (msg == INTERRUPT) relay;
        print(SAVE_ERR1);
        trapped_msg();
        print(SAVE_ERR2);
        if (fp != NULL && msg != CANTCLOSE) close_file(fp);
        abort_command();
    }
}


void try_to_load(FILE *fp, char *name, bool prev_data, bool raw) {
    char run_file[PATH_LENGTH + 1];

    run {
            /** if (!raw) do_load_data(fp,name,raw); **/
            do_load_data(fp, name, raw);
            close_file(fp);
            if (raw) {
                strcpy(run_file, name);
                make_filename(run_file, FORCE_EXTENSION, PREP_EXT);
                if (!redirect_input(run_file, FALSE)) {
                    sprintf(ps, "unable to run file '%s'... skipping initialization\n",
                            run_file);
                    pr();
                    try_to_unload(fp, FALSE, TRUE, FALSE, TRUE); /* saves data */
                } else {
                    sprintf(ps, "running initialization commands from file '%s'...\n",
                            run_file);
                    pr();
                }
            }
        } on_error {
        sprintf(ps, LOAD_ERROR, (prev_data ? PREV_LOST : ""));
        pr();
        if (msg != CANTCLOSE) close_file(fp);
        if (msg == BADDATA) {
            print(BADDATA_reason);
            nl();
        } else trapped_msg();
        if (data_loaded()) print(LOAD_MAYBE);
        abort_command();
    }
}


void mapm_data_info(FILE *fp) {
    if (!data_loaded()) return;
    sprintf(ps, DATA_LOADED, raw.filename, data_info(TRUE));
    fpr(fp);
}


/* Now the real stuff */

command new_load_data(void) {
    FILE *fp = NULL;
    char name[PATH_LENGTH + 1];
    bool prev_data;

    use_uncrunched_args();
    mapm_ready(MAYBE_DATA, 0, 0, NULL);
    get_one_arg(stoken, "", name);
    prev_data = data_loaded();

    if (nullstr(name)) {
        if (!prev_data) print("no data are loaded\n");
        else {
            sprintf(ps, DATA_LOADED, raw.filename, data_info(TRUE));
            pr();
        }
    } else {
        fp = try_to_open(name, FORCE_EXTENSION, DATA_EXT, prev_data);
        if (prev_data) try_to_unload(fp, TRUE, auto_save, TRUE, just_prepared);
        try_to_load(fp, name, prev_data, FALSE);
        just_prepared = FALSE;
    }
}


command new_prepare(void) {
    char name[PATH_LENGTH + 1];
    FILE *fp = NULL;
    bool prev_data;

    use_uncrunched_args();
    mapm_ready(MAYBE_DATA, 0, 0, NULL);
    get_one_arg(stoken, sREQUIRED, name);
    prev_data = data_loaded();

    /* end=len(name)-1;
       if (streq(&name[end-4],DATA_EXT) || streq(&name[end-4],DATA_OLD) ||...
       error(RAW_EXTENSION); */

    fp = try_to_open(name, DEFAULT_EXTENSION, RAW_EXT, prev_data);
    if (prev_data) try_to_unload(fp, TRUE, auto_save, TRUE, just_prepared);
    try_to_load(fp, name, prev_data, TRUE);
    just_prepared = TRUE;
}


command new_save_data(void) {
    char name[PATH_LENGTH + 1], old[PATH_LENGTH + 1];
    bool new_name;

    use_uncrunched_args();
    mapm_ready(ANY_DATA, 0, 0, NULL);
    get_one_arg(stoken, "", name);
    new_name = !nullstr(name);
    /* want to change this so it sets raw.filename only AFTER it writes OK */

    if (new_name && !make_filename(name, FORCE_EXTENSION, DATA_EXT)) {
        sprintf(ps, BAD_FILENAME, name);
        error(ps);
    }

    run {
            if (new_name) {
                strcpy(old, raw.filename);
                strcpy(raw.filename, name);
            }
            try_to_unload(NULL, FALSE, TRUE, FALSE, (new_name || just_prepared));
            just_prepared = FALSE;
        } on_error {
        if (new_name) strcpy(name, old);
        relay_messages;
    }
}


#define CLASS_EXISTS  "class '%s' is already defined"
#define CLASS_CANT_BE "'%s' is already defined"
#define NOT_A_CLASS "'%s' is not a defined class\nEither use 'make class' first or select a defined class name."
#define CANT_MAKE_CLASS "warning: unable to make class named '%s'\n%s\n"

command set_class(void) {
    int classnum, i, n_loci, *locus = NULL;
    char name[TOKLEN + 1];

    mapm_ready(ANY_DATA, 1, LIST_SEQ, NULL);
    get_one_arg(stoken, sREQUIRED, name);
    if (!isa_class(name, &classnum)) {
        sprintf(ps, NOT_A_CLASS, name);
        error(ps);
    }

    run {
            alloc_list_of_all_loci(seq, &locus, &n_loci);
            for (i = 0; i < n_loci; i++) class[locus[i]] = classnum;
            sprintf(ps, "markers in sequence will now be considered in class '%s'\n",
                    class_name[classnum]);
            pr();
        } on_exit {
        unarray(locus, int);
        relay_messages;
    }
}


command make_classes(void) {
    char name[TOKLEN + 1], *errmsg;

    mapm_ready(ANY_DATA, 0, 0, NULL);

    while (stoken(&uncrunched_args, sREQUIRED, name))
        if (!make_new_class(name, &errmsg)) {
            sprintf(ps, CANT_MAKE_CLASS, name, errmsg);
            pr();
        }

    print("classes defined: ");
    print_class_names();
    nl();
}


command set_age(void) {
    int mod, i, n_loci, *locus = NULL;
    char token[TOKLEN + 1];

    mapm_ready(ANY_DATA, 1, LIST_SEQ, NULL);

    if (!stoken(&args, sREQUIRED, token)) usage_error(1);
    else if (streq(token, "new")) mod = TRUE;
    else if (streq(token, "old")) mod = FALSE;
    else usage_error(1);

    run {
            alloc_list_of_all_loci(seq, &locus, &n_loci);
            for (i = 0; i < n_loci; i++) modified[locus[i]] = mod;
            sprintf(ps, "markers in sequence will now be considered %s\n",
                    (mod ? "new" : "old"));
            pr();

        } on_exit {
        unarray(locus, int);
        relay_messages;
    }
}


#define ERROR_PROB_IS \
"markers in sequence now have apriori error probability %.2lf percent\n"

command set_error_rate(void) {
    int i, n_loci, *locus = NULL;
    real rate, prob;

    mapm_ready(ANY_DATA, 1, LIST_SEQ, NULL);

    if (!rtoken(&args, rREQUIRED, &rate) || !rrange(&rate, 0.0, 10.0))
        usage_error(1);
    prob = rate / 100.0;

    run {
            alloc_list_of_all_loci(seq, &locus, &n_loci);
            for (i = 0; i < n_loci; i++) error_rate[locus[i]] = prob;
            sprintf(ps, ERROR_PROB_IS, rate);
            pr();

        } on_exit {
        unarray(locus, int);
        relay_messages;
    }
}


#define NOTE_FORM "%10s: %s\n"

command make_note(void) {
    int n, n_loci, i, *locus = NULL;
    char *name = get_temp_string(), *rest, *errmsg;

    mapm_ready(ANY_DATA, 0, 0, &n_loci);

    /* In here was: if (note[marker-1]==NULL || (strcmp(note[marker-1],"")==0))
       FIX: IS THERE ANY CHANCE THAT A NOTE CAN BE NULL? */

    run {

            if (nullstr(uncrunched_args)) {
                if (n_loci == 0) error("note args error"); /* FIX */
                alloc_list_of_all_loci(seq, &locus, &n_loci);
                print("\nNotes:\n");
                for (i = 0; i < n_loci; i++) {
                    if (nullstr(note[locus[i]]))
                        sprintf(ps, NOTE_FORM, loc2str(locus[i]), "<no note>");
                    else
                        sprintf(ps, NOTE_FORM, loc2str(locus[i]), note[locus[i]]);
                    pr();
                }

            } else if (!split_string(uncrunched_args, &rest, ':')) {
                stoken(&uncrunched_args, sREQUIRED, name);
                if (!is_a_locus(name, &n, &errmsg) || !nullstr(uncrunched_args))
                    error("note args error"); /* FIX */
                print("\nNote:\n");
                if (nullstr(note[n]))
                    sprintf(ps, NOTE_FORM, loc2str(n), "<no note>");
                else

                    sprintf(ps, NOTE_FORM, loc2str(n), note[n]);
                pr();

            } else {
                stoken(&uncrunched_args, sREQUIRED, name);
                if (!is_a_locus(name, &n, NULL) || !nullstr(uncrunched_args))
                    error("note args error"); /* FIX */

                if (len(rest) > MAX_NOTE_LEN) print("warning: note too long\n");
                nstrcpy(note[n], rest, MAX_NOTE_LEN);
            }

        } on_exit {
        unarray(locus, int);
    }
}


/***************************** Sequence Commands *****************************/

command sequence(void) {
    int *locus = NULL, num_loci;
    use_uncrunched_args();
    mapm_ready(ANY_DATA, 0, 0, NULL);

    if (nullstr(args)) {
        print_sequence();
        return;
    }
    run {
            nstrcpy(new_seq, args, MAX_SEQ_LEN);
            set_current_seq(new_seq, FALSE);
            print_sequence();
            if (alloc_list_of_all_loci(seq, &locus, &num_loci)) {
                crunch_locus_list(locus, &num_loci, CRUNCH_WARNINGS, TRUE, MAYBE);
                unarray(locus, int);
            }
        } except_when(BADSEQ) print_badseq();
}


command expand_sequence(void) {
    int *locus = NULL, num_loci;
    use_uncrunched_args();
    mapm_ready(ANY_DATA, 0, 0, NULL);

    run {
            if (nullstr(args)) {
                if (current_chrom == NO_CHROM) strcpy(new_seq, seq_string);
                else sprintf(new_seq, "%s: %s", chrom2str(current_chrom), seq_string);
            } else nstrcpy(new_seq, args, MAX_SEQ_LEN);
            set_current_seq(new_seq, TRUE);
            print_sequence();
            if (alloc_list_of_all_loci(seq, &locus, &num_loci)) {
                crunch_locus_list(locus, &num_loci, CRUNCH_WARNINGS, TRUE, MAYBE);
                unarray(locus, int);
            }
        } except_when(BADSEQ) print_badseq();
}


#define BADEDITEDSEQ "An illegal sequence was specified.\nsequence= %s\n"

command edit_sequence(void) {
    char name[TOKLEN + 1], prompt[TOKLEN + 1], *value, *err;
    int *locus = NULL, num_loci;
    bool set_seq = FALSE;
    mapm_ready(ANY_DATA, 0, 0, NULL);

    /* edit_line() will send an error if it can't be run now */
    if (nullstr(args)) {
        /* if edit_line cannot be called, an error is sent */
        set_seq = TRUE;
        strcpy(prompt, "sequence= ");
        if (current_chrom == NO_CHROM) strcpy(new_seq, seq_string);
        else sprintf(new_seq, "%s: %s", chrom2str(current_chrom), seq_string);
    } else {
        get_one_arg(stoken, sREQUIRED, name);
        if (is_an_old_sequence(name, &value, &err)) {
            set_seq = TRUE;
            strcpy(prompt, "sequence= ");
            strcpy(new_seq, value);
        } else if (is_a_named_sequence(name, &value)) { /* NOT special! */
            set_seq = FALSE;
            sprintf(prompt, "%s= ", name);
            strcpy(new_seq, value);
        } else {
            sprintf(ps, "'%s' is not a user-defined name or an old sequence", name);
            error(ps);
        }
    }
    edit_line(prompt, new_seq, MAX_SEQ_LEN, new_seq);
    if (set_seq) {
        run {
                set_current_seq(new_seq, FALSE);
                if (alloc_list_of_all_loci(seq, &locus, &num_loci)) {
                    crunch_locus_list(locus, &num_loci, CRUNCH_WARNINGS, TRUE, MAYBE);
                    unarray(locus, int);
                }
            } except_when(BADSEQ) print_badseq();
    } else {
        if (!name_sequence(name, new_seq, &err, FALSE)) error(err);
    }
    print("ok\n");
}


command show_seq_history(void) {
    int num_to_do;

    mapm_ready(ANY_DATA, 0, 0, NULL);
    get_one_arg(itoken, 20, &num_to_do);

    print("Previous Sequences:\n");
    print_history_seqs(num_to_do);
}


command let(void) {
    char name[TOKLEN + 1], *seq, *err;
    use_uncrunched_args();
    mapm_ready(ANY_DATA, 0, 0, NULL);

    if (nullstr(args) || !split_arglist(&seq, '=')) usage_error(0);
    if (!stoken(&args, sREQUIRED, name) || !nullstr(args)) usage_error(1);
    if (nullstr(seq)) strcpy(seq, "none");
    strcpy(new_seq, seq);
    if (!name_sequence(name, new_seq, &err, FALSE)) error(err);
    sprintf(ps, "%s= %s\n", name, new_seq);
    pr();
}


command let_expanding(void) {
    char name[TOKLEN + 1], *seq, *err;
    bool set_me = FALSE;

    use_uncrunched_args();
    mapm_ready(ANY_DATA, 0, 0, NULL);

    if (nullstr(args)) usage_error(0);
    if (split_arglist(&seq, '=')) set_me = TRUE;
    if (!stoken(&args, sREQUIRED, name) || !nullstr(args)) usage_error(1);
    if (set_me) {
        if (nullstr(seq)) strcpy(seq, "none");
        if (!name_sequence(name, seq, &err, TRUE)) error(err);
        if (!is_a_sequence(name, &seq, &err)) send(CRASH);
    } else {
        if (!is_a_sequence(name, &seq, &err)) {
            sprintf(ps, "name '%s' is not defined", name);
            error(ps);
        }
        /* but it may be a special sequence: name_sequence will catch this */
    }
    sprintf(ps, "%s= %s\n", name, seq);
    pr();
}


command names(void) {
    mapm_ready(ANY_DATA, 0, 0, NULL);
    nomore_args(0);

    print_special_sequences();
    print_user_sequences();
}


command forget(void) {
    char *errmsg, *name = get_temp_string();

    mapm_ready(ANY_DATA, 0, 0, NULL);
    get_one_arg(stoken, sREQUIRED, name);

    if (!unname_sequence(name, &errmsg)) {
        sprintf(ps, "cannot forget name '%s'\n%s", name, errmsg);
        error(ps);
    } else print("ok\n");
}


command new_delete(void) {
    int i, j, found, num_seq_tokens, k, *locus = NULL, num_loci;
    int *seq_locus = NULL, seq_loci;
    char locus_name[TOKLEN + 1], locus_num[TOKLEN + 1], locus_plus[TOKLEN + 1];

    mapm_ready(ANY_DATA, 1, MAYBE_PERM, NULL);
    if (nullstr(args)) usage_error(0);
    parse_locus_args(&locus, &num_loci); /* error if fails */

    run {
            strcpy(new_seq, seq_string); /* expand_seq_names(new_seq); */
            tokenize_seq(new_seq, seq_tokens, &num_seq_tokens);
            for (k = 0; k < num_loci; k++) {
                strcpy(locus_name, raw.locus_name[locus[k]]);
                sprintf(locus_plus, "%s+", raw.locus_name[locus[k]]);
                sprintf(locus_num, "%d", locus[k] + 1);

                i = 0;
                found = FALSE;
                while (i < num_seq_tokens) {
                    if (xstreq(locus_num, seq_tokens[i]) ||
                        xstreq(locus_plus, seq_tokens[i]) ||
                        xstreq(locus_name, seq_tokens[i])) {
                        for (j = i; j < num_seq_tokens - 1; j++)
                            strcpy(seq_tokens[j], seq_tokens[j + 1]);
                        strcpy(seq_tokens[j], "");
                        num_seq_tokens--;
                        found = TRUE;
                        break;
                    }
                    i++;
                }
                if (!found) {
                    sprintf(ps, "warning: locus %s (%s) not found in the sequence\n",
                            locus_num, locus_name);
                    pr();
                }
            }
            untokenize_seq(new_seq, seq_tokens, num_seq_tokens);
            set_current_seq(new_seq, FALSE);
            if (alloc_list_of_all_loci(seq, &seq_locus, &seq_loci))
                crunch_locus_list(seq_locus, &seq_loci, CRUNCH_WARNINGS, TRUE, MAYBE);
            print_sequence();
        } on_exit {
        if (msg == BADSEQ) print_badseq();
        unarray(locus, int);
        unarray(seq_locus, int);
    }
}


command new_append(void) {
    int *locus = NULL, num_loci;

    mapm_ready(ANY_DATA, 1, MAYBE_PERM, NULL);
    if (nullstr(args)) usage_error(0);
    despace(args);

    run {
            strcpy(new_seq, seq_string);
            strcat(new_seq, " ");
            strcat(new_seq, args);
            set_current_seq(new_seq, FALSE);
            print_sequence();
            if (alloc_list_of_all_loci(seq, &locus, &num_loci)) {
                crunch_locus_list(locus, &num_loci, CRUNCH_WARNINGS, TRUE, MAYBE);
                unarray(locus, int);
            }
        } except_when(BADSEQ) print_badseq();
}


command new_insert(void) {
    int i, j, locus, found, num_seq_tokens, append_at_top;
    char locus_name[TOKLEN + 1], locus_num[TOKLEN + 1], locus_plus[TOKLEN + 1];
    char *appendage, name[TOKLEN + 1], *errmsg;
    int *loci = NULL, num_loci;

    mapm_ready(ANY_DATA, 1, MAYBE_PERM, NULL);
    if (nullstr(args)) usage_error(0);
    if (!split_arglist(&appendage, ':') || nullstr(appendage) ||
        !stoken(&args, "", name) || !nullstr(args))
        usage_error(1);
    if (len(appendage) > TOKLEN) error("sequence to insert is too long");

    if (!nullstr(name)) {
        append_at_top = FALSE;
        if (!is_a_locus(name, &locus, &errmsg)) error(errmsg); /* do better */
        strcpy(locus_name, raw.locus_name[locus]);
        sprintf(locus_num, "%d", locus + 1);
        sprintf(locus_plus, "%s+", raw.locus_name[locus]);
    } else append_at_top = TRUE;

    run {
            tokenize_seq(new_seq, seq_tokens, &num_seq_tokens);
            i = 0;
            found = FALSE;
            if (append_at_top) {
                found = TRUE;
                i = -1;
            } /* see below */
            else
                while (!found && i < num_seq_tokens) {
                    if (streq(locus_num, seq_tokens[i]) ||
                        xstreq(locus_plus, seq_tokens[i]) ||
                        xstreq(locus_name, seq_tokens[i])) {
                        found = TRUE;
                        break;
                    } else i++;
                }
            if (!found) {
                sprintf(ps, "locus %s (%s) not found in the sequence",
                        locus_num, locus_name);
                error(ps);
            }

            for (j = num_seq_tokens - 1; j > i; j--)
                strcpy(seq_tokens[j + 1], seq_tokens[j]);
            /* this is severly broken, if len(appendage) > TOKLEN */
            strcpy(seq_tokens[i + i], appendage);
            num_seq_tokens++;

            untokenize_seq(new_seq, seq_tokens, num_seq_tokens);
            set_current_seq(new_seq, FALSE);
            print_sequence();
            if (alloc_list_of_all_loci(seq, &loci, &num_loci)) {
                crunch_locus_list(loci, &num_loci, CRUNCH_WARNINGS, TRUE, MAYBE);
                unarray(loci, int);
            }
        } except_when(BADSEQ) print_badseq();
}


command translate(void) {
    int i, num_loci, *locus = NULL, source;
    char c;

    mapm_ready(ANY_DATA, MAYBE_SEQ, UNCRUNCHED_LIST, &num_loci);
    run {
            if (!nullstr(args)) {
                parse_locus_args(&locus, &num_loci); /* error if fails */
                if (num_loci == 0) print("no loci\n");
                source = IN_ARGS;
            } else {
                if (!alloc_list_of_all_loci(seq, &locus, &num_loci))
                    error(NEED_SEQ_OR_ARGS);
                /* print_sequence(); */
                source = IN_SEQ;
            }
            crunch_locus_list(locus, &num_loci, CRUNCH_WARNINGS, ANY_CHROMS, source);
            for (i = 0; i < num_loci; i++) {
                c = ' ';
                if (use_haplotypes) {
                    if (haplotype_subordinate(locus[i])) c = '*';
                    else if (haplotyped(locus[i])) c = '+';
                }
                sprintf(ps, "%4d %s%c", locus[i] + 1, raw.locus_name[locus[i]], c);
                pad_to_len(ps, 14);
                pr();
                if (i == num_loci - 1 || i % 5 == 4) nl(); else print("  ");
            }
        } on_exit {
        unarray(locus, int);
        relay_messages;
    }
}
