/******************************************************************************

 #    #    ##       #    #    #           ####
 ##  ##   #  #      #    ##   #          #    #
 # ## #  #    #     #    # #  #          #
 #    #  ######     #    #  # #   ###    #
 #    #  #    #     #    #   ##   ###    #    #
 #    #  #    #     #    #    #   ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#include "mapm.h"

void setup_commands(void);

void make_help_entries(void); /* move these to lib */
//extern char help_filename[];

int main(int argc, char *argv[]) {
    char help_filename[PATH_LENGTH + 1];
    FILE *fp;

#ifdef THREAD
//    quit_if_not_single_invocation(argv[0]);
    threads_init();
#endif

    custom_lib_init();
    get_cmd_line_args(&argc, argv); /* side-effects file_arg vars */
    tty_hello();
    seedrand(RANDOM);

#ifdef HELP_PATH
    strcpy(help_filename,HELP_PATH "mapmaker");
#else
    strcpy(help_filename, "mapmaker");
#endif

#ifdef THREAD
    shell_init("MAPMAKER/EXP", "3.0b (Parallel, 2019)", "1987-1992", help_filename);
#else
    shell_init("MAPMAKER/EXP","3.0b (2019)","1987-1992",help_filename);
#endif
    banner();

    photo_banner_hook = mapm_data_info;
    quit_save_hook = mapm_save_on_exit;

    data_init();
    npt_cmds_init();
    state_init();
    map_init();
    sequence_init();

    setup_commands();
    make_help_entries();
    reset_state();

#ifdef HAVE_CEPH
    prep_init();
#endif

    if (!nullstr(file_arg[PHOTO_FILE_ARG])) {
        run {
                nl();
                if ((append_it && !photo_to_file(file_arg[PHOTO_FILE_ARG], APPEND)) ||
                    (!append_it && !photo_to_file(file_arg[PHOTO_FILE_ARG], WRITE)))
                    send(CANTOPEN);
                sprintf(ps, "photo is on: file is '%s'\n", photo_file);
                pr();
            } on_error { print("error opening photo file\n"); }
    }

    if (!nullstr(file_arg[LOAD_FILE_ARG])) {
        run {
                nl();
                fp = open_file(file_arg[LOAD_FILE_ARG], READ);
                try_to_load(fp, file_arg[LOAD_FILE_ARG], FALSE, prep_it); /*verbose*/
            } on_error { print("error opening load or prep file\n"); }
    }

    if (!nullstr(file_arg[RUN_FILE_ARG])) {
        run {
                redirect_input(file_arg[RUN_FILE_ARG], TRUE); /*verbose*/
            } on_error { print("error opening run file\n"); }
    }

    command_loop();
    /* screen_end(); */
#ifdef THREAD
    threads_shutdown();
#endif
    exit_main();
}


bool mapm_save_on_exit(bool do_it_now) {
    if (!do_it_now) return (auto_save && data_loaded());
    if (auto_save && data_loaded()) do_save_data(raw.filename, FALSE);
    return (TRUE); /* => OK to exit now */
}


#define NO_DATA_ERR \
"No data have been loaded.\nUse the 'load data' command first."
#define IS_DATA_ERR \
"%sYou must use the '%s' command before loading data."
#define IS_DATA_ERR1 \
"Data have already been loaded.\n"
#define WRONG_DATA_ERR \
"Data of the wrong type are loaded.\nThe '%s' command only works with %s data."

#define NO_SEQ_ERR \
"The current sequence is empty.\nUse the 'sequence' command first."
#define SHORT_SEQ_ERR \
"%sThe '%s' command requires a sequence containing at least %d loci."
#define SHORT_SEQ_ERR1 \
"The current sequence is too short.\n"

#define SEQ_PERM_ERR \
"The current sequence specifies only one order.\n\
The '%s' command requires a sequence which specifies more than one order\n\
of loci. %s"
#define SEQ_PERM_WARN \
"warning: The current sequence specifies more than one order.\n\
The '%s' command ignores alternative orders.\n%s"
#define SEQ_FIX_WARN \
"warning: The current sequence specifies fixed distances.\n\
The '%s' command ignores order and distance information.\n%s"
#define SEQ_EXP_EMPTY "After expanding names, the current sequence is empty."
#define SEQ_HELP      "Type 'help sequence' for details."

void mapm_ready(int data_type, /* CEPH, F2, NO_DATA, ANY_DATA, or MAYBE_DATA */ int min_seq_loci, /* 0 indicates no seq is needed */
                bool permable_seq, /* TRUE, FALSE or MAYBE, ignored if min_seq_loci==0 */ int *seq_loci) {
    int loci;

    if (data_type == NO_DATA && data_loaded()) {
        sprintf(ps, IS_DATA_ERR, IS_DATA_ERR1, com);
        error(ps);
    }
    if (data_type != MAYBE_DATA && !data_loaded()) error(NO_DATA_ERR);
    /* CHECK THESE! */
    if ((data_type == CEPH_DATA && raw.data_type != CEPH) ||
        (data_type == F2_DATA && raw.data_type != F2)) {
        sprintf(ps, WRONG_DATA_ERR, com,
                data_type == CEPH ? "CEPH-type data" : "F2-type data");
        error(ps);
    }

    seq = NULL;
    loci = 0;
    run check_current_seq(&loci);
    except_when(BADSEQ) {
        if (min_seq_loci != 0) { /* n>0 or MAYBE_SEQ (-1) means need seq */
            print_badseq();
            abort_command();
        }
    }
    if (seq_loci != NULL) *seq_loci = loci;

    if (min_seq_loci == 0) {
        seq = NULL;
        if (seq_loci != NULL) send(CRASH);
        return;
    }

    if (min_seq_loci > 0) { /* test min num loci in seq */
        if (loci == 0) {
            if (!nullstr(seq_string) && !xstreq(seq_string, "none"))
                error(SEQ_EXP_EMPTY);
            else error(NO_SEQ_ERR);
        }
        if (loci < min_seq_loci) {
            sprintf(ps, SHORT_SEQ_ERR, SHORT_SEQ_ERR1, com, min_seq_loci);
            error(ps);
        }
    }

    /* test permability and crunchiness */
    if (seq != NULL) get_list_of_all_loci(seq, seq_locus, &loci, MAX_SEQ_LOCI);
    switch (permable_seq) {
        case PERM_SEQ:   /* ordered seq - need to use for_all_orders() */
            if (seq != NULL) {
                if (!permutable(seq)) {
                    sprintf(ps, SEQ_PERM_ERR, com, SEQ_HELP);
                    error(ps);
                }
                crunch_locus_list(seq_locus, &loci, ORDER_ERRORS, TRUE, TRUE);
            }
            break;
        case MAYBE_PERM: /* ordered seq - for_all_orders() or get_one_order() */
            /* no permability test, but we need a crunch test */
            if (seq != NULL)
                crunch_locus_list(seq_locus, &loci, ORDER_ERRORS, TRUE, TRUE);
            break;
        case ONE_ORDER:  /* ordered seq - get_one_order() */
            if (seq != NULL) {
                if (permutable(seq)) {
                    sprintf(ps, SEQ_PERM_WARN, com, SEQ_HELP);
                    pr();
                    nl();
                    nl();
                }
                crunch_locus_list(seq_locus, &loci, ORDER_ERRORS, TRUE, TRUE);
            }
            break;
        case LIST_SEQ:   /* not ordered - alloc_ or get_list_of_all_loci() */
            if (seq != NULL) {
                if (permutable(seq)) {
                    sprintf(ps, SEQ_PERM_WARN, com, SEQ_HELP);
                    pr();
                    nl();
                    nl();
                }
                if (has_fixed_dists(seq)) {
                    sprintf(ps, SEQ_FIX_WARN, com, SEQ_HELP);
                    pr();
                    nl();
                    nl();
                }
                crunch_locus_list(seq_locus, &loci, CRUNCH_WARNINGS, TRUE, TRUE);
            }
            break;
        case UNCRUNCHED_LIST: /* not ordered, for assign and haplo commands */
            if (seq != NULL) {
                if (permutable(seq)) {
                    sprintf(ps, SEQ_PERM_WARN, com, SEQ_HELP);
                    pr();
                    nl();
                    nl();
                }
                if (has_fixed_dists(seq)) {
                    sprintf(ps, SEQ_FIX_WARN, com, SEQ_HELP);
                    pr();
                    nl();
                    nl();
                }
            }
            break;
        default:
            send(CRASH);
    }
}


#define WARN_CHROMOSOME \
  "%s: %s markers not assigned to current chromosome"
#define WARN_DUPLICATES \
  "%s: marker(s) listed %smultiple times... %s"
#define WARN_CONVERT \
  "%s: haplotype group(s) %s%sindicated by first locus"
#define WARN_HAPLO_DUPS \
  "%s: haplotype group(s) listed %smultiple times... %s"

bool crunch_locus_list(int *locus, int *num_loci, bool verbose, /* ORDER_ERRORS, or CRUNCH_WARNINGS, or FALSE (silent) */ bool check_assignments,
                       bool in_sequence /* adjusts output: TRUE, FALSE, or MAYBE */) {
    int i, n;
    bool haplos_converted, haplo_dups, other_dups, wrong_chrom;
    bool abort_on_error;
    /* we use the seq_locus and use_locus globals from sequence.c */

    haplos_converted = haplo_dups = other_dups = wrong_chrom = FALSE;
    abort_on_error = (verbose == ORDER_ERRORS);

    for (i = 0; i < raw.num_markers; i++) use_locus[i] = FALSE;
    for (i = 0; i < *num_loci; i++) seq_locus[i] = locus[i];

    if (use_haplotypes)
        for (i = 0; i < *num_loci; i++)
            if (seq_locus[i] != NO_LOCUS && haplotyped(seq_locus[i])) {
                if (haplo_first[seq_locus[i]] == seq_locus[i]) continue;
                if (!use_locus[haplo_first[seq_locus[i]]]) {
                    seq_locus[i] = haplo_first[seq_locus[i]];
                    use_locus[seq_locus[i]] = TRUE; /* now using the haplo_first */
                    haplos_converted = TRUE;
                } else {
                    seq_locus[i] = NO_LOCUS;
                    haplo_dups = TRUE;
                }
            }

    for (i = 0; i < *num_loci; i++) {
        if (!use_locus[seq_locus[i]]) use_locus[seq_locus[i]] = TRUE;
        else {
            seq_locus[i] = NO_LOCUS;
            other_dups = TRUE;
        }
    }

    for (i = 0, n = 0; i < *num_loci; i++)
        if (seq_locus[i] != NO_LOCUS) locus[n++] = seq_locus[i];
    *num_loci = n;

    if (current_chrom != NO_CHROM && check_assignments)
        for (i = 0, n = 0; i < *num_loci; i++)
            if (assigned(seq_locus[i]) && !assigned_to(seq_locus[i], current_chrom))
                wrong_chrom = TRUE;

    if (verbose) { /* TRUE or MAYBE */
        if (wrong_chrom) {
            sprintf(ps, WARN_CHROMOSOME, "error",
                    (in_sequence ? "sequence contains" : "can't use"));
            pr();
            nl();
        }
        if (haplos_converted) {
            sprintf(ps, WARN_CONVERT, "warning", (in_sequence ? "in sequence " : ""),
                    (in_sequence == MAYBE ? "will be " : ""));
            pr();
            nl();
        }
        if (haplo_dups) {
            sprintf(ps, WARN_HAPLO_DUPS, (abort_on_error ? "error" : "warning"),
                    (in_sequence ? "in sequence " : ""),
                    (in_sequence == MAYBE ? "" : "extras deleted"));
            pr();
            nl();
        }
        if (other_dups) {
            sprintf(ps, WARN_DUPLICATES, (abort_on_error ? "error" : "warning"),
                    (in_sequence ? "in sequence " : ""),
                    (in_sequence == MAYBE ? "" : "extras deleted"));
            pr();
            nl();
        }
        if ((abort_on_error && (haplos_converted || haplo_dups)) ||
            wrong_chrom)
            abort_command();
    }
    return (other_dups || haplos_converted || haplo_dups);
}


#define CHROM_NOT_EXISTS "there is no chromosome named '%s'"
#define CHROM_NOT_SET    "no chromosome is selected"
#define CHROM_NOT_ANY    "you must select a chromosome ('any' is not allowed)"

int get_chrom_arg(bool allow_no_chrom) {
    char name[TOKLEN + 1];
    int chrom;

    get_arg(stoken, "", name);
    if (nullstr(name)) {
        if (current_chrom == NO_CHROM) {
            if (!allow_no_chrom) error(CHROM_NOT_SET);
            else chrom = NO_CHROM;
        } else chrom = current_chrom;
    } else if (streq(name, "any") || streq(name, "all")) {
        if (!allow_no_chrom) error(CHROM_NOT_ANY);
        else chrom = NO_CHROM;
    } else if (!isa_chrom(name, &chrom)) {
        sprintf(ps, CHROM_NOT_EXISTS, name);
        error(ps);
    }
    return (chrom);
}


bool input_dist(real *dist) {
    if (*dist < 0.0) return (FALSE);
    else if (*dist <= 0.5) return (TRUE);
    else if (*dist >= 999.0) return (FALSE);

    *dist = ((*mapfunction->dist_to_rec)(*dist / 100.0));
    return (TRUE);
}


void setup_commands(void) {
    two_pt_touched = FALSE;
    three_pt_touched = FALSE;

    /*  command-name	             Abbrev	function	   	type  
	1234567890123456789012345    123	12345678901234...  	123
	------------------------- sp ---  tab   --------------	  tab   --- */

    /* npt commands */
    mc("map", "m", make_map, CMD);
    mc("draw map", "", draw_map, CMD);
    mc("compare", "c", compare, CMD);
    mc("try", "", try, CMD);
    mc("genotypes", "", genotypes, CMD);
    mc("ripple", "", ripple, CMD);
    mc("order", "", order_maker, CMD);
    mc("build", "", greedy, CMD);

    /* state.c */
    mc("print names", "", set_print_names, OPT);
    mc("tolerance", "", set_tolerance, PAR);
    mc("units", "", set_units, PAR);
    mc("centimorgan function", "", set_cm_func, PAR);
    mc("auto save data", "", set_autosave, OPT);
    mc("more mode", "", set_more_mode, OPT);

    mc("default linkage criteria", "", set_default_linkage, PAR);
    mc("use three point", "", set_use_3pt, OPT);
    mc("triple linkage criteria", "", set_3pt_linkage, PAR);
    mc("triple exclusion criteria", "", set_3pt_threshold, PAR);
    mc("triple error detection", "", set_3pt_errors, OPT);
    mc("multipoint criteria", "", set_npt_threshold, PAR);
    mc("informativeness criteria", "", set_inf_threshold, PAR);
    mc("print maps", "", set_print_all_maps, OPT);

    mc("error detection", "", set_use_error_rate, OPT);
    mc("error probability", "", set_error_rate, CMD);
    mc("error thresholds", "", set_error_lod_thresh, PAR);

    /* sequence commands now in sys_cmds.c */
    mc("sequence", "s", sequence, CMD);
    mc("expand sequence", "x", expand_sequence, CMD);
    mc("delete", "d", new_delete, CMD);
    mc("insert", "i", new_insert, CMD);
    mc("append", "a", new_append, CMD);
    mc("flip", "", flip, CMD);
    mc("history", "h", show_seq_history, CMD);
/*  mc("let",                        "l",	let,			CMD);*/
    mc("let", "l", let_expanding, CMD);
    mc("names", "n", names, CMD);
    mc("forget named sequence", "", forget, CMD);
    mc("edit sequence", "e", edit_sequence, CMD);
    mc("translate", "t", translate, CMD);

    /* general commands in sys_cmds.c */
    mc("prepare data", "pd", new_prepare, CMD);
    mc("load data", "ld", new_load_data, CMD);
    mc("save data", "", new_save_data, CMD);
    mc("age", "", set_age, CMD);
    mc("class", "", set_class, CMD);
    mc("make class", "", make_classes, CMD);
    mc("list loci", "ll", list_loci, CMD);

    /* shell commands in shell.h */
    mc("quit", "q", quit, CMD);
    mc("photo", "", do_photo, CMD);
    mc("help", "?", help, CMD);
    mc("about mapmaker", "", about, CMD);
    mc("system", "", system_command, CMD);
    mc("previous commands", "p", show_cmd_history, CMD);
    mc("review output", "", review_output, CMD);
    mc("change directory", "cd", cd_command, CMD);
    mc("run", "", run_from_file, CMD);
    mc("time", "", show_time, CMD);
    mc("remark", "", comment, CMD);
    mc("comment", "", comment, CMD);
    mc("wizard mode", "", set_wizard, CMD);

    /* 2pt commands in cmds_2.c */
    mc("two point", "", two_point, CMD);
    mc("big lods", "", biglods, CMD);
    mc("all lods", "", all_lods, CMD);
    mc("lod table", "", lodtable, CMD);
    mc("near", "", near_locus, CMD);
    mc("pairwise", "", pairwise, CMD);
    mc("group", "", group, CMD);
    mc("suggest subset", "", suggest_subset, CMD);
    mc("join haplotypes", "", haplotype, CMD);
    mc("list haplotypes", "lh", list_haplotypes, CMD);
    mc("restore haplotypes", "", unhaplotype, CMD);

    /* 3pt commands in cmds_2.c */
    mc("three point", "", three_point, CMD);
    mc("forget three point", "", forget_three_point, CMD);

    /* auto commands in auto.c */
    mc("links", "", near_chrom, CMD);
    mc("make chromosome", "", make_chromosome, CMD);
    mc("anchor", "", set_anchors, CMD);
    mc("assign", "", assign, CMD);
    mc("unassign", "", unassign, CMD);
    mc("attach", "", attach, CMD);
    mc("framework", "", set_framework, CMD);
    mc("place", "", place, CMD);
    mc("together", "", place_together, CMD);
    mc("list chromosomes", "lc", list_chroms, CMD);
    mc("list assignments", "la", list_assignments, CMD);
    mc("list status", "ls", list_mapping, CMD);
    mc("draw chromosome", "", draw_chromosome, CMD);
    mc("draw all chromosomes", "", draw_all_chromosomes, CMD);

    /* special commands */
    mc("import data", "", import, CMD);
    mc("export data", "", export, CMD);

/*  mc("fake maps",                  "",	set_fake_maps,		WIZ);*/


#ifdef PUNT_THESE_FOR_NOW
    /* wizard, obsolete, and undeleted CEPH commands */
    mc("likely",              likely,          	CMD);
    mc("revise data",         revisedat,       	CMD);
    mc("new load",            new_load_data,   	CMD);
    mc("new save",            new_save_data,   	CMD);
    mc("use hmm",                   set_use_hmm,                	CMD);

    mc("note",                       "",	make_note,       	CMD);
    mc("sex specific",              set_sex_specific,           CMD);
    mc("segregation distortion",    set_segregation_distortion, CMD);
    mc("print maps",                set_print_maps,             CMD);
    mc("inner tolerance",           set_inner_tolerance,        CMD);
    mc("startrecombs",              set_startrecombs,           CMD);
    mc("inner loop",                set_inner_loop,             CMD);
    mc("print problem size",        set_print_problem_size,     CMD);
    mc("max problem size",          set_max_problem_size,       CMD);
    mc("time stamping",             set_time_stamping,          CMD);
    mc("print dots",                set_print_dots,             CMD);
    mc("seg dist",            seg_dist,        CMD);
#endif
}
