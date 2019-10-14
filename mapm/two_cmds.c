/******************************************************************************

  #####  #    #   ####            ####   #    #  #####    ####            ####
    #    #    #  #    #          #    #  ##  ##  #    #  #               #    #
    #    #    #  #    #          #       # ## #  #    #   ####           #
    #    # ## #  #    #          #       #    #  #    #       #   ###    #
    #    ##  ##  #    #          #    #  #    #  #    #  #    #   ###    #    #
    #    #    #   ####  #######   ####   #    #  #####    ####    ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#include "mapm.h"

/* internal procedures */

void print_2pt_criteria(char *str, real lod, real theta);

void parse_2pt_criteria(char **argptr, real *lod, real *theta);

static void print_biglod_line(int i, int j, real theta, real lod, real dist);

bool print_biglod(int i, int j, real lodbound, real thetabound, bool sex, int chrom);

void print_lodtable(int *locus1, int *locus2, int num_loci1, int num_loci2, int how);

void print_lod_line(int loc, int *toploc, int topnum, int how, int topref, int downref);

#define PRINT_HALF_LODTABLE 1
#define PRINT_FULL_LODTABLE 2

void print_2pt_criteria(char *str, real lod, real theta) {
    sprintf(ps, "%s at min LOD %.2lf, max Distance %s", str, lod, rag(rf2str(theta)));
    pr();
}

void parse_2pt_criteria(char **argptr, real *lod, real *theta) {
    if (!rtoken(argptr, default_lod, lod) || !rrange(lod, 0.0, 1000.0))
        error("bad value for min LOD score");
    if (!rtoken(argptr, default_theta, theta) || !input_dist(theta))
        error("bad value for max distance");
}


command two_point(void) {
    int a, b, num_loci, count, total, *locus = NULL;

    mapm_ready(ANY_DATA, 2, UNCRUNCHED_LIST, NULL);
    nomore_args(0);

    run {
            count = total = 0;
            alloc_list_of_all_loci(seq, &locus, &num_loci);
            crunch_locus_list(locus, &num_loci, CRUNCH_WARNINGS, ANY_CHROMS, IN_SEQ);
            total = (num_loci * (num_loci - 1)) / 2;

            for_all_locus_pairs (locus, num_loci, a, b) {
                    keep_user_amused("pair", count++, total);
                    compute_two_pt(a, b);
                }
            print("two-point data are available.\n");

        } on_exit {
        unarray(locus, int);
        relay_messages;
    }
}


#define GROUP_DIVIDER "-------\n"

command group(void) {
    real lodbound, thetabound;
    int *loci = NULL, num_loci, *linkage_group = NULL, group_size, num_left;
    int *really_unlinked = NULL, num_unlinked, num_linkage_groups, i;

    mapm_ready(ANY_DATA, 2, UNCRUNCHED_LIST, NULL);
    parse_2pt_criteria(&args, &lodbound, &thetabound);
    nomore_args(2);

    run {
            alloc_list_of_all_loci(seq, &loci, &num_loci);
            crunch_locus_list(loci, &num_loci, CRUNCH_WARNINGS, ANY_CHROMS, IN_SEQ);
            array(linkage_group, num_loci, int);
            array(really_unlinked, num_loci, int);
            print_2pt_criteria("Linkage Groups", lodbound, thetabound);
            nl();

            num_groups = 0; /* The global - set to 0 just in case something dies */
            num_linkage_groups = num_unlinked = 0;
            for (i = 0; i < num_loci; i++) my_group[i] = NO_GROUP;
            num_left = num_loci;

            do {
                get_linkage_group(loci, &num_left, linkage_group, &group_size,
                                  lodbound, thetabound);
                if (group_size >= 2) {
                    inv_isort(linkage_group, group_size);
                    if (num_linkage_groups > 0) print(GROUP_DIVIDER); else nl();
                    sprintf(ps, "group%d= ", num_linkage_groups + 1);
                    pr();
                    for (i = 0; i < group_size; i++) {
                        print(rag(loc2str(linkage_group[i])));
                        print(" ");
                        my_group[linkage_group[i]] = num_linkage_groups;
                    }
                    nl();
                    num_linkage_groups++;
                } else {
                    /* it is already in NO_GROUP */
                    really_unlinked[num_unlinked++] = linkage_group[0];
                }
            } while (num_left > 0);

            if (num_unlinked > 0) {
                if (num_linkage_groups > 0) print(GROUP_DIVIDER); else nl();
                print("unlinked= ");
                for (i = 0; i < num_unlinked; i++) {
                    print(rag(loc2str(really_unlinked[i])));
                    print(" ");
                    my_group[really_unlinked[i]] = num_linkage_groups;
                }
                nl();
            }
            num_groups = num_linkage_groups + 1; /* last one is "unlinked" group */

        } on_exit {
        unarray(linkage_group, int);
        unarray(really_unlinked, int);
        unarray(loci, int);
        relay_messages;
    }
}


#define HAP_DEL_NOTE \
"note: Deleting halplotype groups containing %s markers!\n"
#define HAP_OFF_NOTE   "'join haplotypes' is off.\n"
#define HAP_ON_NOTE    "'join haplotypes' is on.\n"
#define HAP_NO_GROUPS  "No linkage groups found.\n"
#define HAP_NO_HAPLOS  "No appropriate haplotype groups found.\n"

command haplotype(void) {
    real lodbound, thetabound;
    int *loci = NULL, num_loci, num_left, *linkage_group = NULL, group_size;
    int *haplo_group = NULL, num_haplo, num_linkage_groups, num_haplo_groups;
    int *changed = NULL, num_changed, *obs1 = NULL, *obs2 = NULL, i;
    bool find_em, any;
    char first[TOKLEN + 1];

    mapm_ready(F2_DATA, 2, UNCRUNCHED_LIST, &num_loci);
    get_arg(stoken, "", first);
    if (streq(first, "on")) {
        find_em = TRUE;
        parse_2pt_criteria(&args, &lodbound, &thetabound);
        print_2pt_criteria("Haplotyping", lodbound, thetabound);
        print("\nHaplotype groups allow NO Obligate Recombinants.\n");
    } else if (streq(first, "off")) {
        find_em = FALSE;
    } else if (nullstr(args)) {
        for (any = FALSE, i = 0; i < raw.num_markers; i++)
            if (haplotyped(i)) {
                any = TRUE;
                break;
            }
        if (any) print(HAP_ON_NOTE); else print(HAP_OFF_NOTE);
        return;
    } else usage_error(1);
    nomore_args(0);

    run {
            array(loci, num_loci, int);
            array(linkage_group, num_loci, int);
            array(haplo_group, num_loci, int);
            array(obs1, raw.data.f2.num_indivs, int);
            array(obs2, raw.data.f2.num_indivs, int);
            array(changed, 2 * raw.num_markers, int);
            get_list_of_all_loci(seq, loci, &num_loci, num_loci); /* no crunching */

            num_linkage_groups = num_haplo_groups = num_changed = 0;
            if (!find_em) {
                num_loci = raw.num_markers;
                for (i = 0; i < raw.num_markers; i++) loci[i] = i;
            }
            if (delete_haplo_groups(loci, num_loci, changed, &num_changed)) {
                sprintf(ps, HAP_DEL_NOTE, (num_loci == raw.num_markers ? "all" : "these"));
                pr();
            }
            if (find_em) {
                nl();
                num_left = num_loci;
                do {
                    get_linkage_group(loci, &num_left, linkage_group, &group_size,
                                      lodbound, thetabound);
                    if (group_size >= 2) {
                        sort_loci(linkage_group, group_size);
                        do {
                            find_haplo_group(linkage_group, &group_size, haplo_group,
                                             &num_haplo, obs1, obs2);
                            if (num_haplo >= 2) {
                                setup_haplo_group(haplo_group, num_haplo);
                                sprintf(ps, "Haplotype Group %2d: %s= [",
                                        ++num_haplo_groups,
                                        locname(haplo_first[haplo_group[0]], TRUE));
                                pr();
                                for (i = 0; i < num_haplo; i++) {
                                    if (i != 0) print(" ");
                                    print(rag(locname(haplo_group[i], FALSE)));
                                    changed[num_changed++] = haplo_group[i];
                                }
                                print("]\n");
                            } /* num_haplo>2 */
                        } while (group_size >= 1);
                        num_linkage_groups++;
                    } /* group_size>2 */
                } while (num_left > 0);
                nl();
                if (num_linkage_groups == 0) print(HAP_NO_GROUPS);
                else if (num_haplo_groups == 0) print(HAP_NO_HAPLOS);
            }
            if (num_changed > 0) {
                if (!find_em) nl();
                sort_loci(changed, num_changed);
                bash_mapping_data(changed, num_changed);
                bash_order_info(changed, num_changed);
            }
            for (any = FALSE, i = 0; i < raw.num_markers; i++)
                if (haplotyped(i)) {
                    any = TRUE;
                    break;
                }
            if (any) print(HAP_ON_NOTE); else print(HAP_OFF_NOTE);

        } on_exit {
        unarray(linkage_group, int);
        unarray(haplo_group, int);
        unarray(loci, int);
        unarray(changed, int);
        unarray(obs1, int);
        unarray(obs2, int);
        relay_messages;
    }
}


command unhaplotype(void) {
    int *loci = NULL, num_loci, *changed = NULL, num_changed;

    mapm_ready(F2_DATA, 0, 0, NULL);
    run {
            if (nullstr(args)) usage_error(1);
            parse_locus_args(&loci, &num_loci); /* error if fails, is uncrunched */
            array(changed, 2 * raw.num_markers, int);

            if (delete_haplo_groups(loci, num_loci, changed, &num_changed)) {
                sprintf(ps, HAP_DEL_NOTE, "these");
                pr();
            } else print("No haplotype groups to delete\n");
            if (num_changed > 0) {
                sort_loci(changed, num_changed);
                bash_mapping_data(changed, num_changed);
                bash_order_info(changed, num_changed);
            }
        } on_exit {
        unarray(loci, int);
        unarray(changed, int);
        relay_messages;
    }
}


command list_haplotypes(void) {
    int num_loci, *locus = NULL;

    mapm_ready(ANY_DATA, MAYBE_SEQ, UNCRUNCHED_LIST, NULL);
    run {
            if (!nullstr(args)) {
                parse_locus_args(&locus, &num_loci); /* error if fails */
                if (num_loci == 0) error("no loci in arguments\n");
            } else {
                if (!alloc_list_of_all_loci(seq, &locus, &num_loci))
                    error(NEED_SEQ_OR_ARGS);
            }
            nl();
            crunch_locus_list(locus, &num_loci, SILENTLY, ANY_CHROMS, 0);
            hold(more_mode) print_haplo_summary(locus, num_loci);

        } on_exit {
        unarray(locus, int);
        relay_messages;
    }
}


command list_loci(void) {
    int source, num_loci, *locus = NULL;

    mapm_ready(ANY_DATA, MAYBE_SEQ, UNCRUNCHED_LIST, NULL);
    run {
            if (!nullstr(args)) {
                parse_locus_args(&locus, &num_loci); /* error if fails */
                if (num_loci == 0) error("no loci in arguments\n");
                source = IN_ARGS;
            } else {
                if (!alloc_list_of_all_loci(seq, &locus, &num_loci))
                    error(NEED_SEQ_OR_ARGS);
                source = IN_SEQ;
            }
            crunch_locus_list(locus, &num_loci, CRUNCH_WARNINGS, ANY_CHROMS, source);
            nl();
            hold(more_mode) print_locus_summary(locus, num_loci, use_haplotypes);

        } on_exit {
        unarray(locus, int);
        relay_messages;
    }
}


#define NOSEX_BIGLODS "\n Marker-1   Marker-2   Theta    LOD     cM\n"
#define SEX_BIGLODS \
  "\n Marker-1  Marker-2    Theta(M) Theta(F)  LOD        cM(M)    cM(F)\n"

command biglods(void) {
    int a, b, n, num_loci, *loci = NULL;
    real lodbound, thetabound;

    mapm_ready(ANY_DATA, 2, UNCRUNCHED_LIST, NULL);
    parse_2pt_criteria(&args, &lodbound, &thetabound);
    nomore_args(2);

    run {
            alloc_list_of_all_loci(seq, &loci, &num_loci);
            crunch_locus_list(loci, &num_loci, CRUNCH_WARNINGS, ANY_CHROMS, IN_SEQ);
            print_2pt_criteria("Linked Marker Pairs", lodbound, thetabound);
            nl();
            if (!sex_specific) print(NOSEX_BIGLODS); else print(SEX_BIGLODS);
            n = 0;
            hold(more_mode) for_all_locus_pairs(loci, num_loci, a, b) if (a != b &&
                                                                          print_biglod(a, b, lodbound, thetabound, sex_specific, NO_CHROM))
                            n++;
            if (n == 0) print("no linked markers in sequence\n");
        } on_exit {
        unarray(loci, int);
        relay_messages;
    }
}


struct two_pt_result {
    int i, j;
    real thetam, thetaf, lod, distm, distf;
};

command all_lods(void) {
    int a, b, n, sort = -1, nr, num_loci, *loci = NULL;
    struct two_pt_result r, *results = NULL;
    char sort_arg[TOKLEN + 1];

    mapm_ready(ANY_DATA, 2, UNCRUNCHED_LIST, NULL);
    if (!nullstr(args)) {
        get_arg(stoken, "", sort_arg);
        if (strcasecmp(sort_arg, "theta") == 0) sort = 0;
        else if (strcasecmp(sort_arg, "lod") == 0) sort = 1;
        else if (strcasecmp(sort_arg, "distance") == 0) sort = 2;
        else usage_error(1);
    }
    nomore_args(1);

    run {
            alloc_list_of_all_loci(seq, &loci, &num_loci);
            crunch_locus_list(loci, &num_loci, CRUNCH_WARNINGS, ANY_CHROMS, IN_SEQ);
            sprintf(ps, "All Marker Pairs\n"
                        "(Marker pairs with both Theta = %.2f and LOD = %.2f are unlinked)\n",
                    UNLINKED_THETA, UNLINKED_LOD);
            pr();

            nr = num_loci * (num_loci - 1) / 2;
            array(results, nr, struct two_pt_result);

            print(NOSEX_BIGLODS);
            n = 0;

            for_all_locus_pairs(loci, num_loci, r.i, r.j) if (r.i != r.j) {
                        get_two_pt(r.i, r.j, &r.lod, &r.thetam);
                        r.distm = cm_dist(r.thetam);
                        results[n++] = r;
                    }

            if (sort != -1) {
                for (a = 0; a < nr - 1; a++) {
                    for (b = a + 1; b < nr; b++) {
                        if ((sort == 0 && ((results[a].i > results[b].i) ||
                                           (results[a].i == results[b].i &&
                                            results[a].thetam > results[b].thetam) ||
                                           (results[a].i == results[b].i &&
                                            results[a].thetam == results[b].thetam &&
                                            results[a].i > results[b].j))) ||
                            (sort == 1 && ((results[a].i > results[b].i) ||
                                           (results[a].i == results[b].i &&
                                            results[a].lod < results[b].lod) ||
                                           (results[a].i == results[b].i &&
                                            results[a].lod == results[b].lod &&
                                            results[a].i > results[b].j))) ||
                            (sort == 2 && ((results[a].i > results[b].i) ||
                                           (results[a].i == results[b].i &&
                                            results[a].distm > results[b].distm) ||
                                           (results[a].i == results[b].i &&
                                            results[a].distm == results[b].distm &&
                                            results[a].i > results[b].j)))) {
                            r = results[a];
                            results[a] = results[b];
                            results[b] = r;
                        }
                    }
                }
            }

            hold(more_mode) for (a = 0; a < nr; a++) {
                    if (a > 0 && results[a].i != results[a - 1].i) nl();
                    print_biglod_line(results[a].i, results[a].j,
                                      results[a].thetam,
                                      results[a].lod,
                                      results[a].distm);
                    nl();
                }

            if (n == 0) print("no markers in sequence\n");
        } on_exit {
        unarray(loci, int);
        unarray(results, two_pt_result);
        relay_messages;
    }
}


command near_locus(void) {
    int a, b, i, j, n, num_loci, *locus = NULL, *trial_locus = NULL, num_trials;
    real lodbound, thetabound;
    char *rest;

    mapm_ready(ANY_DATA, 2, UNCRUNCHED_LIST, NULL);
    if (split_arglist(&rest, ':') && !nullstr(rest)) {
        parse_2pt_criteria(&rest, &lodbound, &thetabound);
        if (!nullstr(rest)) usage_error(0);
    } else {
        lodbound = default_lod;
        thetabound = default_theta;
    }
    if (nullstr(args)) usage_error(0);
    parse_locus_args(&trial_locus, &num_trials);
    crunch_locus_list(trial_locus, &num_trials, CRUNCH_WARNINGS, ANY_CHROMS,
                      IN_ARGS);

    run {
            alloc_list_of_all_loci(seq, &locus, &num_loci);
            crunch_locus_list(locus, &num_loci, CRUNCH_WARNINGS, ANY_CHROMS, IN_SEQ);
            print_2pt_criteria("Linked Markers in Sequence", lodbound, thetabound);
            nl();
            if (!sex_specific) print(NOSEX_BIGLODS); else print(SEX_BIGLODS);

            hold(more_mode) for (j = 0; j < num_trials; j++) {
                    a = trial_locus[j];
                    n = 0;
                    for (i = 0; i < num_loci; i++) {
                        b = locus[i];
                        if (a == b) continue;
                        if (print_biglod(a, b, lodbound, thetabound, sex_specific,
                                         assignment_chrom(locus[i])))
                            n++;
                    }
                    if (n == 0) {
                        sprintf(ps, " %-10s no linked markers\n", loc2str(a));
                        pr();
                    }
                }
        } on_exit {
        unarray(trial_locus, int);
        unarray(locus, int);
        relay_messages;
    }
}


command near_chrom(void) {
    /* should move to auto_cmd.c? */
    int a, b, i, n, num_loci, *locus = NULL, chrom;
    real lodbound, thetabound;
    char str[MAXLINE + 1];

    mapm_ready(ANY_DATA, 2, LIST_SEQ, &num_loci);
    chrom = get_chrom_arg(TRUE);
    parse_2pt_criteria(&args, &lodbound, &thetabound);
    nomore_args(2);

    run {
            alloc_list_of_all_loci(seq, &locus, &num_loci);
            crunch_locus_list(locus, &num_loci, CRUNCH_WARNINGS, ANY_CHROMS, IN_SEQ);
            if (chrom == NO_CHROM) strcpy(str, "Linked Markers in Data Set");
            else sprintf(str, "Linked Markers on Chromosome %s", chrom2str(chrom));
            print_2pt_criteria(str, lodbound, thetabound);
            nl();
            if (!sex_specific) print(NOSEX_BIGLODS); else print(SEX_BIGLODS);

            hold(more_mode) for (i = 0; i < num_loci; i++) {
                    a = locus[i];
                    n = 0;
                    for (b = 0; b < raw.num_markers; b++)
                        if ((!use_haplotypes || !haplotype_subordinate(b)) &&
                            (chrom == NO_CHROM || assigned_to(b, chrom)))
                            if (a != b &&
                                print_biglod(a, b, lodbound, thetabound, sex_specific,
                                             assignment_chrom(b)))
                                n++;
                    if (n == 0) {
                        sprintf(ps, " %-10s no linked markers\n", loc2str(a));
                        pr();
                    }
                }
        } on_exit {
        unarray(locus, int);
        relay_messages;
    }
}


command lodtable(void) {
    int *loci = NULL, num_loci, type, source;
    char name[TOKLEN + 1];

    mapm_ready(ANY_DATA, MAYBE_SEQ, UNCRUNCHED_LIST, &num_loci);
    get_one_arg(stoken, "half", name);
    if (matches(name, "half")) type = PRINT_HALF_LODTABLE;
    else if (matches(name, "full")) type = PRINT_FULL_LODTABLE;
    else usage_error(1);

    run {
            if (!nullstr(args)) {
                parse_locus_args(&loci, &num_loci); /* error if fails */
                if (num_loci == 0) error("no loci in arguments\n");
                source = IN_ARGS;
            } else {
                if (!alloc_list_of_all_loci(seq, &loci, &num_loci))
                    error(NEED_SEQ_OR_ARGS);
                source = IN_SEQ;
            }
            crunch_locus_list(loci, &num_loci, CRUNCH_WARNINGS, ANY_CHROMS, source);
            hold(more_mode) print_lodtable(loci, loci, num_loci, num_loci, type);

        } on_exit {
        unarray(loci, int);
        relay_messages;
    }
}


command pairwise(void) {
    int num_trials, *try_me = NULL, num_loci, *loci = NULL;

    mapm_ready(ANY_DATA, 2, UNCRUNCHED_LIST, &num_loci);
    run {
            if (nullstr(args)) usage_error(0);
            parse_locus_args(&try_me, &num_trials);
            crunch_locus_list(try_me, &num_trials, CRUNCH_WARNINGS, ANY_CHROMS,
                              IN_ARGS);
            alloc_list_of_all_loci(seq, &loci, &num_loci);
            crunch_locus_list(loci, &num_loci, CRUNCH_WARNINGS, ANY_CHROMS, IN_SEQ);

            hold(more_mode) print_lodtable(try_me, loci, num_trials, num_loci,
                                           PRINT_FULL_LODTABLE);
        } on_exit {
        unarray(try_me, int);
        unarray(loci, int);
        relay_messages;
    }
}


void print_biglod_line(int i, int j, real theta, real lod, real dist) {
    sprintf(ps, " %-10s %-10s %s   %s  %s", loc2str(i), loc2str(j),
            rsd(5.2, theta), rsd(5.2, lod), rsd(6.2, dist));
    pr();
}

bool print_biglod(int i, int j, real lodbound, real thetabound, bool sex, int chrom) {
    real lod, theta, thetam, thetaf;

    if (!sex) {
        if (!get_two_pt(i, j, &lod, &theta)) return (FALSE);
        if (lod >= lodbound && theta <= thetabound) {
            sprintf(ps, " %-10s %-10s %s   %s  %s", loc2str(i), loc2str(j),
                    rsd(5.2, theta), rsd(5.2, lod), rsd(6.2, cm_dist(theta)));
            pr();
            if (chrom != NO_CHROM) {
                sprintf(ps, "   (%s)", chrom2str(chrom));
                pr();
            }
            nl();
            return (TRUE);
        } else return (FALSE);
    } else { /* sex */
        if (!get_sex_two_pt(i, j, &lod, &thetam, &thetaf)) return (FALSE);
        if (lod >= lodbound && (thetam <= thetabound || thetaf <= thetabound)) {
            sprintf(ps, " %-10s %-10s  %s    %s   %s    %s    %s",
                    loc2str(i), loc2str(j),
                    rsd(5.2, thetam), rsd(5.2, thetaf), rsd(5.2, lod),
                    rsd(6.2, cm_dist(thetam)), rsd(6.2, cm_dist(thetaf)));
            pr();
            if (chrom != NO_CHROM) {
                sprintf(ps, "   (%s)", chrom2str(chrom));
                pr();
            }
            nl();
            return (TRUE);
        } else return (FALSE);
    }
}


void print_lodtable(int *locus1, int *locus2, int num_loci1, int num_loci2, int how) {
    int across, down, num_across, num_done, across_to_print;
    char *line1, *line2;

    print("Bottom number is LOD score, top number is ");
    if (units == CENTIMORGANS) print("centimorgan distance:\n");
    else print("recombination fraction:\n");

    num_done = 0;
    do {
        num_across = (num_loci1 > 10 ? 10 : num_loci1);
        across_to_print = num_loci1 - (how == PRINT_HALF_LODTABLE ? 1 : 0);

        nl();
        if (!print_names) {
            print("      ");
            for (across = 0; across < num_across; across++) {
                if (across != across_to_print) {
                    sprintf(ps, "%s ", loc2str(locus1[across]));
                    print(ps);
                }
            }
        } else {
            /* hack in names to look nice but also fit as before */
            line1 = get_temp_string();
            line2 = get_temp_string();
            sprintf(line1, "         ");
            sprintf(line2, "            ");
            for (across = 0; across < num_across; across++) {
                if (across != across_to_print) {
                    if (across % 2 == 0) {
                        sprintf(ps, "%s   ", loc2str(locus1[across]));
                        strcat(line1, ps);
                        strcat(line2, "   ");
                    } else {
                        sprintf(ps, "%s", loc2str(locus1[across]));
                        strcat(line2, ps);
                    }
                }
            }
            print(line1);
            nl();
            print(line2);
        }

        nl();
        for (down = 0; down < num_loci2; down++)
            print_lod_line(locus2[down], locus1, num_across, how, num_done, down);
        nl();

        num_loci1 -= num_across;
        num_done += num_across;
        if (num_loci1 > 0) locus1 += num_across;
    } while (num_loci1 > 0);
}


void print_lod_line(int loc, int *toploc, int topnum, int how, int topref, int downref) {
    int across, stop_at;
    real lod, theta, thetam, thetaf;

    if (how == PRINT_HALF_LODTABLE) {
        if (topref >= downref) return;
        stop_at = ((downref - topref) > topnum ? topnum : (downref - topref));
    } else stop_at = topnum;

    if (sex_specific) {
        if (print_names) print("\n        ");
        else print("\n    ");
        for (across = 0; across < stop_at; across++) {
            if (sex_specific) {
                if (get_sex_two_pt(loc, toploc[across], &lod, &thetam, &thetaf))
                    sprintf(ps, " %s", rf2str(thetam));
                else
                    sprintf(ps, "      ");
                print(ps);
            } else {
                /* print for UNLINKED_TWO_PT Similarly, test all return values
                   from get_two_pt and get_sex_two_pt below */
            }
        }
    }
    sprintf(ps, "\n%s", loc2str(loc));
    print(ps);
    for (across = 0; across < stop_at; across++) {
        if (sex_specific) {
            if (get_sex_two_pt(loc, toploc[across], &lod, &thetam, &thetaf))
                sprintf(ps, "%s ", rf2str(thetaf));
            else
                sprintf(ps, "  -   ");
        } else {
            if (get_two_pt(loc, toploc[across], &lod, &theta))
                sprintf(ps, "%s ", rf2str(theta));
            else
                sprintf(ps, "  -   ");
        }
        print(ps);
    }
    if (print_names) print("\n        ");
    else print("\n    ");
    for (across = 0; across < stop_at; across++) {
        if (sex_specific) {
            if (get_sex_two_pt(loc, toploc[across], &lod, &thetam, &thetaf))
                sprintf(ps, " %s", rsd(5.2, lod));
            else
                sprintf(ps, "      ");
        } else {
            if (get_two_pt(loc, toploc[across], &lod, &theta))
                sprintf(ps, " %s", rsd(5.2, lod));
            else
                sprintf(ps, "      ");
        }
        print(ps);
    }
    nl();
}


#define THREE_NO_GRPS  "no linkage groups of 3 or more loci found\n"
#define THREE_NO_TRIPS "no linked triplets found in %d linkage group%s\n"
#define THREE_DO_THIS  "%d linked triplet%s in %d linkage group%s\n"
#define THREE_HEAD  "Triplet criteria: LOD %.2lf, Max-Dist %s, #Linkages %d\n"
#define TRIPLET_SEX FALSE

command three_point(void) {
    int num_loci, *loci = NULL, num_trips, num_groups, num_links;
    int *linkage_group = NULL, group_size, num_unlinked, three_locus[3], foo;
    SEQ_NODE *three_seq;
    MAP *map = NULL;
    real three_like[3], lodbound2, thetabound2;
    real lodbound3, thetabound3;

    mapm_ready(ANY_DATA, 3, LIST_SEQ, &num_loci);
    nomore_args(0);

    /* these are the defaults */
    lodbound2 = default_lod;
    thetabound2 = default_theta;
    lodbound3 = triplet_lod;
    thetabound3 = triplet_theta;
    num_links = triplet_num_links;

    print_2pt_criteria("Linkage Groups", lodbound2, thetabound2);
    nl();
    sprintf(ps, THREE_HEAD, lodbound3, rag(rf2str(thetabound3)), num_links);
    pr();
    if (triplet_error_rate == 0.0) print("'triple error detection' is off.\n");
    else if (triplet_error_rate == LOCUS_ERROR_RATE)
        print("'triple error detection' is on.\n");
    else {
        print("'triple error detection' in three-point analysis is on.\n");
        sprintf(ps, "'error probability' for all loci is fixed at %.2lf%%.\n",
                triplet_error_rate * 100.0);
    }
    /* add: args! (what should they be?) pre-allocating 3pt, use_3pt=off msg,
       smart recalc, keep user amused, print thresholds in banner  */


    run {
            alloc_list_of_all_loci(seq, &loci, &num_loci); /* may delete haps! */
            array(linkage_group, num_loci, int);
            map = allocate_map(3);

            print("counting...");
            flush();
            num_trips = num_groups = 0;
            num_unlinked = num_loci;
            do {
                get_linkage_group(loci, &num_unlinked, linkage_group, &group_size,
                                  lodbound2, thetabound2);
                if (group_size >= 3) {
                    num_groups++;
                    /* compute acceptable three-point orders in the group */
                    sort_loci(linkage_group, group_size);
                    for_all_3pt_seqs (linkage_group, group_size, three_seq) {
                        get_list_of_all_loci(three_seq, three_locus, &foo, 3);
                        if (three_linked(three_locus, lodbound3, thetabound3,
                                         num_links, TRIPLET_SEX))
                            num_trips++;
                    }
                }
            } while (num_unlinked > 0);

            if (num_groups == 0) {
                sprintf(ps, THREE_NO_GRPS);
                pr();
                abort_command();
            } else if (num_trips == 0) {
                sprintf(ps, THREE_NO_TRIPS, num_groups, maybe_s(num_groups));
                pr();
                abort_command();
            } /* else */
            sprintf(ps, THREE_DO_THIS, num_trips, maybe_s(num_trips), num_groups,
                    maybe_s(num_groups));
            pr();

            if (!print_names) {
                print("\n                    log-likelihood differences\n");
                print(" count  markers         a-b-c  b-a-c  a-c-b\n");
                /*       12345:   1234 1234 1234 */
            } else {
                print("\n                                     log-likelihood differences\n");
                print(" count  markers                        a-b-c  b-a-c  a-c-b\n");
                /*       12345:   123456789 123456789 123456789  */
            }

            /* do it again, for real this time... */
            get_list_of_all_loci(seq, loci, &num_unlinked, num_loci);
            crunch_locus_list(loci, &num_loci, SILENTLY, ANY_CHROMS, IN_SEQ);
            num_trips = 0;
            do {
                get_linkage_group(loci, &num_unlinked, linkage_group, &group_size,
                                  lodbound2, thetabound2);
                if (group_size >= 3) {
                    /* compute acceptable three-point orders in the group */
                    inv_isort(linkage_group, group_size);
                    for_all_3pt_seqs (linkage_group, group_size, three_seq) {
                        get_list_of_all_loci(three_seq, three_locus, &foo, 3);
                        if (three_linked(three_locus, lodbound3, thetabound3,
                                         num_links, TRIPLET_SEX)) {
                            compute_3pt(three_seq, TRIPLET_SEX, triplet_error_rate,
                                        three_like, map);
                            if (!print_names) {
                                sprintf(ps, "%5d:  %d %d %d", ++num_trips,
                                        three_locus[0] + 1, three_locus[1] + 1,
                                        three_locus[2] + 1);
                                pad_to_len(ps, 23);
                                pr();
                            } else {
                                sprintf(ps, "%5d:  %s %s %s", ++num_trips,
                                        raw.locus_name[three_locus[0]],
                                        raw.locus_name[three_locus[1]],
                                        raw.locus_name[three_locus[2]]);
                                pad_to_len(ps, 38);
                                pr();
                            }
                            sprintf(ps, "%s %s %s\n", rsd(6.2, three_like[0]),
                                    rsd(6.2, three_like[1]), rsd(6.2, three_like[2]));
                            pr();
                        }
                    }
                }
            } while (num_unlinked > 0);

        } on_exit {
        free_map(map);
        unarray(loci, int);
        unarray(linkage_group, int);
        three_pt_touched = TRUE;
        relay_messages;
    }
}


command forget_three_point(void) {
    mapm_ready(ANY_DATA, 0, 0, NULL);
    nomore_args(0);
    bash_all_three_pt(raw.num_markers);
    print("All three-point data forgotten.\n");
}


#define INF_CRITERIA \
  "Informativeness: min #Individuals %d%s, min Distance %s\n"

command suggest_subset(void) {
    real lodbound, thetabound;
    int *loci = NULL, num_loci, *linkage_group = NULL, group_size, groups_done;
    int *subset = NULL, subset_size, num_unlinked, i, prev;

    mapm_ready(F2_DATA, 3, LIST_SEQ, &num_loci); /* F2 because #indivs */
    parse_2pt_criteria(&args, &lodbound, &thetabound);
    nomore_args(0);

    print_2pt_criteria("Informative Subgroups", lodbound, thetabound);
    nl();
    sprintf(ps, INF_CRITERIA, npt_min_indivs, (npt_codominant ? " (codominant)" : ""),
            rag(rf2str(npt_min_theta)));
    pr();

    run {
            alloc_list_of_all_loci(seq, &loci, &num_loci);
            array(linkage_group, num_loci, int);
            array(subset, num_loci, int);

            num_orders = 0; /* global - deletes all pre-existing orders */
            for (i = 0; i < raw.num_markers; i++) order_next[i] = NO_LOCUS;

            num_unlinked = num_loci;
            groups_done = 0;
            nl();
            do {
                get_linkage_group(loci, &num_unlinked, linkage_group, &group_size,
                                  lodbound, thetabound);
                if (group_size >= 3) {
                    if (groups_done != 0) print(GROUP_DIVIDER);
                    sort_loci(linkage_group, group_size);
                    groups_done++;
                    sprintf(ps, "Linkage group %d: ", groups_done);
                    pr();
                    for (i = 0; i < group_size; i++) {
                        print(rag(loc2str(linkage_group[i])));
                        if (i != group_size - 1) print(" ");
                    }
                    nl();

                    informative_subset(linkage_group, group_size,
                                       npt_min_indivs, npt_min_theta, npt_codominant,
                                       use_haplotypes, subset, &subset_size);

                    if (subset_size < 1)
                        print("No informative markers.\n");
                    else if (subset_size == group_size)
                        print("All markers are informative.\n");
                    else {
                        sprintf(ps, "%d Marker%s in informative subset:\n",
                                subset_size, maybe_s(subset_size));
                        pr();
                    }

                    sprintf(ps, "order%d= ", groups_done);
                    pr();
                    if (subset_size == 0) print("none");
                    for (i = 0; i < subset_size; i++) {
                        print(rag(loc2str(subset[i])));
                        if (i != subset_size - 1) print(" ");
                        if (i == 0) order_first[groups_done - 1] = subset[i];
                        else order_next[prev] = subset[i];
                        prev = subset[i]; /* lint warning is OK */
                    }
                    nl();
                    unorder_first[groups_done - 1] = NO_LOCUS;
                } /* if group size */
            } while (num_unlinked > 0);

            if (groups_done == 0) print("No Linkage Groups Found.\n");
            num_orders = groups_done;

        } on_exit {
        unarray(loci, int);
        unarray(linkage_group, int);
        unarray(subset, int);
        relay_messages;
    }
}
