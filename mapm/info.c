/******************************************************************************

    #    #    #  ######   ####            ####
    #    ##   #  #       #    #          #    #
    #    # #  #  #####   #    #          #
    #    #  # #  #       #    #   ###    #
    #    #   ##  #       #    #   ###    #    #
    #    #    #  #        ####    ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#include "mapm.h"

void expand_two_pt(int num_entries);

/* global vars which happen to be declared here */
char **note;
real *error_rate;
int *modified;
int *my_group, num_groups;
int *haplo_next, *haplo_first; /* [#loci] */
int *order_first, *unorder_first, *order_next, num_orders; /* [#loci] */
int *class;
char **class_name;

bool two_pt_touched;
bool three_pt_touched;

/* internal for 2pt */
TWO_PT_DATA ***two_pt_data = NULL;
int two_pt_allocated = 0, two_pt_used = 0, two_pt_max = 0, two_pt_collect = 0;
TWO_PT_DATA **two_pt_list = NULL;
TWO_PT_DATA *UNLINKED_TWO_PT, the_unlinked_two_pt;

/* internal for 3pt */
typedef struct three_entry {
    int locus2, locus3;
    float delta1, delta2, delta3;
    struct three_entry *next, *next_same2;
} TRIPLE_LIST;

typedef struct {
    struct three_entry *entry, *unused;
} THREE_STRUCT;

THREE_STRUCT *three_pt_data;

TRIPLE_LIST *allocate_three_entries(int num_to_alloc);

void deallocate_triple_list(TRIPLE_LIST *p);

int num_threes_allocated, num_threes_deallocated; /* for testing purposes */

bool replace_triple(int locus1, int locus2, int locus3, real d1, real d2, real d3);

void save_triple(int locus1, int locus2, int locus3, real d1, real d2, real d3);

bool same_2_loop(TRIPLE_LIST *p, int locus3, real *d1, real *d2, real *d3);

bool get_triple(int locus1, int locus2, int locus3, real *d1, real *d2, real *d3);

static void remove_triple_list(TRIPLE_LIST *p);

static void return_to_unused(TRIPLE_LIST *p);


#define UNFILLED 9999.0


/**************** Stuff to handle global 2pt data ****************/

void allocate_two_pt(int num_loci) {
    int i, j;

    two_pt_data = NULL;
    two_pt_list = NULL;
    two_pt_max = 0;

    array(two_pt_data, num_loci, TWO_PT_DATA**);
    for (i = 0; i < num_loci; i++) {
        array(two_pt_data[i], (i + 1), TWO_PT_DATA*);
        for (j = 0; j <= i; j++) two_pt_data[i][j] = NULL;
        two_pt_max += i + 1;
    }
    array(two_pt_list, two_pt_max, TWO_PT_DATA*);

    UNLINKED_TWO_PT = &the_unlinked_two_pt;
    UNLINKED_TWO_PT->used = FALSE;
    for (int i = 0; i < 3; i++) {
        UNLINKED_TWO_PT->lodscore[i] = UNLINKED_LOD;
        UNLINKED_TWO_PT->theta[i] = UNLINKED_THETA;
    }
//    UNLINKED_TWO_PT->lodscore[NOSEX]=   UNLINKED_LOD;
//    UNLINKED_TWO_PT->lodscore[SEXSPEC]= UNLINKED_LOD;
//    UNLINKED_TWO_PT->theta[NOSEX]=      UNLINKED_THETA;
//    UNLINKED_TWO_PT->theta[MALE]=       UNLINKED_THETA;
//    UNLINKED_TWO_PT->theta[FEMALE]=     UNLINKED_THETA;

    two_pt_used = two_pt_allocated = two_pt_collect = 0;
    two_pt_touched = FALSE;
}


void free_two_pt(int num_loci) {
    int i;
    if (two_pt_data == NULL) return;

    two_pt_touched = FALSE;
    for (i = 0; i < two_pt_allocated; i++) {unsingle(two_pt_list[i], TWO_PT_DATA); }
    unarray(two_pt_list, TWO_PT_DATA*);

    for (i = 0; i < num_loci; i++) {unarray(two_pt_data[i], TWO_PT_DATA*); }
    unarray(two_pt_data, TWO_PT_DATA**);
}


TWO_PT_DATA *get_next_two_pt_entry(int a, int b, bool *used_needs_incrementing) {
    TWO_PT_DATA *p = NULL;

    if (two_pt_data[a][b] != NULL && two_pt_data[a][b] != UNLINKED_TWO_PT) {
        p = two_pt_data[a][b];
        p->used = FALSE;
        for (int i = 0; i < 3; i++) {
            p->lodscore[i] = UNLINKED_LOD;
            p->theta[i] = UNLINKED_THETA;
        }
        two_pt_data[a][b] = NULL;
        two_pt_collect = 0;
        *used_needs_incrementing = FALSE;
        return (p);
    }

    if (two_pt_used == two_pt_allocated) {
        while (two_pt_collect < two_pt_used)
            if (!two_pt_list[two_pt_collect]->used) {
                *used_needs_incrementing = FALSE;
                {
                    p = two_pt_list[two_pt_collect];
                    p->used = FALSE;
                    for (int i = 0; i < 3; i++) {
                        p->lodscore[i] = UNLINKED_LOD;
                        p->theta[i] = UNLINKED_THETA;
                    }
                }
                return (two_pt_list[two_pt_collect]);
            } else two_pt_collect++;
        /* none free if we fall to here */
        expand_two_pt(EXPAND_DEFAULT);
    }
    if (two_pt_used >= two_pt_allocated || two_pt_list[two_pt_used] == NULL)
        send(CRASH);

    *used_needs_incrementing = TRUE;
    return (two_pt_list[two_pt_used]); /* note NOT incremented yet! */
}


void expand_two_pt(int num_entries) {
    int i;

    if (num_entries == EXPAND_DEFAULT)
        num_entries = INIT_2PT_SIZE(raw.num_markers);

    for (i = 0; i < num_entries; i++) {
        if (two_pt_allocated >= two_pt_max) break;
        single(two_pt_list[two_pt_allocated], TWO_PT_DATA);
        two_pt_list[two_pt_allocated]->used = FALSE;
        two_pt_allocated++;
    }

}


/* The following procedures are the only external accessors to the 
   TWO_PT_DATA structure. */

void compute_two_pt(int a, int b) {
    /* internal use only */
    TWO_PT_DATA *two_pt = NULL;
    int temp;
    bool inc;

    if (a < b) {
        temp = a;
        a = b;
        b = temp;
    }
    two_pt = get_next_two_pt_entry(a, b, &inc);
    quick_two_pt(a, b, two_pt/*,FALSE*/);
    if (raw.data_type == CEPH) {
        quick_two_pt(a, b, two_pt/*,TRUE*/); /* also w/sex */
    }

    if (two_pt->lodscore[NOSEX] <= UNLINKED_LOD &&
        two_pt->theta[NOSEX] >= UNLINKED_THETA &&
        (raw.data_type != CEPH ||
         (two_pt->lodscore[SEXSPEC] <= UNLINKED_LOD &&
          two_pt->theta[MALE] >= UNLINKED_THETA &&
          two_pt->lodscore[FEMALE] >= UNLINKED_THETA))) {
        two_pt_data[a][b] = UNLINKED_TWO_PT;
    } else {
        two_pt->used = TRUE;
        two_pt_data[a][b] = two_pt;
        if (inc) two_pt_used++;
    }
    two_pt_touched = TRUE;
}


bool get_two_pt(int a, int b, real *lod, real *theta) {
    /* TRUE if not UNLINKED_TWO_PT */
    int temp;

    if (a < b) {
        temp = a;
        a = b;
        b = temp;
    }
    if (two_pt_data[a][b] == NULL) compute_two_pt(a, b);

    if (lod != NULL) *lod = two_pt_data[a][b]->lodscore[NOSEX];
    if (theta != NULL) *theta = two_pt_data[a][b]->theta[NOSEX];
    return (two_pt_data[a][b] != UNLINKED_TWO_PT);
}


bool get_sex_two_pt(int a, int b, real *lod, real *thetam, real *thetaf) {
    /* TRUE if not UNLINKED_TWO_PT */
    int temp;

    if (a < b) {
        temp = a;
        a = b;
        b = temp;
    }
    if (two_pt_data[a][b] == NULL) compute_two_pt(a, b);

    if (lod != NULL) *lod = two_pt_data[a][b]->lodscore[SEXSPEC];
    if (thetam != NULL) *thetam = two_pt_data[a][b]->theta[MALE];
    if (thetaf != NULL) *thetaf = two_pt_data[a][b]->theta[FEMALE];
    return (two_pt_data[a][b] != UNLINKED_TWO_PT);
}


/**************** Stuff for Global 3pt data ****************/

void allocate_three_pt(int num_total) {
    int i;

    run {
            if (three_pt_data != NULL) free_three_pt(raw.num_markers);
            num_threes_allocated = 0;
            three_pt_touched = FALSE;
            single(three_pt_data, THREE_STRUCT);
            three_pt_data->entry = NULL;
            three_pt_data->unused = NULL;
            array(three_pt_data->entry, num_total, TRIPLE_LIST);
            for (i = 0; i < num_total; i++) {
                three_pt_data->entry[i].next = NULL;
                three_pt_data->entry[i].next_same2 = NULL;
                three_pt_data->entry[i].locus2 = NO_LOCUS;
                three_pt_data->entry[i].locus3 = NO_LOCUS;
                three_pt_data->entry[i].delta1 = UNFILLED;
                three_pt_data->entry[i].delta2 = UNFILLED;
                three_pt_data->entry[i].delta3 = UNFILLED;
            }
            three_pt_data->unused =
                    allocate_three_entries(INIT_3PT_SIZE(num_total));
        } when_aborting {
        free_three_pt(num_total);
        relay;
    }
}

TRIPLE_LIST *allocate_three_entries(int num_to_alloc) {
    TRIPLE_LIST *part;

    if (num_to_alloc == 0) return (NULL);
    single(part, TRIPLE_LIST);
    num_threes_allocated++;
    part->next = NULL;
    part->next_same2 = NULL;
    part->locus2 = NO_LOCUS;
    part->locus3 = NO_LOCUS;
    part->delta1 = UNFILLED;
    part->delta2 = UNFILLED;
    part->delta3 = UNFILLED;
    part->next = allocate_three_entries(num_to_alloc - 1);
    part->next_same2 = part->next;
    return (part);
}


void bash_all_three_pt(int num_total) {
    free_three_pt(num_total);
    allocate_three_pt(num_total);
}


void free_three_pt(int num_total) {
    int i;

    num_threes_deallocated = 0; /* for testing purposes */

    for (i = 0; i < num_total; i++) {
        if (three_pt_data->entry[i].next != NULL)
            deallocate_triple_list(three_pt_data->entry[i].next);
    }
    unarray(three_pt_data->entry, TRIPLE_LIST);
    deallocate_triple_list(three_pt_data->unused);
    unsingle(three_pt_data, THREE_STRUCT);
    three_pt_touched = FALSE;
    three_pt_data = NULL;

/*  sprintf(ps,"3pt stats - allocated: %d, deallocated: %d\n",
       num_threes_allocated, num_threes_deallocated); pr(); */
}

void deallocate_triple_list(TRIPLE_LIST *p) {
    if (p->next != NULL) deallocate_triple_list(p->next);
    unsingle(p, TRIPLE_LIST);
    num_threes_deallocated++;
}


/* insert_triple() arranges loci and likelihoods in numeric order (loci)
   and calls save_triple() to store the data in the global three_pt_data */

void insert_triple(int loc1, int loc2, int loc3, real d1, real d2, real d3) {
/* d1 - 123, d2 - 132, d3= 213 */
    int temploc;
    real tempreal;

    if (loc2 < loc1) {
        temploc = loc2;
        tempreal = d2;
        loc2 = loc1;
        d2 = d1;
        loc1 = temploc;
        d1 = tempreal;
    }
    if (loc3 < loc2) {
        temploc = loc3;
        tempreal = d3;
        loc3 = loc2;
        d3 = d2;
        loc2 = temploc;
        d2 = tempreal;
        if (loc2 < loc1) {
            temploc = loc2;
            tempreal = d2;
            loc2 = loc1;
            d2 = d1;
            loc1 = temploc;
            d1 = tempreal;
        }
    }
    save_triple(loc1, loc2, loc3, d1, d2, d3);
}


void save_triple(int locus1, int locus2, int locus3, real d1, real d2, real d3) {
    TRIPLE_LIST *t, *p;
    real r1, r2, r3;

    if (three_pt_data->entry[locus1].locus2 == NO_LOCUS) {
        three_pt_data->entry[locus1].locus2 = (short) locus2;
        three_pt_data->entry[locus1].locus3 = (short) locus3;
        three_pt_data->entry[locus1].delta1 = (float) d1;
        three_pt_data->entry[locus1].delta2 = (float) d2;
        three_pt_data->entry[locus1].delta3 = (float) d3;

    } else if (get_triple(locus1, locus2, locus3, &r1, &r2, &r3)) {
        /* over-write current entry */
        replace_triple(locus1, locus2, locus3, d1, d2, d3);

    } else {    /* must get blank from unused pile */
        if (three_pt_data->unused == NULL)
            three_pt_data->unused =
                    allocate_three_entries(INIT_3PT_SIZE(raw.num_markers));

        t = three_pt_data->unused;
        three_pt_data->unused = three_pt_data->unused->next;
        t->next = NULL;
        t->next_same2 = NULL;
        t->locus2 = (short) locus2;
        t->locus3 = (short) locus3;
        t->delta1 = (float) d1;
        t->delta2 = (float) d2;
        t->delta3 = (float) d3;

        if (three_pt_data->entry[locus1].next == NULL) {
            three_pt_data->entry[locus1].next = t; /* places t in linked list */
            if (three_pt_data->entry[locus1].locus2 == locus2)
                three_pt_data->entry[locus1].next_same2 = t;

        } else {
            for (p = three_pt_data->entry[locus1].next; p->next != NULL;
                 p = p->next) {}
            p->next = t;  /* places t on the end of the linked list ->next */

            /* The following stuff creates the ->next_same2 linked list
               which separately links all triple whose first two markers
               are the same. This aids in the search algorithm in get_triple(),
               speeding it up by about a factor of 20 */

            if (three_pt_data->entry[locus1].locus2 == locus2) {
                if (three_pt_data->entry[locus1].next_same2 == NULL)
                    three_pt_data->entry[locus1].next_same2 = t;
                else {
                    for (p = three_pt_data->entry[locus1].next_same2;
                         p->next_same2 != NULL; p = p->next_same2) {}
                    p->next_same2 = t;
                }

            } else {
                for (p = three_pt_data->entry[locus1].next; p->locus2 != locus2;
                     p = p->next) {}
                if (p != t) {
                    for (; p->next_same2 != NULL; p = p->next_same2) {}
                    p->next_same2 = t;
                }
            }
        }
    }
}


bool replace_triple(int locus1, int locus2, int locus3, real d1, real d2, real d3) {
    TRIPLE_LIST *p, *new_entry;

    if (three_pt_data->entry[locus1].locus2 == locus2 &&
        three_pt_data->entry[locus1].locus3 == locus3) {
        three_pt_data->entry[locus1].delta1 = (float) d1;
        three_pt_data->entry[locus1].delta2 = (float) d2;
        three_pt_data->entry[locus1].delta3 = (float) d3;
        return (TRUE);

    } else {

        for (p = three_pt_data->entry[locus1].next;
             p->next != NULL &&
             !(p->next->locus2 == locus2 && p->next->locus3 == locus3);
             p = p->next) {}

        if (p->next != NULL) {
            if (three_pt_data->unused == NULL)
                allocate_three_entries(INIT_3PT_SIZE(raw.num_markers));

            new_entry = three_pt_data->unused;
            three_pt_data->unused = three_pt_data->unused->next;
            new_entry->locus2 = locus2;
            new_entry->locus3 = locus3;
            new_entry->delta1 = (float) d1;
            new_entry->delta2 = (float) d2;
            new_entry->delta3 = (float) d3;
            new_entry->next = p->next->next;
            new_entry->next_same2 = p->next->next_same2;
            p->next = new_entry;
            /* fix the next_same2 list */
            if (three_pt_data->entry[locus1].locus2 == locus2) {
                if (three_pt_data->entry[locus1].next_same2->locus3 == locus3)
                    three_pt_data->entry[locus1].next_same2 = new_entry;
                else
                    p = three_pt_data->entry[locus1].next_same2;

            } else {
                for (p = three_pt_data->entry[locus1].next; p->locus2 != locus2;
                     p = p->next) {}
            }

            for (; p->next_same2 != NULL && p->next_same2->locus3 != locus3;
                   p = p->next_same2) {}

            if (p->next_same2 != NULL) {
                /* i.e., if it isn't the first in the next_same2 list */
                p->next_same2 = new_entry;
            }

        } else {
            three_pt_data->entry[locus1].next->delta1 = (float) d1;
            three_pt_data->entry[locus1].next->delta2 = (float) d2;
            three_pt_data->entry[locus1].next->delta3 = (float) d3;
        }
    }
    return (TRUE);
}

bool restore_triple(int locus1, int locus2, int locus3, real *d1, real *d2, real *d3) {
    /* d1= 123, d2= 132, d3= 213 -
       these transformations are correct!  */

    if (locus2 < locus1) {
        if (locus3 < locus2)
            return (get_triple(locus3, locus2, locus1, d1, d3, d2));
        else if (locus3 < locus1)
            return (get_triple(locus2, locus3, locus1, d2, d3, d1));
        else
            return (get_triple(locus2, locus1, locus3, d3, d2, d1));

    } else {
        if (locus3 < locus1)
            return (get_triple(locus3, locus1, locus2, d3, d1, d2));
        else if (locus3 < locus2)
            return (get_triple(locus1, locus3, locus2, d2, d1, d3));
        else
            return (get_triple(locus1, locus2, locus3, d1, d2, d3));
    }
}


bool get_triple(int locus1, int locus2, int locus3, real *d1, real *d2, real *d3) {
    TRIPLE_LIST *p;

    if (three_pt_data->entry[locus1].locus2 == locus2) {
        if (three_pt_data->entry[locus1].locus3 == locus3) {
            *d1 = (real) three_pt_data->entry[locus1].delta1;
            *d2 = (real) three_pt_data->entry[locus1].delta2;
            *d3 = (real) three_pt_data->entry[locus1].delta3;
            return (TRUE);
        } else {
            if (!same_2_loop(three_pt_data->entry[locus1].next_same2, locus3,
                             d1, d2, d3))
                return (FALSE);
            else
                return (TRUE);
        }
    } else {
        p = three_pt_data->entry[locus1].next;
        if (p == NULL) return (FALSE);
        while (p->locus2 != locus2) {
            if (p->next == NULL)
                return (FALSE);
            p = p->next;
        }
        if (!same_2_loop(p, locus3, d1, d2, d3))
            return (FALSE);
        else
            return (TRUE);
    }
}


bool same_2_loop(TRIPLE_LIST *p, int locus3, real *d1, real *d2, real *d3) {
    if (p == NULL) return (FALSE);
    while (p->locus3 != locus3) {
        if (p->next_same2 == NULL)
            return (FALSE);
        p = p->next_same2;
    }
    *d1 = (real) p->delta1;
    *d2 = (real) p->delta2;
    *d3 = (real) p->delta3;
    return (TRUE);
}


bool three_linked(int *locus, /* array of 3 locus numbers */ real lodbound, real thetabound, int num_links, bool sex) {
    int x, i, j, count;
    real lod, theta, thetam, thetaf;

    count = 0;
    for (x = 0; x < 3; x++) {
        i = locus[x];
        j = locus[(x + 1) % 3];
        if (!sex) {
            get_two_pt(i, j, &lod, &theta);
            if ((lod >= lodbound) && (theta <= thetabound)) count++;
        } else {
            get_sex_two_pt(i, j, &lod, &thetam, &thetaf);
            if (lod >= lodbound && (thetam <= thetabound || thetaf <= thetabound))
                count++;
        }
    }
    return (count >= num_links);
}


void compute_3pt(SEQ_NODE *seq, bool sex, real trip_err_rate, real *like, MAP *map) {
    int i, k;
    real best;

    k = 0;
    for_all_orders(seq, map) {
        init_for_ctm(map, sex, trip_err_rate != 0.0, MAYBE);
        if (trip_err_rate == LOCUS_ERROR_RATE)
            for (i = 0; i < 3; i++) map->error_rate[i] = error_rate[map->locus[i]];
        else if (trip_err_rate > 0.0)
            for (i = 0; i < 3; i++) map->error_rate[i] = trip_err_rate;
        converge_to_map(map);
        like[k++] = map->log_like;
    }
    for (i = 0, best = VERY_UNLIKELY; i < 3; i++) best = rmaxf(best, like[i]);
    for (i = 0; i < 3; i++) like[i] -= best;

    insert_triple(map->locus[0], map->locus[1], map->locus[2],
                  like[2], like[1], like[0]);
}


/***************************** Haplotype-Groups *****************************/

void setup_haplo_group(int *locus, int num_loci) {
    int first, i;

    first = locus[0];
    for (i = 0; i < num_loci; i++) {  /* build as a linked list */
        haplo_first[locus[i]] = first;
        if (i == num_loci - 1) haplo_next[locus[i]] = NO_LOCUS; /*end*/
        else haplo_next[locus[i]] = locus[i + 1];
    }
}


bool delete_haplo_groups(int *locus, int num_loci, int *old_locus, int *num_old) {
/* delete ALL haplo groups containing any of these loci */
    int i, j, next, any;

    any = FALSE;
    *num_old = 0;
    for (i = 0; i < num_loci; i++)
        if ((j = haplo_first[locus[i]]) != NO_LOCUS) {
            any = TRUE;
            do {
                old_locus[(*num_old)++] = j;
                next = haplo_next[j];
                haplo_next[j] = NO_LOCUS;
                haplo_first[j] = NO_LOCUS;
            } while ((j = next) != NO_LOCUS);
        }
    return (any);
}


/* msg is designed to look like those in chroms.c */
#define INSANE "%s- not the name of its haplotype-group, using %s\n"

bool force_haplo_sanity(int *locus, /* just a ptr to one int, maybe side-effected */ bool verbose) {
    if (haplo_first[*locus] == NO_LOCUS || *locus == haplo_first[*locus])
        return (TRUE);
    if (verbose) {
        sprintf(ps, INSANE, loc2str(*locus), rag(loc2str(haplo_first[*locus])));
        pr();
    }
    *locus = haplo_first[*locus];
    return (FALSE);
}


/**************** Classes ****************/
/* defined here because they are easiest to save/load in the table */

bool isa_class(char *name, int *num) {
    int i, classnum = -1;
    if (streq(name, "")) return (FALSE);
    for (i = 0; i < NUM_CLASSES; i++)
        if (class_name[i][0] != '\0' && xstreq(name, class_name[i])) {
            classnum = i;
            break;
        }
    if (classnum == -1) return (FALSE);
    if (num != NULL) *num = classnum;
    return (TRUE);
}


bool make_new_class(char *name, char **why_not) {
    int i, classnum = -1;

    if (!valid_name(name)) /* checks non-null */
    {
        *why_not = ptr_to("illegal name");
        return (FALSE);
    } else if (!valid_new_name(name)) /* check class names too */
    {
        *why_not = ptr_to("name is already in use");
        return (FALSE);
    }

    for (i = 0; i < NUM_CLASSES; i++)
        if (class_name[i][0] == '\0') {
            classnum = i;
            break;
        }
    if (classnum == -1) {
        *why_not = ptr_to("no more classes can be defined");
        return (FALSE);
    }
    strcpy(class_name[classnum], name);
    return (TRUE);
}


void print_class_names(void) {
    /* let print() auto_wrap... */
    int i;
    for (i = 0; i < NUM_CLASSES; i++) {
        print(class_name[i]);
        print(" ");
    }
}



/**************** Other Stuff, and Alloc/Save/Load/Bash Funcs ****************/


#define BASH_ORDER1 \
"Resetting group and order information for entire data set...\n"
#define BASH_ORDER2 \
"Resetting age, class, and error-prob information for joined loci...\n"
#define BASH_ORDER3 \
"Resetting two-point and three-point information for joined loci...\n"

//void return_to_unused(), remove_triple_list();

void bash_order_info(int *changed, int num_changed) {
    int a, b, i, j;
    TRIPLE_LIST *p, *q, *prev;
    /* this is merciless - might do *something* smarter */

    print(BASH_ORDER1); /* bash all groups/orders */
    num_groups = num_orders = 0;
    for (i = 0; i < raw.num_markers; i++) {
        my_group[i] = NO_GROUP;
        order_next[i] = NO_LOCUS;
    }
    print(BASH_ORDER2); /* bash these age/class/error_probs */
    for (i = 0; i < num_changed; i++) {
        modified[changed[i]] = FALSE;
        class[changed[i]] = NO_CLASS;
        error_rate[changed[i]] = DEFAULT_ERROR_RATE;
    }

    print(BASH_ORDER3);  /* bash two-point data, garbage collecting */
    for (i = 0; i < num_changed; i++)
        for (j = 0; j < raw.num_markers; j++) {
            if (changed[i] > j) {
                a = changed[i];
                b = j;
            } else {
                a = j;
                b = changed[i];
            }
            if (two_pt_data[a][b] == NULL) continue;
            two_pt_data[a][b]->used = FALSE;
            two_pt_data[a][b] = NULL;
            /* while (two_pt_used>0 && !(two_pt_list[two_pt_used-1]->used))
               two_pt_used--; */
            two_pt_collect = 0;
        }

    /* bash three-point data */
    for (i = 0; i < num_changed; i++) {
        three_pt_touched = TRUE;
        a = changed[i];
        /* first get rid of triple_list starting with marker a */
        if (three_pt_data->entry[a].next != NULL)
            remove_triple_list(three_pt_data->entry[a].next);
        three_pt_data->entry[a].next = NULL;
        three_pt_data->entry[a].next_same2 = NULL;
        three_pt_data->entry[a].locus2 = NO_LOCUS;
        three_pt_data->entry[a].locus3 = NO_LOCUS;
        three_pt_data->entry[a].delta1 = UNFILLED;
        three_pt_data->entry[a].delta2 = UNFILLED;
        three_pt_data->entry[a].delta3 = UNFILLED;

        /* now get rid of all others in which a appears */

        for (j = 0; j < a; j++) {
            prev = NULL;
            p = three_pt_data->entry[j].next;
            while (p != NULL) {
                if (p->locus2 == a || p->locus3 == a) {
                    /* go back through list and find next_same2 pointing at
                       this one, set it equal to whatever this next_same2
                       points at */
                    if (three_pt_data->entry[j].next_same2 == p) {
                        three_pt_data->entry[j].next_same2 = p->next_same2;
                    } else {
                        q = three_pt_data->entry[j].next;
                        while (q != p && q != NULL) {
                            if (q->next_same2 == p) {
                                q->next_same2 = p->next_same2;
                                break;
                            }
                            q = q->next;
                        }
                    }
                    /* fix the next pointer */
                    if (prev == NULL) three_pt_data->entry[j].next = p->next;
                    else prev->next = p->next;

                    /* set p for next pass through loop, prev remains the same */
                    q = p->next;

                    /* now p is completely disconnected, return it to unused */
                    return_to_unused(p);

                    p = q;
                } else {
                    prev = p;
                    p = p->next;
                }
            }

            /* now all have been checked except the initial entry[j] */
            if (three_pt_data->entry[j].locus2 == a ||
                three_pt_data->entry[j].locus3 == a) {
                p = three_pt_data->entry[j].next;
                if (p != NULL) {
                    three_pt_data->entry[j].locus2 = p->locus2;
                    three_pt_data->entry[j].locus3 = p->locus3;
                    three_pt_data->entry[j].delta1 = p->delta1;
                    three_pt_data->entry[j].delta2 = p->delta2;
                    three_pt_data->entry[j].delta3 = p->delta3;
                    three_pt_data->entry[j].next = p->next;
                    three_pt_data->entry[j].next_same2 = p->next_same2;
                    /* initial entry[j] has been set to be what was previously
                       entry[j].next, so send the entry[j].next back to unused */
                    return_to_unused(p);
                } else {
                    three_pt_data->entry[j].next = NULL;
                    three_pt_data->entry[j].next_same2 = NULL;
                    three_pt_data->entry[j].locus2 = NO_LOCUS;
                    three_pt_data->entry[j].locus3 = NO_LOCUS;
                    three_pt_data->entry[j].delta1 = UNFILLED;
                    three_pt_data->entry[j].delta2 = UNFILLED;
                    three_pt_data->entry[j].delta3 = UNFILLED;
                }
            }
        }
    }
}

void return_to_unused(TRIPLE_LIST *p) {
    /* Places a previously allocated but now unused elemtn on top of the
       unused stack */

    p->next = NULL;
    p->next_same2 = NULL;
    p->locus2 = p->locus3 = NO_LOCUS;
    p->delta1 = p->delta2 = p->delta3 = UNFILLED;

    p->next = three_pt_data->unused;
    three_pt_data->unused = p;
}

void remove_triple_list(TRIPLE_LIST *p) {
    if (p->next != NULL) remove_triple_list(p->next);
    return_to_unused(p);
}


void allocate_order_data(int num_markers) {
    int i;

    array(my_group, num_markers, int);
    array(haplo_first, num_markers, int);
    array(haplo_next, num_markers, int);
    array(order_first, num_markers, int);
    array(unorder_first, num_markers, int);
    array(order_next, num_markers, int);
    array(class, num_markers, int);
    array(modified, num_markers, int);
    array(error_rate, num_markers, real);
    matrix(note, num_markers, MAX_NOTE_LEN + 1, char);

    num_groups = num_orders = 0;
    for (i = 0; i < num_markers; i++) {
        modified[i] = FALSE;
        class[i] = NO_CLASS;
        my_group[i] = NO_GROUP;
        error_rate[i] = DEFAULT_ERROR_RATE;
        haplo_first[i] = haplo_next[i] = NO_LOCUS;
        order_first[i] = unorder_first[i] = NO_LOCUS;
    }

    matrix(class_name, NUM_CLASSES, NAME_LEN + 1, char);
    strcpy(class_name[0], "no_class");
    /* for (i=1; i<NUM_CLASSES; i++) sprintf(class_name[i],"class%d",i); */
    for (i = 1; i < NUM_CLASSES; i++) class_name[i][0] = '\0';
}


void free_order_data(int num_markers) {
    unarray(my_group, int);
    unarray(haplo_first, int);
    unarray(haplo_next, int);
    unarray(order_first, int);
    unarray(unorder_first, int);
    unarray(order_next, int);
    num_groups = 0;
    num_orders = 0;

    num_groups = 0;
    num_orders = 0;
    unarray(modified, int);
    unarray(error_rate, real);
    unmatrix(note, num_markers, char);

    unarray(class, int);
    unmatrix(class_name, NUM_CLASSES + 1, char);
}


void write_order_data(FILE *fp) {
    int i, locus;

    sprintf(ps, "*OrderInfo: %d %d\n", num_groups, num_orders);
    fpr(fp);
    for (locus = 0; locus < raw.num_markers; locus++) {
        sprintf(ps, "*%-8s %4d   %7.5lf %4d %4d %4d %4d %4d %4d %4d\n",
                raw.locus_name[locus], modified[locus], error_rate[locus],
                my_group[locus], haplo_first[locus], haplo_next[locus],
                order_first[locus], unorder_first[locus], order_next[locus],
                class[locus]);
        fpr(fp);
    }
    fprint(fp, "*Classes:\n");
    for (i = 0; i < NUM_CLASSES; i++) {
        fprint(fp, "*");
        fprint(fp, class_name[i]);
        fnl(fp);
    }
}


void read_order_data(FILE *fp) {
    int i, locus;
    int mod, group, first, next, ord_first, un_first, ord_next, class_num;
    real rate;
    char temp_locus_name[NAME_LEN + 2], word[TOKLEN + 1];

    fgetln(fp);
    if (sscanf(ln, "%s %d %d", word, &num_groups, &num_orders) != 3 ||
        !streq(word, "*OrderInfo:"))
        baddata("expected '*OrderInfo: # #'");

    for (locus = 0; locus < raw.num_markers; locus++) {
        fgetln(fp);

        if (!nstoken(&ln, sREQUIRED, temp_locus_name, NAME_LEN + 1) ||
            temp_locus_name[0] != '*' || len(temp_locus_name) < 2)
            baddata("expected *name");
        else if (!streq(raw.locus_name[locus], &temp_locus_name[1]))
            baddata("locus names don't match");

        if (sscanf(ln, "%d %lf %d %d %d %d %d %d %d", &mod, &rate, &group,
                   &first, &next, &ord_first, &un_first, &ord_next, &class_num) != 9)
            baddata("bad order info line");

        modified[locus] = mod;
        error_rate[locus] = rate;
        my_group[locus] = group;
        haplo_first[locus] = first;
        haplo_next[locus] = next;
        order_first[locus] = ord_first;
        unorder_first[locus] = un_first;
        order_next[locus] = ord_next;
        class[locus] = class_num;
    }
    fgetln(fp);
    if (!streq(ln, "*Classes:")) baddata("bad classes");
    for (i = 0; i < NUM_CLASSES; i++) {
        fgetln(fp);
        strcpy(class_name[i], ln + 1);
    }
}


void read_two_pt(FILE *fp) {
    int a, b, n;
    int i, j, k, n_unlinked = 0, n_miss = 0, n_real = 0, n_to_read;
    bool inc, missing_means_unlinked = FALSE; /* eg missing means missing */
    real theta, lod, thetam, thetaf, lodsex;
    TWO_PT_DATA *new_two;

    getdataln(fp);
    if (sscanf(ln, "%d %d %d", &n_real, &n_unlinked, &n_miss) != 3)
        baddata("expected two-pt count");
    if (n_unlinked > n_miss) missing_means_unlinked = TRUE;
    expand_two_pt(n_real + 2);

    if (missing_means_unlinked) {
        for (i = 0; i < raw.num_markers; i++)
            for (j = 0; j <= i; j++)
                two_pt_data[i][j] = UNLINKED_TWO_PT;
    }
    if (missing_means_unlinked) n_to_read = n_real + n_miss;
    else n_to_read = n_real + n_unlinked;

    run {
            for (k = 0; k < n_to_read; k++) {
                getdataln(fp);
                n = sscanf(ln, "%d %d %lf %lf %lf %lf %lf", &a, &b, &theta, &lod,
                           &thetam, &thetaf, &lodsex);
                new_two = get_next_two_pt_entry(a, b, &inc);
                if (n == 2) {
                    if (missing_means_unlinked) two_pt_data[a][b] = NULL;
                    else two_pt_data[a][b] = UNLINKED_TWO_PT;
                } else if (n == 4) {
                    new_two->theta[NOSEX] = theta;
                    new_two->lodscore[NOSEX] = lod;
                    two_pt_data[a][b] = new_two;
                    new_two->used = TRUE;
                    two_pt_used++;
                } else if (n == 7) {
                    new_two->theta[NOSEX] = theta;
                    new_two->lodscore[NOSEX] = lod;
                    new_two->theta[MALE] = thetam;
                    new_two->theta[FEMALE] = thetaf;
                    new_two->lodscore[SEXSPEC] = lodsex;
                    two_pt_data[a][b] = new_two;
                    new_two->used = TRUE;
                    if (inc) two_pt_used++;
                }
            }
        } except_when(ENDOFILE) {}
}


void write_two_pt(FILE *fp) {
    int i, j, n_unlinked = 0, n_miss = 0, n_real = 0;
    bool missing_means_unlinked = FALSE; /* eg missing means missing */

    for (i = 0; i < raw.num_markers - 1; i++)
        for (j = 0; j <= i; j++) {
            if (two_pt_data[i][j] == NULL) n_miss++;
            else if (two_pt_data[i][j] == UNLINKED_TWO_PT) n_unlinked++;
            else n_real++;
        }

    if (n_unlinked > n_miss) missing_means_unlinked = TRUE;
    sprintf(ps, "%d %d %d \n", n_real, n_unlinked, n_miss);
    fpr(fp);

    for (i = 0; i < raw.num_markers - 1; i++) {
        for (j = 0; j <= i; j++) {
            if (two_pt_data[i][j] == NULL) {
                if (!missing_means_unlinked) continue;
                else sprintf(ps, "%d %d\n", i, j);
                fpr(fp);
            } else if (two_pt_data[i][j] == UNLINKED_TWO_PT) {
                if (missing_means_unlinked) continue;
                sprintf(ps, "%d %d\n", i, j);
                fpr(fp);
            } else if (raw.data_type == F2) {
                sprintf(ps, "%d %d %.3lf %.3lf\n", i, j,
                        two_pt_data[i][j]->theta[NOSEX],
                        two_pt_data[i][j]->lodscore[NOSEX]);
                fpr(fp);
            } else {
                sprintf(ps, "%d %d %.3lf %.3lf %.3lf %.3lf %.3lf\n", i, j, /* lf */
                        two_pt_data[i][j]->theta[NOSEX],
                        two_pt_data[i][j]->lodscore[NOSEX],
                        two_pt_data[i][j]->theta[MALE],
                        two_pt_data[i][j]->theta[FEMALE],
                        two_pt_data[i][j]->lodscore[SEXSPEC]);
                fpr(fp);
            }
        }
    }
}


void read_three_pt(FILE *fp) {
    int i, j, k;
    real d1, d2, d3;

    while ((fscanf(fp, "%d %d %d %lf %lf %lf\n", &i, &j, &k, &d1, &d2, &d3)) == 6) {
        insert_triple(i, j, k, d1, d2, d3);
    }
    three_pt_touched = FALSE;
}


void write_three_pt(FILE *fp) {
    int i;
    TRIPLE_LIST *p;

    for (i = 0; i < raw.num_markers; i++) {
        if (three_pt_data->entry[i].locus2 != NO_LOCUS) {
            sprintf(ps, "%d %d %d %.3lf %.3lf %.3lf\n", i,
                    three_pt_data->entry[i].locus2, three_pt_data->entry[i].locus3,
                    three_pt_data->entry[i].delta1, three_pt_data->entry[i].delta2,
                    three_pt_data->entry[i].delta3);
            fpr(fp);
            p = three_pt_data->entry[i].next;
            while (p != NULL) {
                sprintf(ps, "%d %d %d %.3lf %.3lf %.3lf\n", i, p->locus2, p->locus3,
                        p->delta1, p->delta2, p->delta3);
                fpr(fp);
                p = p->next;
            }
        }
    }
}
