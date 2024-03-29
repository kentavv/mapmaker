/******************************************************************************

  ####   #####   #####   ######  #####    ####            ####
 #    #  #    #  #    #  #       #    #  #               #    #
 #    #  #    #  #    #  #####   #    #   ####           #
 #    #  #####   #    #  #       #####        #   ###    #
 #    #  #   #   #    #  #       #   #   #    #   ###    #    #
  ####   #    #  #####   ######  #    #   ####    ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#include "mapm.h"

bool good_seed(MAP *map, MAP *temp_map, real thresh);

void alloc_3pt_matrix(int num_group);

void free_3pt_matrix(void);

bool is_excluded(int a, int b, int c);

void left_exclusion(bool *excluded, int *locus, int num_loci, int newmarker);

void right_exclusion(bool *excluded, int *locus, int num_loci, int newmarker);

void middle_exclusion(bool *excluded, int *locus, int num_loci, int leftmost, int rightmost, int init_window_size, int newmarker);

bool find_seed_order(bool is_subset, int *locus, int num_loci, int size, int max_tries, real thresh, MAP *map, MAP *temp_map, bool **temp);

bool extend_order_by_1(MAP *placed, PLACEME **unplaced, int *num_unplaced, real npt_thresh, bool *contradiction, PLACE **placements, MAP *temp_map);

void randomize_markers_to_try(PLACEME **unplaced, int num_unplaced);

void rank_markers_to_try(PLACEME **unplaced, int num_unplaced);

int compare_markers_to_try(const void *a, const void *b);

void add_to_order(int how, int *order, int *num_ordered, PLACEME **unplaced, int index, int *num_unplaced);

bool middles_excluded(int *excluded, int num_order);

char ***three_pt_excluded = NULL;
int three_pt_size = 0;
int *three_pt_index = NULL;

#define EXCLUDED   0
#define ALLOWED    1
#define UNCOMPUTED 2

//void alloc_3pt_matrix(); /* args: int num_loci; sets a global, will free
//   and realloc if it's too small, use to prealloc for efficiency */
//void free_3pt_matrix();



/**************** Two-point, Haplotypes ****************/

void get_linkage_group(int *locus, int *num_loci, int *linkage_group, int *group_size, real lodbound, real thetabound) {
    int i, j, linked, *unlinked_list, unlinked_size, oldsize, a, b;
    real lod, theta, thetam, thetaf;

    /* start linkage_group with first marker in the list of loci */
    linkage_group[0] = locus[0];
    *group_size = 1;
    unlinked_list = locus;
    *num_loci -= 1;
    unlinked_size = *num_loci;
    locus = &unlinked_list[1];

    do {
        oldsize = unlinked_size;
        unlinked_size = 0;
        for (i = 0; i < oldsize; i++) { /* for all loci a NOT in the group */
            a = locus[i];
            for (j = 0, linked = FALSE; j < *group_size; j++) { /* for all in group*/
                b = linkage_group[j];
                if (!sex_specific) {
                    get_two_pt(a, b, &lod, &theta);
                    if (lod >= lodbound && theta <= thetabound) {
                        linked = TRUE;
                        break;
                    }
                } else {
                    get_sex_two_pt(a, b, &lod, &thetam, &thetaf);
                    if (lod > lodbound &&
                        (thetam <= thetabound || thetaf <= thetabound)) {
                        linked = TRUE;
                        break;
                    }
                }
            }
            if (linked) {
                linkage_group[*group_size] = a;
                ++*group_size;
            } else {
                /* Effectively, write back into the locus array, relying on the
                   fact that unlinked_size<=*num_loci */
                unlinked_list[unlinked_size] = a;
                ++unlinked_size;
            }
        } /* for i */
        locus = unlinked_list;
        *num_loci = unlinked_size;
    } while (unlinked_size != oldsize && unlinked_size != 0);
}


void find_haplo_group(int *locus, int *num_loci, int *haplo_group, int *num_haplo, int *old_obs, int *new_obs) {
    int i, num_inf, num_dom, best_i = 0, foo;
    int best_score, score;

    /* get the most informative marker */
    best_score = 0;
    for (i = 0; i < *num_loci; i++) {
        f2_genotype(locus[i], FALSE, old_obs);
        num_inf = f2_count_infs(&num_dom, &foo, old_obs);
        score = 2 * num_inf - num_dom;
        if (score > best_score) {
            best_score = score;
            best_i = i;
        }
    }

    /* delete it, and make it the haplo group */
    haplo_group[0] = locus[best_i];
    *num_haplo = 1;
    f2_genotype(locus[best_i], FALSE, old_obs);
    locus[best_i] = locus[--*num_loci];
    if (*num_loci == 0) return;

    do { /* find one to add to the group */
        best_score = 0;
        for (i = 0; i < *num_loci; i++)
            if (merge_genotypes(locus[i], old_obs, new_obs)) { /* if no recs */
                num_inf = f2_count_infs(&num_dom, &foo, new_obs);
                score = 2 * num_inf - num_dom;
                if (score > best_score) {
                    best_score = score;
                    best_i = i;
                }
            }
        if (best_score != 0) { /* if we found any with no recs */
            merge_genotypes(locus[best_i], old_obs, old_obs); /*updates old_obs*/
            haplo_group[*num_haplo] = locus[best_i];
            ++*num_haplo; /* add it */
            locus[best_i] = locus[--*num_loci]; /* delete it */
        }
    } while (*num_loci > 0 && best_score != 0); /* got one and more to try */
}


/**************** Three Point ****************/


void alloc_3pt_matrix(int num_group) {
    if (three_pt_excluded != NULL && num_group <= three_pt_size) return;

    run {
            if (three_pt_excluded != NULL) free_3pt_matrix();
            three_pt_size = num_group;
            three_pt_excluded = alloc_char_3d_matrix(num_group, num_group, num_group);
            array(three_pt_index, raw.num_markers, int);
        } when_aborting {
        free_3pt_matrix();
        relay_messages;
    }
}


void
free_3pt_matrix(void) {
    free_char_3d_matrix(three_pt_excluded, three_pt_size, three_pt_size);
    unarray(three_pt_index, int);
    three_pt_excluded = NULL;
    three_pt_size = 0;
}


void setup_3pt_data(int *locus, int num_loci, real threshold) {
    int i, j, k;
    real d1, d2, d3;

    alloc_3pt_matrix(num_loci);

    /* three_pt_index[] is lookup by raw locus# to three_pt_exclusion index */
    for (i = 0; i < raw.num_markers; i++) three_pt_index[i] = -1;
    for (i = 0; i < num_loci; i++) three_pt_index[locus[i]] = i;

    for (i = 0; i < num_loci; i++)
        for (j = 0; j < num_loci; j++)
            for (k = 0; k < num_loci; k++)
                three_pt_excluded[i][j][k] = UNCOMPUTED;

    for (i = 0; i < num_loci - 2; i++) {
        for (j = i + 1; j < num_loci - 1; j++) {
            for (k = j + 1; k < num_loci; k++) {

                if (!restore_triple(locus[i], locus[j], locus[k], &d1, &d2, &d3))
                    continue;

                /* maybe later compute on the fly ADD keep_user_amused
                   compute_triple(locus[i],locus[j],locus[k],&d1,&d2,&d3); */

                if (d1 < threshold) {
                    three_pt_excluded[i][j][k] = EXCLUDED;
                    three_pt_excluded[k][j][i] = EXCLUDED;
                } else {
                    three_pt_excluded[i][j][k] = ALLOWED;
                    three_pt_excluded[k][j][i] = ALLOWED;
                }

                if (d2 < threshold) {
                    three_pt_excluded[i][k][j] = EXCLUDED;
                    three_pt_excluded[j][k][i] = EXCLUDED;
                } else {
                    three_pt_excluded[i][k][j] = ALLOWED;
                    three_pt_excluded[j][k][i] = ALLOWED;
                }

                if (d3 < threshold) {
                    three_pt_excluded[j][i][k] = EXCLUDED;
                    three_pt_excluded[k][i][j] = EXCLUDED;
                } else {
                    three_pt_excluded[j][i][k] = ALLOWED;
                    three_pt_excluded[k][i][j] = ALLOWED;
                }
            }
        }
    }
}

bool is_excluded(int a, int b, int c) {
    return (three_pt_excluded[three_pt_index[a]][three_pt_index[b]]
            [three_pt_index[c]] == EXCLUDED);
}


void left_exclusion(bool *excluded, int *locus, int num_loci, int newmarker) {
    int middle, end, i;

    for (middle = num_loci - 2; middle >= 0; middle--) {
        for (end = num_loci - 1; end > middle; end--) {
            if (is_excluded(newmarker, locus[middle], locus[end])) {
                for (i = 0; i <= middle; i++)
                    excluded[i] = TRUE;
                return;
            }
        }
    }
}


void right_exclusion(bool *excluded, int *locus, int num_loci, int newmarker) {
    int middle, end, i;

    for (middle = 1; middle <= num_loci - 1; middle++) {
        for (end = 0; end < middle; end++) {
            if (is_excluded(locus[end], locus[middle], newmarker)) {
                for (i = num_loci; i > middle; i--)
                    excluded[i] = TRUE;
                return;
            }
        }
    }
}


void middle_exclusion(bool *excluded, int *locus, int num_loci, int leftmost, int rightmost,  /* presumably, the left and rightmost intervals? */
                      int init_window_size, int newmarker) {
    int first, last, start, end, width, i, j;

    for (width = init_window_size; width >= 1; width--) {
        first = max(1, leftmost - width + 1);
        last = min(rightmost, num_loci - width);
        for (i = first; i <= last; i++) {
            start = i;
            end = i + width - 1;
            if (is_excluded(locus[start - 1], newmarker, locus[end])) {
                for (j = start; j <= end; j++)
                    excluded[j] = TRUE;
                if (start > leftmost && width > 1)
                    middle_exclusion(excluded, locus, num_loci, leftmost, start - 1,
                                     width - 1, newmarker);
                if (end == rightmost) return;
                leftmost = end + 1;
            }
        }
    }
}


/* orig frame:         a - b - c - d - e - f - g
   orig frame index:   0 - 1 - 2 - 3 - 4 - 5 - 6    num_placed=7
   exclusion index:  0 - 1 - 2 - 3 - 4 - 5 - 6 - 7                */

int three_pt_exclusions(int *order, int num_placed, int newmarker, bool *excluded /* side-effected */) {
    int i, j, left, right, num_places;
    for (i = 0; i < num_placed + 1; i++) excluded[i] = FALSE;

    left_exclusion(excluded, order, num_placed, newmarker);
    right_exclusion(excluded, order, num_placed, newmarker);

    left = 1;
    right = num_placed;
    j = 1;
    while (j < num_placed && excluded[j]) j++;
    if (j != num_placed)
        middle_exclusion(excluded, order, num_placed, left, right, right - left,
                         newmarker);
    for (num_places = 0, i = 0; i < num_placed + 1; i++)
        if (!excluded[i]) ++num_places;
    return (num_places);
}


bool three_pt_verify(int *locus, int num_loci, int window_size) {
    int i, j, k;
    if (num_loci < 3 || window_size < 3) return (TRUE); /* crash? */

    for (i = 0; i < num_loci - 2; i++)
        for (k = i + 2; k < num_loci && k < i + window_size; k++)
            for (j = i + 1; j < k; j++)
                if (is_excluded(locus[i], locus[j], locus[k])) return (FALSE);

    return (TRUE);
}


/**************** N-Pt Maker Placer ****************/

void place_locus(MAP *map,          /* map to place locus into, should be ctm-ed! */
                 int locus, int start,
                 int finish, /* The leftmost and rightmost frame markers (i) to use */
                 bool *excluded, /* [interval#], one must be FALSE, NOT set */
                 PLACE **result, /* [interval#]->foo side-effected like=NO_LIKE if untried */
                 int *best_pos,   /* side effected */
                 MAP *best_map, MAP *temp  /* max_loci must be >=finish-start+2 */) {
    int i, j, k, num_allowed, last;
    real best_like, theta;
    real theta1, theta2, dist1, dist2, distframe, distscale;
#ifdef THREAD
    int nqueued, id;
#else
    real lod;
#endif

    printf("Placing locus %s\n", loc2str(locus));

    /* Generate N-point maps from left to right markers with locus 
       inserted into each of the allowed positions, even one */

    for (j = 0, num_allowed = 0; j <= map->num_loci; j++) {
        if (!excluded[j]) num_allowed++;
        result[j]->like = NO_LIKE;
        result[j]->dist = NO_THETA;
        result[j]->net_error = 0.0;
        result[j]->worst_error = 0.0;
        result[j]->zero = FALSE;
    }

    if (sex_specific || num_allowed == 0 || best_map->max_loci < finish - start + 2 ||
        temp->max_loci < finish - start + 2)
        send(CRASH);

    last = finish + 1;
    best_like = VERY_UNLIKELY;
#ifndef THREAD
    for (i=start; i<=last; i++) /* i=the interval # to place in */
      if (!excluded[i]) {
      clean_map(temp); k=0;
      if (i==0) temp->locus[k++]=locus;
      for (j=start; j<=finish; j++) { /* locus #s */
          temp->locus[k++]=map->locus[j];
          if (j+1==i) temp->locus[k++]=locus;
      }
      temp->num_loci= k;
      init_rec_fracs(temp);
      converge_to_map(temp);
      result[i]->like=temp->log_like;

      if (temp->log_like>best_like) { /* keep the best */
          best_like=temp->log_like; *best_pos=i;
          mapcpy(best_map,temp,FALSE);
      }

      /* normalize placement dists by any expansion of this interval
         if i=first or last interval in the FRAMEWORK, then
         dist= the correct rec frac, otherwise the placement dist =
         DL*(Dframe/(DL+DR)) */
      if (i==start) {
          result[i]->dist= theta= temp->rec_frac[0][0];
          if (theta<ZERO_DIST) result[i]->zero=TRUE;
      } else if (i==last) {
          result[i]->dist= theta= temp->rec_frac[i-start-1][MALE];
          if (theta<ZERO_DIST) result[i]->zero=TRUE;
      } else {
          theta1= temp->rec_frac[i-start-1][MALE];
          theta2= temp->rec_frac[i-start][MALE];
          if (theta1<ZERO_DIST || theta2<ZERO_DIST) result[i]->zero=TRUE;
          dist1= (*mapfunction->rec_to_dist)(theta1);
          dist2= (*mapfunction->rec_to_dist)(theta2);
          distframe= (*mapfunction->rec_to_dist)(map->rec_frac[i-1][0]);
          distscale= dist1*(distframe/(dist1+dist2));
          result[i]->dist= (*mapfunction->dist_to_rec)(distscale);
      }

      /* determine the error lods for any particular placement.
         Here is a stupid algorithm which looks only at this locus's
         lods. we realy need to look at flankers too! */
      result[i]->net_error=0.0; result[i]->worst_error=0.0;
      if (temp->allow_errors && (i>start || i<last) &&
          temp->error_rate[i-start]>0.0) {
          for (j=0; j<raw.data.f2.num_indivs; j++) { /* must be F2 */
          lod= temp->error_lod[i-start][j];
          if (lod>result[i]->worst_error) result[i]->worst_error=lod;
          if (lod>=error_lod_thresh)      result[i]->net_error+=lod;
          }
      }
      } /* next interval i */
#else
    for (nqueued = 0, i = start; nqueued < queue_size() && i <= last; i++) {/* i=the interval # to place in */
        if (!excluded[i]) {
            clean_map(temp);
            k = 0;
            if (i == 0) temp->locus[k++] = locus;
            for (j = start; j <= finish; j++) { /* locus #s */
                temp->locus[k++] = map->locus[j];
                if (j + 1 == i) temp->locus[k++] = locus;
            }
            temp->num_loci = k;
            add_converge_request(temp, i);
            nqueued++;
        }
    }
    while (nqueued) {
        int converged = get_converge_result(temp, &id);
        if (converged) {
            result[id]->like = temp->log_like;

            if (temp->log_like > best_like) { /* keep the best */
                best_like = temp->log_like;
                *best_pos = id;
                mapcpy(best_map, temp, FALSE);
            }

            /* normalize placement dists by any expansion of this interval
               if i=first or last interval in the FRAMEWORK, then
               dist= the correct rec frac, otherwise the placement dist =
               DL*(Dframe/(DL+DR)) */
            if (id == start) {
                result[id]->dist = theta = temp->rec_frac[0][0];
                if (theta < ZERO_DIST) result[id]->zero = TRUE;
            } else if (id == last) {
                result[id]->dist = theta = temp->rec_frac[id - start - 1][MALE];
                if (theta < ZERO_DIST) result[id]->zero = TRUE;
            } else {
                theta1 = temp->rec_frac[id - start - 1][MALE];
                theta2 = temp->rec_frac[id - start][MALE];
                if (theta1 < ZERO_DIST || theta2 < ZERO_DIST) result[id]->zero = TRUE;
                dist1 = (*mapfunction->rec_to_dist)(theta1);
                dist2 = (*mapfunction->rec_to_dist)(theta2);
                distframe = (*mapfunction->rec_to_dist)(map->rec_frac[id - 1][0]);
                distscale = dist1 * (distframe / (dist1 + dist2));
                result[id]->dist = (*mapfunction->dist_to_rec)(distscale);
            }

#if 0
            /* determine the error lods for any particular placement.
               Here is a stupid algorithm which looks only at this locus's
               lods. we realy need to look at flankers too! */
            result[id]->net_error=0.0; result[id]->worst_error=0.0;
            if (temp->allow_errors && (id>start || id<last) &&
                temp->error_rate[id-start]>0.0) {
              for (j=0; j<raw.data.f2.num_indivs; j++) { /* must be F2 */
                lod= temp->error_lod[id-start][j];
                if (lod>result[id]->worst_error) result[id]->worst_error=lod;
                if (lod>=error_lod_thresh)      result[id]->net_error+=lod;
              }
            }
#endif
        } else {
            sprintf(ps, "Did not converge: ");
            pr();
            for (j = 0; j < map->num_loci; j++) {
                print(loc2str(map->locus[j]));
            }
            nl();
        }
        nqueued--;
        for (; i <= last; i++) {/* i=the interval # to place in */
            if (!excluded[i]) {
                clean_map(temp);
                k = 0;
                if (i == 0) temp->locus[k++] = locus;
                for (j = start; j <= finish; j++) { /* locus #s */
                    temp->locus[k++] = map->locus[j];
                    if (j + 1 == i) temp->locus[k++] = locus;
                }
                temp->num_loci = k;
                add_converge_request(temp, i);
                nqueued++;
                i++;
                break;
            }
        }
    }
#endif
    if (best_like == VERY_UNLIKELY) send(CRASH);

    /* normalize the likes */
    for (i = start; i <= last; i++)
        if (!excluded[i]) {
            if (i == *best_pos) result[i]->like = 0.0;
            else result[i]->like -= best_like;
        }

    /* if all placements with nearly 0 like are zero dist from a framework 
       locus, remove all but the best (thus place_this will place the marker 
       uniquely) */
    /* here we have a simple (and likely not complete enough) test... 
       are the likelihoods in consecutive intervals really zero, 
       and does each placement have SOME zero dist. */
    if (result[*best_pos]->zero) {
        for (i = *best_pos + 1; !excluded[i] && i <= last &&
                                result[i]->like > ZERO_PLACE && result[i]->zero; i++)
            result[i]->like = ZERO_LIKE;
        for (i = *best_pos - 1; !excluded[i] && i >= start &&
                                result[i]->like > ZERO_PLACE && result[i]->zero; i--)
            result[i]->like = ZERO_LIKE;
    }
}


int npt_exclusions(
/* return #ints ok */
        PLACE **place,  /* results from above */
        int num_loci,  /* in order */
        real like_thresh,
        real one_err_thresh,
        real net_err_thresh,
        bool *excluded,     /* contents side-effected */
        bool *zero_place,   /* side-effected, T if all acceptable places are zeros */
        int *off_ends,     /* side-effected, >0 if all acceptable places are end(s) */
        real *error_lod,    /* side-effected, is greatest for all acceptable places */
        bool *single_error, /* side-effected, T if error_place && err due to summing */
        real *next_best,    /* side-effected, best acceptable but non_zero like (<0) */
        real *best_unallowed /* side-effected, best like >threshold */
) {
    int i, best_i, num, zero, one, ends;
    real second, unallowed, this_like, err, x;
    if (like_thresh >= 0.0) send(CRASH);

    num = 0;
    zero = TRUE;
    err = NO_ERRORS;
    one = FALSE, ends = 2;
    best_i = -1;
    second = NO_LIKE;
    unallowed = NO_LIKE;  /* may stay this way */

    for (i = 0; i <= num_loci; i++) {
        this_like = place[i]->like;
        if (this_like != NO_LIKE && this_like != ZERO_LIKE &&
            this_like > like_thresh) { /* it's acceptable */
            num++;
            if (excluded != NULL) excluded[i] = FALSE;
            if (this_like == 0.0 && best_i == -1) best_i = i;
            else if (this_like > second) second = this_like;
            if ((x = place[i]->worst_error) >= one_err_thresh) {
                err = x;
                one = TRUE;
            } else if ((x = place[i]->net_error) >= net_err_thresh) err = x;

        } else { /* it's unacceptable */
            if (excluded != NULL) excluded[i] = TRUE;
            if (i == 0 || i == num_loci) ends--;
            if (this_like != NO_LIKE) {
                if (this_like != ZERO_LIKE) zero = FALSE;
                if (this_like > unallowed) unallowed = this_like;
            }
        }
    }
    if (best_i == -1) send(CRASH);
    if (place[best_i]->worst_error >= one_err_thresh) one = TRUE; /* override */

    if (zero_place != NULL) *zero_place = (num == 1 && zero);
    if (error_lod != NULL) *error_lod = err;
    if (single_error != NULL) *single_error = one;
    if (off_ends != NULL) *off_ends = (ends > 0 ? TRUE : FALSE);
    if (next_best != NULL) *next_best = second;
    if (best_unallowed != NULL) *best_unallowed = unallowed;
    return (num);
}


/**************** Automatic Sub-Order Finder and Stuff ****************/


void informative_subset(int *locus, int num_loci, int min_infs, real min_theta, bool codom_only, bool haplos, int *subset, int *num) {
    int n_infs, n_dom_obs, n_het_obs, type, i, j, a, b, delete, deleted;
    int a_infs, b_infs;
    real lod, theta;
    *num = 0;

    /* filter for informativeness */
    for (i = 0; i < num_loci; i++) {
        f2_genotype(locus[i], haplos, observations);
        n_infs = f2_count_infs(&n_dom_obs, &n_het_obs, observations);
        type = raw.data.f2.cross_type;

        if (n_infs < min_infs) continue;
        if (type == F3_SELF || type == F2_INTERCROSS) {
            if (codom_only && n_dom_obs != 0) continue;
        }
        subset[(*num)++] = locus[i];
    }
    if (*num < 2) return;

    /* find zeros - a bit kludgey */
    do {
        deleted = FALSE;
        for (i = 0; i < *num - 1; i++) {
            for (j = i + 1; j < *num; j++) {
                a = subset[i];
                b = subset[j];
                get_two_pt(a, b, &lod, &theta);
                if (theta < min_theta) {
                    /* delete least informative marker, to bad we don't save
                       the infomativeness values */
                    f2_genotype(a, haplos, observations);
                    a_infs = f2_count_infs(&n_dom_obs, &n_het_obs, observations);
                    f2_genotype(b, haplos, observations);
                    b_infs = f2_count_infs(&n_dom_obs, &n_het_obs, observations);
                    /* do something smarter with mix/dr markers */
                    if (a_infs > b_infs) delete = (coin_flip() ? i : j);
                    else if (a_infs > b_infs) delete = j;
                    else delete = i;
                    if (delete < *num - 1) subset[delete] = subset[*num - 1];
                    *num -= 1;
                    deleted = TRUE;
                    break;
                }
            }
            if (deleted) break;
        }
    } while (deleted);
}


#define irand(n) (int)(randnum()*(real)n)
#define SEED_FAILED \
  "Failed to find a starting order at log-likelihood threshold %.2lf.\n"
#define SEARCH_INF \
"Searching for a unique starting order containing %d of %d informative loci...\n"
#define SEARCH_ALL \
  "Searching for a starting order containing %d of all %d loci...\n"
#define MAX_SAFETY 100

bool find_seed_order(bool is_subset, int *locus, int num_loci, int size, int max_tries, real thresh, MAP *map, MAP *temp_map,
                     bool **temp /* [num_loci][num_loci] */) {
    int try, pick, match, safety, i, j, k;
    /* size must be >=num_loci; if size==num_loci, there is only one choice.
       if size==num_loci-1, there are only num_loci choices, else there are
       many orders. if size>num_loci it's an error */

    if (num_loci < 3) send(CRASH);
    if (size > num_loci) size = num_loci;
    try = 0;
    temp_map->num_loci = size;
    map->num_loci = 0;

    if (is_subset) sprintf(ps, SEARCH_INF, size, num_loci);
    else sprintf(ps, SEARCH_ALL, size, num_loci);
    pr();

    if (size == num_loci) {
        for (i = 0; i < num_loci; i++) temp_map->locus[i] = locus[i];
        keep_user_amused("subset", try + 1, 0);
        if (good_seed(map, temp_map, thresh)) return (TRUE);
        sprintf(ps, SEED_FAILED, -thresh);
        pr();
        return (FALSE);

    } else if (size == num_loci - 1) {
        for (j = 0; j < num_loci; j++) { /* locus to punt */
            k = 0;
            for (i = 0; i < num_loci; i++)
                if (i != j) temp_map->locus[k++] = locus[i];
            keep_user_amused("subset", ++try, 0);
            if (good_seed(map, temp_map, thresh)) return (TRUE);
            if (try > max_tries) break;
        }
        sprintf(ps, SEED_FAILED, -thresh);
        pr();
        return (FALSE);

    } else /* size << num_loci */ {
        do { /* try a subset */
            safety = 0;
            do { /* find an untried subset */
                /* first use temp as a flag to indicate loci in order */
                for (i = 0; i < num_loci; i++) temp[try][i] = FALSE;
                for (i = 0; i < size; i++) {
                    do pick = irand(num_loci); while (temp[try][pick]);
                    temp_map->locus[i] = locus[pick];
                    temp[try][pick] = TRUE;
                }
                /* have we tried this subset before? */
                sort_loci(temp_map->locus, size);
                for (j = 0, match = FALSE; j < try; j++) {
                    for (match = TRUE, k = 0; k < size; k++)
                        if (temp_map->locus[k] != temp[j][k]) {
                            match = FALSE;
                            break;
                        }
                    if (match) break;
                }
                /* safety is a very lame way to see when we should give up,
                   because we have tried all possible subsets */
                safety++;
            } while (match && safety <= MAX_SAFETY);
            if (match) break;
            for (k = 0; k < size; k++) temp[try][k] = temp_map->locus[k];
            keep_user_amused("subset", ++try, 0);
            /* for (i=0; i<size; i++)
               { sprintf(ps,"%d ",temp_map->locus[i]+1); pr(); } nl(); */
            if (good_seed(map, temp_map, thresh)) return (TRUE);
        } while (try <= max_tries);
        sprintf(ps, SEED_FAILED, -thresh);
        pr();
        return (FALSE);
    }
}


bool good_seed(MAP *map, MAP *temp_map, real thresh) {
    real best2 = VERY_UNLIKELY, best = VERY_UNLIKELY;
#ifdef THREAD
    int nqueued, id, j;
#endif

    make_compare_seq(temp_map->locus, temp_map->num_loci, 0, temp_map->num_loci);
#ifndef THREAD
    for_all_orders(seq,temp_map) {
    if (use_three_pt &&
        !three_pt_verify(temp_map->locus,temp_map->num_loci,
                 three_pt_window)) continue;
    init_rec_fracs(temp_map);
    converge_to_map(temp_map);
    if (temp_map->log_like>best)
      { best2=best; best=temp_map->log_like; mapcpy(map,temp_map,TRUE); }
    else if (temp_map->log_like>best2) { best2=temp_map->log_like; }
    }
#else
    for (nqueued = 0, Oagain = TRUE, reset_seq(seq, temp_map != NULL);
         nqueued < queue_size() && Oagain && clean_map(temp_map) && get_map_order(seq, temp_map);
         Oagain = perm_seq(seq, FALSE, FALSE)) {
        if (use_three_pt &&
            !three_pt_verify(temp_map->locus, temp_map->num_loci,
                             three_pt_window))
            continue;
        add_converge_request(temp_map, 0);
        nqueued++;
    }
    while (nqueued) {
        int converged = get_converge_result(temp_map, &id);
        if (converged) {
            if (temp_map->log_like > best) {
                best2 = best;
                best = temp_map->log_like;
                mapcpy(map, temp_map, TRUE);
            }
            else if (temp_map->log_like > best2) { best2 = temp_map->log_like; }
        } else {
            sprintf(ps, "Did not converge: ");
            pr();
            for (j = 0; j < map->num_loci; j++) {
                print(loc2str(map->locus[j]));
            }
            nl();
        }
        nqueued--;
        for (;
                Oagain && clean_map(temp_map) && get_map_order(seq, temp_map);
                Oagain = perm_seq(seq, FALSE, FALSE)) {
            if (use_three_pt &&
                !three_pt_verify(temp_map->locus, temp_map->num_loci,
                                 three_pt_window))
                continue;
            add_converge_request(temp_map, 0);
            nqueued++;
            Oagain = perm_seq(seq, FALSE, FALSE);
            break;
        }
    }
#endif
    if (best == VERY_UNLIKELY) return (FALSE);
    else if (best2 == VERY_UNLIKELY || best2 - best < thresh) {
        sprintf(ps, "Got one at log-likelihood %.2lf\n", best - best2);
        pr();
        return (TRUE);
    } else return (FALSE);
}



/**************** Automatic Order Maker ****************/

#define BY_NPT 1
#define BY_3PT 0
#define BY_NPT_ENDS (-1)
#define BY_NPT_ERRORS (-2)
#define BY_3PT_ENDS (-3)
#define NO_PLACES 9999

#define ALL_EXCLUDED \
"Three-point analysis excludes %d %smarker%s from all intervals\n"

#define sf_locnum(locus, paren) \
  sprintf(ps,"%s%d%s%s",(paren ? "(":""),(locus)+1, \
     (use_haplotypes && haplotyped(locus) ? "+":""),(paren ? ")":""))

void
extend_order(
        MAP *placed,        /* starting order, must be total long */
        PLACEME **unplaced,
        int *num_unplaced,
        real npt_thresh,        /* both are side-effected, and may have NO_LOCUS */
        bool print_anyway
) {
    int i, j, total = 0;
    bool placed_any, contradiction;
    PLACE **placements = NULL;
    MAP *temp_map = NULL;

    run {
            total = placed->num_loci + *num_unplaced;
            /* if (*num_unplaced==0) return; let it fall through */
            if (placed->num_loci < 2) send(CRASH);

            temp_map = allocate_map(total);
            parray(placements, total, PLACE);
            for (i = 0; i < *num_unplaced; i++) {
                unplaced[i]->best_pos = (-1);
                unplaced[i]->num_places = 0;
                unplaced[i]->best_map->num_loci = 0;
            }

            print("Start:   ");
            for (j = 0; j < placed->num_loci; j++) {
                sf_locnum(placed->locus[j], FALSE);
                pr();
                if (j != placed->num_loci - 1) print(" ");
            }
            nl();

            placed_any = FALSE;
            while (*num_unplaced > 0)
                if (!extend_order_by_1(placed, unplaced, num_unplaced, npt_thresh,
                                       &contradiction, placements, temp_map)) {
                    sprintf(ps, "No unique placements for %d %smarker%s\n",
                            *num_unplaced, (placed_any ? "remaining " : ""),
                            maybe_s(*num_unplaced));
                    pr();
                    if (contradiction == TRUE) {
                        sprintf(ps, "Three-point analysis allows no %s.\n",
                                (placed_any ? "further ordering" : "orders at all"));
                        pr();
                    } else if (contradiction == MAYBE) {
                        print("Maximum number of loci in one map reached.\n");
                    }
                    break;
                } else placed_any = TRUE;
            if (*num_unplaced == 0) {
                sprintf(ps, "Uniquely ordered all %d markers\n\n", total);
                pr();
            }

            if (print_anyway || placed_any) {
                init_not_fixed(placed);
                converge_to_map(placed);
                nl();
                print_long_map(placed, "Map:");
                if (*num_unplaced > 0) {
                    nl();
                    print("Markers placed relative to above map:\n");
                    new_print_placements(placed, unplaced, *num_unplaced);
                }
            }
        } on_exit {
        unparray(placements, total, PLACE);
        free_map(temp_map);
        relay_messages;
    }
}


bool extend_order_by_1(MAP *placed, PLACEME **unplaced, int *num_unplaced, real npt_thresh,
                       bool *contradiction, /* iff return FALSE, is TRUE if 3pt excludes all */
                       PLACE **placements, /* these are all temps */
                       MAP *temp_map) {
    int places, best_places, best_i, i, j;
    int left, right, how;
    real next_best, best_unallowed, error_lod;
    bool zero_place, off_ends, single_error;

    if (placed->num_loci >= MAX_MAP_LOCI) {
        *contradiction = MAYBE;
        return (FALSE);
    }
    randomize_markers_to_try(unplaced, *num_unplaced);
    for (i = 0; i < *num_unplaced; i++) unplaced[i]->best_map->num_loci = 0;

    /* Try to place remaining loci, into order, in their order */
    if (use_three_pt) {
        best_places = NO_PLACES;
        best_i = -1;
        for (i = 0; i < *num_unplaced; i++) { /* sorted */
            unplaced[i]->num_places = places =
                    three_pt_exclusions(placed->locus, placed->num_loci,
                                        unplaced[i]->locus, unplaced[i]->excluded);
            unplaced[i]->best_pos = (-1); /* NO_POS? */
            if (places > 0 && places < best_places) {
                best_places = places;
                best_i = i;
                if (places == 1 && !middles_excluded(unplaced[i]->excluded,
                                                     placed->num_loci)) {
                    /* It doesn't get any better than this! */
                    add_to_order(BY_3PT, placed->locus, &placed->num_loci,
                                 unplaced, best_i, num_unplaced);
                    return (TRUE);
                }
            }
        } /* for i: loop over num_to_place trying 3pt */
        if (best_places == NO_PLACES) {
            *contradiction = TRUE;
            return (FALSE);
        }

        /**** We get to here if all allowed positions are off-end, or if no
              unique placements were found ****/
#ifdef DONT_DO_THIS
        for (i=0; i<*num_unplaced; i++) if (unplaced[i]->num_places==1) {
            /* a troublesome however unique placement */
            add_to_order(BY_3PT_ENDS,placed->locus,&placed->num_loci,unplaced,
                 i,num_unplaced);
            return(TRUE);
        }
#endif
        /* otherwise fall-through to Npt tests below */

    } else { /* !use_three-pt */
        for (i = 0; i < *num_unplaced; i++) {
            /* unplaced[i]->num_places= placed->num_loci+1; */
            for (j = 0; j < placed->num_loci; j++) unplaced[i]->excluded[j] = FALSE;
        }
    }

    /**** Have No UNIQUE Placements, so try N-pt analysis ****/
    for (i = 0; i < placed->num_loci - 1; i++)
        placed->rec_frac[i][MALE] = placed->rec_frac[i][FEMALE] = NOT_FIXED;
    init_rec_fracs(placed);
    converge_to_map(placed);
    rank_markers_to_try(unplaced, *num_unplaced);

    best_places = NO_PLACES;
    best_i = -1;
    for (i = 0; i < *num_unplaced; i++) { /* sorted */
        keep_user_amused("marker", i + 1, *num_unplaced);
        if (unplaced[i]->num_places == 0) continue;
        find_window(placed->locus, placed->num_loci,
                    unplaced[i]->locus, unplaced[i]->excluded,
                    npt_window, &left, &right);
        place_locus(placed, unplaced[i]->locus, left, right,
                    unplaced[i]->excluded, placements, &unplaced[i]->best_pos,
                    unplaced[i]->best_map, temp_map);
        unplaced[i]->num_places = places =
                npt_exclusions(placements, placed->num_loci, npt_thresh,
                               error_single_thresh, error_net_thresh, /* GLOBALS */
                               unplaced[i]->excluded, &zero_place, &off_ends,
                               &error_lod, &single_error, &next_best, &best_unallowed);
        unplaced[i]->off_end = off_ends;
        unplaced[i]->errors = (error_lod != NO_ERRORS);

        /* preferentially add loci which place well */
        if (places < best_places && !unplaced[i]->off_end &&
            !unplaced[i]->errors) {
            best_places = places;
            best_i = i;
            if (places == 1) { /* Greedy: It doesn't get any better than this! */
                add_to_order(i + 1, placed->locus, &placed->num_loci, unplaced,
                             best_i, num_unplaced);
                return (TRUE);
            }
        }
    } /* for i: loop over num_to_place trying 3pt */

    /**** We get to here if all placements are troublesome, or if no unique
          placements were found ****/
    for (i = 0; i < *num_unplaced; i++)
        if (unplaced[i]->num_places == 1) {
            /* a troublesome however unique placement */
            how = (unplaced[i]->off_end ? BY_NPT_ENDS : BY_NPT_ERRORS);
            add_to_order(how, placed->locus, &placed->num_loci, unplaced, i,
                         num_unplaced);
            return (TRUE);
        }

    /**** Else we have no unique placements, we fail ****/
    *contradiction = FALSE;
    return (FALSE);
}


void randomize_markers_to_try(PLACEME **unplaced, int num_unplaced) {
    int i, n_dom, n_het, num;

    for (i = 0; i < num_unplaced; i++) {
        unplaced[i]->num_places = 0;
        f2_genotype(unplaced[i]->locus, use_haplotypes, observations);
        num = f2_count_infs(&n_dom, &n_het, observations);
        unplaced[i]->priority =
                imaxf(num - 2 * n_dom, 0) +
                irand(raw.data.f2.num_indivs / 10);
    }
    qsort(unplaced, num_unplaced, sizeof(PLACEME *), compare_markers_to_try);
}


void rank_markers_to_try(PLACEME **unplaced, int num_unplaced) {
    qsort(unplaced, num_unplaced, sizeof(PLACEME *), compare_markers_to_try);
}


int compare_markers_to_try(const void *a, const void *b) {
    PLACEME *x = *(PLACEME *const *) a;
    PLACEME *y = *(PLACEME *const *) b;

    /* num_places is three_pt criteria */
    if (x->num_places < y->num_places) return (-1);
    else if (x->num_places > y->num_places) return (1);
    else if (x->priority > y->priority) return (-1);
    else if (x->priority < y->priority) return (1);
    else return (0);
}


void add_to_order(int how, int *order, int *num_ordered, PLACEME **unplaced, int index, int *num_unplaced) {
    int i, pos, last, new_marker;
    PLACEME *temp;

    for (pos = 0; pos <= *num_ordered; pos++)
        if (!unplaced[index]->excluded[pos]) break;
    if (pos == *num_ordered + 1) send(CRASH);

    new_marker = unplaced[index]->locus;
    if (pos < *num_ordered)
        for (i = *num_ordered - 1; i >= pos; i--) order[i + 1] = order[i];
    order[pos] = new_marker;
    *num_ordered += 1;

    last = *num_unplaced - 1;
    if (*num_unplaced == 1) *num_unplaced = 0;
    else if (index == last) *num_unplaced -= 1;
    else {
        temp = unplaced[index];
        unplaced[index] = unplaced[last];
        unplaced[last] = temp;
        *num_unplaced -= 1;
    }

#ifndef DEBUGGING
    if (how > 0) {
        sprintf(ps, "Npt-%d:\t ", how);
        pr();
    }
#else
        if (how>0)                        { sprintf(ps,"Npt:\t "); pr(); }
#endif
    else if (how == BY_3PT) print("3pt:     ");
    else if (how == BY_NPT_ENDS) print("Npt-End: ");
    else if (how == BY_NPT_ERRORS) print("Npt-Err: ");
    else if (how == BY_3PT_ENDS) print("3pt-End: ");
    else
        send(CRASH);
    to_column(9);

    for (i = 0; i < *num_ordered; i++) { /* print new order */
        sf_locnum(order[i], (order[i] == new_marker));
        pr();
        if (i != *num_ordered - 1) print(" ");
    }
    nl();
}


/* to keep this all straight:
   orig frame:         a - b - c - d - e - f - g
   orig frame index:   0 - 1 - 2 - 3 - 4 - 5 - 6    num_placed=7
   exclusion index:  0 - 1 - 2 - 3 - 4 - 5 - 6 - 7 
   new frame:              b - c - d - e - f    
   new frame index:        0 - 1 - 2 - 3 - 4        */

void find_window(int *order, /*NOTUSED*/ int num_placed, int new_locus, int *excluded,
                 int min_window,    /* #of loci to left and right, should be odd, e.g. 5 */
                 int *start, int *end  /* side-effected with interval indecies */) {
    int j, step;
    /* simple implementation for now - this should later use informativeness */

    step = max((min_window - 1) / 2, 1);
    for (j = 0; j <= num_placed && excluded[j]; j++) {}
    *start = max(j - step, 0);
    for (j = num_placed; j > 0 && excluded[j]; j--) {}
    *end = min(j + step - 1, num_placed - 1);   /* Is this right? */
}


bool middles_excluded(int *excluded, int num_order) {
    int j;

    for (j = 1; j < num_order; j++) if (!excluded[j]) return (FALSE);
    return (TRUE);
}
