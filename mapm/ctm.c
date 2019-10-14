/******************************************************************************

  ####    #####  #    #           ####
 #    #     #    ##  ##          #    #
 #          #    # ## #          #
 #          #    #    #   ###    #
 #    #     #    #    #   ###    #    #
  ####      #    #    #   ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

/************** Hidden Markov Model based Converge to Map for F2 **************/

#include "mapm.h"

int *observations = NULL; /* [indiv] externed for use as a temp */
int obligate_recs[6][6] = { /* externed in toplevel.h */
        /* -  A  H  B  C  D */
        {0, 0, 0, 0, 0, 0}, /* - */
        {0, 0, 1, 2, 1, 0}, /* A */
        {0, 1, 0, 1, 0, 0}, /* H */
        {0, 2, 1, 0, 0, 1}, /* B */
        {0, 1, 0, 0, 0, 0}, /* C */
        {0, 0, 1, 0, 0, 0}  /* D */
};

typedef struct {
    real recs;   /* for now */
    real norecs;
} SUFF_STATS;

typedef struct {
    int n_indivs, n_loci, n_intervals;
    int max_observations, max_states;
    int **n_poss_states, ***poss_state;     /* [indiv][locus],[indiv][locus][#] */
    real ***trans_prob, ***exp_transitions;  /* [interval][from_state][to_state] */
    real ***implied_recs, ***implied_norecs; /* [interval][from_state][to_state] */
    real ***obs_prob, ***exp_observations;   /* [locus][observation/MD][state] */
    //bool *obs_prob_fixed;  /* [locus] - if all are fixed then exp_obs==NULL */
    int **observation;    /* [indiv][locus] => OBS_ code */
    int *the_observation; /* [indiv] a temp used only in f2_setup_hmm */
    real **apriori_prob;   /* [indiv][state] - for leftmost locus */
    int no_data;          /* index [indiv] into apriori_prob for right locus */
    real *e_step_temp;     /* [state]- used in LCP and RCP calculations */
    /*real ***exp_genotype;*/  /* [indiv][locus][state] for ERROR-KLUDGE */
    int *observations; /* [indiv] externed for use as a temp */

    void (*hmm_set_probs)();

    //bool error_kludge; /* TRUE for ERROR-KLUDGE, e.g. error_rate>0 */
    bool separate_recs;
    /* TRUE if distinct implied_rec/norec entries for each interval */
    SUFF_STATS **sufficient_stats; /* [locus] => ptr to struct */

    /* Internals for HMM E-Step (obs stuff moved to toplevel,h) */
    real **left_prob, **right_prob;  /* [locus#][state#] */
    real ***indiv_transitions;       /* [interval][from_state][to_state] */
    //real ***indiv_observations;      /* [locus][observation][state], maybe NULL */
    int ri_type;  /* RI_SELF, RI_SIB, or RI_NOT */

    /* There are 4 F2 states and 16 hidden F2/F3 states: */
    /*
    int A_AMAM, A_AMAP, A_APAM, A_APAP, H_AMAM, H_AMBP, H_BPAM, H_BPBP;
    int J_BMBM, J_BMAP, J_APBM, J_APAP, B_BMBM, B_BMBP, B_BPBM, B_BPBP;
    */
    /* These are temps for setup_f3_self_hmm(), although f2_recs is also
       used in hmm_f3_self_set_probs() */
    int f3_states, mat_chrom[16], pat_chrom[16], f2_state[16];
    /*real f2_recs[4][4];*/
    real null_like;
} ctm_state_t;

static ctm_state_t *ctm_states = NULL;
static int n_ctm_states = 0;
static MAP *map2 = NULL;

#define HMM_MAX_THETA 0.4995
#define HMM_MIN_THETA 0.0004

/* F2 backcross/intercross, RI states */
#define STATE_A 0
#define STATE_B 1
#define STATE_H 2   /* intercross only */

/**** F3 specific stuff follows: ****/

#define A 0  /* F2 States... */
#define H 1
#define J 2
#define B 3

/* States for one F3 chromosome (grand- maternal or paternal origin) */
#define M 0
#define P 1

/**** The internal routines ****/

static int hmm_converge_to_map(ctm_state_t *ctm_state, MAP *map);

static int f3_state(ctm_state_t *ctm_state,
                    const char *name /* for debugging */,
                    int f2, int mat, int pat, int obs);

static int hmm_e_step(ctm_state_t *ctm_state, real *new_likep);

static void hmm_count_stats(ctm_state_t *ctm_state);

static void hmm_make_new_map(const ctm_state_t *ctm_state, real new_like, MAP *map);

static void setup_hmm(ctm_state_t *ctm_state, MAP *map);

static void setup_bc_like_hmm(ctm_state_t *ctm_state,
                              int *locus /* The locus numbers in order */,
        /*double *error_rate,*/ int cross_type);

static void hmm_bc_like_set_probs(ctm_state_t *ctm_state, MAP *map);

static void setup_f2_hmm(ctm_state_t *ctm_state,
                         int *locus /* The locus numbers in order */ /*,
				       double *error_rate*/);

static void hmm_f2_set_probs(ctm_state_t *ctm_state, MAP *map);

static void setup_f3_self_hmm(ctm_state_t *ctm_state,
                              int *locus /* The locus numbers in the sequence */);

static void hmm_f3_self_set_probs(ctm_state_t *ctm_state, MAP *map);

#if 0
static void count_error_lods(const ctm_state_t *ctm_state, MAP *map);
#endif


/**** External Routines, mostly ****/

void allocate_hmm_temps(int total_loci, int num_indivs, int cross_type) {
    int i, j, num_states = 0, num_observations = 0, num_loci, num_intervals;
    bool special_recs_hack;

#ifdef THREAD
    n_ctm_states = threads_nthreads + 1; /* The first state is always available for non-threaded purposes */
#else
    n_ctm_states = 1;
#endif

    ctm_states = malloc(sizeof(ctm_state_t) * n_ctm_states);

    for (j = 0; j < n_ctm_states; j++) {
//    printf("Setting up HMM temp space for worker thread %d\n", j);
        switch (cross_type) {
            case F2_BACKCROSS:
                num_states = 2;
                num_observations = 3;
                break;
            case F2_INTERCROSS:
                num_states = 3;
                num_observations = 6;
                break;
            case RI_SELF:
                num_states = 2;
                num_observations = 3;
                break;
            case RI_SIB:
                num_states = 2;
                num_observations = 3;
                break;
            case F3_SELF:
                num_states = 16;
                num_observations = 6;
                break;
            default:
                send(CRASH);
        }
        num_loci = MAX_MAP_LOCI;
        num_intervals = num_loci - 1;
        special_recs_hack = (cross_type == F2_INTERCROSS ? TRUE : FALSE);

        parray(ctm_states[j].sufficient_stats, num_loci, SUFF_STATS);
        matrix(ctm_states[j].left_prob, num_loci, num_states, real);
        matrix(ctm_states[j].right_prob, num_loci, num_states, real);
        //ctm_states[j].indiv_observations=NULL;

        array(ctm_states[j].trans_prob, num_intervals, real**);
        array(ctm_states[j].exp_transitions, num_intervals, real**);
        array(ctm_states[j].indiv_transitions, num_intervals, real**);
        array(ctm_states[j].implied_recs, num_intervals, real**);
        array(ctm_states[j].implied_norecs, num_intervals, real**);
        for (i = 0; i < num_intervals; i++) {
            matrix(ctm_states[j].trans_prob[i], num_states, num_states, real);
            matrix(ctm_states[j].exp_transitions[i], num_states, num_states, real);
            matrix(ctm_states[j].indiv_transitions[i], num_states, num_states, real);
            if (special_recs_hack || i == 0) {
                matrix(ctm_states[j].implied_recs[i], num_states, num_states, real);
                matrix(ctm_states[j].implied_norecs[i], num_states, num_states, real);
            }
        }
        if (!special_recs_hack)
            for (i = 1; i < num_intervals; i++) {
                /*  Do the funky pointer KLUDGE */
                ctm_states[j].implied_recs[i] = ctm_states[j].implied_recs[0];
                ctm_states[j].implied_norecs[i] = ctm_states[j].implied_norecs[0];
            }

        //array(ctm_states[j].obs_prob_fixed,num_loci,bool);
        array(ctm_states[j].obs_prob, num_loci, real**);
        for (i = 0; i < num_loci; i++) {
            matrix(ctm_states[j].obs_prob[i], num_observations, num_states, real);
        }

        matrix(ctm_states[j].n_poss_states, num_indivs, num_loci, int);
        array(ctm_states[j].poss_state, num_indivs, int**);
        for (i = 0; i < num_indivs; i++) {
            matrix(ctm_states[j].poss_state[i], num_loci, num_states, int);
        }

        ctm_states[j].exp_observations = NULL;  /* unused at present */

        matrix(ctm_states[j].observation, num_indivs, num_loci, int);
        array(ctm_states[j].the_observation, num_indivs, int);
        matrix(ctm_states[j].apriori_prob, num_indivs + 1, num_states, real);
        array(ctm_states[j].e_step_temp, num_states, real);

#if 0
        array(ctm_states[j].exp_genotype,num_indivs,real**);
        for (i=0; i<num_indivs; i++) {
          matrix(ctm_states[j].exp_genotype[i],num_loci,num_states,real);
        }
#endif
        array(ctm_states[j].observations, num_indivs, int);
    }

    map2 = allocate_map(2); /* should have NULL error_lod stuff */
    array(observations, num_indivs, int); /* extern */
}


void free_hmm_temps(int total_loci, int num_indivs, int cross_type) {
    int i, j, num_states = 0, num_observations = 0, num_loci = 0, num_intervals;
    bool special_recs_hack;

    for (j = 0; j < n_ctm_states; j++) {
        switch (cross_type) {
            case F2_BACKCROSS:
                num_states = 2;
                num_observations = 3;
                break;
            case F2_INTERCROSS:
                num_states = 3;
                num_observations = 6;
                break;
            case RI_SELF:
                num_states = 2;
                num_observations = 3;
                break;
            case RI_SIB:
                num_states = 2;
                num_observations = 3;
                break;
            case F3_SELF:
                num_states = 16;
                num_observations = 6;
                break;
            default:
                send(CRASH);
        }
        num_loci = MAX_MAP_LOCI;
        num_intervals = num_loci - 1;
        special_recs_hack = (cross_type == F2_INTERCROSS ? TRUE : FALSE);

        unparray(ctm_states[j].sufficient_stats, num_loci, SUFF_STATS);
        unmatrix(ctm_states[j].left_prob, num_loci, real);
        unmatrix(ctm_states[j].right_prob, num_loci, real);

        for (i = 0; i < num_intervals; i++) {
            unmatrix(ctm_states[j].trans_prob[i], num_states, real);
            unmatrix(ctm_states[j].exp_transitions[i], num_states, real);
            unmatrix(ctm_states[j].indiv_transitions[i], num_states, real);
            if (special_recs_hack || i == 0) {
                unmatrix(ctm_states[j].implied_recs[i], num_states, real);
                unmatrix(ctm_states[j].implied_norecs[i], num_states, real);
            }
        }
        unarray(ctm_states[j].trans_prob, real**);
        unarray(ctm_states[j].exp_transitions, real**);
        unarray(ctm_states[j].indiv_transitions, real);
        unarray(ctm_states[j].implied_recs, real**);
        unarray(ctm_states[j].implied_norecs, real**);
        /* indiv_observations was never allocated */

        //unarray(ctm_states[j].obs_prob_fixed,bool);
        for (i = 0; i < num_loci; i++) {
            unmatrix(ctm_states[j].obs_prob[i], num_observations, real);
        }
        unarray(ctm_states[j].obs_prob, real**);

        unmatrix(ctm_states[j].n_poss_states, num_indivs, int);
        for (i = 0; i < num_indivs; i++) {
            unmatrix(ctm_states[j].poss_state[i], num_loci, int);
        }
        unarray(ctm_states[j].poss_state, int);

        /* exp_observations was never allocated */

        unmatrix(ctm_states[j].observation, num_indivs, int);
        unarray(ctm_states[j].the_observation, int);
        unmatrix(ctm_states[j].apriori_prob, num_indivs + 1, real);
        unarray(ctm_states[j].e_step_temp, real);

#if 0
        for (i=0; i<num_indivs; i++) {
          unmatrix(ctm_states[j].exp_genotype[i],num_loci,real);
        }
        unarray(ctm_states[j].exp_genotype,real**);
#endif
        unarray(ctm_states[j].observations, int); /* extern */
    }

    free(ctm_states);
    free_map(map2); /* should have NULL error_lod stuff */
    unarray(observations, int); /* extern */
}


int converge_to_map_mt(MAP *map, int index) {
    if (map == NULL || map->num_loci < 2 || raw.data_type != F2) send(CRASH);
    if (map->num_loci > MAX_MAP_LOCI) {
        // Not thread safe...
        strcpy(ps, "too many loci in one map");
        error(ps);
    }

    setup_hmm(&ctm_states[index], map); /* set statics above, and set the initial parameters */
    return hmm_converge_to_map(&ctm_states[index], map);
}

int converge_to_map(MAP *map) {
    return converge_to_map_mt(map, 0);
}


void f2_genotype(int locus, bool haplo /* merge haplotypes */,
                 int *observation /* array of raw.num_indivs observations [raw.num_indivs] */) {
    int j;

    if (raw.data_type != F2) send(CRASH);

    if (haplo && haplotyped(locus))
        locus = haplo_first[locus];  /* go to the head of the list */

    for (j = 0; j < raw.data.f2.num_indivs; j++)
        switch (raw.data.f2.allele[locus][j]) {
            case PARENTAL_TYPE_A:
                observation[j] = OBS_A;
                break;
            case PARENTAL_TYPE_B:
                observation[j] = OBS_B;
                break;
            case HYBRID_TYPE_H:
                observation[j] = OBS_H;
                break;
            case TYPE_NOT_A:
                observation[j] = OBS_NOT_A;
                break;
            case TYPE_NOT_B:
                observation[j] = OBS_NOT_B;
                break;
            default:
                observation[j] = OBS_MISSING;
                break;
        }

    if (haplo && haplo_first[locus] != NO_LOCUS)
        while ((locus = haplo_next[locus]) != NO_LOCUS)
            merge_genotypes(locus, observation, observation);
}


bool merge_genotypes(int locus,
                     const int *observation /* array of raw.num_indivs observations [raw.num_indivs] */,
                     int *new_observation /* array of raw.num_indivs observations [raw.num_indivs] */) {
    int j, the_obs;

    for (j = 0; j < raw.data.f2.num_indivs; j++) {
        switch (raw.data.f2.allele[locus][j]) {
            case PARENTAL_TYPE_A:
                the_obs = OBS_A;
                break;
            case PARENTAL_TYPE_B:
                the_obs = OBS_B;
                break;
            case HYBRID_TYPE_H:
                the_obs = OBS_H;
                break;
            case TYPE_NOT_A:
                the_obs = OBS_NOT_A;
                break;
            case TYPE_NOT_B:
                the_obs = OBS_NOT_B;
                break;
            default:
                the_obs = OBS_MISSING;
                break;
        }
        if (obligate_recs[observation[j]][the_obs])
            return (FALSE);
        if ((observation[j] == OBS_MISSING && the_obs != OBS_MISSING) ||
            (!fully_inf(observation[j]) && fully_inf(the_obs)))
            new_observation[j] = the_obs;
        else new_observation[j] = observation[j];
    }
    return (TRUE);
}


int f2_count_infs(int *num_dom, int *num_het, /* side-effected */
                  const int *observation /* array of raw.num_indivs observations [raw.num_indivs] */) /* return #inf */
{
    int homo, het, dom, missing, j;

    homo = het = dom = missing = 0;
    for (j = 0; j < raw.data.f2.num_indivs; j++)
        switch (observation[j]) {
            case OBS_A:
                homo++;
                break;
            case OBS_B:
                homo++;
                break;
            case OBS_H:
                het++;
                break;
            case OBS_NOT_A:
                dom++;
                break;
            case OBS_NOT_B:
                dom++;
                break;
            default:
                missing++;
                break;
        }
    *num_dom = dom;
    *num_het = het;
    return raw.data.f2.num_indivs - missing;
}


/**** The real Mapmaker itself ****/

static int hmm_converge_to_map(ctm_state_t *ctm_state, MAP *map) {
    int iter = 0;
    real old_like = 0.0, new_like = 0.0;

    for (iter = 0; TRUE; ++iter) {
        (*ctm_state->hmm_set_probs)(ctm_state, map);
        old_like = new_like;
        if (!hmm_e_step(ctm_state, &new_like)) {
            return 0;
        }
        if (ctm_state->n_loci == 2 && iter == 0) ctm_state->null_like = new_like;

        if (fabs(new_like - old_like) < tolerance) break; /* converged! */
        hmm_count_stats(ctm_state);

        hmm_make_new_map(ctm_state, new_like, map);

    } /* iterate */

    return 1;
}


static void setup_hmm(ctm_state_t *ctm_state, MAP *map) {
    int type;

    ctm_state->n_loci = map->num_loci;
    ctm_state->n_intervals = ctm_state->n_loci - 1;
    ctm_state->n_indivs = raw.data.f2.num_indivs;
    ctm_state->no_data = raw.data.f2.num_indivs; /* extra index into apriori_prob */
    //ctm_state->error_kludge=FALSE;
    ctm_state->ri_type = RI_NOT;
    type = raw.data.f2.cross_type;

    switch (type) {
        case F2_INTERCROSS:
            setup_f2_hmm(ctm_state, map->locus/*, map->error_rate*/);
            break;
        case F3_SELF:
            setup_f3_self_hmm(ctm_state, map->locus);
            break;
        case F2_BACKCROSS:
        case RI_SELF:
        case RI_SIB:
            setup_bc_like_hmm(ctm_state, map->locus/*, map->error_rate*/, type);
            break;
        default:
            send(CRASH);
            break;
    }
}


static void
setup_bc_like_hmm(ctm_state_t *ctm_state, int *locus /* The locus numbers in order */, /*double *error_rate,*/ int cross_type) {
    int i, j, k, obs;
    real **recs, **norecs/*, error_prob*/;

    ctm_state->max_states = 2;    /* STATE_A or _B */
    ctm_state->max_observations = 3; /* OBS_A, _H, or _MISSING (yup, this is odd for BC) */
    ctm_state->separate_recs = FALSE;
    ctm_state->hmm_set_probs = hmm_bc_like_set_probs;
    //if (error_rate!=NULL) ctm_state->error_kludge=TRUE;

    switch (cross_type) {
        case F2_BACKCROSS:
            ctm_state->ri_type = RI_NOT;
            break;
        case RI_SELF:
            ctm_state->ri_type = RI_SELF;
            break;
        case RI_SIB:
            ctm_state->ri_type = RI_SIB;
            break;
        default:
            send(CRASH);
            break;
    }

    /*error_prob=0.0;*/
    for (i = 0; i < ctm_state->n_loci; i++) {
        //if (ctm_state->error_kludge) error_prob=error_rate[i]; /* else =0.0 */
        ctm_state->obs_prob[i][OBS_A][STATE_A] = 1.0/* - error_prob*/;
        ctm_state->obs_prob[i][OBS_H][STATE_A] = 0.0/* + error_prob*/;
        ctm_state->obs_prob[i][OBS_A][STATE_B] = 0.0/* + error_prob*/;
        ctm_state->obs_prob[i][OBS_H][STATE_B] = 1.0/* - error_prob*/;
        ctm_state->obs_prob[i][OBS_MISSING][STATE_A] = 1.0;
        ctm_state->obs_prob[i][OBS_MISSING][STATE_B] = 1.0;
        //ctm_state->obs_prob_fixed[i]=TRUE;

        f2_genotype(locus[i], use_haplotypes, ctm_state->the_observation);
        for (j = 0; j < ctm_state->n_indivs; j++) {
            obs = ctm_state->observation[j][i] = ctm_state->the_observation[j];
            k = 0;
            if (obs == OBS_MISSING/* || error_prob>0.0*/) {
                ctm_state->poss_state[j][i][k++] = STATE_A;
                ctm_state->poss_state[j][i][k++] = STATE_B;
            } else if (obs == OBS_A) {
                ctm_state->poss_state[j][i][k++] = STATE_A;
            } else if (obs == OBS_H) {
                ctm_state->poss_state[j][i][k++] = STATE_B;
            } else
                send(CRASH);
            ctm_state->n_poss_states[j][i] = k;
        }
    }

    for (i = 0; i < ctm_state->n_intervals; i++) {
        /* set the implied_recs */
        recs = ctm_state->implied_recs[i];
        norecs = ctm_state->implied_norecs[i];
        recs[STATE_A][STATE_A] = 0.0;
        norecs[STATE_A][STATE_A] = 1.0;
        recs[STATE_A][STATE_B] = 1.0;
        norecs[STATE_A][STATE_B] = 0.0;
        recs[STATE_B][STATE_A] = 1.0;
        norecs[STATE_B][STATE_A] = 0.0;
        recs[STATE_B][STATE_B] = 0.0;
        norecs[STATE_B][STATE_B] = 1.0;
    } /* end loop over intervals (i) */

    for (j = 0; j < ctm_state->n_indivs; j++) {
        obs = ctm_state->observation[j][0];
        if (obs != OBS_MISSING) {
            ctm_state->apriori_prob[j][STATE_A] = 0.50 * ctm_state->obs_prob[0][obs][STATE_A];
            ctm_state->apriori_prob[j][STATE_B] = 0.50 * ctm_state->obs_prob[0][obs][STATE_B];
        } else {
            ctm_state->apriori_prob[j][STATE_A] = 0.50;
            ctm_state->apriori_prob[j][STATE_B] = 0.50;
        }
    }
    ctm_state->apriori_prob[ctm_state->no_data][STATE_A] = 0.50;
    ctm_state->apriori_prob[ctm_state->no_data][STATE_B] = 0.50;
}


static void setup_f2_hmm(ctm_state_t *ctm_state, int *locus /* The locus numbers in order */  /*, double *error_rate*/) {
    int i, j, k, obs;
    real total, **recs, **norecs/*, error_prob*/;

    ctm_state->max_states = 3;    /* STATE_A,_H, or _B */
    ctm_state->max_observations = 6; /* OBS_A, _B, _H, _NOT_A, _NOT_B, or _MISSING */
    ctm_state->separate_recs = TRUE;
    ctm_state->hmm_set_probs = hmm_f2_set_probs;
    //if (error_rate!=NULL) ctm_state->error_kludge=TRUE;

    /*error_prob=0.0;*/
    for (i = 0; i < ctm_state->n_loci; i++) {
        //if (ctm_state->error_kludge) error_prob= error_rate[i]; /* else =0.0 */
        ctm_state->obs_prob[i][OBS_A][STATE_A] = 1.0/* - error_prob*/;
        ctm_state->obs_prob[i][OBS_B][STATE_A] = 0.0/* + error_prob/2.0*/;
        ctm_state->obs_prob[i][OBS_H][STATE_A] = 0.0/* + error_prob/2.0*/;

        /* Fix the NOT_A/B cases */
        ctm_state->obs_prob[i][OBS_NOT_A][STATE_A] = 0.0/* + (2.0 * error_prob)*/; /* KLUDGE */
        ctm_state->obs_prob[i][OBS_NOT_B][STATE_A] = 1.0/* - error_prob*/;

        ctm_state->obs_prob[i][OBS_A][STATE_H] = 0.0/* + error_prob/2.0*/;
        ctm_state->obs_prob[i][OBS_B][STATE_H] = 0.0/* + error_prob/2.0*/;
        ctm_state->obs_prob[i][OBS_H][STATE_H] = 1.0/* - error_prob*/;
        ctm_state->obs_prob[i][OBS_NOT_A][STATE_H] = 1.0/* - error_prob*/;
        ctm_state->obs_prob[i][OBS_NOT_B][STATE_H] = 1.0/* - error_prob*/;

        ctm_state->obs_prob[i][OBS_A][STATE_B] = 0.0/* + error_prob/2.0*/;
        ctm_state->obs_prob[i][OBS_B][STATE_B] = 1.0/* - error_prob*/;
        ctm_state->obs_prob[i][OBS_H][STATE_B] = 0.0/* + error_prob/2.0*/;
        ctm_state->obs_prob[i][OBS_NOT_A][STATE_B] = 1.0/* - (2.0 * error_prob)*/;
        ctm_state->obs_prob[i][OBS_NOT_B][STATE_B] = 0.0/* + error_prob*/;

        ctm_state->obs_prob[i][OBS_MISSING][STATE_A] = 1.0;
        ctm_state->obs_prob[i][OBS_MISSING][STATE_B] = 1.0;
        ctm_state->obs_prob[i][OBS_MISSING][STATE_H] = 1.0;

        //ctm_state->obs_prob_fixed[i]=TRUE;

        /* make a copy of the raw data, translated for speed. Now also
           adds the haplotype info, when available - note order
           of indices for the observation matrix is opposite that
           for the raw struct! */

        f2_genotype(locus[i], use_haplotypes, ctm_state->the_observation);
        for (j = 0; j < ctm_state->n_indivs; j++) {
            obs = ctm_state->observation[j][i] = ctm_state->the_observation[j];
            k = 0;
            if (obs == OBS_MISSING/* || error_prob>0.0*/) {
                ctm_state->poss_state[j][i][k++] = STATE_A;
                ctm_state->poss_state[j][i][k++] = STATE_B;
                ctm_state->poss_state[j][i][k++] = STATE_H;
            } else if (obs == OBS_NOT_A) {
                ctm_state->poss_state[j][i][k++] = STATE_B;
                ctm_state->poss_state[j][i][k++] = STATE_H;
            } else if (obs == OBS_NOT_B) {
                ctm_state->poss_state[j][i][k++] = STATE_A;
                ctm_state->poss_state[j][i][k++] = STATE_H;
            } else if (obs == OBS_A) {
                ctm_state->poss_state[j][i][k++] = STATE_A;
            } else if (obs == OBS_B) {
                ctm_state->poss_state[j][i][k++] = STATE_B;
            } else if (obs == OBS_H) {
                ctm_state->poss_state[j][i][k++] = STATE_H;
            } else
                send(CRASH);
            ctm_state->n_poss_states[j][i] = k;
        }
    }

    for (i = 0; i < ctm_state->n_intervals; i++) {
        /* set most of the implied_recs, only the [H][H] entries vary */
        recs = ctm_state->implied_recs[i];
        norecs = ctm_state->implied_norecs[i];
        recs[STATE_A][STATE_A] = 0.0;
        norecs[STATE_A][STATE_A] = 2.0;
        recs[STATE_A][STATE_B] = 2.0;
        norecs[STATE_A][STATE_B] = 0.0;
        recs[STATE_A][STATE_H] = 1.0;
        norecs[STATE_A][STATE_H] = 1.0;
        recs[STATE_B][STATE_A] = 2.0;
        norecs[STATE_B][STATE_A] = 0.0;
        recs[STATE_B][STATE_B] = 0.0;
        norecs[STATE_B][STATE_B] = 2.0;
        recs[STATE_B][STATE_H] = 1.0;
        norecs[STATE_B][STATE_H] = 1.0;
        recs[STATE_H][STATE_A] = 1.0;
        norecs[STATE_H][STATE_A] = 1.0;
        recs[STATE_H][STATE_B] = 1.0;
        norecs[STATE_H][STATE_B] = 1.0;
        recs[STATE_H][STATE_H] = 0.0;
        norecs[STATE_H][STATE_H] = 0.0;
        /* _H _H entries are a function of theta[i] and are changed
           each iteration by set_probs() */
    } /* end loop over intervals (i) */

    /* Lastly, we fill in apriori state probs for leftmost locus only:
       THIS HAS BEEN CHANGED TO CALCULATE THE UN-NORMALIZED PRIORS! */

    for (j = 0; j < ctm_state->n_indivs; j++) {
        obs = ctm_state->observation[j][0];
        total = 0.0;
        if (obs != OBS_MISSING) { /* NOT MISSING  - general version??? */

#ifdef HMMMM
            total+= ctm_state->apriori_prob[j][STATE_A]= ctm_state->obs_prob[0][obs][STATE_A];
            total+= ctm_state->apriori_prob[j][STATE_B]= ctm_state->obs_prob[0][obs][STATE_B];
            total+= ctm_state->apriori_prob[j][STATE_H]= ctm_state->obs_prob[0][obs][STATE_H];
#else
            total += ctm_state->apriori_prob[j][STATE_A] = 0.25 * ctm_state->obs_prob[0][obs][STATE_A];
            total += ctm_state->apriori_prob[j][STATE_B] = 0.25 * ctm_state->obs_prob[0][obs][STATE_B];
            total += ctm_state->apriori_prob[j][STATE_H] = 0.50 * ctm_state->obs_prob[0][obs][STATE_H];
#endif

/*	    ctm_state->apriori_prob[j][STATE_A]/= total;
	    ctm_state->apriori_prob[j][STATE_B]/= total;
	    ctm_state->apriori_prob[j][STATE_H]/= total; */

        } else { /* OBSERVATION IS MISSING */
            ctm_state->apriori_prob[j][STATE_A] = 0.25;
            ctm_state->apriori_prob[j][STATE_B] = 0.25;
            ctm_state->apriori_prob[j][STATE_H] = 0.50;
        }
    }
    ctm_state->apriori_prob[ctm_state->no_data][STATE_A] = 0.25;
    ctm_state->apriori_prob[ctm_state->no_data][STATE_B] = 0.25;
    ctm_state->apriori_prob[ctm_state->no_data][STATE_H] = 0.50;
}


static void setup_f3_self_hmm(ctm_state_t *ctm_state, int *locus /* The locus numbers in the sequence */) {
    int i, j, k, the_obs;
    real recs;
    /*
    int a[3], b[3], h[3], n[3], aa, ab, ah, ba, bb, bh, ha, hb, hh;
    int A_AMAM, A_AMAP, A_APAM, A_APAP, H_AMAM, H_AMBP, H_BPAM, H_BPBP;
    int J_BMBM, J_BMAP, J_APBM, J_APAP, B_BMBM, B_BMBP, B_BPBM, B_BPBP;
    */
    real f2_recs[4][4];

    ctm_state->max_states = 16;
    ctm_state->f3_states = 0;
    ctm_state->max_observations = 6; /* OBS_A, _B, _H, _NOT_A, _NOT_B, or _MISSING */
    ctm_state->separate_recs = FALSE;
    ctm_state->hmm_set_probs = hmm_f3_self_set_probs;

    /*ctm_state->A_AMAM=*/ f3_state(ctm_state, "A-AmAm", A, M, M, OBS_A); /* 0  */
    /*ctm_state->A_AMAP=*/ f3_state(ctm_state, "A-AmAp", A, M, P, OBS_A); /* 1  */
    /*ctm_state->A_APAM=*/ f3_state(ctm_state, "A-ApAm", A, P, M, OBS_A); /* 2  */
    /*ctm_state->A_APAP=*/ f3_state(ctm_state, "A-ApAp", A, P, P, OBS_A); /* 3  */

    /*ctm_state->H_AMAM=*/ f3_state(ctm_state, "H-AmAm", H, M, M, OBS_A); /* 4  */
    /*ctm_state->H_AMBP=*/ f3_state(ctm_state, "H-AmBp", H, M, P, OBS_H); /* 5  */
    /*ctm_state->H_BPAM=*/ f3_state(ctm_state, "H-BpAm", H, P, M, OBS_H); /* 6  */
    /*ctm_state->H_BPBP=*/ f3_state(ctm_state, "H-BpBp", H, P, P, OBS_B); /* 7  */

    /*ctm_state->J_BMBM=*/ f3_state(ctm_state, "J-BmBm", J, M, M, OBS_B); /* 8  */
    /*ctm_state->J_BMAP=*/ f3_state(ctm_state, "J-BmAp", J, M, P, OBS_H); /* 9  */
    /*ctm_state->J_APBM=*/ f3_state(ctm_state, "J-ApBm", J, P, M, OBS_H); /* 10 */
    /*ctm_state->J_APAP=*/ f3_state(ctm_state, "J-ApAp", J, P, P, OBS_A); /* 11 */

    /*ctm_state->B_BMBM=*/ f3_state(ctm_state, "B-BmBm", B, M, M, OBS_B); /* 12 */
    /*ctm_state->B_BMBP=*/ f3_state(ctm_state, "B-BmBp", B, M, P, OBS_B); /* 13 */
    /*ctm_state->B_BPBM=*/ f3_state(ctm_state, "B-BpBm", B, P, M, OBS_B); /* 14 */
    /*ctm_state->B_BPBP=*/ f3_state(ctm_state, "B-BpBp", B, P, P, OBS_B); /* 15 */

    //for (i=0; i<2; i++) a[i]=b[i]=h[i]=n[i]=0;
    //aa=ab=ah=ba=bb=bh=ha=hb=hh=0;
    for (i = 0; i < ctm_state->n_loci; i++) {
        /* Here we have a KLUDGE, in that the poss_state and n_poss_state
           vars could be used to reduce the computation. However, as yet
           we don't */

        //ctm_state->obs_prob_fixed[i]= TRUE;

        /* Make a copy of the raw data, translated for speed
           KLUDGE update to use f2_genotype() */

        for (j = 0; j < ctm_state->n_indivs; j++) {
            ctm_state->n_poss_states[j][i] = ctm_state->f3_states; /* the slow way */
            for (k = 0; k < ctm_state->f3_states; k++) ctm_state->poss_state[j][i][k] = k;
            switch (raw.data.f2.allele[locus[i]][j]) {
                case PARENTAL_TYPE_A:
                    ctm_state->observation[j][i] = OBS_A;
                    break;
                case PARENTAL_TYPE_B:
                    ctm_state->observation[j][i] = OBS_B;
                    break;
                case HYBRID_TYPE_H:
                    ctm_state->observation[j][i] = OBS_H;
                    break;
                case TYPE_NOT_A:
                    ctm_state->observation[j][i] = OBS_NOT_A;
                    break;
                case TYPE_NOT_B:
                    ctm_state->observation[j][i] = OBS_NOT_B;
                    break;
                case MISSING_DATA:
                    ctm_state->observation[j][i] = OBS_MISSING;
                    break;
                default:
                    send(CRASH);
            }
            /*
            if (i<2) {
            if (observation[j][i]==OBS_A) { a[i]+=1; n[i]+=1; }
            if (observation[j][i]==OBS_B) { b[i]+=1; n[i]+=1; }
            if (observation[j][i]==OBS_H) { h[i]+=1; n[i]+=1; }
            }
            */
            /*
            if (i==0) {
              if (observation[j][i]==OBS_A && observation[j][i+1]==OBS_A) aa++;
              if (observation[j][i]==OBS_A && observation[j][i+1]==OBS_B) ab++;
              if (observation[j][i]==OBS_A && observation[j][i+1]==OBS_H) ah++;
              if (observation[j][i]==OBS_B && observation[j][i+1]==OBS_A) ba++;
              if (observation[j][i]==OBS_B && observation[j][i+1]==OBS_B) bb++;
              if (observation[j][i]==OBS_B && observation[j][i+1]==OBS_H) bh++;
              if (observation[j][i]==OBS_H && observation[j][i+1]==OBS_A) ha++;
              if (observation[j][i]==OBS_H && observation[j][i+1]==OBS_B) hb++;
              if (observation[j][i]==OBS_H && observation[j][i+1]==OBS_H) hh++;
            }
            */
        } /* end loop over indivs (j) */
    } /* end loop over loci (i) */

/*    printf("\n\nL: A=%4d  H=%4d  B=%4d\nR: A=%4d  H=%4d  B=%4d\n\n",
	   a[0],h[0],b[0],a[1],h[1],b[1]); printf(
"AA=%4d  AH=%4d  AB=%4d\nHA=%4d  HH=%4d  HB=%4d\nBA=%4d  BH=%4d  BB=%4d\n\n",
	   aa,ah,ab,ha,hh,hb,ba,bh,bb); */

    /* For F3, the implied recs entries are all constant. Thus, we do
       the same pointer KLUDGE we do for the obs_probs to save space. */

    /* Number of recs in the F2 generation: norecs= 2-recs */
    f2_recs[A][A] = 0.0;
    f2_recs[H][A] = 1.0;
    f2_recs[J][A] = 1.0;
    f2_recs[B][A] = 2.0;
    f2_recs[A][H] = 1.0;
    f2_recs[H][H] = 0.0;
    f2_recs[J][H] = 2.0;
    f2_recs[B][H] = 1.0;
    f2_recs[A][J] = 1.0;
    f2_recs[H][J] = 2.0;
    f2_recs[J][J] = 0.0;
    f2_recs[B][J] = 1.0;
    f2_recs[A][B] = 2.0;
    f2_recs[H][B] = 1.0;
    f2_recs[J][B] = 1.0;
    f2_recs[B][B] = 0.0;

    for (i = 0; i < ctm_state->f3_states; i++) {
        for (j = 0; j < ctm_state->f3_states; j++) {
            ctm_state->implied_recs[0][i][j] = recs = (f2_recs[ctm_state->f2_state[i]][ctm_state->f2_state[j]] +
                                                       (ctm_state->mat_chrom[i] == ctm_state->mat_chrom[j] ? 0.0 : 1.0) +
                                                       (ctm_state->pat_chrom[i] == ctm_state->pat_chrom[j] ? 0.0 : 1.0));
            ctm_state->implied_norecs[0][i][j] = 4.0 - recs;
        }
    }

    /* Lastly, we fill in apriori state probs for leftmost locus only */
    /* In F3 HMM, we enumerate every possible situtation: Thus, 4 F2 genotypes
       each can give rise to 4 F3-Self genotypes. Therefore, the apriori for
       any state is 1/4 * 1/4. These are not normalized. */

    for (j = 0; j < ctm_state->n_indivs; j++) {
        the_obs = ctm_state->observation[j][0];
        if (the_obs == OBS_MISSING) {
            for (k = 0; k < ctm_state->f3_states; k++)
                ctm_state->apriori_prob[j][k] = 0.0625;
        } else { /* not OBS_MISSING */
            for (k = 0; k < ctm_state->f3_states; k++)
                ctm_state->apriori_prob[j][k] = 0.0625 * ctm_state->obs_prob[0][the_obs][k];
        }
    }
    for (k = 0; k < ctm_state->f3_states; k++) ctm_state->apriori_prob[ctm_state->no_data][k] = 0.0625;
}


static int f3_state(ctm_state_t *ctm_state, const char *name /* for debugging */, int f2, int mat, int pat, int obs) {
    int i, j;
    i = ctm_state->f3_states++; /* assign f3 state number, f3_states is a global */
    ctm_state->f2_state[i] = f2;
    ctm_state->mat_chrom[i] = mat;
    ctm_state->pat_chrom[i] = pat;

    /* Assume the KLUDGE allowing us to only have to set obs_prob[0] */
    /* well, actually, let's not and see what happens */

    for (j = 0; j < ctm_state->n_loci; j++) {
        ctm_state->obs_prob[j][OBS_A][i] = (obs == OBS_A ? 1.0 : 0.0);
        ctm_state->obs_prob[j][OBS_B][i] = (obs == OBS_B ? 1.0 : 0.0);
        ctm_state->obs_prob[j][OBS_H][i] = (obs == OBS_H ? 1.0 : 0.0);
        ctm_state->obs_prob[j][OBS_NOT_A][i] = (obs == OBS_H || obs == OBS_B ? 1.0 : 0.0);
        ctm_state->obs_prob[j][OBS_NOT_B][i] = (obs == OBS_H || obs == OBS_H ? 1.0 : 0.0);
        ctm_state->obs_prob[j][OBS_MISSING][i] = 1.0;
    }

    return (i);
}


static void hmm_bc_like_set_probs(ctm_state_t *ctm_state, MAP *map)
/* for now, assume obs_probs are fixed and initialized by setup_hmm() */
{
    int i;
    real rec = 0., norec, theta;
    real **prob;

    /* set the trans_probs (right geno, conditional on left) given theta */
    for (i = 0; i < ctm_state->n_intervals; i++) {
        theta = map->rec_frac[i][0];
        if (ctm_state->ri_type == RI_NOT) rec = theta;
        else if (ctm_state->ri_type == RI_SELF) rec = (2.0 * theta) / (1.0 + 2.0 * theta);
        else if (ctm_state->ri_type == RI_SIB) rec = (4.0 * theta) / (1.0 + 6.0 * theta);
        else
            send(CRASH); /* From Haldane and Wadington (inverse of y(x)) */
        prob = ctm_state->trans_prob[i];
        norec = 1.0 - rec;
        prob[STATE_A][STATE_A] = norec;
        prob[STATE_A][STATE_B] = rec;
        prob[STATE_B][STATE_A] = rec;
        prob[STATE_B][STATE_B] = norec;
    }

    /* For BC/RI we do not need to set any implied_recs/norecs - all entries
       are constant and are initialized by setup_bc1_hmm(). For RI this
       means that we are counting observed genotypic classes (as rec/norec)
       not really the implied_recs. See make_new_map. */
    /* Set obs_probs based on penetrance, etc. - not implemented yet
       exp_observations!=NULL means that obs_probs are not all fixed */
    if (ctm_state->exp_observations != NULL) send(CRASH);
}


static void hmm_f2_set_probs(ctm_state_t *ctm_state, MAP *map) {
/* for now, assume obs_probs are fixed and initialized by setup_hmm() */
    int i;
    real theta, one_minus_theta, no_recs, one_rec, two_recs;
    real **prob, recs;

    /* set the trans_probs (right geno, conditional on left) given theta */
    for (i = 0; i < ctm_state->n_intervals; i++) {
        theta = map->rec_frac[i][0];
        one_minus_theta = 1.0 - theta;
        no_recs = one_minus_theta * one_minus_theta;
        one_rec = theta * one_minus_theta;
        two_recs = theta * theta;
        prob = ctm_state->trans_prob[i];
        prob[STATE_A][STATE_A] = no_recs;
        prob[STATE_A][STATE_B] = two_recs;
        prob[STATE_A][STATE_H] = 2.0 * one_rec;
        prob[STATE_B][STATE_A] = two_recs;
        prob[STATE_B][STATE_B] = no_recs;
        prob[STATE_B][STATE_H] = 2.0 * one_rec;
        prob[STATE_H][STATE_A] = one_rec;
        prob[STATE_H][STATE_B] = one_rec;
        prob[STATE_H][STATE_H] = no_recs + two_recs;

        /* as an optimization we only set implied_recs[H][H] - the
           other entries are constant and are initialized by setup_f2_hmm() */
        recs = (2.0 * two_recs) / (two_recs + no_recs);
        ctm_state->implied_recs[i][STATE_H][STATE_H] = recs;
        ctm_state->implied_norecs[i][STATE_H][STATE_H] = 2.0 - recs;
    }

    /* Set obs_probs based on penetrance, etc. - not implemented yet
       exp_observations!=NULL means that obs_probs are not all fixed */
    if (ctm_state->exp_observations != NULL) send(CRASH);
}


static void hmm_f3_self_set_probs(ctm_state_t *ctm_state, MAP *map) {
/* for now, assume obs_probs are fixed and initialized by setup_hmm() */
    int i, k, j;
    real rec, norec, r0, r1, r2;
    real **prob;
    real f2_recs[4][4];

    /* set the trans_probs (right geno, conditional on left) given theta */
    for (i = 0; i < ctm_state->n_intervals; i++) {
        rec = map->rec_frac[i][0];
        norec = 1.0 - rec;
        r0 = norec * norec;
        r1 = rec * norec;
        r2 = rec * rec;

        /* f2_recs[][] is really the trans_prob for the F2 generation. */
        f2_recs[A][A] = r0;
        f2_recs[H][A] = r1;
        f2_recs[J][A] = r1;
        f2_recs[B][A] = r2;
        f2_recs[A][H] = r1;
        f2_recs[H][H] = r0;
        f2_recs[J][H] = r2;
        f2_recs[B][H] = r1;
        f2_recs[A][J] = r1;
        f2_recs[H][J] = r2;
        f2_recs[J][J] = r0;
        f2_recs[B][J] = r1;
        f2_recs[A][B] = r2;
        f2_recs[H][B] = r1;
        f2_recs[J][B] = r1;
        f2_recs[B][B] = r0;

        prob = ctm_state->trans_prob[i];
        for (j = 0; j < ctm_state->f3_states; j++) {
            for (k = 0; k < ctm_state->f3_states; k++) {
                prob[j][k] = (f2_recs[ctm_state->f2_state[j]][ctm_state->f2_state[k]] *
                              (ctm_state->mat_chrom[j] == ctm_state->mat_chrom[k] ? norec : rec) *
                              (ctm_state->pat_chrom[j] == ctm_state->pat_chrom[k] ? norec : rec));
            }
        }
        /* For F3 we do not need to set any implied_recs/norecs - all entries
           are constant and are initialized by setup_f3_self_hmm() */
    }

    /* Set obs_probs based on penetrance, etc. - not implemented yet
       exp_observations!=NULL means that obs_probs are not all fixed */
    if (ctm_state->exp_observations != NULL) send(CRASH);
}

static void hmm_count_stats(ctm_state_t *ctm_state)
/* for now, assume exp_observations==NULL, and is not used */
{
    int i, k, j;
    real **transitions_exp, **recs_implied = NULL, **norecs_implied = NULL;
    real recs, norecs;

    if (!ctm_state->separate_recs) {
        recs_implied = ctm_state->implied_recs[0];
        norecs_implied = ctm_state->implied_norecs[0];
    }
    for (i = 0; i < ctm_state->n_intervals; i++) {
        if (ctm_state->separate_recs) {
            recs_implied = ctm_state->implied_recs[i];
            norecs_implied = ctm_state->implied_norecs[i];
        }
        transitions_exp = ctm_state->exp_transitions[i];
        recs = norecs = 0.0;
        for (k = 0; k < ctm_state->max_states; k++)
            for (j = 0; j < ctm_state->max_states; j++) {
                recs += transitions_exp[j][k] * recs_implied[j][k];
                norecs += transitions_exp[j][k] * norecs_implied[j][k];
            }
        ctm_state->sufficient_stats[i]->recs = recs;
        ctm_state->sufficient_stats[i]->norecs = norecs;
    }
}


static void hmm_make_new_map(const ctm_state_t *ctm_state, real new_like, MAP *map) {
/* for now, assume no penetrance stuff */
    int i;
    real recs, norecs, sum, theta = 0., ratio;

    for (i = 0; i < ctm_state->n_intervals; i++) {
        if (map->fix_interval[i]) continue;
        sum = recs = ctm_state->sufficient_stats[i]->recs;
        sum += norecs = ctm_state->sufficient_stats[i]->norecs;
        if (sum >= 0.00001) ratio = recs / sum; else ratio = 0.0; /* WHAT TO DO? */
        if (ctm_state->ri_type == RI_NOT) theta = ratio;
        else if (ctm_state->ri_type == RI_SELF) theta = ratio / (2.0 * (1.0 - ratio));
        else if (ctm_state->ri_type == RI_SIB) theta = ratio / (2.0 * (2.0 - 3.0 * ratio));
        else
            send(CRASH); /* From Haldane and Wadington (inverse of y(x)) */

        if (theta > HMM_MAX_THETA) theta = HMM_MAX_THETA;
        else if (theta < HMM_MIN_THETA) theta = HMM_MIN_THETA;
        map->rec_frac[i][0] = theta;
    }
    map->log_like = new_like;
}


static int hmm_e_step(ctm_state_t *ctm_state, real *new_likep) {
/* for now, assume exp_observations==NULL, meaning "don't bother" */
/* return likelihood */
    int indiv, i, the_observation, lstate, rstate, state, leftk, rightk;
    int locus, prev, num_left_states, num_right_states, num_inf;
    int *the_left_state, *the_right_state;
    real sum, total, val, *pheno_given_geno, *the_prob, *prev_prob;
    real indiv_like, total_like;
    real **the_transitions, **my_transitions, *temp, **the_trans_prob;

#ifdef DEBUGGING
    /* debugging */
    int x, y;
    real r, n;
#endif

    /* Sometime unroll array refs in here for a little speed */

    /* For the moment, we don't need all of the indiv_transitions matrix - 
       so we KLUDGE and use entry [0] as BOTH my_transitions, and as temp.
       I will fix this someday  - I THINK I JUST DID IT */

    temp = ctm_state->e_step_temp;

    total_like = 0.0;
    for (i = 0; i < ctm_state->n_intervals; i++)
        for (lstate = 0; lstate < ctm_state->max_states; lstate++) {
            for (rstate = 0; rstate < ctm_state->max_states; rstate++)
                ctm_state->exp_transitions[i][lstate][rstate] = 0.0;
#if 0
            if (ctm_state->error_kludge)
              for (indiv=0; indiv<ctm_state->n_indivs; indiv++)
                ctm_state->exp_genotype[indiv][i][lstate]=0.0;
#endif
        }

    for (indiv = 0; indiv < ctm_state->n_indivs; indiv++) {
        for (locus = 0, num_inf = 0; locus < ctm_state->n_loci; locus++)
            if (ctm_state->observation[indiv][locus] != OBS_MISSING) num_inf++;
        if (num_inf < 2) continue;

        /**** first, compute the left conditioned probs & loglike ****/
        locus = 0; /* leftmost locus */
        total = 0.0;
        /* change this to use n_poss_states */
        for (state = 0; state < ctm_state->max_states; state++)
            total += ctm_state->left_prob[locus][state] = ctm_state->apriori_prob[indiv][state];
        for (state = 0; state < ctm_state->max_states; state++)
            ctm_state->left_prob[locus][state] /= total;
        indiv_like = log10(total);

        /* loop over intervals - prev is marker on left and locus is
           that on right */
        for (prev = 0, locus = 1; locus < ctm_state->n_loci; prev = locus, locus++) {
            num_left_states = ctm_state->n_poss_states[indiv][prev];
            num_right_states = ctm_state->n_poss_states[indiv][locus];
            the_left_state = ctm_state->poss_state[indiv][prev];
            the_right_state = ctm_state->poss_state[indiv][locus];
            the_prob = ctm_state->left_prob[locus];
            prev_prob = ctm_state->left_prob[prev];
            the_trans_prob = ctm_state->trans_prob[prev];
            the_observation = ctm_state->observation[indiv][locus];
            pheno_given_geno = ctm_state->obs_prob[locus][the_observation];
            for (state = 0; state < ctm_state->max_states; state++) the_prob[state] = 0.0;
            total = 0.0;
            for (rightk = 0; rightk < num_right_states; rightk++) {
                rstate = the_right_state[rightk];
                sum = 0.0;
                for (leftk = 0; leftk < num_left_states; leftk++) {
                    lstate = the_left_state[leftk];
                    sum += prev_prob[lstate] * the_trans_prob[lstate][rstate];
                }
                total += temp[rstate] = sum * pheno_given_geno[rstate];
            }
            if (fabs(total) < 0.00001) return 0; /*send(CRASH);*/
            for (rightk = 0; rightk < num_right_states; rightk++) {
                rstate = the_right_state[rightk];
                the_prob[rstate] = temp[rstate] / total;
            }
            indiv_like += log10(total);
        }


        /**** now, the right conditioned probs ****/
        locus = ctm_state->n_loci - 1;
        val = 1.0 / ctm_state->max_states;
        for (state = 0; state < ctm_state->max_states; state++) ctm_state->right_prob[locus][state] = val;

        /* Not the right way: for an F2: 0.25 0.50 0.25, etc
           for (state=0; state<max_states; state++)
           right_prob[locus][state]= apriori_prob[no_data][state]; */

        /* this time, prev is on the right and locus is on the left */
        for (prev = locus, locus = ctm_state->n_loci - 2; locus >= 0; prev = locus, locus--) {
            num_left_states = ctm_state->n_poss_states[indiv][locus];
            num_right_states = ctm_state->n_poss_states[indiv][prev];
            the_left_state = ctm_state->poss_state[indiv][locus];
            the_right_state = ctm_state->poss_state[indiv][prev];
            the_prob = ctm_state->right_prob[locus];
            prev_prob = ctm_state->right_prob[prev];
            the_trans_prob = ctm_state->trans_prob[locus];
            the_observation = ctm_state->observation[indiv][prev]; /* was [locus] */
            pheno_given_geno = ctm_state->obs_prob[prev][the_observation]; /* ditto */
            for (state = 0; state < ctm_state->max_states; state++) the_prob[state] = 0.0;
            total = 0.0;
            for (leftk = 0; leftk < num_left_states; leftk++) {
                lstate = the_left_state[leftk];
                sum = 0.0;
                for (rightk = 0; rightk < num_right_states; rightk++) {
#ifdef NOT_THIS /* CHECK: what is correct PGG below? left or right? */
                    sum+= right_prob[prev][right]
                      * trans_prob[locus][left][right]
                      * pheno_given_geno[left];
#else
                    rstate = the_right_state[rightk];
                    sum += prev_prob[rstate] * the_trans_prob[lstate][rstate]
                           * pheno_given_geno[rstate];
#endif
                }
                total += sum;
                temp[lstate] = sum;
            }
            if (fabs(total) < 0.00001) return 0; /*send(CRASH);*/
            for (leftk = 0; leftk < num_left_states; leftk++) {
                lstate = the_left_state[leftk];
                the_prob[lstate] = temp[lstate] / total;
            }
        }


        /**** and finally, the expected transition matrix for this indiv ****/
        for (i = 0; i < ctm_state->n_intervals; i++) {
            my_transitions = ctm_state->indiv_transitions[i];
            num_left_states = ctm_state->n_poss_states[indiv][i];
            num_right_states = ctm_state->n_poss_states[indiv][i + 1];
            the_left_state = ctm_state->poss_state[indiv][i];
            the_right_state = ctm_state->poss_state[indiv][i + 1];
            prev_prob = ctm_state->left_prob[i];
            the_prob = ctm_state->right_prob[i + 1];
            the_trans_prob = ctm_state->trans_prob[i];
            the_observation = ctm_state->observation[indiv][i + 1];
            pheno_given_geno = ctm_state->obs_prob[i + 1][the_observation];
            total = 0.0;
            for (leftk = 0; leftk < num_left_states; leftk++) {
                lstate = the_left_state[leftk];
                for (rightk = 0; rightk < num_right_states; rightk++) {
                    rstate = the_right_state[rightk];
                    my_transitions[lstate][rstate] =
                            prev_prob[lstate] * the_trans_prob[lstate][rstate]
                            * pheno_given_geno[rstate] * the_prob[rstate];
                    total += my_transitions[lstate][rstate];
                }
            }

            /* normalize and massage these data in with that for all indivs
               also massage into the exp_genotypes global for ERROR KLUDGE
               NB: THE EXP-GENOTYPE FOR THE RIGHTMOST INTERVAL IS NOT SET! */
            the_transitions = ctm_state->exp_transitions[i];
            for (leftk = 0; leftk < num_left_states; leftk++) {
                lstate = the_left_state[leftk];
                for (rightk = 0; rightk < num_right_states; rightk++) {
                    rstate = the_right_state[rightk];
                    my_transitions[lstate][rstate] /= total; /* clean this? */
                    the_transitions[lstate][rstate] +=
                            my_transitions[lstate][rstate];
#if 0
                    if (ctm_state->error_kludge) ctm_state->exp_genotype[indiv][i][lstate]+=
                      my_transitions[lstate][rstate];
#endif
                }
            }
        }

#ifdef DEBUGGING
        test_dump0(left_prob,"left_prob");
        test_dump0(right_prob,"right_prob");
        sprintf(ps,"my_trans-indiv %d ",indiv);
        test_dump(indiv_transitions,max_states,ps);
        r= n= 0.0;
        for (lstate=0; lstate<ctm_state->max_states; lstate++) {
          for (rstate=0; rstate<ctm_state->max_states; rstate++) {
              r+= my_transitions[lstate][rstate]
            * implied_recs[0][lstate][rstate];
              n+= my_transitions[lstate][rstate]
            * implied_norecs[0][lstate][rstate];
          }
        }
        printf("indiv recs=%f  norecs=%f\n",r,n);
#endif

        total_like += indiv_like;
    } /* loop over indivs */

    *new_likep = total_like;
    return 1;
}


int quick_two_pt(int locus0, int locus1, TWO_PT_DATA *two_pt) {
    if (raw.data_type == F2 &&
        (raw.data.f2.cross_type == F2_INTERCROSS ||
         raw.data.f2.cross_type == F2_BACKCROSS)) {
        f2_quick_two_pt(locus0, locus1, two_pt);
    } else { /* bail out - use HMM */
        map2->num_loci = 2;
        map2->locus[0] = locus0;
        map2->locus[1] = locus1;
        map2->rec_frac[0][0] = 0.499;
        map2->fix_interval[0] = FALSE;    /* map2->error_rate=NULL */
        setup_hmm(&ctm_states[0], map2);
        if (!hmm_converge_to_map(&ctm_states[0], map2)) {
            return 0;
        }
        two_pt->theta[NOSEX] = map2->rec_frac[0][0];
        two_pt->lodscore[NOSEX] = map2->log_like - ctm_states[0].null_like;
        if (two_pt->lodscore[NOSEX] < 0.0) two_pt->lodscore[NOSEX] = 0.0;
    }
    return 1;
}
