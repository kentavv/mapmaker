/******************************************************************************
                                         #####   #####
  ####   #    #     #     ####   #    # #     # #     #           ####
 #    #  #    #     #    #    #  #   #        #       #          #    #
 #    #  #    #     #    #       ####    #####   #####           #
 #  # #  #    #     #    #       #  #   #             #   ###    #
 #   #   #    #     #    #    #  #   #  #       #     #   ###    #    #
  ### #   ####      #     ####   #    # #######  #####    ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#include "mapm.h"

char symbols[] = {
        HYBRID_TYPE_H,
        PARENTAL_TYPE_A,
        PARENTAL_TYPE_B,
        TYPE_NOT_A,
        TYPE_NOT_B
};
#define NUM_SYMS    5

/* internal functions */

real quick_back(LOCUS locus, RECVECTOR **rec_frac, real *conv_like, real *unconv_like);

real quick_f2(LOCUS locus, RECVECTOR **rec_frac, real *conv_like, real *unconv_like);

void probinit(real probdist[6][4]);

void calclike(real prob1[4], real prob2[4], real theta, LOCUS locus, real *likelihood, real *numerator);

int lookup(int c);

int changes(int i, int j);

real f2_prob(real theta, int diffs);

real power(real a, int x);

void f2_quick_two_pt(int loc1, int loc2, TWO_PT_DATA *two_pt, bool sexflag) {
    LOCUS locus;
    real conv_like = 0., unconv_like = 0.;
    RECVECTOR **rec_frac;
    if (raw.data_type != F2) send(CRASH);

    locus.size = 2;
    array(locus.Entry, locus.size, int);
    matrix(rec_frac, locus.size - 1, locus.size, RECVECTOR);
    locus.Entry[0] = loc1;
    locus.Entry[1] = loc2;

    if (raw.data.f2.cross_type == F2_INTERCROSS)
        quick_f2(locus, rec_frac, &conv_like, &unconv_like);
    else if (raw.data.f2.cross_type == F2_BACKCROSS)
        quick_back(locus, rec_frac, &conv_like, &unconv_like);
    else
        send(CRASH);

    two_pt->lodscore[NOSEX] = conv_like - unconv_like;
    two_pt->theta[NOSEX] = rec_frac[0][1][MALE];

    unmatrix(rec_frac, locus.size - 1, RECVECTOR);
    unarray(locus.Entry, int);
}


real quick_back(LOCUS locus, RECVECTOR **rec_frac, real *conv_like, real *unconv_like) {
    int i;
    int indiv;
    int x, y;
    real recs;
    real theta, pr_zero, pr_one;
    int pairs[3];  /* pairs: [0]=A-A, B-B; [1]=A-B; [2]=no data */

    for (i = 0; i < 3; i++)
        pairs[i] = 0;
//    theta=startrecombs;
    for (indiv = 0; indiv < raw.data.f2.num_indivs; indiv++) {
        x = raw.data.f2.allele[locus.Entry[0]][indiv];
        y = raw.data.f2.allele[locus.Entry[1]][indiv];
        if (x == y && x != MISSING_DATA)
            pairs[0]++;
        else if (x != y && x != MISSING_DATA && y != MISSING_DATA)
            pairs[1]++;
        else
            pairs[2]++;
    }
    if (pairs[0] + pairs[1] == 0) {
        rec_frac[0][1][0] = .5;
        return (0.0);
    }


    recs = (real) pairs[1];
    theta = recs / (pairs[0] + pairs[1]);
    theta = (theta > .5) ? .5 : theta;

    pr_zero = (1. - theta);
    pr_one = theta;

    *conv_like = 0.0;
    if (pr_one != 0.0)
        *conv_like = pairs[0] * log10(1. * pr_zero) + pairs[1] * log10(1. * pr_one);
    else
        *conv_like = pairs[0] * log10(1. * pr_zero);

    *unconv_like = 0.0;
    pr_zero = .5;
    pr_one = .5;
    *unconv_like = pairs[0] * log10(1. * pr_zero) + pairs[1] * log10(1. * pr_one);

    rec_frac[0][1][0] = theta;
    return (*conv_like - *unconv_like);
}


real quick_f2(LOCUS locus, RECVECTOR **rec_frac, real *conv_like, real *unconv_like) {
    int i, j;
    long meioses;
    int indiv, histogram[NUM_SYMS + 1][NUM_SYMS + 1];
    char x, y;
    real old_log_like, new_log_like, likelihood, numerator;
    real theta, recs;
    real prob[6][4];

    theta = startrecombs;
    probinit(prob);

    for (i = 0; i < NUM_SYMS + 1; i++) {
        for (j = 0; j < NUM_SYMS + 1; j++)
            histogram[i][j] = 0;
    }

    for (indiv = 0; indiv < raw.data.f2.num_indivs; indiv++) {
        x = raw.data.f2.allele[locus.Entry[0]][indiv];
        y = raw.data.f2.allele[locus.Entry[1]][indiv];
        histogram[lookup(x)][lookup(y)] += 1;
    }                            /* for(indiv) */

    meioses = 0;
    for (i = 0; i < NUM_SYMS; i++) {
        for (j = 0; j < NUM_SYMS; j++) {
            meioses += (long) histogram[i][j];
        }
    }
    meioses *= 2;

    if (meioses == 0) {
        rec_frac[0][1][0] = .5;
        return (0.0);        /* (un)converged likelihood is meaningless */
    }

    new_log_like = 0.0;
    do {
        old_log_like = new_log_like;
        new_log_like = 0.0;
        recs = 0.0;

        for (i = 0; i < NUM_SYMS; i++) {
            for (j = 0; j < NUM_SYMS; j++) {
                if (histogram[i][j] != 0) {
                    calclike(prob[i], prob[j], theta, locus, &likelihood, &numerator);
                    new_log_like = histogram[i][j] * log10(likelihood);
                    recs += (real) histogram[i][j] * numerator / likelihood;
                }
            }
        }
        theta = recs / meioses;
        theta = (theta > .5) ? .5 : theta;
    } while (fabs(new_log_like - old_log_like) > tolerance);

    /* Recalculate likelihood using final theta */

    *conv_like = 0.0;
    for (i = 0; i < NUM_SYMS; i++) {
        for (j = 0; j < NUM_SYMS; j++) {
            if (histogram[i][j] != 0) {
                calclike(prob[i], prob[j], theta, locus, &likelihood, &numerator);
                *conv_like += histogram[i][j] * log10(likelihood);
            }
        }
    }
    *unconv_like = 0.0;
    for (i = 0; i < NUM_SYMS; i++) {
        for (j = 0; j < NUM_SYMS; j++) {
            if (histogram[i][j] != 0) {
                calclike(prob[i], prob[j], 0.5, locus, &likelihood, &numerator);
                *unconv_like += histogram[i][j] * log10(likelihood);
            }
        }
    }
    rec_frac[0][1][0] = theta;
    return (*conv_like - *unconv_like);
}


void probinit(real probdist[6][4]) {
/* A */
    probdist[1][0] = 1.0;
    probdist[1][1] = 0.0;
    probdist[1][2] = 0.0;
    probdist[1][3] = 0.0;
/* B */
    probdist[2][0] = 0.0;
    probdist[2][1] = 0.0;
    probdist[2][2] = 0.0;
    probdist[2][3] = 1.0;
/* C */
    probdist[3][0] = 0.0;
    probdist[3][1] = 1.0 / 3.0;
    probdist[3][2] = 1.0 / 3.0;
    probdist[3][3] = 1.0 / 3.0;
/* D */
    probdist[4][0] = 1.0 / 3.0;
    probdist[4][1] = 1.0 / 3.0;
    probdist[4][2] = 1.0 / 3.0;
    probdist[4][3] = 0.0;
/* H */
    probdist[0][0] = 0.0;
    probdist[0][1] = 1.0 / 2.0;
    probdist[0][2] = 1.0 / 2.0;
    probdist[0][3] = 0.0;
/* Missing Data */
    probdist[5][0] = 1.0 / 4.0;
    probdist[5][1] = 1.0 / 4.0;
    probdist[5][2] = 1.0 / 4.0;
    probdist[5][3] = 1.0 / 4.0;
}


void calclike(real prob1[4], real prob2[4], real theta, LOCUS locus, real *likelihood, real *numerator) {
    int i, j, diffs;
    real d1[4], d2[4], p;
    real f2_prob();

    *likelihood = 0.0;
    *numerator = 0.0;

    for (i = 0; i < 4; i++) {
        if (segregation_distortion) {
            d1[i] = raw.data.f2.allelic_distribution[locus.Entry[0]][i];
            d2[i] = raw.data.f2.allelic_distribution[locus.Entry[1]][i];
        } else {
            d1[i] = 0.25;
            d2[i] = 0.25;
        }
    }

    for (i = 0; i < 4; i++) {
        if (prob1[i] > 0.0) {
            for (j = 0; j < 4; j++) {
                if (prob2[j] > 0.0) {
                    diffs = changes(i, j);
                    p = f2_prob(theta, diffs);
                    *numerator += d1[i] * d2[j] * p * diffs;
                    *likelihood += d1[i] * d2[j] * p;
                }
            }
        }
    }
}


int lookup(int c) {
    /*************************************************************\
    * 		Convert a character to it's index in the      *
    * 	symbol[] array.  If it doesn't appear in symbol[],    *
    * 	return NUM_SYMS.				      *
    \*************************************************************/
    int i;

    for (i = 0; i < NUM_SYMS; i++) {
        if (c == symbols[i]) {
            return (i);
        }
    }
    return (NUM_SYMS);
}


int changes(int i, int j) {
    static int ch[4][4], first;

    if (first == 0) {
        ch[0][0] = 0;
        ch[0][1] = 1;
        ch[0][2] = 1;
        ch[0][3] = 2;
        ch[1][0] = 1;
        ch[1][1] = 0;
        ch[1][2] = 2;
        ch[1][3] = 1;
        ch[2][0] = 1;
        ch[2][1] = 2;
        ch[2][2] = 0;
        ch[2][3] = 1;
        ch[3][0] = 2;
        ch[3][1] = 1;
        ch[3][2] = 1;
        ch[3][3] = 0;
        first = 1;
    }
    return (ch[i][j]);
}


real f2_prob(real theta, int diffs) {
    real answer;
    answer = power(theta, diffs) * power(1.0 - theta, 2 - diffs);
    return (answer);
}


real power(real a, int x) {
    real b, result;
    int x1;

    x1 = x;
    result = 1.0;
    for (b = a; x1 != 0; x1 >>= 1, b = b * b) {
        if (x1 & 1) {
            result *= b;
        }
    }
    return (result);
}


#ifdef HAVE_CEPH

void
ceph_quick_two_pt (int loc1, int loc2, TWO_PT_DATA *two_pt, bool sexflag)
{
    LOCUS locus;
    real conv_like, unconv_like, theta;
    RECVECTOR **rec_frac;

    two_pt_touched=TRUE;
    locus.size=2;
    array(locus.Entry, locus.size, int);
    matrix(rec_frac, locus.size-1, locus.size, RECVECTOR);
    locus.Entry[0]=loc1;
    locus.Entry[1]=loc2;

    if (raw.data_type==F2) {

        two_pt->lodscore[NOSEX]= hmm_quick_two(loc1, loc2, &theta);
    two_pt->theta[NOSEX]= theta;

    } else if (raw.data_type==CEPH) {
        if(using()==USE_PHASE_KNOWN) {
            quick_known(locus, rec_frac, &conv_like, &unconv_like, sexflag);
        } else {
            quick_unknown(locus, rec_frac, &conv_like, &unconv_like, sexflag);
    }

    if(sexflag) {
        two_pt->lodscore[SEXSPEC]=conv_like - unconv_like;
        two_pt->theta[MALE]=rec_frac[0][1][MALE];
        two_pt->theta[FEMALE]=rec_frac[0][1][FEMALE];
    } else {
        two_pt->lodscore[NOSEX]=conv_like - unconv_like;
        two_pt->theta[NOSEX]=rec_frac[0][1][MALE];
    }
    }
    unmatrix(rec_frac, locus.size-1, RECVECTOR);
    unarray(locus.Entry, int);
}


real
quick_known (LOCUS locus, RECVECTOR **rec_frac, real *conv_like, real *unconv_like, bool sexflag)
{
    event_vector skip_flag, x;
    FLDB *block1, *block2;
    int fam, kid, male_recs, male_no_recs, female_recs, female_no_recs,
         recs, no_recs, n;
    real male_p, female_p, p;

    recs=0; male_recs=0; female_recs=0;
    no_recs=0; male_no_recs=0; female_no_recs=0;

    for(fam=0; fam<raw.data.ceph.num_families; fam++) {
        block1=&raw.data.ceph.fldb[fam][locus.Entry[0]];
    block2=&raw.data.ceph.fldb[fam][locus.Entry[1]];
    skip_flag=~(block1->no_DNA | block1->double_het
             | block2->no_DNA | block2->double_het);
    if((block1->pat_phase==KNOWN) && ! block1->pat_homo
       && (block2->pat_phase==KNOWN) && ! block2->pat_homo) {
        for(n=0, x=1; n<block1->nkids; n++, x <<= 1) {
            if(x & skip_flag) {
            if((block1->pat & x)==(block2->pat & x)) {
            male_no_recs++;
            } else {
            male_recs++;
            }
        }		/* if (x & skip_flag) */
        }			/* for (n, x) */
    }			/* if( ... && ... && ... && ... ) */
    if((block1->mat_phase==KNOWN) && ! block1->mat_homo
      && (block2->mat_phase==KNOWN) && ! block2->mat_homo) {
        for(n=0, x=1; n<block1->nkids; n++, x <<= 1) {
        if(x & skip_flag) {
            if((block1->mat & x)==(block2->mat & x)) {
            female_no_recs++;
            } else {
            female_recs++;
            }
        }		/* if (x & skip_flag) */
        }			/* for (n, x) */
    }			/* if( ... && ... && ... && ... ) */
    }				/* for (fam) */
    if(male_recs + male_no_recs==0) {
    male_p=.5;
    } else {
    male_p=(float) male_recs / (float) (male_recs + male_no_recs);
    male_p=max(min(male_p, BIGRECOMBS), SMALLRECOMBS);
    }
    if(female_recs + female_no_recs==0) {
    female_p=.5;
    } else {
    female_p=(float)female_recs / (float)(female_recs + female_no_recs);
    female_p=max(min(female_p, BIGRECOMBS), SMALLRECOMBS);
    }

    recs=male_recs + female_recs;
    no_recs=male_no_recs + female_no_recs;

    if(recs + no_recs==0) {
    p=.5;
    } else {
    p=(float) recs / (float) (recs + no_recs);
    p=max(min(p, BIGRECOMBS), SMALLRECOMBS);
    }

    if(sexflag) {
    rec_frac[0][1][MALE]=male_p;
         rec_frac[0][1][FEMALE]=female_p;
    }
    else {
    rec_frac[0][1][0]=p;
    }
    *conv_like=recs * log10(p) + no_recs * log10(1 - p);
    *unconv_like=(recs + no_recs) * log10(.5);
    return(*conv_like - *unconv_like);
}



real
quick_unknown (LOCUS locus, RECVECTOR **rec_frac, real *conv_like, real *unconv_like, bool sexflag)
{
    int num_parents, n, fam, cross, scene;
    bool pat_inf, mat_inf, p_phase, m_phase;
    real cross0_M, cross0_F, cross1_M, cross1_F, cross1DH, cross0_or_2;
    real m_recs, f_recs, m_no_recs, f_no_recs;
    CROSS *crosses;
    int *n_parents, *n_scenes;
    real alloc0_or_2, alloc1DH_M, alloc1DH_F, thetaHH, denom;
    real like[4], likely, male_prob, female_prob;
    real old_log_like, new_log_like;
    FLDB *block1, *block2;
    event_vector all, x, skip_flag, p_rec, m_rec, double_het;

    array(crosses, raw.data.ceph.num_families, CROSS);
    array(n_parents, raw.data.ceph.num_families, int);
    array(n_scenes, raw.data.ceph.num_families, int);
    for(fam=0; fam<raw.data.ceph.num_families; fam++) {
    for(scene=0; scene<4; scene++) {
        for(cross=0; cross<NUM_CROSS; cross++) {
            crosses[fam][scene][cross]=0;
        }
    }
    block1=&raw.data.ceph.fldb[fam][locus.Entry[0]];
    block2=&raw.data.ceph.fldb[fam][locus.Entry[1]];
    pat_inf=! block1->pat_homo && ! block2->pat_homo;
    mat_inf=! block1->mat_homo && ! block2->mat_homo;
    if(! pat_inf && ! mat_inf) {
        n_parents[fam]=0;
        n_scenes[fam]=0;
    } else if(pat_inf && !mat_inf) {
        skip_flag=~(block1->no_DNA | block1->double_het
                | block2->no_DNA | block2->double_het);
        n_parents[fam]=1;
        for(n=0, x=1; n<block1->nkids; n++, x <<= 1) {
        if(x & skip_flag) {
            if((block1->pat & x)==(block2->pat & x)) {
            crosses[fam][0][CROSS0_M]++;
            } else {
            crosses[fam][0][CROSS1_M]++;
            }
        }		/* if(x & skip_flag) */
        }			/* for(n, x) */
        if((block1->pat_phase==KNOWN) && (block2->pat_phase==KNOWN)) {
        n_scenes[fam]=1;
        } else {
        n_scenes[fam]=2;
        crosses[fam][1][CROSS0_M]=crosses[fam][0][CROSS1_M];
        crosses[fam][1][CROSS1_M]=crosses[fam][0][CROSS0_M];
        }
    } else if(! pat_inf && mat_inf) {
        skip_flag=~(block1->no_DNA | block1->double_het
                | block2->no_DNA | block2->double_het);
        n_parents[fam]=1;
        for(n=0, x=1; n<block1->nkids; n++, x <<= 1) {
        if(x & skip_flag) {
            if((block1->mat & x)==(block2->mat & x)) {
            crosses[fam][0][CROSS0_F]++;
            } else {
            crosses[fam][0][CROSS1_F]++;
            }
        }		/* if(x & skip_flag) */
        }			/* for(n, x) */
        if((block1->mat_phase==KNOWN) && (block2->mat_phase==KNOWN)) {
        n_scenes[fam]=1;
        } else {
        n_scenes[fam]=2;
        crosses[fam][1][CROSS0_F]=crosses[fam][0][CROSS1_F];
        crosses[fam][1][CROSS1_F]=crosses[fam][0][CROSS0_F];
        }
    } else {		/* pat_inf && mat_inf */
        n_parents[fam]=2;
        skip_flag=~(block1->no_DNA | block2->no_DNA);
        p_rec=block1->pat ^ block2->pat;
        m_rec=block1->mat ^ block2->mat;
        p_phase=(block1->pat_phase==KNOWN)
               && (block2->pat_phase==KNOWN);
        m_phase=(block1->mat_phase==KNOWN)
               && (block2->mat_phase==KNOWN);
        double_het=block1->double_het | block2->double_het;
        all=0;
        for(n=0, x=1; n<block1->nkids; n++, x <<= 1) {
        all=all | x;
        }
        if(p_phase && m_phase) {
        n_scenes[fam]=1;
        do_scene(p_rec, m_rec, double_het, crosses[fam][0],
             block1->nkids, skip_flag);
        } else if(p_phase) {
        n_scenes[fam]=2;
        do_scene(p_rec, m_rec, double_het, crosses[fam][0],
             block1->nkids, skip_flag);
        do_scene(p_rec, m_rec ^ all, double_het, crosses[fam][1],
             block1->nkids, skip_flag);
        } else if(m_phase) {
        n_scenes[fam]=2;
        do_scene(p_rec, m_rec, double_het, crosses[fam][0],
             block1->nkids, skip_flag);
        do_scene(p_rec ^ all, m_rec, double_het, crosses[fam][1],
             block1->nkids, skip_flag);
        } else {
        n_scenes[fam]=4;
        do_scene(p_rec, m_rec, double_het, crosses[fam][0],
             block1->nkids, skip_flag);
        do_scene(p_rec, m_rec ^ all, double_het, crosses[fam][1],
             block1->nkids, skip_flag);
        do_scene(p_rec ^ all, m_rec, double_het, crosses[fam][2],
             block1->nkids, skip_flag);
        do_scene(p_rec^all, m_rec^all, double_het, crosses[fam][3],
             block1->nkids, skip_flag);
        }
    }			/* if */
    }				/* for (fam) */

    male_prob=startrecombs;
    female_prob=startrecombs;
    new_log_like=0.0;
    do {
        old_log_like=new_log_like;
    new_log_like=0.0;
    m_recs=0;
    f_recs=0;
    m_no_recs=0;
    f_no_recs=0;
    for(fam=0; fam<raw.data.ceph.num_families; fam++){
        num_parents=n_parents[fam];
        if(num_parents==0) {
        continue;
        }
        likely=0.0;
        cross0_M=0;
        cross0_F=0;
        cross1_M=0;
        cross1_F=0;
        cross1DH=0;
        cross0_or_2=0;
        for (scene=0; scene<n_scenes[fam]; scene++) {
            like[scene] =
            power((1.0 - male_prob), crosses[fam][scene][CROSS0_M])
            * power((1.0 - female_prob),
                    crosses[fam][scene][CROSS0_F])
            * power(male_prob, crosses[fam][scene][CROSS1_M])
            * power(female_prob, crosses[fam][scene][CROSS1_F])
            * power((male_prob * (1.0 - female_prob)
                  + (female_prob) * (1.0 - male_prob)) / 2.0,
                crosses[fam][scene][CROSS1DH])
            * power(((male_prob * female_prob
                  + (1.0 - male_prob) * (1.0 - male_prob)) / 2.0),
                crosses[fam][scene][CROSS0_OR_2]);
            likely += like[scene];
        }			/* for (scene) */
        new_log_like += log10(likely);
        for(scene=0; scene<n_scenes[fam]; scene++) {
        cross0_M += (like[scene] / likely)
            * crosses[fam][scene][CROSS0_M];
        cross0_F += (like[scene] / likely)
            * crosses[fam][scene][CROSS0_F];
        cross1_M += (like[scene] / likely)
            * crosses[fam][scene][CROSS1_M];
        cross1_F += (like[scene] / likely)
            * crosses[fam][scene][CROSS1_F];
        cross1DH += (like[scene] / likely)
            * crosses[fam][scene][CROSS1DH];
        cross0_or_2 += (like[scene] / likely)
                 * crosses[fam][scene][CROSS0_OR_2];
        }
        likely /= n_scenes[fam];
            thetaHH=((male_prob * female_prob)
                + (1. - male_prob)*(1. - female_prob));
            alloc0_or_2=(male_prob * female_prob) / thetaHH;
        denom=male_prob * (1.0 - female_prob)
              + female_prob * (1.0 - male_prob);
        alloc1DH_M=male_prob * (1.0 - female_prob) / denom;
        alloc1DH_F=female_prob * (1.0 - male_prob) / denom;
            m_recs += cross1_M + alloc1DH_M * cross1DH
            + alloc0_or_2 * cross0_or_2;
            f_recs += cross1_F + alloc1DH_F * cross1DH
            + alloc0_or_2 * cross0_or_2;
        m_no_recs += cross0_M +(1.0 - alloc1DH_M) * cross1DH
            + (1.0 - alloc0_or_2) * cross0_or_2;
        f_no_recs += cross0_F + (1.0 - alloc1DH_F) * cross1DH
            + (1.0 - alloc0_or_2) * cross0_or_2;
    }			/* for (fam) */
    if(m_recs + f_recs + m_no_recs + f_no_recs==0) {
        male_prob=startrecombs;
        female_prob=startrecombs;
    } else {
        if(sexflag) {
        if(m_recs + m_no_recs==0) {
            male_prob=startrecombs;
        } else {
            male_prob=rmaxf(rminf(m_recs / (m_recs + m_no_recs),
                BIGRECOMBS), SMALLRECOMBS);
        }
        if(f_recs + f_no_recs==0) {
            female_prob=startrecombs;
        } else {
            female_prob=rmaxf(rminf(f_recs / (f_recs + f_no_recs),
                  BIGRECOMBS), SMALLRECOMBS);
        }
        } else {
        male_prob=rmaxf(rmin((m_recs+f_recs)/
                      (m_recs+f_recs+m_no_recs+f_no_recs),
                      BIGRECOMBS),SMALLRECOMBS);
        female_prob=male_prob;
        }
    }
    new_log_like=calc_log_like(crosses, n_parents, n_scenes,
             like, male_prob, female_prob);
    } while(fabs(new_log_like - old_log_like) > tolerance);

    *conv_like=new_log_like;
    *unconv_like=calc_log_like(crosses, n_parents, n_scenes, like, .5, .5);
    rec_frac[0][1][MALE]=male_prob;
    rec_frac[0][1][FEMALE]=female_prob;

    unarray(crosses, CROSS);
    unarray(n_parents, int);
    unarray(n_scenes, int);

    return(*conv_like - *unconv_like);
}



void do_scene(p_rec, m_rec, double_het, cross, nkids, skip_flag)
event_vector p_rec, m_rec, double_het;
long cross[];
int nkids;
event_vector skip_flag;

{
    int n;
    event_vector x;

    for(n=0, x=1; n<nkids; n++, x <<= 1) {
    if(x & skip_flag) {
        if(x & double_het) {
        if(p_rec & m_rec & x) {
            cross[CROSS0_OR_2]++;
        } else if(~p_rec & ~m_rec & x) {
            cross[CROSS0_OR_2]++;
        } else {
            cross[CROSS1DH]++;
        }
        } else {
        if(p_rec & m_rec & x) {
            cross[CROSS1_M]++;
            cross[CROSS1_F]++;
        } else if(p_rec & x) {
            cross[CROSS1_M]++;
            cross[CROSS0_F]++;
        } else if(m_rec & x) {
            cross[CROSS1_F]++;
            cross[CROSS0_M]++;
        } else {
            cross[CROSS0_M]++;
            cross[CROSS0_F]++;
        }
        }
    }			/* if(x & skip_flag) */
    }				/* for(n, x) */
}



real calc_log_like(crosses, n_parents, n_scenes, like, m_p, f_p)
CROSS *crosses;
int *n_parents, *n_scenes;
real *like, m_p, f_p;

{
    int num_parents, fam, scene;
    real likely, new_log_like;

    new_log_like=0.0;
    for (fam=0; fam<raw.data.ceph.num_families; fam++) {
    num_parents=n_parents[fam];
    if (num_parents==0) {
        continue;
    }
    likely=0.0;
    for (scene=0; scene<n_scenes[fam]; scene++) {
        like[scene] =
        power((1.0 - m_p), crosses[fam][scene][CROSS0_M])
        * power((1.0 - f_p), crosses[fam][scene][CROSS0_F])
        * power(m_p, crosses[fam][scene][CROSS1_M])
        * power(f_p, crosses[fam][scene][CROSS1_F])
        * power((m_p * (1.0 - f_p) + (f_p) * (1.0 - m_p)) / 2.0,
            crosses[fam][scene][CROSS1DH])
        * power(((m_p * f_p + (1.0 - m_p) * (1.0 - m_p)) / 2.0),
            crosses[fam][scene][CROSS0_OR_2]);
        likely += like[scene];
    }			/* for (scene) */
    likely /= n_scenes[fam];
    new_log_like += log10 (likely);
    }				/* for (fam) */
    return (new_log_like);
}

#endif


