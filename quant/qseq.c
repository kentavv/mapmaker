/******************************************************************************

  ####    ####   ######   ####            ####
 #    #  #       #       #    #          #    #
 #    #   ####   #####   #    #          #
 #  # #       #  #       #  # #   ###    #
 #   #   #    #  #       #   #    ###    #    #
  ### #   ####   ######   ### #   ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

/********** QSEQ.C - QTL INTERVAL SEQUENCE HANDLERS **********/

//#define INC_LIB
//#define INC_TABLE
//#define INC_SHELL     /* for now need wizard_mode definition, otherwise punt */
//#define INC_CALLQCTM
//#define INC_QTOPLEVEL
//#define INC_QLOWLEVEL
#include "qtl.h"

//void free_qtl_sequence (QTL_SEQUENCE *p);

/***** Global *****/
QTL_SEQUENCE *ints;
char *ints_string;
int *ints_free;
int seq_history_num;
TABLE *named_sequences, *sequence_history;
bool Omore; /* used in for_all_orders iterator macro */

/***** Local *****/
//void seq_init(void);
bool set_qtl_sequence(char *str, char *errmsg, int *errpos);

void free_qtl_sequence(QTL_SEQUENCE *p);

void free_seq_options(QTL_SEQ_OPTION *q);

QTL_SEQUENCE *compile_intervals(char *str);

QTL_SEQUENCE *int_compiler(char **str);

QTL_SEQ_OPTION *mkcontinuous(int trait, real fix_weight);

QTL_SEQ_OPTION *mkinterval(int left, int right);

bool try_right(char **str, QTL_SEQ_OPTION *opt);

bool try_fix_pos(char **str, QTL_SEQ_OPTION *opt);

bool try_range(char **str, QTL_SEQ_OPTION *opt);

void get_genetics_spec(char **str, GENETICS *g);

void enumerate_possibilities(QTL_SEQUENCE *p);

bool reset_state(QTL_SEQUENCE *p, bool wiggle, int *pnum_intervals, int *pcont_vars, int *pnum_orders, int *pwiggle_ints);

//void get_seq_free_genetics(QTL_SEQUENCE *p, bool *free);
bool get_order(QTL_SEQUENCE *p, bool wiggle, QTL_MAP *map);

bool next_order(QTL_SEQUENCE *p, int *perm);

bool next_wiggle(QTL_SEQUENCE *p, int *perm, real cm_step);

bool contig_perm(QTL_SEQUENCE *p, int *perm, bool wiggle);

bool discont_perm(QTL_SEQUENCE *p, int *perm, bool wiggle);

bool wiggle_perm(QTL_SEQUENCE *p, bool *contig, bool *wiggled, real cm_step);

bool test_perm(QTL_SEQUENCE *p, int *perm, bool wiggle, bool perm_rightmost);

bool name_sequence(char *name, char *seq, char **why_not);

bool unname_sequence(char *name, char **why_not);

void add_to_seq_history(char *seq);

char *expand_named_entries(char *str);

void swap_for_dash(char *str);

void unswap_for_dash(char *str);

void name_peaks(WIGGLE_PEAK *peak, char *prefix, bool forget);

//QTL_SEQUENCE *compile_intervals();
//QTL_SEQUENCE *int_compiler();
//QTL_SEQ_OPTION *mkinterval();
//QTL_SEQ_OPTION *mkcontinuous();
//void free_seq_options();
//bool try_fix_pos();
//bool try_range();
//void get_genetics_spec();
//void enumerate_possibilities();
//bool contig_perm();
//bool discont_perm();
//bool wiggle_perm();
//bool test_perm();
//void swap_for_dash();
//void unswap_for_dash();
//
bool try_right(char **str /* may be side-effected */, QTL_SEQ_OPTION *opt /* may be side-effected */);

bool get_order(QTL_SEQUENCE *p, bool wiggle, QTL_MAP *map);

int ichoose(int n, int k); /* Best algorithm for small k */
int add_continuous_var(QTL_MAP *map, int trait, real fix_weight);

char *seqtoken, *seqerr, *seqtoken_ptr;

#define LEFT_ANGLE    "<"
#define RIGHT_ANGLE    ">"
#define LEFT_BRACKET    "["
#define RIGHT_BRACKET    "]"
#define LEFT_SQUIG    "{"
#define RIGHT_SQUIG    "}"
#define PLUS_SIGN       "+"
#define COLON        ":"
#define SURROGATE_DASH    "~"
#define SURROGATE_DASH_char    '~'
#define VERTICAL_BAR    "|"
#define WORD_OF    "of"
#define SEQ_SELF_DELIMITING    "{}[]+:-|~"


void
seq_init(void) {
    ints = NULL;
    array(ints_string, SEQ_LEN + 1, char);
    array(seqtoken, TOKLEN + 1, char);
    seqtoken_ptr = seqtoken;
    array(seqerr, MAXLINE + 1, char);
    sequence_history =
            allocate_table(MAX_HISTORY_SEQS, SEQ_LEN, CANT_EXPAND, INDEX_BY_NUMBER);
    seq_history_num = next_entry_number(sequence_history);
    named_sequences =
            allocate_table(MAX_NAMED_SEQS, SEQ_LEN, EXPANDS_BY(MAX_NAMED_SEQS),
                           INDEX_BY_NAME);
}


bool
set_qtl_sequence(
        char *str, /* str may be uncrunched */
        char *errmsg,
        int *errpos
) {
    QTL_SEQUENCE *new_ints;
    char new_seq[MAXLINE + 1];

    run {
            nstrcpy(new_seq, str, MAXLINE);
            new_ints = compile_intervals(new_seq);
            free_qtl_sequence(ints);
            ints = new_ints;
            add_to_seq_history(str);
            nstrcpy(ints_string, str, SEQ_LEN);
        } except_when(BADSEQ) {
        if (errmsg != NULL) strcpy(errmsg, BADSEQ_errmsg);
        if (errpos != NULL) *errpos = BADSEQ_errpos;
        return (FALSE);
    }
    return (TRUE);
}


void
free_qtl_sequence(QTL_SEQUENCE *p) {
    if (p == NULL || p->dont_free) return;
    free_qtl_sequence(p->next);
    free_seq_options(p->option);
    unarray(p->left, int);
    unarray(p->right, int);
    unarray(p->fix_pos, real);
    unarray(p->count, int);
    unarray(p->contig, bool);
    unsingle(p, QTL_SEQUENCE);
}

void
free_seq_options(QTL_SEQ_OPTION *q) {
    if (q == NULL) return;
    if (q->next != NULL) free_seq_options(q->next);
    unsingle(q, QTL_SEQ_OPTION);
}


/********** COMPILE_INTERVALS **********/

/*** This procedure is used to turn a string representation of the
intervals one is interested in into a heirarchical INTERVALS tree
struct. It returns NULL and prints an error message if
interval_compiler() fails (in this case, nothing is malloced that
isn't freed). Otherwise it mallocs and returns a INTERVALS struct for
the ints string. ***/

#define err_EOL    "Unexpected end of sequence"
#define err_NORIGHTP    "Expected ')'."
#define err_BADPOS    "Expected QTL position (cM or RF) following '+'."
#define err_POSRANGE    "QTL position is greater than interval length."
#define err_BADRIGHT    "Expected a locus name or number following '-' or '|'."
#define err_BADINT    "Zero length or 'backwards' interval."
#define err_BADRANGE    "Zero length or 'backwards' range."
#define err_BADTOKEN    "Expected a locus name or number."
#define err_BADTOKENOR    "Expected a locus name or number, ':', or ']'."
#define err_LEFTANY    "Expected '[', '{', or 'n of' to start."
#define err_LEFTBRACK    "Expected '[' following '%d of'."
#define err_REPEATCONT  "'n of' is not allowed with continuous variables."
#define err_REPEATNUM    "Repeat count must be in range 1...%d (maybe missing '[' or '{')."
#define err_REPEATDEL   "Expected 'of' to follow repeat count (maybe missing '[' or '{')."
#define err_REPEAT    "Repeat count is larger than the number of possible intervals for this QTL."
#define err_INTERX      "Expected 'additive', 'dominant', 'recessive', 'fixed', or 'try'."
#define err_BACKX       "Expected 'free' or 'fixed' (the only options for backcross data)."
#define err_REPEATEST    "Can't have a QTL with both a repeat count and :try."
#define err_FIXFORM    "Expected 'W=<value>' (BC1) or 'A=<value> D=<value>' (F2)."
#define err_NOTLAST "Can't use the last locus in the map as the left locus of an interval."
#define err_NOTRAIT     "Expected trait name(s) and/or number(s) inside {...}"
#define err_TRAITISCV    "Illegal sequence and trait combination.\nThe current trait is used as a continuous variable also."

QTL_SEQUENCE *
compile_intervals(char *str) {
    QTL_SEQUENCE *first, *last, *p;
    char *ptr, *old_self_delim;

    if (nullstr(str)) send(CRASH);
    /* expand_named_entries takes care of any named sequences */
    /* KLUDGE until this code is rewritten */
    ptr = expand_named_entries(str);

    first = last = NULL;
    swap_for_dash(ptr);
    old_self_delim = self_delimiting;
    self_delimiting = ptr_to(SEQ_SELF_DELIMITING);

    run {
            do {
                p = int_compiler(&ptr);
                if (p == NULL) send(CRASH);
                enumerate_possibilities(p);
                if (first == NULL) first = last = p;
                else {
                    last->next = p;
                    last = p;
                }
            } while (!nullstr(ptr));

        } on_exit {
        self_delimiting = old_self_delim;
        BADSEQ_errpos = len(str) - len(ptr); /* makes sense if BADSEQ sent */
        if (BADSEQ_errpos > 0) BADSEQ_errpos--;
        while (white(str[BADSEQ_errpos]) && BADSEQ_errpos > 0) BADSEQ_errpos--;
        unswap_for_dash(str);
        relay_messages;
    }
    return (first);
}


/*** The function int_compiler() really implements compile_intervals(). ***/

/* Note that no ';' should follows calls to FAIL_ or NEXT_TOKEN_ */
#define FAIL_(msg)  { BADSEQ_errmsg=ptr_to(msg); send(BADSEQ); }
#define LOOKAHEAD(chars)  parse_char(str,chars,TRUE,&dummy)
#define NEXT_TOKEN_ \
{ seqtoken=seqtoken_ptr; if (!stoken(str,sREQUIRED,seqtoken)) FAIL_(err_EOL) }

char dummy;

QTL_SEQUENCE *
int_compiler(char **str) {
    int left, right, repeat, trait;
    QTL_SEQUENCE *me;
    QTL_SEQ_OPTION *p, *first, *last;
    GENETICS genetics;
    bool isa_cont_var;

    first = NULL;
    last = NULL;
    me = NULL;
    left = right = NO_LOCUS;
    genetics.backx_weight = DONT_FIX;
    genetics.interx_type = FREE;
    /* 2/21 initialization needed */
    genetics.a = 0.0;
    genetics.b = 0.0;
    genetics.c = 0.0;
    if (nullstr(*str)) return (NULL);

    /* If we have a integer >=1, then it is the number of times to
       repeat this QTL_SEQUENCE interval (we check max_intervals later).
       Otherwise, we expect a left bracket to start this interval off. */

    NEXT_TOKEN_
    if (itoken(&seqtoken, iREQUIRED, &repeat)) {
        if (!irange(&repeat, 1, max_intervals)) {
            sprintf(seqerr, err_REPEATNUM, max_intervals);
            FAIL_(seqerr)
        }
        NEXT_TOKEN_
        if (!streq(lowercase(seqtoken), WORD_OF)) FAIL_(err_REPEATDEL)
        NEXT_TOKEN_
    } else repeat = -1;
    if (streq(seqtoken, LEFT_SQUIG)) isa_cont_var = TRUE;
    else if (streq(seqtoken, LEFT_BRACKET)) isa_cont_var = FALSE;
    else if (repeat > 0) FAIL_(err_LEFTBRACK)
    else FAIL_(err_LEFTANY)

    if (repeat < 0) repeat = 1;
    if (repeat > 1 && isa_cont_var) FAIL_(err_REPEATCONT)

    run {

            /* This is likely to change! */
            if (isa_cont_var)
                do {
                    NEXT_TOKEN_
                    p = NULL; /* So that we can free them when failing */

                    if (wizard_mode && nmatches(seqtoken, "epistasis", 3))
                        p = mkcontinuous(EPISTASIS_TERM, DONT_FIX);
                    else if (valid_trait_str(seqtoken, &trait, seqerr)) {
                        p = mkcontinuous(trait, DONT_FIX); /* No fix-weight for now */

                    } else if (streq(seqtoken, RIGHT_SQUIG)) {
                        if (first == NULL) FAIL_(err_NOTRAIT)
                        else break; /* all done with the do loop */

                    } else FAIL_(seqerr)

                    if (first == NULL) first = last = p;
                    else {
                        last->next = p;
                        last = p;
                    }
                } while (TRUE);

            else
                do { /* !isa_cont_var */
                    NEXT_TOKEN_
                    p = NULL; /* So that we can free them when failing */

                    /* Obsolete code for intervals of the form "<x y>". If we get one
                       it could be followed by fix pos syntax.
                    if (streq(seqtoken,LEFT_ANGLE)) {
                    NEXT_TOKEN_
                    if (streq(seqtoken,"deleted")) break;
                    else if (!valid_locus_str(seqtoken,&left,seqerr)) FAIL_(seqerr)
                    if (left==raw.n_loci-1) FAIL_(err_NOTLAST)
                    right=left+1;
                    NEXT_TOKEN_
                    if (!valid_locus_str(seqtoken,&right,seqerr)) FAIL_(seqerr)
                    if (left>=right) FAIL_(err_BADINT)
                    NEXT_TOKEN_
                    if (!streq(seqtoken,RIGHT_ANGLE)) FAIL_(err_NORIGHTP)
                    p=mkinterval(left,right);
                    try_fix_pos(str,p);
                    } else */

                    /* Otherwise, look for a locus num or name alone... If we get one
                       it could be followed by either '-', '|', or '+'. Both '-'
                       and '+' are not allowed. */
                    if (valid_locus_str(seqtoken, &left, seqerr)) {
                        if (left == raw.n_loci - 1) FAIL_(err_NOTLAST)
                        p = mkinterval(left, left + 1);
                        try_right(str, p);
                        if (!try_fix_pos(str, p)) try_range(str, p);

                        /* Make sure we have gotten at least one interval
                           before we try for a ':' or ']'. */
                    } else if (first == NULL) {
                        FAIL_(err_BADTOKEN)  /* This was commented out. Why? */

                        /* A ']' here then would end the interval.. */
                    } else if (streq(seqtoken, RIGHT_BRACKET)) {
                        break; /* quit the do {} loop */

                        /* The only other possibility is a genetics spec */
                    } else if (streq(seqtoken, COLON)) {
                        get_genetics_spec(str, &genetics);
                        if (raw.data_type == INTERCROSS && repeat > 1 &&
                            genetics.interx_type == TEST_MODELS) FAIL_(err_REPEATEST)
                        /* get_genetics_spec gobbles the ']',
                       so we're done */
                        break;  /* quit the do {} loop */

                        /* Otherwise, it wasn't a ':' or ']', or interval. */
                    } else FAIL_(err_BADTOKENOR)

                    if (first == NULL) first = last = p;
                    else {
                        last->next = p;
                        last = p;
                    }
                } while (TRUE);

            /* If we land here, we're done... */
            single(me, QTL_SEQUENCE);
            me->option = first;
            copy_genetics(&me->genetics, &genetics);
            me->repeat = repeat;
            me->next = NULL;
            me->dont_free = FALSE;
            me->isa_continuous_var = isa_cont_var;
            me->left = NULL;
            me->right = NULL;
            me->fix_pos = NULL;
            me->contig = NULL;
            me->cont_var = NULL;
            me->count = NULL;

        } when_aborting {
        unsingle(me, QTL_SEQUENCE);
        unsingle(p, QTL_SEQ_OPTION);
        while (first != NULL) {
            p = first->next;
            unsingle(first, QTL_SEQ_OPTION);
            first = p;
        }
        relay_messages;
    }
    return (me);
}


QTL_SEQ_OPTION *
mkcontinuous(int trait, real fix_weight) {
    QTL_SEQ_OPTION *p;

    single(p, QTL_SEQ_OPTION);
    p->type = CONT_VAR;
    p->isa.cont_var.num = trait;
    p->isa.cont_var.fix_weight = fix_weight;
    return (p);
}


QTL_SEQ_OPTION *
mkinterval(int left, int right) {
    QTL_SEQ_OPTION *p;

    single(p, QTL_SEQ_OPTION);
    p->type = INTERVAL;
    p->isa.interval.left = left;
    p->isa.interval.right = right;
    p->isa.interval.fix_pos = DONT_FIX;
    return (p);
}


bool try_right(char **str /* may be side-effected */, QTL_SEQ_OPTION *opt /* may be side-effected */) {
    int num;
    char right[TOKLEN + 1];

    if (LOOKAHEAD(VERTICAL_BAR)) {
        if (!stoken(str, sREQUIRED, right)) FAIL_(err_BADRIGHT)
        if (!valid_locus_str(right, &num, seqerr)) FAIL_(seqerr)
        if (num <= opt->isa.interval.left) FAIL_(err_BADINT)
        opt->isa.interval.right = num;
        return (TRUE);
    } else return (FALSE);
}


bool
try_fix_pos(
        char **str, /* may be side-effected */
        QTL_SEQ_OPTION *opt /* may be side-effected */
) {
    real pos, rf;

    if (LOOKAHEAD(PLUS_SIGN)) {
        if (!rtoken(str, rREQUIRED, &pos) || pos < 0.0) FAIL_(err_BADPOS)
        if (pos > 0.5) pos = unhaldane_cm(pos); /* KLUDGE? */
        rf = map_length(opt->isa.interval.left, opt->isa.interval.right);
        if (pos < 0 || pos > (MAX_FRAC_OF_RF * rf)) FAIL_(err_POSRANGE)
        pos = max(pos, MIN_REC_FRAC);
        opt->isa.interval.fix_pos = pos;
        return (TRUE);
    } else return (FALSE);
}


bool
try_range(
        char **str,  /* may be side-effected */
        QTL_SEQ_OPTION *opt /* may be side-effected (turned into a RANGE) */
) {
    int num;
    char right[TOKLEN + 1];

    if (LOOKAHEAD(SURROGATE_DASH)) {
        if (!stoken(str, sREQUIRED, right)) FAIL_(err_BADRIGHT)
        if (!valid_locus_str(right, &num, seqerr)) FAIL_(seqerr)
        if (num <= opt->isa.interval.left) FAIL_(err_BADRANGE)
        opt->isa.range.from = opt->isa.interval.left;
        opt->isa.range.to = num;
        opt->type = RANGE;
        return (TRUE);
    } else return (FALSE);
}


#define INTERX_TYPE_(t) \
if (raw.data_type==INTERCROSS) g->interx_type=t; else FAIL_(err_INTERX)

void
get_genetics_spec(char **str, GENETICS *g) {

    real a, d;
    char c;

    NEXT_TOKEN_
    if (nmatches(seqtoken, "additive", 3)) { INTERX_TYPE_(ADDITIVE)}
    else if (nmatches(seqtoken, "dominant", 3)) { INTERX_TYPE_(DOMINANT)}
    else if (nmatches(seqtoken, "recessive", 3)) { INTERX_TYPE_(RECESSIVE)}
    else if (nmatches(seqtoken, "try", 3)) { INTERX_TYPE_(TEST_MODELS)}
    else if (wizard_mode && nmatches(seqtoken, "f3dominant", 3)) { INTERX_TYPE_(F3DOMINANT)}
    else if (wizard_mode && nmatches(seqtoken, "f3recessive", 3)) { INTERX_TYPE_(F3RECESSIVE)}
    else if (nmatches(seqtoken, "fixed", 3)) {
        if (!parse_char(str, "aAwW", TRUE, &c)) FAIL_(err_FIXFORM)
        if (!parse_char(str, "=", TRUE, &c)) FAIL_(err_FIXFORM)
        if (!rtoken(str, rREQUIRED, &a)) FAIL_(err_FIXFORM)
        if (raw.data_type == INTERCROSS) {
            if (!parse_char(str, "dD", TRUE, &c)) FAIL_(err_FIXFORM)
            if (!parse_char(str, "=", TRUE, &c)) FAIL_(err_FIXFORM)
            if (!rtoken(str, rREQUIRED, &d)) FAIL_(err_FIXFORM)
            g->interx_type = FIXED;
            g->a = a;
            g->b = d;
        } else g->backx_weight = a;
    } else if (raw.data_type == INTERCROSS) FAIL_(err_INTERX)
    else FAIL_(err_BACKX)

    if (!parse_char(str, "]", TRUE, &c)) FAIL_(err_FIXFORM)
}


/********** THINGS FOR BUILDING INTERVALS TREES **********/

void
enumerate_possibilities(QTL_SEQUENCE *p) {
    int i, left, prev_right, contig;
    QTL_SEQ_OPTION *q;
    if (p == NULL) send(CRASH);

    i = 0;
    q = p->option;
    if (q == NULL) send(CRASH);
    do { /* for all options for this interval */
        if (p->isa_continuous_var) {
            if (q->type == CONT_VAR) i++;
            else
                send(CRASH);
        } else { /* !isa_continuous_var */
            if (q->type == INTERVAL) i++;
            else if (q->type == RANGE) i += q->isa.range.to - q->isa.range.from;
            else
                send(CRASH);
        }
    } while ((q = q->next) != NULL);
    p->num_options = i;

    if (p->isa_continuous_var) {
        if (p->repeat != 1) FAIL_(err_REPEATCONT)
        array(p->cont_var, i, int);
        array(p->fix_cont_var_weight, i, real);
        array(p->contig, i, bool);
        array(p->count, p->repeat, int); /* state */

        q = p->option;
        i = 0;
        do { /* once again, for all options for this interval */
            p->cont_var[i] = q->isa.cont_var.num;
            p->fix_cont_var_weight[i] = q->isa.cont_var.fix_weight;
            p->contig[i] = (i == 0 ? TRUE : FALSE);
            i++;
        } while ((q = q->next) != NULL);

    } else { /* !isa_continuous_var */
        if (p->repeat > i) FAIL_(err_REPEAT)
        array(p->left, i, int);
        array(p->right, i, int);
        array(p->fix_pos, i, real);
        array(p->contig, i, bool);
        array(p->count, p->repeat, int); /* state */

        q = p->option;
        i = 0;
        prev_right = -1;
        do { /* once again, for all options for this interval */
            if (q->type == INTERVAL) {
                left = p->left[i] = q->isa.interval.left;
                p->fix_pos[i] = q->isa.interval.fix_pos;
                p->contig[i] = (prev_right < 0 || left == prev_right);
                prev_right = p->right[i] = q->isa.interval.right;
                ++i;
            } else { /* type==RANGE */
                contig = (prev_right < 0 || q->isa.range.from == prev_right);
                for (left = q->isa.range.from; left < q->isa.range.to; left++) {
                    p->left[i] = left;
                    prev_right = p->right[i] = left + 1;
                    p->contig[i] = contig;
                    p->fix_pos[i] = DONT_FIX;
                    ++i;
                    contig = TRUE;
                }
            }
        } while ((q = q->next) != NULL);
    }

    if (i != p->num_options) send(CRASH); /* someone screwed up */
}


/********** THE TREE PERMUTATION BIT **********/


#ifdef TREE
Usage of this stuff: (Typically, commands should use
      qtl_ready() and the for_all_orders() macro, however.)

  reset_state(int_struct,&num_intervals,&n_orders);
  map= alloc_qtl_map(num_intervals);
  do {
      reset_map(map); /* needed because get_order() ADDS ints to map */
      get_order(int_struct,map); /* add ints to map */
      ...
  } while(next_order(map));

get_order(), and reset_map() return bools so that they may be used in
the for_all_orders() macro.
#endif

bool reset_state(p, wiggle, pnum_intervals, pcont_vars, pnum_orders, pwiggle_ints)
        QTL_SEQUENCE *p;
        bool wiggle;
        int *pnum_intervals, *pcont_vars, *pnum_orders, *pwiggle_ints;
{
    int num_orders, num_intervals, num_models, wiggle_ints, i, cont_vars;
    if (p == NULL) send(CRASH);

    num_orders = 1;
    num_intervals = 0;
    wiggle_ints = 0;
    cont_vars = 0;
    do { /* for all interval specs, which may be repeated */
        num_models = 1;
        p->test = 0;
        if (raw.data_type == INTERCROSS && p->genetics.interx_type == TEST_MODELS)
            num_models = 4;
        if (wiggle && p->next == NULL) { /* rightmost interval in a wiggle */
            if (p->repeat > 1 || p->isa_continuous_var) return (FALSE);
            num_orders *= num_models;
            num_intervals += 1;
            wiggle_ints = p->num_options;
            p->count[0] = 0;
            p->interval_cm = haldane_cm(map_length(p->left[0], p->right[0]));
            p->pos_cm = 0.0;
            for (i = 0; i < p->num_options; i++)
                if (p->fix_pos[i] != DONT_FIX) return (FALSE);
        } else {
            num_orders *= (num_models * ichoose(p->num_options, p->repeat));
            if (!p->isa_continuous_var) num_intervals += p->repeat;
            else cont_vars++;  /* FROB - BUG FIX!! */
            for (i = 0; i < p->repeat; i++) p->count[i] = i;
            if (p->isa_continuous_var)
                for (i = 0; i < p->num_options; i++)
                    if (p->cont_var[i] == trait) error(err_TRAITISCV);
            p->pos_cm = p->interval_cm = 0.0;
        }
    } while ((p = p->next) != NULL);

    if (pnum_intervals != NULL) *pnum_intervals = num_intervals;
    if (pcont_vars != NULL) *pcont_vars = cont_vars;
    if (pnum_orders != NULL) *pnum_orders = num_orders;
    if (wiggle && pwiggle_ints != NULL) *pwiggle_ints = wiggle_ints;
    return (TRUE);
}


void
get_seq_free_genetics(
        QTL_SEQUENCE *p,
        bool *free /* [#intervals] should be alloced - elts are side-effected */
) {
    int i, j;

    i = 0;
    do
        for (j = 0; j < p->repeat; j++) {
            if (p->isa_continuous_var) continue; /* don't increment i */
            else if (raw.data_type == INTERCROSS)
                free[i++] = (p->genetics.interx_type == FREE);
            else free[i++] = (p->genetics.backx_weight == DONT_FIX);
        }
    while ((p = p->next) != NULL);
}


bool get_order(QTL_SEQUENCE *p, bool wiggle, QTL_MAP *map) {
    int i, j;
    real fix_pos;
    GENETICS test_genetics;

    if (p == NULL || map == NULL) send(CRASH);

    do
        for (i = 0; i < p->repeat; i++) {
            j = p->count[i];
            if (p->isa_continuous_var) {
                add_continuous_var(map, p->cont_var[j], p->fix_cont_var_weight[j]);
            } else { /* !isa_continuous_var */
                if (wiggle && p->next == NULL) fix_pos = unhaldane_cm(p->pos_cm);
                else fix_pos = p->fix_pos[j];
                if (raw.data_type == INTERCROSS &&
                    p->genetics.interx_type == TEST_MODELS) {
                    test_genetics.interx_type = p->test;
                    test_genetics.a = test_genetics.b = test_genetics.c = 0.0;
                    add_interval(map, p->left[j], p->right[j], fix_pos,
                                 &test_genetics);
                } else /* BACKCROSS or genetics.interx_type!=TEST_MODELS */
                    add_interval(map, p->left[j], p->right[j], fix_pos, &p->genetics);
            }
        }
    while ((p = p->next) != NULL);
    return (TRUE);
}


/* The permutation routines treat cont_var's as a discontiguous permutation,
   simply by the fact that the contig[i] flags are FALSE for i!=0. */

bool
next_order(
        QTL_SEQUENCE *p,
        int *perm /* side-effected to be the interval# permed */
) {
    int permed_int, type;

    permed_int = 0;
    if (contig_perm(p, &permed_int, FALSE)) type = CONTIG;
    else if (test_perm(p, &permed_int, FALSE, TRUE)) type = TEST;
    else if (discont_perm(p, &permed_int, FALSE)) type = SKIP;
    else return (FALSE);

    if (perm != NULL) *perm = type;
    return (TRUE);
}


bool
next_wiggle(
        QTL_SEQUENCE *p,
        int *perm,          /* side-effected to be the interval# permed */
        real cm_step
) {
    int permed_int, last, type, contig, in_wig;
    QTL_SEQUENCE *left, *right;

    if (p->next != NULL) {
        left = p;
        last = 0;
        for (right = p; right->next != NULL; right = right->next) last++;
    } else {
        left = NULL;
        right = p;
        last = 0;
    }

    permed_int = 0;
    if (wiggle_perm(right, &contig, &in_wig, cm_step)) {
        permed_int = last;
        if (in_wig) type = WIGGLE; else type = (contig ? CONTIG : SKIP);
    } else if (test_perm(right, &permed_int, TRUE, TRUE)) type = TEST;
    else if (left == NULL) return (FALSE);
    else if (test_perm(left, &permed_int, TRUE, FALSE)) type = TEST;
    else if (contig_perm(left, &permed_int, TRUE)) type = ORDER;
    else if (discont_perm(left, &permed_int, TRUE)) type = ORDER;
    else return (FALSE);

    if (perm != NULL) *perm = type;
    return (TRUE);
}


bool
contig_perm(
        QTL_SEQUENCE *p,
        int *perm, /* side-effected to be the interval# permed */
        bool wiggle
) {
    int the_perm, i, j, *k, min_k, max_k, contig_end;

    if (wiggle && p->next == NULL) return (FALSE);

    the_perm = *perm + p->repeat;
    if (p->next != NULL && contig_perm(p->next, &the_perm, wiggle)) {
        *perm = the_perm;
        return (TRUE);
    }

    if (p->isa_continuous_var) return (FALSE);

    for (i = p->repeat - 1; i >= 0; i--) { /* repeats right to left */
        max_k = p->num_options - (p->repeat - i);
        contig_end = p->count[i];
        while (contig_end < max_k && p->contig[contig_end + 1]) contig_end++;
        /* contig_end-= i; (p->repeat-i-1); */
        if (p->count[i] < contig_end) { /* can permute, usually... */
            p->count[i] += 1;
            /* except for some wierd fence-post error, which this if handles */
            if (p->repeat > i + 1 && p->count[i] == p->count[i + 1]) continue;
            *perm += i;
            /* contig-reset all intervals to the right of the permuted one */
            do {
                if (i < 0 && wiggle && p->next == NULL) {
                    p->count[0] = 0;
                    p->pos_cm = 0.0;
                    p->interval_cm =
                            haldane_cm(map_length(p->left[0], p->right[0]));
                } else
                    for (j = i + 1; j < p->repeat; j++) {
                        /* back up to the start of the contiguous block - for
                           repeats back up at most to count=count[j-1]. */
                        k = &(p->count[j]);
                        min_k = (j > 0 ? p->count[j - 1] + 1 : 0);
                        while (*k > min_k && p->contig[*k]) --*k;
                    }
                i = -1; /* for looping-all right intervals get all reps reset */
            } while ((p = p->next) != NULL);
            return (TRUE);
        }
    }
    return (FALSE);
}


bool
discont_perm(
        QTL_SEQUENCE *p,
        int *perm, /* side-effected to be the interval# permed */
        bool wiggle
) {
    int the_perm, i, j, *k, min_k, max_k;

    if (p->next == NULL && wiggle) return (FALSE);

    the_perm = *perm + p->repeat;
    if (p->next != NULL && discont_perm(p->next, &the_perm, wiggle)) {
        /* contig-reset to LEFT via recursive returns */
        for (j = 0; j < p->repeat; j++) { /* This works for cont_vars, I think. */
            k = &(p->count[j]);
            min_k = (j > 0 ? p->count[j - 1] + 1 : 0);
            while (*k > min_k && p->contig[*k]) --*k;
            p->test = 0; /* reset genetics permutation state also */
        }
        *perm = the_perm;
        return (TRUE);
    }

    for (i = p->repeat - 1; i >= 0; i--) {  /* intervals right to left */
        /* Permutation is like "while (count[i]<=(#options-(#repeats-i)))" */
        /* This works for cont_vars, I think. */
        max_k = p->num_options - (p->repeat - i);
        while (p->count[i] < max_k && p->contig[p->count[i] + 1]) p->count[i] += 1;
        if (p->count[i] < max_k) {
            p->count[i] += 1;
            *perm += i;
            if (p->contig[p->count[i]]) send(CRASH);
            /* contig-reset all reps to the left of the permuted one */
            for (j = 0; j < i; j++) {
                k = &(p->count[j]);
                min_k = (j > 0 ? p->count[j - 1] + 1 : 0);
                while (*k > j && p->contig[*k]) --*k;
            }
            /* full-reset all intervals to the right of the permuted one */
            do {
                if (i < 0 && wiggle && p->next == NULL) {
                    p->count[0] = 0;
                    p->pos_cm = 0.0;
                    p->interval_cm =
                            haldane_cm(map_length(p->left[0], p->right[0]));
                } else
                    for (j = i + 1; j < p->repeat; j++)
                        p->count[j] = (j > 0 ? p->count[j - 1] + 1 : 0);
                p->test = 0; /* reset genetics permutation state also */
                i = -1; /* for looping-all right intervals get all reps reset */
            } while ((p = p->next) != NULL);
            return (TRUE);
        }
    }
    return (FALSE);
}


bool
wiggle_perm(
        QTL_SEQUENCE *p,
        bool *contig,
        bool *wiggled, /* side-effected */
        real cm_step
) {
    int i;

    if (p->repeat > 1 || p->next != NULL || p->isa_continuous_var) send(CRASH);
    p->pos_cm += cm_step;
    if (p->pos_cm < p->interval_cm) {
        *contig = *wiggled = TRUE;
        return (TRUE);
    }

    else if (p->count[0] < p->num_options - 1) {
        /* don't change this without changing the code in store_wiggle! */
        i = (p->count[0] += 1);
        p->pos_cm = 0.0;
        p->interval_cm = haldane_cm(map_length(p->left[i], p->right[i]));
        *contig = p->contig[i];
        *wiggled = FALSE;
        return (TRUE);
    } else return (FALSE);
}


bool
test_perm(
        QTL_SEQUENCE *p,
        int *perm, /* side-effected to be the interval# permed */
        bool wiggle,
        bool perm_rightmost
) {
    int the_perm, j, *k, min_k;

    if (raw.data_type == BACKCROSS) return (FALSE);
    if (!perm_rightmost && p->next == NULL) return (FALSE);

    the_perm = *perm + p->repeat;
    if (p->next != NULL && test_perm(p->next, &the_perm, wiggle, perm_rightmost)) {
        /* contig-reset to LEFT via recursive returns */
        for (j = 0; j < p->repeat; j++) {
            k = &(p->count[j]);
            min_k = (j > 0 ? p->count[j - 1] + 1 : 0);
            while (*k > min_k && p->contig[*k]) --*k;
        }
        *perm = the_perm;
        return (TRUE);
    }

    if (p->isa_continuous_var) return (FALSE);
    if (p->genetics.interx_type != TEST_MODELS) return (FALSE);
    if (p->repeat > 1) send(CRASH);

    if (p->test < NUM_MODELS - 1) {
        p->test += 1; /* permute model */
        /* contig-reset this interval and those to the RIGHT of it */
        do {
            if (wiggle && p->next == NULL) {
                p->count[0] = 0;
                p->pos_cm = 0.0;
                p->interval_cm = haldane_cm(map_length(p->left[0], p->right[0]));
            } else
                for (j = 0; j < p->repeat; j++) {
                    /* back up to the start of the contiguous block - for
                       repeats back up at most to count=count[j-1]. */
                    k = &(p->count[j]);
                    min_k = (j > 0 ? p->count[j - 1] + 1 : 0);
                    while (*k > min_k && p->contig[*k]) --*k;
                }
        } while ((p = p->next) != NULL);
        return (TRUE);
    }
    return (FALSE);
}


/********** SEQUENCE NAMES **********/


bool
name_sequence(char *name, char *seq, char **why_not) {
    if (name[0] == '*') name++;
    if (!valid_name(name)) {
        *why_not = ptr_to("illegal name");
        return (FALSE);
    }
    put_named_entry(name, seq, context[active_context]->named_sequences);
    return (TRUE);
}

bool
unname_sequence(char *name, char **why_not) {
    int foo, fail;
    char *err = get_temp_string();

    if (!valid_name(name)) {
        *why_not = ptr_to("illegal name");
        return (FALSE);
    }
    if (name[0] == '*') name++;
    if (valid_locus_str(name, &foo, err)) {
        *why_not = ptr_to("name is ambiguous");
        return (FALSE);
    }
    if (!delete_named_entry(name, context[active_context]->named_sequences,
                            &fail)) {
        if (fail == NAME_DOESNT_MATCH)
            *why_not = ptr_to("name is not defined");
        else
            *why_not = ptr_to("name is ambiguous");
        return (FALSE);
    }
    return (TRUE);
}

void
add_to_seq_history(char *seq) {
    int seqnum;

    seqnum = context[active_context]->seq_history_num;
    put_numbered_entry(seq, context[active_context]->sequence_history, &seqnum);
    seqnum++;
    context[active_context]->seq_history_num = seqnum;
}

char *
expand_named_entries(char *str) {
    char *every, *super_exp, *new_str, *tok, *expansion, *dummy;
    char *temp_str;
    int err, num;

    new_str = NULL;
    every = get_temp_string();
    tok = get_temp_string();
    temp_str = get_temp_string();
    run {
            array(new_str, SEQ_LEN + 1, char);
            self_delimiting = ptr_to("[]");
            while (stoken(&str, sREQUIRED, tok)) {
                if (print_mapm_loci) {
                    self_delimiting = ptr_to("[]+:-|");
                    if (itoken(&tok, iREQUIRED, &num)) {
                        sprintf(ps, "%d", num);
                        strcat(new_str, ps);
                        strcat(new_str, " ");
                    } else if (get_named_entry(tok, &dummy, &expansion,
                                               context[active_context]->named_sequences, &err)) {
                        if (streq("all", tok)) {
                            while (stoken(&expansion, sREQUIRED, every)) {
                                if (get_named_entry(every, &dummy, &super_exp,
                                                    context[active_context]->named_sequences,
                                                    &err))
                                    strcat(temp_str, super_exp);
                                while (itoken(&temp_str, iREQUIRED, &num)) {
                                    if (temp_str[0] == '\0') break; /* huh ? */
                                    sprintf(ps, "%d", num);
                                    strcat(new_str, ps);
                                    strcat(new_str, " ");
                                }
                            }
                        } else {
                            while (itoken(&expansion, iREQUIRED, &num)) {
                                if (expansion[0] == '\0') break;
                                sprintf(ps, "%d", num);
                                strcat(new_str, ps);
                                strcat(new_str, " ");
                            }
                        }
                    } else if (tok[0] == '#') {
                        tok = &tok[1];
                        if (itoken(&tok, iREQUIRED, &num)) {
                            --num;
                            if (get_numbered_entry(num, &temp_str,
                                                   context[active_context]->sequence_history)) {
                                /* has obtained sequence reference, now
                                   expand that reference and enter in sequence */
                                temp_str = expand_named_entries(temp_str);
                                strcat(new_str, temp_str);
                                strcat(new_str, " ");
                            }
                        }
                    } else strcat(new_str, tok);
                    strcat(new_str, " ");
                } else {
                    if (get_named_entry(tok, &dummy, &expansion,
                                        context[active_context]->named_sequences, &err))
                        strcat(new_str, expansion);
                    else strcat(new_str, tok);
                    strcat(new_str, " ");
                }
            }
            self_delimiting = ptr_to("");
        } when_aborting {
        unarray(new_str, char);
        relay_messages;
    }
    return (new_str);
}


/* We can't have the '-' symbol be self-delimiting, because then negative
   numbers parse wrong. So, we swap '-' for SURROGATE_DASH which should be
   self-delimiting, and swap back again when compile_sequence() is done. */

void
swap_for_dash( /* internal use only */
        char *str
) {
    int i;
    for (i = 0; str[i] != '\0'; i++) if (str[i] == '-') str[i] = SURROGATE_DASH_char;
}


void
unswap_for_dash( /* internal use only */
        char *str
) {
    int i;
    for (i = 0; str[i] != '\0'; i++) if (str[i] == SURROGATE_DASH_char) str[i] = '-';
}


#ifdef OBSOLETE
int
check_interval (/* make sure all is OK - return T or F */
    int *left,
    int *right,	              /* swap them if ordered wrong */
    real *fix_pos
)
{
    int temp; /* KLUDGED: NEEDS WORK */
    /* check int len? */
    if (!valid_locus_num(left)) return(FALSE);
    if (!valid_locus_num(right)) return(FALSE);
    if (*fix_pos!=DONT_FIX && zero_interval(*left,*right)) return(FALSE);

    if (special_locus(*left)) {
        if (special_locus(*right)) return(FALSE);
        else { temp= *left; *left= *right; *right=temp; } /* swap */
        } else if (zero_interval(*left,*right)) { /* a point */
        if (!valid_interval_num(*left)) return(FALSE);
        ++*right; *fix_pos=0.0;
    } else if (*left>*right && !special_locus(*right)) {
        temp= *left; *left= *right; *right=temp; /* swap */
    }
    if (*fix_pos!=DONT_FIX) {
        if (*fix_pos > map_length(*left,*right))
          *fix_pos=map_length(*left,*right);
    }
    return(TRUE);
}

int
check_range (int *from, int *to, int zero_pos)
{
    int temp;

        if (!valid_locus_num(from) || !valid_locus_num(to) ||
        *from==INF_LOCUS || *from==NO_LOCUS || *from==ANY_LOCUS ||
        *to==INF_LOCUS   || *to==NO_LOCUS   || *to==ANY_LOCUS)
            return(FALSE);
    if (*from>*to) { temp= *from; *from= *to; *to=temp; }
    if (*to==raw.n_loci-1 && !zero_pos) return(FALSE);
       else return(TRUE);
}


void 
sticky_minus (
    char *str /* str must be despaced, and is side-effected */
)
{
    int i, j;
    
    for (i=0, j=0; str[i]!='\0'; i++) 
      if (str[i]=='-' && j>0 && str[j]==' ') j--;
      else if (str[i]==' ' && j>0 && str[j]=='-') continue;
      else str[j++]=str[i];		
    str[j++]='\0'; 
}


#endif


void
name_peaks(WIGGLE_PEAK *peak, char *prefix, bool forget) {
    char *name, *real_name, *seq_str, *seq_to_name, *err;
    int i, fail, count, peaks;

    name = NULL;
    seq_to_name = NULL;
    i = 0;
    count = 0;
    run {
            array(name, NAME_LEN + 1, char);
            array(seq_to_name, SEQ_LEN, char);
            /* forgets old names */
            if (forget) {
                while (TRUE) {
                    sprintf(name, "%s%d", prefix, i + i);
                    if (!delete_named_entry(name,
                                            context[active_context]->named_sequences, &fail))
                        break;
                    else i++;
                }
            } else {  /* else count up to see where to start numbering */
                while (TRUE) {
                    sprintf(name, "%s%d", prefix, count + 1);
                    if (!get_named_entry(name, &real_name, &seq_str,
                                         context[active_context]->named_sequences, &fail))
                        break;
                    else count++;
                }
            }
            /* Now, create the new names */
            for (peaks = 0; peak != NULL; peak = peak->next, peaks++) {
                sprintf(seq_to_name, "%s +%s",
                        interval_str(peak->left, peak->right, FALSE),
                        dist_str(peak->qtl_pos, FALSE));
                if (!streq(genetics_str(&peak->map->constraint
                [peak->map->num_intervals - 1], TRUE), "free")) {
                    sprintf(ps, " :%s", genetics_str(&peak->map->constraint
                    [peak->map->num_intervals - 1], TRUE));
                    strcat(seq_to_name, ps);
                }
                sprintf(name, "%s%d", prefix, count + peaks + 1);
                name_sequence(name, seq_to_name, &err);
            }
        } on_exit {
        unarray(seq_to_name, char);
        unarray(name, char);
    }
}
