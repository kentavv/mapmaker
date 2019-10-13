/******************************************************************************

  ####    #####  #####   #          #    #####            ####
 #          #    #    #  #          #    #    #          #    #
  ####      #    #    #  #          #    #####           #
      #     #    #####   #          #    #    #   ###    #
 #    #     #    #   #   #          #    #    #   ###    #    #
  ####      #    #    #  ######     #    #####    ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

/******************** STRLIB.C - STRING & PARSING STUFF ********************/

//#define INC_LIB
//#define INC_HELP_DEFS
#include "system.h"

//int  get_token();

char Cw, Ct, *self_delimiting, *null_string;

void nstrcpy(char *to, char *from, int n) {
    int i;

    for (i = 0; from[i] != '\0' && i < n; i++) to[i] = from[i];
    to[i] = '\0';
}


void
maxstrcat(char *to, char *from, int n) {
    int i, j;

    if ((i = len(to)) >= n) {
        to[n] = '\0';
        return;
    }
    for (j = 0; from[j] != '\0' && i < n; i++, j++) to[i] = from[j];
    to[i] = '\0';
}


void
nstrcat(char *to, char *from, int n) {
    int i, j;

    i = len(to);
    for (j = 0; from[j] != '\0' && j < n; i++, j++) to[i] = from[j];
    to[i] = '\0';
}


void
strins(char *to, char *from) {
    int i, current_length, num_to_ins;

    current_length = len(to);
    num_to_ins = len(from);
    if (num_to_ins == 0) return;
    for (i = current_length; i >= 0; i--) to[i + num_to_ins] = to[i];
    for (i = 0; i < num_to_ins; i++) to[i] = from[i];
}


void
nstrins(char *to, char *from, int max_to_ins) {
    int i, current_length, num_to_ins;

    current_length = len(to);
    num_to_ins = iminf(len(from), max_to_ins);
    if (num_to_ins == 0) return;
    for (i = current_length; i >= 0; i--) to[i + num_to_ins] = to[i];
    for (i = 0; i < num_to_ins; i++) to[i] = from[i];
}


void
maxstrins(char *to, char *from, int max_total_length) {
    int i, current_length, num_to_ins;

    current_length = len(to);
    if (current_length >= max_total_length) return;
    num_to_ins = iminf(len(from), max_total_length - current_length);
    if (num_to_ins == 0) return;
    for (i = current_length; i >= 0; i--) to[i + num_to_ins] = to[i];
    for (i = 0; i < num_to_ins; i++) to[i] = from[i];
}


void
strdel(char *str, int num_chars) {
    int i;
    char *rest;

    if (len(str) <= num_chars) {
        str[0] = '\0';
        return;
    }
    rest = str + num_chars;
    i = 0;
    do str[i] = rest[i]; while (str[i++] != '\0');
}


char *
mkstrcpy(char *str) {
    char *foo;

    array(foo, istrlen(str) + 1, char);
    strcpy(foo, str);
    return (foo);
}


char *
ptr_to(char *str) { return (str); } /* A Kludge to be sure, but the only way to do this in C */


int
nullstr(char *str) {
    int i;

    if (str == NULL) return (TRUE);
    for (i = 0; str[i] != '\0'; i++)
        if (!white(str[i]) && !trash(str[i])) return (FALSE);
    return (TRUE);
}


/* 2/21 - the name strfind has been changed to strfinder due to
   a conflict with something of the same name in libc.a on
   some UNIX machines. */

int
strfinder(    /* index of chr, or -1 if can't find chr */
        char *str,
        int chr
) {
    int i;

    if (chr == '\0') return (strlen(str));
    else for (i = 0; str[i] != '\0'; i++) if (str[i] == chr) return (i);
    return (-1);
}


bool
nmatches(char *s1, char *s2, int n)
/* match at least n chars in the first token in s1 to s2 "template" 
   for example "pr" matches "print" (n<=2). 
   if s1=="" or s2=="" it returns FALSE, if n==0 it returns TRUE 
   if s1 is longer than s2 it must return FALSE
   it assumes both s1 and s2 are despace()ed and filter()ed, but not 
   necessarily lowercase()ed */
{
    int i;

    if (s1[0] == '\0' || s2[0] == '\0') return (FALSE);
    if (n == 0) return (TRUE);
    for (i = 0; s1[i] != '\0'; i++) {
        if (s2[i] == '\0') return (FALSE);
        if (tolower(s1[i]) != tolower(s2[i])) return (FALSE);
        if (white(s1[i])) break;
    }
    return (i >= n ? TRUE : FALSE);
}


int
xstreq(    /* !!! currently broken !!! */
        char *s1,
        char *s2
) {
    int i, j, preceeding_space;

    i = 0;
    j = 0;
    while (white(s1[i]) || trash(s1[j])) i++;
    while (white(s2[j]) || trash(s2[j])) j++;

    for (preceeding_space = FALSE; s1[i] != '\0'; i++)
        if (trash(s1[i])) continue;
        else if (white(s1[i])) {
            preceeding_space = TRUE;
            continue;
        }
        else { /* char must be matched, checking for space first */
            while (trash(s2[j])) j++;
            if (preceeding_space) {
                if (!white(s2[j])) return (FALSE); /* incl \0 */
                do j++; while (white(s2[j]) || trash(s2[j]));
                preceeding_space = FALSE;
            }
            if (tolower(s1[i]) != tolower(s2[j])) /* incl s2[j]==\0 */
                return (FALSE);
            do j++; while (trash(s2[j]));
            continue;
        }

    /* here, we do not care about preceeding_space */
    if (nullstr(&s2[j])) return (TRUE); else return (FALSE);
}


char *
truncstr(char *str, int length) {
    if (str == NULL) send(CRASH);
    if (len(str) > length) str[length] = '\0';
    return (str);
}


char *
pad_to_len(char *str, int length) {
    int start, i;

    start = len(str);
    if (start > length) return (str);
    for (i = start; i < length; i++) str[i] = ' ';
    str[length] = '\0';
    return (str);
}


char *
append_spaces(char *str, int num) {
    int i, j;

    for (i = len(str), j = 0; j < num; i++, j++) str[i] = ' ';
    str[i] = '\0';
    return (str);
}


char *
despace(char *str)
/* All whitespace is changed to single spaces. Leading and trailing whitespace
   is done away with entirely. */
{
    int i, j, space;

    if (str == NULL) send(CRASH);
    for (i = 0, j = 0, space = TRUE; str[i] != '\0'; i++)
        if (white(str[i])) {
            if (!space) {
                space = TRUE;
                if (j != 0) str[j++] = ' ';
            } /* else space==TRUE and we ignore this char */
        } else {  /* !white */
            space = FALSE;
            str[j++] = str[i];
        }
    if (j > 0 && str[j - 1] == ' ') j--; /* Step back over a trailing space */
    str[j] = '\0';
    return (str);
}


char *
lowercase(char *str) {
    int i;
    char c;

    if (str == NULL) send(CRASH);
    for (i = 0; str[i] != '\0'; i++)
        if ((c = str[i]) >= 'A' && c <= 'Z')
            str[i] = c - 'A' + 'a';
    return (str);
}


char *
_filter(char *str) {
    int i, j;

    for (i = 0, j = 0; str[i] != '\0'; i++)
        if (!trash(str[i])) str[j++] = str[i];
    str[j] = '\0';
    return (str);
}


char *
filter_nonspaces(char *str) {
    int i, j;

    for (i = 0, j = 0; str[i] != '\0'; i++)
        if (white(str[i])) str[j++] = ' ';
        else if (!trash(str[i])) str[j++] = str[i];
    str[j] = '\0';
    return (str);
}


char *
crunch(char *str) {
    despace(str);
    _filter(str);
    lowercase(str);
    return (str);
}


char *
uppercase(char *str) {
    int i;
    char c;

    if (str == NULL) send(CRASH);
    for (i = 0; str[i] != '\0'; i++)
        if ((c = str[i]) >= 'a' && c <= 'z')
            str[i] = c - 'a' + 'A';
    return (str);
}


/* get_token() is an internal routine for token parsers below...
   If it fails, then p_rest points to the '\0' at the end of the string and 
   tok=="". The global self_delimiting indicates what chars always split into 
   tokens, whether surrounded by whitespace or not. 
   For example, parsing "(1+ 12)" with self_delimiting=="()+" 
   will give tokens "(","1","+","12", and ")". Length indicates the max # of
   characters that will be parsed, if this is exceeded the rest are thrown 
   out and *truncated=TRUE. */

#define self(c) strin(self_delimiting,c)

bool
get_token(
        char *str,
        char *tok,
        char **p_rest,
        int length,
        bool *truncated /* side-effected */
) {
    int i, j;

    *truncated = FALSE;
    i = 0;
    while (white(str[i])) i++;
    if (str[i] == '\0') {
        tok[0] = '\0';
        *p_rest = str + i;
        return (FALSE);
    }
    if (self(str[i])) { /* must be !white() here */
        tok[0] = str[i];
        tok[1] = '\0';
        *p_rest = str + i + 1;
        return (TRUE);
    }
    for (j = 0; str[i] != '\0' && !white(str[i]) && !self(str[i]); i++)
        if (j >= length) {
            /* gobble to the next delimiter */
            while (str[i] != '\0' && !self(str[i]) && !white(str[i])) i++;
            tok[j] = '\0';
            *p_rest = str + i;
            *truncated = TRUE;
            return (TRUE);
        } else {
            tok[j++] = str[i];
            continue;
        }

    tok[j] = '\0';
    while (str[i] != '\0' && white(str[i])) i++;
    *p_rest = str + i;
    return (TRUE);
}


bool
split_string(
        char *str,
        char **rest,
        int divider /* character that the string should be split on */
) {
    int i, j;

    *rest = NULL;
    i = 0;
    for (i = 0; str[i] != '\0'; i++)
        if (str[i] == divider) {
            j = i + 1;
            while (white(str[j])) j++;
            *rest = str + j;
            i--;
            while (i > 0 && white(str[i])) i--;
            str[i + 1] = '\0';
            return (TRUE);
        }
    return (FALSE);
}


/**** see strlib.h for info on itoken(), ltoken(), rtoken() and stoken() ****/

char *tok, *x; /* local */

int
itoken(char **p_str, int default_val, int *p_val) {
    char *rest;
    long foo;
    bool toolong;

    if (!get_token(*p_str, tok, &rest, TOKLEN, &toolong)) {
        *p_val = default_val;
        *p_str = rest; /* will be the end of str */
        if (default_val == iREQUIRED) return (FALSE); else return (TRUE);
    } else if (toolong || sscanf(tok, "%ld%s", &foo, x) != 1 ||
               foo > MAXINT || foo < -MAXINT) {
        /* a token exists, but it's bad, so we leave *p_str alone */
        *p_val = default_val;
        return (FALSE);
    } else { /* all is OK: set *p_str and val */
        *p_val = (int) foo;
        *p_str = rest;
        return (TRUE);
    }
}


int
ltoken(char **p_str, long default_val, long *p_val) {
    char *rest;
    long foo;
    bool toolong;

    if (!get_token(*p_str, tok, &rest, TOKLEN, &toolong)) {
        *p_val = default_val;
        *p_str = rest; /* will be the end of str */
        if (default_val == lREQUIRED) return (FALSE); else return (TRUE);
    } else if (toolong || sscanf(tok, "%ld%s", &foo, x) != 1) {
        /* a token exists, but it's bad, so we leave *p_str alone */
        *p_val = default_val;
        return (FALSE);
    } else { /* all is OK: set *p_str and val */
        *p_val = foo;
        *p_str = rest;
        return (TRUE);
    }
}


int
rtoken(char **p_str, real default_val, real *p_val) {
    char *rest;
    double foo;
    bool toolong;

    if (!get_token(*p_str, tok, &rest, TOKLEN, &toolong)) {
        *p_val = default_val;
        *p_str = rest;
        if (default_val == rREQUIRED) return (FALSE); else return (TRUE);
    } else if (toolong || sscanf(tok, "%lf%s", &foo, x) != 1 ||
               foo > MAXREAL || foo < -MAXREAL ||
               (foo > ((double) 0) && foo < MINREAL) ||
               (foo < ((double) 0) && foo > MINREAL)) { /* leave *p_str */
        *p_val = default_val;
        return (FALSE);
    } else { /* token is OK */
        *p_val = (real) foo;
        *p_str = rest;
        return (TRUE);
    }
}


/* only to be driven by macros in strlib.c */
int
stok(char **p_str, char *default_val, char *val, int num_chars, int allow_truncation, char *ok_chars) {
    char *rest;
    bool toolong;
    int i;

    if (!get_token(*p_str, val, &rest, num_chars, &toolong)) {
        *p_str = rest; /* will be the '\0' at the end of str */
        if (default_val == sREQUIRED) {
            val[0] = '\0';
            return (FALSE);
        }
        else {
            nstrcpy(val, default_val, TOKLEN);
            return (TRUE);
        }
    } else {
        if (!allow_truncation && toolong) {
            if (default_val == NULL)
                val = NULL;
            else
                nstrcpy(val, default_val, num_chars);
            return (FALSE);
        }
        if (ok_chars != NULL)
            for (i = 0; i <= len(*p_str); i++)
                if (!strin(ok_chars, (*p_str)[i])) {
                    nstrcpy(val, default_val, num_chars);
                    return (FALSE);
                }
    }
    *p_str = rest; /* The next token - val is set */
    return (TRUE);
}

int
stoken(char **p, char *def, char *val) { return (stok(p, def, val, TOKLEN, TRUE, NULL)); }


int
parse_char( /* see strlib.h for info */
        char **p_str,
        char *ok_list,
        int skip_whitespace,
        char *p_ch
) {
    int i;

    i = 0;
    while ((*p_ch = (*p_str)[i++]) != '\0') {
        if (*p_ch == ' ' || *p_ch == '\n' || *p_ch == '\t') {
            if (skip_whitespace) continue;
            else {
                *p_ch = '\0';
                return (FALSE);
            } /* leave p_str */
        } else if (nullstr(ok_list) || strchr(ok_list, *p_ch)) {
            *p_str += i;
            return (TRUE);
        } else return (FALSE); /* leave p_str; p_ch=the bad char */
    }
    *p_str += i - 1;  /* hit end of the string, set p_str to point to \0 */
    return (FALSE); /* p_ch=\0 */
}


void
parse_whitespace(char **p_str) { while (**p_str != '\0' && white(**p_str)) ++*p_str; }


int
count_tokens(char *str) {
    int i, n;

    if (str == NULL) send(CRASH);
    i = n = 0;
    do {
        while (white(str[i]) || strin(self_delimiting, str[i]))
            if (str[i] == '\0') return (n); else i++;
        n++;
        while (!white(str[i]) && !strin(self_delimiting, str[i]))
            if (str[i] == '\0') return (n); else i++;
    } while (TRUE);
    return (-1); /* never reached */
}


bool
is_a_token(char *str) {
    int length, i;

    if ((length = len(str)) == 0) return (FALSE);
    if (length == 1 && strin(self_delimiting, str[0])) return (TRUE);
    for (i = 0; i < length; i++) if (white(str[i])) return (FALSE);
    return (TRUE);
}


int field(char **p_str, int length, char *val)   /* OBSOLETE */
{
    if (len(*p_str) < length) return (FALSE);
    nstrcpy(val, *p_str, length);
    *p_str += length;
    despace(val);
    return (TRUE);
}


/**** see strlib.h for info on irange(), lrange(), and rrange() ****/

int
irange(int *p_var, int min_val, int max_val) {
    if (*p_var < min_val) {
        *p_var = min_val;
        return (FALSE);
    }
    else if (*p_var > max_val) {
        *p_var = max_val;
        return (FALSE);
    }
    else return (TRUE);
}


int
lrange(long *p_var, long min_val, long max_val) {
    if (*p_var < min_val) {
        *p_var = min_val;
        return (FALSE);
    }
    else if (*p_var > max_val) {
        *p_var = max_val;
        return (FALSE);
    }
    else return (TRUE);
}


int
rrange(real *p_var, real min_val, real max_val) {
    if (*p_var < min_val) {
        *p_var = min_val;
        return (FALSE);
    }
    else if (*p_var > max_val) {
        *p_var = max_val;
        return (FALSE);
    }
    else return (TRUE);
}


char *
binary(
        int num,
        int bits,
        char *str  /* side effected */
) {
    int i;

    for (i = bits; i >= 0; i--, num >>= 1)
        if (num & 1) str[i] = '1'; else str[i] = '0';
    str[bits] = '\0';
    return (str);
}


static char **tempstr;        /* internal use only */
static int tempstr_n;         /* internal use only */

char *
get_temp_string(void) {
    char *str;

    if (tempstr_n == NUM_TEMP_STRINGS - 1) str = tempstr[tempstr_n = 0];
    else str = tempstr[++tempstr_n];
    str[0] = '\0';
    return (str);
}


/***************** The print_info struct used by rs(), etc. *****************/

typedef struct {
    int total_spaces;
    int decimal_places;
    int leading_spaces;
    real max_neg_number;
    char *pr_format_str;
    char *prn_pos_format_str;
    char *prn_neg_format_str;
    char *prd_format_str;
} PRINT_INFO;

PRINT_INFO **print_info; /* [(int)(format*10)] => ptr to a PRINT_INFO */


static PRINT_INFO *
lookup_print_info(real format_num) {
    int int_format;
    PRINT_INFO *p;

    int_format = (int) (format_num * 10 + 0.25);
    if (int_format > 99 || int_format < 0) send(CRASH); /* illegal format */
    p = print_info[int_format];
    if (p == NULL) send(CRASH); /* illegal format */
    return (p);
}


void
init_print_info(void) {
    int int_format, total_spaces, decimal_places;
    PRINT_INFO *p;

    array(print_info, 100, PRINT_INFO*);
    int_format = 0;

    for (total_spaces = 0; total_spaces <= 9; total_spaces++) {
        for (decimal_places = 0; decimal_places <= 9; decimal_places++) {
            if (total_spaces < 1 || decimal_places > total_spaces - 2) {
                print_info[int_format++] = NULL;
                continue;
            }
            single(print_info[int_format], PRINT_INFO);
            p = print_info[int_format++];
            array(p->pr_format_str, 9, char);
            array(p->prn_pos_format_str, 9, char);
            array(p->prd_format_str, 9, char);

            p->total_spaces = total_spaces;
            p->decimal_places = decimal_places;
            p->leading_spaces = total_spaces - decimal_places - 1;
            p->max_neg_number = -0.5 * 1.0 / exp10((real) decimal_places);

            sprintf(p->pr_format_str, "%%-%d.%dlf", total_spaces, decimal_places);
            sprintf(p->prn_pos_format_str, " %%-%d.%dlf",
                    total_spaces - 1, decimal_places);
            p->prn_neg_format_str = p->pr_format_str;
            sprintf(p->prd_format_str, "%%%d.%dlf", total_spaces, decimal_places);
        }
    }
}


/**************** The funky real-number printing routines ********************/

char *
rs(real format, real number) {
    int i, max_spaces;
    char *str;
    PRINT_INFO *info;

    str = get_temp_string();
    info = lookup_print_info(format);

    if (number >= 0.0 || number <= info->max_neg_number)
        sprintf(str, info->pr_format_str, number);
    else sprintf(str, info->pr_format_str, 0.0); /* to avoid printing -0.0 */

    if ((i = strfinder(str, '.')) >= (max_spaces = info->total_spaces) - 1) {
        if (i > max_spaces) /* the int portion is too big */
            for (i = 0; i < max_spaces; i++) str[i] = '*';
        else if (i == max_spaces - 1) str[i] = ' '; /* don't want  "123." */
        /* else if i==max_spaces then '.' will be erased below */
    }
    str[max_spaces] = '\0'; /* Cuts string to proper length, just in case */
    return (str);
}


char *
rsn(real format, real number) {
    int i, max_spaces;
    char *str;
    PRINT_INFO *info;

    str = get_temp_string();
    info = lookup_print_info(format);

    if (number >= 0.0)
        sprintf(str, info->prn_pos_format_str, number);
    else if (number <= info->max_neg_number)
        sprintf(str, info->prn_neg_format_str, number);
    else
        sprintf(str, info->prn_pos_format_str, 0.0); /* to avoid printing -0.0 */

    if ((i = strfinder(str, '.')) >= (max_spaces = info->total_spaces) - 1) {
        if (i > max_spaces) /* the int portion is too big */
            for (i = 0; i < max_spaces; i++) str[i] = '*';
        else if (i == max_spaces - 1) str[i] = ' '; /* don't want  "123." */
        /* else if i==max_spaces then '.' will be erased below */
    }
    str[max_spaces] = '\0'; /* Cuts string to proper length, just in case */
    return (str);
}


char *
rsd(real format, real number) {
    int i, max_spaces;
    char *str;
    PRINT_INFO *info;

    str = get_temp_string();
    info = lookup_print_info(format);

    if (number >= 0.0 || number <= info->max_neg_number)
        sprintf(str, info->prd_format_str, number);
    else
        sprintf(str, info->prd_format_str, 0.0); /* to avoid printing -0.0 */

    if ((i = strfinder(str, '.')) >= (max_spaces = info->total_spaces) - 1) {
        if (i > max_spaces) /* the int portion is too big */
            for (i = 0; i < max_spaces; i++) str[i] = '*';
        else if (i == max_spaces - 1) str[i] = ' '; /* don't want  "123." */
        /* else if i==max_spaces then '.' will be erased below */
    }
    str[max_spaces] = '\0'; /* Cuts string to proper length, just in case */
    return (str);
}


void
str_init(void) {
    array(tok, TOKLEN + 1, char);
    tok[0] = '\0';
    array(x, TOKLEN + 1, char);
    x[0] = '\0';
    self_delimiting = mkstrcpy("");
    null_string = mkstrcpy("");

    matrix(tempstr, NUM_TEMP_STRINGS, TEMP_STRING_LEN, char);
    tempstr_n = 0;
    init_print_info();
}



