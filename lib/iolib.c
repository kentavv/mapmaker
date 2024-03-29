/******************************************************************************

    #     ####   #          #    #####            ####
    #    #    #  #          #    #    #          #    #
    #    #    #  #          #    #####           #
    #    #    #  #          #    #    #   ###    #
    #    #    #  #          #    #    #   ###    #    #
    #     ####   ######     #    #####    ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#include "system.h"
#include "shell.h"
#include "iolib.h"

static void dump_memory_lines(int start, int num) /* internal */;

/***** globals- see descriptions in iolib.h *****/
char *ps, *ps_, *ln, *ln_;
bool logging, more_mode, more, ignore_eof, scrollback;
int tty_lines, term;
bool screen, interactive;

/* Local declarations - externed in iolib.h only for syscode.c  */
char *out_file, *in_file;
char *photo_file;
FILE *in, *out, *photo;

char *linebuf;  /* stdout buffer: one screen line max */
int cursor, buf, printed_lines;
bool more_break_pending, memorize;
int lines_written;
char *lnptr;

int temp_print_spaces;
#define temp_print_up (temp_print_spaces>0)

char **memory, **held;
int memory_end, memory_wrapped, holding, hold_count;

FILE **in_fp; /* the "stack" of input files for redirecting input */
int redirs;

/* a list associating fp's with names, etc */
FILE_INFO **files;


/***** ERROR HANDLING *****/

#define try_filenum(i) if (files[i]->fp==fp) {    \
    if (name!=NULL) *name=files[i]->name;    \
    if (modechar!=NULL) *modechar=files[i]->modechar; \
    return(i); }

int lookup_fp(FILE *fp, char **name, char *modechar /* side-effected if non-null */) {
    if (fp == NULL) {
        if (name != NULL) *name = null_string;
        if (modechar != NULL) *modechar = '?';
        return -1;
    }

    try_filenum(3); /* MUST HAVE 0...MAXFILES-1 OF THESE! */
    try_filenum(4);
    try_filenum(5);
    try_filenum(6);
    try_filenum(7);
    try_filenum(8);
    try_filenum(9);
    try_filenum(0); /* These are in, out, and photo */
    try_filenum(1);
    try_filenum(2);

    if (name != NULL) *name = null_string;
    if (modechar != NULL) *modechar = '?';
    return -1;
}


void ioerror(char *errmsg, FILE *fp, char *ioline) {
    IOERROR_errmsg = errmsg;
    if (ioline != NULL) {
        nstrcpy(IOERROR_linecopy, ioline, 64);
    } else {
        ioline[0] = '\0';
    }
    lookup_fp(fp, &IOERROR_filename, &IOERROR_modechar);
    send(IOERROR);
}


/********** FILE HANDLING STUFF ***********/

bool make_filename_in_dir(char *str, /* str is side-effected: it should be PATH_LENGTH+1 chars long */
                          bool force_ext,
                          char *ext, /* may have the preceeding '.' or not */
                          int add_dir_mode, /* values defined in iolib.h */
                          char *dir /* assume no trailing divider, except VMS dirs, which have [..] */) {
    int end, dir_chars, has_dir, i, j;
    char cwd[PATH_LENGTH + 1];

    if (!make_filename(str, force_ext, ext)) {
        return FALSE;
    }

    end = len(str) - 1;

    for (i = end; i >= 0; i--) {
        if (strin(PATH_DIR_SEPARATORS, str[i])) break;
    }

    /* has_dir= (i<0); --> this is opposite of what it should be MJD 1-21-93 */

    has_dir = (i >= 0);
    dir_chars = i + 1;

    if (has_dir && add_dir_mode == DEFAULT_DIR) return TRUE;
    if (dir == NULL) { /* eg dir==CURRENT_DIR */
        dir = cwd;
        if (!get_directory(dir)) {
            dir[0] = '\0';
            if (!has_dir) {
                return TRUE; /* Nothing we can do, but we're OK */
            }
        }
    } else if (streq(dir, HOME_DIR)) {
        dir = cwd;
        if (!get_home_directory(dir)) {
            dir[0] = '\0';
            if (!has_dir) {
                return TRUE;
            }
        }
    } else if (streq(dir, CODE_DIR)) {
        dir = cwd;
        if (!get_code_directory(dir)) {
            dir[0] = '\0';
            if (!has_dir) {
                return TRUE;
            }
        }
    }

    /* DO WE REALLY WANT count_tokens() HERE? */
    despace(dir);
    if (count_tokens(dir) != 1) {
        return FALSE;
    }
    if (PATH_UPPERCASE) {
        uppercase(dir);
    }
    for (i = 0, j = 0; i <= len(dir); i++)
        if (strin(PATH_DIR_CHARS, dir[i])) {
            dir[j++] = dir[i];
        }

    if (has_dir) strdel(str, dir_chars);
    maxstrcat(dir, PATH_DIR_FILE_INSERT, PATH_LENGTH);
    maxstrins(str, dir, PATH_LENGTH);
    return TRUE;
}

#define NO_DOT (-1)

bool make_filename(char *str, /* str is side-effected: it should be PATH_LENGTH+1 chars long */
                   bool force_ext,
                   char *ext_ /* may have the preceeding '.' or not */) {
    int i, j, first_dot, next_dot, last_dot, root_start, end;
    bool has_ext, two_ext;
    char *ext_buf = strdup(ext_);
    char *ext = ext_buf;

    if (nullstr(str) || ext == NULL) {
        send(CRASH);
    }

    /* Elimiate extra whitespace, despace, truncate, filter, etc. */
    despace(str);
    truncstr(str, PATH_LENGTH);
    if (PATH_UPPERCASE) {
        uppercase(str);
    }
    if (count_tokens(str) != 1) {
        free(ext_buf);
        return FALSE;
    }
    for (i = 0, j = 0; i <= len(str); i++) {
        if (strin(PATH_OK_CHARS, str[i])) str[j++] = str[i];
    }

    /* Also for the extension... */
    despace(ext);
    if (ext[0] == '.') {
        ext++;
    }
    if (PATH_UPPERCASE) {
        uppercase(ext);
    }
    if (count_tokens(ext) != 1) {
        free(ext_buf);
        return FALSE;
    }
    for (i = 0, j = 0; i <= len(ext); i++) {
        if (strin(PATH_OK_CHARS, ext[i])) ext[j++] = ext[i];
    }

    /* Find the dot for the extension. If PATH_SINGLE_EXTENSION, find 
       the leftmost dot to the right of a PATH_DIR_SEPARATOR, otherwise 
       find the rightmost dot of a PATH_DIR_SEPARATOR */

    end = len(str) - 1;

    for (i = end; i >= 0; i--)
        if (strin(PATH_DIR_SEPARATORS, str[i])) break;
    root_start = i + 1;
    if (root_start == end + 1) {
        return FALSE;
    }

    for (i = root_start, first_dot = NO_DOT; i <= end; i++) {
        if (str[i] == '.') {
            first_dot = i;
            break;
        }
    }

    for (i = first_dot + 1, next_dot = NO_DOT; i <= end; i++) {
        if (str[i] == '.') {
            next_dot = i;
            break;
        }
    }

    for (i = end, last_dot = NO_DOT; i >= root_start; i--) {
        if (str[i] == '.') {
            last_dot = i;
            break;
        }
    }

    has_ext = first_dot != NO_DOT;
    two_ext = has_ext && next_dot != NO_DOT;

    if (PATH_SINGLE_EXTENSION && two_ext) {
        /* truncate */
        str[next_dot] = '\0';
        last_dot = first_dot;
        next_dot = NO_DOT;
    }

    if (!has_ext && !nullstr(ext)) {
        truncstr(str, PATH_LENGTH - 2); /* then add a dot */
        strcat(str, ".");
        maxstrcat(str + len(str), ext, PATH_LENGTH);
    }

    if (force_ext && has_ext) {
        str[last_dot] = '\0';
        truncstr(str, PATH_LENGTH - 2);
        if (!nullstr(ext)) {
            strcat(str, ".");
        }
        maxstrcat(str + len(str), ext, PATH_LENGTH);
    }

    free(ext_buf);
    return TRUE;
}


FILE *open_file(char *name,  /* It's best to use make_filename() on name first */ char *mode)  /* use a #define in iolib.h for mode */ {
    FILE *fp;
    int i, n;

    for (i = 3, n = -1; i < MAXFILES; i++) {
        if (files[i]->fp == NULL) n = i;
    }
    if (n == -1) {
        ioerror("can't open: too many files", (FILE *) NULL, name);
    }

    if ((fp = fopen(name, mode)) == NULL) { /* Can't open */
        nstrcpy(CANTOPEN_path, name, PATH_LENGTH);
        CANTOPEN_modechar = mode[0];
        send(CANTOPEN);
    }

    files[n]->fp = fp;
    nstrcpy(files[n]->name, name, PATH_LENGTH);
    files[n]->modechar = mode[0];
    return fp;
}


void close_file(FILE *fp) {
    /* ignore NULL fp */
    int n, i;

    if (fp == NULL) return;

    for (i = 3, n = -1; i < MAXFILES; i++) if (files[i]->fp == fp) n = i;
    if (n == -1) send(CANTCLOSE);
    if (!xclose(fp)) send(CANTCLOSE);
    files[n]->fp = NULL;
    strcpy(files[n]->name, "");
}


void do_fwrite(FILE *fp, char *str_) {
    char *str = strdup(str_);
    _filter(str);
    lib_puts(fp, str);
    free(str);
} /* DO NOTHING FANCY */


void fprint(FILE *fp, char *str_) {
    char *str = strdup(str_);
    _filter(str);
    lib_puts(fp, str);
    free(str);
} /* JUST FOR NOW - GET FANCIER LATER */


void finput(FILE *fp, char *str, int length) {
    if (fp == stdin) {
        input("? ", str, length);
        return;
    }
    if (!file_gets(fp, str, length)) {
        send(ENDOFILE);
    }
    _filter(str);
}


void fgetln(FILE *fp) {
    /* ln is side-effected. */
    if (fp == stdin) {
        getln("? ");
        return;
    }
    ln = lnptr;
    if (!file_gets(fp, ln, MAXLINE)) send(ENDOFILE);
    _filter(ln);
}


void fgetdataln(FILE *fp, int *count) {
    /* ln is side-effected */
    char *p;

    if (fp == stdin) {
        getln("? ");
        if (count != NULL) ++*count;
        return;
    }
    ln = lnptr;
    do {
        if (!file_gets(fp, ln, MAXLINE)) {
            send(ENDOFILE);
        }
        if (count != NULL) {
            ++*count;
        }
        for (p = ln; white(*p) || trash(*p); ++p) {
        }
    } while (*p == '\0' || *p == '#');
    _filter(ln);
}


bool end_of_text(FILE *fp) {
    int c;

    while ((c = getc(fp)) != EOF) {
        if (!white(c) && !trash(c)) {
            ungetc(c, fp);
            return FALSE;
        }
    }

    return TRUE;
}


/*********************** OUTPUT ROUTINES... local ***********************/

void flush_linebuf(void) {
    bool endoline;
    if (buf == 0) {
        return;
    }
    endoline = (linebuf[buf - 1] == '\n');
    linebuf[buf] = '\0';
    buf = 0;

    if (holding) {
        maxstrcat(held[hold_count], linebuf, LINE + 1);
    } else {
        kill_temp_print();
        dump_to_screen(linebuf);
    }

    if (logging) {
        lib_puts(photo, linebuf);
    }
    mem_puts(linebuf);

    if (endoline) {
        cursor = 0;
        if (holding) {
            hold_count++;
            held[hold_count][0] = '\0';
            if (hold_count == MAX_HOLD_LINES || hold_count == tty_lines - 1) {
                if (!dump_held_lines()) {
                    held[0][0] = '\0';
                    hold_count = 0;
                    holding = 0;
                    more = FALSE;
                    more_break_pending = FALSE;
                    send(INTERRUPT);
                }
            }
        }
    }
}


bool dump_held_lines(void) {
/* When one does hold {}, we flush_and_force_nl(TRUE). THUS, held lines will
   always start dumping into the left column. Thus clear_screen works right. */
    int i;

    if (held[hold_count][0] != '\0') {
        hold_count++;
        held[hold_count][len(held[hold_count]) - 1] = '\n';
    }
    if (hold_count == 0) {
        return TRUE;
    }

    kill_temp_print();
    if (more_break_pending) {
        if (!really_do_more()) {
            hold_count = 0;
            return FALSE;
        }
    } else if (more) {
        if (printed_lines > tty_lines / 2 &&
            printed_lines + hold_count >= tty_lines - 2) {
            if (!really_do_more()) {
                hold_count = 0;
                return FALSE;
            }
        }
    }

    /* dump_to_screen may add new 'more' breaks if needed */
    for (i = 0; i <= hold_count; i++) {
        if (!dump_to_screen(held[i])) return FALSE;
    }
    hold_count = 0;
    held[0][0] = '\0';
    return TRUE;
}


bool dump_to_screen(char *str) {
    if (more_break_pending) {
        if (!really_do_more()) {
            return FALSE;
        }
    }
    lib_puts(out, str);

    if (len(str) > 0 && str[len(str) - 1] == '\n') {
        if (++lines_written > 1000000) send(CRASH);
        if (more) {
            if (printed_lines >= tty_lines - 3) more_break_pending = TRUE;
            else printed_lines++;
        } else printed_lines = 0;
    }
    return TRUE;
}


bool really_do_more(void) {
    bool continue_flag = TRUE;
    int i;

    if (!interactive || redirecting_input) return TRUE;

    run {
            do_highlight(TRUE);
            check_tty_lines();
            lib_puts(out, "Hit RETURN for more or Q to quit...");
            fflush(out);
            do_highlight(FALSE);
            if (!tty_gets(ln_, MAXLINE - 10) && !ignore_eof) {
                nl();
                send(ENDOINPUT);
            }
            do_delete_previous_line();
            for (i = 0; ln_[i] != '\0'; i++) {
                if (ln_[i] == 'q' || ln_[i] == 'Q') send(INTERRUPT);
            }
            printed_lines = 0;
            more_break_pending = FALSE;
        } except_when(INTERRUPT) {
        do_delete_previous_line();
        continue_flag = FALSE;
    }
    return continue_flag;
}


void mem_puts(char *str) {
    int last;
    if (str[0] == '\0' || !memorize) return;

    last = len(str) - 1;
    maxstrcat(memory[memory_end], str, LINE + 1);
    if (str[last] == '\n') {
        if (memory_end == MEMORY_LINES - 1) {
            memory_end = 0;
            memory_wrapped = TRUE;
        } else {
            memory_end++;
        }
        memory[memory_end][0] = '\0';
    }
}


void flush_and_force_nl(bool nl_on_screen_also) {
    kill_temp_print();
    flush();
    if (cursor != 0) {
        if (nl_on_screen_also) {
            nl();
            flush();
        } else {
            if (logging) {
                lib_puts(photo, "\n");
            }
            mem_puts("\n");
            cursor = 0;
            if (more) {
                if (printed_lines >= tty_lines - 3) {
                    more_break_pending = TRUE;
                } else {
                    printed_lines++;
                }
            } else {
                printed_lines = 0;
            }
        }
    }
}


void kill_temp_print(void) {
    if (temp_print_up) {
        lib_puts(out, "\n");
        temp_print_spaces = 0;
    }
}


bool lib_clear_screen(void) {
    check_tty_lines();
    if (more_break_pending) {
        if (!really_do_more()) {
            send(INTERRUPT);
        }
    }
    if (do_clear_screen()) {
        printed_lines = -1;
        cursor = 0;
        return TRUE;
    } else {
        return FALSE;
    }
}


/**************** application level output routines ****************/

int photo_to_file(char *new_log_name, /* best to run through make_filename() first */ char *new_log_mode /* use a #define as for open_file() */) {
/* Returns FALSE if it fails, in which case nothing is changed. If 
   nullstr(new_log_name) then photoing is stopped. */
    FILE *fp;
    flush_and_force_nl(TRUE);

    if (nullstr(new_log_name)) {
        if (log_open) fclose(photo);
        logging = FALSE;
        photo_file[0] = '\0';
        return TRUE;
    } else {
        if ((fp = fopen(new_log_name, new_log_mode)) == NULL) {
            return FALSE;
        } else {
            if (log_open) {
                fclose(photo);
            }
            logging = TRUE;
            photo = fp;
            strcpy(photo_file, new_log_name);
            return TRUE;
        }
    }
}


void print(const char *str) {
/* Sends IOERROR if an error occurs */
/* Linebuf should never have an embedded '\n' (only one at the end, maybe) 
   and should only have as many chars as can fit on the rest of the line. 
   Print interprets tabs and newlines, and punts any other control chars. */

    int i, j, k, found;
    char c, save[21];

    kill_temp_print();

    if (str == NULL) {
        send(CRASH);
    }

    for (i = 0; str[i] != '\0'; i++) {
        c = str[i];

        if (c == '\n' || (c == ' ' && cursor > LINE) ||
            (c == '\t' && cursor > LASTTAB)) {
            linebuf[buf++] = '\n';
            flush_linebuf();
        } else if (cursor > LINE) {
            linebuf[buf] = '\0';
            /* Line would wrap, so we find a nice place to break it */
            for (found = FALSE, j = buf - 1, k = (buf > 20 ? buf - 20 : 0); j >= k; j--)
                if (white(linebuf[j])) {
                    found = TRUE;
                    break;
                }
            if (found) {
                strcpy(save, &linebuf[j + 1]);
                linebuf[j++] = '\n';
                buf = j;
                flush_linebuf(); /* now buf==0, cursor==0 */
                for (k = 0; save[k] != '\0'; k++) {
                    linebuf[buf++] = save[k];
                    cursor++;
                }
            } else { /* didn't find a whitespace */
                linebuf[buf++] = '\n';
                flush_linebuf();
            }
            if (!trash(c)) {
                linebuf[buf++] = c;
                cursor++;
            }


        } else if (c == '\t' && TRANSLATE_TABS) {
            do {
                linebuf[buf++] = ' ';
                cursor++;
            } while (cursor % 8 != 0);

        } else if (!trash(c)) { /* KLUDGE: WHAT CTRL CHARS TO HANDLE?*/
            linebuf[buf++] = c;
            cursor++;
        }
    }
}


bool temp_print(char *simple_str, char *fancy_str) {
    char *str;

    /* If the last print command was not temp_print(), we must first get the
       screen up to date... */
    if (!temp_print_up) flush();

    if (term == TERM_UNKNOWN) str = simple_str; else str = fancy_str;

    /* NOTE: we don't need to fflush in here for the do_ functions, as they
       in syscode fflush(out), and the curses functions update()... */

    if (nullstr(str)) { /* kill old temp print thing */
        if (temp_print_up) {
            lib_puts(out, "\n");
            fflush(out);
            temp_print_spaces = 0;
        } /* else do nothing */
        return TRUE;

    } else if (!temp_print_up) { /* new temp print thing */
        /* do_highlight does nothing if TERM_UNKNOWN */
        if (term == TERM_UNKNOWN) {
            lib_puts(out, str);
            fflush(out);
        } else {
            do_highlight(TRUE);
            lib_puts(out, " ");
            lib_puts(out, str);
            do_highlight(FALSE);
        }
        temp_print_spaces = len(str);
        return TRUE;

    } else { /* there is an existing temp print str, which we must deal with */
        if (temp_print_spaces == 79) {
            lib_puts(out, "\n");
            temp_print_spaces = 0;
        }
        lib_puts(out, ".");
        fflush(out);
        temp_print_spaces++;
        return TRUE;
    }
}


bool do_hold(bool start, bool more_on) {
    if (start == TRUE && !holding) {
        flush_and_force_nl(TRUE);
        check_tty_lines();
        hold_count = 0;
        holding++;
        more = more_on;
        printed_lines = 0;
        more_break_pending = FALSE;
    } else if (start == FALSE) {
        check_tty_lines();
        dump_held_lines();
        holding = 0;
        more = FALSE;
    } /* else start==MAYBE */
    return FALSE;
}


void flush(void) {
    dump_held_lines();
    flush_linebuf();
}


void space(int n) {
    while (n-- > 0) print(" ");
}


bool to_column(int num) {
    if (num > LINE - 1) {
        nl();
        return TRUE;
    }
    if (num < cursor) return FALSE;
    while (cursor < num) {
        linebuf[buf++] = ' ';
        cursor++;
    }
    return TRUE;
}


int at_column(void) {
    return cursor;
}


bool maybe_clear_screen(void) {
    flush_and_force_nl(FALSE);
    if ((!scrollback) && lib_clear_screen()) return TRUE;
    else return FALSE;
}


void review_memory(void) {
    int i, old_photo, start, num, old_memory;

    if (memory_end == 0 && !memory_wrapped) {
        return;
    }
    if (more_break_pending) {
        if (!really_do_more()) {
            send(INTERRUPT);
        }
    }
    flush_and_force_nl(TRUE);

    temp_logging(FALSE, &old_photo);
    do_hold(TRUE, TRUE);
    old_memory = memorize;
    memorize = FALSE;

    start = i = (memory_wrapped ? memory_end + 1 : 0); /* start */
    num = 0;
    while (TRUE) { /* find next group of lines */
        if (i == memory_end) {
            break; /* end - was m_e-1 */
        }
        if (i == MEMORY_LINES - 1) {
            i = 0;
        } else {
            i++;
        }
        num++;
        if (num == tty_lines - 1) {
            dump_memory_lines(start, num);
            start = i;
            num = 0;
        }
    }
    dump_memory_lines(start, num);

    do_hold(FALSE, FALSE);
    prev_logging(old_photo);
    memorize = old_memory;
}


static void dump_memory_lines(int start, int num) {
    /* internal */
    int j;

    for (j = start; num > 0; num--) {
        print(memory[j]);
        if (j == MEMORY_LINES - 1) {
            j = 0;
        } else {
            j++;
        }
    }
}


/*** INPUT ROUTINES ***/

int redirect_input(char *fname, bool verbose) {
    /* fname=NULL to interrupt */
    bool retoin = FALSE;
    FILE *fp;

    if (fname == NULL) {
        if (redirs != 0 && verbose) {
            print("\n\n\t...Input file interrupted...\n");
        }
        while (redirs > 0) {
            close_file(in_fp[redirs--]);
        }
        return TRUE;
    }
    if (redirs >= MAXREDIRECTS) {
        if (verbose) {
            print("error: Too many open input files\n");
        }
        return FALSE;
    }
    run {
            fp = open_file(fname, READ);
            in_fp[++redirs] = fp;
            retoin = TRUE;
        } except_when(CANTOPEN) retoin = FALSE;

    if (verbose) {
        if (!retoin) {
            sprintf(ps, "error: Unable to open input file '%s'\n", fname);
        } else {
            sprintf(ps, "\n\t...Running commands from input file '%s'...\n", fname);
        }
        pr();
    }
    return retoin;
}


/* Input(): If the prompt wraps (over, say, the 75th column), everything 
   may become confused!  Str must be at least length+2 long to hold
   "\n\0", however, the trailing '\n' is stripped from the returned
   string, which is also despace()ed and filter()ed. */

void input(char *prompt, char *str, int length) {
    bool eof = FALSE;
    more_break_pending = FALSE;

    if (redirs == 0) { /* it's not redirected via redirect_input() */
        flush();
        eof = !do_gnu_readline(prompt, str, length);  /* see syscode.c */
        if (eof) {
            str[0] = '\0';
            if (interactive && ignore_eof) {
                send(INTERRUPT);
            } else {
                send(ENDOINPUT);
            }
        } else {
            _filter(str);
            despace(str);
            if (use_gnu_readline) {
                mem_puts(prompt);
                if (logging) {
                    lib_puts(photo, prompt);
                }
            }
            mem_puts(str);
            mem_puts("\n");
            if (logging) {
                lib_puts(photo, str);
                lib_puts(photo, "\n");
            }
        }
    } else { /* we have redirected input */
        print(prompt);
        if (!file_gets(in_fp[redirs], str, length)) {
            /* if (cursor!=0) nl(); */
            print("\n\t...end of input file...\n\n");  /* move to shell? */
            close_file(in_fp[redirs--]);
            input(prompt, str, length); /* try again */
        }
        _filter(str);
        despace(str);
        print(str);
        nl();
    }
}


/* edit_line(): essentially the same rules as input() */

void edit_line(char *prompt, char *str, int length, char *initial) {
    bool eof;
    more_break_pending = FALSE;
    eof = FALSE;

    if (interactive && redirs == 0 && use_gnu_readline) {
        flush();
        eof = !do_gnu_edit(prompt, str, length, initial);
        if (eof) {
            str[0] = '\0';
            if (ignore_eof) {
                send(INTERRUPT);
            } else {
                send(ENDOINPUT);
            }
        } else {
            _filter(str);
            despace(str);
            mem_puts(prompt);
            mem_puts(str);
            mem_puts("\n");
            if (logging) {
                lib_puts(photo, prompt);
                lib_puts(photo, str);
                lib_puts(photo, "\n");
            }
        }
    } else { /* we have redirected input */
        error("Unable to edit sequence from file input");
    }
}


void getln(char *prompt /* may signal IOERROR or ENDOINPUT */) /* ln is side-effected, is filtered, despaced, lowercased */ {
    ln = lnptr;
    input(prompt, ln, MAXLINE - 2);
    lowercase(ln);
}


/*** TEMPORARILY CHANGE TERMINAL I/O MODES ***/

bool temp_logging(bool new, bool *old) {
    if (new && !log_open) return FALSE;
    flush();
    if (cursor != 0) {
        nl();
    }
    *old = logging;
    logging = new;
    return TRUE;
}


void prev_logging(int old) {
    if (old && !log_open) {
        send(CRASH);
    }
    flush();
    if (cursor != 0) {
        nl();
    }
    logging = old;
}


void io_init(void)  /* this is run from lib_init() */ {
    int i;

    array(ps, MAXLINE + 3, char);
    array(ps_, MAXLINE + 3, char);
    array(linebuf, LINE + 3, char);
    array(ln, MAXLINE + 3, char);
    lnptr = ln;
    array(ln_, MAXLINE + 3, char);

    parray(files, MAXFILES, FILE_INFO);
    for (i = 0; i < MAXFILES; i++) {
        array(files[i]->name, PATH_LENGTH + 1, char);
        files[i]->fp = NULL;
        files[i]->modechar = '?';
    }

    in_file = files[0]->name;
    strcpy(in_file, "stdin");
    array(in_fp, MAXREDIRECTS + 1, FILE*);
    redirs = 0;
    for (i = 1; i <= MAXREDIRECTS; i++) {
        in_fp[i] = NULL;
    }
    in_fp[0] = in = files[0]->fp = stdin;

    out_file = files[1]->name;
    strcpy(out_file, "stdout");
    out = files[1]->fp = stdout;

    photo_file = files[2]->name;

    matrix(held, MAX_HOLD_LINES + 1, LINE + 3, char);
    holding = FALSE;
    hold_count = 0;
    for (i = 0; i < MAX_HOLD_LINES + 1; i++) {
        held[i][0] = '\0';
    }

    matrix(memory, MEMORY_LINES, LINE + 3, char);
    memory_end = 0;
    memory_wrapped = FALSE;
    for (i = 0; i < MEMORY_LINES; i++) {
        memory[i][0] = '\0';
    }

    more_break_pending = FALSE;
    memorize = TRUE;
    cursor = 0;
    printed_lines = 0;
    lines_written = 0;
    buf = 0;
    linebuf[0] = '\0';
    holding = 0;
    temp_print_spaces = 0;

    /* These may all be changed by tty_init() or wimp_io_init() */
    more = FALSE;
    more_mode = DEFAULT_MORE_MODE;
    logging = FALSE;
    ignore_eof = TRUE;
    interactive = TRUE;
    use_gnu_readline = FALSE;
    scrollback = FALSE;
    screen = TRUE;
    tty_errors = 0;
    file_errors = 0;
    puts_errors = 0;
    tty_lines = 24;
    term = TERM_UNKNOWN;
}
