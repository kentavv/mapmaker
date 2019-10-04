/******************************************************************************

    #     ####   #          #    #####           #    #
    #    #    #  #          #    #    #          #    #
    #    #    #  #          #    #####           ######
    #    #    #  #          #    #    #   ###    #    #
    #    #    #  #          #    #    #   ###    #    #
    #     ####   ######     #    #####    ###    #    #

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

/************* IOLIB.H - terminal and text file I/O functions *************
OVERVIEW:

Terminal I/O is in general assumed to be "glass-TTY" type line
oriented blocking I/O with most relevent pre-processing handled by the
operating system or this library. (In other words, line editing
[delete keys and and so forth], buffering, echoing characters, and
type-ahead are handled.) Input will not be recieved until an entire
line has been read, and calls to input functions will hang until a
line is ready. Assume that 78 columns of output are possible, that output
lines scroll off the screen, and that at least 24 lines can fit on a
screen. Print()ed strings may be buffered until a '\n', flush(), or
input()/getln() call. These functions should be reasonably well-behaved 
if the program's input or output has been redirected, piped, etc. 

The text file I/O routines operate very similarly. Random access or
binary files are not in general a good idea as both the code for such
files and the files themselves are largely unportable. It is also not
recomended that anything but a disk file be opened using the functions
here (eg do not do open_file("/dev/tty",READ)).

Fancy output modes, including "photo-ing", MORE-processing, scrollback
and input/output redirection are provided by this library. These are
implemented in such a way that most application code can almost
completely ignore them: if the library I/O routines are used
exclusively (eg print()), something reasonable will happen whan thay
are called. The few routines which are exceptions return a boolean on
success or failure which can be tested. 

Even fancier stuff, like better line editing, split screen, and
asychronous scrollback has been or will be hacked into this library,
using curses, X-windows, or the Macintosh toolbox. This now works
under Unix, and VMS, and should work on PCs and Macs soon.

The print string ps may be used to hold one screen line for output
using sprintf() and print() [or pr()]. Always print ps IMMEDIATELY
after it is set, otherwise it may get side-effected in your own
procedure calls!  Ps is large (2K chars) for paranoia's sake, as
sprintf() can spew out some very large strings of characterss in error
conditions. Assume its really only one (or maybe 2-3) lines. ps_ is similar, 
but for internal use only by the functions in this library.

Ln is a global input string side-effected by getln() and fgetln().
The pointer ln itself may be side-effected, (eg: by the token parsing
routines) as it is reset during each call to getln() or fgetln().
However, the string space used to hold the input line is reused, and
any pointers into the ln string which you saved may point to trash
after one of these calls. Ln is alloced for 2K chars, leaving lots of
room for paranoia (again, assume it's much smaller).

Note that the print() and fprint() functions do an awful lot of
processing on the output string, including trying to wrap lines at
word breaks, converting tabs to spaces (every 8 chars), buffering
output (perhaps multiple lines!) and more. The do_fwrite() function does
not do this. Similarly, the getln() function despace()es,
lowercase()es, and filter()s its input, while fgetln(), input(), and
finput() functions do not [yes, it's a somewhat odd inconsistancy, but
it's handy].

It is recomended that only the I/O functions mentioned in here be
used, and that most other I/O functions never be used directly. This
is both because these functions have been tailored for portability and
robustness (crash-proofing), and because they provide a powerful
line-oriented interface which takes as much advantage of the operating
system as it can. Certainly never intermix code that uses these with
code that uses Unix-ish I/O directly!

Iocheck() someday may be a function which should be called
occasionally to check for the break key, pending output for MORE, etc
on computers which don't support asynchronous terminal I/O. This is
best done inside some loop in any computationally ugly piece of code
such that it gets called at least once every few seconds. Iocheck()
returns TRUE if an entire line of input is waiting to be read with
input() or getln(). Do not place iocheck() calls inside
computationally critical loops, as iocheck() may take some time itself
(though likely doesn't).  NOTE: iocheck() IS NOT YET IMPLEMENTED - it
exists but does nothing, and it may never...
***************************************************************************/


extern char *ps, *ln;  /* these global strings are malloced by io_init() */
void iocheck();   /* no args; CURRENTLY A NOP */


/************************ Terminal output routines ***********************/ 
void print(const char *);	 /* args: char *string; does lots of processing */
#define nl() print("\n") 
#define pr() print(ps) 
void flush();	 /* no args. forces everything print()ed to be output */ 

bool temp_print(); /* args: char *str; The string is printed out and flush()ed,
   although it will be erased upon the next call to print(), input(), getln(),
   etc. A following call to temp_print() will erase the previous string
   temp_print()ed, and temp_print(NULL) or temp_print("") explicitly
   clears the last temp_print()ed string. Strings for temp_print may not 
   wrap, and may not contain tabs or newlines. Return FALSE if this is not
   possible on the user's terminal type. */

/* Nice screen control functions - the boolean functions do their stuff and 
   return TRUE if the operation was possible, or return FALSE otherwise. */

bool boing(); /* make a noise, if possible */

bool highlight();  /* args: bool val; sets highlighting (reverse video) mode-
   NOTE: this mode may or may not be effective for only one screen line of 
   output! This needs to be handled as yet... */

bool maybe_clear_screen();  /* no args; clears the screen, if it can,
   although not if the terminal is believed to support scrollback or
   not if curses is enabled. Returns TRUE if it did, FALSE if not. */

bool clear_screen(); /* no args; does it anyway and returns TRUE if possible */

bool to_column(); /* args: int num; (leftmost=0), returns FALSE if the cursor 
   is already right of the specified column */ 
int at_column(); /* no args; returns the current column (leftmost=0) */
void space();     /* args: int n; prints n spaces, regardless */

void do_more();       /* no args: does a "Hit return for more" thing */
extern int tty_lines; /* You may look at this, but please do not set it. */



/********************* Terminal input routines ****************************/
void input();   /* args: char *prompt, *string; int max_input_chars; 
   String must have room for max_input_chars+1 chars, and only
   max_input_chars-1 chars can usually be read, as a '\n' may be read in
   at the end and then deleted from the string. */

void getln();   /* args: char *prompt; side-effects global string ln */

/* NOTE: input() filter()s non-printing characters from the returned
string, but otherwise returns it verbatim, while getln() filter()s,
lowercase()s, and despace()s it! */

bool redirect_input(); /* args: char *filename; bool verbose; */
/* Semantics are kind of like photo_to_file(). Begin taking standard
   input from the specified file instead of where-ever it was comming
   from.  Return FALSE if we fail (e.g. can't open file), in which case
   nothing is changed. Note that input redirections may be layered inside
   one-another, up to some reasonable depth. Use redirect_input(NULL) to
   bag ALL redirects.  Redirect_input(TTY) should redirect to the
   terminal (THIS IS NOT IMPLEMENTED YET!). */ 
#define TTY ""

/* You may look at, but not set, these variables... */
#define redirecting_input (redirs>0)
extern bool interactive;        /* TRUE if input is currently a terminal */


/********************* Fancy Terminal I/O Modes **************************/ 
extern bool more, more_mode, ignore_eof;
/* User code could, but probably shouldn't set these variables.
   Generally, more=FALSE, except in hold() blocks when more_mode=TRUE
   Ignore_eof is on whether the standard input is a terminal or not. 
   All are set by tty_init() and get_cmd_line_args() */

bool do_hold(bool start, bool more_on);

#define hold(with_more_on) \
  for (do_hold(TRUE,with_more_on); holding>0; do_hold(FALSE,FALSE))
extern int holding;
/* This is the preferred interface to more, clearing screen, etc */
#define unhold() do_hold(FALSE,FALSE)

bool photo_to_file(); /* args: char *name, *mode; Returns TRUE if it succeeds, 
   FALSE if it fails (e.g. unable to open file). If a photo file is already 
   open, it is closed and a new one is opened. This function will not affect 
   the current photo mode if it fails. Use name==NULL to turn photo-ing off. */

#define log_open (photo_file[0]!='\0') /* Is photo-ing enabled? This will
   return TRUE even if photo-ing has been temporarily disabled by 
   temp_logging(). To see if photoing is actually happening, look at
   the logging variable. */

extern char *photo_file; /* valid only if log_open is TRUE */
extern bool logging;     /* you may examine this, but don't set it directly */

/* Here are some new functions to temporarily change fancy modes. These each 
   require a pointer to an int to hold the previous state. For example, to 
   temporarily inhibit logging to a file, use: 

   	temp_logging(FALSE,&save);
	...
	prev_logging(save);

   If logging is already off, these calls will have no effect on the code.
   DO NOT do anything silly to set these modes while inside of a temp mode 
   change, as the results are unpredictable! For example, do not call 
   photo_to_file() inside calls to temp_logging() and prev_logging();
   Also, do not set the variable more inside calls to temp_more_mode() and 
   prev_more_mode(). In general, this means make sure that these calls do 
   not span any code which would toggle these things! */

bool temp_logging();      /* args logging, *save_state; turn logging on/off: 
   To turn it on, a log-file must already be opened with photo_to_file()!
   If not (the only failure case), FALSE is returned, otherwise TRUE is. 
   However, even in this case, the routine behaves well (that is, you can 
   ignore its return value if you don't care). */
void prev_logging();      /* args state; restore logging */

bool temp_more_mode(); /* args new_more_mode,*save_state; turn more on/off */
void prev_more_mode(); /* args state; */

void review_memory(); /* flips back through previous output in a nice way */


/******* Curses and WIMP (Windows, Icons, Mouse, and Pointers) Support *******/

/* Use one (and only one) of these functions after custom_lib_init()
   to enable and control really fancy screen management (see system.h for
   details). These may invoke customized code for really funky WIMP stuff. */

bool screen_init();  /* args: int *argc_ptr; char *argv[]; both may be side-
   effected. Sets up funky screen I/O if possible - this may clear the
   screen if it is appropriate to do so, among other things. This returns
   TRUE if screen management (either curses or wimp) was enabled. */

bool split_screen_init(); 
/* args: int *argc_ptr; char *argv[]; int top_lines; void (*update_func)(); 
         both argc_ptr and argv may be side-effected. 
	 update_func() gets args: char **text; int top_lines, top_columns;

   Like screen_init(), except the window is broken into two parts - a
   top half which may be used to display status information, and a bottom
   half which scrolls like an ordinary terminal. The top half will have
   the specified # of lines, and will either be in reverse video or will
   be separated from the bottom half somehow. The procedure
   update_func may be called asynchronously, and should side-effect
   the matrix of characters (text) to reflect what the top part should
   display. */

bool update_top();  /* Forces the top part of a split_screen to be updated. 
   FALSE is returned if split_screen() has not successfully set things up. */
void screen_end();  /* no args: undoes screen_init() or split_screen() */

bool string_editor(); /* args: char *prompt, *str; int max_num_chars; */


/**************************** Text file routines  ****************************/
/* IMPORTANT: It is NOT OK to use the text file I/O routines on terminals or
   other such devices: these routines are geared for real life text disk
   files! (The obvious exception is logging and input redirection, for which
   you should use the rather carefuly crafted routines described above!) */

bool make_filename(); /* args: char *name; int mode; char *extension; 
   This turns name into a perfectly valid file name, and returns TRUE if this
   was possible. The name may include a directory specification - otherwise it
   will be opened in the current working directory. name must be able to hold
   PATH_LENGTH+1 chars. The argument mode should be one of... */
#define DEFAULT_EXTENSION  0  /* if name has no ext, one is added */
#define FORCE_EXTENSION	   1  /* even if name has an ext, it's changed */

bool make_filename_in_dir(); 
/* args: char *name; int ext_mode; char *ext; int dir_mode; char *dir; 
   Like make_filename, but you can also specify a directory. You may use: */
#define FORCE_DIR	   0
#define DEFAULT_DIR  	   1
#define HOME_DIR           "*h*"
#define CODE_DIR           "*c*"
#define CURRENT_DIR        NULL

FILE *open_file();  /* args: char *name, *mode; on failure sends CANTOPEN */
#define WRITE  "w"
#define READ   "r"
#define APPEND "a"
void close_file();  /* args: FILE *fp; */

bool end_of_file(); /* args: FILE *fp; returns TRUE if finput() or fgetln()  
   will be able to grab a new line from the file. */
bool end_of_text(); /* args: FILE *fp; Like an end-of-file test, except it also
  returns TRUE if rest of the file is white. As a side-effect, it will step
  through the file up to the first non-white char. */

void fprint();    /* args: FILE *fp; char *string; processes output */
void do_fwrite(); /* args: FILE *fp; char *string; no processing happens */
#define fwrite(fp,str) (do_fwrite(fp,str))  /* UNIX already has an fwrite() */

#define fpr(fp) fprint(fp,ps)
#define fwp(fp) fwrite(fp,ps)
#define fnl(fp) fwrite(fp,"\n")

void finput();  /* args: FILE *fp; char *str; int max_input_chars; */
void fgetln();  /* args: FILE *fp; side-effects global char *ln; 
   Finput() and fgetln() both return a filter()ed line, and on end-of-file, the
   ENDOFILE message is sent. Str must have room for max_input_chars+1 chars,
   and only max_input_chars-1 chars can usually be read, as a '\n' may be read
   in at the end and then deleted from the string. */
   
void fgetdataln(); /* args: FILE *fp; int *count; side-effects global ln;
   Like fgetln(), although this skips null (white) and comment lines (those
   beginning with a '#' in the leftmost position. Also, each time any line is
   read from the file (data, null, or comment), *count is incremented. */
		     
#define frewind(fp) fseek(fp,0L,0)
#define fflush(fp)  do_fflush(fp)  /* redeclare the C library function */
void do_fflush(); /* never call this directly */

bool rename_file(); 
/* args: char *old, new; in syscode.c renames old file with new name */

bool fgoto_line(); /* in syscode.c - replacement for fseek() */
/* args: FILE *fp; long index; */


/***** Local declarations provided for syscode.c: for internal use only! *****/

extern char *out_file, *in_file, *photo_file;
extern char out_modechar, in_modechar, photo_modechar; 
extern FILE *in, *out, *photo;

#define MAX_FILE_ARGS  3
#define LOAD_FILE_ARG  0
#define RUN_FILE_ARG   1
#define PHOTO_FILE_ARG 2
extern char **file_arg;
extern int num_file_args, dos_output, prep_it, append_it;

extern char *ps_, *ln_;     /* input and output strings for the library only */
extern char *linebuf;       /* tty output buffer */
extern char *gotln, *lnptr; /* input line for getln(), fgetln() */

extern FILE **in_fp;   /* the "stack" of input files for redirecting input */
extern int redirs;
#define MAXREDIRECTS 	4
#define MAXFILES 	10

/* a list associating fp's with names etc */
typedef struct { FILE *fp; char modechar; char *name; } FILE_INFO;
extern FILE_INFO **files;
int lookup_fp();
void ioerror();

void tty_init();  /* no args: in syscode.c, called by lib_init() */
void io_init();   /* no args: in iolib.c, called by lib_init() */
void tty_hello(); /* no args */

extern int cursor, buf, printed_lines, supress_more, lines_written;
extern int term, tty_lines;
extern bool screen, scrollback;
extern int tty_errors, file_errors, puts_errors;

#define MAXLINE  2050
#define LINE     79
#define LASTTAB  72

#define TERM_UNKNOWN      0 /* term types */
#define ANSI              1
#define HP_TERM           2
#define CURSES            8
#define SCROLLING_ANSI    3
#define NONSCROLLING_ANSI 4
#define PC_CONSOLE        5
#define MAC_WINDOW        6
#define WIMP              7

bool boing();
bool do_clear_screen();
bool do_highlight();
bool do_delete_previous_line();
bool do_cursor_left(); /* args int spaces; char *str_to_then_print; */
#define FAR_LEFT -1
bool check_tty_lines();

void lib_puts();
bool tty_gets();
bool file_gets();
void ioerror();
void flush_and_force_nl();
extern bool supress_more, more_break_pending;

extern bool curses, split;
extern bool tried_curses, tried_split, have_drawn_top;
/* These functions exist if HAVE_CURSES is defined. */
bool curses_init();
bool curses_split();
void curses_end();
void curses_error(); 
void curses_clr_scrn();
void curses_del_prev_ln();
void curses_del_this_ln();
void curses_del_to_eol();
void curses_goto_far_left();
void curses_cursor_left();
void curses_set_highlight();
bool curses_gets();
void curses_puts();
void curses_scrollup();
void curses_flush();
void curses_boing();
void curses_refresh();
void curses_update_top();
void curses_draw_top();

extern bool wimp, tried_wimp;
/* These functions exist if HAVE_WIMP is defined. */
void do_text_wimp_init();
void do_split_wimp_init();
void do_custom_wimp_init(); /* Only MAY exist... */

extern bool use_gnu_readline;
bool do_gnu_readline();
bool do_gnu_edit();
bool gnu_copyright();

void edit_line(); /* make real decl above */

