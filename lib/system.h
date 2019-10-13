/******************************************************************************

  ####    #   #   ####    #####  ######  #    #          #    #
 #         # #   #          #    #       ##  ##          #    #
  ####      #     ####      #    #####   # ## #          ######
      #     #         #     #    #       #    #   ###    #    #
 #    #     #    #    #     #    #       #    #   ###    #    #
  ####      #     ####      #    ######  #    #   ###    #    #

******************************************************************************/
/* This file is part of MAPMAKER, Copyright 1987-1992 Whitehead Institute.
   See the READ.ME file for license agreement and non-warranty information. */

#include <ctype.h>
#include <errno.h>
#include <malloc.h>
#include <math.h>
#include <setjmp.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

typedef double  real;
typedef int     flag;
typedef int     bool;
#define TRUE 	(1)
#define FALSE	(0)
#define MAYBE	(-1)

#include "eqn.h"
#include "iolib.h"
#include "mathlib.h"
#include "memlib.h"
#include "misclib.h"
#include "msglib.h"
#include "shell.h"
#include "stats.h"
#include "strlib.h"
#include "table.h"


/***************************************************************************
This file, along with the syscode.c file, is used to provide most, if
not all system specific, definitions allowing our library files to
port from one system to another.  Note that these declarations are
provided for the shell library code, NOT NECESSARILY FOR THE APPLICATION
CODE! See the apropriate ...lib.h file for a list of the functions you
should use.

With the increasing encroachment of ANSI C, we can now get away with
only only one system.h files (with #ifdefs, of course) for all
systems. An ANSI compiler doesn't help too much, because Sun/DEC (and
others) include files are just not ANSI. At present, fixes are in
place for:

* SunOS 4.1.x with the bundled cc (ANSI-C, just say no)
* DEC Ultrix 4.2 on the nice but now obsolete DECStation line (cc, not acc)
* Apple's A/UX 3.0 on high end Macintoshes
* WATCOM C/386 9.0 on MicroSloth's MS-DOG, optionally with M$-Windoze 3.1

Previous versions of this library have been ported in various flavors
to Mac ThinkC and VMS as well as other flavors of Unix - hopefully
these ports will happen again for this version soon. If you port this
code, we suggest you use the following flags:

   _SYS_SUNOS	Specifics for SunOS 4.1.x (only tried SPARC)
   _SYS_SOLARIS	Specifics for Solaris 2.x (SPARC and Intel versions identical?)
   _SYS_AUX	Specifics for Apple's A/UX 3.0
   _SYS_MACHTEN	Specifics for Tennon's MachTen on Macintosh (tried & failed!) 
   _SYS_ULTRIX	Specifics for DEC RISC/ULTRIX 4.x (for MIPS based DECStations)
   _SYS_OSF	Specifics for OSF/1 for DEC Alpha (like, if it ever ships)
   _SYS_AIX	Specifics for the inferior but marketable AIX (RISC or other?)
   _SYS_HPUX	Specifics for HP-UX (8.x?) for HP 9000/700 series.
   _SYS_UNIX	Basic Unix semantics, defined if any one of the above is

   _SYS_WATCOM  Specifics for WATCOM C/386 9.0 with apropriate libraries
   _SYS_DOS	Basic DOS semantics, defined if the above is

   _SYS_THINKC  Specifics for ThinkC (5.0?)
   _SYS_MPW     Specifics for MPW C (ver???)
   _SYS_MAC	Basic MAC semantics, should be defined if one of the above is

   Things I'm not quite sure what to do with...
   _SYS_NT	Windoze/NT (Not There) on Intel (what about Alpha/MIPS/etc?)
   _SYS_OS2	OS/2 2.x (gives me the willies just to think about it)
   _SYS_NEXT	NeXTStep 3.0, wiTh ThE fUNkey CaPItaLizAtION (040 or Intel?)
   _SYS_VMS	Basic VMS semantics (what options are there?)

Of course most of these flags are not handled yet. Happy porting :-)

Among other things ,this is the only file which understands which of
the system's #include files need to be used, as well as what order
they should be included in! (This is handled at the end of this file).
You should include standard files by #defining one or more of the
following and then #including "system.h":

   INC_IO      File and TTY I/O routines, including stdio.
   INC_MATH    Simple math routines. 
   INC_MSG     Message sending and trapping routines. 
   INC_MEM     Memory allocation/free routines. 
   INC_STR     String functions, including C's and this library's.
   INC_LIB     All five of the above.

   INC_MISC    Misc. functions (time, subshells, array dumps, sorting, etc).  
   INC_SHELL   The ubiquitous MAPMAKERish shell.
   INC_TABLE   A useful data struct, used by the shell and other things. 
*************************************************************************/

//#ifdef  _SYS_SUNOS
//#define _SYS_UNIX
//#else
//#ifdef  _SYS_ULTRIX
//#define _SYS_UNIX
//#else
//#ifdef  _SYS_AUX
//#define _SYS_UNIX  /* Note: NOT _SYS_MAC */
//#endif
//#endif
//#endif
//#ifdef _SYS_WATCOM
//#define _SYS_DOS
//#endif



/************************ File name syntax ****************************
Setup the following for the OS's path types (these are used by the 
make_filename() procedure in iolib.c). Note that the lengths of various 
elements of a path are not checked, only its total length.
*************************************************************************/
#define PATH_LENGTH		200

//#ifndef _SYS_DOS /* e.g. _SYS_UNIX or some POSIX like thing */
#define HELP_EXT                ".help"
#define ARG_CHAR                '-' /* Usual char for command line switches */
#define PATH_REQUIRE_EXTENSION 	FALSE
#define PATH_SINGLE_EXTENSION  	FALSE
#define PATH_UPPERCASE 		FALSE
#define PATH_DIR_SEPARATORS 	"/"    /* rightmost separating char */
#define PATH_DIR_FILE_INSERT    "/"    /* insert between a dir and filename */
#define PATH_OK_CHARS \
"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_.~#/"
#define PATH_DIR_CHARS \
"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_.~#/"

//#else /* is _SYS_DOS */
//#define HELP_EXT                ".hlp"
//#define ARG_CHAR                '/' /* Usual char for command line switches */
//#define PATH_REQUIRE_EXTENSION 	TRUE
//#define PATH_SINGLE_EXTENSION  	TRUE
//#define PATH_UPPERCASE 		TRUE
//#define PATH_DIR_SEPARATORS 	"\\"    /* rightmost separating char */
//#define PATH_DIR_FILE_INSERT    "\\"    /* insert between a dir and filename */
//#define PATH_OK_CHARS "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_.$#/\\"
//#define PATH_DIR_CHARS "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_.$#/\\"
//#endif


/**************************** Subshells *********************************
To help portability, we use the system() function to spawn subshells
on all operating systems. (However, system specific code may be placed
in syscode.c to do it another way). Under Unix, if TRY_SHELL_VAR is
TRUE, getenv() will be used to get the value of SHELL to run by
calling system(getenv("SHELL")). If the SHELL variable is not set we
try calling system(TRY_SHELL_CMD). Naturally, these presume
HAVE_GETENV is defined.

On MS-DOS we will first try system(getenv("COMSPEC")). Otherwise,
TRY_SHELL_CMD should be "command.com". System() is not supported on the
Mac, and on VMS, SHELL_CMD should be "SPAWN".
**************************************************************************/
/* #define NO_SYSTEM_FUNC *//* ANSI/POSIX - All should have system() now */

//#ifdef _SYS_DOS
//#define TRY_COMSPEC_VAR
//#define TRY_SHELL_CMD "command.com"
//#define SHELL_MESSAGE "Running DOS prompt: Type 'exit' to return...\n"
//
//#else /* _SYS_UNIX */
#define TRY_SHELL_VAR
#define TRY_SHELL_CMD "/bin/csh"  /* csh is the default on Sun, DEC, A/UX */
#define SHELL_MESSAGE \
"Running Unix shell: Type 'exit' or hit Control-D to return...\n"
//#endif


/************************ Error handling routines ************************* 
Errno lookups: strerror(), which is allegedly ANSI standard, should
return an error message (char*) when given a correct system errno
(generally from ferror()). Don't use perror(). For systems without
strerror(), use sys_errlist[] and sys_nerr (this includes Lightspeed,
System V, HP-UX, and SunOS). If sys_nerr does not exist, we ignore it
and hope that sysnerr[errno] is valid... On some systems (Sun)
sys_errlist does not appear in any #include file, although errno
(extern int) is defined (usu. in math.h?).

Also, define HAVE_MATHERR if matherr() trapping is supported, as it
is on System V (See msglib.h). SIGHANDLE should be the return value
type of a signal handling function used as an argument to signal().
***************************************************************************/
#define HAVE_MATHERR  /* Doesn't hurt to not define it... */

//#ifdef _SYS_WATCOM  /* strerror is defined, as it should be */
//#define SIGHANDLE int
//#define SIGQUIT SIGABRT /* bizzare */
//#endif
//#ifdef _SYS_SUNOS /* the @#!$% map pages lie */
//#define SIGHANDLE void
//#define strerror(num) (num<sys_nerr ? sys_errlist[num]:(char*)NULL)
//extern char *sys_errlist[]; extern int sys_nerr;
//#endif
//#ifdef _SYS_AUX
//#define SIGHANDLE void
//#define strerror(num) (num<sys_nerr ? sys_errlist[num]:(char*)NULL)
//extern char *sys_errlist[]; extern int sys_nerr;
//#endif
//#ifdef _SYS_ULTRIX
//#define SIGHANDLE void
//#endif



/// ************************* Random Numbers **********************************
//There are at least 3 standard random number generators in C-libraries:
//rand() (ANSI, fairly ubiquitous), drand48() (System V), and random() (BSD).
//Rand() is a lousy random number generator: don't use it unless no
//other is available. See the code for these in syscode.c.  The limits sont
//seem to be real portable, so be careful.
//***************************************************************************/
//#ifdef _SYS_AUX
//#define USE_DRAND48
//#else
//#ifdef _SYS_SUNOS
//#define USE_RANDOM
//#else
//#ifdef _SYS_ULTRIX
//#define USE_RANDOM
//#else
//#define USE_SRAND
//#endif
//#endif
//#endif


// ***************************** C-Library **************************************
//Here are some casts to make various incompatible C-library functions
//compatible accross machines. NONE OF THESE FUNCTIONS SHOULD BE CALLED
//DIRECTLY! RATHER, THEY SHOULD BE CALLED ONLY BY THE HELPERS ROUTINES
//WHICH DRIVE THEM!  Define the constants appropriately so that the
//following declarations are correct on the target machine:
//
//      CALLOC_PTR_TO *calloc(num,sizeof(...))
//          CALLOC_NUM_TYPE num;
//	  SIZEOF_TYPE size;
//
//      TIME_TYPE time();
//      char *ctime(ptr)
//          TIME_TYPE *ptr;
//
//      qsort(data,array_length,width,compare)  \*** portability hell ***
//          QSORT_DATA_PTR_TO *data; \* an array \*
//          QSORT_LENGTH array_length;
//	  SIZEOF_TYPE width;
//          int (*compare)();
//      int compare(a,b)
//          QSORT_COMPARE_PTR_TO(type) *a, *b;
//
//where type is the type of thing that pointers into data (eg a and b) really
//point to. For example, to compare long integers, the code should read (Note:
//better interfaces are provided in the library - see misclib.h for them!)
//
//      int long_compare(a,b)
//      QSORT_COMPARE_PTR_TO(long) *a;
//      QSORT_COMPARE_PTR_TO(long) *b;
//      {
//          long x,y;
//	  x= *(long*)a; y= *(long*)b;
//	  if (x<y) return(-1); else if (x==y) return(0); else return(1);
//      }
//      ...
//      qsort(QSORT_CAST(data),(QSORT_LENGTH)length,sizeof(long),long_compare)
//
//Qsort()'s return value, if it has one, is ignored. We assume that qsort never
//fails (perhaps this is stupid).
//
//Define (and use!) exit_main as the statement which main() should execute
//when the program normally terminates. abnormal_exit should be an exit()
//call for abnormal situations (mostly for the HP 800 exit() bus-error bug!)
//
//Define HAVE_CHDIR, HAVE_GETCWD, and/or HAVE_GETENV if those functions exist.
//Check out their usage (arguments, etc) in syscode.c to be safe!
//***************************************************************************/
#define HAVE_GETENV   /* Don't all ANSI compilers have these now? */
#define HAVE_GETCWD
#define HAVE_CHDIR
#define normal_exit()     exit(0)
#define abnormal_exit()   exit(1)	
#define exit_main()       return(0)
//#define TIME_TYPE 	  time_t /* for ctime - yes even on SunOS! */
//#define QSORT_CAST(x)      ((QSORT_DATA_PTR_TO*)(x))
//#define QSORT_COMPARE_PTR_TO(type) type

//#ifdef _SYS_ULTRIX     /* way old and quite kludgey C specs, BSD 4.2 da! */
//#define CALLOC_PTR_TO      char
//#define CALLOC_NUM_TYPE	   size_t /* unsigned in man page is wrong? */
//#define SIZEOF_TYPE	   size_t
//#define QSORT_DATA_PTR_TO  char
//#define QSORT_LENGTH       int   /* actually width is an int, not a unsigned */
//#endif
//
//#ifdef _SYS_AUX     /* just like ULTRIX? */
//#define CALLOC_PTR_TO      char
//#define CALLOC_NUM_TYPE	   size_t /* unsigned in man page is wrong? */
//#define SIZEOF_TYPE	   size_t
//#define QSORT_DATA_PTR_TO  char
//#define QSORT_LENGTH       int   /* actually width is an int, not a unsigned */
//#endif

//#ifdef _SYS_SUNOS     /* semi-ANSI is neither better nor worse */
//#define CALLOC_PTR_TO      void   /* the man pages lie? */
//#define CALLOC_NUM_TYPE	   size_t
//#define SIZEOF_TYPE	   size_t
//#define QSORT_DATA_PTR_TO  char
//#define QSORT_LENGTH       int   /* actually width is an int, not a unsigned */
//int system(); /* why is this in no include file */
//int fseek();  /* ditto */
///* Fucking Suns: exit() is really of type int, but it's void in the lint-lib */
//#endif

//#ifdef _SYS_WATCOM /* should be good for any ANSI system */
//#define CALLOC_PTR_TO      void
//#define CALLOC_NUM_TYPE	   size_t
//#define SIZEOF_TYPE	   size_t
//#define QSORT_DATA_PTR_TO  void
//#define QSORT_LENGTH       size_t
//#endif



/**************************** Terminal I/O ********************************
If TRY_GETENV_TERM is set, the system first tries to recognize the
string returned by getenv("TERM"). The system can deal with the
following terminal types, in terms of very, very basic screen I/O:

						Default  Escape	 Assume
 TERM variable	 Terminal			#Lines	 Codes   Scrollback?
 ------------	 -------			------	 -----	 -----------
 ansi*		 simple DEC VT-100 compatible	24	 ANSI	 Yes (new)
 dec*, vt*		   "			"        "       "
 xterm*	 	 X-Windows vt-100 emulator	"	 "	 "
 mac*	 	 A/UX CommandShell window	"	 "	 "
 sun*		 Sun workstation console or	36	 ANSI	 Yes (new)
 		 SunView shelltool/cmdtool
 pc*		 PC with ANSI.SYS loaded	25	 ANSI	 No
 hp*		 most Hewlett-Packards		24	 HP	 Yes
 300* 		 HP 9000 console		46	 HP	 Yes
 anything else	 unknown terminal type		24	 none	 No

[Note: here * means any characters can follow, including possibly none.]

If the TERM shell variable is set, then on Unix we may try to lookup
the screen length using the termcap functions, if HAVE_TERMCAP is
defined.  On any system we then try to use the LINES environment
variable (LINES over-rides termcap). We do not try to look up anything
else using termcap (too ugly, bleach). If HAVE_WINSIZE is defined, the
termios (POSIX) ioctl(fileno(stdout),TIOCGWINSZ,struct winsize*) is
used to get the #lines, and it WILL be used appropriately to see if
the window's size changed (we do not have a SIGWINCH handler however).
This value will supercede the defaults above. As yet no attempt is
made to determine the screen width (see note below).

There was a TRY_TERMINFO flag, to use terminfo instead of termcap, but
it may be broken. All Unixes seem to have termcap these days. Someone
wanna fix it?

IMPORTANT NOTE FOR WINDOWING SYSTEMS: If HAVE_WINSIZE does not work,
then it will be impossible for the shell to tell if the #lines changes
after it has started (e.g. if a terminal emulator window is resized).

In ALL cases, we really want to use an 80-character wide window or
screen! We assume 79 chars max per output line with tabs every 8
spaces, to maximize portability. (Really avoid system dependent
terminal widths, or otherwise transcript files must be printed on a
wider printer, other system's editors and so forth must be able to
deal, E-mailing files can be a problem, and... ick!). It's probably a
good idea to always keep TRANSLATE_TABS set.

If DEFAULT_SCROLLBACK is defined (to either TRUE or FALSE) the default
assumptions about scrollback are ignored. If scrollback is assumed, then
the screen is not cleared. More-mode is now always diabled at first.
***************************************************************************/
#define TRY_GETENV_TERM   /* Might as well always leave these defined? */
#define TRY_GETENV_LINES  /* see syscode.c */
#define TRY_ISATTY        
/* #define DEFAULT_SCROLLBACK TRUE */ /* define ONLY to over-ride */ 
#define DEFAULT_MORE_MODE TRUE /* is over-ridden with -nomore */
#define TRANSLATE_TABS TRUE
#define MEMORY_LINES   150
#define MAX_HOLD_LINES 150

//#ifdef _SYS_UNIX
#define DEFAULT_TERM_TYPE SCROLLING_ANSI 
#define TRY_TERMCAP
#define TRY_WINSIZE
//
//#else /* _SYS_DOS */
//#define DEFAULT_TERM_TYPE PC_CONSOLE
//#endif




/************************** Stdio interface *******************************
The TTY/Text File I/O library now assumes that the following functions
exist and work as in ANSI (or Unix System V, if not in K&R). Note that
we use a very small subset of the USUAL stdio/Unix-I/O library, for
portability's sake. If a function is not listed here, then we probably
found a portability problem somewhere, so don't use it! (Actually, as
always, the application code should only use the functions described
in iolib.h).

K&R did not specify what the output routine putc() does on errors,
nor does it specify its return value type. The same is true of fclose().
In practice, on different systems they do different things. Define 
xputc() and xclose() here to return TRUE if OK or FALSE if error.

We assume that fgets() provides standard I/O preprocessing (delete key 
and so forth). For details of the curses interface, see cursesio.c.

   printf(...)	  	     not used in general, but only in special cases
   fprintf(...)	  	     ditto
   scanf(...)	  	     ditto
   fscanf(...)	  	     ditto
   sprintf(...)		     used everywhere (aliased to "sf" in strlib.h) 
   int sscanf(...)    	     returns # tokens parsed
   fflush(FILE *fp)	     can't fail if fp is OK

   int feof(FILE *fp)	     returns T/F AFTER a call to xgetc() gives EOF
   int ferror(FILE *fp)      returns errno>=0 or 0 indicating "no error"
   int getc(FILE *fp)        returns c or EOF 
   ungetc(int c; FILE *fp)   it can't fail

   FILE *fopen(char *nam,*mode)         NULL if fail
   fseek(FILE *fp; long num; int dir)   it can't fail if fp is OK
   int fgets(FILE *fp; char *s; int n)  may return EOF char

Note that functions listed without return types here may or may not be of 
type void: on many systems they are of type int! However, the fuction's
return value is always (and should always be) discarded, if it exists.

We provide a workaround for fseek(), which seems to have bugs on VMS.
***************************************************************************/
#define xputc(c,fp) (putc(c,fp)!= EOF)
#define xclose(fp)  (fclose(fp)!= EOF)
#define file_seek(fp,char_num) fseek(fp,((long)(char_num)),0)
/* #define REPLACE_FSEEK *//* Needed on VMS? */


/**************** GLOBAL DECLARATIONS FOR USERS OF LIB *****************
Good working assumptions about data types, etc:

  * Doubles are big (say e+/-30, at least): avoid single precision (float)
    or other floating point types like the plauge.
  * A short is 2 bytes and a long is 4 bytes, both are signed.
  * An int is either a short or a long. On everything we use now, its long.
  * A char can take 7 (NOT 8) bits (0-127), positive values only.
  * unsigned ints exist, but there are a number of reasons not to use them 
    (it is easy for arguments to get passed incorrectly).
  * Never-ever assume that shifts wrap around correctly (or don't wrap 
    correctly), or that overflow/underflow ever produces predictable results
  * God knows what a pointer is. Always cast pointer types correctly,
    and avoid recasting pointers in general. 
  * Always make sure that arguments passed=arguments declared EXACTLY! Many
    bugs arrise if you don't. Never assume that parameters in a function call 
    will be automaticly be cast correctly. Be explicit.
  * Pass structs/unions/arrays as args by pointer ONLY. Passing by value is not
    even close to portable. Variable length argument lists are not portable.
    Returning structs/unions/arrays by value is also unportable.
  * Never-ever-ever try to side-effect constants. The biggest culprit here
    is passing a string argument to a function which then try to change 
    characters in the string. This causes bus errors on most systems.
  * void exists, but void* doesn't, nor does enum. Use void, not int, for
    procedures that return nothing
  * Really avoid other implicit int declarations (ex arguments w/o a type 
    are assumed to be ints - YUK!)
  * Only assume that the result of a boolean expression (ex x<y) is zero
    (FALSE) or nonzero (TRUE). Thus, expressions like "if (flag==TRUE) ..."
    are a very bad idea. Instead say "if (flag) ...".
  * Otherwise, what the original K&R version 1 says goes.

Note that the following definitions are specs and should never be
changed! They exist only to clarify code, not to allow any sort of
data-type abstraction!
******************************************************************************/

//#ifndef NO_INCLUDE_HELPERS

/* One of the lib_init() routines must be called before ANYTHING else
is done, or the program will crash in wierd ways! In general, user's
_init routines should come right after the lib_init() call and should
just initialize variables, malloc things, etc. They should not do
anything important, and particularly should not generate any I/O. */

//void lib_init(); /* no args - assumes line type tty I/O */
//
//void lib_inits(); /* args: int *argc_pointer; char *argv[]; both are
//   side-effected. This is more or less the same as lib_init() except that it
//   checks the environment and may start up either curses or wimp I/O via
//   screen_init(); */
//
//void custom_lib_init(); /* no args: This is like lib_init() except
//   that tty_init() is not called. A screen_init() routine must be invoked
//   before ANY I/O is attempted, or else... (see iolib.h) */
//
//void get_cmd_line_args();  /* args: int *argc_ptr; char **argv;
//   side-effects file_arg[] and nulls the args it parses, so another
//   arg sucker can go at it afterwards */

//#endif


/*************************** Library Include files **************************/
//#ifdef  INC_LIB
//#define INC_MATH
//#define INC_MSG
//#define INC_IO
//#define INC_STR
//#define INC_MEM
//#endif

//#ifndef _SYS_WATCOM
//#endif

//#ifdef INC_IO
//#endif

//#ifdef INC_MATH
//#ifdef _SYS_WATCOM
//#undef real
//#endif
//#ifdef _SYS_WATCOM
//#define real double
//#endif
//#endif

//#ifdef INC_MSG
//#endif

//#ifdef INC_MEM
//#endif

//#ifdef INC_STR
//#endif

/* The following are NOT included by INC_LIB */

//#ifdef INC_MISC
//#endif

//#ifdef INC_SHELL
//#endif

//#ifdef INC_TABLE
//#endif

//#ifdef INC_EQN
//#endif

//#ifdef INC_HISTO
//void make_histo(); /* that's all? */
//#endif

//#ifdef INC_STATS
//#endif

/********************* Defs internal to the library *********************
Here are definitions needed ONLY to compile the helpers files themselves.
None of this stuff should be used by the user code: much better (that is,
robust and portable) interfaces are provided by the library. */
/*************************************************************************/

//#ifdef INC_HELP_DEFS

//#ifdef TRY_WINSIZE
//#endif

///* The BSD random number functions... seemingly not declared anywhere */
//#ifdef USE_RANDOM
//long random();
////int srandom();
//#endif

///* The HPUX (and System V?) random number functions... ditto */
//#ifdef USE_DRAND48
//double drand48();
//void srand48();
//#endif

/* Library declarations only to be used by the helpers library code itself */
extern char *ps_, *ln_; /* defined in iolib.c */
//void dummy_math_calls();

//#endif /* for #ifdef INC_HELP_DEFS */

#undef fflush  /* a special version is used by the application code */

/*************************************************************************/
