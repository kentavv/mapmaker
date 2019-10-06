/******************************************************************************

  ####    #   #   ####    ####    ####   #####   ######           ####
 #         # #   #       #    #  #    #  #    #  #               #    #
  ####      #     ####   #       #    #  #    #  #####           #
      #     #         #  #       #    #  #    #  #        ###    #
 #    #     #    #    #  #    #  #    #  #    #  #        ###    #    #
  ####      #     ####    ####    ####   #####   ######   ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

/**************** SYSCODE.C - SYSTEM SPECIFIC CODE ***************************/
/* Most of this file is doccumented in system.h - it contains most of the 
   system dependent code for the helpers library. */

#define INC_LIB
#define INC_EQN
#define INC_HELP_DEFS
#include "system.h"
#include "syscode.h"

#ifdef _GNU_READLINE
#include "readline/readline.h"
#include "readline/history.h"
#endif

/*********************** C-Library Extensions ********************************/

/***** Time functions *****/
/* Note that time() and ctime() seem to be the only portable time functions.
   However, time() returns different types with different C compilers! */

static TIME_TYPE old_stamp, new_stamp;   /* For local use only! */

real usertime(bool do_reset) /* return time in seconds, or -1.0 if fail */
{
    real rtime;
    new_stamp= time((TIME_TYPE *)NULL); rtime= (real)(new_stamp - old_stamp);
    if (do_reset) old_stamp= new_stamp;
    return(rtime);
}

char *time_string(void)    /* return ptr to "" if fail */
{ 
    TIME_TYPE the_time;  /* note that asctime() does not always exist */
    char *str;
    int end;

    the_time=time((TIME_TYPE *)NULL); 
    str=ctime(&the_time); 
    if (str==NULL) return(ptr_to(""));
    end=len(str)-1;
    if (str[end]=='\n') str[end]='\0';
    return(str);
}


/***** subprocess functions *****/
 

bool shell_command(char *cmd)
{ 
    bool success;
    success=FALSE;

#ifdef NO_SYSTEM_FUNC
    return(FALSE);
#else
/* NEED WIMP HOOK */
#ifdef HAVE_CURSES
    bool had_curses= curses;
    if (curses) curses_end();  /* DO SOMETHING */
#endif
#ifndef VMS
    if (system(cmd)==0) success=TRUE;
#endif
#ifdef VMS
    if (system(cmd)!=0) success=TRUE;
#endif
#ifdef HAVE_CURSES
    if (success && had_curses) curses_refresh();
#endif
    return(success);
#endif
}


bool subshell(void)
{
        char *shell_name, cmd[120];
	bool success=FALSE;
	
#ifdef NO_SYSTEM_FUNC
    return(FALSE);
#else
/* NEED WIMP HOOK */
#ifdef HAVE_CURSES
	bool had_curses= curses;
	if (curses) curses_end();  /* DO SOMETHING */
#endif
#ifdef TRY_SHELL_VAR
	if ((shell_name=getenv("SHELL"))!=NULL && !nullstr(shell_name)) {
	    nstrcpy(cmd,shell_name,100); /* maxstrcat(cmd," -i",110); why? */
	    if (system(cmd)==0) success=TRUE;
	}
#endif
#ifdef TRY_COMSPEC_VAR
	if (!success && (shell_name=getenv("COMSPEC"))!=NULL && 
	    !nullstr(shell_name)) {
	    if (system(shell_name)==0) success=TRUE;
	}
#endif
#ifdef TRY_SHELL_CMD
	if (!success && !nullstr(TRY_SHELL_CMD)) {
	    if (system(TRY_SHELL_CMD)==0) success=TRUE;
	}
#endif
#ifdef HAVE_CURSES
	if (success && had_curses) curses_refresh();
#endif
	return(success);
#endif
}

		    
/***** get/set directories *****/

bool change_directory(char *dir)
{
    if (dir==NULL) send(CRASH);
#ifdef HAVE_CHDIR
    if (chdir(dir)==0) return(TRUE); 
#endif
    return(FALSE); 
}

bool get_directory(char *buf)
{ 
    if (buf==NULL) send(CRASH);
#ifdef HAVE_GETCWD
    if (getcwd(buf,PATH_LENGTH-2)!=NULL) return(TRUE); 
#endif
    return(FALSE); 
}

bool get_home_directory(char *buf)
{ 
    char *dir;
    if (buf==NULL) send(CRASH);
#ifdef HAVE_GETENV
    if ((dir=getenv("HOME"))!=NULL) 
      { nstrcpy(buf,dir,PATH_LENGTH); return(TRUE); }
#endif
    return(FALSE);
}

bool get_code_directory(char *buf)
{ 
    char *dir;
    if (buf==NULL) send(CRASH);
#ifdef HAVE_GETENV
    if ((dir=getenv("MAPM_LIB"))!=NULL)
      { nstrcpy(buf,dir,PATH_LENGTH); return(TRUE); }
#endif
#ifdef _CODE_DIR /* compiled in default */
    if (!nullstr(_CODE_DIR)) 
      { nstrcpy(buf,_CODE_DIR,PATH_LENGTH); return(TRUE); }
#endif
    return(FALSE);
}

bool rename_file(char *original_name, char *new_name)
{
#ifdef _SYS_DOS
    sprintf(ps,"copy %s %s",original_name,new_name);
    if (system(ps)==0) return(TRUE);
    else return(FALSE);
#else
    if (rename(original_name,new_name)==-1) return(FALSE);
    else return(TRUE);
#endif
}


bool fgoto_line(FILE *fp, long index)
{
#ifdef REPLACE_FSEEK
    long fseekvalue= 0L;
    frewind(fp);
    run while (fseekvalue < index-1) {
	fgetln(help_file);
	fseekvalue+=len(ln)+1;
    } except_when(ENDOFILE) return(FALSE);
    return(TRUE);
#endif
    return(fseek(fp,index,0)!=-1);
}

/***** random number functions *****/

long mkseed(long x)
{ if (x==RANDOM) return((long)time(NULL)); else return(x); }

//#ifdef USE_RANDOM
void do_seedrand(long x) { srandom((int)mkseed(x)); }
real randnum() { return(((real)random())/2147483648.0); }
//#else
//#ifdef USE_DRAND48
//void do_seedrand(long x) long x; { srand48(mkseed(x)); }
//real randnum(void) { return(drand48()); }
//#else /* USE_SRAND */
//void do_seedrand(x) long x; { srand((int)mkseed(x)); }
//real randnum() { return(((real)rand())/((real)(RAND_MAX+1))); }
//#endif
//#endif


/***** message and signal handling *****/

void untrapped_msg(void) /* DO NOT ASSUME THAT MSGNAME IS SET! */
{
    /* if (msg!=IOERROR) flush(); most are disk errors */
    if (msg<1 || msg>MSGS) 
      { fprintf(stderr, "Untrapped error %d (?)\n",msg); exit(1); }
    fprintf(stderr,"Untrapped error %d (%s)\n",msg,mname[msg]);
    (*(mstrmsg[msg]))(ps_); fprintf(stderr,ps_); fprintf(stderr,"\n");
}

void trapped_msg(void) /* DO NOT ASSUME THAT MSGNAME IS SET! */
{
    /* if (msg!=IOERROR) flush(); most are disk errors */
    if (msg<1 || msg>MSGS) 
      { fprintf(stderr,"Error %d (?)\n",msg); exit(1); }
    fprintf(stderr,"Error %d (%s)\n",msg,mname[msg]);
    (*(mstrmsg[msg]))(ps_); fprintf(stderr,ps_); fprintf(stderr,"\n");
}

#define SHUTDOWN1 "NOTE: In some extreme cases you may have to shut down and "
#define SHUTDOWN2 "restart the\nprogram in order to resume proper operation. "
#define SHUTDOWN3 "Hit <return> to continue..."

void verbose_untrapped_msg(void) /* DO NOT ASSUME THAT MSGNAME IS SET! */
{
  fprintf(stderr,"*** Drats! An unhandled internal error occured. ***\n");
  if (msg<1 || msg>MSGS) { fprintf(stderr,"error #%d (?)\n",msg); exit(1); }
    else fprintf(stderr,"Error message #%d (%s) was sent.\n",msg,mname[msg]);
  (*(mstrmsg[msg]))(ps_); fprintf(stderr,ps_); 
  if (!nullstr(ps_)) fprintf(stderr,"\n");
  fprintf(stderr,SHUTDOWN1); fprintf(stderr,SHUTDOWN2); 
//  fprintf(stderr,SHUTDOWN3); fgets(ps_,MAXLINE,stdin);
  fprintf(stderr,"\n");
}

void do_trap(void)
{ 
  if (msg<1 || msg>MSGS) 
      { fprintf(stderr,"Illegal trap message\n"); exit(1); }
  fprintf(stderr,"Trapped message %d (%s)\n",msg,mname[msg]);
  (*(mstrmsg[msg]))(ps_); fprintf(stderr,ps_); 
}


static int signals;

void sigcounter(void)
{ signals++; if (signals>MAX_BAD_SIGNALS) send(CRASH); }


void signal_trap_init(void)
{  
  signals=0;

  signal(SIGQUIT, handle_quit);		/* ANSI Signals - for Microsoft C */
  signal(SIGINT,  handle_interrupt);
  signal(SIGFPE,  handle_matherror);
  signal(SIGILL,  handle_weird_signal);
  signal(SIGSEGV, handle_buserror);
  signal(SIGTERM, handle_quit);
}



/********************************** I/O *************************************/

int old_term, old_lines, old_scrollback, old_more, dos_output;
bool tried_curses, tried_wimp, tried_split, have_drawn_top;
int tty_errors, file_errors, puts_errors;
int curses, split, wimp, use_gnu_readline; /* externed global bools */
char **file_arg;
int prep_it, append_it;


bool do_gnu_readline(char *prompt, char *str, int num)
{
#ifndef _GNU_READLINE
    send(CRASH);
    return(FALSE);
#else
    char *result=NULL;

    result= readline(prompt);
    if (result==NULL) return(FALSE); /* EOF */

    nstrcpy(str,result,num-2);
    add_history(str);
    free((char*)result);
    return(TRUE);
#endif
}


#ifdef _GNU_READLINE
static char *default_text = NULL;

int set_rl_default(void)
{
    if(default_text)  {
        int n = strlen(default_text);
        rl_extend_line_buffer(n + 1);
        strcpy(rl_line_buffer, default_text);
        rl_point = rl_end = n;
    }
    return 0;
}
#endif

	
bool do_gnu_edit(char *prompt, char *str, int num, char *initial /* initial may be = str */)
{
#ifndef _GNU_READLINE
    send(CRASH);
    return(FALSE);
#else
    char *result=NULL, *hist_entry=NULL;
    Function *old_handler;

    old_handler = rl_startup_hook;
    rl_startup_hook = set_rl_default;
    default_text = initial;
    result = readline(prompt);
    rl_startup_hook = old_handler;

    if (result==NULL) return(FALSE); /* EOF */

    nstrcpy(str,result,num-2);
    add_history(str);
    free((char*)result);
    return(TRUE);
#endif
}

	
bool gnu_copyright(char *str /* side-effected, so it must be big enough */)
{
#ifndef _GNU_READLINE
    return(FALSE);
#else
    if (!use_gnu_readline) return(FALSE);
    sprintf(str,"GNU Readline Copyright 1988-1989, Free Software Foundation");
    return(TRUE);
#endif
}


bool tty_gets(char *str, int num)
{
    int n, i;

    in_tty_gets = TRUE; hit_interrupt = FALSE;

    if (num<2) send(CRASH);
    for (i=0; i<num+2; i++) str[i]='\0';      /* To get ok string if EOF */
    
    if (fgets(str,num+1,in)==NULL) {         /* EOF or error */
	str[num]='\0';                       /* in case fgets() is weird */
	if (++tty_errors>MAX_IO_FAILURES) send(CRASH);
	if ((n=ferror(in))!=0) ioerror(strerror(n),in,str); /* error */
	if (feof(in)) return(FALSE);         /* else must be EOF */
	ioerror("fgets() failed",in,str);    /* I dunno? */
    }

    if(hit_interrupt) { str[0] = '\n'; str[1] = '\0'; }
    in_tty_gets = FALSE;

    cursor=0; printed_lines=0;
    if (str[0]=='\0' || str[len(str)-1]!='\n') { /* no \n => truncated */
	if (++tty_errors>MAX_IO_FAILURES) send(CRASH);
	else ioerror("input line too long",in,str);
    }
    tty_errors= 0;
    return(TRUE);
}


bool file_gets(FILE *fp /* must be opened with file_open() */, char *str  /* must be num+2 chars long, but use num+3 in case of weirdness */, int num    /* num chars, not including the '\n' or '\0', will be read */)
{
    int i, c, n;
    
    for (i=0; (c=fgetc(fp))!='\n'; i++)
      if (c==EOF) { /* error or EOF */
	  str[i]='\0'; 
	  if (++file_errors>MAX_IO_FAILURES) send(CRASH);
	  else if (feof(fp)) { if (i==0) return(FALSE); else return(TRUE); }
	  else if ((n=ferror(fp))!=0) { ioerror(strerror(n),fp,str); }
	  else ioerror("fgetc() failed",fp,str);
      } else if (i==num-1) { 
	  str[i]='\0'; 
	  while ((c=fgetc(fp))!='\n') if (c==EOF) break;
	  if (++file_errors>MAX_IO_FAILURES) send(CRASH);
	  ioerror("input line too long",fp,str);
      } else { /* all is OK */
	  str[i]=c; continue;
      }
    str[i]='\0'; 
    file_errors= 0;
    return(TRUE);
}


void lib_puts(FILE *fp, char *str)
{
    int i, n;
    char c;
    
    if (fp==NULL || str==NULL) send(CRASH);
#ifdef HAVE_CURSES
    if (curses && fp==out) { curses_puts(str); } else
#endif
    {
	for (i=0; (c=str[i])!='\0'; i++) {
	    if (c=='\n' && dos_output && fp!=out && !xputc('\015',fp)) {
		if (++puts_errors>=MAX_IO_FAILURES) send(CRASH);
		if ((n=ferror(fp))!=0) ioerror(strerror(n),fp,str);
		else ioerror("xputc() failed",fp,str);
	    }
	    if (!xputc(c,fp)) {
		if (++puts_errors>=MAX_IO_FAILURES) send(CRASH);
		if ((n=ferror(fp))!=0) ioerror(strerror(n),fp,str);
		else ioerror("xputc() failed",fp,str);
	    }
	}
	fflush(out); 
    }
    puts_errors= 0;
}


void iocheck(void) { return; }


void tty_init(void)
{ 
    char *tty_type, *num_lines, copy[10], bp[1025];
    int x;

/* THIS HAS NO WIMP HOOK - IT SHOULD NEVER BE CALLED IF WIMP I/O IS USED! */

#ifdef TRY_ISATTY
    if (!isatty(fileno(stdin)))  { interactive=FALSE; ignore_eof=FALSE; }
    if (!isatty(fileno(stdout))) { screen=FALSE; }
#endif
    if (!screen) return; /* term=TERM_UNKNOWN, more_mode=scrollback=FALSE */

/* If we DO assume that a tty type has scrollback, we should NOT clear
its screen, NOR should more_mode be on by default. If it does NOT have
scrollback, we MAY clear its screen and might turn on more_mode (more
mode will be way ugly w/o cursor motion however). In some sense,
scrollback=TRUE is the conservative option.*/

#ifdef TRY_GETENV_TERM
    if ((tty_type=getenv("TERM"))!=NULL) {
	nstrcpy(copy,tty_type,9); crunch(copy);
	if (nstreq(copy,"hp",2)) 
	  { term=HP_TERM; scrollback=TRUE; tty_lines=24; }
	else if (nstreq(copy,"300h",4)) /* hp series 300 console */
	  { term=HP_TERM; scrollback= TRUE; tty_lines=46; }
	else if (nstreq(copy,"ansi",4) || nstreq(copy,"vt",2) ||
		 nstreq(copy,"dec",3)  || nstreq(copy,"mac",3))
	  { term=ANSI; scrollback=TRUE; tty_lines=24; }
	else if (nstreq(copy,"xterm",5)) /* X-Windows terminal emulator */
	  { term=ANSI; scrollback=TRUE; tty_lines=24; }
	else if (nstreq(copy,"sun",3)) /* sun console/cmdtool/shelltool */
	  { term=ANSI; scrollback=TRUE; tty_lines=34; }
	else if (nstreq(copy,"pc",2)) /* pc console */
	  { term=ANSI; scrollback=FALSE; tty_lines=25; }
#ifdef _GNU_READLINE
    /* the cmd-line override '-simple' sets this back to FALSE */
    if (interactive) use_gnu_readline=TRUE; /* else was set to FALSE */
#endif
#ifdef TRY_TERMCAP
    if (!nullstr(tty_type)) {
	tgetent(bp,tty_type);
	if ((x=tgetnum("li"))>0) tty_lines=x;
    }
#endif
    } else /* can't getenv("TERM") */
#endif /* TRY_GETENV_TERM */

    if (term==TERM_UNKNOWN) {
	/* either no TERM environment variable or its value is unknown */
	use_gnu_readline=FALSE; /* it will fail anyway w/o TERM var */
	if (DEFAULT_TERM_TYPE==HP_TERM) /* all HP_TERMs scrollback? */
	  { scrollback=TRUE; term=HP_TERM; tty_lines=24; }
	else if (DEFAULT_TERM_TYPE==SCROLLING_ANSI)
	  { scrollback=TRUE; term=ANSI; tty_lines=24; }
	else if (DEFAULT_TERM_TYPE==NONSCROLLING_ANSI)
	  { scrollback=FALSE; term=ANSI; tty_lines=24; }
	else if (DEFAULT_TERM_TYPE==PC_CONSOLE)
	  { scrollback=FALSE; term=ANSI; tty_lines=25; }
	else if (DEFAULT_TERM_TYPE==MAC_WINDOW)
	  { scrollback=FALSE; term=MAC_WINDOW;  tty_lines=24; }
	else /* default is TERM_UNKNOWN terminal type */
	  { scrollback=TRUE; term=TERM_UNKNOWN; tty_lines=24; }
    }

#ifdef TRY_GETENV_LINES
    if ((num_lines=getenv("LINES"))!=NULL && sscanf(num_lines,"%d",&x)==1) 
      tty_lines=x;
#endif
    check_tty_lines(); /* ioctl will always over-ride */
#ifdef DEFAULT_SCROLLBACK
    scrollback=DEFAULT_SCROLLBACK;  /* override the settings above */ 
#endif
    /* if (term!=TERM_UNKNOWN) screen=TRUE; else screen=FALSE; */    
}


bool check_tty_lines(void) /* return TRUE and set tty_lines if changed */
{
/* maybe add some weird PC thing here to get #lines */
#ifdef TRY_WINSIZE
    struct winsize thesize;
    if (screen && ioctl(fileno(stdout),TIOCGWINSZ,&thesize)==0)
      tty_lines= thesize.ws_row;
#endif
    return FALSE;
}


/* Use lib_puts(out,...) (not print()) for these screen handling routines! 
   flush() will be executed immediately beforehand.... */

int save_cursor;
char Tcmd[100];

/* These have been tested on a Xterm and vt220 */
#define ansi_tty_init()      lib_puts(out,"\033[0m\n")
#define ansi_clr_scrn()      lib_puts(out,"\033[1;1H\033[2J")
#define ansi_highlight(on)   lib_puts(out,on ? "\033[7m":"\033[0m")
#define ansi_del_prev_ln()   lib_puts(out,"\033[99D\033[K\033[1A\033[K")
#define ansi_boing()         lib_puts(out,"\007")
void ansi_cursor_left(int i, char *s)
{ if(i<0) sprintf(Tcmd, "\033[99D\033[K%s", s); else sprintf(Tcmd, "\033[%dD\033[K%s", i, s);
  lib_puts(out,Tcmd); }

#define hp_tty_init()        lib_puts(out,"\n\033&d@\n")
#define hp_clr_scrn()        lib_puts(out,"\033H\033J")
#define hp_highlight(on)     lib_puts(out,on ? "\033&dB":"\033&d@")
#define hp_del_prev_ln()     lib_puts(out,"\033&a0C\033K\033A\033K")
#define hp_boing()           lib_puts(out,"\007")
void hp_cursor_left(int i, char *s)
{ if(i<0) sprintf(Tcmd, "\033&a0C\033K%s", s); else sprintf(Tcmd, "\033&a-%dC\033K%s", i, s);
  lib_puts(out,Tcmd); }

/* These should be filled in for the Mac's ThinkC "console package" */
#define mac_tty_init()   	printf("\n") /* printf() opens window */
#define mac_clr_scrn()      	{}
#define mac_highlight(on)   	{}
#define mac_del_prev_ln()   	{}
#define mac_boing()       	{}
void mac_cursor_left(void)		{}

void tty_hello(void)
{
    if (term==HP_TERM) hp_tty_init();
    else if (term==ANSI) ansi_tty_init();
    else lib_puts(out,"\n");
    if (!scrollback) { do_clear_screen(); lib_puts(out,"\n"); }

#ifdef _GNU_READLINE
    if (use_gnu_readline) rl_bind_key('\t', rl_insert); /* completion off */
#endif
}


bool do_clear_screen(void)
{ 
/* NEEDS WIMP AND MAC CONSOLE HOOK */
    if (term==HP_TERM)     { hp_clr_scrn();   fflush(out); return(TRUE); }
    else if (term==ANSI)   { ansi_clr_scrn(); fflush(out); return(TRUE); }
#ifdef HAVE_CURSES
    else if (term==CURSES) { curses_clr_scrn(); return(TRUE); }
#endif    
    else return(FALSE);
}


bool do_delete_previous_line(void) /* Needed for the "Hit RETURN for more" thing */
{
/* NEEDS WIMP AND MAC CONSOLE HOOK */
    if (term==HP_TERM)     { hp_del_prev_ln();   fflush(out); return(TRUE); }
    else if (term==ANSI)   { ansi_del_prev_ln(); fflush(out); return(TRUE); }
#ifdef HAVE_CURSES
    else if (term==CURSES) { curses_del_prev_ln(); return(TRUE); }
#endif
    else return(FALSE);
}


bool do_highlight(bool reverse)
{
/* NEEDS WIMP AND MAC CONSOLE HOOK */
#ifdef HAVE_CURSES
    if (term==CURSES)    { curses_set_highlight(reverse); return(TRUE); }
#endif
    if (term==HP_TERM)   { hp_highlight(reverse);  fflush(out); return(TRUE); }
    else if (term==ANSI) { ansi_highlight(reverse);fflush(out); return(TRUE); }
    else return(FALSE);
}


bool do_cursor_left(int num_spaces /* might be FAR_LEFT */, char *str_to_print)
{
/* NEEDS WIMP AND MAC CONSOLE HOOK */
    if (term==HP_TERM)     { hp_cursor_left(num_spaces,str_to_print); }
    else if (term==ANSI)   { ansi_cursor_left(num_spaces,str_to_print); }
#ifdef HAVE_CURSES
    else if (term==CURSES) { curses_cursor_left(num_spaces,str_to_print); }
#endif
    else return(FALSE);
    return(TRUE);
}


bool boing(void)
{
/* NEEDS WIMP AND MAC CONSOLE HOOK */
    if (term==HP_TERM)     { hp_boing(); return(TRUE); } 
    else if (term==ANSI)   { ansi_boing(); return(TRUE); } 
#ifdef HAVE_CURSES
    else if (term==CURSES) { curses_boing(); return(TRUE); } 
#endif
    else return(FALSE);
}


/****************************** TOPLEVEL STUFF ******************************/

void misc_init(void) /* init this file */
{ 
    int i;
    matrix(file_arg,MAX_FILE_ARGS,PATH_LENGTH+1,char);
    for (i=0; i<MAX_FILE_ARGS; i++) file_arg[i][0]='\0';
    old_stamp=new_stamp=time(NULL); dos_output=FALSE; 
}


void custom_lib_init(void)
{
/* These init routines shouldn't really DO much - they really are meant for
   initializing variables, mallocing structures, and so forth. Serious work
   (for ex kicking off fancy screen I/O) should happen elsewhere. */
       
    mem_init();  /* memlib.c */
    msg_init();  /* msglib.c */
    io_init();   /* iolib.c - just allocates things etc for the I/O lib */
    str_init();  /* strlib.c */
    math_init(); /* mathlib.c */
    eqn_init();  /* eqn.c */
    misc_init(); /* syscode.c (this file, dummy) */
    tty_init();  /* also in this file */
}

void lib_init(void)
{ 
    custom_lib_init();
    tty_hello();
}

void lib_inits(int *argc_ptr, char *argv[])
{
    custom_lib_init();
    /* if (!screen_init(argc_ptr,argv)) */ 
    get_cmd_line_args(argc_ptr,argv);
    tty_hello();
}


/* This is a little specialized for MAPMAKER right now. Generalize soon... */
#define ERROR_BADARG \
  "%s error: unrecognized option '%s'\ntype '%s -help' for help\n"
#define ERROR_ONEFILE   "%s error: only one %s argument allowed\n"
#define ERROR_BADFILE   "%s error: bad %s file name '%s'\n"
#define ERROR_NOFILE    "%s error: can't open %s file '%s'\n"
#define ERROR_EMPTYFILE "%s error: %s file '%s' is empty\n"
#define ERROR_NOTSCREEN "%s error: output is not a terminal, can't use %s\n"
#define HELP_ARG_MSG "\
%s optional arguments: [%csimple] [%cnomore] [%cclear] [%chelp] \n\
  [%cload file] [%cprep file] [%crun file] [%cphoto file] [%cout file]\n\
run %s and type 'help' for help with commands and other information \n"

#define NOTSCREEN_ \
{ fprintf(stderr,ERROR_NOTSCREEN,argv[0],argv[i]); abnormal_exit(); }

void get_cmd_line_args(int *argc_ptr, char *argv[])
{
    int i;

    for (i=1; i<*argc_ptr; i++) {
	if (argv[i][0]!=ARG_CHAR) {
	    fprintf(stderr,ERROR_BADARG,argv[0],argv[i],argv[0]); 
	    abnormal_exit();

	} else if (matches(argv[i]+1,"simple")) {
	    if (!screen) NOTSCREEN_
	      term=TERM_UNKNOWN; use_gnu_readline=FALSE;
	    argv[i][0]='\0';
	} else if (matches(argv[i]+1,"nomore")) {
	    if (!screen) NOTSCREEN_
	      more_mode=FALSE;
	    argv[i][0]='\0';
	} else if (matches(argv[i]+1,"clear")) {
	    if (!screen) NOTSCREEN_
	      if (term==PC_CONSOLE) scrollback=TRUE; else scrollback=FALSE;
	    argv[i][0]='\0';

	} else if (matches(argv[i]+1,"load")) {
	    check_file_arg(*argc_ptr-i,argv[i+1],file_arg[LOAD_FILE_ARG],
			   "load","data",argv[0],READ);
	    prep_it=FALSE;
	    argv[i++][0]='\0'; argv[i][0]='\0';
	} else if (matches(argv[i]+1,"prep")) {
	    check_file_arg(*argc_ptr-i,argv[i+1],file_arg[LOAD_FILE_ARG],
			   "prep","raw",argv[0],READ);
	    prep_it=TRUE;
	    argv[i++][0]='\0'; argv[i][0]='\0';
	} else if (matches(argv[i]+1,"run")) {
	    check_file_arg(*argc_ptr-i,argv[i+1],file_arg[RUN_FILE_ARG],
			   "run","in",argv[0],READ);
	    argv[i++][0]='\0'; argv[i][0]='\0';
	} else if (matches(argv[i]+1,"photo")) {
	    check_file_arg(*argc_ptr-i,argv[i+1],file_arg[PHOTO_FILE_ARG],
			   "photo","out",argv[0],APPEND);
	    append_it=TRUE;
	    argv[i++][0]='\0'; argv[i][0]='\0';
	} else if (matches(argv[i]+1,"out")) {
	    check_file_arg(*argc_ptr-i,argv[i+1],file_arg[PHOTO_FILE_ARG],
			   "photo","out",argv[0],WRITE);
	    append_it=FALSE;
	    argv[i++][0]='\0'; argv[i][0]='\0';
	} else if (matches(argv[i]+1,"help")) {
	    fprintf(stderr,HELP_ARG_MSG,argv[0],ARG_CHAR,ARG_CHAR,ARG_CHAR,
		    ARG_CHAR,ARG_CHAR,ARG_CHAR,ARG_CHAR,ARG_CHAR,ARG_CHAR,
		    argv[0]);
	    normal_exit();

	} else { 
	    fprintf(stderr,ERROR_BADARG,argv[0],argv[i],argv[0]); 
	    abnormal_exit();
	}
    }
}


bool check_file_arg(int num, char *arg, char *name, char *type, char *def_ext, char *prog, char *mode)
{
    char file[PATH_LENGTH+1];
    FILE *fp;

    if (!nullstr(name)) {
	fprintf(stderr,ERROR_ONEFILE,prog,type); 
	abnormal_exit();
    }
    nstrcpy(file,arg,PATH_LENGTH);
    if (!make_filename(file,DEFAULT_EXTENSION,def_ext)) {
	fprintf(stderr,ERROR_BADFILE,prog,type,file); 
	abnormal_exit();
    }
    run {
	fp=open_file(file,mode);
	close_file(fp);
	strcpy(name,file);
    } except {
	when IOERROR: 
	when CANTOPEN:
	when CANTCLOSE:
	  fprintf(stderr,ERROR_NOFILE,prog,type,file);
	  abnormal_exit();
	when ENDOFILE:
	  fprintf(stderr,ERROR_EMPTYFILE,prog,type,file);
	  abnormal_exit();
	default: 
	  relay; 
    }
    return FALSE;
}


bool update_top(void)
{ return(FALSE); }




//#ifdef PUNT_FOR_NOW /*******************************************************/
//
//bool screen_init(int *argc_ptr, char *argv[]) /* side-effect wimp, curses, and split */
//{
//    bool try_curses, try_wimp;
//    get_screen_preferences(argc_ptr,argv,&try_curses,&try_wimp);
//#ifdef HAVE_WIMP
//    if (!tried_wimp && try_wimp)
//      if (do_wimp_init(argc_ptr,argv)) term=WIMP; /* wimp=TRUE */
//    tried_wimp= TRUE;
//#endif
//    tty_init();
//#ifdef HAVE_CURSES
//    if(!wimp && !tried_curses && try_curses)
//	if (curses_init(&tty_lines)) term=CURSES; /* sets curses=TRUE */
//    tried_curses= TRUE;
//#endif
//    maybe_clear_screen(); nl(); boing();
//    return(wimp || curses);
//}
//
//
//bool split_screen_init(int *argc_ptr, char *argv[], int top_lines, void (*top_update_function)())
//{
//    bool try_curses, try_wimp;
//    get_screen_preferences(argc_ptr,argv,&try_curses,&try_wimp);
//#ifdef HAVE_WIMP
//    if (!tried_wimp && try_wimp)
//      if (do_split_wimp_init(top_lines,top_update_function,argc_ptr,argv))
//	term=WIMP; /* sets wimp=TRUE */
//    tried_wimp= TRUE;
//#endif
//    tty_init();
//#ifdef HAVE_CURSES
//    if (!wimp && !tried_curses && try_curses)
//      if (curses_init(&tty_lines)) {
//	  if (curses_split(top_lines,CURSES_REVERSE_TOP,top_update_function,
//			   &tty_lines)) term=CURSES; /* sets curses=TRUE */
//	  else curses_end();
//      }
//    tried_curses= TRUE;
//#endif
//    maybe_clear_screen(); update_top(); nl(); boing();
//    return(wimp || curses);
//}
//
//
//void get_screen_preferences(int *argc_ptr, char *argv[], int *try_curses, int *try_wimp)
//{
//    int i, j;
//
//    *try_curses= *try_wimp= MAYBE;
//    for (i=0; i<*argc_ptr; i++) {
//	if (nmatches(argv[i],"-window",2)) *try_wimp=TRUE;
//	else if (nmatches(argv[i],"+window",2)) *try_wimp=FALSE;
//	else if (nmatches(argv[i],"-nowindow",4)) *try_wimp=FALSE;
//	else if (nmatches(argv[i],"-screen",2)) *try_curses=TRUE;
//	else if (nmatches(argv[i],"+screen",2)) *try_curses=FALSE;
//	else if (nmatches(argv[i],"-noscreen",4)) *try_curses=FALSE;
//	else if (nmatches(argv[i],"-line",2)) *try_curses= *try_wimp= FALSE;
//	else continue;
//	/* matched argument, so delete it from the list */
//	for (j=i+1; j<*argc_ptr; j++) argv[j-1]= argv[j];
//	--*argc_ptr;
//    }
//    if (*try_curses==MAYBE) *try_curses=DEFAULT_TRY_CURSES;
//    if (*try_wimp==MAYBE) *try_wimp=DEFAULT_TRY_WIMP;
//}
//
//
///* Custom windows can work similarly to the canned types (text &
//split), except that the widow code is very application specific: it
//may respond to mouse clicks and so forth in special ways and have
//other special stuff, which will probably have to interect with the
//application code via global variables. Such interfaces could be
//implemented using a customized do_wimp_init() routine. To make life
//easy and so that things work similarly, we require (1) that the window
//have a scrolling text region, (2) that it uses the menu struct (see
//shell.c), (3) that it only need to know about state variable changes
//that already call update_top() (or maybe_ok()), and (4) that it
//effects state changes only when inhibit_menus is FALSE (see shell.c).
//Like all init routines, do_custom_wimp_init() should do very little -
//do_wimp_start() in shell.c should do this. */
//
//
//bool update_top(void)
//{
///* NEEDS WIMP HOOK */
//#ifdef HAVE_CURSES
//    if (curses && split) {
//	if (!have_drawn_top) curses_draw_top(); else curses_update_top();
//	have_drawn_top=TRUE; return(TRUE);
//    }
//#endif
//    return(FALSE);
//}
//
//
//void screen_end(void) {
///* NEEDS WIMP HOOK */
///* ALSO NEEDS a "Hit return" hook for quiting under a window system */
//#ifdef HAVE_CURSES
//    if (curses) curses_end(FALSE);
//#endif
//}
//
//
//bool string_editor(char *prompt, char *str, int num /* max #chars in str */)
//{
//    bool successful;
//
///* NEEDS WIMP HOOK */
//    successful = FALSE;
//#ifdef HAVE_CURSES
//    if(curses)
//      successful = curses_string_editor(str,num);
//#endif
//    return(successful);
//}
//
//
//#ifdef HAVE_CURSES
//void curses_error(char *where)
//{
//    curses_end(TRUE);
//    printed_lines=4; more_break_pending=FALSE;
//    term= old_term;  tty_lines=old_lines;
//    more= old_more;  scrollback= old_scrollback;
//    clear_screen();
//
//    if (where==NULL) where=ptr_to("?");
//    fprintf(stderr,"\nwarning: curses failed in %-50s\n",where);
//    fprintf(stderr,"attempting to continue in line mode...\n\n");
//    flush();
//}
//#endif
//
//#endif /* PUNT_FOR_NOW */
