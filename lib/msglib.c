/******************************************************************************

 #    #   ####    ####   #          #    #####            ####
 ##  ##  #       #    #  #          #    #    #          #    #
 # ## #   ####   #       #          #    #####           #
 #    #       #  #  ###  #          #    #    #   ###    #
 #    #  #    #  #    #  #          #    #    #   ###    #    #
 #    #   ####    ####   ######     #    #####    ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

/***** MESSAGE (EXCEPTION) HANDLING ROUTINES ******/

#define INC_IO
#define INC_MSG
#define INC_MEM
#define INC_STR
#define INC_MATH
#define INC_HELP_DEFS
#include "system.h"

int msg; 	/* user-accessible - valid for last message sent */
char *msgname;  /* see list of lib messages in helpers.h - HAS MOVED? */
char *msgstr;
bool exiting1, exiting0;
bool in_tty_gets, hit_interrupt;

/* message vars - valid for specific messages */
char *MATHERROR_type;
real  MATHERROR_arg1, MATHERROR_arg2;
int   NOMEMORY_num_cells, NOMEMORY_cell_size;
char *CANTOPEN_path, CANTOPEN_modechar;
char *IOERROR_errmsg, *IOERROR_filename, IOERROR_modechar, *IOERROR_linecopy;
char *SYSERROR_errmsg;
char *QUIT_errmsg;
char *BADEQN_errmsg;
int   BADEQN_errpos;


void strmsg_default(char *str)
 { str[0]='\0'; }
void strmsg_MATHERROR(str) char *str;
 { sprintf(str, "error type=%-60s\narg1=%lf\narg2=%lf",
           MATHERROR_type, MATHERROR_arg1, MATHERROR_arg2); }
void strmsg_NOMEMORY(str) char *str;
 { sprintf(str, "cell size=%d\nnum cells=%d", NOMEMORY_cell_size, NOMEMORY_num_cells);}
void strmsg_CANTOPEN(str) char *str;
 { sprintf(str, "filename=%-60s\nmode=%c", CANTOPEN_path, CANTOPEN_modechar); }
void strmsg_IOERROR(str) char *str;
 { sprintf(str, "filename=%-60s\nline=%-60s\nerror=%-60s",
           IOERROR_filename, truncstr(IOERROR_linecopy,60), IOERROR_errmsg); }
void strmsg_SYSERROR(str) char *str;
 { sprintf(str, "errmsg=%-60s", SYSERROR_errmsg); }

jmp_buf stk[TRAP_DEPTH];     /* stack of message traps */
int lvl;       	             /* trap stack pointer */
int last; 		     /* last msg num assigned */
char **mname;                /* [msg] => name string */
void (*(maction[MSGS]))(int);   /* [msg] => function which does whatever */
void (*(mstrmsg[MSGS]))(char *);   /* [msg] => function which sets its (char*) arg */
void sighandle(); 

void msg_init()
{
	msgname=NULL; mname= NULL; lvl=0; 
	matrix(mname,MSGS,TOKLEN+1,char);
	array(msgstr,MAXLINE+1,char);

	for(zz=0;zz<MSGS;zz++) {
		strcpy(mname[zz],"???");
		maction[zz]= default_action; 
		mstrmsg[zz]= strmsg_default;
	}

	array(MATHERROR_type,LINE,char);
	array(CANTOPEN_path,PATH_LENGTH,char);
	array(IOERROR_errmsg,LINE,char);
	array(IOERROR_filename,PATH_LENGTH,char);
	array(IOERROR_linecopy,LINE,char);
	array(SYSERROR_errmsg,LINE,char);
	array(QUIT_errmsg,LINE,char);
	
	setmsg(RERUN,"attempt to RERUN failed",sender,strmsg_default);
	setmsg(PUNT,"attempt to PUNT failed",sender,strmsg_default); 
	
	setmsg(CRASH,"internal software error",sender,strmsg_default);
	setmsg(QUIT,"QUIT signal",sender,strmsg_default);
	setmsg(NOMEMORY,"request for memory failed",sender,strmsg_NOMEMORY);
	setmsg(ENDOINPUT,"end of input",sender,strmsg_default);
	setmsg(INTERRUPT,"interrupt signal",sender,strmsg_default);
	setmsg(SOFTABORT,"software abort message",sender,strmsg_default);
	setmsg(MATHERROR,"math error",sender,strmsg_MATHERROR);
	setmsg(IOERROR,"I/O error",sender,strmsg_IOERROR);

	setmsg(ENDOFILE,"end of file",sender,strmsg_default);
	setmsg(CANTOPEN,"unable to open file - maybe no write permission",
	       sender,strmsg_CANTOPEN);	
	setmsg(SYSERROR,"segmentation violation or bus error",sender,
	      strmsg_SYSERROR);
	setmsg(BADEQN,"error in equation",sender,strmsg_default);
	setmsg(CANTCLOSE,"unable to close file - disk may be full",
	       sender,strmsg_default);

	last= USER_MESSAGE(0);
	signal_trap_init(); /* in syscode.c */
}


void setmsg(var,nam,action,disp) /* for a pre-assigned msg number */
int var;
char *nam;
void (*action)();
void (*disp)();
{ if (nam!=NULL) nstrcpy(mname[var],nam,MSGNAMLEN); 
  if (action!=NULL) maction[var]=action; 
  if (disp!=NULL) mstrmsg[var]=disp; }


int lvl_plus_plus() 
{ if (lvl<TRAP_DEPTH) return(lvl++);
  msg=CRASH; untrapped_msg(); abnormal_exit(); return(0); }
  


void sender(int num)
{
  if (num==0) return;

  if (num>0 && num<MSGS) msg=num; else msg=CRASH;
  msgname=mname[msg];
  if (--lvl>=0) { longjmp(stk[lvl],msg); return; }
  else { untrapped_msg(); abnormal_exit(); } /* see syscode.c */
}


void punter(/*int num*/)
{ if (--lvl>=0) longjmp(stk[lvl],QUIT); else {msg=CRASH; untrapped_msg();} }

void trapper(num) int num;
{ void do_trap(); msg=num; do_trap(); } /* do_trap() is in syscode.c */

void default_action(int num)
{ msg=num; untrapped_msg(); abnormal_exit(); }


SIGHANDLE handle_interrupt(n) int n;
{ signal(n,handle_interrupt); 
  if(in_tty_gets) hit_interrupt = TRUE; 
  else send(INTERRUPT); }

SIGHANDLE handle_quit(n) int n;
{ sprintf(QUIT_errmsg,"received QUIT signal (#%d) from system",n); send(QUIT);}

SIGHANDLE handle_matherror(n) int n;
{ signal(n,handle_matherror); sigcounter(); MATHERROR_arg1=MATHERROR_arg2=0.0; 
  strcpy(MATHERROR_type,"unknown math error (SIGFPE)"); send(MATHERROR); }

SIGHANDLE handle_buserror(n) int n;
{ signal(n,handle_buserror); sigcounter();
  strcpy(SYSERROR_errmsg,"segmentation-violation or bus-error"); 
  send(SYSERROR); }

SIGHANDLE handle_weird_signal(n) int n;
{ signal(n,handle_weird_signal); sigcounter();
  sprintf(SYSERROR_errmsg,"system signal #%d",n); send(SYSERROR); }


#ifdef MATHERR_DEFINED

int matherr(ex) 	 	/* System V std math error trap */
struct exception *ex;
{ switch (ex->type) {
    case DOMAIN:    strcpy(MATHERROR_type,"DOMAIN");      break;
    case SING:	    strcpy(MATHERROR_type,"SINGULARITY"); break;
    case OVERFLOW:  strcpy(MATHERROR_type,"OVERFLOW");    break;
    case UNDERFLOW: strcpy(MATHERROR_type,"UNDERFLOW");   break;
    case TLOSS:	    strcpy(MATHERROR_type,"Total-Loss of significance"); break;
    case PLOSS:	strcpy(MATHERROR_type,"Partial-Loss of significance");   break;
    default:	strcpy(MATHERROR_type,"Unknown");         break;
  }
  MATHERROR_arg1= ex->arg1; MATHERROR_arg2= ex->arg2;
  send(MATHERROR);
  return(0);  /*  never reached */
}

#endif

bool stack_check(var)
int *var;
{    
     if (*var==-1) { *var=lvl; return(TRUE); }
     if (lvl== *var) return(TRUE); 
     else return(FALSE);
}



