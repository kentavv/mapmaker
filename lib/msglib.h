#ifndef _MSGLIB_H_
#define _MSGLIB_H_

/******************************************************************************

 #    #   ####    ####   #          #    #####           #    #
 ##  ##  #       #    #  #          #    #    #          #    #
 # ## #   ####   #       #          #    #####           ######
 #    #       #  #  ###  #          #    #    #   ###    #    #
 #    #  #    #  #    #  #          #    #    #   ###    #    #
 #    #   ####    ####   ######     #    #####    ###    #    #

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

/*** Predefined Messages for sending. ***/
#define RERUN		1
#define PUNT		2
#define QUIT		3
#define MATHERROR	4
#define NOMEMORY	5
#define ENDOFILE	6
#define ENDOINPUT	7
#define IOERROR		8
#define CANTOPEN	9
#define CRASH		10
#define INTERRUPT	11
#define SYSERROR	12
#define ABORT		13
#define SOFTABORT	14
#define BADEQN 		15
#define CANTCLOSE	16

/* The first message number accessible to the user code is 20. */
#define USER_MESSAGE(n) (20+n)

/* Message variables for the above */
extern char *MATHERROR_type;    /* allocated LINE long */		
extern real  MATHERROR_arg1;
extern real  MATHERROR_arg2;
extern int   NOMEMORY_num_cells;
extern int   NOMEMORY_cell_size;
extern char *CANTOPEN_path;     /* allocated PATH_LENGTH long */
extern char  CANTOPEN_modechar; 
extern char *IOERROR_errmsg;    /* allocated LINE long */
extern char *IOERROR_linecopy;  /* set the ptr, ths is not allocated */
extern char *IOERROR_filename;  /* allocated PATH_LENGTH long */
extern char  IOERROR_modechar; 
extern char *SYSERROR_errmsg;   /* allocated LINE long */
extern char *SOFTABORT_msgstr;  /* a ptr, not allocated */
extern char *BADEQN_errmsg;
extern int   BADEQN_errpos;

/* user accessible stuff */
extern int msg;		/* the message # */		
extern char **mname;    /* [msg] => name for this msg # */
extern char *msgstr;    /* use for msg vars - is MAXLINE long */


void strmsg_default(char *str);
void strmsg_MATHERROR(char *str);
void strmsg_NOMEMORY(char *str);
void strmsg_CANTOPEN(char *str);
void strmsg_IOERROR(char *str);
void strmsg_SYSERROR(char *str);
void msg_init(void);
void setmsg(int var, char *nam, void (*action)(int), void (*disp)(char *));
int lvl_plus_plus(void);
void sender(int num);
void punter(void);
void trapper(int num);
void default_action(int num);
void handle_interrupt(int n);
void handle_quit(int n);
void handle_matherror(int n);
void handle_buserror(int n);
void handle_weird_signal(int n);
//int matherr(struct exception *ex);
bool stack_check(int *var);





//void setmsg();	/* sets up a message */
///* args int msg; char *msgname; void (*send_proc)(),(*message_proc)();
//   message_proc may be NULL */

//void untrapped_msg();   /* Deal with unhandled messages - in syscode.c */
//void trapped_msg();
//void verbose_untrapped_msg();

///* Possible send_procs - happen when we call send(MSG); args: message msg; */
//void sender(int num);
//void trapper();
//void punter();
//void default_action(int num); /* DON'T USE */

///* Possible message_procs in msglib.c. These produce an error message string
//   for untrapped messages. */
//void strmsg_default(char *);
//void strmsg_NOMEMORY();
//void strmsg_MATHERROR();
//void strmsg_CANTOPEN();
//void strmsg_SYSERROR();

/* Stuff for internal use only! - see msglib.c */
extern jmp_buf stk[];   		
extern int lvl;
extern bool exiting1;
extern void (*(maction[]))(int);
extern void (*(mstrmsg[]))(char *);
#define MSGS 	    50
#define MSGNAMLEN   39
#define TRAP_DEPTH  100
#define MAX_IO_FAILURES 20
#define MAX_BAD_SIGNALS 20
//void sigcounter();    /* count weird signals & CRASH if too many */
//int lvl_plus_plus();  /* return lvl for push onto stk[], chacking overflow */

/*****************************************************************************/
/* Trap Syntax: (Braces are USUALLY not required for single statements)
	
		run {
			blah;	 * Signal() statements can be in this
			blah;	 * block or in procedures called from
			blah;	 * it- break statements will work in here
			blah;	 * return is ABSOLUTELY NOT ALLOWED in here!
		} except {       *
		  when FOO:      * As this is really a switch statement,
			blah;	 * be sure each clause ends with a break,
			blah;	 * relay, or return instruction;
			break;   *
		  when BAR:
		        ...
		  relay_others;	 * This must exist for unhandled 
		}  		 * signals to be passed on!
*/


/* These are the basic syntactic elements...
   The while condition here is never true! (e.g. like while(FALSE)) We use the
   while so that break statements will work inside run {...} blocks. */

#define run 	if ((msg=setjmp(stk[lvl_plus_plus()]))==0) { do
#define except	while(--lvl>10000); } else switch(msg)	
  
#define send(num) 	(*(maction[num]))(num)
#define when 		case
#define relay   	sender(msg)
#define relay_others 	default: sender(msg)


/* An abbreviation to catch only one signal- all others are relayed. 
   Syntax: 

	run { 
	... 
	} except_when(FOO) {	* braces are not required for 1 statement here.
	...		 	* break can be used in an except_when clause
	}			* to stop a loop surrounding the run... thing
*/

#define except_when(num) \
  while(--lvl>10000); } else if (msg!=num) { sender(msg); } else

/* Trapping() provides nicer syntax for the trivial do-nothing traps.
   Syntax:

	run { 
	... 
        } trapping(FOO);   	* the same as except_when(FOO) {}
*/

#define trapping(m) except_when(m) 


/* Another alternative to the except ... when notation shown above. 
   In this case, one block  of code is used to handle all messages, although
   this block may look at the msg variable to tell what happened.
   Note that the break statement should  work correctly in an 
   when_aborting { ... } block. 
   Syntax:

	run { 
	    ...
	} when_aborting { 
	    if (msg==FOO) { 
	    ...
	    } else if (msg==BAR) { 
	    ...
	    } else relay; 		* remember to add this if you don't
	}				* just want to fall out.
*/

#define when_aborting	 while(--lvl>10000); } else 
#define on_error	 while(--lvl>10000); } else 


/* Yet another way to do traps. On_exit is similar to the
   when_aborting construct shown above, except that the on_exit clause is
   always executed, whether there was a message or whether it was just
   "fallen into". As with when_aborting, the variable msg will be set
   within this code correctly (msg==0 if no msg). Note that the break
   statement should also work correctly in an on_exit { ... } block.

   For example, the following code will never leave a chunk of malloced space
   around. If array() fails, NOMEMORY is sent to the caller, after an error
   message is printed. This code behaves correctly, freeing all malloced
   space, even if the first array() call succeeds and the second fails!

	p=q=NULL;  * This must be here as the pointers may not be initialized, 
	...        * which might cause us to try to unarray a bad pointer!
	run {
	    array(p, ... );
	    array(q, ... );
	    ...
	} on_exit {
	    unarray(p,...); 
	    unarray(q,...);
	    relay_messages; 
	}
*/

#define on_exit			\
     while(--lvl>10000);	\
  } 				\
  for (exiting1=TRUE; exiting1; exiting1=FALSE)

#define relay_messages sender(msg) /* nicer syntax than relay for this */

/* We use stack_check() to check to see if the run stack is where it
   was last time this was called with the given variable (an int). Use
   init_stack_check() to set it sometime before calling stack check().
   Remember to call stack_check() from some toplevel loop which should
   not be building up run...except calls.  If stack_check() fails, the
   program should be aborted! */

//bool stack_check();       /* args: int *var; return TRUE if all is OK */
//void init_stack_check();  /* args: int *var; */


/******************** Support for Unix style signals ********************/

///* Signal trap functions for signal_trap_init (which is in syscode.c) */
//void handle_quit();
//void handle_matherror();
//void handle_weird_signal();
//void handle_interrupt();
//void handle_buserror();
//
//void msg_init();


/***** THE FOLLOWING IS AS YET UNUSED, AND MAYBE BEST LEFT THAT WAY *****/

/* After any run...except... thing, msg will be set to either 0 or a
valid message number. Therefore, we can use it as a state variable to
keep track of whether we aborted or fell through to after the except { }
clause..  However use this with extreme care (eg: only use it in a
function AFTER executing a run...except... thing). */

#define relay_trapped_msg   sender(msg)

extern bool in_tty_gets;
extern bool hit_interrupt;

#endif