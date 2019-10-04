/******************************************************************************

  ####   #    #  ######  #       #                ####
 #       #    #  #       #       #               #    #
  ####   ######  #####   #       #               #
      #  #    #  #       #       #        ###    #
 #    #  #    #  #       #       #        ###    #    #
  ####   #    #  ######  ######  ######   ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#define INC_LIB
#define INC_MISC
#define INC_HELP_DEFS
#include "system.h"
#include "shell.h"
#include "table.h"

char *com, *args;       /* available as globals to the command procedures */
char *uncrunched_args;
int com_num, num_args;
char *((*prompt)());
char *default_prompt();
bool (*quit_save_hook)();
bool more_mode;

/***** defs for internal use only *****/
#define UNHELPFUL	(-1)        /* as a help_entry */
#define HELPLESS	((long) -1) /* as a help_key */
#define MAX_SUCCESSIVE_FAILURES 15

#define MAX_NUM_MENUS    20
#define MAX_MENU_ENTRIES 20
#define MENU_ENTRY_LEN   20

COMMAND **cmd;		/* [command#] => ptr to a COMMAND struct */
MENU **menu;            /* [menu#] => ptr to a menu struct */
int num_menus;
int command_num; 	/* num commands in cmd */
int help_entries;	/* number of entries in help_file */
char **topic_name;   	/* [topic_num] set by mktopic() */
int *topic_code;   	/* [topic_num] set by mktopic() */
long *topic_help_key;   /* [topic_num] set by mktopic() */
char **tokens;		/* [token_num] used as a temp by the parser */
char **remaining;  	/* ditto */
bool *matched; 		/* [command_num] parser temp */
int num_matched;
bool wizard_mode;
char *save_args_ptr;    /* the arg list, saved so the ptr can be bashed */
char *save_uncrunched_args;
TABLE *cmd_history;
int autosave;           /* saves state on exit if TRUE */
FILE *help_file;        /* file that contains all the help summaries */
bool user_is_amused;    /* two state vars for keep_user_amused */

bool shell_uses_wimp_cmds, help_uses_wimp_help; 
bool inhibit_menus;
/* If shell_uses_wimp_cmds==TRUE, a command entered to the prompt will
invoke its WIMP version, if one is available, otherwise its normal
text version will be run. help_uses_wimp_help is similar. The WIMP 
code should only try to run a command when inhibit_menus==FALSE. */

int cmd_history_num;
bool photo_update_top_hook;
void (*photo_banner_hook)();
char *the_program, *the_version, *the_copyright;

int alias_match();
int abbrev_match();
int com_matches();
int abbrev_matches();
int try_to_match();

void expand_history_references();
char *centering();

#define NO_HELP_FILE "\
Can't find help file - detailed help information is not available.\n\
See installation instructions for details.\n"

#define NO_HELP_KEY "No additional information available.\n"

#define SURROGATE_ABOUT "\
This program is freely redistributable under certain conditions and is\n\
licensed free of charge for non-commercial applications. No warranty of any\n\
type is provided. See the License Agreement for details.\n"

#define TYPE_HELP_PLEASE "\
Type 'help' for help.\n\
Type 'about' for license, non-warranty, and support information.\n"

void banner() 
{ 
    char line[81];

print(
 "************************************************************************\n");
print(
 "* Welcome to:                                                          *\n");

print(centering("",FALSE));
print(centering(the_program,FALSE));
if (!nullstr(the_version)) 
  { sprintf(line, "(version %s)", the_version); print(centering(line, FALSE)); }
print(centering("",FALSE));
sprintf(line, "Copyright %s, Whitehead Institute for Biomedical Research",
        the_copyright); print(centering(line,FALSE));
if (gnu_copyright(line)) print(centering(line,TRUE));

print(
 "************************************************************************\n");

nl();
print(TYPE_HELP_PLEASE);
if (help_file==NULL) print(NO_HELP_FILE);
}


char *get_version(version_filename)
char *version_filename;
{
    FILE *fp;
    char str[TOKLEN+1], *p, *version, name[PATH_LENGTH+1];

    if (nullstr(version_filename)) return(NULL);
    fp=NULL; version=NULL;

    run {
	strcpy(name,version_filename);
	if (!make_filename_in_dir(name,FORCE_EXTENSION,".v",
				  FORCE_DIR,CODE_DIR)) send(CANTOPEN); 
	fp=open_file(name,READ);
	finput(fp,ln_,MAXLINE); p=ln_; /* don't hack &ln_ */
	close_file(fp);
	if (stoken(&p,sREQUIRED,str)) version=mkstrcpy(str);
    } when_aborting {
	close_file(fp); /* Don't bother relaying messages */
    }
    return(version);
}
    

void photo_banner()
{
    char ver[81];
/*                                         Thu Apr 13 05:31:11 EDT 1989
 123456789012345678901234567890123456789012345678901234567890123456789012 */
lib_puts(photo,
 "************************************************************************\n");
sprintf(ps, "* Output from:                                %24s *\n", time_string());
lib_puts(photo,ps);

lib_puts(photo,centering("",FALSE));
lib_puts(photo,centering(the_program,FALSE));
sprintf(ver, "(version %s)", the_version); lib_puts(photo, centering(ver, FALSE));
lib_puts(photo,centering("",FALSE));
lib_puts(photo,
 "************************************************************************\n");
lib_puts(photo,"\n");
if (photo_banner_hook!=NULL) (*photo_banner_hook)(photo);
}


int prev_spaces=2;

char *centering(str,lineup)  /* into global ps */
char *str;
bool lineup;
{ 
    int spaces, i; 
    if (!lineup) spaces= 36-(len(str)/2); else spaces=prev_spaces;
    prev_spaces= spaces;
    strcpy(ps,
"*                                                                      *\n");
    for (i=0; str[i]!='\0'; i++) ps[spaces+i]=str[i];
    return(ps);
}


void shell_init(program,version,copyright,help_filename)
char *program, *version, *copyright, *help_filename; 
/* help_filename must be side-effectable */
{
    int i;

    array(cmd, MAX_COMMANDS, COMMAND*);
    array(matched, MAX_COMMANDS, bool);
    array(topic_name, MAX_COM_TOPICS+1, char*); 
    array(topic_code, MAX_COM_TOPICS+1, int); 
    array(topic_help_key, MAX_COM_TOPICS+1, long); 
    matrix(tokens, MAX_COM_TOKENS, TOKLEN+1, char);
    array(remaining, MAX_COM_TOKENS+1, char*);
    array(args, ARGS_LEN+1, char); 
    array(uncrunched_args, ARGS_LEN+1, char); 
    save_args_ptr= args; save_uncrunched_args= uncrunched_args;
    
    the_program= mkstrcpy(program);
    the_version= mkstrcpy(version);
    the_copyright= mkstrcpy(copyright);

    for (i=1; i<MAX_COM_TOPICS+1; i++) 
      { topic_name[i]=NULL; topic_code[i]=0; topic_help_key[i]=HELPLESS; }
    topic_name[0]= mkstrcpy("TOPIC ZERO");
    topic_code[0]=0; topic_help_key[0]=HELPLESS;

    cmd_history= allocate_table(21,250,CANT_EXPAND,INDEX_BY_NUMBER);
    cmd_history_num= next_entry_number(cmd_history);

    command_num=0;
    help_entries=0;
    wizard_mode=FALSE;
    photo_update_top_hook=TRUE;
    photo_banner_hook=NULL;
    prompt=default_prompt;
    inhibit_menus=TRUE;
    quit_save_hook=NULL;
    user_is_amused=FALSE;

    parray(menu, MAX_NUM_MENUS, MENU);
    for (i=0; i<MAX_NUM_MENUS; i++) 
      array(menu[i]->title,MENU_ENTRY_LEN+1,char);
    shell_uses_wimp_cmds= FALSE;
    help_uses_wimp_help=  FALSE;

    help_file=NULL;
    if (make_filename_in_dir(help_filename,FORCE_EXTENSION,HELP_EXT,
			     DEFAULT_DIR,CODE_DIR)) {
	run help_file=open_file(help_filename,READ);
	  except_when(CANTOPEN) {}
    }
}


char *default_prompt(s)
char *s; 
{ sprintf(s, "\n%d> ", cmd_history_num + 1); return(s); }


void mktopic(num,nam,code,description_index)
int num;
char *nam;
int code;
long description_index;
{ 
    if (num<=0 || num>MAX_COM_TOPICS || len(nam)>MAX_TOPIC_LEN) send(CRASH); 
    topic_name[num]=mkstrcpy(nam); topic_code[num]=code;
    topic_help_key[num]= description_index;
}


int mkcommand(name,abbrev,func,code) 
char *name, *abbrev;
void (*func)();
int code;
{
    char *str, **toks, c;
    int i, n_tokens;
    
    single(cmd[command_num], COMMAND);
    str= mkstrcpy(name); crunch(str); truncstr(str,MAX_COM_NAME_LEN);
    cmd[command_num]->name= mkstrcpy(str); /* str itself => toks, below */
    
    array(toks,MAX_COM_TOKENS,char*);
    for (i=0; i<MAX_COM_TOKENS; i++) toks[i]=NULL;
    toks[0]= str;

    for (i=0, n_tokens=1; str[i]!='\0'; i++)
      if ((c=str[i])==' ' || c=='-' || c=='_') {
	  str[i]= '\0';
	  if (n_tokens<MAX_COM_TOKENS) toks[n_tokens++]= &str[i+1];
	  else { cmd[command_num]->name[i]= '\0'; break; }
      }
	
    cmd[command_num]->tokens= toks;
    cmd[command_num]->num_tokens= n_tokens;
    cmd[command_num]->procedure= func;
    cmd[command_num]->code= code;

    nstrcpy(cmd[command_num]->abbreviation,abbrev,MAX_ABBREV_LEN);
    crunch(cmd[command_num]->abbreviation);

    cmd[command_num]->topic= 0;
    cmd[command_num]->help_key= HELPLESS;
    cmd[command_num]->help_entry= UNHELPFUL;
    cmd[command_num]->cmd_help= NULL;
    cmd[command_num]->args_help=NULL;
    cmd[command_num]->num_args= -1;
    cmd[command_num]->num_args_prefix= 0;

    cmd[command_num]->wimp_procedure= NULL;
    cmd[command_num]->status_function= NULL;
    cmd[command_num]->wimp_help= NULL;
    cmd[command_num]->wimp_menu_num= -1;
    cmd[command_num]->menu_entry= NULL;
    cmd[command_num]->wimp_shortcut= '\0';

    return(command_num++);
}

		
void mkhelp(cmd_name,abbrev,description_index,num_args_prefix,num_args,
	    code,topic,description,arguments,defaults)
char *cmd_name, *abbrev;
long description_index;
int num_args_prefix, num_args;
int topic, code;
char *description, *arguments, *defaults;
{
    int i, old_wiz;
    char *rest;

    old_wiz=wizard_mode; wizard_mode=TRUE;
    i=parser(cmd_name,&rest,TRUE);
    wizard_mode= old_wiz;

    if (i<0 || !nullstr(rest)) {
	if (code!=HLP) {
	    sprintf(ps, "warning: attempt to make help for non-command '%s'\n",
                cmd_name); pr();
	    return;
	} else i=mkcommand(cmd_name,abbrev,NULL,HLP);
    }

    if (!streq(cmd[i]->name,cmd_name)) {
	sprintf(ps, "warning: names disagree for command '%s'\n", cmd_name);
	pr();
    }
    if (code!=cmd[i]->code) {
	sprintf(ps, "warning: type codes disagree for command '%s'\n", cmd_name);
	pr();
    }
    if (!streq(cmd[i]->abbreviation,abbrev)) {
	sprintf(ps, "warning: abbreviations disagree for command '%s'\n", cmd_name);
	pr();
    }

    array(cmd[i]->cmd_help, MAX_COM_HELP_LEN+1,char);
    array(cmd[i]->def_help, MAX_ARG_HELP_LEN+1,char);
    array(cmd[i]->args_help,MAX_ARG_HELP_LEN+1,char);
    
    nstrcpy(cmd[i]->cmd_help,description,MAX_COM_HELP_LEN);
    nstrcpy(cmd[i]->args_help,arguments,MAX_ARG_HELP_LEN); 
    nstrcpy(cmd[i]->def_help,defaults,MAX_ARG_HELP_LEN);
    
    cmd[i]->code=  code;
    cmd[i]->topic= topic;
    cmd[i]->help_key= description_index;
    cmd[i]->help_entry= help_entries;
    cmd[i]->num_args= num_args;
    cmd[i]->num_args_prefix= num_args_prefix;

    ++help_entries;
}


bool valid_name(str) /* checks the syntax of names */
char *str;
{ 
  int i;

  if (!is_a_token(str)) return(FALSE);
  if (strin(NAME_TAG_CHARS,str[0])) str++;
  if (!strin(NAME_FIRST_CHARS,str[0])) return(FALSE);
  for(i=1; str[i]!='\0'; i++) if (!strin(NAME_CHARS,str[i])) return(FALSE);
  return(TRUE);
}


void null_command() { send(CRASH); }



/**************************** The Command Parser ****************************/

int parser(line,rest,help_ok) 
char *line;
char **rest; /* side-effected */
bool help_ok;
{ 
    int i, n_tokens, last_match;
    char *foo;
    extern char **tokens, **remaining; /* temps */
    extern int *matched, num_matched;
	
    for (i=0; i<command_num; i++) matched[i]= FALSE;
    num_matched= 0;
    if (rest==(char**)NULL) rest= &foo;

    remaining[0]= line;
    for (n_tokens=0; n_tokens<MAX_COM_TOKENS; n_tokens++) 
      if (stoken(&line,sREQUIRED,tokens[n_tokens])) {
	  remaining[n_tokens+1]= line;
	  if (remaining[n_tokens+1][0]==' ')
	    remaining[n_tokens+1]+=1;
      } else break;

    for (i=0; i<command_num; i++) /* try to exactly match an abbreviation */
      if (allowed_cmd(i,help_ok) && xstreq(tokens[0],cmd[i]->abbreviation)) 
	{ matched[i]=TRUE; num_matched=1; *rest=remaining[1]; return(i); }

    try_to_match(tokens,n_tokens,3,3,matched,&num_matched,&last_match,help_ok);
    if (num_matched==1) { *rest=remaining[3]; return(last_match); }

    try_to_match(tokens,n_tokens,2,2,matched,&num_matched,&last_match,help_ok);
    if (num_matched==1) { *rest=remaining[2]; return(last_match); }
    try_to_match(tokens,n_tokens,2,3,matched,&num_matched,&last_match,help_ok);
    if (num_matched==1) { *rest=remaining[2]; return(last_match); }

    try_to_match(tokens,n_tokens,1,1,matched,&num_matched,&last_match,help_ok);
    if (num_matched==1) { *rest=remaining[1]; return(last_match); }
    try_to_match(tokens,n_tokens,1,2,matched,&num_matched,&last_match,help_ok);
    if (num_matched==1) { *rest=remaining[1]; return(last_match); }
    try_to_match(tokens,n_tokens,1,3,matched,&num_matched,&last_match,help_ok);
    if (num_matched==1) { *rest=remaining[1]; return(last_match); }

    *rest= remaining[0]; 
    return(-1);
}


void print_parser_results(rest,help_ok) /* uses the global state */
char *rest;
bool help_ok;
{
    int j;

    if (num_matched==1) { 
	for (j=0; j<command_num; j++) 
	  if (matched[j]) sprintf(ps, "matched command '%s'\n", cmd[j]->name);
	pr();

    } else if (num_matched>1) {
	sprintf(ps, "ambiguous command%s: could be any of the following:\n",
            (help_ok ? " or help topic":"")); pr();
	for (j=0; j<command_num; j++) 
	  if (matched[j]) { print("\t"); print(cmd[j]->name); nl();
	}

    } else { /* num_matched==0 */
	sprintf(ps, "unrecognized command%s: '", (help_ok ? " or help topic" : ""));
	pr();
	if (len(rest)>LINE-25) strcpy(&rest[LINE-28],"..."); 
	print(rest); print("'\n");
    }
}


/* Try to match n_to_try the input tokens against all n_word_commands, and
   set com_match, n_matched and last_match accordingly. Return TRUE only if a
   unique match was found. */
  
int try_to_match(token,n_tokens,n_to_try,n_word_command,
	com_match,n_matched,last_match,allow_help_only_stuff)
char **token;
int n_tokens, n_to_try, n_word_command;
int *com_match, *n_matched, *last_match, allow_help_only_stuff;
{
	int exact, i;
	
	if (n_tokens<n_to_try) return(FALSE);
	for (i=0; i<command_num; i++) 
	    if (allowed_cmd(i,allow_help_only_stuff) &&
		cmd[i]->num_tokens==n_word_command &&
		com_matches(token,cmd[i]->tokens,n_to_try,&exact)) {
		        *last_match= i;
			if (!com_match[i]) 
				{ com_match[i]= TRUE; ++*n_matched; }
			if (exact && n_to_try==n_word_command) return(TRUE);
		}
	if (*n_matched==1) return(TRUE); else return(FALSE);
}


int com_matches(in_tokens,com_tokens,num_to_match,exact)
char **in_tokens, **com_tokens;
int num_to_match, *exact;
{
    int i;
	
    for (i=0, *exact=TRUE; i<num_to_match; i++) 
      if (!matches(in_tokens[i],com_tokens[i])) return(FALSE);
      else if (istrlen(in_tokens[i])!=istrlen(com_tokens[i])) *exact=FALSE;

    return(TRUE);
}


/********************************* The Shell *********************************/


#define HIST_RANGE "command history number is out of range"

void expand_history_references(line) /* sends an error if need be */
char *line; /* line IS side-effected */
{
    char *str, *save;
    int num, first, i;

    /* The cmd_history num will NOT have been incremented between the call
       to prompt and here! Also remember that 1+ that number is printed. */

    save=line;
    if (!itoken(&line,iREQUIRED,&num)) return;
    for (i=0; save+i!=line; i++) save[i]=' ';

    if (get_numbered_entry(num-1,&str,cmd_history)) {
	/* maxstrins(line," ",MAXLINE); */
	maxstrins(line,str,MAXLINE); 
	despace(line);
	print("=");
	if (cmd_history_num>9) print("=");
	if (cmd_history_num>99) print("=");
	print("> ");
	print(line); nl();
	expand_history_references(line);

    } else {
	/* THIS IS A KLUDGE, and results in a 'statement not reached' error. */
	for(Te=cmd_history->list; Te!=NULL; Te=Te->next)
	  { first=Te->id.num; break; } /*NOTREACHED*/ 
	if (cmd_history_num>1)
	  sprintf(ps, "%s\nUse a number from %d to %d.", HIST_RANGE,
	     first+1, cmd_history_num);
	else if (cmd_history_num==1) 
	  sprintf(ps, "%s\nThe only valid command number yet is 1.", HIST_RANGE);
	else if (cmd_history_num==0) 
	  sprintf(ps, "%s\nNo commands have yet been entered.", HIST_RANGE);
	error(ps); return; /* never returns from here */
    }
}


#define nopunt  failures++; break;
#define punt    done=TRUE;  break;
#define toleft  if (cursor!=0) nl();

void command_loop()
{
    char *rest; 
    int done, failures, foo, prev_lvl= -1;
    void (*func)();
	
    done=FALSE; failures=0;
    do {
	/* run { */
	if ((msg=setjmp(stk[lvl_plus_plus()]))==0) { 
	    do {
		if (prev_lvl<0) prev_lvl=lvl;
		  else if (prev_lvl!=lvl) {
		      fprintf(stderr,"Internal Error: stack level changed\n");
		      prev_lvl= lvl;
		  }
		if (log_open) fflush(photo);

		inhibit_menus= FALSE;
		cmd_history_num= next_entry_number(cmd_history);
		uncrunched_args= save_uncrunched_args;
		input((*prompt)(ps),uncrunched_args,ARGS_LEN-1); 
		despace(uncrunched_args); /* should already be filtered */

		if (!nullstr(uncrunched_args)) {
		    expand_history_references(uncrunched_args); 
		    put_numbered_entry(uncrunched_args,cmd_history,&foo);

		    if ((com_num=parser(uncrunched_args,&rest,FALSE))>=0) { 
			uncrunched_args=rest; args=save_args_ptr; 
			strcpy(args,uncrunched_args); lowercase(args);
			if (wimp && shell_uses_wimp_cmds &&
			    cmd[com_num]->wimp_procedure!=NULL) 
		          func=cmd[com_num]->wimp_procedure; 
			else func=cmd[com_num]->procedure;
			if (func==NULL) send(CRASH);
			com=cmd[com_num]->name; 
			num_args=0; inhibit_menus=TRUE;
			(*func)(); /* Run the command */
			failures=0;
		    } else print_parser_results(rest,FALSE); /* error msg */
		}
	    } while(--lvl>10000); 

	    /* } except { */
	} else switch(msg) {
	    when INTERRUPT:  
	    	print("\n*** break ***\n");
	    	if (redirecting_input) {
		    print("\n\t...run cancelled...\n");
		    redirect_input(NULL,FALSE);
		}
	    	nopunt;

	    when NOMEMORY:   print("\n*** out of memory ***\n"); nopunt;
	    when MATHERROR:  print("\n*** floating point error ***"); nopunt;
	    when ENDOINPUT:  toleft; print("\n\t...end of input...\n\n"); punt;
	    when QUIT:       toleft; print("\n\t...goodbye...\n\n"); punt;
	    when SOFTABORT:  nopunt;

	    when IOERROR:
	    when CRASH:
	    when SYSERROR:
	    when CANTOPEN:
	    when ENDOFILE:
	    default:	      verbose_untrapped_msg(); nopunt;
	}

	if (failures>MAX_SUCCESSIVE_FAILURES) {
	    fprintf(stderr,"\n*** too many successive errors: aborting ***\n");
	    done= TRUE;
	}
    } while(!done);
}



/************************ Useful things for commands ************************/


void abort_command() { send(SOFTABORT); }


void error(const char *errmsg) /* guaranteed not to use ps */
{ print("error: "); print(errmsg); nl(); abort_command(); }


void maybe_set_bool(var)
bool *var;
{
    char temp[TOKLEN+1];
    
    if (stoken(&args,sREQUIRED,temp)) {
        if (streq(temp,"on")) *var= TRUE;
	else if(streq(temp,"off")) *var=FALSE;
	else set_usage_error("either 'on' or 'off'");
    }
    sprintf(ps, "'%s' is %s.\n", com, *var ? "on" : "off"); pr();
}


void maybe_set_real(var,lbound,hbound,fmt)
real *var, lbound, hbound;
real fmt;
{
    real temp;
	
    if (!nullstr(args)) {
	if (!rtoken(&args,rREQUIRED,&temp) || !rrange(&temp,lbound,hbound)) {
	    sprintf(ps, "a real number from %s to %s", rs(fmt, lbound), rs(fmt, hbound));
	    set_usage_error(ps); 
	} else *var= temp;
    }
    sprintf(ps, "'%s' is %s\n", com, rs(fmt, *var)); pr();
}


void maybe_set_long(var,lbound,hbound)
long *var, lbound, hbound;
{
	long temp;
	
	if (!nullstr(args)) {
	    if (!ltoken(&args,lREQUIRED,&temp) ||
		!lrange(&temp,lbound,hbound)) {
		   sprintf(ps, "an integer number from %ld to %ld", lbound, hbound);
		   set_usage_error(ps);
	    } else *var= temp;
	}
	sprintf(ps, "'%s' is %ld\n", com, *var); pr();
}


void maybe_set_int(var,lbound,hbound)
int *var, lbound, hbound;
{
	int temp;
	
	if (!nullstr(args)) {
	    if (!itoken(&args,iREQUIRED,&temp) ||
		!irange(&temp,lbound,hbound)) {
		   sprintf(ps, "an integer number from %d to %d", lbound, hbound);
		   set_usage_error(ps);
	    } else *var= temp;
	}
	sprintf(ps, "'%s' is %d\n", com, *var); pr();
}


void set_usage_error(com_args) /* guaranteed not to use ps */
char *com_args;
{ print("error: illegal value for '"); print(com); print("'\n");
  if (cmd[com_num]->args_help!=NULL) 
    { print("correct value: "); print(cmd[com_num]->args_help); print("\n"); }
  print("type '"); print(com); 
  print("' alone to display current value\n"); 
  print("type 'help "); print(com); print("' for details\n"); 
  abort_command(); 
}


void usage_error(num) 
int num; /* num args given, maybe <0 */
{ 
    if (cmd[com_num]->num_args==0) 
      sprintf(ps, "error: The '%s' command takes no arguments.\n", com);
    else 
      sprintf(ps, "error: Missing%s argument(s) for the '%s' command.\n",
              (num!=0 ? " or invalid":""), com);
    pr();

    if (cmd[com_num]->num_args!=0) {
	if (!nullstr(cmd[com_num]->args_help))
	  { sprintf(ps, "expected:  %s\n", cmd[com_num]->args_help); pr(); }
	if (!nullstr(cmd[com_num]->def_help))
	  { sprintf(ps, "default%s %s\n", (cmd[com_num]->num_args == 1 ? ": " : "s:"),
                cmd[com_num]->def_help); pr(); }
    }

    if (cmd[com_num]->help_key!=HELPLESS)
      { sprintf(ps, "type 'help %s' for details.\n", com); pr(); }

    abort_command(); 
}      


#define TOOMANY \
"error: Too many arguments for the '%s' command\n(%s%d argument%s expected).\n"

void nomore_args(n)
int n; /* for now, a dummy arg */
{ 
    int mode, num; 
    num= cmd[com_num]->num_args;
    mode=cmd[com_num]->num_args_prefix;
    if (nullstr(args)) return; /* all is OK, otherwise... */
    
    if (num==0) 
      sprintf(ps, "error: The '%s' command takes no arguments.\n", com);
    else if (num>0)
      sprintf(ps, TOOMANY, com, (mode == EXACTLY ? "" : "up to "), num, maybe_s(num));
    else /* num<0 */ sprintf(ps, "error: Too many arguments.\n");
    pr();
    
    if (!nullstr(cmd[com_num]->args_help))
      { print("expected: "); print(cmd[com_num]->args_help); nl(); }
    
    if (!nullstr(cmd[com_num]->def_help))
      { sprintf(ps, "default%s: %s\n", maybe_s(num), cmd[com_num]->args_help); pr(); }
    
    if (cmd[com_num]->help_key!=HELPLESS)
      { sprintf(ps, "Type 'help %s' for details.\n", com); pr(); }
    
    abort_command();
}


void more_args(num)
int num; /* num args given, maybe <0 */
{
    int mode, want; 
    want= cmd[com_num]->num_args; 
    mode=cmd[com_num]->num_args_prefix;

    if (!nullstr(args)) usage_error(-1); /* instead */

    sprintf(ps, "error: Missing argument%s.\n", maybe_s(want - num)); pr();
    sprintf(ps, "The '%s' command requires %s%d argument%s.\n", com,
       mode==EXACTLY ? "" : "at least ", want, maybe_s(want)); pr();
    
    if (!nullstr(cmd[com_num]->args_help))
      { print("expected: "); print(cmd[com_num]->args_help); nl(); }
    
    if (!nullstr(cmd[com_num]->def_help))
      { sprintf(ps, "default%s: %s\n", maybe_s(num), cmd[com_num]->args_help); pr(); }
    
    if (cmd[com_num]->help_key!=HELPLESS)
      { sprintf(ps, "Type 'help %s' for details.\n", com); pr(); }

    abort_command();
}


void input_error(val,def) 
char *val, *def;
{ 
    print("error: you have given an invalid response\n"); 
    if (!nullstr(val))
      { print("expected: "); print(val); nl(); }
    if (!nullstr(def))
      { print("default:  "); print(def); nl(); }

  if (cmd[com_num]->help_key!=HELPLESS)
    { sprintf(ps, "try typing 'help %s' for details\n", com); pr(); }

  abort_command(); 
}      


#define EXTRA_INPUT   "too many values given in input line"

void expect_nomore_input(str,mode,num)
char *str;
int mode, num;
{ if (nullstr(str)) return;
  else if (num>0)
    sprintf(ps, "error: %s\n %s %d value%s were expected\n", EXTRA_INPUT,
            (mode==EXACTLY ? "only" : "up to"),
            num, maybe_s(num));
  else sprintf(ps, "%s\n", EXTRA_INPUT);
  pr(); print("try 'help "); print(com); print("' for details\n"); 
  abort_command();
}


void expect_more_input(str,mode,num)
char *str;
int mode, num;
{ if (!nullstr(str)) return;
  else if (num>0)
    sprintf(ps, "error: missing input\n%s %d value%s expected\n",
            (mode==EXACTLY ? "":"at least"), num, (num>1 ? "s were":" was"));
  else sprintf(ps, "error: missing input\n");
  pr(); print("try 'help "); print(com); print("' for details\n"); 
  abort_command();
}


bool split_arglist(rest,divider)
char **rest;
char divider; /* character that arglist should be split on */
{
    int i, j;

    *rest=NULL; i=0;
    for (i=0; args[i]!='\0'; i++) 
      if (args[i]==divider) {
	  j=i+1; while (white(args[j])) j++; *rest=args+j;
	  i--;   while (white(args[i]) && i>0) i--;  args[i+1]='\0';
	  return(TRUE);
      }
    return(FALSE);
}


bool split_uncrunched_args(rest,divider)
char **rest;
char divider; /* character that arglist should be split on */
{
    int i, j;

    *rest=NULL; i=0;
    for (i=0; uncrunched_args[i]!='\0'; i++) 
      if (uncrunched_args[i]==divider) {
	  j=i+1; while (white(uncrunched_args[j])) j++; 
	  *rest=uncrunched_args+j;
	  i--;   while (white(uncrunched_args[i]) && i>0) i--; 
	  uncrunched_args[i+1]='\0';
	  return(TRUE);
      }
    return(FALSE);
}


void maybe_ok(str)
char *str;
{
    if (update_top()) {	
	if (photo!=NULL) lib_puts(photo,str); 
	lib_puts(out,"ok");
    } else print(str);
    nl();
}
    

void keep_user_amused(thing,iter,max_iter)
char *thing; /* nullstr(str) means we are done, len<<TOKLEN, like one word */
int iter, max_iter;
{
    char simple[TOKLEN+1], fancy[TOKLEN+1];

    if (iter==0) { 
	usertime(TRUE);
	/* only need to do the simples the first time */
	if (max_iter==0) sprintf(simple, "computing %s%s", thing, "s");
	  else sprintf(simple, "computing %d %s%s", max_iter, thing, maybe_s(max_iter));
    } else strcpy(simple,"foo");

    if (max_iter==0) sprintf(fancy, "%s %d", thing, iter);
      else sprintf(fancy, "%s %d of %d", thing, iter, max_iter);

    if (iter==0 || usertime(FALSE)>1.0)
      { temp_print(simple,fancy); usertime(TRUE); }
}


/**************************** Some Basic Commands ****************************/

command show_cmd_history()
{
    int i, num_to_print, first_to_print, printed_any; 
    char *cmd_str;

    get_one_arg(itoken,1000,&num_to_print);
    if (num_to_print<1) usage_error(1);

    printed_any=FALSE;
    first_to_print= imaxf(cmd_history_num-num_to_print,0);
    for(Te=cmd_history->list; Te!=NULL; Te=Te->next) {
	i=Te->id.num;  cmd_str=Te->string; 
	if (i>=first_to_print && i<cmd_history_num) {
	    if(!printed_any){print("Previous commands:\n"); printed_any=TRUE;} 
	    sprintf(ps, "%3d  ", i + 1); pr(); print(cmd_str); nl();
	}
    }
    if (!printed_any) print("No commands have yet been entered.\n");
}


command quit()
{
    char token[TOKLEN+1];

    if (interactive && !redirecting_input) {
	/* test auto_save && data_loaded */
	if (quit_save_hook==NULL || !(*quit_save_hook)(FALSE)) send(QUIT);
	run getln("save data before quitting? [yes] ");
	except { 
	    when ENDOINPUT: ln[0]='y'; break;
	    when INTERRUPT: abort_command(); break;
	    default: 	    relay; break;
	}
	if (!stoken(&ln,sREQUIRED,token) || !matches(token,"no")) {
	    if (!((*quit_save_hook)(TRUE))) return; /* if saving fails */
	      else send(QUIT); 
	} else send(QUIT);
    } else {
	if (quit_save_hook!=NULL && !((*quit_save_hook)(TRUE))) return;
	send(QUIT);
    }
}


command really_quit()
{ send(QUIT); }


command run_from_file()
{
    char file_name[PATH_LENGTH+1];

    use_uncrunched_args();
    get_one_arg(stoken,sREQUIRED,file_name);

    if (!make_filename(file_name,DEFAULT_EXTENSION,"in")) 
      error("bad input file name");
    else redirect_input(file_name,TRUE); /* verbose -> messages */
}


command do_photo()
{
    char file_name[TOKLEN+1];

    use_uncrunched_args();
    get_one_arg(stoken,"",file_name);

    if (streq(file_name,"off")) { 
	photo_to_file("","");
	sprintf(ps, "'%s' is off.\n", com); pr();
	return;

    } else if (streq(file_name,"on")) {
	usage_error(1);

    } else if (!nullstr(file_name)) {
	if (!make_filename(file_name,DEFAULT_EXTENSION,"out")) 
	  sprintf(ps, "error: Bad photo file name");
	if (!photo_to_file(file_name,"a")) 
	  sprintf(ps, "error: Unable to open photo file '%s'\n", file_name);
	else {
	    photo_banner(); 
	    if (photo_update_top_hook && update_top()) sprintf(ps, "ok\n");
	    else sprintf(ps, "'%s' is on: file is '%s'\n", com, photo_file);
	}
	pr(); return;

    } else { /* no args */
	if (log_open) {
	    fflush(photo);
	    sprintf(ps, "'%s' is on: file is '%s'\n", com, photo_file);
	} else sprintf(ps, "'%s' is off\n", com);
	pr(); return;
    }
}


command set_more()   { maybe_set_bool(&more); }
command set_wizard() { maybe_set_bool(&wizard_mode); }
command set_verbose_mem() { maybe_set_bool(&verbose_mem); }


#define s_or_space_colon (cmd[i]->num_args==1 ? ": ":"s:")
#define NOT_A_TOPIC "There is no help topic number %d.\n"
#define MOREHELP0 \
"Type 'help' for a list of commands and topics.\n"
#define MOREHELP1 \
"Type 'help <topic-number>' for more information about a particular topic.\n"
#define MOREHELP2 \
"Type 'help <command-name>' for more help with a particular command.\n"
#define HELP_LEFT       26
#define HELP_FORMAT     "%-26s.%-52s\n"
#define HELP_DIVIDER \
"=============================================================================\
\n"

command help()
{
    int j, i, k, n;
    char *name=get_temp_string(), *rest=NULL, *str=get_temp_string(); 
    bool got_any, got_it;
	
    crunch(args);
    maybe_clear_screen();

    hold (more_mode) {

	if (nullstr(args)) { /**** Print Help Table ****/
	    print(HELP_DIVIDER); 
	    sprintf(ps, "%s Commands and Options:\n", the_program); pr();
	    for (j=1; j<MAX_COM_TOPICS; j++) { /* NOT j<=MAX_COM_TOPICS! */
		if (nullstr(topic_name[j]) || wizard_only(topic_code[j]))
		  continue;
		strcpy(name,topic_name[j]); uppercase(name);
		sprintf(ps, "\n(%d) %s:\n", j, name); pr();
		for (n=0; n<help_entries; n++) {
		    for (i=0, got_it=FALSE; i<command_num; i++) {
			if (cmd[i]->help_entry!=n) continue;
			if (cmd[i]->topic!=j) break;
			if (!allowed_cmd(i,TRUE)) break;
			got_it=TRUE; break;
		    }
		    if (!got_it) continue; /* for n */

		    if (!nullstr(cmd[i]->abbreviation)) 
		      sprintf(str, "%s (%s)", cmd[i]->name, cmd[i]->abbreviation);
		    else sprintf(str, "%s", cmd[i]->name);
		    if (!nullstr(cmd[i]->cmd_help)) {
			for (k=len(str); k<HELP_LEFT; k++) str[k]='.'; 
			str[HELP_LEFT]='\0';
			sprintf(ps, HELP_FORMAT, str, cmd[i]->cmd_help); pr();
		    } else { /* no cmd_help */
			print(str); nl();
		    }
		}
	    }
	    for (i=0, got_any=FALSE; i<command_num; i++) {
		if (!allowed_cmd(i,TRUE)) continue;
		if (cmd[i]->help_entry==UNHELPFUL) {
		    if (!got_any) 
		      { print("\nOTHER COMMANDS:\n"); got_any=TRUE; }
		    if (!nullstr(cmd[i]->abbreviation))
		      sprintf(ps, "%s (%s)\n", cmd[i]->name, cmd[i]->abbreviation);
		    else sprintf(ps, "%s\n", cmd[i]->name);
		    pr();
		}
	    }
	    nl();
	    print(MOREHELP1); print(MOREHELP2);
	    if (help_file==NULL) { print("\nNote: "); print(NO_HELP_FILE); }
	    print(HELP_DIVIDER); 
	    
	} else if (itoken(&args,iREQUIRED,&i)) { /**** a topic num ****/
	    if (i<=0 || i>MAX_COM_TOPICS || nullstr(topic_name[i])) {
		unhold(); sprintf(ps, NOT_A_TOPIC, i); pr();
		print(MOREHELP0); print(MOREHELP1); print(MOREHELP2);
	    } else {
		print(HELP_DIVIDER);
		strcpy(name,topic_name[i]); uppercase(name);
		sprintf(ps, "(%d) %s\n\n", i, name); pr();
		if (help_file==NULL) print(NO_HELP_FILE);
		else if (topic_help_key[i]==HELPLESS) print(NO_HELP_KEY);
		else {
		    fgoto_line(help_file,topic_help_key[i]);
		    fgetln(help_file);
		    got_any=FALSE;
		    while (ln[0]!='@') {
			got_any=TRUE;
			print(ln); nl();
			fgetln(help_file);
		    }
		    if (!got_any) print(NO_HELP_KEY);
		}
		print(HELP_DIVIDER);
	    }

	} else { /**** a command name? ****/
	    if ((i=parser(args,&rest,TRUE))<0) {
		unhold();
		print_parser_results(rest,TRUE); /* not a command */
		print(MOREHELP0); print(MOREHELP1); print(MOREHELP2);
	    } else {
		/* print out short help */
		print(HELP_DIVIDER);
		nstrcpy(name,cmd[i]->name,MAX_COM_NAME_LEN); uppercase(name);
		if (isa_option(cmd[i]->code))         strcpy(str,"Option");
		else if (isa_parameter(cmd[i]->code)) strcpy(str,"Parameter");
		else if (!help_only(cmd[i]->code))    strcpy(str,"Command");
		else                                 strcpy(str,"Information");
		if (!nullstr(cmd[i]->abbreviation))
		  sprintf(ps, "%s %s (abbreviation: '%s')\n", name, str,
                  cmd[i]->abbreviation);
		  else sprintf(ps, "%s %s\n", name, str);
		pr(); nl();

		if (!nullstr(cmd[i]->cmd_help)) 
		  { sprintf(ps, "Description: %s\n", cmd[i]->cmd_help); pr(); }
		if (cmd[i]->topic!=0) 
		  { sprintf(ps, "Help Topic:  (%d) %s\n", cmd[i]->topic,
                    topic_name[cmd[i]->topic]); pr(); }
		if (cmd[i]->num_args==0) print("No Arguments\n");
		else if (!nullstr(cmd[i]->args_help)) {
		    sprintf(ps, "Argument%s   %s\n", s_or_space_colon,
                    cmd[i]->args_help); pr();
		    if (!nullstr(cmd[i]->def_help)) {
			sprintf(ps, "Default%s    %s\n", s_or_space_colon,
                    cmd[i]->def_help); pr();
		    }
		}

		nl();
		if (help_file==NULL) print(NO_HELP_FILE);
		else if (cmd[i]->help_key==HELPLESS) print(NO_HELP_KEY);
		else {
		    fgoto_line(help_file,cmd[i]->help_key);
		    fgetln(help_file);
		    got_any=FALSE;
		    while (ln[0]!='@') {
			got_any=TRUE;
			print(ln); nl();
			fgetln(help_file);
		    }
		    if (!got_any) print(NO_HELP_KEY);
		}
		print(HELP_DIVIDER);
	    }
	}
    } /* hold */
}


command about()
{
    int n;
	
    print(HELP_DIVIDER);
    sprintf(ps, "%s version %s, Copyright %s Whitehead Institute\n",
            the_program, the_version, the_copyright); pr();
    if (gnu_copyright(ps)) { pr(); nl(); }
    print("All rights reserved.\n"); nl();

    n=cmd[com_num]->help_key;
    if (help_file==NULL)
      { print(NO_HELP_FILE); nl(); print(SURROGATE_ABOUT); }
    else if (n==HELPLESS)
      { print(NO_HELP_KEY); nl(); print(SURROGATE_ABOUT); }
    else {
	fgoto_line(help_file,n);
	fgetln(help_file);
	while (ln[0]!='@') {
	    print(ln); nl();
	    fgetln(help_file);
	}
    }
    nl(); print(MOREHELP0); print(MOREHELP1); print(MOREHELP2);
    print(HELP_DIVIDER);
}


command review_output() { review_memory(); print("ok\n"); }

command show_time() 
{ 
    char *str; 

    if (!nullstr(str=time_string())) { 
	print(str); nl();
    } else print("Time to go home...\n"); 
}


command cd_command()
{
    char dir_name[PATH_LENGTH+1];
    
    use_uncrunched_args();
    get_one_arg(stoken,"",dir_name);
    
    if (nullstr(dir_name)) {
	if (get_directory(dir_name)) 
	  { sprintf(ps, "The current directory is '%s'\n", dir_name); pr(); }
	else error("Can't get current directory name");

    } else { /* !nullstr(dir_name) */
	if (change_directory(dir_name)) {
	    if (get_directory(dir_name)) {
		sprintf(ps, "The current directory is now '%s'\n", dir_name); pr();
	    } else print("ok\n");
	} else { 
	    sprintf(ps, "Can't change to directory '%s'", dir_name);
	    error(ps); 
	}
    }
}


command system_command()
{
    if (!nullstr(args)) { /* run shell command */
	if (shell_command(args)) {
	    if (curses) 
	      print("\n        ...System Output Omitted...\n");
	    else if (logging) 
	      lib_puts(photo,"\n        ...System Output Omitted...\n");
	    print("\nOK\n");
	} else {
	    sprintf(ps, "Unable to run command '%s'", truncstr(args, 40));
	    error(ps);
	}

    } else { /* run a subshell */
	print(SHELL_MESSAGE); nl(); flush(); /* for curses */
	if (subshell()) {
	    if (curses) 
	      print("        ...System Interaction Omitted...\n");
	    else if (logging) 
	      lib_puts(photo,"        ...System Interaction Omitted...\n");
	    lib_puts(out,"\n");  /* C shell does not print one on ctrl-D */
	    print("\nBack in "); print(the_program); print("...\n");
	} else error("Unable to run subshell\n");
    }
}


#define GIMME_A_COMMENT \
"Enter your comment. End it with a period ('.') on a line by itself.\n"

command comment()
{
    if (nullstr(args)) {
	print(GIMME_A_COMMENT);
	do { getln(""); despace(ln); } while(!streq(ln,"."));
    } else print("ok\n"); 
}




/********** WIMP (Windows, Icons, Mouse, and Pointers) Support **********/

#ifdef NOT_PUBLIC
void wimp_start() /* Call from main() */
{
#ifdef HAVE_WIMP
    if (wimp) do_wimp_start(menus,num_menus);
#endif
}


void mkwimp(name,menu_entry,menu_num,wimp_function,status_function,
	  wimp_help_function,shortcut)
char *name, *menu_entry;
int menu_num;
void (*wimp_function)(), (*status_function)(), (*wimp_help_function)();
char shortcut;
{
    int i;

    if ((i=parser(name,(char**)NULL,TRUE))<0) send(CRASH);

    cmd[i]->menu_entry= mkstrcpy(menu_entry);
    cmd[i]->wimp_menu_num= menu_num;
    cmd[i]->wimp_procedure= wimp_function;
    cmd[i]->status_function= status_function;
    cmd[i]->wimp_help= wimp_help_function;
    cmd[i]->wimp_shortcut= shortcut;

    if (menu_num>=MAX_NUM_MENUS || nullstr(menu[menu_num]->title) ||
	menu[menu_num]->num_entries>=MAX_MENU_ENTRIES) send(CRASH);
    menu[menu_num]->entry[(menu[menu_num]->num_entries)++]= cmd[i];
}


void mkmenu(num,title)
int num;
char *title;
{
    if (num>=MAX_NUM_MENUS) send(CRASH);
    menu[num]->title=mkstrcpy(title);
}


COMMAND the_divider; 

void mkdivider(menu_num)
int menu_num;
{
    if (menu_num>=MAX_NUM_MENUS || nullstr(menu[menu_num]->title) ||
	menu[menu_num]->num_entries>=MAX_MENU_ENTRIES) send(CRASH);
    menu[menu_num]->entry[(menu[menu_num]->num_entries)++]= &the_divider;
}    
#endif
