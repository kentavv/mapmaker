#ifndef _SHELL_H_
#define _SHELL_H_

/******************************************************************************

  ####   #    #  ######  #       #               #    #
 #       #    #  #       #       #               #    #
  ####   ######  #####   #       #               ######
      #  #    #  #       #       #        ###    #    #
 #    #  #    #  #       #       #        ###    #    #
  ####   #    #  ######  ######  ######   ###    #    #

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

extern bool wizard_mode;
//void null_command();
extern FILE *help_file;
typedef void command;

extern bool (*quit_save_hook)(bool);
//extern bool more_mode;



void banner(void);

char *get_version(char *version_filename);

void photo_banner(void);

char *centering(char *str, bool lineup);

void shell_init(char *program, char *version, char *copyright, char *help_filename);

char *default_prompt(char *s);

void mktopic(int num, char *nam, int code, long description_index);

int mkcommand(char *name, char *abbrev, void (*func)(void), int code);

void mkhelp(char *cmd_name, char *abbrev, long description_index, int num_args_prefix, int num_args, int code, int topic, char *description,
            char *arguments, char *defaults);

bool valid_name(char *str);

void null_command(void);

int parser(char *line, char **rest, bool help_ok);

void print_parser_results(char *rest, bool help_ok);

int try_to_match(char **token, int n_tokens, int n_to_try, int n_word_command, int *com_match, int *n_matched, int *last_match,
                 int allow_help_only_stuff);

int com_matches(char **in_tokens, char **com_tokens, int num_to_match, int *exact);

void expand_history_references(char *line);

void command_loop(void);

void abort_command(void);

void error(const char *errmsg);

void maybe_set_bool(bool *var);

void maybe_set_real(real *var, real lbound, real hbound, real fmt);

void maybe_set_long(long *var, long lbound, long hbound);

void maybe_set_int(int *var, int lbound, int hbound);

void set_usage_error(char *com_args);

void usage_error(int num);

void nomore_args(int n);

void more_args(int num);

void input_error(char *val, char *def);

void expect_nomore_input(char *str, int mode, int num);

void expect_more_input(char *str, int mode, int num);

bool split_arglist(char **rest, int divider);

bool split_uncrunched_args(char **rest, int divider);

void maybe_ok(char *str);

void keep_user_amused(char *thing, int iter, int max_iter);

command show_cmd_history(void);

command quit(void);

command really_quit(void);

command run_from_file(void);

command do_photo(void);

command set_more(void);

command set_wizard(void);

command set_verbose_mem(void);

command help(void);

command about(void);

command review_output(void);

command show_time(void);

command cd_command(void);

command system_command(void);

command comment(void);

void wimp_start(void);

void mkwimp(char *name, char *menu_entry, int menu_num, void (*wimp_function)(void), void (*status_function)(void), void (*wimp_help_function)(void),
            int shortcut);

void mkmenu(int num, char *title);

void mkdivider(int menu_num);



//int parser(); /* args: char *line, **rest; bool allow_help; */
///* Parser() returns a command table index, or -1 if it fails. *rest is
//   side-effected to point to the unparsed portion of the line, which, on
//   failure, is the entire line (**rest may be NULL). The globals matched[]
//   and num_matched are also set, with num_matched= the number of command
//   names matched (0-max), NOT INCLUDING ABBREVIATIONS, and
//   matched[i]=TRUE if command #i's name was matched. */
//
//void print_parser_results(); /* args: char *rest; also uses global state */

/* These globals are set any time a command is run... */
extern int com_num;  /* command # in the command table */
extern int num_args; /* set to 0 when the command is called */
extern char *com, *args;  /* the command name (full), and it's arguments. Note
  that the args pointer may be side-effected, eg by stoken...(it gets reset) */
extern char *uncrunched_args; /* OK to hack this pointer too */
extern int cmd_history_num;   /* in shell.c -- looks nice in prompt() */

/***** Here are many nice things for writing commands *****/

//void maybe_set_bool();  /* args: bool *var; */
//void maybe_set_real();  /* args: real *var, low_bound, high_bound, lf_format;*/
//void maybe_set_long();  /* args: long *var, low_bound, high_bound; */
//void maybe_set_int();   /* args: long *var, low_bound, high_bound; */
//
///* These all call abort_command(), which sends SOFTABORT */
//void error(const char *errmsg);        /* args: char *msg; prints "error: <msg>\n" and aborts */
//void usage_error();  /* args: int args; 1..n; -1 => unknown; use arg_num! */
//void more_args();    /* args: int args; is #of args present */
//void nomore_args();  /* args: int args; is #of args present */
//void set_usage_error(); /* args: char *msg */

#define ARGS_LEN       1900

#define get_arg(parse_func, default, val_ptr)        \
 if (parse_func(&args,default,val_ptr)) num_args++;    \
 else { if (nullstr(args)) usage_error(num_args);    \
        more_args(++num_args); }

#define get_one_arg(parse_func, default, val_ptr)    \
 { get_arg(parse_func,default,val_ptr); nomore_args(0); }

//bool split_arglist(); /* args: char **rest, divider;
//  side-effects (shortens) the global variable args, splitting the arg list
//  into two strings around a divider (which is removed) */
//bool split_uncrunched_args(); /* the same for uncrunched_args */

#define use_uncrunched_args() strcpy(args,uncrunched_args)

///* older and/or lower level stuff */
//void abort_command();       /* no args - sends SOFTABORT */
//void expect_more_args();    /* args: num_args_given */
//void input_error();         /* args: char *correct_input, *defaults; */
//void expect_nomore_input(); /* args: char *rest; int qualif, vals_wanted; */
//void expect_more_input();   /* args: char *rest; int qualifier, vals_wanted; */

/* these are the qualifiers: EXACTLY and UPTO apply to _nomore_, while EXACTLY
   and ATLEAST apply to _more_ */
#define EXACTLY 0
#define UPTO 1
#define ATLEAST 2

///* existing system commands */
//command quit();
//command really_quit();
//command run_from_file();
//command do_photo();
//command set_more();
//command set_verbose_memory();
//command help();
//command about();
//command show_cmd_history();
//command review_output();
//command show_time();
//command cd_command();
//command system_command();
//command comment();
//command set_wizard();

extern bool photo_update_top_hook;  /* make photo call update_top() */
extern void (*photo_banner_hook)(FILE *); /* make photo_banner call this w/arg fp */

//void maybe_ok(); /* args: char *str; prints the str then a '\n',
//   although only "ok\n" is printed to the screen if update_top() succeeds.
//   Use if the acknowledgement for a command might just be shown up top. */
//
//void keep_user_amused(); /* args: char *str; int iter_num, max_iter_num;
//   keeps user amused - use str==NULL or "" to end amusements. DO NOT
//   PRINT ANYTHING UNTIL THE AMUSEMENTS ARE OVER! */
//
///* This is how to access the shell... */
//void shell_init();   /* args: char *program, *version, *copyright, *help; */
//void banner();       /* no args */
//void command_loop(); /* no args; invokes the shell */
//char *get_version(); /* args: char *version_file_name; (in source dir) */
//
///* To set up the command table... */
//int  mkcommand(); /* args: char *name, *abbrev; void (*f)(); int topic, code;*/
//void mktopic();   /* args: int num (>0); char *name; int code; (see below) */
//void mkhelp();    /* args: char *cmd_name; long help_file_index;
//   int qualifier, num_args; char *line, *args, *defaults; */
//
///* Support for a WIMP interface... */
//void mkwimp(); /* args: char *name, *menu_entry; int menu_num;
//   void (*wimp_function)(), (*status_function)(), (*wimp_help_function)();
//   char shortcut; */
//void mkmenu();     /* args: int menu_num; char *menu_title; */
//void mkdivider();  /* args: int menu_num; */
//void wimp_start(); /* no args */
extern int inhibit_menus; /* FALSE when it's OK to execute a command */

#define mc mkcommand
#define mh mkhelp
#define mt mktopic
#define mm mkmenu
#define md mkdivider

/* mktopic() and mkcommand() codes (DO NOT CHANGE THESE). Note that
setting a topic to WIZard disables the printing of that topic and all
of its commands when help lists all available commands [for
non-wizards]. Individual commands must be WIZed to disable their use,
or to disable the printing of help information for them alone. */

#define CMD 0
#define OPT 1
#define PAR 2
#define HLP 3
#define TOP 4  /* TOP or WIZ for topics */
#define WIZ 10
#define wizard_only(code)   (code>=10)
#define isa_command(code)   (code%10==0)
#define isa_option(code)    (code%10==1)
#define isa_parameter(code) (code%10==2)
#define help_only(code)     (code%10==3)
#define allowed_cmd(num, allow_help)                \
  ((allow_help  || !help_only(cmd[num]->code)) &&  \
   (wizard_mode || !wizard_only(cmd[num]->code)))

#define MAX_COM_TOKENS     3
#define MAX_COM_NAME_LEN 25   /* might be more? 26 in help output */
#define MAX_COM_HELP_LEN 52
#define MAX_ARG_HELP_LEN 78
#define MAX_ABBREV_LEN   3
#define MAX_TOPIC_LEN     70
#define MAX_COMMANDS     200
#define MAX_COM_TOPICS     20

/* If you want a custom prompt, set prompt to point to the function... */
extern char *((*prompt)()); /* ptr to a func which returns a ptr to a string - 
   it gets as an arg char *str; and it should side-effect AND return str. */

typedef struct {
    int num_tokens;
    char *name;
    char abbreviation[MAX_ABBREV_LEN + 1];
    char **tokens;         /* [token#] -> string */
    void (*procedure)(void);     /* ptr to a void function of no arguments */
    char *cmd_help;
    char *args_help;
    char *def_help;
    int num_args;
    int num_args_prefix; /* ATLEAST, EXACTLY, or UPTO */
    long help_key;   /* index into the help file */
    int help_entry; /* entry number in the help file, 1...N */
    int code, topic;

    void (*wimp_procedure)(void);

    void (*status_function)(void);

    void (*wimp_help)(void);

    int wimp_menu_num;
    char *menu_entry;
    char wimp_shortcut;
} COMMAND;

typedef struct {
    char *title;
    int num_entries;
    COMMAND **entry;  /* [entry#] => ptr to an entry in the cmd struct */
} MENU;

extern COMMAND the_divider; /* dummy thing - used to put dividers in menus */

#endif
