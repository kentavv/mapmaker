#ifndef _SYSCODE_H_
#define _SYSCODE_H_

real usertime(bool do_reset); /* return time in seconds, or -1.0 if fail */
char *time_string(void);    /* return ptr to "" if fail */
bool shell_command(char *cmd);

bool subshell(void);

bool change_directory(char *dir);

bool get_directory(char *buf);

bool get_home_directory(char *buf);

bool get_code_directory(char *buf);

bool rename_file(char *original_name, char *new_name);

bool fgoto_line(FILE *fp, long index);

long mkseed(long x);

void do_seedrand(long x);

real randnum(void);

void untrapped_msg(void); /* DO NOT ASSUME THAT MSGNAME IS SET! */
void trapped_msg(void); /* DO NOT ASSUME THAT MSGNAME IS SET! */
void verbose_untrapped_msg(void); /* DO NOT ASSUME THAT MSGNAME IS SET! */
void do_trap(void);

void sigcounter(void);

void signal_trap_init(void);

//void get_screen_preferences();
//void tty_hello();
bool do_gnu_readline(char *prompt, char *str, int num);

int set_rl_default(void);

bool do_gnu_edit(char *prompt, char *str, int num, char *initial /* initial may be = str */);

bool gnu_copyright(char *str /* side-effected, so it must be big enough */);

bool tty_gets(char *str, int num);

bool file_gets(FILE *fp /* must be opened with file_open() */, char *str  /* must be num+2 chars long, but use num+3 in case of weirdness */,
               int num    /* num chars, not including the '\n' or '\0', will be read */);

void lib_puts(FILE *fp, char *str);

//void iocheck(void);
void tty_init(void);

bool check_tty_lines(void); /* return TRUE and set tty_lines if changed */
void ansi_cursor_left(int i, char *s);

void hp_cursor_left(int i, char *s);

//void mac_cursor_left(void);
void tty_hello(void);

bool do_clear_screen(void);

bool do_delete_previous_line(void); /* Needed for the "Hit RETURN for more" thing */
bool do_highlight(bool reverse);

bool do_cursor_left(int num_spaces /* might be FAR_LEFT */, char *str_to_print);

//bool boing(void);
void misc_init(void); /* init this file */
void custom_lib_init(void);

void lib_init(void);

//void lib_inits(int *argc_ptr, char *argv[]);
void get_cmd_line_args(int *argc_ptr, char *argv[]);

bool check_file_arg(int num, char *arg, char *name, char *type, char *def_ext, char *prog, char *mode);

bool update_top(void);
//bool screen_init(int *argc_ptr, char *argv[]); /* side-effect wimp, curses, and split */
//bool split_screen_init(int *argc_ptr, char *argv[], int top_lines, void (*top_update_function)());
//void get_screen_preferences(int *argc_ptr, char *argv[], int *try_curses, int *try_wimp);
//bool update_top(void);
//void screen_end(void);
//bool string_editor(char *prompt, char *str, int num /* max #chars in str */);
//void curses_error(char *where);

#endif