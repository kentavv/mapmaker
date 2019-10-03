/******************************************************************************

  #####    ##    #####   #       ######          #    #
    #     #  #   #    #  #       #               #    #
    #    #    #  #####   #       #####           ######
    #    ######  #    #  #       #        ###    #    #
    #    #    #  #    #  #       #        ###    #    #
    #    #    #  #####   ######  ######   ###    #    #

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

/* Tables are a generally useful structure, allowing one to make a list of
   strings, where each is assigned either a numeric or alphanumeric key. 
   Tables can be made to either expand themselves as needed, or to bash old 
   entries. Numbered tables are kept sorted, and named tables may be sorted. */

typedef struct named_entry {
    char *string;
    union { int num; char *name; } id;
    struct named_entry *next; /* make a linked list */
} TABLE_ENTRY;
    
typedef struct { 
    struct named_entry *list, *unused;
    bool named_entries;
    int string_length, expands_by;
    int next_entry_num; /* used if !named_entries */
} TABLE;

TABLE *allocate_table();
/* args: int initial_num_entries, string_length, expands_by, index_by_name; */
#define EXPANDS_BY(num_entries) (num_entries)
#define CANT_EXPAND 0
#define INDEX_BY_NUMBER FALSE
#define INDEX_BY_NAME   TRUE

void free_table();       /* args: TABLE *list; */

void put_named_entry();    /* args: char *name, *string; TABLE *p; */
void put_numbered_entry(); /* args: char *string; TABLE *p; int *num; */
/* string is saved in table under the index name or num. Names must be
single tokens both despace()ed and filter()ed. Numbered strings are
assigned numbers sequentially, with numbers starting at 0, indicated by
*num. This routine does not check for duplicate names or some other
bad things. If the table is out of space for new entries, either more 
space is added, or the oldest entry is bashed, depending on the setting
of expands_by when the list was allocated. */

bool get_named_entry();    /* args: char *name; char **string, **full_name; 
			            TABLE *p; flag *fail;*/
bool delete_named_entry(); /* args: char *name; TABLE *p; flag *fail; */
/* The string with the specified name is looked up the table using matches()
to find an unambiguous match for the (possibly abbreviated) name. If no
such string can be found, FALSE is returned and *fail is set to one of the
constants below. Otherwise, delete_named_entry() deletes the string from the 
table, while get_named_entry() side-effects *string and *full_name 
accordingly. */
#define NAME_DOESNT_MATCH 1
#define NAME_IS_AMBIGUOUS 2

bool get_numbered_entry();    /* args: int num; char **string; TABLE *p; */
bool delete_numbered_entry(); /* args: int num; TABLE *p; UNIMPLEMENTED! */
/* The string with the specified number is looked up the table. If no
such string can be found, FALSE is returned. Otherwise,
delete_numbered_entry() deletes the string from the table, while
get_numbered_entry() side-effects *string to point to it. */

bool table_full();          /* args: TABLE *p; */
  /* TRUE iff the list is full AND it can't be expanded. */
bool table_empty();         /* args: TABLE *p; */
int  count_table_entries(); /* args: TABLE *p; */
#define next_entry_number(table) ((table)->next_entry_num)

bool valid_name(); /* args: char *name; */

/* iterator macros:

   for_all_numbered_entries()  args: TABLE *p; char *string; int num; 
   for_all_named_entries()     args: TABLE *p; char *string, *name;

Iterate over all entries in the table, setting num to the number or
name to the name, and string to the string (note that these are NOT
pointers to the side-effected variables, as is necessary in C
functions). The loop always goes in a FIFO (first-in, first-out)
fashion. Do not put one for_all...() loop inside another! */

#define for_named_entries(Table,Str,Name) 				\
 for (Te=Table->list, Name=Te->id.name, Str=Te->string; 		\
      Te!=NULL; 							\
      Te=Te->next, Name=Te->id.name, Str=Te->string)

#define for_numbered_entries(Table,Str,Num) 				\
 for (Te=Table->list, Num=Te->id.num, Str=Te->string; 			\
      Te!=NULL; 							\
      Te=Te->next, Num=Te->id.num, Str=Te->string)

extern TABLE_ENTRY *Te;

void save_table(), load_table();

extern TABLE *cmd_history;

#define NAME_TAG_CHARS    "*"
#define NUMBER_TAG_CHARS  "#"

#define NAME_FIRST_CHARS \
	"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
#define NAME_CHARS \
 	"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._"
