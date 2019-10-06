#ifndef _EQN_H_
#define _EQN_H_

/*****************************************************************************

 ######   ####   #    #          #    #
 #       #    #  ##   #          #    #
 #####   #    #  # #  #          ######
 #       #  # #  #  # #   ###    #    #
 #       #   #   #   ##   ###    #    #
 ######   ### #  #    #   ###    #    #

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

/* Eqn is broken up into 2 major parts.  Make_equation takes an 
inputted equation and converts it to an array of structures in reverse
polish notation, and evaluate_equation evlautes the equation over a set
of points ***/


/* This is the data structure that I propose to use.  It is basically 
made up of two parts - is_a which is a flag that can have three 
possible values, SYMBOL, VARIABLE, or NUMBER, which indicate what the 
union is storing. I will use an array of pointers to these structures
of MAX_EQN_SIZE 500.  When sending this array to the function 
make_equation, a ** pointer to the array of pointers is used.*/
                    
#define MAX_EQN_SIZE 500 /* this setting means the largest equation we
			    can handle is about 100 characters      */

typedef struct {                  
    int is_a;  
    union {         
       	int variable;
       	int symbol;
      	real number;
    } val;
} EQUATION;

/***
 *
 * @param original_equation
 * @param variable_lookup
 * @return
 *
*************** The make equation subroutine **********************
* args:  char *original_eqn, char **symbol_table, int table_len
  inputted_eqn is a char pointer to the equation as inputted by the
  user, **symbol_table is a pointer to an array of chars that contains
  all variables in the expression, and table_len is the length of the
  symbol table

 The subroutine  **make_equation takes the inputted equation, and
   converts it to an array of structures.  It returns a pointer to the
   array of pointers.

  Error handling - if an error occurs in make_equation, a string is
    set that will tell what the error is will be sent to the routine send,
    which is a back-checking routine.  If the occurs, I will also free all
    arrays and pointers that have been malloced, and reinitialize all
    internal variables
*/
EQUATION **make_equation(char *original_equation, int (*variable_lookup)());
void parse_equation(char *original_equation, EQUATION **parsed_eqn, int (*variable_lookup)(), int *new_size, int *the_index);
void add_number(int mark, EQUATION **parsed_eqn, char *parsed_token, int *new_size, int *the_index);
void add_to_parsed_eqn(int mark, EQUATION **parsed_eqn, int parsed_token, int *new_size, int *the_index);
void add_parenthesis(int *i, int par, int mark, EQUATION **parsed_eqn);
void check_sizeof_array(int *size);
void postfix(EQUATION **parsed_eqn, EQUATION **postfixed);

/***
 *
 * @param postfixed
 * @param value_find
 * @return
 ******************* evaluating the equation ************************

 args: ELEMENT **parsed_equation, real *variable_values
  parsed_equation is the pointer to the array of pointers to the
  structures, and variable_value contains the value of each
  variable in the equation.

 Error handling - if a math error occurs, this will send the message
   MATHERR, if any other type of error occurs, then it probably a coding
   error that will have to be looked at, so the message CRASH will be sent
*/
real evaluate_equation (EQUATION **postfixed, real (*value_find)());

real pop_stack(void) /* This function returns the top value from the stack */;
real push_stack(real value_to_push);
void eqn_init(void);   /* Takes no arguments, it mallocs pasred and temp_eqn */
int variable_lookup(char *item);
real value_lookup(int index);

#define SYMBOL -1           /* these are flags          */
#define VARIABLE 2          /* to be used in the        */
#define NUMBER 1            /* data structures EQUATION */

#define PLUS 1
#define MINUS 2
#define MULTIPLY 3
#define DIVIDE 4
#define EXP 5
#define LOG 6
#define LN 7
#define SIN 8
#define ASIN 9
#define COS 10
#define ACOS 11
#define TAN 12
#define ATAN 13
#define LEFT_P 14
#define RIGHT_P 15

/* Various important global variables */

extern int stack_pointer;
extern real val[];
extern int eqnlen;
extern int table_size;
extern char **variable_table;
extern real *value_table;

#endif
