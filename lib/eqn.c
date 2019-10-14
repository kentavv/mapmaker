/******************************************************************************

 ######   ####   #    #           ####                                  
 #       #    #  ##   #          #    #                                        
 #####   #    #  # #  #          #                                     
 #       #  # #  #  # #   ###    #                                             
 #       #   #   #   ##   ###    #    #                                        
 #######  ### #  #    #   ###     ####                 

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#include "system.h"

#define VALUE (-1)
#define OP 0
#define LEFT_PAREN 1
#define MISSING_PHENO (-100000.0)

static int stack_pointer;

#define SIZEOF_STACK 100
real val[SIZEOF_STACK];
static int eqnlen;
int table_size;
char **variable_table;
real *value_table;


#define SYMBOL (-1)           /* these are flags          */
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

static void parse_equation(char *original_equation, EQUATION **parsed_eqn, int (*variable_lookup)(char *), int *new_size, int *the_index);

static void add_number(int mark, EQUATION **parsed_eqn, char *parsed_token, int *new_size, int *the_index);

static void add_to_parsed_eqn(int mark, EQUATION **parsed_eqn, int parsed_token, int *new_size, int *the_index);

static void add_parenthesis(const int *i, int par, int mark, EQUATION **parsed_eqn);

static void check_sizeof_array(const int *size);

static void postfix(EQUATION **parsed_eqn, EQUATION **postfixed);

static real pop_stack(void) /* This function returns the top value from the stack */;

static real push_stack(real value_to_push);


EQUATION **make_equation(char *original_equation, int (*variable_lookup)(char *)) {
    EQUATION **parsed_eqn, **postfixed = NULL;
    int i;
    int the_index, new_size;
    /* I'll now call the function parse_equation, that will take
       the inputted equation, parse it into individual tokens, and
       then place the tokens into the array of structures.      */

    run {
            the_index = 0;
            new_size = 0;
            array(parsed_eqn, MAX_EQN_SIZE, EQUATION*);
            for (i = 0; i < MAX_EQN_SIZE; i++) {
                single(parsed_eqn[i], EQUATION);
            }
            parse_equation(original_equation, parsed_eqn,
                           variable_lookup, &new_size, &the_index);

            /* Now that the equation has been parsed, I will parray, then
               call postfix, which will take the parsed_eqn, and place it
               in reverse polish notation, for easy evaluation */

            array(postfixed, MAX_EQN_SIZE, EQUATION*);
            for (i = 0; i < new_size + 1; i++) single(postfixed[i], EQUATION);

            postfix(parsed_eqn, postfixed);
        } on_exit {
        relay_messages;
    }

    return postfixed;
}


/******************** The parse_eqn algorithim *************************/
/* This algorithim uses a state table. The state variable STATE has three
   expected values: VALUE which means a number or varaible is expected, OP
   which means some operation is expected, and LEFT_PAREN which means a left
   parenthesis is expected.  A complete state table follows:
   
   Given State | input          | output              | new state
   ------------------------------------------------------------------
   VALUE       | *,/,),^        | C R A S H ! ! ! - Error message
   ------------------------------------------------------------------
   VALUE       | number,variable| same                | OP
   ------------------------------------------------------------------
   VALUE       | +,-            | U+, U-              | VALUE
   ------------------------------------------------------------------
   VALUE       | (              | ((((                | VALUE
   ------------------------------------------------------------------
   OP          | num,var,ufunct.| C R A S H ! ! ! - Error message
   ------------------------------------------------------------------
   OP          | +,-            | )))+(((,)))-(((     | VALUE
   ------------------------------------------------------------------
   OP          | *,/            | ))*((,))/((         | VALUE
   ------------------------------------------------------------------
   OP          | ^              | )^(                 | VALUE
   ------------------------------------------------------------------
   OP          | )              | ))))                | OP
   ------------------------------------------------------------------
   LEFT_PAREN  | (              | ((((                | VALUE
   ------------------------------------------------------------------
   LEFT_PAREN  | anything else  | C R A S H ! ! ! - Error message
   ------------------------------------------------------------------
   
   It can be plainly seen that this algorithim adds many parenthesis
   to the equation.  The point of this will be seen in the second algorithim,
   postfix equation.  This state table will also make sure the equation is in
   proper form, i.e. nothing like 7 +* 2 or (^4 + 5).    ******************/


void parse_equation(char *original_equation, EQUATION **parsed_eqn, int (*variable_lookup)(char *), int *new_size, int *the_index) {
    int i, mark, state;
    char parsed_token[TOKLEN + 1];
    real tok;
    char *token;

    token = get_temp_string();
    self_delimiting = "()*+-/^";
    state = VALUE; /* Sets initial state in state machine to expect a VALUE */

    eqnlen = len(original_equation);
    while (stoken(&original_equation, sREQUIRED, parsed_token)) {
        /* stoken parses the equation into individual tokens */
        if (state == VALUE) {
            strcpy(token, parsed_token);
            if (streq(parsed_token, "log")) {
                mark = SYMBOL;
                add_to_parsed_eqn(mark, parsed_eqn, LOG, new_size, the_index);
                state = LEFT_PAREN;
            } else if (streq(parsed_token, "sin")) {
                mark = SYMBOL;
                add_to_parsed_eqn(mark, parsed_eqn, SIN, new_size, the_index);
                state = LEFT_PAREN;
            } else if (streq(parsed_token, "ln")) {
                mark = SYMBOL;
                add_to_parsed_eqn(mark, parsed_eqn, LN, new_size, the_index);
                state = LEFT_PAREN;
            } else if (streq(parsed_token, "cos")) {
                mark = SYMBOL;
                add_to_parsed_eqn(mark, parsed_eqn, COS, new_size, the_index);
                state = LEFT_PAREN;
            } else if (streq(parsed_token, "tan")) {
                mark = SYMBOL;
                add_to_parsed_eqn(mark, parsed_eqn, TAN, new_size, the_index);
                state = LEFT_PAREN;
            } else if (streq(parsed_token, "atan")) {
                mark = SYMBOL;
                add_to_parsed_eqn(mark, parsed_eqn, ATAN, new_size, the_index);
                state = LEFT_PAREN;
            } else if (streq(parsed_token, "asin")) {
                mark = SYMBOL;
                add_to_parsed_eqn(mark, parsed_eqn, ASIN, new_size, the_index);
                state = LEFT_PAREN;
            } else if (streq(parsed_token, "acos")) {
                mark = SYMBOL;
                add_to_parsed_eqn(mark, parsed_eqn, ACOS, new_size, the_index);
                state = LEFT_PAREN;
            } else if (streq(parsed_token, "(")) {
                mark = SYMBOL;
                add_to_parsed_eqn(mark, parsed_eqn, LEFT_P, new_size, the_index);
                state = VALUE;
            } else if (rtoken(&token, rREQUIRED, &tok)) {
                mark = NUMBER;
                sprintf(ps, "%lf", tok);
                add_number(mark, parsed_eqn, ps, new_size, the_index);
                state = OP;
            } else if (('a' <= *parsed_token && *parsed_token <= 'z') ||
                       ('0' <= *parsed_token && *parsed_token <= '9')) {
                /* need to execute the lookup function  ******/
                i = (*variable_lookup)(parsed_token);
                mark = VARIABLE;
                add_to_parsed_eqn(mark, parsed_eqn, i, new_size, the_index);
                state = OP;
            } else {
                BADEQN_errpos = eqnlen - len(original_equation) - 1;
                BADEQN_errmsg = "Equation error - VALUE expected";
                send (BADEQN);
            }
        } else if (state == OP) {
            if (streq(parsed_token, "*")) {
                mark = SYMBOL;
                add_to_parsed_eqn(mark, parsed_eqn, MULTIPLY, new_size,
                                  the_index);
                state = VALUE;
            } else if (streq(parsed_token, "+")) {
                mark = SYMBOL;
                add_to_parsed_eqn(mark, parsed_eqn, PLUS, new_size, the_index);
                state = VALUE;
            } else if (streq(parsed_token, "^")) {
                mark = SYMBOL;
                add_to_parsed_eqn(mark, parsed_eqn, EXP, new_size, the_index);
                state = VALUE;
            } else if (streq(parsed_token, "-")) {
                mark = SYMBOL;
                add_to_parsed_eqn(mark, parsed_eqn, MINUS, new_size, the_index);
                state = VALUE;
            } else if (streq(parsed_token, "/")) {
                mark = SYMBOL;
                add_to_parsed_eqn(mark, parsed_eqn, DIVIDE, new_size, the_index);
                state = VALUE;
            } else if (streq(parsed_token, ")")) {
                mark = SYMBOL;
                add_to_parsed_eqn(mark, parsed_eqn, RIGHT_P, new_size, the_index);
                state = OP;
            } else {
                BADEQN_errpos = eqnlen - len(original_equation) - 1;
                BADEQN_errmsg = "Equation error-operation expected";
                send (BADEQN);
            }
        } else if (state == LEFT_PAREN) {
            if (streq(parsed_token, "(")) {
                mark = SYMBOL;
                add_to_parsed_eqn(mark, parsed_eqn, LEFT_P, new_size, the_index);
                state = VALUE;
            } else {
                BADEQN_errpos = eqnlen - len(original_equation) - 1;
                BADEQN_errmsg = "Error - left parenthesis expected";
                send (BADEQN);
            }
        }
    }
}


void add_number(int mark, EQUATION **parsed_eqn, char *parsed_token, int *new_size, int *the_index) {
    real new_number; /* is used to convert tokens into reals */

    if (mark == NUMBER) /* only need to add the number */ {
        parsed_eqn[*the_index]->is_a = mark; /* sets is_a */
        rtoken(&parsed_token, rREQUIRED, &new_number); /* converts the */
        parsed_eqn[*the_index]->val.number = new_number; /* token to real */
        (*the_index)++; /* increment array */
        check_sizeof_array(the_index);
        (*new_size)++;
    }
}


void add_to_parsed_eqn(int mark, EQUATION **parsed_eqn, int parsed_token, int *new_size, int *the_index) {
    /* This is is function that adds the structure to the array */
    int t;

    if (mark == VARIABLE) {
        parsed_eqn[*the_index]->is_a = mark;
        parsed_eqn[*the_index]->val.variable = parsed_token;
        *the_index = (*the_index + 1);
        check_sizeof_array(the_index);
        *new_size = (*new_size + 1);
    } else if (mark == SYMBOL) {
        if (parsed_token == LEFT_P || parsed_token == RIGHT_P) {
            for (t = 0; t < 4; t++, *the_index = (*the_index + 1)) {
                add_parenthesis(the_index, parsed_token, mark, parsed_eqn);
            }
            /* where I add so many parenthesis, I **/
            /* wrote a subroutine to do it for me **/
            check_sizeof_array(the_index);
        } else if (parsed_token == PLUS || parsed_token == MINUS) {
            /* Need to add )))+/-((( */
            for (t = 0; t < 3; t++, *the_index = (*the_index + 1)) {
                add_parenthesis(the_index, RIGHT_P, mark, parsed_eqn);
            }
            parsed_eqn[*the_index]->is_a = mark;
            parsed_eqn[*the_index]->val.symbol = parsed_token;
            *the_index = (*the_index + 1);
            for (t = 0; t < 3; t++, *the_index = (*the_index + 1)) {
                add_parenthesis(the_index, LEFT_P, mark, parsed_eqn);
            }
            check_sizeof_array(the_index);
            *new_size = (*new_size + 1);
        } else if (parsed_token == MULTIPLY || parsed_token == DIVIDE) {
            /* Need to add ))* or /(( */
            for (t = 0; t < 2; t++, *the_index = (*the_index + 1)) {
                add_parenthesis(the_index, RIGHT_P, mark, parsed_eqn);
            }
            parsed_eqn[*the_index]->is_a = mark;
            parsed_eqn[*the_index]->val.symbol = parsed_token;
            *the_index = (*the_index + 1);
            for (t = 0; t < 2; t++, *the_index = (*the_index + 1)) {
                add_parenthesis(the_index, LEFT_P, mark, parsed_eqn);
            }
            check_sizeof_array(the_index);
            *new_size = (*new_size + 1);
        } else if (parsed_token == EXP) {
            /* Need to add )^( */
            add_parenthesis(the_index, RIGHT_P, mark, parsed_eqn);
            *the_index = (*the_index + 1);
            parsed_eqn[*the_index]->is_a = mark;
            parsed_eqn[*the_index]->val.symbol = parsed_token;
            *the_index = (*the_index + 1);
            add_parenthesis(the_index, LEFT_P, mark, parsed_eqn);
            *the_index = (*the_index + 1);
            check_sizeof_array(the_index);
            *new_size = (*new_size + 1);
        } else if (parsed_token == LOG || parsed_token == SIN ||
                   parsed_token == LN || parsed_token == COS ||
                   parsed_token == TAN || parsed_token == ATAN ||
                   parsed_token == ASIN || parsed_token == ACOS) {
            parsed_eqn[*the_index]->is_a = mark;
            parsed_eqn[*the_index]->val.symbol = parsed_token;
            *the_index = (*the_index + 1);
            check_sizeof_array(the_index);
            *new_size = (*new_size + 1);
        }
    }
}


void add_parenthesis(const int *i, int par, int mark, EQUATION **parsed_eqn) {
    parsed_eqn[*i]->is_a = mark;
    parsed_eqn[*i]->val.symbol = par;
}


void check_sizeof_array(const int *size) {
    if (*size > MAX_EQN_SIZE) {
        BADEQN_errpos = -1;
        BADEQN_errmsg = "Error - Equation is too long!";
        send (BADEQN);
    }
}


/******************  The postfix algorithm  ************************/
/**   The principle behind postfix is to take the array of data 
  structures that was created in parse equation and to put it in reverse
  polish notation.  As was mentioned before, this is done using stacks.
  The mark of each structure will be read, and then the char will be 
  accessed.  If it is anything but a right parentheses, the structure
  will be pushed onto the stack.  If it is a right parentheses, then 
  structures will be popped off the stack and added to the new output 
  array until the end or a left parenthesis is reached.  No parentheses
  are added to the output array.  At the end, all other structures left 
  on the stack are popped and placed on the output array.            **/

void postfix(EQUATION **parsed_eqn, EQUATION **postfixed) {
    int i, parsed_index = 0, temp_index = 0, post_index = 0;
    EQUATION **temp_eqn;

    array(temp_eqn, MAX_EQN_SIZE, EQUATION*);
    for (i = 0; i < MAX_EQN_SIZE; i++) {
        single(temp_eqn[i], EQUATION);
    }
    while (parsed_eqn[parsed_index]->is_a != 0) {
        if (parsed_eqn[parsed_index]->is_a == NUMBER ||
            parsed_eqn[parsed_index]->is_a == VARIABLE) {
            temp_eqn[temp_index] = parsed_eqn[parsed_index];
            temp_index++;
        } else {
            if (parsed_eqn[parsed_index]->val.symbol != RIGHT_P &&
                parsed_eqn[parsed_index]->is_a == SYMBOL) {
                temp_eqn[temp_index] = parsed_eqn[parsed_index];
                temp_index++;
            } else {
                if (temp_index >= 0) temp_index--;
                while (1) {
                    if (temp_index < 0) break;
                    if (temp_eqn[temp_index]->val.symbol == LEFT_P &&
                        temp_eqn[temp_index]->is_a == SYMBOL) {
                        break;
                    } else {
                        postfixed[post_index] = temp_eqn[temp_index];
                        post_index++;
                        temp_index--;
                    }
                }
            }
        }
        if (temp_index < 0) temp_index = 0;
        parsed_index++;
    }
    temp_index--;
    while (temp_index >= 0) {
        if (temp_eqn[temp_index]->val.symbol == LEFT_P &&
            temp_eqn[temp_index]->is_a == SYMBOL)
            temp_index--;
        else {
            postfixed[post_index] = temp_eqn[temp_index];
            post_index++;
            temp_index--;
        }
    }
    postfixed[post_index]->is_a = 0;
}


/********************* EVALUATING THE EQUATION **********************/

real evaluate_equation(EQUATION **postfixed, real (*value_find)(int)) {
    int i, missing = FALSE;
    real number_to_use, c, exponent, divisor, var_push, subtract;
    i = 0;
    stack_pointer = 0;

    run {
            while (postfixed[i]->is_a != 0) {
                if (postfixed[i]->is_a == NUMBER) {
                    push_stack(postfixed[i]->val.number);
                } else if (postfixed[i]->is_a == VARIABLE) {
                    var_push = value_lookup(postfixed[i]->val.variable);
                    if (var_push == MISSING_PHENO) {
                        missing = TRUE;
                        break;
                    }
                    push_stack(var_push);
                } else {
                    switch (postfixed[i]->val.symbol) {
                        case PLUS:
                            push_stack(pop_stack() + pop_stack());
                            break;
                        case MINUS:
                            subtract = pop_stack();
                            push_stack(pop_stack() - subtract);
                            break;
                        case MULTIPLY:
                            push_stack(pop_stack() * pop_stack());
                            break;
                        case DIVIDE:
                            divisor = pop_stack();
                            if (divisor != 0.0) {
                                push_stack(pop_stack() / divisor);
                            } else {
                                BADEQN_errpos = -1;
                                BADEQN_errmsg = "Error - division by 0 attempted\n";
                                send(BADEQN);
                            }
                            break;
                        case EXP:
                            exponent = pop_stack();
                            push_stack(pow(pop_stack(), exponent));
                            break;
                        case LOG:
                            number_to_use = pop_stack();
                            if (number_to_use > 0.0) {
                                push_stack(log10(number_to_use));
                            } else {
                                BADEQN_errpos = -1;
                                BADEQN_errmsg = "Error - can't take log of 0 or negative number\n";
                                send(BADEQN);
                            }
                            break;
                        case LN:
                            number_to_use = pop_stack();
                            if (number_to_use > 0.0)
                                push_stack(log(number_to_use));
                            else {
                                BADEQN_errpos = -1;
                                BADEQN_errmsg = "Error - can't take ln of 0 or negative number\n";
                                send(BADEQN);
                            }
                            break;
                        case SIN:
                            push_stack(sin(pop_stack()));
                            break;
                        case ASIN:
                            number_to_use = pop_stack();
                            if ((number_to_use < 1.0 && number_to_use > -1.0)) {
                                push_stack(asin(number_to_use));
                            } else {
                                BADEQN_errpos = -1;
                                BADEQN_errmsg = "Error - can't take asin of number greater than 1\n";
                                send(BADEQN);
                            }
                            break;
                        case COS:
                            push_stack(cos(pop_stack()));
                            break;
                        case ACOS:
                            number_to_use = pop_stack();
                            if ((number_to_use < 1.0 && number_to_use > -1.0)) {
                                push_stack(acos(number_to_use));
                            } else {
                                BADEQN_errpos = -1;
                                BADEQN_errmsg = "Error - can't take acos of number greater than 1\n";
                                send(BADEQN);
                            }
                            break;
                        case TAN:
                            push_stack(tan(pop_stack()));
                            break;
                        case ATAN:
                            push_stack(atan(pop_stack()));
                            break;
                    }
                }
                i++;
            }
        } on_exit {
        relay_messages;
    }

    if (missing == TRUE) {
        return (MISSING_PHENO);
    } else {
        c = pop_stack();
        return (c);
    }
}


real pop_stack(void) {
    /* This function returns the top value from the stack */
    if (stack_pointer <= 0) send (CRASH);
    return (val[--stack_pointer]);
}


real push_stack(real value_to_push) {
    if (stack_pointer > SIZEOF_STACK) {
        send (CRASH);
    }
    return (val[stack_pointer++] = value_to_push);
}


void eqn_init(void) {
    matrix(variable_table, 200, 200, char);
    array(value_table, 200, real);
}


int variable_lookup(char *item) {
    int i;

    for (i = 0; i < table_size; i++) {
        if (streq(item, variable_table[i])) {
            return (i);
        }
    }
    BADEQN_errpos = -1;
    BADEQN_errmsg = "Error - variable name not known";
    send (BADEQN);
    return 0; /* not reached */
}


real value_lookup(int index) {
    if (index < table_size) {
        return (value_table[index]);
    } else {
        BADEQN_errpos = -1;
        BADEQN_errmsg = "Error - no value for variable";
        send (BADEQN);
    }
    return 0; /* not reached */
}
