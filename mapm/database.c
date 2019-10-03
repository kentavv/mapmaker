/******************************************************************************

 #####     ##     #####    ##    #####     ##     ####   ######           ####
 #    #   #  #      #     #  #   #    #   #  #   #       #               #    #
 #    #  #    #     #    #    #  #####   #    #   ####   #####           #
 #    #  ######     #    ######  #    #  ######       #  #        ###    #
 #    #  #    #     #    #    #  #    #  #    #  #    #  #        ###    #    #
 #####   #    #     #    #    #  #####   #    #   ####   ######   ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#define INC_LIB
#define INC_SHELL
#include "mapm.h"

/* undefine the following if you want to enable a database hookup */

#define NO_DATABASE



command import()
{
    char *out_name, *str, first_token[TOKLEN+1];
    int length, prev_data, filenum;
    FILE *fp=NULL;

#ifndef NO_DATABASE

    /* command initialization */
    mapm_ready(MAYBE_DATA, 0, 0, NULL);

    /* set values to NULL before memory allocation */
    raw.locus_name= NULL;
    raw.data.f2.allele= NULL;
    raw.data.f2.allelic_distribution= NULL;
    modified= NULL;
    class= NULL;

    /* get command line argument = ".data" file name 
       to dump imported database data into */
    out_name = get_temp_string();
    get_one_arg(stoken, sREQUIRED, out_name);

    /* force file extension to be ".data" */
    make_filename(out_name, FORCE_EXTENSION, ".data");

    /* make a function call to a procedure which connects to the database */
    get_from_database();

    raw.filenumber = new_magic_number();
    strcpy(raw.filename, out_name);

    print("data from database successfully loaded\n");

    /* writes all the data files */
    try_to_unload(NULL,TRUE,FALSE,TRUE);

    allocate_order_data(raw.num_markers);
    allocate_mapping_data(raw.num_markers);
    allocate_seq_stuff(raw.num_markers);

    allocate_two_pt(raw.num_markers);
    allocate_three_pt(raw.num_markers);
    reset_state()
#else
    printf("no database is currently attached to MAPMAKER\n");
#endif

}


void get_from_database()
{
    int num_indivs, num_markers, cross_type;

    /* must allocate MAPMAKER data structures as follows */

    /* read in num_markers, num_indivs, cross_type from your database *
       cross type must be one of F2_INTERCROSS, F3_SELF, F2_BACKCROSS, 
       RI_SIB, RI_SELF */
    
    /* allocate and set these portions of the data structure */

    allocate_f2_data(num_markers, num_indivs);
    raw.num_markers = num_markers;
    raw.data.f2.num_indivs = num_indivs;
    raw.data_type = F2;
    raw.data.f2.cross_type = cross_type;
	
    /* now fill in the data structure (raw), obtaining the data from your database and
       filling it in the structure as is done in the procedure new_read_f2_raw() in 
       the file reader.c */



}


command export()
{
    /* this command can be designed to send mapping results back to your 
       database if this functionality is required */

#ifndef NO_DATABASE

    mapm_ready(F2,0,0,NULL);

    /* dump out any mapping results to your database */

#else

    printf("no database is currently attached to MAPMAKER\n");

#endif

}
