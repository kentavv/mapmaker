/******************************************************************************
  
   ####    ####   #    #  #####    ####            ####
  #    #  #    #  ##  ##  #    #  #               #    #
  #    #  #       # ## #  #    #   ####           #
  #  # #  #       #    #  #    #       #   ###    #
  #   #   #    #  #    #  #    #  #    #   ###    #    #
   ### #   ####   #    #  #####    ####    ###     ####
  
******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#define INC_LIB
#define INC_SHELL
#define INC_EQN
#define INC_STATS
#define INC_CALLQCTM
#define INC_QTOPLEVEL
#define INC_QLOWLEVEL
#include "qtl.h"

bool isa_test_wiggle(int wiggle_num);

/***** GLOBAL *****/
int print_long_maps;
char default_intercross_chars[10], default_backcross_chars[10];
char geno_chars[10];
bool auto_save;
bool fix_weight_kludge;
real min_trait_val, max_trait_val;

/***** INTERNAL *****/
char *temp, errmsg[MAXLINE+1];
int get_pair_entry();
void load_qtl_files(), save_qtl_files();
bool save_on_exit();
void print_ps_wiggle_order(), print_ps_multi_wiggle();

/***** QTL Commands in this file *****/
command set_intervals();	
command show_seq_history();
command set_trait();		
command set_genetics();
command load_data();		
command qtl_map();		
command singles(); /* FROB */
command qtl_like();		
command compare();		
command pheno_histogram();	
command wiggle();	
command test_wiggle();	
command make_trait();
command predict();
command forget_trait();
command forget_compare();
command forget_all_compares();
command forget_all_scans();
command forget_scan();
command set_units();	
command set_print_iter();	
command set_print_long();	
command set_dbg_qctm();	
command set_max_ints();
command set_brute();		
command set_print_bf();	
command set_print_nam();	
command set_print_scans();
command set_print_rm();	
command set_tolerance();	
command set_bag_qctm();	
command set_autosave();	
command set_more_mode();
command translate();		
command let();
command names();
command forget();
command sequence_editor();
command show_peaks();
command show_wiggle();
command draw_wiggle();
command show_trait();
command list_traits();
command list_wiggles();
command save_status();
command show_trait(); 
command show_test_wiggle();
command show_best();
command show_compare();
command list_compares();
command new_show_map();
command tester();
command dump_traits();
command dump_genome();
command dump_scan();
command tweak_weight();
command set_min_trait();
command set_max_trait();

/* QTL Command Topics */
#define RAWDATA  1
#define SELECT   2
#define MAPPING  3
#define DATAVIEW 4
#define OPTION   5
#define SYSTEM   6
#define ALIAS    7
#define WIZARD   8
#define HELPONLY 9

DATA *data;

/***** Initialize all user variables and the command tables *****/
void cmd_init() {
    
    array(temp,MAXLINE+1,char);

/********** NO LONGER USE THIS ***********************************************
    mktopic(SELECT,  "Interval (sequence) and Trait Selection Command",  CMD);
    mktopic(MAPPING, "QTL-Map Making Function",                          CMD);
    mktopic(DATAVIEW,"QTL-Map Data Viewing Function",                    CMD);
    mktopic(OPTION,  "MAPMAKER/QTL Option",                              OPT);
    mktopic(SYSTEM,  "System Command",                                   CMD);
    mktopic(RAWDATA, "Raw Data Handling Command",                        CMD);
    mktopic(ALIAS,   "Command Aliase",                                   WIZ);
    mktopic(WIZARD,  "Wizard Command",                                   WIZ);
    mktopic(HELPONLY,"Other Help Topics",                                HLP);
*****************************************************************************/
    
/* COMMANDS:   command-name              abbrev c-function       topic   code*/
/*             1234567890123456789012345   123  1234567890123456 12345678 123*/
/*            --------------------------- ----- ---------------- -------- ---*/
/*  mkcommand("1234567890123456789012345","123",1234567890123456,12345678,1);*/
    
    /* system commands in shell.c... */
    mkcommand("help",                     "?",  help,            CMD);
    mkcommand("quit",                     "q",  quit,            CMD);
    mkcommand("run",                      "r",  run_from_file,   CMD);
    mkcommand("photo",                    "",   do_photo,        CMD);
    mkcommand("change directory",         "cd", cd_command,      CMD);
    mkcommand("system",                   "!",  system_command,  CMD);
    mkcommand("previous commands",        "",   show_cmd_history,CMD);
    mkcommand("review output",            "",   review_output,   CMD);
    mkcommand("time",                     "",   show_time,       CMD);
    mkcommand("wizard mode",              "",   set_wizard,      CMD);
    mkcommand("comment",                  "",   comment,         CMD);
    mkcommand("about mapmaker/qtl",       "",   about,           CMD);

    /* QTL specific commands in this file... */
    mkcommand("scan",                     "",   wiggle,          CMD);
    mkcommand("map",                      "m",  qtl_map,         CMD);
    mkcommand("like",                     "",   qtl_like,        WIZ);
    mkcommand("compare",                  "",   compare,         CMD);

    mkcommand("list scans",               "",   list_wiggles,    CMD);
    mkcommand("show scan",                "",   show_wiggle,     CMD);
    mkcommand("draw scan",                "",   draw_wiggle,     CMD);
    mkcommand("show peaks",               "",   show_peaks,      CMD);
    mkcommand("show trys",                "",   show_test_wiggle,CMD);
    mkcommand("forget scan",              "",   forget_scan,     CMD);
    mkcommand("forget all scans",         "",   forget_all_scans,CMD);

    mkcommand("list compares",            "",   list_compares,   CMD);
    mkcommand("show best maps",           "",   show_best,       CMD);
    mkcommand("show compare",             "",   show_compare,    CMD);

/***** ALIASES - ELIMINATE 
    mkcommand("wiggle",                   "",   wiggle,          CMD);
    mkcommand("list wiggles",             "",   list_wiggles,    CMD);
    mkcommand("show wiggle",              "",   show_wiggle,     CMD);
    mkcommand("draw wiggle",              "",   draw_wiggle,     CMD);
*****/

    mkcommand("sequence",                 "s",  set_intervals,   CMD);
    mkcommand("history",                  "h",  show_seq_history,CMD);
    mkcommand("trait",                    "t",  set_trait,       CMD);
    mkcommand("let",                      "",   let,             CMD);
    mkcommand("names",                    "",   names,           CMD);
    mkcommand("forget named sequence",    "",   forget,          CMD);
    mkcommand("edit sequence",            "",   sequence_editor, CMD);


/*    mkcommand("prepare data",             "",   prep_data,       CMD); */
    mkcommand("load data",                "ld", load_data,       CMD);
    mkcommand("save data",                "",   save_status,     CMD);
    mkcommand("auto save data",           "",   set_autosave,    OPT);
    mkcommand("list loci",                "",   translate,       CMD);
    mkcommand("list traits",              "",   list_traits,     CMD);
    mkcommand("show linkage map",         "",   new_show_map,    CMD);
    mkcommand("show trait",               "",   show_trait,      CMD);
    mkcommand("make trait",               "",   make_trait,      CMD);
    mkcommand("forget trait",             "",   forget_trait,    CMD);

    mkcommand("forget compare",           "",   forget_compare,  CMD);
    mkcommand("forget all compares",      "",forget_all_compares,CMD);

    mkcommand("units",                    "",   set_units,       PAR);
    mkcommand("print names",              "",   set_print_nam,   OPT);
    mkcommand("print scans",              "",   set_print_scans, WIZ);
    mkcommand("more mode",                "",   set_more_mode,   OPT);

    mkcommand("print rec mat",            "",   set_print_rm,    WIZ);
    mkcommand("brute force",              "",   set_brute,       WIZ);
    mkcommand("bag qctm",                 "",   set_bag_qctm,    WIZ);
    mkcommand("debug qctm",               "",   set_dbg_qctm,    WIZ);
    mkcommand("print brute force",        "",   set_print_bf,    WIZ);
    mkcommand("print bf",                 "",   set_print_bf,    WIZ);
    mkcommand("print long maps",          "plm",set_print_long,  WIZ);
    mkcommand("print iterations",         "pi", set_print_iter,  WIZ);
    mkcommand("tolerance",                "",   set_tolerance,   WIZ);

    mkcommand("test",                     "",   tester,          WIZ);
    mkcommand("predict trait",            "",   predict,         WIZ);
    mkcommand("dump traits",              "",   dump_traits,     WIZ);
    mkcommand("dump genome",              "",   dump_genome,     WIZ);
    mkcommand("dump scan",                "",   dump_scan,       WIZ);
    mkcommand("tweak weight",             "",   tweak_weight,    WIZ);
    mkcommand("min trait",                "",   set_min_trait,   WIZ);
    mkcommand("max trait",                "",   set_max_trait,   WIZ);

    /* FROB */
    mkcommand("singles",                  "",   singles,         WIZ);

    /* QTL User state variables */	/* source file defined in */
    units= CENTIMORGANS;     		/* qtop.c */
    print_names= FALSE;			/* qtop.c */
    print_mapm_loci = TRUE;             /* qtop.c */
    print_scans = TRUE;                 /* qtop.c */
    print_maps= FALSE;			/* qtop.c */
    print_long_maps= FALSE;		/* qcmds.c - here */
    trait= NOTRAIT;			/* qcmds.c */
    auto_save= TRUE;                    /* qcmds.c */
    more_mode= TRUE;                    /* qcmds.c */
    
    print_iter= FALSE;			/* qctm.c */
    print_rec_mat= FALSE;		/* qctm.c */
    debug_qctm= FALSE;			/* qctm.c */
    bag_qctm= FALSE;			/* qctm.c */
    pos_tolerance= 0.001;		/* qctm.c */
    like_tolerance= 0.001;		/* qctm.c */
    mat_tolerance= 0.00001;		/* qctm.c */
    brute_force= TRUE;			/* qctm.c */
    print_brute_force= FALSE;		/* qctm.c */
    
    /* ints, ints_string, and seq_history are set up in seq_init() */
    
    /* Other global state variables */  /* source file */
    max_intervals= -5;			/* qctm.c *** Negative -> default val*/
    max_continuous_vars= -5;		/* ditto */
    wizard_mode= FALSE;              	/* shell.c */
    quit_save_hook= save_on_exit;       /* shell.c - function to call */
    fix_weight_kludge= FALSE;
    min_trait_val= -VERY_BIG; max_trait_val=VERY_BIG;
} 


void help_init()
{
/* #include "qtl.code" - outdated help mechanism */
}


/********** These are the QTL commands ***********/

command set_print_iter() { maybe_set_bool(&print_iter); }
command set_print_rm()	 { maybe_set_bool(&print_rec_mat); }
command set_print_long() { maybe_set_bool(&print_long_maps); }
command set_print_nam()	 { maybe_set_bool(&print_names); }
command set_print_scans(){ maybe_set_bool(&print_scans); }
command set_brute()	 { maybe_set_bool(&brute_force); }
command set_print_bf()	 { maybe_set_bool(&print_brute_force); }
command set_dbg_qctm()	 { maybe_set_bool(&debug_qctm); }
command set_bag_qctm()	 { maybe_set_bool(&bag_qctm); }
command set_tolerance()	 { maybe_set_real(&like_tolerance,0.0,VERY_BIG,8.6); }
command set_autosave()   { maybe_set_bool(&auto_save); }
command set_more_mode()  { maybe_set_bool(&more_mode); }

command set_min_trait()	 { maybe_set_real(&min_trait_val,0.0,VERY_BIG,8.6); }
command set_max_trait()	 { maybe_set_real(&max_trait_val,0.0,VERY_BIG,8.6); }

command set_max_ints() /* BUG? e can only do this before qctm_init! */
{ maybe_set_int(&max_intervals,1,7); } /* limit is discussed in qctm.c */

command set_units()
{
    if (nullstr(args)) {
	if(units == RECFRACS)
	  print("the 'units' are currently set to recombination-fractions.\n");
	else 
	  print("the 'units' are currently set to haldane centimorgans.\n");

    } else { /* have args */
	crunch(args);
	if (matches(args,"recombination fractions") || 
	    matches(args,"rec-fracs") ||
	    matches(args,"recombination-fractions") ||
	    matches(args,"rec fracs") ||
	    matches(args,"rf")) units = RECFRACS;
	else if (matches(args,"centimorgans") || 
		 matches(args,"cm") ||
		 matches(args,"haldane centimorgans") ||
		 matches(args,"haldane cm")) units= CENTIMORGANS;
	else usage_error(1);
	if(units == RECFRACS)
	  print("the 'units' are now set to recombination-fractions.\n");
	else 
	  print("the 'units' are now set to haldane centimorgans.\n");
    }
}


command set_intervals()
{
    int errpos;
    
    qtl_ready(ANY_DATA,NOSEQ,NOTRAIT,NOQCTM);
    if (nullstr(args)) {
	sprintf(ps, "The sequence is%s%s'\n",
            (len(ints_string)<55 ? " '":":\n'"), ints_string); pr();

    } else if (set_qtl_sequence(args,errmsg,&errpos)) {
	sprintf(ps, "The sequence is now%s%s'",
            (len(ints_string)<55 ? " '":":\n'"), ints_string); 
	maybe_ok(ps);

    } else {
	print("Error in sequence '"); print(args); print("'\n ");
	space((errpos+20) % LINE); print("^\n"); print(errmsg); nl();
    }
}


command show_seq_history()
{
    int i, num_to_print, first_to_print, printed_any; 
    char *cmd_str;
    TABLE *seq_history;

    get_one_arg(itoken,1000,&num_to_print);
    if (num_to_print<1) usage_error(1);

    printed_any=FALSE;
    first_to_print= imaxf(context[active_context]->seq_history_num-num_to_print,0);
    seq_history= context[active_context]->sequence_history;
    for(Te=seq_history->list; Te!=NULL; Te=Te->next) {
	i=Te->id.num;  cmd_str=Te->string; 
	if (i>=first_to_print && i<context[active_context]->seq_history_num) {
	    if(!printed_any){print("Previous sequences:\n"); printed_any=TRUE;} 
	    sprintf(ps, "%3d  ", i + 1); pr(); print(cmd_str); nl(); 
	}
    }
    if (!printed_any) print("No sequences have yet been entered.\n");
}


command set_trait()
{
    qtl_ready(ANY_DATA,NOSEQ,NOTRAIT,NOQCTM);
    run {
	if (nullstr(args)) {
	    if (trait<0) print("The current trait is not set.\n");
	    else { print("The current trait is: "); print(trait_str()); nl(); }

	} else { /* there was an argument */
	    set_trait_spec(args);
	    sprintf(ps, "The current trait is now: %s", trait_str());
	    maybe_ok(ps);
	} 
    }
    except_when(BADTRAIT) {
	sprintf(ps, "Bad trait name or number specified.\n%s", BADTRAIT_errmsg); 
	error(ps);
    }
}


command translate()
{
    int printed, i;
    
    qtl_ready(ANY_DATA,NOSEQ,NOTRAIT,NOQCTM);
    nomore_args(0);
    
    printed= 0;
    print(BIG_DIVIDER); print("LOCI:\n");
    for (i=0; i<raw.n_loci; i++) {
	if ((printed++)%5==0) nl();
	sprintf(ps, "%3d %-10s ", raw.original_locus[i], raw.locus_name[i]); pr();
    }
    nl();

/*    nl(); print("\nTRAITS:\n"); printed=0;
    for (i=0; i<raw.n_traits; i++) {
	if ((printed++)%5==0) nl();
	if (nullstr(raw.trait_name[i])) sprintf(ps,"%3d <deleted>  ",i+1);
	  else sprintf(ps,"%3d %-10s ",i+1,raw.trait_name[i]); 
	pr();
    } nl(); */
    print(BIG_DIVIDER);
    
    /* Old code - what to do if there are args? 
    while (stoken(&args,sREQUIRED,temp)) { 
	if (valid_locus_str(temp,&i,&errmsg)) {
	    sprintf(ps," locus %d %-10s\n",i+1,raw.locus_name[i]); print(ps);
	} else if (valid_trait_str(temp,&i,&errmsg)) {
	    sprintf(ps," trait %d %-10s\n",i+1,raw.trait_name[i]); print(ps);
	} else {
	    sprintf(ps," Error: \"%s\" doesn't match, or is ambiguous.\n",
	       temp); print(ps);
	}
    } */
}


command qtl_map()
{
    int perm;
    
    qtl_ready(ANY_DATA,SEQ,TRAIT,QCTM);
    nl(); 
    
    for_all_orders(ints,map,perm) {
	make_qtl_map(map);
	print_map_divider();
	print_trait(1); nl();
	print_qtl_map(map,free_genetics);
    }
    print_map_divider();
}


command singles() /* FROB */
{
    int perm, last;
    real threshold, scale;
    
    qtl_ready(ANY_DATA,SEQ,TRAIT,QCTM);
    get_arg(rtoken,2.0,&threshold);
    get_arg(rtoken,0.25,&scale);
    last= num_intervals-1;

    nl(); 
    print_trait(num_orders); nl(); 
    nl();
    print_short_title();
    
    for_all_orders(ints,map,perm) {
	map->fix_pos[last]= 0.0;
	make_qtl_map(map);
	print_short_qtl_map(map,threshold,scale);
    }
    nl();
}


command qtl_like()
{
/*
    int perm;
    real like;
    */

/*    qtl_ready(ANY_DATA,SEQ,TRAIT,QCTM);
    nl(); 
    
    for_all_orders(ints,map,perm) {
	 like= qtl_map_like(map); 
	sprintf(ps,"like= %-7.2lf\n",like); pr();
    } */
}


command compare()
{
    int perm, i, comp;
    
    qtl_ready(ANY_DATA,SEQ,TRAIT,QCTM);
    run {
	comp=allocate_compare_struct(trait,ints,ints_string,num_intervals,
				     num_orders);
	
	for_all_orders(ints,map,perm) {
	    if (perm==FIRST || perm>CONTIG) 
	      { i=0; nl(); print_tiny_map(map,-1,0.0); } /* print tiny title */
	    
	    make_qtl_map(map);
	    print_tiny_map(map,i++,0.0);
	    store_compare_map(comp,map,(perm<=CONTIG));
	}
    } on_exit {
	if(msg == NOMEMORY || msg == INTERRUPT || msg == CRASH)
	  bash_compare_struct(num_compares-1);
    }
}


command list_compares()
{
    if (first_compare==num_compares) 
      print("No compare results have been saved.\n");
    else { print("\nSaved compare results:\n\n"); print_saved_compares(); }
}

#define SAVED_COMP_BEST  "Sorted compare results number%s %d.%d"
#define SAVED_COMP_COMPS  "-%d.%d"
#define SAVED_COMP_TRAIT  " for trait %d (%s).\n"
#define SAVED_COMP_NUM    "Saved compare number %d.%d:\n"
#define SAVED_COMP_THRESH "LOD Threshold: %-4.2lf  Falloff: %-4.2lf\n"

command show_best()
{
    char arg[TOKLEN+1];
    int compare, contig, i;
    real threshold, falloff;
    COMPARE_OPERATION *op;

    qtl_ready(ANY_DATA,NOSEQ,NOTRAIT,NOQCTM);
    get_arg(stoken,"",arg);
    get_arg(rtoken,2.0,&threshold);
    get_arg(rtoken,-1.0,&falloff);
    nomore_args(num_args);
    if (threshold<0.0 || falloff==0.0 || falloff>=threshold) 
      usage_error(num_args);

    get_compare_nums(arg,&compare,&contig);
    if (compare<0) { 
	compare= num_compares-1; 
	contig= compares[compare]->num_contigs-1;
    } else if (contig<0 && compares[compare]->num_contigs==1) contig=0;
    op=compares[compare]; nl();

    if (contig<0) {
	sprintf(ps, SAVED_COMP_BEST, maybe_s(op->num_contigs), compare + 1, 1); pr();
	if (op->num_contigs>1)
	  { sprintf(ps, SAVED_COMP_COMPS, compare + 1, op->num_contigs); pr(); }
	sprintf(ps, SAVED_COMP_TRAIT, op->trait + 1, raw.trait_name[op->trait]); pr();
	print_old_seq(op->seq_string); 
	sprintf(ps, SAVED_COMP_THRESH, threshold, falloff); pr();
	for (i=0; i<op->num_contigs; i++) {
	    nl(); sprintf(ps, SAVED_COMP_NUM, compare + 1, i + 1); pr(); 
	    print_best_saved_maps(compare,i,threshold,falloff);
	}

    } else {
	sprintf(ps, SAVED_COMP_BEST, "", compare + 1, contig + 1); pr();
	sprintf(ps, SAVED_COMP_TRAIT, op->trait + 1, raw.trait_name[op->trait]); pr();
	print_old_seq(op->seq_string); 
	sprintf(ps, SAVED_COMP_THRESH, threshold, falloff); pr(); 
	nl(); print_best_saved_maps(compare,contig,threshold,falloff);
    }
}


#define SAVED_COMP_WHICH  "Compare results number%s %d.%d"

command show_compare()
{
    char arg[TOKLEN+1];
    int compare, contig, i;
    COMPARE_OPERATION *op;

    qtl_ready(ANY_DATA,NOSEQ,NOTRAIT,NOQCTM);
    get_arg(stoken,"",arg);
    nomore_args(num_args);

    get_compare_nums(arg,&compare,&contig);
    if (compare<0) { 
	compare= num_compares-1; 
	contig= compares[compare]->num_contigs-1;
    } else if (contig<0 && compares[compare]->num_contigs==1) contig=0;
    op=compares[compare]; nl();

    if (contig<0) {
	sprintf(ps, SAVED_COMP_WHICH, maybe_s(op->num_contigs), compare + 1, 1); pr();
	if (op->num_contigs>1)
	  { sprintf(ps, SAVED_COMP_COMPS, compare + 1, op->num_contigs); pr(); }
	sprintf(ps, SAVED_COMP_TRAIT, op->trait + 1, raw.trait_name[op->trait]); pr();
	print_old_seq(op->seq_string); 
	for (i=0; i<op->num_contigs; i++) {
	    nl(); sprintf(ps, SAVED_COMP_NUM, compare + 1, i + 1); pr(); 
	    print_saved_maps(compare,i);
	}

    } else {
	sprintf(ps, SAVED_COMP_WHICH, "", compare + 1, contig + 1); pr();
	sprintf(ps, SAVED_COMP_TRAIT, op->trait + 1, raw.trait_name[op->trait]); pr();
	print_old_seq(op->seq_string); 
	nl(); print_saved_maps(compare,contig);
    }
}


#define WIGGLE_STORED  "Results have been stored as scan number %d.\n"
#define WIGGLES_STORED "Results have been stored as scan numbers %d.1-%d.%d.\n"
#define THRESHOLD_AND_SCALE "LOD threshold: %-4.2lf  Scale: %-4.2lf per '*'\n"

command wiggle()
{
    int perm, wiggle;
    real inc, threshold, scale;

    
    qtl_ready(ANY_DATA,WIGSEQ,TRAIT,QCTM);
    get_arg(rtoken,2.0,&inc);
    get_arg(rtoken,2.0,&threshold);
    get_arg(rtoken,0.25,&scale);
    nomore_args(num_args);
    if (inc<0.0 || threshold<0.0 || scale<0.0) usage_error(num_args); 
    
    if (print_scans) {
	nl(); print_trait(2); print_seq();
	sprintf(ps, THRESHOLD_AND_SCALE, threshold, scale); pr();
    }
    run { 
	wiggle=allocate_wiggle_struct(trait,ints,ints_string,num_intervals,
				      num_orders,num_ints_to_wiggle);
	
	for_wiggle_orders(ints,map,inc,perm) { 
	    
	    if (print_scans) {
		if (perm==FIRST || perm==TEST || perm==ORDER) {
		    if (perm!=FIRST) { print_wiggle_interval(NULL); }
		    nl(); print_wiggle_left_seq(map); 
		    print_wiggle_genetics(&map->constraint[num_intervals-1]); 
		    nl();
		    print_wiggle_title();   
		    print_wiggle_interval(map);
		    store_wiggle_interval(wiggle,map,TRUE,FALSE,inc);
		}
		else if (perm==CONTIG || perm==SKIP) {
		    if (perm==SKIP) { print_wiggle_interval(NULL); nl(); }
		    print_wiggle_interval(map); 
		    store_wiggle_interval(wiggle,map,FALSE,(perm==CONTIG),inc);
		}
	    }
	    else {
		if (perm==FIRST || perm==TEST || perm==ORDER)
		  store_wiggle_interval(wiggle,map,TRUE,FALSE,inc);
		else if (perm==CONTIG || perm==SKIP)
		    store_wiggle_interval(wiggle,map,FALSE,(perm==CONTIG),inc);
	    }
	    make_qtl_map(map);
	    if (print_scans) print_wiggle_map(map,threshold,scale); 
	    store_wiggle_point(wiggle,map);
	} 
    
	if (print_scans)
	  print_wiggle_interval(NULL); nl(); /* line on the bottom */
	if (wiggle>=0) {
	    if (num_orders==1) sprintf(ps, WIGGLE_STORED, wiggle + 1); 
	    else sprintf(ps, WIGGLES_STORED, wiggle + 1, wiggle + 1, num_orders);
	    pr();
	}
    } on_exit {
	if(msg == NOMEMORY || msg == INTERRUPT || msg == CRASH)
	  bash_wiggle_struct(num_wiggles-1); /* the current one */
    }
}



#define SAVED_WIGGLE_WIG    "Scan results number%s %d.%d"
#define SAVED_WIGGLE_WIGS   "-%d.%d"
#define SAVED_WIGGLE_TRAIT  " for trait %d (%s).\n"

command list_wiggles()
{
    char arg[TOKLEN+1];
    int wiggle, order;
    WIGGLE_OPERATION *op;

    qtl_ready(ANY_DATA,NOSEQ,NOTRAIT,NOQCTM);
    get_arg(stoken,"",arg); nomore_args(num_args);

    if (first_wiggle==num_wiggles) print("No scan results have been saved.\n");
    if (nullstr(arg)) 
      { print("\nSaved scan results:\n\n"); print_saved_wiggles(); return; }

    get_wiggle_nums(arg,&wiggle,&order);
    if (wiggle<0) { wiggle=num_wiggles-1; order=wiggles[wiggle]->num_orders-1;}
    else if (order<0 && wiggles[wiggle]->num_orders==1) order=0;
    op=wiggles[wiggle]; nl();

    sprintf(ps, SAVED_WIGGLE_WIG, maybe_s(op->num_orders), wiggle + 1, 1); pr();
    if (op->num_orders>1) 
      { sprintf(ps, SAVED_WIGGLE_WIGS, wiggle + 1, op->num_orders); pr(); }
    sprintf(ps, SAVED_WIGGLE_TRAIT, op->trait + 1, raw.trait_name[op->trait]); pr();
    print_old_seq(op->seq_string); nl(); 
    print_saved_wiggle(wiggle);
}


command show_wiggle()
{
    char arg[TOKLEN+1];
    int wiggle, order, last;
    real threshold, scale;
    WIGGLE_OPERATION *op;

    qtl_ready(ANY_DATA,NOSEQ,NOTRAIT,NOQCTM);
    get_arg(stoken,"",arg);
    get_arg(rtoken,2.0,&threshold);
    get_arg(rtoken,0.25,&scale);
    nomore_args(num_args);

    if (threshold<0.0 || scale<=0.0) usage_error(2);
    get_wiggle_nums(arg,&wiggle,&order);
    if (wiggle<0) { wiggle=num_wiggles-1; order=wiggles[wiggle]->num_orders-1;}
    else if (order<0 && wiggles[wiggle]->num_orders==1) order=0;
    op=wiggles[wiggle]; nl();

    if (order<0) { 
	sprintf(ps, SAVED_WIGGLE_WIG, maybe_s(op->num_orders), wiggle + 1, 1); pr();
	if (op->num_orders>1) 
	  { sprintf(ps, SAVED_WIGGLE_WIGS, wiggle + 1, op->num_orders); pr(); }
	sprintf(ps, SAVED_WIGGLE_TRAIT, op->trait + 1, raw.trait_name[op->trait]); pr();
	print_old_seq(op->seq_string); nl(); 
	print_saved_wiggle(wiggle);
    } else {
	sprintf(ps, SAVED_WIGGLE_WIG, "", wiggle + 1, order + 1); pr();
	sprintf(ps, SAVED_WIGGLE_TRAIT, op->trait + 1, raw.trait_name[op->trait]); pr();
	print_old_seq(op->seq_string); 
	sprintf(ps, THRESHOLD_AND_SCALE, threshold, scale); pr(); 
	nl(); print_wiggle_left_seq(op->data[order][0]->map);
	last= op->data[order][0]->map->num_intervals-1;
	print_wiggle_genetics(&op->data[order][0]->map->constraint[last]);
	nl(); print_saved_wiggle_order(wiggle,order,threshold,scale);
    }
}
    

#define WHICH_WIGGLE  \
 "%d is an ambiguous scan number.\nUse a number from %d.1-%d.%d."
#define PEAKS_TITLE  "LOD score peaks for scan %d.%d of trait %d (%s).\n"
#define PEAKS_PARAMS "Peak Threshold: %-4.2lf  Falloff: %-4.2lf\n"
#define PEAKS_NAME   "Peak name %s1 has been set.\n"
#define PEAKS_NAMES  "Peak names %s1-%s%d have been set.\n"
#define PEAKS_NONE   "No peaks with a LOD score above the threshold were found!\n(Existing names remain unchanged).\n"

command show_peaks()
{	
    char arg[TOKLEN+1],*name;
    int wiggle, order, last, peak;
    real threshold, falloff;
    WIGGLE_PEAK *peaks=NULL, *p;

    name = NULL;  
    qtl_ready(ANY_DATA,NOSEQ,NOTRAIT,QCTM); /* find_wiggle_peaks calls QCTM */
    get_arg(stoken,"",arg);
    get_arg(rtoken,2.0,&threshold);
    get_arg(rtoken,-2.0,&falloff);
    nomore_args(num_args);

    if (threshold<0.0 || falloff==0.0 || falloff>threshold) 
      usage_error(num_args);
    get_wiggle_nums(arg,&wiggle,&order);
    if (wiggle<0) { wiggle=num_wiggles-1; order=wiggles[wiggle]->num_orders-1;}
    else if (order<0 && wiggles[wiggle]->num_orders==1) order=0;
    if (order<0) { sprintf(ps, WHICH_WIGGLE, wiggle + 1, wiggle + 1, wiggle + 1,
		      wiggles[wiggle]->num_orders-1); error(ps); }

    run {
	array(name, NAME_LEN+1, char);
	peaks=find_wiggle_peaks(wiggle,order,threshold,-2.0,falloff,2.0,TRUE); 
	nl(); 
	sprintf(ps, PEAKS_TITLE, wiggle + 1, order + 1, wiggles[wiggle]->trait + 1,
            raw.trait_name[wiggles[wiggle]->trait]); pr();
	print_old_seq(wiggles[wiggle]->seq_string);
	print_wiggle_left_seq(wiggles[wiggle]->data[order][0]->map);
	last= wiggles[wiggle]->data[order][0]->map->num_intervals-1;
	print_wiggle_genetics(
	  &wiggles[wiggle]->data[order][0]->map->constraint[last]);
	sprintf(ps, PEAKS_PARAMS, threshold, falloff); pr(); nl(); 

	for (peak=0, p=peaks; p!=NULL; p=p->next, peak++) 
	  print_peak(p,peak);
	if (peak>0) {
	    print_map_divider(); nl();
	    nstrcpy(name,raw.trait_name[wiggles[wiggle]->trait],7);
	    strcat(name,".");
	    name_peaks(peaks,name,TRUE);
	    if (peak==1) sprintf(ps, PEAKS_NAME, name);
	    else sprintf(ps, PEAKS_NAMES, name, name, peak + 1);
	} else print(PEAKS_NONE);

    } when_aborting { free_wiggle_peaks(peaks); unarray(name, char) }
}


command new_show_map()
{
    int perm,i,last,first = TRUE;
    real inc,dist;

    qtl_ready(ANY_DATA,WIGSEQ,TRAIT,QCTM);

    inc = 999.0; last= -10; /* just one per interval */

    print("\nlinkage maps:\n");
    print(MAP_DIVIDER);
    for_wiggle_orders(ints,map,inc,perm) {
	if(map->left[0] != last && !first)
	  print(MAP_DIVIDER);
	first = FALSE;
	for(i = 0; i < map->num_intervals; i++) {
	    dist = map_length(map->left[i], map->right[i]);
	    if(dist <= .490) {
		sprintf(ps, " %s   %5.1lf cM   %4.1lf %%\n",
                interval_str(map->left[i],map->right[i],TRUE),
                haldane_cm(dist),dist*100.0);
		pr();
	    }
	    else {
		print(MAP_DIVIDER);
		print("MAP:\n\n");
	    }
	}
	last = map->right[i-1];
    }
    print(MAP_DIVIDER);
    nl();
}


#define TEST_TITLE  "Test genetics results for trait %d (%s).\n" 
#define TEST_PARAMS "Scan numbers: %d.%d-%d.%d  Threshold: %-4.2lf\n"
#define TEST_ISNT   \
"Scan number %d.x did not use test genetics in the scanned interval."

command show_test_wiggle()
{	
    char arg[TOKLEN+1];
    int wiggle, order;
    real threshold;

    qtl_ready(INTERCROSS,NOSEQ,NOTRAIT,NOQCTM);
    get_arg(stoken,"",arg);
    get_arg(rtoken,2.0,&threshold);
    nomore_args(num_args);

    if (threshold<0.0) usage_error(2);
    get_wiggle_nums(arg,&wiggle,&order);
    if (wiggle<0) { wiggle=num_wiggles-1; order=wiggles[wiggle]->num_orders-1;}
    if (!isa_test_wiggle(wiggle)) error(TEST_ISNT);
    if (order<0 && wiggles[wiggle]->num_orders==NUM_MODELS) order=0;
    if (order<0) { sprintf(ps, WHICH_WIGGLE, wiggle + 1, wiggle + 1, wiggle + 1,
		      wiggles[wiggle]->num_orders-1); error(ps); }

	order= (order/NUM_MODELS)*NUM_MODELS; /* first order for this test */
	nl(); 
	sprintf(ps, TEST_TITLE, wiggles[wiggle]->trait + 1,
            raw.trait_name[wiggles[wiggle]->trait]); pr();
	print_old_seq(wiggles[wiggle]->seq_string);
	print_wiggle_left_seq(wiggles[wiggle]->data[order][0]->map);
	sprintf(ps, TEST_PARAMS, wiggle + 1, order + 1, wiggle + 1, order + NUM_MODELS,
            threshold); pr(); nl(); 
	print_test_wiggle_order(wiggle,order,threshold);

}


#define TRAIT_BADNAME "'%s' is not a valid name for a new trait.\n%s"
#define TRAITS_FULL "You have reached maximum number of traits.\nUse 'forget tr"
#define TRAIT_MADE "New trait number %d (%s) had been added to the data set.\n"
command make_trait()
{
    int i,j,k,adjusted_array_size;
    real *normal_array=NULL;
    NORMAL_TEST *normal_check=NULL;
    char *eqn, *name_to_check, *name;
    int trait_redone, trait_index;
    EQUATION **postfixed;

    /* I would like to use despace(uncrunched_args), but split_arglist 
       doesn't like that. Hmmmm... */
    name_to_check = get_temp_string();
    array(postfixed, MAX_EQN_SIZE, EQUATION*);
    qtl_ready(ANY_DATA,NOSEQ,NOTRAIT,NOQCTM);
    if (nullstr(args)) usage_error(0); 
    name=args; 
    if (!split_arglist(&eqn,'=')) usage_error(1);
    if (!nstoken(&name,sREQUIRED,name_to_check,NAME_LEN)) usage_error(2);

    if (!valid_new_trait_name(name_to_check,errmsg)) 
	{ sprintf(ps, TRAIT_BADNAME, name_to_check, errmsg); error(ps); }

    trait_redone = FALSE;
    for (i=0;i<raw.n_traits;i++) {
	if (nullstr(raw.trait_name[i])) {
	    trait_index = i;
	    trait_redone = TRUE;
	    break;
	}
    }
    if (!trait_redone) {
	if (raw.n_traits==raw.max_traits) error(TRAITS_FULL);
	else trait_index = raw.n_traits;
    }

    run { 
	/* The variable value and lookup stuff should be changed to read 
	   right out of the raw struct. */
	eqn_init();
	table_size = raw.n_traits; /* check bounds? */
	for (k=0; k<raw.n_traits; k++) 
	  nstrcpy(variable_table[k],raw.trait_name[k],NAME_LEN);
	postfixed = make_equation(eqn,variable_lookup);

	for (i=0; i<raw.n_indivs; i++) { /* check bounds? */
	    for (k=0; k<raw.n_traits; k++) value_table[k]=raw.trait[i][k];
	    raw.trait[i][trait_index]=
	      evaluate_equation(postfixed,value_lookup);
	}
	strcpy(raw.trait_name[trait_index],name_to_check);
	strcpy(raw.trait_eqn[trait_index],eqn);
	if (!trait_redone) raw.n_traits++;

	nl(); sprintf(ps, TRAIT_MADE, trait_index + 1, name_to_check); pr(); nl();
	array(normal_array,raw.n_indivs,real); 
	adjusted_array_size = raw.n_indivs;
	for (i=0,j=0; i<raw.n_indivs; i++) {
	    if(raw.trait[i][trait_index] == MISSING_PHENO) {
		adjusted_array_size--;
	    }
	    else {
		normal_array[j] = raw.trait[i][trait_index];
		j++;
	    }
	}
	single(normal_check, NORMAL_TEST);
	normal_check = check_normalcy(normal_array,adjusted_array_size);
	print_normal(normal_check, -1.0);
	nl(); print_rhisto(normal_array,adjusted_array_size);
    } on_exit {
	if (msg==BADEQN) { /* msg should be 0 if none sent - CHECK THIS OUT */
	    print("error in equation: "); print(eqn); nl();
	    if (BADEQN_errpos<60) { space(BADEQN_errpos+19); print("^\n"); }
	    print(BADEQN_errmsg); nl();
	}
	unarray(normal_array,real); 
	/* free normal_check???? */
	/* free postfixed - it wasn't alloced here... make it a local! */
    }
}


command predict()
{
    char *name= get_temp_string();
    int trait_redone, trait_index=0, i, perm;

    qtl_ready(INTERCROSS,SEQ,TRAIT,QCTM);
    if (nullstr(args)) usage_error(0); 
    if (!nstoken(&args,sREQUIRED,name,NAME_LEN)) usage_error(2);
    if (!valid_new_trait_name(name,errmsg)) 
	{ sprintf(ps, TRAIT_BADNAME, name, errmsg); error(ps); }

    for (i=0, trait_redone = FALSE; i<raw.n_traits; i++) 
      if (nullstr(raw.trait_name[i])) 
	{ trait_index= i; trait_redone= TRUE; break; }
    if (!trait_redone) {
	if (raw.n_traits==raw.max_traits) error(TRAITS_FULL);
	else trait_index= raw.n_traits;
    }

    if (num_orders!=1) 
      error("The current sequence spcifies more than one model.\n");
    for_all_orders(ints,map,perm) {
	get_fixed_qtl_weights(map);
	
	for (i=0; i<raw.n_indivs; i++) { 
	    raw.trait[i][trait_index]= model_prediction(map,i);
	}
    }
    nl(); sprintf(ps, TRAIT_MADE, trait_index + 1, name); pr(); nl();
    if (!trait_redone) raw.n_traits++;
    strcpy(raw.trait_name[trait_index],name);
    strcpy(raw.trait_eqn[trait_index],ints_string);
}


command dump_scan()
{
    char arg[TOKLEN+1];
    char *name=get_temp_string();
    int i, j;
    FILE *fp=NULL;
    WIGGLE_OPERATION *op;
    WIGGLE_INTERVAL **data;
    WIGGLE_POINT *point;
    int mark,position,order,wiggle;
    real threshold,scale;

    qtl_ready(ANY_DATA,NOSEQ,NOTRAIT,NOQCTM);
    get_arg(stoken,"",arg);
    get_arg(rtoken,2.0,&threshold);
    get_arg(rtoken,0.25,&scale);
    get_wiggle_nums(arg,&wiggle,&order);
    if(wiggle<0) { wiggle=num_wiggles-1;order=wiggles[wiggle]->num_orders-1;}
    else if (order<0 && wiggles[wiggle]->num_orders==1) order = 0;

    sprintf(name, "scan.%d", wiggle);
    sprintf(ps, "dumping to '%s'\n", name); pr();
    fp= open_file(name,WRITE);
    op= wiggles[wiggle];
    data = op->data[order];
    position = 0;
    mark= 0;
    for (i=0;i<op->num_wiggled_intervals;i++) {
	if (i>0 && !data[i]->contig) { mark=0;
				       position=(int)(position/1000)+1;
				       position=(position*1000); }

			       
	for (j=0;j<data[i]->num_points; j++) {
	    point=data[i]->point[j];
	    fprintf(fp,"%lf\t %6.3lf\t %7.3lf\n",position+mark*threshold,
		    point->qtl_weight,point->lod_score);
	    mark++;
	}
    }
    close_file(fp);
}

command dump_traits()
{
    int t[8], num, i, j;
    FILE *fp=NULL;
    
    qtl_ready(ANY_DATA,NOSEQ,NOTRAIT,NOQCTM);
    
    num=0;
    do {
	if (num==8) error("too many traits");
	if (!itoken(&args,iREQUIRED,&t[num]) || 
	    !valid_trait_num(--(t[num++])))  error("bad trait");
    } while (!nullstr(args));

    run {
	print("dumping to 'traits.dump'\n");
	fp= open_file("traits.dump",WRITE);
	fprintf(fp,"*indiv\t");
	for (j=0; j<num; j++) fprintf(fp,"%7s\t",raw.trait_name[t[j]]); 
	fprintf(fp,"\n");
    
	for (i=0; i<raw.n_indivs; i++) {
	    fprintf(fp,"%-4d\t",i+1);
	    for (j=0; j<num; j++) {
		if (raw.trait[i][t[j]] != MISSING_PHENO) 
		  fprintf(fp,"%7.4lf\t",raw.trait[i][t[j]]);
		else fprintf(fp,"\t");
	    }
	    fprintf(fp,"\n");
	}
	
    } on_exit {
	close_file(fp);
    }
}	


command dump_genome()
{
    FILE *fp=NULL;
    int k,l;
    int tot_aa=0,tot_bb=0,tot_ab=0,total_num=0;
    run {
	print("dumping to 'genome.dump'\n");
	fp = open_file("genome.dump",WRITE);
	fprintf(fp,"*loci\t %%AA\t %%AB  %%BB\t n_indivs");
	for(k=0;k<raw.n_loci;k++) {
	    fprintf(fp,"\n");
	    for(l=0;l<raw.n_indivs;l++) {
		if (raw.locus[k][l]=='A') { tot_aa++; total_num++; }
		else if (raw.locus[k][l]=='B') {tot_bb++; total_num++; }
		else if (raw.locus[k][l]=='H') {tot_ab++; total_num++; }
	    }
	    fprintf(fp,"%s\t %.2f\t %.2f\t %.2f\t %d",raw.locus_name[k],
		    100.*tot_aa/total_num,100.*tot_bb/total_num,100.*tot_ab/total_num,
		    total_num);
	    tot_aa=tot_bb=tot_ab=total_num=0;
	}
	fprintf(fp,"\n");
	fprintf(fp,"\n");	    
	
    } on_exit {
	close_file(fp);
	
    }
}

command list_traits()
{
    int i;
    
    print(BIG_DIVIDER); print("TRAITS:\n\n");
    for(i=0; i< raw.n_traits; i++) {
	if (nullstr(raw.trait_name[i])) sprintf(ps, "%3d <deleted>  ", i + 1);
	else sprintf(ps, "%3d %-10s %s ", i + 1, raw.trait_name[i], raw.trait_eqn[i]);
	pr();nl();
    }
    print(BIG_DIVIDER);
}


#define SHOW_WHAT "No trait specified.\nEither supply a trait name or number or use the 'trait' command."

command show_trait()
{
    real *normal_array=NULL;
    NORMAL_TEST *normal_check=NULL;
    int i, t,j, adjusted_array_size;
    char arg[TOKLEN+1], *errmsg=get_temp_string();
    
    qtl_ready(ANY_DATA,NOSEQ,NOTRAIT,NOQCTM);
    get_one_arg(stoken,"",arg); 
    if (!nullstr(arg)) { if (!valid_trait_str(arg,&t,errmsg)) error(errmsg); }
      else { if (trait<0) error(SHOW_WHAT); else t=trait; }

    run {
	sprintf(ps, "\nTrait %d (%s)", t + 1, raw.trait_name[t]); pr();
	if (!nullstr(raw.trait_eqn[t])) 
	  { sprintf(ps, "='%s'", raw.trait_eqn[t]); pr(); }
	print(":\n\n");
	array(normal_array,raw.n_indivs,real);
	adjusted_array_size = raw.n_indivs;
	for (i=0,j=0; i<raw.n_indivs; i++) {
	    if(raw.trait[i][t] == MISSING_PHENO) {
		adjusted_array_size--;
	    }
	    else {
		normal_array[j] = raw.trait[i][t];
		j++;
	    }
	}
	single(normal_check, NORMAL_TEST);
	normal_check = check_normalcy(normal_array,adjusted_array_size);
	print_normal(normal_check, -1.0);
	nl(); print_rhisto(normal_array,adjusted_array_size);
    } on_exit {
	unarray(normal_array,real); 
	/* free normal_check???? */
    }	
}


#define TRAIT_NOTMADE "Can't delete specified trait.\nTrait %d (%s) was made using the 'make trait' command."
#define REALLY_KILL "Trait %d (%s)"

command forget_trait()
{
    char trait_name[TOKLEN+1], c;
    int t;
    
    qtl_ready(ANY_DATA,NOSEQ,NOTRAIT,NOQCTM);
    get_one_arg(stoken,sREQUIRED,trait_name); 
    if (!valid_trait_str(trait_name,&t,errmsg)) error(errmsg);
    if (raw.n_traits==1) error("There is only one trait in the data set!\n.");

    sprintf(ps, "Deleting trait %d (%s)", t + 1, raw.trait_name[t]); pr();
    if (!nullstr(raw.trait_eqn[t])) { sprintf(ps, "='%s'", raw.trait_eqn[t]); pr(); }
    print("\n");
    getln("Are you sure you want to delete it? [no] ");
    if (!parse_char(&ln,"y",TRUE,&c)) return;

    raw.trait_name[t][0]='\0'; raw.trait_eqn[t][0]='\0';
    if (t==raw.n_traits-1) raw.n_traits--; 
    print("Trait Deleted.\n");
    if (t==trait) { trait= -1; print("The current trait is now not set!\n"); }
    /* need to blow away saved wiggles for this trait */
}
#define COMP_LIST_NAMES "%3d.%s  %-10s %s\n"
#define COMP_LIST_NUMS  "%3d.%s  %3d   %s\n"
#define COMP_LIST_TITLE "  NUM  TRAIT%s  SEQUENCE\n"
#define WIG_LIST_TITLE "  NUM  TRAIT%s  SEQUENCE\n"
#define WIG_LIST_NAMES "%3d.%s  %-10s  %s\n"
#define WIG_LIST_NUMS  "%3d.%s  %3d    %s\n"

command forget_compare()
{
    int i,comp_number;
    char c;

    qtl_ready(ANY_DATA,NOSEQ,NOTRAIT,NOQCTM);
    get_one_arg(itoken,iREQUIRED,&comp_number); 
    if (comp_number > num_compares) error(errmsg);
    sprintf (ps, "Deleting compare %d \n", comp_number);pr();
    i=comp_number-1;
    sprintf(ps, COMP_LIST_TITLE, (print_names ? "     " : "")); pr();
    if(print_names)
      sprintf(ps, COMP_LIST_NAMES, i + 1, (compares[i]->num_contigs == 1 
				 ? "1" : "x"),
              raw.trait_name[compares[i]->trait],
              compares[i]->seq_string);
    else
      sprintf(ps, COMP_LIST_NUMS, i + 1, (compares[i]->num_contigs == 1 
				? "1" : "x"),
	 compares[i]->trait+1, compares[i]->seq_string);
    pr();
    print("\n");
    getln("Are you sure you want to delete it? [no] ");
    if (!parse_char(&ln,"y",TRUE,&c)) return;
    bash_compare_struct(i);
}


command forget_all_compares()
{
    int i;
    char c;
    
    qtl_ready(ANY_DATA,NOSEQ,NOTRAIT,NOQCTM);
    sprintf (ps, "\nDeleting ALL compares!\n");pr();
    getln("\nAre you sure you want to delete them? [no] ");
    if (!parse_char(&ln,"y",TRUE,&c)) return;
    for (i=0;i<num_compares;i++) {
	if (compares[i]->data!=NULL) 
	  bash_compare_struct(i);
    }
    num_compares = 0;
}


command forget_scan()
{ 
    char c;
    int scan_number, t;
    WIGGLE_OPERATION *op;

    qtl_ready(ANY_DATA,NOSEQ,NOTRAIT,NOQCTM);
    get_one_arg(itoken,iREQUIRED,&scan_number);
    if (scan_number > num_wiggles) error(errmsg);
    sprintf(ps, "Deleting scan %d \n", scan_number);pr();
    t=scan_number-1;
    sprintf(ps, WIG_LIST_TITLE, (print_names ? "      " : ""));pr();
    op=wiggles[t];
    if (print_names)
      sprintf(ps, WIG_LIST_NAMES, t + 1, (op->num_orders == 1 ? "1" : "x"),
              raw.trait_name[op->trait], op->seq_string);
    else sprintf(ps, WIG_LIST_NUMS, t + 1, (op->num_orders == 1 ? "1" : "x"),
	    op->trait+1, op->seq_string);
    pr();
    print("\n");
    getln("Are you sure you want to delete it? [no] ");
    if (!parse_char(&ln,"y",TRUE,&c)) return;
    bash_wiggle_struct(t);
    print("Scan results deleted.\n");
}

command forget_all_scans()
{
    char c;
    int i;
    
    qtl_ready(ANY_DATA,NOSEQ,NOTRAIT,NOQCTM);
    sprintf (ps, "Deleting ALL saved scan results!\n");pr();
    getln("\nAre you sure you want to delete them? [no] ");
    if (!parse_char(&ln,"y",TRUE,&c)) return;
    for (i=0;i<num_wiggles;i++) {
	if (wiggles[i]->data != NULL)
	  bash_wiggle_struct(i);
    }
    num_wiggles = 0;
    print("All scan results deleted.\n");
}


command let()
{
    char *name,*seqnce,*err;

    name = args;
    split_arglist(&seqnce,'=');

    if(!name_sequence(name,seqnce,&err))
        error(err);
}

command names()
{
    /* This is a KLUDGE for now until we write a macro which ports. */
    for (Te=context[active_context]->named_sequences->list; Te!=NULL;
	 Te=Te->next) {
	sprintf(ps, "%-12s= %s\n", Te->id.name, Te->string);
	print(ps);
    }
}

command forget()
{
    char *name;
    int fail;
    
    name = get_temp_string();
    if(stoken(&args,sREQUIRED,name)) {
	if(!delete_named_entry(name,context[active_context]->named_sequences,&fail)) {
	    if(fail == NAME_DOESNT_MATCH)
	      sprintf(ps, "%s is not a defined name", name);
	    else
	      sprintf(ps, "%s is an ambiguous name", name);
	    error(ps);
	}
    }
}


#define BADEDITEDSEQ \
"An illegal sequence was specified.\nThe interval list remains '%s'.\n"

command sequence_editor()
{
    char prompt[TOKLEN+1], *new_seq;
    int errpos;

    qtl_ready(ANY_DATA,NOSEQ,NOTRAIT,NOQCTM);
    new_seq=NULL;

    run {
        /* edit line will send an error if it can't be run now */
        array(new_seq, SEQ_LEN+1, char); new_seq[0]='\0';

	strcpy(prompt,"sequence= ");
	edit_line(prompt,new_seq,SEQ_LEN,ints_string);
	if (!set_qtl_sequence(new_seq,errmsg,&errpos)) {
	    print("Error in sequence '"); print(new_seq); print("'\n ");
	    space((errpos+20) % LINE); print("^\n"); print(errmsg); nl();
	    send(BADSEQ);
	}
	print("ok\n");
	unarray(new_seq,char);
    } except_when(BADSEQ) { unarray(new_seq,char); }
}


#define LOADED_DATA \
"data files '%s' and '%s' are %sloaded.\n(%d %s progeny, %d loci, %d trait%s)\n"
#define ALREADY_LOADED_DATA \
  "data are already loaded\nshut down and restart the program to load new data"
#define NO_LOADED_DATA \
  "no data are loaded\ntype \"load data <filename>\" to load data\n"
  
command load_data() 
{
    FILE *fpa=NULL, *fpb=NULL, *fpc=NULL;
    int num_of_file=0;
    char *dfile= get_temp_string(), *tfile= get_temp_string(), *mfile= get_temp_string();

    if (nullstr(args)) {
	if (!data_loaded()) print(NO_LOADED_DATA); else {
	    strcpy(tfile,raw.file);
	    make_filename(tfile,FORCE_EXTENSION,TRAIT_EXT);
	    sprintf(ps, LOADED_DATA, raw.file, tfile, "", raw.n_indivs,
                (raw.data_type==BACKCROSS ? "backcross":"intercross"),
                raw.n_loci, raw.n_traits, maybe_s(raw.n_traits)); pr();
	}
	
    } else { /* !nullstr */
	
	nstoken(&args,sREQUIRED,dfile,PATH_LENGTH); tfile[0]='\0';
	nomore_args(num_args);
	run {
	    if (!make_filename(dfile,FORCE_EXTENSION,DATA_EXT)) send(CANTOPEN);
	    fpa= open_file(dfile,READ);
	    fgetln(fpa);

	    if (streq(ln,"prepared data f2 backcross")) {
		strcpy(geno_chars,default_backcross_chars);
		raw.data_type=BACKCROSS;

	    } else if (streq(ln,"prepared data f2 intercross")) {
		strcpy(geno_chars,default_intercross_chars);
		raw.data_type=INTERCROSS; raw.f3=FALSE;

	    } else if (streq(ln,"prepared data f3")) {
		strcpy(geno_chars,default_intercross_chars);
		raw.data_type=INTERCROSS; raw.f3=FALSE; /* KLUDGE!!! */
		
	    } else { 
		error("unrecognized header line in data file");
	    }

	    strcpy(tfile,dfile);
	    if (!make_filename(tfile,FORCE_EXTENSION,TRAIT_EXT))
	      send(CANTOPEN);
	    fpb = open_file(tfile,READ);
	    strcpy(mfile,dfile);
	    if (!make_filename(mfile,FORCE_EXTENSION,MAPS_EXT))
	      send(CANTOPEN);
	    fpc = open_file(mfile,READ);

	    read_data(fpa,fpb,fpc,dfile,num_of_file);

	    crunch_data();
	    allocate_qtl_struct(raw.max_traits*2,raw.max_traits*2);

	    sprintf(ps, LOADED_DATA, dfile, tfile, "", raw.n_indivs,
                (raw.data_type==BACKCROSS ? "backcross":"intercross"),
	       raw.n_loci-dum_loc, raw.n_traits, maybe_s(raw.n_traits)); pr();
	    if (raw.n_traits==1) trait=0; else trait= NOTRAIT;
	    update_top();
	    if(!altered_chroms)
	      load_qtl_files();

	} on_exit {
	    if (msg==CANTOPEN) {
		sprintf(ps, "error: unable to open data file\n"); pr();
	    } else if (msg==BADDATA) {
		sprintf(ps, "error: unable to load data from file\nline %d:",
                BADDATA_line_num); pr();
		print(BADDATA_error); nl();
		strcpy(raw.file,"");
	    } 
	    if (msg!=CANTOPEN && msg!=0) { 
		print("data file not loaded\n"); 
		strcpy(raw.file,"");
		close_file(fpa); close_file(fpb); close_file(fpc);
		relay_messages;
	    } 
	}
    }
}

/********** NO MORE PREP IN QTL


command prep_data()
{
    char *in_name,*out_name,*str;
    FILE *fp;
    int length;

    in_name = NULL; out_name = NULL;
    
    run {
	str = get_temp_string();
	if(!stoken(&args,sREQUIRED,str)) 
	  input("Input file to be prepared: ",str,25);
	in_name = str;

	length = len(str);
	if(streq(&str[length-5],DATA_EXT))
	  error("raw file cannot have '.data' extension");
	if(streq(&str[length-7],TRAIT_EXT))
	  error("raw file cannot have '.trait' extension");

	out_name = get_temp_string();
	strcpy(out_name,in_name);
	make_filename(out_name,FORCE_EXTENSION,DATA_EXT);

 fp = open_file(in_name,READ);
	getdataln(fp);  crunch(ln);
	close_file(fp);
	if(streq(ln,"data type f2 intercross") ||
	   streq(ln,"data type f3 intercross") ||
	   streq(ln,"data type f2 backcross") ||
	   streq(ln,"data type ri sib") ||
	   streq(ln,"data type ri self") ||
	   streq(ln,"data type f2") ||
	   streq(ln,"data type backcross") ||
	   streq(ln,"data type quant intercross") ||
	   streq(ln,"data type quant backcross") ) {
	    f2_prep(in_name,out_name);
	}
	else {
	    sf(ps,"'%s' is not an allowable data type",ln);
	    error(ps);
	}
    } except_when(CANTOPEN) {
	sf(ps,"Can't open file '%s'.",in_name);
	error(ps);
	
    }
}

**********/

command save_status()
{
    char *name, *name2, *name3;
    FILE *fp;

    if(nullstr(raw.file)) {
	print("No data have been loaded, none can be saved.\n");  return; }
    name = mkstrcpy(raw.file);
    name3 = mkstrcpy(raw.file);
    name2 = "dummy.qtls";
    make_filename(name, FORCE_EXTENSION, QTL_EXT);

    /* The following accounts for the case where there is no file 'name'
       (i.e. it makes one)  */
    run { 
	fp = open_file(name,READ);
	close_file(fp);
    } except_when(CANTOPEN) {
	fp = open_file(name,WRITE);
	close_file(fp);
    }

    make_filename(name3, FORCE_EXTENSION, QTL_OLD);
    run {
	fp = open_file(name2, WRITE);
	sprintf(ps, "Now saving %s...\n", name);  pr();
	save_qtl_files(fp);
	close_file(fp);
	if (rename_file(name,name3)) rename_file(name2,name);
    } except {
	when CANTOPEN:
	    sprintf(ps, "Can't open %s.\n", name);  pr(); /* fall through */
	default:
	    rename_file(name3,name);
	    if(msg == INTERRUPT) send(INTERRUPT);
	    print("\nAn error occurred while the data files were being rewritten.\n");
	 print("You may have to turn the 'auto save' option 'off' before quitting.\n");
	    if(redirecting_input) send(QUIT);
	    
    }

    name2="dummy.traits";
    make_filename(name, FORCE_EXTENSION, TRAIT_EXT);
    make_filename(name3, FORCE_EXTENSION, TRAIT_OLD);
    run {
	fp = open_file(name2, WRITE);
	sprintf(ps, "Now saving %s...\n", name);  pr();
	save_traitfile(fp);
	close_file(fp);
	if (rename_file(name,name3)) rename_file(name2,name);
    } except {
	when CANTOPEN:
	    sprintf(ps, "Can't open fle '%s'.\n", name);  pr(); /* fall through */
	default:
	    rename_file(name3,name);
	    if(msg == INTERRUPT) send(INTERRUPT);
	    print("\nAn error occurred while the data files were being rewritten.\n");
	 print("You may have to turn the 'auto save' option 'off' before quitting.\n");
	    if (redirecting_input) send(QUIT);
    }
}


bool save_on_exit(do_it_now)
bool do_it_now;
{
    if(!do_it_now) return(auto_save && data_loaded());
    if(auto_save && data_loaded()) save_status();
    return(TRUE);
}

void save_qtl_files(fp)
FILE *fp;
{
    int i;
    
    fprintf(fp,"%d\n",raw.filenumber);
    fprintf(fp,"%d %d\n",num_wiggles,first_wiggle);
    for(i = 0; i < num_wiggles; i++) {
	save_wiggle(fp,i);
    }
    
    fprint(fp,"#Compares\n");
    
    fprintf(fp,"%d %d\n",num_compares,first_compare);
    for(i = 0; i < num_compares; i++) {
	save_compare(fp,i);
    }
}

#define QTL_LD_ERROR1 \
"An error ocurred while loading the qtl map data in file '%s'.\n"
#define QTL_LD_ERROR2 \
"In some extreme cases, this can prevent correct operation of this program.\n"
#define QTL_LD_ERROR3 \
"If problems occur, delete this file and try reloading the data set.\n"
#define QTL_LD_OLD1 \
"You have prepared a new version of your data files since the last qtl map\n"
#define QTL_LD_OLD2 \
"data were saved in the '.qtls' file. These data will not be loaded.\n" 

void load_qtl_files()
{
    char name[PATH_LENGTH+1];
    FILE *fp=NULL;
    int i, n_wigs, n_comps, filenum;

    run {
	strcpy(name,raw.file);
	make_filename(name,FORCE_EXTENSION,QTL_EXT);
	fp= open_file(name,READ);

	if(fscanf(fp,"%d\n",&filenum) != 1) send(IOERROR); /* KLUDGE */
	else if(filenum != raw.filenumber) {
	    print(QTL_LD_OLD1); print(QTL_LD_OLD2);
	}
	else if(fscanf(fp,"%d %d\n",&n_wigs,&first_wiggle)!= 2) send(IOERROR);

	else {
	    for (i = 0; i < n_wigs; i++) {
		load_wiggle(fp);
	    }

	    fgetdataln(fp,NULL);
	    
	    if (sscanf(ln,"%d %d\n",&n_comps,&first_compare) != 2)
	      send(IOERROR);
    
	    for(i = 0; i < n_comps; i++) {
		load_compare(fp);
	    }
	
	    sprintf(ps, "QTL map data in file '%s' have been loaded.\n", name);
	    pr();
	}
	close_file(fp);

    } except {
	when CANTOPEN: 
	  print("Unable to load any saved QTL map data.\n");
	  msg=0; break; /* No .qtls file - WHY IS THIS MSG=0 NEEDED??? */
	default:
	  sprintf(ps, QTL_LD_ERROR1, name); pr();
	  print(QTL_LD_ERROR2); print(QTL_LD_ERROR3);
	  close_file(fp);
    }
}



command tester()
{
    real theta, f2_sum, f3_sum, a, b, c, x, y, z, left_rf, right_rf;
    int qtl, left, right;

    getln("Theta, L_Pos, R_pos: ");
    sscanf(ln,"%lf %lf %lf",&theta,&left_rf,&right_rf);
    /* left_rf= rmaxf(MIN_REC_FRAC,min(pos,MAX_FRAC_OF_RF*theta));
    right_rf= (theta - left_rf)/(1 - 2*left_rf); */

    sprintf(ps, "left_rf=%lf right_rf=%lf\n", left_rf, right_rf); pr();
    for (left=0; left<4; left++) 
      for (right=0; right<4; right++) {	
	f2_sum=f3_sum= 0.0;
	x=y=z= 0.0;
	raw.f3=TRUE; c=transition_prob(INTERCROSS,left,right,theta);
	/* raw.f3=FALSE; z=transition_prob(INTERCROSS,left,right,theta); */
	sprintf(ps, "left=%d  right=%d  F2-prob=%lf  F3-prob=%lf\n", left, right, z, c); pr();
	for (qtl=0; qtl<4; qtl++) {
	    raw.f3=TRUE; f3_sum+= 
	      ((a=transition_prob(INTERCROSS,left,qtl,left_rf)) *
	       (b=transition_prob(INTERCROSS,qtl,right,right_rf))) /c;
	    /* raw.f3=FALSE; f2_sum+= 
	      ((x=transition_prob(INTERCROSS,left,qtl,left_rf)) *
	       (y=transition_prob(INTERCROSS,qtl,right,right_rf))) /z; */
	    sprintf(ps, "qtl= %d  F2: L=%lf R=%lf  F3: L=%lf R=%lf\n", qtl, x, y, a, b); pr();
	}		
	sprintf(ps, "==== F2 sum= %lf  F3 sum=%lf\n", f2_sum, f3_sum); pr();
    }
}


command tweak_weight()
{
    real start,end,step,weight;
    int perm;

    qtl_ready(ANY_DATA,SEQ,TRAIT,QCTM);
    get_arg(rtoken,0.0,&start);
    get_arg(rtoken,1.0,&end);
    get_arg(rtoken,0.1,&step);
    nomore_args(num_args);

   run {
	fix_weight_kludge = TRUE;
	for_all_orders(ints,map,perm) {
	    nl(); print_wiggle_title(); 
	    print_wiggle_interval(map); 
	    for (weight=start; weight<=end; weight+=step) {
		map->qtl_weight[0]= weight;
		make_qtl_map(map);
		print_wiggle_map(map,2.0,0.25); 
	    }
	    print_wiggle_interval(NULL); 
	}
    } on_exit { fix_weight_kludge = FALSE; }
    fix_weight_kludge = FALSE; 
}


command draw_wiggle()
{
    char arg[TOKLEN+1];
    int wiggle, order;
    real threshold, scale;

    qtl_ready(ANY_DATA,NOSEQ,NOTRAIT,NOQCTM);
    get_arg(stoken,"",arg);
    get_arg(rtoken,2.0,&threshold);
    get_arg(rtoken,0.25,&scale);
    nomore_args(num_args);

    if (threshold<0.0 || scale<=0.0) usage_error(2);
    get_wiggle_nums(arg,&wiggle,&order);
    if (wiggle<0) { wiggle=num_wiggles-1; order=wiggles[wiggle]->num_orders-1;}
    else if (order<0 && wiggles[wiggle]->num_orders==1) order=0;
    nl();

    if (order<0) { 
        print_ps_multi_wiggle(wiggle, threshold);
    } else {
	print_ps_wiggle_order(wiggle, order, threshold);
    }
}
    

