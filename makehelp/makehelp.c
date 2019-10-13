/******************************************************************************

 #    #    ##    #    #  ######  #    #  ######  #       #####            ####
 ##  ##   #  #   #   #   #       #    #  #       #       #    #          #    #
 # ## #  #    #  ####    #####   ######  #####   #       #    #          #
 #    #  ######  #  #    #       #    #  #       #       #####    ###    #
 #    #  #    #  #   #   #       #    #  #       #       #        ###    #    #
 #    #  #    #  #    #  ######  #    #  ######  ######  #        ###     ####

POSTSCRIPT VERSION

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

/* Makes help files (.code and .help) from a template. This is a
   quick-and-dirty program with lots of globals and other ickiness. Sorry. */

/* The syntax of the input file is as follows in this example:

   @title MY MANUAL
   @topic Some Commands
   $Here is a summary of these commands

   @cmd do this (dt) 	name needs to be uniquely parsable by shell, options
			are  @cmd, @opt, @param, @info, @topic, or @end
   $one line text	
   <1 			<n UPTO n args; =n EXACTLY n args; >n at least n args
   %<arg1> <arg2>	
   ^default1 default2	default values

   Text describing command - one blank line is REQUIRED?

   @cmd next command 
   # comment line ignored

   !system output
   !more system output
   ...

   @end 
*/

#include <stdio.h>
#include <string.h>

//#define INC_LIB
//#define INC_SHELL   /* for various definitions */
#include "system.h"

typedef char STRING[MAXLINE+1];

/*** stuff for PostScript formatted manual ***/
#define FONT_SIZE   10
#define LEFT_MARGIN 80   /* 1.12 inch */
#define TOP_MARGIN  62   /* = page ht - margin  - font size*PAGE_LINES */
#define FOOT_MARGIN 50
/*********************************************/

#ifndef _TEXT_MANUAL
#define PAGE_LINES  68   /* NOT including FOOTER */
#define FOOTER "%-30s    %20s                    Page %s\n"
#define MAN_TAB 36
#else
#define PAGE_LINES  60
#define FOOTER "    %-30s %20s              Page %s\n"
#define MAN_TAB "    "
#endif

#define PAGE_NUM "%d"
/*1234567890123456789012345678901234567*9012345678901234567890123456789012345*/
/*    MAPMAKER/EXP V3.0 Manual  12345678901234567890 Commands  ------Page xx*/
#define CONTENTS_LEFT    70

#define LINE_BREAK_LEN   1l

#define SPACE "    "
#define HLP_TAB "    "


/**** other defs ****/

static FILE *file=NULL, *code=NULL, *hlp=NULL;
static STRING file_name, code_name, hlp_name, final_hlp_name, code_failed;
static long pos;
static int n_args, prefix, topic;
static char **entry;
static int entry_type[MAX_COMMANDS+1];
static STRING cmd_description[MAX_COMMANDS+1];
static STRING str, type, name, arguments, defaults, description, title;
static STRING abbreviation, sequence, results, section[MAX_COM_TOPICS];
static long position[MAX_COM_TOPICS];

/**** defs for man_stuff ****/
static int lines, page=0, start_page=0, ps_page=0;
static int entries, pending, entry_page[MAX_COMMANDS+1];
static FILE *man=NULL;
static STRING save, man_name, chapter_title;
#define make_man (!nullstr(man_name))

static void nextstr(void);
static void parse_error(char *msg, int punt /* 0= no, 1= this-entry, 2= entirely */);
static void close_files(char *name);
static void parse_entry(int kind, char *name, char *abbrev);
static void write_mkhelp(char *cmd,char *abbrev,long pos, int prefix, int num_args, char *desc, char *args, char *defs,int topic,int kind);
static void write_topics_and_end(void);
static void man_new_page(void);
static void man_write_line(char *line);
static void man_write_done(void);
static void man_new_entry(int kind, char *name, bool has_description) /* only CMD or OPT */;
static void man_new_topic(char *name, char *description, bool has_description);
static void man_write_contents(void);
static void man_write_title(void);
static void ps_file_start(FILE *fp);
static void ps_file_end(FILE *fp);
static void ps_page_start(FILE *fp, int pagenum);
static void ps_page_end(FILE *fp);
static char *ps_string(char *str);


void nextstr(void)
{
    do finput(file,str,MAXLINE); /* crunches */
    while (str[0]=='#');
}


void parse_error(char *msg, int punt /* 0= no, 1= this-entry, 2= entirely */)
{
    sprintf(ps, "%s  ", msg); pr();
    if (punt>0) do nextstr(); while (str[0]!='@');
    if (punt==2) send(QUIT);
}


void close_files(char *name)
{
    close_file(file);

    do_fwrite(hlp,"@end\n");
    close_file(hlp);

    write_topics_and_end();
    close_file(code);

    man_write_contents();

#ifndef _TEXT_MANUAL
    /* POSTSCRIPT postscript */
    ps_page_end(man);
    ps_file_end(man);
#endif

    close_file(man);

    sprintf(ps, "%s done\n", name); pr();
}


void parse_entry(int kind, char *name, char *abbrev)
{
    bool rest=FALSE;
    int i;
    STRING line, key;

    /* get stuff for cmd/opt */
    description[0]= arguments[0]= defaults[0]= '\0';
    n_args=-1; prefix=0;

    while (!rest) {
	switch(str[0]) {

	    case '@':
	      rest=TRUE; break; /* the while loop */

	    case '$': 		/* $ - one line description */
	      if (len(str+1)>MAX_COM_HELP_LEN)
	        parse_error("description too long",0);
	      nstrcpy(description,str+1,MAX_COM_HELP_LEN);
	      nextstr(); break;
	    
	    case '<':       /* <,=,> */
	    case '=': 
	    case '>':
	      if (white(str[1]) || sscanf(str+1,"%d",&n_args)!=1 ||
		  n_args>9 || n_args<0) parse_error("bad digit",0);
	      else {
		  if (str[0]=='<') prefix= UPTO;
		  else if (str[0]=='=') prefix= EXACTLY;
		  else if (str[0]=='<') prefix= ATLEAST;
		  for (i=2; str[i]!='\0'; i++) if (!white(str[i])) break;
		  if (len(str+i)>MAX_ARG_HELP_LEN)
		    parse_error("arguments too long",0);
		  else strcpy(arguments,str+i);
	      }
	      nextstr(); break;

	    case '%':  		/* % - defaults */
	      if (len(str+1)>MAX_ARG_HELP_LEN)
	        parse_error("defaults too long",0);
	      strcpy(defaults,str+1);
	      nextstr(); break;

#ifdef UNUSED_THIS
	    case '&': 		/* & - sequence */
	      if (len(str+1)>MAX_ARG_HELP_LEN)
	        parse_error("sequence too long",0);
	      strcpy(sequence,str+1);
	      nextstr(); break;
#endif
	    case '*': 		/* * - results */
	      if (len(str+1)>MAX_ARG_HELP_LEN)
	        parse_error("results too long",0);
	      strcpy(results,str+1);
	      nextstr(); break;
	    
	    case '#': 		/* # comment - ignore */
	      nextstr(); break;

	    default:
	      rest=TRUE; break;
	  }
    }
    /* now rest=TRUE, meaning there is more text */
    while (nullstr(str)) nextstr();
    if (str[0]=='@') { parse_error("null help text",0); rest=FALSE; }
    strcpy(key,name); crunch(key);

    /* check for bogasity - WANT MORE? */
    if (n_args>0  &&  nullstr(arguments)) parse_error("no args text",0);
    if (n_args==0 && !nullstr(arguments)) parse_error("extra args text",0);
    if (n_args==0 && !nullstr(defaults)) parse_error("extra default text",0);
    if (kind==TOP && (!nullstr(arguments) || !nullstr(defaults) || 
		      !nullstr(results) || !nullstr(sequence) || n_args>0))
      parse_error("improper entries for topic",0);
    if (nullstr(description)) parse_error("no description",0);

    /* save info and write help file key */
    strcpy(entry[entries],name);
    entry_type[entries]=kind;

    if (kind==TOP) {
	sprintf(ps, "%cTOPIC %s\n", '@', key); do_fwrite(hlp, ps);
	pos+= (long)(len(name)+1+6)+LINE_BREAK_LEN; /* after writing above! */
	topic++; strcpy(section[topic],description); position[topic]=pos;
    } else {
	/* write code, doc file entries */
	sprintf(ps, "%c%s\n", '@', key); do_fwrite(hlp, ps);
	pos+= (long)(len(name)+1)+LINE_BREAK_LEN; /* after writing above! */
	write_mkhelp(key,abbrev,pos,prefix,n_args,description,arguments,
		     defaults,topic,kind);
	strcpy(cmd_description[entries],description);
    }
    if (kind==TOP) man_new_topic(name,description,str[0]!='@');
      else man_new_entry(kind,name,str[0]!='@');

    /* write long text, if any */
    if (rest) 
      while (str[0]!='@') {
	  if (str[0]=='!') sprintf(line,"%s%s",HLP_TAB,str+1);
	    else if (str[0]=='^') { strcpy(line,str+1); uppercase(line); }
	    else strcpy(line,str);
	  nextstr();
	  if (!(nullstr(line) && str[0]=='@')) { /* punt last line if blank */
	      do_fwrite(hlp,line); fnl(hlp);
	      pos+= (long)len(line) + LINE_BREAK_LEN;
	      man_write_line(line);
	  }
      }
    man_write_done();
    entries++;
}


void write_mkhelp(char *cmd,char *abbrev,long pos, int prefix, int num_args, char *desc, char *args, char *defs,int topic,int kind)
{
    STRING temp;

    switch (kind) {
	case TOP: strcpy(temp,"TOP"); break;
	case CMD: strcpy(temp,"CMD"); break;
	case OPT: strcpy(temp,"OPT"); break;
	case PAR: strcpy(temp,"PAR"); break;
	case HLP: strcpy(temp,"HLP"); break;
    }
    sprintf(ps, " mkhelp(\"%s\",\"%s\",%ldl,%s,%d,%s,%d,\n",
            cmd, abbrev, pos,
            (prefix==EXACTLY ? "EXACTLY":(prefix==UPTO ? "UPTO":"ATLEAST")),
            num_args, temp, topic); fpr(code);
    sprintf(ps, "        \"%s\",\n", desc); fpr(code);
    sprintf(ps, "        \"%s\",\n", args); fpr(code);
    sprintf(ps, "        \"%s\");\n", defs); fpr(code);
    /* need to add sequence, results, etc */
}


void write_topics_and_end(void)
{
    int i, s;
    STRING temp;
    if (!make_man) return;

    fnl(code);
    for (i=0, s=1; i<entries; i++) if (entry_type[i]==TOP) {
	strcpy(temp,section[s]); /* uppercase(temp); */
	sprintf(ps, " mktopic(%d,\"%s\",TOP,%ldl);\n", s, temp, position[s]);
	fpr(code); s++;
    }
    fprint(code,"}\n");
}


	    
/*****************************************************************************/


void man_new_page(void)
{
    STRING p, temp, foot;
    if (!make_man) return;

#ifndef _TEXT_MANUAL
    /******* POSTSCRIPT VERSION *******/
    if (page==0) 
      { page=start_page; lines=PAGE_LINES; 
	ps_page=1; ps_page_start(man,ps_page); return; }

    sprintf(p, PAGE_NUM, abs(page)); if (page > 0) page++; else page--;
    ps_page++;
    sprintf(temp, "%s", title);
    sprintf(foot, FOOTER, temp, chapter_title, p);
    
    fprintf(man,"GS %d %d moveto /Courier FF 9 SF F (%s)S GR\n",
	    LEFT_MARGIN,FOOT_MARGIN,foot);
    lines=PAGE_LINES;
    
    ps_page_end(man);

    ps_page_start(man,ps_page);

#else
    /********** ASCII VERSION *********/
    if (page==0) 
      { fprintf(man,"\n"); page=start_page; lines=PAGE_LINES; return; }

    while (lines>0) { fprintf(man,"\n"); lines--; }
    sprintf(p,PAGE_NUM,abs(page)); if (page>0) page++; else page--;
    fprintf(man,"\n\n\n");
    sprintf(temp,"%s",title);
    fprintf(man,FOOTER,temp,chapter_title,p);
    lines= PAGE_LINES;
    fprintf(man,"\n");
#endif

}


void man_write_line(char *line)
{
    if (!make_man) return;

#ifndef _TEXT_MANUAL
    /******** POSTSCRIPT VERSION ********/
    if (pending) {
        if (lines <= 1) man_new_page();
	if (save[0]=='!') { 
	  fprintf(man,"GS %d %d moveto /Courier FF 10 SF F (%s)S GR\n",
		  LEFT_MARGIN+MAN_TAB, (lines*FONT_SIZE)+TOP_MARGIN, 
		  ps_string(save+1));
        } else {
	  fprintf(man,"GS %d %d moveto /Courier FF 10 SF F (%s)S GR\n",
		  LEFT_MARGIN, (lines*FONT_SIZE)+TOP_MARGIN, ps_string(save));
	}
	lines--;
    }
#else
    /********** ASCII VERSION **********/
    if (pending) {
	if (lines<=1) man_new_page();
	fprintf(man,SPACE); 
	if (save[0]=='!') { fprintf(man,MAN_TAB); fprintf(man,save+1); }
	else fprintf(man,save);
	fprintf(man,"\n"); lines--;
    }
#endif

    strcpy(save,line); pending=TRUE; 
}


void man_write_done(void)
{
    if (!make_man) return;

#ifndef _TEXT_MANUAL
    /******** POSTSCRIPT VERSION ********/
    if (save[0]=='!') { 
        fprintf(man,"GS %d %d moveto /Courier FF 10 SF F (%s)S GR\n",
		LEFT_MARGIN+MAN_TAB, (lines*FONT_SIZE)+TOP_MARGIN, 
		ps_string(save+1));
    } else {
        fprintf(man,"GS %d %d moveto /Courier FF 10 SF F (%s)S GR\n",
		LEFT_MARGIN, (lines*FONT_SIZE)+TOP_MARGIN, ps_string(save));
    }
    lines--; pending=FALSE;

#else
    /********** ASCII VERSION **********/
    fprintf(man,SPACE); 
    if (save[0]=='!') { fprintf(man,MAN_TAB); fprintf(man,save+1); }
    else fprintf(man,save);
    fprintf(man,"\n"); lines--;
    pending=FALSE;
#endif

}


void man_new_entry(int kind, char *name, bool has_description) /* only CMD or OPT */
{
    STRING temp0, temp1, temp2, templine;
    int n= n_args; /* -1 => plural */
    int need=5;  /* name + space + args + defs + (space or no_desc line) */
    bool blank=FALSE;
    if (!make_man) return;

    if (n>0 && !nullstr(arguments) && !nullstr(defaults)) need++;
    if (has_description) need+=2;
    if (lines!=PAGE_LINES) need+=3; /* if not 1st on page */

    if (need>lines) { man_new_page(); }
    else {
      lines-=3;
#ifdef _TEXT_MANUAL
      fprintf(man,"\n\n\n");
#endif
    }
    entry_page[entries]= page; 

    if (kind==CMD) strcpy(temp1,"Command");
    else if (kind==OPT) strcpy(temp1,"Option");
    else if (kind==PAR) strcpy(temp1,"Parameter");
    else if (kind==HLP) strcpy(temp1,"Information");
    else strcpy(temp1,"???");
    strcpy(temp0,name); uppercase(temp0);
    if (nullstr(abbreviation)) temp2[0]='\0';
    else sprintf(temp2, " (abbreviation '%s')", abbreviation);

#ifndef _TEXT_MANUAL
    sprintf(templine, "%s %s%s", temp0, temp1, temp2);
    fprintf(man,"GS %d %d moveto /Courier-Bold FF 10 SF F (%s)S GR\n",
	    LEFT_MARGIN,(lines*FONT_SIZE)+TOP_MARGIN,ps_string(templine));
    fprintf(man,"GS .5 LW %d %d moveto /Courier-Bold FF 10 SF F (%s)US GR\n",
	    LEFT_MARGIN,(lines*FONT_SIZE)+TOP_MARGIN-1,ps_string(templine));
    lines-=2;

    if (!nullstr(description)) {
        sprintf(templine, "Summary:     %s", description);
        fprintf(man,"GS %d %d moveto /Courier FF 10 SF F (%s)S GR\n",
		LEFT_MARGIN,(lines*FONT_SIZE)+TOP_MARGIN,
		ps_string(templine));
	lines--;
	blank=TRUE;
    }

    if (n==0) { 
        fprintf(man,"GS %d %d moveto /Courier FF 10 SF F (%s)S GR\n",
		LEFT_MARGIN,(lines*FONT_SIZE)+TOP_MARGIN,"No Arguments");
	lines--; blank=TRUE; 
    }
    else if (!nullstr(arguments)) {
        sprintf(templine, "Argument%s:%s   %s", maybe_s(n), maybe_sp(n), arguments);
	fprintf(man,"GS %d %d moveto /Courier FF 10 SF F (%s)S GR\n",
		LEFT_MARGIN,(lines*FONT_SIZE)+TOP_MARGIN,ps_string(templine));
	lines--;
	blank=TRUE;

	if (!nullstr(defaults)) {
	  sprintf(templine, "Default%s:%s    %s", maybe_s(n), maybe_sp(n), defaults);
	  fprintf(man,"GS %d %d moveto /Courier FF 10 SF F (%s)S GR\n",
		  LEFT_MARGIN,(lines*FONT_SIZE)+TOP_MARGIN,
		  ps_string(templine));
	  lines--;
	}
    }
    if (blank) lines--;

#else
    fprintf(man,"%s%s %s%s\n",SPACE,temp0,temp1,temp2); lines-=1;
    fprintf(man,"\n",temp0); lines-=1;

    if (!nullstr(description)) {
	fprintf(man,"%sDescription: %s\n",SPACE,description); lines--; 
	blank=TRUE;
    }

    if (n==0) { fprintf(man,"%sNo Arguments\n",SPACE); lines--; blank=TRUE; }
    else if (!nullstr(arguments)) {
	fprintf(man,"%sArgument%s:%s   %s\n",SPACE,maybe_s(n),maybe_sp(n),
	       arguments); lines--;
	blank=TRUE;
	if (!nullstr(defaults)) {
	    fprintf(man,"%sDefault%s:%s    %s\n",SPACE,maybe_s(n),maybe_sp(n),
		   defaults); lines--;
	}
    }
    if (blank) { fprintf(man,"\n"); lines--; }
#endif

}


void man_new_topic(char *name, char *description, bool has_description)
{
    STRING upcase,templine;

    if (!make_man) return;

    man_new_page();
    strcpy(chapter_title,name);
    entry_page[entries]= page;

    strcpy(upcase,description); uppercase(upcase);

#ifndef _TEXT_MANUAL
    sprintf(templine, "(%d) %s", topic, upcase);
    fprintf(man,"GS %d %d moveto /Courier-Bold FF 10 SF F (%s)S GR\n",
	    LEFT_MARGIN,(lines*FONT_SIZE)+TOP_MARGIN,ps_string(templine));
    fprintf(man,"GS .5 LW %d %d moveto /Courier-Bold FF 10 SF F (%s)US GR\n",
	    LEFT_MARGIN,(lines*FONT_SIZE)+TOP_MARGIN-1,ps_string(templine));
    lines--;
    if (has_description) { lines--; }
#else
    fprintf(man,"%s(%d) %s\n",SPACE,topic,upcase); lines--;
    if (has_description) { fprintf(man,"\n"); lines--; }
#endif
}


void man_write_contents(void)
{
    int i, s, k;
    STRING temp, upcase, templine;
    if (!make_man) return;

    man_new_page();
    strcpy(upcase,title); uppercase(upcase);

#ifndef _TEXT_MANUAL
    sprintf(templine, "%s COMMAND REFERENCE:", upcase);
    fprintf(man,"GS %d %d moveto /Courier-Bold FF 10 SF F (%s)S GR\n",
	    LEFT_MARGIN,(lines*FONT_SIZE)+TOP_MARGIN,ps_string(templine));
#else
    fprintf(man,"%s%s COMMAND REFERENCE:\n\n",SPACE,upcase);
#endif

    lines-=2;
    strcpy(chapter_title,"Command Reference"); page= -2;
    for (i=0, s=1; i<entries; i++) {
	switch(entry_type[i]) {
	    case TOP:
	      man_write_line("");
	      strcpy(upcase,section[s]); uppercase(upcase);
	      sprintf(temp, " %s\\(%d\\) %s ", (s >= 10 ? "" : " "), s, upcase); s++; break;
	    case CMD: 
	      strcpy(upcase,entry[i]); uppercase(upcase);
	      sprintf(temp, "      %s Command ", upcase); break;
	    case OPT:
	      strcpy(upcase,entry[i]); uppercase(upcase);
	      sprintf(temp, "      %s Option ", upcase); break;
	    case PAR:
	      strcpy(upcase,entry[i]); uppercase(upcase);
	      sprintf(temp, "      %s Parameter ", upcase); break;
	    case HLP:
	      strcpy(upcase,entry[i]); uppercase(upcase);
	      sprintf(temp, "      %s Information ", upcase); break;
	}
	for (k=len(temp); k<CONTENTS_LEFT; k++) temp[k]='.';
	sprintf(temp+k," %d",entry_page[i]);
	man_write_line(temp);
    }

    man_write_line("");
    man_new_page();
    strcpy(upcase,title); uppercase(upcase);

#ifndef _TEXT_MANUAL
    sprintf(templine, "%s QUICK REFERENCE:", upcase);
    fprintf(man,"GS %d %d moveto /Courier-Bold FF 10 SF F (%s)S GR\n",
	    LEFT_MARGIN,(lines*FONT_SIZE)+TOP_MARGIN,ps_string(templine));
    lines-=1;
#else
    fprintf(man,"%s%s QUICK REFERENCE:\n\n",SPACE,upcase);
    lines-=2;
#endif
    strcpy(chapter_title,"Quick Reference");
    for (i=0, s=1; i<entries; i++) {
	if (entry_type[i]==TOP) {
	    man_write_line("");
	    strcpy(upcase,section[s]); uppercase(upcase);
	    sprintf(temp, "%s\\(%d\\) %s ", (s >= 10 ? "" : ""), s, upcase); s++;
	    man_write_line(temp);
	} else {
	    sprintf(temp, "%s%s", entry[i], (entry_type[i] == HLP ? "*" : ""));
	    for (k=len(temp); k<MAX_COM_NAME_LEN+1; k++) temp[k]='.';
	    sprintf(temp+k,"%s",cmd_description[i]);
	    man_write_line(temp);
	}
    }
    man_write_line("");
    man_write_line("");
    sprintf(temp, "* = reference information only - not a command");
    man_write_line(temp);
    man_write_line("");
    man_new_page();
    man_write_done();
}


void man_write_title(void)
{
    if (!make_man) return;

#ifndef _TEXT_MANUAL
    ps_file_start(man);
#endif
}

/********** POSTSCRIPT STUFF **********/

void ps_file_start(FILE *fp)
{
    fprintf(fp,"%%!PS-Adobe-3.0\n");
    fprintf(fp,"%%%%Creator: MAPMAKER\n");
    fprintf(fp,"%%%%LanguageLevel: 1\n");
    fprintf(fp,"%%%%PageOrder: Special\n");

    fprintf(fp,"%%%%BeginProlog\n");
    fprintf(fp,"%%%%BeginResource: procset Map_Painter_prolog\n");
    fprintf(fp,"/Map_Painter_prolog 100 dict def\n");
    fprintf(fp,"Map_Painter_prolog begin\n");
    fprintf(fp,"/CF {dup 0 eq\n");
    fprintf(fp,"       {pop /Times-Bold}\n");
    fprintf(fp,"       {1 eq {/Times-Roman} {/Times-Italic} ifelse}\n");
    fprintf(fp,"     ifelse} def\n");
    fprintf(fp,"/F {setfont} def\n");
    fprintf(fp,"/FF {findfont} def\n");
    fprintf(fp,"/GM {restore} def\n");
    fprintf(fp,"/GR {grestore} def\n");
    fprintf(fp,"/GS {gsave} def\n");
    fprintf(fp,"/LC {currentpoint pop dup /l exch def add /r exch def\n");
    fprintf(fp,"     counttomark 3 idiv {\n");
    fprintf(fp,"       currentpoint pop r ge {l currentpoint exch pop LLD sub moveto} if\n");
    fprintf(fp,"       LD\n");
    fprintf(fp,"     } repeat pop} def\n");
    fprintf(fp,"/LD {/f exch def /n exch def /p exch def\n");
    fprintf(fp,"     p length 0 ne {\n");
    fprintf(fp,"       /Times-Roman findfont LFS 1 sub scalefont setfont p show\n");
    fprintf(fp,"       2.5 0 rmoveto\n");
    fprintf(fp,"     } if\n");
    fprintf(fp,"     f CF findfont LFS scalefont setfont\n");
    fprintf(fp,"     n show 6 0 rmoveto} def\n");
    fprintf(fp,"/LM {save} def\n");
    fprintf(fp,"/LW {setlinewidth} def\n");
    fprintf(fp,"/S {show} def\n");
    fprintf(fp,"/SF {scalefont} def\n");
    fprintf(fp,"/SR {dup stringwidth pop neg 0 rmoveto show} def\n");
    fprintf(fp,"/TR {translate} def\n");
    fprintf(fp,"/XY {moveto} def\n");
    fprintf(fp,"/US {stringwidth pop 0 rlineto stroke} def\n");
    fprintf(fp,"end\n");
    fprintf(fp,"%%%%EndResource\n");
    fprintf(fp,"%%%%EndProlog\n");

    fprintf(fp,"/LFS 8 def /LLD 7.5 def\n");
    fprintf(fp,"%%%%BeginSetup\n");
    fprintf(fp,"Map_Painter_prolog begin\n");
    fprintf(fp,"gsave\n");
    fprintf(fp,"%%%%EndSetup\n");
}

void ps_file_end(FILE *fp)
{
    fprintf(fp,"%%%%Trailer\n");
    fprintf(fp,"grestore\n");
    fprintf(fp,"end %% Map_Painter_prolog\n");
    fprintf(fp,"%%%%EOF\n");
}

void ps_page_start(FILE *fp, int pagenum)
{
    fprintf(fp,"%%%%Page: ? %d\n",pagenum);
    fprintf(fp,"%%%%BeginPageSetup\n");
    fprintf(fp,"LM\n");
    fprintf(fp,"%%%%EndPageSetup\n");
}

void ps_page_end(FILE *fp)
{
    fprintf(fp,"GM showpage\n");
}

char newstr[MAXLINE];

char *ps_string(char *str)
{
    int sl, i, j;
    
    sl=strlen(str);

    for (i=0,j=0; i<sl; i++) {
      if (str[i]=='(' || str[i]==')' || str[i]==92)
          newstr[j++] = 92;
      newstr[j++] = str[i];
    }
    newstr[j] = '\0';
    return(newstr);
}


int main(int argc, char *argv[])
{
    int i, j;

    lib_init();

    if (argc!=5) {
        sprintf(ps, "usage: %s source code help doc dir\n", argv[0]);
        fprint(stderr,ps);
        abnormal_exit();
    }

    run {
            strcpy(file_name,argv[1]);
            strcpy(code_name,argv[2]);
            strcpy(hlp_name, argv[3]);
            strcpy(man_name, argv[4]);

            strcpy(code_failed,code_name);
            make_filename(code_failed,FORCE_EXTENSION,"failed");
            strcpy(final_hlp_name,hlp_name);
            make_filename_in_dir(final_hlp_name,FORCE_EXTENSION,"help",
                                 FORCE_DIR,argv[5]);

            file= open_file(file_name,READ);
            code= open_file(code_name,WRITE);
            hlp=  open_file(hlp_name, WRITE);
            man=  open_file(man_name,WRITE);

            topic= 0;
            entries= 0;
            matrix(entry,MAX_COMMANDS,MAXLINE+1,char);

            /* start help file 12345678901234567890123456789012345 */
            do_fwrite(hlp,       "#MAPMAKER help file - do not edit!\n");
            pos= 34l + LINE_BREAK_LEN;

            /* code file */
            do_fwrite(code,"/* MAPMAKER help code file - do not edit! */ \n\n");
//            do_fwrite(code,"#define INC_LIB \n#define INC_SHELL \n");
            do_fwrite(code,"#include \"system.h\" \n\n");
            /* sprintf(ps,"char help_filename[]= \"%s\";\n\n",final_hlp_name);
               do_fwrite(code,ps); */
            do_fwrite(code,"void make_help_entries(void)\n{\n");

            /* man file */
            man_write_title();

            /* get title */
            while (nullstr(str)) nextstr();
            if (str[0]!='@' || sscanf(str+1,"%s",type)!=1 || !streq(type,"title"))
                parse_error("need to start with a title",2);
            i=0; while (str[i++]!=' ') {}
            strcpy(title,str+i);

            /* get starting page# */
            do nextstr(); while (nullstr(str));
            if (str[0]!='@' || sscanf(str+1,"%s %d",type,&start_page)!=2 ||
                !streq(type,"page"))
                parse_error("need a page number",2);

            nextstr();
            while (TRUE) {
                while (nullstr(str)) nextstr();

                if (str[0]!='@' || sscanf(str+1,"%s",type)!=1) {
                    sprintf(ps, "error:  bad header:%s\n", str); pr();
                    close_files(argv[0]);
                    return(1);
                } else if (streq(type,"end")) {
                    close_files(argv[0]);
                    return(0);
                } else if (topic==0 && !streq(type,"topic"))
                    parse_error("need to start with a topic",2);

                i=0; while (str[i++]!=' ') {}
                j=i; while (str[j]!='(' && str[j]!='\0') j++;
                if (str[j]=='(') {
                    str[j]='\0';
                    if (sscanf(str+j+1,"%s",abbreviation)!=1 ||
                        len(abbreviation)>4 ||
                        abbreviation[len(abbreviation)-1]!=')')
                        parse_error("bad abbreviation",1);
                    else abbreviation[len(abbreviation)-1]='\0'; /* end ')' */
                } else abbreviation[0]='\0';

                strcpy(name,str+i); despace(name);
                nextstr();
                sprintf(ps, "\t%s...  ", name); pr(); flush();

                if      (streq(type,"cmd"))   parse_entry(CMD,name,abbreviation);
                else if (streq(type,"opt"))   parse_entry(OPT,name,abbreviation);
                else if (streq(type,"param")) parse_entry(PAR,name,abbreviation);
                else if (streq(type,"info"))  parse_entry(HLP,name,abbreviation);
                else if (streq(type,"topic")) parse_entry(TOP,name,abbreviation);
                else 	     		  parse_error("unknown type",1);
                nl();
            }

        } on_exit {
        if (msg==ENDOFILE) print("error:  unexpected end of file\n");
        else print("error: makehelp failed");
        close_files(argv[0]);
        rename_file(code_name,code_failed);
        return(1);
    }
    return(1); /* not reached */
}
