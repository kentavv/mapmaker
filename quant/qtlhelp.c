/* MAPMAKER help code file - do not edit! */ 

#define INC_LIB 
#define INC_SHELL 
#include "system.h" 

void make_help_entries()
{
 mkhelp("about mapmaker/qtl","",573l,EXACTLY,-1,CMD,1,
        "License and Contact Information for MAPMAKER/QTL",
        "",
        "");
 mkhelp("release notes","",1584l,EXACTLY,-1,HLP,1,
        "Information on the Release Notes",
        "",
        "");
 mkhelp("starting mapmaker/qtl","",2413l,EXACTLY,-1,HLP,1,
        "How to Start MAPMAKER/QTL and Available Options",
        "",
        "");
 mkhelp("tutorial","",7023l,EXACTLY,-1,HLP,1,
        "Information on the MAPMAKER/QTL Tutorial",
        "",
        "");
 mkhelp("entering commands","",7439l,EXACTLY,-1,HLP,1,
        "How to Type Commands into MAPMAKER/QTL",
        "",
        "");
 mkhelp("abbreviating commands","",8696l,EXACTLY,-1,HLP,1,
        "How to Abbreviate Commands You Type",
        "",
        "");
 mkhelp("keyboard editing","",9815l,EXACTLY,-1,HLP,1,
        "Editing Commands and Sequences Using the Keyboard",
        "",
        "");
 mkhelp("postscript output","",11358l,EXACTLY,-1,HLP,1,
        "What to Do with PostScript Graphic Output",
        "",
        "");
 mkhelp("help","?",13078l,UPTO,1,CMD,2,
        "Read On-Line Help Information",
        "<command name or topic number>",
        "with no arguments, display a list of all commands and topics");
 mkhelp("photo","",13594l,UPTO,1,CMD,2,
        "Begin Saving MAPMAKER/QTL Output to a Text File",
        "<file-name>",
        "with no arguments, display the current 'photo' status");
 mkhelp("load data","ld",15051l,EXACTLY,1,CMD,2,
        "Load Data Files into MAPMAKER/QTL",
        "<file-name>",
        "");
 mkhelp("save data","",16450l,EXACTLY,0,CMD,2,
        "Save MAPMAKER/QTL Data and Status Info to Disk",
        "",
        "");
 mkhelp("quit","q",16733l,EXACTLY,0,CMD,2,
        "Quit from a MAPMAKER/QTL Session",
        "",
        "");
 mkhelp("sequence","s",17829l,UPTO,1,CMD,3,
        "Set the Sequence of Intervals to Look for QTLs in",
        "<sequence>",
        "with no arguments, displays the current sequence");
 mkhelp("history","h",26776l,EXACTLY,1,CMD,3,
        "List Previous Sequences",
        "<number of previous sequences to display>",
        "20");
 mkhelp("list loci","",27302l,EXACTLY,0,CMD,3,
        "Display Trait and Locus Names and Numbers",
        "",
        "");
 mkhelp("edit sequence","",27567l,EXACTLY,0,CMD,3,
        "Interactive Keyboard Driven Sequence Editor",
        "",
        "");
 mkhelp("let","",28049l,EXACTLY,2,CMD,3,
        "Name a Sequence or Portion of a Sequence",
        "<name> = <sequence>",
        "");
 mkhelp("names","",29323l,EXACTLY,0,CMD,3,
        "Display Values of All Names Set Using 'Let'",
        "",
        "");
 mkhelp("forget named sequence","",29611l,EXACTLY,1,CMD,3,
        "Forget the Assignment of a Name (Undo 'Let')",
        "<name>",
        "");
 mkhelp("trait","t",30411l,UPTO,1,CMD,4,
        "Set the Trait for QTL Mapping",
        "<trait name or number>",
        "with no arguments, displays the current trait");
 mkhelp("make trait","",30935l,EXACTLY,2,CMD,4,
        "Add a New Trait (as a Function of Existing Traits)",
        "<name> = <equation>",
        "");
 mkhelp("show trait","",31604l,UPTO,1,CMD,4,
        "Display Statistics and Histogram for a Given Trait",
        "<trait>",
        "the current trait");
 mkhelp("list traits","",32026l,EXACTLY,0,CMD,4,
        "Display a List of Traits Currently in Data Set",
        "",
        "");
 mkhelp("forget trait","",32258l,EXACTLY,1,CMD,4,
        "Remove a Trait from the Data Set",
        "<trait>",
        "");
 mkhelp("map","m",33373l,EXACTLY,0,CMD,5,
        "Display QTL-Maps for Specified Intervals",
        "",
        "");
 mkhelp("scan","",39494l,UPTO,3,CMD,5,
        "Search for QTLs, Stepping Through the Right Interval",
        "<cM-step> <LOD-threshold> <scale>",
        "1.0 2.0 0.25");
 mkhelp("list scans","",45664l,UPTO,1,CMD,5,
        "List Scan Results Which Have Been Saved",
        "<scan number>",
        "");
 mkhelp("show scan","",47707l,UPTO,3,CMD,5,
        "Display Saved Scan Data for a Particular Scan",
        "<scan-number> <LOD-threshold> <scale>",
        "<last-scan> 2.0 0.25");
 mkhelp("draw scan","",48486l,UPTO,1,CMD,5,
        "Draws Saved Scan Data in PostScript",
        "<scan-number> <LOD-threshold>",
        "<last-scan> 2.0");
 mkhelp("forget scan","",49378l,EXACTLY,1,CMD,5,
        "Delete Particular Scan Results from Memory",
        "<scan-number>",
        "");
 mkhelp("forget all scans","",49887l,EXACTLY,0,CMD,5,
        "Delete All Scan Results from Memory",
        "",
        "");
 mkhelp("show peaks","",50245l,UPTO,3,CMD,5,
        "Display LOD Peaks and QTL Maps for a Particular Scan",
        "<scan-number> <LOD-threshold> <falloff>",
        "<last-scan> 2.0 -1.0");
 mkhelp("show trys","",51794l,UPTO,2,CMD,5,
        "Display LOD Scores for a Scan Using ':try' Genetics",
        "<scan-number> <LOD-threshold>",
        "<last-scan> 2.0");
 mkhelp("compare","",55088l,EXACTLY,0,CMD,5,
        "Compare QTL-Maps for Specified Interval Combinations",
        "",
        "");
 mkhelp("list compares","",55933l,EXACTLY,0,CMD,5,
        "List Compare Results Which Have Been Saved",
        "",
        "");
 mkhelp("show compare","",56324l,UPTO,3,CMD,5,
        "Display Saved Data from a Particular 'Compare'",
        "<compare-number> <LOD-threshold> <falloff>",
        "<last-compare> 2.0 -1.0");
 mkhelp("forget compare","",56446l,EXACTLY,1,CMD,5,
        "Delete Particular Compare Results from Memory",
        "<compare-number>",
        "");
 mkhelp("forget all compares","",56644l,EXACTLY,0,CMD,5,
        "Delete All Compare Results from Memory",
        "",
        "");
 mkhelp("show best maps","",56756l,UPTO,3,CMD,5,
        "Display Best Results from a Particular 'Compare'",
        "<compare-number> <LOD-threshold> <falloff>",
        "<last-compare> 2.0 -1.0");
 mkhelp("show linkage map","",57057l,EXACTLY,1,CMD,5,
        "Display Linkage Map Locus Order and Distances",
        "<locus name or number, or 'let' name>",
        "Displays linkage map for all chromosomes");
 mkhelp("print names","",57437l,UPTO,1,OPT,6,
        "Set Displaying of Locus Names or Numbers",
        "<'on' or 'off'>",
        "with no arguments, displays the current value");
 mkhelp("units","",57571l,UPTO,1,PAR,6,
        "Display Recombination Fractions or cM Distances",
        "<'cm' or 'rf'>",
        "with no arguments, displays the current value");
 mkhelp("auto save data","",58521l,UPTO,1,OPT,6,
        "Force MAPMAKER/QTL to Save Data When You Quit",
        "<'on' or 'off'>",
        "with no arguments, displays the current value");
 mkhelp("more mode","",59408l,UPTO,1,OPT,6,
        "Enable Pauses After Each Screen of Output",
        "<'on' or 'off'>",
        "with no arguments, displays the current value");
 mkhelp("run","r",59965l,EXACTLY,1,CMD,7,
        "Accept Commands from an Input File",
        "<file-name>",
        "");
 mkhelp("change directory","cd",61320l,UPTO,1,CMD,7,
        "Changes Default Directory Used by MAPMAKER/QTL",
        "<directory-name>",
        "with no arguments, display the current default directory");
 mkhelp("system","!",61776l,UPTO,1,CMD,7,
        "Run a system command, or go to a system prompt",
        "<command>",
        "");
 mkhelp("previous commands","",63412l,UPTO,1,CMD,7,
        "Display previous MAPMAKER/QTL commands",
        "<number of commands to display>",
        "<all previous commands>");
 mkhelp("review output","",63865l,EXACTLY,0,CMD,7,
        "Display previous MAPMAKER/QTL output",
        "",
        "");
 mkhelp("time","",64223l,EXACTLY,0,CMD,7,
        "Display current time (useful for timing commands)",
        "",
        "");
 mkhelp("comment","",64414l,UPTO,1,CMD,7,
        "Enter a comment into the photo file (does nothing)",
        "<comment>",
        "with no args, allows you to enter a long comment");
 mkhelp("wizard mode","",65039l,EXACTLY,1,CMD,7,
        "Danger Will Robinson, Danger!",
        "<on or off>",
        "");

 mktopic(1,"General Information on MAPMAKER/QTL Version 1.1",TOP,62l);
 mktopic(2,"Basic MAPMAKER/QTL Commands",TOP,12649l);
 mktopic(3,"Sequence Command and Related Features",TOP,16974l);
 mktopic(4,"Trait Command and Related Features",TOP,29874l);
 mktopic(5,"MAPMAKER/QTL mapping commands",TOP,32668l);
 mktopic(6,"MAPMAKER/QTL Parameters and Options",TOP,57340l);
 mktopic(7,"Miscellaneous Commands",TOP,59890l);
}
