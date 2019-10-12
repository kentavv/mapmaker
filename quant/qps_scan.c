/******************************************************************************

  ####   #####    ####            ####    ####     ##    #    #           ####
 #    #  #    #  #               #       #    #   #  #   ##   #          #    #
 #    #  #    #   ####            ####   #       #    #  # #  #          #
 #  # #  #####        #               #  #       ######  #  # #   ###    #
 #   #   #       #    #          #    #  #    #  #    #  #   ##   ###    #    #
  ### #  #        ####  #######   ####    ####   #    #  #    #   ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#define INC_LIB
#define INC_SHELL
#define INC_CALLQCTM
#define INC_QTOPLEVEL
#define INC_QLOWLEVEL
#include "qtl.h"

void ps_file_start(), ps_file_end(), ps_page_start(), ps_page_end();
void do_bezier(), draw_axes(), draw_x();

#define XZERO 0.0
#define YZERO 0.0

#define SOLID_LINE    ""
#define THICK_LINE    "2 LW"
#define THIN_LINE     ".5 LW"
#define EVEN_DASH     "[4 4] 0 setdash"
#define BIG_DASH      "[6 2] 0 setdash"
#define SMALL_DASH    "[2 6] 0 setdash"
#define DOTTED        "[1 2] 0 setdash"

real xscale, yscale, best;
#define MAX_CHROM_LEN 1000

void print_ps_wiggle_order(wiggle, order, threshold)
int wiggle, order;
real threshold;
{
    int i, j, count=0, num_notches=0, pagenum=1;
    double *xval, *yval, *notch, highest, longest, current_len;
    WIGGLE_OPERATION *op;
    WIGGLE_INTERVAL **data;
    WIGGLE_POINT *point;
    FILE *fp;
    char **label, *trait_str, *filename;

    best = 0.0;
    array(xval, MAX_CHROM_LEN, real);
    array(yval, MAX_CHROM_LEN, real);
    array(notch, MAX_CHROM_LOC, real);
    matrix(label, MAX_CHROM_LOC, NAME_LEN+1, char);
    filename = get_temp_string();
    sprintf(filename, "scan%d_%d.ps", wiggle + 1, order + 1);
    fp = fopen(filename,"w");

    ps_file_start(fp);
    ps_page_start(fp, pagenum);

    op = wiggles[wiggle]; data = op->data[order];

    /***** Read through the data and calculate scale based on *****/
    /*****      longest chromosome and highest lodscore.      *****/

    longest=0.0; highest=0.0; current_len=0.0;
    for (i=0; i<op->num_wiggled_intervals; i++) {
        if (i>0 && !data[i]->contig) {
	    /* chromosome done */
	    if(current_len > longest) longest = current_len;
	    current_len = 0.0;
	}
	for (j=0; j<data[i]->num_points; j++) {
	    point=data[i]->point[j];
	    if(point->lod_score > highest) highest = point->lod_score;
	    current_len += 2.0;
	}
    }
    if(current_len > longest) longest = current_len;
    xscale = (double) 700.0/longest;
    yscale = rminf(50.0, 400.0/highest);

    trait_str = get_temp_string();
    sprintf(trait_str, "LOD score - Trait %d (%s)", op->trait + 1,
            raw.trait_name[op->trait]);
    if(!print_names) sprintf(label[0], "%d",
                             raw.original_locus[data[0]->map->left[data[0]->map->num_intervals-1]]);
    else             sprintf(label[0], "%s",
                             raw.locus_name[data[0]->map->left[data[0]->map->num_intervals-1]]);
    for (i=0; i<op->num_wiggled_intervals; i++) {
        if (i>0) { 
	    notch[num_notches++] = count*2.0; 
	    if (!print_names) sprintf(label[num_notches], "%d",
                                  raw.original_locus[data[i-1]->map->right[data[i-1]->map->num_intervals-1]]);
	    else              sprintf(label[num_notches], "%s",
                                  raw.locus_name[data[i-1]->map->right[data[i-1]->map->num_intervals-1]]);
	}
        if (i>0 && !data[i]->contig) {
	    /* a page is done */
	    draw_axes(fp, notch, num_notches, label, trait_str, threshold);
	    do_bezier(fp, xval, yval, count, 0.0, 0.0, SOLID_LINE);
	    ps_page_end(fp);
	    pagenum++;
	    ps_page_start(fp, pagenum);
	    if(!print_names) sprintf(label[0], "%d",
                                 raw.original_locus[data[i]->map->left[data[i]->map->num_intervals-1]]);
	    else             sprintf(label[0], "%s",
                                 raw.locus_name[data[i]->map->left[data[i]->map->num_intervals-1]]);
	    count = 0; num_notches = 0;
	    best = 0.0;
	}
	for(j=0; j<data[i]->num_points; j++) {
	    point=data[i]->point[j];
	    xval[count] = count*2.0;
	    yval[count] = point->lod_score;
	    if(point->lod_score>best) best = point->lod_score;
	    count++;
	}
    }

    /* dump last page */
    notch[num_notches++] = count*2.0; 
    if(!print_names) sprintf(label[num_notches], "%d",
                             raw.original_locus[data[i-1]->map->right[data[i-1]->map->num_intervals-1]]);
    else             sprintf(label[num_notches], "%s",
                             raw.locus_name[data[i-1]->map->right[data[i-1]->map->num_intervals-1]]);

    draw_axes(fp, notch, num_notches, label, trait_str, threshold);
    do_bezier(fp, xval, yval, count, 0.0, 0.0, SOLID_LINE);

    ps_page_end(fp);
    ps_file_end(fp);

    fclose(fp);
    sprintf(ps, "scan %d.%d saved in PostScript file '%s'\n",
       wiggle+1,order+1, filename); pr();
    unarray(xval, real);
    unarray(yval, real);
    unarray(notch, real);
    unmatrix(label, MAX_CHROM_LOC, char);
}


void 
draw_axes (FILE *fp, double *xnotch, int num_notches, char **label, char *y_name, double dotted_val)
{
    int i;
    double prev, next, current;

    fprintf(fp,"500 50 translate\n");
    fprintf(fp,"90 rotate\n");
    fprintf(fp,"%.2lf %.2lf moveto\n", XZERO, YZERO);
    fprintf(fp,"GS 0.5 LW 700 0 rlineto stroke GR\n");
    fprintf(fp,"GS 0.5 LW 0 450 rlineto stroke GR\n");
    fprintf(fp,"GS 0 460 rmoveto /Times-Roman FF 14 SF F (%s)S GR\n",y_name);

    for(i=1; i<=(int)(best+1.0); i++) {
        fprintf(fp,"%.2lf %.2lf moveto\n", XZERO-4, YZERO+((double)i*yscale));
	fprintf(fp,"GS 0.5 LW 9 0 rlineto stroke GR\n");
	fprintf(fp,"GS -3 -3 rmoveto /Times-Roman FF 9 SF F (%d.0)SR GR\n",i);
    }
    /* draw dotted line - y = dotted_val */
    fprintf(fp,"%.2lf %.2lf moveto\n", XZERO, YZERO+(dotted_val*yscale));
    fprintf(fp,"GS 0.5 LW [4 4] 0 setdash 700 0 rlineto stroke GR\n");

    prev = 0.0; 
    for (i=0; i<num_notches; i++) {
        current = XZERO+(xnotch[i]*xscale);
	next = (i < num_notches-1) ? XZERO+(xnotch[i+1]*xscale) : 999.9;

	if(next-current < 5.0) {
	    if(current-prev > 6.0) current -= 1.0;
	}
        fprintf(fp,"%.2lf %.2lf moveto\n",current,YZERO-4);
	fprintf(fp,"GS 0.5 LW 0 2 rmoveto 0 4 rlineto stroke GR\n");
     fprintf(fp,"GS -90 rotate 0 -3 rmoveto /Times-Roman FF 8 SF F (%s)S GR\n",
		label[i+1]);
	prev = current;
    }

    fprintf(fp,"%.2lf %.2lf moveto\n", XZERO, YZERO-4);
    fprintf(fp,"GS -90 rotate 0 -3 rmoveto /Times-Roman FF 8 SF F (%s)S GR\n",
	    label[0]);
    draw_x(fp);
}


void 
do_bezier (FILE *fp, double *xval, double *yval, int num_points, double s0, double sn, char *line_type)
{
    int i;
    double *x0, *x1, *x2, *x3, *y0, *y1, *y2, *y3, *slope, s1, s2;

    array(x0, num_points, real);
    array(x1, num_points, real);
    array(x2, num_points, real);
    array(x3, num_points, real);
    array(y0, num_points, real);
    array(y1, num_points, real);
    array(y2, num_points, real);
    array(y3, num_points, real);
    array(slope, num_points, real);

    /* Bezier functionality draws an appropriate curve between points (x0,y0)
       and (x3,y3) based on Bezier cubic control points (x1,y1) and (x2,y2)
       (x1,y1) lies on the tangent from (x0,y0) and (x2,y2) lies on the tangent
       from (x3,y3) with x0 < x1 < x2 < x3 for our cases since 
       x3 = x0 + scan_interval_length (usually 2 cM) */


    /* (x0,y0) and (x3,y3) are the known data points, scaled to size of page */
 
    for(i = 0; i < num_points-1; i++) {
        x0[i] = xval[i]*xscale + XZERO;
	y0[i] = yval[i]*yscale + YZERO;
	x3[i] = xval[i+1]*xscale + XZERO;
	y3[i] = yval[i+1]*yscale + YZERO;
    }

    /* estimate the slope of the tangent at each data point...given three 
       points on a curve (x1,y1) (x2,y2) (x3,y3), a good guess for the 
       slope of the tangent at (x2,y2) is the average of the slopes of the 
       lines connecting (x1,y1)(x2,y2) and (x2,y2)(x3,y3) - since this method 
       does not work for the first or last point on the curve for obvious 
       reasons, those values are handed into the function as s0 and sn...
       usually values of 0.0 will be just fine */

    slope[0] = s0;
    for(i = 1; i < num_points-1; i++) {
	s1 = (y3[i-1]-y0[i-1])/(x3[i-1]-x0[i-1]);
	s2 = (y3[i]-y0[i])/(x3[i]-x0[i]);
	slope[i] = (s1+s2)/2.0;
	if(slope[i] > 10.0) slope[i] = 10.0;
	if(slope[i] < -10.0) slope[i] = -10.0;
    }
    slope[i] = sn;


    /* calculate points (x1,y1) and (x2,y2) now that we have approximate 
       tangent lines for our real data points */

    for(i = 0; i < num_points-1; i++) {
        x1[i] = (x3[i]-x0[i])/3.0;
	y1[i] = y0[i] + x1[i]*slope[i];
	x2[i] = x1[i];
	y2[i] = y3[i] - x2[i]*slope[i+1];
	x1[i] += x0[i];
	x2[i] += x1[i];
    }

    /* Postscript conveniently draws a Bezier cubic section from the current 
       position to position (x3,y3) based on the calculated control points
       (x1,y1) and (x2,y2) */

    fprintf(fp,"%.2lf %.2lf moveto\n",x0[0],y0[0]);

    for(i = 0; i < num_points-1; i++) {
        fprintf(fp,"GS %s %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf curveto stroke GR\n",
		line_type, x1[i], y1[i], x2[i], y2[i], x3[i], y3[i]);
        fprintf(fp,"%.2lf %.2lf moveto\n", x3[i], y3[i]);
	draw_x(fp);
    }
    fprintf(fp,"%.2lf %.2lf moveto\n", x3[i-1], y3[i-1]);
    draw_x(fp);

    unarray(x0, real);
    unarray(x1, real);
    unarray(x2, real);
    unarray(x3, real);
    unarray(y0, real);
    unarray(y1, real);
    unarray(y2, real);
    unarray(y3, real);
    unarray(slope, real);
}


void 
draw_x (void)
{
#ifdef DRAWX
    fprintf(fp,"GS -3 -3 rlineto stroke GR\n");
    fprintf(fp,"GS 3 3 rlineto stroke GR\n");
    fprintf(fp,"GS 3 -3 rlineto stroke GR\n");
    fprintf(fp,"GS -3 3 rlineto stroke GR\n");
#endif
}



void 
ps_file_start (FILE *fp)
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
    fprintf(fp,"/VSdict 4 dict def\n");
    fprintf(fp,"/VS { VSdict begin\n");
    fprintf(fp,"      /thestring exch def\n");
    fprintf(fp,"      /lineskip exch def\n");
    fprintf(fp,"      thestring\n");
    fprintf(fp,"       {\n");
    fprintf(fp,"        /charcode exch def\n");
    fprintf(fp,"        /thechar ( ) dup 0 charcode put def\n");
    fprintf(fp,"        0 lineskip neg rmoveto\n");
    fprintf(fp,"        gsave\n");
    fprintf(fp,"         thechar stringwidth pop 2 div neg 0 rmoveto\n");
    fprintf(fp,"         thechar show\n");
    fprintf(fp,"        grestore\n");
    fprintf(fp,"      } forall\n");
    fprintf(fp,"     end\n");
    fprintf(fp,"   } def\n");

    fprintf(fp,"end\n");
    fprintf(fp,"%%%%EndResource\n");
    fprintf(fp,"%%%%EndProlog\n");

    fprintf(fp,"/LFS 8 def /LLD 7.5 def\n");
    fprintf(fp,"%%%%BeginSetup\n");
    fprintf(fp,"Map_Painter_prolog begin\n");
    fprintf(fp,"gsave\n");
    fprintf(fp,"%%%%EndSetup\n");
}

void 
ps_file_end (FILE *fp)
{
    fprintf(fp,"%%%%Trailer\n");
    fprintf(fp,"grestore\n");
    fprintf(fp,"end %% Map_Painter_prolog\n");
    fprintf(fp,"%%%%EOF\n");
}

void 
ps_page_start (FILE *fp, int pagenum)
{
    fprintf(fp,"%%%%Page: ? %d\n",pagenum);
    fprintf(fp,"%%%%BeginPageSetup\n");
    fprintf(fp,"LM\n");
    fprintf(fp,"%%%%EndPageSetup\n");
}

void 
ps_page_end (FILE *fp)
{
    fprintf(fp,"GM showpage\n");
}




void 
print_ps_multi_wiggle (int wiggle, real threshold)
{
    int i, j, count=0, num_notches=0, pagenum=1, order;
    double **xval, **yval, *notch, highest, longest, current_len;
    WIGGLE_OPERATION *op;
    WIGGLE_POINT *point;
    FILE *fp;
    char **label, *trait_str, *filename, *line_choice();

    op = wiggles[wiggle]; 

    if(op->num_orders > 4) error("too many wiggle orders to display on one graph");

    best = 0.0;
    matrix(xval, op->num_orders, MAX_CHROM_LEN, real);
    matrix(yval, op->num_orders, MAX_CHROM_LEN, real);
    array(notch, MAX_CHROM_LOC, real);
    matrix(label, MAX_CHROM_LOC, NAME_LEN+1, char);
    filename = get_temp_string();
    sprintf(filename, "scan%d_x.ps", wiggle + 1);
    fp = fopen(filename,"w");

    ps_file_start(fp);
    ps_page_start(fp, pagenum);
    longest=0.0; highest=0.0; current_len=0.0; 
    for (i=0; i<op->num_wiggled_intervals; i++) {
        if (i>0 && !op->data[0][i]->contig) {
	    if(current_len > longest) longest = current_len;
	    current_len = 0.0;
	}
	for (j=0; j<op->data[0][i]->num_points; j++) {
	    point=op->data[0][i]->point[j];
	    if(point->lod_score > highest) highest = point->lod_score;
	    current_len += 2.0;
	}
    }

    if(current_len > longest) longest = current_len;
    xscale = (double) 700.0/longest;
    yscale = rminf(50.0, 400.0/highest);

    trait_str = get_temp_string();
    sprintf(trait_str, "LOD score - Trait %d (%s)", op->trait + 1,
            raw.trait_name[op->trait]);
    if (!print_names) sprintf(label[0], "%d",
                              raw.original_locus[op->data[0][0]->map->left[op->data[0][0]->map->num_intervals-1]]);
    else              sprintf(label[0], "%s",
                              raw.locus_name[op->data[0][0]->map->left[op->data[0][0]->map->num_intervals-1]]);
    for (i=0; i<op->num_wiggled_intervals; i++) {
        if (i>0) { 
	    notch[num_notches++] = count*2.0; 
	    if (!print_names) sprintf(label[num_notches], "%d",
                                  raw.original_locus[op->data[0][i-1]->map->right[op->data[0][i-1]->map->num_intervals-1]]);
	    else              sprintf(label[num_notches], "%s",
                                  raw.locus_name[op->data[0][i-1]->map->right[op->data[0][i-1]->map->num_intervals-1]]);
	}
        if (i>0 && !op->data[0][i]->contig) {

	    draw_axes(fp, notch, num_notches, label, trait_str, threshold);
	    for (order=0; order<op->num_orders; order++)
	      do_bezier(fp, xval[order], yval[order], count, 0.0, 0.0, line_choice(order));

	    fprintf(fp,"0 -40 moveto\n");
	    fprintf(fp,"GS /Times-Roman FF 9 SF F (GENETICS:)S GR\n");
	    for (order=0; order<op->num_orders; order++) {
	      fprintf(fp,"0 %d moveto\n",-50-(order*10));
	      fprintf(fp,"GS /Times-Roman FF 9 SF F (%s)S GR\n",
	 genetics_str(&op->data[order][0]->map->constraint[op->num_intervals-1],FALSE));
	      fprintf(fp,"%d %d moveto\n",50,-50-(order*10));
	      fprintf(fp,"GS %s 20 0 rlineto stroke GR\n",line_choice(order));
	    }
	    ps_page_end(fp);
	    pagenum++;
	    ps_page_start(fp, pagenum);
	    if (!print_names) sprintf(label[0], "%d",
                                  raw.original_locus[op->data[0][i]->map->left[op->data[0][i]->map->num_intervals-1]]);
	    else              sprintf(label[0], "%s",
                                  raw.locus_name[op->data[0][i]->map->left[op->data[0][i]->map->num_intervals-1]]);
	    count = 0; num_notches = 0;
	    best = 0.0;
	}
	for (j=0; j<op->data[0][i]->num_points; j++) {
	  for (order=0; order<op->num_orders; order++) {
	    point=op->data[order][i]->point[j];
	    xval[order][count] = count*2.0;
	    yval[order][count] = point->lod_score;
	    if(point->lod_score>best) best = point->lod_score;
	  }
	  count++;
	}
    }

    notch[num_notches++] = count*2.0; 
    if (!print_names) sprintf(label[num_notches], "%d",
                              raw.original_locus[op->data[0][i-1]->map->right[op->data[0][i-1]->map->num_intervals-1]]);
    else              sprintf(label[num_notches], "%s",
                              raw.locus_name[op->data[0][i-1]->map->right[op->data[0][i-1]->map->num_intervals-1]]);
    draw_axes(fp, notch, num_notches, label, trait_str, threshold);
    for (order=0; order<op->num_orders; order++) 
      do_bezier(fp, xval[order], yval[order], count, 0.0, 0.0, 
		line_choice(order));
      
    fprintf(fp,"0 -40 moveto\n");
    fprintf(fp,"GS /Times-Roman FF 9 SF F (GENETICS:)S GR\n");
    for (order=0; order<op->num_orders; order++) {
        fprintf(fp,"0 %d moveto\n",-50-(order*10));
	fprintf(fp,"GS /Times-Roman FF 9 SF F (%s)S GR\n",
	genetics_str(&op->data[order][0]->map->constraint[op->num_intervals-1],
		     FALSE));
	fprintf(fp,"%d %d moveto\n",50,-50-(order*10));
	fprintf(fp,"GS %s 20 0 rlineto stroke GR\n",line_choice(order));
    }

    ps_page_end(fp);
    ps_file_end(fp);

    fclose(fp);
    sprintf(ps, "scan %d.x saved in PostScript file '%s'\n", wiggle + 1, filename); pr();
    unarray(xval, real);
    unarray(yval, real);
    unarray(notch, real);
    unmatrix(label, MAX_CHROM_LOC, char);
}

char *
line_choice (int order)
{
    switch(order) {
        case 0: return(THICK_LINE);
	case 1: return(THIN_LINE);
	case 2: return(DOTTED);
	case 3: return(EVEN_DASH);
	default: return(SOLID_LINE);
    }
}

