/******************************************************************************

 #####    ####           #    #    ##    #####    ####            ####
 #    #  #               ##  ##   #  #   #    #  #               #    #
 #    #   ####           # ## #  #    #  #    #   ####           #
 #####        #          #    #  ######  #####        #   ###    #
 #       #    #          #    #  #    #  #       #    #   ###    #    #
 #        ####  #######  #    #  #    #  #        ####    ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

//#define INC_LIB
//#define INC_SHELL
#include "mapm.h"
//#include "toplevel.h"
//#include "lowlevel.h"

#define LONGEST_MAP 300.0  /* longest map allowed for one page output */

void ps_file_start(FILE *fp);

void ps_file_end(FILE *fp);

void ps_page_start(FILE *fp, int pagenum);

void ps_page_end(FILE *fp);

real ps_scale(real map_length);

void ps_dump_chrom(FILE *fp, int chrom, real scale);

char *ps_loc_str(int i);

char *ps_frame_str(int i);

void
ps_file_start(FILE *fp) {
    fprintf(fp, "%%!PS-Adobe-3.0\n");
    fprintf(fp, "%%%%Creator: MAPMAKER\n");
    fprintf(fp, "%%%%LanguageLevel: 1\n");
    fprintf(fp, "%%%%PageOrder: Special\n");

    fprintf(fp, "%%%%BeginProlog\n");
    fprintf(fp, "%%%%BeginResource: procset Map_Painter_prolog\n");
    fprintf(fp, "/Map_Painter_prolog 100 dict def\n");
    fprintf(fp, "Map_Painter_prolog begin\n");
    fprintf(fp, "/CF {dup 0 eq\n");
    fprintf(fp, "       {pop /Times-Bold}\n");
    fprintf(fp, "       {1 eq {/Times-Roman} {/Times-Italic} ifelse}\n");
    fprintf(fp, "     ifelse} def\n");
    fprintf(fp, "/F {setfont} def\n");
    fprintf(fp, "/FF {findfont} def\n");
    fprintf(fp, "/GM {restore} def\n");
    fprintf(fp, "/GR {grestore} def\n");
    fprintf(fp, "/GS {gsave} def\n");
    fprintf(fp, "/LC {currentpoint pop dup /l exch def add /r exch def\n");
    fprintf(fp, "     counttomark 3 idiv {\n");
    fprintf(fp, "       currentpoint pop r ge {l currentpoint exch pop LLD sub moveto} if\n");
    fprintf(fp, "       LD\n");
    fprintf(fp, "     } repeat pop} def\n");
    fprintf(fp, "/LD {/f exch def /n exch def /p exch def\n");
    fprintf(fp, "     p length 0 ne {\n");
    fprintf(fp, "       /Times-Roman findfont LFS 1 sub scalefont setfont p show\n");
    fprintf(fp, "       2.5 0 rmoveto\n");
    fprintf(fp, "     } if\n");
    fprintf(fp, "     f CF findfont LFS scalefont setfont\n");
    fprintf(fp, "     n show 6 0 rmoveto} def\n");
    fprintf(fp, "/LM {save} def\n");
    fprintf(fp, "/LW {setlinewidth} def\n");
    fprintf(fp, "/S {show} def\n");
    fprintf(fp, "/SF {scalefont} def\n");
    fprintf(fp, "/SR {dup stringwidth pop neg 0 rmoveto show} def\n");
    fprintf(fp, "/TR {translate} def\n");
    fprintf(fp, "/XY {moveto} def\n");
    fprintf(fp, "end\n");
    fprintf(fp, "%%%%EndResource\n");
    fprintf(fp, "%%%%EndProlog\n");

    fprintf(fp, "/LFS 8 def /LLD 7.5 def\n");
    fprintf(fp, "%%%%BeginSetup\n");
    fprintf(fp, "Map_Painter_prolog begin\n");
    fprintf(fp, "gsave\n");
    fprintf(fp, "%%%%EndSetup\n");
}

void
ps_file_end(FILE *fp) {
    fprintf(fp, "%%%%Trailer\n");
    fprintf(fp, "grestore\n");
    fprintf(fp, "end %% Map_Painter_prolog\n");
    fprintf(fp, "%%%%EOF\n");
}

void
ps_page_start(FILE *fp, int pagenum) {
    fprintf(fp, "%%%%Page: ? %d\n", pagenum);
    fprintf(fp, "%%%%BeginPageSetup\n");
    fprintf(fp, "LM\n");
    fprintf(fp, "%%%%EndPageSetup\n");
}

void
ps_page_end(FILE *fp) {
    fprintf(fp, "GM showpage\n");
}

#define LAST_INTERVAL (-999.0)

void
print_ps_map(FILE *fp, MAP *map) {
    int i;
    double map_length = 0.0, interval_length, ps_length, scale;
    char *loc_str;

    loc_str = get_temp_string();

    ps_file_start(fp);
    ps_page_start(fp, 1);

    /* move to center of page, .5 inch from top */
    fprintf(fp, "4.0 72 mul 10.5 72 mul moveto\n");

    for (i = 0; i < map->num_loci - 1; i++)
        map_length += cm_dist(map->rec_frac[i][0]);

    scale = ps_scale(map_length);

    for (i = 0; i < map->num_loci - 1; i++) {
        sprintf(loc_str, "%s", ps_loc_str(map->locus[i]));
        interval_length = cm_dist(map->rec_frac[i][0]);
        while (interval_length < 1.0) {
            i++;
            sprintf(ps, "  %s", ps_loc_str(map->locus[i]));
            strcat(loc_str, ps);
            if (i == map->num_loci - 1) {
                interval_length = LAST_INTERVAL;
                break;
            } else {
                interval_length += cm_dist(map->rec_frac[i][0]);
            }
        }
        fprintf(fp, "GS .5 LW -4 0 rmoveto 8 0 rlineto stroke GR\n");
        fprintf(fp, "GS 10 -3 rmoveto 0 CF FF 9 SF F (%s)S GR\n", loc_str);

        if (interval_length != LAST_INTERVAL) {
            ps_length = -1.0 * interval_length * scale;
            fprintf(fp, "GS 2 LW 0 %.2lf rlineto stroke GR\n", ps_length);
            fprintf(fp, "GS -10 %.2lf rmoveto /Times-Roman FF 9 SF F (%.1lf)SR GR\n",
                    (ps_length / 2.0) - 3.0, interval_length);
            fprintf(fp, "0 %.2lf rmoveto\n", ps_length);
        }
    }
    if (i == map->num_loci - 1) {
        fprintf(fp, "GS .5 LW -4 0 rmoveto 8 0 rlineto stroke GR\n");
        fprintf(fp, "GS 10 -3 rmoveto 0 CF FF 9 SF F (%s)S GR\n",
                ps_loc_str(map->locus[i]));
    }

    ps_page_end(fp);
    ps_file_end(fp);
}

void
print_ps_chrom(FILE *fp, int chrom) {
    int i;
    real map_length, scale;

    if (chrom == NO_CHROM) {
        print_all_ps_chroms(fp);
        return;
    }

    map_length = 0.0;
    for (i = 0; i < chromosome->map_list[chrom]->num_loci - 1; i++)
        map_length += cm_dist(chromosome->map_list[chrom]->rec_frac[i][0]);

    scale = ps_scale(map_length);

    ps_file_start(fp);
    ps_page_start(fp, 1);

    ps_dump_chrom(fp, chrom, scale);

    ps_page_end(fp);
    ps_file_end(fp);
}

void
print_all_ps_chroms(FILE *fp) {
    int i, j;
    real map_length, best, scale;

    best = 0.0;
    for (i = 0; i < chromosome->num_maps; i++) {
        map_length = 0.0;
        for (j = 0; j < chromosome->map_list[i]->num_loci - 1; j++)
            map_length += cm_dist(chromosome->map_list[i]->rec_frac[j][0]);
        if (map_length > best) best = map_length;
    }
    scale = ps_scale(best);

    ps_file_start(fp);

    for (i = 0; i < chromosome->num_maps; i++) {
        ps_page_start(fp, i + 1);
        ps_dump_chrom(fp, i, scale);
        ps_page_end(fp);
    }

    ps_file_end(fp);
}

real
ps_scale(real map_length) {
    real scale;

    if (map_length > LONGEST_MAP) error("map too long for output");
    scale = 720.0 / map_length;
    if (scale > 12.0) scale = 12.0;

    return (scale);
}

#define PLACE_STRING 200

void
ps_dump_chrom(FILE *fp, int chrom, real scale) {
    int i, j, k, l, num_crunched, interval = 0, marker;
    double interval_length, ps_length, dist = 0.;
    char *loc_str, **placed_markers = NULL;
    MAP *frame;

    frame = chromosome->map_list[chrom];
    matrix(placed_markers, frame->num_loci + 1, PLACE_STRING, char);
    for (i = 0; i < frame->num_loci + 1; i++) strcpy(placed_markers[i], "");

    for (i = 0; i < raw.num_markers; i++) {

        if (haplotyped(i)) marker = haplo_first[i];
        else marker = i;

        if (assigned_to(marker, chrom) && placed_locus(marker)) {
            for (j = 0; j < placement[marker]->num_intervals; j++) {
                if (placement[marker]->like_ratio[j] == 0.00) {
                    interval = placement[marker]->interval[j];
                    dist = placement[marker]->distance[j];
                    break;
                }
            }
            if (interval == 0) {
                if (dist <= ZERO_DIST) {
                    sprintf(ps, "()(%s)%c ", ps_loc_str(i),
                            (placement[marker]->status == M_UNIQUE) ? '1' : '2');
                    strcat(placed_markers[1], ps);
                } else {
                    sprintf(ps, "(%.1lf)(%s)%c ", cm_dist(dist), ps_loc_str(i),
                            (placement[marker]->status == M_UNIQUE) ? '1' : '2');
                    strcat(placed_markers[0], ps);
                }
            } else {
                if (dist <= ZERO_DIST) {
                    sprintf(ps, "()(%s)%c ", ps_loc_str(i),
                            (placement[marker]->status == M_UNIQUE) ? '1' : '2');
                    strcat(placed_markers[interval], ps);

                } else if (interval != frame->num_loci &&
                           frame->rec_frac[interval - 1][0] - dist <= ZERO_DIST) {
                    sprintf(ps, "()(%s)%c ", ps_loc_str(i),
                            (placement[marker]->status == M_UNIQUE) ? '1' : '2');
                    strcat(placed_markers[interval + 1], ps);
                } else {
                    sprintf(ps, "(%.1lf)(%s)%c ", cm_dist(dist), ps_loc_str(i),
                            (placement[marker]->status == M_UNIQUE) ? '1' : '2');
                    strcat(placed_markers[interval], ps);
                }
            }
        }
    }

    if (!nullstr(frame->map_name)) {
        fprintf(fp, "72 11 72 mul 62 sub moveto\n");
        fprintf(fp, "/Times-Roman FF 14 SF F (Chromosome %s)S\n",
                frame->map_name);
    }

    /* move to center of page, .5 inch from top */
    fprintf(fp, "4.0 72 mul 10.5 72 mul moveto\n");

/********** make the centromere ball at top
    fprintf(fp,"GS newpath 4.0 72 mul 10.75 72 mul 10 0 360 arc closepath fill GR\n");
**********/

    fprintf(fp, "GS 2 LW 0 15 rlineto stroke GR\n");
    fprintf(fp, "GS 5 12 mul 7 rmoveto mark %s 16 12 mul LC GR\n", placed_markers[0]);

    for (i = 0; i < frame->num_loci - 1; i++) {
        loc_str = get_temp_string();
        sprintf(loc_str, "%s", ps_frame_str(frame->locus[i]));
        interval_length = cm_dist(frame->rec_frac[i][0]);
        num_crunched = 0;
        while (interval_length < 1.0) {
            i++;
            num_crunched++;
            if (i == frame->num_loci - 1) {
                interval_length = LAST_INTERVAL;
                break;
            } else {
                interval_length += cm_dist(frame->rec_frac[i][0]);
            }
        }
        for (k = 0; k < num_crunched; k++) {
            strcat(placed_markers[i + 1], placed_markers[i - k]);
        }

        for (k = 0; k < num_crunched; k++) {
            dist = 0.00;
            for (l = k; l < num_crunched; l++) dist += cm_dist(frame->rec_frac[i - (l + 1)][0]);

            if (dist <= ZERO_DIST) {
                sprintf(ps, "()(%s)0 ", ps_loc_str(frame->locus[i - k]));
                strcat(placed_markers[i + 1], ps);
            } else {
                sprintf(ps, "(%.1lf)(%s)0 ", dist, ps_loc_str(frame->locus[i - k]));
                strcat(placed_markers[i + 1], ps);
            }
        }

        if (haplotyped(frame->locus[i])) {
            for (j = haplo_first[frame->locus[i]]; j != NO_CHROM; j = haplo_next[j]) {
                if (j != haplo_first[frame->locus[i]]) {
                    sprintf(ps, "()(%s)0 ", ps_loc_str(j));
                    strcat(placed_markers[i + 1], ps);
                }
            }
        }

        fprintf(fp, "GS .5 LW -4 0 rmoveto 8 0 rlineto stroke GR\n");
        fprintf(fp, "GS 10 -3 rmoveto 0 CF FF 9 SF F (%s)S GR\n", loc_str);
        /* do something with placement markers */
        fprintf(fp, "GS 5 12 mul -3 rmoveto mark %s 16 12 mul LC GR\n",
                placed_markers[i + 1]);

        if (interval_length != LAST_INTERVAL) {
            ps_length = -1.0 * interval_length * scale;
            fprintf(fp, "GS 2 LW 0 %.2lf rlineto stroke GR\n", ps_length);
            fprintf(fp, "GS -10 %.2lf rmoveto /Times-Roman FF 9 SF F (%.1lf)SR GR\n",
                    (ps_length / 2.0) - 3.0, interval_length);
            fprintf(fp, "0 %.2lf rmoveto\n", ps_length);
        }
    }
    if (i == frame->num_loci - 1) {
        fprintf(fp, "GS .5 LW -4 0 rmoveto 8 0 rlineto stroke GR\n");
        fprintf(fp, "GS 10 -3 rmoveto 0 CF FF 9 SF F (%s)S GR\n",
                ps_frame_str(frame->locus[i]));

        if (haplotyped(frame->locus[i])) {
            for (j = haplo_first[frame->locus[i]]; j != NO_CHROM; j = haplo_next[j]) {
                if (j != haplo_first[frame->locus[i]]) {
                    sprintf(ps, "()(%s)0 ", ps_loc_str(j));
                    strcat(placed_markers[i + 1], ps);
                }
            }
        }

        fprintf(fp, "GS 5 12 mul -3 rmoveto mark %s 16 12 mul LC GR\n",
                placed_markers[i + 1]);
    }
    unmatrix(placed_markers, frame->num_loci + 1, char);
}


char *
ps_loc_str(int i) {
    char *tempstr;

    if (print_names) return (raw.locus_name[i]);
    else {
        tempstr = get_temp_string();
        sprintf(tempstr, "%d", i + 1);
        return (tempstr);
    }
}

char *
ps_frame_str(int i) {
    int mark;
    char *tempstr;

    if (haplotyped(i)) mark = haplo_first[i];
    else mark = i;

    if (print_names) return (raw.locus_name[mark]);
    else {
        tempstr = get_temp_string();
        sprintf(tempstr, "%d", mark + 1);
        return (tempstr);
    }
}





