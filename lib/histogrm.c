/******************************************************************************

 #    #     #     ####    #####   ####    ####   #####   #    #           ####
 #    #     #    #          #    #    #  #    #  #    #  ##  ##          #    #
 ######     #     ####      #    #    #  #       #    #  # ## #          #
 #    #     #         #     #    #    #  #  ###  #####   #    #   ###    #
 #    #     #    #    #     #    #    #  #    #  #   #   #    #   ###    #    #
 #    #     #     ####      #     ####    ####   #    #  #    #   ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#define INC_LIB
#define INC_MISC
#define INC_HISTO
#include "system.h"


void print_histogram();

void make_histo(dist, dist_size)
real *dist;
int dist_size;
{
    int i, intervals;
    int more_than_2_25=0,more_than_2=0,more_than_1_75=0,more_than_1_5=0;
    int more_than_1_25=0,more_than_1=0,more_than_75=0,more_than_5=0;
    int more_than_25=0,more_than_0=0,less_than_25=0,less_than_5=0;
    int less_than_75=0,less_than_1=0,less_than_1_25=0,less_than_1_5=0;
    int less_than_1_75=0,less_than_2=0,less_than_2_25=0,less_than_2_5=0;
    real probab,mean,sigma,sub_total=0;

    mean = rmean(dist,dist_size);
    for (i=0;i<dist_size;i++) sub_total += (dist[i]-mean)*(dist[i]-mean);
    sigma = sub_total/dist_size;

    /* this first loop asertains the size of each interval for the histogram */
    for (i=0;i<dist_size;i++) {
	if (dist[i]>=mean-2.5*sigma&&dist[i]<mean-2.25*sigma) 
	  more_than_2_25++;
	if (dist[i]>=mean-2.25*sigma&&dist[i]<mean-2*sigma)
	  more_than_2++;
	if (dist[i]>=mean-2*sigma&&dist[i]<mean-1.75*sigma)
	  more_than_1_75++;
	if (dist[i]>=mean-1.75*sigma&&dist[i]<mean-1.5*sigma)
	  more_than_1_5++;
	if (dist[i]>=mean-1.5*sigma&&dist[i]<mean-1.25*sigma)
	  more_than_1_25++;
	if (dist[i]>=mean-1.25*sigma&&dist[i]<mean-sigma)
	  more_than_1++;
	if (dist[i]>=mean-sigma&&dist[i]<mean-.75*sigma)
	  more_than_75++;
	if (dist[i]>=mean-.75*sigma&&dist[i]<mean-.5*sigma)
	  more_than_5++;
	if (dist[i]>=mean-.5*sigma&&dist[i]<mean-.25*sigma)
	  more_than_25++;
	if (dist[i]>=mean-.25*sigma&&dist[i]<mean)
	  more_than_0++;
	if (dist[i]>=mean&&dist[i]<mean+.25*sigma)
	  less_than_25++;
	if (dist[i]>=mean+.25*sigma&&dist[i]<mean+.5*sigma)
	  less_than_5++;
	if (dist[i]>=mean+.5*sigma&&dist[i]<mean+.75*sigma)
	  less_than_75++;
	if (dist[i]>=mean+.75*sigma&&dist[i]<mean+sigma)
	  less_than_1++;
	if (dist[i]>=mean+sigma&&dist[i]<mean+1.25*sigma)
	  less_than_1_25++;
	if (dist[i]>=mean+1.25*sigma&&dist[i]<mean+1.5*sigma)
	  less_than_1_5++;
	if (dist[i]>=mean+1.5*sigma&&dist[i]<mean+1.75*sigma)
	  less_than_1_75++;
	if (dist[i]>=mean+1.75*sigma&&dist[i]<mean+2*sigma)
	  less_than_2++;
	if (dist[i]>=mean+2*sigma&&dist[i]<mean+2.25*sigma)
	  less_than_2_25++;
	if (dist[i]>=mean+2.25*sigma&&dist[i]<mean+2.5*sigma)
	  less_than_2_5++;
    }

    /* now we will print out the histogram, with a * for each .02 in prob */

    probab=(real)(more_than_2_25)/(real)(dist_size);
    intervals=(int)(probab/.004);
    print_histogram(intervals,"-2.50");
    probab=(real)(more_than_2)/(real)(dist_size);
    intervals=(int)(probab/.004);
    print_histogram(intervals,"-2.25");
    probab=(real)(more_than_1_75)/(real)(dist_size);
    intervals=(int)(probab/.004);
    print_histogram(intervals,"-2.00");
    probab=(real)(more_than_1_5)/(real)(dist_size);
    intervals=(int)(probab/.004);
    print_histogram(intervals,"-1.75");
    probab=(real)(more_than_1_25)/(real)(dist_size);
    intervals=(int)(probab/.004);
    print_histogram(intervals,"-1.50");
    probab=(real)(more_than_1)/(real)(dist_size);
    intervals=(int)(probab/.004);
    print_histogram(intervals,"-1.25");
    probab=(real)(more_than_75)/(real)(dist_size);
    intervals=(int)(probab/.004);
    print_histogram(intervals,"-1.00");
    probab=(real)(more_than_5)/(real)(dist_size);
    intervals=(int)(probab/.004);
    print_histogram(intervals,"-0.75");
    probab=(real)(more_than_25)/(real)(dist_size);
    intervals=(int)(probab/.004);
    print_histogram(intervals,"-0.50");
    probab=(real)(more_than_0)/(real)(dist_size);
    intervals=(int)(probab/.004);
    print_histogram(intervals,"-0.25");
    probab=(real)(less_than_25)/(real)(dist_size);
    intervals=(int)(probab/.004);
    print_histogram(intervals,"+0.25");
    probab=(real)(less_than_5)/(real)(dist_size);
    intervals=(int)(probab/.004);
    print_histogram(intervals,"+0.50");
    probab=(real)(less_than_75)/(real)(dist_size);
    intervals=(int)(probab/.004);
    print_histogram(intervals,"+0.75");
    probab=(real)(less_than_1)/(real)(dist_size);
    intervals=(int)(probab/.004);
    print_histogram(intervals,"+1.00");
    probab=(real)(less_than_1_25)/(real)(dist_size);
    intervals=(int)(probab/.004);
    print_histogram(intervals,"+1.25");
    probab=(real)(less_than_1_5)/(real)(dist_size);
    intervals=(int)(probab/.004);
    print_histogram(intervals,"+1.50");
    probab=(real)(less_than_1_75)/(real)(dist_size);
    intervals=(int)(probab/.004);
    print_histogram(intervals,"+1.75");
    probab=(real)(less_than_2)/(real)(dist_size);
    intervals=(int)(probab/.004);
    print_histogram(intervals,"+2.00");
    probab=(real)(less_than_2_25)/(real)(dist_size);
    intervals=(int)(probab/.004);
    print_histogram(intervals,"+2.25");
    probab=(real)(less_than_2_5)/(real)(dist_size);
    intervals=(int)(probab/.004);
    print_histogram(intervals,"+2.50");
}


void print_histogram(intervals, title)
int intervals;
char title[10];
{
    int i;

    sprintf(ps, "%s  |", title);pr();
    for(i=0;i<intervals;i++) print("*");
    nl();
    return;
}

