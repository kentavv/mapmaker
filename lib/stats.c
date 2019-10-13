/******************************************************************************

  ####    #####    ##     #####   ####            ####
 #          #     #  #      #    #               #    #
  ####      #    #    #     #     ####           #
      #     #    ######     #         #   ###    #
 #    #     #    #    #     #    #    #   ###    #    #
  ####      #    #    #     #     ####    ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

//#define INC_LIB
//#define INC_MISC
//#define INC_STATS
//#define INC_HELP_DEFS
#include "system.h"

NORMAL_TEST *check_normalcy(real *dist, int dist_size)
{
    NORMAL_TEST *normal_results;

    single(normal_results, NORMAL_TEST);
    skew(normal_results, dist, dist_size);
    quartile(normal_results, dist, dist_size);
    within(normal_results, dist, dist_size);
    return(normal_results);
}

void skew(NORMAL_TEST *normal_results, real *dist, int dist_size)
{
    real mean=0,second_moment=0,second_total=0,fourth_total=0;
    real third_moment=0,fourth_moment=0,third_total=0;
    real total=0,to_be=0,skewness=0,kurt=0;
    int i;

    for(i=0;i<dist_size;i++) 
      total += dist[i];
    mean = total/dist_size;
    normal_results->mean = mean;

    /* calculate second, third, and fourth moments */
    for(i=0;i<dist_size;i++) {
	to_be = dist[i] - mean;
	second_total += (to_be*to_be);
	third_total += (to_be*to_be*to_be);
	fourth_total += (to_be*to_be*to_be*to_be);
    }
    second_moment = second_total/dist_size;
    normal_results->sigma = sqrt(second_moment);
    third_moment = third_total/dist_size;
    fourth_moment = fourth_total/dist_size;

    skewness = third_moment/pow(second_moment, 1.5);
    kurt = (fourth_moment/pow(second_moment,2.0)) - 3;

    normal_results->skew = skewness;
    normal_results->kurt = kurt;
}


void quartile(NORMAL_TEST *normal_results, real *dist, int dist_size)
{
    real first_quart=0,third_quart=0,quart_ratio=0;

    rsort(dist,dist_size);
    first_quart = dist[dist_size/4];
    third_quart = dist[dist_size - dist_size/4];
    quart_ratio = (third_quart-first_quart)/(1.3490*normal_results->sigma);
    normal_results->quart_ratio = quart_ratio;
}

void within(NORMAL_TEST *normal_results, real *dist, int dist_size)
{
    int i;
    int one_fourth=0,one_half=0,one=0,two=0,three=0;

    for (i=0;i<dist_size;i++) {
	if (dist[i] > normal_results->mean - .25*normal_results->sigma &&
	    dist[i] < normal_results->mean + .25*normal_results->sigma)
	  one_fourth++;
	if (dist[i] > normal_results->mean - .5*normal_results->sigma &&
	    dist[i] < normal_results->mean + .5*normal_results->sigma)
	  one_half++;
	if (dist[i] > normal_results->mean - normal_results->sigma &&
	      dist[i] < normal_results->mean + normal_results->sigma)
	  one ++;
	if (dist[i] > normal_results->mean - 2*normal_results->sigma &&
	    dist[i] < normal_results->mean + 2*normal_results->sigma)
	  two ++;
	if (dist[i] > normal_results->mean - 3*normal_results->sigma &&
	    dist[i] < normal_results->mean + 3*normal_results->sigma)
	  three ++;
    }
    
    normal_results->within_one_fourth = (real)(one_fourth)/(real)(dist_size);
    normal_results->within_one_half = (real)(one_half)/(real)(dist_size);
    normal_results->within_one = (real)(one)/(real)(dist_size);
    normal_results->within_two = (real)(two)/(real)(dist_size);
    normal_results->within_three = (real)(three)/(real)(dist_size);
}




void print_normal(NORMAL_TEST *to_be_printed, real lamda)
{
    

    if (lamda < 0) {
	print("-------------------------------------------------------------");
	print("----------------");
	nl();
	print("distribution:");
	print("                     ");
	print(" quartile |   fraction within n deviations:");
	nl();
        print("mean   ");
	print("sigma   ");
	print("skewness  ");
	print("kurtosis ");
	print(" ratio    |  ");
	print(" 1/4 ");
	print("   1/2 ");
	print("  1   ");
	print("  2   ");
	print("  3   ");
	nl();
	sprintf(ps, "%s ", rs(6.2, to_be_printed->mean));pr();
	sprintf(ps, "%s ", rs(6.2, to_be_printed->sigma));pr();
	sprintf(ps, " %s   ", rs(6.2, to_be_printed->skew));pr();
	sprintf(ps, " %s   ", rs(6.2, to_be_printed->kurt));pr();
	sprintf(ps, " %s     |  ", rs(4.2, to_be_printed->quart_ratio));pr();
	sprintf(ps, " %s  ", rs(4.2, to_be_printed->within_one_fourth));pr();
	sprintf(ps, " %s ", rs(4.2, to_be_printed->within_one_half));pr();
	sprintf(ps, " %s ", rs(4.2, to_be_printed->within_one));pr();
	sprintf(ps, " %s ", rs(4.2, to_be_printed->within_two));pr();
	sprintf(ps, " %s   ", rs(4.2, to_be_printed->within_three));pr();
	/*nl();*/
	print("-------------------------------------------------------------");
	print("----------------");
	nl();
    }
    else {
	sprintf(ps, "%lf  ", lamda);pr();
	sprintf(ps, "%.2lf ", to_be_printed->mean);pr();
	sprintf(ps, " %.2lf  ", to_be_printed->sigma);pr();
	sprintf(ps, " %.2lf ", to_be_printed->skew);pr();
	sprintf(ps, " %.2lf ", to_be_printed->kurt);pr();
	sprintf(ps, " %.2lf   ", to_be_printed->quart_ratio);pr();
	sprintf(ps, " %.2lf ", to_be_printed->within_one_fourth);pr();
	sprintf(ps, " %.2lf ", to_be_printed->within_one_half);pr();
	sprintf(ps, " %.2lf ", to_be_printed->within_one);pr();
	sprintf(ps, " %.2lf ", to_be_printed->within_two);pr();
	sprintf(ps, " %.2lf ", to_be_printed->within_three);pr();
	nl();
	print("-------------------------------------------------------------");
	print("----------------");
	nl();
    }
}
    

void box_cox(real start, real stop, real step, real *dist, int dist_size)
{
    int i;
    real lamda,*new_dist;
    NORMAL_TEST *new_test;

    lamda = start;

    single(new_test,NORMAL_TEST);
    array(new_dist,dist_size,real);
    

    while (lamda <= stop) {
	if (lamda <= 0)
	  return;
	for (i=0;i<dist_size;i++) {
	    if (dist[i] < 0) {
		print ("I'm sorry box_cox doesn't work on negative numbers\n");
		return;
	    } else
	  new_dist[i] = (pow(dist[i], lamda) - 1)/lamda;
    }
	new_test = check_normalcy(new_dist,dist_size);
	print_normal(new_test,lamda);
	lamda += step;
    }

}



static void print_histogram(int intervals, real title);

void print_rhisto(real *dist, int dist_size)
{
    int i, intervals;
    int biggest, more_than_2=0, more_than_1_5=0;
    int more_than_1=0,more_than_5=0,more_than_0=0,less_than_5=0;
    int less_than_1=0,less_than_1_5=0,less_than_2=0,less_than_2_5=0;
    real mean,sigma,sub_total=0.0;

    mean = rmean(dist,dist_size);
    for (i=0;i<dist_size;i++) sub_total += (dist[i]-mean)*(dist[i]-mean);
    sigma = sqrt(sub_total/dist_size);

    /* this first loop asertains the size of each interval for the histogram */
    for (i=0;i<dist_size;i++) {
	if (dist[i]>=mean-2.5*sigma&&dist[i]<mean-2*sigma)
	  more_than_2++;
	if (dist[i]>=mean-2*sigma&&dist[i]<mean-1.5*sigma)
	  more_than_1_5++;
	if (dist[i]>=mean-1.5*sigma&&dist[i]<mean-sigma)
	  more_than_1++;
	if (dist[i]>=mean-sigma&&dist[i]<mean-.5*sigma)
	  more_than_5++;
	if (dist[i]>=mean-.5*sigma&&dist[i]<mean)
	  more_than_0++;
	if (dist[i]>=mean&&dist[i]<mean+.5*sigma)
	  less_than_5++;
	if (dist[i]>=mean+.5*sigma&&dist[i]<mean+sigma)
	  less_than_1++;
	if (dist[i]>=mean+sigma&&dist[i]<mean+1.5*sigma)
	  less_than_1_5++;
	if (dist[i]>=mean+1.5*sigma&&dist[i]<mean+2*sigma)
	  less_than_2++;
	if (dist[i]>=mean+2*sigma&&dist[i]<mean+2.5*sigma)
	  less_than_2_5++;
    }
    /* We must now find which interval has the most values in it */
    
    biggest = more_than_2;
    if (more_than_1_5 > biggest) biggest = more_than_1_5;
    if (more_than_1 > biggest) biggest = more_than_1;
    if (more_than_5 > biggest) biggest = more_than_5;
    if (more_than_0 > biggest) biggest = more_than_0;
    if (less_than_5 > biggest) biggest = less_than_5;
    if (less_than_1 > biggest) biggest = less_than_1;
    if (less_than_1_5 > biggest) biggest = less_than_1_5;
    if (less_than_2 > biggest) biggest = less_than_2;
    if (less_than_2_5 > biggest) biggest = less_than_2_5;
    
    /* now we will print out the histogram, with a * for each .02 in prob */

    intervals=(int)(more_than_2*60/biggest);
    print_histogram(intervals,mean-2*sigma);
    intervals=(int)(more_than_1_5*60/biggest);
    print_histogram(intervals,mean-1.5*sigma);
    intervals=(int)(more_than_1*60/biggest);
    print_histogram(intervals,mean-sigma);
    intervals=(int)(more_than_5*60/biggest);
    print_histogram(intervals,mean-.5*sigma);
    intervals=(int)(more_than_0*60/biggest);
    print_histogram(intervals,mean);
    intervals=(int)(less_than_5*60/biggest);
    print_histogram(intervals,mean+.5*sigma);
    intervals=(int)(less_than_1*60/biggest);
    print_histogram(intervals,mean+sigma);
    intervals=(int)(less_than_1_5*60/biggest);
    print_histogram(intervals,mean+1.5*sigma);
    intervals=(int)(less_than_2*60/biggest);
    print_histogram(intervals,mean+2*sigma);
    intervals=(int)(less_than_2_5*60/biggest);
    print_histogram(intervals,mean+2.5*sigma);
    
}


void print_histogram(int intervals, real title)
{
    int i;

    sprintf(ps, "%s  |", rsd(6.2, title));pr();
    for(i=0;i<intervals;i++) print("*");
    nl();
}
