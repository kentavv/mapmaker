/******************************************************************************

 #    #    ##     #####  #    #  #          #    #####            ####
 ##  ##   #  #      #    #    #  #          #    #    #          #    #
 # ## #  #    #     #    ######  #          #    #####           #
 #    #  ######     #    #    #  #          #    #    #   ###    #
 #    #  #    #     #    #    #  #          #    #    #   ###    #    #
 #    #  #    #     #    #    #  ######     #    #####    ###     ####

******************************************************************************/
/* This file is part of MAPMAKER 3.0b, Copyright 1987-1992, Whitehead Institute
   for Biomedical Research. All rights reserved. See READ.ME for license. */

#define INC_MATH
#define INC_MEM
#define INC_MSG
#define INC_STR 
#define INC_IO

#define INC_MISC
#define INC_HELP_DEFS
#include "system.h"

real sq(r) real r; { return(r*r); }
real rmaxf(r,s) real r,s; { return(r>=s ? r : s); }
real rminf(r,s) real r,s; { return(r<=s ? r : s); }

long lpow2(i)	/* Return 2 to the i-th power for 0<=i<=LONGBITS (31) */
int i;
{
	long result; int j;
	
	if (i<0 || i>LONGBITS) send(CRASH);
	for(result=(long)1, j=0; j<i; j++) result <<= 1;
	return(result);
}

int ipow2(i)	/* Return 2 to the i-th power 0<=i<=INTBITS (15) */
int i;
{
	int result; int j;
	
	if (i<0 || i>INTBITS) send(CRASH);
	for(result=1, j=0; j<i; j++) result <<= 1;
	return(result);
}


long lpow(x,i)	/* Return x to the i-th power */
int x,i;
{
	long result; int j;
	
	if (i<0 || i>LONGBITS) send(CRASH);
	for(result=(long)1, j=0; j<i; j++) result*=(long)x;
	return(result);
}

int ipow(x,i)	/* Return 2 to the i-th power 0<=i<=INTBITS (15) */
int i;
{
	int result; int j;
	
	if (i<0 || i>INTBITS) send(CRASH);
	for(result=1, j=0; j<i; j++) result*=x;
	return(result);
}

int ichoose(n,k) /* Best algorithm for small k */
int n, k;
{
    long a, b, Ln, Lk, L1, Lx;

    if (k>n || n<1 || k<1) send(CRASH);
    Ln= (long)n; Lk= (long)k; Lx= (long)(n-k+1); L1= (long)1;

    for (a=L1; Ln>=Lx; Ln--) a*=Ln; /* is product of n...n-k+1, or n!/(n-k)! */
    for (b=L1; Lk>L1; Lk--) b*=Lk;
    return((int)(a/b));
}

int imaxf(r,s)  int r,s;  { return(r>=s ? r : s); }
int iminf(r,s)  int r,s;  { return(r<=s ? r : s); }
long lmaxf(r,s) long r,s; { return(r>=s ? r : s); }
long lminf(r,s) long r,s; { return(r<=s ? r : s); }


/*** Compare Functions for Array Sort Routines ***/

int icomp(x,y)
QSORT_COMPARE_PTR_TO(int) *x, *y;
{ return( (*(int*)y) - (*(int*)x) ); }

int lcomp(x,y)
QSORT_COMPARE_PTR_TO(long) *x, *y;
{ long a, b; a= *(long*)x; b= *(long*)y; 
  if (a<b) return(-1); else if (a==b) return(0); else return(1); }

int rcomp(x,y)
QSORT_COMPARE_PTR_TO(real) *x, *y;
{ real a, b; a= *(real*)x; b= *(real*)y; 
  if (a<b) return(-1); else if (a==b) return(0); else return(1); }

int scomp(x,y)
QSORT_COMPARE_PTR_TO(char) **x, **y; /* does this work? */
{ return((int)strcmp(*(char**)x,*(char**)y)); }

int inv_icomp(x,y)
QSORT_COMPARE_PTR_TO(int) *x, *y;
{ return( (*(int*)x) - (*(int*)y) ); }

int inv_rcomp(x,y)
QSORT_COMPARE_PTR_TO(real) *x, *y;
{ real a, b; a= *(real*)x; b= *(real*)y; 
  if (a>b) return(-1); else if (a==b) return(0); else return(1); }

/* this needs work */
int rhistogram(data,length,min_num_buckets,
	scale_quantization,scale_limit_quantization)   /* return T or F */
real *data; 
int length, min_num_buckets;
real scale_quantization;	/* for present, ignored */
real scale_limit_quantization; 	/* for present, ignored */
{
    int i, j, index, stars, max_stars;
    int *bucket, num_buckets, largest_bucket; 
    real max_data, scale, min_data;
    real stars_per;
    
    num_buckets= min_num_buckets;
    for (i=0, max_data= -VERY_BIG, min_data= VERY_BIG; i<length; i++) {
	if (data[i]>max_data) max_data=data[i];
	if (data[i]<min_data) min_data=data[i];
    }
    scale= (max_data - min_data)/REAL(num_buckets); /* range per bucket */

    run { array(bucket, num_buckets, int); } 
        except_when (NOMEMORY) return(FALSE);
    for (j=0; j<num_buckets; j++) bucket[j]=0;

    for (i=0; i<length; i++) {
	if (data[i]==max_data) index=num_buckets-1; else 
	  index= (int)((data[i]-min_data)/scale);	
	if (!irange(&index,0,num_buckets-1)) send(CRASH);
	++(bucket[index]);
    }

    max_stars= LINE-10;
    for (j=0, largest_bucket=0; j<num_buckets; j++) 
      if (bucket[j]>largest_bucket) largest_bucket=bucket[j];
    if (largest_bucket<=max_stars) stars_per=REAL(max_stars/largest_bucket);
    else stars_per=REAL(max_stars)/REAL(largest_bucket);

    for (j=0; j<num_buckets; j++) {
	sf(ps,"%6.3lf | ",((real)j)*scale+(scale/2.0)+min_data); print(ps);
	stars= INT(REAL(bucket[j]) * stars_per);
	irange(&stars,0,max_stars);
	if (bucket[j]>0 && stars<1) print(":");
	else for (i=0; i<INT(stars); i++) print("*");
	nl();
    }
    unarray(bucket, int);
    return(TRUE);
}


real rmean(data,length)
real *data;
int length;
{
	int i;
	real total;
	
	if (length==0) return(0.0);
	for (i=0, total=0.0; i<length; i++) total+= data[i];
	return(total/REAL(length));
}


real rmaxin(data,length)
real *data;
int length;
{
	int i;
	real largest;
	
	if (length==0) return(0.0);
	for (i=0, largest= -VERY_BIG; i<length; i++) 
	    if (data[i]>largest) largest=data[i];
	return(largest);
}


real rmedian(data,length)
real *data;
int length;
{
	real *copy, median;
	
	if (length==0) return(0.0);
	array(copy, length, real);
	rcopy(copy, data, length);
	rsort(copy,length);
	median= rmiddle(copy,length);
	unarray(copy, real);
	return(median);
}


real rmiddle(data,length)
real *data;
int length;
{
	if (length==0) return(0.0);
	if (length<3)  return(data[0]);
	return(data[(length-1)/2]);
}
	

void rcopy(to,from,length)
real *to, *from;
int length;
{ int i; for (i=0; i<length; i++) to[i]=from[i]; }


void dummy_math_calls()
{
    real x,y,z; x=y=z=0.5;
 
    x=log10(y);
    x=pow(y,z);
    x=log(y);
    x=exp(y);

    x=floor(y);
    x=ceil(y);
    x=fabs(y);

    x=sin(y);
    x=cos(y);
    x=tan(y);
    x=sinh(y);
    x=cosh(y);
    x=tanh(y);
    x=asin(y);
    x=acos(y);
    x=atan(y);
}


real two_sigma_sq; 

real normal_func(deviation) /* return P(deviation) */
real deviation;
{ return(exp(-sq(deviation)/two_sigma_sq)); }

DISTRIBUTION *make_normal_dist(mu,sigma,inc,limit)
real mu, sigma, inc, limit;      /* inc and limit are fractions of sigma */
{
    two_sigma_sq= 2.0 * sq(sigma);    
    return(make_distribution(inc*sigma,inc,limit*sigma,mu,normal_func));
}



DISTRIBUTION *make_distribution(sigma,inc,limit,mean,prob_func)
real sigma, inc, limit, mean;   /* inc and limit here are offsets from mean */
real (*prob_func)();		/* a ptr to a function */
{
    real d, total_prob, cum_prob, end;
    real *prob, *deviation;
    int length, i, j;
    DISTRIBUTION *dist;

    prob=NULL; deviation=NULL; dist=NULL;
    run {
	end= limit + VERY_SMALL;
	length= ((int)(2.0*limit/inc + 0.99999)) + 1;
	
	array(prob, length, real);
	array(deviation, length, real);
	single(dist, DISTRIBUTION);
	
	total_prob= 0.0;
	for (d= -limit, i=0; d<end; d+=inc, i++) {
	    if (i>=length) send(CRASH);
	    deviation[i]= d+mean;
	    total_prob+= prob[i]= (*prob_func)(d);
	}
	
	cum_prob= 0.0;
	for (j=0; j<i; j++) {
	    cum_prob+= prob[j];
	    prob[j]= cum_prob/total_prob;
	}
	
	dist->entries= i;
	dist->start= deviation[0];
	dist->increment= inc;
	dist->mean= mean;
	dist->deviation= deviation;
	dist->prob= prob;

    } when_aborting {
	unsingle(dist, DISTRIBUTION);
	unarray(prob, real);
	unarray(deviation, real);
	relay;
	return((DISTRIBUTION*)NULL); /* never reached */
    }
    return(dist);
}


real pick_from_distribution(dist,prob)
DISTRIBUTION *dist;
real *prob;  /* may be NULL */
{
    real p, *probs;
    int i, last;
    
    p= randnum();
    probs= dist->prob; last= dist->entries-1;
    for (i=0; probs[i]<p; i++) if (i==last) break;
    if (prob!=NULL) *prob=p;
    return(dist->deviation[i]);
}


void eliminate(); /* defined below */

void mat_invert(m,size,m_inverse)  /* Invert square matrix by Gauss' method */
real **m;
int size;
real **m_inverse;		
/* m_inverse should be 2*size columns (2nd index) by size rows (first index)
   its left square will be left with the result (ignore the right side) */
{
    int row, col, c, twice_size;
    real diff, value;
    
    twice_size= 2 * size;
    for (row=0; row<size; row++)
      for (col=0; col<size; col++) {
	  m_inverse[row][col]= m[row][col];
	  m_inverse[row][col+size]= (row==col ? 1.0 : 0.0);
      }
    
    for (col=0; col<size; col++) {
	if (m_inverse[col][col]==0.0) {
	    for (row=0; row<size && m_inverse[row][col]==0.0; row++);
	    if (row==size) send(CRASH);
	    for (c=0; c<twice_size; c++)
	      m_inverse[col][c]+= m_inverse[row][c];
	}
	for (row=0; row<size; row++)
	  if (row!=col) eliminate(row,col,col,m_inverse,twice_size);
    }
    
    for (row=0; row<size; row++) {
	value= m_inverse[row][row];
	for (col=0; col<twice_size; col++)
	  m_inverse[row][col]/= value;
    }
    
    for (row=0; row<size; row++)
      for (col=0; col<size; col++)
	m_inverse[row][col]= m_inverse[row][col+size];
    
/*    if ((diff=test_mult(m,m_inverse,size))>mat_tolerance) 
      { MATINV_diff= diff; send(MATINV); } OLD */
}


void eliminate(row,col,row_to_sub,m,n_cols)
int row, col, row_to_sub;
real **m;
int n_cols;
{
    int j;
    real value;
    
    /* Get m[row][col]=0 by subtracting the row with a 1 in col
       from m[row], where the row is scaled appropriately */
    
    if (m[row][col]==0.0) return;
    if (m[row_to_sub][col]==0.0) send(CRASH); /* used to send(SINGMAT) */
                                              /* but SINGMAT was undefined */
    value= m[row][col]/m[row_to_sub][col];
    for (j=0; j<n_cols; j++) m[row_to_sub][j]*= value;
    for (j=0; j<n_cols; j++) m[row][j]-= (m[row_to_sub][j]);
}


void mat_mult(m,m2,size,result) /* NEEDS TESTING */
real **m, **m2;
int size;
real **result;
{ 
    int i, j, k;
    real entry;
    
    for (i=0; i<size; i++)
      for (j=0; j<size; j++) {
	  for (k=0, entry= 0.0; k<size; k++) entry+= m[i][k] * m2[k][j];
	  result[i][j]= entry;
      }
}


void array_times_matrix(a,b,rows,columns,c)
real *a, **b;
int rows,columns;
real *c;
{
    int i, j;
    
    for (i=0; i <columns; i++) 
      for (j=0, c[i] = 0.0; j < rows; j++)  c[i] += a[j] * b[j][i];
}



