/*
 *  test_hash_speed.c
 *  pfr
 *
 *  Created by Eric Anderson on 3/11/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */



#include <stdio.h>
#include <stdlib.h>

#define UN_EXTERN

#include "ECA_MemAlloc.h"
#include "MathStatRand.h"
#include "ECA_Opt.h"
#include "ranlib.h"
#include "uthash.h"


/* this lets us make a hash, so we can retrieve the d of an individual, given his name */
struct str_hash {
	double d;  
	char Name[15];
	UT_hash_handle hh;         /* makes this structure hashable */
};


int main(int argc, char *argv[])
{
	int i,j,k,t,N=100,REPS=10000000,n1,n2,n3;
	double ***P;
	double d=1.0,sum=0.0;
	struct str_hash *Phash=NULL;
	struct str_hash *IN,*OUT;
	char str[20];
	
if(0) {
	/* HERE IS OUR CODE FOR DOING IT ALL IN TERMS OF A 3-D ARRAY */
	P=(double ***)calloc(N,sizeof(double **));
	for(i=0;i<N;i++)  {
		P[i]=(double **)calloc(N,sizeof(double *));
		for(j=0;j<N;j++)  {
			P[i][j]=(double *)calloc(N,sizeof(double));
			for(k=0;k<N;k++)  {
				P[i][j][k]=d;
				d+=.23453;
			}
		}
	}
	
	printf("Done Filling Matrix\n");
	
	
	
	for(t=0;t<REPS;t++)  {
		n1 = UniformRV(0,N-1);
		n2 = UniformRV(0,N-1);
		n3 = UniformRV(0,N-1);
		
		sum += P[n1][n2][n3];
		
		if(t%10000==0) printf("Done with rep %d.  Tot=%f\n",t,sum);
	}	
return(0);
}

	/* HERE IS OUR CODE FOR DOING IT ALL IN TERMS OF A HASH */
	for(i=0;i<N;i++)  {
		for(j=0;j<N;j++)  {
			for(k=0;k<N;k++)  {
				IN=(struct str_hash *)malloc(sizeof(struct str_hash));
				sprintf(IN->Name,"%d,%d,%d",i,j,k);
				IN->d=d;
				HASH_ADD_STR(Phash,Name,IN);
				d+=.23453;
			}
		}
	}
	
	
	printf("Done Filling Hash \n");
	
	for(t=0;t<REPS;t++)  {
		n1 = UniformRV(0,N-1);
		n2 = UniformRV(0,N-1);
		n3 = UniformRV(0,N-1);
		
		sprintf(str,"%d,%d,%d",n1,n2,n3);
		HASH_FIND_STR(Phash,str,OUT);
		sum += OUT->d;
		
		if(t%10000==0) printf("Done with rep %d.  Tot=%f\n",t,sum);
	}
	

	return(0);
}