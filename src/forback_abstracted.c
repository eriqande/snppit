/*
 *  forback_abstracted.c
 *  pfr
 *
 *  Created by Eric C. Anderson on 3/9/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "forback_abstracted.h"
#include "MathStatRand.h"


/* this function is the type of function that will be the main one for 
doing the forward calculation.  I am going to do this one here for the case of
d=1.  that is, each A-state is in a group characterized by a v-vector that is 
(0) or (1).  

It takes as arguments:

	forabs_univ *FU  : all the settings for the problem that are not the A-states
	
	double **A  :  a vector of a states subscripted by [t][j], where t=1,...,T and j = 0,...,NA-1
	
	double **PX  : an output parameter which is a two-dimensional array where we store the probability of
		  each X state given the A states that belong to it.  (A simple sum, but we will want to keep 
		  this thing stored for later use on the backward step.)  Note that memory must be allocated to this already!
		  
	double *Over : an output parameter which is an array for the prob that the state is in S^uparrow, for each
		  time t.  This must have memory allocated to it before callingn this function.
		  
It returns:

	A d+1-dimensional array.  (So, if d=1, it is essentially just a list of vectors, each element of which is the 
	probability of the S variable for a given time t.  
	

*/
double **ForwardStep_One_d(forabs_univ *FU, double **PA, double **PX, double *Over  ) 
{
	double **ret;
	int x,t;
	int s[MAX_FORABS_DIMS+1], snew[MAX_FORABS_DIMS+1];
	int T=FU->T;
	int d=FU->d;
	int NX=FU->NX;
	int **v=FU->v;
	int *Smax=FU->Smax;
	int Test;
	double Prob;
	
	
	/* first, we have to compute the X-probs */
	ComputeXProbs( FU, PA, PX );
	
	/* then initialize Over */
	for(t=0;t<=T;t++)  {
		Over[t] = 0.0;
	}
	
	/* now, allocate memory to ret.  All this gets initialized to zero automatically by calloc */
	ret=(double **)calloc(FU->T+1, sizeof(double *));
	for(t=0;t<=FU->T;t++) {
		ret[t]=(double *)calloc(   ECA_MIN(FU->Smax[0],t)  + 1  , sizeof(double));
	}
	
	
	/* now, do the forward step */
	for(t=0;t<T;t++)  {
		
		/* cycle over all the states at time t */
		for(s[0]=0;  s[0]<=ECA_MIN(FU->Smax[0],t);  s[0]++) {
			
			/* cycle over X states */
			for(x=0;x<NX;x++)  {
				Test = S_status_plus( s, v[x], t, Smax, snew, d);
				Prob = ret[t][ s[0] ] * PX[t+1][x];
				if(Test==0) { /* if adding v to s makes an allowable s for the next time step */
					ret[t+1][ snew[0] ] += Prob;
				 } 
				 else if(Test==1) {  /* otherwise, if adding v to s makes s_t+1 go into the S^upparow set, then add it to over */
					Over[t+1] += Prob;
				 }
				 else {
					fprintf(stderr,"Error! Weird dude, S_status_plus returned %d in ForwardStep_One_d\n",Test);
				 }
			}
		}
	}
	
	return(ret);
}



double ***ForwardStep_Two_d(forabs_univ *FU, double **PA, double **PX, double *Over  ) 
{
	double ***ret;
	int x,t;
	int s[MAX_FORABS_DIMS+1], snew[MAX_FORABS_DIMS+1];
	int T=FU->T;
	int d=FU->d;
	int NX=FU->NX;
	int **v=FU->v;
	int *Smax=FU->Smax;
	int Test;
	double Prob;
	
	
	/* first, we have to compute the X-probs */
	ComputeXProbs( FU, PA, PX );
	
	/* then initialize Over */
	for(t=0;t<=T;t++)  {
		Over[t] = 0.0;
	}
	
	/* now, allocate memory to ret.  All this gets initialized to zero automatically by calloc */
	ret=(double ***)calloc(FU->T+1, sizeof(double **));
	for(t=0;t<=FU->T;t++) {
		ret[t]=(double **)calloc(   ECA_MIN(FU->Smax[0],t)  + 1  , sizeof(double *));
		for(s[0]=0; s[0]<=ECA_MIN(FU->Smax[0],t); s[0]++)  {
			ret[t][ s[0] ] = (double *)calloc(  ECA_MIN(FU->Smax[1],t)  + 1, sizeof(double));
		}
	}
	
	
	/* now, do the forward step */
	for(t=0;t<T;t++)  {
		
		/* cycle over all the states at time t */
		for(s[0]=0;  s[0]<=ECA_MIN(FU->Smax[0],t);  s[0]++) {
		for(s[1]=0;  s[1]<=ECA_MIN(FU->Smax[1],t);  s[1]++) {
			
			/* cycle over X states */
			for(x=0;x<NX;x++)  {
				Test = S_status_plus( s, v[x], t, Smax, snew, d);
				Prob = ret[t][ s[0] ][ s[1] ] * PX[t+1][x];
				if(Test==0) { /* if adding v to s makes an allowable s for the next time step */
					ret[t+1][ snew[0] ][ snew[1] ] += Prob;
				 } 
				 else if(Test==1) {  /* otherwise, if adding v to s makes s_t+1 go into the S^upparow set, then add it to over */
					Over[t+1] += Prob;
				 }
				 else {
					fprintf(stderr,"Error! Weird dude, S_status_plus returned %d in ForwardStep_One_d\n",Test);
				 }
			}
		}}
	}
	
	return(ret);
}



/* this is a simple function that computes the P_X^{(t)}(k) for each X-state */
void ComputeXProbs( forabs_univ *FU, double **PA, double **PX)
{
	int t,j,k;
	int T=FU->T;
	int NX=FU->NX;
	int *CX=FU->CX;
	int **X=FU->X;
	
	/* cycle over time steps (loci) */
	for(t=1;t<=T;t++)  {
		
		/* initialize all X Probs to zero */
		for(j=0;j<NX;j++)  {
			PX[t][j] = 0.0;
		}
		
		/* now, add stuff up into each X state */
		for(j=0;j<NX;j++)  { /* cycle over the X states */
			for(k=0;k<CX[j];k++)  { /* cycle over the A states in this X-state */
				PX[t][j] += PA[t][ X[j][k] ];
			}
		}		
	}
}


/* these functions takes a current S state, s, which is a vector of length d, 
   and it adds (or substracts)  the vector v from that and then tests 
   whether the result is a valid point in S^\downarrow.  if it is it returns 0.  
   If one of the components of the resulting vector is less than 0, it returns -1.  If
   one of the components is greater than t or the corresponding component in smax, it returns 1. 
   
   note that t should be the index of the S variable that you are adding or subtracting something
   from!!!!   In other words, if you are testing what happens when you add v[x] to S_t, then
   t that you pass to the function will be t. The t and smax comparisons will be made based on
   T+1.  Same if you want to know what happens when you subtract v[x] from S_t: you pass t to the function,
   but it will make the relevant comparison to t-1.
   
   newS is an output parameter that is a single array of length d giving the S vector that results from adding
   or subtracting the v vector to s.
   
*/	
int S_status_plus( int *s, int *v, int t, int *Smax, int *newS, int d)
{
	int i,r;
	
	for(i=0;i<d;i++)  {
		r = s[i]+v[i];
		
		if( r > ECA_MIN(Smax[0],t+1) ) {
			return(1);
		}
		if( r < 0 ) {   /* this really shouldn't need to be in here.  I may take it out for efficiency */
			return(-1);
		}
		newS[i]=r;
	}
	return(0);
}


int S_status_minus( int *s, int *v, int t, int *Smax, int *newS, int d)
{
	int i,r;
	
	for(i=0;i<d;i++)  {
		r = s[i]-v[i];
		
		if( r > ECA_MIN(Smax[0],t-1) ) { 
			return(1);
		}
		if( r < 0 ) {   
			return(-1);
		}
		newS[i]=r;
	}
	return(0);
}






