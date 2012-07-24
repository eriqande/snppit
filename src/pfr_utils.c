/*
 *  pfr_utils.c
 *  pfr
 *
 *  Created by Eric C. Anderson on 8/2/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <limits.h>

#include "ECA_MemAlloc.h"
#include "MathStatRand.h"
#include "ranlib.h"
#include "ECA_Opt.h"
#include "MCTypesEtc.h"
#include "pfr_forback.h"
#include "pfr_utils.h"



/***************************************************************

Following functions all deal with setting up the value of 
different state and probability variables

***************************************************************/ 

int XNumOfaType(analtype aType) 
{
	switch(aType) {
		case(MATRI):
		case(PATRI):
			return(3);
		case(TRI):
			return(2);
		case(MAPATRI): 
			return(5);
		default:
			fprintf(stderr,"Unknown analtype %d in XNumOfaType. Exiting...\n",aType);
			exit(1);
	}
}


/*! This function returns the number of different genotype configurations possible corresponding 
to the type of analysis.  It is used to set the NumA variable of the pfr_forback_data struct right in the 
beginning of run time. */
int NumAofaType(analtype aType)
{
	switch(aType) {
		case(MATRI):
		case(PATRI):
		case(TRI):
		case(MAPATRI): 
			return(27);
		default:
			fprintf(stderr,"Unknown analtype %d in NumAOfaType. Exiting...\n",aType);
			exit(1);
	}

}
int XdimOfaType(analtype aType) 
{
	switch(aType) {
		case(MATRI):
		case(PATRI):
			return(2);
		case(TRI):
			return(1);
		case(MAPATRI): 
			return(3);
		default:
			fprintf(stderr,"Unknown analtype %d in XdimOfaType. Exiting...\n",aType);
			exit(1);
	}
}

/*! Allocate memory to and set values of XKeys given analysis type */
int **SetXKeys(analtype aType) {
	int **ret;
	int dim,num,i;

	dim = XdimOfaType(aType);
	num = XNumOfaType(aType);	
	ret = (int **)calloc(num,sizeof(int *));
	for(i=0;i<num;i++)  {
		ret[i] = (int *)calloc(dim,sizeof(int));
	}
	
	if(aType==TRI) {
		ret[0][0] = 0;
		ret[1][0] = 1;
	}
	else if(aType==MATRI || aType==PATRI) {
		ret[0][0] = 0;  /* this is the way I named these in my notes PFR 7/254/06.  It isn't the most rational, but */
		ret[0][1] = 0;  /* I am sort of stuck with it now. i hope that I am not fighting it the whole way. */
		
		ret[1][0] = 1;
		ret[1][1] = 1;
		
		ret[2][0] = 0;
		ret[2][1] = 1;
	}
	else if(aType==MAPATRI) {
		ret[0][0] = 0; 
		ret[0][1] = 0;
		ret[0][2] = 0;
		
		ret[1][0] = 1;
		ret[1][1] = 0;
		ret[1][2] = 1;
		
		ret[2][0] = 0;
		ret[2][1] = 1;
		ret[2][2] = 1;
		
		ret[3][0] = 1;
		ret[3][1] = 1;
		ret[3][2] = 1;
		
		ret[4][0] = 0;
		ret[4][1] = 0;
		ret[4][2] = 1;
	}
	else {
		fprintf(stderr,"Error! Unknown analtype %d in SetKeys. Exiting....\n",aType);
		exit(1);
	}
	
	return(ret);
}


/* functions to return states given 0-26 */
int Xstate(int a, analtype aType)
{
	if(aType==MAPATRI) {
		switch(a) {
			case(0): return(0); break; 
			case(1): return(0); break; 
			case(2): return(2); break; 
			case(3): return(0); break; 
			case(4): return(0); break; 
			case(5): return(2); break; 
			case(6): return(1); break; 
			case(7): return(1); break; 
			case(8): return(3); break; 
			case(9): return(4); break; 
			case(10): return(0); break; 
			case(11): return(0); break; 
			case(12): return(0); break; 
			case(13): return(0); break; 
			case(14): return(0); break; 
			case(15): return(0); break; 
			case(16): return(0); break; 
			case(17): return(4); break; 
			case(18): return(3); break; 
			case(19): return(1); break; 
			case(20): return(1); break; 
			case(21): return(2); break; 
			case(22): return(0); break; 
			case(23): return(0); break; 
			case(24): return(2); break; 
			case(25): return(0); break; 
			case(26): return(0); break; 
			default:
				fprintf(stderr,"Error in Xstate().  a=%d.  not in [0,26]. Exiting...\n",a);
				exit(1);
		}
	}
	else if(aType==MATRI) {
		switch(a) {
			case(0): return(0); break; 
			case(1): return(0); break; 
			case(2): return(1); break; 
			case(3): return(0); break; 
			case(4): return(0); break; 
			case(5): return(1); break; 
			case(6): return(2); break; 
			case(7): return(2); break; 
			case(8): return(1); break; 
			case(9): return(2); break; 
			case(10): return(0); break; 
			case(11): return(0); break; 
			case(12): return(0); break; 
			case(13): return(0); break; 
			case(14): return(0); break; 
			case(15): return(0); break; 
			case(16): return(0); break; 
			case(17): return(2); break; 
			case(18): return(1); break; 
			case(19): return(2); break; 
			case(20): return(2); break; 
			case(21): return(1); break; 
			case(22): return(0); break; 
			case(23): return(0); break; 
			case(24): return(1); break; 
			case(25): return(0); break; 
			case(26): return(0); break; 
			default:
				fprintf(stderr,"Error in XState().  a=%d.  not in [0,26]. Exiting...\n",a);
				exit(1);
		}
	}
	else if(aType==PATRI) {
		switch(a) {
			case(0): return(0); break; 
			case(1): return(0); break; 
			case(2): return(2); break; 
			case(3): return(0); break; 
			case(4): return(0); break; 
			case(5): return(2); break; 
			case(6): return(1); break; 
			case(7): return(1); break; 
			case(8): return(1); break; 
			case(9): return(2); break; 
			case(10): return(0); break; 
			case(11): return(0); break; 
			case(12): return(0); break; 
			case(13): return(0); break; 
			case(14): return(0); break; 
			case(15): return(0); break; 
			case(16): return(0); break; 
			case(17): return(2); break; 
			case(18): return(1); break; 
			case(19): return(1); break; 
			case(20): return(1); break; 
			case(21): return(2); break; 
			case(22): return(0); break; 
			case(23): return(0); break; 
			case(24): return(2); break; 
			case(25): return(0); break; 
			case(26): return(0); break; 
			default:
				fprintf(stderr,"Error in XState().  a=%d.  not in [0,26]. Exiting...\n",a);
				exit(1);
		}
	}
	else if(aType==TRI) {
		switch(a) {
			case(0): return(0); break; 
			case(1): return(0); break; 
			case(2): return(1); break; 
			case(3): return(0); break; 
			case(4): return(0); break; 
			case(5): return(1); break; 
			case(6): return(1); break; 
			case(7): return(1); break; 
			case(8): return(1); break; 
			case(9): return(1); break; 
			case(10): return(0); break; 
			case(11): return(0); break; 
			case(12): return(0); break; 
			case(13): return(0); break; 
			case(14): return(0); break; 
			case(15): return(0); break; 
			case(16): return(0); break; 
			case(17): return(1); break; 
			case(18): return(1); break; 
			case(19): return(1); break; 
			case(20): return(1); break; 
			case(21): return(1); break; 
			case(22): return(0); break; 
			case(23): return(0); break; 
			case(24): return(1); break; 
			case(25): return(0); break; 
			case(26): return(0); break; 
			default:
				fprintf(stderr,"Error in XState().  a=%d.  not in [0,26]. Exiting...\n",a);
				exit(1);
		}
	}
	else {
		fprintf(stderr,"Error, unknown analysis type %d in Xstate. Exiting...\n",aType);
		exit(1);
	}
	exit(1);
	return(-99999.99);
}


/*! returns the index of the population having name Str.  If Str is not a known population
then it returns -1 */
int IdxFromPopName(char *Str, pfr_forback_data *D) 
{
	int i, ret = -1;

	for(i=0;i<D->NumPops;i++) {
		if(strcmp(Str,D->Pops[i]->Name)==0) {
			return(i);
		}
	}
	
	return(ret);
}

