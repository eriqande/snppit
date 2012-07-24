/*
 *  forback_abstracted.h
 *  pfr
 *
 *  Created by Eric C. Anderson on 3/9/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */


/* Here we will define some data structures for dealing with the forback_calculation
   in a suitably abstracted way so that it can be defined for all sorts of different
   situations.
   
   I have some writings on this in:
   
   /Users/eriq/Documents/work/prj/PFR/doc/forback_abstracted.tex

	That might be useful to consult, as it probably influences the choice of variable names
	etc.

*/


#define MAX_FORABS_DIMS 10



/* these are parameters that are constant across all the different sets of
	probabilities of A-states we might have */
typedef struct {
	int T;  /* number of time steps (i.e. number of loci) */
	int NA;  /* number of states of each A_t.  */
	int NX; /* total number of sets X_0,...,X_{NX-1}, that partition the integers from 0 to NA-1 */
	
	int *CX; /* the cardinality of each of the NX sets X_k (i.e. the number of elements in each one */
	int **X; /* the actual integers (which are indexes of the A-states) contained in each set X */
	
	int d; /* the length of the v-vectors of 0's and 1's */
	int **v;  /* an NX by d matrix of 0's and 1's.  v[i][j] is the j-th element of the v vector corresponding
	             to set X_i (all base-0 subscripted) */
	
	int *Smax; /* array of length d giving the s^max vector */
	
} forabs_univ;


/* function prototypes */
void ComputeXProbs( forabs_univ *FU, double **PA, double **PX);
int S_status_plus( int *s, int *v, int t, int *Smax, int *newS, int d);
int S_status_minus( int *s, int *v, int t, int *Smax, int *newS, int d);


double **ForwardStep_One_d(forabs_univ *FU, double **PA, double **PX, double *Over  );