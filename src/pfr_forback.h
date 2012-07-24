/*
 *  pfr_forback.h
 *  pfr
 *
 *  Created by Eric C. Anderson on 8/2/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
 
 
#define MAX_POP_NAME 100
#define MAX_NUM_POPS 2000
#define MAX_NUM_MCOPTS 2000
#define MAX_NUM_LAMOPTS 2000
#define MAX_XSTATES 5


/*! these enums give the type of analysis.  The type of analysis determines a lot of 
different parameters and variable values.*/
typedef enum {
	TRI,  /*!< only trios have been screened */
	MATRI, /*!< first mother pairs, then trios were screened */
	PATRI,  /*!< first father pairs, then trios were screened */
	MAPATRI /*!< mother pairs, then father pairs, then trios were screened */
} analtype;



typedef struct {
	double *PA;  /*!< frequencies of the 26 trio states which the user will input (Prob of A, where A is the genotypic state of trio) */
	double *PX; /*!< freqs of the different mismatch states (X states).  Length 5 if MaPaTri, 3 if MaTri or PaTri, and 2 if Tri.  That length stored in XNum. */
	double **PAinX; /*!< freqs of the different A's normalized to sum to one within within each of the X-states.  Size is XNum by NumAinX[i] */
} pfrfb_locus;



typedef struct {
	int Idx; /*!< The index of the population */
	char Name[MAX_POP_NAME]; /*!< the name of the population */
	pfrfb_locus **Loci;  /*!< array of pointers to locus information for this population */
	double ****dim3PY; /*! for MAPATRI, this stores probabilities subscripted by locus,Y1,Y2,Y3 */
	
	/* eventually, this is where I will have a pointer to the results of a forward calculation */
} pfr_pop_struct;


/*! This structure holds parameters relating to the types of monte carlo exercises 
we would like to do. */
typedef struct {
	int LN; /*!< Index of population for numerator of lambda */
	int LD; /*!< Index of population to serve as denominator of lambda */
	
	int T; /*!< Index of the "True" population, i.e. that which we want to simulate genotypes from */
	int IS; /*!< Index of the population that should be the importance sampling distribution.  If -1, then we just will do vanilla Monte Carlo */
	
	double LC;  /*!< Lambda_c -- i.e. the critical value of lambda. */					
	int N;  /*!< The number of importance sampling reps */
	
	int FullDist; /*!< 0 ==> don't print out the lambdas.  1 ==> print out the lambdas */
} mc_struct;


/*! This is a structure for holding information about how to compute lambda
in the pfr_bits world.  It is just like the mc_struct, but we only retain
the parts that we need---since there is no simulating, we just need LN, LD, and
LC */
typedef struct {
	int LN; /*!< Index of population for numerator of lambda */
	int LD; /*!< Index of population to serve as denominator of lambda */
	double LC;  /*!< Lambda_c -- i.e. the critical value of lambda. */					
	
} lam_struct;


typedef struct {
	int L; /*!< number of loci */
	analtype aType;  /*!< what type of process has been followed in selecting 
						the trios for further analysis.  (i.e., did you look for compatible fathers, then compatible mothers,
						and then construct trios from them, etc.) */
	int Xdim; /*!< dimension of the X and Y variables.  This gets set according to the aType. It is 1,2, or 3. */
	int XNum; /*!< number of X states.  Depends on aType.  It is 5, 3, or 2. */
	int NumA;  /*!< The total number of A-states.  Typically 27 for no-missing data trios */ 
	int* NumAinX; /*!< for each X-state, this gives the number of A-states in it.  Depends on aType */
	int **IdxAinX; /*!< Index of the A-state corresponding to the elements in PAinX. Also depends on aType */
	int **Xkeys; /*!< array of length XNum. each element is an array of length Xdim which indicates incompatibility at a locus.
					For example, if aType==MAPATRI, then XKey[0] = 0,0,0, which means that there were no incompat loci
					in the pa-kid, ma-kid, or overall trio respectively.  All these values are set in function SetXKeys(). */
	int *Maxes; /*!< array of length Xdim giving the max value of each state to keep track of */
	 
	int NumPops; /*!< The number of different distributions of A-states we will be dealing with.  Typically
						these are going to arise from different breeding populations.  Hence the name */
	pfr_pop_struct **Pops;  /*!< an array of pointers to the different breeding pops */
	
	
	
	int NumMC;  /*!< The number of different sets of Monte Carlo experiments to do */
	mc_struct **MC;  /*!< array of pointers to structs holding information about what type of MC experiment to do */
	
} pfr_forback_data;


