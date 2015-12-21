



/* here we define some structures that will hold information about the indivdiuals
in each offspring collection */

typedef struct {
	int NumPops;  /* just the number of parental pops, copied over here */
	int NumOffs;  /* the number of offspring in the offspring collection under consideration */
	int OC_idx;   /* the index of the collection.  This can be good to have for debugging purposes */
	
	double *pi;  /* mixing proportion of the different populations.  Length=NumPops */
	double *pi_prior;  /* parallel vector of prior parameters for the Dirichlet rv */
	
	double *phi;  /* if you are an offspring from pop x, this records the probability 
					  that your parent pair is in the data base.  Length=NumPops  */
	double *phi_prior;
	
	double **theta; /* fraction of trios of different types that are NOT parental.
					   indexed by pop and then by pedigree type.  Length = NumPops 
					   by NumSpedPeds */
	double **theta_prior;  /* will use the same prior for all pops, but want to have
							the flexibility to change that later, so we will store a 
							separate direchlet parameter vector for each */
	
	int *setZ;  /* vector of 0's and 1's of length NumPops.  1 means the pop is in the
				   set of OffCollPossPops, 0 means it isn't.  Essentially, this is just
				   a row of OffCollPossPops.  */
	
	int *setW;  /* a vector of length NumPops.  >0 denotes that there is
				at least one pair of parents from that population that is non-excluded
				from parenthood (i.e. is in S-down) in at least one individuals in the
				offspring collection under consideration */
	
	int *NumPairs;  /* for each offspring in the collection, the number of non-excluded
						putative parent pairs.  Length=NumOffs */
	
	int **W;	/* array of indicators telling us the population of origin of different pairs.
				   Subscripted by:
					"idx"		the index of the offspring in the collection (0,...,NumOffs-1)
					"p"			the index of the pair of putative parents non-excluded to indiv idx.
								this will be (0,...,NumPairs[idx]-1).
					The values of this array will be the population indices (0 to NumPops-1) */
	
	int ****V; /* big array of zeros and ones that tell us about putative parental individuals
				  shared between non-excluded pairs of each offspring individual in the 
				  offspring collection under consideration.  Subscripted by:
						"idx"    the index of the offspring in the collection (0,...,NumOffs-1)
						"p1"	 the index of the first pair in the comparison (0,...,NumPairs[idx]-1)
						"p2"	 the index of the second pair in the comparison (0,...,NumPairs[idx]-1)
						"sex"    0=MALE, 1=FEMALE (a la the SexEnum).  Will figure out how to deal with 
								 sex-unknown individuals at some point, but not at this time. (0,1) 
				 So, for example:  V[idx][p1][p2][0] == 1 means that the putative father in parent pair p1 and p2
				 of offspring idx in the collection under consideration is the same individual.  Note that this
				 big array will be symmetrical in p1 and p2.    */
	
	double ****P;  /* probabilities of the genotypes of all the members of a trio made up of indiv "idx" and pair
					 "p", conditional on being in S-down, and given trio relationship "r".  Subscripted by:
						"idx"				the index of the offspring (0,...,NumOffs-1)
						"p"					the index of the pair (0,...,NumPairs[i]-1)
						"cross_pop?"		a 1 if the offspring and the parents are in different populations 
											and a 0 if they are in the same. (0,1)
						"ps_p"				if cross_pop==0 this is the index of the trio type in the ped_spec array.
											if cross_pop==1 this is the index of the possible population of origin of 
											the individual.  (0,...,NumPops-1) but nonzero only for those indices also
											in setZW.   &*/
	
	double **PairlessP;  /* probability of the genotypes of the offspring, irrespective of whether they have non-excluded parent pairs.
							this is just the genotype prob of each offspring, and it will be used to do the mixture stuff for 
							offspring with no non-excluded parents.  Subscripted by:
								"idx"		the index of the offspring  (0,...,NumOffs-1)
								"pop"		the index of the population (0,...,NumPops-1) 
							*/
	
	int **U;  /* indicators telling us whether the putative father or putative mother of each pedigree type is 
				 equal to self.  Subscripted by:
					"ped"		the index of the specified pedigree (0,...,NumSpecPeds-1)
					"sex"		MALE=0  Female=1.    (0,1)  
				these values will not change during the course of the program, but I want to keep them local
				in this struct anyway.  */
	
	
} trio_mixture_offcoll_struct;



/* function prototypes */
int **ReturnPedSpecSelfIndicators(void);
trio_mixture_offcoll_struct *PrepareTrioMixtureOffCollStruct(pbt_high_level_vars *HLV, int OC_idx);
double ProbOfSingleLocusGenotypeDataOfTrio(char kid, char pa, char ma, double *Aprobs);
double CondProbOfMultiLocusGenotypeDataOfTrio(pfr_offspring *kid,		
											  int kid_pop,					/* the population of the kid.  Usually not known, but given to compute necessary probs */
											  pfr_parent *pa,
											  pfr_parent *ma,
											  int the_ped,					/* note that this is ignored if kid_pop is not the same as the pop of the ma and pa */
											  pbt_high_level_vars *HLV );
double SimpleIndivGenotypeProb(char *G, int pop, pbt_high_level_vars *HLV);
void SuperNaiveSimpleTrioEM(trio_mixture_offcoll_struct *T);

