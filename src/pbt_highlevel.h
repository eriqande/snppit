/* here are some macros to define to get more verbose output
 
 VERBOSE_SINGLE_PARENT_COMPAT_WITH_OFFSPRING				print the parent pairs compatible with each offspring
 VERBOSE_NUM_NON_EXC_PAIRS_BY_POP							execute all the print statements at the end of RecordNonExcParentPairsFromPop
 
*/

#ifdef VERBOSE_ALL
#define VERBOSE_SINGLE_PARENT_COMPAT_WITH_OFFSPRING
#define VERBOSE_NUM_NON_EXC_PAIRS_BY_POP
#endif

//#define VERBOSE_SINGLE_PARENT_COMPAT_WITH_OFFSPRING

//#define VERBOSE_NUM_NON_EXC_PAIRS_BY_POP




/* these are data structures and prototypes for the highest level functions
and for the options that the user will actually pass to the programs, 
etc.
*/



#define MAX_FILE_PATH_LENGTH 5000
#define MAX_NUM_NONEXC_PARENT_PAIRS 100000


/* to return the index of the collection in PurePopTrioColls given then PedSpecIdx and 
 the population.  This calculation requires the Number of PedSpecIdx's which is also an
 argument to the macro */
#define PPTC_IDX(NumPedSpec,  PedSpecIdx,  PopIdx)  ( PopIdx * NumPedSpec + PedSpecIdx )
#define CPTC_IDX(NumPops,  PopIdxParent,  PopIdxKid)  ( PopIdxParent * NumPops +  PopIdxKid)


/* at the end of all the parentage assignment, we are going to want to have a list of indiduals with 
 the max-LOD parent from eac individual.  We are going to store that in this structure.  We will want
 to be able to access it in a number of different ways, which is why it has this sort of structure.
*/
typedef struct {
	int n;  /* number of offspring with max-lod, non-excluded parents from this population */
	pfr_offspring **inds;  /* array of pointers to all the individuals with max-lod nonexcluded parents from this population */
	
	int NumNonExcCnts; /* amongst all these offspring, the number of unique numbers of non-excluded pairs to each kid */
	int *NonExcCnts;  /* array recording the number of non-excluded pairs. length of NumNonExcCnts.  The values in this should
	                    be in sorted order. */
	int *NumKidsWithNonExc;  /* the number of offspring with max lod parents from this population that have a certain number of non-excluded parents.
								this thing is parallel to NonExcCnts. */
	pfr_offspring ***inds_by_cnt; /* suscripted by i,j where i is the index of the number of non-excluded parents (the number of 
					non-exc parents is found in NonExcCnts[i]), and j is between 0 and NonExcCnts[i]-1. */
	
	double **BetaArraysByHundredths;  /* for keeping track of how often a true parental trio would have lower LOD than the max of the 
										simulated nonparental ones.  Indexed first by number of non excluded parent pairs, and then by 
										0,...,99 which corresponds to alpha cutoffs of 0, .01, .02, ... ,99 */  
	
	double *ParentPValueDensityHistogramEstimate;  /* the 0 element of this is the estimated density of the portion of the distribution that is Pvalues greater
	 than all of those observed in the non parentals.  This corresponds at to the range of P-values between 0 and 1/NumMcReps.
	 The 1 element is for that portion of the distribution between 1/NumMC_Reps and .01; the 2 element is for that
	 portion between .01 and .02, and so forth.  The 100 element is for that portion between .99 and 1.0. */
	
} inds_with_max_lod_parents_from_this_pop;



typedef struct {
	char *DataFileName;
	int *big_smax;  /* the smax to use when computing false negative rates.  This will typically
						 be higher than the smax applied later on.  This might be 5 5 8 for example */
	int d_from_big_smax;  /* the number of components of smax as given on the command line */
	
	double MI_fnr_target;
	
	int MaxAllowMissingLociInParents;  /* the maximum number of missing loci allowed for a parent.  Default is 10 */
	
	double **Pi;
	
	int DryRun;
	
	int PszForAll;  /* user may specify an average "chinook" popsize to coerce upon all the pops.  Initialize to -1 */

} pbt_user_opts;


/* this is a little structure for capturing information about individuals with max-LOD non-excluded parents from a particular population.  It will be made sortable. */
typedef struct {
	pfr_offspring *ind; /* a pointer to the offspring associated with this data */
	double p_value;   /* the p-value associated with this offspring's max LOD parent pair */
	double *parental_cdf;
	double approx_p_dens;  /* the approximate density value corresponding to the p-value */
	double *parental_p_value_dens;
} pbt_fdr_container;
	
	
typedef struct {
	pfr_geno_data *PFR;
	pbt_user_opts *PBUO;
	int *smax_to_use;  /* the smax that we decided to use (because it gives a reasonably
						low false negative rate)  */
	
	
	
	/* these are for storing the Trio Forward-backward structures */
	/* here is an array to store the results for only parental trios 
	 at the largest s-max the user wants to consider (smax_initial---typically 5 5 7 or something
	 like that, unless genotyping error rates are totally out of hand). It is 
	 indexed by population. */
	FB_Vars *ParentalTrioBigSmax;
	
	/* here is a pointer for storing the collections for trios that consist entirely of 
	 individuals from a single population. Inside it, the Colls field is indexed by [ped_spec][population],
	 but that has to get translated into a single array of pointers, so it is 
	 ped_spec_idx * num_pops + pop_idx.  */
	FB_Vars *PurePopTrioColls;
	
	/* and here is a pointer to store the collections that are the cross-bred trios,
	 i.e. the parents are from pop i and the kid is from pop j.  The Colls within these are indexed
	 by [i][j], and, of course, the pedigree specified for all of them is of type 
	 C_U_U.  The 2-D index gets translated into a 1-D index as i*num_pops + j.  So,
	 [i][i] is parents and youths from the same population i, and the "off-diagonals" 
	 are the cross-bred trios.   */
	FB_Vars *CrossPopTrioColls;
	
	
	
	/* here is a thing to keep track of the groups of individuals with non-excluded max-LOD parents from each 
	 of the different populations */
	inds_with_max_lod_parents_from_this_pop **KidsWithMaxLOD_Parents;  /* indexed by population */
	
	
	/* some file pointers for output */
	FILE *MaxLodNonExpPar_File;
	FILE *SingleParentCompats_File;
	FILE *TrioPosteriors_File;
	FILE *BasicSummaries_File;
	FILE *FDR_Summaries_File;
	
	
} pbt_high_level_vars;



/* this is just a simple little struct that I can use to sort simulated values 
   and retain the parent pair number of pedigree type that they originated from 
*/
typedef struct  {
	double LogL;  /* the log likelihood ratio of parental over non-parental */
	double ParPosterior;  /* the Posterior prob of parental */
	double PostVec[NUM_SPEC_PEDS];
	int ped;  /* the index of the type of trio simulated from */
	int parpair;  /* the index of the non-excluded parent pair from which the missing data mask was obtained */
	
	double PrPar;  /* the prob of *simulating* the genotype of the trio under the parental hypothesis.  */
	double PrNonPar; /* the prob of *simulating* the genotype of the trio under the non-parental hypothesis */
	
} pbt_simmed_values_struct;


/* prototypes */
pbt_user_opts *GetPBT_UserOpts(int argc, char *argv[]);
void SelectAnSmax3(pbt_high_level_vars *HLV);
double SumOfS_states3(FB_Vars *A, int *s,  int C);
void ComputeCrossPopTrioColls(pbt_high_level_vars *HLV);
void AssignMatchingSingleParents(pbt_high_level_vars *HLV, int OC_idx);
void AssignMatchingParentPairs(pbt_high_level_vars *HLV, int OC_idx, int MaxAllowable);
void CalculateTrioPosteriorProbs(pbt_high_level_vars *HLV, int OC_idx);
void AssessFPR_and_FNR_ByBackwardSimulation(pbt_high_level_vars *HLV, int OC_idx, FB_Vars *CPFB);
inds_with_max_lod_parents_from_this_pop *RecordNonExcParentPairsFromPop(int Pop, pbt_high_level_vars *HLV);
void MarginalSimForFDR_Calcs(pbt_high_level_vars *HLV, int ThePop, FB_Vars *CPFB);
void GrabRandomGenos(pbt_high_level_vars *HLV, int ThePop, char **kid_geno, char **pa_geno, char **ma_geno);
int cmp_simmed_values(const void *a, const void *b);
int cmp_fdr_containers(const void *a, const void *b);
int cmp_doubles(const void *va, const void *vb);
double *EmpiricalCrossCDF100(double *x, double *y, int lx, int ly);
void DoFDR_ForAPop(int ThePop, pbt_high_level_vars *HLV);
double EstimateFractionParentalByEM_UsingApproxDensOfPvalues(pbt_fdr_container *f, int N, double NullDens, double Precis, int MaxIts);
double *Return_ArchTypeChinookPiVec(int AVE, int *outx, int *outpsz);
void NegotiatePiVectors(pbt_user_opts *PBUO, pfr_geno_data *PFR);
