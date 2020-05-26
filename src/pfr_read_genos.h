/* Here are some things to define to get some verbose output for debugging
 
	VERBOSE_FIRST_PASS_THROUGH_DATA     print out all the stuff that can be spit out on the first pass through the data file 
	VERBOSE_SECOND_PASS_THROUGH_DATA    just what it sounds like.  This also controls verbosity output in functions called from
										within CollectDataOnSecondPass, like CollectExtraColsOfPopInds, CollectExtraColsOfPopInds,
										SetUpPFR_Parent, etc.
 
	VERBOSE_PRINT_FIRST_PASS_SUMMARY	after doing the first pass, print out the summary of the first pass.
 
*/
#ifdef VERBOSE_ALL
#define VERBOSE_FIRST_PASS_THROUGH_DATA
#define VERBOSE_SECOND_PASS_THROUGH_DATA
#define VERBOSE_PRINT_FIRST_PASS_SUMMARY
#endif

//#define VERBOSE_FIRST_PASS_THROUGH_DATA
//#define VERBOSE_SECOND_PASS_THROUGH_DATA
//#define VERBOSE_PRINT_FIRST_PASS_SUMMARY



 /* max number of times an individual can be regenotyped */
#define MAX_REGENO_TIMES  10 

#define MAX_INDIV_NAME_LENGTH 100
#define MAX_POP_COLL_NAME_LENGTH 100

#define MAXPOPS  2000
#define MAXOFFCOLLS 2000

#define MAX_DISTINCT_ALLELES 5

#define MAX_NUM_VALS_IN_RANGE  500
#define EXCOLMAX 20

#define POPCOLSTART 5
#define OFFCOLSTART 8
#define NUMKEYWORDS 11
typedef enum {
	PREAMBLE, /* this is not a keyword, but is does let me know if we are still in the preamble as a "regime" */
	POP,
	OFFSPRING,
	MISSING_ALLELE,
	CHINOOK_AVE_POP_SIZE,
	POPCOLUMN_SEX,
	POPCOLUMN_REPRO_YEARS,
	POPCOLUMN_SPAWN_GROUP,
	OFFSPRINGCOLUMN_BORN_YEAR,
	OFFSPRINGCOLUMN_SAMPLE_YEAR,
	OFFSPRINGCOLUMN_AGE_AT_SAMPLING,
} KeywordStatus;

typedef enum {
	MALE,
	FEMALE,
	SEX_UNKNOWN
} SexEnum;


/* for telling us whether a name of a collection refers to a parental collection or an offpring collection */
typedef enum {
	POP_NAME,
	OFF_COLL_NAME
} coll_type_name_enum;

/* for referring to the alleles 0, 1, or missing */
typedef enum {
	ALLELE_0,
	ALLELE_1,
	ALLELE_MISSING
} AllelicTypeEnum;




struct indiv_name_hash_struct {
    char name[MAX_INDIV_NAME_LENGTH];             /* key */

	int AbsIdx;  /* a little something to record the absolute index of this individual */
	
	int parent_pop;  /* if it is found in a POP collection, this tells us which one */
	int offspring_coll;
	
	/* and these tell us the index in those collections*/
	int idx_in_pop;   /* note that  the position of the individual will be [parent_pop][Sex][idx_in_pop] */
	int idx_in_offcoll;  /* the position of this indiv will be [offspring_coll][idx_in_offcol] */
	/* these two are for counting how many times the genotype is given */
	int cnt_par;  /* number of occurrences amongst the parents */
	int cnt_offs; /* number of occurrence amonst the offspring.  */  
	
	SexEnum  Sex;  /* we are going to check to make sure that the sex is the same each time */
	
	
	
    UT_hash_handle hh;         /* makes this structure hashable */
};


struct pop_and_coll_name_hash_struct {
    char name[MAX_POP_COLL_NAME_LENGTH];             /* key */
	coll_type_name_enum CollType; /* tells us whether the name is for an offspring collection or a parental (POP) group */
    int idx;  /* the idx given the pop/coll name */
    UT_hash_handle hh;         /* makes this structure hashable */
};




/* here are a few things for keeping track of spawning group names.  Spawning groups can be specified 
 in a hierarchical fashion so that there can be some missing data in the lower levels, but you can stil
 identify possible mates at an upper level.  For example: you know that fish X spawned on day Y but you
 don't know the exact pair.  Then, eventually you will be able to specify the spawner group as:
 Y|? or something like that.  Or  11/13/09$GroupSix  or something like that.  At first I am not going
 to bother parsing the multiple possible hierarchies, but I will set things up so that we can do that
 eventually.  
 
 the important part is that we will be hashing strings to ints in the context of different hierarchical 
 levels.  Note that ints specifying spawner groups in any particular population might not be in series
 because we just want to get an integer (any integer) that corresponds to a string.
 */
#define MAX_SPAWN_GROUP_LEVELS  10 
#define MAX_SPAWN_GROUP_NAME_LENGTH  500
#define MAX_NUM_SPAWN_GROUP_NAMES_WITHIN_LEVELS 5000
struct spawn_group_name_hash {
	char string[MAX_SPAWN_GROUP_NAME_LENGTH];
	int *ID;  /* the integer that corresponds to the string, which depends on the level (so is indexed by level) */
	UT_hash_handle hh;
};


/* Here I define a few things for handling individuals that can be in multiple spawner groups.
 * I am going to give up on the whole hierarchical levels of spawning groups, as that ends up not
 * being as useful as having individuals capable of being in multiple spawner groups.
 */
#define MAX_NUM_SPAWN_GROUP_PER_INDIV 50




/* this is just for simply keeping track of things like the year in which different individuals reproduced
 or might have been born, etc */
typedef struct {
	int Lo;  /* The lowest year in the range */
	int Len; /* the length of bits.  Note that Lo+Len-1 is the maximum year in the range */
	char *bits;  /* an array of 0's or 1's.  if bits[i] is 1 then the year Lo+i occurs in the range */  
} pfr_year_range;



typedef struct {
	char *Name; /* the identifier given to this individual in the data file */
	int Pop;  /* the index of the population it belongs to */
	int Idx;  /* the index of the individual within the population (also depends on sex)  */
	int AbsIdx;  /* the absolute index of this individual's name (so, we will hash an individual's name and keep an
	 index for all POP and OFFSPRING individuals that way so we can deal with duplicates in the data set.
	 So, we can compare an individuals identifier by comparing AbsIdx's which should be way faster than using
	 strcmp on their string names. */
	
	SexEnum Sex; 
	pfr_year_range *ReproYears;
	int NumSpawnGroups;  /* the number of different spawning groups an individual belongs to */
	int* SpawnGroup;  /* denotes the different spawning groups (recorded as integers) that an individual can belong
	                     to*/
	
	int NumGenos;  /* the total number of genotyped samples of this individual in this context.  */
	char **geno;  /* indexed by i,l where i is the regenotyping number of l is the locus. Possible values are 
	 0 -> no copies of the ``1'' allele 
	 1 -> one copy of the 1 allele
	 2 -> two copies of the 1 allele 
	 3 -> missing data   */
	
	int *NumMissingLoci;  /* number of missing loci in the individual. Subscripted by the ReGeno number */
	
} pfr_parent;



/* this is a little structure to hold some information about particular parent pairs that form a 
 non-excluded trio with the offspring individual that will be carrying this structure around.
 */
typedef struct {
	pfr_parent *mamapapa[2];  /* subscripted  by MALE or FEMALE.  */
	
	int pop;  /* the population they come from */
	
	/* records the number of Mendelian Incompatibilities between the mother, the father, and in the trio as a whole */
	int ma_MI;
	int pa_MI;
	int trio_MI;
	
	int *A_64_idx;  /* an array to store the genotype of the trio formed by this parent pair and the youth--an int from  0 to 63 which allows for missing data and is
					   congruent with the AProbs_M arrays) at each locus. */
	
	int *X_states_64;  /* an array to hold the X states of each locus in this trio formed by this parent pair and youth */
	
	double *TrioPosteriors;  /* an array of length NUM_SPEC_PEDS that will hold the posterior prob for each PED type */
	double LogL;  /* holds the log-likelihood ratio of parental versus non-parental.  Gets calculated at the same time as TrioPosteriors */
	
	double PrPar;  /* prob of the trio genotype under the parental hypothesis */
	double PrNonPar; /* prob of the trio genotype under the non-parental hypothesis */
	
	double MaxPosterior;  /* the posterior for the trio relationship that has the highest posterior for this pair */
	int MaxPosteriorRelat;  /* the indes of the relationship that has highest posterior probability */ 
	
} ParentPairMatch;



/*
 This is a little structure that an individual will carry around with him/herself
 which will store pointers to individuals and pairs of individuals that satisfy
 the Mendelian incompatibility criterion we have set up.
 
 */
typedef	struct {
	pfr_parent **Parents[3];  /* A 2-d array of pointers to parents.  Subscripted by SexEnum (MALE,FEMALE,SEX_UNKNOWN) and then
	 by the idx of the parent within that sex */
	int NumParents[3];  /* number of parents in each sex category */
	int *Num_MI_Parents[3];  /* number of loci with Mendelian incompatibilities with each parent. Subscripted by sex and then the idx of the parent. */
	
	ParentPairMatch *ParPairs;
	int NumParPairs;  /* this will reflect the number of unexcluded parents that go forward for further evaluations.  Note that if pairs get
                       tossed due to low LogL, this number will reflect that still. */
  int NumParPairsMendelian;  /* this records the raw number or parent pairs with few enough Mendelian incompatibilities. */
  int NumParPairsMend_LogL;       /* This records the number of parent pairs left after applying the MinLogL criterion to the Mendelian compatible ones */
  int NumParPairsMend_LogL_and_Rank; /* This records the number of parent pairs left after Mendelian, LogL, and max-number parent pairs allowed. */ 

	
} indiv_matchers;





typedef struct {
	char *Name; /* the identifier given to this individual in the data file */
	int OffColl;  /* the index of the offspring collection it belongs to */
	int Idx;  /* the index of the individual within the offspring collection. */
	int AbsIdx;  /* the absolute index of this individual's name (so, we will hash an individual's name and keep an
	 index for all POP and OFFSPRING individuals together so we can deal with duplicates in the data set.
	 So, we can compare an individuals identifier by comparing AbsIdx's which should be way faster than using
	 strcmp on their string names. */
	
	int SampledYear;  /* year it was sampled */
	
	int NumBornYears;  /* number of differnt years when it might have been born */
	int *BornYears;  /* the years when it might have been born.  This is either supplied by the user or computed from 
						SampledYear and the AgesAtSampling */
	int NumAgesAtSampling;  /* how many possible ages could it be at sampling */
	int *AgesAtSampling;  /* may get data read into it, but this should be immediately translated into BornYears which should
									actually then get used by the program.  */	
	
	int NumGenos;  /* the number of regenotyped samples of this individual */
	char **geno;  /* indexed by i,l where i is the regenotyping number and l is the locus. Possible values are 
	 0 -> no copies of the ``1'' allele 
	 1 -> one copy of the 1 allele
	 2 -> two copies of the 1 allele 
	 3 -> missing data   */
	
	int *NumMissingLoci;  /* number of missing loci in the individual. Subscripted by the ReGeno number */	
	
	indiv_matchers *MendComps;  /* pointer to struct that holds info about whom this individual is mendelian-compatible with (given our smax threshold) */
	
	ival *MargAlphas;  /* for gathering information about how many times simulated LODS exceed this individual's LOD from marginal simulation */
	
	double *PValueParentalCDF;  /* the i-th element of this records the fraction of simulated parental trios conditional on this offspring that had a LogL which 
	                               was less than the LogL at which i% of the non-parental max LogL's are smaller.  And the 0 element records the fraction of 
									simulated parental pairs with LOD less than all the simulated nonparental pairs.  */
	double p_value;  
	double FDR; /* after the FDR calculations, this will hold the FDR value for this individual */
	
} pfr_offspring;


typedef struct {
	int NumLoci;
	char MissingAlleleString[50];
	char **LocusNames;
	int NumPops;
	int NumOffColls;
	char **PopNames;
	char **OffCollNames;
	
	int **NumInPops; /* subscripted by Pop Idx, then by Sex */
	int *TotInPops; /* the total number of indivs in each pop */
	int *TotRowsInPops;  /* the total number of rows in each pop */
	pfr_parent ***Pops;  /* subscripted first by Pop Idx, then by Sex (0,1,2) and then by idx within pop-sex */
	
	int **OffCollPossPops;  /* array of 0's and 1's that tell us whether the offspring in a collection could be children of the corresponding parent Pops.  This 
							is subscripted by OffColl then by POP */
	int *NumInOffColls;
	int *TotRowsInColls;  /* for the total number of rows in each offspring collection */
	pfr_offspring **Offs;  /* subscripted by OffColl then by idx with OffColl */
	
	int AbsNum;  /* for the absolute number of individuals in the data set */
	
	/* this hashes indiv names so we know if we have seen them yet or not. */
	struct indiv_name_hash_struct *IndivNameHash;
	
	/* this hashes the collection names so we can check for uniqueness */
	struct pop_and_coll_name_hash_struct* PopCollNameHash;
	
	
	/* these are for before we might drop some loci (we aren't dropping any at this point, though) */
	int **AlleNamesRaw;  /* array that is subscripted by [Locus][allele] which gives the int that is each allele at each locus */
	int *NumAlleRaw;  /* number of alleles at each locus.  This should be 2!  Set it on the first pass and discard those that 
					  are monomorphic or which have more than 3 alleles.  Ooops. I currently don't discard monomorphic loci.  */  
	
	
	double ***AlleleCounts;  /* subscripted by POP idx, then Locus, then AllelicTypeEnum.  For alleles 0 and 1 it counts the number of occurrences of each 
							allele, and for 3 (=missing data) it counts the number of loci at which there are missing data. */
	double ***AlleFreqs; /* subscripted by Pop idx, then Locus, then 0 or 1 (and now can be subscripted by 2, which holds the frequency of missing data at this locus in the population).
	                   For allele freqs these are just the posterior means, assuming a Beta(1/2,1/2) prior,
						of the allele freqs in each POP. The missing data are discarded for this calculation, and the alle freqs at a locus in a pop
						sum to one. */
	
	/* this is some stuff to record about the input file */
	KeywordStatus *ExtraPopCols;  /* this is the type of each of the extra columns in order */
	KeywordStatus *ExtraOffCols;
	int *KeywordFlags;  /* gets a zero if the keyword has not used and a >0 otherwise */
	int NumExtraPopCols;
	int NumExtraOffCols;
	
	
	/* stuff for turning spawner group names into integers */
	char ***SpawnGroupNames;  /* indexed by the level in the hierarchy and then number within that level */
	int *TopSpawnGroupInt;  /* gets pre-fix incremented with each new spawn group name at a level */
	struct spawn_group_name_hash *SpawnGroupNameHash;
	
	double *GtypErrRates;  /* assumed GtypErrorRates assumed for the analysis */
	
	int *AveragePopSizes;  /* for storing the average pop sizes in the data file when the CHINOOK_AVE_POP_SIZE method is used */
	
} pfr_geno_data;







pfr_geno_data *FirstPassThroughData(const char *FileName);
int *pfr_StringToArrayOfInts(char *Str, int MaxNums, int *length);
void PrintSummaryOfInputData(FILE *s, pfr_geno_data *P);
void SummarizeAllelicTypes(FILE *s, pfr_geno_data *P);
void CountAlleles(pfr_geno_data *P);
void PrintSummaryOfAllelicCounts(pfr_geno_data *P);
void PrintSummaryOfAllelicCountsAndFreqs(FILE *s, pfr_geno_data *P);
void ComputeAlleleFreqsFromCounts(FILE *s, pfr_geno_data *P);
void SummarizeLocusNameAndGtypRates(FILE *s, pfr_geno_data *P);
int ReadAndCheckPopSize(FILE *in, char *S);
void CollectDataOnSecondPass(pfr_geno_data *ret, const char *FileName);
void TokenizeSpawnerGroupString(char *str, int *num_toks_out, char tokens[MAX_NUM_SPAWN_GROUP_PER_INDIV][MAX_SPAWN_GROUP_NAME_LENGTH]);
