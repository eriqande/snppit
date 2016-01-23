
/* these are some things to define either here or at compile time with the -D flag
   to give verbose output in different parts of the program for debugging and
   checking, etc
 
   VERBOSE_READING_FB_INPUT
 
*/

#ifdef VERBOSE_ALL
#define VERBOSE_READING_FB_INPUT
#endif

//#define VERBOSE_READING_FB_INPUT



#define MAX_COLLECTIONS 10000
#define MAX_COLLECTION_NAME_LENGTH 5000
#define MAX_LOCI 1000
#define MAX_X_STATES 400
#define MAX_A_STATES 400
#define MANY_MORE  300  /* weird little thing to make sure I don't run out of MULT_USE_OPTION instances in some places */



typedef struct {
	int Nx; /*!< number of X-states */
	int d; /*!< length of the v-vectors */
	int **v_vecs; /*!< array that is Nx by d giving values of v vectors */
	int L; /*!< the number of loci */
	int *smax; /*!< the s_max vector */
	
	int Na; /*!< the total number of A-states.  Assume they are indexed 0,...,Na-1 */
	int *AinX; /*!< array of length Na telling which X state each A belongs to */
	int *NaInX; /*!< array of length Nx telling us how many A states are in each X state */
	int **aInX2A; /*!< subcripted by [i][j] it gives the index of the A state that is 
	 the j-th A-state associated with the i-th X state*/
	int HasFakeProbs;  /*!< set to 1 if we have fake probs (as if we ignored missing data (i.e gave it a prob of 1.0)) to compute */
	
	int NY;  /*!< the number of genotype states underlying each A-state (3 in the case of trios) */
	int **Ystates; /* indexed by i,j.  Ystates[i][j] is the state of the j-th Y state underlying the i-the A state */
	int *ymaxes;  /* keep track of the max value taken in each of the y-states so we can make a quick hash to produce the A state from all the Y states
	              at least with no missing data at this point */
	int *AgivenY; /* pick an element of this using ReturnVdex to get to A from the Y's */
	
	int NumMissModes;  /* the number of patterns of missing data. Equal to 8 when dealing with trios.  */
	int **MissModes;   /* 2-d array indexed by i,j, with i =0,...,NumMissModes-1 and  j= 0,...,NY-1.  A 1 in any
						entry indicates that the j-th genotype is missing in the i-th missing data mode.  For trios,
						j=0 is kid; j=1 is pa;  j=2 is ma. */
	
	int MissModeFrom0sAnd1s[2][2][2];  /* stick a 1 or a 0 in if the y, pa, ma is missing data or not, and the value in the
											array element should be the MissMode.  i.e.  MissModeFrom0sAnd1s[0][0][0] = 0, etc.
										I just hard wired this because I am tired of all the extensible crap that isn't ever going to 
										be used anyway.  */
	int **MissModeCensoring;  /* a 2-d array indexed by i,j with i=0,...,NumMissModes-1 and  j= 0,...,NY-1.  A 1 in any
					entry indicates that being in the i-th missing data mode makes it impossible to observe the X-state 
					indexed by i.  I don't know if I am really going to use this much yet.  */
	
	int Na_M;  /* number of A-states when considering that there might be missing data.  With trios it is 64 and will just get hard-wired as such. */
	int **Ystates_M;  /* when the possible values of each individuals Y is increased by 1 to allow missing data (with SNPs missing data == 3) this
						records the possible states.  These are going to just be hard-wired in, not read in.  0 0 0 = 0,   0 0 1 = 2,   0 0 3 = 4, etc.*/
	int *ymaxes_M;  /* every element is one larger than in ymaxes.  Used to quickly hash values using ReturnVDex */
	int *AgivenY_M;  /* pick an element of this using ReturnVdex to get to the A-with-missing-data given a Y-with-missing data. This will probably just end up being the
						identity array, but I want to have this in here to offer more future flexibility. */
	
	int **AtoA_64;  /* this is a little array indexed by miss mode [mm] say and then by uncensored A state [a].   mm = 0,...,NumMissModes and 
					a = 1,...,Na.  It tells us how a pattern of missing data yields a 64-state A from a corresponding 27 state A. */
	
} RunPars;


/*! now here is the stuff for the forwards backwards algorithm */
struct fb_cell{
	int idx;  /* position of the cell in 0,...,K_t-1 */
	double p; /* probability*/
	double fake_p; /* to store the "value" which is like probability, but includes the fake-probs when appropriate for missing data */
	int *s; /* stores the value of s that this cell refers to */
	struct fb_cell **up; /* array of length d of pointers each corresponding to an X-value pointing to cells in next time step */
	struct fb_cell **down; /* like up, only each element points to the appropriate cell in the previous time step */
	double *dp; /* array.  each element holds the probability of the "down" cell times the probability of getting here from there, normalized to sum to one */
} ;

typedef struct {
	int *Kt; /* the number of states at each time step between 0 and L-1 */
	struct fb_cell **FB;  /*!< This is where the actual cells will get put */
	double S_up;
	double S_down;
	double *FinalStepNormalizedProbs; /*!< probs in FB[L-1] normalized to sum to 1. */
} fb_struct;


typedef struct {
	char *CollName;
	
	double **Xprobs; /*< array that is L by Nx which gives the probs of the X-states in this collection */
	double **Xfakeprobs; /*< array that is L by Nx which gives the probs of the X-states when missing data gets a prob of one.  
	 I am doing this in the hopes that it will allow me to treat missing data reasonably well.  But
	 currently feel I am on tenuous ground in terms of theory, hence the silly name "fakeprobs" */
	double **Aprobs; /*< array that is L by Na which gives the probs of the A-states */
	double ***AprobsGX; /* array that is L by Nx by NaInX that gives the probs of each a state with 
	 X given it is in X.  (i.e. these sum to 1 over the NaInX of each */
	
	double **Aprobs_M; /*< array that is L by Na_M which gives the probs of the A-states when missing data are possible and summed out */
	
	fb_struct *FBS;
	
} Collection;


/*
 This is a little struct to describe the genotypes of a trio in a way that 
 will be useful with regard to the fb_struct.  
 */
typedef struct {
	int anchor;  /* the index of the fb_cell that corresponds to the S vector at the
	 last locus (the one with index L-1) */
	int *path;  /* this has length L and it records the path taken through the down pointers
	 from L-1,...,1.  Each element is between 0 and Nx-1 inclusive.  It corresponds
	 to the values simluated from the dp array in each cell visited in the 
	 backward path.  And then the 0th element should hold the X state of the 0-th
	 locus.  */
	
	int *agx;  /* the index of the A state within the X-state at each time step from 
	 0 up to L-1.  This index corresponds to a probability in the AprobsGX array */
	
	int *a;  /* the absolute index of the A state */
	
} trio_genotype;


typedef struct {
	RunPars *RP;
	
	int NumColls;
	Collection **Colls;  /* array of pointers.  To access them use the macro PPTC_IDX or CPTC_IDX */
	
	int *vdex;  /* we will assume that the vdex is the same for all the collections.  This means that
					in each collection there must be non-zero probability of every genotypic state.  That 
					means that no allele frequency can be zero and no gtyping error rate can be zero!  */
	
} FB_Vars;




/* function prototypes */
void SumSV(int *Sprev, int *v, int *Sout, int d);
void AllocFB_innards( struct fb_cell *f, int d, int Nx);
void FreeFB_innards( struct fb_cell *f);
void AllocTrioGenotypeInnards(trio_genotype *GT, int L);
void FreeTrioGenotypeInnards(trio_genotype *GT, int L);
int ReturnVDex(int d, int *smax, int *s);
int InSDown(int *x, int *s, int d);
fb_struct *ForwardStep(RunPars *R, double **X, double **FX, int *vdex, int VDexSet);
void FreeFB_Vars_innards(RunPars *RP, fb_struct *F);
void PrepForBackwardStep(fb_struct *FBS, RunPars *RP, double **X);
double SimBackward(fb_struct *FBS, double ***A, RunPars *RP, trio_genotype *GT);
double CondProbOfGeno(Collection *C, trio_genotype *GT, RunPars *RP, int Condition);
FB_Vars *Get_pbt_C_fb_Opts(int argc, char *argv[]);
void PrintDetailedCollectionSummary(Collection *C, RunPars *RP);
Collection *AllocToCollection(RunPars *RP);
void ConditionAprobsOnYobs(double **AprobsIn, int **Yobs, RunPars *RP, double **AprobsOut);
void FillYsFromAs(int *A, int **Y, int L, int **Ystates, int NY, int *Mask);
int GiveAFromY(int *y, RunPars *RP);
int GiveAFromY_M(int *y, RunPars *RP);
void ReturnAProbsM_HardWiredForTrios(RunPars *RP, double *AP, double *ret);
int *A_64_ArrayFrom3IndivCharVecs(char *y, char *pa, char *ma, RunPars *RP);
double UnconditionalGenoProb_A64(double **Aprobs_M, int *A, int L);
double LogUnconditionalGenoProb_A64(double **Aprobs_M, int *A, int L);
double *TrioPostProbs(int *A_array, FB_Vars *PP,int SinglePopCollWithPedsOrdered, double *Pi, int Pop, double *LogL, double *UseThisAsOutput, double *PrPar, double *PrNonPar);
int cmp_parent_pair(const void *a, const void *b);
int ReturnXstateOfA_64(int A);
void XStuffFromAStuff(double **Aprobs, RunPars *RP, Collection *Out);
