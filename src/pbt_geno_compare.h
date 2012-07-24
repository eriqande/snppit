/*
	Define these macros to get more verbose output.
 
	VERBOSE_PRINT_COMPAT_TRIOS						print the compatible trios and the their genotypes before computing posteriors for all of them
 
 
*/
#ifdef VERBOSE_ALL
#define VERBOSE_PRINT_COMPAT_TRIOS
#endif

//#define VERBOSE_PRINT_COMPAT_TRIOS


#define MAX_TEMP_PARS 1000000
#define MAX_TEMP_PAIRS 1000000



/* some function prototypes */
int PairCompatByMendel(char *Par, char *Kid, int MaxAllowable, int L, int *NumInc);
int BirthTimingWorks(
					 int NBY,					/* number of possible birth years of the kid */
					 int *BY,					/* array of possible birth years of the kid */
					 pfr_year_range *RY			/* struct holding the possible reproductive years of the parent */
					 );


pfr_parent **ReturnMatchingParents(
								   pfr_offspring *Y,			/* all the info about the youth */
								   int ReGenoNumKid,				/* which genotype (from various regenotypings) to use for offspring.  Typically will be zero */
								   int ReGenoNumPar,				/* which genotype to use for the parent.  Typically zero */
								   pbt_high_level_vars *HLV,			/* all the other info we might need */
								   SexEnum ParSex,				/* the sex of the parents we are looking for (MALE/FEMALE/SEX_UNKNOWN) */
								   int *NumRet,					/* output parameter---number of parents in returned pfr_parent * array. */
								   int **NumMI						/* another output parameter, for the number of mendelian incompatibilities at each parent in the array. 
																	we have to pass this in as an int ** because we will allocate memory to it in this function. */
								   );

int MateTimingWorksOut(pfr_year_range *A, pfr_year_range *B, int *intsct);
int TrioIsIncompat(char m, char f, char y);
ParentPairMatch *ReturnMatchingTrios(
									 pfr_offspring *Y,						/* pointer to the offspring */
									 pbt_high_level_vars *HLV,				/* pointer to everything */
									 int MaxAllowable,							/* max allowable number of mendelian incompatibilities at the trio */
									 int ReGenoMa,							/* re-genotyping numbers for ma, pa, kid---typically 0 */
									 int ReGenoPa,
									 int ReGenoKid,
									 int *NumTrios							/* output parameters telling us how many trios were compatible */
									 );


