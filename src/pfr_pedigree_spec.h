#define NUM_SPEC_PEDS 18

#ifdef UN_EXTERN
#define GLOB
#else 
#define GLOB extern
#endif


FB_Vars *CreateCompleteFB_Vars(
							   int ParentPopIdx,		/* if >=0, compute forback tables assuming parents come from this single population, if -1 do all the pops */
							   int KidPopIdx,			/* if >=0, compute forback tables assuming kids come from this single population, if -1 do all the pops, if -2 then kids' pop = parent's pop always */
							   int PedSpecIdx,			/* if>=0, just compute forback tables for this one relationship. If -1 compute for all relationships.  */
							   int *smax,				/* the smax value that is desired */
							   int smax_d,				/* number of components of this smax vector */
							   pfr_geno_data *PFR,		/* holds all the data about number of loci and allele freqs and such */
							   char *SSP_InputFileName,	/* file into which to spew the snpSumPed input command-file */
							   char *SSP_OutFileName		/* file into which to spew the snpSumPed output */
							   );




void Print_snpSumPed_CommandFile(char *FileName,
								 int PedIdx,
								 double **p,
								 double **youth_p,  /* if different than parental p.  Otherwise pass in as NULL */
								 double *mu,
								 int *MissMeth,
								 int DoPream,
								 char *MaPopPopName,
								 char *KidPopName,
								 char *RunName,
								 int NumLoc,
								 const char *SmaxStr
								 );




const char *SpecPedIdx2Str(int P);


/* a global variable to use as the path to snpSumPed */
/* char gSnpSumPedPath[5000]; */
/* int gHasSnpSumPedPath; */












