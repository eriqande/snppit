typedef struct {
	char infile_name[4000];
	
} trio_prior_clopts;

typedef struct {
	int MA; /* max age of individuals */
	
	int LoYear;  /* the earliest year for which the pop size data will be supplied */
	int HiYear;  /* the latest year from which the pop size data will be supplied */
	int NumYears; /* HiYear - LoYear + 1 */
	
	
	/* these could be int's, but much of the time we will be doing floating point ops on them anyway, so save them as doubles to 
	 avoid casting all the time */
	double **Nf; /* Number of females at each year.  Subscripted by years (from 0 to NumYears-1) and by ages (from ages 0 to MA) */
	double **Nm; /* Number of males each year.  Subscripted like Nf */
	
	double **Nef; /* effective number of females.  Basically just Nf times the appropriate Lambda */
	double **Nem;
	
	double *TotNf;  /* sum of females spawners in a years over all age classes */
	double *TotNm;
	
	double *Lambda_f;  /* Ne to N ratio for females in different years */
	double *Lambda_m;  /* Ne to N ratio for females in different years */
	
	double **Rate_f; /* age specific relative reproductive rates as a function of year. Subscripted by year and then ages */
	double **Rate_m; 
	
	double *FSR;  /* the "Full-Sib Rate" i.e. the conditional probability that a pair shares a father given that they share a mother. Indexed by year */
	
} trio_prior_pars;





/* function prototypes */
trio_prior_clopts *GetTrioPriorOpts(int argc, char *argv[]);
trio_prior_pars *ReadTrioPriorInputFile(trio_prior_clopts *TPC);
void PrintTrioPriorPars(trio_prior_pars *T);
double HalfSibProb(double *Ne, double *A, int MA);
void ProbOfR(int tyc, int tcms, int tcfs, int tmc, int tfc, int tcmc,  int tcfc,  trio_prior_pars *T, double *P);
void ProbOfR_marg(int tyc, int tcms, int tcfs,   trio_prior_pars *T, double *P);
void GenerateAllOutput(trio_prior_pars *T);
double TranslateProbsToPFR_Categories(double *P, int i);