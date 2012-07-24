#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <stdarg.h>
#include "MathStatRand.h"
#include "ECA_MemAlloc.h"
#include "ranlib.h"
#include "MCTypesEtc.h"
//#include "ECA_Opt.h"
#include "ECA_print.h"

/* Some verbosity constants */
#define RECCYC_VERB 0x10000000
#define RECSUMOVER_VERB 0x20000000
#define PROBPEDSTATE_VERB 0x40000000
#define REPORTPED_VERB 0x80000000

/* this is what we use to get the standard "pre-xml and pre-eca-opt" output */
#define STANDARD_VERB 0x01000000

/* this gives us the ECA_OPT output */
#define ECA_OPT_VERB 0x02000000

/* this is a flag saying that we want the ECA_OPT output, but we want it bear---without the ECA_OPT_OUTPUT tags */
#define ECA_OPT_BARE_VERB 0x04000000


#define XML_VERB 0x08000000



typedef struct ind_ptr {
	int ID;
	struct ind_ptr *ma;
	struct ind_ptr *pa;
	struct ind_ptr *pair;  /* pointer to ind with which it has pairwise relatedness not shown in the pedigree */
	double K[3];  /* for the k-coefficients of relatedness */
	int InPair;  /* gets 1 if this individual is part of a pairwise relatedness pair */
	int g;  /* genotype---0,1, or 2 */
	int fixed;  /* 1 if it is in the "fixed" subset */
	int alive;  /* 1 if it occurs in the pedigree, 0 otherwise */
} ind;


/* for parameter values that the user may supply to the program */
typedef struct {
	double p; /* minor allele freq's at a single locus.  This gets used as a sort of holding variable */
	double *Ps; /* an array of the allele freqs to look at */
	int NumPs;  /* the number of allele freqs to look at */
	double *mu;  /* single-gene-copy genotyping error rate */
	int NumMu;
	double *YouthPs;
	int NumYouthPs;
	char *File;
	int L; /* number of lines in the pedfile*/
	int Nfixed; /* number of fixed genotype indivs */
	int *fixed;  /* list of  fixed genotypes */
	int MaxN; /* largest value that an individuals integer ID could possibly take */
	ind *ped;
	char RunIdName[1000];
	int PrintXMLPream;
	char CollectionGroupName[5000];
	int PrintXMLOutput;
	int PrintXMLClose;
	int Smax[6]; 
	
	int MissMeth;  /* 0=don't keep track of missing loci;  1=count all loci with at least one missing member of trio;  3=keep track of missing loci in kids, moms and dads. */
	int NumMiss;
	double **MissRates;  /* rate at which data are missing amongst the youths, putative fathers, and putative mothers, respectively
	 at a single locus.  Subscripted by locus and then by 0, 1, or 2 for amongst kids, fathers, or mothers, respectively. */  
} sum_ped_pars;


/* function prototypes */
double PairWiseProbFunc(int g1, int g2, double p, double K[3]);
ind *AllocPedStruct(int MaxN);
void ProbOfObservedStates(double mu, double TrueProbs[3][3][3], double ObsProb[3][3][3]);
void ProbOfObservedStatesWithMiss(double *mu, double TrueProbs[3][3][3], double ObsProb[4][4][4], int MissingIsUnity);
double ProbOfPedState(ind *ped, int MaxN, double p);
void RecursSumOver(ind *ped, int i, int MaxN, double p);
void RecursCycleOverFixed(ind *ped, int Nfixed, int *fixed, int i, int *states, double p, int MaxN);
void ReportPedState(ind *ped, int MaxN);
sum_ped_pars *getSumPedOpts(int argc, char **argv);
double CorrectTrioProb(int m,int f, int y);
int IndicIncomp(int m, int f, int y);
double ObsProbs(int a, int b, double mu);
double MissObsProb(int a, int b, double d, int MissingIsUnity);
double OneRandomParent(int m, int y, double p);
double PofG(int G, double p);
void Print_ECA_Opt_LocusProbs(const char *Ind, double ObsProbs[3][3][3], double MissObsProbs[4][4][4], double MissFakeProbs[4][4][4], sum_ped_pars *pars);
int snpSumPed(int argc, char **argv, FILE *outfile);

