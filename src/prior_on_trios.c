#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <limits.h>



#include "ECA_MemAlloc.h"
#include "MathStatRand.h"
#include "ranlib.h"
#include "ECA_Opt3.h"
#include "MCTypesEtc.h"
#include "prior_on_trios.h"



trio_prior_clopts *GetTrioPriorOpts(int argc, char *argv[]) 
{
	trio_prior_clopts *ret = (trio_prior_clopts *)malloc(sizeof(trio_prior_clopts));
	
	int infile_f=0;
	DECLARE_ECA_OPT_VARS  
	
	/* This is a good place to set some default values for variables in your program */
	
	/* some information about the program, author, etc */
	SET_PROGRAM_NAME("prior_on_trios");  /* Note:  QUOTED string */
	SET_PROGRAM_SHORT_DESCRIPTION("Compute expected fraction of different trio types"); /* QUOTED string */
	SET_PROGRAM_LONG_DESCRIPTION(This is a simple program to compute the expected fraction 
								 of trios falling into different relationship categories given 
								 an age-structured semelparous population and counts of the number of
								 spawners in the past two generations or so.);  /* UN-QUOTED string! */
	SET_PROGRAM_AUTHOR_STRING("Eric C. Anderson"); /* QUOTED string */
	SET_VERSION("Version 00");  /* QUOTED string */
	SET_VERSION_HISTORY("Version 00 project launched on July 21, 2009\n"); /*QUOTED string */
	
	
	BEGIN_OPT_LOOP  
	
	
	OPEN_SUBSET(Input Options,  /* Name of Subset for output in --help, --help-full, and --help-nroff */
				Input Options, /* Name of subset for lumping them together in GuiLiner */
				Blah blah blah /* Description of subset.  Not really used currently */
	)   /* NOTE.  These are all UNQUOTED strings.  Beware of commas. */
	
	
	
	if(REQUIRED_OPTION(
			Input file,
			infile_f,
			f,
			,
			F,
			path to input file,
			F is the path to the input file. This file is a simple white-space-delimited file.  The first 
					   item in the file must be an integer giving the maximum age of a reproducing individual in the population.
					   The second item must be the year of the earliest year for which you will be supplying counts of reproductive
					   individuals (for example\054 1994 or 2001).  The third item is the latest year for which the file will be
					   providing data.  Data must be given for consecutive years between the earliest and the latest\054 inclusive.
					   Following these three items are a series of blocks\054 one for each year.
					   The first element in each block is just the year.  Then two doubles which are the Ne/N ratios for males and females
					   respectively in this year.  Then a single double which gives the prob that two randomly selected offspring shares 
					   a father given that they
					   shares a mother. with Following this is a row of MA doubles giving the relative reproductive successes of 
					   males different ages classes from 1 up to MA.  Then a row which is the same thing for females.  Then there is a series
					   of doubles giving the numbers of males of
					   age 1 2 3 and so forth up to the maximum age.  Then a series of doubles giving the numbers of females.  Note that these
					   can be ints but they will get read as doubles by the program.  Comments can be included in the file by enclosing them
					   within & symbols.
	   )) {
		if(ARGS_EQ(1)) {
			GET_STR(ret->infile_name);
		}
	}
	
	
	CLOSE_SUBSET;  /* done with the input options */

	END_OPT_LOOP   /* Macro that loses the Main ECA_Opt Loop and does some error checking 
	 behind the scenes */
	
	return(ret);
}





trio_prior_pars *ReadTrioPriorInputFile(trio_prior_clopts *TPC)
{
	int i,j,tempyear,k;
	trio_prior_pars *ret;
	FILE *in;
	
	if(  (in=fopen(TPC->infile_name,"r"))==NULL ) {
		fprintf(stderr, "Error! Failed to open file %s.  Exiting\n",TPC->infile_name);
		exit(1);
	}
	
	ret = (trio_prior_pars *)malloc(sizeof(trio_prior_clopts));
	
	eat_comments(in,'&');
	fscanf(in," %d",&(ret->MA));
	
	eat_comments(in,'&');
	fscanf(in," %d",&(ret->LoYear));
	
	eat_comments(in,'&');
	fscanf(in," %d",&(ret->HiYear));
	
	
	/* compute the number of years */
	ret->NumYears = ret->HiYear - ret->LoYear + 1;
	
	ret->Nf = (double **)calloc(ret->NumYears,sizeof(double *));
	ret->Nm = (double **)calloc(ret->NumYears,sizeof(double *));
	ret->Rate_f = (double **)calloc(ret->NumYears,sizeof(double *));
	ret->Rate_m = (double **)calloc(ret->NumYears,sizeof(double *));
	ret->Nef = (double **)calloc(ret->NumYears,sizeof(double *));
	ret->Nem = (double **)calloc(ret->NumYears,sizeof(double *));
	
	ret->TotNf = (double *)calloc(ret->NumYears,sizeof(double));
	ret->TotNm = (double *)calloc(ret->NumYears,sizeof(double));
	ret->Lambda_f = (double *)calloc(ret->NumYears,sizeof(double));
	ret->Lambda_m = (double *)calloc(ret->NumYears,sizeof(double));
	ret->FSR = (double *)calloc(ret->NumYears,sizeof(double));
	
	for(i=ret->NumYears-1;i>=0;i--)  {
		k = ret->NumYears-(i+1);
		
		/* get the year and check it */
		eat_comments(in,'&');
		fscanf(in," %d",&tempyear);
		
		if(tempyear != ret->LoYear+k) {
			fprintf(stderr,"Error!  Year mismatch.  Reading the data for year block %d was expecting year %d but read %d. Exiting\n",i,ret->LoYear+i,tempyear);
			exit(1);
		}
		
		
		/* get the lambdas */
		eat_comments(in,'&');
		fscanf(in," %lf",&(ret->Lambda_m[i]));
		
		eat_comments(in,'&');
		fscanf(in," %lf",&(ret->Lambda_f[i]));
		
		
		
		/* get the father identity given mother identity conditional probability */
		eat_comments(in,'&');
		fscanf(in," %lf",&(ret->FSR[i]));
		
		
		
		/* get the reproductive rates for males and females */
		ret->Rate_m[i] = (double *)calloc(ret->MA+1,sizeof(double));
		for(j=1;j<=ret->MA;j++)  {
			eat_comments(in,'&');
			fscanf(in," %lf",&(ret->Rate_m[i][j]));
		}
		ret->Rate_f[i] = (double *)calloc(ret->MA+1,sizeof(double));
		for(j=1;j<=ret->MA;j++)  {
			eat_comments(in,'&');
			fscanf(in," %lf",&(ret->Rate_f[i][j]));
		}
		
		
		
		
		/* get the Nm's and compute Nem's while at it */
		ret->Nm[i] = (double *)calloc(ret->MA+1,sizeof(double));
		ret->Nem[i] = (double *)calloc(ret->MA+1,sizeof(double));
		
		for(j=1;j<=ret->MA;j++)  {
			eat_comments(in,'&');
			fscanf(in," %lf",&(ret->Nm[i][j]));
			ret->TotNm[i] += ret->Nm[i][j];
			ret->Nem[i][j] = ret->Nm[i][j] * ret->Lambda_m[i];	
		}
		
		
		/* get the Nfs and compute the Nef's while at it */
		ret->Nf[i] = (double *)calloc(ret->MA+1,sizeof(double));
		ret->Nef[i] = (double *)calloc(ret->MA+1,sizeof(double));

		for(j=1;j<=ret->MA;j++)  {
			eat_comments(in,'&');
			fscanf(in," %lf",&(ret->Nf[i][j]));
			ret->TotNf[i] += ret->Nf[i][j];
			ret->Nef[i][j] = ret->Nf[i][j] * ret->Lambda_f[i];
		}
	}
	return(ret);
}

void PrintTrioPriorPars(trio_prior_pars *T)
{
	int i,j;
	
	printf("TRIO_PRIOR_PARS: MaxAge   : %d\n",T->MA);
	printf("TRIO_PRIOR_PARS: LoYear   : %d\n",T->LoYear);
	printf("TRIO_PRIOR_PARS: HiYear   : %d\n",T->HiYear);
	printf("TRIO_PRIOR_PARS: NumYears : %d\n",T->NumYears);
	
	for(i=0;i<T->NumYears;i++)  {
		printf("TRIO_PRIOR_PARS: ArrayElement : %d :  Year : %d  : Lambda_m= %f  Lambda_f= %f\n",i,T->HiYear-i,T->Lambda_m[i],T->Lambda_f[i]);
		printf("TRIO_PRIOR_PARS: ArrayElement : %d :  Year : %d  : ConditionalFullSibRate= %f\n",i,T->HiYear-i,T->FSR[i]);
		printf("TRIO_PRIOR_PARS: ArrayElement : %d :  Year : %d  : Males   : ",i,T->HiYear-i);
		for(j=1;j<=T->MA;j++)  {
			printf(" %.1f ",T->Nm[i][j]);
		}
		printf("  (Total: %.1f )  : Ne's: ",T->TotNm[i]);
		for(j=1;j<=T->MA;j++)  {
			printf(" %.1f ",T->Nem[i][j]);
		}
		printf("  : Rates: ");
		for(j=1;j<=T->MA;j++)  {
			printf(" %.1f ",T->Rate_m[i][j]);
		}
		printf("\n");
		printf("TRIO_PRIOR_PARS: ArrayElement : %d :  Year : %d  : Females : ",i,T->HiYear-i);
		for(j=1;j<=T->MA;j++)  {
			printf(" %.1f ",T->Nf[i][j]);
		}
		printf("  (Total: %.1f )  : Ne's: ",T->TotNf[i]);
		for(j=1;j<=T->MA;j++)  {
			printf(" %.1f ",T->Nef[i][j]);
		}
		printf("  : Rates: ");
		for(j=1;j<=T->MA;j++)  {
			printf(" %.1f ",T->Rate_f[i][j]);
		}
		printf("\n");
	}
}



/* given potential parents in a given year having effective counts of different ages in N
 and age-specific repro rates in A, return the probability 
 that two randomly chosen individuals have the same parent.
 
 MA is the max age of individuals.
 
*/
double HalfSibProb(double *Ne, double *A, int MA)
{
	int i;
	double ret=0.0;
	double P[100];  /* for the probability that an offspring comes from a parent in a certain age class */
	double normo=0.0;
	
	
	/* compute the probs that offspring are from parents in the different age classes. */
	for(i=0;i<=MA;i++)  {
		P[i]=Ne[i]*A[i];
		normo+=P[i];
	}
	for(i=0;i<=MA;i++)  {
		P[i] /= normo;
	}
	
	
	/* now, sum over the prob of identity within each age-group, weighted by the 
	   probability that both offspring have parents in the same age-group */
	for(i=0;i<=MA;i++)  {
		if(Ne[i]>0.001) {
			ret += P[i]*P[i] * 1.0/Ne[i];
		}
	}
	
	return(ret);
	
}









/*
	Given certain times in the ancestry:
 tcfs : time the candidate female parent (mother) spawned
 tyc   : time the youth was conceived (i.e. the time the true mother and father spawned)
 tcms : time the candidate male parent spawned
 tcfc : time the candidate female was conceived
 tfc  : time the youth's true female parent was conceived
 tmc  : time the youths male parent was conceived
 tcmc : time the candidate male's parent was conceived
 
 Note that all these times are numbers of years counting backwards from the latest year for which 
 population size data are available.
 
 And given the population sizes and parameters in T
 
 This function puts probababilities of the 22 trio states into the P vector.
 
 This is not at all pretty.  I just kind of brute-force-hard-wire this thing.
 
*/
void ProbOfR(int tyc, int tcms, int tcfs, int tmc, int tfc, int tcmc,  int tcfc,  trio_prior_pars *T, double *P)
{
	int i;
	/* these are for holding probabilities like female parent and candidate parent are Self (f_cfSe) (or Sib (Si) or Unrelated (U)) */
	/* See file:///Users/eriq/Documents/work/prj/PFR/prior_on_trios/prior_on_trios.tex for an enumeration of these things. */
	double f_cfSe, f_cfSi, f_cfU,
		   f_cmSe, f_cmSi, f_cmU,
		   m_cfSe, m_cfSi, m_cfU,
		   m_cmSe, m_cmSi, m_cmU,
		   cm_yB, cm_yH, 
		   cf_yB, cf_yH,
		   cm_y_cfB;
	
	double *TotNf = T->TotNf,
			*TotNm = T->TotNm,
			**Nef = T->Nef,
			**Nem = T->Nem,
			**Rate_m = T->Rate_m,
			**Rate_f = T->Rate_f,
			*FSR = T->FSR;	
	int MA = T->MA;
	double temp;
	double NoHorB;
	
	
	/* now we just write down what these are.  Note that we first check to see if each candidate is Self, if not
	 then we consider the other options. */
	f_cfSe = (tyc == tcfs) ?  1.0/TotNf[tyc] : 0.0;
	
	printf("f_cfSe= %.5e\n",f_cfSe);
	
	
	/* given f_cfSe is not the case, here is the prob that f_cfSi */
	f_cfSi = (tcfc == tfc) ? HalfSibProb(Nef[tfc], Rate_f[tfc],MA) * FSR[tfc]  :  0.0;
	printf("f_cfSi= %.5e\n",f_cfSi);
	
	f_cfU = (1.0 - f_cfSe) * (1.0 - f_cfSi);  /* note that Unrelated means neither self or full sib.  We disregard the occurrence of 
									candidate's being half sibs of the true parents because that doesn't seem to make much difference.  It
									would be easy enough to change that later should we want to.  */
	printf("f_cfU= %.5e\n",f_cfU);
	
	
	
	f_cmSe = 0.0;  /* candidate male can't be a female in reality */
	f_cmSi = (tcmc == tfc) ?  HalfSibProb(Nef[tfc], Rate_f[tfc],MA) * FSR[tfc]  :  0.0;
	f_cmU = (1.0 - f_cmSi);
	
	
	m_cmSe = (tyc == tcms) ?  1.0/TotNm[tyc] : 0.0;
	m_cmSi = (tcmc == tmc) ?  HalfSibProb(Nef[tmc], Rate_f[tmc],MA) * FSR[tmc]  :  0.0;
	m_cmU = (1.0 - m_cmSe)*(1.0 - m_cmSi);
	
	
	m_cfSe = 0.0;
	m_cfSi = (tcfc == tmc) ?  HalfSibProb(Nef[tmc], Rate_f[tmc],MA) * FSR[tmc]  :  0.0;
	m_cfU = (1.0 - m_cfSe)*(1.0 - m_cfSi); 
	
	
	
	
		
	
	/* note that under our model, an individual is a half sib of the youth if they were both conceived in the same year,
	 they share a mother and they don't share a father, plus the prob that they were both conceived in the same year times,
	 they don't share a mother, but the do share a father.  this of course all has to be weighted by the probability that both the 
	 cf and cm are not full sibs of the y. */
	temp = 0.0;
	if(tcmc == tyc || tcfc == tyc) {
		temp = HalfSibProb(Nef[tyc],Rate_f[tyc],MA);
	}
	
	cm_y_cfB = (tcmc == tyc) && (tcfc == tyc)  ?  temp*FSR[tyc] * temp*FSR[tyc]  :  0.0;
	
	/* These are the marginal probs of these evbents */
	cm_yB = (tcmc == tyc) ?  temp * FSR[tyc]  : 0.0 ;
	cf_yB = (tcfc == tyc) ?  temp * FSR[tyc]  : 0.0 ;
	
	
	/* these will represent half sibness and no relationship through the other candidate */
	cm_yH = (tcmc == tyc) ?  (temp * (1.0-FSR[tyc]) +  (1.0-temp)*HalfSibProb(Nem[tyc],Rate_m[tyc],MA))  : 0.0;
	cf_yH = (tcfc == tyc) ?  (temp * (1.0-FSR[tyc]) +  (1.0-temp)*HalfSibProb(Nem[tyc],Rate_m[tyc],MA))  : 0.0;

	
	
	/* now we are ready to combine those into the proabablities of different trios  */
	NoHorB = 1.0 - cm_y_cfB 
		- (cm_yB * (1.0 - cf_yB - cf_yH)) 
		- (cm_yH * (1.0 - cf_yB - cf_yH)) 
		- (cf_yB * (1.0 - cm_yB - cm_yH))
		- (cf_yH * (1.0 - cm_yB - cm_yH));
	
	
	
	
	
	P[12] = cf_yB * (1.0 - cm_yB - cm_yH) * m_cmSe;
	P[13] = cm_yB * (1.0 - cf_yB - cf_yH) * f_cfSe;
	P[14] = cm_yH * (1.0 - cf_yB - cf_yH) * f_cfSe;
	P[15] = cf_yH * (1.0 - cm_yB - cm_yH) * m_cmSe;
	P[16] = cm_yB * (1.0 - cf_yB - cf_yH) * f_cfSi;
	P[17] = cm_yB * (1.0 - cf_yB - cf_yH) * m_cfSi;
	P[18] = cf_yB * (1.0 - cm_yB - cm_yH) * f_cmSi;
	P[19] = cf_yB * (1.0 - cm_yB - cm_yH) * m_cmSi;
	P[20] = cm_yB * (1.0 - cf_yB - cf_yH) * m_cfU * (1.0-f_cfSe)*(1.0-f_cfSi);
	P[21] = cf_yB * (1.0 - cm_yB - cm_yH) * m_cmU * f_cmU;
	P[22] = cm_y_cfB;

	
	NoHorB = 1.0;
	for(i=12;i<=22;i++)  {
		NoHorB -= P[i];
	}
	
	P[0] = NoHorB * m_cmSe * f_cfSe;
	P[1] = NoHorB * f_cfSe * (1.0-m_cmSe) * m_cmSi;  // we assume the possibility that cm is a Si of both m and f is too small to 
													// worry about (though we could stick it in there if we wanted, and just keep in 
													// mind that the C_Se_Si class includes that case
	P[2] = NoHorB * m_cmSe * (1.0-f_cfSe) * f_cfSi;
	P[3] = NoHorB * f_cfSe * (1.0-m_cmSe)*(1.0-m_cmSi-m_cfSi);
	P[4] = NoHorB * m_cmSe * (1.0-f_cfSe)*(1.0-f_cfSi-m_cfSi);
	P[5] = NoHorB * (1.0-f_cfSe) * f_cfSi * (1.0-m_cmSe) * m_cmSi;
	P[6] = NoHorB *  (1.0-m_cmSe) * f_cmSi * (1.0-f_cfSe) * m_cfSi;
	P[7] = NoHorB * (1.0-f_cfSe) *f_cfSi * (1.0-m_cmSe)*(1.0-m_cmSi);
	P[8] = NoHorB * (1.0-m_cmSe) * f_cmSi * (1.0-f_cfSe)*(1.0-f_cfSi);
	P[9] = NoHorB * (1.0-m_cmSe) * m_cmSi * (1.0-f_cfSe)*(1.0-f_cfSi);
	P[10] = NoHorB * (1.0-f_cfSe) * m_cfSi * (1.0-m_cmSe)*(1.0-f_cmSi);
	P[11] = NoHorB * (1.0-f_cfSe)*(1.0-f_cfSi-m_cfSi) * (1.0-m_cmSe)*(1.0-m_cmSi-f_cmSi);
		
	//printf("NoHorB= %.4e\n",NoHorB);
	
}




/*
 Pretending that we know tyc -- the year the youth was conceived -- and knowing tcfs and tcms (which we
 should know, since we know the year that fish spawn) this sums over all the years when the parents and candidate parents 
 could have been conceived and given us what we want in P.  This is to be used when we don't know the ages of cf and cm
 themselves, but we have estimates of the fraction of fish of different ages in the year they were spawned.  
 
 As above, the years are all given counting backward with 0 being the latest year for which 
 population data are available. 
 
 Note that if we had trait data on the candidate male and female then we could adjust the 
 Ptcfc and Ptcmc values as need be if the trait gave us any information about the age of the 
 candidate male or female.  
 
*/
void ProbOfR_marg(int tyc, int tcms, int tcfs,   trio_prior_pars *T, double *P)
{
	int i;
	int tcfc, tfc, tmc, tcmc;
	int cfc, fc, mc, cmc;
	double Ptcfc, Ptfc, Ptmc, Ptcmc;
	int MA = T->MA;
	double tempA[30];
	
	
	/* initialize P to accumulate a sum */
	for(i=0;i<30;i++) {P[i] = 0.0;}  
	
	
	/* first cycle over all possible years when the true father might have been conceived */
	for(mc=1;mc<=MA;mc++)  {
		tmc = tyc + mc;
		Ptmc = T->Nm[tyc][mc]/T->TotNm[tyc];  /* this is the probability that he was conceived at tmc given what we know about the
												 age distribution of males at time tyc */
		

		/* then cycle over all possible years when the true mother might have been conceived */
		for(fc=1;fc<=MA;fc++)  {
			tfc = tyc+fc;   /* this is the actual year when she was conceived */
			Ptfc = T->Nf[tyc][fc]/T->TotNf[tyc];  /* this is the probability that she was conceived then given what we know about the
													  age distribution of females at time tyc */
		
		
				
			/* and now we cycle over all the possible years when the candidate father might have been conceived */
			for(cmc=1;cmc<=MA;cmc++)  {
				tcmc = tcms + cmc;
				Ptcmc = T->Nm[tcms][cmc]/T->TotNm[tcms];
				
				
				/* and finally we cycle over all the possible years when the candidate female might have been conceived */
				for(cfc=1;cfc<=MA;cfc++)  {
					tcfc = tcfs + cfc;
					Ptcfc = T->Nf[tcfs][cfc]/T->TotNf[tcfs];
					
					/* now we compute the probs of all the different trio relationships and then add them into P
					   while weighting by the product of the P terms.  But only do this if the probability of all the 
					   years of conception is greater than 0. */
					if(Ptmc * Ptfc * Ptcmc * Ptcfc > 0.0000000000001) {
						ProbOfR(tyc, tcms, tcfs, tmc, tfc, tcmc, tcfc, T, tempA );
						printf("INTERMEDIATE: tmc= %d  tfc= %d  tcmc= %d  tcfc= %d   Prob of this:  %.4e\n",tmc, tfc, tcmc, tcfc, Ptmc * Ptfc * Ptcmc * Ptcfc);
						for(i=0;i<=22;i++)  {
							printf("INTERMEDIATE_PROBS:  %d  %.4e\n",i,tempA[i]);
							P[i] += tempA[i] * Ptmc * Ptfc * Ptcmc * Ptcfc; 
						}
					}
				}
			}
		}
	}
}


/* this is a little function that takes the probabilities in the Probs array and combines them as need be
 to get results that correspond to the categories in the pfr program. See
 file://localhost/Users/eriq/Documents/work/prj/PFR/prior_on_trios/prior_on_trios.tex
 for a description of these.
 
 i is the PFR relationship category and P is the array of Probs that have been computed for the
 different categories in this program.
 
*/
double TranslateProbsToPFR_Categories(double *P, int i) 
{
	switch(i) {
		case(0): 
			return(P[0]);
			break;
		case(1): 
			return(P[2]);
			break;
		case(2): 
			return(P[1]);
			break;
		case(3): 
			return(P[4]);
			break;
		case(4): 
			return(P[3]);
			break;
		case(5): 
			return(P[5]+P[6]);
			break;
		case(6): 
			return(P[9]+P[10]);
			break;
		case(7): 
			return(P[7]+P[8]);
			break;
		case(8): 
			return(P[11]);
			break;
		case(9): 
			return(P[12]);
			break;
		case(10): 
			return(P[13]);
			break;
		case(11): 
			return(P[14]);
			break;
		case(12): 
			return(P[15]);
			break;
		case(13): 
			return(P[16]+P[17]);
			break;
		case(14): 
			return(P[18]+P[19]);
			break;
		case(15): 
			return(P[20]);
			break;
		case(16): 
			return(P[21]);
			break;
		case(17): 
			return(P[22]);
			break;
		default:
			fprintf(stderr, "Unrecognized PFR category i= %d in TranslateProbsToPFR_Categories. Exiting....\n",i);
			exit(1);
	}			
}


/* given all the trio prior pars, this function cycles over all possible spawning years of candidate males and females,  within each of those
 it then cycles over the possible conception years of youths that might possibly be compared to these males and females and it computes the 
 array of probabilities for each such tcms/tcfs and tyc compbo.  It also prints that in a format that can be read in later where each value 
 is tagged by the respective years.
 
 These are the limits we set on the years we will cycle over which are dictated by the subscripts we will encounter in ProbOfR_marg,
 and subsequently in ProbOfR().

 In particular, it appears that tyc, tcms, and tcfs are used as subscripts in  ProbOfR_marg, so we want to be sure that none of these
 ever exceeds NumYears-1
 
 But, then, additionally, in ProbOfR, these variables are used as subscripts:  tfc and tmc.  The max values these can take depend on other 
 variables when values are getting assigned to them in ProbOfR:
	Max of tfc = tyc + MA
	Max of tmc = tyc + MA
 
 So, it would appear that we want to restrict tyc to be <= NumYears-1 - MA.
 
 And, on the low end for tyc, we clearly want to choose 0.  OK.
 
 So we will cycle over tcms and tcfs from 0 to NumYears-1, and within each of those we can cycle over tyc from 0 to NumYears-1 - MA.
 
 
 */
void GenerateAllOutput(trio_prior_pars *T)
{
	int tcps;  /* time candidate parent spawned.  cycle over this and set tcms and tcfs to it. */
	int tyc;  /* time youth conceived */
	int tcms, tcfs;
	int NumYears = T->NumYears;
	int MA = T->MA;
	double *Probs = (double *)calloc(30, sizeof(double));
	int i;
	
	for(tcps=0;tcps<NumYears;tcps++)  {
		tcms=tcfs=tcps;
		for(tyc=0;tyc<NumYears-MA;tyc++) {
			ProbOfR_marg(tyc, tcms, tcfs, T, Probs);
			printf("FINAL_OUTPUT:  tyc= %d   tcms= %d   tcfs= %d  ",tyc, tcms, tcfs);
			printf("  %d   %d  ",T->HiYear-tyc, T->HiYear-tcms); /* this just prints the actual years of the youth conception and of the candidate parent spawning */
			for(i=0;i<23;i++) {
				printf("%.16e ",Probs[i]);
			}
			printf("\n");
			
			/* here we print them according to their PFR categories */
			printf("PFR_RELAT_OUTPUT:   %d   %d  ",T->HiYear-tyc, T->HiYear-tcms);
			for(i=0;i<18;i++)  {
				printf(" %.10e",TranslateProbsToPFR_Categories(Probs,i));
			}
			printf("\n");
			
		}
	}

}

int main(int argc, char *argv[])
{
	int i;
	trio_prior_clopts *TPC;
	trio_prior_pars *TPP;
	int tyc=0, tcps=0;
	int tcms, tcfs;
	
	tcms=tcfs=tcps; 
	
	double *Probs = (double *)calloc(30, sizeof(double));
	
	TPC = GetTrioPriorOpts(argc,argv);
	TPP = ReadTrioPriorInputFile(TPC);
	PrintTrioPriorPars(TPP);
	
	
	GenerateAllOutput(TPP);
	
	return(0);
	
	
	ProbOfR_marg(tyc, tcms, tcfs, TPP, Probs);
	
	printf("FINAL_PROBS:  tyc= %d   tcms= %d   tcfs= %d\n",tyc, tcms, tcfs);
	for(i=0;i<23;i++) {
		printf("FINAL_PROBS: %d  %.5e\n",i,Probs[i]);
	}
	
	return(0);
	
}
