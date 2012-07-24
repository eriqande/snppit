

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <limits.h>

#define UN_EXTERN


#include "ECA_MemAlloc.h"
#include "MathStatRand.h"
#include "ranlib.h"
#include "ECA_Opt.h"
#include "MCTypesEtc.h"
#include "pfr_forback.h"
#include "pfr_utils.h"







pfr_forback_data *GetPFR_ForBackOptions(int argc, char *argv[])
{
	int i,j,*tempa;
	int Popf=0,
		Af=0,
		Mf=0,
		mcoptsf=0;
	pfr_forback_data *ret = (pfr_forback_data *)malloc(sizeof(pfr_forback_data));
	DECLARE_ECA_OPT_VARS;
	
	/* set some default values and do some allocation, etc */
	ret->Pops = (pfr_pop_struct **)calloc(MAX_NUM_POPS,sizeof(pfr_pop_struct *));
	ret->MC = (mc_struct **)calloc(MAX_NUM_MCOPTS,sizeof(mc_struct *));
	ret->L = 0;
	ret->NumPops = 0;
	ret->NumMC = 0;
	
	


	/* some defaults */
	SET_OPT_WIDTH(25);
	SET_ARG_WIDTH(34);
	SET_VERSION("\npfr_forback -- calculate probs of Mendelian incompats and simulate SNP genotypes conditional on same --\n\nVERSION: 1.0 Beta\nAUTHOR: Eric C. Anderson (eric.anderson@noaa.gov)\nDATE: 24 July 2006\nCOPYRIGHT: None -- US Govt Employee Work\n\n");
	SET_DESCRIBE("\npfr_forback -- calculate probs of Mendelian incompats and simulate SNP genotypes conditional on same --\n");

	BEGIN_OPT_LOOP 	 
		if(REQUIRED_OPTION(
			Af,
			a,
			analysis-type,
			S,
			name of analysis,
			S should be one of: mapatri patri matri tri.  mapatri means that you screened father pairs and mother pairs
			and then screened the full trios composed of compatible pairs.  patri means you screened father pairs and then 
			trios that they could be members in.  matri means you screened mother pairs and then 
			trios that they could be members in.  tri means you just screened trios and no pairs of any sort. ) ) 
		{

			if(ARGS_EQ(1)) { char tempStr[100];
				GET_STR(tempStr);
				if(strcmp(tempStr,"mapatri")==0) {
					ret->aType = MAPATRI;
				}
				else if(strcmp(tempStr,"patri")==0) {
					ret->aType = PATRI;
				}
				else if(strcmp(tempStr,"matri")==0) {
					ret->aType = MATRI;
				}
				else if(strcmp(tempStr,"tri")==0) {
					ret->aType = TRI;
				}
				else {
					fprintf(stderr,"String \"%s\" is an invalid option to -a/--analysis-type.  Expecting mapatri, patri, matri, or tri. Exiting...\n",tempStr);
					exit(1);
				}
				printf("READ_DATA: Arg \"%s\" to -a/--analysis-type sets variable aType to %d \n",tempStr,ret->aType);
				/* Here we do a minor bit of setting a few more aType related variables */
				ret->Xdim = XdimOfaType(ret->aType);
				ret->XNum = XNumOfaType(ret->aType);
				ret->Xkeys = SetXKeys(ret->aType);
				ret->NumA = NumAofaType(ret->aType);
			}
		}
		if(MULT_USE_REQ_OPTION(
			Popf,
			p,
			pop,
			S Loc1 R1 ... RN Loc2 R1 ...,
			supply information about genotype freqs in a pop,
			Use this option to input genotype frequencies for all the loci in the study.  The input format is somewhat complex.
				S is the name given to this population.  It might be something like PopA_NullHypothesis or PopA_Parental.  The point
				is that what is mean by a population can include a hypothesis about relatedness.  After the name of the population the
				probabilities of the different combinations of possible genotypes among the individuals of interest are put on the command line.
				For analysis of trios there are 27 such frequencies for each locus---i.e. RN is R27 in the argument list.  If one includes missing data as a genotypic state then
				there are 64 such states.  If you are doing an individual matching analysis there are only 9 states.  The locus information must
				go in order.  For the X-th locus the frequencies must be preceded by LocX. For example: Loc1 or Loc122.  This tag is
				case sensitive.  The total number of loci will be calculated by the program from the length of the command line.  The number
				of loci must be identical for each population.  , 
				MAX_NUM_POPS))
		{
			if(ALREADY_HAS(Af,-a/--analysis-type)) {
				if(ARGS_GEQ(ret->NumA+2)) { int tempL; int tempArgs; char tempPop[MAX_POP_NAME]; char tempLocName[50]; char tempLocFlag[500];
					tempArgs = COUNT_ARGS;
					GET_STR(tempPop);
					if( (tempArgs - 1) % (ret->NumA + 1) != 0) {
						fprintf(stderr,"Error!  You have %d args to option -p/--pop with (pop-name \"%s\").  The number of args minus one should be perfecty divisible by %d. Exiting...\n",tempArgs,tempPop,ret->NumA+1); 
						exit(1);
					}
					tempL = (tempArgs - 1) / (ret->NumA+1);
					if(ret->L==0) {
						ret->L = tempL;
						printf("READ_DATA: Setting Number of Loci to %d\n",ret->L);
					}
					else if(ret->L != tempL) {
						fprintf(stderr,"Error! Previous calls to -p/--pop set values for %d loci.  The call with pop-name %s is setting values for %d loci.  Exiting...\n",ret->L,tempPop,tempL);
						exit(1);
					}
					/* down here and we are good to go */
					/* allocate some memory */
					ret->Pops[ret->NumPops] = (pfr_pop_struct *)malloc(sizeof(pfr_pop_struct));
					ret->Pops[ret->NumPops]->Loci = (pfrfb_locus **)calloc(ret->L,sizeof(pfrfb_locus *));
					
					/* check to make sure the popname is unique compared to previous names, then transfer the pop name */
					if(IdxFromPopName(tempPop,ret)>-1) {
						fprintf(stderr,"Error! While reading population with index %d: PopName %s has been used previously.  PopNames must be unique.  Exiting...\n",ret->NumPops,tempPop);
						exit(1);
					}	
					sprintf(ret->Pops[ret->NumPops]->Name,tempPop);
					
					/* then cycle over loci and get all the info */
					for(i=0;i<ret->L;i++)  {
						/* pick up the correct locus flag */
						sprintf(tempLocName,"Loc%d",i+1);
						GET_STR(tempLocFlag);
						if(strcmp(tempLocFlag,tempLocName) != 0) {
							fprintf(stderr,"Error! While reading data from -p/--pop option for pop-name.  Expecting to chew flag %s, but got %s instead.  Exiting...\n",tempLocName,tempLocFlag);
							exit(1);  
						}
						
						/* allocate necessary memory to the locus */
						ret->Pops[ret->NumPops]->Loci[i] = (pfrfb_locus *)malloc(sizeof(pfrfb_locus));
						ret->Pops[ret->NumPops]->Loci[i]->PA = (double *)calloc(ret->NumA,sizeof(double));
						
						/* then get the trio genotype frequencies */
						for(j=0;j<ret->NumA;j++) {
							ret->Pops[ret->NumPops]->Loci[i]->PA[j] = GET_DUB;
						}
					}
					ret->NumPops++;
				}
			} /* closes if ALREADY_HAS */
		}
		
		
		if(REQUIRED_OPTION(
			Mf,
			m,
			max-misses,
			J1 ...,
			max value of different types of states to keep track of,
			This option allows you to specify the maximum value of particular states that you are keeping track of. 
			It must come after the -a/--analysis-type option since that option tells us what the different states are
			and how many there will be.  The different arguments refer to the different states.  For example if we use
			-a mapatri then the program will expect three integer arguments to -m/--max-misses.  Namely J1 is the max
			number of Mendelian incompatibilities between the father and the kid alone.  J2 is the max number of Mendelian 
			incompatibilities between the mother and the kid alone. J3 is the max number of incompatibilities in the 
			trio as a whole.  For matri J1 is between mother and kid and J2 is between trio as a whole.  For patri J1 is 
			between father and kid and J2 is between trio as a whole.  For tri J1 is max number of incompatibilities in the
			trio as a whole. ))
		{
			if(ALREADY_HAS(Af, -a/--analysis-type)) {
				if(ARGS_EQ(ret->Xdim)) { int i;
					ret->Maxes = (int *)calloc(ret->Xdim,sizeof(int));
					for(i=0;i<ret->Xdim;i++)  {
						ret->Maxes[i] = GET_INT;
					}
				}
			}
		}
		
		if(MULT_USE_OPTION(
			mcoptsf,
			M,
			monte-carlo,
			S1 S2 S3 S4 R1 J1 J2,
			specify a set of Monte Carlo reps to be done,
			This option may be used multiple times---each time it specifies a set of Monte Carlo simulations to do.
			S1 is the name of the population that should be the numerator of the log-likelihood ratio. S2 is the denominator.  S3 is the
			population---i.e. the true relationship T.  S4 is the name of the population from which genotypes should be simulated 
			for an importance sampling approach---i.e. it defines the importance sampling distribution which in turn will---with S3---define
			the importance weights.  If S4 is the word NULL then vanilla Monte Carlo will be pursued instead of importance sampling.  R1 is the critical value of Lambda.  
			If it is non-zero then just the probability that lambda is less than R1 will be printed.  If it is zero then each replicate will be printed so the 
			distribution of lambdas or of the importance weights can be inspected.  J1 is the number of monte carlo samples to take.  This option should only
			be given after all the -p/--pop options have been given.  J2 is a flag that says whether or not the distribution of Lambdas should be printed.  J2==0 means
			do not print.  In this case the program will just return the estimated probability that Lambda < Lambda_c.  J2==1 means do print.  , 
			MAX_NUM_MCOPTS))
		{
			if(ARGS_EQ(7)) {
				char tempstr[MAX_POP_NAME];
				
				/* allocate memory and do some variable initialization */
				ret->MC[ret->NumMC] = (mc_struct *)malloc(sizeof(mc_struct));
				ret->MC[ret->NumMC]->N=0;
				ret->MC[ret->NumMC]->LC=0.0;
				ret->MC[ret->NumMC]->IS=-1;
				
				
				/* get logl numerator */
				GET_STR(tempstr);
				if(IdxFromPopName(tempstr,ret)==-1) {
					fprintf(stderr,"Error! Unrecognized population name %s used for argument S1 in -M/--monte-carlo.  Exiting...\n",tempstr);
					exit(1);
				}
				else {
					ret->MC[ret->NumMC]->LN = IdxFromPopName(tempstr,ret);
				}
				/* get logl denomninator */
				GET_STR(tempstr);
				if(IdxFromPopName(tempstr,ret)==-1) {
					fprintf(stderr,"Error! Unrecognized population name %s used for argument S2 in -M/--monte-carlo.  Exiting...\n",tempstr);
					exit(1);
				}
				else {
					ret->MC[ret->NumMC]->LD = IdxFromPopName(tempstr,ret);
				}
				/* get true relationship population */
				GET_STR(tempstr);
				if(IdxFromPopName(tempstr,ret)==-1) {
					fprintf(stderr,"Error! Unrecognized population name %s used for argument S3 in -M/--monte-carlo.  Exiting...\n",tempstr);
					exit(1);
				}
				else {
					ret->MC[ret->NumMC]->T = IdxFromPopName(tempstr,ret);
				}
				/* get the importance sampling population */
				GET_STR(tempstr);
				if(strcmp("NULL",tempstr)==0) {
					ret->MC[ret->NumMC]->IS = -1;
				}
				else {
					if(IdxFromPopName(tempstr,ret)==-1) {
						fprintf(stderr,"Error! Unrecognized population name %s used for argument S3 in -M/--monte-carlo.  Exiting...\n",tempstr);
						exit(1);
					}
					else {
						ret->MC[ret->NumMC]->IS = IdxFromPopName(tempstr,ret);
					}
				}
				/* then get lambda_c and N */
				ret->MC[ret->NumMC]->LC = GET_DUB;
				ret->MC[ret->NumMC]->N = GET_INT;
				ret->MC[ret->NumMC]->FullDist = GET_INT;

				
				/* increment the counter of the number of MC experiments */
				ret->NumMC++;
			}
			
		}
		
	END_OPT_LOOP
	
	
		ret->NumAinX = (int *)calloc(ret->XNum,sizeof(int));
	/* figure out how many A's are in each X state */
	for(i=0;i<ret->NumA;i++)  {
		ret->NumAinX[ Xstate(i,ret->aType) ]++;
	}
	/* allocate memory to the indexes of those states */
	ret->IdxAinX = (int **)calloc(ret->XNum,sizeof(int *));
	for(i=0;i<ret->XNum;i++)  {
		ret->IdxAinX[i] = (int *)calloc(ret->NumAinX[i], sizeof(int));
	}
	
	/* allocate to temporary counters */
	tempa = (int *)calloc(ret->XNum,sizeof(int));
	for(i=0;i<ret->NumA;i++) {
		ret->IdxAinX[Xstate(i,ret->aType)][ tempa[ Xstate(i,ret->aType) ]++ ] = i;
	}
	free(tempa);
	
	
	
	return(ret);
	
}	

/*! Prints out a summary of all the parameters that get set according to the analysis type.
This is primarily here for error checking */
void PrintATypeSummary(pfr_forback_data *D)
{
	int i,j;
	
	printf("SUMMARY: Loci : %d\n",D->L);
	printf("SUMMARY: NumPops : %d\n",D->NumPops);
	printf("SUMMARY: aType : %d\n",D->aType);
	printf("SUMMARY: Xdim : %d\n",D->Xdim);
	printf("SUMMARY: XNum : %d\n",D->XNum);
	printf("SUMMARY: NumAinX : ");
	for(i=0;i<D->XNum;i++)  {
		printf(" %d",D->NumAinX[i]);
	}
	printf("\n");
	for(i=0;i<D->XNum;i++)  {
		printf("SUMMARY:  XKeys :  %d  : ",i);
		for(j=0;j<D->Xdim;j++)  {
			printf(" %d",D->Xkeys[i][j]);
		}
		printf("\n");
	}
	for(i=0;i<D->XNum;i++)  {
		printf("SUMMARY:  IdxAinX :  %d  : ",i);
		for(j=0;j<D->NumAinX[i];j++)  {
			printf(" %d",D->IdxAinX[i][j]);
		}
		printf("\n");
	}
}


/*! Summarizes the information in a locus's PA field and allocates memory to the PX
	and PAinX fields and puts the appropriate info in there.  The appropriate info depends
	on aType, but aType dictates the NumAinX and the IdxAinX and the XNum fields, and those
	are all we need to get these summaries.  */
void CompleteLocusInfo(pfrfb_locus *L, int XNum, int *NumAinX, int **IdxAinX)
{
	int i,j;
		
	/* allocate memory */
	L->PX = (double *)calloc(XNum,sizeof(double));
	L->PAinX = (double **)calloc(XNum,sizeof(double *));
	for(i=0;i<XNum;i++)  {
		L->PAinX[i] = (double *)calloc(NumAinX[i],sizeof(double));
	}
	
	/* now put probabilities where they need to be, and sum them into L->PX */
	for(i=0;i<XNum;i++)  {
		for(j=0;j<NumAinX[i];j++)  {
			L->PX[i] += L->PA[ IdxAinX[i][j] ];
			L->PAinX[i][j] = L->PA[ IdxAinX[i][j] ];
		}
	}
	
	/* and then normalize the PAinX elements */
	for(i=0;i<XNum;i++)  {
		for(j=0;j<NumAinX[i];j++)  {
			L->PAinX[i][j] /= L->PX[i];
		}
	}
}


void FillAllLocusInfo(pfr_forback_data *D)
{
	int i,j;
	
	/* cycle over pops */
	for(i=0;i<D->NumPops;i++)  {
		for(j=0;j<D->L;j++)  {
			CompleteLocusInfo(D->Pops[i]->Loci[j], D->XNum, D->NumAinX, D->IdxAinX);
		}
	}
}

void PrintPopLocSummaries(pfr_forback_data *D)
{
	int i,j,k,l;
	
	for(i=0;i<D->NumPops;i++)  {
		for(j=0;j<D->L;j++) {
			printf("POP_LOC_SUMM : Pop= %d  NameOfPop= %s  Loc= %d\n",i,D->Pops[i]->Name,j);
			printf("POP_LOC_SUMM : Prob of XStates : ");
			for(k=0;k<D->XNum;k++)  {
				printf("  %f",D->Pops[i]->Loci[j]->PX[k]);
			}
			printf("\n");
			for(k=0;k<D->XNum;k++)  {
				printf("POP_LOC_SUMM : Prob of A states given X= %d :",k);
				for(l=0;l<D->NumAinX[k];l++)  {
					printf("  %f",D->Pops[i]->Loci[j]->PAinX[k][l]);
				}
				printf("\n");
			}
		}
	}
}


/*! This function is a pile of crap.  A vestige from when I tried to deal with things differently */
int BigHash(int *x, int *m, int l, analtype aType) 
{
	int A,B1,B2,C,v1,v2;
	
	switch(aType) {
		case(TRI):
			return(x[0]);
			break;
			
		case(PATRI):
		case(MATRI):
			return(x[0]*ECA_MIN(l,m[1]+1) - (x[0]*(x[0]-1))/2 + x[1]);
			break;
		case(MAPATRI):
			v2 = ECA_MIN(l,m[2]+1);
			v1 = ECA_MIN(l,m[1]+1);
			A = x[2] - ECA_MAX(x[0],x[1]);
			B1 = (x[0]+1)*(v2-x[0]+1);
			if(x[0]>x[1]-2) {
				B2=0;
			}
			else {
				B2 = (x[1]-1-x[0])*(v2+1) - ( (x[1]-1)*(x[1]))/2 + ( (x[0]+1)*(x[0]+2))/2;
			}
			C = v2*( (x[0]*(x[0]-1))/2)  - ( (x[0]-1)*x[0]*(2*(x[0]-1)+1) )/6 + ( (v2+1)*v1 - (v1*(v1-1))/2) * x[0] + (
					( (x[0]-1)*x[0]*(2*(x[0]-1)+1))/6  + 3 * ( x[0]*(x[0]-1))/2 + 2*x[0] )/2;
			return(A+B1+B2+C);
		default:
			fprintf(stderr, "Error! Dropped out of switch in BigHash unexpectedly.  Exiting...\n");
			exit(1);
	}
	
	
	return(-9999);
}


void MaPaTriForward(pfr_forback_data *D, pfr_pop_struct *Pop)
{
	int i,l,y0,y1,y2,t;
	int x0,x1,x2,x;
	int L = D->L;
	int Xdim = D->Xdim;
	int XNum = D->XNum;
	int Maxes[MAX_XSTATES];  /* put this in the stack.  might be faster */
	double ****P;
	pfrfb_locus *Loc;
	int **Xkeys = D->Xkeys;

#define MAX_VAL(X,Y) ECA_MIN(X,Maxes[Y]+1)

	/* allocate memory to D->dim3PY */
	Pop->dim3PY = (double ****)calloc(L+1,sizeof(double***));
	P = Pop->dim3PY;
	
	/* transfer the Maxes */
	for(i=0;i<Xdim;i++)  {
		Maxes[i] = D->Maxes[i];
	}
	
	/* allocate to the starting state */
/*	P[0] = (double ***)calloc(1,sizeof(double **));
	P[0][0] = (double **)calloc(1,sizeof(double *));
	P[0][0][0] = (double *)calloc(1,sizeof(double));
	P[0][0][0][0] = 1.0;  /* all mass on this starting state */
	
	
	/* take a pass over all t and allocate necessary memory */
	/* here we cycle over "time" t.  At time t=1, we are dealing with locus 1, but that locus has
	subscript 0.  Hence the l=t-1 later on.   Note that at some point I can try allocating these in such
	a way that I don't waste memory on values that I will never reach, but I won't worry about that for now. */
	for(t=0;t<=L;t++)  {
		P[t] = (double ***)calloc( MAX_VAL(t,0)+1, sizeof(double **));
		for(y0=0;y0<=MAX_VAL(t,0);y0++)  {
			P[t][y0] = (double **)calloc(MAX_VAL(t,1)+1, sizeof(double *));
			for(y1=0;y1<=MAX_VAL(t,1)+1;y1++)  {
				P[t][y0][y1] = (double *)calloc(MAX_VAL(t,2)+1, sizeof(double));
			}
		}
	}
	
/*printf("Done allocating memory \n");*/
	
	/* now, let us cycle over t=1 to L and compute the necessary probabilities iteratively */
	P[0][0][0][0] = 1.0;  /* Starting condition: all mass on this starting state */
	for(t=1;t<=L;t++)  {

		l=t-1;
		Loc = Pop->Loci[l];
		/* cycle over all places we could have been in the previous time step (t-1) */
		for(x0=0;x0<=MAX_VAL(t-1,0);x0++)  {
			for(x1=0;x1<=MAX_VAL(t-1,1);x1++) {
				for(x2=ECA_MAX(x0,x1);x2<=MAX_VAL(t-1,2);x2++)  {
				
/*printf("STARTING POINT: t= %d  y[t-1] =  %d %d %d  Prob= %.20e\n",t,x0,x1,x2,P[t-1][x0][x1][x2]); */
				

					/* then from such a place in t-1, cycle over all the places we could get to in t. */
					for(x=0;x<XNum;x++)  {
						/* We'll name that place using y0, y1, and y2. */
						y0 = ECA_MIN(x0 + Xkeys[x][0], Maxes[0]+1);
						y1 = ECA_MIN(x1 + Xkeys[x][1],Maxes[1]+1);
						y2 = ECA_MIN(x2 + Xkeys[x][2],Maxes[2]+1);
						
						/* the probability of ending up there is the prob of starting in the starting place and taking a certain route */
						P[t][y0][y1][y2] += P[t-1][x0][x1][x2] * Loc->PX[x];
/* printf("\tx-state= %d-%d-%d  leading to: %d-%d-%d   adding:  %.20e  resulting in: %.20e\n",Xkeys[x][0], Xkeys[x][1], Xkeys[x][2], y0, y1, y2, P[t-1][x0][x1][x2] * Loc->PX[x],
																					P[t][y0][y1][y2]); */
					}
				}
			}
		}
		
	}
	
#undef MAX_VAL(X,Y)
} 



void PrintMaPaTriForwardResults(pfr_forback_data *D, pfr_pop_struct *Pop)
{
	int i,t;
	int x0,x1,x2;
	int L = D->L;
	int Xdim = D->Xdim;
	int Maxes[MAX_XSTATES];  /* put this in the stack.  might be faster */
	double ****P;
	double **Marg;  /* to hold the marginal probs for each dimension */
	double sum;
	double TotIn=0.0;  /* TotIn is for recording the total prob that a trio is non-excluded */ 
	
	#define MAX_VAL(X,Y) ECA_MIN(X,Maxes[Y]+1)
	
	/* transfer the Maxes */
	for(i=0;i<Xdim;i++)  {
		Maxes[i] = D->Maxes[i];
	}
	
	/* allocate memory */
	Marg=(double **)calloc(D->Xdim,sizeof(double *));
	for(i=0;i<D->Xdim;i++)  {
		Marg[i] = (double *)calloc(MAX_VAL(L,i)+1,sizeof(double));
	}
	
	/* set P to address of the forback results */
	P = Pop->dim3PY;


	
	/* set t to the final locus */
	t=L;
	/* cycle over _all_ values of the number of Mendelian incompatibilities of various types */
	for(x0=0;x0<=MAX_VAL(t,0);x0++)  {
		for(x1=0;x1<=MAX_VAL(t,1);x1++) {
			for(x2=ECA_MAX(x0,x1);x2<=MAX_VAL(t,2);x2++)  {
				printf("FINAL_PROBS: pop= %s   t= %d    state=  %d %d %d  Prob= %.20e\n",Pop->Name,t,x0,x1,x2,P[t][x0][x1][x2]);
				Marg[0][x0] += P[t][x0][x1][x2];
				Marg[1][x1] += P[t][x0][x1][x2];
				Marg[2][x2] += P[t][x0][x1][x2];
			}
		}
	}
	
	/* now we print the marginals */
	printf("MARG_PROB:PopName\tX\tx\tPX.eq.x\tPX.leq.x\tPX.gt.x\n");
	sum=0.0;
	for(i=0;i<MAX_VAL(t,0);i++)  {
		sum += Marg[0][i];
		printf("MARG_PROB:%s\tx0\t%d\t%.8e\t%.8e\t%.8e\n",Pop->Name,i,Marg[0][i],sum,1.0-sum);
		
	}
	printf("MARG_CHK: Pop= %s  Prob  x0  >=  %d   is  %.10e\n",Pop->Name,i,Marg[0][i]);  /* this last one is the MAX_VAL(t,0).  It should be 1-sum at this point
																		so I include it in there just to be able to check that everything is
																		working correctly */

	sum=0.0;
	for(i=0;i<MAX_VAL(t,1);i++)  {
		sum += Marg[1][i];
		printf("MARG_PROB:%s\tx1\t%d\t%.8e\t%.8e\t%.8e\n",Pop->Name,i,Marg[1][i],sum,1.0-sum);
		
	}
	printf("MARG_CHK: Pop= %s  Prob  x1  >=  %d   is  %.10e\n",Pop->Name,i,Marg[1][i]);
	
	sum=0.0;
	for(i=0;i<MAX_VAL(t,2);i++)  {
		sum += Marg[2][i];
		printf("MARG_PROB:%s\tx2\t%d\t%.8e\t%.8e\t%.8e\n",Pop->Name,i,Marg[2][i],sum,1.0-sum);
		
	}
	printf("MARG_CHK: Pop= %s  Prob  x2  >=  %d   is  %.10e\n",Pop->Name,i,Marg[2][i]);


	/* at the very end we print the total probs */
	/* cycle over all _allowed_ values of the number of Mendelian incompatibilities */
	for(x0=0;x0<MAX_VAL(t,0);x0++)  {
		for(x1=0;x1<MAX_VAL(t,1);x1++) {
			for(x2=ECA_MAX(x0,x1);x2<MAX_VAL(t,2);x2++)  {
				TotIn += P[t][x0][x1][x2];
			}
		}
	}
	printf("TOTAL_PROB_OF_TRIO_INCLUSION: %s  %.30e\n",Pop->Name,TotIn);
	printf("TOTAL_PROB_OF_TRIO_EXCLUSION: %s  %.30e\n",Pop->Name,1.0-TotIn);
	
#undef MAX_VAL(X,Y)
} 


double  MaPaTriBackward_Once(pfr_forback_data *D, mc_struct *MC)
{
	int i,y0,y1,y2,t;
	int x0,x1,x2,x,drawn_x,drawn_A;
	int L = D->L;
	int Xdim = D->Xdim;
	int XNum = D->XNum;
	int Maxes[MAX_XSTATES];  /* put this in the stack.  might be faster */
	double Probs[MAX_XSTATES];  /* this is where we will normalize the probs of the backward steps, so we can sample from them. */
	pfrfb_locus *Loc;
	int **Xkeys = D->Xkeys;
	double sum,rando;
	double ****T, ****IS;
	double Weight, Lambda, LNv=1.0, LDv=1.0, Tv=1.0, ISv=1.0;  /* for accumulating products of values over loci */
	double full_sum;
	
	/* assign values to T and IS */
	T = D->Pops[ MC->T ]->dim3PY;
	if(MC->IS > -1) {
			IS = D->Pops[ MC->IS ]->dim3PY;
	}
	else {
		IS = T;
	}
	
			
	/* transfer the Maxes */
	for(i=0;i<Xdim;i++)  {
		Maxes[i] = D->Maxes[i];
	}
	
	/* here we have to sample our "first" state which is the one at time T.  We cycle over all possible states once
	to get the sum, then again to sample from it: */
	t = L;
	full_sum=0.0;
	for(x0=0;x0<=Maxes[0];x0++)  {
		for(x1=0;x1<=Maxes[1];x1++) {
			for(x2=ECA_MAX(x0,x1);x2<=Maxes[2];x2++)  {
				full_sum += IS[t][x0][x1][x2];
			}
		}
	}
	
/*printf("FULLSUM is %e\n",full_sum);*/


	y0=y1=y2=-1;
	rando = ranf() * full_sum;

/*printf("RANDO_IS is %e\n",rando); */
	
	sum=0.0;
	for(x0=0;x0<=Maxes[0];x0++)  {
		for(x1=0;x1<=Maxes[1];x1++) {
			for(x2=ECA_MAX(x0,x1);x2<=Maxes[2];x2++)  {
				sum += IS[t][x0][x1][x2];
				if(sum>rando) {
					y0=x0;
					y1=x1;
					y2=x2;
					x0=Maxes[0]+1;  /* ensure that we get out of all the loops */
					x1=Maxes[1]+1;
					x2=Maxes[2]+1;
				}
			}
		}
	}
	if(y0==-1 || y1==-1 || y2==-1) { /* dealing with case where ranf() returns 1.0 identically */
		y0 = Maxes[0];
		y1 = Maxes[1];
		y2 = Maxes[2];
	}
	/* now, y0, y1, y2 is the current state, and for each t from T-1 to 1 those variable will be updated */

	if(MC->IS > -1)  { /* if we are going to be doing importance sampling, we have to compute the Prob of the realization under T and IS */
		ISv *= IS[t][y0][y1][y2] / full_sum;
		
		/* and then to compute the T part we have to sum over all allowable values in T to normalize properly */
		full_sum = 0.0;
		for(x0=0;x0<=Maxes[0];x0++)  {
			for(x1=0;x1<=Maxes[1];x1++) {
				for(x2=ECA_MAX(x0,x1);x2<=Maxes[2];x2++)  {
					full_sum += T[t][x0][x1][x2];
				}
			}
		}
		Tv *= T[t][y0][y1][y2] / full_sum;

	}
	

	/* now we cycle backward over t, and each time we select at random the state that it could have come from.  Note that the 
	route it took to get there is the Xstate of the individual, and we can use that to randomly sample a genotypic configuration. */
	for(t=L-1;t>=0;t--) {
		if(MC->IS > -1) {
			Loc = D->Pops[ MC->IS ]->Loci[t];
		}
		else {
			Loc = D->Pops[ MC->T ]->Loci[t];
		}
		sum = 0.0;
		for(x=0;x<XNum;x++)  {  /* cycle over possible X states of Locus with index t */
			x0 = y0-Xkeys[x][0];
			x1 = y1-Xkeys[x][1];
			x2 = y2-Xkeys[x][2];
			
			if(x0>=0 && x1>=0 && x2>=0 && x0<=t && x1<=t && x2<=t) {
				Probs[x] = Loc->PX[x] * IS[t][ x0 ][ x1 ][ x2 ];  /* this is the prob that it came from y0-x[0], etc in this step */
			}
			else {
				Probs[x] = 0.0;
			}
			sum += Probs[x];
		}
		
		/* then sample the x from the Probs array */
		drawn_x = IntFromDoublesRV(Probs,sum,XNum);
		
		/* if we are importance sampling, we have to multiply the prob of this onto ISv, and compute the prob under T. */
		if(MC->IS > -1) {
			ISv *= Probs[drawn_x] / sum;
			
			/* then compute the prob under T */
			Loc = D->Pops[ MC->T ]->Loci[t];
			sum = 0.0;
			for(x=0;x<XNum;x++)  {  /* cycle over possible X states of Locus with index t */
				x0 = y0-Xkeys[x][0];
				x1 = y1-Xkeys[x][1];
				x2 = y2-Xkeys[x][2];
				
				if(x0>=0 && x1>=0 && x2>=0 && x0<=t && x1<=t && x2<=t) {
					Probs[x] = Loc->PX[x] * T[t][ x0 ][ x1 ][ x2 ];  /* this is the prob that it came from y0-x[0], etc in this step */
				}
				else {
					Probs[x]=0.0;
				}
				sum += Probs[x];
			}
			
			Tv *= Probs[drawn_x] / sum;
		}
		
		/* now, down here, we have to sample the genotypic configuration given the drawn_x, and compute the numerator and the denominator of
		lambda for each.  */
		if(MC->IS > -1) {
			Loc = D->Pops[ MC->IS ]->Loci[t];
		}
		else {
			Loc = D->Pops[ MC->T ]->Loci[t];
		}
		
		drawn_A = IntFromProbsRV(Loc->PAinX[drawn_x],0,D->NumAinX[drawn_x]);
		
		/* compute the importance weights if need be */
		if(MC->IS > -1) {
			ISv *= D->Pops[ MC->IS ]->Loci[t]->PAinX[drawn_x][drawn_A];
			Tv *=  D->Pops[ MC->T ]->Loci[t]->PAinX[drawn_x][drawn_A];
		}
		
		/* multiply terms onto the numerator and denominator of the likelihood ratio */
		LNv *= D->Pops[ MC->LN ]->Loci[t]->PA[ D->IdxAinX[drawn_x][drawn_A] ];
		LDv *= D->Pops[ MC->LD ]->Loci[t]->PA[ D->IdxAinX[drawn_x][drawn_A] ];
		
		
		/* at the end of all this, we now set y0, y1, and y2 to the current values  */
		y0 -= Xkeys[drawn_x][0];
		y1 -= Xkeys[drawn_x][1];
		y2 -= Xkeys[drawn_x][2];
		
	}
	
	/* and now that we have gotten down here, we just have to figure out how to return the result.
		Basically:  if MC->FullDist == 1 then we  print Lambda and the importance weight.  Otherwise we don't.  
		As far as returning the value---if Lambda > MC->LC then we return the importance weight.  Otherwise zero. */
	Lambda = log(LNv) - log(LDv);
	if(MC->IS > -1) {
		Weight = Tv / ISv;
	}
	else {
		Weight = 1.0;
	}
	if(MC->FullDist) {
		printf("MC_IMP   :   %.10f   T    %.10e\n",Lambda,Weight);
	}	
	
	if(Lambda > MC->LC) {
		return(Weight);
	}
	else {
		return(0.0);
	}
	
} 


void Do_MC_Experiment(pfr_forback_data *D, mc_struct *MC) 
{
	int i;
	char temp[MAX_POP_NAME];
	dval *Term;
	int NumUnInformative = 0;
	
	/* Allocate to the dval.  This also initializes it */
	Term = AllocDval(0,0,0);
	
	
	/* first thing. Print out a banner */
	if(MC->IS==-1) {
		sprintf(temp,"NULL");
	}
	else {
		sprintf(temp,"%s",D->Pops[ MC->IS ]->Name);
	}
	printf("STARTING_MC_SET: LN= %d %s  LD= %d %s  T= %d %s  IS= %d %s  NumReps= %d  LambdaC= %f\n",MC->LN,D->Pops[ MC->LN ]->Name,
		MC->LD,D->Pops[ MC->LD ]->Name,  MC->T,D->Pops[ MC->T ]->Name,  MC->IS,temp, MC->N,MC->LC);
		
		
	/* do the Monte Carlo iterations */
	for(i=0;i<MC->N;i++)  {
		Term->v = MaPaTriBackward_Once(D,MC);
		if(Term->v == 0.0) {  /* if it is identically zero, then the value of lambda was below lambda crit */
			NumUnInformative++;
		}
		IncrementDval(Term);
	}
	
	/* now print the results. */
	printf("COMPLETED_MC_SET: LN= %d %s  LD= %d %s  T= %d %s  IS= %d %s  NumReps= %d  NumInformative= %d  LambdaC= %f  Ave= %e  Var= %e  SE= %e\n",MC->LN,D->Pops[ MC->LN]->Name,
		MC->LD,D->Pops[ MC->LD ]->Name,  MC->T,D->Pops[ MC->T ]->Name,  MC->IS,temp, MC->N, MC->N-NumUnInformative, MC->LC, Term->Ave,  Term->Var,  sqrt(Term->Var) ); 
}

int main (int argc, char *argv[]) 
{
	int i;
	pfr_forback_data *Dat;

	Dat = GetPFR_ForBackOptions(argc,argv);
	
	PrintATypeSummary(Dat);
	FillAllLocusInfo(Dat);
	PrintPopLocSummaries(Dat);
	
	/* do the forward calculations for all the populations */
	for(i=0;i<Dat->NumPops;i++) {
		MaPaTriForward(Dat, Dat->Pops[i]);
		PrintMaPaTriForwardResults(Dat, Dat->Pops[i]);
	}
	
	/* seed the random number generator */
	SeedFromFile("pfr_forback_seeds");
	
	/* now, cycle over the MC experiments and do each */
	for(i=0;i<Dat->NumMC;i++)  {
		Do_MC_Experiment(Dat, Dat->MC[i]);
	}
	
	SeedToFile("pfr_forback_seeds");
	
    return(0);
}
