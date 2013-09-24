

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
#include "pbt_C_fb.h"
#include "uthash.h"
#include "pfr_read_genos.h"
#include "pfr_pedigree_spec.h"
#include "pbt_highlevel.h"



/*************************************************************************************/

/* given a starting S-value (which could be null if it is all zeroes)
 add v to the elements and return the result in the output s
*/
void SumSV(int *Sprev, int *v, int *Sout, int d)
{
	int i;
	
	if(Sprev==NULL) {
		for(i=0;i<d;i++)  {
			Sout[i] = v[i];
		}
	}
	else {
		for(i=0;i<d;i++)  {
			Sout[i] = v[i] + Sprev[i];
		}
	}
}

/* allocate the innards of an fb_cell */
void AllocFB_innards( struct fb_cell *f, int d, int Nx)
{
	f->s = (int *)calloc(d,sizeof(int));
	f->up = (struct fb_cell **)calloc(Nx, sizeof(struct fb_cell *));
	f->down = (struct fb_cell **)calloc(Nx, sizeof(struct fb_cell *));
	f->dp = (double *)calloc(Nx,sizeof(double));
	f->p = 0.0;
	f->fake_p = 0.0;
}



/* allocate the innards of an fb_cell */
void FreeFB_innards( struct fb_cell *f)
{
	free(f->s);
	free(f->up);
	free(f->down);
	free(f->dp);
}




void AllocTrioGenotypeInnards(trio_genotype *GT, int L)
{
	GT->path = (int *)calloc(L,sizeof(int));
	GT->agx = (int *)calloc(L,sizeof(int));
	GT->a = (int *)calloc(L,sizeof(int));
}

void FreeTrioGenotypeInnards(trio_genotype *GT, int L)
{
	free(GT->path);
	free(GT->agx);
	free(GT->a);
}


/* 
 Given d and smax, find the index in the big vdex array that gives 
 you a particular value of x.  This can also be used to 
 figure out an index for an a-state given the underlying y-states.
*/
int ReturnVDex(int d, int *smax, int *s)
{
	int i;
	int prod=1;
	int sum=s[d-1];
	
	for(i=d-2;i>=0;i--) {
		prod *= (smax[i+1]+1);
		sum += s[i]*prod;
	}
	return(sum);
}

/* return 1 if every element of x is <= the
 corresponding element in s.  Each of those vectors
 is of length d.
*/
int InSDown(int *x, int *s, int d)
{
	int i;
	int ret = 1;
	
	for(i=0;i<d;i++) {
		ret *= (x[i]<=s[i]);
	}
	
	return(ret);
}

/* this is the function that does the forward step.  vdex is an array 
   that is used to index the different values of s below smax that have already been
   seen.  It is accessed with a macro that depends on knowing what smax is and how
   many components it has.   The first time this is called, vdex should be 
	set to a bunch of -1's. 
 
 VDexSet should be set to 0 for the first time this function is executed and then
 set to 1 afterwards.  If it is set to one, then there shouldn't be any new 
 s values reached.  Actually, this is not really working now. Set to 0 always, and
 makd sure that vdex is always set to -1's.
 
 I may fix this later.  Basically I want to use it to ensure that the structure
 in each fb_struct is identical between collections.  For now that is ensured so
 long as there is nonzero probability for every X-state.
 
 I could easily change this so that I passed in VdexIn and if that were NULL
 I would allocate memory to it and the point vdex as a local variable to that.
 If it were not null then I could just point vdex to it.  And then I should probably 
 also pass Kt in as well and have that be common across all collections.  But that
 could be a shrieking horror to change at this point, so I won't worry about it for now.
 
 It turns out that if we are computing conditional probs on the basis of genotypes and then
 dividing by the prob of SDown, it doesn't really matter.
 
 Note that FX are the "Fake Probs"
 
 */
fb_struct *ForwardStep(RunPars *R, double **X, double **FX, int *vdex, int VDexSet)
{
	int t,i,n,j,idx;
	int temp_v[MAX_X_STATES];
	fb_struct *ret;
	double S_up = 0.0;
	
	/* some memory allocation */
	ret = (fb_struct *)malloc(sizeof(fb_struct));
	ret->Kt = (int *)calloc(R->L,sizeof(int));
	ret->FB = (struct fb_cell **)calloc(R->L,sizeof(struct fb_cell *));
	
	
	/* now initialize for t=0 */
	n=0;
	
	
	/************************* INITIALIZE AT LOCUS (TIME) 0 ****************************/
	/* first time through we just count up the new S-cells */
	for(i=0;i<R->Nx;i++) {
		if( InSDown(R->v_vecs[i],R->smax,R->d) ) {
			if( vdex[ ReturnVDex(R->d,R->smax,R->v_vecs[i])] == -1 ) {
				if(VDexSet) {
					fprintf(stderr,"Error! vdex element == -1 while initializing t=0 and VDexSet != 0. Exiting\n");
					exit(1);
				}
				else {
					n++;
				}
			}
		}
	}
	
	/* now, here n holds the number of new fb cells.  So allocate to them and and then re-run and fill them */
	ret->FB[0] = (struct fb_cell *)calloc(n,sizeof(struct fb_cell));
	n=0;
	for(i=0;i<R->Nx;i++) {
		if( InSDown(R->v_vecs[i],R->smax,R->d) ) {
			if( vdex[ ReturnVDex(R->d,R->smax,R->v_vecs[i])] == -1 ) {
				if(VDexSet) {
					fprintf(stderr,"Error! vdex element == -1 while initializing t=0 and VDexSet != 0. Exiting\n");
					exit(1);
				}
				else {
					AllocFB_innards(&(ret->FB[0][n]),R->d, R->Nx);
					ret->FB[0][n].idx = n;
					ret->FB[0][n].p = X[0][i];
					if(R->HasFakeProbs) {
						ret->FB[0][n].fake_p = FX[0][i];
					}
					else {
						ret->FB[0][n].fake_p = X[0][i];
					}
			
					SumSV(NULL,R->v_vecs[i],ret->FB[0][n].s,R->d);
					vdex [ ReturnVDex(R->d,R->smax,R->v_vecs[i]) ] = n;
					/*printf("INIT_FB: just initialized cell %d at time %d.  prob= %e  S= %d %d %d \n",n,0,ret->FB[0][n].p,
						   ret->FB[0][n].s[0],ret->FB[0][n].s[1],ret->FB[0][n].s[2]); */
					n++;
				}
			}
		}
		else {
			S_up += X[0][i];
			//printf("S_UP_PROB_FB: Just added %e to S_up while initializing at locus 0. S_up is now %e\n",X[0][i],S_up);
			
		}
	}
	
	/* finally, record the number of entries here at time 0 */
	ret->Kt[0] = n;
	
	/********** DONE INITIALIZING *******************/
	
	
	/************* NOW CYCLE OVER LOCI FROM 1 to L-1 ******************/
	for(t=1;t<R->L;t++) {  /* cycle over loci */
		
		/* first time through we just count up the new S-cells */
		n=0;
		for(j=0;j<ret->Kt[t-1];j++) {  /* cycle over the states in the previous time step */
			for(i=0;i<R->Nx;i++) {
				
				SumSV(ret->FB[t-1][j].s, R->v_vecs[i], temp_v, R->d);  /* now temp_v holds the resulting v-vector */
				if( InSDown(temp_v,R->smax,R->d) ) {
					if( vdex[ ReturnVDex(R->d,R->smax,temp_v)] == -1 ) {
						if(VDexSet) {
							fprintf(stderr,"Error! vdex element == -1 while initializing t=0 and VDexSet != 0. Exiting\n");
							exit(1);
						}
						else {
							vdex[ ReturnVDex(R->d,R->smax,temp_v)] = ret->Kt[t-1] + n; /* we set it here so that it doesn't get double counted */
							n++;
						}
					}
				}
			}
		}
		
		/* now we set ret->Kt[t] and allocate to ret->FB[t] */
		ret->Kt[t] = ret->Kt[t-1] + n;
		//printf("COUNT_FB: At time %d we have counted %d fb_cells in s_down\n",t,ret->Kt[t]);
		ret->FB[t] = (struct fb_cell *)calloc(ret->Kt[t], sizeof(struct fb_cell));
		for(j=0;j<ret->Kt[t];j++)  {
			AllocFB_innards(&(ret->FB[t][j]), R->d, R->Nx);
		}
		
		/* and here we cycle over everything again and accumulate some sums of probabilities and 
		 point some pointers around */
		for(j=0;j<ret->Kt[t-1];j++) {  /* cycle over the states in the previous time step */
			for(i=0;i<R->Nx;i++) {
				SumSV(ret->FB[t-1][j].s, R->v_vecs[i], temp_v, R->d);  /* now temp_v holds the resulting v-vector */
				if( InSDown(temp_v,R->smax,R->d) ) {
					if( vdex[ ReturnVDex(R->d,R->smax,temp_v)] == -1 ) {
						fprintf(stderr,"Error! vdex element == -1 after pre-counting states.  Shouldn't happen. Exiting\n");
						exit(1);
					}
					else {
						idx = vdex[ ReturnVDex(R->d,R->smax,temp_v)]; /* get the index of the state that this creates */
						ret->FB[t][idx].idx = idx;
						ret->FB[t][idx].p += ret->FB[t-1][j].p * X[t][i];
						if(R->HasFakeProbs) {
							ret->FB[t][idx].fake_p += ret->FB[t-1][j].fake_p * FX[t][i];
						}
						else {
							ret->FB[t][idx].fake_p += ret->FB[t-1][j].fake_p * X[t][i];
						}
						SumSV(ret->FB[t-1][j].s, R->v_vecs[i], ret->FB[t][idx].s, R->d);
						/*printf("CYCLE_FB: just added prob %e to cell %d at time %d.  Reaching forward from cell %d in previous step.  PrevS= %d %d %d   CurrS= %d %d %d\n",ret->FB[t-1][j].p * X[t][i],idx,t,j,
							ret->FB[t-1][j].s[0],ret->FB[t-1][j].s[1],ret->FB[t-1][j].s[2],
							  ret->FB[t][idx].s[0],ret->FB[t][idx].s[1],ret->FB[t][idx].s[2] ); */
						/* now set some pointers */
						ret->FB[t-1][j].up[i] = &(ret->FB[t][idx]);
						ret->FB[t][idx].down[i] = &(ret->FB[t-1][j]);
					}
				}
				else {  /* otherwise we add stuff up to S_up */
					S_up += ret->FB[t-1][j].p * X[t][i];
					/*printf("S_UP_PROB_FB: Just added %e to S_up while cycling at locus %d  j= %d and i= %d   previous_S= %d %d %d    S_up is now %e.  One minus it is %e\n",
						   ret->FB[t-1][j].p * X[t][i],  t,j,i,
						   ret->FB[t-1][j].s[0], ret->FB[t-1][j].s[1], ret->FB[t-1][j].s[2],
						   S_up,1.0-S_up); */

				}
				
			}
		}
	}
	
	ret->S_up =S_up;
	ret->S_down = 1.0 - S_up;
	
	return(ret);
}


/* this just frees up all the memory allocated to the innards on an FB_Vars struct so that after we
 do that we can safely free up the pointer itself */
void FreeFB_Vars_innards(RunPars *RP, fb_struct *F) 
{
	int i,j, L = RP->L;
	
	for(i=0;i<L;i++)  {
		for(j=0;j<F->Kt[i];j++)  {
			FreeFB_innards( &(F->FB[i][j]));	
		}
		free(F->FB[i]);
	}
	free(F->FB);
	free(F->Kt);
}


/* use this to allocate space to a collection that we can then use to put values conditioned on individual genotypes into.
 
 Note that, in anticipation of changing the dimensions of some of the variables when we have missing data, I am going to allocate
 more space than needed by default for some (i.e. 64 for some of the things, etc.  Since we only make a few of these, this should be just fine.
 
 */
Collection *AllocToCollection(RunPars *RP)
{
	int i,j;
	Collection *ret;
	
	ret = (Collection *)malloc(sizeof(Collection));
	ret->CollName = (char *)calloc(MAX_COLLECTION_NAME_LENGTH,sizeof(char));
	
	/* now allocate memory to the elements of the collection */
	ret->Xprobs = (double **)calloc(RP->L,sizeof(double *));
	ret->Xfakeprobs = (double **)calloc(RP->L,sizeof(double *));
	ret->Aprobs_M = (double **)calloc(RP->L,sizeof(double *));
	if(RP->Na) {  /* if we are also specifying the A's */
		ret->Aprobs = (double **)calloc(RP->L,sizeof(double *));
		ret->AprobsGX = (double ***)calloc(RP->L,sizeof(double **));
	}
	for(i=0;i<RP->L;i++) {
		ret->Xprobs[i] = (double *)calloc(RP->Nx,sizeof(double));
		ret->Xfakeprobs[i] = (double *)calloc(RP->Nx,sizeof(double));
		ret->Aprobs_M[i] = (double *)calloc(64,sizeof(double));
		if(RP->Na) {  /* if we are also specifying the A's */
			ret->Aprobs[i] = (double *)calloc(/*RP->Na*/ 64,sizeof(double));
			ret->AprobsGX[i] = (double **)calloc(RP->Nx,sizeof(double *));
			for(j=0;j<RP->Nx;j++) {
				ret->AprobsGX[i][j] = (double *)calloc(/*RP->NaInX[j]*/ 64,sizeof(double));
			}
		}
	}
	
	return(ret);
	
}



/*
 Given the probs of all the A states at all the loci and also all the info in the RunPars struct RP,
 this function fills out all the corresponding values (Xprobs, etc) in the output Collection, Out. 
 
 Note that memory has to be allocated to Out first.
 
 */
void XStuffFromAStuff(double **Aprobs, RunPars *RP, Collection *Out)  {
	int i,j,l;
	int Na = RP->Na;
	int *AinX = RP->AinX;
	int **aInX2A = RP->aInX2A;
	int L=RP->L;
	int Nx = RP->Nx;
	int *NaInX = RP->NaInX;
	
	/* cycle over all the loci and do this */
	for(l=0;l<L;l++)  {
		/* First thing, let's copy the Aprobs over */ 
		for(i=0;i<Na;i++)  {
			Out->Aprobs[l][i] = Aprobs[l][i];
		}
		
		/* then do the X probs */
		/* first we initialize to zero to get a sum */
		for(i=0;i<Nx;i++)  {
			Out->Xprobs[l][i] = 0.0;
		}
		/* then cycle over A states and sum them into the X states */
		for(i=0;i<Na;i++)  {
			Out->Xprobs[l][ AinX[i] ] += Aprobs[l][i];
		}
		
		/* and now we want to put the normalized probs into the AprobsGX field */
		for(i=0;i<Nx;i++) {
			for(j=0;j<NaInX[i];j++) {
				Out->AprobsGX[l][i][j] = Aprobs[l][ aInX2A[i][j] ] / Out->Xprobs[l][i];
			}
		}
	}
}



/* 
 Given a bunch of A probs at L loci in AprobsIn AND some observations on the Y-states underlying the 
 A-states, (the **Yobs),  and information about what those Ystates mean (in the RunPars struct RP),
 this function writes new values to AprobsOut in which it has conditioned on the observed Y states
 and renormalized so that they sum to one.
 
 If nothing is observed at a Yobs element it should be -1.  Otherwise it is 0, 1, or 2 for SNP genotypes. 
 
 Yobs[l][j] is the state of the j-th y at locus l.  
*/
void ConditionAprobsOnYobs(double **AprobsIn, int **Yobs, RunPars *RP, double **AprobsOut)
{
	int j,l,k;
	int Na = RP->Na;
	int L=RP->L;
	int **Ystates = RP->Ystates;
	int NY = RP->NY;
	double mult,normo;
	
	for(l=0;l<L;l++)  {
		normo=0.0;
		for(k=0;k<Na;k++)  {
			mult=1.0;
			for(j=0;j<NY;j++)  {
				if(Yobs[l][j]>=0 && Yobs[l][j] != Ystates[k][j]) {
					mult=1e-12;  /* looks like zeroes might cause problems, so must put a really small value in there. */
					break;
				}
			}
			AprobsOut[l][k] = mult * AprobsIn[l][k];
			normo+=AprobsOut[l][k];
		}
		/* now, normalize those dudes */
		for(k=0;k<Na;k++)  {
			AprobsOut[l][k] /= normo;
		}
	}
	
}




/*
 This function takes an array of L A-states in A and it returns the corresponding
 values of the Y-states into Y.  Ystates holds the key to what Y's each A implies.
 Usually that is in the RunPars.  NY is the number of Y states for each A state.  
 
 Mask is an array of length NY that tells whether the corresponding Y value for each
 A state should be returned (if it is a 0) or if the value -1 should be returned in its place (if 
 it is 1---i.e. mask==1 means don't report the Y-state).  This is useful if you want to
 create a Y that you will use to only condition on the kid of the trio for example.  
 
 
 You have to allocate space to Y outside of this function.  It should be of dimension
 [L][NY].
 
*/
void FillYsFromAs(int *A, int **Y, int L, int **Ystates, int NY, int *Mask)
{
	int i,j;

	for(j=0;j<L;j++)  {
		for(i=0;i<NY;i++)  {
			if(Mask[i]==0) {
				Y[j][i] = Ystates[A[j]][i];
			}
			else {
				Y[j][i] = -1;
			}
		}
	}
	
}




/* given a vector of NY Y's, and the info in RunPars RP, this 
 returns the index of the corresponding A state ONLY VALID WHEN
 THERE IS NO MISSING DATA!
*/
int GiveAFromY(int *y, RunPars *RP) {
	return(RP->AgivenY[ ReturnVDex(RP->NY,RP->ymaxes, y) ]);
}



/* given a vector of NY Y's and the info in RunPars RP, for the	
   case where there could be missing data, this returns the missing-data-A-state
   given the Y */
int GiveAFromY_M(int *y, RunPars *RP) {
	return(RP->AgivenY_M[ ReturnVDex(RP->NY,RP->ymaxes_M, y) ]);
}



/* here is a function to allocate to and return a whole array of these 
 A-states-with-missing data given the observed data on individuals in a trio.  
 y is the youth,  ma and pa are self explanatory. */
int *A_64_ArrayFrom3IndivCharVecs(char *y, char *pa, char *ma, RunPars *RP)
{
	int i,g[3];
	int L = RP->L;
	int *ret = (int *)calloc(L,sizeof(int));
	
	for(i=0;i<L;i++)  {
		g[0]=y[i];
		g[1]=pa[i];
		g[2]=ma[i];
		
		ret[i] = GiveAFromY_M(g,RP);
	}
	
	return(ret);
}




/* given the observed genotypes at a trio, A, (expressed as 0 to 63) at each locus and the associated
 probabilities for each in Aprobs_M, at L loci, return the prob of the trio genotype at all L loci 
*/
double UnconditionalGenoProb_A64(double **Aprobs_M, int *A, int L) 
{
	int i;
	double ret=1.0;
	
	
	/*printf("COMPING_POSTS_INSIDE_UnconditionalGenoProb_A64:  ");
	for(k=0;k<L;k++)  {
		printf("%d ",A[k]);
	}
	printf("\n");*/
	
	
	for(i=0;i<L;i++)  {
		//printf("PARTIAL_OVER_LOCI loc= %d     A[i]=  %d    Aprobs_M= %f\n",i,A[i], Aprobs_M[i][ A[i] ]);
		ret *= Aprobs_M[i][ A[i] ];
		//printf("PARTIAL_OVER_LOCI_AFTER_RET_PROD:  loc= %d     A[i]=  %d    Aprobs_M= %f\n",i,A[i], Aprobs_M[i][ A[i] ]);
	}
	
	return(ret);
}

/* given an A-state array having values for L loci (each can be 0 to 63 and hence can deal with missing data),
 and an FB_Vars struct that holds the PurePopTrioColls, 
 and a vector of length NUM_SPEC_PEDS that holds the prior probability of the trio falling into
 any of the pedigree relationships, and the Population from which these individuals are supposed
 to all be from, 
 compute the posterior probability of each pedigree category given that youth, pa, and ma are all from the
 same population.  This is the prior for that category times the trio prob given that category divided by the 
 sum of the prior times probability for all the categories.
 
 This also returns the log of the likelihood ratio of parental over non-parental (where the non-parental distribution 
 is a mixture of the non-parental categories) as the output variable LogL.
 
 Note that Pi has to be properly normalized and sum to one for the LogL to be correct. 
 
 Note that this only works right if the 0 category is the parental relationship.
 
 Note that if UseThisAsOutput is NULL, then this function will allocate new memory and return it.
 If it is not null, then it will just put the posterior probs into the space in UseThisAsOutput (and
 will also return a pointer to it).
 
 Two more output variables that are useful for my importance sampling scheme:  PrPar is the genotype probability of the trio
 under the parental hypothesis.  PrNonPar is the genotype probability of the trio under the non-parental hypothesis.
 
 SinglePopCollWithPedsOrdered is a flag.  If 1 it means that the collections in PP are those for a single population and the are simply
 indexed according to the pedigree index.  If it is 0, it means that PP is for all the pops and then Pop will be used 
 to get the right index for each collection out of it.  Basically we use 1 when we want to compute the prob of something 
 being simulated from CPFB and 0 when we want the prob based on, say HLV->PurePopTrioColls.
 */
double *TrioPostProbs(int *A_array, FB_Vars *PP, int SinglePopCollWithPedsOrdered, double *Pi, int Pop, double *LogL, double *UseThisAsOutput, double *PrPar, double *PrNonPar)
{
	int i,x;
	double normo=0.0;
	int L = PP->RP->L;
	double *ret;
	double LogLnumer, LogLdenom;
	
	
	if(UseThisAsOutput==NULL) {
		ret = (double *)calloc(NUM_SPEC_PEDS,sizeof(double));
	}
	else {
		ret = UseThisAsOutput;
	}
	for(i=0;i<NUM_SPEC_PEDS;i++)  {
		
		if(SinglePopCollWithPedsOrdered) {
			x=i;
		}
		else {
			x = PPTC_IDX(NUM_SPEC_PEDS,  i,  Pop);  /* this is the index of the Collection we want */
		}
		if(Pi[i] > 1e-50) {  /* if the prior prob is very small (like 0) we don't even worry about it */
			
						
			ret[i] =  Pi[i] *  UnconditionalGenoProb_A64( PP->Colls[x]->Aprobs_M, A_array, L);
		}
		else {
			ret[i] = 0.0;
		}
		normo += ret[i];
	}
	
	/* at this point, we can get the likelihood ratio by getting the numerator as the raw parental prob as ret[0]/Pi[0] and then getting
	 the denominator by subtracting ret[0] from normo, and then dividing that by (1-Pi[0]). */
	LogLnumer = ret[0]/Pi[0];
	LogLdenom = (normo - ret[0]) / (1.0-Pi[0]);
	
	*PrPar = LogLnumer;
	*PrNonPar = LogLdenom;
	
	
	*LogL = log(LogLnumer/LogLdenom);
	
	
	
	/* now normalize everything */
	for(i=0;i<NUM_SPEC_PEDS;i++)  {
		ret[i] /= normo;
	}
	
	return(ret);
}


/* this is just hardwired for the order in which we will put in the states of trios when we include missing data in them. I am just going to use a big switch statement using 
 the things I determined in:  http://users.soe.ucsc.edu/~eriq/dokuwiki/doku.php?id=project:pfr:enumerating_the_v_vecs_for_the_64_a_statesm
 */
int ReturnXstateOfA_64(int A)
{
	switch(A) {
		case(0):
		case(1):
		case(3):
		case(4):
		case(5):
		case(7):
		case(12):
		case(13):
		case(15):
		case(17):
		case(18):
		case(19):
		case(20):
		case(21):
		case(22):
		case(23):
		case(24):
		case(25):
		case(27):
		case(28):
		case(29):
		case(30):
		case(31):
		case(37):
		case(38):
		case(39):
		case(41):
		case(42):
		case(43):
		case(45):
		case(46):
		case(47):
		case(48):
		case(49):
		case(50):
		case(51):
		case(52):
		case(53):
		case(54):
		case(55):
		case(56):
		case(57):
		case(58):
		case(59):
		case(60):
		case(61):
		case(62):
		case(63):
			return(0);
			break;
		case(16):
		case(26):
			return(1);
			break;
		case(2):
		case(6):
		case(14):
		case(36):
		case(40):
		case(44):
			return(2);
			break;
		case(8):
		case(9):
		case(11):
		case(33):
		case(34):
		case(35):
			return(3);
			break;
		case(10):
		case(32):
			return(4);
			break;
		default:
			fprintf(stderr,"Error.  A64 state %d unrecognized in ReturnXstateOfA_64().  Exiting.\n",A);
			exit(1);
			break;
	}
}






/* this takes the Aprobs array (hardwired for the case of trios) and fills the
   Aprobs_M array in ret.  It assumes that missing data are denoted by 3 and assumes that
   memory is already allocated to ret.

 */
void ReturnAProbsM_HardWiredForTrios(RunPars *RP, double *AP, double *ret) 
{
	int l,y[3],ylo[3],yhi[3],yin[3];
	
	int a,am;
	
	/*double *ret = (double *)calloc(64,sizeof(double)); */
	
	/* cycle over all possible 64 states */
	for(y[0]=0;y[0]<4;y[0]++) {
		for(y[1]=0;y[1]<4;y[1]++) {
			for(y[2]=0;y[2]<4;y[2]++) {
				
				am = GiveAFromY_M(y,RP);  /* this is the index of the A state (with missing data) */
				
				/* now, determine which A-states must be summed over to account for the missing data */
				for(l=0;l<3;l++)  {
					if(y[l]==3) {
						ylo[l]=0;
						yhi[l]=2;
					}
					else {
						ylo[l]=y[l];
						yhi[l]=y[l];
					}
				}
				
				/* now cycle over all the A states that contribute to the AProbsM */
				for(yin[0]=ylo[0]; yin[0]<=yhi[0]; yin[0]++) {
					for(yin[1]=ylo[1]; yin[1]<=yhi[1]; yin[1]++) {
						for(yin[2]=ylo[2]; yin[2]<=yhi[2]; yin[2]++) {
							a = GiveAFromY(yin,RP);  /* this is the index of the A state with no missing data */
							
							/* collect a sum.  It is already initialized to 0.0 from calloc */
							ret[am] += AP[a];
						}
					}
				}
			}
		}
	}
}


/* compare two ParentPairMatch structs on the basis of their posterior prob of parentage.  Used for sorting.
 This is intended to sort ones with higher posterior probability first. */
int cmp_parent_pair(const void *a, const void *b)
{
	ParentPairMatch *ia = (ParentPairMatch *)a;
	ParentPairMatch *ib = (ParentPairMatch *)b;
	
	if(ia->TrioPosteriors[0] > ib->TrioPosteriors[0]) {
		return(-1);
	}
	else if(ia->TrioPosteriors[0] < ib->TrioPosteriors[0]) {
		return(1);
	}
	return(0);
}






/*
 This normalizes a lot of probabilities in the dp arrays of each fb_cell,
 and it normailzes the probs of all the S-states at the last locus and 
 stores those in FinalStepNormalizedProbs.
*/
void PrepForBackwardStep(fb_struct *FBS, RunPars *RP, double **X)
{
	int i,j,t;
	double normo;
	int L=RP->L,
		Nx=RP->Nx;
	double *d;
	
	/* allocate memory to FinalStepNormalizedProbs */
	FBS->FinalStepNormalizedProbs = (double *)calloc(FBS->Kt[L-1],sizeof(double));
	
	/* then store probs and normalize */
	normo = 0.0;
	for(i=0;i<FBS->Kt[L-1];i++)  {
		FBS->FinalStepNormalizedProbs[i] = FBS->FB[L-1][i].p;
		normo += FBS->FinalStepNormalizedProbs[i];
	}
	for(i=0;i<FBS->Kt[L-1];i++)  {
		FBS->FinalStepNormalizedProbs[i] /= normo;
	}
	
	/* here we print something out to check that normo = 1-Sup. */
	//fprintf(stderr,"PROGRESS_CHECK!  Sum of cells at last time step= %e (1-that= %e )  and 1-Sup  is %e   (Sup is %e)\n",normo,1.0-normo, 1.0-FBS->S_up,FBS->S_up);
	FBS->S_down = 1.0 - FBS->S_up;
	
	/* now we go and normalize all the dp probs */
	for(t=L-1;t>=1;t--)  {
		for(i=0;i<FBS->Kt[t];i++)  {
			d = FBS->FB[t][i].dp;
			normo = 0.0;
			for(j=0;j<Nx;j++)  {
				if(FBS->FB[t][i].down[j]==NULL)  {
					d[j] = 0.0;
				}
				else {
					d[j] = FBS->FB[t][i].down[j]->p * X[t][j];
					normo += d[j];
				}
			}
			/* now cycle again to normalize */
			for(j=0;j<Nx;j++)  {
				d[j] /= normo;
			}
		}
	}
	
}



/* Simulate a backward step.  Requires that memory is allocated to GT already.
 A are the AprobsGX's.  In the process, it will also compute the probability
 of this realization. */
double SimBackward(fb_struct *FBS, double ***A, RunPars *RP, trio_genotype *GT)
{
	int t;
	int L=RP->L;
	struct fb_cell **FB = FBS->FB;
	int *Kt = FBS->Kt;
	int Nx=RP->Nx;
	int *NaInX=RP->NaInX;
	int **aInX2A = RP->aInX2A;
	struct fb_cell *cur;
	double ret=1.0;
	int rando,rando2;
	
	
	/* draw the starting point, and record that in the probability */
	GT->anchor = IntFromProbsRV(FBS->FinalStepNormalizedProbs, 0, Kt[L-1]);
	cur = &(FB[L-1][GT->anchor]);
	ret *= FBS->FinalStepNormalizedProbs[ GT->anchor ];
	
	/* now take a path from anchor back through the FB struct */
	for(t=L-1;t>=0;t--)  {
		
		if(t>0) {  /* for all cases with t>0 */
			/* choose the X-state that takes you to the previous step and record the probability thereof */
			rando = IntFromProbsRV(cur->dp, 0, Nx);
			GT->path[t] = rando;
			ret *= cur->dp[rando];
			
			
			/* now choose the agx associated with that and multiply the probability of that in too */
			rando2 = IntFromProbsRV(A[t][rando], 0, NaInX[rando]);
			GT->agx[t] = rando2;
			ret *= A[t][rando][rando2];
			
			/* now set cur to the cell pointed from the appropriate element in down */
			cur = cur->down[rando];
		}
		else {  
			/* for the first locus, now that we are here, we know that path[0] is just the idx of cur. The
			 probability of that occurring should be accounted for by the dp from [1] that got us here,
			 so we shouldn't need to add it in here. */
			GT->path[t] = cur->idx;
			
			/* now choose an a */
			rando2 = IntFromProbsRV(A[t][cur->idx], 0, NaInX[cur->idx]);
			GT->agx[t] = rando2;
			ret *= A[t][cur->idx][rando2];
		}
		/* here, at the end we can set the GT->a variable */
		GT->a[t] = aInX2A[ GT->path[t] ][ GT->agx[t] ];
		
	}
	
	return(ret);
}


/* assuming that the forwards algorithm has been run on each collection, and that
   a trio genotype has been simulated into GT.  This function returns the probability that
   the trio genotype would have occurred in collection C, conditional on having 
   number of Mendelian incompatibilities as specified by smax 
 
 Condition is a flag.  If it is 0 then we don't condition on the fact that it was in S-down.
 Otherwise we do.
*/
double CondProbOfGeno(Collection *C, trio_genotype *GT, RunPars *RP, int Condition)
{
	int t;
	int L=RP->L;
	double **A=C->Aprobs;
	int *a=GT->a;
	double prod=1.0;
	
	/* get the product of genotype probs */
	for(t=0;t<L;t++) {
		prod *= A[t][a[t]];
	}
	
	if(Condition!=0) {
		/* divide by S_down */
		prod /= C->FBS->S_down;
	}
	return(prod);
}

/* here are some macros for checking that preamble input is in there already */
#define PREAM_DONE ALREADY_HAS(Nx_f,--Nx) && ALREADY_HAS(d_f,-d) && ALREADY_HAS(L_f,-L) && ALREADY_HAS(smax_f,--smax) \
	&& ALREADY_HAS(v_vecs_f,--v-vecs) && ALREADY_HAS(Na_f,--Na) && ALREADY_HAS(AinX_f,--AinX) \
	&& ALREADY_HAS(NaInX_f,--NaInX)
FB_Vars *Get_pbt_C_fb_Opts(int argc, char *argv[])
{
	FB_Vars *ret;
	int Nx_f=0,
		d_f=0,
		L_f=0,
		smax_f=0,
		v_vecs_f=0,
		Na_f=0,
		AinX_f=0,
		NaInX_f=0,
		NY_f=0,
		Ystates_f=0,
		num_miss_modes_f=0,
		miss_data_modes_f=0,
		miss_mode_vvec_cens_f=0,
		new_coll_f=0,
		new_loc_f=0,
	    Xprob_f=0,
		Xfakeprob_f=0,
		Aprobs_f=0,
		AprobsGX_f=0;
	int loc=-1; /* to count the number of locs read in for each collection. */
	int Xpr=0, Xfpr=0, Apr=0, AprGX=0; /* for counting the number of XProbs, Xfakeprobs, AProbs and AProbsGX collected for each locus */
	int temp[MAX_X_STATES];
	int i,j;

	DECLARE_ECA_OPT_VARS  
	
	/* This is a good place to set some default values for variables in your program */
	
	/* some information about the program, author, etc */
	SET_PROGRAM_NAME("pbt_C_fb");  /* Note:  QUOTED string */
	SET_PROGRAM_SHORT_DESCRIPTION("the forwards-bacwards algorithm for PBT implemented in C"); /* QUOTED string */
	SET_PROGRAM_LONG_DESCRIPTION(xxx\054 xxx);  /* UN-QUOTED string! */
	SET_PROGRAM_AUTHOR_STRING("Eric C. Anderson"); /* QUOTED string */
	SET_VERSION("Version XX");  /* QUOTED string */
	SET_VERSION_HISTORY("Version 00 started Dec. 5, 2009\n"); /*QUOTED string */
	
	
	
	/* some memory allocation */
	ret = (FB_Vars *)malloc(sizeof(FB_Vars));
	ret->RP = (RunPars *)malloc(sizeof(RunPars));
	ret->Colls = (Collection **)calloc(MAX_COLLECTIONS,sizeof(Collection *));
	
	/* some defaults and initializations  */
	ret->RP->Na = 0;
	ret->RP->NaInX = NULL;
	ret->RP->AinX = NULL;
	ret->NumColls = -1;
	ret->RP->HasFakeProbs = 0;
	ret->RP->NY=0;
	ret->RP->Ystates=NULL;
	
	BEGIN_OPT_LOOP  
	
	
	/* use this to start a block of related options.  Close the block with the CLOSE_SUBSET macro (below)  */
	OPEN_SUBSET(Preamble Options,  /* Name of Subset for output in --help, --help-full, and --help-nroff */
				Preamble Options, /* Name of subset for lumping them together in GuiLiner */
				Blah blah blah /* Description of subset.  Not really used currently */
	)   /* NOTE.  These are all UNQUOTED strings.  Beware of commas. */
	
	if(REQUIRED_OPTION(
					   Number of X states,
					   Nx_f,	
					   ,
					   Nx,
					   J,
					   Number of X states,
					   )) {
		if(ARGS_EQ(1)) {
			ret->RP->Nx = GET_INT;
			#ifdef VERBOSE_READING_FB_INPUT
				printf("READING_INPUT: Nx : %d\n",ret->RP->Nx);
			#endif
		}
	}
	
	if(REQUIRED_OPTION(
					   Length of v-vectors,
					   d_f,	
					   d,
					   ,
					   J,
					   Length of the v-vectors,
					   )) {
		if(ARGS_EQ(1)) {
			ret->RP->d = GET_INT;
			#ifdef VERBOSE_READING_FB_INPUT
				printf("READING_INPUT: d : %d\n",ret->RP->d);
			#endif
		}
	}
	
	if(REQUIRED_OPTION(
					   Number of Loci,
					   L_f,	
					   L,
					   ,
					   J,
					   Number of loci,
					   )) {
		if(ARGS_EQ(1)) {
			ret->RP->L = GET_INT;
			#ifdef VERBOSE_READING_FB_INPUT
				printf("READING_INPUT: L : %d\n",ret->RP->L);
			#endif
		}
	}
	
	if(REQUIRED_OPTION(
					   s_max vector,
					   smax_f,	
					   ,
					   smax,
					   J1 ... Jd,
					   s_max vector,
					   )) {
		if(ALREADY_HAS(d_f,-d) && ARGS_EQ(ret->RP->d)) { int i;
			ret->RP->smax = (int *)calloc(ret->RP->d,sizeof(int));
			#ifdef VERBOSE_READING_FB_INPUT
				printf("READING_INPUT: smax : ");
			#endif
			for(i=0;i<ret->RP->d;i++)  {
				ret->RP->smax[i] = GET_INT;
				#ifdef VERBOSE_READING_FB_INPUT
					printf(" %d", ret->RP->smax[i]);
				#endif
			}
			#ifdef VERBOSE_READING_FB_INPUT
				printf("\n");
			#endif
		}
	}
	
	if(REQUIRED_OPTION(
					   V vectors,
					   v_vecs_f,	
					   ,
					   v-vecs,
					   J1 ... J(Nx x d),
					   v-vectors all in series.,
					   Here you give the v-vectors in the order that you want to put them. If you have Nx X states and each 
					   v vector is of length d\054 then this options will expect Nx*d arguments.
					   )) {
		if(ALREADY_HAS(Nx_f,--Nx) && ALREADY_HAS(d_f,-d) && ARGS_EQ(ret->RP->Nx * ret->RP->d)) { int i; int j;
			ret->RP->v_vecs = (int **)calloc(ret->RP->Nx,sizeof(int *));
			for(i=0;i<ret->RP->Nx;i++)  {
				ret->RP->v_vecs[i] = (int *)calloc(ret->RP->d,sizeof(int));
				#ifdef VERBOSE_READING_FB_INPUT
					printf("READING_INPUT: v-vector %d : ",i);
				#endif
				for(j=0;j<ret->RP->d;j++) {
					ret->RP->v_vecs[i][j] = GET_INT;
					#ifdef VERBOSE_READING_FB_INPUT
						printf(" %d", ret->RP->v_vecs[i][j]);
					#endif
				}
				#ifdef VERBOSE_READING_FB_INPUT
					printf("\n");
				#endif
			}			
		}
	}
	
	
	if(OPTION(
			   Number of A states,
			   Na_f,	
			   ,
			   Na,
			   J1,
			   Number of A states total,
			   )) {
	
		if(ARGS_EQ(1)) {
			ret->RP->Na = GET_INT;
			#ifdef VERBOSE_READING_FB_INPUT
				printf("READING_INPUT: Na : %d\n",ret->RP->Na);
			#endif
		}
	}
	
	if(OPTION(
			  X state membership of A states,
			  AinX_f,	
			  ,
			  AinX,
			  J1 ... JNa,
			  The X-states to which each A state belongs,
			  )) {
		
		if(ALREADY_HAS(Na_f,--Na)  &&  ARGS_EQ(ret->RP->Na)) { int i;
			#ifdef VERBOSE_READING_FB_INPUT
				printf("READING_INPUT: AinX : ");
			#endif
			ret->RP->AinX = (int *)calloc(ret->RP->Na, sizeof(int));
			for(i=0;i<ret->RP->Na;i++) {
				ret->RP->AinX[i] = GET_INT;
				#ifdef VERBOSE_READING_FB_INPUT
					printf(" %d",ret->RP->AinX[i]);
				#endif
			}
			#ifdef VERBOSE_READING_FB_INPUT
				printf("\n");
			#endif
		}
	}
	
	
	if(OPTION(
			  Number of A states in each X state,
			  NaInX_f,	
			  ,
			  NaInX,
			  J1 ... JNx,
			  Number of A states in each X state.  Hey Eric\054 it should be possible to just get this
			  information directly from the AinX input and make this option unnecessary.,
			  )) {
		
		if(ALREADY_HAS(Na_f,--Na) && ALREADY_HAS(Nx_f,--Nx) && ARGS_EQ(ret->RP->Nx)) { int i; int sum;
			#ifdef VERBOSE_READING_FB_INPUT
				printf("READING_INPUT: NaInX : ");
			#endif
			ret->RP->NaInX = (int *)calloc(ret->RP->Nx, sizeof(int));
			sum=0;
			for(i=0;i<ret->RP->Nx;i++) {
				ret->RP->NaInX[i] = GET_INT;
				#ifdef VERBOSE_READING_FB_INPUT
					printf(" %d",ret->RP->NaInX[i]);
				#endif
				sum+=ret->RP->NaInX[i];
			}
			#ifdef VERBOSE_READING_FB_INPUT
				printf("\n");
			#endif
			if(sum!=ret->RP->Na) {
				fprintf(stderr,"Error! The sum over NaInX is %d but we expect it to be the same as Na = %d.  Exiting",sum,ret->RP->Na);
				exit(1);
			}
		}
	}
	
	if(OPTION(
			  Number of Y states underlying each A state,
			  NY_f,
			  ,
			  NY,
			  J,
			  Number of Y states underlying each A state,
			  Number of Y states underlying each A state
			  )) {
		if(ARGS_EQ(1)) {
			ret->RP->NY = GET_INT;
		}
	}
	if(OPTION(
			  Define Y states underlying each A state,
			  Ystates_f,
			  ,
			  Ystates,
			  J1 ... JNY*NA,
			  Define Y states underlying each A state,
			  Define Y states underlying each A state.  For each A state in order just list the underlying Y states.  There 
			  should be NY * NA of these elements.  Requires giving the Na and NY options first.
			  )) {
		if(ALREADY_HAS(Na_f,Na) && ALREADY_HAS(NY_f,NY) ) { int i; int j; int k;  int l;
			if(ARGS_EQ(ret->RP->Na * ret->RP->NY)) {
				ret->RP->Ystates = (int **)calloc(ret->RP->Na, sizeof(int *));
				ret->RP->ymaxes = (int *)calloc(ret->RP->NY,sizeof(int));
				for(i=0;i<ret->RP->Na;i++)  {
					ret->RP->Ystates[i] = (int *)calloc(ret->RP->NY,sizeof(int));
					for(j=0;j<ret->RP->NY;j++)  {
						ret->RP->Ystates[i][j] = GET_INT;
						if(ret->RP->Ystates[i][j] > ret->RP->ymaxes[j]) {
							ret->RP->ymaxes[j] = ret->RP->Ystates[i][j];
						}
					}
				}
				/* now that is done, we can make our little hash to get to a from y */
				ret->RP->AgivenY = (int *)calloc( ReturnVDex(ret->RP->NY, ret->RP->ymaxes, ret->RP->ymaxes)+1, sizeof(int));
				for(i=0;i<ret->RP->Na;i++)  {
					ret->RP->AgivenY[ ReturnVDex(ret->RP->NY, ret->RP->ymaxes, ret->RP->Ystates[i]) ] = i;
				}
				
				/* and now we are also going to set up all the things for dealing with missing data at these genotypes.  this is pretty much
				 just hard-wired into here, assuming that we are dealing with trios and that the missing data can occur at any member of the trios.
				 Recall that missing data is 3 while the observed genotypes are 0, 1, or 2. */
				ret->RP->ymaxes_M = (int *)calloc(ret->RP->NY,sizeof(int));
				for(i=0;i<ret->RP->NY;i++) { ret->RP->ymaxes_M[i] = ret->RP->ymaxes[i]+1; }
				ret->RP->Na_M = 64;
				ret->RP->AgivenY_M = (int *)calloc( ReturnVDex(ret->RP->NY, ret->RP->ymaxes_M, ret->RP->ymaxes_M)+1, sizeof(int));
				ret->RP->Ystates_M = (int **)calloc(ret->RP->Na_M, sizeof(int *));
				l=0;
				for(i=0;i<4;i++) { for(j=0;j<4;j++) { for(k=0;k<4;k++) {
					ret->RP->Ystates_M[l] = (int *)calloc(ret->RP->NY, sizeof(int));  /* this is pretty much hard-wired for the case that we have trios and SNP data */
					ret->RP->Ystates_M[l][0] = i;
					ret->RP->Ystates_M[l][1] = j;
					ret->RP->Ystates_M[l][2] = k;
					ret->RP->AgivenY_M[ ReturnVDex(ret->RP->NY, ret->RP->ymaxes_M, ret->RP->Ystates_M[l]) ] = l++;
				}}}
				
				
				
			}
		}
	}
	
	if(OPTION(
		Number of Missing Data Modes,
		num_miss_modes_f,
		,
		num-miss-modes,
		J,
		Number of possible patterns of missing data amongst the NY genotypes,
		Number of possible patterns of missing data amongst the NY genotypes.  For trios\054 NY=3 and
			  the number of missing data patterns/modes is 8.))
	   {
		   if(ARGS_EQ(1)) {
			   ret->RP->NumMissModes = GET_INT;
		   }
	   }
	   
	if(OPTION(
		Missing Data Modes,
		miss_data_modes_f,
		,
		miss-data-modes,
		J1 ... J(NumMissModes*NY),
		specify the modes/patterns of missing data possible,
		Use this option to specify the modes of missing data possible at a trio.  If used\054 it must be issued after the --num-miss-modes option and the --NY option.  You can 
			  think of having NumMissModes rows\054 each of length NY (the number of Y-states underlying each A-state).  Each row specifies which of the
			  Y-state genotypes is missing.  That is what J1 up to J(NumMissModes*NY) is.   This can only be issued after the --NY option\054 the 
			  --num-miss-modes option\054 and the --Ystates option have been issued. )) {
		if(ALREADY_HAS(NY_f,--NY) && ALREADY_HAS(num_miss_modes_f,--num-miss-modes)  && ALREADY_HAS(Ystates_f,--Ystates) ) {
			if(ARGS_EQ(ret->RP->NumMissModes * ret->RP->NY)) { int i; int j; int k; int l; int m; int n;
				ret->RP->MissModes = (int **)calloc(ret->RP->NumMissModes,sizeof(int *));
				for(i=0;i<ret->RP->NumMissModes;i++)  {
					ret->RP->MissModes[i] = (int *)calloc(ret->RP->NY,sizeof(int));
					for(j=0;j<ret->RP->NY;j++)  {
						ret->RP->MissModes[i][j] = GET_INT;
					}
				}
				
				/* right here we define an array that will get us back to the Miss Mode index given the pattern of missing data. */
				for(i=0;i<ret->RP->NumMissModes;i++)  {
					ret->RP->MissModeFrom0sAnd1s[ ret->RP->MissModes[i][0]  ]
												[ ret->RP->MissModes[i][1]  ]
												[ ret->RP->MissModes[i][2]  ] = i;
				}			
				
				/* and right here we want to compute the AtoA_64 vector, and I am just going to hard-wire it for trios  */
				ret->RP->AtoA_64 = (int **)calloc(ret->RP->NumMissModes,sizeof(int *));
				for(m=0;m<ret->RP->NumMissModes;m++)  { int y[3];
					ret->RP->AtoA_64[m] = (int *)calloc(ret->RP->Na,sizeof(int));
					l=0;
					for(i=0;i<3;i++) { for(j=0;j<3;j++) { for(k=0;k<3;k++) {
						y[0]=i; y[1]=j; y[2]=k;
						/*printf("CHECKING_MISS_STUFF: MissMode= %d   l= %d  y27= %d %d %d     MissPattern= ",m,l,y[0],y[1],y[2]); */
						for(n=0;n<3;n++) {
							if(ret->RP->MissModes[m][n]==1) {
								y[n]=3;
							}
							/*printf("%d ",ret->RP->MissModes[m][n]); */
						}
						/* printf("      y64= %d %d %d ",y[0],y[1],y[2]); */
						ret->RP->AtoA_64[m][l++] = ret->RP->AgivenY_M[ ReturnVDex(ret->RP->NY, ret->RP->ymaxes_M, y) ];
						/* printf("   AtoA64= %d\n",ret->RP->AtoA_64[m][l-1]); */
					}}}
				}
			}
		}
	}
	
	
	if(OPTION(
		Missing data censoring of X-states,
		miss_mode_vvec_cens_f,
		,
		miss-mode-vvec-cens,
		J1 ... J(NumMissModes*Nx),
		specify which X states are undetectable given each missing data mode/pattern,
		specify which X states are undetectable given each missing data mode/pattern.  There should be 
	    NumMissModes rows of integers each with Nx columns.  If the i\054j-th element is equal to 1 it means
		that having data with missing-data-mode i makes it impossible to observe mendelian incompatibility of the
		sort described by X-state j.  A 0 entry means the missing data does not interfere with being able to observe the 
	    Mendelian incompatibilities of the X-state. Requires prior issuance of the --Nx and --num-miss-modes options. )) {
			if(ALREADY_HAS(Nx_f,--Nx) && ALREADY_HAS(num_miss_modes_f,--num-miss-modes)) {
				if(ARGS_EQ(ret->RP->NumMissModes * ret->RP->Nx)) { int i; int j;
					ret->RP->MissModeCensoring = (int **)calloc(ret->RP->NumMissModes,sizeof(int *));
					for(i=0;i<ret->RP->NumMissModes;i++)  {
						ret->RP->MissModeCensoring[i] = (int *)calloc(ret->RP->Nx,sizeof(int));
						for(j=0;j<ret->RP->Nx;j++)  {
							ret->RP->MissModeCensoring[i][j] = GET_INT;
						}
					}
				}
			}
			
	}
	
	CLOSE_SUBSET;  /* done with the preamble options */
	
	
	
	
	OPEN_SUBSET(Collection Options,  /* Name of Subset for output in --help, --help-full, and --help-nroff */
				Collection Options, /* Name of subset for lumping them together in GuiLiner */
				Blah blah blah /* Description of subset.  Not really used currently */
				)   /* NOTE.  These are all UNQUOTED strings.  Beware of commas. */
	
	if(MULT_USE_REQ_OPTION(
						   Collection Name,
						   new_coll_f,
						   ,
						   coll-name,
						   S,
						   Begin a new collection and give it a name,
						   Issuing this option causes the program to prepare to read the probabilities associated with a new collection, 
						   MAX_COLLECTIONS)) {
		if(PREAM_DONE && ARGS_EQ(1)) { int i;
			
			ret->NumColls++;
			ret->Colls[ret->NumColls] = (Collection *)malloc(sizeof(Collection));
			ret->Colls[ret->NumColls]->CollName = (char *)calloc(MAX_COLLECTION_NAME_LENGTH,sizeof(char));
			GET_STR(ret->Colls[ret->NumColls]->CollName);
			/* here we do some error checking */
			if(ret->NumColls>0 && (loc+1)!=ret->RP->L) {
				fprintf(stderr,"Error! Trying to initialize collection %s but we only read %d loci in collection %s. Exiting...\n",loc+1,ret->Colls[ret->NumColls-1]->CollName);
				exit(1);
			}
			/* here we have to check to make sure that the last locus we read in in the previous collection had everything it needed to have in it */
			/* now we check that we got everything in the last locus */
			else if(loc>0 && Xpr+1 != ret->RP->Nx) {
				fprintf(stderr,"Error! Read incorrect number of --Xprob commands for locus %d. Read %d and expected %d. Exiting...\n",loc, Xpr+1, ret->RP->Nx);
				exit(1);
			}
			else if(loc>0 && (Xfpr+1 != ret->RP->Nx   &&  (Xfpr>-1 || Xfakeprob_f>0) )  ) {
				fprintf(stderr,"Error! Read incorrect number of --X-fake-prob commands for locus %d. Read %d and expected 0 or %d.  Recall that if you issue even just one X-fake-prob command for any X state you have to supply them to all.  Exiting...\n",loc, Xfpr+1, ret->RP->Nx);
				exit(1);
			}
			
			if(ret->RP->Na) {
				if(loc>0 && Apr != 1) {
					fprintf(stderr,"Error! Read incorrect number of --Aprob commands for locus %d. Read %d and expected %d. Exiting...\n",loc, Apr, 1);
					exit(1);
				}
				if(loc>0 && AprGX+1 != ret->RP->Nx) {
					fprintf(stderr,"Error! Read incorrect number of --AprobGX commands for locus %d. Read %d and expected %d. Exiting...\n",loc, AprGX+1, ret->RP->Nx);
					exit(1);
				}
			}
			
			loc=-1; /* reinitialize to keep track of the current locus */
			#ifdef VERBOSE_READING_FB_INPUT
				printf("READING_INPUT: NewCollection : %d %s\n",ret->NumColls,ret->Colls[ret->NumColls]->CollName);
			#endif
			
			/* now allocate memory to the elements of the collection */
			ret->Colls[ret->NumColls]->Xprobs = (double **)calloc(ret->RP->L,sizeof(double *));
			ret->Colls[ret->NumColls]->Xfakeprobs = (double **)calloc(ret->RP->L,sizeof(double *));
			if(ret->RP->Na) {  /* if we are also specifying the A's */
				ret->Colls[ret->NumColls]->Aprobs = (double **)calloc(ret->RP->L,sizeof(double *));
				ret->Colls[ret->NumColls]->Aprobs_M = (double **)calloc(ret->RP->L,sizeof(double *));
				ret->Colls[ret->NumColls]->AprobsGX = (double ***)calloc(ret->RP->L,sizeof(double **));
			}
			for(i=0;i<ret->RP->L;i++) {
				ret->Colls[ret->NumColls]->Xprobs[i] = (double *)calloc(ret->RP->Nx,sizeof(double));
				ret->Colls[ret->NumColls]->Xfakeprobs[i] = (double *)calloc(ret->RP->Nx,sizeof(double));
				if(ret->RP->Na) {  /* if we are also specifying the A's */
					ret->Colls[ret->NumColls]->Aprobs[i] = (double *)calloc(ret->RP->Na,sizeof(double));
					ret->Colls[ret->NumColls]->Aprobs_M[i] = (double *)calloc(64,sizeof(double));  /* this is just hard-wired for the trios case */
					/*ret->Colls[ret->NumColls]->Aprobs_M[i] = (double *)calloc(ret->RP->Na_M,sizeof(double));  Not necessary here, we do it with a function down below . */
					ret->Colls[ret->NumColls]->AprobsGX[i] = (double **)calloc(ret->RP->Nx,sizeof(double *));
					if(NaInX_f) { int j;
						for(j=0;j<ret->RP->Nx;j++) {
							ret->Colls[ret->NumColls]->AprobsGX[i][j] = (double *)calloc(ret->RP->NaInX[j],sizeof(double));
						}
					}
					else {
						fprintf(stderr,"Error!  Allocating memory for Collection %s.  Since you specifid the number of A states, you must also specify the NaInX option. Exiting.\n",
								ret->Colls[ret->NumColls]->CollName);
						exit(1);
					}
				}
			}
		}
	}
	
	
	if(MULT_USE_REQ_OPTION(
						  New Locus,
						  new_loc_f,
						  ,
						  locus,
						  ,
						  Signify the beginning of a new locus,
						  Issuing this option causes the program to prepare to read the probabilities associated with a new locus, 
						   MAX_LOCI*MAX_COLLECTIONS))  {
		if(PREAM_DONE && ALREADY_HAS(new_coll_f,--new-coll) && ARGS_EQ(0)) {
			loc++; 
			
			/* now we check that we got everything in the last locus */
			if(loc>0 && Xpr+1 != ret->RP->Nx) {
				fprintf(stderr,"Error! Read incorrect number of --Xprob commands for locus %d. Read %d and expected %d. Exiting...\n",loc-1, Xpr+1, ret->RP->Nx);
				exit(1);
			}
			if(loc>0 && (Xfpr+1 != ret->RP->Nx   &&  (Xfpr>-1 || Xfakeprob_f>0) )  ) {
				fprintf(stderr,"Error! Read incorrect number of --X-fake-prob commands for locus %d. Read %d and expected 0 or %d.  Recall that if you issue even just one X-fake-prob command for any X state you have to supply them to all.  Exiting...\n",loc, Xfpr+1, ret->RP->Nx);
				exit(1);
			}
			
			if(ret->RP->Na) {
				if(loc>0 && Apr != 1) {
					fprintf(stderr,"Error! Read incorrect number of --Aprob commands for locus %d. Read %d and expected %d. Exiting...\n",loc-1, Apr, 1);
					exit(1);
				}
				if(loc>0 && AprGX+1 != ret->RP->Nx) {
					fprintf(stderr,"Error! Read incorrect number of --AprobGX commands for locus %d. Read %d and expected %d. Exiting...\n",loc-1, AprGX+1, ret->RP->Nx);
					exit(1);
				}
			}
			
			/* and now we set out counters */
			Xpr = -1;
			Xfpr = -1;
			Apr = 0;
			AprGX = -1;
			
			#ifdef VERBOSE_READING_FB_INPUT
				printf("READING_INPUT: NewLocus : %d\n",loc);
			#endif
		}
	}
	
	if(MULT_USE_OPTION(
						   A Probs,
						   Aprobs_f,
						   ,
						   Aprobs,
						   R1 ... RNa,
						   Specify the probs of the A states at the current locus. ,
						   Specify the probs of the A states at the current locus.  For each locus you must specify the probs of all of the A states\054 
						     in the correct order\054 with a single applications of this option.  This is required 
						     before issuing the next locus command\054 if you gave the --Na option.  Ultimately I am going to have to change the dependence
						     here so that this is only necessary if you specified in the preamble that you were going to give Aprobs for this particular locus. This
							 can easily be done by specifying Na for every locus.  Later..., 
						   MAX_LOCI * MAX_A_STATES * MANY_MORE))  {
		if(PREAM_DONE && ALREADY_HAS(new_coll_f,--new-coll) && ARGS_EQ(ret->RP->Na)) { int i;
			for(i=0;i<ret->RP->Na;i++) {
				ret->Colls[ret->NumColls]->Aprobs[loc][i] = GET_DUB;
				#ifdef VERBOSE_READING_FB_INPUT
					printf("READING_INPUT: Aprobs  %d : %.16f\n",i,ret->Colls[ret->NumColls]->Aprobs[loc][i]);
				#endif
			}
			/* and now that we have those, we ought to crunch them into the Aprobs_M as well */
			ReturnAProbsM_HardWiredForTrios(ret->RP, ret->Colls[ret->NumColls]->Aprobs[loc],ret->Colls[ret->NumColls]->Aprobs_M[loc]);
			
			Apr++;;
		}
		
	}
	if(MULT_USE_REQ_OPTION(
						   X Prob,
						   Xprob_f,
						   ,
						   Xprob,
						   R,
						   Specify the prob on an X state. ,
						   Specify the prob on an X state.  For each locus you must specify all of the X states with separate applications of this option\054 in the correct order\054
						    before issuing the next locus command, 
						   MAX_LOCI * MAX_X_STATES * MANY_MORE))  {
		if(PREAM_DONE && ALREADY_HAS(new_coll_f,--new-coll) && ARGS_EQ(1)) {
			ret->Colls[ret->NumColls]->Xprobs[loc][++Xpr] = GET_DUB;
			#ifdef VERBOSE_READING_FB_INPUT
				printf("READING_INPUT: Xprob %d : %.16f\n",Xpr,ret->Colls[ret->NumColls]->Xprobs[loc][Xpr]);
			#endif
		}
		
	}
	
	if(MULT_USE_OPTION(
						   X Fake Prob,
						   Xfakeprob_f,
						   ,
						   X-fake-prob,
						   R,
						   Specify the prob of an X state when missing data are ignored ,
						   Specify the prob on an X state when missing data are given a probability of one (effectively they are ignored).
						   For each locus you must specify all of the fake X states with separate applications of this option\054 in the correct order\054
						   before issuing the next locus command.  This is not a required option.  But\054 if you give it for any X-state at any locus
						   you must provide it for every X-state in every locus., 
						   MAX_LOCI * MAX_X_STATES * MANY_MORE))  {
		if(PREAM_DONE && ALREADY_HAS(new_coll_f,--new-coll) && ARGS_EQ(1)) {
			ret->Colls[ret->NumColls]->Xfakeprobs[loc][++Xfpr] = GET_DUB;
			ret->RP->HasFakeProbs = 1; 
			#ifdef VERBOSE_READING_FB_INPUT
				printf("READING_INPUT: X-fake-prob %d : %.16f\n",Xfpr,ret->Colls[ret->NumColls]->Xfakeprobs[loc][Xfpr]);
			#endif
		}
		
	}
	
	if(MULT_USE_REQ_OPTION(
						   A Probs Given X,
						   AprobsGX_f,
						   ,
						   AprobsGX,
						   R1 ... R_NaInX,
						   Specify the probs of the NaInX A states at a locus ,
						   Specify the probs of the NaInX A states at a locus.  For each locus you must specify the probs of the A states given each X state with separate applications of this option
						     for each X state. Issue them in the correct order., 
						   MAX_LOCI * MAX_X_STATES * MANY_MORE))  {
		if(PREAM_DONE && ALREADY_HAS(new_coll_f,--new-coll) && ARGS_EQ(ret->RP->NaInX[++AprGX])) { int i;
			for(i=0;i<ret->RP->NaInX[AprGX];i++) {
				ret->Colls[ret->NumColls]->AprobsGX[loc][AprGX][i] = GET_DUB;
				#ifdef VERBOSE_READING_FB_INPUT
					printf("READING_INPUT: AprobsGX  %d  %d : %.16f\n",AprGX,i,ret->Colls[ret->NumColls]->AprobsGX[loc][AprGX][i]);
				#endif
			}
		}
		
	}	
	CLOSE_SUBSET;  /* done with the preamble options */
	
	
	END_OPT_LOOP   /* Macro that loses the Main ECA_Opt Loop and does some error checking 
	 behind the scenes */
	
	/* here at the very end we have to check that the correct number of options were issued for the last locus and the last collection */
	/* now we check that we got everything in the last locus */
	else if(loc>0 && Xpr+1 != ret->RP->Nx) {
		fprintf(stderr,"Error! Read incorrect number of --Xprob commands for locus %d. Read %d and expected %d. Exiting...\n",loc, Xpr+1, ret->RP->Nx);
		exit(1);
	}
	if(ret->RP->Na) {
		if(loc>0 && Apr != 1) {
			fprintf(stderr,"Error! Read incorrect number of --Aprob commands for locus %d. Read %d and expected %d. Exiting...\n",loc, Apr, 1);
			exit(1);
		}
		if(loc>0 && AprGX+1 != ret->RP->Nx) {
			fprintf(stderr,"Error! Read incorrect number of --AprobGX commands for locus %d. Read %d and expected %d. Exiting...\n",loc, AprGX+1, ret->RP->Nx);
			exit(1);
		}
	}
	/* check the last collection */
	ret->NumColls++;  /* this also sets the number of collections to the correct number */
	if(ret->NumColls>0 && (loc+1)!=ret->RP->L) {
		fprintf(stderr,"Error! Trying to initialize collection %s but we only read %d loci in collection %s. Exiting...\n",loc+1,ret->Colls[ret->NumColls-1]->CollName);
		exit(1);
	}
	#ifdef VERBOSE_READING_FB_INPUT
		printf("READING_INPUT: Number of Collections : %d\n",ret->NumColls);
		printf("READING_INPUT: Done with reading command line options with ECA_Opt\n");
	#endif
	
	/* now, down at the very end here we will allocate memory to and fill the aInX2A array.  This is a bit silly.  I should do this
	 up above and only have the AinX option and do everything from that one piece of info, but, oh well.   Print it out at the end. */
	for(i=0;i<MAX_X_STATES;i++) temp[i] = 0;
	ret->RP->aInX2A = (int **)calloc(ret->RP->Nx,sizeof(int *));
	for(j=0;j<ret->RP->Nx;j++)  {
		ret->RP->aInX2A[j] = (int *)calloc(ret->RP->NaInX[j],sizeof(int));
	}
	for(i=0;i<ret->RP->Na;i++) {
		ret->RP->aInX2A[ ret->RP->AinX[i] ][ temp[ret->RP->AinX[i]]++ ] = i;
	}
	for(j=0;j<ret->RP->Nx;j++)  {
		for(i=0;i<ret->RP->NaInX[j];i++) {
			; //printf("FILLING_aInX2A:  %d  %d  :  %d\n",j,i,ret->RP->aInX2A[j][i]);
		}
	}

	return(ret);
}


/* print out a detailed summary of the collection after running the forwards algorithm AND preparing for the backward algorithm */
void PrintDetailedCollectionSummary(Collection *C, RunPars *RP)
{
	int i,j;

	printf("COLLECTION_SUMMARY: %s  Sup and Sdown %e  %e \n",C->CollName, C->FBS->S_up, C->FBS->S_down);
	for(i=0;i<C->FBS->Kt[RP->L-1];i++)  {
		printf("COLLECTION_SUMMARY %s    s-val (%d) :  ",C->CollName,i);
		for(j=0;j<RP->d;j++) {
			printf("%d ",C->FBS->FB[RP->L-1][i].s[j]);
		} 
		printf("  %e   %e    %e\n",C->FBS->FB[RP->L-1][i].p,  C->FBS->FB[RP->L-1][i].fake_p, C->FBS->FinalStepNormalizedProbs[i]);
	}
	
}
							   
