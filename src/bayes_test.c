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
#include "ECA_Opt3.h"
#include "MCTypesEtc.h"

/* this is going to be a super quick program.
 
 From stdin it will read a single set of 27 doubles that are the joint probs of trios under the 
 parental hypothesis.
 
 Then from those it will compute the marginal prob for any individual's genotype at that locus, and it will
 compute the conditional probability of an offspring given its parents.
 
 Then it will create a bunch of individuals, some fraction of which, at random will have children according
 to a poisson process, and then we will also make a bunch of unrelated "offspring," and we will see how well
 we can estimate the fraction of "offspring" that really have parents in the initial set.
 
 */

// a structure to hold the parameters for the program
typedef struct {
	int L; /* number of loci */
	int N; /* number of potential parent pairs */
	double PoisMu; /* the mean of the Poisson distributed number of offspring each potential parent will have */
	double SampFrac;  /* the fraction of our offspring pool that will have parents in the parent pair data base */
	int InputGenos; /* if 1, then it expects genotypes to be fed in on stdin */
	int D;  /* number of potential moms coming in on stdin */
	int S;  /* number of potential dads coming in on stdin */
	int K; /* number of kids coming in on stdin */
	
	int **y_ma;
	int **y_pa;
	int **y_kid;
	
	
	
	double *CTP;  /* conditional trio probs 1,...,27 */
	double *MP; /* Marginal probs 1,...,3 */
	double *SPC;  /* single parent conditional probs 1,...,9 */
	
	double ***StoredTP;  /* store the trio probs for every kid, pa, ma.  Indexed by [kid][pa][ma]  Obviously means we can't have really large data sets. */
 	double **StoredMaP; /* to store the pair probs with moms.  Indexed by [kid][ma] */
	double **StoredPaP; /* ditto for dads */
	
	double *SumTP;  /* sum over pairs of parents of the conditional prob of each kid.  indexed by [k]. */
	double *SumPaP;
	double *SumMaP;

	
	double **pi; /* a 2 by 2 array for pi00, pi01, pi10, and pi11 */  
	
} bt_vars;



/* here is a macro for the index of the genotype of a trio */
#define GENO(kid,pa,ma)  (kid*9 + pa*3 + ma)


/* here is a macro the index of the genotype of a parent-offspring pair */
#define SGENO(kid,pa) (kid*3 + pa)

/* some globals */
bt_vars *gBT;


/* for kid a vector of L ints that are 0, 1, or 2 (the number of 1 alleles)  and for marginal probabilities
 of each genotype given in marg_prob, thiks returns the probability of Kid's genotype just occurring at
 random in the population */
double MargProb(int *Kid, double *marg_prob, int L)
{
	int i;
	double ret=1.0;
	for(i=0;i<L;i++)  {
		ret *= marg_prob[Kid[i]];
	}
	return(ret);
}


/* for kid a vector of L ints that are 0, 1, or 2 (the number of 1 alleles)  and for conditional probabilities
 of each genotype given in CondProb, this returns the probability of Kid's genotype conditional on Pa and Ma. */
double CondProb(int *Kid, int *Pa, int *Ma, int L, double *CondProbs) 
{
	int i;
	double ret=1.0;
	
	for(i=0;i<L;i++)  {
		ret *= CondProbs[ GENO(Kid[i],Pa[i],Ma[i]) ];
	}
	return(ret);
}


/* for kid a vector of L ints 0, 1, or 2, and for the 9 single parent conditional probs in SPC,
 this returns the conditional prob of the kids genotype given its one parent Par */
double SP_CondProb(int *Kid, int *Par, int L, double *SCP) 
{
	int i;
	double ret=1.0;
	
	for(i=0;i<L;i++) {
		ret *= SCP [ SGENO(Kid[i],Par[i]) ];
	}
	return(ret);
	
}



bt_vars *GetBTOpts(int argc, char *argv[])
{
	
	int Lf=0,
	Nf=0,
	PoisMuf=0,
	SampFracf=0,
	InputGenosf=0;
	bt_vars *ret;
	
	DECLARE_ECA_OPT_VARS  
	
	
	
	/* This is a good place to set some default values for variables in your program */
	ret = (bt_vars *)malloc(sizeof(bt_vars));
	ret->InputGenos = 0;
	
	
	/* some information about the program, author, etc */
	SET_PROGRAM_NAME("bayes_test");  /* Note:  QUOTED string */
	SET_PROGRAM_SHORT_DESCRIPTION("yyy"); /* QUOTED string */
	SET_PROGRAM_LONG_DESCRIPTION(xxx\054 xxx);  /* UN-QUOTED string! */
	SET_PROGRAM_AUTHOR_STRING("Eric C. Anderson"); /* QUOTED string */
	SET_VERSION("Version XX");  /* QUOTED string */
	SET_VERSION_HISTORY("Version XX written Nov. 8, 2007\nGrew out of simpler programs.\nblah, blah..."); /*QUOTED string */
	
	
	BEGIN_OPT_LOOP  
	
	
	/* use this to start a block of related options.  Close the block with the CLOSE_SUBSET macro (below)  */
	OPEN_SUBSET(SubsetName1,  /* Name of Subset for output in --help, --help-full, and --help-nroff */
				SubsetName2, /* Name of subset for lumping them together in GuiLiner */
				Blah blah blah /* Description of subset.  Not really used currently */
				)   /* NOTE.  These are all UNQUOTED strings.  Beware of commas. */
	
	if(OPTION(
			  Take data on stdin,
			  InputGenosf,
			  s,
			  ,
			  J1 J2 J3,
			  Take input of moms and dads and kids from stdin,
			  Use this option to let the program know you want to feed it a bunch of genotype data.
			  J1 is the number of mas coming and J2 is the number of pas and J3 is the number of kids.
			  They come on stdin after the 27 doubles giving the trio probs for the allele frequency (note 
																									  that we assume the freq is the same for all loci---because this is primarily for a simple simulation).
			  The format for the individuals is quite simple: NO indiv identifier just a single string of L
	          0s 1s or 2s each giving the number of 1 alleles at the locus.   Note that you feed them in as moms first then
			  dads and then kids.)) {
		if(ARGS_EQ(3)) {
			ret->D = GET_INT;
			ret->S = GET_INT;
			ret->K = GET_INT;
			ret->InputGenos=1;
		}
		
	}
	
	CLOSE_SUBSET;  /* done with the greeting-related options */
	
	END_OPT_LOOP   /* Macro that loses the Main ECA_Opt Loop and does some error checking 
	 behind the scenes */
	
	return(ret);
}


/* once the marg and cond probs have been read, this is a simple function to read in the 
 genos that are input on the stdin */
void ReadGenosFromStdin(bt_vars *BT)
{
	int i,j;
	char c;
	
	char tempstr[5000]; /* for holding the strings temporarily */
	
	
	BT->y_ma=(int **)calloc(BT->D,sizeof(int *));
	BT->y_pa=(int **)calloc(BT->S,sizeof(int *));
	BT->y_kid=(int **)calloc(BT->K,sizeof(int *));
	
	
	fscanf(stdin," %s",tempstr);
	BT->L=strlen(tempstr);
	
	/* transer genotype to ma[0] */
	printf("ma 0  ");
	BT->y_ma[0] = (int *)calloc(BT->L,sizeof(int));
	for(j=0;j<BT->L;j++) {
		BT->y_ma[0][j] = tempstr[j] - '0';
		printf("%d",BT->y_ma[0][j]);
	}
	printf("\n");
	
	for(i=1;i<BT->D;i++) {
		BT->y_ma[i] = (int *)calloc(BT->L,sizeof(int));
		fscanf(stdin," %s",tempstr);
		printf("ma %d  ",i);
		for(j=0;j<BT->L;j++) {
			BT->y_ma[i][j] = tempstr[j] - '0';
			printf("%d",BT->y_ma[i][j]);
		}
		printf("\n");
	}
	for(i=0;i<BT->S;i++) {
		BT->y_pa[i] = (int *)calloc(BT->L,sizeof(int));
		fscanf(stdin," %s",tempstr);
		printf("pa %d  ",i);
		for(j=0;j<BT->L;j++) {
			BT->y_pa[i][j] = tempstr[j] - '0';
			printf("%d",BT->y_pa[i][j]);
		}
		printf("\n");	
	}
	for(i=0;i<BT->K;i++) {
		BT->y_kid[i] = (int *)calloc(BT->L,sizeof(int));
		fscanf(stdin," %s",tempstr);
		printf("kid %d  ",i);
		for(j=0;j<BT->L;j++) {
			BT->y_kid[i][j] = tempstr[j] - '0';
			printf("%d",BT->y_kid[i][j]);
		}
		printf("\n");	
	}
}


/* After reading all the input genos, compute all the trio probs and all the pair probs that can be formed by the sample */
void StoreProbs(bt_vars *BT)
{
	int k,s,d;
	
	BT->StoredTP = (double ***)calloc(BT->K,sizeof(double **));
	BT->StoredMaP = (double **)calloc(BT->K,sizeof(double *));
	BT->StoredPaP = (double **)calloc(BT->K,sizeof(double *));
	
	for(k=0;k<BT->K;k++)  {
		BT->StoredTP[k] = (double **)calloc(BT->S,sizeof(double *));
		BT->StoredPaP[k] = (double *)calloc(BT->S,sizeof(double ));
		for(s=0;s<BT->S;s++) {
			BT->StoredTP[k][s] = (double *)calloc(BT->D,sizeof(double));
			BT->StoredPaP[k][s] = SP_CondProb(BT->y_kid[k], BT->y_pa[s], BT->L, BT->SPC);
			//printf("PAIR_PROB_PA: %d %d  %e\n",k,s,BT->StoredPaP[k][s]);
			for(d=0;d<BT->D;d++) {
				BT->StoredTP[k][s][d] = CondProb(BT->y_kid[k], BT->y_pa[s], BT->y_ma[d], BT->L, BT->CTP);
				//printf("TRIO_PROB_PA_MA: %d %d %d  %e\n",k,s,d,BT->StoredTP[k][s][d]);
			}
		}
		/* and after all that we do the ma pair probs */
		BT->StoredMaP[k] = (double *)calloc(BT->D,sizeof(double ));
		for(d=0;d<BT->D;d++) {
			BT->StoredMaP[k][d] = SP_CondProb(BT->y_kid[k], BT->y_ma[d], BT->L, BT->SPC);
			//printf("PAIR_PROB_MA: %d %d  %e\n",k,d,BT->StoredMaP[k][d]);
		}
	}
}


/* store the sums (one for each kid, summed over individual parents or parent pairs)
 of the conditional probabilities */
void StoreSums(bt_vars *BT)
{
	int k,s,d;
	
	//printf("About to Store Sums\n");
	
	/* allocate the memory */
	BT->SumTP = (double *)calloc(BT->K,sizeof(double));
	BT->SumPaP = (double *)calloc(BT->K,sizeof(double));
	BT->SumMaP = (double *)calloc(BT->K,sizeof(double));
	
	for(k=0;k<BT->K;k++) {
		printf("k=%d\n",k);
		/* first sum the pairs */
		for(s=0;s<BT->S;s++) {
			for(d=0;d<BT->D;d++) {
				BT->SumTP[k] += BT->StoredTP[k][s][d];
			}
		}
		
		//printf("Done with pairs\n");
		/* sum the the pa's */
		for(s=0;s<BT->S;s++) {
			BT->SumPaP[k] += BT->StoredPaP[k][s];
		}
		
		//printf("Done with pas\n");
		
		/* sum the ma's */
		for(d=0;d<BT->D;d++) {
			BT->SumMaP[k] += BT->StoredMaP[k][d];
		}		
	}
	//printf("Done Storing Sums\n");
}

/* estimate all four elements of Pi by an EM algorithm */
void EstimatePiByEM(bt_vars *BT)
{
	int i,j,k,s,d;
	
	/* p is temp for holding values for a trio; q is for the sum over all trios for a kid; and t is for the sum over all the kids */
	double p[2][2], q[2][2], t[2][2];
	
	double dK=(double)BT->K;
	double dS=(double)BT->S;
	double dD=(double)BT->D;
	double **pi, **new_pi;
	double ***STP=BT->StoredTP;
	double **SP=BT->StoredPaP;
	double **SM=BT->StoredMaP;
	double normo,summed;
	double both, dad, mom, diff=1;;  /* prob that both parents are in the dataset, or that dad is or that mom is (not exclusive) */
	
	/* allocate to pi */
	BT->pi = (double **)calloc(2,sizeof(double *));
	new_pi = (double **)calloc(2,sizeof(double *));
	for(i=0;i<2;i++) {
		BT->pi[i] = (double *)calloc(2,sizeof(double));
		new_pi[i] = (double *)calloc(2,sizeof(double));
	}
	/* then initialize to .25's... That seems as reasonable a starting point as any. */
	for(i=0;i<2;i++) {
		new_pi[i][0] = .25;
		new_pi[i][1] = .25;
	}
	pi = BT->pi;
	
	
	while(diff>.0001) {
		
		/* transfer new values to pi */
		for(i=0;i<2;i++) for(j=0;j<2;j++) pi[i][j] = new_pi[i][j];
		
		/* initialize some values */
		t[0][0]=t[1][0]=t[0][1]=t[1][1] = 0.0;
		
		/* then cycle over kiddos and pairs.  */
		for(k=0;k<BT->K;k++) {
			q[0][0]=q[1][0]=q[0][1]=q[1][1] = 0.0;
			summed=0.0;
			for(s=0;s<BT->S;s++)  {
				for(d=0;d<BT->D;d++)  {
					p[0][0] = pi[0][0] * MargProb(BT->y_kid[k], BT->MP, BT->L);
					p[1][0] = (pi[1][0]/dS) * SP[k][s];
					p[0][1] = (pi[0][1]/dD) * SM[k][d];
					p[1][1] = (pi[1][1] / (dS*dD)) * STP[k][s][d];
					
					normo = p[0][0] + p[1][0] + p[0][1] + p[1][1]; 
					
					for(i=0;i<2;i++) for(j=0;j<2;j++) {q[i][j] += p[i][j] / normo;}
					summed += 1.0;
					
				}
			}
			for(i=0;i<2;i++) for(j=0;j<2;j++) t[i][j] += q[i][j] / summed;
		}
		
		both = t[1][1] * dS * dD / dK;
		dad = (t[1][0]+t[1][1]) * dS /dK;
		mom = (t[0][1]+t[1][1]) * dD /dK;
		
		new_pi[1][1] = both;
		new_pi[1][0] = ECA_MAX(0, dad-both);
		new_pi[0][1] = ECA_MAX(0, mom-both);
		new_pi[0][0] = 1.0 - new_pi[1][1] - new_pi[1][0] - new_pi[0][1];
		
		diff = ECA_ABS(new_pi[0][0]-pi[0][0]) + ECA_ABS(new_pi[1][1]-pi[1][1]) + ECA_ABS(new_pi[1][0]-pi[1][0]) + ECA_ABS(new_pi[0][1]-pi[0][1]);
		
		printf("EM Progress: 11= %f    10= %f      01= %f      00= %f    diff= %f\n",new_pi[1][1], new_pi[1][0], new_pi[0][1], new_pi[0][0], diff);
	}
	
	//printf("Estimate of fraction of kids with both parents: %e\n",t[1][1] * dS * dD / dK);
	//printf("Estimate of fraction with dad in the data: %e\n", (t[1][0]+t[1][1]) * dS /dK);
	//printf("Estimate of fraction with mom in the data: %e\n", (t[0][1]+t[1][1]) * dD /dK);
	
	
}


/* function to return the log-likelihood of pi.  Note that the 
 elements pi11, pi10, pi01, will have to be put into elements 1,2,3 
 of the vector q for this to work with powell's method, the way it 
 is coded with base 1 subscripting by the Numerical Recipes folks. */
float PiNegLogLikelihood(float q[])
{
	int k,s,d;
	double ret=0.0;
	double p00,p01,p10,p11;
	bt_vars *BT = gBT;  /* put a local pointer on the stack */
	
	double dS=(double)BT->S;
	double dD=(double)BT->D;
	double ***STP=BT->StoredTP;
	double **SP=BT->StoredPaP;
	double **SM=BT->StoredMaP;
	double margprob;
	
	/* first, return a very large value if any of the q's are out of range */
	if( q[1]<0.00001 || q[1] >.999999 ||
	   q[2]<0.00001 || q[2] >.999999 ||
	   q[3]<0.00001 || q[3] >.999999 ||
	   q[1] + q[2] + q[3] > .99999 ) {
		return(999.999e26);
	}
	
	p11 = (double)q[1];
	p10 = (double)q[2];
	p01 = (double)q[3];
	
	p00 = 1.0 - p11 - p01 - p10;
	
	for(k=0;k<BT->K;k++)  {
		ret += log(p00 * MargProb(BT->y_kid[k], BT->MP, BT->L) +
				   (p10/dS) * BT->SumPaP[k] + 
				   (p01/dD) * BT->SumMaP[k] +
				   (p11/(dS*dD)) * BT->SumTP[k]);
				    
	}
	
	return((float)(-ret));			
	
	
}

/* stuff for testing Powell's method from NR */
float testN(float x[])
{
	return( exp( (x[1]-2.0)*(x[1]-2.0)/2.0 ) * exp( (x[2]+0.5)*(x[2]+0.5)/2.0 )  );
}


#define NRANSI
#include "nrutil.h"
#define ITMAX 200

void powell(float p[], float **xi, int n, float ftol, int *iter, float *fret, float (*func)(float []))
{
	void linmin(float p[], float xi[], int n, float *fret,
				float (*func)(float []));
	int i,ibig,j;
	float del,fp,fptt,t,*pt,*ptt,*xit;
	
	
	printf("past the declarations\n");
	
	pt=vector(1,n);
	ptt=vector(1,n);
	xit=vector(1,n);
	*fret=(*func)(p);
	
	printf("past the setup\n");
	
	for (j=1;j<=n;j++) pt[j]=p[j];
	for (*iter=1;;++(*iter)) {
		fp=(*fret);
		ibig=0;
		del=0.0;
		for (i=1;i<=n;i++) {
			printf("about to start xiting\n");
			for (j=1;j<=n;j++) xit[j]=xi[j][i];
			fptt=(*fret);
			printf("about to call linmin with p = %f %f %f  and fret = %f\n",p[1],p[2],p[3],*fret);
			
			linmin(p,xit,n,fret,func);
			if (fabs(fptt-(*fret)) > del) {
				del=fabs(fptt-(*fret));
				ibig=i;
			}
		}
		if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) {
			free_vector(xit,1,n);
			free_vector(ptt,1,n);
			free_vector(pt,1,n);
			return;
		}
		if (*iter == ITMAX) nrerror("powell exceeding maximum iterations.");
		for (j=1;j<=n;j++) {
			ptt[j]=2.0*p[j]-pt[j];
			xit[j]=p[j]-pt[j];
			pt[j]=p[j];
		}
		fptt=(*func)(ptt);
		if (fptt < fp) {
			t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
			if (t < 0.0) {
				linmin(p,xit,n,fret,func);
				for (j=1;j<=n;j++) {
					xi[j][ibig]=xi[j][n];
					xi[j][n]=xit[j];
				}
			}
		}
	}
}
#undef ITMAX
#undef NRANSI






int main(int argc, char *argv[])
{
	int i,j,k,ma,pa,kid,n;
	double par_trio_probs[27], cond_trio_probs[27];
	double marg_prob[3],single_par_joint_prob[9], single_par_cond_prob[9];
	int N=100;  /* number of potential parent pairs */
	int L=40; /* number of loci */
	int ***Pairs;  /* indexed by i=1,..,N and then by 0 for pa and 1 for ma, then by l=1,..,L */
	float PoisMu=.4;  /* the mean of the Poisson distributed number of offspring each potential parent will have */
	int *NumOffs; /* number of offspring each potential parent has */
	int TotNumOffs=0;
	double SampFrac=.1;  /* the fraction of our offspring pool that will have parents in the parent pair data base */
	int TotSamp;  /* total number sampled in the pool of possible offspring */
	int **SamKids;
	double *CondProbSum;
	double pi,maxpi,maxlogl;
	double sum;
	bt_vars *BTO;
	
	
	/* just testing powell's method */
	float p[3];
	float **xi;
	int nn=2, iter;
	float ftol=.0000001, fret;
	
	p[1] = -1.3;
	p[2] = 0.0;
	
	xi = (float **)calloc(4,sizeof(float *));
	for(i=0;i<4;i++) xi[i] = (float *)calloc(4,sizeof(float));
	xi[1][1]=1.0;
	xi[2][2]=1.0;
	xi[3][3]=1.0;
	
	
	printf("Going for it!\n");
	
	powell(p, xi, nn, ftol,  &iter, &fret, testN);
	
	printf("%f  %f \n",p[1],p[2]);
	//return(0);
	
	
	
	BTO = GetBTOpts(argc, argv);
	
	
	SeedFromFile("bayes_test_seeds");
	
	/* read in the probs */
	for(i=0;i<27;i++)  {
		scanf(" %lf",&par_trio_probs[i]);
		printf("%f\n",par_trio_probs[i]);
	}
	
	/* verify that we have our indexing right */
	printf("\n\n");
	for(i=0;i<3;i++) for(j=0;j<3;j++) for(k=0;k<3;k++) printf("JOINT:  %d %d %d   %f\n",i,j,k, par_trio_probs[GENO(i,j,k)]);
	
	
	
	/* get the marginal prob of the 3 possible genotypes. Do it by holding ma constant and summing over all the others.*/
	for(ma=0;ma<3;ma++) {
		marg_prob[ma]=0.0;
		for(kid=0;kid<3;kid++) for(pa=0;pa<3;pa++) marg_prob[ma]+=par_trio_probs[ GENO(kid,pa,ma) ];
	}
	for(i=0;i<3;i++) {
		printf("MARG:  %d  %f\n",i,marg_prob[i]);
	}
	
	
	
	/* get the joint prob of kid and pa.  This is done by summing over ma */
	for(kid=0;kid<3;kid++) {
		for(pa=0;pa<3;pa++) {
			single_par_joint_prob[ SGENO(kid,pa) ] = 0.0;
			for(ma=0;ma<3;ma++) {
				single_par_joint_prob[ SGENO(kid,pa) ] += par_trio_probs[ GENO(kid,pa,ma) ];
			}
		}
	}
	/* verify and at the same time get the conditional pair prob by dividing out prob of pa */
	for(kid=0;kid<3;kid++) { 
		for(pa=0;pa<3;pa++) { 
			printf("PAIR_JOINT: %d  %d  %f\n",kid,pa,single_par_joint_prob[ SGENO(kid,pa) ]);
			single_par_cond_prob[ SGENO(kid,pa) ] = single_par_joint_prob[ SGENO(kid,pa) ] / marg_prob[pa];
		}
	}
	
	/* finally, print to verify the single_par_cond_prob */
	for(kid=0;kid<3;kid++) { 
		for(pa=0;pa<3;pa++) { 
			printf("PAIR_COND: %d  %d  %f\n",kid,pa,single_par_cond_prob[ SGENO(kid,pa) ]);
		}
	}
	
	
	/* now compute the vector of conditional probabilities.  Do this by dividing each element of the joint prob vector by 
	 the product of the prob of the two parents */
	for(ma=0;ma<3;ma++) for(kid=0;kid<3;kid++) for(pa=0;pa<3;pa++)  {
		cond_trio_probs[ GENO(kid,pa,ma) ] =  par_trio_probs[ GENO(kid,pa,ma) ] / (marg_prob[pa] * marg_prob[ma]);
		printf("COND: %d %d %d  %f\n",kid,pa,ma,cond_trio_probs[ GENO(kid,pa,ma) ]);
	}
	
	
	if(BTO->InputGenos==1) { float q[] = {0.0, .25, .25, .25 };
		/* transfer some data */
		BTO->CTP = cond_trio_probs;
		BTO->MP = marg_prob;
		BTO->SPC = single_par_cond_prob;
		gBT = BTO;  /* make a global pointer to it */
		
		ReadGenosFromStdin(BTO);
		StoreProbs(BTO);
		StoreSums(BTO);
		
		printf("%f %f %f    NegLogLike: %f\n",q[1],q[2],q[3], PiNegLogLikelihood( q));
		powell(q, xi, 3, ftol,  &iter, &fret, PiNegLogLikelihood);
		printf("Min at:  %f %f %f    NegLogLike: %f\n",q[1],q[2],q[3], fret);
		return(0);
		//EstimatePiByEM(BTO);
	}
	else {  /* here we do just the simulation of things */
		
		/* now we are ready to simulate some individual parent pairs */
		Pairs = (int ***)calloc(N,sizeof(int **));
		NumOffs = (int *)calloc(N,sizeof(int));
		for(i=0;i<N;i++) {
			Pairs[i] = (int **)calloc(2,sizeof(int *));
			for(k=0;k<2;k++) {
				Pairs[i][k] = (int *)calloc(L,sizeof(int));
				/* give them some genotypes */
				for(j=0;j<L;j++)  {
					Pairs[i][k][j] = IntFromProbsRV(marg_prob,0,3);
				}		
			}
			/* and simultate how many offspring this pair will have */
			NumOffs[i] = ignpoi(PoisMu);
			TotNumOffs+=NumOffs[i];
		}
		printf("TotNumOffs: %d\nNumOffs:",TotNumOffs);
		for(i=0;i<N;i++) printf(" %d",NumOffs[i]);
		printf("\n");
		
		
		
		/* now simulate the offspring for these */
		TotSamp = (int)((double)TotNumOffs/SampFrac);
		SamKids = (int **)calloc(TotSamp,sizeof(int *));
		n=0;
		for(i=0;i<N;i++)  {
			for(j=0;j<NumOffs[i];j++)  {
				SamKids[n] = (int *)calloc(L,sizeof(int));
				for(k=0;k<L;k++)  { double sum; double rando;
					pa=Pairs[i][0][k];
					ma=Pairs[i][1][k];
					/* now sample the offspring genotype */
					sum=0.0;
					rando=ranf();
					for(kid=0;kid<3;kid++)  {
						sum+=cond_trio_probs[ GENO(kid,pa,ma) ];
						if(sum>rando) break;
					}
					if(kid==3) { /* in this case the cond prob summed less than one perhaps */
						kid=2;
					}
					SamKids[n][k] = kid;
					printf("KIDPAMA n= %d    %d %d %d\n",n,kid,pa,ma);
					
				}
				n++;
			}
		}
		/* now put the rest in there */
		for(;n<TotSamp;n++)  {
			SamKids[n] = (int *)calloc(L,sizeof(int));
			for(j=0;j<L;j++)  {
				SamKids[n][j] = IntFromProbsRV(marg_prob,0,3);
				printf("MARG_KID: n= %d    %d\n",n,SamKids[n][j]);
			}	
		}
		
		printf("TOTSAMP: %d\n",TotSamp);
		
		
		
		
		/* now we just need to compute the sums for each */
		CondProbSum = (double *)calloc(TotSamp,sizeof(double));
		for(i=0;i<TotSamp;i++) {
			for(j=0;j<N;j++)  {
				CondProbSum[i] += CondProb(SamKids[i], Pairs[j][0], Pairs[j][1],  L, cond_trio_probs);
			}
		}
		
		
		/* now that we have that, we just have to compute the likelihood for pi */
		maxlogl=-99999999999999999.9;
		maxpi=0.0;
		for(pi=.01;pi<=.99;pi+=.01) {
			sum = 0.0;
			for(i=0;i<TotSamp;i++)  {
				sum += log(CondProbSum[i] * pi / (double)N + (1.0-pi) * MargProb(SamKids[i], marg_prob, L));
			}
			printf("PI:   %f     %f\n",pi,sum);
			if(sum>maxlogl) {
				maxpi = pi;
				maxlogl = sum;
			}
		}
		
		printf("MLE_PI:  %f   %f\n",maxpi,maxlogl);
		
		SeedToFile("bayes_test_seeds");
	}
	return(0);
	
	
	
}

























#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#define NR_END 1
#define FREE_ARG char*

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;
	
	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;
	
	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

unsigned char *cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	unsigned char *v;
	
	v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) nrerror("allocation failure in cvector()");
	return v-nl+NR_END;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;
	
	v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;
	
	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;
	
	/* allocate pointers to rows */
	m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;
	
	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	/* return pointer to array of pointers to rows */
	return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;
	
	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;
	
	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;
	
	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;
	
	
	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	/* return pointer to array of pointers to rows */
	return m;
}

float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
				  long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	float **m;
	
	/* allocate array of pointers to rows */
	m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in submatrix()");
	m += NR_END;
	m -= newrl;
	
	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;
	
	/* return pointer to array of pointers to rows */
	return m;
}

float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
 declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
 and ncol=nch-ncl+1. The routine should be called with the address
 &a[0][0] as the first argument. */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;
	
	/* allocate pointers to rows */
	m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;
	
	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	float ***t;
	
	/* allocate pointers to pointers to rows */
	t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;
	
	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;
	
	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;
	
	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}
	
	/* return pointer to array of pointers to rows */
	return t;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(unsigned char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
				   long ndl, long ndh)
/* free a float f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

#else /* ANSI */
/* traditional - K&R */

#include <stdio.h>
#define NR_END 1
#define FREE_ARG char*

void nrerror(error_text)
char error_text[];
/* Numerical Recipes standard error handler */
{
	void exit();
	
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

float *vector(nl,nh)
long nh,nl;
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;
	
	v=(float *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;
	
	v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

unsigned char *cvector(nl,nh)
long nh,nl;
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	unsigned char *v;
	
	v=(unsigned char *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) nrerror("allocation failure in cvector()");
	return v-nl+NR_END;
}

unsigned long *lvector(nl,nh)
long nh,nl;
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;
	
	v=(unsigned long *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}

double *dvector(nl,nh)
long nh,nl;
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;
	
	v=(double *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

float **matrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;
	
	/* allocate pointers to rows */
	m=(float **) malloc((unsigned int)((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;
	
	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	/* return pointer to array of pointers to rows */
	return m;
}

double **dmatrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;
	
	/* allocate pointers to rows */
	m=(double **) malloc((unsigned int)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;
	
	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;
	
	/* allocate pointers to rows */
	m=(int **) malloc((unsigned int)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;
	
	
	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	/* return pointer to array of pointers to rows */
	return m;
}

float **submatrix(a,oldrl,oldrh,oldcl,oldch,newrl,newcl)
float **a;
long newcl,newrl,oldch,oldcl,oldrh,oldrl;
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	float **m;
	
	/* allocate array of pointers to rows */
	m=(float **) malloc((unsigned int) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in submatrix()");
	m += NR_END;
	m -= newrl;
	
	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;
	
	/* return pointer to array of pointers to rows */
	return m;
}

float **convert_matrix(a,nrl,nrh,ncl,nch)
float *a;
long nch,ncl,nrh,nrl;
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
 declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
 and ncol=nch-ncl+1. The routine should be called with the address
 &a[0][0] as the first argument. */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;
	
	/* allocate pointers to rows */
	m=(float **) malloc((unsigned int) ((nrow+NR_END)*sizeof(float*)));
	if (!m)	nrerror("allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;
	
	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

float ***f3tensor(nrl,nrh,ncl,nch,ndl,ndh)
long nch,ncl,ndh,ndl,nrh,nrl;
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	float ***t;
	
	/* allocate pointers to pointers to rows */
	t=(float ***) malloc((unsigned int)((nrow+NR_END)*sizeof(float**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;
	
	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(float **) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(float*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;
	
	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(float *) malloc((unsigned int)((nrow*ncol*ndep+NR_END)*sizeof(float)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;
	
	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}
	
	/* return pointer to array of pointers to rows */
	return t;
}

void free_vector(v,nl,nh)
float *v;
long nh,nl;
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(v,nl,nh)
int *v;
long nh,nl;
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(v,nl,nh)
long nh,nl;
unsigned char *v;
/* free an unsigned char vector allocated with cvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(v,nl,nh)
long nh,nl;
unsigned long *v;
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(v,nl,nh)
double *v;
long nh,nl;
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(m,nrl,nrh,ncl,nch)
float **m;
long nch,ncl,nrh,nrl;
/* free a float matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_dmatrix(m,nrl,nrh,ncl,nch)
double **m;
long nch,ncl,nrh,nrl;
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
long nch,ncl,nrh,nrl;
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_submatrix(b,nrl,nrh,ncl,nch)
float **b;
long nch,ncl,nrh,nrl;
/* free a submatrix allocated by submatrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(b,nrl,nrh,ncl,nch)
float **b;
long nch,ncl,nrh,nrl;
/* free a matrix allocated by convert_matrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_f3tensor(t,nrl,nrh,ncl,nch,ndl,ndh)
float ***t;
long nch,ncl,ndh,ndl,nrh,nrl;
/* free a float f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

#endif /* ANSI */

#define NRANSI
#include "nrutil.h"
#define TOL 2.0e-4

int ncom;
float *pcom,*xicom,(*nrfunc)(float []);

void linmin(float p[], float xi[], int n, float *fret, float (*func)(float []))
{
	float brent(float ax, float bx, float cx,
				float (*f)(float), float tol, float *xmin);
	float f1dim(float x);
	void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb,
				float *fc, float (*func)(float));
	int j;
	float xx,xmin,fx,fb,fa,bx,ax;
	
	ncom=n;
	pcom=vector(1,n);
	xicom=vector(1,n);
	nrfunc=func;
	for (j=1;j<=n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
	*fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
	for (j=1;j<=n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	free_vector(xicom,1,n);
	free_vector(pcom,1,n);
}
#undef TOL
#undef NRANSI








#include <math.h>
#define NRANSI
#include "nrutil.h"
#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

float brent(float ax, float bx, float cx, float (*f)(float), float tol,
			float *xmin)
{
	int iter;
	float a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	float e=0.0;
	
	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*f)(x);
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=(*f)(u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	nrerror("Too many iterations in brent");
	*xmin=x;
	return fx;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef NRANSI













#include <math.h>
#define NRANSI
#include "nrutil.h"
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb, float *fc,
			float (*func)(float))
{
	float ulim,u,r,q,fu,dum;
	
	*fa=(*func)(*ax);
	*fb=(*func)(*bx);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*func)(*cx);
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
		(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,(*func)(u))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=(*func)(u);
		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT
#undef NRANSI















#define NRANSI
#include "nrutil.h"

extern int ncom;
extern float *pcom,*xicom,(*nrfunc)(float []);

float f1dim(float x)
{
	int j;
	float f,*xt;
	
	xt=vector(1,ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=(*nrfunc)(xt);
	free_vector(xt,1,ncom);
	return f;
}
#undef NRANSI

