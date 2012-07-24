

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
#include "pbt_geno_compare.h"
#include "pbt_trio_mixture.h"





/* this function returns the array U for the trio_mixture_offcoll_struct.  It is pretty
 lame how I have this hard-wired to relate to the following specification of pedigrees:
 
 "C_Se_Se",       0      
 "C_Se_Si",       1      
 "C_Si_Se",       2      
 "C_Se_U",       3      
 "C_U_Se",       4      
 "C_Si_Si",       5      
 "C_Si_U",       6      
 "C_U_Si",       7      
 "C_U_U",       8      
 "Se_B",       9      
 "B_Se",       10      
 "H_Se",       11      
 "Se_H",       12      
 "B_Si",       13      
 "Si_B",       14      
 "B_U",       15      
 "U_B"       16   
 "B_B"		  17

 So, if I change the ped spec in Print_snpSumPed_CommandFile() then I have to change things here too.
 
 But, I just want to get this implemented to see if it has a prayer of working!
 
 */
int **ReturnPedSpecSelfIndicators(void)  
{
	int i;
	int **ret;
	
	ret = (int **)calloc(NUM_SPEC_PEDS,sizeof(int *));
	for(i=0;i<NUM_SPEC_PEDS;i++)  {
		ret[i] = (int *)calloc(2,sizeof(int));
	}
	
	/* now set the values according to the "Se"'s above.  i.e. C_Se_Si means the father is Self and the mother is Sib */
	ret[0][0]=1;	ret[0][1]=1;
	ret[1][0]=1;
	ret[2][1]=1;
	ret[3][0]=1;
	ret[4][1]=1;
	ret[9][0]=1;
	ret[10][1]=1;
	ret[11][1]=1;
	ret[12][0]=1;
	
	return(ret);
}


/*
 Prepare all the info for this trio mixture operation.  OC_idx is the index of the offspring collection
 under consideration.  This assumes that Compat Trios have already been found for that collection.
*/
trio_mixture_offcoll_struct *PrepareTrioMixtureOffCollStruct(pbt_high_level_vars *HLV, int OC_idx)
{
	int i,j,k,s;
	pfr_offspring *Offs=HLV->PFR->Offs[OC_idx];
	int NumOffs = HLV->PFR->NumInOffColls[OC_idx];
	trio_mixture_offcoll_struct *ret=(trio_mixture_offcoll_struct *)malloc(sizeof(trio_mixture_offcoll_struct));
	int NumPops=HLV->PFR->NumPops;
	
	
	
	/************ memory allocation and some necessary copying of values *************/
	ret->pi = (double *)calloc(NumPops,sizeof(double));
	ret->pi_prior = (double *)calloc(NumPops,sizeof(double));
	ret->phi = (double *)calloc(NumPops,sizeof(double));
	ret->phi_prior = (double *)calloc(NumPops,sizeof(double));
	ret->theta = (double **)calloc(NumPops,sizeof(double *));
	ret->theta_prior = (double **)calloc(NumPops,sizeof(double *));
	for(i=0;i<NumPops;i++)  {
		ret->theta[i] = (double *)calloc(NUM_SPEC_PEDS,sizeof(double));
		ret->theta_prior[i] = (double *)calloc(NUM_SPEC_PEDS,sizeof(double));
	}
	ret->setZ = (int *)calloc(NumPops,sizeof(int));
	ret->setW = (int *)calloc(NumPops,sizeof(int));
	ret->NumPairs = (int *)calloc(NumOffs,sizeof(int));
	ret->W = (int **)calloc(NumOffs,sizeof(int *));
	for(i=0;i<NumOffs;i++) {
		ret->NumPairs[i] = Offs[i].MendComps->NumParPairs;							/* copying values for NumPairs here */
		if(ret->NumPairs[i]>0) {
			ret->W[i] = (int *)calloc(ret->NumPairs[i],sizeof(int));
		}
		else {
			ret->W[i]=NULL;
		}
	}

	ret->V = (int ****)calloc(NumOffs,sizeof(int ***));
	for(i=0;i<NumOffs;i++)  {
		if(ret->NumPairs[i]>0) {
			ret->V[i] = (int ***)calloc(ret->NumPairs[i],sizeof(int **));
			for(j=0;j<ret->NumPairs[i];j++)  {
				ret->V[i][j] = (int **)calloc(ret->NumPairs[i],sizeof(int *));
				for(k=0;k<ret->NumPairs[i];k++)  {
					ret->V[i][j][k] = (int *)calloc(2,sizeof(int));
				}
			}
		}
		else {
			ret->V[i] = NULL;
		}
	}

	ret->PairlessP = (double **)calloc(NumOffs,sizeof(double *));
	for(i=0;i<NumOffs;i++)  {
		ret->PairlessP[i] = (double *)calloc(NumPops,sizeof(double));
	}
	
	ret->P = (double ****)calloc(NumOffs,sizeof(double ***));
	for(i=0;i<NumOffs;i++)  {
		if(ret->NumPairs[i]>0) {
			ret->P[i] = (double ***)calloc(ret->NumPairs[i],sizeof(double **));
			for(j=0;j<ret->NumPairs[i];j++) {
				ret->P[i][j] = (double **)calloc(2,sizeof(double *));
				ret->P[i][j][0] = (double *)calloc(NUM_SPEC_PEDS,sizeof(double));
				ret->P[i][j][1] = (double *)calloc(NumPops,sizeof(double));
			}
		}
		else  {
			ret->P[i]=NULL;
		} 
	}
	
	
	
	/******** setting of values, do first the ones that we will need later ************/
	/* copy setZ */
	for(i=0;i<NumPops;i++)  {
		ret->setZ[i] = HLV->PFR->OffCollPossPops[OC_idx][i];
	}
	
	/* fill W and setW */
	for(i=0;i<NumOffs;i++)  {
		for(j=0;j<ret->NumPairs[i];j++)  {
			/* here, assume that the pop is the same for putative father and mother...It certainly should be! */
			ret->W[i][j] = Offs[i].MendComps->ParPairs[j].mamapapa[0]->Pop;
			
			/* now record that in setW too */
			ret->setW[ ret->W[i][j] ] = 1;
		}
	}
	
	/* get U */
	ret->U = ReturnPedSpecSelfIndicators();
	/* now, get the V's */
	for(i=0;i<NumOffs;i++)  {
		for(j=0;j<ret->NumPairs[i];j++)  {
			for(k=0;k<ret->NumPairs[i];k++)  {
				for(s=0;s<2;s++)  { pfr_parent *par1; pfr_parent *par2;
					par1 = Offs[i].MendComps->ParPairs[j].mamapapa[s];
					par2 = Offs[i].MendComps->ParPairs[k].mamapapa[s];
					
					if(par1->AbsIdx == par2->AbsIdx) {
						ret->V[i][j][k][s] = 1;
						//printf("Setting_V: OC_idx= %d  ind= %s   i= %d    j= %d    k= %d   s= %d   V= %d    par1= %s   par2= %s\n",OC_idx,Offs[i].Name,i,j,k,s,ret->V[i][j][k][s],par1->Name,par2->Name);
					}
				}
			}
		}
	}
	
	/* now, let's set the values in the P array and also in PairlessP */
	printf("P_Matrix: / PairlessP:  Setting P for Offspring Collection %d\n",OC_idx);
	for(i=0;i<NumOffs;i++)  { pfr_offspring *kid;
		kid = &(Offs[i]);
//if(OC_idx==2 && i==NumOffs-1) { printf("***********	printing ret->P[0][0][0][0] = %e assigned kid pointer...\n",ret->P[0][0][0][0]);}
		/* cycle over pops to get PairlessP's */
		for(j=0;j<NumPops;j++)  {
//if(OC_idx==2 && i==NumOffs-1) { printf("***********	printing ret->P[0][0][0][0] = %e about to assign to pairless p[%d]...\n",ret->P[0][0][0][0],j);}
			ret->PairlessP[i][j] = SimpleIndivGenotypeProb(kid->geno[0],j,HLV);
//if(OC_idx==2 && i==NumOffs-1) { printf("***********	printing ret->P[0][0][0][0] = %e done assigning to pairless p[%d]...\n",ret->P[0][0][0][0],j);}
			printf("PairlessP: Offspring: %s   KidPop:  %d   Prob:  %e\n",kid->Name,j,ret->PairlessP[i][j]);
		}
//if(OC_idx==2 && i==NumOffs-1) { printf("***********	printing ret->P[0][0][0][0] = %e after pairless p...\n",ret->P[0][0][0][0]);}
		for(j=0;j<ret->NumPairs[i];j++)  { pfr_parent *ma; pfr_parent *pa; int mapa_pop;
			pa = Offs[i].MendComps->ParPairs[j].mamapapa[MALE];
			ma = Offs[i].MendComps->ParPairs[j].mamapapa[FEMALE];
			mapa_pop = ma->Pop;
			
//if(OC_idx==2 && i==NumOffs-1) { printf("***********	printing ret->P[0][0][0][0] = %e down past mapa_pop assignment...\n",ret->P[0][0][0][0]);}
			
			for(k=0;k<NUM_SPEC_PEDS;k++)  {
				ret->P[i][j][0][k] = CondProbOfMultiLocusGenotypeDataOfTrio(kid,mapa_pop,pa,ma,k,HLV);
				printf("P_Matrix_SamePops: Trio: %s %s %s   Pedigree: %d   Prob:  %e\n",kid->Name,pa->Name,ma->Name,k,ret->P[i][j][0][k]);
			}
			for(k=0;k<NumPops;k++)  {
				ret->P[i][j][1][k] = CondProbOfMultiLocusGenotypeDataOfTrio(kid,k,pa,ma,8,HLV);
				printf("P_Matrix_CrossPops: Trio: %s %s %s   KidPop:  %d   Prob:  %e\n",kid->Name,pa->Name,ma->Name,k,ret->P[i][j][1][k]);
			}
/*			if(OC_idx==2 && i==0) {
				printf("*********************** GOT HERE ************************  P[0][0] = %x    P[0][0][0][0]=%e\n",ret->P[0][0],ret->P[0][0][0][0]);
				WatchIt=ret->P[0][0];
				onit=1;
			}
*/
		}
//if(OC_idx==2) {printf("***********	printing ret->P[0][0][0][0] = %e for fun...\n",ret->P[0][0][0][0]);}
	}
	
	/* finally, set the values of NumOffs and NumPop */
	ret->NumPops = NumPops;
	ret->NumOffs = NumOffs;
	ret->OC_idx = OC_idx;
	
	
		
	return(ret);
	
}


/*
 Given probs of A states 0 up to 26 for a trio. These are like so:
	idx:   kid	pa	ma
	0:		0	0	0
	1:		0	0	1
	2:		0	0	2
	3:		0	1	0
 and so forth.
 
 And given the genotypes of the individuals, with missing data denoted
 by 3.
 
 This function returns the probability of the genotype of the trio,
 summing over the missing data as necessary.
 
*/
double ProbOfSingleLocusGenotypeDataOfTrio(char kid, char pa, char ma, double *Aprobs)
{
	char i,j,k;
	int l=0;  /* for cycling over kid, pa, ma, respectively */
	double sum=0.0;
	
	for(i=0;i<=2;i++)  {
		for(j=0;j<=2;j++)  {
			for(k=0;k<=2;k++) {
				if(  (kid==i || kid==3) && 
				     (pa==j  ||  pa==3) &&
				   (ma==k  ||  ma==3) ) {
					sum+=Aprobs[l];
				}
				l++;
			}
		}
	}
	return(sum);
}



/* 
 This is a relatively high level function.  It figures out the pop of the parents
 involved and chooses to use PurePopTrioColls or CrossPopTrioColls as appropriate.  
 
 Then it returns the conditional prob of the trio's genotype given pedigree type
 the_ped and given that it is in S-down.

 This just uses the first of however many regenotyped versions of the individual there are. 
 
*/
double CondProbOfMultiLocusGenotypeDataOfTrio(pfr_offspring *kid,		
										  int kid_pop,					/* the population of the kid.  Usually not known, but given to compute necessary probs */
										  pfr_parent *pa,
										  pfr_parent *ma,
										  int the_ped,					/* note that this is ignored if kid_pop is not the same as the pop of the ma and pa */
										  pbt_high_level_vars *HLV )
{
	int j;
	int ma_pop = ma->Pop;
	int pa_pop = pa->Pop;
	int pick;
	double *Aprobs;
	double prod=1.0;
	
	Collection *C;
	
	
	
	if(ma_pop != pa_pop) {
		fprintf(stderr,"Error!  Ma and Pa are from different populations in CondProbOfMultiLocusGenotypeDataOfTrio.  ma_pop= %d and pa_pop= %d.  ma_name= %s   pa_name= %s   Exiting.\n",
				ma_pop,pa_pop,ma->Name,pa->Name);
		exit(1);
	}
	
	if(ma_pop==kid_pop)  {
		pick = PPTC_IDX(NUM_SPEC_PEDS,  the_ped,  kid_pop);
		C = HLV->PurePopTrioColls->Colls[pick];
	}
	else {
		pick = CPTC_IDX(HLV->PFR->NumPops,  ma_pop,  kid_pop);
		C = HLV->CrossPopTrioColls->Colls[pick];
	}
	
	
	for(j=0;j<HLV->PFR->NumLoci;j++)  {
		Aprobs = C->Aprobs[j];
		prod *= ProbOfSingleLocusGenotypeDataOfTrio(kid->geno[0][j], pa->geno[0][j], ma->geno[0][j], Aprobs);
	}
	
	/* return the prob, divided by the prob of being in S_down, to make it the conditional prob.
	 Note that I could probably do some funny stuff here to try a quick and dirty way of dealing with
	 missing data (i.e. adjust the S_down a little bit */
	return(prod / C->FBS->S_down);
	
}
										  

/* this just returns the prob of an individual's genotype.  Pass in its genotype and the pop you want the 
 probability from and it will give you that prob.  HLV is there for the allele freqs, etc.
 
 Note that this uses the allele freqs which are the posterior means given a Beta(.5,.5) prior. 
 
 Missing data are ignored  (i.e. given a value of 1)  
 
 */
double SimpleIndivGenotypeProb(char *G, int pop, pbt_high_level_vars *HLV)
{
	int i;
	double **F = HLV->PFR->AlleFreqs[pop];
	double prod=1.0;
	int L=HLV->PFR->NumLoci;
	
	for(i=0;i<L;i++)  {
		switch(G[i]) {
			case(0):
				prod *= F[i][0]*F[i][0];
				break;
			case(1):
				prod *= 2.0 * F[i][0] * F[i][1];
				break;
			case(2):
				prod *= F[i][1] * F[i][1];
				break;
			default:  /* missing data = do nothing */
				break;
		}
	}

	return(prod);
}


/*
 Treat every trio as an independent observation and use an EM-like algorithm to
 estimate the fraction of the individuals from each of the populations, and also the 
 fractions of different trio types.  This is just fast and loose, but might 
 be good for initializing the composite likelihood approach taken later.
*/
void SuperNaiveSimpleTrioEM(trio_mixture_offcoll_struct *T)
{
	int i,j,k,iter;
	int NumOffs = T->NumOffs;
	int NumPops = T->NumPops;
	int *setZ = T->setZ;
	int **W = T->W;
	double ****P = T->P;
	
	double *pi = (double *)calloc(NumPops,sizeof(double));
	double *new_pi = (double *)calloc(NumPops,sizeof(double));
	
	double **t = (double **)calloc(NumPops,sizeof(double *));
	double **new_t = (double **)calloc(NumPops,sizeof(double *));
	double normo;
	int tset;
	
	
	for(i=0;i<NumPops;i++) {
		t[i] = (double *)calloc(NUM_SPEC_PEDS,sizeof(double));
		new_t[i] = (double *)calloc(NUM_SPEC_PEDS,sizeof(double));
	}
	
	
	/* start them up with uniform freqs for pops in setZ (for the mixture) and setW for the ped specs */
	normo=0.0; 
	for(i=0;i<NumPops;i++)  {
	 	pi[i] = (double)setZ[i];
		normo += pi[i];
		if(setZ[i]) {
			for(j=0;j<NUM_SPEC_PEDS;j++)  {
				t[i][j] = 1.0/NUM_SPEC_PEDS;
			}
		}
	}
	for(i=0;i<NumPops;i++)  {
	 	pi[i] /= normo;
	}
	
	
	for(iter=0;iter<100;iter++)  {
		
		/* initialize the "new" versions of pi and t to accumulate a sum */
		for(i=0;i<NumPops;i++)  {
			new_pi[i]=0.0;
			for(j=0;j<NUM_SPEC_PEDS;j++)  {
				new_t[i][j] = 0.0;
			}
		}
		
		/* now, cycle over the offspring, and within those over the parent pairs, and start normalizing and average stuff */
		tset=0;
		for(i=0;i<NumOffs;i++)  {
			if(T->NumPairs[i]==0) {
				normo=0.0;
				for(j=0;j<NumPops;j++)  {
					normo += T->PairlessP[i][j] * pi[j];
				}
				for(j=0;j<NumPops;j++)  {
					new_pi[j] += (T->PairlessP[i][j] * pi[j]) / normo;
				}
			}
			else {
				for(j=0;j<T->NumPairs[i];j++)  { int mapop; double norm2;
					tset=1;
					normo = 0.0;
					mapop = W[i][j];
					/* first get the normalizing constant */
					norm2 = 0.0;
					for(k=0;k<NUM_SPEC_PEDS;k++)  {
						normo += P[i][j][0][k] * t[mapop][k] * pi[mapop];
						norm2 += P[i][j][0][k] * t[mapop][k] * pi[mapop];
					}
					for(k=0;k<NumPops;k++) {
						normo += P[i][j][1][k] * pi[k];
					}
					
					/* then record the values in new_t and new_pi */
					for(k=0;k<NUM_SPEC_PEDS;k++)  {
						new_t[mapop][k] += P[i][j][0][k] * t[mapop][j] * pi[mapop] / normo;
					}
					for(k=0;k<NumPops;k++) {
						if(k==mapop) {
							new_pi[k] += norm2/normo;
						}
						else {
							new_pi[k] += P[i][j][1][k] * pi[k]/normo;
						}
					}				
				}
			}
		}
		
		/* now cycle over new_pi and new_t to normalize */
		if(tset) {
			normo=0.0;
			for(i=0;i<NumPops;i++)  {
				for(j=0;j<NUM_SPEC_PEDS;j++)  {
					normo += new_t[i][j];
				}
			}
			for(i=0;i<NumPops;i++)  {
				for(j=0;j<NUM_SPEC_PEDS;j++)  {
					new_t[i][j] /=normo;
				}
			}
		}
		normo=0.0;
		for(i=0;i<NumPops;i++)  {
			normo += new_pi[i];
		}
		for(i=0;i<NumPops;i++)  {
			new_pi[i] /=normo;
		}
		
		/* print those dudes out */
		for(i=0;i<NumPops;i++)  {
			if(setZ[i]) {
				printf("SuperNaiveEMProgress: Iter= %d   OC_idx= %d    Pop= %d  pi= %f    t = ",iter,T->OC_idx,i,new_pi[i]);
				for(j=0;j<NUM_SPEC_PEDS;j++)  {
					printf(" %.4f",new_t[i][j]);
				}
				printf("\n");
			}
		}
		
		/* then swap the new with the originals */
		if(tset) {
			for(i=0;i<NumPops;i++)  {
				for(j=0;j<NUM_SPEC_PEDS;j++)  {
					t[i][j] = new_t[i][j];
				}
			}
		}
		for(i=0;i<NumPops;i++)  {
			pi[i] = new_pi[i];
		}
	}
	
}


void PrintCollectionMix(trio_mixture_offcoll_struct *ret)
{
	int i,j,k;
	/* now, let's set the values in the P array and also in PairlessP */
	printf("TESTING!!! P_Matrix: / PairlessP:  Printing P for Offspring Collection %d\n",ret->OC_idx);
	for(i=0;i<ret->NumOffs;i++)  {
		
		/* cycle over pops to get PairlessP's */
		for(j=0;j<ret->NumPops;j++)  {
			printf("TESTING!!! PairlessP: Offspring: %d   KidPop:  %d   Prob:  %e\n",i,j,ret->PairlessP[i][j]);
		}
		for(j=0;j<ret->NumPairs[i];j++)  { 
			
			for(k=0;k<NUM_SPEC_PEDS;k++)  {
				printf("TESTING!!! P_Matrix_SamePops: Trio: %s %s %s   Pedigree: %d   Prob:  %e\n",i,j,j,k,ret->P[i][j][0][k]);
			}
			for(k=0;k<ret->NumPops;k++)  {
				printf("TESTING!!! P_Matrix_CrossPops: Trio: %s %s %s   KidPop:  %d   Prob:  %e\n",i,j,j,k,ret->P[i][j][1][k]);
			}
		}
	}
	
}