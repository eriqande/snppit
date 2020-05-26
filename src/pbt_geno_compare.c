

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
#include "pbt_C_fb.h"
#include "uthash.h"
#include "pfr_read_genos.h"
#include "pfr_pedigree_spec.h"
#include "pbt_highlevel.h"
#include "pbt_geno_compare.h"



/*
 This function returns a 1 if the number of mendelian incompatibilities between the L-locus
 genotypes Par and Kid is less than or equal to MaxAllowable.  In this case, it returns the
 actual number of incompatible loci in the output parameter NumInc.
 
 If there are more than MaxAllowable incompatibilities, this function returns 0 and the value 
 of NumInc is not defined.
*/
int PairCompatByMendel(char *Par, char *Kid, int MaxAllowable, int L, int *NumInc)
{
	int i;
	char *p=Par,*k=Kid;
	int cnt=0;
	
	for(i=0;i<L;i++,p++,k++) {
		cnt += (*p<3 && *k<3) && ( (*p==0 && *k==2) || (*p==2 && *k==0));
		//printf("i= %d  k= %d  p= %d  cnt= %d\n",i,*k,*p,cnt);
		if(cnt>MaxAllowable) {
			return(0);
		}
	}
	
	*NumInc = cnt;
	return(1);
	
}


/*
 This function returns a 1 if the number of mendelian incompatibilities between the L-locus
 genotypes Pa, Ma, and Kid is less than or equal to MaxAllowable.  In this case, it returns the
 actual number of incompatible loci in the output parameter NumInc.
 
 If there are more than MaxAllowable incompatibilities, this function returns 0 and the value 
 of NumInc is not defined.
 
*/
int TrioCompatByMendel(char *Pa, char *Ma, char *Kid, int MaxAllowable, int L, int *NumInc)
{
	int i;
	char *p=Pa,*k=Kid,*m=Ma;
	int cnt=0;
	
	for(i=0;i<L;i++,p++,k++,m++) {
		cnt += TrioIsIncompat(*m, *p, *k);
		//printf("i= %d  k= %d  p= %d  m= %d  cnt= %d\n",i,*k,*p,*m,cnt);
		if(cnt>MaxAllowable) {
			return(0);
		}
	}
	
	*NumInc = cnt;
	return(1);
	
}





/*
 Returns a one if the possible birth years (of an individual) and the possible reproduction years (of a parent)
 intersect.  Otherwise return 0.
*/
int BirthTimingWorks(
					 int NBY,					/* number of possible birth years of the kid */
					 int *BY,					/* array of possible birth years of the kid */
					 pfr_year_range *RY			/* struct holding the possible reproductive years of the parent */
					 )
{
	int i;
	int Lo;
	int Len;
	char *bits;
	int ret=0;
	int *b;
	
	
	/* If there are no possible birth years for an offspring, we assume it could have been born at any time,
	 AND
		If RY is NULL then we assume that the parent could have been reproductive in any year,
	 SO
		we just return 1 in that case
	 */
	if(NBY==0 || RY==NULL)  {
		return(1);
	}
	
	Lo = RY->Lo;
	Len = RY->Len;
	bits = RY->bits;
	b=BY;
	
	for(i=0;i<NBY;i++,b++)  {
		if( *b-Lo<Len  && *b-Lo>=0 && bits[*b-Lo] ) {
			return(1);
		}	
	}
	
	return(ret);
}


/*
	Returns a one if the two pfr_year_ranges intersect.  0 otherwise.  
	intsct is an output variable to return the first year in which they intersect. 
	It returns -1 if they are null, etc.
*/
int MateTimingWorksOut(pfr_year_range *A, pfr_year_range *B, int *intsct)
{
	int i,year;
	int ALo;
	int ALen;
	char *Abits;
	int BLo;
	int BLen;
	char *Bbits;
	
	
	*intsct = -1; 
	
	/* IF either A or B is NULL it means that mating years were not specified for that individual and we assume that it could 
	 have mated at any time, and hence we want to return a 1 here */
	if(A==NULL || B==NULL) {
		return(1);
	}
	
	
	ALo = A->Lo;
	ALen = A->Len;
	Abits = A->bits;
	BLo = B->Lo;
	BLen = B->Len;
	Bbits = B->bits;
	
	
	/* we will cycle over the values in A, and return 1 if we intersect with any in B */
	for(i=0;i<ALen;i++)  {
		if(Abits[i])  {
			year = ALo + i;
			if( year-BLo<BLen  &&  Bbits[year-BLo] ) {
				*intsct = year;
				return(1);
			}	
		}
	}
	
	return(0);
}






/*
  This is a function that compares the genotype of the youth, Y, to the genotypes of all possible parents
 (using the possible pops and the ages and repro-years).  The ParSex decides whether we are looking for mothers,
 fathers, or SEX_UNKNOWNS.  
 
*/
pfr_parent **ReturnMatchingParents(
									pfr_offspring *Y,			/* all the info about the youth */
									int ReGenoNumKid,				/* which genotype (from various regenotypings) to use for offspring.  Typically will be zero */
									int ReGenoNumPar,				/* which genotype to use for the parent.  Typically zero */
								    pbt_high_level_vars *HLV,			/* all the other info we might need */
									SexEnum ParSex,				/* the sex of the parents we are looking for (MALE/FEMALE/SEX_UNKNOWN) */
									int *NumRet,					/* output parameter---number of parents in returned pfr_parent * array. */
								    int **NumMI						/* another output parameter, for the number of mendelian incompatibilities at each parent in the array. 
																		we have to pass this in as an int ** because we will allocate memory to it in this function. */
								   )
{
	int i,j,NumP;
	
	int Cnt=0;  /* keeping count of how many there are */
	pfr_parent **TempPars;
	int *NumIncompats;
	pfr_parent **ret;
	int *retMI;
	
	
	/* get some local variables on the stack */
	char *G = Y->geno[ReGenoNumKid];
	int NBY = Y->NumBornYears;
	int *BY = Y->BornYears;
	pfr_parent ***Parray = HLV->PFR->Pops;	/* local  array of possible parents */
	pfr_parent *P;
	int *OCPP = HLV->PFR->OffCollPossPops[Y->OffColl];		/* array giving the possible parental pops for this individual  */
	int NumPops = HLV->PFR->NumPops;
	int MaxMiss = HLV->smax_to_use[0];  /* NOTE: HERE I AM JUST ASSUMING THAT smax[0]==smax[1].  Might have to alter that if we run into
									situations where we are dealing with hybrid populations where the males are from one pop and the females from another. */
	int L = HLV->PFR->NumLoci;
	int YAbsIdx = Y->AbsIdx;
	
	
	/* allocate some temporary memory that we will free at the end of the function */
	TempPars = (pfr_parent **)calloc(MAX_TEMP_PARS, sizeof(pfr_parent *));
	NumIncompats = (int *)calloc(MAX_TEMP_PARS,sizeof(int));

	
	for(i=0;i<NumPops;i++)  {  /* cycle over parental pops */
		if(OCPP[i])  {  /* make sure that this is a population from which this offspring could have descended */
			NumP = HLV->PFR->NumInPops[i][ParSex];    /* get the number of parents and check to make sure it is greater than 0 */
			if(NumP>0) {
				P = Parray[i][ParSex];
			}
			for(j=0;j<NumP;j++,P++) {  int tempNumInc; /* cycle over the parents of sex ParSex in this pop  */
				if( YAbsIdx != P->AbsIdx &&  /* make sure that P is not Y in the parent pool */
				    P->NumMissingLoci[ReGenoNumPar] <= HLV->PBUO->MaxAllowMissingLociInParents && /* don't consider the parent if it has too much missing data */
					BirthTimingWorks(NBY, BY, P->ReproYears)  &&   
					PairCompatByMendel(P->geno[ReGenoNumPar], G, MaxMiss, L, &tempNumInc)  ) {
					TempPars[Cnt] = P;
					NumIncompats[Cnt] = tempNumInc;
					
					Cnt++;
				}
			}
		}
	}
	
	
	if(Cnt==0) {
		ret = NULL;
		retMI = NULL;
	}
	else {
		ret = (pfr_parent **)calloc(Cnt, sizeof(pfr_parent *));
		retMI = (int *)calloc(Cnt,sizeof(int));
		for(i=0;i<Cnt;i++)  {
			ret[i] = TempPars[i];
			retMI[i] = NumIncompats[i];
		}
	}

	/* assign values to the output parameters */
	*NumRet = Cnt;
	*NumMI = retMI;
	
	
	free(TempPars);
	free(NumIncompats);
	return(ret);
}
	

/*
 Given genotypes scored as 0,1, and 2, (number of copies of the 
 minor allele) for mother, father, and kid, OR a 3 for missing data, 
 this is the indicator function as to whether the trio is 
 incompatible (1) or compatible (0), at the locus.  
 
 I am sure I could make this more efficient, but I don't think that will
 be necessary. 
 
 */
int TrioIsIncompat(char m, char f, char y) {
	/* deal first with missing data situations */
	if(y==3) return(0);
	if( (m==3 && f<3) ) {
		if( (f==0 && y==2) || (f==2 && y==0) )  { return(1); }
		else {return(0);}
	}  
	if( (f==3 && m<3)  ) {
		if( (m==0 && y==2) || (m==2 && y==0) ) { return(1); }
		else { return(0); }
	} 
	if( f==3 && m==3 ) return(0);
	if(m+f==0 && y>0) return(1);
	if(m+f==1 && y>1) return(1);
	if(m+f==2 && m!=f && y!=1) return(1);
	if(m+f==3 && y<1) return(1);
	if(m+f==4 && y<2) return(1);
	return(0);
}



/*
Return a 1 if the two individuals could belong to the same spawning group and a 
 0 otherwise.  THIS IS CURRENTLY JUST A QUICK HACK ASSUMING THAT THERE IS A 
 SINGLE HIERARCHICAL LEVEL OF SPAWNING GROUP.  I WILL NEED TO EXPAND THIS IN THE 
 FUTURE.  NOTE ALSO THAT I SHOULD ALLOW INDIVIDUALS TO BE INCLUDED IN MULTIPLE 
 SPAWNING GROUPS (FOR REPEAT SPAWNERS, ETC). 
*/
int SpawnGroupMatches( pfr_parent *Ma, pfr_parent *Pa)
{
	
	int i,j;
  int nsgMa,nsgPa;
  
  nsgMa = Ma->NumSpawnGroups;
  nsgPa = Pa->NumSpawnGroups;
  
  /* note that we should only ever have 0's for the SpawnGroup in the first position
   * of the SpawnGroup array.  So test for it here, and throw an error if we see it later.
   * A spawning group of 0 means missing data, so the fish matches any spawn group.
   */
  if(Ma->SpawnGroup[0]==0 || Pa->SpawnGroup[0]==0)  {
    return(1);
  }
  
  for(i=0;i<nsgMa;i++) {
    for(j=0;j<nsgPa;j++) {
      
      /* Throw an error if we see other zeroes in there. */
      if(Ma->SpawnGroup[i]==0 || Pa->SpawnGroup[j]==0)  {
        fprintf(stderr, "Error! Missing Spawn Group in token > 0 for fish %s or %s\n", Ma->Name, Pa->Name);
        exit(1);
      }
      if( Ma->SpawnGroup[i] == Pa->SpawnGroup[j] ) {
        return(1);
      }
    }
  }
  
	return(0);
}


/*
Returns a 1 if a pair of individuals could reasonably be mates.  They must:
	1. be from the same population
	2. have reproduced in at least one year together  
	3. must be in the same spawning group
*/
int MateTimeAndPlaceWorksOut( pfr_parent *Ma, pfr_parent *Pa)
{
	int dummy;
	
	if( Ma->Pop != Pa->Pop ) {
		return(0);
	}
	if(!MateTimingWorksOut(Ma->ReproYears, Pa->ReproYears,&dummy)) {
		return(0);
	}
	if(!SpawnGroupMatches(Ma,Pa)) {
		return(0);
	}
	
	return(1);
	
}


/* assuming that we have already identified the possible single parents for the individual, this function 
 checks the possible pairs of those single parents taking account of sex, and population, and repro years, and 
 spawner group.  If it passes all those tests, then it counts the number of incompatibilities in the trio as a whole
 and compares that to the MaxAllowable */
ParentPairMatch *ReturnMatchingTrios(
									pfr_offspring *Y,						/* pointer to the offspring */
									pbt_high_level_vars *HLV,				/* pointer to everything */
									int MaxAllowable,							/* max allowable number of mendelian incompatibilities at the trio */
									int ReGenoMa,							/* re-genotyping numbers for ma, pa, kid---typically 0 */
									int ReGenoPa,
									int ReGenoKid,
									int *NumTrios							/* output parameters telling us how many trios were compatible */
									) 
{
	int i,j;
	
	int Cnt=0;  /* keeping count of how many there are */
	
	ParentPairMatch *TempPairs;
	ParentPairMatch *ret;
	pfr_parent *ma,*umom;
	pfr_parent *pa,*upa;
	int NumInc;
	
	TempPairs = (ParentPairMatch *)calloc(MAX_TEMP_PAIRS, sizeof(ParentPairMatch));
	
	
	/* get some local variables on the stack */
	char *G = Y->geno[ReGenoKid];
	int L = HLV->PFR->NumLoci;
	pfr_parent **Mas = Y->MendComps->Parents[FEMALE];
	int NumMas = Y->MendComps->NumParents[FEMALE];
	pfr_parent **Pas = Y->MendComps->Parents[MALE];
	int NumPas = Y->MendComps->NumParents[MALE];
	pfr_parent **Unks = Y->MendComps->Parents[SEX_UNKNOWN];
	int NumUnks = Y->MendComps->NumParents[SEX_UNKNOWN];
	int **Num_MI_Parents=Y->MendComps->Num_MI_Parents;

	
	/* form all possible pairs of parents.  If sex of both parents is known, use that to determine if they can be
	 a pair or not.  If sex of one parent is known, coerce the sex of the other parent (since it is unknown) to be 
	 the opposite.  If both parents are unknown, assign the one with the lowest AbsIdx to be the MALE */
	
	if(NumPas>0) for(i=0;    i<NumPas;     i++)  {  /* cycle over possible dads */
		pa=Pas[i];
		if(NumMas>0) for(j=0;    j<NumMas;   j++) { /* cycle over moms within dads */
			ma=Mas[j];
			if(MateTimeAndPlaceWorksOut(ma,pa) &&
			   TrioCompatByMendel(pa->geno[ReGenoPa], ma->geno[ReGenoMa], G, MaxAllowable, L, &NumInc)) {
				TempPairs[Cnt].mamapapa[MALE]=pa;
				TempPairs[Cnt].mamapapa[FEMALE]=ma;
				TempPairs[Cnt].pop = pa->Pop;
				TempPairs[Cnt].ma_MI = Num_MI_Parents[FEMALE][j];
				TempPairs[Cnt].pa_MI = Num_MI_Parents[MALE][i];
				TempPairs[Cnt].trio_MI = NumInc;
				Cnt++;
			}
		}
		
		
		if(NumUnks>0) for(j=0;    j<NumUnks;   j++) {  /* cycle over unknowns as moms within dads */
			ma=Unks[j];
			if(MateTimeAndPlaceWorksOut(ma,pa) &&
			   TrioCompatByMendel(pa->geno[ReGenoPa], ma->geno[ReGenoMa], G, MaxAllowable, L, &NumInc)) {
				TempPairs[Cnt].mamapapa[MALE]=pa;
				TempPairs[Cnt].mamapapa[FEMALE]=ma;
				TempPairs[Cnt].pop = pa->Pop;
				TempPairs[Cnt].ma_MI = Num_MI_Parents[SEX_UNKNOWN][j];
				TempPairs[Cnt].pa_MI = Num_MI_Parents[MALE][i];
				TempPairs[Cnt].trio_MI = NumInc;
				Cnt++;
			}
		}
	}
	
	/* now we cycle over mothers and within them cycle over unknowns as possible dads */
	
	if(NumMas>0) for(j=0;    j<NumMas;   j++) { /* cycle over moms */
		ma=Mas[j];
		if(NumUnks>0) for(i=0;    i<NumUnks;   i++) {  /* cycle over unknowns as dads within moms */
			pa=Unks[i];
			if(MateTimeAndPlaceWorksOut(ma,pa) &&
			   TrioCompatByMendel(pa->geno[ReGenoPa], ma->geno[ReGenoMa], G, MaxAllowable, L, &NumInc)) {
				TempPairs[Cnt].mamapapa[MALE]=pa;
				TempPairs[Cnt].mamapapa[FEMALE]=ma;
				TempPairs[Cnt].pop = pa->Pop;
				TempPairs[Cnt].ma_MI = Num_MI_Parents[FEMALE][j];
				TempPairs[Cnt].pa_MI = Num_MI_Parents[SEX_UNKNOWN][i];
				TempPairs[Cnt].trio_MI = NumInc;
				Cnt++;
			}
		}
	}
	
	
	/* finally, we cycle over all ordered pairs of the SEX_UNKNOWN individuals */
	if(NumUnks>1) for(i=0;    i<NumUnks;   i++) {  /* cycle over first member of unknown pairs  which is father if rev = 0 (mother if rev=1) */
		if(NumUnks>1) for(j=i+1; j<NumUnks; j++)  { int rev; /* cycle over second member of unknown pairs which is mother if rev=0 (father if rev=1) */
			upa = Unks[i];
			umom = Unks[j];
			if(upa->AbsIdx < umom->AbsIdx) {
				pa=upa;
				ma=umom;
				rev=0;
			}
			else {
				pa=umom;
				ma=upa;
				rev=1;
			}
			if(MateTimeAndPlaceWorksOut(ma,pa) &&
			   TrioCompatByMendel(pa->geno[ReGenoPa], ma->geno[ReGenoMa], G, MaxAllowable, L, &NumInc)) {
				TempPairs[Cnt].mamapapa[MALE]=pa;
				TempPairs[Cnt].mamapapa[FEMALE]=ma;
				TempPairs[Cnt].pop = pa->Pop;
				TempPairs[Cnt].trio_MI = NumInc;
				if(rev==0) {
					TempPairs[Cnt].ma_MI = Num_MI_Parents[SEX_UNKNOWN][j];
					TempPairs[Cnt].pa_MI = Num_MI_Parents[SEX_UNKNOWN][i];
				}
				if(rev==1) {
					TempPairs[Cnt].ma_MI = Num_MI_Parents[SEX_UNKNOWN][i];
					TempPairs[Cnt].pa_MI = Num_MI_Parents[SEX_UNKNOWN][j];
				}
				Cnt++;
			}
			
		}
	}
	
	/* At this point, all we have to do is allocate space to ret and then copy the contents in TempPairs into it! */
	if(Cnt>0) {
		ret = (ParentPairMatch *)calloc(Cnt, sizeof(ParentPairMatch));
		#ifdef VERBOSE_PRINT_COMPAT_TRIOS
			printf("COMPAT_TRIOS: ind= %s  Num=%d  : ",Y->Name,Cnt);
		#endif
		for(i=0;i<Cnt;i++)  { int k;
			ret[i].mamapapa[MALE]      =   TempPairs[i].mamapapa[MALE];
			ret[i].mamapapa[FEMALE]    =   TempPairs[i].mamapapa[FEMALE];
			ret[i].pop                 =   TempPairs[i].pop;
			ret[i].ma_MI               =   TempPairs[i].ma_MI;
			ret[i].pa_MI               =   TempPairs[i].pa_MI;
			ret[i].trio_MI             =   TempPairs[i].trio_MI;
			
			#ifdef VERBOSE_PRINT_COMPAT_TRIOS
				printf("%s--x--%s  [%d] (%d %d %d)   ",ret[i].mamapapa[MALE]->Name,  ret[i].mamapapa[FEMALE]->Name, ret[i].pop, ret[i].ma_MI, ret[i].pa_MI,ret[i].trio_MI);
			#endif
			
			/* here we want to record the A_64_idx of each locus at this trios.  Note that, typically, the ParentalTrioBigSmax is the only FBVArs struct that
			 will have been filled out by this point.  */
			ret[i].A_64_idx = A_64_ArrayFrom3IndivCharVecs(Y->geno[ReGenoKid], ret[i].mamapapa[MALE]->geno[ReGenoPa], ret[i].mamapapa[FEMALE]->geno[ReGenoMa], HLV->ParentalTrioBigSmax->RP);
			
			/* and now that we have done that we shall make an array parallel to A_64_idx that has the X-states of all the loci */
			ret[i].X_states_64 = (int *)calloc( HLV->ParentalTrioBigSmax->RP->L, sizeof(int) );
			for(k=0;k<HLV->ParentalTrioBigSmax->RP->L;k++)  {
				ret[i].X_states_64[k] = ReturnXstateOfA_64(ret[i].A_64_idx[k]);
			}
						
		}
		
		#ifdef VERBOSE_PRINT_COMPAT_TRIOS
			printf("\n");
		#endif
		
		
		
		#ifdef VERBOSE_PRINT_COMPAT_TRIOS
			/* here as a test, we print out some lines showing the genotypes */
			for(i=0;i<Cnt;i++) { int j;
				printf("\t\tTRIO_GENOTYPES of %s with %s--x--%s   :   ",Y->Name, ret[i].mamapapa[MALE]->Name,  ret[i].mamapapa[FEMALE]->Name);
				for(j=0;j<L;j++)  {
					printf("%d ",ret[i].A_64_idx[j]);
				}
				printf("\n");
			}
		#endif
		
	}
	else {
		/*printf("COMPAT_TRIOS: ind= %s  Num=%d  \n: ",Y->Name,0); */
		ret=NULL;
	}
	
	*NumTrios=Cnt;
	
	free(TempPairs);
	return(ret);
}
