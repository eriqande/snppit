
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <limits.h>


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





pbt_user_opts *GetPBT_UserOpts(int argc, char *argv[])
{
	pbt_user_opts *ret;
	
	int file_f=0,
		big_smax_f=0,
		mi_fnr_f=0,
		forback_help_f=0,
		trio_cat_prob_f=0,	
		max_allow_par_miss_f=0,
		dry_run_f=0,
		psz_for_all_f=0,
		min_logl_f=0,
		max_par_pair_f=0;
	DECLARE_ECA_OPT_VARS  
	
	/* This is a good place to set some default values for variables in your program */
	
	
	/* some information about the program, author, etc */
	SET_PROGRAM_NAME("snppit");  /* Note:  QUOTED string */
	SET_PROGRAM_SHORT_DESCRIPTION("fast and accurate parentage with SNPs"); /* QUOTED string */
	SET_PROGRAM_LONG_DESCRIPTION(Fast and accurate inference of parent pairs with SNPs);  /* UN-QUOTED string! */
	SET_PROGRAM_AUTHOR_STRING("Eric C. Anderson"); /* QUOTED string */
	SET_VERSION("Version 1.0.1");  /* QUOTED string */
	SET_VERSION_HISTORY("Version 1.0.1 June 20, 2019 --- includes Logl and Rank filtering for non-excluded pairs\nVersion 1.0.0 written Nov. 8, 2007 --- initial release."); /*QUOTED string */
	
	/* setting some defaults */
	ret = (pbt_user_opts *)malloc(sizeof(pbt_user_opts));
	ret->big_smax = (int *)calloc(3,sizeof(int));
	ret->big_smax[0] = 10;
	ret->big_smax[1] = 10;
	ret->big_smax[2] = 15;
	ret->d_from_big_smax = 3;
	ret->MI_fnr_target = .005;
	ret->MaxAllowMissingLociInParents = 10;
	ret->Pi =  NULL;
	gHasSnpSumPedPath = 0;
	ret->DryRun = 0;
	ret->PszForAll = -1;
	ret->Pi = (double **)calloc(MAXPOPS,sizeof(double *));
	ret->EnforceParPairMinLogL = 0;
	ret->ParPairMinLogL = 0.0;
	ret->MaxNumRetainedParPairs = -1;
	
	BEGIN_OPT_LOOP  
	
	
	/* use this to start a block of related options.  Close the block with the CLOSE_SUBSET macro (below)  */
	OPEN_SUBSET(Main Options,  /* Name of Subset for output in --help, --help-full, and --help-nroff */
				MainOptions, /* Name of subset for lumping them together in GuiLiner */
				Blah blah blah /* Description of subset.  Not really used currently */
	)   /* NOTE.  These are all UNQUOTED strings.  Beware of commas. */
	
	if(REQUIRED_OPTION(
					Data File Name,
					file_f,
					f,
					,
					S,
					Path to the data file with genotypes and metadata,
					Path to the data file with genotypes and metadata)) {
		if(ARGS_EQ(1)) { char temp[MAX_FILE_PATH_LENGTH];
			GET_STR(temp);
			ret->DataFileName = (char *)calloc( strlen(temp)+1,sizeof(char));
			sprintf(ret->DataFileName,"%s",temp);
		}
	}
	if(OPTION(Big Initial Smax,
		   big_smax_f,
		   s,
		   big-smax,
		   J1 J2 J3,
		   smax vector to use on first pass when checking false negative rates,
		    smax vector to use on first pass when checking false negative rates. This will
			typically be larger than the final smax that will be used in the analysis.  This 
			is optional.  The default is 10 10 15. But you might want to make it bigger if you have
			very high genotyping error rates and hundreds of SNPs.)) {
		if(ARGS_EQ(3)) { int i;
			for(i=0;i<3;i++)  {
				ret->big_smax[i] = GET_INT;
			}
		}
	}
	
	if(OPTION(Max Allowed Missing Loci in Parents,
			  max_allow_par_miss_f,
			  ,
			  max-par-miss,
			  J,
			  The maximum allowable number of missing loci before a candidate parent is discarded from consideraion,
			  The maximum allowable number of missing (i.e. unsuccessfully genotyped) loci before a candidate parent is discarded from consideraion.  If a candidate parent
			  has more than J missing loci then it will be discarded from the analysis.  Default is 10.)) {
		if(ARGS_EQ(1)) {
			ret->MaxAllowMissingLociInParents = GET_INT;
		}
		
	}
	if(OPTION(
			  Coerce All Pops to Common Size,
			  psz_for_all_f,
			  ,
			  psz-for-all,
			  J,
			  Assume chinook and apply popsize J to all parental populations,
			  If this option is invoked then the population size of all parental populations in the data set will be coerced to size J unless the 
			  population has a pi vector set using the --Pi option.  This option overrides the population sizes given in the data file using the 
			  CHINOOK_AVE_POP_SIZE method. )) {
		if(ARGS_EQ(1)) {
			ret->PszForAll = GET_INT;
		}
		
	}
	if(MULT_USE_OPTION(Trio Category Prior,
			trio_cat_prob_f,
			,
			Pi,
			G R_0 R_2 ... R_17,
			Prior probabilities for 18 trio categories\054 from 0 to 17 in pop in range in G (advanced option),
			Prior probabilities for 18 trio categories\054 from 0 to 17.  I need to write out here what those categories are.  For now\054 just see the function
					   Print_snpSumPed_CommandFile.  These will be used as average Pi values for spawner groups with indices in the range specified by G.  Note that
						   populations are counted starting with 1.  Internally they are represented counting from 0 but we account for that. So count them
						   starting with one!  This can be used multiple times.  Later invocations having the same values in the range G as a previous invocation will
						   overwrite the old one. Note that these values\054 if given on the command line\054 will override any values that might be obtained by
					       specifying the population sizes in the data file using the CHINOOK_AVE_POP_SIZE method.  Values given here will also override the value given 
						   with the --psz-for-all option.  See the output file snppit_output_PopSizesAnPiVectors.txt to see how the population size and/or
					       the pi-vectors were set in any particular run.,
			MAXPOPS)) {
	    if(ARGS_EQ(19)) { int i; int j; int *range; char temp[2000]; double *tempd;
			
			/* first get the range */
			GET_STR(temp);
			range = StringToIntArray(temp, 0, MAXPOPS, 1);
			
			/* then get the actual values */
			tempd = (double *)calloc(18,sizeof(double));
			for(i=0;i<18;i++) {
				tempd[i] = GET_DUB;
			}
			/* now copy them over to their final home as indicated by range */
			for(i=1; i<MAXPOPS; i++)  {
				if(range[i]) {
					ret->Pi[i-1] = (double *)calloc(18,sizeof(double));
					for(j=0;j<18;j++)  {
						ret->Pi[i-1][j] = tempd[j];  /* note that we turn these back into base-0 subscripting here */
					}
				}
			}
			free(range);
			free(tempd);
		 }
	   }
	
	
	if(OPTION(Mend Inc Target FNR,
			  mi_fnr_f,
			  ,
			  mi-fnr,
			  R,
			  false negative rate for trios that you do not want to exceed using mendelian incompatibility,
			   false negative rate for trios that you do not want to exceed using mendelian incompatibility.  This is optional.  The default is
			   0.005.
			  )) {
		if(ARGS_EQ(1)) {
			ret->MI_fnr_target = GET_DUB;
		}
	}
	
	if(OPTION(
			  Print help for forback module,
			  forback_help_f,
			  ,
			  for-back-help,
			  S,
			  print help info for for-back module. S= help\054 help-full\054 help-nroff\054 or help-xml,
			  print help info for for-back module. S= help\054 help-full\054 help-nroff\054 or help-xml.  Pass those
			  strings in without any quotes or any leading dashes.  This will cause the help for the for-back module
			  to be printed to stdout\054 after which the program will terminate.)) {
		if(ARGS_EQ(1)) { char tempstr[1000]; char *tempargv[2]; int tempargc; char tempstr2[1000]; char tempstr3[1000]; 
			GET_STR(tempstr);
			if(!(strcmp(tempstr,"help")==0  ||  /* make sure that S is one of the correct four possibilities */
			strcmp(tempstr,"help-full")==0 ||
			strcmp(tempstr,"help-nroff")==0 ||
			strcmp(tempstr,"help-xml")==0)) {
			   fprintf(stderr,"Error!  Unrecognized string option \"%s\" to option --for-back-help. Exiting\n",tempstr);
			   exit(1);
			}
			   tempargv[0] = tempstr2;
			   tempargv[1] = tempstr3;
			   sprintf(tempargv[0],"pbt_C_forback");
			   sprintf(tempargv[1],"--%s",tempstr);
			   tempargc=2;
			   Get_pbt_C_fb_Opts(tempargc, tempargv);
			   exit(1);
		}
	}
	if(OPTION(
			  Read Data Then Quit,
			  dry_run_f,
			  ,
			  dry-run,
			  ,
			  Simply read data\054 print summary of data\054 then quit. ,
			  Simply read data\054 print summary of data\054 then quit.  Use this to make sure that your data are being read in correctly before committing to a full-blown run.)) {
		if(ARGS_EQ(0)) {
			ret->DryRun = 1;
		}
	}
	if(OPTION(
			  Set Min LogL Threshold,
			  min_logl_f,
			  ,
			  min-logl,
			  R,
			  Discard parent pairs with a Log Likelihood less than R. ,
			  Sometimes there might not be a lot of power to exclude parent pairs on the basis of Mendelian incompatibility\054
			  in which case you end up with a huge number of unexcluded parent pairs.  This can really bog down the TrioPosteriors
			  calculation and the backward simulation step.  One way around the problem might be to enforce a minimum Log Likelihood
			  ratio for the parent pair.  Any pair with a log likelihood ratio less than R will be discarded.  As a consequence\054
			  the backward simulation might not be accurate\054 so you might end up with lower FDRs than you should get for some individuals.
			  That said\054 if there are so many unexcluded parent pairs it probably will not make too much of a difference\054 because
			  the offspring should not get assigned to any of those parents with high confidence.
			  )) {
		if(ARGS_EQ(1)) {
			ret->ParPairMinLogL = GET_DUB;
			ret->EnforceParPairMinLogL = 1;
		}
	}
	if(OPTION(
			  Set Max Number of Parent Pairs,
			  max_par_pair_f,
			  ,
			  max-par-pair,
			  J,
			  Retain only the top J parent pairs for each individual. ,
			  Sometimes there might not be a lot of power to exclude parent pairs on the basis of Mendelian incompatibility\054
			  in which case you end up with a huge number of unexcluded parent pairs. This is an experimental option to only
			  retain the top J parent pairs. This will skew
			  the backward simulation might not be accurate\054 so you might end up with lower FDRs than you should get for some individuals.
			  But this might be a way of dealing with a handful of individuals that have too many unexcluded parents.  Obviously you
			  should not go too low on this!
			  )) {
		if(ARGS_EQ(1)) {
			ret->MaxNumRetainedParPairs = GET_INT;
		}
	}
	
	CLOSE_SUBSET;  /* done with the greeting-related options */
	
	END_OPT_LOOP   /* Macro that loses the Main ECA_Opt Loop and does some error checking 
	 behind the scenes */
	
	return(ret);
}





/* do a run of the forward algorithm with just the parental trios (and all populations) using the big_smax then summarize that 
 in a good way to select the smax that we will use.  
*/
void SelectAnSmax3(pbt_high_level_vars *HLV)
{
	FB_Vars *A;
	int x[3],C;
	int *big_smax = HLV->PBUO->big_smax;
	double Res;
	double target = 1.0-HLV->PBUO->MI_fnr_target;
	int *chosen;
	int gotit;
	char tag;
	FILE *out; 
	
	
	chosen = (int *)calloc(3,sizeof(int));
	
	
	/* do the forback algorithm on parental trios for all pops */
	HLV->ParentalTrioBigSmax = CreateCompleteFB_Vars(-1,-2,0, HLV->PBUO->big_smax, HLV->PBUO->d_from_big_smax, HLV->PFR,  "BigSmax_Input","BigSmax_Output.txt");
	
	/* now that we have done that, we want to find the smax with the smallest values of smax[0] and smax[1] that provides a 
	 probability of being in sdown, given parental, that is less than HLV->PBUO->FNR_thresh  (the user-set false negative
	 rate threshold) for all the populations.  Note that we choose small values of smax[0] and smax[1] (which we constrain to
	 be equal) which means we might have large-ish values of smax[2], but that is fine.  
	 
	 We choose a single smax for all collections, which means that we will take the largest smax[0], smax[1], and smax[2] required
	 across the different populations.
	 
	 */
	
	A = HLV->ParentalTrioBigSmax;
	for(C=0; C<A->NumColls;  C++) {  /* cycle over collections */
		gotit=0;
		for(x[0]=0,x[1]=0; x[0]<=big_smax[0] && x[1]<=big_smax[1];  x[0]++,x[1]++) {
			for(x[2]=x[0]; x[2]<=big_smax[2];  x[2]++)  {
				tag=' ';
				Res =  SumOfS_states3(A, x, C);
				if(Res > target && gotit==0) {
					chosen[0] = chosen[1] = ECA_MAX(chosen[0],x[0]);
					chosen[2] = ECA_MAX(chosen[2],x[2]);
					tag='*';
					gotit=1;
				}
				//printf("Collection %d   %d-%d-%d    %.10e  %c\n",C, x[0],x[1],x[2], Res,tag);
			}
		}
	}
	
	
	if( (out=fopen("snppit_output_ChosenSMAXes.txt","w"))==NULL) {
		fprintf(stderr,"Error! Unable to open file \"snppit_output_ChosenSMAXes.txt\" to write to in from function NegotiatePiVectors(). Exiting!\n");
		exit(1);
	}
	
	fprintf(out,"\n\n\nCHOSEN_SMAX:  Chosen smax is %d-%d-%d which gives FNRs of:\n",chosen[0],chosen[1],chosen[2]);
	for(C=0; C<A->NumColls;  C++) {
		fprintf(out,"CHOSEN_SMAX:\t\tCollection %d (%s) : FNR= %f\n",C,A->Colls[C]->CollName,  1-SumOfS_states3(A,chosen,C));
	}
	fprintf(out,"\n\n\n");
	
	HLV->smax_to_use = chosen;
}


/*
 Given an FB_Vars struct that has been filled out, with an smax greater than or equal to the s argument of this function,
 this function returns the probability that at the last locus, the s vector is less than or equal to the s argument
 for collection C.   This function assumes smax has 3 components.
 
*/
double SumOfS_states3(FB_Vars *A, int *s,  int C)
{
	int i,k,abort=0;
	int d = A->RP->d;
	int *smax = A->RP->smax;
	int L = A->RP->L;
	double sum=0.0;
	int x[3];
	struct fb_cell **f;
	int *vdex = A->vdex;
	
	/* first do some checking */
	if(d!=3) {
		fprintf(stderr,"Error!  SumOfS_states3 was called with smax have %d dimensions.  Not 3. Exiting\n",d);
		exit(1);
	}
	for(i=0;i<d;i++) {
		if(s[i]>smax[i]) {
			fprintf(stderr,"Error.  Component %d of s in SumOfS_states is %d which exceeds the correspondng component of smax which is %d.  Prepare to abort!\n",i,s[i],smax[i]);
			abort=1;
		}
	}
	if(abort) {
		exit(1);
	}
	
	/* set f to point to the fb_struct of the correct collection */
	f = A->Colls[C]->FBS->FB;
	
	
	/* now, cycle over the s-states using the x vector */
	for(x[0]=0; x[0]<=s[0];  x[0]++)  {
		for(x[1]=0; x[1]<=s[1];  x[1]++)  {
			for(x[2]=ECA_MAX(x[0],x[1]); x[2]<=s[2];  x[2]++)  {
				k = vdex[ReturnVDex(d, smax, x)];
				//printf("\t\tx=( %d %d %d )  k = %d    cell_s= ( %d %d %d )  p=%.10e\n",x[0],x[1],x[2], k, f[L-1][k].s[0], f[L-1][k].s[1], f[L-1][k].s[2], f[L-1][k].p);
				sum += f[L-1][k].p;
			}
		}
	}
	
	return(sum);
}


/* this assumes that data have already been input, etc, and that HLV->smax_to_use is already set.
    Typically that will get set in SelectAnSmax3 */
 void ComputePurePopTrioColls(pbt_high_level_vars *HLV) 
{
	int i,j,k;
	
	HLV->PurePopTrioColls = CreateCompleteFB_Vars(-1,-2,-1, HLV->smax_to_use,HLV->PBUO->d_from_big_smax, HLV->PFR,  "PurePop_Input","PurePop_Output.txt");
	
	//printf("Computed PurePopTrioColls:\n");
	for(i=0;i<HLV->PFR->NumPops;i++)  {
			for(j=0;j<NUM_SPEC_PEDS;j++) {
			fflush(stdout);
			k = PPTC_IDX(NUM_SPEC_PEDS,  j,  i);
			//printf("\t\tPopIdx= %d  PedIdx= %d  Num= %d  Name= %s   Sup= %.10e   Sdown= %.10e\n",i,j,k,  HLV->PurePopTrioColls->Colls[k]->CollName,
				   //HLV->PurePopTrioColls->Colls[k]->FBS->S_up, HLV->PurePopTrioColls->Colls[k]->FBS->S_down );
		}
	}
	
}



/* this assumes that data have already been input, etc, and that HLV->smax_to_use is already set.
 Typically that will get set in SelectAnSmax3 */
void ComputeCrossPopTrioColls(pbt_high_level_vars *HLV) 
{
	int i,j,k;
	
	HLV->CrossPopTrioColls = CreateCompleteFB_Vars(-1,-1,8, HLV->smax_to_use,HLV->PBUO->d_from_big_smax, HLV->PFR,  "CrossPop_Input","CrossPop_Output.txt");
	
	//printf("Computed PurePopTrioColls:\n");
	for(i=0;i<HLV->PFR->NumPops;i++)  {
		for(j=0;j<HLV->PFR->NumPops;j++) {
			k = CPTC_IDX(HLV->PFR->NumPops, i, j);
			//printf("\t\tParentPopIdx= %d  KidPedIdx= %d  Num= %d  Name= %s   Sup= %.10e   Sdown= %.10e\n",i,j,k,  HLV->CrossPopTrioColls->Colls[k]->CollName,
			//	   HLV->CrossPopTrioColls->Colls[k]->FBS->S_up, HLV->CrossPopTrioColls->Colls[k]->FBS->S_down );
		}
	}
}


/*
 For all the offspring in OffColl OC_idx:
 Find the parents that have fewer than the use_this_smax number of incompatibilities
 with single parents from all the SexEnum categories (MALE/FEMALE/SEX_UNKNOWN)
 

*/
void AssignMatchingSingleParents(pbt_high_level_vars *HLV, int OC_idx)
{
	int i,sex;
	pfr_offspring *Offs = HLV->PFR->Offs[OC_idx];
	int NumOffs = HLV->PFR->NumInOffColls[OC_idx];
	
	
	for(i=0;i<NumOffs;i++) {  /* cycle over all the individuals */
	
		for(sex=0;sex<3;sex++)  {  /* cycle over the different sexes of the parents */
	
			Offs[i].MendComps->Parents[sex] =  ReturnMatchingParents(
									   &(Offs[i]),			/* all the info about the youth */
									   0,				/* which genotype (from various regenotypings) to use for offspring.  Typically will be zero */
									   0,				/* which genotype to use for the parent.  Typically zero */
									   HLV,				/* all the other info we might need */
									   sex,				/* the sex of the parents we are looking for (MALE/FEMALE/SEX_UNKNOWN) */
									   &(Offs[i].MendComps->NumParents[sex]),					/* output parameter---number of parents in returned pfr_parent * array. */
									   &(Offs[i].MendComps->Num_MI_Parents[sex])						/* another output parameter, for the number of mendelian incompatibilities at each parent in the array. 
	 we have to pass this in as an int ** because we will allocate memory to it in this function. */
										   );
			
			#ifdef VERBOSE_SINGLE_PARENT_COMPAT_WITH_OFFSPRING
				printf("Offspring= %s  Sex= %d  NumberCompatibleParentsOfSex= %d ",Offs[i].Name, sex, Offs[i].MendComps->NumParents[sex]);
				//for(j=0;j<Offs[i].MendComps->NumParents[sex];j++) {
				//	printf(" %s [%d]   ",Offs[i].MendComps->Parents[sex][j]->Name,Offs[i].MendComps->Num_MI_Parents[sex][j]);
				//}
				printf("\n");
			#endif
		}
		if(i%100==0) {
			if(i==0) printf("0");
			else printf("\n%d",i);
		}
		else {
			printf(".");
		}
		fflush(stdout);
	}
	printf("\n");
}




/*
 For all the offspring in OffColl OC_idx:
 Find the parent pairs that have fewer than MaxAllowable number of incompatibilities
 with each offspring.
 */
void AssignMatchingParentPairs(pbt_high_level_vars *HLV, int OC_idx, int MaxAllowable)
{
	int i;
	pfr_offspring *Offs = HLV->PFR->Offs[OC_idx];
	int NumOffs = HLV->PFR->NumInOffColls[OC_idx];
	
	
	for(i=0;i<NumOffs;i++) {  /* cycle over all the individuals */
		
		Offs[i].MendComps->ParPairs = ReturnMatchingTrios(
														   &(Offs[i]),						/* pointer to the offspring */
														   HLV,				/* pointer to everything */
														   MaxAllowable,							/* max allowable number of mendelian incompatibilities at the trio */
														   0,							/* re-genotyping numbers for ma, pa, kid---typically 0 */
														   0,
														   0,
														   &(Offs[i].MendComps->NumParPairs)							/* output parameters telling us how many trios were compatible */
														  );
		if(i%100==0) {
			if(i==0) printf("0");
			else printf("\n%d",i);
		}
		else {
			printf(".");
		}
		fflush(stdout);
		
	}
	
	printf("\n");
}




/* For al the offspring in OffColl OC_idx:
   If there are non-excluded parent pairs for the offspring, then compute the 
	posterior probabilities that they form a trio of the different pedigree
	types.  When this is done, sort them on the basis of the posterior of the 
	parentage relationship 
*/
void CalculateTrioPosteriorProbs(pbt_high_level_vars *HLV, int OC_idx)
{
	int i,j,k;
	pfr_offspring *Offs = HLV->PFR->Offs[OC_idx];
	int NumOffs = HLV->PFR->NumInOffColls[OC_idx];
	int the_pop;
	int gotta_toss;
	
	double *FixedPri;  /* Here is what I used for a few simulations once: {1.212474e-04 ,6.272982e-05 ,4.857895e-05 ,4.293684e-03 ,4.003158e-03 ,7.843700e-05 ,9.025607e-03 ,9.207556e-03 ,9.719246e-01 ,0.000000e+00 ,0.000000e+00 ,0.000000e+00 ,0.000000e+00 ,1.538702e-06 ,1.408175e-07 ,8.823070e-04 ,3.363333e-04 ,2.608895e-06};*/

	
	
	
	for(i=0;i<NumOffs;i++) {  /* cycle over all the individuals */
		if(Offs[i].MendComps->NumParPairs>0) {
			for(j=0;j<Offs[i].MendComps->NumParPairs;j++)  {   int MaxPostRelat; double MaxPost;
				/* get the pop of the pair from the male parent */
				the_pop = Offs[i].MendComps->ParPairs[j].pop;
				
				if(HLV->PBUO->Pi[the_pop]==NULL) {
					fprintf(stderr,"ERROR!  We need a Pi prior vector for population %d (counting from 1).  Problem discovered in CalculateTrioPosteriorProbs\n",the_pop+1); 
					exit(1);
				}
				else {	
					FixedPri = HLV->PBUO->Pi[the_pop];
				}
				
				/* right here I can add a call to a function to compute the Pi priors for a given trio conditional on the 
				   spawner year of the parents and an assumed age distribution of the offspring.  But oh my that will be a mess!  And I don't have
				 time to worry about that at this point.  So, I am just going to use an average figure for each population. */
				Offs[i].MendComps->ParPairs[j].TrioPosteriors = TrioPostProbs(Offs[i].MendComps->ParPairs[j].A_64_idx, HLV->PurePopTrioColls, 0, FixedPri,
																			  the_pop, &(Offs[i].MendComps->ParPairs[j].LogL), NULL,
																			  &(Offs[i].MendComps->ParPairs[j].PrPar), &(Offs[i].MendComps->ParPairs[j].PrNonPar));
				
				/* let's record the maximum posterior and the relationship class associated with it */
				MaxPostRelat=0; 
				MaxPost=Offs[i].MendComps->ParPairs[j].TrioPosteriors[0];
				for(k=1;k<NUM_SPEC_PEDS;k++)  {
					if(Offs[i].MendComps->ParPairs[j].TrioPosteriors[k] > MaxPost) {
						MaxPost = Offs[i].MendComps->ParPairs[j].TrioPosteriors[k];
						MaxPostRelat = k;
					}
				}
				Offs[i].MendComps->ParPairs[j].MaxPosterior = MaxPost;
				Offs[i].MendComps->ParPairs[j].MaxPosteriorRelat = MaxPostRelat;
				
			
			}
			/* now sort that array if there are 2 or more ParentPairs nonexcluded */
			if(Offs[i].MendComps->NumParPairs>1) {
				qsort( (void *)(Offs[i].MendComps->ParPairs), (size_t)(Offs[i].MendComps->NumParPairs), sizeof(ParentPairMatch),cmp_parent_pair);
			}
			
			/* and now we can do the thresholding of things.  We cycle these things in order until we find a LogL less than a certain amount
			 (if we are enforcing that).
			 */
			Offs[i].MendComps->NumParPairsMendelian = Offs[i].MendComps->NumParPairs; 
			if(HLV->PBUO->EnforceParPairMinLogL != 0) {
			  gotta_toss = 0;
			  for(j=0;j<Offs[i].MendComps->NumParPairs;j++)  {
			    if(Offs[i].MendComps->ParPairs[j].LogL < HLV->PBUO->ParPairMinLogL) {
			      gotta_toss = 1;
			      break;
			    }
			  }
			  if(gotta_toss) {
			    Offs[i].MendComps->NumParPairs = j;
			  }
			}
			Offs[i].MendComps->NumParPairsMend_LogL = Offs[i].MendComps->NumParPairs;
			
			/* and now we will enforce the Max number of parent pairs criterion, too */
      if(HLV->PBUO->MaxNumRetainedParPairs > 0 && Offs[i].MendComps->NumParPairs > HLV->PBUO->MaxNumRetainedParPairs) {
        Offs[i].MendComps->NumParPairs = HLV->PBUO->MaxNumRetainedParPairs;
      }
      Offs[i].MendComps->NumParPairsMend_LogL_and_Rank = Offs[i].MendComps->NumParPairs;
			
			
			
			/* now, let's print that stuff to file too */
			for(j=0;j<Offs[i].MendComps->NumParPairs;j++)  {
				fprintf(HLV->TrioPosteriors_File,"%s\t%s\t%s\t%s\t%d\t%f",HLV->PFR->OffCollNames[OC_idx],Offs[i].Name,Offs[i].MendComps->ParPairs[j].mamapapa[MALE]->Name, Offs[i].MendComps->ParPairs[j].mamapapa[FEMALE]->Name, j+1, Offs[i].MendComps->ParPairs[j].LogL);
				for(k=0;k<NUM_SPEC_PEDS;k++)  {
					fprintf(HLV->TrioPosteriors_File,"\t%.4f", Offs[i].MendComps->ParPairs[j].TrioPosteriors[k]);
				}
				fprintf(HLV->TrioPosteriors_File,"\t%d\t%d\t%d",Offs[i].NumMissingLoci[0],
						Offs[i].MendComps->ParPairs[0].mamapapa[MALE]->NumMissingLoci[0],
						Offs[i].MendComps->ParPairs[0].mamapapa[FEMALE]->NumMissingLoci[0] );
				fprintf(HLV->TrioPosteriors_File,"\t%d\t%d\t%d\n",Offs[i].MendComps->ParPairs[j].pa_MI,Offs[i].MendComps->ParPairs[j].ma_MI,Offs[i].MendComps->ParPairs[j].trio_MI );
			}
		}
		if(i%100==0) {
			if(i==0) printf("0");
			else printf("\n%d",i);
		}
		else {
			printf(".");
		}
		fflush(stdout);
	}
	printf("\n");
}




/* For al the offspring in OffColl OC_idx:
 If there are non-excluded parent pairs for the offspring, then do the backward simulation
 of genotypes from either parental of non parental pops and record the 
 false positive and false negative rates given the LogL of the highest posterior
 prob pair.
 
 CPFB is some allocated memory to use for some calculations. 
 
 
 */
void AssessFPR_and_FNR_ByBackwardSimulation(pbt_high_level_vars *HLV, int OC_idx, FB_Vars *CPFB)
{
	int i,j,k,p,r;
	pfr_offspring *Offs = HLV->PFR->Offs[OC_idx];
	int NumOffs = HLV->PFR->NumInOffColls[OC_idx];
	int **Yobs;
	trio_genotype *GT;
	RunPars *RP = HLV->PurePopTrioColls->RP;
	int L = RP->L;
	int NY = HLV->PurePopTrioColls->RP->NY;
	double *FixedPri;  /* Here is what I used for some spip simulation tests: {1.212474e-04 ,6.272982e-05 ,4.857895e-05 ,4.293684e-03 ,4.003158e-03 ,7.843700e-05 ,9.025607e-03 ,9.207556e-03 ,9.719246e-01 ,0.000000e+00 ,0.000000e+00 ,0.000000e+00 ,0.000000e+00 ,1.538702e-06 ,1.408175e-07 ,8.823070e-04 ,3.363333e-04 ,2.608895e-06}; */
	int ThePop;  /* the population of the maximum a posteriori putative parent pair */
	int *vdex;  /* just going to have an internal one here, rather than using the ones in the various FB_Vars.  */
	double	CondPi[NUM_SPEC_PEDS],  /* this is for the probability of all the categories given that they are in s-down and given the offspring genotype */ 
			CondPiNonParental[NUM_SPEC_PEDS], /* this is for the prob of the categories given they are in S-down, the offspring geno, and that they are not Parental */
			normo;
	char *ma_geno,*pa_geno,*kid_geno;
	int A64[MAX_LOCI];
	int **AtoA_64 = RP->AtoA_64;
	double *dummy=(double *)calloc(NUM_SPEC_PEDS,sizeof(double));
	double NumExceed;
	pbt_simmed_values_struct *SimmedVals;  /* for storing and sorting simulated values */
	double NonParNormoC;  /* For storing the sum of pi * conditional_pi over all non-parental categories---useful for importance sampling.  */
	int imp_sampling=1;  /* a flag that says whether we are importance sampling or not */
	double imp_weight;
	
	int NumMC_Reps=1000;
	double **TopTwoLogls;  /* to store the top two logls from the non-parental simulated pairs on each rep */
	double *ParentalLogls;  /* to store the logls of the pairs simulated under the hypothesis that they are parental (and in S-down) */
	double *NonParentalLogls; /* to store the logls of the max-lod pairs simulated under the non-parental (but in S-down) hypothesis */
	double ParLogLExceedsObs;
	double ParLogSmallerThanNonPar;
	
	ParentalLogls = (double *)calloc(NumMC_Reps,sizeof(double));
	NonParentalLogls = (double *)calloc(NumMC_Reps,sizeof(double));
	TopTwoLogls = (double **)calloc(NumMC_Reps,sizeof(double *));
	
	for(i=0;i<NumMC_Reps;i++)  {
		TopTwoLogls[i] = (double *)calloc(2,sizeof(double));
	}
	
	
	
		
	
	/* some memory allocation */
	SimmedVals = (pbt_simmed_values_struct *)calloc(MAX_TEMP_PAIRS,sizeof(pbt_simmed_values_struct));
	vdex = (int *)calloc( ReturnVDex(RP->d, RP->smax,RP->smax)+1,sizeof(int));
	/* allocate to a Yobs */
	Yobs = (int **)calloc(HLV->PurePopTrioColls->RP->L,sizeof(int *));
	for(i=0;i<HLV->PurePopTrioColls->RP->L;i++)  {
		Yobs[i] = (int *)calloc(NY,sizeof(int));
	}
	/* allocate to a genotype struct */
	GT = (trio_genotype *)malloc(sizeof(trio_genotype));
	AllocTrioGenotypeInnards(GT,HLV->PurePopTrioColls->RP->L);
	
	
	
	for(i=0;i<NumOffs;i++) {  /* cycle over all the individuals */
		if(Offs[i].MendComps->NumParPairs>0) {
			
			/* ThePop is the population of the max posterior parent pair.  We will just do our simulations assuming they are all from this pop.  */
			ThePop = Offs[i].MendComps->ParPairs[0].pop;
			
			if(HLV->PBUO->Pi[ThePop]==NULL)  {
				fprintf(stderr,"ERROR! User needs to specify Pi for population %d.  Problem found in AssessFPR_and_FNR_ByBackwardSimulation\n",ThePop+1);
				exit(1);
			}
			else {
				FixedPri = HLV->PBUO->Pi[ThePop];
			}
			
			
			
			/* first we fill Yobs so as to condition on the offspring of the youth */
			for(j=0;j<L;j++)  {
				/* set everything to -1's first */
				for(k=0;k<NY;k++)  {
					Yobs[j][k] = -1;
				}
				/* then, if the offspring is observed at this locus put his Y value into. Just go ahead and use the first re-genotyping...  */
				if(Offs[i].geno[0][j]<3) {
					Yobs[j][0] = Offs[i].geno[0][j];
				}
			}
			
			
			/* now we cycle over all the different pedigree categories and condition on the offspring, prepare for and then do the forward step
			   and then prep for the backward step */
			for(p=0;p<NUM_SPEC_PEDS;p++)  { int PedPick;
				
				/* figure out the index of the appropriate collection */
				PedPick = PPTC_IDX(NUM_SPEC_PEDS,p,ThePop);
				
				/* now condition on that Yobs for pedigree p.  Put the result in the Aprobs for the appropriate CPFB Collection */
				ConditionAprobsOnYobs(HLV->PurePopTrioColls->Colls[PedPick]->Aprobs, Yobs, HLV->PurePopTrioColls->RP, CPFB->Colls[p]->Aprobs);
				
				
				if(imp_sampling==1)  {  int ll; /* to compute the importance weights, we have to fill out the Aprobs_M array on this CPFB thing */
					for(ll=0;ll<L;ll++)  {
						ReturnAProbsM_HardWiredForTrios(RP, CPFB->Colls[p]->Aprobs[ll], CPFB->Colls[p]->Aprobs_M[ll]);
					}
				}
				
				
				/* now compute the X-probs in there given those conditioned Aprobs */
				XStuffFromAStuff(CPFB->Colls[p]->Aprobs, RP, CPFB->Colls[p]);
				
				/* then, get ready to do a forward step by initializing vdex to -1's*/
				for(j=0;j<=ReturnVDex(RP->d, RP->smax,RP->smax);j++) {
					vdex[j]=-1;
				}
				
				/* then do this forward step and prepare for the backward steps  */
				CPFB->Colls[p]->FBS = ForwardStep(RP, CPFB->Colls[p]->Xprobs, CPFB->Colls[p]->Xfakeprobs, vdex, 0);
				PrepForBackwardStep(CPFB->Colls[p]->FBS, RP, CPFB->Colls[p]->Xprobs);
				
								
			}
			
			/* now we compute the expected fraction of pedigree types conditional on the 
			 offspring genotype and on being in S-down */
			normo = 0.0;
			NonParNormoC = 0.0;
			for(j=0;j<NUM_SPEC_PEDS;j++)  {
				CondPi[j] = FixedPri[j] *  CPFB->Colls[j]->FBS->S_down;
				normo += CondPi[j];
				if(j>0) {
					NonParNormoC += CondPi[j];
				}
			}
			CondPiNonParental[0] = 0.0;  
			for(j=0;j<NUM_SPEC_PEDS;j++)  {
				CondPi[j]/=normo;
				if(j>0) {
					CondPiNonParental[j] = CondPi[j] / (1.0 - CondPi[0]);
				}
			}
			
			//printf("DONE_WITH_FORWARD_STEPS_FOR_OFFSPRING_AND_PREPPED_FOR_BACKWARD  %s\n",Offs[i].Name);
//			for(j=0;j<NUM_SPEC_PEDS;j++)  {
//				printf("%f   %f    %f\n",CondPi[j],CondPiNonParental[j], FixedPri[j]);
//			}
			
			
			/* now we just cycle over the reps and for each one, cycle over the number of non excluded parent pairs, and for each one
			 draw a non-parental trio type and simulate from it.  */
			for(imp_sampling=0;imp_sampling<1;imp_sampling++) {
				if(imp_sampling==0) {
					imp_weight = 1.0;
				}
				NumExceed = 0.0;
				for(r=0;r<NumMC_Reps;r++)  { int ped;  double LogL;  int pp; double PrPar; double PrNonPar;
					
					for(pp=0;pp<Offs[i].MendComps->NumParPairs;pp++)  { int  tt;
						
						
						if(imp_sampling) {
							ped = 0;
						}
						else {
							ped = IntFromProbsRV(CondPiNonParental, 0, NUM_SPEC_PEDS);
						}
						SimBackward(CPFB->Colls[ped]->FBS, CPFB->Colls[ped]->AprobsGX,RP,GT);
						
						/* now apply the missing data pattern in a parent pair (cycling over them successively via pp) to the
						   simulated trio genotype and convert them all into an A64 vector. */
						kid_geno = Offs[i].geno[0];
						pa_geno = Offs[i].MendComps->ParPairs[pp].mamapapa[MALE]->geno[0];
						ma_geno = Offs[i].MendComps->ParPairs[pp].mamapapa[FEMALE]->geno[0];
						for(j=0;j<RP->L;j++)  { int m;
							m = RP->MissModeFrom0sAnd1s	[ kid_geno[j]==3 ][ pa_geno[j]==3 ][ ma_geno[j]==3 ];
							A64[j] = AtoA_64[m][ GT->a[j] ];
						}
						
						/* and now we compute the posterior probabilities of that genotype.  Note that PrPar and PrNonPar are included in the call
						 but we don't use those values, here.  We grab them later if we are doing some importance sampling.  */
						TrioPostProbs(A64,HLV->PurePopTrioColls, 0, FixedPri, ThePop, &LogL, dummy,&PrPar,&PrNonPar);
						/*printf("\tr= %d   post= %f   simpost= %f   LogL= %f   ped=  %d    Exceeds= %d\n",r,Offs[i].MendComps->ParPairs[0].TrioPosteriors[0],
							   dummy[0],LogL,ped, 
							   dummy[0]>Offs[i].MendComps->ParPairs[0].TrioPosteriors[0]); */
						
						/* and then transfer that information into the simmed_vals array for sorting */
						SimmedVals[pp].LogL = LogL;
						SimmedVals[pp].ParPosterior = dummy[0];
						for(tt=0;tt<NUM_SPEC_PEDS;tt++)  {
							SimmedVals[pp].PostVec[tt] = dummy[tt];
						}
						SimmedVals[pp].ped = ped;
						SimmedVals[pp].parpair = pp;
											
						
						
						/* here, if we are doing the importance sampling thing we also have to compute the prob of having simulated that
						   genotype so we do TrioPostProbs given the conditioned FBS we are simulating from.  We pass 0 for the third arg to 
							announce that. */
						if(imp_sampling) {
							TrioPostProbs(A64,CPFB, 1, FixedPri, ThePop, &LogL, dummy,&PrPar,&PrNonPar);
							SimmedVals[pp].PrPar = PrPar;
							SimmedVals[pp].PrNonPar = PrNonPar;					
						}
					}
					
					/* now sort that SimmedVals array if there are 2 or more ParentPairs nonexcluded */
					if(Offs[i].MendComps->NumParPairs>1) {
						qsort( (void *)(SimmedVals), (size_t)(Offs[i].MendComps->NumParPairs), sizeof(pbt_simmed_values_struct),cmp_simmed_values);
	/*					printf("SORTED_SIMMED_VALS: ");
						for(pp=0;pp<Offs[i].MendComps->NumParPairs;pp++)  {
							printf("%d  [%d %d]  %f  %f     | ",pp,SimmedVals[pp].ped,SimmedVals[pp].parpair,SimmedVals[pp].LogL,SimmedVals[pp].ParPosterior);
						}
						printf("\n");
	*/
					}
					
					if(imp_sampling)  {  /* here we compute the importance weight */
						imp_weight = (SimmedVals[0].PrNonPar / NonParNormoC)  /  SimmedVals[0].PrPar;  /* this is not quite correct yet.  I still have to adjust PrPar to account for being in Sdown */
						//printf("\tIMP_WEIGHT: %s  r= %d   %e\n", Offs[i].Name,r, (SimmedVals[0].LogL > Offs[i].MendComps->ParPairs[0].LogL) * imp_weight); 
					}
					
					NumExceed += (SimmedVals[0].LogL > Offs[i].MendComps->ParPairs[0].LogL) * imp_weight;
					
					/* also, store the top two for use in computing prob that if the parent were there it would have higher LOD */
					TopTwoLogls[r][0] = SimmedVals[0].LogL;
					TopTwoLogls[r][1] = SimmedVals[1].LogL;
					NonParentalLogls[r] = SimmedVals[0].LogL; 
					
					
									
					//if(r>0 && r%2000==0) printf("\tBackSim rep %d done\n",r);
				}
				Offs[i].p_value=(double)NumExceed/(double)NumMC_Reps;
				
				
				/*******************************************************************************************************************************************************************************************/
				/* NOW DO THIS SIMULATIONS ASSUMING PARENTAL ********************************************************/
				for(r=0;r<NumMC_Reps;r++)  { int ped;  double LogL;  int pp; double PrPar; double PrNonPar;
					
					ped = 0;  /* always simulate from parental here */
					pp = UniformRV(0,Offs[i].MendComps->NumParPairs-1);  /* simulate the missing data pattern from a randomly-chosen non-excluded pair. */ 
					
						
						SimBackward(CPFB->Colls[ped]->FBS, CPFB->Colls[ped]->AprobsGX,RP,GT);
						
						/* now apply the missing data pattern in the pp parent pair to the
						 simulated trio genotype and convert them all into an A64 vector.  Just assume that we are using RegenoNumber 0 */
						kid_geno = Offs[i].geno[0];
						pa_geno = Offs[i].MendComps->ParPairs[pp].mamapapa[MALE]->geno[0];
						ma_geno = Offs[i].MendComps->ParPairs[pp].mamapapa[FEMALE]->geno[0];
						for(j=0;j<RP->L;j++)  { int m;
							m = RP->MissModeFrom0sAnd1s	[ kid_geno[j]==3 ][ pa_geno[j]==3 ][ ma_geno[j]==3 ];
							A64[j] = AtoA_64[m][ GT->a[j] ];
						}
						
						/* and now we compute the posterior probabilities of that genotype.  Note that PrPar and PrNonPar are included in the call
						 but we don't use those values, here.  We grab them later if we are doing some importance sampling.  */
						TrioPostProbs(A64,HLV->PurePopTrioColls, 0, FixedPri, ThePop, &LogL, dummy,&PrPar,&PrNonPar);
						
						/* and then transfer that information into the ParentalLogls array for later use. */
						ParentalLogls[r] = LogL;
				}
				
				/* OK, now we have that sample of parental pair LogL's and we compare them to the non-parental ones to determine both the FNR, and the prob that,
				 if the true parents are amongst the non-excluded, they will have a LogL higher than any of the non-parental ones.  (I THINK THIS MIGHT BE BOGUS) */
				{
					double TotComparisons;
					int t;
					double pl;
					double NP;
					
					ParLogLExceedsObs = 0.0;
					ParLogSmallerThanNonPar = 0.0;
					TotComparisons = 0.0;
					
					NP = Offs[i].MendComps->NumParPairs;
					if(NP<1) {
						ParLogSmallerThanNonPar = 0.0;
					}
					for(r=0;r<NumMC_Reps;r++)  {
						pl = ParentalLogls[r];
						if(pl > Offs[i].MendComps->ParPairs[0].LogL) {  ParLogLExceedsObs += 1.0; }
						if(NP>1) {
							for(t=0;t<NumMC_Reps;t++)  {
								/* note:  1.0-(1.0/NP) is the prob that the highest nonparental Logl is still in there, even if one of the pairs was replaced by a parental */
								ParLogSmallerThanNonPar += (  (1.0-(1.0/NP)) * (pl < TopTwoLogls[t][0])  +  (1.0/NP) * (pl < TopTwoLogls[t][1]));
								TotComparisons += 1.0;
							}
						}
						
					}
					ParLogLExceedsObs /= (double)NumMC_Reps;
					if(NP>1) {  ParLogSmallerThanNonPar /= TotComparisons;  }
					
				}
													 
													 
				/* Here is a different block of code in which we provide for each individuals vector of LogL values corresponding to different p-values */
				{
					Offs[i].PValueParentalCDF = EmpiricalCrossCDF100(NonParentalLogls, ParentalLogls, NumMC_Reps, NumMC_Reps);
				}
													
				
				/*************************************************************************************************************************************************/
								
				/*fprintf(HLV->MaxLodNonExpPar_File,"%s\t%s\t%s\t%.6f\t%f\t%f\t%f\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t",
						Offs[i].Name,
						Offs[i].MendComps->ParPairs[0].mamapapa[MALE]->Name,
						Offs[i].MendComps->ParPairs[0].mamapapa[FEMALE]->Name,
						(double)NumExceed/(double)NumMC_Reps,
						Offs[i].MendComps->ParPairs[0].LogL,
						Offs[i].MendComps->ParPairs[0].TrioPosteriors[0],
						Offs[i].MendComps->ParPairs[0].MaxPosterior,
						SpecPedIdx2Str(Offs[i].MendComps->ParPairs[0].MaxPosteriorRelat),
						Offs[i].MendComps->NumParents[MALE],
						Offs[i].MendComps->NumParents[FEMALE],
						Offs[i].MendComps->NumParents[SEX_UNKNOWN],
						Offs[i].MendComps->NumParPairs,
						Offs[i].NumMissingLoci[0],
						Offs[i].MendComps->ParPairs[0].mamapapa[MALE]->NumMissingLoci[0],
						Offs[i].MendComps->ParPairs[0].mamapapa[FEMALE]->NumMissingLoci[0],
						Offs[i].MendComps->ParPairs[0].pa_MI,
						Offs[i].MendComps->ParPairs[0].ma_MI,
						Offs[i].MendComps->ParPairs[0].trio_MI
					   );
				 */
				/* now provide a comma separated list of locus indexes that have Mendelian incompatibilities in the trio. */
				/*p=0;
				for(j=0;j<L;j++)  {
					if(Offs[i].MendComps->ParPairs[0].X_states_64[j] > 0) {
						if(p>0) {
							fprintf(HLV->MaxLodNonExpPar_File,",");
						}
						p++;
						fprintf(HLV->MaxLodNonExpPar_File,"%d",j+1);
					}
				}
				if(p==0) { fprintf(HLV->MaxLodNonExpPar_File,"NA"); }
				fprintf(HLV->MaxLodNonExpPar_File,"\n");
				*/
				
				
				
				
				/* only do the importance sampling if there were zero exceedances using the vanilla method.
				 Note that this might bias things upward a little bit and may not be what I really want to do. */
				//if( imp_sampling==0 && NumExceed>0  ) { 
				//	imp_sampling++;
				//}  			
			
			}
			
			
			/* down here we deallocate all the memory we allocated for those forward steps */
			for(p=0;p<NUM_SPEC_PEDS;p++)  {
				FreeFB_Vars_innards(RP, CPFB->Colls[p]->FBS);
				free(CPFB->Colls[p]->FBS->FinalStepNormalizedProbs);  /* this is another thing that gets put in there when prepping for the backward step */
				free(CPFB->Colls[p]->FBS);
			}
		}
		/* here we still are going to want to report them in the output, even if they have no non-excluded parent pairs */
		else {
			/*fprintf(HLV->MaxLodNonExpPar_File,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t\t\n",
					Offs[i].Name,
					"---",
					"---",
					"---",
					"---",
					"---",
					"---",
					"---",
					Offs[i].MendComps->NumParents[MALE],
					Offs[i].MendComps->NumParents[FEMALE],
					Offs[i].MendComps->NumParents[SEX_UNKNOWN],
					Offs[i].MendComps->NumParPairs,
					Offs[i].NumMissingLoci[0],
					"---",
					"---",
					"---",
					"---",
					"---"
				   );
			 */
			;
			
		}
		if(i%100==0) {
			if(i==0) printf("0"); 
			else printf("\n%d",i);
		}
		else {
			printf(".");
		}
		fflush(stdout);
	}
	
	
	/* down here at the end we will deallocate some memory */
	for(i=0;i<HLV->PurePopTrioColls->RP->L;i++)  {
		free(Yobs[i]);
	}
	free(Yobs);
	free(vdex);
	FreeTrioGenotypeInnards(GT,L);
	free(GT);
	free(SimmedVals);
	for(i=0;i<NumMC_Reps;i++)  {
		free(TopTwoLogls[i]);
	}
	free(TopTwoLogls);
	free(ParentalLogls);
	free(NonParentalLogls);
	
	
	printf("\n");
	
}


/* record info in the inds_with_max_lod_parents_from_this_pop for population Pop */
inds_with_max_lod_parents_from_this_pop *RecordNonExcParentPairsFromPop(int Pop, pbt_high_level_vars *HLV)
{
	inds_with_max_lod_parents_from_this_pop *ret = (inds_with_max_lod_parents_from_this_pop *)malloc(sizeof(inds_with_max_lod_parents_from_this_pop));
	int i,j,n;
	int *NumNonExc = (int *)calloc(MAX_NUM_NONEXC_PARENT_PAIRS,sizeof(int));
	int *FinalIdx = (int *)calloc(MAX_NUM_NONEXC_PARENT_PAIRS,sizeof(int));
	int temp_non_exc, temp_idx;
	pfr_offspring *off;
	
	
	ret->n=0;
	
	/* first cycle over all the offspring in all the collections and count up what we have here */
	for(i=0;i<HLV->PFR->NumOffColls;i++)  {
		for(j=0;j<HLV->PFR->NumInOffColls[i];j++) {
			off = &(HLV->PFR->Offs[i][j]);
			if(off->MendComps->NumParPairs>0 && off->MendComps->ParPairs[0].pop == Pop) {
				temp_non_exc = off->MendComps->NumParPairs;
				NumNonExc[temp_non_exc]++;
				ret->n++;
			}
		}
	}
	
	
	/* now, see how many classes that fits everyone into */
	ret->NumNonExcCnts=0;
	for(i=0;i<MAX_NUM_NONEXC_PARENT_PAIRS;i++)  {
		if(NumNonExc[i]>0) {
			ret->NumNonExcCnts++;
		}
	}
	/* then allocate space for those classes */
	ret->NonExcCnts = (int *)calloc(ret->NumNonExcCnts,sizeof(int));
	ret->NumKidsWithNonExc = (int *)calloc(ret->NumNonExcCnts,sizeof(int));
	
	/* then cycle over the possible values of number of non-excluded pairs and put those all in there */
	n=0;
	for(i=0;i<MAX_NUM_NONEXC_PARENT_PAIRS;i++)  {
		if(NumNonExc[i]>0) {
			ret->NonExcCnts[n] = i;
			ret->NumKidsWithNonExc[n] = NumNonExc[i];
			FinalIdx[i] = n;
			n++;
		}
	}
	
	
	/* now allocate memory to the pointers to the kids */
	ret->inds = (pfr_offspring **)calloc(ret->n,sizeof(pfr_offspring *));
	ret->inds_by_cnt = (pfr_offspring ***)calloc(ret->NumNonExcCnts,sizeof(pfr_offspring **));
	for(i=0;i<ret->NumNonExcCnts;i++)  {
		ret->inds_by_cnt[i] = (pfr_offspring **)calloc(ret->NumKidsWithNonExc[i],sizeof(pfr_offspring *));
	}
	
	
	
	/* allocate memory to keep track of how often a randomly simulated parental trio has lower LOD than one of the max simulated
	 ones from a non-parental trio */
	/*ret->BetasAtEachNonExcCategory = DvalVector(0, ret->NumNonExcCnts-1, 0.0, 0.0, -1.0);*/
	
	ret->BetaArraysByHundredths = (double **)calloc(ret->NumNonExcCnts,sizeof(int *));
	for(i=0;i<ret->NumNonExcCnts;i++)  {
		ret->BetaArraysByHundredths[i] = (double *)calloc(100,sizeof(int));
	}
	
	/* and now we are ready to cycle over all the individuals in the offspring collection again and set pointers where they
	   need to be, etc */
	/* first initialize stuff to collect a sum */
	for(i=0;i<MAX_NUM_NONEXC_PARENT_PAIRS;i++)  {
		NumNonExc[i] = 0;
	}
	/* then have at it */
	n=0;
	for(i=0;i<HLV->PFR->NumOffColls;i++)  {
		for(j=0;j<HLV->PFR->NumInOffColls[i];j++) {
			off = &(HLV->PFR->Offs[i][j]);
			if(off->MendComps->NumParPairs>0 && off->MendComps->ParPairs[0].pop == Pop) {
				temp_non_exc = off->MendComps->NumParPairs;
				temp_idx = FinalIdx[temp_non_exc];
				
				/* now we can make the necessary assignments */
				ret->inds[n] = off;
				ret->inds_by_cnt[temp_idx][NumNonExc[temp_idx]++] = off;
				n++;
			}
		}
	}
										 
	
	free(NumNonExc);
	free(FinalIdx);
	
#ifdef VERBOSE_NUM_NON_EXC_PAIRS_BY_POP
	/* here we can print out a little report about this */
	printf("NUM_NON_EXC_PAIRS_BY_POP:  Pop= %d  Name= %s\n",Pop, HLV->PFR->PopNames[Pop]);
	printf("NUM_NON_EXC_PAIRS_BY_POP:  NumNonExcCats= %d  NumNonExcOffs=%d \n",ret->NumNonExcCnts,ret->n);
	for(i=0;i<ret->NumNonExcCnts;i++)  {
		printf("NUM_NON_EXC_PAIRS_BY_POP:  Idx= %d  NumNonExc= %d   NumKidsWithNonExc= %d\n",i, ret->NonExcCnts[i], ret->NumKidsWithNonExc[i]);
		for(j=0;j<ret->NumKidsWithNonExc[i];j++)  {
			printf("NUM_NON_EXC_PAIRS_BY_POP:              Idx= %d    Ind= %s    NumNonExc= %d \n",j,ret->inds_by_cnt[i][j]->Name, ret->inds_by_cnt[i][j]->MendComps->NumParPairs);
		}
	}
#endif
	
	return(ret);
	
}






/* this was copied from the Conditional simulation function and then highly modified. 
 For a particular population with index Pop
 If there are non-excluded parent pairs for the offspring, then do the backward simulation
 of genotypes from either parental of non parental pops and record the 
 false positive and false negative rates given the LogL of the highest posterior
 prob pair.
 
 CPFB is some allocated memory to use for some calculations. 
 
 
 
*/
void MarginalSimForFDR_Calcs(pbt_high_level_vars *HLV, int ThePop, FB_Vars *CPFB)
{
	int i,j,k,p,r;
	int **Yobs;
	trio_genotype *GT;
	RunPars *RP = HLV->PurePopTrioColls->RP;
	int L = RP->L;
	int NY = HLV->PurePopTrioColls->RP->NY;
	double *FixedPri;  /* Here is what I used for some spip simulation tests: {1.212474e-04 ,6.272982e-05 ,4.857895e-05 ,4.293684e-03 ,4.003158e-03 ,7.843700e-05 ,9.025607e-03 ,9.207556e-03 ,9.719246e-01 ,0.000000e+00 ,0.000000e+00 ,0.000000e+00 ,0.000000e+00 ,1.538702e-06 ,1.408175e-07 ,8.823070e-04 ,3.363333e-04 ,2.608895e-06}; */
	int *vdex;  /* just going to have an internal one here, rather than using the ones in the various FB_Vars.  */
	double	CondPi[NUM_SPEC_PEDS],  /* this is for the probability of all the categories given that they are in s-down and given the offspring genotype */ 
	CondPiNonParental[NUM_SPEC_PEDS], /* this is for the prob of the categories given they are in S-down, the offspring geno, and that they are not Parental */
	normo;
	char *ma_geno,*pa_geno,*kid_geno;
	int A64[MAX_LOCI];
	int **AtoA_64 = RP->AtoA_64;
	double *dummy=(double *)calloc(NUM_SPEC_PEDS,sizeof(double));
	double NonParNormoC;  /* For storing the sum of pi * conditional_pi over all non-parental categories---useful for importance sampling.  */
	
	int NumMC_Reps=1000;
	double **TopTwoLogls;  /* to store the top two logls from the non-parental simulated pairs on each rep */
	double *ParentalPostPs;  /* to store the posterior probs of the pairs simulated under the hypothesis that they are parental (and in S-down) */
	double **SimmedMaxPostPsNonPar; /* to store the max posterior probs of being parental under the non-parental hypothesis.  This is indexed by the number of non-excluded pairs */
	
	
	SimmedMaxPostPsNonPar = (double **)calloc(HLV->KidsWithMaxLOD_Parents[ThePop]->NumNonExcCnts,sizeof(double *));
	for(i=0;i<HLV->KidsWithMaxLOD_Parents[ThePop]->NumNonExcCnts;i++)  {
		SimmedMaxPostPsNonPar[i] = (double *)calloc(NumMC_Reps,sizeof(double));
	}
	ParentalPostPs = (double *)calloc(NumMC_Reps,sizeof(double));
	TopTwoLogls = (double **)calloc(NumMC_Reps,sizeof(double *));
	for(i=0;i<NumMC_Reps;i++)  {
		TopTwoLogls[i] = (double *)calloc(2,sizeof(double));
	}
	
	
	if(HLV->PBUO->Pi[ThePop]==NULL)  {
		fprintf(stderr,"ERROR! User needs to specify Pi for population %d.  Problem found in MarginalSimForFDR_Calcs\n",ThePop+1);
		exit(1);
	}
	else {
		FixedPri = HLV->PBUO->Pi[ThePop];
	}
	
	
	
	/* some memory allocation */
	vdex = (int *)calloc( ReturnVDex(RP->d, RP->smax,RP->smax)+1,sizeof(int));
	/* allocate to a Yobs */
	Yobs = (int **)calloc(HLV->PurePopTrioColls->RP->L,sizeof(int *));
	for(i=0;i<HLV->PurePopTrioColls->RP->L;i++)  {
		Yobs[i] = (int *)calloc(NY,sizeof(int));
	}
	/* allocate to a genotype struct */
	GT = (trio_genotype *)malloc(sizeof(trio_genotype));
	AllocTrioGenotypeInnards(GT,HLV->PurePopTrioColls->RP->L);
	
	
	
			
			
			
	/* first we fill Yobs so as to condition on absolutely nothing! */
	for(j=0;j<L;j++)  {
		/* set everything to -1's first */
		for(k=0;k<NY;k++)  {
			Yobs[j][k] = -1;
		}
	}
			
			
	/* now we cycle over all the different pedigree categories and condition on the offspring, prepare for and then do the forward step
	 and then prep for the backward step */
	for(p=0;p<NUM_SPEC_PEDS;p++)  { int PedPick;
		
		/* figure out the index of the appropriate collection */
		PedPick = PPTC_IDX(NUM_SPEC_PEDS,p,ThePop);
		
		/* now condition on that Yobs for pedigree p.  Put the result in the Aprobs for the appropriate CPFB Collection */
		ConditionAprobsOnYobs(HLV->PurePopTrioColls->Colls[PedPick]->Aprobs, Yobs, HLV->PurePopTrioColls->RP, CPFB->Colls[p]->Aprobs);
		
			
		/* now compute the X-probs in there given those conditioned Aprobs */
		XStuffFromAStuff(CPFB->Colls[p]->Aprobs, RP, CPFB->Colls[p]);
		
		/* then, get ready to do a forward step by initializing vdex to -1's*/
		for(j=0;j<=ReturnVDex(RP->d, RP->smax,RP->smax);j++) {
			vdex[j]=-1;
		}
		
		/* then do this forward step and prepare for the backward steps  */
		CPFB->Colls[p]->FBS = ForwardStep(RP, CPFB->Colls[p]->Xprobs, CPFB->Colls[p]->Xfakeprobs, vdex, 0);
		PrepForBackwardStep(CPFB->Colls[p]->FBS, RP, CPFB->Colls[p]->Xprobs);		
	}
			
	/* now we compute the expected fraction of pedigree types conditional on the 
	 offspring genotype and on being in S-down */
	normo = 0.0;
	NonParNormoC = 0.0;
	for(j=0;j<NUM_SPEC_PEDS;j++)  {
		CondPi[j] = FixedPri[j] *  CPFB->Colls[j]->FBS->S_down;
		normo += CondPi[j];
		if(j>0) {
			NonParNormoC += CondPi[j];
		}
	}
	CondPiNonParental[0] = 0.0;  
	for(j=0;j<NUM_SPEC_PEDS;j++)  {
		CondPi[j]/=normo;
		if(j>0) {
			CondPiNonParental[j] = CondPi[j] / (1.0 - CondPi[0]);
		}
	}
	
	/* allocate memory to the MargAlpha fields in all the relevant offspring individuals and init them all to zeroes */
	for(i=0;i<HLV->KidsWithMaxLOD_Parents[ThePop]->n;i++) {
		HLV->KidsWithMaxLOD_Parents[ThePop]->inds[i]->MargAlphas = AllocIval(0, 1, -1.0);
	}
	
	
	/* before we start into this, we simulate NumMC_Reps truly parental trios and record their PosteriorProbs in a big array */
	for(r=0;r<NumMC_Reps;r++)  { int ped;  double LogL;   double PrPar; double PrNonPar;
		ped=0;
		SimBackward(CPFB->Colls[ped]->FBS, CPFB->Colls[ped]->AprobsGX,RP,GT);
		
		/* now apply the missing data pattern in a randomly chosen set of parents from the population and from a kid randomly
		 chosen from those with max-LOD, non-excluded individuals in the population */
		GrabRandomGenos(HLV, ThePop, &kid_geno, &pa_geno, &ma_geno);
		for(j=0;j<RP->L;j++)  { int m;
			m = RP->MissModeFrom0sAnd1s	[ kid_geno[j]==3 ][ pa_geno[j]==3 ][ ma_geno[j]==3 ];
			A64[j] = AtoA_64[m][ GT->a[j] ];
		}
		
		
		/* and now we compute the posterior probabilities of that genotype.  Note that PrPar and PrNonPar are included in the call
		 but we don't use those values, here.  We grab them later if we are doing some importance sampling.  */
		TrioPostProbs(A64,HLV->PurePopTrioColls, 0, FixedPri, ThePop, &LogL, dummy,&PrPar,&PrNonPar);
		
		ParentalPostPs[r] = dummy[0];
	}
	
	
			
	/* now we just cycle over the reps and within that keep cycling over the number
	 of nonexluded parents.  for each one
	 draw a non-parental trio type and simulate from it.  Keep track of the max, and occasionally, after 
	 an appropriate number of simulated trios, check to say where the observed values on the individuals
	 with that same number of non-excluded parents falls.  */
	for(r=0;r<NumMC_Reps;r++)  { int ped;  double LogL;   double PrPar; double PrNonPar; int level; double MaxPrPar;
		
		level = 0;
		MaxPrPar = -99999999999.99999;
		for(i=0;i<HLV->KidsWithMaxLOD_Parents[ThePop]->NonExcCnts[ HLV->KidsWithMaxLOD_Parents[ThePop]->NumNonExcCnts - 1 ]; i++) {
			ped = IntFromProbsRV(CondPiNonParental, 0, NUM_SPEC_PEDS);
			SimBackward(CPFB->Colls[ped]->FBS, CPFB->Colls[ped]->AprobsGX,RP,GT);
			
			/* now apply the missing data pattern in a randomly chosen set of parents from the population and from a kid randomly
			 chosen from those with max-LOD, non-excluded individuals in the population */
			GrabRandomGenos(HLV, ThePop, &kid_geno, &pa_geno, &ma_geno);
			for(j=0;j<RP->L;j++)  { int m;
				m = RP->MissModeFrom0sAnd1s	[ kid_geno[j]==3 ][ pa_geno[j]==3 ][ ma_geno[j]==3 ];
				A64[j] = AtoA_64[m][ GT->a[j] ];
			}
			
			
			/* and now we compute the posterior probabilities of that genotype.  Note that PrPar and PrNonPar are included in the call
			 but we don't use those values, here.  We grab them later if we are doing some importance sampling.  */
			TrioPostProbs(A64,HLV->PurePopTrioColls, 0, FixedPri, ThePop, &LogL, dummy,&PrPar,&PrNonPar);
			
			if(dummy[0] > MaxPrPar) {
				MaxPrPar = dummy[0];
			}
			
			/* here we do the occasional comparision.  Note, we are assuming that the parental hypothesis is Pedigree 0 (in dummy and the TrioPosterior in the ParPairs)  */
			if(i==HLV->KidsWithMaxLOD_Parents[ThePop]->NonExcCnts[level]-1) {
				SimmedMaxPostPsNonPar[level][r] = MaxPrPar;
				for(k=0;k<HLV->KidsWithMaxLOD_Parents[ThePop]->NumKidsWithNonExc[level];k++)  {
					HLV->KidsWithMaxLOD_Parents[ThePop]->inds_by_cnt[level][k]->MargAlphas->v =   (HLV->KidsWithMaxLOD_Parents[ThePop]->inds_by_cnt[level][k]->MendComps->ParPairs[0].TrioPosteriors[0] < MaxPrPar);
					IncrementIval(HLV->KidsWithMaxLOD_Parents[ThePop]->inds_by_cnt[level][k]->MargAlphas);
				}
				/*for(k=0;k<NumMC_Reps;k++)  {
					HLV->KidsWithMaxLOD_Parents[ThePop]->BetasAtEachNonExcCategory[level]->v = (ParentalPostPs[k] < MaxPrPar);
					IncrementDval(HLV->KidsWithMaxLOD_Parents[ThePop]->BetasAtEachNonExcCategory[level]);
				}
				*/
				level++;  /* move on to the next one */
			}
		}
		
		
		printf("Done with Overall Rep %d.  MaxPrPar = %f  \n",r,MaxPrPar);
	}
				

	/* now, for fun, print out the MargAlphas of these guys */
	for(i=0;i<HLV->KidsWithMaxLOD_Parents[ThePop]->n;i++)  { pfr_offspring *kiddo;
		kiddo = HLV->KidsWithMaxLOD_Parents[ThePop]->inds[i];
		printf("MARG_ALPHAS: Pop= %d  Ind= %s  Pa= %s   Ma= %s  Alpha= %f   N= %d   NumNonExc= %d\n", ThePop, kiddo->Name,
			   kiddo->MendComps->ParPairs[0].mamapapa[MALE]->Name,
			   kiddo->MendComps->ParPairs[0].mamapapa[FEMALE]->Name,
			   kiddo->MargAlphas->Ave, kiddo->MargAlphas->NumAved, kiddo->MendComps->NumParPairs);  
	}
								

	/* and here we compare the parental values at each different level */
	/* first sort those things */
	for(i=0;i<HLV->KidsWithMaxLOD_Parents[ThePop]->NumNonExcCnts;i++)  {
		qsort(SimmedMaxPostPsNonPar[i],NumMC_Reps,sizeof(double),cmp_doubles);
	}
	{	double inc;
		int x; int q;
		double summand;
		
		summand = 1.0/(double)NumMC_Reps;
		for(q=0,inc=0.0;inc<1.0;inc+=.01,q++) {
			x = (int)(NumMC_Reps * inc);
			for(r=0;r<NumMC_Reps;r++) {
				for(i=0;i<HLV->KidsWithMaxLOD_Parents[ThePop]->NumNonExcCnts;i++)  {
					if(SimmedMaxPostPsNonPar[i][x] > ParentalPostPs[r]) {
							HLV->KidsWithMaxLOD_Parents[ThePop]->BetaArraysByHundredths[i][q] +=  summand;
					}
				}
			}
			
		}
		
		/* and here we report it at the .05 levels */
		for(i=0;i<HLV->KidsWithMaxLOD_Parents[ThePop]->NumNonExcCnts;i++)  {
			printf("BETA_AT_ALPHA=.05:     NumNonExc= %d    Beta= %f\n",HLV->KidsWithMaxLOD_Parents[ThePop]->NonExcCnts[i],HLV->KidsWithMaxLOD_Parents[ThePop]->BetaArraysByHundredths[i][5]);
		}
			
	}
			
			
	/* down here we deallocate all the memory we allocated for those forward steps */
	for(p=0;p<NUM_SPEC_PEDS;p++)  {
		FreeFB_Vars_innards(RP, CPFB->Colls[p]->FBS);
		free(CPFB->Colls[p]->FBS->FinalStepNormalizedProbs);  /* this is another thing that gets put in there when prepping for the backward step */
		free(CPFB->Colls[p]->FBS);
	}
			
	
	/* down here at the end we will deallocate some memory */
	for(i=0;i<HLV->PurePopTrioColls->RP->L;i++)  {
		free(Yobs[i]);
	}
	free(Yobs);
	free(vdex);
	FreeTrioGenotypeInnards(GT,L);
	free(GT);
	for(i=0;i<NumMC_Reps;i++)  {
		free(TopTwoLogls[i]);
	}
	free(TopTwoLogls);
	free(ParentalPostPs);
	for(i=0;i<HLV->KidsWithMaxLOD_Parents[ThePop]->NumNonExcCnts;i++)  {
		free(SimmedMaxPostPsNonPar);
	}
	free(SimmedMaxPostPsNonPar);
	
}



/* this function grabs a random kid genotype from amongst the kids with max-LOD non-exluded parents
 in population ThePop, and it grabs a ma and pa from ThePop, too.  This will be used to insert missing
 data into the simulations.  kid_geno, pa_geno, and ma_geno are all output parameters.  It is going to not
 worry about drawing from the non-sexed pile of individuals, hoping that they will be few and far between.
 I could change that later as need be.
 
 Also, I just use the first re-genotyping of all of these.
 */
void GrabRandomGenos(pbt_high_level_vars *HLV, int ThePop, char **kid_geno, char **pa_geno, char **ma_geno)
{
	int padx, madx, kidx;
	
	/* start with ma and pa */
	padx = UniformRV(0,HLV->PFR->NumInPops[ThePop][MALE]-1);
	(*pa_geno) = HLV->PFR->Pops[ThePop][MALE][padx].geno[0];
	
	madx = UniformRV(0,HLV->PFR->NumInPops[ThePop][FEMALE]-1);
	(*ma_geno) = HLV->PFR->Pops[ThePop][FEMALE][madx].geno[0];
	
	/* then go ahead and grab your kid from amongst all the ones that have max-LOD non-excluded parents in ThePop */
	kidx = UniformRV(0,HLV->KidsWithMaxLOD_Parents[ThePop]->n-1);
	(*kid_geno) = HLV->KidsWithMaxLOD_Parents[ThePop]->inds[kidx]->geno[0];
}


/* compare two SimmedValue structs on the basis of their LogL.  Used for sorting.
 This is intended to sort ones with higher LogL first. */
int cmp_simmed_values(const void *a, const void *b)
{
	pbt_simmed_values_struct *ia = (pbt_simmed_values_struct *)a;
	pbt_simmed_values_struct *ib = (pbt_simmed_values_struct *)b;
	
	if(ia->LogL > ib->LogL) {
		return(-1);
	}
	else if(ia->LogL < ib->LogL) {
		return(1);
	}
	return(0);
}



int cmp_fdr_containers(const void *a, const void *b)
{
	pbt_fdr_container *ia = (pbt_fdr_container *)a;
	pbt_fdr_container *ib = (pbt_fdr_container *)b;
	
	if(ia->p_value < ib->p_value) {
		return(-1);
	}
	else if(ia->p_value > ib->p_value) {
		return(1);
	}
	return(0);
}


int cmp_doubles(const void *va, const void *vb) 
{
	double *a = (double *)va;
	double *b = (double *)vb;
	
	if (*a==*b)
		return(0);
	else if (*a < *b)
		return(-1);
	else
		return(1);
}


/* given lx values in x and ly values in y, this function sorts them (yes, it will modify them) and 
 then it determines the quantiles of x at 0, .01, .02, ..., 1.0.  It then returns the fraction 
 of y values that are less than the quantiles of x at 0, .01, .02, ... 1.0.
 
 This does a little quick and dirty interpolation.
*/
double *EmpiricalCrossCDF100(double *x, double *y, int lx, int ly)
{
	double *ret=(double *)calloc(101, sizeof(double));
	int i,j;
	double d,lo,hi,xval;
	int lox, hix;
	double frac;
	
	qsort(x,lx,sizeof(double),cmp_doubles);
	qsort(y,ly,sizeof(double),cmp_doubles);
	
	/*for(i=0;i<lx;i++) {
		printf("x  %d  %f\n",i,x[i]);
	}
	for(i=0;i<ly;i++) {
		printf("y  %d  %f\n",i,y[i]);
	}
	*/
	
	ret[0]=0.0;
	/* here, set ret[0] to the fraction of y's that are less than all the x's */
	for(j=0;j<ly;j++) {
		if(y[j]<x[0]) {
			ret[0]+=1.0;
		}
		else {
			break;
		}
	}
	ret[0] /= ly;
	
	/* now do the rest of the values */
	for(i=1,d=.01;  i<=100;  i++,d+=.01) {
		lo = lx*d;
		lox = (int)lo;
		hix = lox+1;
		
		if(hix<lx) {
			xval = x[lox] + (x[hix]-x[lox])/(double)(hix-lox) * (lo-(double)lox);
		}
		else if(lox<lx) {
			xval = x[lox];
		}
		else {
			xval = x[lx-1];
		}
		
		
		
		for(j=0;y[j]<=xval && j<ly ;j++) {
			;
		}
		/* here, j indexes the first value in y greater than xval */
		if(j>0) {
			j--;
		}
		//printf("xval i= %d   %f   j= %d\n",i,xval,j);
		
		/* now, the simple estimate of the fraction of y less than xval is (ly-j)/ly. but we can probabably do a little interpolation. */
		if(j>=ly-1)  {
			frac = 1.0;
		}
		else if(j==0) {
			frac = 0.0;
		}
		else {
			lo = y[j];
			hi = y[j+1];
			
			frac = ((double)j + (xval-lo)/(hi-lo)) /(double)ly;
		}
		
		ret[i]=frac;
	}
	return(ret);
}



/* this compiles all the offspring with max-LOD non-excluded parents in ThePop and it sorts them and estimates 
 the fraction of true parental trios amongst those, and then prints them all out in sorted order along with the 
 estimated False Discovery Rate and the estimated False Negative Rate.
 */
void DoFDR_ForAPop(int ThePop, pbt_high_level_vars *HLV) 
{
	int i,j,k,N=0,maxi,pv_int;
	double FDR,slope,newslope;
	double NumEstParental, NumEstNonParental;
	pbt_fdr_container *fdrc;
	double par_cdf[1001], par_dens[1001];
	double par_lower_than_nonpar[1001];  /* this is a quick little thing to attempt an estimate of the number of offspring having parents in S-down.  So, this computes that marginal probaiblity 
											that a kid in S-down has parents with lower LOD than the max-LOD nonparental pair */
	double marg_prob_par_lower_than_nonpar;
	

	//printf("Entered DoFDR_ForAPop for ThePop= %d\n",ThePop);

	for(i=0;i<HLV->PFR->NumOffColls;i++)  {
		for(j=0;j<HLV->PFR->NumInOffColls[i];j++)  {
			if(HLV->PFR->Offs[i][j].MendComps->NumParPairs>0 && HLV->PFR->Offs[i][j].MendComps->ParPairs[0].pop==ThePop) {
				N++;
			}
		}
	}
	
	//printf("Counted N=%d\n",N);

	if(N==0) {  /* just get the hell out if no one is assigned to these pops.  Otherwise PC's crash and burn
		     when calloc-ing length 0 */
	  //printf("About to return from DoFDR_ForApop.  N= %d\n",N);
	  return;	
	}
	
	//printf("Oops! Down here and N= %d\n",N);

	fdrc = (pbt_fdr_container *)calloc(N,sizeof(pbt_fdr_container));
	N=0;
	for(i=0;i<HLV->PFR->NumOffColls;i++)  {
		for(j=0;j<HLV->PFR->NumInOffColls[i];j++)  {
			if(HLV->PFR->Offs[i][j].MendComps->NumParPairs>0 && HLV->PFR->Offs[i][j].MendComps->ParPairs[0].pop==ThePop) {
				fdrc[N].ind = &(HLV->PFR->Offs[i][j]);
				fdrc[N].p_value = fdrc[N].ind->p_value;
				
#ifndef THIS_IS_NOT_DEFINED
				fdrc[N].parental_cdf = fdrc[N].ind->PValueParentalCDF;
				
				/* I am going to experiment here with approximating the densities of parental and non-parental conditional on offspring genotype here */
				par_dens[0] = fdrc[N].parental_cdf[0];
				for(k=1;k<101;k++) {
					par_dens[k] = fdrc[N].parental_cdf[k] - fdrc[N].parental_cdf[k-1];
				}
				fdrc[N].parental_p_value_dens = (double *)calloc(101,sizeof(double));
				fdrc[N].parental_p_value_dens[0] = 99999999999.99;
				for(k=1;k<101;k++) {
					fdrc[N].parental_p_value_dens[k] = par_dens[101-k];
				}
				if(fdrc[N].p_value < .001) { /* currently hard-wired to NumMC_Reps=1000 */
					fdrc[N].approx_p_dens = fdrc[N].parental_p_value_dens[0];
				}
				else { int temp;
					temp = (int)( (fdrc[N].p_value + .01) * 100.0);
					fdrc[N].approx_p_dens = fdrc[N].parental_p_value_dens[temp];
					if(fdrc[N].approx_p_dens<0.0) {
						fdrc[N].approx_p_dens = .000001;
					}
				}
#endif
				//printf("KID_APPROX_P_DENS:  p_value= %f    approx_dens= %f\n",fdrc[N].p_value,fdrc[N].approx_p_dens);
				
				N++;
			}
		}
	}
	
	//printf("Filled fdrc\n");
	
	/* now sort that according to p-value */
	qsort( (void *)(fdrc), (size_t)(N), sizeof(pbt_fdr_container),cmp_fdr_containers);
	
	//printf("Q-sorted that dude. N=%d\n",N);

#ifndef THIS_IS_NOT_DEFINED
	/* get the average values for the parental_cdfs */
	for(j=0;j<101;j++) {
		par_cdf[j]=0.0;
		par_lower_than_nonpar[j]=0.0;
	}
	for(i=0;i<N;i++) {
		for(j=0;j<101;j++) {
			par_cdf[j] += fdrc[i].parental_cdf[j];
			par_lower_than_nonpar[j]+= .01 * fdrc[i].parental_cdf[j];  /* the .01 arises because each element in the parental_cdf represents a little unit of probabilty mass of .01 corresponding to the 
																		fraction of nonparental pairs falling into that bin */
			//printf("BOO!:  PAR_CDF  i=  %d   j=  %d   val= %f \n",i,j,fdrc[i].parental_cdf[j]);
		}
	}
	marg_prob_par_lower_than_nonpar = 0.0;
	for(j=0;j<101;j++) {
		par_cdf[j] /= (double)N;
		par_lower_than_nonpar[j] /= (double)N;
		marg_prob_par_lower_than_nonpar += par_lower_than_nonpar[j];
		//printf("HERE_IS_THE_MARG_PAR_CDF:  %d    %f\n",j,par_cdf[j]);
	}
	

	/* now, store the histogram estimate of the density into ParentPValueDensityHistogramEstimate */
	HLV->KidsWithMaxLOD_Parents[ThePop]->ParentPValueDensityHistogramEstimate = (double *)calloc(101,sizeof(double));
	par_dens[0] = par_cdf[0];
	for(j=1;j<101;j++) {
		par_dens[j] = par_cdf[j] - par_cdf[j-1];
	}
	HLV->KidsWithMaxLOD_Parents[ThePop]->ParentPValueDensityHistogramEstimate[0]=99999999.9;
	for(j=1;j<101;j++) {
		HLV->KidsWithMaxLOD_Parents[ThePop]->ParentPValueDensityHistogramEstimate[j] = par_dens[101-j];
	}
	for(j=0;j<101;j++) {
		; //printf("ParentPValueDensityHistogramEstimate: %d  %f\n",j,HLV->KidsWithMaxLOD_Parents[ThePop]->ParentPValueDensityHistogramEstimate[j]);
	}
	
	/* now, we want to assign these approximate densities to the different individuals on the basis of their p-values */
	for(i=0;i<N;i++)  {
		if(fdrc[i].p_value < .001) { /* currently hard-wired to NumMC_Reps=1000 */
			fdrc[i].approx_p_dens = HLV->KidsWithMaxLOD_Parents[ThePop]->ParentPValueDensityHistogramEstimate[0];
		}
		else { int temp;
			temp = (int)( (fdrc[i].p_value + .01) * 100.0);
			fdrc[i].approx_p_dens = HLV->KidsWithMaxLOD_Parents[ThePop]->ParentPValueDensityHistogramEstimate[temp];
		}
		//printf("KID_APPROX_P_DENS:  p_value= %f    approx_dens= %f\n",fdrc[i].p_value,fdrc[i].approx_p_dens);
	}

	
	/* here, try to estimate the fraction of parentals by EM */
	//printf("ESTIMATED_NUM_PARENTALS_BY_EM: %f      N= %d\n",N*EstimateFractionParentalByEM_UsingApproxDensOfPvalues(fdrc, N, 0.01, .00001, 200), N);
#endif
	
	/* here we need to estimate the number of nonparental cases (M0) */
	slope=(1.0-fdrc[0].p_value)/(double)N;  /* start with this */
	if(N>1) {
	  for(i=1;i<N;i++) {
	    newslope = (1.0-fdrc[i].p_value)/(double)(N-i);
	    //printf("i=%d  N=%d  p-value=%f  slope=%f  newslope=%f\n",i,N,fdrc[i].p_value,slope,newslope);
	    if(newslope<slope) {
	      
	      slope=(newslope+slope)/2.0; /* here, instead of setting slope to newslope, I take the average of the two adjacent slopes.  This is 
					     a huge kluge to deal with situations where the p-values go along at about .02, .03, .05, and then they
					     jump immediately to 1.0. */
	      maxi = i;
	      
	      //printf("newslope<slope.  Set slope to average of two: %f  and maxi=%d\n",slope,maxi);
	      break;
	    }
	    else {
	      slope=newslope;
	      maxi = i;
	    }
	  }
	  /* here, slope is the slope of the first line that has a lesser slope than the previous one and maxi is the index of that point. */
	  //printf("About to try assigning to NumEst(Non)Parental\n");
	  NumEstParental = (maxi+1) - fdrc[maxi].p_value/slope;
	  NumEstNonParental = N - NumEstParental;
	}
	else {
	  NumEstNonParental = 1.0;  /* this is just a hack, so the FDR comes out to be the p.value down below */
	}
	//printf("Assigned Values to NumEst(Non)Parental\n");
	
	/* now print those results */
	for(i=0;i<N;i++)  {
		
		
		FDR = fdrc[i].p_value *  NumEstNonParental / (double)(i+1);
		pv_int = (int)((1.0-fdrc[i].p_value)*100.0);
		fdrc[i].ind->FDR = FDR;  /* record the FDR value in the individual */
		
		/* have to take the average over all individuals of the FNR's here.  Wait! Do it for each of the 100 values in those arrays first (outside of this loop)
		 and then use the results.  */
		
		
		fprintf(HLV->FDR_Summaries_File,"%s\t%d\t%s\t%s\t%s\t%f\t%.2f\t%f\n",
			   HLV->PFR->PopNames[ThePop],
			   i+1,
			   /*NumEstParental,*/
			   fdrc[i].ind->Name,
			   fdrc[i].ind->MendComps->ParPairs[0].mamapapa[MALE]->Name,
			   fdrc[i].ind->MendComps->ParPairs[0].mamapapa[FEMALE]->Name,
			   FDR,
			   FDR*(i+1), 
			   fdrc[i].p_value
			   /*par_cdf[pv_int],
			   par_cdf[pv_int] * NumEstParental,
			   pv_int,
			   NumEstParental / (1.0-marg_prob_par_lower_than_nonpar)*/
			   );
		
	}
	
	
	free(fdrc);
}



/* given N individuals with non-excluded parents from a certain population.  Estimate the fraction
 of those that are truly parental given the approximate density values in teh pbt_fdr_containers, and 
 the approximate density of the Null Distribution (which in our case will typically be .01).  Iterate until
 a change of less than Precis occurs.  However, don't do more than MaxIts iterations.  */
double EstimateFractionParentalByEM_UsingApproxDensOfPvalues(pbt_fdr_container *f, int N, double NullDens, double Precis, int MaxIts)
{
	int i,j;
	double p=.5,q=.5;  /* p=frac of parentals.  Start assuming that it is half and half */
	double sum;
	
	for(j=0;j<MaxIts;j++)  {
		for(i=0;i<N;i++)  {
			sum+= p*f[i].approx_p_dens / (p*f[i].approx_p_dens + q*NullDens);
		}
		sum /= (double)N;
		
		if( ECA_ABS(p-sum) < Precis) {
			return(sum);
		}
		p=sum;
		q=1.0-p;
		//printf("EM_PROGRESS   j= %d   p= %f\n",j,p);
	}
	
	return(p);
}








/* this is a silly function that I can use to return a reasonable value of a pi-vector for a value AVE of the "average" run size of a chinook population. */
double *Return_ArchTypeChinookPiVec(int AVE, int *outx, int *outpsz)
{
	int i,x=-1;
	double *ret;
	
	
	
	/* here we statically initialize these monster arrays */
	int NumInArrays=245;
	int psz[]={10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,500,550,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500,1550,1600,1650,1700,1750,1800,1850,1900,1950,2000,2050,2100,2150,2200,2250,2300,2350,2400,2450,2500,2550,2600,2650,2700,2750,2800,2850,2900,2950,3000,3100,3200,3300,3400,3500,3600,3700,3800,3900,4000,4100,4200,4300,4400,4500,4600,4700,4800,4900,5000,5100,5200,5300,5400,5500,5600,5700,5800,5900,6000,6200,6400,6600,6800,7000,7200,7400,7600,7800,8000,8200,8400,8600,8800,9000,9200,9400,9600,9800,10000,10500,11000,11500,12000,12500,13000,13500,14000,14500,15000,15500,16000,16500,17000,17500,18000,18500,19000,19500,20000,20500,21000,21500,22000,22500,23000,23500,24000,24500,25000,25500,26000,26500,27000,27500,28000,28500,29000,29500,30000,31000,32000,33000,34000,35000,36000,37000,38000,39000,40000,41000,42000,43000,44000,45000,46000,47000,48000,49000,50000,51000,52000,53000,54000,55000,56000,57000,58000,59000,60000,61000,62000,63000,64000,65000,66000,67000,68000,69000,70000,71000,72000,73000,74000,75000,76000,77000,78000,79000,80000,81000,82000,83000,84000,85000,86000,87000,88000,89000,90000,91000,92000,93000,94000,95000,96000,97000,98000,99000,100000};
	double pivec[][18]={
		{1.685957e-02,2.235759e-02,2.268928e-02,7.333640e-04,-7.050800e-03,1.114592e-01,2.107540e-01,2.062806e-01,4.255883e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.914812e-02,2.286543e-02,2.313388e-02},
		{2.695106e-03,3.779132e-03,4.228146e-03,1.076819e-02,1.052405e-02,2.378956e-02,1.339345e-01,1.338830e-01,6.646894e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.733766e-04,0.000000e+00,1.727543e-02,1.099622e-02,4.345458e-03},
		{1.006946e-03,1.430586e-03,1.530362e-03,9.239473e-03,8.912462e-03,8.114361e-03,8.834950e-02,8.840303e-02,7.776914e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.720801e-04,2.259771e-05,1.271723e-02,7.620548e-03,1.587770e-03},
		{5.302330e-04,7.750115e-04,8.296662e-04,7.299230e-03,7.322652e-03,4.434650e-03,6.726292e-02,6.759944e-02,8.312178e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.104238e-04,1.237980e-05,9.389796e-03,6.119343e-03,8.388102e-04},
		{3.369412e-04,4.977446e-04,5.349365e-04,6.047238e-03,6.214422e-03,2.884920e-03,5.488097e-02,5.578235e-02,8.612384e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,9.124529e-05,9.575418e-06,8.069634e-03,5.299395e-03,5.610095e-04},
		{2.209739e-04,3.262485e-04,3.423652e-04,5.119101e-03,5.250567e-03,1.859229e-03,4.486859e-02,4.573373e-02,8.864655e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,6.151921e-05,6.302962e-06,6.663234e-03,4.294515e-03,3.541323e-04},
		{1.639098e-04,2.406016e-04,2.545367e-04,4.490452e-03,4.641887e-03,1.349348e-03,3.869366e-02,3.949508e-02,9.018033e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,5.236103e-05,5.147339e-06,5.864809e-03,3.812809e-03,2.568006e-04},
		{1.212598e-04,1.796071e-04,1.874860e-04,3.970492e-03,4.037949e-03,1.005310e-03,3.362834e-02,3.427626e-02,9.150077e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.252619e-05,3.592936e-06,4.930149e-03,3.271638e-03,1.861145e-04},
		{9.396047e-05,1.389213e-04,1.453530e-04,3.484826e-03,3.661122e-03,7.959816e-04,3.006629e-02,3.066579e-02,9.239494e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.845645e-05,3.209444e-06,4.525450e-03,2.961293e-03,1.500481e-04},
		{7.777410e-05,1.159898e-04,1.225207e-04,3.200534e-03,3.352041e-03,6.633475e-04,2.736304e-02,2.801340e-02,9.306645e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.238327e-05,2.505384e-06,4.114650e-03,2.719190e-03,1.240676e-04},
		{6.109101e-05,9.046463e-05,9.401138e-05,2.862350e-03,3.019820e-03,5.105841e-04,2.433082e-02,2.490519e-02,9.382962e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.993749e-05,2.347410e-06,3.697453e-03,2.440166e-03,9.549688e-05},
		{5.262744e-05,7.815035e-05,8.161900e-05,2.690097e-03,2.796505e-03,4.396984e-04,2.249864e-02,2.303330e-02,9.428833e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.585677e-05,1.741530e-06,3.453152e-03,2.258046e-03,8.397775e-05},
		{4.410843e-05,6.588169e-05,6.829171e-05,2.467878e-03,2.588077e-03,3.709920e-04,2.073846e-02,2.129489e-02,9.473618e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.378218e-05,1.517633e-06,3.152029e-03,2.072326e-03,6.844053e-05},
		{3.816024e-05,5.692601e-05,5.947384e-05,2.311304e-03,2.410971e-03,3.208789e-04,1.936531e-02,1.979408e-02,9.509688e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.179595e-05,1.417427e-06,2.918458e-03,1.950669e-03,5.960938e-05},
		{3.344667e-05,5.000534e-05,5.248765e-05,2.170370e-03,2.267927e-03,2.843896e-04,1.817580e-02,1.859052e-02,9.539399e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,9.940572e-06,1.136285e-06,2.754227e-03,1.854187e-03,5.372620e-05},
		{2.927398e-05,4.367904e-05,4.577460e-05,2.043245e-03,2.124408e-03,2.459183e-04,1.695955e-02,1.734192e-02,9.570473e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,8.785855e-06,9.903285e-07,2.568562e-03,1.699892e-03,4.537593e-05},
		{2.564664e-05,3.838548e-05,3.996929e-05,1.914462e-03,2.002763e-03,2.158920e-04,1.592641e-02,1.633106e-02,9.595960e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,7.832563e-06,9.307885e-07,2.424283e-03,1.615740e-03,4.038906e-05},
		{2.244885e-05,3.354662e-05,3.479955e-05,1.796907e-03,1.880539e-03,1.890417e-04,1.492468e-02,1.529184e-02,9.621511e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,6.983177e-06,7.848283e-07,2.276422e-03,1.513206e-03,3.511869e-05},
		{2.048573e-05,3.062575e-05,3.203655e-05,1.716484e-03,1.803649e-03,1.729034e-04,1.428028e-02,1.464035e-02,9.637810e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,6.369440e-06,7.523937e-07,2.179970e-03,1.447018e-03,3.201317e-05},
		{1.845025e-05,2.780054e-05,2.895270e-05,1.630248e-03,1.717252e-03,1.583664e-04,1.365740e-02,1.400125e-02,9.653910e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,5.802333e-06,6.689800e-07,2.079778e-03,1.386323e-03,2.925188e-05},
		{1.521773e-05,2.276687e-05,2.371810e-05,1.485042e-03,1.567976e-03,1.286445e-04,1.236730e-02,1.268108e-02,9.686271e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.945880e-06,5.824071e-07,1.898794e-03,1.260191e-03,2.387473e-05},
		{1.267126e-05,1.900099e-05,1.969286e-05,1.366792e-03,1.428552e-03,1.067517e-04,1.126905e-02,1.155102e-02,9.714131e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.019930e-06,4.569356e-07,1.727964e-03,1.149927e-03,1.980228e-05},
		{1.081638e-05,1.627530e-05,1.685040e-05,1.265499e-03,1.324036e-03,9.143226e-05,1.043089e-02,1.070673e-02,9.735293e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.463226e-06,3.908335e-07,1.601900e-03,1.061568e-03,1.683944e-05},
		{9.346671e-06,1.401620e-05,1.468395e-05,1.174475e-03,1.237777e-03,7.950760e-05,9.731509e-03,9.972695e-03,9.753302e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.901197e-06,3.466595e-07,1.484769e-03,9.993361e-04,1.476064e-05},
		{8.127958e-06,1.222211e-05,1.274316e-05,1.098988e-03,1.154812e-03,6.921936e-05,9.080370e-03,9.312302e-03,9.769770e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.535142e-06,3.018949e-07,1.384542e-03,9.318194e-04,1.279307e-05},
		{7.118446e-06,1.072059e-05,1.116739e-05,1.032707e-03,1.080612e-03,6.049918e-05,8.490937e-03,8.709483e-03,9.784665e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.195239e-06,2.527571e-07,1.296817e-03,8.702307e-04,1.116191e-05},
		{6.302311e-06,9.449412e-06,9.844544e-06,9.723801e-04,1.019437e-03,5.307681e-05,7.973799e-03,8.179889e-03,9.797649e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.987990e-06,2.319591e-07,1.224107e-03,8.188295e-04,9.823415e-06},
		{5.585509e-06,8.377068e-06,8.687721e-06,9.164102e-04,9.616742e-04,4.712473e-05,7.516695e-03,7.707113e-03,9.809362e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.738294e-06,2.077912e-07,1.148946e-03,7.717162e-04,8.711531e-06},
		{5.020238e-06,7.561179e-06,7.857474e-06,8.678056e-04,9.150218e-04,4.273949e-05,7.159507e-03,7.351698e-03,9.818300e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.612931e-06,1.896739e-07,1.100623e-03,7.380677e-04,7.911850e-06},
		{4.541145e-06,6.836988e-06,7.115230e-06,8.269251e-04,8.703870e-04,3.862459e-05,6.803013e-03,6.984978e-03,9.827350e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.440429e-06,1.654822e-07,1.046763e-03,6.991806e-04,7.131969e-06},
		{4.100500e-06,6.177783e-06,6.416612e-06,7.871978e-04,8.270854e-04,3.489328e-05,6.475170e-03,6.641799e-03,9.835833e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.303465e-06,1.544802e-07,9.887148e-04,6.663311e-04,6.403023e-06},
		{3.730893e-06,5.612159e-06,5.829173e-06,7.517344e-04,7.898995e-04,3.160248e-05,6.168976e-03,6.329016e-03,9.843500e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.198053e-06,1.403025e-07,9.480063e-04,6.347274e-04,5.830688e-06},
		{3.425267e-06,5.157184e-06,5.370550e-06,7.190575e-04,7.589976e-04,2.915927e-05,5.919798e-03,6.077670e-03,9.849759e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.106472e-06,1.293558e-07,9.127002e-04,6.104043e-04,5.392132e-06},
		{3.139363e-06,4.724709e-06,4.908011e-06,6.904366e-04,7.254603e-04,2.663160e-05,5.661633e-03,5.808628e-03,9.856378e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,9.997528e-07,1.181912e-07,8.685773e-04,5.841639e-04,4.914282e-06},
		{2.895218e-06,4.357643e-06,4.530714e-06,6.627769e-04,6.979765e-04,2.458030e-05,5.439332e-03,5.586320e-03,9.861925e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,9.347869e-07,1.090222e-07,8.384579e-04,5.611323e-04,4.538521e-06},
		{2.384029e-06,3.593205e-06,3.733012e-06,6.032025e-04,6.337244e-04,2.026777e-05,4.947633e-03,5.079558e-03,9.874504e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,7.704427e-07,8.966318e-08,7.592309e-04,5.085831e-04,3.719683e-06},
		{2.000584e-06,3.015056e-06,3.132997e-06,5.531098e-04,5.814212e-04,1.702947e-05,4.534279e-03,4.655361e-03,9.884970e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,6.438419e-07,7.490559e-08,6.970443e-04,4.669096e-04,3.134039e-06},
		{1.696973e-06,2.554784e-06,2.652586e-06,5.093375e-04,5.367521e-04,1.441674e-05,4.177233e-03,4.289164e-03,9.894011e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,5.538193e-07,6.520583e-08,6.428980e-04,4.309397e-04,2.652691e-06},
		{1.470845e-06,2.215549e-06,2.304729e-06,4.752604e-04,4.993122e-04,1.249621e-05,3.887153e-03,3.989529e-03,9.901403e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.713303e-07,5.498274e-08,5.973057e-04,4.002299e-04,2.294914e-06},
		{1.278447e-06,1.926322e-06,2.001173e-06,4.436432e-04,4.657589e-04,1.086023e-05,3.625024e-03,3.720856e-03,9.908048e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.086516e-07,4.779957e-08,5.562433e-04,3.741786e-04,1.994779e-06},
		{1.122712e-06,1.693683e-06,1.756703e-06,4.157732e-04,4.370814e-04,9.543801e-06,3.399338e-03,3.491496e-03,9.913737e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.641661e-07,4.230056e-08,5.237363e-04,3.505267e-04,1.752718e-06},
		{9.925998e-07,1.496178e-06,1.554387e-06,3.909577e-04,4.114993e-04,8.439466e-06,3.198971e-03,3.284025e-03,9.918851e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.198751e-07,3.781174e-08,4.919919e-04,3.300755e-04,1.551345e-06},
		{8.852518e-07,1.333354e-06,1.383606e-06,3.695041e-04,3.886916e-04,7.514784e-06,3.017979e-03,3.098226e-03,9.923432e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.835344e-07,3.358315e-08,4.639633e-04,3.118319e-04,1.382571e-06},
		{7.969862e-07,1.201277e-06,1.250347e-06,3.506242e-04,3.690870e-04,6.783577e-06,2.867518e-03,2.943251e-03,9.927255e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.565925e-07,3.025128e-08,4.414038e-04,2.967358e-04,1.249816e-06},
		{7.181472e-07,1.083125e-06,1.126016e-06,3.326895e-04,3.508160e-04,6.117546e-06,2.723033e-03,2.796756e-03,9.930896e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.332761e-07,2.737745e-08,4.200342e-04,2.817337e-04,1.126215e-06},
		{6.506972e-07,9.816415e-07,1.019256e-06,3.170288e-04,3.338305e-04,5.539529e-06,2.591337e-03,2.661666e-03,9.934242e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.092842e-07,2.468747e-08,3.987390e-04,2.683549e-04,1.019379e-06},
		{5.928156e-07,8.944482e-07,9.294909e-07,3.026945e-04,3.187581e-04,5.054595e-06,2.476408e-03,2.543484e-03,9.937177e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.912969e-07,2.246001e-08,3.807670e-04,2.557811e-04,9.271415e-07},
		{5.426154e-07,8.182137e-07,8.502807e-07,2.897590e-04,3.049989e-04,4.615569e-06,2.367309e-03,2.431304e-03,9.939928e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.768707e-07,2.062703e-08,3.648973e-04,2.448880e-04,8.483212e-07},
		{4.978447e-07,7.507146e-07,7.797568e-07,2.776790e-04,2.921955e-04,4.234964e-06,2.267323e-03,2.328698e-03,9.942462e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.612569e-07,1.886253e-08,3.496922e-04,2.345216e-04,7.797801e-07},
		{4.584664e-07,6.916295e-07,7.182669e-07,2.665089e-04,2.805345e-04,3.901182e-06,2.176741e-03,2.235840e-03,9.944765e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.483422e-07,1.748695e-08,3.352939e-04,2.251863e-04,7.173924e-07},
		{4.239860e-07,6.393863e-07,6.641386e-07,2.561891e-04,2.700155e-04,3.608846e-06,2.094325e-03,2.151037e-03,9.946857e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.382748e-07,1.628064e-08,3.229077e-04,2.166499e-04,6.631772e-07},
		{3.918505e-07,5.913072e-07,6.139972e-07,2.463837e-04,2.596392e-04,3.338268e-06,2.013848e-03,2.068706e-03,9.948896e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.273708e-07,1.499016e-08,3.102998e-04,2.086220e-04,6.137087e-07},
		{3.654353e-07,5.516003e-07,5.730553e-07,2.381722e-04,2.505534e-04,3.112717e-06,1.944040e-03,1.996603e-03,9.950675e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.177648e-07,1.379739e-08,2.992968e-04,2.011535e-04,5.716471e-07},
		{3.408045e-07,5.140885e-07,5.343473e-07,2.298966e-04,2.421729e-04,2.902837e-06,1.878002e-03,1.928609e-03,9.952350e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.106471e-07,1.295112e-08,2.893951e-04,1.943795e-04,5.330469e-07},
		{3.180230e-07,4.800127e-07,4.983731e-07,2.221869e-04,2.339151e-04,2.709406e-06,1.814671e-03,1.863578e-03,9.953962e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.030349e-07,1.215368e-08,2.792123e-04,1.878557e-04,4.969840e-07},
		{2.978560e-07,4.494520e-07,4.669118e-07,2.150138e-04,2.264898e-04,2.537532e-06,1.756041e-03,1.803538e-03,9.955445e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,9.615461e-08,1.135022e-08,2.704577e-04,1.817575e-04,4.661897e-07},
		{2.795938e-07,4.218309e-07,4.380933e-07,2.084044e-04,2.194189e-04,2.378505e-06,1.700500e-03,1.746479e-03,9.956846e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,9.089720e-08,1.067658e-08,2.623466e-04,1.761345e-04,4.373971e-07},
		{2.627966e-07,3.965318e-07,4.118355e-07,2.021128e-04,2.127139e-04,2.237841e-06,1.649819e-03,1.694396e-03,9.958141e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,8.553751e-08,1.002719e-08,2.542470e-04,1.706753e-04,4.107959e-07},
		{2.475290e-07,3.735419e-07,3.878638e-07,1.962232e-04,2.064338e-04,2.107167e-06,1.600934e-03,1.644237e-03,9.959382e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,8.012616e-08,9.413634e-09,2.464720e-04,1.656285e-04,3.866121e-07},
		{2.334456e-07,3.523114e-07,3.656965e-07,1.904771e-04,2.006200e-04,1.987775e-06,1.554852e-03,1.597424e-03,9.960543e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,7.595105e-08,8.900645e-09,2.396673e-04,1.608919e-04,3.647256e-07},
		{2.203767e-07,3.325118e-07,3.451540e-07,1.851248e-04,1.949196e-04,1.876109e-06,1.510724e-03,1.551726e-03,9.961666e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,7.158576e-08,8.390792e-09,2.328577e-04,1.563958e-04,3.445154e-07},
		{2.092428e-07,3.157575e-07,3.280825e-07,1.803680e-04,1.899766e-04,1.782523e-06,1.472141e-03,1.512256e-03,9.962642e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,6.800995e-08,7.962788e-09,2.269945e-04,1.525376e-04,3.274843e-07},
		{1.981963e-07,2.990481e-07,3.108540e-07,1.755345e-04,1.849518e-04,1.688807e-06,1.433282e-03,1.472030e-03,9.963634e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,6.421945e-08,7.576880e-09,2.208163e-04,1.484938e-04,3.102800e-07},
		{1.880604e-07,2.838859e-07,2.947722e-07,1.709999e-04,1.801891e-04,1.602685e-06,1.396319e-03,1.434448e-03,9.964566e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,6.124334e-08,7.209542e-09,2.153612e-04,1.446548e-04,2.945230e-07},
		{1.787218e-07,2.698015e-07,2.801368e-07,1.667241e-04,1.756848e-04,1.522975e-06,1.361252e-03,1.398426e-03,9.965456e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,5.815620e-08,6.845694e-09,2.099135e-04,1.410709e-04,2.798072e-07},
		{1.704437e-07,2.572783e-07,2.673685e-07,1.628991e-04,1.714997e-04,1.452089e-06,1.328918e-03,1.364957e-03,9.966280e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,5.513169e-08,6.469093e-09,2.048140e-04,1.376491e-04,2.667079e-07},
		{1.621306e-07,2.446868e-07,2.540466e-07,1.588781e-04,1.672903e-04,1.380416e-06,1.296054e-03,1.331162e-03,9.967114e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,5.262024e-08,6.201430e-09,1.997383e-04,1.343078e-04,2.535223e-07},
		{1.547337e-07,2.334691e-07,2.425772e-07,1.551761e-04,1.635029e-04,1.317841e-06,1.266194e-03,1.300556e-03,9.967869e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,5.019666e-08,5.908421e-09,1.952805e-04,1.312477e-04,2.422297e-07},
		{1.476699e-07,2.229204e-07,2.315255e-07,1.516676e-04,1.596907e-04,1.258228e-06,1.237575e-03,1.271071e-03,9.968601e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.790412e-08,5.634669e-09,1.906442e-04,1.281218e-04,2.308759e-07},
		{1.411162e-07,2.129709e-07,2.211503e-07,1.482290e-04,1.561592e-04,1.202478e-06,1.209916e-03,1.242758e-03,9.969301e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.592790e-08,5.398782e-09,1.865250e-04,1.252914e-04,2.207711e-07},
		{1.350925e-07,2.038975e-07,2.117439e-07,1.450754e-04,1.527646e-04,1.150520e-06,1.183590e-03,1.215677e-03,9.969969e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.394837e-08,5.180550e-09,1.823919e-04,1.225936e-04,2.111679e-07},
		{1.294142e-07,1.953114e-07,2.028742e-07,1.419825e-04,1.495623e-04,1.102183e-06,1.158379e-03,1.189937e-03,9.970606e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.207630e-08,4.954887e-09,1.786011e-04,1.199675e-04,2.022947e-07},
		{1.240043e-07,1.871658e-07,1.942651e-07,1.390341e-04,1.463711e-04,1.055717e-06,1.133796e-03,1.164615e-03,9.971230e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.035777e-08,4.735538e-09,1.748526e-04,1.174497e-04,1.938369e-07},
		{1.190540e-07,1.796610e-07,1.865167e-07,1.362207e-04,1.434502e-04,1.013481e-06,1.110846e-03,1.141037e-03,9.971812e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.870161e-08,4.541620e-09,1.712726e-04,1.150795e-04,1.860599e-07},
		{1.143528e-07,1.725506e-07,1.792104e-07,1.335214e-04,1.405854e-04,9.731023e-07,1.088536e-03,1.118147e-03,9.972377e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.715088e-08,4.369003e-09,1.678620e-04,1.127892e-04,1.786861e-07},
		{1.098293e-07,1.657947e-07,1.720717e-07,1.308523e-04,1.377980e-04,9.352382e-07,1.067200e-03,1.096339e-03,9.972919e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.575288e-08,4.207630e-09,1.645166e-04,1.105353e-04,1.715687e-07},
		{1.056505e-07,1.594676e-07,1.655583e-07,1.283479e-04,1.351659e-04,8.995569e-07,1.046722e-03,1.075161e-03,9.973439e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.438643e-08,4.044972e-09,1.614138e-04,1.084546e-04,1.651448e-07},
		{1.017279e-07,1.535129e-07,1.594101e-07,1.259325e-04,1.326585e-04,8.660154e-07,1.026962e-03,1.054911e-03,9.973940e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.306559e-08,3.888028e-09,1.583583e-04,1.063977e-04,1.589652e-07},
		{9.774562e-08,1.475222e-07,1.531493e-07,1.234353e-04,1.300573e-04,8.324702e-07,1.007056e-03,1.034515e-03,9.974446e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.189999e-08,3.755856e-09,1.553067e-04,1.043478e-04,1.527950e-07},
		{9.449565e-08,1.426763e-07,1.481938e-07,1.213916e-04,1.278535e-04,8.054664e-07,9.903743e-04,1.017353e-03,9.974870e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.075520e-08,3.616466e-09,1.527025e-04,1.026376e-04,1.478761e-07},
		{9.116465e-08,1.376395e-07,1.429653e-07,1.192288e-04,1.255960e-04,7.770412e-07,9.726760e-04,9.992291e-04,9.975319e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.963059e-08,3.483088e-09,1.499632e-04,1.007760e-04,1.426042e-07},
		{8.797584e-08,1.328115e-07,1.379489e-07,1.171461e-04,1.233750e-04,7.497504e-07,9.554948e-04,9.814687e-04,9.975755e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.857916e-08,3.360112e-09,1.473146e-04,9.905901e-05,1.377147e-07},
		{8.493154e-08,1.282612e-07,1.331774e-07,1.150864e-04,1.212526e-04,7.242397e-07,9.391522e-04,9.647592e-04,9.976170e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.764952e-08,3.251011e-09,1.447956e-04,9.730140e-05,1.328902e-07},
		{8.209244e-08,1.239433e-07,1.287606e-07,1.131538e-04,1.192102e-04,6.997020e-07,9.231472e-04,9.483218e-04,9.976575e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.671358e-08,3.146453e-09,1.423522e-04,9.567479e-05,1.284424e-07},
		{7.938051e-08,1.198459e-07,1.244652e-07,1.112833e-04,1.172172e-04,6.765301e-07,9.076866e-04,9.324050e-04,9.976968e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.579204e-08,3.039910e-09,1.399072e-04,9.407941e-05,1.241704e-07},
		{7.433408e-08,1.122300e-07,1.165431e-07,1.076937e-04,1.134515e-04,6.334066e-07,8.782964e-04,9.023008e-04,9.977712e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.417199e-08,2.837223e-09,1.354854e-04,9.102151e-05,1.162681e-07},
		{6.975612e-08,1.053042e-07,1.093427e-07,1.043275e-04,1.099111e-04,5.941952e-07,8.507573e-04,8.740291e-04,9.978410e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.272767e-08,2.672027e-09,1.312646e-04,8.819603e-05,1.090899e-07},
		{6.555975e-08,9.898367e-08,1.027731e-07,1.011693e-04,1.065474e-04,5.586581e-07,8.249993e-04,8.475111e-04,9.979066e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.134217e-08,2.507414e-09,1.272240e-04,8.548460e-05,1.025220e-07},
		{6.178464e-08,9.326923e-08,9.686369e-08,9.822566e-05,1.034303e-04,5.262670e-07,8.007447e-04,8.225206e-04,9.979682e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.008896e-08,2.360816e-09,1.234614e-04,8.295446e-05,9.654423e-08},
		{5.827946e-08,8.800221e-08,9.133727e-08,9.540443e-05,1.004666e-04,4.964363e-07,7.777527e-04,7.990212e-04,9.980264e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.899341e-08,2.229402e-09,1.199771e-04,8.058565e-05,9.107474e-08},
		{5.506813e-08,8.311825e-08,8.627650e-08,9.273947e-05,9.767117e-05,4.688921e-07,7.558984e-04,7.764659e-04,9.980819e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.788822e-08,2.109261e-09,1.165386e-04,7.835024e-05,8.606260e-08},
		{5.215979e-08,7.875261e-08,8.178492e-08,9.023371e-05,9.509121e-05,4.445085e-07,7.359326e-04,7.560407e-04,9.981325e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.700122e-08,1.999974e-09,1.135547e-04,7.630175e-05,8.160626e-08},
		{4.944133e-08,7.465461e-08,7.752219e-08,8.787178e-05,9.256782e-05,4.214507e-07,7.166438e-04,7.362225e-04,9.981816e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.608316e-08,1.892098e-09,1.104751e-04,7.426116e-05,7.729724e-08},
		{4.693581e-08,7.085721e-08,7.356961e-08,8.563203e-05,9.018678e-05,3.998097e-07,6.980163e-04,7.170784e-04,9.982287e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.528256e-08,1.795749e-09,1.076981e-04,7.235699e-05,7.341455e-08},
		{4.461577e-08,6.735438e-08,6.993399e-08,8.347281e-05,8.795652e-05,3.801084e-07,6.807055e-04,6.992928e-04,9.982727e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.455740e-08,1.713978e-09,1.050152e-04,7.056087e-05,6.975683e-08},
		{4.246409e-08,6.412463e-08,6.658514e-08,8.146475e-05,8.578610e-05,3.618574e-07,6.640216e-04,6.821405e-04,9.983150e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.379906e-08,1.621972e-09,1.023887e-04,6.883390e-05,6.639540e-08},
		{4.045520e-08,6.108940e-08,6.342016e-08,7.950820e-05,8.374638e-05,3.447596e-07,6.482363e-04,6.659005e-04,9.983552e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.316859e-08,1.551478e-09,9.992964e-05,6.720576e-05,6.323254e-08},
		{3.860204e-08,5.828125e-08,6.051156e-08,7.767003e-05,8.181121e-05,3.287848e-07,6.330534e-04,6.503200e-04,9.983935e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.258006e-08,1.479850e-09,9.768188e-05,6.564309e-05,6.035972e-08},
		{3.686171e-08,5.565788e-08,5.778362e-08,7.591363e-05,7.993648e-05,3.140361e-07,6.187281e-04,6.356030e-04,9.984300e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.199233e-08,1.410795e-09,9.539435e-05,6.412713e-05,5.760330e-08},
		{3.521713e-08,5.317055e-08,5.519467e-08,7.419287e-05,7.814879e-05,3.000257e-07,6.047712e-04,6.212847e-04,9.984654e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.147191e-08,1.348277e-09,9.328334e-05,6.269413e-05,5.505265e-08},
		{3.373207e-08,5.092793e-08,5.289723e-08,7.260317e-05,7.649588e-05,2.874765e-07,5.919670e-04,6.081029e-04,9.984979e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.097615e-08,1.293424e-09,9.128916e-05,6.138459e-05,5.276204e-08},
		{3.230204e-08,4.877881e-08,5.064160e-08,7.105021e-05,7.486259e-05,2.752966e-07,5.793066e-04,5.951620e-04,9.985299e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.052887e-08,1.239536e-09,8.937004e-05,6.007778e-05,5.052734e-08},
		{3.097340e-08,4.676646e-08,4.855136e-08,6.958571e-05,7.329538e-05,2.638742e-07,5.671646e-04,5.826325e-04,9.985608e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.008397e-08,1.187777e-09,8.747294e-05,5.881651e-05,4.842411e-08},
		{2.971814e-08,4.487594e-08,4.659581e-08,6.816632e-05,7.179835e-05,2.532613e-07,5.556741e-04,5.708209e-04,9.985900e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,9.672639e-09,1.138344e-09,8.568297e-05,5.759868e-05,4.645869e-08},
		{2.854073e-08,4.309419e-08,4.474366e-08,6.680202e-05,7.036398e-05,2.431737e-07,5.445225e-04,5.593742e-04,9.986183e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,9.300870e-09,1.095666e-09,8.397534e-05,5.645419e-05,4.461388e-08},
		{2.742655e-08,4.141195e-08,4.298644e-08,6.549164e-05,6.897621e-05,2.336401e-07,5.337384e-04,5.483123e-04,9.986456e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,8.939207e-09,1.050810e-09,8.234066e-05,5.533813e-05,4.287233e-08},
		{2.639027e-08,3.984167e-08,4.136914e-08,6.424224e-05,6.766396e-05,2.247603e-07,5.234919e-04,5.377901e-04,9.986716e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,8.593912e-09,1.011271e-09,8.075694e-05,5.428134e-05,4.124535e-08},
		{2.539570e-08,3.834547e-08,3.980820e-08,6.302004e-05,6.638288e-05,2.163462e-07,5.136245e-04,5.276437e-04,9.986966e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,8.277065e-09,9.741001e-10,7.923414e-05,5.325695e-05,3.969748e-08},
		{2.443641e-08,3.689503e-08,3.830155e-08,6.181214e-05,6.512651e-05,2.081992e-07,5.038854e-04,5.176621e-04,9.987214e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,7.974997e-09,9.389639e-10,7.773366e-05,5.224718e-05,3.819875e-08},
		{2.359228e-08,3.562682e-08,3.699497e-08,6.073902e-05,6.398835e-05,2.010874e-07,4.951290e-04,5.086737e-04,9.987436e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,7.683493e-09,9.039082e-10,7.637162e-05,5.133617e-05,3.689193e-08},
		{2.275031e-08,3.435971e-08,3.567290e-08,5.964536e-05,6.284148e-05,1.939561e-07,4.862927e-04,4.995870e-04,9.987660e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,7.415229e-09,8.724477e-10,7.500831e-05,5.041899e-05,3.558052e-08},
		{2.196193e-08,3.316427e-08,3.443558e-08,5.860730e-05,6.173991e-05,1.871811e-07,4.777178e-04,4.907685e-04,9.987878e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,7.150783e-09,8.423864e-10,7.367485e-05,4.953985e-05,3.434282e-08},
		{2.121043e-08,3.202979e-08,3.325580e-08,5.759625e-05,6.067904e-05,1.807585e-07,4.694501e-04,4.822989e-04,9.988087e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,6.908839e-09,8.121011e-10,7.242998e-05,4.867871e-05,3.316657e-08},
		{2.049718e-08,3.095040e-08,3.213382e-08,5.661922e-05,5.965136e-05,1.746475e-07,4.614658e-04,4.741017e-04,9.988289e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,6.683392e-09,7.862718e-10,7.120596e-05,4.785872e-05,3.204862e-08},
		{1.981456e-08,2.992190e-08,3.106496e-08,5.567594e-05,5.864653e-05,1.688650e-07,4.537809e-04,4.661894e-04,9.988485e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,6.457697e-09,7.594247e-10,6.999999e-05,4.704894e-05,3.098050e-08},
		{1.855701e-08,2.802461e-08,2.908915e-08,5.388410e-05,5.675570e-05,1.581143e-07,4.391088e-04,4.511299e-04,9.988857e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,6.051452e-09,7.112563e-10,6.774841e-05,4.552572e-05,2.900326e-08},
		{1.741780e-08,2.630238e-08,2.731065e-08,5.219484e-05,5.499932e-05,1.484452e-07,4.254587e-04,4.371000e-04,9.989203e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,5.680276e-09,6.684543e-10,6.565170e-05,4.412611e-05,2.724096e-08},
		{1.637626e-08,2.472777e-08,2.567253e-08,5.062130e-05,5.332218e-05,1.395312e-07,4.125045e-04,4.237875e-04,9.989532e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,5.337240e-09,6.276566e-10,6.365054e-05,4.277743e-05,2.560847e-08},
		{1.542652e-08,2.329753e-08,2.418897e-08,4.913545e-05,5.175255e-05,1.314686e-07,4.003909e-04,4.113366e-04,9.989840e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,5.021546e-09,5.906355e-10,6.175788e-05,4.152103e-05,2.411808e-08},
		{1.455702e-08,2.198170e-08,2.282148e-08,4.772896e-05,5.027882e-05,1.240206e-07,3.889189e-04,3.995477e-04,9.990131e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.746618e-09,5.585571e-10,6.001363e-05,4.033827e-05,2.275941e-08},
		{1.375249e-08,2.076669e-08,2.155749e-08,4.639285e-05,4.887124e-05,1.171827e-07,3.780572e-04,3.883957e-04,9.990407e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.482986e-09,5.272320e-10,5.832401e-05,3.920505e-05,2.149813e-08},
		{1.302396e-08,1.966892e-08,2.041942e-08,4.514426e-05,4.756551e-05,1.110015e-07,3.679465e-04,3.780231e-04,9.990663e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.247099e-09,5.000170e-10,5.676979e-05,3.816568e-05,2.036740e-08},
		{1.234745e-08,1.864695e-08,1.936011e-08,4.396249e-05,4.630975e-05,1.052320e-07,3.582688e-04,3.680545e-04,9.990909e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.022693e-09,4.734922e-10,5.525931e-05,3.714991e-05,1.930238e-08},
		{1.172071e-08,1.769942e-08,1.837285e-08,4.283428e-05,4.511906e-05,9.986464e-08,3.490205e-04,3.585638e-04,9.991143e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.821884e-09,4.495136e-10,5.384960e-05,3.619623e-05,1.832153e-08},
		{1.114211e-08,1.682562e-08,1.746719e-08,4.176302e-05,4.399480e-05,9.493623e-08,3.403060e-04,3.496083e-04,9.991364e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.633110e-09,4.276322e-10,5.250338e-05,3.529409e-05,1.741700e-08},
		{1.060832e-08,1.602082e-08,1.663446e-08,4.074890e-05,4.293011e-05,9.041922e-08,3.320868e-04,3.411778e-04,9.991573e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.457273e-09,4.068325e-10,5.123098e-05,3.444012e-05,1.658666e-08},
		{1.010812e-08,1.526520e-08,1.584923e-08,3.977822e-05,4.190621e-05,8.615260e-08,3.241627e-04,3.330285e-04,9.991774e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.293656e-09,3.879349e-10,5.000320e-05,3.362222e-05,1.580470e-08},
		{9.643201e-09,1.456243e-08,1.511862e-08,3.885229e-05,4.093343e-05,8.217406e-08,3.165967e-04,3.252691e-04,9.991966e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.145188e-09,3.701094e-10,4.885356e-05,3.283933e-05,1.507662e-08},
		{9.209815e-09,1.390778e-08,1.444003e-08,3.797378e-05,4.000000e-05,7.847882e-08,3.094053e-04,3.178613e-04,9.992149e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.001378e-09,3.531226e-10,4.773080e-05,3.208449e-05,1.439457e-08},
		{8.802619e-09,1.329202e-08,1.379790e-08,3.712474e-05,3.910819e-05,7.499777e-08,3.024713e-04,3.107387e-04,9.992325e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.867448e-09,3.378092e-10,4.665925e-05,3.137145e-05,1.375893e-08},
		{8.425447e-09,1.272443e-08,1.321087e-08,3.631935e-05,3.826353e-05,7.181637e-08,2.959838e-04,3.040866e-04,9.992489e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.746020e-09,3.232062e-10,4.565595e-05,3.069470e-05,1.317135e-08},
		{8.070538e-09,1.218732e-08,1.265270e-08,3.554511e-05,3.745190e-05,6.877531e-08,2.896675e-04,2.975953e-04,9.992649e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.633626e-09,3.100800e-10,4.469559e-05,3.004349e-05,1.261685e-08},
		{7.736879e-09,1.168473e-08,1.213020e-08,3.480605e-05,3.666701e-05,6.593976e-08,2.836222e-04,2.913772e-04,9.992803e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.522120e-09,2.970368e-10,4.374803e-05,2.941816e-05,1.209366e-08},
		{7.424498e-09,1.121246e-08,1.164012e-08,3.409921e-05,3.591749e-05,6.326827e-08,2.778250e-04,2.854239e-04,9.992950e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.419626e-09,2.847779e-10,4.285604e-05,2.881268e-05,1.160423e-08},
		{7.131195e-09,1.076903e-08,1.118215e-08,3.341484e-05,3.520598e-05,6.077853e-08,2.722976e-04,2.797430e-04,9.993090e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.324152e-09,2.737499e-10,4.200738e-05,2.824566e-05,1.115056e-08},
		{6.466698e-09,9.765902e-09,1.013761e-08,3.182437e-05,3.352395e-05,5.510363e-08,2.592859e-04,2.663804e-04,9.993420e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.109120e-09,2.481323e-10,4.000671e-05,2.689354e-05,1.010861e-08},
		{5.892689e-09,8.900096e-09,9.239731e-09,3.037751e-05,3.200540e-05,5.023179e-08,2.475527e-04,2.543296e-04,9.993718e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.921870e-09,2.261976e-10,3.819272e-05,2.567610e-05,9.213695e-09},
		{5.391736e-09,8.142499e-09,8.453785e-09,2.906079e-05,3.061303e-05,4.594690e-08,2.367653e-04,2.432408e-04,9.993992e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.757701e-09,2.068347e-10,3.652781e-05,2.455544e-05,8.427358e-09},
		{4.951367e-09,7.477460e-09,7.762790e-09,2.784905e-05,2.933791e-05,4.219410e-08,2.268919e-04,2.331043e-04,9.994242e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.614705e-09,1.899895e-10,3.501164e-05,2.353509e-05,7.741468e-09},
		{4.562844e-09,6.891025e-09,7.153801e-09,2.673592e-05,2.816312e-05,3.888433e-08,2.178187e-04,2.237794e-04,9.994473e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.487473e-09,1.750821e-10,3.360271e-05,2.259202e-05,7.131785e-09},
		{4.218520e-09,6.371180e-09,6.614464e-09,2.570694e-05,2.708163e-05,3.595405e-08,2.094517e-04,2.151814e-04,9.994685e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.375349e-09,1.619048e-10,3.231264e-05,2.172454e-05,6.594540e-09},
		{3.909825e-09,5.904593e-09,6.129613e-09,2.474754e-05,2.607410e-05,3.331972e-08,2.016400e-04,2.071620e-04,9.994883e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.276000e-09,1.502342e-10,3.111276e-05,2.091648e-05,6.111800e-09},
		{3.637561e-09,5.493598e-09,5.703221e-09,2.387099e-05,2.514990e-05,3.100002e-08,1.944859e-04,1.998150e-04,9.995065e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.186695e-09,1.396732e-10,3.001157e-05,2.017584e-05,5.686771e-09},
		{3.390955e-09,5.121287e-09,5.316938e-09,2.304733e-05,2.428389e-05,2.890127e-08,1.877878e-04,1.929316e-04,9.995235e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.106306e-09,1.302262e-10,2.897827e-05,1.948113e-05,5.301878e-09},
		{3.168281e-09,4.785132e-09,4.967519e-09,2.227966e-05,2.347189e-05,2.700327e-08,1.815216e-04,1.864896e-04,9.995394e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.033240e-09,1.216747e-10,2.800390e-05,1.883043e-05,4.952481e-09},
		{2.967178e-09,4.481405e-09,4.652269e-09,2.156045e-05,2.271635e-05,2.528968e-08,1.756685e-04,1.804820e-04,9.995542e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,9.679942e-10,1.139669e-10,2.710636e-05,1.822458e-05,4.639185e-09},
		{2.784762e-09,4.205539e-09,4.366098e-09,2.088871e-05,2.200595e-05,2.372963e-08,1.701668e-04,1.748262e-04,9.995682e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,9.081847e-10,1.069090e-10,2.625698e-05,1.765275e-05,4.352851e-09},
		{2.618560e-09,3.954837e-09,4.105837e-09,2.025547e-05,2.134001e-05,2.231864e-08,1.650267e-04,1.695454e-04,9.995812e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,8.539003e-10,1.005560e-10,2.546121e-05,1.711988e-05,4.093808e-09},
		{2.466591e-09,3.725398e-09,3.867182e-09,1.965988e-05,2.071134e-05,2.102108e-08,1.601620e-04,1.645505e-04,9.995936e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,8.047589e-10,9.469132e-11,2.471407e-05,1.661365e-05,3.855515e-09},
		{2.327704e-09,3.515423e-09,3.649522e-09,1.909735e-05,2.012136e-05,1.983791e-08,1.555939e-04,1.598556e-04,9.996052e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,7.596239e-10,8.943710e-11,2.400944e-05,1.614088e-05,3.638810e-09},
		{2.199725e-09,3.322181e-09,3.448731e-09,1.856606e-05,1.955979e-05,1.874686e-08,1.512529e-04,1.553954e-04,9.996162e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,7.175610e-10,8.444691e-11,2.333838e-05,1.569051e-05,3.438613e-09},
		{2.082755e-09,3.145569e-09,3.265566e-09,1.806589e-05,1.903276e-05,1.775035e-08,1.471801e-04,1.512095e-04,9.996265e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,6.794538e-10,8.000961e-11,2.270864e-05,1.526815e-05,3.255728e-09},
		{1.974748e-09,2.982534e-09,3.096430e-09,1.759101e-05,1.853327e-05,1.683157e-08,1.433158e-04,1.472423e-04,9.996363e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,6.440852e-10,7.581580e-11,2.211327e-05,1.486729e-05,3.087214e-09},
		{1.874553e-09,2.831190e-09,2.939148e-09,1.713960e-05,1.805690e-05,1.597690e-08,1.396330e-04,1.434580e-04,9.996457e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,6.115142e-10,7.197269e-11,2.154542e-05,1.448503e-05,2.930474e-09},
		{1.782063e-09,2.691524e-09,2.794232e-09,1.671117e-05,1.760628e-05,1.518950e-08,1.361489e-04,1.398789e-04,9.996545e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,5.812305e-10,6.842288e-11,2.100609e-05,1.412334e-05,2.785778e-09},
		{1.696167e-09,2.561717e-09,2.659410e-09,1.630391e-05,1.717673e-05,1.445492e-08,1.328170e-04,1.364543e-04,9.996630e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,5.533904e-10,6.514331e-11,2.049569e-05,1.377929e-05,2.651592e-09},
		{1.616345e-09,2.441177e-09,2.534224e-09,1.591577e-05,1.676771e-05,1.377525e-08,1.296569e-04,1.332072e-04,9.996710e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,5.272145e-10,6.207949e-11,2.000579e-05,1.345156e-05,2.526735e-09},
		{1.541942e-09,2.328802e-09,2.417557e-09,1.554533e-05,1.637756e-05,1.314090e-08,1.266391e-04,1.301071e-04,9.996786e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,5.030591e-10,5.922299e-11,1.954108e-05,1.313789e-05,2.410336e-09},
		{1.472762e-09,2.224375e-09,2.309228e-09,1.519243e-05,1.600635e-05,1.255244e-08,1.237676e-04,1.271592e-04,9.996859e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.804023e-10,5.653738e-11,1.909860e-05,1.283999e-05,2.302402e-09},
		{1.407890e-09,2.126281e-09,2.207259e-09,1.485452e-05,1.564963e-05,1.199805e-08,1.210064e-04,1.243198e-04,9.996929e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.591326e-10,5.407039e-11,1.867071e-05,1.255360e-05,2.200720e-09},
		{1.347428e-09,2.035114e-09,2.112751e-09,1.453236e-05,1.530977e-05,1.148419e-08,1.183843e-04,1.216271e-04,9.996996e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.394068e-10,5.171974e-11,1.826628e-05,1.228154e-05,2.106376e-09},
		{1.290722e-09,1.949381e-09,2.023870e-09,1.422269e-05,1.498499e-05,1.100106e-08,1.158698e-04,1.190428e-04,9.997060e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.209927e-10,4.957533e-11,1.787847e-05,1.202125e-05,2.017886e-09},
		{1.237388e-09,1.868854e-09,1.940027e-09,1.392645e-05,1.467166e-05,1.054552e-08,1.134472e-04,1.165547e-04,9.997121e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.037188e-10,4.751712e-11,1.750588e-05,1.176942e-05,1.934267e-09},
		{1.187451e-09,1.793526e-09,1.861909e-09,1.364211e-05,1.437320e-05,1.012166e-08,1.111424e-04,1.141874e-04,9.997180e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.874200e-10,4.560696e-11,1.714947e-05,1.153016e-05,1.856420e-09},
		{1.140467e-09,1.722465e-09,1.788194e-09,1.337006e-05,1.408549e-05,9.719767e-09,1.089145e-04,1.118971e-04,9.997236e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.720116e-10,4.378895e-11,1.680552e-05,1.129860e-05,1.782674e-09},
		{1.096145e-09,1.655523e-09,1.718644e-09,1.310767e-05,1.380933e-05,9.342004e-09,1.067770e-04,1.097026e-04,9.997290e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.576090e-10,4.209147e-11,1.647724e-05,1.107762e-05,1.713644e-09},
		{1.054360e-09,1.592446e-09,1.653144e-09,1.285575e-05,1.354341e-05,8.985994e-09,1.047241e-04,1.075926e-04,9.997343e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.439155e-10,4.048646e-11,1.615839e-05,1.086416e-05,1.648076e-09},
		{1.014944e-09,1.532931e-09,1.591400e-09,1.261299e-05,1.328822e-05,8.650519e-09,1.027507e-04,1.055646e-04,9.997393e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.310695e-10,3.897647e-11,1.585402e-05,1.065948e-05,1.586573e-09},
		{9.774562e-10,1.476264e-09,1.532517e-09,1.237755e-05,1.304096e-05,8.330584e-09,1.008343e-04,1.035972e-04,9.997441e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.189999e-10,3.755856e-11,1.555960e-05,1.046118e-05,1.527950e-09},
		{9.425023e-10,1.423494e-09,1.477771e-09,1.215433e-05,1.280561e-05,8.032741e-09,9.901298e-05,1.017267e-04,9.997487e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.075332e-10,3.620236e-11,1.527923e-05,1.027258e-05,1.473386e-09},
		{9.091375e-10,1.373117e-09,1.425507e-09,1.193711e-05,1.257723e-05,7.748770e-09,9.724710e-05,9.991193e-05,9.997532e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.966541e-10,3.492363e-11,1.500680e-05,1.008939e-05,1.421317e-09},
		{8.774666e-10,1.325303e-09,1.375809e-09,1.172782e-05,1.235584e-05,7.478809e-09,9.553917e-05,9.815604e-05,9.997576e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.862581e-10,3.370688e-11,1.474118e-05,9.911993e-06,1.371626e-09},
		{8.474736e-10,1.280003e-09,1.328790e-09,1.152544e-05,1.214324e-05,7.223234e-09,9.389261e-05,9.646598e-05,9.997617e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.765233e-10,3.255693e-11,1.448864e-05,9.741528e-06,1.324903e-09},
		{8.190135e-10,1.236959e-09,1.284140e-09,1.133066e-05,1.193724e-05,6.979816e-09,9.229765e-05,9.482632e-05,9.997658e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.671885e-10,3.145528e-11,1.424238e-05,9.575724e-06,1.280232e-09},
		{7.919457e-10,1.196126e-09,1.241752e-09,1.114172e-05,1.173853e-05,6.749999e-09,9.076428e-05,9.325092e-05,9.997697e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.583431e-10,3.042010e-11,1.400496e-05,9.416729e-06,1.238045e-09},
		{7.416560e-10,1.120146e-09,1.162853e-09,1.078206e-05,1.136007e-05,6.321059e-09,8.783555e-05,9.024231e-05,9.997771e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.420377e-10,2.849718e-11,1.355415e-05,9.112710e-06,1.159379e-09},
		{6.960101e-10,1.051219e-09,1.091300e-09,1.044536e-05,1.100471e-05,5.932020e-09,8.508958e-05,8.742044e-05,9.997841e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.270961e-10,2.674022e-11,1.312957e-05,8.827828e-06,1.087993e-09},
		{6.544568e-10,9.884719e-10,1.026151e-09,1.012886e-05,1.067125e-05,5.578055e-09,8.251128e-05,8.477226e-05,9.997906e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.135331e-10,2.513568e-11,1.273213e-05,8.560241e-06,1.023080e-09},
		{6.165358e-10,9.311853e-10,9.666837e-10,9.831068e-06,1.035760e-05,5.254498e-09,8.008278e-05,8.227683e-05,9.997968e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.011746e-10,2.368335e-11,1.235810e-05,8.308751e-06,9.637979e-10},
		{5.817850e-10,8.787009e-10,9.121860e-10,9.550064e-06,1.006154e-05,4.958393e-09,7.779449e-05,7.992580e-05,9.998026e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.898316e-10,2.234899e-11,1.200450e-05,8.071153e-06,9.094377e-10},
		{5.499052e-10,8.305355e-10,8.621709e-10,9.284799e-06,9.782008e-06,4.686581e-09,7.563183e-05,7.770371e-05,9.998081e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.793820e-10,2.112329e-11,1.167029e-05,7.846781e-06,8.595872e-10},
		{5.206101e-10,7.863033e-10,8.163173e-10,9.033951e-06,9.518128e-06,4.437294e-09,7.359272e-05,7.560882e-05,9.998133e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.698424e-10,1.999911e-11,1.135576e-05,7.635432e-06,8.138747e-10},
		{4.935566e-10,7.454733e-10,7.738893e-10,8.796173e-06,9.267571e-06,4.206922e-09,7.165689e-05,7.362054e-05,9.998182e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.610454e-10,1.895929e-11,1.105721e-05,7.434313e-06,7.715787e-10},
		{4.685721e-10,7.077098e-10,7.346866e-10,8.570847e-06,9.029834e-06,3.993580e-09,6.981663e-05,7.172991e-05,9.998228e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.528867e-10,1.799705e-11,1.077381e-05,7.243520e-06,7.325122e-10},
		{4.454275e-10,6.727667e-10,6.984180e-10,8.356543e-06,8.804071e-06,3.796481e-09,6.807250e-05,6.993740e-05,9.998273e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.453206e-10,1.710891e-11,1.050383e-05,7.062367e-06,6.962944e-10},
		{4.239744e-10,6.403566e-10,6.647671e-10,8.152713e-06,8.589589e-06,3.613525e-09,6.641160e-05,6.823195e-05,9.998315e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.383492e-10,1.628718e-11,1.024840e-05,6.890425e-06,6.627753e-10},
		{4.040090e-10,6.102123e-10,6.334649e-10,7.958571e-06,8.384861e-06,3.443470e-09,6.483051e-05,6.660676e-05,9.998355e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.318170e-10,1.552097e-11,1.000350e-05,6.726296e-06,6.315363e-10},
		{3.854460e-10,5.821559e-10,6.043532e-10,7.773664e-06,8.189963e-06,3.284998e-09,6.332145e-05,6.505662e-05,9.998393e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.257600e-10,1.480579e-11,9.771271e-06,6.569737e-06,6.025119e-10},
		{3.681165e-10,5.560019e-10,5.771769e-10,7.596952e-06,8.003769e-06,3.137447e-09,6.188300e-05,6.357917e-05,9.998430e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.201241e-10,1.413979e-11,9.549412e-06,6.420292e-06,5.754224e-10},
		{3.519129e-10,5.315176e-10,5.517671e-10,7.427851e-06,7.825724e-06,2.999347e-09,6.050599e-05,6.216409e-05,9.998465e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.148279e-10,1.351725e-11,9.336755e-06,6.277563e-06,5.501084e-10},
		{3.368134e-10,5.087201e-10,5.281210e-10,7.266706e-06,7.656062e-06,2.870793e-09,5.919454e-05,6.081694e-05,9.998498e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.098941e-10,1.293751e-11,9.134302e-06,6.141508e-06,5.265250e-10},
		{3.226238e-10,4.872887e-10,5.058670e-10,7.112025e-06,7.493088e-06,2.749870e-09,5.793491e-05,5.952274e-05,9.998530e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.052630e-10,1.239271e-11,8.939670e-06,6.010718e-06,5.043282e-10},
		{3.093193e-10,4.671887e-10,4.849919e-10,6.963910e-06,7.336926e-06,2.636324e-09,5.672628e-05,5.828061e-05,9.998561e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.009242e-10,1.188317e-11,8.753406e-06,5.885624e-06,4.835317e-10},
		{2.968256e-10,4.483228e-10,4.654135e-10,6.821796e-06,7.187336e-06,2.529916e-09,5.556942e-05,5.709255e-05,9.998590e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,9.684910e-11,1.140065e-11,8.575184e-06,5.765471e-06,4.640125e-10},
		{2.850674e-10,4.305666e-10,4.469810e-10,6.685446e-06,7.043436e-06,2.429695e-09,5.445758e-05,5.594995e-05,9.998618e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,9.300110e-11,1.094831e-11,8.403262e-06,5.650111e-06,4.456233e-10},
		{2.739890e-10,4.138275e-10,4.295913e-10,6.554260e-06,6.905265e-06,2.335188e-09,5.338880e-05,5.485186e-05,9.998645e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,8.940666e-11,1.052492e-11,8.238632e-06,5.539225e-06,4.282964e-10},
		{2.635630e-10,3.980791e-10,4.132575e-10,6.428359e-06,6.772615e-06,2.246352e-09,5.236325e-05,5.379802e-05,9.998671e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,8.599480e-11,1.012366e-11,8.080157e-06,5.432680e-06,4.119888e-10},
		{2.537041e-10,3.831916e-10,3.977943e-10,6.307042e-06,6.644737e-06,2.162333e-09,5.137493e-05,5.278276e-05,9.998696e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,8.277661e-11,9.745327e-12,7.927506e-06,5.330194e-06,3.965772e-10},
		{2.443641e-10,3.690805e-10,3.831435e-10,6.189723e-06,6.521461e-06,2.082728e-09,5.042075e-05,5.180265e-05,9.998721e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,7.974997e-11,9.389639e-12,7.780604e-06,5.231323e-06,3.819875e-10},
		{2.355922e-10,3.558363e-10,3.694035e-10,6.077584e-06,6.403394e-06,2.008022e-09,4.950757e-05,5.086450e-05,9.998744e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,7.688069e-11,9.051291e-12,7.639876e-06,5.136667e-06,3.682969e-10},
		{2.272475e-10,3.432351e-10,3.563151e-10,5.969034e-06,6.288961e-06,1.936903e-09,4.862320e-05,4.995601e-05,9.998766e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,7.415612e-11,8.730959e-12,7.503266e-06,5.044951e-06,3.552507e-10},
		{2.193485e-10,3.313030e-10,3.439336e-10,5.864433e-06,6.178638e-06,1.869583e-09,4.777055e-05,4.907972e-05,9.998788e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,7.156933e-11,8.426926e-12,7.371424e-06,4.956416e-06,3.428955e-10},
		{2.118473e-10,3.199693e-10,3.321647e-10,5.763245e-06,6.072159e-06,1.805601e-09,4.694668e-05,4.823336e-05,9.998809e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,6.913681e-11,8.140066e-12,7.244599e-06,4.870884e-06,3.311621e-10},
		{2.047246e-10,3.092134e-10,3.209992e-10,5.665625e-06,5.969129e-06,1.744892e-09,4.615065e-05,4.741531e-05,9.998829e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,6.680507e-11,7.865889e-12,7.121509e-06,4.788289e-06,3.200230e-10},
		{1.979559e-10,2.989921e-10,3.103869e-10,5.571201e-06,5.869643e-06,1.687233e-09,4.538151e-05,4.662532e-05,9.998848e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,6.459544e-11,7.604473e-12,7.002941e-06,4.708456e-06,3.094494e-10},
		{1.915209e-10,2.892704e-10,3.002952e-10,5.479900e-06,5.773475e-06,1.632324e-09,4.463702e-05,4.586030e-05,9.998867e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,6.249792e-11,7.357952e-12,6.888271e-06,4.631351e-06,2.993896e-10},
		{1.853886e-10,2.800086e-10,2.906782e-10,5.391473e-06,5.680308e-06,1.580070e-09,4.391699e-05,4.512052e-05,9.998886e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,6.049605e-11,7.122444e-12,6.777004e-06,4.556582e-06,2.897972e-10},
		{1.795488e-10,2.711852e-10,2.815159e-10,5.305894e-06,5.590123e-06,1.530276e-09,4.321933e-05,4.440367e-05,9.998903e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,5.858147e-11,6.897876e-12,6.669183e-06,4.484194e-06,2.806654e-10},
		{1.739873e-10,2.627876e-10,2.728107e-10,5.223016e-06,5.502938e-06,1.482948e-09,4.254568e-05,4.371158e-05,9.998920e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,5.676996e-11,6.684418e-12,6.565265e-06,4.414367e-06,2.719873e-10},
		{1.686726e-10,2.547663e-10,2.644757e-10,5.142644e-06,5.418246e-06,1.437691e-09,4.189141e-05,4.303958e-05,9.998937e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,5.504145e-11,6.480106e-12,6.464357e-06,4.346388e-06,2.636785e-10},
		{1.636005e-10,2.470999e-10,2.565171e-10,5.064796e-06,5.336110e-06,1.394381e-09,4.125569e-05,4.238645e-05,9.998953e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,5.338496e-11,6.284721e-12,6.366452e-06,4.280472e-06,2.557480e-10},
		{1.587517e-10,2.397790e-10,2.489187e-10,4.989188e-06,5.256456e-06,1.353087e-09,4.064038e-05,4.175407e-05,9.998969e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,5.179957e-11,6.098606e-12,6.271202e-06,4.216558e-06,2.481605e-10},
		{1.541192e-10,2.327804e-10,2.416522e-10,4.915810e-06,5.179243e-06,1.313580e-09,4.004251e-05,4.114015e-05,9.998984e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,5.029399e-11,5.921133e-12,6.179251e-06,4.154659e-06,2.409226e-10},
		{1.496808e-10,2.260790e-10,2.346935e-10,4.844552e-06,5.104098e-06,1.275775e-09,3.946223e-05,4.054369e-05,9.998999e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.884169e-11,5.750778e-12,6.089365e-06,4.094418e-06,2.339781e-10},
		{1.454371e-10,2.196649e-10,2.280381e-10,4.775408e-06,5.031216e-06,1.239545e-09,3.889797e-05,3.996405e-05,9.999013e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.745678e-11,5.587258e-12,6.002525e-06,4.035874e-06,2.273427e-10},
		{1.413673e-10,2.135225e-10,2.216557e-10,4.708132e-06,4.960328e-06,1.204892e-09,3.835036e-05,3.940160e-05,9.999027e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.613306e-11,5.430815e-12,5.918063e-06,3.978977e-06,2.209807e-10},
		{1.374618e-10,2.076210e-10,2.155311e-10,4.642630e-06,4.891360e-06,1.171605e-09,3.781699e-05,3.885348e-05,9.999040e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.485641e-11,5.280737e-12,5.835691e-06,3.923689e-06,2.148793e-10},
		{1.337302e-10,2.019870e-10,2.096874e-10,4.579160e-06,4.824535e-06,1.139835e-09,3.730049e-05,3.832291e-05,9.999053e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.363699e-11,5.137438e-12,5.755954e-06,3.870101e-06,2.090512e-10},
		{1.301380e-10,1.965613e-10,2.040536e-10,4.517249e-06,4.759304e-06,1.109222e-09,3.679636e-05,3.780495e-05,9.999066e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.246440e-11,4.999490e-12,5.678054e-06,3.817753e-06,2.034319e-10},
		{1.266904e-10,1.913527e-10,1.986441e-10,4.457039e-06,4.695821e-06,1.079799e-09,3.630507e-05,3.730004e-05,9.999079e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.133993e-11,4.867433e-12,5.602336e-06,3.766901e-06,1.980430e-10},
		{1.233792e-10,1.863526e-10,1.934552e-10,4.398394e-06,4.634091e-06,1.051597e-09,3.582771e-05,3.680980e-05,9.999091e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.025977e-11,4.739531e-12,5.528793e-06,3.717319e-06,1.928699e-10},
		{1.201944e-10,1.815431e-10,1.884626e-10,4.341307e-06,4.573842e-06,1.024451e-09,3.536222e-05,3.633143e-05,9.999103e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.921735e-11,4.616981e-12,5.456814e-06,3.669020e-06,1.878891e-10},
		{1.171298e-10,1.769126e-10,1.836521e-10,4.285606e-06,4.515173e-06,9.983073e-10,3.490844e-05,3.586520e-05,9.999114e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.822291e-11,4.499837e-12,5.386912e-06,3.621937e-06,1.830960e-10},
		{1.141864e-10,1.724664e-10,1.790408e-10,4.231417e-06,4.458079e-06,9.732270e-10,3.446705e-05,3.541164e-05,9.999125e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.725959e-11,4.386535e-12,5.318705e-06,3.576079e-06,1.784921e-10},
		{1.113476e-10,1.681797e-10,1.745883e-10,4.178512e-06,4.402310e-06,9.490338e-10,3.403607e-05,3.496890e-05,9.999136e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.633283e-11,4.277578e-12,5.252134e-06,3.531383e-06,1.740544e-10},
		{1.086062e-10,1.640379e-10,1.702881e-10,4.126692e-06,4.347858e-06,9.256689e-10,3.361463e-05,3.453600e-05,9.999147e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.544443e-11,4.173173e-12,5.187248e-06,3.487712e-06,1.697722e-10},
		{1.059835e-10,1.600779e-10,1.661799e-10,4.076545e-06,4.295060e-06,9.033322e-10,3.320628e-05,3.411649e-05,9.999157e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.458652e-11,4.072006e-12,5.124319e-06,3.445383e-06,1.656789e-10},
		{1.034434e-10,1.562420e-10,1.621957e-10,4.027418e-06,4.243272e-06,8.816831e-10,3.280607e-05,3.370536e-05,9.999168e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.375711e-11,3.974489e-12,5.062496e-06,3.403876e-06,1.617077e-10},
		{1.009966e-10,1.525459e-10,1.583604e-10,3.979526e-06,4.192760e-06,8.608301e-10,3.241570e-05,3.330416e-05,9.999177e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.295577e-11,3.880306e-12,5.002133e-06,3.363341e-06,1.578807e-10},
		{9.863337e-11,1.489753e-10,1.546527e-10,3.932672e-06,4.143455e-06,8.406728e-10,3.203423e-05,3.291228e-05,9.999187e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.218937e-11,3.789928e-12,4.943408e-06,3.323736e-06,1.541842e-10},
		{9.635221e-11,1.455305e-10,1.510768e-10,3.886971e-06,4.095221e-06,8.212295e-10,3.166158e-05,3.252932e-05,9.999197e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.144257e-11,3.702117e-12,4.885785e-06,3.285070e-06,1.506167e-10},
		{9.414954e-11,1.422042e-10,1.476232e-10,3.842294e-06,4.048148e-06,8.024657e-10,3.129768e-05,3.215555e-05,9.999206e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.072346e-11,3.617039e-12,4.829685e-06,3.247296e-06,1.471758e-10},
		{9.202266e-11,1.389910e-10,1.442878e-10,3.798645e-06,4.002172e-06,7.843173e-10,3.094177e-05,3.178982e-05,9.999215e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,3.003016e-11,3.535553e-12,4.774865e-06,3.210433e-06,1.438509e-10},
		{8.996505e-11,1.358832e-10,1.410608e-10,3.755942e-06,3.957180e-06,7.667835e-10,3.059406e-05,3.143258e-05,9.999224e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.935842e-11,3.456525e-12,4.721134e-06,3.174327e-06,1.406323e-10},
		{8.797657e-11,1.328788e-10,1.379410e-10,3.714208e-06,3.913200e-06,7.498280e-10,3.025385e-05,3.108301e-05,9.999232e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.870644e-11,3.380052e-12,4.668553e-06,3.139026e-06,1.375229e-10},
		{8.605560e-11,1.299782e-10,1.349342e-10,3.673405e-06,3.870276e-06,7.334807e-10,2.992220e-05,3.074228e-05,9.999241e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.808067e-11,3.306318e-12,4.617392e-06,3.104648e-06,1.345254e-10},
		{8.419416e-11,1.271688e-10,1.320149e-10,3.633467e-06,3.828191e-06,7.176298e-10,2.959710e-05,3.040836e-05,9.999249e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.747523e-11,3.234754e-12,4.567250e-06,3.070868e-06,1.316154e-10},
		{8.239348e-11,1.244470e-10,1.291894e-10,3.594433e-06,3.787006e-06,7.022537e-10,2.927835e-05,3.008088e-05,9.999257e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.688711e-11,3.165382e-12,4.518159e-06,3.037816e-06,1.288000e-10},
		{8.064916e-11,1.218134e-10,1.264560e-10,3.556186e-06,3.746712e-06,6.873986e-10,2.896711e-05,2.976100e-05,9.999265e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.631675e-11,3.098422e-12,4.469977e-06,3.005485e-06,1.260705e-10},
		{7.896111e-11,1.192631e-10,1.238081e-10,3.518748e-06,3.707317e-06,6.730022e-10,2.866208e-05,2.944778e-05,9.999273e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.576809e-11,3.033747e-12,4.423064e-06,2.973904e-06,1.234331e-10},
		{7.732336e-11,1.167903e-10,1.212403e-10,3.482086e-06,3.668656e-06,6.590523e-10,2.836354e-05,2.914091e-05,9.999280e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.523217e-11,2.970883e-12,4.376819e-06,2.942911e-06,1.208706e-10},
		{7.573819e-11,1.143944e-10,1.187543e-10,3.446221e-06,3.630851e-06,6.455191e-10,2.807086e-05,2.884026e-05,9.999288e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.471482e-11,2.909798e-12,4.331775e-06,2.912544e-06,1.183920e-10},
		{7.419971e-11,1.120724e-10,1.163416e-10,3.411046e-06,3.593787e-06,6.324191e-10,2.778455e-05,2.854619e-05,9.999295e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.421441e-11,2.850655e-12,4.287614e-06,2.882795e-06,1.159870e-10},
		{7.270605e-11,1.098153e-10,1.139992e-10,3.376532e-06,3.557446e-06,6.196883e-10,2.750352e-05,2.825738e-05,9.999302e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.372613e-11,2.793252e-12,4.244211e-06,2.853663e-06,1.136535e-10},
		{7.126261e-11,1.076360e-10,1.117387e-10,3.342834e-06,3.521966e-06,6.073990e-10,2.722930e-05,2.797570e-05,9.999309e-01,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.325441e-11,2.737816e-12,4.201877e-06,2.825211e-06,1.113991e-10}};
	
	for(i=0;i<NumInArrays;i++) {
		if(psz[i]>AVE) {
			if(i==0) x=0;
			else x=i-1;
			break;
		}
	}	
	if(x==-1) x=NumInArrays-1;  /* set it to the pi vector for the largest pop size we investigated */
	
	ret = (double *)calloc(NUM_SPEC_PEDS,sizeof(double));
	for(i=0;i<NUM_SPEC_PEDS;i++)  {
		ret[i] = pivec[x][i];
	}
	*outx = x;
	*outpsz = psz[x];
	return(ret);
}




/* after parsing the command line and reading all the data, this function
 decides what to include for the Pi vectors.  Recall that the final resting
 place for these variables is in the PBUO->Pi arrays.  Here is how we deal 
 with multiple definitions:
   -If the PBUO->Pi[i] vector was set on the command line (using --Pi), it remains untouched
   -If not, then:
     -If the --psz-for-all option was invoked, then we use that, even if it has a value given from the data file
     -if the --psz-for-all option was not invoked then:
       -we use the data from the data file if given
       -if value is not given in data file, then we use the standard default of whatever it will be (probably 1000). 
 */
void NegotiatePiVectors(pbt_user_opts *PBUO, pfr_geno_data *PFR) 
{
	int i,j,N = PFR->NumPops;
	int x,psz;
	int DefaultChinookPopSizeValue = 1000;
	FILE *out;
	
	if( (out=fopen("snppit_output_PopSizesAnPiVectors.txt","w"))==NULL) {
		fprintf(stderr,"Error! Unable to open file \"snppit_output_PopSizesAnPiVectors.txt\" to write to in from function NegotiatePiVectors(). Exiting!\n");
		exit(1);
	}
	fprintf(out, "PopName    PopSizeSetBy    GivenPopSize    PopSizeUsedForPi     PiVectorElements....\n");
	
	for(i=0;i<N;i++)  {
		fprintf(out,"%s   ",PFR->PopNames[i]);
		if(PBUO->Pi[i]) {
			fprintf(out,"--Pi    XXXX   XXXX  ");
		}
		else {
			if(PBUO->PszForAll>0) {
				PBUO->Pi[i] = Return_ArchTypeChinookPiVec(PBUO->PszForAll, &x, &psz);
				fprintf(out,"--psz-for-all   %d   %d  ",PBUO->PszForAll,psz);
			}
			else if (PFR->AveragePopSizes==NULL || PFR->AveragePopSizes[i]<=0) {
				PBUO->Pi[i] = Return_ArchTypeChinookPiVec(DefaultChinookPopSizeValue, &x, &psz);
				fprintf(out,"Default   %d   %d  ",DefaultChinookPopSizeValue,psz);
			}
			else {  
 				PBUO->Pi[i] = Return_ArchTypeChinookPiVec(PFR->AveragePopSizes[i], &x, &psz);
				fprintf(out,"InDataFile   %d   %d  ",PFR->AveragePopSizes[i],psz);
			}
			
		}
		for(j=0;j<NUM_SPEC_PEDS;j++)  {
			fprintf(out, "%.3e  ", PBUO->Pi[i][j]);
		}
		fprintf(out,"\n");
	}
	fclose(out);
}





void PrintFinalIndivReportWithFDRs(pbt_high_level_vars *HLV)
{
	int k,j,i,p;
	pfr_offspring *Offs;
	FILE *out;
	int intsct;

	
	if( (out=fopen("snppit_output_ParentageAssignments.txt","w"))==NULL)  {
		fprintf(stderr,"Error! Failed to open file \"snppit_output_ParentageAssignments.txt\" to write to it in function PrintFinalIndivReportWithFDRs().  Exiting\n");
		exit(1);
	}
	
	fprintf(out,"OffspCollection\tKid\tPa\tMa\tPopName\tSpawnYear\tFDR\tPvalue\tLOD\tP.Pr.C_Se_Se\tP.Pr.Max\tMaxP.Pr.Relat\tTotPaNonExc\tTotMaNonExc\tTotUnkNonExc\tTotPairsMendCompat\tTotPairsMendAndLogL\tTotParsMendLoglAndRank\tTotPairsNonExc\tKidMiss\tPaMiss\tMaMiss\tMI.Kid.Pa\tMI.Kid.Ma\tMI.Trio\tMendIncLoci\n");

	for(j=0;j<HLV->PFR->NumOffColls;j++)  {
		Offs = HLV->PFR->Offs[j];
		
		for(i=0;i<HLV->PFR->NumInOffColls[j];i++)  {
			if(Offs[i].MendComps->NumParPairs>0) {
				/* get the first intersecting year (note, with iteroparous individuals we will have to modify this */
				MateTimingWorksOut(Offs[i].MendComps->ParPairs[0].mamapapa[MALE]->ReproYears, Offs[i].MendComps->ParPairs[0].mamapapa[FEMALE]->ReproYears, &intsct);
				fprintf(out,"%s\t%s\t%s\t%s\t%s\t%d\t%.6f\t%.6f\t%f\t%f\t%f\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t",
						HLV->PFR->OffCollNames[j],
						Offs[i].Name,
						Offs[i].MendComps->ParPairs[0].mamapapa[MALE]->Name,
						Offs[i].MendComps->ParPairs[0].mamapapa[FEMALE]->Name,
						HLV->PFR->PopNames[Offs[i].MendComps->ParPairs[0].mamapapa[FEMALE]->Pop],
						intsct,
						Offs[i].FDR,
						Offs[i].p_value,
						Offs[i].MendComps->ParPairs[0].LogL,
						Offs[i].MendComps->ParPairs[0].TrioPosteriors[0],
						Offs[i].MendComps->ParPairs[0].MaxPosterior,
						SpecPedIdx2Str(Offs[i].MendComps->ParPairs[0].MaxPosteriorRelat),
						Offs[i].MendComps->NumParents[MALE],
						Offs[i].MendComps->NumParents[FEMALE],
						Offs[i].MendComps->NumParents[SEX_UNKNOWN],
            Offs[i].MendComps->NumParPairsMendelian,
            Offs[i].MendComps->NumParPairsMend_LogL,
            Offs[i].MendComps->NumParPairsMend_LogL_and_Rank,
						Offs[i].MendComps->NumParPairs,
						Offs[i].NumMissingLoci[0],
						Offs[i].MendComps->ParPairs[0].mamapapa[MALE]->NumMissingLoci[0],
						Offs[i].MendComps->ParPairs[0].mamapapa[FEMALE]->NumMissingLoci[0],
						Offs[i].MendComps->ParPairs[0].pa_MI,
						Offs[i].MendComps->ParPairs[0].ma_MI,
						Offs[i].MendComps->ParPairs[0].trio_MI
						);
				p=0;
				for(k=0;k<HLV->PFR->NumLoci;k++)  {
					if(Offs[i].MendComps->ParPairs[0].X_states_64[k] > 0) {
						if(p>0) {
							fprintf(out,",");
						}
						p++;
						fprintf(out,"%s",HLV->PFR->LocusNames[k]);
					}
				}
				if(p==0) { fprintf(out,"---"); }
				fprintf(out,"\n");
			}
			else {
				fprintf(out,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t---\n",
						HLV->PFR->OffCollNames[j],
						Offs[i].Name,
						"---",
						"---",
						"---",
						"---",
						"---",
						"---",
						"---",
						"---",
						"---",
						"---",
						Offs[i].MendComps->NumParents[MALE],
						Offs[i].MendComps->NumParents[FEMALE],
						Offs[i].MendComps->NumParents[SEX_UNKNOWN],
            Offs[i].MendComps->NumParPairsMendelian,
            Offs[i].MendComps->NumParPairsMend_LogL,
            Offs[i].MendComps->NumParPairsMend_LogL_and_Rank,
						Offs[i].MendComps->NumParPairs,
						Offs[i].NumMissingLoci[0],
						"---",
						"---",
						"---",
						"---",
						"---"
						);
			}
			
		}
		
		
	}
	
	
}


