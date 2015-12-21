
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
#include "pbt_trio_mixture.h"




int main(int argc, char *argv[])
{
	int i;
	pbt_high_level_vars *HLV = (pbt_high_level_vars *)malloc(sizeof(pbt_high_level_vars));
	FB_Vars  *CPFB;  /* "Compat Pairs FB_vars"  this is going to be heavily re-used space, each time we are dealing with
									an offspring and all of his non-excluded parent pairs.  It will contain 
									a Colls array that is linear with the Specified pedigrees.  */
	 
	
	
	
	HLV->PBUO = GetPBT_UserOpts(argc, argv);
	
	
	/* open up a file stream for the basic summary data  */
	if( (HLV->BasicSummaries_File = fopen("snppit_output_BasicDataSummary.txt","w"))==NULL) {
		fprintf(stderr,"Error! Failed to open file snppit_output_BasicDataSummary.txt to write to it.  You may have it open and locked in another application. Exiting...\n");
		exit(1);
	}
	printf("\n\n");
	
	
	HLV->PFR = FirstPassThroughData(HLV->PBUO->DataFileName);
	
	#ifdef VERBOSE_PRINT_FIRST_PASS_SUMMARY
		PrintFirstPassSummaryOfPopsCollsAndIndivs(HLV->PFR);
	#endif
	CollectDataOnSecondPass(HLV->PFR, HLV->PBUO->DataFileName);
	NegotiatePiVectors(HLV->PBUO, HLV->PFR);
	PrintSummaryOfInputData(HLV->BasicSummaries_File,HLV->PFR);
	SummarizeLocusNameAndGtypRates(HLV->BasicSummaries_File, HLV->PFR);
	SummarizeAllelicTypes(HLV->BasicSummaries_File,HLV->PFR);
	CountAlleles(HLV->PFR);
	ComputeAlleleFreqsFromCounts(HLV->BasicSummaries_File, HLV->PFR);
	PrintSummaryOfAllelicCountsAndFreqs(HLV->BasicSummaries_File,HLV->PFR);
	fflush(stdout);
	fclose(HLV->BasicSummaries_File);
	
	
	
	
	
	if(HLV->PBUO->DryRun>0) {
		printf("\n\nData have been read in and summaries compiled on this dry run.\n");
		printf("Please check the data summaries in file \"snppit_output_BasicDataSummary.txt\"\n");
		printf("to confirm that the program is running correctly. If it all looks good,\n");
		printf("then try a full-blown run by removing the --dry-run option from the\n");
		printf("command line.\n\n");
		
		return(0);
	}
	else {
		printf("\n\n\nDATA HAVE BEEN READ.  SUMMARIES APPEAR IN:  snppit_output_BasicDataSummary.txt\n\n\n");
	}
	

	
	
	
	/* now compute the parental trios forwards probs using the Big Smax and select the smax to use 
	   for the future analyses */
	printf("COMPUTING AN APPROPRIATE S-MAX\n");
	SelectAnSmax3(HLV);
	
	
	
	printf("\n\n");
	for(i=0;i<HLV->PFR->NumOffColls;i++) {
		#ifdef VERBOSE_SINGLE_PARENT_COMPAT_WITH_OFFSPRING
			printf("EXCLUDING SINGLE PARENTS.  COLLECTION %d  %s   with %d indivs in collection. \nDone with individual index:\n",i+1,HLV->PFR->OffCollNames[i],HLV->PFR->NumInOffColls[i]);
		#endif
		AssignMatchingSingleParents(HLV, i);
	}
	
	
	
	printf("\n\n");
	for(i=0;i<HLV->PFR->NumOffColls;i++) {
		printf("FINDING NON EXCLUDED PARENT PAIRS.  COLLECTION %d  %s   with %d indivs in collection. \nDone with individual index:\n",i+1,HLV->PFR->OffCollNames[i], HLV->PFR->NumInOffColls[i]);
		fflush(stdout);
		AssignMatchingParentPairs(HLV, i, HLV->smax_to_use[2]);
		
	}
	
	
	printf("\n\n");
	printf("COMPUTING THE FORWARD STEP AND PREPARING FOR BACKWARD STEP FOR ALL POPULATIONS\n");
	fflush(stdout);
	ComputePurePopTrioColls(HLV);  /* this should do the forward step AND prepare for the backward step in all these */
	//ComputeCrossPopTrioColls(HLV); 
	
	
	/* this can only be done AFTER computing all the PurePopTrioColls */
	/* open up a file stream for the posteriors */
	if( (HLV->TrioPosteriors_File = fopen("snppit_output_TrioPosteriors.txt","w"))==NULL) {
		fprintf(stderr,"Error! Failed to open file snppit_output_TrioPosteriors.txt to write to it.  You may have it open and locked in another application. Exiting...\n");
		exit(1);
	}
	fprintf(HLV->TrioPosteriors_File,"OffspCollection\tKid\tPa\tMa\tRank\tLOD");
	{int k;
		for(k=0;k<NUM_SPEC_PEDS;k++)  {
			fprintf(HLV->TrioPosteriors_File,"\tP.Pr.%s",SpecPedIdx2Str(k));
		}
	}
	fprintf(HLV->TrioPosteriors_File,"\tKidMiss\tPaMiss\tMaMiss\tMI.Kid.Pa\tMI.Kid.Ma\tMI.Trio\n");
	
	
	
	printf("\n\n");
	for(i=0;i<HLV->PFR->NumOffColls;i++) {
		printf("COMPUTING POSTERIORS:  COLLECTION %d  %s      with %d indivs in collection. \nDone with individual index:\n",i+1,HLV->PFR->OffCollNames[i],HLV->PFR->NumInOffColls[i]);
		fflush(stdout);
		CalculateTrioPosteriorProbs(HLV, i);
	}
	fclose(HLV->TrioPosteriors_File);
	
	
	HLV->KidsWithMaxLOD_Parents = (inds_with_max_lod_parents_from_this_pop **)calloc(HLV->PFR->NumPops,sizeof(inds_with_max_lod_parents_from_this_pop *));
	for(i=0;i<HLV->PFR->NumPops;i++)  {
		HLV->KidsWithMaxLOD_Parents[i] = RecordNonExcParentPairsFromPop(i,HLV);
	}
	
	
			
	/* allocate space to the areas where we will do FB algorithm successively, for each offspring with compatible parents */
	CPFB = (FB_Vars *)malloc(sizeof(FB_Vars));
	CPFB->RP = HLV->PurePopTrioColls->RP;
	CPFB->NumColls = NUM_SPEC_PEDS;
	CPFB->Colls = (Collection **)calloc(CPFB->NumColls, sizeof(Collection *));
	for(i=0;i<CPFB->NumColls;i++)  {
		CPFB->Colls[i] = AllocToCollection(CPFB->RP);
	}
	
	
	
	/* open up a file stream where we will store the max LOD parents */
	/*if( (HLV->MaxLodNonExpPar_File = fopen("snppit_output_MaxLodNonExParents.txt","w"))==NULL) {
		fprintf(stderr,"Error! Failed to open file snppit_output_MaxLodNonExParents.txt to write to it.  You may have it open and locked in another application. Exiting...\n");
		exit(1);
	}
	 */
	/*fprintf(HLV->MaxLodNonExpPar_File,"Kid\tPa\tMa\tPvalue\tLOD\tP.Pr.C_Se_Se\tP.Pr.Max\tMaxP.Pr.Relat\tTotPaNonExc\tTotMaNonExc\tTotUnkNonExc\tTotPairsNonExc\tKidMiss\tPaMiss\tMaMiss\tMI.Kid.Pa\tMI.Kid.Ma\tMI.Trio\tMendIncLoci\n");
	*/
	printf("\n\n");
	
	SeedFromFile("snppit_seeds");
	
	printf("\n\n");
	for(i=0;i<HLV->PFR->NumOffColls;i++) {
		printf("COMPUTING P-VALUES BY SIMULATION:  COLLECTION %d  %s    with %d indivs in collection\nDone with individual index:\n",i+1,HLV->PFR->OffCollNames[i],HLV->PFR->NumInOffColls[i]);
		AssessFPR_and_FNR_ByBackwardSimulation(HLV, i, CPFB);
	}
	/*fclose(HLV->MaxLodNonExpPar_File);*/
	
	
	if( (HLV->FDR_Summaries_File = fopen("snppit_output_FDR_Summary.txt","w"))==NULL) {
		fprintf(stderr,"Error! Failed to open file \"snppit_output_FDR_Summary.txt\" to write to it.  You may have it open and locked in another application. Exiting...\n");
		exit(1);
	}
	printf("\n\n");
	printf("PERFORMING FALSE DISCOVERY RATE CORRECTIONS\n");
	fprintf(HLV->FDR_Summaries_File,"PopName\tRankInFDR\tKid\tPa\tMa\tFDR\tFDC.est.to.pop\tPvalue\n");
	for(i=0;i<HLV->PFR->NumPops;i++) {
		DoFDR_ForAPop(i,HLV);
		//printf("Done with FDR for Pop= %d\n",i);
	}
	fclose(HLV->FDR_Summaries_File);
	
	printf("\n\n");
	SeedToFile("snppit_seeds");
	printf("\n\n");
	
	
	printf("PRINTING FINAL PARENTAGE REPORT\n");
	PrintFinalIndivReportWithFDRs(HLV);
	
	printf("\n\n");
	printf("SNPPIT PROGRAM EXECUTION COMPLETED.\n");
	printf("\nOutput is in the following files:\n");
	printf("\tsnppit_output_ParentageAssignments.txt -- Main output file that gives false discovery rates for all offspring with the most likely parents\n");
	printf("\tsnppit_output_BasicDataSummary.txt -- Basic information about the data that got read in.\n"); 
	printf("\tsnppit_output_ChosenSMAXes.txt -- Information about the smax vectors used in the analysis.\n"); 
	printf("\tsnppit_output_FDR_Summary.txt -- Offspring assigned to parents in each population, ranked by false discovery rate.\n"); 
	printf("\tsnppit_output_PopSizesAnPiVectors.txt -- Information about the sizes of the populations and the expected fraction of different trios thereby implied.\n"); 
	printf("\tsnppit_output_TrioPosteriors.txt -- Posterior probs for all non-excluded parent pairs of every offspring in the data file.\n"); 
	printf("\n\n");
	printf("Questions, etc.? Send them to eric.anderson@noaa.gov\n\n");
	
	
	return(0);		
}


