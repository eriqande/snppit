#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "MCTypesEtc.h"
#include "pbt_C_fb.h"
#include "uthash.h"
#include "pfr_read_genos.h"
#include "pfr_pedigree_spec.h"
#include "snp_sumped.h"


/* this is a function that automates the task of filling up the Collection arrays
   with forwards-backwards information for trios.  It does this with some system calls to
   snpSumPed which puts output into a file that is then read by 
   Get_pbt_C_fb_Opts() and which creates an FB_Vars struct.
*/
FB_Vars *CreateCompleteFB_Vars(
									   int ParentPopIdx,		/* if >=0, compute forback tables assuming parents come from this single population, if -1 do all the pops */
									   int KidPopIdx,			/* if >=0, compute forback tables assuming kids come from this single population, if -1 do all the pops, if -2 then kids' pop = parent's pop always */
									   int PedSpecIdx,			/* if>=0, just compute forback tables for this one relationship. If -1 compute for all relationships.  */
									   int *smax,				/* the smax value that is desired */
									   int smax_d,				/* number of components of this smax vector */
									   pfr_geno_data *PFR,		/* holds all the data about number of loci and allele freqs and such */
									   char *SSP_InputFileName,	/* file into which to spew the snpSumPed input command-file */
									   char *SSP_OutFileName		/* file into which to spew the snpSumPed output */
									   )
{
	int i,j,k;
	int HiParentIdx,HiKidIdx,HiPedIdx;
	int LoParentIdx,LoKidIdx,LoPedIdx;
	int TimesThrough=0;
	FB_Vars *A;
	char *argv[3];
	char CommFile[]="--command-file";
	int argc = 3;
	int *vdex;
	char smax_string[200];
	char *SNPSUMPED;
	/* we use this stuff now that we call snpSumPed internally */
	char qopt[]="-q";
	char *sargv[4];
	int sargc = 4;
	FILE *ssp_outp;
	
	/* open up a file to write the snpsumped output to multiple times */
	if( (ssp_outp = fopen(SSP_OutFileName,"w"))==NULL) {
		fprintf(stderr,"Error.  Unable to open file %s to write to it in CreateCompleteFB_Vars.\n",SSP_OutFileName);
		exit(1);
	}
	
	
	/* this is what we used to initialize SNPSUMPED to: "/Users/eriq/Documents/xp_dev_svn_checkouts/snpSumPed/OSX/bin/Release/snpSumPed" */
	
	if(gHasSnpSumPedPath) {
		SNPSUMPED = gSnpSumPedPath;
	}
	
	/* figure out the limits for the for loops */
	if(ParentPopIdx>=0) {
		LoParentIdx=HiParentIdx=ParentPopIdx;
	}
	else {
		LoParentIdx = 0;
		HiParentIdx = PFR->NumPops-1;
	}
	if(KidPopIdx>=0) {
		LoKidIdx = HiKidIdx = KidPopIdx;
	}
	else {  /* note, we set it here for the -1 case, but will modify it within the for loop below if it is -2 */
		LoKidIdx = 0;
		HiKidIdx = PFR->NumPops-1;
	}
	if(PedSpecIdx>=0) {
		LoPedIdx=HiPedIdx=PedSpecIdx;
	}
	else {
		LoPedIdx = 0;
		HiPedIdx = NUM_SPEC_PEDS-1;
	}
	
	
	/* set the smax string */
	sprintf(smax_string," --s-max ");
	for(i=0;i<smax_d;i++) {
		sprintf(smax_string,"%s %d ",smax_string,smax[i]);
	}
	
	
	printf("Compiling trio type probabilities for %d parental collections\n",PFR->NumPops); 
	fflush(stdout);
	for(i=LoParentIdx;  i<=HiParentIdx;   i++)   {
		
		

		
		if(KidPopIdx==-2) {
			LoKidIdx = i;
			HiKidIdx = i;
		}
		
		for(k=LoKidIdx;    k<=HiKidIdx;  k++) {
			
			for(j=LoPedIdx;     j<=HiPedIdx;     j++) { double **youthp;
				
				youthp=NULL;
				
				if(i != k) {
					if(j != 8) {
						fprintf(stderr,"Error! Whoa there pardner! You are asking for a cross-pop trio with i=%d and k=%d, but you request PedIdx=%d and not C_U_U which is PedIdx 8. Exiting\n",i,k,j);
						exit(1);
					}
					else {
						youthp = PFR->AlleFreqs[k];
					}
				}
				
				Print_snpSumPed_CommandFile(SSP_InputFileName,			/* name of command file for snpSumPed */
											j,								/*PedIdx		*/
											PFR->AlleFreqs[i],				/* p			*/
											youthp,							/* youth p		*/
											PFR->GtypErrRates,				/* mu			*/
											0,								/* MissMeth		*/
											TimesThrough==0,				/* DoPream		*/
											PFR->PopNames[i],				/* MaPopPopName */
											PFR->PopNames[k],				/* KidPopName	*/
											"Boingo",						/* RunName		*/
											PFR->NumLoci,					/* NumLoc		*/
											smax_string						/* SmaxStr		*/
											);
				
				
				/*  THIS WAS OLD STUFF I USED WHEN I CALLED ALL THIS STUFF USING SYSTEM CALLS
				if(TimesThrough==0) {
					sprintf(SysCall," %s -q --command-file %s > %s ",SNPSUMPED, SSP_InputFileName,SSP_OutFileName);
					if(system(SysCall)) {
						fprintf(stderr,"Error! System call: \"%s\" exited with failure status. i=%d  j=%d  k=%d.  Exiting.\n",SysCall,i,j,k);
					}
				}
				else {
					sprintf(SysCall," %s -q --command-file %s >> %s ",SNPSUMPED, SSP_InputFileName,SSP_OutFileName);
					if(system(SysCall)) {
						fprintf(stderr,"Error! System call: \"%s\" exited with failure status. i=%d  j=%d  k=%d.  Exiting.\n",SysCall,i,j,k);
					}
				}
				*/
				sargv[0]=NULL;
				sargv[1]=qopt;
				sargv[2]=CommFile;
				sargv[3]=SSP_InputFileName;
				
				snpSumPed(sargc, sargv, ssp_outp);
				
				
				TimesThrough++;
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
	
	fclose(ssp_outp);
	
	/* now that we have created the output file that we want, we read make the FB_Vars with it. */
	argv[0]=NULL;
	argv[1]=CommFile;
	argv[2]=SSP_OutFileName;
	
	/* read all the data in */
	A = Get_pbt_C_fb_Opts(argc,argv);
	
	/*	NOTE!  I SHOULD DO SOME MEMORY FREEING AT SOME POINT, AND CAN PROBABLY DO THAT BY WRITING SOME MORE CODE IN ECA_OPT, BECUASE
	 AS IT IS NOW ARGV ENDS UP BEING HUGE AFTER READING ALL THIS */
	
	
	/* then do the forward algorithm */
	/* do the forward step and prep for the backward steps on all collections */
	vdex = (int *)calloc( ReturnVDex(A->RP->d, A->RP->smax,A->RP->smax)+1,sizeof(int));
	printf("\nPerforming Forward Step on %d Collections of Trio Probabilities\n",A->NumColls);
	for(j=0;j<A->NumColls;j++)  {
		for(i=0;i<=ReturnVDex(A->RP->d, A->RP->smax,A->RP->smax);i++) {
			vdex[i]=-1;
		}
		A->Colls[j]->FBS = ForwardStep(A->RP, A->Colls[j]->Xprobs, A->Colls[j]->Xfakeprobs, vdex, 0);	
		PrepForBackwardStep(A->Colls[j]->FBS, A->RP, A->Colls[j]->Xprobs);
		
		if(j%100==0) {
			if(j==0) printf("0");
			else printf("\n%d",j);
		}
		else {
			printf(".");
		}
		fflush(stdout);
		
	}
	
	printf("\n");
	A->vdex = vdex;	
	
	return(A);
}






/* this function will write information to a file which can then be run by snpSumPed */
void Print_snpSumPed_CommandFile(char *FileName,
								int PedIdx,
								double **p,
								double **youth_p,  /* if different than parental p.  Otherwise pass in as NULL */
								double *mu,
								int *MissMeth,
								int DoPream,
								char *MaPopPopName,
								char *KidPopName,
								char *RunName,
								int NumLoc,
								const char *SmaxStr
								)
{
	
	/* the first specifier relates to the father, the second to the mother in these names
	 
	 Note that 2= putative father  3=putative mother  1=youth  */
	
	
	int i;
	FILE *in;
	char *PedSpecs[] = { \
		" & PedSpec_C_Se_Se: (0)  C-type which is a true parental trio & \
		--ped 1 2 3   \
		--ped 2 0 0  \
		--ped 3 0 0  ", \
		" & PedSpec_C_Se_Si: (1)  C-type with father self and the mother sib & \
		--ped 1 4 5   \
		--ped 2 0 0  \
		--ped 3 0 0  \
		--ped 4 0 0  \
		--ped 5 0 0  \
		--kappa 2 4 0 0 1  \
		--kappa 3 5 .25 .5 .25  ", \
		" & PedSpec_C_Si_Se: (2)  C-type with mother self and the father sib & \
		--ped 1 4 5   \
		--ped 2 0 0  \
		--ped 3 0 0  \
		--ped 4 0 0  \
		--ped 5 0 0  \
		--kappa 3 5 0 0 1  \
		--kappa 2 4 .25 .5 .25  ", \
		" & PedSpec_C_Se_U: (3)  C-type with father self and the mother unrelated & \
		--ped 1 4 5   \
		--ped 2 0 0  \
		--ped 3 0 0  \
		--ped 4 0 0  \
		--ped 5 0 0  \
		--kappa 2 4 0 0 1  ", \
		" & PedSpec_C_U_Se: (4)  C-type with mother self and the father unrelated & \
		--ped 1 4 5   \
		--ped 2 0 0  \
		--ped 3 0 0  \
		--ped 4 0 0  \
		--ped 5 0 0  \
		--kappa 3 5 0 0 1  ", \
		" & PedSpec_C_Si_Si: (5)  C-type with one putative parent sib and the other sib & \
		--ped 1 4 5   \
		--ped 2 0 0  \
		--ped 3 0 0  \
		--ped 4 0 0  \
		--ped 5 0 0  \
		--kappa 2 4 .25 .5 .25  \
		--kappa 3 5 .25 .5 .25  ", \
		" & PedSpec_C_Si_U: (6)  C-type with one putative father a sib and the other unrelated & \
		--ped 1 4 5   \
		--ped 2 0 0  \
		--ped 3 0 0  \
		--ped 4 0 0  \
		--ped 5 0 0  \
		--kappa 2 4 .25 .5 .25  ", \
		" & PedSpec_C_U_Si: (7)  C-type with putative mother sib and the other unrelated & \
		--ped 1 4 5   \
		--ped 2 0 0  \
		--ped 3 0 0  \
		--ped 4 0 0  \
		--ped 5 0 0  \
		--kappa 3 5 .25 .5 .25  ", \
		" & PedSpec_C_U_U: (8)  C-type with unrelated putative parents & \
		--ped 1 0 0   \
		--ped 2 0 0  \
		--ped 3 0 0  ", /* note that you have to specify this with a --ped 1 0 0 if it is to with the CrossPop trios */ \
		" & PedSpec_Se_B: (9)  putative mother is a full sib of the offspring and putative father is self & \
		--ped 1 4 5  \
		--ped 3 4 5  \
		--ped 2 0 0  \
		--ped 4 0 0  \
		--ped 5 0 0  \
		--kappa 2 4 0 0 1  ", \
		" & PedSpec_B_Se: (10)  putative father is a full sib of the offspring and putative mother is self & \
		--ped 1 4 5  \
		--ped 2 4 5  \
		--ped 3 0 0  \
		--ped 4 0 0  \
		--ped 5 0 0  \
		--kappa 3 5 0 0 1  ", \
		" & PedSpec_H_Se: (11)  father is a half sib and the mother is self & \
		--ped 1 4 5  \
		--ped 2 4 6  \
		--ped 3 0 0  \
		--ped 4 0 0  \
		--ped 5 0 0  \
		--ped 6 0 0  \
		--kappa 3 5 0 0 1  ", \
		" & PedSpec_Se_H: (12)  mother is a half sib and the father is self & \
		--ped 1 4 5  \
		--ped 3 4 6  \
		--ped 2 0 0  \
		--ped 4 0 0  \
		--ped 5 0 0  \
		--ped 6 0 0  \
		--kappa 2 4 0 0 1  ", \
		" & PedSpec_B_Si: (13)  putative father is a full sib of the offspring and the putative mother is a sib of the true mother (or of the true father) & \
		--ped 1 4 5  \
		--ped 2 4 5  \
		--ped 3 0 0  \
		--ped 4 0 0  \
		--ped 5 0 0  \
		--kappa 3 5 .25 .5 .25  ", \
		" & PedSpec_Si_B: (14)  putative mather is a full sib of the offspring and the putative father is a sib of the true father (or true mother) & \
		--ped 1 4 5  \
		--ped 3 4 5  \
		--ped 2 0 0  \
		--ped 4 0 0  \
		--ped 5 0 0  \
		--kappa 2 4 .25 .5 .25  ", \
		" & PedSpec_B_U: (15)  putative father is a full sib of the offspring and the mother is unrelated & \
		--ped 1 4 5  \
		--ped 2 4 5  \
		--ped 3 0 0  \
		--ped 4 0 0  \
		--ped 5 0 0  ", \
		" & PedSpec_U_B: (16)  putative mother is a full sib of the offspring and the father is unrelated & \
		--ped 1 4 5  \
		--ped 3 4 5  \
		--ped 2 0 0  \
		--ped 4 0 0  \
	    --ped 5 0 0  ", \
		" & PedSpec_B_B: (17)  putative mother is a full sib of the offspring and the putative father is a full sib of the offspring too & \
		--ped 1 4 5  \
		--ped 2 4 5  \
		--ped 3 4 5  \
		--ped 4 0 0  \
		--ped 5 0 0  "};
	
	
	/* NOTE!! If you change these, you also have to change the function that returns the string name of each pedigree.  Really stupid I know!  I should
	 just incorporate it all into a single place somewhere.  Oh Well. */
	char *PedSpecNames[] = { \
		"C_Se_Se",     /*  0  */    \
		"C_Se_Si",     /*  1  */    \
		"C_Si_Se",     /*  2  */    \
		"C_Se_U",     /*  3  */    \
		"C_U_Se",     /*  4  */    \
		"C_Si_Si",     /*  5  */    \
		"C_Si_U",     /*  6  */    \
		"C_U_Si",     /*  7  */    \
		"C_U_U",     /*  8  */    \
		"Se_B",     /*  9  */    \
		"B_Se",     /*  10  */    \
		"H_Se",     /*  11  */    \
		"Se_H",     /*  12  */    \
		"B_Si",     /*  13  */    \
		"Si_B",     /*  14  */    \
		"B_U",     /*  15  */    \
		"U_B",     /*  16  */    \
		"B_B"      /*  17  */    \
	};
	
	if((in=fopen(FileName,"w"))==NULL) {
		fprintf(stderr,"Error! Unable to open file %s to write to it in Print_snpSumPed_CommandFile. Exiting...\n",FileName);
		exit(1);
	}
	
	fprintf(in,"--miss-meth %d\n",MissMeth); 
	fprintf(in,"-f 1 2 3\n");
	fprintf(in,"--xml-output %sx%s__%s\n",MaPopPopName,KidPopName,PedSpecNames[PedIdx]);
	if(DoPream) {
		fprintf(in,"--xml-pream %s\n%s\n",RunName,SmaxStr);
		
	}
	fprintf(in,"-p ");
	for(i=0;i<NumLoc;i++) {
		fprintf(in," %f",p[i][0]);
	}
	fprintf(in,"\n");
	if(youth_p) {
		fprintf(in," --youth-p ");
		for(i=0;i<NumLoc;i++) {
			fprintf(in," %f",youth_p[i][0]);
		}
		fprintf(in,"\n");
	}
	fprintf(in,"-e ");
	for(i=0;i<NumLoc;i++) {
		fprintf(in," %f",mu[i]);
	}
	fprintf(in,"\n");
	fprintf(in,"%s\n",PedSpecs[PedIdx]);
	
	
	fclose(in);
	
}
								

const char *SpecPedIdx2Str(int P)
{
	switch(P) {
		case(0):
			return("C_Se_Se");
			break;
		case(1):
			return("C_Se_Si");
			break;
		case(2):
			return("C_Si_Se");
			break;
		case(3):
			return("C_Se_U");
			break;
		case(4):
			return("C_U_Se");
			break;
		case(5):
			return("C_Si_Si");
			break;
		case(6):
			return("C_Si_U");
			break;
		case(7):
			return("C_U_Si");
			break;
		case(8):
			return("C_U_U");
			break;
		case(9):
			return("Se_F");
			break;
		case(10):
			return("F_Se");
			break;
		case(11):
			return("H_Se");
			break;
		case(12):
			return("Se_H");
			break;
		case(13):
			return("F_Si");
			break;
		case(14):
			return("Si_F");
			break;
		case(15):
			return("F_U");
			break;
		case(16):
			return("U_F");
			break;
		case(17):
			return("F_F");
			break;
		default:
			fprintf(stderr,"Unrecognized Spec Ped Idx %d in SpecPedIdx2Str().  Exiting...\n");
			exit(1);
			break;
	}
}

