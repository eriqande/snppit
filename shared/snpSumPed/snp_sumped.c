
/*



	snp_sumped SOURCE CODE WRITTEN BY ERIC C ANDERSON
	
	eric.anderson@noaa.gov
	
	Public domain source code developed by Federal employee

	Sofware for calculating the power to assign individuals to parents
	using binary markers, like SNPs.

*/

#define UN_EXTERN

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
#include "ECA_Opt3.h"
#include "ECA_print.h"
#include "snp_sumped.h"

/* global variables */
double gSumOverFixed, gSumOverAll;
double gTrioGProbs[3][3][3];  /* for the prob of the genotypes of the fixed indiv's
						      if there are only three such fixed ones */
double gTrioObsProbs[3][3][3];  /* prob of the trios observed probs given genotyping error */
double gTrioMissObsProbs[4][4][4]; /* prob of the trios when you include missing as a data state */
double gTrioMissFakeProbs[4][4][4];  /* "prob" of trios when you just ignore the missing data at whichever
										individuals have missing data. In some sense this is kind of like marginalizing
										out the missing data at an individual and assuming that they are all independent...
										or is it? This might be totally hare-brained. */
unsigned long int gVerb;  /* controls the verbosity via bitmasking. */  
int gOutputEcaOptToFile;
FILE *gOutfile;
double gYouthP;  /* this is part of a little hack to let me compute the prob of unrelated trios in which the two parents are from the 
					population with the frequency p as given on the command line and the youth (indiv with the 1 index) is from a 
					different population with frequencies given by the youth-p option */
int gHasYouthP=0;  /* for denoting whether or not the youth-p option was given. */








/*
	This simply returns the index of the X-state of the trio.
	X-states are as follows:
		
	Idx :  kid-pa  kid-ma  whole-trio
    0 :   0  0  0  (i.e. the kid is compatible with pa and ma)
	1 :   0  0  1
    2 :   0  1  1
    3 :   1  0  1
	4 :   1  1  1
*/
int XVectorOfGenos(int kid, int pa, int ma)
{
	int pa_inc, ma_inc, tot_inc;
	
	pa_inc = 0;
	pa_inc = ( (kid==2 && pa==0) || (kid==0 && pa==2) );
	
	ma_inc = 0;
	ma_inc = ( (kid==2 && ma==0) || (kid==0 && ma==2) );
	
	tot_inc = 0;
	/* tot_inc is 1 only for that special case of parents compatible with the kid individually, but the whole
	   trio is not compatible. */
	tot_inc = (kid==1 && ( (pa==0 && ma==0) || (pa==2) && (ma==2)));
	
	if(pa_inc + ma_inc + tot_inc == 0) return(0);
	if(tot_inc) return(1);
	if(ma_inc && !pa_inc) return(2);
	if(pa_inc && !ma_inc) return(3);
	if(pa_inc && ma_inc) return(4);
	
	fprintf(stderr,"Got to where we shouldn't in XVectorOfGenos. Exiting.\n");
	exit(1);
	
	return(-9999);
	
}


/*
 This simply returns the index of the X-state of the trio.
 
 The inputs are the genotypes of the individuals where 0-2 count the
 number of 1 alleles and 3 is missing data.
 
 X-states are as follows:
 
 Idx :  kid-pa  kid-ma  whole-trio  at-least-one-miss

 0:   0 0 0 0
 1:   0 0 1 0
 2:   0 1 1 0
 3:   1 0 1 0
 4:   1 1 1 0
 5:   0 0 0 1
 6:   0 1 1 1
 7:   1 0 1 1

 This is pretty straightforward.  If there are  more then two indivs with missing data in the trio
 then there can be no incompats.  With one, there can never be a 001 incompat.
 
 
 */
int XVectorOfGenosWithMiss1(int kid, int pa, int ma)
{
	int pa_inc, ma_inc, tot_inc;
	int kidmiss, mamiss, pamiss;
	
	pa_inc = 0;
	pa_inc = ( (kid==2 && pa==0) || (kid==0 && pa==2) );
	
	ma_inc = 0;
	ma_inc = ( (kid==2 && ma==0) || (kid==0 && ma==2) );
	
	tot_inc = 0;
	/* tot_inc is 1 only for that special case of parents compatible with the kid individually, but the whole
	 trio is not compatible. */
	tot_inc = (kid==1 && ( (pa==0 && ma==0) || (pa==2) && (ma==2)));
	
	kidmiss = (kid==3);
	mamiss = (ma==3);
	pamiss = (pa==3);
	
	if(kidmiss+mamiss+pamiss==0) {
		if(pa_inc + ma_inc + tot_inc == 0) return(0);
		if(tot_inc) return(1);
		if(ma_inc && !pa_inc) return(2);
		if(pa_inc && !ma_inc) return(3);
		if(pa_inc && ma_inc) return(4);
	}
	else if(kidmiss+mamiss+pamiss>=2) {
		return(5);
	}
	else {  /* in this case there is only one individual missing data at the locus */
		if(kidmiss) {
			return(5);
		}
		else if(pamiss && ma_inc) {
			return(6);
		}
		else if(mamiss && pa_inc) {
			return(7);
		}
		else {  /* down here, we are missing an individual but have no incompatibilities. */
			return(5);
		}
	}
	
	fprintf(stderr,"Got to where we shouldn't in XVectorOfGenosWithMiss1.  kid=%d  pa=%d  ma=%d    pa_inc=%d  ma_inc=%d  tot_inc=%d   kidmiss=%d  pamiss=%d   mamiss=%d  Exiting.\n",
			kid,ma,pa,pa_inc,ma_inc,tot_inc,kidmiss,pamiss,mamiss);
	exit(1);
	
	return(-9999);
}	


/*
 This simply returns the index of the X-state of the trio.
 
 The inputs are the genotypes of the individuals where 0-2 count the
 number of 1 alleles and 3 is missing data.
 
 X-states are as follows:
 
 Idx :  kid-pa  kid-ma  whole-trio  kidd-miss   ma-miss    pa-miss
 
 0:  000000
 1:  001000
 2:  101000
 3:  011000
 4:  111000
 5:  000100    
 6:  000010    
 7:  101010    
 8:  000110    
 9:  000001    
 10:  011001   
 11:  000101
 12:  000011
 13:  000111
 
 I'm sure I could have coded this up more elegantly or efficiently, but I just wanted to bang it out without much
 thought, even if it is totally ugly and stupid.
 
 This is the "3" version because it is used when counting number of missing loci at all 3 individuals in the trio.
 
 */
int XVectorOfGenosWithMiss3(int kid, int pa, int ma)
{
	int pa_inc, ma_inc, tot_inc;
	int kidmiss, mamiss, pamiss;
	
	pa_inc = 0;
	pa_inc = ( (kid==2 && pa==0) || (kid==0 && pa==2) );
	
	ma_inc = 0;
	ma_inc = ( (kid==2 && ma==0) || (kid==0 && ma==2) );
	
	tot_inc = 0;
	/* tot_inc is 1 only for that special case of parents compatible with the kid individually, but the whole
	 trio is not compatible. */
	tot_inc = (kid==1 && ( (pa==0 && ma==0) || (pa==2) && (ma==2)));
	
	kidmiss = (kid==3);
	mamiss = (ma==3);
	pamiss = (pa==3);
	
	if(kidmiss+mamiss+pamiss==0) {
		if(pa_inc + ma_inc + tot_inc == 0) return(0);
		if(tot_inc) return(1);
		if(ma_inc && !pa_inc) return(2);
		if(pa_inc && !ma_inc) return(3);
		if(pa_inc && ma_inc) return(4);
	}
	if(kidmiss && !(mamiss || pamiss) ) {
		if(pa_inc + ma_inc + tot_inc == 0) return(5);
		else {
			fprintf(stderr,"Error in XVectorOfGenosWithMiss. Incompatibiles with kidmiss = 1\n");
			exit(1);
		}
	}
	if(mamiss && !(kidmiss || pamiss) ) {
		if(pa_inc + ma_inc + tot_inc == 0) return(6);
		else if(pa_inc) return(7);
		else {
			fprintf(stderr,"Error in XVectorOfGenosWithMiss when mamiss && !(kidmiss || pamiss) \n");
			exit(1);
		}
	}
	if(mamiss && kidmiss && !pamiss) {
		if(pa_inc + ma_inc + tot_inc == 0) return(8);
		else {
			fprintf(stderr,"Error in XVectorOfGenosWithMiss when mamiss && kidmiss \n");
			exit(1);
		}
	}
	if(pamiss && !(kidmiss || mamiss) ) {
		if(pa_inc + ma_inc + tot_inc == 0) return(9);
		else if(ma_inc) return(10);
		else {
			fprintf(stderr,"Error in XVectorOfGenosWithMiss when pamiss && !(kidmiss || mamiss) \n");
			exit(1);
		}
	}
	if(kidmiss && pamiss && !mamiss) {
		if(pa_inc + ma_inc + tot_inc == 0) return(11);
		else {
			fprintf(stderr,"Error in XVectorOfGenosWithMiss when kidmiss && pamiss && !mamiss \n");
			exit(1);
		}
	}
	if(!kidmiss && pamiss && mamiss) {
		if(pa_inc + ma_inc + tot_inc == 0) return(12);
		else {
			fprintf(stderr,"Error in XVectorOfGenosWithMiss when kidmiss && !pamiss && mamiss\n");
			exit(1);
		}
	}
	if(kidmiss && pamiss && mamiss) {
		if(pa_inc + ma_inc + tot_inc == 0) return(13);
		else {
			fprintf(stderr,"Error in XVectorOfGenosWithMiss when kidmiss && pamiss && mamiss\n");
			exit(1);
		}
	}
	
	
	fprintf(stderr,"Got to where we shouldn't in XVectorOfGenos.  kid=%d  pa=%d  ma=%d    pa_inc=%d  ma_inc=%d  tot_inc=%d   kidmiss=%d  pamiss=%d   mamiss=%d  Exiting.\n",
				kid,ma,pa,pa_inc,ma_inc,tot_inc,kidmiss,pamiss,mamiss);
	exit(1);
	
	return(-9999);
	
}








const char *PrintXStateFromIdx(int X)
{
	switch(X) {
		case(0):
			return("0 0 0");
			break;
		case(1):
			return("0 0 1");
			break;
		case(2):
			return("0 1 1");
			break;
		case(3):
			return("1 0 1");
			break;
		case(4):
			return("1 1 1");
			break;
		default:
			fprintf(stderr,"Unrecognized X in PrintXStateFromIdx.  Exiting.\n");
			exit(1);
	}
	return("  DISASTER  ");
}


/*
returns a string with the v-vector for the different X states when the fourth component
 of the v-vector indicates whether or not there is at least one individual in the trio
 who is missing data at this vector.
*/
const char *PrintXStateFromIdxWithMiss1(int X)
{
	switch(X) {
		case(0):
			return("0 0 0 0");
			break;
		case(1):
			return("0 0 1 0");
			break;
		case(2):
			return("0 1 1 0");
			break;
		case(3):
			return("1 0 1 0");
			break;
		case(4):
			return("1 1 1 0");
			break;
		case(5):
			return("0 0 0 1");
			break;
		case(6):
			return("0 1 1 1");
			break;
		case(7):
			return("1 0 1 1");
			break;
		default:
			fprintf(stderr,"Unrecognized X in PrintXStateFromIdxWithMiss1.  Exiting.\n");
			exit(1);
	}
	return("  DISASTER  ");
}



/* this prints out a string with the v-vector for the case where you also
 have missing data indicators.  So, here the first three bits are the
 incompatibity indicators for ma pa and both and the next three are indicators
 for missing data at kid pa ma.  I used the script:
 /Users/eriq/Documents/svn_code_checkouts/snpSumPed/script/EnumerateV-VecsWithMiss.awk
 to help enumerate and print these states.  
 
 This is the "3" version because it keeps track of the number of missing loci in 
 all three members of the trio.
 
 */
const char *PrintXStateFromIdxWithMiss3(int X)
{
	switch(X) {
		case(0):
			return("0 0 0 0 0 0");
			break;
		case(1):
			return("0 0 1 0 0 0");
			break;
		case(2):
			return("1 0 1 0 0 0");
			break;
		case(3):
			return("0 1 1 0 0 0");
			break;
		case(4):
			return("1 1 1 0 0 0");
			break;
		case(5):
			return("0 0 0 1 0 0");
			break;
		case(6):
			return("0 0 0 0 1 0");
			break;
		case(7):
			return("1 0 1 0 1 0");
			break;
		case(8):
			return("0 0 0 1 1 0");
			break;
		case(9):
			return("0 0 0 0 0 1");
			break;
		case(10):
			return("0 1 1 0 0 1");
			break;
		case(11):
			return("0 0 0 1 0 1");
			break;
		case(12):
			return("0 0 0 0 1 1");
			break;
		case(13):
			return("0 0 0 1 1 1");
			break;			
		default:
			fprintf(stderr,"Unrecognized X in PrintXStateFromIdxWithMiss3.  Exiting.\n");
			exit(1);
	}
	return("  DISASTER  ");
}



int V_vec_comp(int X, int i) 
{
	switch(X) {
		case(0):
			return(0);
		case(1):
			switch(i) {
				case(0):
				case(1):
					return(0);
				case(2):
					return(1);
				default:
					fprintf(stderr,"Error X=%d and i=%d in V_vec_comp. Exiting\n",X,i); 
					exit(1);
			}
		case(2):
			switch(i) {
				case(0):
					return(0);
				case(1):	
				case(2):
					return(1);
				default:
					fprintf(stderr,"Error X=%d and i=%d in V_vec_comp. Exiting\n",X,i); 
					exit(1);
			}
		case(3):
			switch(i) {
				case(1):
					return(0);
				case(0):	
				case(2):
					return(1);
				default:
					fprintf(stderr,"Error X=%d and i=%d in V_vec_comp. Exiting\n",X,i); 
					exit(1);
			}
		case(4):
			return(1);
			break;
		default:
			fprintf(stderr,"Unrecognized X in PrintXStateFromIdx.  Exiting.\n");
			exit(1);
	}
	return(-9999999);
}


void PrintXML_Preamble(const char *Ind, sum_ped_pars *pars)
{
	int i,j,k,X,numA;
	char Indent[5000];
	
	sprintf(Indent,"XML_OUTPUT:%s",Ind);
	
	eprintf(gVerb & XML_VERB, "%s<run id=\"%s\" xsi:noNamespaceSchemaLocation=\"FourLocusSimple.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" >\n",Indent,pars->RunIdName);
	eprintf(gVerb & XML_VERB, "%s  <N>5</N>\n",Indent);
	eprintf(gVerb & XML_VERB, "%s  <d>3</d>\n",Indent);
	eprintf(gVerb & XML_VERB, "%s  <T>%d</T>\n",Indent,pars->NumPs);
	for(i=0;i<3;i++) {
		eprintf(gVerb & XML_VERB, "%s  <smax_item>%d</smax_item>\n",Indent,pars->Smax[i]);
	}
	
	
	for(X=0;X<5;X++) {
		eprintf(gVerb & XML_VERB, "%s  <X-desc id=\"%d\">\n",Indent,X);
		for(i=0;i<3;i++) {
			eprintf(gVerb & XML_VERB, "%s    <v-vector-item>%d</v-vector-item>\n",Indent, V_vec_comp(X,i));
		}
		numA=0;
		for(i=0;i<3;i++)  {
			for(j=0;j<3;j++)  {
				for(k=0;k<3;k++)  {
					if(XVectorOfGenos(i,j,k)==X) {
						eprintf(gVerb & XML_VERB, "%s    <A-desc pid=\"%d\" id=\"%d\">\n",Indent,X,numA);
						eprintf(gVerb & XML_VERB, "%s      <genotype> %d %d %d </genotype>\n",Indent,i,j,k);
						eprintf(gVerb & XML_VERB, "%s    </A-desc>\n",Indent,numA);
						numA++;
					}
				}
			}
		}
		eprintf(gVerb & XML_VERB, "%s  </X-desc>\n",Indent);
	}
}


void PrintECA_OPT_Preamble(const char *Ind, sum_ped_pars *pars)
{
	int i,j,k,X;
	char Indent[5000];
	int NaInX[100];
	int Nx=5,d=3,Na=27,As=3;
	
	if(pars->MissMeth==1) {
		Nx=8;
		d=4;
		Na=64;
		As=4;  /* this is the number of allele states.  We count missing as an allele state for the purposes here */
	}
	
	if(pars->MissMeth==3) {
		Nx=14;
		d=6;
		Na=64;
		As=4;  /* this is the number of allele states.  We count missing as an allele state for the purposes here */
	}
	
	
	/* if we want this to be bare output, we don't put the ECA_OPT_OUTPUT tag on it */
	if(gVerb & ECA_OPT_BARE_VERB) {
			sprintf(Indent,"%s",Ind);
	}
	else {
			sprintf(Indent,"ECA_OPT_OUTPUT:%s",Ind);
	}
	
	if(gOutputEcaOptToFile==0) {
		eprintf(gVerb & ECA_OPT_VERB, "%s--Nx %d\n",Indent,Nx);
		eprintf(gVerb & ECA_OPT_VERB, "%s-d %d\n",Indent,d);
		eprintf(gVerb & ECA_OPT_VERB, "%s-L %d\n",Indent,pars->NumPs);
		eprintf(gVerb & ECA_OPT_VERB, "%s--smax ",Indent);
	}
	else {
		efprintf(gVerb & ECA_OPT_VERB, gOutfile, "%s--Nx %d\n",Indent,Nx);
		efprintf(gVerb & ECA_OPT_VERB, gOutfile, "%s-d %d\n",Indent,d);
		efprintf(gVerb & ECA_OPT_VERB, gOutfile, "%s-L %d\n",Indent,pars->NumPs);
		efprintf(gVerb & ECA_OPT_VERB, gOutfile, "%s--smax ",Indent);
	}
	for(i=0;i<d;i++) {
		if(gOutputEcaOptToFile==0) {
			eprintf(gVerb & ECA_OPT_VERB, "%d ",pars->Smax[i]);
		}
		else {
			efprintf(gVerb & ECA_OPT_VERB, gOutfile, "%d ",pars->Smax[i]);
		}
	}
	
	if(gOutputEcaOptToFile==0) {
		eprintf(gVerb & ECA_OPT_VERB, "\n");
		eprintf(gVerb & ECA_OPT_VERB, "%s--v-vecs\n",Indent);
	}
	else {
		efprintf(gVerb & ECA_OPT_VERB, gOutfile, "\n");
		efprintf(gVerb & ECA_OPT_VERB, gOutfile, "%s--v-vecs\n",Indent);
	}
	
	
	
	for(X=0;X<Nx;X++) {
		if(pars->MissMeth==1) {
			if(gOutputEcaOptToFile==0) eprintf(gVerb & ECA_OPT_VERB, "%s %s\n",Indent,PrintXStateFromIdxWithMiss1(X));
			else efprintf(gVerb & ECA_OPT_VERB, gOutfile, "%s %s\n",Indent,PrintXStateFromIdxWithMiss1(X));
		}
		else if(pars->MissMeth==3) {
			if(gOutputEcaOptToFile==0) eprintf(gVerb & ECA_OPT_VERB, "%s %s\n",Indent,PrintXStateFromIdxWithMiss3(X));
			else efprintf(gVerb & ECA_OPT_VERB, gOutfile, "%s %s\n",Indent,PrintXStateFromIdxWithMiss3(X));
		}
		else if(pars->MissMeth==0) {
			if(gOutputEcaOptToFile==0) eprintf(gVerb & ECA_OPT_VERB, "%s %s\n",Indent,PrintXStateFromIdx(X));
			else efprintf(gVerb & ECA_OPT_VERB, gOutfile, "%s %s\n",Indent,PrintXStateFromIdx(X));
		}
		else {
			fprintf(stderr,"Error! Unrecognized MissMeth %d in PrintECA_OPT_Preamble. Exiting\n",pars->MissMeth);
			exit(1);
		}
		NaInX[X]=0;  /* do this initialization to collect a sum later */
	}
	
	if(gOutputEcaOptToFile==0) eprintf(gVerb & ECA_OPT_VERB, "%s--Na %d\n",Indent,Na);
	else efprintf(gVerb & ECA_OPT_VERB, gOutfile, "%s--Na %d\n",Indent,Na);
	
	/* now, print the AinX's */
	if(gOutputEcaOptToFile==0) eprintf(gVerb & ECA_OPT_VERB, "%s--AinX",Indent);
	else  efprintf(gVerb & ECA_OPT_VERB, gOutfile, "%s--AinX",Indent);
	
	for(i=0;i<As;i++)  {
		for(j=0;j<As;j++)  {
			for(k=0;k<As;k++)  {
				if(pars->MissMeth==3) {
					X=XVectorOfGenosWithMiss3(i,j,k);
				}
				else if(pars->MissMeth==1) {
					X=XVectorOfGenosWithMiss1(i,j,k);
				}
				else {
					X=XVectorOfGenos(i,j,k);
				}
				if(gOutputEcaOptToFile==0) eprintf(gVerb & ECA_OPT_VERB," %d",X);
				else efprintf(gVerb & ECA_OPT_VERB, gOutfile," %d",X);
				
				NaInX[X]++;
			}
		}
	}
	
	if(gOutputEcaOptToFile==0) {
		eprintf(gVerb & ECA_OPT_VERB, "\n");
		eprintf(gVerb & ECA_OPT_VERB, "%s--NaInX",Indent);
	}
	else {
		efprintf(gVerb & ECA_OPT_VERB, gOutfile, "\n");
		efprintf(gVerb & ECA_OPT_VERB, gOutfile, "%s--NaInX",Indent);
	}
	
	
	
	for(X=0;X<Nx;X++) {
		if(gOutputEcaOptToFile==0) eprintf(gVerb & ECA_OPT_VERB, " %d",NaInX[X]);
		else efprintf(gVerb & ECA_OPT_VERB, gOutfile, " %d",NaInX[X]);
	}
	
	if(gOutputEcaOptToFile==0) {
		eprintf(gVerb & ECA_OPT_VERB, "\n");
		eprintf(gVerb & ECA_OPT_VERB, "%s--NY 3\n",Indent);
		eprintf(gVerb & ECA_OPT_VERB, "%s--Ystates\n",Indent);
	}
	else {
		efprintf(gVerb & ECA_OPT_VERB, gOutfile, "\n");
		efprintf(gVerb & ECA_OPT_VERB, gOutfile, "%s--NY 3\n",Indent);
		efprintf(gVerb & ECA_OPT_VERB, gOutfile, "%s--Ystates\n",Indent);
	}
	
	for(i=0;i<As;i++)  {
		for(j=0;j<As;j++)  {
			for(k=0;k<As;k++)  {
				if(gOutputEcaOptToFile==0) eprintf(gVerb & ECA_OPT_VERB, "%s        %d %d %d\n",Indent,i,j,k);
				else efprintf(gVerb & ECA_OPT_VERB, gOutfile, "%s        %d %d %d\n",Indent,i,j,k);
			}
		}
	}
	
	if(gOutputEcaOptToFile==0) {
		eprintf(gVerb & ECA_OPT_VERB, "\n");
		eprintf(gVerb & ECA_OPT_VERB, "%s--num-miss-modes   8\n",Indent);
		eprintf(gVerb & ECA_OPT_VERB, "%s--miss-data-modes   \n",Indent);
		eprintf(gVerb & ECA_OPT_VERB, "%s    0     0    0\n\
%s    0     0    1\n\
%s    0     1    0\n\
%s    1     0    0\n\
%s    0     1    1\n\
%s    1     1    0\n\
%s    1     0    1\n\
%s    1     1    1\n\n",Indent,Indent,Indent,Indent,Indent,Indent,Indent,Indent);
	eprintf(gVerb & ECA_OPT_VERB, "%s--miss-mode-vvec-cens\n",Indent);
	eprintf(gVerb & ECA_OPT_VERB, "%s    0   0   0   0   0\n\
%s    0   1   1   0   1\n\
%s    0   1   0   1   1\n\
%s    0   1   1   1   1\n\
%s    0   1   1   1   1\n\
%s    0   1   1   1   1\n\
%s    0   1   1   1   1\n\
%s    0   1   1   1   1\n\n",Indent,Indent,Indent,Indent,Indent,Indent,Indent,Indent);
 	}
	else {
	efprintf(gVerb & ECA_OPT_VERB, gOutfile, "\n");
	efprintf(gVerb & ECA_OPT_VERB, gOutfile, "%s--num-miss-modes   8\n",Indent);
	efprintf(gVerb & ECA_OPT_VERB, gOutfile, "%s--miss-data-modes   \n",Indent);
	efprintf(gVerb & ECA_OPT_VERB, gOutfile, "%s    0     0    0\n\
%s    0     0    1\n\
%s    0     1    0\n\
%s    1     0    0\n\
%s    0     1    1\n\
%s    1     1    0\n\
%s    1     0    1\n\
%s    1     1    1\n\n",Indent,Indent,Indent,Indent,Indent,Indent,Indent,Indent);
	efprintf(gVerb & ECA_OPT_VERB, gOutfile, "%s--miss-mode-vvec-cens\n",Indent);
	efprintf(gVerb & ECA_OPT_VERB, gOutfile, "%s    0   0   0   0   0\n\
%s    0   1   1   0   1\n\
%s    0   1   0   1   1\n\
%s    0   1   1   1   1\n\
%s    0   1   1   1   1\n\
%s    0   1   1   1   1\n\
%s    0   1   1   1   1\n\
%s    0   1   1   1   1\n\n",Indent,Indent,Indent,Indent,Indent,Indent,Indent,Indent);
	}
}







void PrintXML_LocusProbs(const char *Ind, double Pr[3][3][3], char *LocName)
{
	int i,j,k,X,numA;
	char Indent[5000];
	double PX=0.0;
	
	sprintf(Indent,"XML_OUTPUT:%s",Ind);
	
	eprintf(gVerb & XML_VERB, "%s<locus name=\"%s\">\n",Indent,LocName);
	for(X=0;X<5;X++) {
		eprintf(gVerb & XML_VERB, "%s  <X id=\"%d\">\n",Indent,X);
		
		/* first cycle over all A states to compute the prob of this X state */
		PX=0.0;
		for(i=0;i<3;i++)  {
			for(j=0;j<3;j++)  {
				for(k=0;k<3;k++)  {
					if(XVectorOfGenos(i,j,k)==X) {
						PX += Pr[i][j][k];
					}
				}
			}
		}
		eprintf(gVerb & XML_VERB, "%s    <PX>%.16f</PX>\n",Indent,PX);
		
		/* now, cycle again to print out the probabilities of the A-states */
		numA=0;
		for(i=0;i<3;i++)  {
			for(j=0;j<3;j++)  {
				for(k=0;k<3;k++)  {
					if(XVectorOfGenos(i,j,k)==X) {
						eprintf(gVerb & XML_VERB, "%s    <A pid=\"%d\" id=\"%d\">\n",Indent,X,numA);
						eprintf(gVerb & XML_VERB, "%s      <PA>%.16f</PA>\n",Indent,Pr[i][j][k]);
						eprintf(gVerb & XML_VERB, "%s    </A>\n",Indent,numA);
						numA++;
					}
				}
			}
		}
		eprintf(gVerb & XML_VERB, "%s  </X>\n",Indent);
	}
	eprintf(gVerb & XML_VERB, "%s</locus>\n",Indent);
}


void Print_ECA_Opt_LocusProbs(const char *Ind, double ObsProbs[3][3][3], double MissObsProbs[4][4][4], double MissFakeProbs[4][4][4], sum_ped_pars *pars)
{
	int i,j,k,X;
	char Indent[5000];
	double PX=0.0;
	int Nx=5,d=3,Na=27,As=3;
	double Pr[4][4][4];
	
	
	if(pars->MissMeth==1) {
		Nx=8;
		d=4;
		Na=64;
		As=4;  /* this is the number of allele states.  We count missing as an allele state for the purposes here */
	}
	
	if(pars->MissMeth==3) {
		Nx=14;
		d=6;
		Na=64;
		As=4;  /* this is the number of allele states.  We count missing as an allele state for the purposes here */
	}
	
		
	/* some weird foolish stuff because I am using static 3-d arrays and passing both in, blah, blah */
	for(i=0;i<As;i++)  {
		for(j=0;j<As;j++)  {
			for(k=0;k<As;k++)  {
				if(pars->MissMeth>0) {
					Pr[i][j][k] = MissObsProbs[i][j][k];
				}
				else {
					Pr[i][j][k] = ObsProbs[i][j][k];
				}
			}
		}
	}
	
	
	
	/* if we want this to be bare output, we don't put the ECA_OPT_OUTPUT tag on it */
	if(gVerb & ECA_OPT_BARE_VERB) {
		sprintf(Indent,"%s",Ind);
	}
	else {
		sprintf(Indent,"ECA_OPT_OUTPUT:%s",Ind);
	}
	
	
	if(gOutputEcaOptToFile==0) {
		eprintf(gVerb & ECA_OPT_VERB, "%s--locus\n",Indent);
		eprintf(gVerb & ECA_OPT_VERB, "%s--Aprobs",Indent);
	}
	else  {
		efprintf(gVerb & ECA_OPT_VERB, gOutfile, "%s--locus\n",Indent);
		efprintf(gVerb & ECA_OPT_VERB, gOutfile, "%s--Aprobs",Indent);
	}
	
	/* first cycle over all A states and report their probs */
	
	for(i=0;i<As;i++)  {
		for(j=0;j<As;j++)  {
			for(k=0;k<As;k++)  {
				if(gOutputEcaOptToFile==0) eprintf(gVerb & ECA_OPT_VERB, " %.16f",Pr[i][j][k]);
				else efprintf(gVerb & ECA_OPT_VERB, gOutfile, " %.16f",Pr[i][j][k]);
			}
		}
	}
	
	
	if(gOutputEcaOptToFile==0) eprintf(gVerb & ECA_OPT_VERB, "\n");
	else efprintf(gVerb & ECA_OPT_VERB, gOutfile, "\n");
	
	
	
	/* finally cycle over the x's and and print out the prob of each X and the probs of the A's within those */
	for(X=0;X<Nx;X++) {
		/* first cycle over all A states to compute the prob of this X state */
		PX=0.0;
		for(i=0;i<As;i++)  {
			for(j=0;j<As;j++)  {
				for(k=0;k<As;k++)  {
					if(pars->MissMeth==0 && XVectorOfGenos(i,j,k)==X) {
						PX += Pr[i][j][k];
					}
					else if(pars->MissMeth==1 && XVectorOfGenosWithMiss1(i,j,k)==X) {
						PX += Pr[i][j][k];
					}
					else if(pars->MissMeth==3 && XVectorOfGenosWithMiss3(i,j,k)==X) {
						PX += Pr[i][j][k];
					}
				}
			}
		}
		
		if(gOutputEcaOptToFile==0) eprintf(gVerb & ECA_OPT_VERB, "%s--Xprob %.16f\n",Indent,PX);
		else efprintf(gVerb & ECA_OPT_VERB, gOutfile, "%s--Xprob %.16f\n",Indent,PX);
		
		/* now cycle again to get the fake-x-prob */
		if(pars->MissMeth>0 && MissFakeProbs != NULL)  { double PFX;
			PFX=0.0;
			for(i=0;i<As;i++)  {
				for(j=0;j<As;j++)  {
					for(k=0;k<As;k++)  {
						if(pars->MissMeth==0 && XVectorOfGenos(i,j,k)==X) {
							PFX += MissFakeProbs[i][j][k];
						}
						else if(pars->MissMeth==1 && XVectorOfGenosWithMiss1(i,j,k)==X) {
							PFX += MissFakeProbs[i][j][k];
						}
						else if(pars->MissMeth==3 && XVectorOfGenosWithMiss3(i,j,k)==X) {
							PFX += MissFakeProbs[i][j][k];
						}
					}
				}
			}
			if(gOutputEcaOptToFile==0) eprintf(gVerb & ECA_OPT_VERB, "%s--X-fake-prob %.16f\n",Indent,PFX);
			else efprintf(gVerb & ECA_OPT_VERB, gOutfile, "%s--X-fake-prob %.16f\n",Indent,PFX);
		}
		
		/* now, cycle again to print out the probabilities of the A-states */
		if(gOutputEcaOptToFile==0) eprintf(gVerb & ECA_OPT_VERB, "%s--AprobsGX",Indent);
		else efprintf(gVerb & ECA_OPT_VERB, gOutfile, "%s--AprobsGX",Indent);
		
		
		for(i=0;i<As;i++)  {
			for(j=0;j<As;j++)  {
				for(k=0;k<As;k++)  {
					if(pars->MissMeth==0 && XVectorOfGenos(i,j,k)==X) {
						if(gOutputEcaOptToFile==0) eprintf(gVerb & ECA_OPT_VERB, " %.16f",Pr[i][j][k]/PX);
						else efprintf(gVerb & ECA_OPT_VERB, gOutfile, " %.16f",Pr[i][j][k]/PX);
					}
					else if(pars->MissMeth==1 && XVectorOfGenosWithMiss1(i,j,k)==X) {
						if(gOutputEcaOptToFile==0) eprintf(gVerb & ECA_OPT_VERB, " %.16f",Pr[i][j][k]/PX);
						else efprintf(gVerb & ECA_OPT_VERB, gOutfile, " %.16f",Pr[i][j][k]/PX);
					}
					else if(pars->MissMeth==3 && XVectorOfGenosWithMiss3(i,j,k)==X) {
						if(gOutputEcaOptToFile==0) eprintf(gVerb & ECA_OPT_VERB, " %.16f",Pr[i][j][k]/PX);
						else efprintf(gVerb & ECA_OPT_VERB, gOutfile, " %.16f",Pr[i][j][k]/PX);
					}
				}
			}
		}
		if(gOutputEcaOptToFile==0) eprintf(gVerb & ECA_OPT_VERB, "\n");
		else  efprintf(gVerb & ECA_OPT_VERB, gOutfile, "\n");
	}
	
}



/* this is like the old main, so you can pass in argc and argv to it.  If outfile is null then 
 it prints to stdout.  Otherwise if outfile is not, and you are using the -q and --xml-output 
 options, then the eca_opt formatted output gets sent to the outfile and nothing goes to stdout.
*/
int snpSumPed(int argc, char **argv, FILE *outfile)
{
	int i,j,k,mus,a;
	sum_ped_pars *pars;
	int *states;
	double ProbIncompat=0.0; 
	
	
	if(outfile!=NULL) {
		gOutputEcaOptToFile = 1;
		gOutfile = outfile;
	}
	else {
		gOutputEcaOptToFile = 0;
	}
	
	pars = getSumPedOpts(argc, argv);
	
	
	
	if(pars->PrintXMLPream) {
		PrintXML_Preamble("",pars);
		PrintECA_OPT_Preamble("",pars);
	}
	if(pars->PrintXMLOutput) {
		eprintf(gVerb & XML_VERB, "XML_OUTPUT:  <collection-group name=\"%s\">\n",pars->CollectionGroupName);
		
		if( !(gVerb & ECA_OPT_BARE_VERB) ) {
			if(gOutputEcaOptToFile==0) eprintf(gVerb & ECA_OPT_VERB, "ECA_OPT_OUTPUT:    --coll-name %s\n",pars->CollectionGroupName);
			else efprintf(gVerb & ECA_OPT_VERB, gOutfile, "ECA_OPT_OUTPUT:    --coll-name %s\n",pars->CollectionGroupName);
		}
		else if(gOutputEcaOptToFile==0) {
			eprintf(gVerb & ECA_OPT_VERB, "    --coll-name %s\n",pars->CollectionGroupName);
		}
		else {
			efprintf(gVerb & ECA_OPT_VERB, gOutfile,  "    --coll-name %s\n",pars->CollectionGroupName);
		}
	}
	
	/*************************************************************************************************/
	/******** NASTY BLOCK FOR DEALING WITH A SINGLE allele frequency *********/
	
	for(a=0;a<pars->NumPs;a++) {
		/* here we set the pars->p variable successively to each one in the Ps array */
		pars->p = pars->Ps[a];
		if(gHasYouthP) {
			gYouthP = pars->YouthPs[a];
		}
		
		/* initialize some of our global variables and others  */
		ProbIncompat = 0.0;
		for(i=0;i<3;i++)  {
			for(j=0;j<3;j++)  {
				for(k=0;k<3;k++)  {
					gTrioGProbs[i][j][k] = 0.0;
					gTrioObsProbs[i][j][k] = 0.0;
				}
			}
		}
			
		/* tag the fixed individuals in the ped array */
		for(i=0;i<pars->Nfixed;i++)  {
			pars->ped[pars->fixed[i]].fixed = 1;
		}
		
		/* report some values */
		eprintf(gVerb & STANDARD_VERB, "ALLELE_FREQ : %f\n",pars->p);
		if(gHasYouthP) {
			eprintf(gVerb & STANDARD_VERB, "YOUTH_ALLELE_FREQ : %f\n",gYouthP);
		}
		eprintf(gVerb & STANDARD_VERB, "GTYPING_ERROR_RATES : ");
		for(i=0;i<pars->NumMu;i++) {
			eprintf(gVerb & STANDARD_VERB, " %f",pars->mu[i]);
		}
		eprintf(gVerb & STANDARD_VERB, "\n");
		
		/* ReportPedState(pars->ped, pars->MaxN); 
		RecursSumOver(pars->ped,1,pars->MaxN, pars->p);
		return(0); */
		
		states = (int *)ECA_CALLOC(pars->Nfixed,sizeof(int));
		gSumOverAll = 0.0;
		RecursCycleOverFixed(pars->ped, pars->Nfixed, pars->fixed, 0, states, pars->p, pars->MaxN);
		
		eprintf(gVerb & STANDARD_VERB, "SUM_OVER_ALL_STATES_FIXED_ANO_OTHERWISE : %f\n",gSumOverAll);
		
		
		/* now we print out the probs of the true trio states */
		for(i=0;i<3;i++)  {
			for(j=0;j<3;j++)  {
				for(k=0;k<3;k++)  {
					eprintf(gVerb & STANDARD_VERB, "TRIO_TRUE_STATES : %d %d %d    %f     X-State-> %d:  %s\n",i,j,k,gTrioGProbs[i][j][k],
							XVectorOfGenos(i,j,k),PrintXStateFromIdx(XVectorOfGenos(i,j,k)) );
				}
			}
		}
		
		
		/*set the mus variable to the mu for the locus a */
		mus = a;
		
		/* now compute the prob of the observed trios given gtyping error */
		ProbOfObservedStates(pars->mu[mus], gTrioGProbs, gTrioObsProbs);
		
		/* now compute the prob of the observed trios given gtyping error and missing data rates, if missing data rates were
		 supplied */
		if(pars->MissMeth>0) {
			ProbOfObservedStatesWithMiss(pars->MissRates[mus],gTrioObsProbs,gTrioMissObsProbs,0);
			ProbOfObservedStatesWithMiss(pars->MissRates[mus],gTrioObsProbs,gTrioMissFakeProbs,1);
		}
		
		/*  and report them all, and while doing that, compute the prob that the locus is incompatible in the trio */
		for(i=0;i<3;i++)  {
			for(j=0;j<3;j++)  {
				for(k=0;k<3;k++)  {
					eprintf(gVerb & STANDARD_VERB, "TRIO_OBS_STATES :  %f  : %d %d %d    %f\n",pars->mu[mus],i,j,k,gTrioObsProbs[i][j][k]);
					ProbIncompat += (double)IndicIncomp(k,j,i) * gTrioObsProbs[i][j][k];
				}
			}
		}
		
		/* report the stuff with missing data, too, if appropriate */
		if(pars->MissMeth>0) {
			for(i=0;i<4;i++)  {
				for(j=0;j<4;j++)  {
					for(k=0;k<4;k++)  {
						eprintf(gVerb & STANDARD_VERB, "TRIO_MISS_OBS_STATES :  %f  %f   %f     :    %d %d %d    %f    ",pars->MissRates[mus][0],pars->MissRates[mus][1],pars->MissRates[mus][2],
							   i,j,k,gTrioMissObsProbs[i][j][k]);
						if(pars->MissMeth==3) {
							   eprintf(gVerb & STANDARD_VERB, " X-State3-> %d:  %s \n",
									  XVectorOfGenosWithMiss3(i,j,k),PrintXStateFromIdxWithMiss3(XVectorOfGenosWithMiss3(i,j,k)));
						}
						else if(pars->MissMeth==1) {
							eprintf(gVerb & STANDARD_VERB, " X-State1-> %d:  %s\n",
							   XVectorOfGenosWithMiss1(i,j,k),PrintXStateFromIdxWithMiss1(XVectorOfGenosWithMiss1(i,j,k)));
						}
					}
				}
			}
		}
		
		
		/* report the "fake probs" with missing data, too, if appropriate */
		if(pars->MissMeth>0) {
			for(i=0;i<4;i++)  {
				for(j=0;j<4;j++)  {
					for(k=0;k<4;k++)  {
						eprintf(gVerb & STANDARD_VERB, "TRIO_MISS_FAKE_PROB_STATES :  %f  %f   %f     :    %d %d %d    %f    ",pars->MissRates[mus][0],pars->MissRates[mus][1],pars->MissRates[mus][2],
							   i,j,k,gTrioMissFakeProbs[i][j][k]);
						if(pars->MissMeth==3) {
							eprintf(gVerb & STANDARD_VERB, " X-State3-> %d:  %s \n",
								   XVectorOfGenosWithMiss3(i,j,k),PrintXStateFromIdxWithMiss3(XVectorOfGenosWithMiss3(i,j,k)));
						}
						else if(pars->MissMeth==1) {
							eprintf(gVerb & STANDARD_VERB, " X-State1-> %d:  %s\n",
								   XVectorOfGenosWithMiss1(i,j,k),PrintXStateFromIdxWithMiss1(XVectorOfGenosWithMiss1(i,j,k)));
						}
					}
				}
			}
		}
		
		
		/* finally, produce some output that we can put right into the command line for the importance sampling
			program that will do the power analysis  */  
		eprintf(gVerb & STANDARD_VERB, "OUTPUT_FOR_ENSUING_COMMAND_LINE : p %f mu %f probs ",pars->p,pars->mu[mus]);
		for(i=0;i<3;i++)  {
			for(j=0;j<3;j++)  {
				for(k=0;k<3;k++)  {
					eprintf(gVerb & STANDARD_VERB, " %.16f",gTrioObsProbs[i][j][k]);
				}
			}
		}
		eprintf(gVerb & STANDARD_VERB, "\n");
		
		/* finally, print out the prob that the locus is incompatible */
		eprintf(gVerb & STANDARD_VERB, "PROB_LOC_IS_INCOMPAT : %f\n",ProbIncompat);
		
		/* and here, if it is requested, print the XML output for each locus */
		if(pars->PrintXMLOutput) { char temp[50];
			sprintf(temp,"Locus_%d",a+1);
			PrintXML_LocusProbs("    ",gTrioObsProbs,temp);
			Print_ECA_Opt_LocusProbs("        ",gTrioObsProbs,gTrioMissObsProbs,gTrioMissFakeProbs,pars);
			
		}
		
		
	} /************************** DONE CYCLING OVER ALLE FREQS ************/
	
	if(pars->PrintXMLOutput) {
		eprintf(gVerb & XML_VERB, "XML_OUTPUT:  </collection-group>\n");
	}
	if(pars->PrintXMLClose) {
		eprintf(gVerb & XML_VERB, "XML_OUTPUT:</run>\n");
	}

	
	
	free(states);
	free(pars->File);
	free(pars->ped);
	free(pars->Ps);
	free(pars->fixed);
	if(pars->NumMu>0) {
		free(pars->mu);
	}
	if(pars->NumMiss>0) {
		for(i=0;i<pars->NumMiss;i++)  {
			free(pars->MissRates[i]);
		}
		free(pars->MissRates);
	}
	if(gHasYouthP == 1) {
		free(pars->YouthPs);
	}
	/* lastly free pars */
	free(pars);
	
	return(0);
}


/* takes in an array of true genotype probabilities, and, by summing over those, fills
out each entry in ObsProb---the array of the probabilities of observing all the different
genotypes */
void ProbOfObservedStates(double mu, double TrueProbs[3][3][3], double ObsProb[3][3][3])
{
	int i,j,k,ii,jj,kk;
	double sum;
	
	for(i=0;i<3;i++)  {
		for(j=0;j<3;j++)  {
			for(k=0;k<3;k++)  {
				sum=0.0;
				for(ii=0;ii<3;ii++)  {
				for(jj=0;jj<3;jj++)  {
				for(kk=0;kk<3;kk++)  {
					sum += ObsProbs(ii, i, mu) * ObsProbs(jj, j, mu) * ObsProbs(kk, k, mu) * TrueProbs[ii][jj][kk];
				}}}
				/* assign that to the output variable */
				ObsProb[i][j][k] = sum;
			}
		}
	}
	
}

/* takes in an array of underlying genotype probabilities which are called TrueProbs, but might
 include genotyping error, and, by summing over those, fills
 out each entry in ObsProb---the array of the probabilities of observing all the different
 genotypes when "missing" is included as state 3.  The prob of missingness is given for the first, second
 and third members of the trio separately in the array mu.  Note that memory must be allocated to 
 ObsProb outside of this function.
 
 Note that kid is the i subscript, pa is the j subscript, ma is the k subscript.
 
 So mu[0] = miss prob amongst the youths
    mu[1] = miss prob amongst the fathers
    mu[2] = miss prob amongst the mothers
 
 This applies to all missing data modes whether we are ultimately keeping track of the number of missing
 loci in each individual on just the number of loci at which at least one member is missing a locus
 in the trio.
 
 MissingIsUnity is this strange flag used to compute what the probabilities would end up being if we just
 ignored missing loci (and hence gave them "probability" of 1).  MissingIsUnity=0 gives us the actual probability
 and MissingIsUnity==1 gives us the other weird quantity that I hope relate somehow to what happens when you
 just ignore the missing data.
 
 */
void ProbOfObservedStatesWithMiss(double *mu, double TrueProbs[3][3][3], double ObsProb[4][4][4], int MissingIsUnity)
{
	int i,j,k,ii,jj,kk;
	double sum;
	
	for(i=0;i<4;i++)  {
		for(j=0;j<4;j++)  {
			for(k=0;k<4;k++)  {
				sum=0.0;
				for(ii=0;ii<3;ii++)  {
					for(jj=0;jj<3;jj++)  {
						for(kk=0;kk<3;kk++)  {
							sum += MissObsProb(ii, i, mu[0],MissingIsUnity) * MissObsProb(jj, j, mu[1],MissingIsUnity) * MissObsProb(kk, k, mu[2],MissingIsUnity) * TrueProbs[ii][jj][kk];
						}}}
				/* assign that to the output variable */
				ObsProb[i][j][k] = sum;
			}
		}
	}
	
}


/* compute the probability of genotype g2 given genotype g1 and
freq p, and K-coefficients K specifying the pairwise relatedness between
g1 and g2.  */
double PairWiseProbFunc(int g1, int g2, double p, double K[3])
{
	double result,q=1.0-p;
	
	switch(g1) {
		case 0:
			switch(g2) {
				case 0:
					result = K[0]*PofG(g2,p) + K[1]*p + K[2];
					break;
				case 1:
					result = K[0]*PofG(g2,p) + K[1]*q;
					break;
				case 2:
					result = K[0]*PofG(g2,p);
					break;
				default:
					fprintf(stderr,"Error! g2=%d which is <0 or >2 in PairWiseProbFunc().  Exiting!\n",g2);
					exit(1);
				}
				break;
		case 1:
			switch(g2) {
				case 0:
					result = K[0]*PofG(g2,p) + K[1]*.5*p;
					break;
				case 1:
					result = K[0]*PofG(g2,p) + K[1]*.5 + K[2];
					break;
				case 2:
					result = K[0]*PofG(g2,p) + K[1]*.5*q;
					break;
				default:
					fprintf(stderr,"Error! g2=%d which is <0 or >2 in PairWiseProbFunc().  Exiting!\n",g2);
					exit(1);
				}
				break;
		case 2: 
			switch(g2) {
				case 0:
					result = K[0]*PofG(g2,p);
					break;
				case 1:
					result = K[0]*PofG(g2,p) + K[1]*p;
					break;
				case 2:
					result = K[0]*PofG(g2,p) + K[1]*q + K[2];
					break;
				default:
					fprintf(stderr,"Error! g2=%d which is <0 or >2 in PairWiseProbFunc().  Exiting!\n",g2);
					exit(1);
				}
				break;
		default:
			fprintf(stderr,"Error! g1=%d which is <0 or >2 in PairWiseProbFunc().  Exiting!\n",g1);
			exit(1);
			
	}	
	
	return(result);
}

/* computes the probability of the current pedigree state.  Does this merely by a product
over nodes of the prob of the node given its parents */
double ProbOfPedState(ind *ped, int MaxN, double p)
{
	int i;
	double temp,prod=1.0;
	
	for(i=1;i<=MaxN;i++)  {
		if(ped[i].alive==1) {
			if(ped[i].ma != NULL && ped[i].pa != NULL)  { 	
				temp = CorrectTrioProb(ped[i].ma->g, ped[i].pa->g, ped[i].g);
				eprintf(gVerb & PROBPEDSTATE_VERB, "child ");
				
			}
			else if(ped[i].ma == NULL && ped[i].pa != NULL) {
				temp = OneRandomParent(ped[i].pa->g, ped[i].g, p);
				eprintf(gVerb & PROBPEDSTATE_VERB, "Missing Mom ");
			}
			else if(ped[i].ma != NULL && ped[i].pa == NULL) {
				eprintf(gVerb & PROBPEDSTATE_VERB, "MIssing Dad ");
				temp = OneRandomParent(ped[i].ma->g, ped[i].g, p);
			}
			else {  /* both parents must be NULL */
				if(ped[i].pair != NULL) {  /* if it is the origination of the pointer in a pairwise related pair */
					eprintf(gVerb & PROBPEDSTATE_VERB, "PairMember ");
					temp = PairWiseProbFunc(ped[i].pair->g,ped[i].g,p,ped[i].K);
				} 
				else {
					eprintf(gVerb & PROBPEDSTATE_VERB, "Founder ");
					if(gHasYouthP && ped[i].ID==1) {  /* here is a kluge to get the youth from a different population */
 						temp = PofG(ped[i].g, gYouthP);
					}
					else {
						temp = PofG(ped[i].g,p);
					}
				}
			}
			prod *= temp;
			
			eprintf(gVerb & PROBPEDSTATE_VERB, "INPROBPEDSTATE : %d     %f     %f\n",i,temp,prod);
		}
	}
	return(prod);
}


/* function to recursively cycle over the the genotypes of the fixed individuals in the pedigree */
void RecursCycleOverFixed(ind *ped, int Nfixed, int *fixed, int i, int *states, double p, int MaxN)
{
	int j,k;
	
	for(j=0;j<3;j++) {
	
		states[i]=j;
		
		if(i<Nfixed-1) {
			RecursCycleOverFixed(ped,Nfixed,fixed,i+1,states, p, MaxN);
		}
		
		else  {
			/* set the fixed state in the pedigree to what it should be */
			for(k=0;k<Nfixed;k++)  {
				ped[fixed[k]].g = states[k];
			}	
			/* now compute sum over all other states of the pedigree, keeping the fixed
				states fixed */
			gSumOverFixed = 0;
			RecursSumOver(ped,1,MaxN,p);
			
			/* now report this */
			eprintf(gVerb & RECCYC_VERB, "FIXED_STATES : ");
			for(k=0;k<Nfixed;k++)  {
				eprintf(gVerb & RECCYC_VERB, "%d (%d)     ",fixed[k],states[k]);
			}	
			eprintf(gVerb & RECCYC_VERB, "\n");
			eprintf(gVerb & RECCYC_VERB, "FIXED_ST_PROBS : %.10f\n",gSumOverFixed);
			
			/* down here we will store this state so long as we are dealing with just a trio of fixed individuals */
			if(Nfixed==3)  {
				gTrioGProbs[ states[0] ][ states[1] ][ states[2] ] = gSumOverFixed;
			}
			
			
			
			gSumOverAll += gSumOverFixed;
		}
	}
}

/* function to recursively cycle over the genotypes of non-fixed individuals alive in the pedigree
  and, given a set of genotype states for the fixed individuals, compute the probability of all
  the genotypes on the pedigree.  Note that the genotypes of the fixed individuals must already be
  set in the tree. */
void RecursSumOver(ind *ped, int i, int MaxN, double p)
{
	int j;
	double temp;
	
	if(ped[i].alive==0)  {
		if(i==MaxN) {
			eprintf(gVerb & RECSUMOVER_VERB, "PEDSTATE\n");
			if(gVerb & REPORTPED_VERB) ReportPedState(ped,MaxN);
			temp = ProbOfPedState(ped,MaxN,p);
			gSumOverFixed +=  temp;
			eprintf(gVerb & RECSUMOVER_VERB, "PEDPROB : %f\n",temp);
			eprintf(gVerb & RECSUMOVER_VERB, "\n\n\n");
		}
		else {
			RecursSumOver(ped,i+1,MaxN,p);
		}
		return;
	}
	if(ped[i].alive==1 && ped[i].fixed==1 ) {
		if(i==MaxN) {
			eprintf(gVerb & RECSUMOVER_VERB, "PEDSTATE\n");
			if(gVerb & REPORTPED_VERB) ReportPedState(ped,MaxN);
			temp = ProbOfPedState(ped,MaxN,p);
			gSumOverFixed +=  temp;
			eprintf(gVerb & RECSUMOVER_VERB, "PEDPROB : %f\n",temp);
			eprintf(gVerb & RECSUMOVER_VERB, "\n\n\n");
		}
		else {
			RecursSumOver(ped,i+1,MaxN,p);
		}
		return;
	}
	else if(ped[i].alive==1 && ped[i].fixed==0) {
		for(j=0;j<3;j++)  {
			ped[i].g = j;
			if(i==MaxN) {
				eprintf(gVerb & RECSUMOVER_VERB, "PEDSTATE\n");
				if(gVerb & REPORTPED_VERB) ReportPedState(ped,MaxN);
				temp = ProbOfPedState(ped,MaxN,p);
				gSumOverFixed +=  temp;
				eprintf(gVerb & RECSUMOVER_VERB, "PEDPROB : %f\n",temp);
				eprintf(gVerb & RECSUMOVER_VERB, "\n\n\n");
			}
			else {
				RecursSumOver(ped,i+1,MaxN,p);
			}
		}
		return;
	}
}


/* allocate space to and initialize pedigree stuff for M spaces. Typically we just set MaxN to 300, and that
will be more than enough! */
ind *AllocPedStruct(int MaxN)
{
	int i;
	ind *ped = (ind *)ECA_CALLOC(MaxN+1,sizeof(ind));
	
	/* initialize values */
	for(i=0;i<=MaxN;i++)  {
		ped[i].ID = i;
		ped[i].ma = NULL;
		ped[i].pa = NULL;
		ped[i].pair = NULL;
		ped[i].InPair = 0;
		ped[i].g = 0;
		ped[i].fixed = 0;
		ped[i].alive = 0;
	}
	
	return(ped);	
}


void ReportPedState(ind *ped, int MaxN)
{
	int i;
	
	eprintf(gVerb & REPORTPED_VERB, "Kid   Ma    Pa  : reported as   ID (genotype) [alive] {fixed}\n");
	for(i=0;i<=MaxN;i++) {
		if(ped[i].alive==1 ) {
			if(ped[i].pa!=NULL && ped[i].ma != NULL) {
				eprintf(gVerb & REPORTPED_VERB, "%d (%d) [%d] {%d}       %d (%d) [%d] {%d}       %d (%d) [%d] {%d} \n",
					ped[i].ID, ped[i].g, ped[i].alive, ped[i].fixed,
					ped[i].pa->ID, ped[i].pa->g, ped[i].pa->alive, ped[i].pa->fixed,
					ped[i].ma->ID, ped[i].ma->g, ped[i].ma->alive, ped[i].ma->fixed);
			}
			else if(ped[i].pa==NULL && ped[i].ma != NULL) {
				eprintf(gVerb & REPORTPED_VERB, "%d (%d) [%d] {%d}       0 (0) [0] {0}       %d (%d) [%d] {%d} \n",
					ped[i].ID, ped[i].g, ped[i].alive, ped[i].fixed,
					ped[i].ma->ID, ped[i].ma->g, ped[i].ma->alive, ped[i].ma->fixed);
			}
			else if(ped[i].pa!=NULL && ped[i].ma == NULL) {
				eprintf(gVerb & REPORTPED_VERB, "%d (%d) [%d] {%d}       %d (%d) [%d] {%d}       0 (0) [0] {0} \n",
					ped[i].ID, ped[i].g, ped[i].alive, ped[i].fixed,
					ped[i].pa->ID, ped[i].pa->g, ped[i].pa->alive, ped[i].pa->fixed);
			}
			else {
				eprintf(gVerb & REPORTPED_VERB, "%d (%d) [%d] {%d}       0 (0) [0] {0}       0 (0) [0] {0} \n",
					ped[i].ID, ped[i].g, ped[i].alive, ped[i].fixed);
			}
		}
	}
}


/* returns the prob of the genotype (0,1, or 2) given p */
double PofG(int G, double p)
{
	double q=1.0-p;
	
	switch(G) {
		case 0:
			return(p*p);
			break;
		case 1:
			return(2.0*p*q);
			break;
		case 2:
			return(q*q);
			break;
		default:
			return(-99999999.99);
	}	
}


/*
	Given true genotypes scored as 0,1, and 2, (number of copies of the 
	minor allele) for mother, father, and kid, this function returns the 
	probabilities based on mendelian segregation.  
*/
double CorrectTrioProb(int m,int f, int y) {
  if(m+f==0 && y>0) return(0.0);
  if(m+f==1 && y>1) return(0.0);
  if(m+f==2 && m!=f && y!=1) return(0.0);
  if(m+f==3 && y<1) return(0.0);
  if(m+f==4 && y<2) return(0.0);

  if(m+f==0 && y==0) return(1.0);
  if(m+f==1 && y<2) return(.5);
  if(m+f==3 && y>0) return(.5);
  if(m+f==4 && y==2) return(1.0);

  if(m+f==2 && m!=f && y==1) return(1.0);
  if(m+f==2 && m==f && y==1) return(.5);
  if(m+f==2 && m==f && y!=1) return(.25);
  
  fprintf(stderr,"Error! No Return Out of ObsProbs.  m=%d, f=%d, y=%d\nExiting...\n",m,f,y);
  exit(1);
  
  return(-99999999.99);
}


/* this gives the probability of y given m, and a random f */
double OneRandomParent(int m, int y, double p) {
  	double q = 1.0 - p;
  	
  	if( (m==0 && y==2) || (m==2 && y==0) ) return(0.0);
  	if( (m==0 && y==1) || (m==2 && y==2) ) return(q);
  	if( (m==0 && y==0) || (m==2 && y==1) ) return(p);
  	if(m==1 && y==0) return(.5*p);
  	if(m==1 && y==2) return(.5*q);
  	if(m==1 && y==1) return(.5);
  	
 	fprintf(stderr,"Error.  Got through all conditionals in OneRandomParent().\tExiting.\n\n");
 	exit(1);
 	
 	return(-99999999.99);
}



/*
	Given true genotypes scored as 0,1, and 2, (number of copies of the 
	minor allele) for mother, father, and kid, this is the indicator function
	as to whether the trio is incompatible (1) or compatible (0), at the locus.
*/
int IndicIncomp(int m, int f, int y) {
  if(m+f==0 && y>0) return(1);
  if(m+f==1 && y>1) return(1);
  if(m+f==2 && m!=f && y!=1) return(1);
  if(m+f==3 && y<1) return(1);
  if(m+f==4 && y<2) return(1);
  return(0);
}


/*
Given true genotypes "a" scored as 0,1, and 2, (number of copies of the 
	minor allele) for an individual, and observed genotypes "b" and a genotyping error
	rate mu, this returns the probability of the observed genotype given the true.
*/
double ObsProbs(int a, int b, double mu)
{
	if(a==0 && b==0) return( pow(1.0-mu,2) );
	if(a==0 && b==1) return( 2*mu*(1.0-mu) );
	if(a==0 && b==2) return( mu*mu );
	if(a==1 && b==0) return( mu*(1.0-mu) );
	if(a==1 && b==1) return( mu*mu + pow(1.0-mu,2) );
	if(a==1 && b==2) return( mu*(1.0-mu) );
	if(a==2 && b==0) return( mu*mu );
	if(a==2 && b==1) return( 2*mu*(1.0-mu) );
	if(a==2 && b==2) return( pow(1.0-mu,2) );

	fprintf(stderr,"Error! No Return Out of ObsProbs.  a=%d, b=%d, mu=%f\nExiting...\n",a,b,mu);
	exit(1);
	
	return(-999999999.99);
}


/*
 Given an underlying genotype "a" scored as 0,1, and 2, (number of copies of the 
 minor allele) observed on an individual (this can include genotyping 
 error already), and an observed genotype "b" scored as 0, 1, 2, and 3, where 
 3 is the state of "missing", and also given the probability of missingness, d,
 this returns the probability of observing the genotype (with missingness) given
 the underlying genotype (without missingness)
 
 The MissingIsUnity flag is a weird thing for computing values of genotypes that would
 be used if the missing loci were just ignored (Hence given the value of 1 so that they
 don't change any probabilities).  
 
*/
double MissObsProb(int a, int b, double d, int MissingIsUnity)
{
	if(a<0 || b<0 || a>2 || b>3)  {
		fprintf(stderr,"ERROR!  a= %d  or b=%d  out of range in MissObsProb.  Exiting...\n",a,b);
		exit(1);
	}
	if(!MissingIsUnity) {
		if(a==b) return(1-d);
		if(b==3) return(d);
	}
	else if(MissingIsUnity) {
		if(a==b) return(1.0);
		if(b==3) return(1.0);
	}
	return(0.0);
}


sum_ped_pars *getSumPedOpts(int argc, char **argv)
{

	sum_ped_pars *temp = (sum_ped_pars *)ECA_MALLOC(sizeof(sum_ped_pars));
	int only_eca_opt_outputF=0,
		freqsF=0,
		muF = 0,
		dF = 0,
		fixedF=0, 
		pedF=0,
		pairKF=0,
		kappaF=0,
		verbF=0,
		xml_pream_f=0,
		xml_output_f=0,
		close_xml_f=0,
	    smax_f=0,
		missmeth_f=0,
		youth_p_f=0;
	
	
	DECLARE_ECA_OPT_VARS;
	
	
	
	SET_OPT_WIDTH(20);
	SET_ARG_WIDTH(26);
	SET_VERSION("\nSNP_SumPed -- probability of pedigree subsets using binary loci\n\nVERSION: 1.0 beta\nAUTHOR: Eric C. Anderson (eric.anderson@noaa.gov)\nDATE: 04 JANUARY 2005\nCOPYRIGHT: No Copyright---U.S. Federal Government Work\n\n")
	SET_VERSION_HISTORY("VERSION: 1.0\nAUTHOR: Eric C. Anderson (eric.anderson@noaa.gov)\nRELEASE DATE: 30 MAY 2006\nCOPYRIGHT: No Copyright---U.S. Federal Government Work\n\n\nVERSION: 1.0 beta\nAUTHOR: Eric C. Anderson (eric.anderson@noaa.gov)\nDATE: 04 JANUARY 2005\nCOPYRIGHT: No Copyright---U.S. Federal Government Work\n\n")

	
	SET_PROGRAM_SHORT_DESCRIPTION("probability of pedigree subsets using binary loci");
	
	
	/* set some default flags for some of the variables */
	temp->p = 0.0;
	temp->mu = NULL;
	temp->File = (char *)ECA_CALLOC(5000,sizeof(char));
	temp->MaxN = 0;
	temp->PrintXMLPream = 0;
	temp->PrintXMLOutput = 0;
	temp->PrintXMLClose = 0;
	temp->Smax[0]=2; temp->Smax[1]=2; temp->Smax[2]=3;  /* default smax values */
	temp->Smax[3]=10; temp->Smax[4]=10; temp->Smax[5]=10;  /* default smax values */
	temp->MissMeth=0;
	temp->MissRates=NULL;
	temp->NumMiss=0;
	temp->NumMu = 0;
	gVerb=(STANDARD_VERB | ECA_OPT_VERB | XML_VERB);  /* by default it is primed to print lots of this stuff */
	
	
	/* allocate space to the pedigree */
	temp->ped = AllocPedStruct(300);
	

	BEGIN_OPT_LOOP
	
		if(__OptCommentLevel > 0) {
			printf("\n   ****  Command Line Switches for inputting variables  ****\n\n");
		}	
	if(OPTION(Only produce eca opt output,
			  only_eca_opt_outputF,
			  q,
			  eca-pbt-output,
			  ,
		      only produce output in ECA_OPT format for the snppit program,
			  only produce output in ECA_OPT format for the snppit program.  This suppressed COMMFILE insertion statements
			  and all other things like line tags and what not to make a clean file for pbtC (i.e. SNPPIT).  Note that to be effective
			  this must be issued before any --command-file statements. )) {
		if(ARGS_EQ(0)) {
			gVerb = (ECA_OPT_VERB | ECA_OPT_BARE_VERB);
			__OptIncludeCommFileVerbose = 0;
		}
	}
	
	
		if (  REQUIRED_OPTION(Allele Frequencies, 
			freqsF,
			p,
			minor-freq,
			[R1] [R2] ...,
			Freqs of the 0 allele at all the loci.,
			The frequencies of the 0 allele at all the loci.  Allele 0 not be the least frequent allele in the population. Required.) )	{
			if(ARGS_GEQ(1)) { int i;
				temp->NumPs = COUNT_ARGS;
				temp->Ps = (double *)ECA_CALLOC(temp->NumPs,sizeof(double));
				for(i=0;i<temp->NumPs;i++)  { 
					temp->Ps[i] = GET_DUB;
				}
			}
		}
		if( OPTION(Genotyping Error Rates,
				muF, 
				e ,
				gtyp-error, 
				[R1] [R2]. . . ,
				per-gene-copy genotyping error rate for all loci,
				Per-gene-copy genotyping error rate for all loci.  Default is zero. Optional.  Default is zero.
				) ) { int i;
			if(ALREADY_HAS(freqsF,-p/--minor-freq )) {
				temp->NumMu = COUNT_ARGS;
				if(temp->NumMu != temp->NumPs) {
					fprintf(stderr,"Error. Expecting %d mu's to be parallel to the freqs.  Reading %d instead. Exiting.\n",temp->NumPs,temp->NumMu); 
					exit(1);
				}
				temp->mu = (double *)ECA_CALLOC(temp->NumMu,sizeof(double));
				for(i=0;i<temp->NumMu;i++)  { 
					temp->mu[i] = GET_DUB;
				}
			}
		}
	
	
	
	if( OPTION(Missing Data Rates,
			   dF, 
			   d ,
			   miss-rates, 
			   [R1] [R2]. . . ,
			   missing data rates\054 one for each locus,
			   Missing data rates\054 one for each locus.  Optional.  Default is NULL.  If given then they should not all be zero\054
			   because this will be more inefficient than just not issuing the command.  Note that currently you just specify a single missing data rate for youths and mothers and
			   fathers\054 however the program is set up to store and use different rates for these different populations if that
			   feature is to be added in the future.
			   ) ) { int i;
		if(ALREADY_HAS(freqsF,-p/--minor-freq )) {
			temp->NumMiss = COUNT_ARGS;
			if(temp->NumMiss != temp->NumPs) {
				fprintf(stderr,"Error. Expecting %d mu's to be parallel to the freqs.  Reading %d instead. Exiting.\n",temp->NumPs,temp->NumMu); 
				exit(1);
			}
			temp->MissRates = (double **)ECA_CALLOC(temp->NumMiss,sizeof(double *));
			for(i=0;i<temp->NumMiss;i++)  { int j; double xxx;
				temp->MissRates[i] = (double *)calloc(3,sizeof(double));
				xxx = GET_DUB;
				for(j=0;j<3;j++)  {
					temp->MissRates[i][j] = xxx;
				}
			}
		}
	}
		
		/* DOWN HERE, I SCUTTLED THE WAY I HAD IT SET UP TO BE ABLE TO COMPUTE PROBABIILITIES ON MORE THAN 
		THREE SUBSETS IN THE INTEREST OF IT BEING ABLE TO USE FOR PARENTAGE SORTS OF PROBLEMS.  THE FUNCTIONALITY
		COULD PROBABLY BE REVIVED TO SOME EXTENT. */
		if( REQUIRED_OPTION(Fixed Genotype Indices,
				fixedF, 
				f ,
				fixed-geno, 
				[J-kid] [J-pa] [J-ma] ,
				Indices of youth/father/mother ,
				Use this option to list the indexes of the individuals that you want to specify as the members of a trio consisting of
				a putative youth a putative father and a putative mother.  The first integer should be the index 
				of the putative youth; the second the index of the putative father;  the third the index of the 
				putative mother.  The program will then
				go through all 27 possible assignments of genotypes to those three individuals and report the probability
				for each one.  Required.
				) ) { int i;
				if(ARGS_GT(0)) {
					temp->Nfixed = COUNT_ARGS;
					temp->fixed = ECA_CALLOC(temp->Nfixed,sizeof(int));
					for(i=0;i<temp->Nfixed;i++) {
						temp->fixed[i] = GET_INT;
					}
				}
		}
		if(MULT_USE_REQ_OPTION(Pedigree Command, 
				pedF,
				,
				ped,
				[J-kid] [J-pa] [J-ma],
				pedigree entry.  Individual names are positive integers.,
				Pedigree entry.  The names of the individuals must be positive integers.
				If an individual is a founder then it has both J-pa and J-ma set equal to 
				zero.  If just one parent is regarded as randomly drawn from the population then
				that parent gets a 0---i.e. either J-pa or J-ma is zero.  All the --ped statements must
				be given before any --pair statements can be given. If an individual occurs in the pedigree
				somewhere as a parent but not as a child then it is assumed to be a founder---the program
				automatically assigns 0s for both parents. Required.,
				100
			)) {
				if(ARGS_EQ(3)) { int kid, ma, pa;
					if(pairKF==0) {
						kid = GET_INT;
						ma = GET_INT;
						pa = GET_INT;
						temp->ped[kid].alive = 1;
			
						if(pa==0)
							temp->ped[kid].pa = NULL;
						else  {
							temp->ped[kid].pa = &(temp->ped[pa]);
							temp->ped[kid].pa->alive = 1;
						}
						if(ma==0)
							temp->ped[kid].ma = NULL;
						else  {
							temp->ped[kid].ma = &(temp->ped[ma]);
							temp->ped[kid].ma->alive = 1;
						}
						
						/* here we keep track of the highest number of any of them */
						if(kid > temp->MaxN) temp->MaxN = kid;
						if(ma > temp->MaxN) temp->MaxN = ma;
						if(pa > temp->MaxN) temp->MaxN = pa;
					}
				}
				else {
					fprintf(stderr,"Error!  The --ped option was used after a --pair option.  Not allowed!\n");
					OPT_ERROR
				}
		}
		if(MULT_USE_OPTION(Kappa Coefficients, 
				kappaF,
				,
				kappa,
				[J1] [J2] [R0] [R1] [R2],
				pairwise relatedness between J1 and J2,
				This allows the user to input a value for the pairwise relatedness coefficients between individuals
				J1 and J2 in the pedigree.  J1 and J2 are the integer identifiers of the two individuals which MUST
				be founders in the pedigree---i.e. they can have no non-zero parents in the pedigree.  
				R0 is the probability that J1 and J2 share zero gene copies IBD. R1 is the probability that J1 and J2 
				share exactly one gene copy IBD. R2 is the probability that
				J1 and J2 share exactly 2 gene copies identical by descent at a locus.     These
				relate to the Cotterman coefficients as  R0 = k_0; R1 = 2k_1; R2 = k_2.  R0+R1+R2 must sum to 1. Each time --kappa is used the pair J1 and
				J2 must be entirely distinct with the pair in every other use of --pair-k or --kappa., 
				100
			
			)) {
				if(ARGS_EQ(5)) { int j1, j2; double k0, k1, k2;
					j1 = GET_INT;
					j2 = GET_INT;
					k0 = GET_DUB;
					k1 = GET_DUB;
					k2 = GET_DUB;
					
					if(temp->ped[j1].ma != NULL || 
							temp->ped[j1].pa != NULL || 
							temp->ped[j2].ma != NULL || 
							temp->ped[j2].pa != NULL) {
						fprintf(stderr, "Error! --kappa option given with non-founders in pair at: \"--kappa %d %d %f %f %f\"\n",
								j1, j2, k2, k1, k0);
						OPT_ERROR;
					}
					else if(temp->ped[j1].InPair +  temp->ped[j2].InPair > 0) {
						fprintf(stderr, "Error! --kappa option given with same pair member as previously at: \"--kappa %d %d %f %f %f\"\n",
								j1, j2, k2, k1, k0);
						OPT_ERROR;
					}
					else if(k1+k2+k0 != 1.0)  {
						fprintf(stderr, "Error! Sum of R2, R1, and R0 not equal to 1.0 at: \"--kappa %d %d %f %f %f\"\n",
								j1, j2, k2, k1, k0);
						OPT_ERROR;
					}
					else {
						temp->ped[j1].InPair = 1;
						temp->ped[j2].InPair = 1;
						temp->ped[j1].pair = &(temp->ped[j2]);
						
						temp->ped[j1].K[0] = k0;
						temp->ped[j1].K[1] = k1;
						temp->ped[j1].K[2] = k2;
					}
					
				}
		}
		if(MULT_USE_OPTION(Old Style Kappa Command, 
				pairKF,
				,
				pair-k,
				[J1] [J2] [R2] [R1] [R0],
				OLD STYLE pairwise relatedness between J1 and J2,
				This allows the user to input a value for the pairwise relatedness between individuals
				J1 and J2 in the pedigree.  J1 and J2 are the integer identifiers of the two individuals which MUST
				be founders in the pedigree---i.e. they can have no non-zero parents in the pedigree.  R2 is the probability that
				J1 and J2 share exactly 2 gene copies identical by descent at a locus.  R1 is the probability that J1 and J2 
				share exactly one gene copy IBD.  R0 is the probability that J1 and J2 share zero gene copies IBD.  These
				relate to the Cotterman coefficients as R2 = k_2; R1 = 2k_1; R0 = k_0.  R2+R1+R0 must sum to 1. Each time --pair-k is used the pair J1 and
				J2 must be entirely distinct with the pair in every other use of --pair-k OR --kappa.  This --pair-k was the option that was originally in place
				for inputting pairwise relatedness into the program.  However the kappa option is now preferred because it follows the notation
				of the Genetics article more closely. , 
				100
			
			)) {
				if(ARGS_EQ(5)) { int j1, j2; double k0, k1, k2;
					j1 = GET_INT;
					j2 = GET_INT;
					k2 = GET_DUB;
					k1 = GET_DUB;
					k0 = GET_DUB;
					
					if(temp->ped[j1].ma != NULL || 
							temp->ped[j1].pa != NULL || 
							temp->ped[j2].ma != NULL || 
							temp->ped[j2].pa != NULL) {
						fprintf(stderr, "Error! --pair option given with non-founders in pair at: \"--pair %d %d %f %f %f\"\n",
								j1, j2, k2, k1, k0);
						OPT_ERROR;
					}
					else if(temp->ped[j1].InPair +  temp->ped[j2].InPair > 0) {
						fprintf(stderr, "Error! --pair option given with same pair member as previously at: \"--pair %d %d %f %f %f\"\n",
								j1, j2, k2, k1, k0);
						OPT_ERROR;
					}
					else if(k1+k2+k0 != 1.0)  {
						fprintf(stderr, "Error! Sum of R2, R1, and R0 not equal to 1.0 at: \"--pair %d %d %f %f %f\"\n",
								j1, j2, k2, k1, k0);
						OPT_ERROR;
					}
					else {
						temp->ped[j1].InPair = 1;
						temp->ped[j2].InPair = 1;
						temp->ped[j1].pair = &(temp->ped[j2]);
						
						temp->ped[j1].K[0] = k0;
						temp->ped[j1].K[1] = k1;
						temp->ped[j1].K[2] = k2;
					}
					
				}
		}
		
		if(OPTION(Verbosity Control,
			verbF,
			,
			verbose,
			,
			spout out ridiculous amounts of output,
			Giving this option results in the program spitting out all sorts of intermediate output which is 
			useful for debugging but probably is not of interest to the general user.)) {
				if(ARGS_EQ(0)) {
					gVerb = 2;
				}
		}
	
	if(REQUIRED_OPTION(Missing Data Handling Method, 
			  missmeth_f,
			  ,
			  miss-meth,
			  J1,
			  integers specifying how to treat missing loci,
			  J1 is an integer that specifies how we are going to try to keep track of missing loci at different individuals in the trios.
				J1=0 means we will not be counting any missing loci.  Accordingly the v-vectors will all have length of 3 (as must smax).
				J1=1 means that we will count the number of loci at which at least one individual in the trio has a missing locus.  the v-vectors
				thus have four components and the fourth indicates whether there are any individuals missing data at the locus.
				J1=3 means that we will keep track of the number of missing loci at kids\054 moms and dads.  This may be exceptionally memory-intensive.
				The v-vectors in this mode have 6 components.  The final three count the number of missing loci in kids\054 moms\054 and dads\054
				respectively.  This option is required and must be issued before issuing the --s-max option.
	   )) { 
		if(ARGS_EQ(1)) { int flag;
			flag = GET_INT;
			switch(flag) {
				case(0):
					temp->MissMeth=0;
					break;
				case(1):
					temp->MissMeth=1;
					break;
				case(3):
					temp->MissMeth=3;
					break;
				default:
					fprintf(stderr,"Error!  Unrecognized argument %d to option --miss-meth.  Exiting.\n",flag);
					exit(1);
			}
			
		}
	}
	
	if( OPTION(Youth Allele Freqs,
			   youth_p_f, 
			    ,
			   youth-p, 
			   [R1] [R2]. . . ,
			   allele frequencies of the minor allele in the youth population if different that parents,
			   This option can be used if the youth in a trio came from a different population than 
			   the parents.  This will only have an effect when the pedigree is such that the youth 
			   is considered a founder in the pedigree---mostly when the trio is unrelated.  This option can
			   only be given after the -p/--minor-freq option and must have the same number of arguments as that
			   option.  Note that for this to work you have to have a --ped 1 0 0 statement.  You cannot specify 1 as
			   the offspring of any indivduals\054 even if they are unrelated.
			   ) ) { int i;
		if(ALREADY_HAS(freqsF,-p/--minor-freq )) {
			gHasYouthP = 1;
			temp->NumYouthPs = COUNT_ARGS;
			if(temp->NumYouthPs != temp->NumPs) {
				fprintf(stderr,"Error. Expecting %d youth p's to be parallel to the freqs.  Reading %d instead. Exiting.\n",temp->NumPs,temp->NumYouthPs); 
				exit(1);
			}
			temp->YouthPs = (double *)ECA_CALLOC(temp->NumYouthPs,sizeof(double));
			for(i=0;i<temp->NumYouthPs;i++)  { 
				temp->YouthPs[i] = GET_DUB;
			}
		}
	}
	
	
	
	if(OPTION(SMAX Option,
			  smax_f,
			  ,
			  s-max,
			  J1 J2 J3 ...,
			  integers specifying smax.  eg. 2 2 3,
			  If this option is not issued\054 the default is 2 2 3 10 10 10.  The number of arguments should be commensurate with
			  whether or not there is considered to be missing data.  If so then 2 2 3 10 10 10 might be appropriate. So\054 for now
			  we should have only 3 or 6 arguments.)) {		
		if(ALREADY_HAS(missmeth_f, --miss-meth)) {  int numargs;
		   numargs = COUNT_ARGS;
			if( (temp->MissMeth==0 && numargs!=3)  ||
				(temp->MissMeth==1 && numargs!=4) ||
				(temp->MissMeth==3 && numargs!=6) ) {
				fprintf(stderr,"Error reading data for smax.  Should have either 3, 4 or 6 args depending on miss-meth. It has %d with MissMeth=%d.   Exiting\n",numargs,temp->MissMeth);
				exit(1);
			}
			else { int i;
				for(i=0;i<numargs;i++) {
					temp->Smax[i] = GET_INT;
				}
			}
		}
	}
	
	if(OPTION(PrintXML Preamble,
			  xml_pream_f,
			  ,
			  xml-pream,
			  S,
			  print out the xml preamble for this.  S is the run-id,
			  Giving this option results in the program spitting out the xml preamble part---number of X states\054 number
			  of components of X\054 all the v-vectors\054 etc.  After doing so\054 it will also print out the xml-output for 
			  whatever collection group is listed. )) {
		if(ARGS_EQ(1)) {
			GET_STR(temp->RunIdName);
			temp->PrintXMLPream=1;
		}
	}
	
	if(OPTION(Produce XML Output,
			  xml_output_f,
			  ,
			  xml-output,
			  S,
			  print output to xml format. S is the name of the collection group.,
			   )) {
		if(ARGS_EQ(1)) {
			GET_STR(temp->CollectionGroupName);
			temp->PrintXMLOutput=1;
		}
	}
	
	
	if(OPTION(Print Closing XML tag,
			  close_xml_f,
			  ,
			  xml-close,
			  ,
			  put the enclosing </run> tag at the very end,
			  Issue this option when this is the last time you will run snpSumPed for a particular run. )) {
		if(ARGS_EQ(0)) {
			temp->PrintXMLClose=1;
		}
	}
	

		  		
	END_OPT_LOOP
	
		
	return(temp);
}