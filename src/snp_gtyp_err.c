/*
 *  snp_gtyp_err.c
 *  pfr
 *
 *  Created by Eric C. Anderson on 8/2/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

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


typedef struct {
	double CE;  /* the constant per-allele error rate */
	double *VE;  /* if error rates are variable, then they are stored in here */
	int NumL;  /* is this is 0, then we use CE, otherwise we use VE */
} gtyperr_pars;



gtyperr_pars *GetGtypErrOptions(int argc, char *argv[])
{
	int ceF=0;
	gtyperr_pars *ret = (gtyperr_pars *)malloc(sizeof(gtyperr_pars));
	
	DECLARE_ECA_OPT_VARS;
	
	
	/* some defaults */
	SET_OPT_WIDTH(25);
	SET_ARG_WIDTH(34);
	SET_VERSION("\nsnp_gtyp_err --\n\nVERSION: 1.0 Beta\nAUTHOR: Eric C. Anderson (eric.anderson@noaa.gov)\nDATE: 24 July 2006\nCOPYRIGHT: None -- US Govt Employee Work\n\n");
	SET_DESCRIBE("\nsnp_gtyp_err -- make genotyping errors in pfr data --\n");

	/* some default values */
	ret->NumL = 0;

	BEGIN_OPT_LOOP 	 
	
		if(OPTION(
			ceF,
			c,
			const-error,
			R1,
			genotyping error rate common to all loci,
			))
		{
			if(ARGS_EQ(1)) {
				ret->CE = GET_DUB;
			}
		}
	
	END_OPT_LOOP

	return(ret);
}



void PrintMutatedX(char x, gtyperr_pars *GEP, int l)
{
	double r;
	long NumErrs;
	
	if(GEP->NumL==0) {
		r = GEP->CE;
	}
	else {
		r = GEP->VE[l];
	}
	
	/* now choose how many errors there will be at the locus.  This is a binomial r.v.
	of two trials with success prob r */
	NumErrs = ignbin(2,(float)r);
	
	if(NumErrs>0) {  /* you get in here because the first gene copy has an error */
		switch(x) {
			case('0'):
				if(NumErrs==2) {
					printf("2");
				}
				else {
					printf("1");
				}
				break;
			case('2'):
				if(NumErrs==2) {
					printf("0");
				}
				else  {
					printf("1");
				}
				break;
			case('1'):
				if(NumErrs==2) { 
					printf("1");  /* no change if both have errors */
				}
				else if(ranf()<.5) {
					printf("2");
				}
				else {
					printf("0");
				}
				break;
			default:
				fprintf(stderr,"Unknown SNP type %c in PrintMutatedX().  Exiting...",x);
				exit(1);
		}
	}
	else {
		printf("%c",x);
	}
}

int main(int argc, char *argv[])
{
	int i,j;
	int N,L;
	char x,tempstr[2000];
	gtyperr_pars *GEP;
	
	
	
	GEP = GetGtypErrOptions(argc,argv);

	
	SeedFromFile("snp_gtyp_err_seeds");
	
	scanf(" %d",&N);
	printf("%d\n",N);
	
	for(i=0;i<N;i++) {
		scanf(" %s",tempstr);
		printf("%s ",tempstr);
		
		scanf(" %s",tempstr);
		L=strlen(tempstr);
		for(j=0;j<L;j++)  {
			x = tempstr[j];
			PrintMutatedX(x,GEP,j);
		}
		printf("\n");
	}
	
	SeedToFile("snp_gtyp_err_seeds");
	
	return(0);
}