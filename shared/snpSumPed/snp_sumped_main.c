#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <stdarg.h>
#include "MathStatRand.h"
#include "snp_sumped.h"




int main(int argc, char *argv[]) 
{
	
	FILE *dummy = NULL;
	
	/* dummy = fopen("FPRINTF_version.txt","w"); */
	
	snpSumPed(argc, argv, dummy);
	
	return(0);
}