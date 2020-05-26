
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
#include "uthash.h"
#include "pfr_read_genos.h"
#include "pbt_C_fb.h"
#include "pfr_pedigree_spec.h"
#include "pbt_highlevel.h"
#include "pbt_geno_compare.h"



/* here is a global variable to provide the ID of an individual where an error or warning occurs */
char gCurrentIndivID[MAX_INDIV_NAME_LENGTH];

char *KeywordEnumAsString(KeywordStatus K)
{
	switch (K) {
		case(PREAMBLE): return("PREAMBLE");
		case(POP): return("POP");
		case(OFFSPRING): return("OFFSPRING");
		case(MISSING_ALLELE): return("MISSING_ALLELE");
		case(POPCOLUMN_SEX): return("POPCOLUMN_SEX");
		case(POPCOLUMN_REPRO_YEARS): return("POPCOLUMN_REPRO_YEARS");
		case(POPCOLUMN_SPAWN_GROUP): return("POPCOLUMN_SPAWN_GROUP");
		case(OFFSPRINGCOLUMN_BORN_YEAR): return("OFFSPRINGCOLUMN_BORN_YEAR");
		case(OFFSPRINGCOLUMN_SAMPLE_YEAR): return("OFFSPRINGCOLUMN_SAMPLE_YEAR");
		case(OFFSPRINGCOLUMN_AGE_AT_SAMPLING): return("OFFSPRINGCOLUMN_AGE_AT_SAMPLING");
		default:
			return("UNKNOWN_KEYWORD!!!");
			break;
		return("UNKNOWN_KEYWORD!!!");
	}
}
/* read a string and return a 1 if it is a keyword.  Tell us what keyword it is in K.
 otherwise return a 0
 */
int CheckForKeyword(const char *S, KeywordStatus *K)
{	
	if(strcmp(S,"POP")==0) {  *K = POP;   return(1);  }
	if(strcmp(S,"OFFSPRING")==0) {  *K = OFFSPRING;   return(1);  }
	if(strcmp(S,"CHINOOK_AVE_POP_SIZE")==0) {  *K = CHINOOK_AVE_POP_SIZE;  return(1);  }
	if(strcmp(S,"POPCOLUMN_SEX")==0) {  *K = POPCOLUMN_SEX;   return(1);  }
	if(strcmp(S,"POPCOLUMN_REPRO_YEARS")==0) {  *K = POPCOLUMN_REPRO_YEARS;   return(1);  }
	if(strcmp(S,"POPCOLUMN_SPAWN_GROUP")==0) {  *K = POPCOLUMN_SPAWN_GROUP;   return(1);  }
	if(strcmp(S,"OFFSPRINGCOLUMN_BORN_YEAR")==0) {  *K = OFFSPRINGCOLUMN_BORN_YEAR;   return(1);  }
	if(strcmp(S,"OFFSPRINGCOLUMN_SAMPLE_YEAR")==0) {  *K = OFFSPRINGCOLUMN_SAMPLE_YEAR;   return(1);  }
	if(strcmp(S,"OFFSPRINGCOLUMN_AGE_AT_SAMPLING")==0) {  *K = OFFSPRINGCOLUMN_AGE_AT_SAMPLING;   return(1);  }
	if(strcmp(S,"MISSING_ALLELE")==0) {  *K = MISSING_ALLELE;   return(1);  }
	
	return(0);
}


SexEnum SexStrToEnum(const char *S)
{
	if(strcmp(S,"?")==0) {
		return(SEX_UNKNOWN);
	}
	if(strcmp(S,"M")==0) {
		return(MALE);
	}
	if(strcmp(S,"F")==0) {
		return(FEMALE);
	}
	
	fprintf(stderr,"ERROR! Unrecognized Sex identifier \"%s\".  Last read individual = %s  ...  Exiting.\n",S, gCurrentIndivID);
	exit(1);
	
	return(-999);
}

/* do what needs to be done after a PREAMBLE column-specifying keyword */
void HandlePreambleKeywords(KeywordStatus kw_stat,  KeywordStatus *ExtraPopCols, KeywordStatus *ExtraOffCols, int *NumExtraPopCols, int *NumExtraOffCols, int *KeywordFlags)
{
	if(kw_stat >= POPCOLSTART && kw_stat < OFFCOLSTART) { /* here we just read a pop col keyword */
		ExtraPopCols[(*NumExtraPopCols)++] = kw_stat;
	}
	else if(kw_stat >= OFFCOLSTART) {
		ExtraOffCols[(*NumExtraOffCols)++] = kw_stat;
	}
	KeywordFlags[kw_stat]++;
}


/* check to make sure that this pop name is unique and give it an idx */
void CheckAndInsertPopCollName(char *Name, coll_type_name_enum C, pfr_geno_data *P)
{
	struct pop_and_coll_name_hash_struct *s=NULL;
	char temp[100],temp2[100];
	
	if(C==POP_NAME) {
		sprintf(temp,"POP");
	}
	else if(C==OFF_COLL_NAME) {
		sprintf(temp,"OFFSPRING");	
	}
	else {
		fprintf(stderr,"Error! unrecognized enum %d for argument C in CheckAndInsertPopCollName.  Last read individual = %s  ...  Exiting\n",C, gCurrentIndivID);
		exit(1);
	}
	
	HASH_FIND_STR(P->PopCollNameHash,Name,s);
	if(s) {
		if(s->CollType==POP_NAME) {
			sprintf(temp2,"POP");
		}
		else if(s->CollType==OFF_COLL_NAME) {
			sprintf(temp2,"OFFSPRING");	
		}
		fprintf(stderr,"Error!  Just read the %s name %s, but it was previously used for an earlier %s name.  Last read individual = %s  ...  Exiting.\n",
          temp2,Name,temp, gCurrentIndivID);
		exit(1);
	}
	else {
		s = (struct pop_and_coll_name_hash_struct *)malloc(sizeof(struct pop_and_coll_name_hash_struct));
		s->CollType = C;
		sprintf(s->name,"%s",Name);
		if(s->CollType==POP_NAME) {
			s->idx = 1+P->NumPops++;
		}
		else {
			s->idx = 1+P->NumOffColls++;
		}
		HASH_ADD_STR( P->PopCollNameHash, name, s );
	}
}


/*
 Check the hash list for string.  If it already exists at the level, simply
 return its index.  If it doesn't exist at this level yet, increment the counter.
 If it doesn't exist in the hash at all then insert it into it and increment the
 appropriate counter, and, of course, return that counter.
 
 Also, record that string in the appropriate place in SGNames. 
 
 If the string is a question mark, then we just return a 0.
*/
int HandleSpawnGroupString(char *string, int level, pfr_geno_data *P)
{
	int AddNameToLevel=0,len;
	struct spawn_group_name_hash* s;
	
	if(level>MAX_SPAWN_GROUP_LEVELS-1) {
		fprintf(stderr,"Error! Spawn Group Name %s occurs at level %d which is greater than MAX_SPAWN_GROUP_LEVELS-1. Last read individual = %s  ...   Exiting!\n",
          string,level, gCurrentIndivID);
	}
	
	if(strcmp(string,"?")==0) {
		return(0);
	}
	
	
	/* first check to see if it is in the hash yet */
	HASH_FIND_STR(P->SpawnGroupNameHash,string,s);
	
	if(s==NULL) {  /* insert it into the hash */
		s = (struct spawn_group_name_hash *)malloc(sizeof(struct spawn_group_name_hash));
		s->ID = (int *)calloc(MAX_SPAWN_GROUP_LEVELS,sizeof(int));
		sprintf(s->string,"%s",string);
		AddNameToLevel=1;
		HASH_ADD_STR(P->SpawnGroupNameHash,string,s);
	}
	else {
		if(P->TopSpawnGroupInt[level]>0) {
			return(s->ID[level]);
		}
		else {
			AddNameToLevel=1;
		}
	}
	
	if(AddNameToLevel==1) {
		P->TopSpawnGroupInt[level]++;
		s->ID[level] = P->TopSpawnGroupInt[level];
		len = strlen(string);
		P->SpawnGroupNames[level][s->ID[level]] = (char *)calloc(len+1,sizeof(char));
		sprintf(P->SpawnGroupNames[level][s->ID[level]],"%s",string);
	}
	
	return(s->ID[level]);
}

/* check to make sure that this pop name is unique and give it an idx */
void CheckAndInsertIndivName(char *Name, KeywordStatus Regime, pfr_geno_data *P,int CurrPop, SexEnum Sex)
{
	struct indiv_name_hash_struct *s=NULL;
	
	if( !(Regime==POP || Regime==OFFSPRING ) ) {
		fprintf(stderr,"Error! unrecognized enum %d for argument Regime in CheckAndInsertIndivName. Last read individual = %s  ...  Exiting\n",Regime, gCurrentIndivID);
		exit(1);
	}
	
	HASH_FIND_STR(P->IndivNameHash,Name,s);
	
	
	/* if we already have that indvidual we just count this extra occurrence and make sure some metadata are consistent */
	if(s) {
	   if(Regime==POP) {
			/* check to make sure that the current pop is the same */
			if(s->parent_pop != CurrPop && s->parent_pop != -1) { /* if s->parent_pop is -1, this means it hasn't been seen before */
				fprintf(stderr,"Error!  We previously observed invividual %s in a POP with index %d. Now he is showing up in a different population.  Last read individual = %s  ... Exiting!\n",
						Name,s->parent_pop, gCurrentIndivID);
				exit(1);
			}
		   if(s->Sex != Sex) {
			   fprintf(stderr,"Warning!  We previously observed invividual %s with a SexEnum of %d. Now he is showing up with a SexEnum of %d.  It will retain its original Sex.\n",
					   Name,s->Sex,Sex);
		   }
		   /* here we have to deal with the case where the individual has only been read as an offspring so far */
		   if(s->cnt_par==0) {
			   s->parent_pop = P->NumPops;
			   s->Sex = Sex;
			   s->idx_in_pop = P->NumInPops[P->NumPops][s->Sex]++;
			   P->TotInPops[P->NumPops]++;			   
		   }
		   s->cnt_par++;
	   }
	   else {
		   if(s->offspring_coll != CurrPop  && s->offspring_coll!=-1 ) {
			   fprintf(stderr,"Error!  We previously observed invividual %s in a OFFSPRING collection with index %d. Now he is showing up in a different collection. Last read individual = %s  ...  Exiting!\n",
					   Name,s->offspring_coll, gCurrentIndivID);
			   exit(1);
		   }
		   else {
			   if(s->cnt_offs==0) {  /* this is necessary in case the individual already appeared in a POP collection */
				   s->offspring_coll = P->NumOffColls;
				   s->idx_in_offcoll = P->NumInOffColls[P->NumOffColls]++;
			   }
			   else {
				   s->offspring_coll = CurrPop;
			   }
		   }
		   s->cnt_offs++;
	   }
	}
	
	/* otherwise we put it in the hash table and figure out its position in the other data structures we will have */
	else {
		s = (struct indiv_name_hash_struct *)malloc(sizeof(struct indiv_name_hash_struct));
		sprintf(s->name,"%s",Name);
		s->parent_pop=-1;
		s->offspring_coll=-1;
		s->cnt_par = 0;
		s->cnt_offs = 0;
		if(Regime==POP) {
			s->parent_pop = P->NumPops;
			s->cnt_par++;
			s->Sex = Sex;
			s->idx_in_pop = P->NumInPops[P->NumPops][s->Sex]++;
			P->TotInPops[P->NumPops]++;
		}
		else {
			s->offspring_coll = P->NumOffColls;
			s->cnt_offs++;
			s->idx_in_offcoll = P->NumInOffColls[P->NumOffColls]++;
		}
		s->AbsIdx = P->AbsNum++;
		HASH_ADD_STR( P->IndivNameHash, name, s );
	}
}



/* 
 Return the subscript of the individual given its name and whether it is in a POP of and OFFSPRING array
 
 This also returns a pointer to the actual element in the hash table.
*/
struct indiv_name_hash_struct *SubscriptsOfName(char *Name, coll_type_name_enum CE, pfr_geno_data *P, int *a, int *b, int *c)
{
	struct indiv_name_hash_struct *s=NULL;
	
	HASH_FIND_STR(P->IndivNameHash,Name,s);
	if(s==NULL) {
		return(s);
	}
	
	
	if(CE==POP_NAME)  {
		*a = s->parent_pop;
		*b = s->Sex;
		*c = s->idx_in_pop;
	}
	else if(CE==OFF_COLL_NAME) {
		*a = s->offspring_coll;
		*b = s->idx_in_offcoll;
	}
	else {
		fprintf(stderr,"Error! Unrecognized coll_type_name_enum %d in SubscriptsOfName. Last read individual = %s  ... Exiting.\n",CE, gCurrentIndivID);
		exit(1);
	}
	return(s);
}



/* 
 given two string alleles A and B at locus i,
 
 check to see if they are missing (i.e. match the sting Missing)
 
 track the number of distinct alleles observed (up to 2) at a locus
 and return whether we are going to call the locus missing or 0,1,2
 
*/
int ReturnGeno(int i,int **AlleNamesRaw, int *NumAlleRaw, char *A, char *B, char *Missing) 
{
	
	int j,k,l;
	int a[2];
	int b[2];
	int got_it;
	
	if(strcmp(A,Missing)==0 && strcmp(B,Missing)==0) {
		return(3);  /* return "missing" for this locus */
	}
	else if((strcmp(A,Missing)==0 && strcmp(B,Missing)!=0)  || 
			(strcmp(A,Missing)!=0 && strcmp(B,Missing)==0)) {
		fprintf(stderr,"Warning: Locus %d at individual %s has only one missing allele. (A= \"%s\"  B= \"%s\".  Coerced both to missing.  No need to abort.\n",i,gCurrentIndivID,A,B);
	}
	else {
		a[0] = atoi(A);
		a[1] = atoi(B);
		
		for(j=0;j<2;j++)  {
			got_it=0;
			for(k=0;k<NumAlleRaw[i];k++)  {
				if(AlleNamesRaw[i][k]==a[j]) {
					got_it=1;
					b[j] = k;
				}
			}
			if(got_it==0) {
				AlleNamesRaw[i][NumAlleRaw[i]] = a[j];
				b[j] = NumAlleRaw[i]++;
			}
		}
		if(b[0]>1 || b[1]>1)  {
			fprintf(stderr,"Error at locus %d.  We've detected more than 2 alleles: ",i);
			for(l=0;l<NumAlleRaw[i];l++) fprintf(stderr,"%d ",AlleNamesRaw[i][l]);
			fprintf(stderr,"  Last read individual = %s  ...  Exiting.\n", gCurrentIndivID);
			exit(1);
		}
	}
			
	return(b[0]+b[1]);
}

/* make a first pass through the data and count everything up (i.e. number of males and female an unknowns in each population, etc)
 and figure out which alleles will be called 0's and which ones called 1's 
*/
pfr_geno_data *FirstPassThroughData(const char *FileName)
{
	FILE *in;
	pfr_geno_data *ret;
	char tempstr[1000], tempstr2[1000],tempIndivID[1000];
	KeywordStatus kw_stat;
	int HitPop=0,loccount=0;
	KeywordStatus Regime=PREAMBLE;  /* to record whether we are in the preamble, in the POPs section, or in the OFFSPRING sections */ 

	int i;
	char tempPopName[200], tempOffColName[200];
	coll_type_name_enum CurrentCollType;
	char PopSizeTempStr[600];
	
	sprintf(gCurrentIndivID, "Undefined");
	
	if( (in=fopen(FileName,"r"))==NULL) {
		fprintf(stderr,"Failed in FirstPassThroughData trying to read file %s. Exiting.\n",FileName);
		exit(1);
	}
	
	ret = (pfr_geno_data *)malloc(sizeof(pfr_geno_data));
	
	/* initialize some values */
	sprintf(ret->MissingAlleleString,"-9");
	ret->NumPops=-1;
	ret->NumOffColls=-1;
	ret->AbsNum=0;
	ret->IndivNameHash = NULL;
	ret->PopCollNameHash = NULL;
	ret->SpawnGroupNameHash = NULL;
	ret->NumInPops = (int **)calloc(MAXPOPS,sizeof(int *));
	for(i=0;i<MAXPOPS;i++) {
		ret->NumInPops[i] = (int *)calloc(3,sizeof(int));
	}
	ret->NumInOffColls = (int *)calloc(MAXPOPS,sizeof(int));
	ret->TotInPops = (int *)calloc(MAXPOPS,sizeof(int));
	ret->TotRowsInPops = (int *)calloc(MAXPOPS,sizeof(int));
	ret->TotRowsInColls = (int *)calloc(MAXPOPS,sizeof(int));
	ret->ExtraPopCols = (KeywordStatus *)calloc(EXCOLMAX,sizeof(KeywordStatus));
	ret->ExtraOffCols = (KeywordStatus *)calloc(EXCOLMAX,sizeof(KeywordStatus));
	ret->NumExtraPopCols = 0;
	ret->NumExtraOffCols = 0;
	ret->KeywordFlags = (int *)calloc(NUMKEYWORDS,sizeof(int));
	ret->SpawnGroupNames = (char ***)calloc(MAX_SPAWN_GROUP_LEVELS,sizeof(char **));
	for(i=0;i<MAX_SPAWN_GROUP_LEVELS;i++) {
		ret->SpawnGroupNames[i] = (char **)calloc(MAX_NUM_SPAWN_GROUP_NAMES_WITHIN_LEVELS,sizeof(char *));
		ret->SpawnGroupNames[i][0] = (char *)calloc(2,sizeof(char));  /* set the 0 ones to "?\0" */
		sprintf(ret->SpawnGroupNames[i][0],"?");
	}
	ret->TopSpawnGroupInt = (int *)calloc(MAX_SPAWN_GROUP_LEVELS,sizeof(int));
	
	/* get the number of loci */
	fscanf(in," %s",tempstr);
	if(strcmp(tempstr,"NUMLOCI")!=0) {
		fprintf(stderr,"Error! Expecting file %s to start with NUMLOCI keyword, but it started with %s. Exiting\n",FileName,tempstr);
		exit(1);
	}
	fscanf(in," %d",&(ret->NumLoci));
	
	/* now allocate space to locus names and such */
	ret->LocusNames = (char **)calloc(ret->NumLoci,sizeof(char *));
	ret->GtypErrRates = (double *)calloc(ret->NumLoci,sizeof(double));
	ret->AlleNamesRaw = (int **)calloc(ret->NumLoci,sizeof(int *));
	for(i=0;i<ret->NumLoci;i++)  {
		ret->AlleNamesRaw[i] = (int *)calloc(MAX_DISTINCT_ALLELES, sizeof(int));
	}
	ret->NumAlleRaw = (int *)calloc(ret->NumLoci,sizeof(int));
	
	
	/* now, start eating up locus names, but checking to see if they might be keywords.  POP is one of our keywords.  If we hit that, then we will exit
	 this loop.  */
	do { 
			fscanf(in," %s",tempstr);
			if(CheckForKeyword(tempstr,&kw_stat)) {
				if(kw_stat==POP) {
					if(loccount != ret->NumLoci) {
						fprintf(stderr,"Reached the first POP keyword, but have only recorded %d locus names.  Was expecting %d.  Are all locus names followed by a genotyping error rate column?  They should be.  Last recorded locus name = %s. Exiting\n",
								loccount,ret->NumLoci,ret->LocusNames[loccount-1]);
						exit(1);
					}
					else {
						Regime=POP;
						CurrentCollType=POP_NAME;
						HitPop=1;
						fscanf(in,"%s ", tempPopName);  /* eat the population name */
						CheckAndInsertPopCollName(tempPopName, CurrentCollType, ret);
						/* and if we are expecting the popsizes, get them here */
						if(ret->KeywordFlags[CHINOOK_AVE_POP_SIZE]) { int dummy;
							/* just check to make sure that the next token is a pop size in the format of:  ave_sz_1550 */
							dummy = ReadAndCheckPopSize(in,PopSizeTempStr);
						}
					}
				}
				else if(kw_stat==OFFSPRING)
				{
					fprintf(stderr,"Error! Hit the OFFSPRING keyword before any POPs.  Last read individual = %s  ... Exiting!\n", gCurrentIndivID);
					exit(1);
				}
				else if(kw_stat==MISSING_ALLELE) {
					fscanf(in," %s",(ret->MissingAlleleString));
				}
				else {
					HandlePreambleKeywords(kw_stat,  ret->ExtraPopCols, ret->ExtraOffCols, &(ret->NumExtraPopCols), &(ret->NumExtraOffCols),ret->KeywordFlags);
				}
			}
			else { 
				int len;
				len = strlen(tempstr);
				ret->LocusNames[loccount] = (char *)calloc(len+1,sizeof(char));
				sprintf(ret->LocusNames[loccount],"%s",tempstr);
				fscanf(in," %s",tempstr);
				ret->GtypErrRates[loccount]=atof(tempstr);  /* this is a total hack here---no ability to check for keywords, etc. but I couldn't figure out how to do this mo' better
															   now that it has been so long since I originally coded this stuff up */
				loccount++;
				}
	} while (HitPop==0);
	
	
	/* this is just output to check some stuff */
	#ifdef VERBOSE_FIRST_PASS_THROUGH_DATA
		for(i=0;i<loccount;i++) {
			printf("LOCUS: %s\n",ret->LocusNames[i]);
		}
		printf("POPCOLS\n");
		for(i=0;i<ret->NumExtraPopCols;i++) {
			printf("POPCOLS: %s\n",KeywordEnumAsString(ret->ExtraPopCols[i]));
		}
		printf("OFFCOLS\n");
		for(i=0;i<ret->NumExtraOffCols;i++) {
			printf("OFFCOLS: %s\n",KeywordEnumAsString(ret->ExtraOffCols[i]));
		}
	#endif
	
	
	/* now we start reading through all the POP indivs and offspring and counting them up.  */
	while (!feof(in)) {
		HitPop=0;  /* here this will tell us that we just hit a POP or an OFFSPRING and so should skip some instructions */
		if(Regime==POP) {
			
			/* read the ID */
			fscanf(in," %s",tempIndivID);
			sprintf(gCurrentIndivID,"%s",tempIndivID);
			
			#ifdef VERBOSE_FIRST_PASS_THROUGH_DATA
				printf("\nREADING_FIRST_PASS_PARENTS ");
				printf("%s ",tempIndivID);
			#endif
			if(CheckForKeyword(tempIndivID,&kw_stat)) {
				if(kw_stat==POP) {
					Regime=POP;
					CurrentCollType=POP_NAME;
					HitPop=1;
					fscanf(in,"%s ", tempPopName);  /* eat the population name */
					CheckAndInsertPopCollName(tempPopName, CurrentCollType, ret);
					/* and if we are expecting the popsizes, get them here */
					if(ret->KeywordFlags[CHINOOK_AVE_POP_SIZE]) { int dummy;
						/* just check to make sure that the next token is a pop size in the format of:  ave_sz_1550 */
						dummy = ReadAndCheckPopSize(in,PopSizeTempStr);
					}
					
				}
				else if(kw_stat==OFFSPRING) {
					Regime=OFFSPRING;
					CurrentCollType=OFF_COLL_NAME;
					HitPop=1;
					fscanf(in,"%s ", tempOffColName);  /* eat the population name */
					CheckAndInsertPopCollName(tempOffColName, CurrentCollType, ret);
					fscanf(in,"%s ", tempstr);  /* eat the possible parental pops*/
				}
				else {
					fprintf(stderr,"ERROR!  Was expecting to read an indiv ID and instead read the keyword %s.  At indiv %d in POP=%s   Last read individual = %s  ...  Exiting\n",
						tempstr,ret->TotInPops[ret->NumPops],tempPopName, gCurrentIndivID);
					exit(1);
				}
			}
			
			if(!HitPop) { SexEnum tempSex;
				tempSex = SEX_UNKNOWN;
				/* read the extra columns */
				for(i=0;i<ret->NumExtraPopCols;i++) {
					fscanf(in," %s",tempstr);
					#ifdef VERBOSE_FIRST_PASS_THROUGH_DATA
						printf("%s ",tempstr);
					#endif
					if(CheckForKeyword(tempstr,&kw_stat)) {
						fprintf(stderr,"ERROR!  Was expecting to read a value in a %s column and instead read the keyword %s.  At indiv %d in POP=%s   Last read individual = %s  ... Exiting\n",
								KeywordEnumAsString(ret->ExtraPopCols[i]),tempstr,ret->TotInPops[ret->NumPops],tempPopName, gCurrentIndivID);
						exit(1);
					}
					if(ret->ExtraPopCols[i]==POPCOLUMN_SEX) {
						tempSex = SexStrToEnum(tempstr);
					}
				}
				/* HERE MUST HASH THE FISH ID TO KNOW IF IT IS NEW OR IS A DUPLICATE */
				CheckAndInsertIndivName(tempIndivID,Regime,ret,ret->NumPops, tempSex);
				
				/* we also do know that we have read a valid row, so */
				ret->TotRowsInPops[ret->NumPops]++;
				
				/* now read the loci */
				for(i=0;i<ret->NumLoci;i++) {  int boing;
					fscanf(in," %s %s ",tempstr,tempstr2);
					#ifdef VERBOSE_FIRST_PASS_THROUGH_DATA
						printf(" %s %s ",tempstr,tempstr2);
					#endif
					if(CheckForKeyword(tempstr,&kw_stat)) {
						fprintf(stderr,"ERROR!  Was expecting to read a value in for Locus %d and instead read the keyword %s.  At indiv %d in POP=%s  Last read individual = %s  ...  Exiting\n",
								i,tempstr,ret->TotInPops[ret->NumPops],tempPopName, gCurrentIndivID);
						exit(1);
					}
					if(CheckForKeyword(tempstr2,&kw_stat)) {
						fprintf(stderr,"ERROR!  Was expecting to read a value in for Locus %d and instead read the keyword %s.  At indiv %d in POP=%s   Last read individual = %s  ... Exiting\n",
								i,tempstr2,ret->TotInPops[ret->NumPops],tempPopName, gCurrentIndivID);
						exit(1);
					}
					/* and figure out what the alleles are */
					boing = ReturnGeno(i,ret->AlleNamesRaw, ret->NumAlleRaw, tempstr, tempstr2, ret->MissingAlleleString);
					#ifdef VERBOSE_FIRST_PASS_THROUGH_DATA
						printf(" =[%d] ",boing);
					#endif
				}
				
			}
		}
		else if(Regime==OFFSPRING) {
			/* read the ID */
			fscanf(in," %s",tempIndivID);
			
			#ifdef VERBOSE_FIRST_PASS_THROUGH_DATA
				printf("\nREADING_FIRST_PASS_OFFSPRING ");
				printf("%s ",tempIndivID);
			#endif
			if(CheckForKeyword(tempIndivID,&kw_stat)) {
				if(kw_stat==POP) {
					Regime=POP;
					CurrentCollType=POP_NAME;
					HitPop=1;
					fscanf(in,"%s ", tempPopName);  /* eat the population name */
					CheckAndInsertPopCollName(tempPopName, CurrentCollType, ret);
				}
				else if(kw_stat==OFFSPRING) {
					Regime=OFFSPRING;
					CurrentCollType=OFF_COLL_NAME;
					HitPop=1;
					fscanf(in,"%s ", tempOffColName);  /* eat the population name */
					CheckAndInsertPopCollName(tempOffColName, CurrentCollType, ret);
					fscanf(in,"%s ", tempstr);  /* eat the possible parental pops*/
				}
				else {
					fprintf(stderr,"ERROR!  Was expecting to read an indiv ID and instead read the keyword %s.  At indiv %d in POP=%s   Last read individual = %s  ... Exiting\n",
							tempIndivID,ret->NumInOffColls[ret->NumOffColls],tempOffColName, gCurrentIndivID);
					exit(1);
				}
			}
			
			if(!HitPop) {
				/* read the extra columns */
				for(i=0;i<ret->NumExtraOffCols;i++) {
					fscanf(in," %s",tempstr);
					#ifdef VERBOSE_FIRST_PASS_THROUGH_DATA
						printf("%s ",tempstr);
					#endif
					if(CheckForKeyword(tempstr,&kw_stat)) {
						fprintf(stderr,"ERROR!  Was expecting to read a value in a %s column and instead read the keyword %s.  At indiv %d in OFFSPRING collection=%s  Last read individual = %s  ...  Exiting\n",
								KeywordEnumAsString(ret->ExtraOffCols[i]),tempstr,ret->NumInOffColls[ret->NumOffColls],tempOffColName, gCurrentIndivID);
						exit(1);
					}
				}
				
				/* DUDE! YOU MUST HASH THE FISH ID TO KNOW IF IT IS NEW OR IS A DUPLICATE */
				CheckAndInsertIndivName(tempIndivID,Regime,ret,ret->NumOffColls, SEX_UNKNOWN);
				ret->TotRowsInColls[ret->NumOffColls]++;
				
				/* now read the loci */
				for(i=0;i<ret->NumLoci;i++) { int boing;
					fscanf(in," %s %s ",tempstr,tempstr2);
					#ifdef VERBOSE_FIRST_PASS_THROUGH_DATA
						printf(" %s %s ",tempstr,tempstr2);
					#endif
					if(CheckForKeyword(tempstr,&kw_stat)) {
						fprintf(stderr,"ERROR!  Was expecting to read a value in for Locus %d and instead read the keyword %s.  At indiv %d in OFFSPRING collection=%s   Last read individual = %s  ... Exiting\n",
								i,tempstr,ret->NumInOffColls[ret->NumOffColls],tempOffColName, gCurrentIndivID);
						exit(1);
					}
					if(CheckForKeyword(tempstr2,&kw_stat)) {
						fprintf(stderr,"ERROR!  Was expecting to read a value in for Locus %d and instead read the keyword %s.  At indiv %d in OFFSPRING collection=%s  Last read individual = %s  ...  Exiting\n",
								i,tempstr2,ret->NumInOffColls[ret->NumOffColls],tempOffColName, gCurrentIndivID);
						exit(1);
					}
					
					/* and figure out what the alleles are */
					boing = ReturnGeno(i,ret->AlleNamesRaw, ret->NumAlleRaw, tempstr, tempstr2, ret->MissingAlleleString);
					#ifdef VERBOSE_FIRST_PASS_THROUGH_DATA
						printf(" =[%d] ",boing);
					#endif
				}
				
			}
			
		}
	}
	
	ret->NumPops++;
	ret->NumOffColls++;
	
	fclose(in);
	return(ret);
}



pfr_parent *SetUpPFR_Parent(char *Name, pfr_geno_data *P)
{
	int a,b,c,i;
	coll_type_name_enum CE = POP_NAME;
	pfr_parent *ret;
	struct indiv_name_hash_struct *s;
	
	if((s=SubscriptsOfName(Name, CE, P, &a, &b, &c))==NULL) {
		fprintf(stderr,"Error!  Indiv Name %s not found in hash table during second pass through data! Exiting\n",Name);
		exit(1);
	}
	#ifdef VERBOSE_FIRST_PASS_THROUGH_DATA
		printf(" ([%d][%d][%d]) ",a,b,c);
	#endif
	ret = &(P->Pops[a][b][c]);  /* this sets the pointer to the place in memory where this indiv will live */
	
	/* here we put a bunch of values in there and allocate some memory as need be */
	ret->Idx = c;
	ret->Pop = a;
	ret->Sex = s->Sex;
	ret->AbsIdx = s->AbsIdx;
	ret->Name = (char *)calloc(strlen(Name)+1,sizeof(char));
	sprintf(ret->Name,"%s",Name);
	if(ret->NumGenos==0) {
		ret->geno = (char **)calloc(s->cnt_par,sizeof(char *));
		for(i=0;i<s->cnt_par;i++)  {
			ret->geno[i] = (char *)calloc(P->NumLoci,sizeof(char));
		}
	}
	ret->NumMissingLoci = (int *)calloc(MAX_REGENO_TIMES,sizeof(int));
	return(ret);
}



pfr_offspring *SetUpPFR_Offspring(char *Name, pfr_geno_data *P)
{
	int a,b,c,i;
	coll_type_name_enum CE = OFF_COLL_NAME;
	pfr_offspring *ret;
	struct indiv_name_hash_struct *s;
	
	if((s=SubscriptsOfName(Name, CE, P, &a, &b, &c))==NULL) {
		fprintf(stderr,"Error!  Indiv Name %s not found in hash table during second pass through data! Exiting\n",Name);
		exit(1);
	}
	#ifdef VERBOSE_FIRST_PASS_THROUGH_DATA
		printf(" ([%d][%d]) ",a,b);
	#endif
	ret = &(P->Offs[a][b]);  /* this sets the pointer to the place in memory where this indiv will live */
	
	/* here we put a bunch of values in there and allocate some memory as need be */
	ret->Idx = b;
	ret->OffColl = a;
	ret->AbsIdx = s->AbsIdx;
	ret->Name = (char *)calloc(strlen(Name)+1,sizeof(char));
	sprintf(ret->Name,"%s",Name);
	if(ret->NumGenos==0) {
		ret->geno = (char **)calloc(s->cnt_offs,sizeof(char *));
		for(i=0;i<s->cnt_offs;i++)  {
			ret->geno[i] = (char *)calloc(P->NumLoci,sizeof(char));
		}
	}
	ret->MendComps = (indiv_matchers *)malloc(sizeof(indiv_matchers));
	for(i=0;i<3;i++)  {
		ret->MendComps->Parents[i]=NULL;
		ret->MendComps->NumParents[i] = 0;
		ret->MendComps->Num_MI_Parents[i]=NULL;
		ret->MendComps->ParPairs=NULL;
		ret->MendComps->NumParPairs = 0;
		ret->MendComps->NumParPairsMendelian = 0;
		ret->MendComps->NumParPairsMend_LogL = 0;
		ret->MendComps->NumParPairsMend_LogL_and_Rank = 0;
	}
	ret->NumMissingLoci = (int *)calloc(MAX_REGENO_TIMES,sizeof(int));
	
	return(ret);
}

/* 
 Given a string of nonegative numbers in ascending order separated 
 by ,'s or -'s, (like "1-5,7,10-15"), this function will
 return an array of lenth length of those numbers in increasing order.
 
 MaxNums is the maximum number of integers allowed in the range.
*/
int *pfr_StringToArrayOfInts(char *Str, int MaxNums, int *length)
{
	char *temp,*thedash;
	int lo=-1,hi=-1,i,now=0,old=-1;
	int *temp_ints = (int *)calloc(MaxNums,sizeof(int));
	int *result;
	int idx=0;
	
	temp = (char *)ECA_CALLOC(strlen(Str)+10,sizeof(char));
	
	sprintf(temp,"%s",Str);
	
	temp=strtok(temp,",");
	
	while(temp) {
		
		if(strchr(temp,'-'))  {
			thedash = strchr(temp,'-');
			hi = atoi(thedash+1);
			*thedash = '\0';
			lo = atoi(temp);
			if(lo>=hi) {
				fprintf(stderr,"Error in StringToIntArray().  Number string %s includes a range that is non-increasing.  Last read individual = %s  ... \n",Str, gCurrentIndivID);
				exit(1);
			}
			for(i=lo;i<=hi;i++)  {
				temp_ints[idx++] = i;
			}
		}
		else {	
			now = atoi(temp);
			if(now<=hi || now<=old) {
				fprintf(stderr,"Error in StringToIntArray().  Number string %s includes a singleton that is non-increasing. Last read individual = %s  ... \n",Str, gCurrentIndivID);
				exit(1);
			}
			temp_ints[idx++] = now;
			old=now;
		}
		temp=strtok(NULL,",");
	}
	
	/* now transfer all that */
	result = (int *)calloc(idx,sizeof(int));
	for(i=0;i<idx;i++)  {
		result[i] = temp_ints[i];
	}
	free(temp_ints);
	*length = idx;
	return(result);
}


/* tokenizes the spawner group string on commas to allow individuals to
 * belong to multiple spawner groups
 */
void TokenizeSpawnerGroupString(char *str, int *num_toks_out, char tokens[MAX_NUM_SPAWN_GROUP_PER_INDIV][MAX_SPAWN_GROUP_NAME_LENGTH]) {
  char *token;
  int num_toks = 0;
  const char s[2] = ","; /* the delimiter */
	
	/* get the first token */
	token = strtok(str, s);
	
	/* walk through other tokens */
	while( token != NULL ) {
	  sprintf(tokens[num_toks], "%s", token);
	  
	  token = strtok(NULL, s);
	  
	  num_toks++;
	}
	
	*num_toks_out = num_toks;
	
}



/* 
 Converts a string to a pfr_year_range.  
*/
pfr_year_range *CreatePFR_YearRange(char *S)
{
	int *temp,len,i;
	pfr_year_range *ret;
	
	/* first get back all the values in it */
	temp = pfr_StringToArrayOfInts(S, MAX_NUM_VALS_IN_RANGE, &len);
	
	ret = (pfr_year_range *)malloc(sizeof(pfr_year_range));
	ret->Lo = temp[0];
	ret->Len = temp[len-1]-temp[0]+1;
	ret->bits = (char *)calloc(ret->Len, sizeof(char));
	for(i=0;i<len;i++) {
		ret->bits[ temp[i] - ret->Lo ]++;
	}
	
	return(ret);
}


/* given a string S slurped up in the context of column type K, handle this and put it into
 the offspring pointed to by R 
*/
void CollectExtraColsOfOffspring(KeywordStatus K, char *S, pfr_offspring *R)
{
	switch (K) {
		case(OFFSPRINGCOLUMN_BORN_YEAR):
			R->BornYears = pfr_StringToArrayOfInts(S, MAX_NUM_VALS_IN_RANGE, &(R->NumBornYears));
			#ifdef VERBOSE_FIRST_PASS_THROUGH_DATA
		{ int i;
				printf(" ("); 
				for(i=0;i<R->NumBornYears;i++)  {
					printf("%d ",R->BornYears[i]);
				}
				printf(") ");
				}
			#endif
			break;
		case(OFFSPRINGCOLUMN_SAMPLE_YEAR):
			R->SampledYear = atoi(S);
			#ifdef VERBOSE_FIRST_PASS_THROUGH_DATA
				printf(" (%d) ",R->SampledYear);
			#endif
			break;
		case(OFFSPRINGCOLUMN_AGE_AT_SAMPLING):
			R->AgesAtSampling = pfr_StringToArrayOfInts(S, MAX_NUM_VALS_IN_RANGE, &(R->NumAgesAtSampling));
			#ifdef VERBOSE_FIRST_PASS_THROUGH_DATA
			{	int i;
				printf(" ("); 
				for(i=0;i<R->NumAgesAtSampling;i++) {
					printf("%d ",R->AgesAtSampling[i]);
				}
				printf(") ");
			}
			#endif
			break;
		default:
			fprintf(stderr,"Unrecognized Offspring Column Type %d in CollectExtraColsOfOffspring(). Exiting.\n",K);
			exit(1);
			break;
	}
}


/* given a string S slurped up in the context of column type K, handle this and put it into
 the parent pointed to by R 
 */
void CollectExtraColsOfPopInds(KeywordStatus K, char *S, pfr_parent *R, pfr_geno_data *P)
{
  
  int NumSG_Tokens, i;
  char SG_Tokens[MAX_SPAWN_GROUP_NAME_LENGTH][MAX_SPAWN_GROUP_NAME_LENGTH];
  
  
	switch (K) {
		case(POPCOLUMN_SEX):
			switch(S[0]) {
				case('F'):
					R->Sex = MALE;
					break;
				case('M'):
					R->Sex = FEMALE;
					break;
				case('?'):
					R->Sex = SEX_UNKNOWN;
					break;
				default:
					fprintf(stderr,"Unknown sex identifier %c having read string %s in CollectExtraColsOfPopInds(). Exiting.\n",S[0],S);
					exit(1);
					break;
			}
			#ifdef VERBOSE_SECOND_PASS_THROUGH_DATA
				printf(" (%d) ",R->Sex);
			#endif
			break;
		case(POPCOLUMN_REPRO_YEARS):
			R->ReproYears = CreatePFR_YearRange(S);
			#ifdef VERBOSE_SECOND_PASS_THROUGH_DATA
			{	int i;
				printf(" (lo=%d Len=%d | ",R->ReproYears->Lo,R->ReproYears->Len);  
				for(i=0;i<R->ReproYears->Len;i++) { 
					printf("%d ",R->ReproYears->bits[i]); printf(") ");
				}
			}
			#endif
			break;
		case(POPCOLUMN_SPAWN_GROUP):
			/* Tokenize the spawner group string on commas and store that in our temp character array.
			 * Then allocate memory for the integer equivalents of those spawner groups and assign them.
			 */ 
			TokenizeSpawnerGroupString(S, &NumSG_Tokens, SG_Tokens);
		  
			R->SpawnGroup = (int *)calloc(NumSG_Tokens,sizeof(int));
			for(i=0;i<NumSG_Tokens;i++) {
			  R->SpawnGroup[i] = HandleSpawnGroupString(SG_Tokens[i], 0, P);
			  printf("SpawnerGroup Assignment %d for indiv %s is %d from token %s\n", i, R->Name, R->SpawnGroup[i], SG_Tokens[i]);
			}
			R->NumSpawnGroups = NumSG_Tokens;
			
			#ifdef VERBOSE_SECOND_PASS_THROUGH_DATA 
				printf(" (%d<-->%s) ",R->SpawnGroup[0],P->SpawnGroupNames[0][R->SpawnGroup[0]]);
			#endif
			break;
		default:
			fprintf(stderr,"Unrecognized Offspring Column Type %d in CollectExtraColsOfOffspring(). Exiting.\n",K);
			exit(1);
			break;
	}
}



/* 
If BornYears is not given but Sampling Time and Sampling Age are, this function computes born years 
Returns NumBornYears as an output variable.  Check to make sure that it has everything it needs to compute
 this stuff before using this function.
*/
int *ComputeBornYears(int SampleYear, int *SampleAges, int NumSampleAges, int *NumBornYears)
{
	int i,lo,hi,*temp,*ret=NULL,cnt=0,n=0;
	
	lo = SampleYear - SampleAges[NumSampleAges-1];
	hi = SampleYear - SampleAges[0];
	
	temp=(int *)calloc(hi-lo+1,sizeof(int));
	for(i=0;i<NumSampleAges;i++)  {
		temp[ SampleYear-SampleAges[i]-lo]++;
	}
	for(i=0;i<hi-lo+1;i++)  {
		if(temp[i]) cnt++;
	}
	ret=(int *)calloc(cnt,sizeof(int));
	for(i=0;i<hi-lo+1;i++)  {
		if(temp[i]) ret[n++] = i+lo;
	}
	*NumBornYears = n;
	free(temp);
	return(ret);
}



/* parse the string of possible pops (pps) and fill out the appropriate spots in the
 OffCollPossPops array.  CollIdx is the index of the offspring collection whose data 
 is currently being read.
 
 Recall that pps should either be a "?", signifying all pops are possible, or a comma-separated
 string of possible pops.
 */
void HandleOffCollPossPops(pfr_geno_data *P,char *pps, int CollIdx)
{
	int i;
	char *p;
	struct pop_and_coll_name_hash_struct *s;
	
	#ifdef VERBOSE_SECOND_PASS_THROUGH_DATA
		printf("OFF_COLL_POSS_POPS: OffColl= %d (%s) Str= \"%s\".  Results in array: ",CollIdx, P->OffCollNames[CollIdx],pps);
	#endif
	
	if(strcmp(pps,"?")==0) {
		for(i=0;i<P->NumPops;i++)  {
			P->OffCollPossPops[CollIdx][i] = 1;
			#ifdef VERBOSE_SECOND_PASS_THROUGH_DATA
				printf(" %d", P->OffCollPossPops[CollIdx][i]);
			#endif
		}
		#ifdef VERBOSE_SECOND_PASS_THROUGH_DATA
			printf("\n");
		#endif
		return;
	}
	
	else {
		p = strtok(pps,",");
		if(p==NULL) {
			fprintf(stderr,"Error! Bad news dude.  p is null in HandleOffCollPossPops when it shouldn't be.  pps= \"%s\",  CollIdx= %d.   Last read individual = %s  ... Exiting\n",
           pps,CollIdx, gCurrentIndivID);
			exit(1);
		}
		while(p) {
			HASH_FIND_STR(P->PopCollNameHash,p,s);
			if(s && s->CollType==POP_NAME) {
				i = s->idx;
				P->OffCollPossPops[CollIdx][i]++;
			}
			else {
				fprintf(stderr,"Error!  While parsing OffCollsPossPop string \"%s\", working on token \"%s\".  That token is not a recognized Population Name. Last read individual = %s  ...  Exiting\n",
						pps,p, gCurrentIndivID);
				exit(1);
			}
			p = strtok(NULL,",");
		}
		
		for(i=0;i<P->NumPops;i++)  {
			#ifdef VERBOSE_SECOND_PASS_THROUGH_DATA
				printf(" %d", P->OffCollPossPops[CollIdx][i]);
			#endif
		}
		#ifdef VERBOSE_SECOND_PASS_THROUGH_DATA
			printf("\n");
		#endif
	}
	
}




/* ret should have been set by the first pass through the data.
 
 This allocates more memory where needed and then zips through the data one more
 time doing less error checking, but it puts all the data where it needs to be. 
 */
void CollectDataOnSecondPass(pfr_geno_data *ret, const char *FileName)
{
	int i,j,k;
	FILE *in;
	char tempstr[1000], tempstr2[1000],tempIndivID[1000],tempPopName[1000],tempParentalPoss[MAXPOPS*MAX_POP_COLL_NAME_LENGTH+1];
	int HitPop=0;
	coll_type_name_enum CE;
	pfr_offspring *pfo;
	pfr_parent *pfp;
	char PopSizeTempStr[600];
	
	/* initialize some variables */
	ret->AveragePopSizes = NULL;
	
	
	/* first allocate memory to where we will put all the geno structs, etc */
	ret->Pops = (pfr_parent ***)calloc(ret->NumPops,sizeof(pfr_parent **));
	for(i=0;i<ret->NumPops;i++)  {
		ret->Pops[i] = (pfr_parent **)calloc(3,sizeof(pfr_parent *));
		for(j=0;j<3;j++)  {
			#ifdef VERBOSE_SECOND_PASS_THROUGH_DATA
				printf("SECOND_PASS:  Alloc: Pop %d  Sex %d  NumInPops= %d\n",i,j,ret->NumInPops[i][j]);
			#endif
			if(ret->NumInPops[i][j]>0) {
				ret->Pops[i][j] = (pfr_parent *)calloc(ret->NumInPops[i][j],sizeof(pfr_parent));
				
				for(k=0;k<ret->NumInPops[i][j];k++)  {
					ret->Pops[i][j][k].NumGenos = 0;
				  
				  /* We allocate to variable number of spawner groups while reading, now, but leave this
				   * here to ensure that something is allocated to it whether or not the spawner
				   * groups are being used.
				   */
					ret->Pops[i][j][k].SpawnGroup = (int *)calloc(MAX_SPAWN_GROUP_LEVELS,sizeof(int));
					
					/* we also set the number of spawn groups to one by default so that our
					 * code will evaluate pairs even if no spawner groups info is given.
					 */
					ret->Pops[i][j][k].NumSpawnGroups = 1;
				}
				
			}
		}
	}
	ret->Offs = (pfr_offspring **)calloc(ret->NumOffColls,sizeof(pfr_offspring *));
	for(i=0;i<ret->NumOffColls;i++)  {
		ret->Offs[i] = (pfr_offspring *)calloc(ret->NumInOffColls[i],sizeof(pfr_offspring));
		for(j=0;j<ret->NumInOffColls[i];j++)  {
			ret->Offs[i][j].NumGenos = 0;
		}
	}
	/* then allocate some memory for pop names, OffCollPossPops, etc */
	ret->PopNames = (char **)calloc(ret->NumPops,sizeof(char *));
	ret->OffCollNames = (char **)calloc(ret->NumOffColls,sizeof(char *));
	ret->OffCollPossPops = (int **)calloc(ret->NumOffColls,sizeof(int *));
	for(i=0;i<ret->NumOffColls;i++)  {
		ret->OffCollPossPops[i] = (int *)calloc(ret->NumPops,sizeof(int));
	}
	
	/* open the file */
	if( (in=fopen(FileName,"r"))==NULL) {
		fprintf(stderr,"Failed in FirstPassThroughData trying to read file %s. Exiting.\n",FileName);
		exit(1);
	}
	
	/* now we eat words until we hit the first POP */
	do {
		fscanf(in," %s", tempstr);
		if(strcmp(tempstr,"POP")==0) {
			HitPop=1;
		}
	} while (HitPop==0);
	
	
	/* now we cycle over the POPs and within each of those we cycle over the number of rows, and within each of those we cycle over the columns */
	/* currently I am just scanning and printing to make sure that this is working correctly */
	CE=POP_NAME;
	for(i=0;i<ret->NumPops;i++)  { 
		/* first get the name of the population */
		fscanf(in," %s",tempPopName);
		ret->PopNames[i] = (char *)calloc(strlen(tempPopName)+1,sizeof(char));
		#ifdef VERBOSE_SECOND_PASS_THROUGH_DATA
			printf("SECOND_PASS:  %s\n",tempPopName);
		#endif
		sprintf(ret->PopNames[i],"%s",tempPopName);
		/* and if we are expecting the popsizes, get them here */
		if(ret->KeywordFlags[CHINOOK_AVE_POP_SIZE]) {
			if(ret->AveragePopSizes == NULL) {
				ret->AveragePopSizes = (int *)calloc(MAXPOPS,sizeof(int));
			}
			/* just check to make sure that the next token is a pop size in the format of:  ave_sz_1550 */
			ret->AveragePopSizes[i] = ReadAndCheckPopSize(in,PopSizeTempStr);
		}
		

		for(j=0;j<ret->TotRowsInPops[i];j++)  {
			pfp=NULL;
			fscanf(in," %s",tempIndivID);
			#ifdef VERBOSE_SECOND_PASS_THROUGH_DATA
				printf("SECOND_PASS:  %s   ",tempIndivID);
			#endif
			
			
			/* now set it up if this is the first time we have seen it on the first pass */
			pfp=SetUpPFR_Parent(tempIndivID, ret);
			
			
			for(k=0;k<ret->NumExtraPopCols;k++)  {
				fscanf(in," %s",tempstr);
				#ifdef VERBOSE_SECOND_PASS_THROUGH_DATA
					printf(" %s ",tempstr);
				#endif
				if(pfp->NumGenos==0) {  /* only collect on first appearance of genotype */
					CollectExtraColsOfPopInds(ret->ExtraPopCols[k],tempstr,pfp,ret);
				}
			}
			
			for(k=0;k<ret->NumLoci;k++)  { char boing;
				fscanf(in," %s %s ",tempstr,tempstr2);
				boing = (char)ReturnGeno(k,ret->AlleNamesRaw, ret->NumAlleRaw, tempstr, tempstr2, ret->MissingAlleleString);
				#ifdef VERBOSE_SECOND_PASS_THROUGH_DATA
					printf(" %s %s ",tempstr,tempstr2);
					printf(" =[%d]  ",boing);
				#endif
				pfp->geno[pfp->NumGenos][k] = boing;
				pfp->NumMissingLoci[pfp->NumGenos] += (boing==3);  /* if boing is 3 it means the locus is missing */
			}
			#ifdef VERBOSE_SECOND_PASS_THROUGH_DATA
				printf("NumMissingLoci= %d ",pfp->NumMissingLoci[pfp->NumGenos]);
				printf("\n");
			#endif
			/* here we increment p->NumGenos */
			pfp->NumGenos++;
			
			
		}
		
		/* if you still have POPs left over, eat the POP keyword */
		if(i<ret->NumPops-1) {
			fscanf(in," %s",tempstr);
			#ifdef VERBOSE_SECOND_PASS_THROUGH_DATA
				printf("SECOND_PASS:  Just ate a %s!\n",tempstr);
			#endif
		}
	}
	
	
	
	/* Here we cycle over the OFFSPRINGS and within each of those we cycle over the number of rows, and within each of those we cycle over the columns */
	/* currently I am just scanning and printing to make sure that this is working correctly */
	CE=OFF_COLL_NAME;
	for(i=0;i<ret->NumOffColls;i++)  {
		
		/* first eat the OFFSPRING keyword */
		fscanf(in," %s",tempstr);
		if(strcmp(tempstr,"OFFSPRING")!=0) {
			fprintf(stderr,"Error!  Expected to read \"OFFSPRING\" on second pass through data on i=%d.  Exiting\n",i);
			exit(1);
		}
		
		
		
		/* first get the name of the collection */
		fscanf(in," %s",tempPopName);
		ret->OffCollNames[i] = (char *)calloc(strlen(tempPopName)+1,sizeof(char));
		sprintf(ret->OffCollNames[i],"%s",tempPopName);
		
		/* then get the name of the parent POPs that could be parental to these OFFSRPING */
		fscanf(in," %s",tempParentalPoss);
		
		#ifdef VERBOSE_SECOND_PASS_THROUGH_DATA
			printf("SECOND_PASS:  %s\n",tempPopName);
			printf("SECOND_PASS:  %s\n",tempParentalPoss);
		#endif
		HandleOffCollPossPops(ret,tempParentalPoss,i);
		
		/* I need to put in a function here that processes these possible parent pops */
		
		for(j=0;j<ret->TotRowsInColls[i];j++)  {
			fscanf(in," %s",tempIndivID);
			#ifdef VERBOSE_SECOND_PASS_THROUGH_DATA
				printf("SECOND_PASS:  %s   ",tempIndivID);
			#endif
			
			/* now set it up if this is the first time we have seen it on the first pass */
			pfo=SetUpPFR_Offspring(tempIndivID, ret);
			
			
			for(k=0;k<ret->NumExtraOffCols;k++)  {
				fscanf(in," %s",tempstr);
				#ifdef VERBOSE_SECOND_PASS_THROUGH_DATA
					printf(" %s ",tempstr);
				#endif
				if(pfo->NumGenos==0) {  /* only get these extra column values for the first appearance of the individual */
					CollectExtraColsOfOffspring(ret->ExtraOffCols[k], tempstr, pfo);
				}
			}
			
			/* if we weren't given "BORN YEARS" but we do have the year sampled and its age, compute the born years from that */
			if( pfo->NumGenos==0 && !ret->KeywordFlags[OFFSPRINGCOLUMN_BORN_YEAR] && (ret->KeywordFlags[OFFSPRINGCOLUMN_SAMPLE_YEAR] && ret->KeywordFlags[OFFSPRINGCOLUMN_AGE_AT_SAMPLING]) ) {
				pfo->BornYears = ComputeBornYears(pfo->SampledYear, pfo->AgesAtSampling, pfo->NumAgesAtSampling, &(pfo->NumBornYears));
				#ifdef VERBOSE_SECOND_PASS_THROUGH_DATA
				{ int x;
					printf(" (CBY: ");
					for(x=0;x<pfo->NumBornYears;x++)  {
						printf("%d ",pfo->BornYears[x]);
					}
					printf(") ");
				}
				#endif
			}
			
			for(k=0;k<ret->NumLoci;k++)  { char boing;
				fscanf(in," %s %s ",tempstr,tempstr2);
				boing = (char)ReturnGeno(k,ret->AlleNamesRaw, ret->NumAlleRaw, tempstr, tempstr2, ret->MissingAlleleString);
				
				#ifdef VERBOSE_SECOND_PASS_THROUGH_DATA 
					printf(" %s %s ",tempstr,tempstr2);
					printf(" =[%d]  ",boing);
				#endif
				
				pfo->geno[pfo->NumGenos][k] = boing;
				pfo->NumMissingLoci[pfo->NumGenos] += (boing==3);  /* if boing is 3 it means the locus is missing */
			}
			#ifdef VERBOSE_SECOND_PASS_THROUGH_DATA
				printf("NumMissingLoci= %d ",pfo->NumMissingLoci[pfo->NumGenos]);
				printf("\n");
			#endif
			
			/* here we increment pfo->NumGenos */
			pfo->NumGenos++;
			
		}
	}
	
	
		
	fclose(in);
}

void PrintFirstPassSummaryOfPopsCollsAndIndivs(pfr_geno_data *PFR)
{
	
	struct pop_and_coll_name_hash_struct *s;
	struct indiv_name_hash_struct *I;

	
	printf("\n\nPOPS_AND_COLLECTIONS\n");
	for(s=PFR->PopCollNameHash; s != NULL; s=s->hh.next) {
        printf("popcoll enum %d    id %d      name %s\n", s->CollType,s->idx, s->name);
    }
	printf("\n\nINDIVS\n");
	for(I=PFR->IndivNameHash; I != NULL; I=I->hh.next) {
        printf("indiv  xPOP %d     xOFF %d      AbsIdx %d     ParentPop %d    OffsColl %d    name %s     PopIdx: [%d][%d][%d]    OffCollIdx:  [%d][%d]\n",
			   I->cnt_par,I->cnt_offs,I->AbsIdx,I->parent_pop,I->offspring_coll,I->name,
			   I->parent_pop,I->Sex,I->idx_in_pop,
			   I->offspring_coll,I->idx_in_offcoll);
    }
}

void PrintSummaryOfInputData(FILE *s, pfr_geno_data *P)
{
	int i;
	char str[]="BASIC_DATA_SUMMARY:";
	
	fprintf(s,"%s NumberOfPOP_Collections %d\n",str,P->NumPops);
	for(i=0;i<P->NumPops;i++)  {
		fprintf(s,"%s POP_Collection %d    Name= %s   NumMales= %d     NumFemales= %d     NumUnknownSex= %d   Total= %d\n",
			   str,i,P->PopNames[i],P->NumInPops[i][MALE],P->NumInPops[i][FEMALE],P->NumInPops[i][SEX_UNKNOWN],
				P->NumInPops[i][MALE]+P->NumInPops[i][FEMALE]+P->NumInPops[i][SEX_UNKNOWN]);
	}
	fprintf(s,"%s NumberOfOffspring_Collections %d\n",str,P->NumOffColls);
	for(i=0;i<P->NumOffColls;i++)  {
		fprintf(s,"%s OFFSPRING_Collection %d    Name= %s   NumIndivs= %d\n",
			   str,i,P->OffCollNames[i],P->NumInOffColls[i]);
	}
}

/* counts allele freqs and missing data in P and records it in P */
void CountAlleles(pfr_geno_data *P)
{
	int i,j,k,l,m;
	double ng;
	
	/*  allocate data and add stuff up */
	P->AlleleCounts = (double ***)calloc(P->NumPops,sizeof(double **));
	for(i=0;i<P->NumPops;i++)  {
		P->AlleleCounts[i] = (double **)calloc(P->NumLoci,sizeof(double *));
		for(j=0;j<P->NumLoci;j++)  {
			P->AlleleCounts[i][j] = (double *)calloc(3,sizeof(double));
			
			/* cycle over the males females and sex unknown indivs and add up the counts */
			for(k=0;k<3;k++)  {
				for(l=0;l<P->NumInPops[i][k];l++)  {
					ng = (double)P->Pops[i][k][l].NumGenos;
					for(m=0;m<P->Pops[i][k][l].NumGenos;m++)  {  /* note: take the total weight of all duplicately genotyped loci in an individual to be 1.0 */
						switch(P->Pops[i][k][l].geno[m][j]) {
							case(0):
								P->AlleleCounts[i][j][ALLELE_0] += 2.0/ng;
								break;
							case(1):
								P->AlleleCounts[i][j][ALLELE_0] += 1.0/ng;
								P->AlleleCounts[i][j][ALLELE_1] += 1.0/ng;
								break;
							case(2):
								P->AlleleCounts[i][j][ALLELE_1] += 2.0/ng;
								break;
							case(3):
								P->AlleleCounts[i][j][ALLELE_MISSING] += 1.0/ng;  /* note here that this counts the number of missing loci, *not* the number of missing alleles */
								break;
							default:
								fprintf(stderr,"Error! Unrecognized genotype %d at individual %s (duplicate m= %d) locus %d while counting allele frequencies.  Exiting\n",
										P->Pops[i][k][l].geno[m][j], P->Pops[i][k][l].Name, m, j);
								exit(1);
								
								
						}
					}
				}
			}			
		}
	}
	
	
}


/* tells us how the allelic identifiers correspond to the internal representation of them */
void SummarizeAllelicTypes(FILE *s, pfr_geno_data *P)
{
	int i,j;
	char str[]="ALLELE_NAMES_SUMMARY:";
	
	
	for(i=0;i<P->NumLoci;i++)  {
		fprintf(s,"%s  %s   ",str,P->LocusNames[i]);
		for(j=0;j<P->NumAlleRaw[i];j++)  {
			fprintf(s,"%d= %d    ",j,P->AlleNamesRaw[i][j]);
		}
		if(P->NumAlleRaw[i]==1) {
			fprintf(s," 1= xxx   ");
		}
		if(P->NumAlleRaw[i]<1 || P->NumAlleRaw[i]>2) {
			fprintf(stderr,"Error! NumAlleRaw[%d] is %d which is <1 or >2.  Exiting!\n",i,P->NumAlleRaw[i]);
			exit(1);
		}
		fprintf(s,"\n");
	}
}



/* tells us how the allelic identifiers correspond to the internal representation of them */
void SummarizeLocusNameAndGtypRates(FILE *s, pfr_geno_data *P)
{
	int i;
	char str[]="LOCUS_NAMES_AND_GTYP_ERR_RATES_SUMMARY:";
	
	
	
	for(i=0;i<P->NumLoci;i++)  {
		fprintf(s,"%s  %s    number= %d     gtyp_err_rate=  %f\n",str,P->LocusNames[i],i+1,P->GtypErrRates[i]);
	}
}





void PrintSummaryOfAllelicCountsAndFreqs(FILE *s, pfr_geno_data *P)
{
	int h,i;
	char str[]="ALLELE_COUNTSANDFREQS_SUMMARY:";
	
	for(h=0;h<P->NumPops;h++)  {
		for(i=0;i<P->NumLoci;i++)  {
			fprintf(s,"%s  %s   %s   ",str,P->PopNames[h],P->LocusNames[i]);
			fprintf(s,"  0=%d  --> %.2f  (%f)  ",P->AlleNamesRaw[i][0], P->AlleleCounts[h][i][ALLELE_0],P->AlleFreqs[h][i][ALLELE_0]);
			if(P->NumAlleRaw[i]==1) {
				fprintf(s,"  1=xxx  --> %.2f (%f)  ", P->AlleleCounts[h][i][ALLELE_1],P->AlleFreqs[h][i][ALLELE_1]);
			}
			else {
				fprintf(s,"  1=%d  --> %.2f  (%f)  ",P->AlleNamesRaw[i][1], P->AlleleCounts[h][i][ALLELE_1],P->AlleFreqs[h][i][ALLELE_1]);
			}
			fprintf(s,"Missing -->  %.2f",P->AlleleCounts[h][i][ALLELE_MISSING]);
			fprintf(s,"\n");
		}
	}
}

/* this just computes the posterior mean allele freqs assuming a beta(1/2.1/2) prior and the allele counts.
 Note that you must have already counted the alleles first! */
void ComputeAlleleFreqsFromCounts(FILE *s, pfr_geno_data *P)
{
	int i,j;
	double tot;
	
	/*  allocate data first */
	P->AlleFreqs = (double ***)calloc(P->NumPops,sizeof(double **));
	for(i=0;i<P->NumPops;i++)  {
		P->AlleFreqs[i] = (double **)calloc(P->NumLoci,sizeof(double *));
		for(j=0;j<P->NumLoci;j++)  {
			P->AlleFreqs[i][j] = (double *)calloc(3,sizeof(double));
			
			tot = P->AlleleCounts[i][j][ALLELE_0] + P->AlleleCounts[i][j][ALLELE_1] + 1.0;
			P->AlleFreqs[i][j][ALLELE_0] = (P->AlleleCounts[i][j][ALLELE_0]+.5)/tot;
			P->AlleFreqs[i][j][ALLELE_1] = (P->AlleleCounts[i][j][ALLELE_1]+.5)/tot;
			
			/* down here we want to have frequencies of missing data in each population, assuming a beta 1/2,1/2 prior on missing. */
			tot = 0.5 * tot + P->AlleleCounts[i][j][ALLELE_MISSING] + 0.5;
			P->AlleFreqs[i][j][ALLELE_MISSING] = (P->AlleleCounts[i][j][ALLELE_MISSING]+0.5)/tot;
			fprintf(s,"MISSING_DATA_REPORT: Pop= %d  Locus= %d  MissFreq= %f\n",i,j,P->AlleFreqs[i][j][ALLELE_MISSING]);
			
		}
	}
}


/* a simple function that reads the next word off of stream in and 
 then tests to make sure that it starts with "ave_sz_".  If it does
 not it exits with an error message.  Otherwise it returns the integer
 that follows ave_sz_ .   You must have already allocated space to S */
int ReadAndCheckPopSize(FILE *in, char *S)
{
	int ret = 0;
	char temp[500] = "";
	
	fscanf(in," %s",S);
	strncpy(temp,S,7);  /* now temp should be "ave_sz_" */
	 
	if(strcmp(temp,"ave_sz_")!=0) {
		fprintf(stderr,"Error reading in a popsize when using the CHINOOK_AVE_POP_SIZE method.  Expected to read a population size XXX in the format: ave_sz_XXX, but instead read: \"%s\".  Exiting...\n",S);
		exit(1);
	}
	else {
		strcpy(temp,S+7);
		ret = atoi(temp);
		if(ret<=0) {
			fprintf(stderr,"Error reading in a popsize when using the CHINOOK_AVE_POP_SIZE method.  Read in \"%s\", which gives an average popsize <= 0.  Exiting...\n",S);
			exit(1);
		}
	}
		
	return(ret);
	
}
