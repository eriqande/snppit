

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
#include "pfr_forback.h"
#include "pfr_utils.h"

/*!
Maximum length of an individual's identifier.
*/
#define MAX_ID_LENGTH 100
#define MAX_FILE_PATH_LENGTH 5000
#define MAX_NUM_COLLECTIONS 5000


/* two different names for the same region in memory */
typedef union {
	vector signed int    v;
	signed int		  elem[4];
} sint_vector;


typedef union {
	vector unsigned char    v;
	unsigned char		  elem[16];
} uchar_vector;

typedef union {
	vector unsigned int    v;
	unsigned int		  elem[4];
} uint_vector;

typedef enum {
	Ma,
	Pa,
	Kid
} categ;

typedef struct {
	char *id; 
	int pop;  /*!< The zero-based index of the individual's population. (Determined by the order in which population files are read into the program)   */
	categ cat;  /*!< tells whether the individual is a mother, father, or kid */
	int num;  /*!< The index of the individual in the pop/category */
	uint_vector a; /*!< 128 bits holding the type of one of the two alleles   */
	uint_vector s; /*!< 128 bits: 1 means homozygous (i.e. s is for "same") and 0 is heterozygous  */
	uint_vector g; /*!< 128 bits telling us if we have got data or not.  a 1 means "got it" and 0 means data missing at the locus   */
	int *G; /*!< array of chars.  Each will be -1, 0, 1, 2 as an int.  This is inefficient on storage, but I will use this till
					I figure out a fast way to convert a,s,g to an element in a genotype freq array. */
} ind_struct;


typedef struct {
	char MaFile[MAX_FILE_PATH_LENGTH];
	char PaFile[MAX_FILE_PATH_LENGTH];
	char KidFile[MAX_FILE_PATH_LENGTH];
	
	int ListOnlyMode;

	pfr_forback_data* D;  /*!< This is here to hold freqs of genotype configurations  */
	
	lam_struct **Lambdas; /*!< Array of pointers to lambdas that we will want to compute */
	int NumLams;  /*!< Length of Lambdas array */
	
} pfr_opts;

typedef struct {
	int NumMa;
	int NumPa;
	int NumKids;
	
	ind_struct *Moms;
	ind_struct *Dads;
	ind_struct *Kids;
	
} geno_data;

int gNumLoci = -1; /*!< used for making sure all the genotype lengths are correct */

/*! Returns a vector of 16 unsigned chars, each one showing the number of 1-bits in each corresponding
8-bit segment of v.  This code is from a posting on the web by Ian Ollman: 
http://www.simdtech.org/altivec/archive/msg?list_name=altivec&monthdir=200404&msg=msg00006.html */
vector unsigned char bitcount( vector unsigned char v )
{
       const vector unsigned char table1 = (vector unsigned char)
               ( 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 );
        const vector unsigned char table2 = vec_add( table1, vec_splat_u8(1));

       vector unsigned char hiBits = vec_sr( v, vec_splat_u8(5 ) );

       return vec_add( vec_perm( table1, table2, v),
                       vec_perm( table1, table2, hiBits ) );
}


/*! this is a simple function to read 32 characters in a string, and from those,
pack bits into an unsigned int
\param str The string---should be composed of 0, 1, 2, or n or N or ? (n,N, and ? denotes missing)
\param a output parameter.  It is the indicator of the type of the first allele
\param s output parameter.  This holds bits which are 0 if heterozygous and 1 if homozyg.
\param g output parameter.  This holds bits which are 1 if that locus is genotyped, 0 if not.  
 */
void StringTo32Bits(char *str, unsigned int *a, unsigned int *s, unsigned int *g)
{
	int i;
	int Len = strlen(str);
	unsigned int one = 1;
	
	/* initialize */
	*a = 0;
	*s = 0;
	*g = 0;
	
	for(i=0;i<32;i++)  {
		/* shift everybody */
		(*a) <<= 1; 
		(*s) <<= 1;
		(*g) <<= 1;

		/* then give '1' values in certain cases */
		if(i<Len) {
			switch(str[i]) {
				case('0'):
					(*s) |= one;
					(*g) |= one;
					break;
				case('1'): 
					(*g) |= one;
					break;
				case('2'):
					(*a) |= one;  /* shift a one onto a */
					(*s) |= one;
					(*g) |= one;
					break;
				case('n'):
				case('?'):
				case('N'):
					break;
				default:
					fprintf(stderr,"Error! Unrecognized character \"%c\" in string \"%s\" at position %d (starting from 0). Exiting...\n",str[i],str,i);
					exit(1);
			}
					
		}
		
	}
}




/*! print out the 0's and 1's in a 32-bit unsigned int.  The bits that get printed "further to the left" are the 
"more significant" ones.
\param x the unsigned int whose bits we want to print  */
void Print32Bits(unsigned int x) {
	unsigned int mask = 0x80000000;
	int count;
	
	for (count=0; count<32; count++,mask>>=1) {
		if( (x & mask) == 0) {
			printf("0");
		}
		else {
			printf("1");
		}
	}
}




/*! This computes the log-likelihood (Lambda) for a trio of individuals (y,m,f) = (youth, mother, father)
	given frequencies of the 27 different mapatri trio genotype configurations which get passed in as 
	pointers to pointers of pfrfb_locus structs (could probably make a more efficient structure for it here)
	
	L is the number of loci.  
	
	Note that for each individual, I amm imagining passing in the G field of the ind_struct.
	
	NOTE!  THIS IS NOT CURRENTLY DEVELOPED TO DEAL WITH MISSING DATA!
	
*/	
double Lambda_MaPaTri(int *y, int *m, int *f, pfrfb_locus **Num, pfrfb_locus **Denom, int L)
{
	int i;
	double LN=1.0, LD=1.0;  /* numerator and denominators */
	int pick;
	double ln, ld;
	
	for(i=0;i<L;i++)  {
		pick = 9*y[i] + 3*f[i] + m[i];
		ln = Num[i]->PA[pick];
		ld = Denom[i]->PA[pick];

		LN *= ln;
		LD *= ld;
		
//printf("locus %d  pick= %d ln= %f  ld= %f  LN= %e   LD= %e \n",i+1,pick,ln,ld,LN,LD);
	}	
	return( log(LN)-log(LD) );
}


/*! Convert a char '0', '1', '2', or 'n' to an int 0, 1, 2, or -1, respectively */
int GenoCharToInt(char G)
{
	switch(G) {
		case('0'):
			return(0);
			break;
		case('1'):
			return(1);
			break;
		case('2'):
			return(2);
			break;
		case('n'):
			return(-1);
			break;
		default:
			fprintf(stderr,"Error! Unknown argument %c to GenoCharToInt().  Exiting\n",G);
			exit(1);
	}
}

/*! Print out all the information about an individual */
void PrintInd(ind_struct x)
{
	int i;
	
	
	printf("IND: %s     %d %d %d\n",x.id,x.pop,x.cat,x.num);
	printf("a_var:");
	for(i=0;i<4;i++)  {
		printf(" ");
		Print32Bits(x.a.elem[i]);
	}
	printf("\n");
	
	printf("s_var:");
	for(i=0;i<4;i++)  {
		printf(" ");
		Print32Bits(x.s.elem[i]);
	}
	printf("\n");
	
	printf("g_var:");
	for(i=0;i<4;i++)  {
		printf(" ");
		Print32Bits(x.g.elem[i]);
	}
	printf("\n");
	
	printf("AsInt: ");
	for(i=0;i<gNumLoci;i++)  {
		printf("%d",x.G[i]);
	}
	printf("\n");
}

/*!
This function opens the file "File" to read in a series of genotypes.  The format of the 
file should be such that it starts with an integer giving the number of individuals in the file.  
Then for each individual there should be two strings.  The first is an identifier for that individual
and the second is a string of 012n giving its genotype.  

This function can be called in ListOnlyMode which means merely go through the individuals and
list their identifiers, so that we can hash out a list of indexes of individuals that may have
mated with one another (using, for example, awk).  

\param File the name of the file.
\param ListOnlyMode 1 means only list the names and the indices. 0 means collect data into structures.
\param pop The zero-based index of the population these individuals are coming from.
\param cat The category of these individuals --- i.e. are they ma's, pa's, or kids.
\param Num Output parameter for number of individuals read.
\return pointer to the array which gets memory allocated to it, and stuff put in it.  
*/
ind_struct *ReadGenosFromFile(char *File, int ListOnlyMode, int pop, categ cat, int *Num)
{
	int i,j,N;
	char *munch=(char *)ECA_CALLOC(140,sizeof(char));
	char tname[MAX_ID_LENGTH];
	FILE *in;
	int Glength;
	ind_struct *ret;
	
	if( (in=fopen(File,"r"))==NULL) {
		fprintf(stderr, "Error! Unable to open file \"%s\" to read in data.  Exiting\n",File);
		exit(1);
	}
	
	fscanf(in," %d",&N);
	*Num = N;
	printf("Preparing to read %d individuals from file %s\n",N,File);
	
	if(!ListOnlyMode)  {
		ret = (ind_struct *)ECA_CALLOC(N,sizeof(ind_struct));
	}
	
	for(i=0;i<N;i++)  {
		fscanf(in," %s %s",tname,munch);
		if(ListOnlyMode)  {
			printf("LIST: %s  %d  %d  %d\n",tname,pop,cat,i);
		}
		else {  /* in this case we store it to memory */
			/* get the name, pop, cat, and num */
			ret[i].id=(char *)ECA_CALLOC(MAX_ID_LENGTH,sizeof(char));
			sprintf(ret[i].id,"%s",tname);
			ret[i].pop = pop;
			ret[i].cat = cat;
			ret[i].num = i;
			
			/* check that the genotypes are all the same length */
			Glength = strlen(munch);
			if(gNumLoci<0) {
				gNumLoci = Glength;
			}
			if(gNumLoci != Glength) {
				fprintf(stderr,"Error! Indiv %s has genotype of %d loci, but those before it had genotypes of %d. (Or %d loci were specified in a use of the -p/--pop option). Exiting.\n",tname,Glength,gNumLoci,gNumLoci);
				exit(1);
			}
			if(Glength>128) {
				fprintf(stderr,"Error! Indiv %s has genotype of %d loci, which is greater than 128. Exiting.\n",tname,Glength);
				exit(1);
			}
			
			/* cycle over the 4 32-bit unsigned ints that make up our vector register, and stick them in there.  */
			for(j=0;j<4;j++)  {
				if(j <  ceil(Glength/32.0) ) {
					StringTo32Bits(munch+j*32, &(ret[i].a.elem[j]), &(ret[i].s.elem[j]),&(ret[i].g.elem[j]));
				}
				else {  /* if we have read past the end of the string, we just fill all with 0's */
					ret[i].a.elem[j] = 0;
					ret[i].s.elem[j] = 0;
					ret[i].g.elem[j] = 0;
				}
			}
			
			/* now put the genotype in the individual in an integer representation (held as chars, since there are very 
			few values that go in each element */
			ret[i].G = (int *)calloc(gNumLoci,sizeof(int));
			//printf("ID: %s  ",ret[i].id);
			for(j=0;j<gNumLoci;j++)  {
				ret[i].G[j] = GenoCharToInt(munch[j]);
				//printf("%d",ret[i].G[j]);
			}
			//printf("\n");
			
		}
		
/*PrintInd(ret[i]);*/
	}
	
	
	return(ret);
	
}


geno_data GetGenoData(pfr_opts PF) 
{
	geno_data gd;

	gd.Moms = ReadGenosFromFile(PF.MaFile, PF.ListOnlyMode, 0, Ma, &(gd.NumMa));
	gd.Dads = ReadGenosFromFile(PF.PaFile, PF.ListOnlyMode, 0, Ma, &(gd.NumPa));
	gd.Kids = ReadGenosFromFile(PF.KidFile, PF.ListOnlyMode, 0, Ma, &(gd.NumKids));
	
	return(gd);

}



/*! This function counts the number of mismatches between a m, f, and y
\todo I think I can avoid some of the intermediate variables here */
int CountMismTrio(ind_struct m, ind_struct f, ind_struct y)
{
	sint_vector S;
	vector signed int zero;
	
	/* initialize zero */
	zero = vec_splat_s32(0);
	
	/* do the Boolean operations */
	S.v = vec_sums( (vector signed char) vec_sum4s(bitcount( (vector unsigned char)vec_and(  vec_and( m.g.v , vec_and(f.g.v,y.g.v)) ,  
			vec_or(   vec_and(y.s.v,vec_or(vec_and( vec_xor(m.a.v,y.a.v) , m.s.v) , vec_and( vec_xor(f.a.v,y.a.v) , f.s.v)))  , 
					vec_andc(vec_and(vec_and(f.a.v,m.a.v),vec_and(f.s.v,m.s.v)),y.s.v)
				  )
			    ) ),zero) , zero) ;
				
	/* then count the bits of D */
	//C.v = bitcount( (vector unsigned char)(D.v) ) ;
		
		
	/* and here we count up the number of one bits over all parts */
	//S.v = vec_sums( (vector signed char) vec_sum4s( C.v, zero ) , zero );
		
	return(S.elem[3]);
}

/*! This function counts the number of mismatches between a single putative parent and putative offspring */
int CountMismPair(ind_struct f, ind_struct y) 
{
	
	sint_vector S;
	vector signed int zero;
	
	/* initialize zero */
	zero = vec_splat_s32(0);
	
	S.v = vec_sums( (vector signed char)vec_sum4s( 
			bitcount( (vector signed char) 
				vec_and(  vec_and(f.g.v,y.g.v) , vec_and( vec_xor(f.a.v,y.a.v), vec_and(f.s.v,y.s.v) ) )
			) ,  zero), zero);
			
	return(S.elem[3]);
}


pfr_opts GetPFROptions(int argc, char *argv[])
{
	int mothers_f=0,
		fathers_f=0,
		kids_f=0,
		Popf=0,
		Af=0,
		Mf=0,
		lamoptsf=0;
	pfr_opts PF;
	pfr_forback_data *ret = (pfr_forback_data *)malloc(sizeof(pfr_forback_data));
	DECLARE_ECA_OPT_VARS;
	
	/* set some default values and do some allocation, etc */
	ret->Pops = (pfr_pop_struct **)calloc(MAX_NUM_POPS,sizeof(pfr_pop_struct *));
	ret->MC = (mc_struct **)calloc(MAX_NUM_MCOPTS,sizeof(mc_struct *));
	ret->L = 0;
	ret->NumPops = 0;
	ret->NumMC = 0;

	
	

	/* some defaults and mem alloc for PF */
	PF.ListOnlyMode = 0;
	PF.Lambdas = (lam_struct **)calloc(MAX_NUM_LAMOPTS,sizeof(lam_struct));

	SET_OPT_WIDTH(25);
	SET_ARG_WIDTH(34);
	SET_VERSION("\npfr -- vectorized bitwise comparison of SNP genotypes --\n\nVERSION: 1.0 Beta\nAUTHOR: Eric C. Anderson (eric.anderson@noaa.gov)\nDATE: 2 January 2006\nCOPYRIGHT: None -- US Govt Employee Work\n\n");
	SET_DESCRIBE("\npfr -- vectorized bitwise comparison of SNP genotypes --\n");

	BEGIN_OPT_LOOP 	 
		if(REQUIRED_OPTION(
			Af,
			a,
			analysis-type,
			S,
			name of analysis,
			S should be one of: mapatri patri matri tri.  mapatri means that you screened father pairs and mother pairs
			and then screened the full trios composed of compatible pairs.  patri means you screened father pairs and then 
			trios that they could be members in.  matri means you screened mother pairs and then 
			trios that they could be members in.  tri means you just screened trios and no pairs of any sort. ) ) 
		{

			if(ARGS_EQ(1)) { char tempStr[100];
				GET_STR(tempStr);
				if(strcmp(tempStr,"mapatri")==0) {
					ret->aType = MAPATRI;
				}
				else if(strcmp(tempStr,"patri")==0) {
					ret->aType = PATRI;
				}
				else if(strcmp(tempStr,"matri")==0) {
					ret->aType = MATRI;
				}
				else if(strcmp(tempStr,"tri")==0) {
					ret->aType = TRI;
				}
				else {
					fprintf(stderr,"String \"%s\" is an invalid option to -a/--analysis-type.  Expecting mapatri, patri, matri, or tri. Exiting...\n",tempStr);
					exit(1);
				}
				printf("READ_DATA: Arg \"%s\" to -a/--analysis-type sets variable aType to %d \n",tempStr,ret->aType);
				/* Here we do a minor bit of setting a few more aType related variables */
				ret->Xdim = XdimOfaType(ret->aType);
				ret->XNum = XNumOfaType(ret->aType);
				ret->Xkeys = SetXKeys(ret->aType);
				ret->NumA = NumAofaType(ret->aType);
			}
		}
		
		if(MULT_USE_REQ_OPTION(
			Popf,
			p,
			pop,
			S Loc1 R1 ... RN Loc2 R1 ...,
			supply information about genotype freqs in a pop,
			Use this option to input genotype frequencies for all the loci in the study.  The input format is somewhat complex.
				S is the name given to this population.  It might be something like PopA_NullHypothesis or PopA_Parental.  The point
				is that what is mean by a population can include a hypothesis about relatedness.  After the name of the population the
				probabilities of the different combinations of possible genotypes among the individuals of interest are put on the command line.
				For analysis of trios there are 27 such frequencies for each locus---i.e. RN is R27 in the argument list.  If one includes missing data as a genotypic state then
				there are 64 such states.  If you are doing an individual matching analysis there are only 9 states.  The locus information must
				go in order.  For the X-th locus the frequencies must be preceded by LocX. For example: Loc1 or Loc122.  This tag is
				case sensitive.  The total number of loci will be calculated by the program from the length of the command line.  The number
				of loci must be identical for each population.  , 
				MAX_NUM_POPS))
		{ int i,j;
			if(ALREADY_HAS(Af,-a/--analysis-type)) {
				if(ARGS_GEQ(ret->NumA+2)) { int tempL; int tempArgs; char tempPop[MAX_POP_NAME]; char tempLocName[50]; char tempLocFlag[500];
					tempArgs = COUNT_ARGS;
					GET_STR(tempPop);
					if( (tempArgs - 1) % (ret->NumA + 1) != 0) {
						fprintf(stderr,"Error!  You have %d args to option -p/--pop with (pop-name \"%s\").  The number of args minus one should be perfecty divisible by %d. Exiting...\n",tempArgs,tempPop,ret->NumA+1); 
						exit(1);
					}
					tempL = (tempArgs - 1) / (ret->NumA+1);
					if(ret->L==0) {
						ret->L = tempL;
						printf("READ_DATA: Setting Number of Loci to %d\n",ret->L);
					}
					else if(ret->L != tempL) {
						fprintf(stderr,"Error! Previous calls to -p/--pop set values for %d loci.  The call with pop-name %s is setting values for %d loci.  Exiting...\n",ret->L,tempPop,tempL);
						exit(1);
					}
					if(gNumLoci>=0) {
						if(ret->L != gNumLoci) {
							fprintf(stderr,"Error!  -p/--pop option implies %d loci, but one of the -m, -f, or -y options already established there are %d loci. Exiting...\n",ret->L,gNumLoci);
							exit(1);
						}
					}
					else { /* in this case we will set the gNumLoci here */
						gNumLoci = ret->L;
					}
					/* down here and we are good to go */
					/* allocate some memory */
					ret->Pops[ret->NumPops] = (pfr_pop_struct *)malloc(sizeof(pfr_pop_struct));
					ret->Pops[ret->NumPops]->Loci = (pfrfb_locus **)calloc(ret->L,sizeof(pfrfb_locus *));
					
					/* check to make sure the popname is unique compared to previous names, then transfer the pop name */
					if(IdxFromPopName(tempPop,ret)>-1) {
						fprintf(stderr,"Error! While reading population with index %d: PopName %s has been used previously.  PopNames must be unique.  Exiting...\n",ret->NumPops,tempPop);
						exit(1);
					}	
					sprintf(ret->Pops[ret->NumPops]->Name,tempPop);
					
					/* then cycle over loci and get all the info */
					for(i=0;i<ret->L;i++)  {
						/* pick up the correct locus flag */
						sprintf(tempLocName,"Loc%d",i+1);
						GET_STR(tempLocFlag);
						if(strcmp(tempLocFlag,tempLocName) != 0) {
							fprintf(stderr,"Error! While reading data from -p/--pop option for pop-name.  Expecting to chew flag %s, but got %s instead.  Exiting...\n",tempLocName,tempLocFlag);
							exit(1);  
						}
						
						/* allocate necessary memory to the locus */
						ret->Pops[ret->NumPops]->Loci[i] = (pfrfb_locus *)malloc(sizeof(pfrfb_locus));
						ret->Pops[ret->NumPops]->Loci[i]->PA = (double *)calloc(ret->NumA,sizeof(double));
						
						/* then get the trio genotype frequencies */
						for(j=0;j<ret->NumA;j++) {
							ret->Pops[ret->NumPops]->Loci[i]->PA[j] = GET_DUB;
						}
					}
					ret->NumPops++;
				}
			} /* closes if ALREADY_HAS */
		}

		if(REQUIRED_OPTION(
			Mf,
			M,
			max-misses,
			J1 ...,
			max value of different types of states to keep track of,
			This option allows you to specify the maximum value of particular states that you are keeping track of. 
			It must come after the -a/--analysis-type option since that option tells us what the different states are
			and how many there will be.  The different arguments refer to the different states.  For example if we use
			-a mapatri then the program will expect three integer arguments to -m/--max-misses.  Namely J1 is the max
			number of Mendelian incompatibilities between the father and the kid alone.  J2 is the max number of Mendelian 
			incompatibilities between the mother and the kid alone. J3 is the max number of incompatibilities in the 
			trio as a whole.  For matri J1 is between mother and kid and J2 is between trio as a whole.  For patri J1 is 
			between father and kid and J2 is between trio as a whole.  For tri J1 is max number of incompatibilities in the
			trio as a whole. ))
		{
			if(ALREADY_HAS(Af, -a/--analysis-type)) {
				if(ARGS_EQ(ret->Xdim)) { int i;
					ret->Maxes = (int *)calloc(ret->Xdim,sizeof(int));
					for(i=0;i<ret->Xdim;i++)  {
						ret->Maxes[i] = GET_INT;
					}
				}
			}
		}

		if(MULT_USE_OPTION(
			lamoptsf,
			L,
			lambdas,
			S1 S2 R1,
			specify how lambda will be calculated and used,
			This option may be used multiple times---each time it specifies a Lambda calculation to do.
			S1 is the name of the population that should be the numerator of the log-likelihood ratio. S2 is the denominator.
			R1 is the critical value of Lambda.  
			If the computed lambda is greater than R1 then that trio will be printed to standard output.  This option
			should only be issued after all the -p/--pop options have been given. , 
			MAX_NUM_LAMOPTS))
		{
			if(ARGS_EQ(3)) {
				char tempstr[MAX_POP_NAME];
				
				/* allocate memory and do some variable initialization */
				PF.Lambdas[PF.NumLams] = (lam_struct *)malloc(sizeof(lam_struct));
				PF.Lambdas[PF.NumLams]->LC=0.0;
				
				
				/* get logl numerator */
				GET_STR(tempstr);
				if(IdxFromPopName(tempstr,ret)==-1) {
					fprintf(stderr,"Error! Unrecognized population name %s used for argument S1 in -L/--lambdas.  Exiting...\n",tempstr);
					exit(1);
				}
				else {
					PF.Lambdas[PF.NumLams]->LN = IdxFromPopName(tempstr,ret);
				}
				
				/* get logl denomninator */
				GET_STR(tempstr);
				if(IdxFromPopName(tempstr,ret)==-1) {
					fprintf(stderr,"Error! Unrecognized population name %s used for argument S2 in -M/--monte-carlo.  Exiting...\n",tempstr);
					exit(1);
				}
				else {
					PF.Lambdas[PF.NumLams]->LD = IdxFromPopName(tempstr,ret);
				}
				
				/* then get lambda_c and N */
				PF.Lambdas[PF.NumLams]->LC = GET_DUB;
				
				/* increment the counter of the number of MC experiments */
				PF.NumLams++;
			}
			
		}
				
		if( REQUIRED_OPTION(
				mothers_f,
				m,
				mothers,
				S,
				name of file holding genotypes of potential mothers,
				Name of file holding genotypes of potential mothers. The file itself
				must start with the number of individuals in it.  Then for each individual a string identifier is
				given followed by a single string giving the genotype like 00112122001n0021.  0 denotes homozygous for
				the 0 allele.  1 denotes heterozygous and 2 denotes homozygous for the 1 allele.  n denotes missing data.
		) ) { 
			GET_STR(PF.MaFile);
			}
			
			
			if( REQUIRED_OPTION(
				fathers_f,
				f,
				fathers,
				S,
				name of file holding genotypes of potential fathers,
				Name of file holding genotypes of potential fathers. The file itself
				must start with the number of individuals in it.  Then for each individual a string identifier is
				given followed by a single string giving the genotype like 00112122001n0021.  0 denotes homozygous for
				the 0 allele.  1 denotes heterozygous and 2 denotes homozygous for the 1 allele.  n denotes missing data.
		) ) {  
			GET_STR(PF.PaFile);
			}
			
			if( REQUIRED_OPTION(
				kids_f,
				k,
				kids,
				S,
				name of file holding genotypes of potential offspring,
				Name of file holding genotypes of potential offsrpring. The file itself
				must start with the number of individuals in it.  Then for each individual a string identifier is
				given followed by a single string giving the genotype like 00112122001n0021.  0 denotes homozygous for
				the 0 allele.  1 denotes heterozygous and 2 denotes homozygous for the 1 allele.  n denotes missing data.
		) ) {  
			GET_STR(PF.KidFile);
			}
			
					
	END_OPT_LOOP
	
	/* at the end, assign ret to PF.D */
	PF.D = ret;
	
	return(PF);
	
}	


//Initialize a prefetch constant for use with vec_dst(), vec_dstt(), vec_dstst 
//or vec_dststt 
inline unsigned int GetPrefetchConstant( int blockSizeInVectors, 
int blockCount, 
int blockStride ) 
{ 
//ASSERT( blockSizeInVectors > 0 && blockSizeInVectors <= 32 ); 
//ASSERT( blockCount > 0 && blockCount <= 256 ); 
//ASSERT( blockStride > MIN_SHRT && blockStride <= MAX_SHRT ); 
return ((blockSizeInVectors << 24) & 0x1F000000) | 
 ((blockCount << 16) && 0x00FF0000) | 
(blockStride & 0xFFFF); 
} 


/* note that Lam is the index of the lambda spec to use */
void DoMaPaTri_RandoMate(pfr_opts PF, geno_data GD, int Lam)
{
	int m0, m1, m2;  /* the maxes */
	int *moms, Nmoms;  /* an array of indices of compatible mothers */
	int *dads, Ndads;  /* an array of indices of compatible fathers */
	int Ntrios;
	int kid, ma, pa, mism;
	double Lambda, Crit;
	pfrfb_locus **Num, **Denom;
	
	
	/* set the maxes */
	m0 = PF.D->Maxes[0];
	m1 = PF.D->Maxes[1];
	m2 = PF.D->Maxes[2];
	
	/* allocate to moms and dads */
	moms = (int *)calloc(GD.NumMa,sizeof(int));
	dads = (int *)calloc(GD.NumPa,sizeof(int));
	
	/* if Lam>=0 get local addresses of the Num and Denom */
	if(Lam>=0) {
		Num = PF.D->Pops[ PF.Lambdas[Lam]->LN ]->Loci;
		Denom = PF.D->Pops[ PF.Lambdas[Lam]->LD ]->Loci;
		Crit = PF.Lambdas[Lam]->LC;
	}
	
	/* cycle over kids */
	for(kid=0;kid<GD.NumKids;kid++) {
		printf("MAPATRI_KID: %s ",GD.Kids[kid].id);  /* print the kid's name out */
		
		/* get the compatible moms */
		Nmoms=0;
		for(ma=0;ma<GD.NumMa;ma++)  {
			mism = CountMismPair(GD.Moms[ma],GD.Kids[kid]);
			if(mism<=m1) {
				moms[Nmoms++] = ma;
			}
		}
		printf(" CompatMoms= %d ",Nmoms);
		
		/* get the compatible dads */
		Ndads = 0;
		for(pa=0;pa<GD.NumPa;pa++)  {
			mism = CountMismPair(GD.Dads[pa],GD.Kids[kid]);
			if(mism<=m0) {
				dads[Ndads++] = pa;
			}
		}
		printf(" CompatDads= %d \n",Ndads);
		
		
		/* now cycle over all pairs of compatible moms and dads and compute Lambda for the ones we're interested in */
		Ntrios=0;
		for(ma=0;ma<Nmoms;ma++)  {
			for(pa=0;pa<Ndads;pa++)  {
				mism = CountMismTrio(GD.Dads[ dads[pa] ],  GD.Moms[ moms[ma] ],  GD.Kids[kid]);
				
				if(mism <= m2) {
					Ntrios++;
					if(Lam>=0) {
						Lambda = Lambda_MaPaTri(GD.Kids[kid].G, GD.Moms[ moms[ma] ].G, GD.Dads[ dads[pa] ].G, Num, Denom, gNumLoci);
						//printf("LAMBDA: %f\n",Lambda);
						if(Lambda>Crit) {
							printf("MATCHING_TRIO: kid= %s  pa= %s  ma= %s  mism= %d  lambda= %f\n",  
								GD.Kids[kid].id, GD.Moms[ moms[ma] ].id, GD.Dads[ dads[pa] ].id, mism, Lambda );
						}
					}
				}
			}
		}
		printf("MAPATRI_TRIO_COUNT: %s  trios<m= %d\n",GD.Kids[kid].id,Ntrios); 
		
	}
	
	
}



int main (int argc, char *argv[]) 
{
//	int ma,pa,kid;
//	uint_vector *boing;
//	uchar_vector count;
//	vector signed int zero;
//	sint_vector intcounts;
//	unsigned int a,s,g;
	pfr_opts PF;
	geno_data GD;
//	int mism;
//	unsigned int cword;


/*******************************************************
	THIS IS OLD UNUSED STUFF. COMMENTED OUT WITH THE #IFDEF
******************************************************/
#ifdef FALSE	
	//cword = GetPrefetchConstant(10,128,1000);
	
	
	zero = vec_splat_s32(0);
	
	boing = (uint_vector *)calloc(100,sizeof(uint_vector));
	
	for(j=0;j<1;j++) {
		boing[j].v = vec_splat_u32(15);
		
		for(i=0;i<4;i++) {
			printf("%d ",boing[j].elem[i]);
		}	
		printf("\n");
		count.v = bitcount( (vector unsigned char)(boing[j].v) ) ;
		for(i=0;i<16;i++) {
			printf("%d ",count.elem[i]);
		}	
		printf("\n");
		
		/* and here we count up the number of one bits over all parts */
		intcounts.v = vec_sums( (vector signed char) vec_sum4s( count.v, zero ) , zero );
		
		for(i=0;i<4;i++) {
			printf("%d ",intcounts.elem[i]);
		}
		printf("\n");
		
	}
	
	// insert code here...
    printf("Hello, World!\n");
	
	StringTo32Bits("0000111N?222nnnn000011112222nnn?",&a,&s,&g);
	printf("a: ");
	Print32Bits(a);
	printf("\n");
	
	printf("s: ");
	Print32Bits(s);
	printf("\n");

	printf("g: ");
	Print32Bits(g);
	printf("\n");

	Print32Bits(1);
	printf("\n");
#endif
/*************************************************************
END UNUSED STUFF 
********************************************/

	/* Get options and then data */
	PF = GetPFROptions(argc,argv);
	GD = GetGenoData(PF);
	
	
	if(0) {
		pfrfb_locus **Num; pfrfb_locus **Denom; 
		int NumLam; int kid; int ma; int pa; double Lambda;
		Num = PF.D->Pops[ PF.Lambdas[0]->LN ]->Loci;
		Denom = PF.D->Pops[ PF.Lambdas[0]->LD ]->Loci;
		NumLam=0;
		for(kid=0;kid<GD.NumKids;kid++) {
		for(ma=0;ma<GD.NumMa;ma++)  {
		for(pa=0;pa<GD.NumPa;pa++)  {
			Lambda = Lambda_MaPaTri(GD.Kids[kid].G, GD.Moms[ ma ].G, GD.Dads[ pa ].G, Num, Denom, gNumLoci);
			if(NumLam % 1000000 == 0) {
				printf("%d  %f\n",NumLam,Lambda);
			}
			NumLam++;
		}
		}printf("Done with kid %d\n",kid);}
	return(0);
	}
	
	DoMaPaTri_RandoMate(PF, GD, 0);


	return 0;
	
#ifdef FALSE	
	/* for now the default is just going to be a mapatri analysis */
	for(kid=0;kid<GD.NumKids;kid++) {
		
		for(ma=0;ma<GD.NumMa;ma++)  {
		
			for(pa=0;pa<GD.NumPa;pa++)  {
			
				//Prefetch the Dads.  Note, it is not worth it here!
				//vec_dst(&(GD.Dads[pa+1]),cword,0);
				
				mism = CountMismTrio(GD.Dads[pa],GD.Moms[ma],GD.Kids[kid]);
				
				printf("MISM_NUM_OF_PFR : (%s) %s %s : %d\n",GD.Moms[ma].id, GD.Dads[pa].id, GD.Kids[kid].id,mism);
			}
		}
	}
#endif

}
