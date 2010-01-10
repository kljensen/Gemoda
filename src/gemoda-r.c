/*
   (c) Massachusetts Institute of Technology, 2003 
 */

#include "bitSet.h"
#include "convll.h"
#include "FastaSeqIO/fastaSeqIO.h"
#include <unistd.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
// #include <time.h>
/*#include "efence.h"*/
#include "realIo.h"
#include "realCompare.h"

/**************************************************************
 * Displays hopefully useful output for the user
 * **************************************************************/

void
usage(char **argv)
{
	fprintf(stdout,
		"Usage: %s -i <Fasta sequence file> "
		"-l <word size> \n\t-k <support> -g <threshold> -m <matrix name> [-z] \n\t[-c <cluster method [0|1]>] [-p <unique support>] \n\n\n"
		"Required flags and input:\n\n"
		"-i <Fasta sequence file>:\n\t"
		"File containing all sequences to be searched, in Fasta format.\n\n"
		"-l <word size>:\n\t"
		"Minimum length of motifs; also the sliding window length\n\t"
		"over which all motifs must meet the similarity criterion\n\n"
		"-k <support>:\n\t"
		"Minimum number of motif occurrences.\n\n"
		"-g <threshold>:\n\t"
		"Similarity threshold.  Two windows, when scored with the\n\t"
		" similarity matrix defined by the -m flag, must have at least\n\t"
		" this score in order to be deemed 'connected'.  This criterion\n\t"
		" must be met over all sliding windows of length l.\n\n"
		"-c <cluster method [0|1]>:\n\t"
		"The clustering method to be used after evaluating the "
		"\n\tsimilarity of the unique words in the input.  Note that the "
		"\n\tclustering method will have a significant impact on both the "
		"\n\tresults that one obtains and the computation time.\n\n\t"
		"0: clique-finding\n\t\t"
		"Uses established methods to find all maximal cliques in the "
		"\n\t\tdata.  This will give the most thorough results (that are "
		"\n\t\tprovably exhaustive), but will also give less-significant "
		"\n\t\tresults in addition to the most interesting and most\n\t"
		"significant ones.  The results are deterministic but may take some "
		"\n\t\ttime on data sets with high similarity or if the similarity "
		"\n\t\tthreshold is set extremely low.\n\t"
		"1: single-linkage clustering\n\t\t"
		"Uses a single-linkage-type clustering where all nodes that "
		"\n\t\tare connected are put in the same cluster.  This method is "
		"\n\t\talso deterministic and will be faster than clique-finding, "
		"\n\t\tbut it loses guarantees of exhaustiveness in searching the "
		"\n\t\tdata set.\n\n",
		
		"-p <unique support>:\n\t"
		"A pruning parameter that requires the motif to occur in "
		"\n\tat least <unique support> different input sequences.  Note "
		"\n\tthat this parameter must be less than or equal to the total "
		"\n\tsupport parameter set by the -k flag.\n\n",
		
			argv[0]);
	fprintf(stdout, "\n");
}

/*
   Functions defined in newConv.c 
 */
extern cll_t * convolve(bitGraph_t *bg, int support, int R, int *indexToSeq, int p, int clusterMethod, int **offsetToIndex, int numberOfSequences, 
		int noConvolve, FILE *OUTPUT_FILE);
extern bitGraph_t * pruneBitGraph(bitGraph_t *bg, int *indexToSeq, 
		int **offsetToIndex, int numOfSeqs, int p);

int
countExtraParams(char *s){
	int i = 0;
	int numParams = 1;
	for (i = 0; i < strlen(s); i++) {
		if(s[i]==','){
			numParams++;
		}
	}
	return numParams;
}
// Borrowed from the gemoda-p code, there it used to parse filenames, here
//  we are parsing comma-separated lists of doubles
double*
parseExtraParams(char *s,int numParams){
	int i = 0, j = 0, k = 0;
	int startLength = 0;
	double *extraParams = NULL;
	char *paramString = NULL;
	
	extraParams = (double *) malloc(numParams * sizeof(double));
	if(extraParams == NULL) {
		fprintf(stderr, "Can't allocate extra params!\n");
		exit(0);
	}
	j = 0;
	k = 0;
	startLength = strlen(s);
	for(i = 0; i < startLength; i++) {
		if(s[i] == ',') {
			// We've found an end.  So point the pointer to
			//  the beginning of the previous string.
			paramString = &s[k];
			// Terminate the string where the comma used to be
			s[i] = '\0';
			// Update the location for the next string beginning
			k = i+1;
			// Convert to a double and update the param number.
			extraParams[j] = atof(paramString);
			j++;
		}
	}
	// Don't forget to do the last one, which isn't comma-terminated.
	paramString = &s[k];
	extraParams[j] = atof(paramString);
	return(extraParams);
}

int
main(int argc, char **argv)
{
	int inputOption = 0;
	char *sequenceFile = NULL;
	FILE *SEQUENCE_FILE = NULL;
	char *outputFile = NULL;
	FILE *OUTPUT_FILE = NULL;
	int L = 0;
	int status = 0;
	double g = 0;
	int sup = 2;
	int R = 1;
	int P = 0;
	int compFunc = 0;
	double *extraParams = NULL;
	int numExtraParams = 0;
	int i = 0, j = 0;
	/*int j, k, i, l;*/
	int noConvolve = 0;
	int samp = 1;
	int supportDim = 0, lengthDim = 0;
	bitGraph_t *oam = NULL;
	unsigned int **d = NULL;
	int oamSize = 0;

	cll_t *allCliques = NULL;
	/*cll_t *curCliq = NULL;*/
	/*int curSeq;*/
	/*int curPos;*/
	int clusterMethod = 0;
	int joelOutput = 0;

	// gemoda-r new stuff
	rdh_t *data = NULL;

	/*
	   Get command-line options 
	 */
	while ((inputOption = getopt(argc, argv, "p:m:e:i:o:l:g:k:c:njs:")) != EOF) {
		switch (inputOption) {
		// Comparison metric
		case 'm':
			compFunc = atoi(optarg);
			break;
		// Input file
		case 'i':
			sequenceFile = optarg;
			break;
		// Output file
		case 'o':
			outputFile = (char*) malloc((strlen(optarg)+1)*sizeof(char));
			if (outputFile == NULL) {
				fprintf(stderr,"Error allocating memory for options.\n");
				exit(EXIT_FAILURE);
			} else {
				strcpy(outputFile,optarg);
			}
			break;
		// Minimum motif length
		case 'l':
			L = atoi(optarg);
			break;
		// Minimum motif similarity score
		case 'g':
			g = atof(optarg);
			status++;
			break;
		// Minimum support (number of motif occurrences)
		case 'k':
			sup = atoi(optarg);
			break;

/***************************************************************
 * Recursive initial pruning: an option for clique finding.
 *   It takes all nodes with less than the minimum
 *   number of support and removes all of their nodes, and does this 
 *   recursively so that nodes that are connected to many sparsely connected
 *   nodes will be removed and not left in the 
 * This option is deprecated as it is at worst no-gain and at best useful.
 *   It will be on by default for clique-finding, but can be turned 
 *   back off with some
 *   minor tweaking.  For almost all cases in which it does not speed
 *   up computations, it will have a trivial time to perform.  Thus, if 
 *   clique-finding is turned on, then R is set to 1 by default.
		case 'r':
			R = 1;
			break;
************************************************************************/
		// Optional pruning parameter to require at motif occurrences
		//   in at least P distinct input sequences
		
		case 'p':
			P = atoi(optarg);
			break;
			
		// Clustering method.
		case 'c':
			clusterMethod = atoi(optarg);
			break;
		// Extra parameters for comparison function
		case 'e':
			numExtraParams = countExtraParams(optarg);
			extraParams = parseExtraParams(optarg,numExtraParams);
			break;
		case 'n':
			noConvolve = 1;
			break;
		case 'j':
			joelOutput = 1;
			break;
		case 's':
			samp = atoi(optarg);
			break;
		// Catch-all.
		case '?':
			fprintf(stderr, "Unknown option `-%c'.\n", optopt);
			usage(argv);
			return EXIT_SUCCESS;
		default:
			usage(argv);
			return EXIT_SUCCESS;
		}
	}
	// Require an input file, a nonzero length, and a similarity threshold
	//  to be set.
	if (sequenceFile == NULL || L == 0 || status < 1) {
		usage(argv);
		return EXIT_SUCCESS;
	}
	// Open the sequence file
	if ((SEQUENCE_FILE = fopen(sequenceFile, "r")) == NULL) {
		fprintf(stderr, "Couldn't open file %s; %s\n", sequenceFile, strerror(errno));
		exit(EXIT_FAILURE);
	}
	// Open the output file
	if (outputFile != NULL) {
		if ((OUTPUT_FILE = fopen(outputFile, "w")) == NULL) {
			fprintf(stderr, "Couldn't open file %s; %s\n", outputFile, strerror(errno));
			exit(EXIT_FAILURE);
		}
	} else {
		OUTPUT_FILE = stdout;
	}

	

	// Verbosity in output helps to distinguish output files.
	fprintf(OUTPUT_FILE,"Input file = %s\n", sequenceFile);
	fprintf(OUTPUT_FILE,"l = %d, k = %d, g = %f\n", L, sup, g);
	if (P > 1) {
		fprintf(OUTPUT_FILE,"Minimum # of sequences with motif = %d\n", P);
	}
	if (R > 0) {
		fprintf(OUTPUT_FILE,"Recursive pruning is ON.\n");
	}

	data = readRealData(SEQUENCE_FILE);
	fclose(SEQUENCE_FILE);
//	printf("size = %d,indexSize = %d\n",data->size,data->indexSize);
//	printf("size1 = %d,size2 = %d\n",data->seq[0]->size1,data->seq[0]->size2);
//	for(i = 0; i < 2; i++) {
//		for(j = 0; j < 3; j++) {
//	printf("%lf,%lf,%lf\n",gsl_matrix_get(data->seq[i],j,0),
//			gsl_matrix_get(data->seq[i],j,1),
//			gsl_matrix_get(data->seq[i],j,2));}}
	oam = realComparison(data, L, g,compFunc,extraParams);
//	printf("oam->size = %d\n", oam->size);
	if ((samp > 0) && (clusterMethod == 0)) {
		// We are currently using one gap per sequence, as done in 
		//  realCompare.c's call to initRdhIndex in realComparison.
		//  Note that this is data->size, NOT oam->size.
		d = getStatMat(oam,sup,L,&supportDim,&lengthDim,data->size,samp,OUTPUT_FILE);
	} else {
		d = NULL;
		supportDim = 0;
	}

	allCliques = convolve(oam, sup, R, data->indexToSeq, P, clusterMethod,data->offsetToIndex, data->size,noConvolve,OUTPUT_FILE);

	oamSize = oam->size;
	// Do some early memory cleanup since this is so big.
	deleteBitGraph(oam);

	if ((samp > 0) && (clusterMethod == 0)) {
		cumDMatrix(d,allCliques,supportDim,lengthDim,oamSize,data->size);
		calcStatAllCliqs(d,allCliques,oamSize - data->size);
		allCliques = sortByStats(allCliques);
	}

	if (joelOutput == 0) {
		outputRealPats(data, allCliques, L, OUTPUT_FILE,d);
	}
	else {
	outputRealPatsWCentroid(data,allCliques,L,OUTPUT_FILE,extraParams,compFunc);
	}

	freeD(d,supportDim);
	freeRdh(data);
	allCliques = popAllCll(allCliques);
	fclose(OUTPUT_FILE);

	return 0;
}
