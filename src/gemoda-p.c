/*
   (c) Massachusetts Institute of Technology, 2003 
 */

/*#include "pdbprot.h"*/
#include "pdb.h"
#include "protAlign.h"
#include "bitSet.h"
#include "convll.h"

#include <unistd.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>


void
usage(char **argv)
{
	fprintf(stdout,
		"Usage: %s -i <comma separated list of processed pdb files>\n\t"
		"-l <word size> -g <threshold> -k <support> [-s <uniq support>]\n\t [-c <cluster method (0|1)>] [-x <output format (0|1)>]\n\n\n"
		"Required flags and input:\n\n"
		"-i <comma-separated list of processed pdb files>:\n\t"
		"A list of files to be searched for motifs.  These should\n\t"
		"be files created from pdb files using the gemoda-pdb-preprocess\n\t"
		"script in the 'scripts' directory.  Each file contains\n\t"
		"one chain from a protein structure; multiple subunits of\n\t"
		"a larger protein may be searched by including them as \n\t"
		"separate files.\n\n"
		"-l <word size>:\n\t"
		"Minimum length of motifs; also the sliding window length\n\t"
		"over which all motifs must meet the similarity criterion\n\n"
		"-k <support>:\n\t"
		"Minimum number of motif occurrences.\n\n"
		"-g <threshhold>:\n\t"
		"Similarity threshold.  Two windows, when scored with the\n\t"
		"similarity function, must have at most this distance threshold\n\t"
		"to be deemed 'connected'.  This criterion must be met over\n\t"
		"all sliding windows of length L.  For a RMSD similarity\n\t"
		"metric, this is the maximum RMSD that two windows can have\n\t"
		"and still be 'connected'.\n\n\n"
		"Optional flags and input:\n\n"
		"-x <output format (0|1)>:\n\t"
		"Allows the user to choose what format the output file\n\t"
		"should be in.  Default: 1.\n\n\t"
		"0: PDB format\n\t\t"
		"Traditional format for all PDB files, including remarks\n\t\t"
		"describing what is contained in the files.\n\t\t"
		"(currently being implemented)\n\n\t"
		"1: XML format\n\t\t"
		"An xml format suitable for use with the gemoda-p pattern\n\t\t"
		"viewer (GPV).  This is a simple python-based viewer that\n\t\t"
		"allows the user to view any given motif one at a time\n\t\t"
		"with all occurrences of that motif superimposed and\n\t\t"
		"represented by stick structures of different colors.\n\n"
		"-s <unique support>:\n\t"
		"A pruning parameter that requires the motif to occur in\n\t"
		"at least <unique support> different input sequences.  Note\n\t"
		"that this parameter must be less than or equal to the total\n\t"
		"support parameter set by the -k flag.\n\n"
		"-c <cluster method [0|1]>:\n\t"
		"The clustering method to be used after evaluating the\n\t"
		"similarity of the unique windows in the input.  Note that\n\t"
		"the clustering method will have a significant impact on\n\t"
		"both the results that are obtained and the computation\n\t"
		"time.\n\n\t"
		"0: clique-finding\n\t\t"
		"The default clustering method, uses established\n\t\t"
		"algorithms to find all maximal cliques in the data.\n\t\t"
		"This will give the most thorough results (that are\n\t\t"
		"provably exhaustive), but will also give less-significant\n\t\t"
		"results in addition to the more interesting and\n\t\t"
		"signicant motifs.  The results are deterministic but\n\t\t"
		"may take some time on data sets with high similarity or\n\t\t"
		"if the similarity threshold is set extremely low.\n\t"
		"1: single-linkage clustering\n\t\t"
		"Uses a single-linkage-type clustering where all nodes\n\t\t"
		"that are connected are put in the same cluster.  This\n\t\t"
		"method is also deterministic and will be much faster than\n\t\t"
		"clique-finding, but it loses guarantees of exhaustiveness\n\t\t"
		"in searching the data set.\n\n",
		argv[0]);
	fprintf(stdout, "\n");
}

/*
   Functions defined in newConv.c 
 */
extern cll_t * convolve(bitGraph_t *bg, int support, int R, int *indexToSeq, int P, int clusterMethod, int **offsetToIndex, int numberOfSequences,
		int noConvolve);




//  parsing comma separated list
int
countFileNames(char *s){
	int i;
	int numFiles=1;
	for(i=0 ; i<strlen(s) ;i++){
		if(s[i]==','){
			numFiles++;
		}
	}
	return numFiles;
}

//  parsing comma separated list
char **
parseFileNames(char *s,int numFiles){
	int i,j,k;
	int startLength;
	char **files = NULL;

	files = (char **)malloc( numFiles * sizeof(char *));
	if(files == NULL){
		fprintf(stderr, "can't allocate file names!");
		exit(0);
	}
	j=0;
	k=0;
	startLength=strlen(s);
	for(i=0 ; i<startLength ;i++){
		if(s[i]==','){
			files[j] = &s[k];
			j++;
			s[i] = '\0';
			k=i+1;
		}
	}
	files[j] = &s[k];
	return files;
}

int
main(int argc, char **argv)
{
	int inputOption = 0;
	char *outputFile = NULL;
	FILE *OUTPUT_FILE = NULL;
	int L = -1;
	int status = 0;
	double g = 0.0;
	int sup = -1;
	char **inputFiles = NULL;
	int P = 0, R = 0;

	int i;
	bitGraph_t *bg = NULL;
	//int *prevRowArray;

	int numFiles=0;
	Pdbentry_ **pdbs = NULL;
	caPoint_t *cap = NULL;
	int numWindows = 0;
	int numChains = 0;
	cll_t *allCliques = NULL;
	int us = -1;
	int clusterMethod = 0;
	int noConvolve = 0;

	// outputFormat: 0=pdb, 1=xml
	int outputFormat = 1;



	/*
	   Get use command-line options 
	 */
	while ((inputOption = getopt(argc, argv, "i:o:l:g:k:s:c:xn")) != EOF) {
		switch (inputOption) {
		case 'i':
			numFiles = countFileNames(optarg);
			inputFiles = parseFileNames(optarg, numFiles);
			break;
		case 'o':
			outputFile = optarg;
			break;
		case 'l':
			L = atoi(optarg);
			break;
		case 'g':
			g = atof(optarg);
			status++;
			break;
		case 'k':
			sup = atoi(optarg);
			break;
		case 's':
			us = atoi(optarg);
			break;
		case 'x': 
			outputFormat = 1;
			break;
		case 'c':
			clusterMethod = atoi(optarg);
			break;
		case 'n':
			noConvolve = 1;
			break;
		case '?':
			fprintf(stderr, "Unknown option `-%c'.\n", optopt);
			usage(argv);
			return EXIT_SUCCESS;
		default:
			usage(argv);
			return EXIT_SUCCESS;
		}
	}
	if (inputFiles == NULL || L == -1 || status < 1 || sup == -1) {
		usage(argv);
		return EXIT_SUCCESS;
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

	//  read in each pdb file
	pdbs = (Pdbentry_ **)malloc(numFiles * sizeof(Pdbentry_ *));
	if (pdbs == NULL) {
		fprintf(stderr,"Memory error - gemoda-p - %s\n",strerror(errno));
		fflush(stderr);
		exit(0);
	}
	
	for( i=0 ; i<numFiles ; i++){
		//  get alpha carbons ONLY!
		pdbs[i] = get_pdb(inputFiles[i]); 
	}


	//  get pointers to translate from the triplet
	//  (pdb file, chain number, residue) -> single array index
	cap = getCaPoint(pdbs, numFiles, L, &numWindows, &numChains);
	printf("there are %d windows!!\n", numWindows);
	printf("there are %d chains!!\n", numChains);

	bg = alignPdbs(cap, pdbs, L, g, numWindows);

	//  now we're convolving on offsets
	/*printBitGraph(oam);*/
	
	
	allCliques = convolve(bg, sup, R, NULL, P, clusterMethod,NULL,numChains,noConvolve);
	/*printCll(allCliques);*/
	outputProtPats(cap, pdbs, allCliques, L, us, OUTPUT_FILE, outputFormat);

	allCliques = popAllCll(allCliques);

	/*deleteBitGraph(bg);*/
	/*deleteBitGraph(oam);*/

	for( i=0 ; i<numFiles ; i++){
		free_pdb(pdbs[i]);
	}
	free(pdbs);
	free(inputFiles);
	free(cap);
	deleteBitGraph(bg);


	fclose(OUTPUT_FILE);
	exit(0);
}
