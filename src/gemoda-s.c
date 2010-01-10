/*
   (c) Massachusetts Institute of Technology, 2003 
 */

#include "bitSet.h"
#include "spat.h"
#include "convll.h"
#include "matdata.h"
#include "FastaSeqIO/fastaSeqIO.h"
//#include "probability.h"
#include <unistd.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include "patStats.h"
// #include <time.h>

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
		"-m <matrix name>:\n\t"
		"Name of the similarity matrix to be used to compare windows.\n\t"
		"Use -z to see a list of matrices installed by default.\n\n\n"
		"Optional flags and input:\n\n"
		"-z:\n\t"
		"Lists all of the similarity matrices available with the\n\t"
		"initial installation of Gemoda.  Note that this overrides\n\t"
		"all other options and will only give this output.\n\n"
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
		"\n\t\tdata set.\n\n"
		"-p <unique support>:\n\t"
		"A pruning parameter that requires the motif to occur in "
		"\n\tat least <unique support> different input sequences.  Note "
		"\n\tthat this parameter must be less than or equal to the total "
		"\n\tsupport parameter set by the -k flag.\n\n",
			argv[0]);
	fprintf(stdout, "\n");
}

void
matrixlist(void) 
{
	fprintf(stdout, "\nThe following similarity matrices are installed "
			"with the default Gemoda installation.\n  Most of these "
			"were obtained from publically available BLAST distributions. \n\n"
			"dna_idmat:\n\t"
			"Identity matrix for DNA: returns 1 when A,C,G,T are "
			"compared to \n\tthemselves, 0 otherwise.\n\n"
			"identity_aa:\n\t"
			"Identity matrix for amino acids: returns 1 when any \n\t"
			"letter but J,O,U are compared to themselves, and 0 "
			"otherwise.\n\n"
			"idmat:\n\t"
			"Similar to identity_aa, but it returns 10 in place "
			"of 1.\n\n"
			"est_idmat:\n\t"
			"Similar to idmat, but it returns -10 in place of 0. "
			"\n\n"
			"pam100:\n"
			"pam110:\n"
			"pam120:\n"
			"pam130:\n"
			"pam140:\n"
			"pam150:\n"
			"pam160:\n"
			"pam190:\n"
			"pam200:\n"
			"pam210:\n"
			"pam220:\n"
			"pam230:\n"
			"pam240:\n"
			"pam250:\n"
			"pam260:\n"
			"pam280:\n"
			"pam290:\n"
			"pam300:\n"
			"pam310:\n"
			"pam320:\n"
			"pam330:\n"
			"pam340:\n"
			"pam360:\n"
			"pam370:\n"
			"pam380:\n"
			"pam390:\n"
			"pam400:\n"
			"pam430:\n"
			"pam440:\n"
			"pam450:\n"
			"pam460:\n"
			"pam490:\n"
			"pam500:\n\t"
			"PAM matrices for various evolutionary distances.\n\n"
			"blosum30:\n"
			"blosum35:\n"
			"blosum40:\n"
			"blosum45:\n"
			"blosum50:\n"
			"blosum55:\n"
			"blosum60:\n"
			"blosum62:\n"
			"blosum65:\n"
			"blosum70:\n"
			"blosum75:\n"
			"blosum80:\n"
			"blosum85:\n"
			"blosum90:\n"
			"blosum100:\n\t"
			"BLOSUM matrices for various evolutionary distances.\n\n"
			"blosumn:\n\t"
			"BLOSUM matrix of unknown origin.\n\n"
			"dayhoff:\n\t"
			"'Vanilla-flavored' pam250, very similar to pam250.\n\n"
			"phat_t75_b73:\n"
			"phat_t80_b78:\n"
			"phat_t85_b82:\n\t"
			"BLOSUM-clustered scoring matrix with target frequency\n\t"
			"PHDhtm clustering = {75,80,85}percent and background frequency\n\t"
			"Persson-Argos clustering = {73,78,82}percent.\n\t"
			"From Ng, Henikoff, & Henikoff, Bioinformatics 16: 760.\n\n"
			"coil_mat:\n"
			"alpha_mat:\n"
			"beta_mat:\n\t"
			"Three structure-specific matrices described by Luthy,\n\t"
			"McLachlan, and Eisenberg in Proteins 10, 229-239, obtained from AAindex.\n\n");
	fprintf(stdout,"\n");
}

/*
 * Functions defined in matrices.c
 */
extern void getMatrixByName(char name[], int mat[][MATRIX_SIZE]);

/*
   Functions defined in align.c 
 */
extern bitGraph_t *alignWordsMat_bit(sPat_t * words, int wc, int mat[][MATRIX_SIZE], int threshold);

/*
   Functions defined in words.c 
 */
extern sPat_t *countWords2(fSeq_t * seq, int numSeq, int L, int *numWords);

/*
   Functions defined in newConv.c 
 */
extern cll_t * convolve(bitGraph_t *bg, int support, int R, int *indexToSeq, int p, int clusterMethod, int **offsetToIndex, int numberOfSequences,
		int noConvolve,FILE *OUTPUT_FILE);
extern bitGraph_t * pruneBitGraph(bitGraph_t *bg, int *indexToSeq, 
		int **offsetToIndex, int numOfSeqs, int p);


int
main(int argc, char **argv)
{
	int inputOption = 0;
	char *sequenceFile = NULL;
	char *outputFile = NULL;
	char *matName = NULL;
	FILE *SEQUENCE_FILE = NULL;
	FILE *OUTPUT_FILE = NULL;
	int L = 0;
	int numberOfSequences = 0;
	fSeq_t *mySequences = NULL;
	fSeq_t *(*seqReadFunct) () = &ReadFSeqs;
	sPat_t *words = NULL;
	int wc;
	int status = 0;
	int g = 0;
	int sup = 2;
	int R = 1;
	int P = 0;
	int (*mat)[MATRIX_SIZE] = NULL;
	int noConvolve = 0;

	int j, k, i, l;
	bitGraph_t *bg = NULL;
	bitGraph_t *oam = NULL;

	// new
	int **offsetToIndex = NULL;
	int *indexToSeq = NULL;
	int *indexToPos = NULL;
	int numberOfOffsets = 0;
	int pos1, pos2;
	//int *prevRowArray;
	sOffset_t *offset1, *offset2;
	cll_t *allCliques = NULL;
	cll_t *curCliq = NULL;
	int curSeq;
	int curPos;
	int clusterMethod = 0;

	// patStats
	int samp = 1;
	unsigned int **d = NULL;
	int supportDim = 0, lengthDim = 0;
	int oamSize = 0;
/*
	// for probability calculations
	int** augmat = NULL;
	int augcount = 0;
	int matmaxscore = 0;
//	probWorkspace_t* work = NULL;
	double* letterfreqs = NULL;
	double** scorefreqs = NULL;
*/
	/*
	   Get command-line options 
	 */
	while ((inputOption = getopt(argc, argv, "i:o:l:g:k:m:p:zc:ns:")) != EOF) {
		switch (inputOption) {
		// Input file
		case 'i':
			sequenceFile = optarg;
			seqReadFunct = &ReadFSeqs;
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
			g = atoi(optarg);
			status++;
			break;
		// Minimum support (number of motif occurrences)
		case 'k':
			sup = atoi(optarg);
			break;
		// Similarity matrix used to find similarity score
		case 'm':
			getMatrixByName(optarg, &mat);
			matName = (char*) malloc(strlen(optarg)*sizeof(char));
			if (matName == NULL) {
				fprintf(stderr,"Error allocating memory for options.\n");
				exit(EXIT_FAILURE);
			} else {
				strcpy(matName,optarg);
			}
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
		case 'n':
			noConvolve = 1;
			break;
		case 's':
			samp = atoi(optarg);
			break;
		// Catch-all.
		case '?':
			fprintf(stderr, "Unknown option `-%c'.\n", optopt);
			usage(argv);
			return EXIT_SUCCESS;
		case 'z':
			matrixlist();
			return EXIT_SUCCESS;
		default:
			usage(argv);
			return EXIT_SUCCESS;
		}
	}
	// Require a similarity matrix
	if (mat == NULL) {
		usage(argv);
		return EXIT_SUCCESS;
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
	// Allocate some sequences
	mySequences = seqReadFunct(SEQUENCE_FILE, &numberOfSequences);
	if (mySequences == NULL) {
		fprintf(stderr, "\nError reading your sequences/text.");
		fprintf(stderr, "\nCheck the format/size of the file.");
		fprintf(stderr, "\nERROR:  %s\n", strerror(errno));
		return EXIT_FAILURE;
	}
	// Close the input files
	fclose(SEQUENCE_FILE);

	// Verbosity in output helps to distinguish output files.
	fprintf(OUTPUT_FILE,"\nMatrix used = %s\n", matName);
	fprintf(OUTPUT_FILE,"Input file = %s\n", sequenceFile);
	fprintf(OUTPUT_FILE,"l = %d, k = %d, g = %d\n", L, sup, g);
	if (P > 1) {
		fprintf(OUTPUT_FILE,"Minimum # of sequences with motif = %d\n", P);
	}
	if (R > 0) {
		fprintf(OUTPUT_FILE,"Recursive pruning is ON.\n");
	}

	// Find the unique words in the input.
	words = countWords2(mySequences, numberOfSequences, L, &wc);

	/*fprintf(stderr, "Counted %d words\n", wc);*/
	/*fflush(stderr);*/
	
	// Align the words that we just found by applying the similarity
	//  matrix to each pair of them.  Note that
	//  bg is the adjacency matrix of words, but we
	//  need an adjacency matrix of offsets instead.  
	bg = alignWordsMat_bit(words, wc, mat, g);

	fprintf(OUTPUT_FILE,"\nAligned!  Creating offset matrix...\n");
	fflush(NULL);

	
	// Create an intermediate translation matrix
	//  to store the offset number of each sequence number/position.
	//
	//  Note that this matrix is better called "Index to offset", and
	//   the other matrices are better called "offset to Seq" and
	//   "offset to Pos"
	offsetToIndex = (int **) malloc (numberOfSequences * sizeof(int *));
	if (offsetToIndex == NULL) {
		fprintf(stderr,"Unable to allocate memory - offsetToIndex in gemoda.c\n%s\n",
				strerror(errno));
		fflush(stderr);
		exit(0);
	}
	for ( i=0 ; i<numberOfSequences ; i++){
		// MPS 5/23/05: Added in "-L+2" to make there only be one
		//  blank between sequences.
		offsetToIndex[i] = malloc((strlen(mySequences[i].seq)-L+2) * sizeof(int));
		if (offsetToIndex[i] == NULL) {
			fprintf(stderr,"Unable to allocate memory - offsetToIndex[%d] in gemoda.c\n%s\n",
					i,strerror(errno));
			fflush(stderr);
			exit(0);
		}
		// MPS 5/23/05: Added in "-L+2" to make there only be one
		//  blank between sequences.
		for ( j=0 ; j<(strlen(mySequences[i].seq)-L+2) ; j++){
			offsetToIndex[i][j] = numberOfOffsets;
			numberOfOffsets++;
		}
	}
	// Now create translation matrices such that we can get the sequence
	//  or position number of a given offset.
	indexToSeq = (int *) malloc (numberOfOffsets * sizeof(int));
	if (indexToSeq == NULL) {
		fprintf(stderr,"Unable to allocate memory - indexToSeq in gemoda.c\n%s\n",
				strerror(errno));
		fflush(stderr);
		exit(0);
	}
	indexToPos = (int *) malloc (numberOfOffsets * sizeof(int));
	if (indexToPos == NULL) {
		fprintf(stderr,"Unable to allocate memory - indexToPos in gemoda.c\n%s\n",
				strerror(errno));
		fflush(stderr);
		exit(0);
	}
	k=0;
	for ( i=0 ; i<numberOfSequences ; i++){
		// MPS 5/23/05: Added in "-L+2" to make there only be one
		//  blank between sequences.
		for ( j=0 ; j<(strlen(mySequences[i].seq)-L+2) ; j++){
			indexToSeq[k] = i;
			indexToPos[k] = j;
			k++;
		}
	}
	//  Now make an offset adjacency matrix!  
	//
	oam = newBitGraph(numberOfOffsets);
	// Go through each unique word
	for (i = 0; i < wc; i++) {
		offset1 = words[i].offset;
		// Go through each occurrence
		for (k = 0; k < words[i].support; k++) {
			// Use the offsetToIndex translation to get the offset
			//  of the first occurrence
			pos1 = offsetToIndex[ offset1[k].seq ][ offset1[k].pos ];
			// And go through each word in the first offset to 
			//  find words that meet the similarity threshold
			for (j = 0; j < wc; j++) {
				if (bitGraphCheckBit(bg,i,j)) {
					offset2 = words[j].offset;
					// And find all of their occurrences,
					// using offsetToIndex to get the
					// offsets, and then setting those
					// locations in the offset adjacency
					// matrix true.
					for (l = 0; l < words[j].support; l++) {
						pos2 = offsetToIndex[ offset2[l].seq ][ offset2[l].pos ];
						bitGraphSetTrueSym(oam, pos1, pos2);
					}
				}
			}
		}
	}

	fprintf(OUTPUT_FILE,"Offset matrix created...");
	deleteBitGraph(bg);

	if ((samp > 0) && (clusterMethod == 0)) {
		fprintf(OUTPUT_FILE," taking preliminary statistics.\n");
		fflush(NULL);
		d = getStatMat(oam,sup,L,&supportDim,&lengthDim,numberOfSequences,samp,OUTPUT_FILE);
		fprintf(OUTPUT_FILE,"Now filtering...\n");
		fflush(NULL);
	} else {
		fprintf(OUTPUT_FILE," now filtering.\n");
		fflush(NULL);
		d = NULL;
		supportDim = 0;
	}
	//  Now we're convolving on offsets
	allCliques = convolve(oam, sup,R,indexToSeq,P, clusterMethod,offsetToIndex,numberOfSequences,noConvolve,OUTPUT_FILE);

	// Do some early memory cleanup to limit usage
	oamSize = oam->size;
	deleteBitGraph(oam);

	fprintf(OUTPUT_FILE,"Convolved!  Now making output...\n");
	fflush(NULL);

	if ((samp > 0) && (clusterMethod == 0)) {
		cumDMatrix(d,allCliques,supportDim,lengthDim,oamSize,numberOfSequences);
		calcStatAllCliqs(d,allCliques,numberOfOffsets-numberOfSequences);
		allCliques = sortByStats(allCliques);
	}

	//  walk over the cliques and give some output in the format:
	//   pattern <pattern id num>: len=<motif length> sup=<motif instances>
	//    <sequence num> <position num> <motif instance>
	//    ...
	curCliq = allCliques;

// NOTE that down from here is code for Mark's previous probability attempt... it's not being used,
//  but is being put in for the commit since there are other things to commit, but is being commented
//  out so that it won't have any impact on the code.
/*	// This code creates data structures used inside the loop so that
	//  they don't have to be created and destroyed each time a p-value
	//  calculation is made.
	augmat = create_augmented_matrix(mat);
	augcount = augment_matrix(augmat);
	matmaxscore = find_max_score(augmat);

	// Note that maxsize<x> is the maximum SIZE of the array, so that the
	//  maximum array index is the maximum score
	work = makeProbWorkspace(matmaxscore*L + 1,matmaxscore+1,matmaxscore*L+1);
	letterfreqs = (double *) malloc(MATRIX_SIZE * sizeof(double));

	scorefreqs = (double **) malloc(L * sizeof(double*));
	for (i = 0; i < L; i++) {
		scorefreqs[i] = (double *) malloc((matmaxscore+1) * sizeof(double));
	}

	for (i = 0; i < MATRIX_SIZE; i++) {
		letterfreqs[0] = 0;
	}

	letterfreqs[0] = .25;
	letterfreqs[4] = .25;
	letterfreqs[7] = .25;
	letterfreqs[16] = .25;
	printf("First prob = %le\n\n", prob_calc_motif(augmat,matmaxscore,
				curCliq,L,letterfreqs,work,
				g,indexToSeq,indexToPos,mySequences,
				scorefreqs));
*/
// NOTE end of probability code, except for a bit down a little ways.
	i = 0;
	while(curCliq != NULL){
		fprintf(OUTPUT_FILE, "pattern %d:\tlen=%d\tsup=%d", i, curCliq->length+L, curCliq->set->size);
		if (d != NULL) {
			fprintf(OUTPUT_FILE, "\tsignif=%le\n",curCliq->stat);
		} else {
			fprintf(OUTPUT_FILE,"\n");
		}
		// INSERT code to calculate and output motif probability
		//  <--------- HERE
		/*
		polya = totaledge_prob(augmat,maxscore,basestring,stringlength,L,
				letterfreqs,polya,maxsizea,polyb,maxsizeb,
				polyc,maxsizec);
		for(i = 0; i < sizea; i++) {
			j = (sizec - 1) - 1;
			polya[j] = polya[j] + polya[j+1];
		}
		pval = polya[g + augcount*(curCliq->length + L)];
*/
		for(j=0 ; j<curCliq->set->size ; j++){
			pos1 = curCliq->set->members[j];
			curSeq = indexToSeq[pos1];
			curPos = indexToPos[pos1];
			fprintf(OUTPUT_FILE, "   %d\t%d\t", curSeq, curPos);
			for(k=curPos ; k<curPos+curCliq->length+L ; k++){
				fprintf(OUTPUT_FILE, "%c", mySequences[curSeq].seq[k]);
			}
			fprintf(OUTPUT_FILE, "\n");
		}
		fprintf(OUTPUT_FILE, "\n\n");
		curCliq = curCliq->next;
		i++;
	}

	// And do some memory cleanup
// And cleanup of probability stuff...
/*
	free(letterfreqs);
	delete_augmented_matrix(augmat);
*/

	allCliques = popAllCll(allCliques);
	free(indexToSeq);
	indexToSeq = NULL;
	free(indexToPos);
	indexToPos = NULL;

	for ( i=0 ; i<numberOfSequences ; i++){
		free(offsetToIndex[i] );
		offsetToIndex[i] = NULL;
	}
//Free'ing added by MPS, 6/4
	for (i=0 ; i<wc ; i++){
		free(words[i].offset);
	}
	free(words);
// End free'ing added by MPS
	
	free(offsetToIndex);
	offsetToIndex = NULL;
	// -------------------------------------------

	// Free up fastaSequences
	FreeFSeqs(mySequences, numberOfSequences);
	fclose(OUTPUT_FILE);
	return 0;
}
