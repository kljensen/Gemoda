#include "realIo.h"
#include "realCompare.h"
#include "patStats.h"

/*#include "efence.h"*/

//  Turn the substring of s starting at
//  char s[begin] and ending at s[end]
//  int a double.
//  INPUT:
//  a string s
//  integer begin
//  integer end
//  OUTPUT
//  a double
//  NOTE:
//  Throws an error and dies if there's a problem
//  making the double from the substring.  No room
//  for ill-formated data files.
double
wordToDouble(char *s, int begin, int end){
	char *str = NULL;
	char *endptr;
	double val;
	int size;
	int memsize;

	//  Check for a sane substring
	if(end-begin <= 0){
		fprintf(stderr, "\nInvalid argument to wordToDouble!\n");
		fflush(stderr);
		exit(0);
	}

	//  Get the required string size
	memsize = end - begin + 2; // An extra space in mem for null-termination
	size = end -begin +1;

	//  Get memory for a temporary string
	str = (char *) malloc ( memsize * sizeof(char));
	if ( str == NULL){
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(0);
	}
	// Make sure the string ends with a null char
	str[size] = '\0';

	//  Copy the word into str
	str = strncpy(str, s+begin, size);

	//  Set endptr to str as initial value
	endptr = str;
	val = strtod(str, &endptr);

	//  endptr should point to the last char
	//  used in the conversion if strtod worked
	if(val == 0 && endptr == str){
		fprintf(stderr, "\nError making double from string: %s\n", str);
		fflush(stderr);
		exit(0);
	}
	free(str);
	return val;
}

// Count the number of fields (delimited by 'sep')
// in a single string.  I was going to use strsep
// in string.h for this; however, I don't like that
// it changes the input string, which makes free-ing
// the string later more tricky.
//
// Ignores consecutive seperators
int countFields(char *s, char sep){
	int i;
	int begin=0;
	int end=0;
	int status=0; // 0 = in sep, 1 = in word
	int fieldCount=0;
	double val;

	if (s==NULL){
		fprintf(stderr, "Passed NULL string to countFields -- error!");
		fflush(stderr);
		exit(0);
	}

	//  Loop over the length of the string
	for( i=0 ; i<strlen(s) ; i++ ){
		//  The previous state was space
		if(status == 0){
			//  We hit a word
			if(s[i] != sep){
				begin = i;
				status = 1;
			}else{ // We hit more space
				continue;
			}
		}else{ // The previous state was word
			if(s[i] != sep){
				continue;
			}else{ // We hit a space
				end = i-1;
				status = 0;

				//  being and end now delimit a word,
				//  turn that word into a double
				val = wordToDouble(s, begin, end);
				fieldCount++;
			}
		}
	}
	//  At the end, if we were in a word, we have
	//  one more field
	if(status == 1){  // We're in a word
		val = wordToDouble(s, begin, strlen(s));
		fieldCount++;
	}

	
	return fieldCount;
}

//  Check that each sequence has the same dimensionality
//  and that, within a sequence, each dimension has the
//  same number of entries.
//
//  NOTE: THIS ROUTINE ALTERS *nunSeq_p and *dim_p!!!!
//
//  Also, you MUST call this routine before calling 
//  parseRealData.  Otherwise, parseRealData is garunteed
//  to barf if the data turns out to be ill-formatted.
int checkRealDataFormat(char **buf, int nl, char sep, int *numSeq_p, int *dim_p){
	int i;
	int thisDim=0;
	int status=1;
	int width;
	int fieldCount = 0;  // number of positions in a single sequence
	int numSeq = 0; // number of sequences 
	int dim = 0;  // The dimensionality of the sequences
	
// NOTE this is not checking the dimensionality of the last sequence...
//    that's bad.  We can fix that though.
	//  Check the dimensionality of each sequence
	for ( i=0 ; i<nl ; i++){
		if(buf[i][0] == '>'){
			//  If this is only the second sequence we've seen,
			//  record the dimensionality of the first sequence
			//  as the dim to insist upon from here on out
			if(numSeq == 1){
				dim = thisDim;

			//  For other sequences, we need to check to make sure
			//  that they've got the same dimensions as previous
			//  sequences
			}else if(numSeq > 1){
				//  If the dimensions are wrong, quit with status=0
				if(thisDim != dim ){
					status = 0;
					break;
				}
			}
			numSeq++;
			width=0;
			thisDim=0;
		}else{

			//  Field count can be different for each sequence but
			//  must be the same for each dimension in a single sequence
			fieldCount = countFields(buf[i], sep);

			//  If this is the first row of this sequence,
			//  then store the number of fields
			if(thisDim == 0){
				width = fieldCount;

			//  If it's not the first row, make sure it has the
			//  same number of fields as previous rows in this
			//  sequence
			}else{
				if(fieldCount != width){
					status = 0;
					break;
				}
			}
			thisDim++;
		}
	}

	//  Pass back the numSeq and dim
	*numSeq_p = numSeq;
	*dim_p = thisDim;
	return status;
}

//  Count the number of fields in each sequence and
//  return the sum of these
int countTotalFields(char **buf, int nl, char sep){
	int i=0;
	int totalFields=0;
	int seqNo=0;
	while(i<nl){
		// Hit a new sequence
		if(buf[i][0] == '>'){
			seqNo++;
			//  Assume that the sequence has at least
			//  one row (should have called checkRealDataFormat!
			//  and that each row has the same number of fields
			totalFields += countFields(buf[i+1], sep);
		}
		i++;
	}
	return totalFields;
}

rdh_t *
initRdh(int x){
	int i;
	rdh_t *data = NULL;

	//  Allocate space for our structure
	data = (rdh_t *) malloc ( sizeof(rdh_t) );
	if(data == NULL){
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(0);
	}

	data->size = x;

	//  Index has to be initialized later, once
	//  we know the word size.
	data->indexSize = 0;
	data->indexToSeq = NULL;
	data->indexToPos = NULL;
	/*data->indexSize = y;*/

	data->label = (char **) malloc ( data->size * sizeof(char *) );
	if(data->label == NULL){
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(0);
	}
	data->seq = (gsl_matrix **) malloc ( data->size * sizeof(gsl_matrix *) );
	if(data->seq == NULL){
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(0);
	}
	for( i=0 ; i<data->size ; i++){
		data->label[i] = NULL;
		data->seq[i] = NULL;
	}
	return data;
}

int
getRdhSeqLength(rdh_t *data, int seqNo){
	if (data==NULL|| data->seq == NULL || data->seq[seqNo] == NULL){
		fprintf(stderr, "Passed bad data to getRdhSeqLength -- error!");
		fflush(stderr);
		exit(0);
	}
	return data->seq[seqNo]->size1;
}

//  seqgap is added in case we later want to do
//  wierd shit with the convolution that would 
//  require large gaps between sequences in the
//  bitGraph_t.
int
initRdhIndex(rdh_t *data, int wordSize, int seqGap){
	int i,j,k;
	int numWindows=0;
	int thisNumWindows;
	int numSeq;
	int seqLen = 0;

	//  The number of sequences
	numSeq = data->size;
	
	//   Allocate offsetToIndex's outer structure
	data->offsetToIndex = (int **) malloc (numSeq * sizeof(int *));
	if(data->offsetToIndex == NULL) {
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(0);
	}


	// For each sequence
	for ( i=0 ; i<numSeq ; i++){
		// How many windows are in this sequence
		seqLen = getRdhSeqLength(data, i);
		numWindows += seqLen - wordSize + 1;
		// And also use this to further allocate offsetToIndex
		data->offsetToIndex[i] = (int *) malloc ((seqLen - wordSize + 1) * sizeof(int));
		if(data->offsetToIndex[i] == NULL) {
			fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
			fflush(stderr);
			exit(0);
		}
	}

	//  One index for each word plus seqGap between each sequence
	//  and a gap at the end
	data->indexSize = numWindows + numSeq*seqGap;

	//  Allocate indexToSeq 
	//    NOTE that it should be size of int, not int *... I think we got
	//    fortunate in the previous revision because they are the same
	//    size
	data->indexToSeq = (int *) malloc ( data->indexSize * sizeof(int) );
	if(data->indexToSeq == NULL){
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(0);
	}
	//  Allocate indexToPos 
	//  See above for int vs. int* argument.
	data->indexToPos = (int *) malloc ( data->indexSize * sizeof(int) );
	if(data->indexToPos == NULL){
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(0);
	}
	
	//  Fill in the values
	k=0;
	for ( i=0 ; i<numSeq ; i++){
		//  How many windows are in this sequence?
		thisNumWindows =  getRdhSeqLength(data, i) - wordSize + 1;

		// For each window, make an entry in the indexToSeq
		// and indexToPos and offsetToIndex
		for (j=0 ; j<thisNumWindows ; j++){
			data->indexToSeq[k] = i;
			data->indexToPos[k] = j;
			data->offsetToIndex[i][j] = k;
			k++;
		}
		//  Add gaps between sequences in the index.
		//  Usually seqGap is just 1;
		for (j=0 ; j<seqGap ; j++){
			//  -1 means no sequence and no position
			data->indexToSeq[k] = -1;
			data->indexToPos[k] = -1;
			k++;
		}
	}
	return 0;
}

rdh_t *
freeRdh(rdh_t *data){
	int i;
	if(data != NULL){
		if(data->indexToPos != NULL){
			free(data->indexToPos);
			data->indexToPos = NULL;
		}
		if(data->indexToSeq != NULL){
			free(data->indexToSeq);
			data->indexToSeq = NULL;
		}
		if(data->offsetToIndex != NULL) {
			for (i = 0; i < data->size; i++) {
				free(data->offsetToIndex[i]);
				data->offsetToIndex[i] = NULL;
			}
			free(data->offsetToIndex);
			data->offsetToIndex = NULL;
		}
		for( i=0 ; i<data->size ; i++){
			if(data->seq[i] != NULL){
				gsl_matrix_free(data->seq[i]);
				data->seq[i] = NULL;
			}
			if(data->label[i] != NULL){
				free(data->label[i]);
				data->label[i] = NULL;
			}
		}
		if(data->seq != NULL){
			free(data->seq);
			data->seq = NULL;
		}
		if(data->label != NULL){
			free(data->label);
			data->label = NULL;
		}
		free(data);
		data = NULL;
	}
	return data;
}

int
getRdhDim(rdh_t *data){
	if (data==NULL|| data->seq == NULL || data->seq[0] == NULL){
		fprintf(stderr, "Passed bad data to getRdhSeqLength -- error!");
		fflush(stderr);
		exit(0);
	}
	return data->seq[0]->size2;
}

int
setRdhLabel(rdh_t *data, int seqNo, char *s){
	if ( data->seq == NULL || data->label == NULL){
		fprintf(stderr, "Passed bad data to setRdhLabel -- error!");
		fflush(stderr);
		exit(0);
	}
	data->label[seqNo] = strdup(s);
	if(data->label[seqNo] == NULL){
		fprintf(stderr, "\nMemory Error allocating label!\n%s\n", strerror(errno));
		fflush(stderr);
		exit(0);
	}
	return 0;
}

int
setRdhValue(rdh_t *data, int seqNo, int posNo, int dimNo, double val){
	if (data==NULL|| data->seq == NULL || data->seq[seqNo] == NULL ||
				posNo>getRdhSeqLength(data, seqNo) ||
				dimNo>getRdhDim(data)){
		fprintf(stderr, "Passed bad data to setRdhValue -- error!");
		fflush(stderr);
		exit(0);
	}
	gsl_matrix_set( data->seq[seqNo], posNo, dimNo, val);
	return 0;
}


int
setRdhIndex(rdh_t *data, int seqNo, int posNo, int index){
	if ( data==NULL || data->indexToSeq == NULL ||
			data->indexToPos==NULL || index > data->indexSize){
		fprintf(stderr, "Passed bad data to getRdhValue -- error!");
		fflush(stderr);
		exit(0);
	}
	/*printf("Setting index %d -> %d, %d\n", index, seqNo, posNo);*/
	/*fflush(stdout);*/
	data->indexToSeq[index] = seqNo;
	data->indexToPos[index] = posNo;
	return 0;
}

// Alters seq and pos!
int
getRdhIndexSeqPos(rdh_t *data, int index, int *seq, int *pos){
	if ( data==NULL || data->indexToSeq == NULL ||
			data->indexToPos==NULL || index > data->indexSize){
		fprintf(stderr, "Passed bad data to getRdhIndexSeqPos -- error!");
		fflush(stderr);
		exit(0);
	}
	/*printf("Setting index %d -> %d, %d\n", index, seqNo, posNo);*/
	/*fflush(stdout);*/
	*seq = data->indexToSeq[index];
	*pos = data->indexToPos[index];
	return 0;
}

double
getRdhValue(rdh_t *data, int seqNo, int posNo, int dimNo){
	if ( data==NULL || data->seq == NULL || data->seq[seqNo] == NULL ||
				posNo>getRdhSeqLength(data, seqNo) ||
				dimNo>getRdhDim(data)){
		fprintf(stderr, "Passed bad data to getRdhValue -- error!");
		fflush(stderr);
		exit(0);
	}
	return gsl_matrix_get( data->seq[seqNo], posNo, dimNo);
}

char *
getRdhLabel(rdh_t *data, int seqNo){
	if ( data==NULL || data->label == NULL || data->label[seqNo] == NULL){
		fprintf(stderr, "Passed bad data to getRdhLabel -- error!");
		fflush(stderr);
		exit(0);
	}
	return data->label[seqNo];
}


int
printRdhSeq(rdh_t *data, int seqNo, FILE *FH){
	int i,j;
	int len;
	int dim;
	len = getRdhSeqLength(data, seqNo);
	dim = getRdhDim(data);
	fprintf(FH, "%s\n", getRdhLabel(data, seqNo));
	for( i=0 ; i<len ; i++){
		for( j=0 ; j<dim ; j++){
			fprintf(FH, "%3.1f ", getRdhValue(data, seqNo, i, j));
		}
		fprintf(FH, "\n");
	}
	return 0;
}

//  Similar code to countFields
int 
setRdhColFromString (rdh_t *data, int seqNo, int colNo, char *s, char sep){
	int i;
	int begin=0;
	int end=0;
	int status=0; // 0 = in sep, 1 = in word
	int fieldCount=0;
	double val;


	//  Make sure the string is not null and
	//  the rdh_t gsl_matrix array is not null
	//  and the selected gsl_matrix is not null
	if (s == NULL || data->seq == NULL || data->seq[seqNo] == NULL){
		fprintf(stderr, "Passed bad data to setRdhColFromString -- error!");
		fflush(stderr);
		exit(0);
	}

	//  Loop over the length of the string
	for( i=0 ; i<strlen(s) ; i++ ){
		//  The previous state was space
		if(status == 0){
			//  We hit a word
			if(s[i] != sep){
				begin = i;
				status = 1;
			}else{ // We hit more space
				continue;
			}
		}else{ // The previous state was word
			if(s[i] != sep){
				continue;
			}else{ // We hit a space
				end = i-1;
				status = 0;
				val = wordToDouble(s, begin, end);

				//  Go to the gsl_matrix object data->seq[seqNo]
				//  and set the (fieldCount, colNo) = val;
				setRdhValue(data, seqNo, fieldCount, colNo, val);
				fieldCount++;
			}
		}
	}
	//  At the end, if we were in a word, we have
	//  one more field
	if(status == 1){  // We're in a word
		val = wordToDouble(s, begin, strlen(s));
		// Added in, MPS 5/3/05 ---
		//  And don't forget to set the RdhValue!
		setRdhValue(data,seqNo,fieldCount,colNo,val);
		fieldCount++;
	}

	
	return fieldCount;
}

int 
initRdhGslMat(rdh_t *data, int seqNo, int x, int y){
	data->seq[seqNo] = gsl_matrix_alloc(x, y);
	if(data->seq[seqNo] == NULL){
		return 0;
	}else{
		return 1;
	}
}

//  This is the ONLY routine that should be used to
//  alter the rdh_t structure as we're reading in 
//  sequences.  This routine uses a few static variables
//  to keep track of the state of the rdh_t object.
int 
pushOnRdhSeq(rdh_t *data, char **buf, int startLine, int dim, char sep){
	int i,j,k;
	int numFields;


	// NOTE THAT THESE ARE STATIC VARIABLES!!!!!
	// That is, they retain their last value on
	// each call to this function!
	static int seqNo=0;
	/*static int indexNo=0;*/



	i = startLine;

	//  Assume that the sequence has at least
	//  one row (should have called checkRealDataFormat!
	numFields = countFields(buf[i+1], sep);

	//  Initialize the gsl_matrix object for this 
	//  sequence in 'data'
	//
	//  NOTE THAT WE STORE THE TRANSPOSE OF WHAT'S IN
	//  THE INPUT FILE -- x,y = position x, dimension y
	initRdhGslMat(data, seqNo, numFields, dim);

	// Set the sequence label
	setRdhLabel(data, seqNo, buf[i]);

	//  Read in 'dim' rows
	for( j=i+1, k=0 ; j<i+1+dim ; j++, k++){
		/*printf("%d\n", countFields(buf[j], sep));*/

		//  Set the k-th dimension of this sequence
		//  STILL NOTE THE TRANSPOSE!
		setRdhColFromString (data, seqNo, k, buf[j], sep);

	}
	/*
	for ( l=0 ; l<numFields ; l++ ){
		setRdhIndex(data, seqNo, l, indexNo);
		indexNo++;
	}
	*/
	seqNo++;
	// Augment indexNo once more to have a -1 between each sequence!
	/*indexNo++;  */
	return 0;
}

//  Parse the data
rdh_t *
parseRealData(char **buf, int nl, char sep, int numSeq, int dim){
	int i;
	int seqNo=-1;
	int totalNumFields;
	rdh_t *data = NULL;

	totalNumFields = countTotalFields(buf, nl, sep);

	/*data = initRdh(numSeq, totalNumFields + numSeq - 1);*/
	data = initRdh(numSeq);


	//  We're going to add an empty index between
	//  windows that correspond to different 
	//  sequences

	// Fast forward to the first sequence
	i=0;
	while(i<nl){
		// Hit a new sequence
		if(buf[i][0] == '>'){
			seqNo++;  // Note that seqNo started at -1!

			pushOnRdhSeq(data, buf, i, dim, sep);
			i+=dim+1;
		}else{
			i++;
		}
	}

	/*printRdhSeq(data, 0, stdout);*/

	return data;
}

rdh_t *readRealData(FILE * INPUT){
	char **buf = NULL;
	int nl;
	int i;
	char sep = ' ';
	int numSeq = 0;
	int dimensions = 0;
	int status = 1;
	rdh_t *data = NULL;

	// Read the entire INPUT file and put it's
	// contents into 'buf'.  This function also
	// alters the contents of the location pointed
	// to by &nl.  Now nl is the number of lines
	// in the file (or the size of the buff array.
	buf = ReadFile(INPUT, &nl);
	if(buf == NULL){
		return NULL;
	}
	status = checkRealDataFormat(buf, nl, sep, &numSeq, &dimensions);
	if(numSeq <= 0 || dimensions <= 0 || status == 0){
		fprintf(stderr, "Data file is poorly formatted or no sequences read!\n");
		fprintf(stderr, "Each sequence needs to be the same dimensionality!  QUITTING!\n");
		fprintf(stderr, "numSeq = %d, dimensions = %d, status = %d\n",numSeq,dimensions,status);
		exit(EXIT_FAILURE);
	}

	//  From here on, we assume that the sequence file is well-formatted
	//  to make the code more simple.
	data = parseRealData(buf, nl, sep, numSeq, dimensions);

	//  Free up our buffer
	for ( i=0 ; i<nl ; i++ ){
		if(buf[i] != NULL){
			free(buf[i]);
		}
	}
	if(buf != NULL){
		free(buf);
	}
	return data;
}

int
outputRealPats(rdh_t *data, cll_t *allPats, int L, FILE *OUTPUT_FILE, 
		int **d){
	int i,j,pos1;
	int curSeq, curPos;
	cll_t *curCliq = NULL;
	
	curCliq = allPats;
	i=0;
	while(curCliq != NULL){
		fprintf(OUTPUT_FILE,"pattern %d:\tlen=%d\tsup=%d\t", 
			i, curCliq->length+L, curCliq->set->size);
		if (d != NULL) {
			fprintf(OUTPUT_FILE,"\tsignif=%le\n",curCliq->stat);
		} else {
			fprintf(OUTPUT_FILE,"\n");
		}
		for(j=0 ; j<curCliq->set->size ; j++){
			pos1 = curCliq->set->members[j];
			getRdhIndexSeqPos(data, pos1, &curSeq, &curPos);

			fprintf(OUTPUT_FILE, "   %d\t%d\t", curSeq, curPos);
			fprintf(OUTPUT_FILE, "%lf\t",
				gsl_matrix_get(data->seq[curSeq],curPos,0));
			/*
			for(k=curPos ; k<curPos+curCliq->length+L ; k++){
				fprintf(OUTPUT_FILE, "%c", mySequences[curSeq].seq[k]);
			}
			*/
			fprintf(OUTPUT_FILE, "\n");
		}
		fprintf(OUTPUT_FILE, "\n\n");
		curCliq = curCliq->next;
		i++;
	}
	return 0;
}

// We'll settle for taking the integer value of the comparison function
//   and calling the getCompFunc in here so we don't have to worry
//   about passing function pointers as arguments... this is already 
//   confusing enough.
//
int
findCliqueCentroid(rdh_t *data, cll_t *curCliq, int L,int compFunc,
		double *extraParams, int *candidates) {
	double (*comparisonFunc) (rdh_t*,int,int,int,double*) = NULL;
	int i = 0, j = 0, indmin = -1, counter = 0;
	double sim = 0, min = 0, flagmin = 0;
	double *cliqueAdjMat = NULL;

	cliqueAdjMat = (double *) malloc(curCliq->set->size * sizeof(double));
	if (cliqueAdjMat == NULL) {
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(0);
	}

	for (i = 0; i < curCliq->set->size; i++) {
		cliqueAdjMat[i] = 0;
	}
	
	// We'll accumulate our comparison function values... except here
	//  we're really assuming that we're using a match factor, with
	//  value less than one, so that we can subtract it from one to
	//  get a distance, and then find the centroid by identifying the
	//  node with the smallest cumulative Euclidean distance to all
	//  nodes.  
	// Note that we only need to compare each unique pair, and can apply
	//  the results from each comparison to each member of the pair,
	//  hence the somewhat odd indices of initiation for the for loops.
	comparisonFunc = getCompFunc(compFunc);
	for (i = 0; i < curCliq->set->size; i++) {
		for (j = i+1; j < curCliq->set->size; j++) {
			sim = comparisonFunc(data,curCliq->set->members[i],
					curCliq->set->members[j],L,extraParams);
//			printf("i = %d, j = %d, L = %d, extra = %lf, sim = %lf\n",i,j,L,extraParams[0],sim);
			cliqueAdjMat[i] += pow(1-sim,2);
			cliqueAdjMat[j] += pow(1-sim,2);
		}
	}

	// Now we find the minimum Euclidean distance.
	
	min = cliqueAdjMat[0];
	indmin = 0;
	for (i = 1; i < curCliq->set->size; i++) {
//		printf("index %d product = %lf\n",i,cliqueAdjMat[i]);
		if (cliqueAdjMat[i] < min) {
			indmin = i;
			min = cliqueAdjMat[i];
			flagmin = 0;
		} else if (cliqueAdjMat[i] == min) {
			flagmin = 1;
		}
	}

	// If we had a duplicate on the minimum, we locate all duplicates.
	if (flagmin == 1) {
		counter = 0;
		for (i = 0; i < curCliq->set->size; i++) {
			if(cliqueAdjMat[i] == min) {
				counter++;
				candidates[counter] = i;
			}
		}
		// Store the number of candidates at the array's beginning
		candidates[0] = counter;
		free(cliqueAdjMat);
		return(-1);
	} else {
		free(cliqueAdjMat);
		return(indmin);
	}
}
// In order to make the centroid decision slightly less dependent on 
//  input order, we decide to choose from the tied candidates the one
//  whose relative position in the sequence is highest.  There is no basis
//  in theory for this, it is done so that a consistent choice is made.  Only
//  rarely will two spectra be tied for being a centroid and have the same
//  sequence number... in that case, we pretty much have to default
//  to the sequence number, which is what would be done without this function.
// Note that now though we are less sensitive to the order of input of the
//  sequences, we are now more sensitive to the context surrounding a given
//  spectrum.  That is, if it is put in the beginning of the "sequence",
//  it is more likely to be chosen.  
// This choice can only be justified insofar as if multiple choices are tied,
//  then they are the same cumulative distance to the clique, and so *any*
//  should be allowed to be chosen equally.  There should be little difference
//  in terms of tangible results.  This just makes the semantics consistent.

int makeAlternateCentroid(rdh_t *data, cll_t *curCliq, int *candidates) {
	int indmin, min, i;
	int curSeq, curPos;
	int numCandidates = candidates[0];

	indmin = candidates[1];
	getRdhIndexSeqPos(data,curCliq->set->members[indmin],&curSeq,&curPos);
	min = curPos;

	// We use less-than-or-equal here because we're starting at 1, 
	//  so we want 1 to end.  The length of candidates is one more than
	//  the maxSup, so we know we can reach candidates[maxSup] without
	//  a segfault.
	for (i = 2; i <= numCandidates; i++) {
		getRdhIndexSeqPos(data,curCliq->set->members[candidates[i]],
				&curSeq,&curPos);
		if (curPos < min) {
			indmin = candidates[i];
			min = curPos;
		}
	}

	return(indmin);
}
		
int
outputRealPatsWCentroid(rdh_t *data, cll_t *allPats, int L, FILE *OUTPUT_FILE,
		double *extraParams,int compFunc){
	int i,j,k,pos1,centroid;
	int curSeq, curPos;
	int maxSup = 0;
	cll_t *curCliq = NULL;
	double mfToCentroid = 0;
	double (*comparisonFunc) (rdh_t*,int,int,int,double*) = NULL;
	int *candidates = NULL;
	
	curCliq = allPats;
	while(curCliq != NULL) {
		if(curCliq->set->size > maxSup) {
			maxSup = curCliq->set->size;
		}
		curCliq = curCliq->next;
	}
	candidates = (int *) malloc((maxSup+1) * sizeof(int));
	if (candidates == NULL) {
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(0);
	}
	for (i = 0; i <= maxSup; i++) {
		candidates[i] = 0;
	}
	
	comparisonFunc = getCompFunc(compFunc);
	curCliq = allPats;
	i=0;
	while(curCliq != NULL){
		fprintf(OUTPUT_FILE, "pattern %d:\tlen=%d\tsup=%d\n", 
				i, curCliq->length+L, curCliq->set->size);
		centroid = findCliqueCentroid(data,curCliq,L,compFunc,
				extraParams,candidates);
		if (centroid < 0) {
			centroid = makeAlternateCentroid(data,curCliq,
					candidates);
//			fprintf(OUTPUT_FILE, "WARNING: No single node in"
//			" cluster has non-zero similarity to all other\n nodes"
//			" in cluster; centroid set to first node.\n");
//			centroid = 0;
		}
		for(j=0 ; j<curCliq->set->size ; j++){
			pos1 = curCliq->set->members[j];
			getRdhIndexSeqPos(data, pos1, &curSeq, &curPos);
			fprintf(OUTPUT_FILE, "   %d\t%d\t", curSeq, curPos);
//			fprintf(OUTPUT_FILE, "%lf\t",
//				gsl_matrix_get(data->seq[curSeq],curPos,0));
			mfToCentroid = comparisonFunc(data,
					curCliq->set->members[j],
					curCliq->set->members[centroid],
					L,extraParams);
			fprintf(OUTPUT_FILE, "%lf\t",
					mfToCentroid);
			/*
			for(k=curPos ; k<curPos+curCliq->length+L ; k++){
				fprintf(OUTPUT_FILE, "%c", mySequences[curSeq].seq[k]);
			}
			*/
			fprintf(OUTPUT_FILE, "\n");
		}
		fprintf(OUTPUT_FILE, "\n\n");
		curCliq = curCliq->next;
		i++;
		for (k = 0; k <= maxSup; k++) {
			candidates[k] = 0;
		}
	}
	free(candidates);
	return 0;
}

