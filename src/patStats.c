#include <math.h>
#include "patStats.h"

int getLargestSupport(cll_t *cliqs){
	int size=0;
	cll_t *curCliq = NULL;

	curCliq = cliqs;
	while(curCliq != NULL){
		if(curCliq->set->size > size){
			size = curCliq->set->size;
		}
		curCliq = curCliq->next;
	}
	return size;
}

int getLargestLength(cll_t *cliqs){
	int len=0;
	cll_t *curCliq = NULL;

	curCliq = cliqs;
	while(curCliq != NULL){
		if(curCliq->length > len){
			len = curCliq->length;
		}
		curCliq = curCliq->next;
	}

	// We return (len + 1) because the length of the shortest streak
	//  is one, but is stored in the cluster data structure as being
	//  zero (number of extensions that have been made).
	return (len + 1);
}

int measureDiagonal(const bitGraph_t *bg, const int i, const int j){
	int len = 0;
	while(bitGraphCheckBit(bg, i+len, j+len) != 0){
		len++;
	}
	return len;
}

// dimToChange is 1 for the first dimension (support), 2 for the second 
//  dimension (length).  newVal is the new value for the dimension to be
//  changed, not including the "1" that should be added... so it should
//  just be some integer times the initial support.

unsigned int**
increaseMem(unsigned int **d,int dimToChange,int currSupport,int currLength,int newVal){ 
	int i = 0, j = 0;
	
	if (dimToChange == 1) {
		d = (unsigned int **) realloc(d, (newVal+1) * sizeof(unsigned int*));
		if (d == NULL) {
			fprintf(stderr,"\nMemory error --- couldn't allocate array!"
					"\n%s\n",strerror(errno));
			fflush(stderr);
			exit(0);
		}

		for (i = currSupport + 1; i < newVal + 1; i++) {
			d[i] = (unsigned int *) malloc((currLength + 1) * sizeof(unsigned int));
			if (d[i] == NULL) {
				fprintf(stderr,"\nMemory error --- couldn't allocate array!"
						"\n%s\n",strerror(errno));
				fflush(stderr);
				exit(0);
			}
			for (j = 0; j < currLength+1; j++) {
				d[i][j] = 0;
			}
		}
		return d;
	} else if (dimToChange == 2) {
		for (i = 0; i < currSupport + 1; i++) {
			d[i] = (unsigned int *) realloc(d[i],(newVal + 1) * sizeof(unsigned int));
			if (d[i] == NULL) {
				fprintf(stderr,"\nMemory error --- couldn't allocate array!"
						"\n%s\n",strerror(errno));
				fflush(stderr);
				exit(0);
			}
			for (j = currLength+1; j < newVal + 1; j++) {
				d[i][j] = 0;
			}
		}
		return d;
	} else {
		fprintf(stderr,"Invalid arguments to increaseMem!\n\n");
		fflush(stderr);
		exit(0);
	}
}



	// OK, here is something that is a little bit "hackish" but that we
	//  have to do.  Since our initial matrix is being pruned and filtered
	//  before being clustered, but we need to calculate stats based
	//  on the original matrix, we need to get information from the 
	//  matrix before pruning, so we're using this function.  We could
	//  just make a copy of that matrix, but it's far too big, and that
	//  would cause an unneccessary constraint on memory, limiting the
	//  size of problems we can address.  But we need to define just how
	//  big our d matrix is before we can use it.  We could go through
	//  and compute the longest streak beforehand, and then redo 
	//  everything, but we've already found the first step of finding
	//  all of the streaks to be fairly expensive (KLJ).  So instead 
	//  what we'll do is use the user's parameters as a benchmark and
	//  expand from there.  We'll assume that most of the time, the 
	//  biggest streak (number of extensions) will be less than 50 times
	//  the length given as input by the user, and the biggest support
	//  will be less than 50 times the minimum number of support given
	//  by the user.  This seems perhaps overly conservative, but 
	//  otherwise is reasonable.  We then realize that even on a 64-bit
	//  computer, if the user gives L=50 and K=50, we'll still use
	//  less than 48 MB of memory... and if L=50 and K=50, it is 
	//  extremely likely that doubling the adjacency matrix would
	//  have been a much worse option.  Scaling back to more common values
	//  of L~20 and K~20, the memory used shoots down to ~9MB, which is
	//  definitely acceptable.
	// Now, if for some reason our initial allocation wasn't enough,
	//  then we'll have to go through and realloc all of our memory again.
	//  Somewhat time-consuming, but hopefully not done too often.
	//  Each time we find we try to put something in an index that
	//  doesn't exist, we'll reallocate our memory, adding twice as much
	//  in the dimension that was violated.
	// It is important to us that we get back the final dimensions of this
	//  matrix, since in the support dimension we'll have to sum across
	//  all values, and in the length dimension we'll have to be sure
	//  we're not at the edge of a matrix during our d manipulations
	//  later on.

unsigned int** 
oldGetStatMat(bitGraph_t *bg, int support, int length, int *supportDim, 
		int *lengthDim, int numBlanks){
	int *Q = NULL;
	unsigned int **d = NULL;
	int i,j;
	int x,y;
	bitSet_t *X = NULL;
	int currSupport;
	int currLength;
	int multiplier = 50;
	time_t probStart, probEnd;
	int timeNeeded = 0;

	currSupport = support*multiplier;
	currLength = length*multiplier;

	X = newBitSet(bg->size);
//	printf("Made bitSet of size %d\n", bg->size);
	
	Q = (int *) malloc ( bg->size * sizeof(int) );
	if (Q == NULL) {
		fprintf(stderr,"\nMemory error --- couldn't allocate array!"
				"\n%s\n",strerror(errno));
		fflush(stderr);
		exit(0);
	}
	for ( i=0 ; i<bg->size ; i++){
		Q[i]=0;
	}


	d = (unsigned int **) malloc ( (currSupport+1) * sizeof(unsigned int *) );
	if (d == NULL) {
		fprintf(stderr,"\nMemory error --- couldn't allocate array!"
				"\n%s\n",strerror(errno));
		fflush(stderr);
		exit(0);
	}

	for ( i=0 ; i<currSupport+1 ; i++){
		d[i] = (unsigned int *) malloc ( (currLength+1) * sizeof(unsigned int) );
		if (d[i] == NULL) {
			fprintf(stderr,"\nMemory error --- couldn't allocate array!"
					"\n%s\n",strerror(errno));
			fflush(stderr);
			exit(0);
		}
		for ( j=0 ; j<currLength+1 ; j++){
			d[i][j]=0;
		}
	}


	time(&probStart);
	for ( i=0 ; i<bg->size ; i++){
		if(i == 200) {
			time(&probEnd);
			timeNeeded = ((double) (probEnd-probStart)) / 
					((double) 60)*((double) bg->size)/
					((double) 200);
			if (timeNeeded > 2) {
				printf("Max total time to calculate probability:\n");
				printf("\t%d minutes\n",timeNeeded);
				printf("Actual time will be less than this, but at",
					"least half of it.\n");
				printf("To bypass excessive probability calculations,",
					"cancel and use the '-d' flag.\n");
				fflush(NULL);
			}
		}

		for ( j=bg->size-1 ; j>i ; j--){
			bitGraphRowIntersection(bg, i, j, X);
			x = countSet(X);
			if(Q[j-1] != 0){
				y = Q[j-1] - 1;
				Q[j] = Q[j-1] - 1;
			}else{
				y = measureDiagonal(bg,i,j);
				Q[j] = y;
			}
			while(x > currSupport){
				d = increaseMem(d,1,currSupport,currLength,
					currSupport + support*multiplier);
				currSupport += support*multiplier;
			}
			while(y > currLength){
				d = increaseMem(d,2,currSupport,currLength,
					currLength + length*multiplier);
				currLength += length*multiplier;
			}
			d[x][y]++;
			/*
			if(x != 0){
				printf("%d:\t%d %d\n", j, x, y);
				fflush(stdout);
			}
			*/
		}
		/*
		printf("done\n");
		fflush(stdout);
		*/
	}

	// We know that the "blanks", inserted to delimit unique sequences
	//  and prevent convolution through them, will skew our statistics,
	//  so we subtract them.  We know that they will never be similar to
	//  any others, so will only add to the d[0][0] number.  Furthermore,
	//  we know how many they add.  Since d never hits the main diagonal
	//  and only does the upper half of the matrix, the first one 
	//  contributes bgsize - 1 to d[0][0], the next bgsize - 2, etc.
	for (i = 0; i < numBlanks; i++) {
		d[0][0] -= bg->size - 1 - i;
	}

	deleteBitSet(X);
	free(Q);
	*supportDim = currSupport;
	*lengthDim = currLength;
	return(d);
	
}

unsigned int** getStatMat(bitGraph_t *bg, int support, int length, int *supportDim, 
		int *lengthDim, int numBlanks, int s, FILE *OUTPUT_FILE){
	int *Q = NULL;
	unsigned int **d = NULL;
	int i,j,k;
	int x,y;
	bitSet_t *X = NULL;
	int currSupport;
	int currLength;
	int multiplier = 50;
	int diagonal = 0;
	time_t probStart, probEnd;
	int timeNeeded = 0;
	int sampleCounter = 1;
//	int visitCounter = 0, uniqCounter = 0;


	currSupport = support*multiplier;
	currLength = length*multiplier;

	X = newBitSet(bg->size);
//	printf("Made bitSet of size %d\n", bg->size);
	
	Q = (int *) malloc ( bg->size * sizeof(int) );
	if (Q == NULL) {
		fprintf(stderr,"\nMemory error --- couldn't allocate array!"
				"\n%s\n",strerror(errno));
		fflush(stderr);
		exit(0);
	}
	for ( i=0 ; i<bg->size ; i++){
		Q[i]=0;
	}


	d = (unsigned int **) malloc ( (currSupport+1) * sizeof(unsigned int *) );
	if (d == NULL) {
		fprintf(stderr,"\nMemory error --- couldn't allocate array!"
				"\n%s\n",strerror(errno));
		fflush(stderr);
		exit(0);
	}

	for ( i=0 ; i<currSupport+1 ; i++){
		d[i] = (unsigned int *) malloc ( (currLength+1) * sizeof(unsigned int) );
		if (d[i] == NULL) {
			fprintf(stderr,"\nMemory error --- couldn't allocate array!"
					"\n%s\n",strerror(errno));
			fflush(stderr);
			exit(0);
		}
		for ( j=0 ; j<currLength+1 ; j++){
			d[i][j]=0;
		}
	}
//	printf("size=%d\n",bg->size);
	time(&probStart);
	for ( i=0 ; i<bg->size ; i++){
		if(i == 200) {
			time(&probEnd);
			timeNeeded = ((double) (probEnd-probStart)) / 
					((double) 60)*((double) bg->size)/
					((double) 200);
			if (timeNeeded > 2) {
				fprintf(OUTPUT_FILE,"Max total time to calculate probability:\n");
				fprintf(OUTPUT_FILE,"\t%d minutes\n",timeNeeded);
				fprintf(OUTPUT_FILE,"Actual time will be less than this, " 
					"but at least half of it.\n");
				fprintf(OUTPUT_FILE,"To bypass excessive probability calculations,"
					" cancel and use a different value\n"
					" for the '-s' flag (samples every "
					"'s' points).\n");
				fflush(NULL);
			}
		}
		j = nextBitBitSet(bg->graph[i],0);
		while(j >= 0) {
			k = nextBitBitSet(bg->graph[i],j+1);
			while(k >= 0) {
				if(checkBit(bg->graph[j],k) == 0) {
					if (sampleCounter == s) {
						bitGraphRowIntersection(bg,j,k,X);
//					visitCounter++;
						if (nextBitBitSet(X,0) >= i) {
//						uniqCounter++;
							x = countSet(X);
							while(x > currSupport) {
								d = increaseMem(d,1,
									currSupport,
									currLength,
									currSupport + 
									support*multiplier);
								currSupport += support*
									multiplier;
							}
							d[x][0] += 1;
						}
						sampleCounter = 0;
					}
					sampleCounter++;
				}
				k = nextBitBitSet(bg->graph[i],k+1);
			}
			if (j <= i) {
				j = nextBitBitSet(bg->graph[i],j+1);
				continue;
			}
			bitGraphRowIntersection(bg, i, j, X);
			x = countSet(X);
			// Note, now we're using "diagonals" rather than
			//  location in a horizontal array.  So you always
			//  start from the main diagonal at 0 and move out.
			diagonal = j - i;
			// We change this to greater-than-one because
			//  after Q[diagonal] is reduced to one, it isn't 
			//  visited again until we reach a new streak, (because
			//  the next bit in the diagonal is a zero), and at
			//  that point we want to start with a new diagonal
			//  measure.
			if(Q[diagonal] > 1){
				y = Q[diagonal] - 1;
				Q[diagonal]--;
			}else{
				y = measureDiagonal(bg,i,j);
				Q[diagonal] = y;
			}
			while(x > currSupport){
				d = increaseMem(d,1,currSupport,currLength,
					currSupport + support*multiplier);
				currSupport += support*multiplier;
			}
			while(y > currLength){
				d = increaseMem(d,2,currSupport,currLength,
					currLength + length*multiplier);
				currLength += length*multiplier;
			}
			d[x][y]++;
			j = nextBitBitSet(bg->graph[i],j+1);
			/*
			if(x != 0){
				printf("%d:\t%d %d\n", j, x, y);
				fflush(stdout);
			}
			*/
		}
		/*
		printf("done\n");
		fflush(stdout);
		*/
	}

	// We need to rescale by the sampling factor for all i>0 in d[i][0].
	//
	for (i = 1; i < currSupport; i++) {
		d[i][0] *= s;
	}

	// Now we only need to assign the correct value for d[0][0]...
	//  but rather than figuring that out, we will just assign it in the
	//  cumulative function, since there it is merely the number of unique
	//  non-self comparisons and is easy to calculate.

	deleteBitSet(X);
	free(Q);
	*supportDim = currSupport;
	*lengthDim = currLength;
	return(d);
	
}

int
cumDMatrix(unsigned int **d, cll_t *cliqs, int currSupport, int currLength, int bgSize,
		int numSeqs) {
	int maxSup = 0;
	int maxLen = 0;
	int i, j;
	int numWins = 0;
/*	
	for (i = 0; i <= currSupport; i++) {
		printf("support = %d:\t",i);
		for (j = 0; j <= currLength; j++) {
			printf("%d\t",d[i][j]);
		}
		printf("\n");
	}
*/
	maxSup = getLargestSupport(cliqs);
	maxLen = getLargestLength(cliqs);
/********* COMMENTED OUT
	// First we note that the number of unique streaks of a given
	//  support is defined by d[support][1], where as 1 increases,
	//  the value of d decreases because only unique streaks are
	//  counted.
	// We also note that the number of disjoint node-pairs with a given
	//  number of other nodes in common is defined by d[support][0].
	// So, in order to properly account for all "unique" comparisons 
	//  (which is equal to (# streaks + # disjoint node-pairs), we must
	//  add d[support][1] to d[support][0].
	
	for (i = 0; i < currSupport + 1; i++) {
		d[i][0] += d[i][1];
	}
	********************/

	// We no longer need to do that, since now we sum across both
	//  the support and the length dimensions.  Now, d[support][0] will
	//  necessarily include d[support][1] being added to it.  We don't 
	//  want to add this anymore, otherwise we would be underestimating
	//  the probability of making that first connection.  For instance,
	//  if there were no nodes with 20 in common that weren't also 
	//  connected, and no nodes whatsoever with more than 20 in common,
	//  we'd want the p[20][0] to be 1, which would be 
	//  d[20][1]/d[20][0].  When summing across length directions,
	//  this happens naturally, whereas before we needed to do it 
	//  artificially as per above.  If we did above, we'd have the
	//  probability of each node being 1/2 instead of 1.
	
	// Rather than storing doubles and doing lots of multiplications,
	//  we're going to limit the number of operations done in the actual
	//  probability calculation by only storing cumulative sums in d.
	// Now remember, what we're storing at each location is the 
	//  number of nodes with [i] or more nodes in common (including
	//  each other and selves) that can be extended [j] times (with
	//  their initial similarity counting as 1).
	//
	// We go up to the last possible index in the length direction, which
	//  means going up to [maxLen].  We know that this is legitimate
	//  because maxLen is less than or equal to the longest possible 
	//  diagonal, and the longest possible diagonal will be less
	//  than or equal to currLength.  Since we have allotted 
	//  (currLength + 1) integers, we know we're OK to access [currLength].
	for (j = 0; j < currLength + 1; j++) {
		// We start at currSupport - 1, because currSupport will
		//  clearly not be changed, and this makes it a much easier
		//  loop to read.
		for (i = currSupport - 1; i >= 0; i--) {
			d[i][j] += d[i+1][j];
		}
	}
	
	for (i = 0; i < currSupport + 1; i++) {
		for (j = currLength - 1; j >= 0; j--) {
			d[i][j] += d[i][j+1];
		}
	}

	// Now we need to forcibly set d[0][0] to its correct value... it's 
	//  just the total number of comparisons, not including comparisons
	//  to delimiter 0's meant to separate sequences.  The number of 
	//  windows is equal to the number of offsets minus the number
	//  of sequences (assuming one delimiter per sequence).  We don't count
	//  the main diagonal, so the first row has one less, and we want to
	//  sum over all the subsequent rows in the upper half of the matrix.
	//  So it's (numWins - 1)*(numWins - 1 + 1)/2 to sum that up.
	
	numWins = bgSize - numSeqs;
	d[0][0] = numWins*(numWins - 1)/2;
/*
	for (i = 0; i <= maxSup; i++) {
		printf("support = %d:\t",i);
		for (j = 0; j <= maxLen; j++) {
			printf("%d\t",d[i][j]);
		}
		printf("\n");
	}
*/

	return 1;
}

// We also want to throw in the number of trials so that the "stat" will be
//  closer to an expected value than a probability, and so that the impact
//  of support will be lessened.  We also don't want the total (prob * trials)
//  to ever overflow or underflow, so we'll throw the steps in piecemeal to
//  try to keep the numbers reasonable.
//  Note that what we want in the end is (numberOfWindows choose Support),
//  which is [ numWindows*(numWindows-1)*...*(numWindows-support+1)/support! ]
//  Of course, this assumes independent trials, which is an approximation.

double
calcStatCliq(unsigned int **d, cll_t *cliq, int numWindows) {
	double stat = 0;
	int i = 0;
	int supChooseTwo = 0;
	double interimP = 0;
	int support = cliq->set->size;
	int length = cliq->length;
	double numTrials = 0;
	
	if (support < 2) {
		fprintf(stderr, "Support for cluster less than 2... exiting.\n");
		fflush(stderr);
		exit(0);
	}

	// OK, so support is at least two.  So we make the connections all
	//  on the first level, knowing that each node being connected has
	//  at least zero in common.  There are [(size of cluster) - 1] of
	//  these connections to be made.
	// And we know we can call for d[0][1] because if the second index
	//  were out of bounds, then there would be no similarities, and
	//  there would be no reason to call this function.
	
	interimP = ((double) d[0][1]) / ((double) d[0][0]);
	stat = pow(interimP,support-1);
	stat *= ((double) numWindows*(numWindows-1))/((double) 2);
	// Now we actually calculate the probability... the first connection
	//  has to be made no matter what, and after that we multiply for 
	//  every connection after the first one.  So we descend iteratively
	//  until we have made all connections, terminating after we've made
	//  the single i = (n - 2) connection.  There is no i = (n - 1) 
	//  connection.
	for(i = 1; i < support - 1; i++) {
		interimP = ((double) d[i][1]) / ((double) d[i][0]);
		stat *= pow(interimP,support-i-1);
		stat *= ((double) (numWindows - (i+1)))/((double) (i+2));
	}
	supChooseTwo = (support*(support - 1))/2;
	
	// Remember that length = (numwindows - 1), or alternatively,
	//  the number of extensions... normally we'd want to have the last
	//  p be p[support][numwindows - 1], which corresponds to 
	//  alteredD[support][numwindows]/alteredD[support][numwindows-1],
	//  so that means we want our last d to be d[support][numwindows].
	// Here, we note that the calculation of p's would be continuously
	//  re-normalizing, so multiplying all p's is the same as dividing
	//  the last d by the initial d.
	interimP = ((double) d[support][length+1]) / ((double) d[support][1]);
	stat *= pow(interimP,supChooseTwo);

	return stat;
}

int
calcStatAllCliqs(unsigned int **d, cll_t *allCliqs, int numWindows) {
	cll_t *curr = NULL;

	curr = allCliqs;

	while (curr != NULL) {
		curr->stat = calcStatCliq(d,curr,numWindows);
		curr = curr->next;
	}

	return(0);
}
	
int
freeD(unsigned int **d, int supportDim) {
	int i = 0;
	
	if (d == 0) {
		return 0;
	} else {
		// Still, it's supportDim + 1, because we have an extra
		//  one for the "0" support.
		for(i = 0; i < supportDim + 1; i++) {
			free(d[i]);
		}
		free(d);
		return 0;
	}
}

int
statCompare(const cll_t **first, const cll_t **second){
	double difference = (*first)->stat - (*second)->stat;

	if (difference < 0) {
		return(-1);
	} else if (difference > 0) {
		return(1);
	} else {
		return(0);
	}
}


// OK, we can't sort a linked list, so instead we'll create a data type
//  that can access each node of the linked list.
cll_t *
sortByStats(cll_t *allCliqs) {
	cll_t *curCliq = NULL;
	cll_t **arrayOfCliqs = NULL;
	int numOfCliqs = 0;
	int i = 0;

	curCliq = allCliqs;
	if (curCliq != NULL) {
		numOfCliqs = curCliq->id + 1;
	} else {
		return(NULL);
	}

	arrayOfCliqs = (cll_t **) malloc (numOfCliqs * sizeof(cll_t *));
	
	for (i = 0; i < numOfCliqs; i++) {
		arrayOfCliqs[i] = curCliq;
		curCliq = curCliq->next;
	}
	
	qsort(arrayOfCliqs,numOfCliqs,sizeof(cll_t *),statCompare);
	
	for (i = 0; i < numOfCliqs - 1; i++) {
		arrayOfCliqs[i]->next = arrayOfCliqs[i+1];
	}
	arrayOfCliqs[numOfCliqs - 1]->next = NULL;

	return(arrayOfCliqs[0]);
}
