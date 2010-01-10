
/*
   (c) Massachusetts of Technology, 2003 
 */
# include "fastaSeqIO.h"
# include <stdlib.h>
# include <string.h>
# include <errno.h>

#define BUFFER 100000
#define BIG_BUFFER 1000000

int
printFSeqSubSeq(fSeq_t *seq, int start, int stop){
	int i;
	for(i=start; i<stop; i++){
		putchar(seq->seq[i]);
	}
	return 0;
}

// Measure a line in a file, returning the file to
// the caller at the original location.
long
measureLine(FILE * INPUT)
{
	long start;
	long count = 0;
	int myChar;
	start = ftell(INPUT);
	myChar = fgetc(INPUT);
	count++;
	while (myChar != '\n' && myChar != EOF) {
		count++;
		myChar = fgetc(INPUT);
	}
	fseek(INPUT, start, SEEK_SET);
	return count;
}


// Count the number of lines that start with a '>'
long
CountFSeqs(FILE * INPUT)
{
	long start;
	long count = 0;
	int myChar;
	int newLine = 1;
	start = ftell(INPUT);
	myChar = fgetc(INPUT);
	while (myChar != EOF) {
		if (newLine == 1 && myChar == '>') {
			count++;
		}
		if (myChar == '\n') {
			newLine = 1;
		} else {
			newLine = 0;
		}
		myChar = fgetc(INPUT);
	}
	fseek(INPUT, start, SEEK_SET);
	return count;
}

// Count the number of lines 
long
countLines(FILE * INPUT)
{
	long start;
	long count = 1;
	int myChar;
	int status = 0;
	start = ftell(INPUT);
	myChar = fgetc(INPUT);
	while (myChar != EOF) {
		if (myChar == '\n') {
			count++;
			status = 1;
		} else {
			status = 0;
		}
		myChar = fgetc(INPUT);
	}
	if (status == 1) {
		count--;
	}
	fseek(INPUT, start, SEEK_SET);
	return count;
}

int
initAofFSeqs(fSeq_t * aos, int numSeq)
{
	int i;
	for (i = 0; i < numSeq; i++) {
		aos[i].seq = NULL;
		aos[i].label = NULL;
	}
	return 1;
}

char **
ReadFile(FILE * INPUT, int *n)
{
	char **buf = NULL;
	long nl;
	long tls = 0;
	int i=0;

	nl = countLines(INPUT);
	if( nl == 0){
		fprintf(stderr, "\nNo sequences! Error!\n\n");
		fflush(stderr);
		return NULL;
	}
	buf = (char **) malloc ( (int)(nl+1) * sizeof(char *));
	if ( buf == NULL){
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(0);
	}

	//  measure the first line
	tls = measureLine(INPUT) + 1;
	if(tls != 0){
		buf[i] = (char *) malloc ( tls * sizeof(char));
		if ( buf[i] == NULL){
			fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
			fflush(stderr);
			exit(0);
		}
	}
	fgets(buf[i], tls, INPUT);
	do{
		if(buf[i][ strlen(buf[i])-1 ] == '\n'){
			buf[i][ strlen(buf[i])-1 ] = '\0';
		}
		tls = measureLine(INPUT) + 1;
		if(tls != 0){
			i++;
			buf[i] = (char *) malloc ( tls * sizeof(char) );
			if ( buf[i] == NULL){
				fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
				fflush(stderr);
				exit(0);
			}
		}
	}while( fgets(buf[i], tls, INPUT) != NULL );
	free(buf[i]);
	buf = (char **) realloc ( buf, i * sizeof(char *) );
	if ( buf == NULL){
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		return NULL;
	}
	//  I think that 'i' might actually be the # of lines
	//  plus one here?  somehow line 131 isn't being freed,
	//  or at least 2 bytes of it.
	*n = i;
	return buf;
}


fSeq_t *
ReadTxtSeqs(FILE * INPUT, int *numberOfSequences){
	int i;
	int nl;
	char **buf = NULL;
	fSeq_t *aos;

	buf = ReadFile(INPUT, &nl);
	if(buf == NULL){
		return NULL;
	}
	aos = (fSeq_t *) malloc ( nl * sizeof(fSeq_t));
	if( aos == NULL){
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(0);
	}
	initAofFSeqs(aos, nl);
	for ( i=0 ; i<nl ; i++ ){
		aos[i].seq = buf[i];
	}
	free(buf);
	*numberOfSequences = nl;
	return (aos);
}


fSeq_t *
ReadFSeqs(FILE * INPUT, int *numberOfSequences){
	int i,j,k;
	int nl, ns=0;
	char **buf = NULL;
	fSeq_t *aos;
	sSize_t *ss;
	sSize_t *ll;

	buf = ReadFile(INPUT, &nl);
	if(buf == NULL){
		return NULL;
	}

	//  Count how many sequences we have
	for( j=0 ; j<nl ; j++){
		//  Test for null lines
		if(strlen(buf[j]) >= 1){
			if(buf[j][0] == '>'){
				ns++;
			}
		}
	}
	ss = (sSize_t *) malloc ( ns * sizeof(sSize_t) );
	if(ss == NULL){
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(0);
	}
	ll = (sSize_t *) malloc ( ns * sizeof(sSize_t) );
	if(ll == NULL){
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(0);
	}

	//  find the first sequence
	k=0;
	while( strlen(buf[k]) < 1 || buf[k][0] != '>'){
		k++;
	}

	//  record how large each sequence is
	i = -1;
	for( j=k ; j<nl ; j++){
		if(strlen(buf[j]) >= 1){
			if(buf[j][0] == '>'){
				i++;
				ll[i].start = j;
				ll[i].stop = j;
				ll[i].size = strlen( buf[j] );;
				ss[i].start = j+1;
				ss[i].size = 0;
			}else{
				ss[i].stop = j;
				ss[i].size += strlen( buf[j] );;
			}
		}
	}

	aos = (fSeq_t *) malloc ( ns * sizeof(fSeq_t));
	if( aos == NULL){
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(0);
	}
	initAofFSeqs(aos, ns);

	for ( i=0 ; i<ns ; i++ ){
		if( ll[i].size > 0 ){
			aos[i].label = (char *) malloc ( (ll[i].size+1) * sizeof(char) );
			if( aos[i].label == NULL){
				fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
				fflush(stderr);
				exit(0);
			}
			aos[i].label[0] = '\0';
			for ( j=ll[i].start ; j<=ll[i].stop ; j++ ){
				/*printf("(%d): label.size = %d ; buf.size = %d\n", i,ll[i].size, strlen(buf[j]));*/

				//  both instances of strcat here are using
				//  .label/.seq's that are NULL and that is
				//  throwing a memory error in valgrind
				aos[i].label = strcat ( aos[i].label, buf[j] );
			}
		}
		if( ss[i].size > 0 ){
			aos[i].seq = (char *) malloc ( (ss[i].size+1) * sizeof(char) );
			if( aos[i].seq == NULL){
				fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
				fflush(stderr);
				exit(0);
			}
			aos[i].seq[0] = '\0';
			for ( j=ss[i].start ; j<=ss[i].stop ; j++ ){
				/*printf("(%d): seq.size = %d ; buf.size = %d\n", i,ss[i].size, strlen(buf[j] ) );*/
				aos[i].seq = strcat ( aos[i].seq, buf[j] );
			}
		}
	}
	free(ll);
	free(ss);

	for ( i=0 ; i<nl ; i++ ){
		free(buf[i]);
	}
	free(buf);

	*numberOfSequences = ns;
	return aos;
}

int
FreeFSeqs(fSeq_t * arrayOfSequences, int numberOfSequences)
{
    int i;
    for (i = 0; i < numberOfSequences; i++) {
        if (arrayOfSequences[i].label != NULL) {
            free(arrayOfSequences[i].label);
        }
        arrayOfSequences[i].label = NULL;

        if (arrayOfSequences[i].seq != NULL) {
            free(arrayOfSequences[i].seq);
        }
        arrayOfSequences[i].seq = NULL;
    }
    if (arrayOfSequences != NULL) {
        free(arrayOfSequences);
    }
    arrayOfSequences = NULL;
    return EXIT_SUCCESS;
}

// write an array of fastaSequence objects to a file 
// starts from 'start' and stops at 'stop'.  trying
// to print a sequence that isn't there will segfault
// 
int
WriteFSeqA(FILE * MY_FILE, fSeq_t * arrayOfSequences, int start, int stop)
{
    int i;
    for (i = start; i <= stop; i++) {
        fprintf(MY_FILE, "%s\n", arrayOfSequences[i].label);
        fprintf(MY_FILE, "%s\n", arrayOfSequences[i].seq);
    }
    return EXIT_SUCCESS;
}
