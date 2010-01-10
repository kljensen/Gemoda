#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "spat.h"
#include "FastaSeqIO/fastaSeqIO.h"


// Prime number generator: returns first prime number
// equal or less than n
int
sieve3(long n)
{
	int i, p, j;
	int *a;
	a = (int *) malloc((n + 1) * sizeof(int));
	if (a == NULL) {
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(0);
	}
	a[0] = 0;
	a[1] = 0;
	for (i = 2; i < n; i++) {
		a[i] = 1;
	}
	p = 2;
	do {
		j = 2 * p;
		do {
			a[j] = 0;
			j = j + p;
		} while (j <= n);
		p = p + 1;
	} while (p * p < 2 * n);
	for (i = n; i > 2; i--) {
		if (a[i]) {
			free(a);
			return i;
		}
	}
	free(a);
	return 0;
}

// A hashing function
unsigned long
hash1(unsigned char *str)
{
	unsigned long hash = 5381;
	int c;

	while ((c = *str++))
		hash = ((hash << 5) + hash) + c;	/* hash * 33 + c */

	return hash;
}


// A hashing function
int
hashpjw(char *s)
{
	char *p;
	unsigned int h, g;

	h = 0;
	for (p = s; *p != '\0'; p++) {
		h = (h << 4) + *p;
		if ((g = h & 0xF0000000)) {
			h ^= g >> 24;
			h ^= g;
		}
	}
	return h;
}



// Type for a hash table entry
typedef struct {
	char *key;
	int L;
	int data;
	int idx;
} sHashEntry_t;


// Type for a hash table
typedef struct {
	int *hashSize;
	int *iHashSize;
	int totalSize;
	sHashEntry_t **hash;
} sHash_t;

// Allocates the memory for a sHash table and
// initializes some of the elements
sHash_t
initSHash(int n)
{
	int i = 0;
	int step = 0;
	sHash_t this;

	this.totalSize = n;
	this.hashSize = (int *) malloc(n * sizeof(int));
	if (this.hashSize == NULL) {
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(0);
	}
	this.iHashSize = (int *) malloc(n * sizeof(int));
	if (this.iHashSize == NULL) {
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(0);
	}
	this.hash = (sHashEntry_t **) malloc(n * sizeof(sHashEntry_t *));
	if (this.hash == NULL) {
		fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
		fflush(stderr);
		exit(0);
	}
	for (i = 0; i < n; i++) {
		this.hash[i] = NULL;
		this.hashSize[i] = 0;
		this.iHashSize[i] = step;
	}
	return this;
}

#define SHASH_MAX_KEY_SIZE 1000

// This function has two purposes.  It searches for entries
// in the hash table and it puts new entries in.
sHashEntry_t *
searchSHash(sHashEntry_t * newEntry, sHash_t * thisHash, int create)
{
	char string[SHASH_MAX_KEY_SIZE];
	unsigned long (*hashFunction) () = &hash1;
	int i, thisIndex;
	int status = 0;

	// A string to store the key
	strncpy(string, newEntry->key, newEntry->L);
	string[newEntry->L] = '\0';

	// The index that this key hashes to
	thisIndex = hashFunction((unsigned char *) string) % thisHash->totalSize;

	// For each member that has this index, check to see
	// if the key is the same
	for (i = 0; i < thisHash->hashSize[thisIndex]; i++) {
		if (strncmp(thisHash->hash[thisIndex][i].key, string, newEntry->L) == 0) {

			// We found a match
			/*
			   printf("\t%s already in hash table!\n"); 
			 */
			status = 1;
			return &(thisHash->hash[thisIndex][i]);
			break;

		}
	}

	// If we didn't find the key and we're told to create it,
	// then allocate new memory for the hashEntry and put it in
	if (status == 0 && create != 0) {

		// Allocate space for the new entry at this index
		if(thisHash->iHashSize[thisIndex] == 0){
			thisHash->hash[thisIndex] = (sHashEntry_t *) malloc( sizeof(sHashEntry_t) );
		}else{
			thisHash->hash[thisIndex] = (sHashEntry_t *) realloc( thisHash->hash[thisIndex], 
						(thisHash->iHashSize[thisIndex] + 1) * sizeof(sHashEntry_t));
		}
		if (thisHash->hash[thisIndex] == NULL) {
			fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
			fflush(stderr);
			exit(0);
		}
		// Increase our record of the size
		i = thisHash->hashSize[thisIndex];
		thisHash->hash[thisIndex][i] = *newEntry;
		thisHash->iHashSize[thisIndex]++;
		thisHash->hashSize[thisIndex]++;

		/*
		   printf("Put %s into key %d, member %d (new size = %d)\n", string, thisIndex, i, 
		 */
		/*
		   thisHash->iHashSize[thisIndex]); 
		 */

		// Return a pointer to this entry
		return &(thisHash->hash[thisIndex][i]);
	}
	return NULL;
}


// Destroy a hash table, freeing the memory
int
destroySHash(sHash_t * thisHash)
{
	int i;
	free(thisHash->iHashSize);
	free(thisHash->hashSize);
	for (i = 0; i < thisHash->totalSize; i++) {
		if(thisHash->hash[i] != NULL){
			free(thisHash->hash[i]);
			thisHash->hash[i] = NULL;
		}
	}
	if(thisHash->hash != NULL){
		free(thisHash->hash);
		thisHash->hash = NULL;
	}
	return 0;
}


// Print the hash out
int
printSHash(sHash_t * thisHash, FILE * FH)
{
	int i, j;
	char string[SHASH_MAX_KEY_SIZE];

	for (i = 0; i < thisHash->totalSize; i++) {
		for (j = 0; j < thisHash->hashSize[i]; j++) {

			strncpy(string, thisHash->hash[i][j].key, thisHash->hash[i][j].L);
			string[thisHash->hash[i][j].L] = '\0';
			fprintf(FH, "%s %d\n", string, thisHash->hash[i][j].data);

		}
	}
	return 0;
}

int
printSPats(sPat_t * a, int n)
{
	char *s = NULL;
	int i, j;
	int size = 0;
	for (i = 0; i < n; i++) {
		if (a[i].length > size) {
			s = (char *) realloc(s, a[i].length * sizeof(char));
		}
		strncpy(s, a[i].string, a[i].length);
		s[a[i].length] = '\0';
		printf("%d:  %s\n", i, s);
		for (j = 0; j < a[i].support; j++) {
			printf("\t%d %d -> (%d, %d)\n", a[i].offset[j].seq, a[i].offset[j].pos,
				   a[i].offset[j].prev, a[i].offset[j].next);
		}
		printf("\n");
	}
	free(s);
	return 0;
}

int
destroySPatA(sPat_t * words, int wc)
{
	int i;
	for (i = 0; i < wc; i++) {
		if (words[i].offset != NULL) {
			free(words[i].offset);
			words[i].offset = NULL;
		}
	}
	free(words);
	words = NULL;
	return 0;
}

// Counts words of size L in the input FastA sequences
// and returns an array of sPat_t objects 
sPat_t *
countWords2(fSeq_t * seq, int numSeq, int L, int *numWords)
{
	int i, j;
	int totalChars = 0;
	int hashSize;
	sHashEntry_t newEntry;
	sHashEntry_t *ep;
	sHash_t wordHash;
	sPat_t *words = NULL;
	int wc = 0;
	int prev = -1;
	int l;


	// Count the total number of characters.  This
	// is the upper limit on how many words we can have
	for (i = 0; i < numSeq; i++) {
		totalChars += strlen(seq[i].seq);
	}

	// Get a prime number for the size of the hash table
	hashSize = sieve3((long) (2 * totalChars));
	wordHash = initSHash(hashSize);

	// Chop up each sequence and hash out the words of size L
	for (i = 0; i < numSeq; i++) {
		prev = -1;

		//  skip sequences that are too short to have
		//  a pattern
		if(strlen(seq[i].seq) < L){
			continue;
		}
		for (j = 0; j < strlen(seq[i].seq) - L + 1; j++) {

			// Make a hash table entry for this word
			newEntry.key = &(seq[i].seq[j]);
			newEntry.data = 1;
			newEntry.idx = wc;
			newEntry.L = L;

			// Check to see if it's already in the hash table
			ep = searchSHash(&newEntry, &wordHash, 0);
			if (ep == NULL) {

				// If it's not, create an entry for it
				ep = searchSHash(&newEntry, &wordHash, 1);

				// Increase the size of our word array
				words = (sPat_t *) realloc(words, (wc + 1) * sizeof(sPat_t));
				if (words == NULL) {
					fprintf(stderr, "Error!\n");
					fflush(stderr);
				}
				// Add the new word
				words[wc].string = &(seq[i].seq[j]);
				words[wc].length = L;
				words[wc].support = 1;
				words[wc].offset = (sOffset_t *) malloc(1 * sizeof(sOffset_t));
				if (words[wc].offset == NULL) {
					fprintf(stderr, "\nMemory Error\n%s\n", strerror(errno));
					fflush(stderr);
					exit(0);
				}
				words[wc].offset[0].seq = i;
				words[wc].offset[0].pos = j;
				words[wc].offset[0].prev = prev;
				words[wc].offset[0].next = -1;

				if (prev != -1) {
					words[prev].offset[words[prev].support - 1].next = wc;
				}
				prev = wc;
				wc++;

			} else {

				// If it is, increase the count for this word
				ep->data++;

				// add a new offset to the word array
				l = words[ep->idx].support;
				words[ep->idx].offset =
					(sOffset_t *) realloc(words[ep->idx].offset, (l + 1) * sizeof(sOffset_t));
				words[ep->idx].offset[l].seq = i;
				words[ep->idx].offset[l].pos = j;
				words[ep->idx].offset[l].prev = prev;
				words[ep->idx].offset[l].next = -1;

				// Update the next/prev
				if (prev != -1) {
					words[prev].offset[words[prev].support - 1].next = ep->idx;
				}
				prev = ep->idx;

				// Have to put this down here for cases when we create
				// a word and it is immeadiately followed by itself!!
				words[ep->idx].support += 1;
			}
		}
	}


	destroySHash(&wordHash);
	*numWords = wc;
	return words;
}

