/*
 *    (c) Massachusetts Institute of Technology, 2003 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "FastaSeqIO/fastaSeqIO.h"
/*#include "matrices.h"*/
#include "spat.h"
#include "bitSet.h"
#include "matdata.h"

/*  This file defines functions that are used to
 *  create our similarity graph via the comparison
 *  of small windows in our dataset.  Usually, this
 *  is via alignment (hence align.c).
 */

#define ALIGN_ALPHABET 256
/*#define ALIGN_MATRIX(i,j) idmat[i][j]*/
/*#define ALIGN_MATRIX(i,j) blosum62mt2[i][j]*/

//  aaOrder is defined int matrices.h
extern const int aaOrder[];


//  Given two pointers to strings, compare the strings
//  for a length L, using the specified scoring matrix.
int
alignMat(char *s1, char *s2, int L, int mat[][MATRIX_SIZE])
{
	int i;
	int points = 0;
	int x, y;

	//  Go over each character in the L-length window
	for (i = 0; i < L; i++) {
		//  The integer corresponding to the character in
		//  the first string, so that we can look it up 
		//  in one of our scoring matricies.
		x = aaOrder[(int) s1[i]];
		//  And for the second character
		y = aaOrder[(int) s2[i]];
		//  If the characters aren't going to be in the scoring
		//  matrix, they get a -1 value...which we'll give zero
		//  points to here.
		if (x != -1 && y != -1) {
			//  Otherwise, they get a score that is looked up
			//  in the scoring matrix
			points += mat[x][y];
		}
	}
	return points;
}

//  This uses the function above. Here, we have an
//  array of words (sPat_t objects) and we compare (align) them all.
//  If their score is above 'threshold' then we will set a bit to
//  'true' in a bitGraph_t that we create.  A bitGraph_t is essentially
//  an adjacency matrix, where each member of the matrix contains only
//  a single bit: are the words equal, true or false?
bitGraph_t *
alignWordsMat_bit(sPat_t * words, int wc, int mat[][MATRIX_SIZE], int threshold)
{
	bitGraph_t *sg = NULL;
	int score;
	int i,j;


	//  Assign a new bitGraph_t object, with (wc x wc) possible
	//  true/false values
	sg = newBitGraph(wc);
	for (i = 0; i < wc; i++) {
		for (j = i; j < wc; j++) {
			//  Get the score for the alignment of word i and word j
			score = alignMat(words[i].string, words[j].string, words[i].length, mat);
			//  If that score is greater than threshold, set
			//  a bit to 'true' in our bitGraph_t object
			if(score >= threshold){
				//  We use 'bitGraphSetTrueSym' because, if i=j,
				//  then j=i for most applications.  However, this
				//  can be relaxed for masochists.
				bitGraphSetTrueSym(sg, i, j);
			}
		}
	}

	//  Return a pointer to this new bitGraph_t object
	return sg;
}
