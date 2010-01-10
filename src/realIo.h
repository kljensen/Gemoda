/*
   (c) Massachusetts of Technology, 2003
 */
#ifndef __REALIO_H_
#define __REALIO_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <gsl/gsl_matrix.h>
#include "FastaSeqIO/fastaSeqIO.h"
#include "convll.h"

typedef struct{
	int size;
	int indexSize;
	char **label;
	gsl_matrix **seq;
	int *indexToSeq;
	int *indexToPos;
	int **offsetToIndex;
} rdh_t;  // "Real data holder"  -- Lame, I know.

rdh_t
*readRealData(FILE * INPUT);

rdh_t *
freeRdh(rdh_t *data);

int
initRdhIndex(rdh_t *data, int wordSize, int seqGap);
int
getRdhIndexSeqPos(rdh_t *data, int index, int *seq, int *pos);
int
getRdhDim(rdh_t *data);
int
outputRealPats(rdh_t *data, cll_t *allPats, int L, FILE *OUTPUT_FILE,
		int **d);
int
outputRealPatsWCentroid(rdh_t *data, cll_t *allPats, int L, FILE *OUTPUT_FILE,
		double *extraParams, int compFunc);

#endif				/* __FASTA_SEQ_IO_H_ */
