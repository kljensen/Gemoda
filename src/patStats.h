/*
   (c) Massachusetts of Technology, 2003
 */
#ifndef __PATSTATS_H_
#define __PATSTATS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "bitSet.h"
#include "convll.h"
#include <time.h>

unsigned int**
getStatMat(bitGraph_t *bg, int support, int length, int *supportDim,
		int *lengthDim, int numBlanks, int s, FILE *OUTPUT_FILE);

int
cumDMatrix(unsigned int **d, cll_t *cliqs, int currSupport, int currLength,
		int bgSize, int numSeqs);

int
calcStatAllCliqs(unsigned int **d, cll_t *allCliqs,int numWindows);

cll_t*
sortByStats(cll_t *allCliqs);

int
freeD(unsigned int **d, int supportDim);
#endif				
