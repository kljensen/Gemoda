
#include <stdio.h>
#include <string.h>
#include "matdata.h"
#include "matrixmap.h"

#define DEFAULT_MATRIX blosum62
/*#define DEFAULT_MATRIX idmat*/
/*#define DEFAULT_MATRIX blosum62mt2*/
/*#define DEFAULT_MATRIX dna_idmat*/

/***********************************************************************
 * A simple function to take the matrix name argument given as input to
 *  gemoda and return the physical memory location of that matrix by using the
 *  matrix_map construct.
 * Input: a string containing the matrix name a pointer to a two-dimensional
 *  array.
 * Output: None, though the value of the pointer given as input is changed
 *  to reflect the location of the matrix
 * ********************************************************************/

void
getMatrixByName(char name[], const int (**matp)[MATRIX_SIZE]) {
	int i;
	for (i = 0; matrix_map[i].name != NULL; i++) {
		if (strcmp(name, matrix_map[i].name) == 0) {
			break;
		}
	}
	if (matrix_map[i].name != NULL) {
		*matp = (matrix_map[i].mat);
	} else {
		*matp = (DEFAULT_MATRIX);
	}
}
