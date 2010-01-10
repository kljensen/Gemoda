/*
   (c) Massachusetts of Technology, 2003
 */
#ifndef __PROT_ALIGN_H_
#define __PROT_ALIGN_H_

/*
 *  This file contains the functions neccessary for aligning
 *  small windows of protein structures and is used exclusively
 *  by the gemoda-p program.
 *
 *
 *  The basic concept here is that, gemoda-p reads in the x,y,z
 *  coordinates of alpha-carbon atoms in a set of proteins.
 *  To compare two windows of length L, it first makes each of 
 *  their centroids be zero, then uses an SVD technique that 
 *  minimizes the RMSD between the windows.  The resulting RMSD
 *  is then compared to the threshold, and if it qualifies then
 *  this similarity is indicated in the bitgraph.
 *
 */

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include "pdb.h"
#include "bitSet.h"
#include "convll.h"



typedef struct {
	int pdbNo;
	int chainNo;
	int caNo;
} caPoint_t;

caPoint_t *getCaPoint(Pdbentry_ **pdbs, int num, int window, int *numWindows, int *numChains);

bitGraph_t *
alignPdbs(caPoint_t *cap, Pdbentry_ **pdbs, int L, double g, int numWindows);


int moveToCentroid(gsl_matrix *m);

/*gsl_matrix **/
int
pdb2mat(Pdbentry_ *s, int chainNum, int resStart, int resStop, gsl_matrix *m);

gsl_matrix *
sumOuterProd(gsl_matrix *a, gsl_matrix *b);

int
gsl_matrix_pretty_fprintf(FILE *fp, const gsl_matrix *a, const char *format);

int
rotateMats(gsl_matrix *a, gsl_matrix *b);

double
gsl_matrix_rmsd(gsl_matrix *a, gsl_matrix *b);

int
outputProtPats(caPoint_t *cap, Pdbentry_ **pdbs, cll_t *allPats,int L, int us, FILE *FH, int format);

#endif /*__PROT_ALIGN_H_*/
