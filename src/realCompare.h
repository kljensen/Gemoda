/*
   (c) Massachusetts of Technology, 2003
 */
#ifndef __REALCOMPARE_H_
#define __REALCOMPARE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <gsl/gsl_matrix.h>
#include "realIo.h"
#include "bitSet.h"
#include "protAlign.h"  // for moveToCentroid and for #include <math.h>

double
rmsdCompare(rdh_t *data, int win1, int win2, int L, double *extraParams);

double
generalMatchFactor(rdh_t *data, int win1, int win2, int L, double *extraParams);

double
massSpecCompareWElut(rdh_t *data, int win1, int win2, int L, double *extraParams);

bitGraph_t *
realComparison(rdh_t *data, int l, double g, int compFunc, double* extraParams);

double (*getCompFunc(int compFunc))(rdh_t*,int,int,int,double*);

#endif				
