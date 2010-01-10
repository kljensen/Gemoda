
#include "realCompare.h"

//  Calculate the rmsd between two windows, with optional
//  translation and rotation.
double
rmsdCompare(rdh_t *data, int win1, int win2, int L, double *extraParams){
	int trans = 1;
	int rot = 1;
	int dim;
	double result=0;
	int seq1,pos1;
	int seq2,pos2;
	gsl_matrix_view view1;
	gsl_matrix_view view2;
	gsl_matrix *mat1;
	gsl_matrix *mat2;
	gsl_matrix *mat1copy;
	gsl_matrix *mat2copy;

	// The "rint" function is in math.h and rounds a number to the 
	//  nearest integer.  It raises an "inexact exception" if the
	//  number initially wasn't an integer.
	if (extraParams != NULL) {
		trans = rint(extraParams[0]);
		rot = rint(extraParams[1]);
	}
	dim = getRdhDim(data);

	//  Find out which seq,pos pairs these two
	//  windows correspond to
	getRdhIndexSeqPos(data, win1, &seq1, &pos1);
	getRdhIndexSeqPos(data, win2, &seq2, &pos2);

	//  Get a reference to a submatrix.  That is,
	//  'chop out' the window.
	view1 = gsl_matrix_submatrix(data->seq[seq1], pos1, 0, L, dim);
	view2 = gsl_matrix_submatrix(data->seq[seq2], pos2, 0, L, dim);

	//  This just makes it easier to handle the views
	mat1 = &view1.matrix;
	mat2 = &view2.matrix;

	//  Create copies of the windows, because our comparison
	//  will require altering the matrices
	mat1copy = gsl_matrix_alloc(mat1->size1, mat1->size2);
	mat2copy = gsl_matrix_alloc(mat2->size1, mat2->size2);


	gsl_matrix_memcpy(mat1copy, mat1);
	gsl_matrix_memcpy(mat2copy, mat2);



	/*
	printf("matrix1:\n");
	gsl_matrix_pretty_fprintf(stdout, mat1copy, "%f ");
	printf("\nmatrix2:\n");
	gsl_matrix_pretty_fprintf(stdout, mat2copy, "%f ");
	*/

	//  Are we going to do a translation?
	if(trans == 1){
		moveToCentroid(mat1copy);
		moveToCentroid(mat2copy);
	}
	//  Are we going to do a rotation?
	if(rot == 1){
		//  Rotate mat2copy to have a minimal
		//  rmsd with mat1copy
		rotateMats(mat1copy,mat2copy);
	}

	//  Compute the rmsd between mat2copy and mat2copy
	result = gsl_matrix_rmsd(mat1copy, mat2copy);

	gsl_matrix_free(mat1copy);
	gsl_matrix_free(mat2copy);
	return result;
}
double
generalMatchFactor(rdh_t *data, int win1, int win2, int L, double *extraParams){
	int i,j;
	double numerator=0.0;
	/*double denominator=0.0;*/
	double xsum;
	double ysum;
	double ldenom=0.0;
	double rdenom=0.0;
	int dim;
	int seq1,pos1;
	int seq2,pos2;
	gsl_matrix_view view1;
	gsl_matrix_view view2;
	gsl_matrix *mat1;
	gsl_matrix *mat2;

	dim = getRdhDim(data);

	//  Find out which seq,pos pairs these two
	//  windows correspond to
	getRdhIndexSeqPos(data, win1, &seq1, &pos1);
	getRdhIndexSeqPos(data, win2, &seq2, &pos2);

	//  Get a reference to a submatrix.  That is,
	//  'chop out' the window.
	view1 = gsl_matrix_submatrix(data->seq[seq1], pos1, 0, L, dim);
	view2 = gsl_matrix_submatrix(data->seq[seq2], pos2, 0, L, dim);

	//  Some error checking here would be nice!
	//  Did we get the matrices we wanted?

	//  This just makes it easier to handle the views
	mat1 = &view1.matrix;
	mat2 = &view2.matrix;

	//  Loop over each position
	for ( i=0 ; i<mat1->size1 ; i++){
		xsum=0.0;
		ysum=0.0;

		//  Loop over each dimension at each position
		for ( j=0 ; j<dim ; j++){
			xsum += gsl_matrix_get(mat1, i, j);
			ysum += gsl_matrix_get(mat2, i, j);
		}

		numerator += (i+1) * sqrt(xsum * ysum);
		ldenom += (i+1) * xsum;
		rdenom += (i+1) * ysum;
	}

	return pow(numerator, 2.0)/(ldenom*rdenom);
}
double
massSpecCompareWElut(rdh_t *data, int win1, int win2, int L, double *extraParams){
	int i,j;
	double numerator=0.0;
	/*double denominator=0.0;*/
	double xsum;
	double ysum;
	double cum;
	double ldenom=0.0;
	double rdenom=0.0;
	int dim;
	int seq1,pos1;
	int seq2,pos2;
	double weight = 2.0;
	gsl_matrix_view view1;
	gsl_matrix_view view2;
	gsl_matrix *mat1;
	gsl_matrix *mat2;
	double maxElut = -1;
	
	if (extraParams != NULL) {
		maxElut = extraParams[0];
	}

	dim = getRdhDim(data);

	//  Find out which seq,pos pairs these two
	//  windows correspond to
	getRdhIndexSeqPos(data, win1, &seq1, &pos1);
	getRdhIndexSeqPos(data, win2, &seq2, &pos2);

	//  Get a reference to a submatrix.  That is,
	//  'chop out' the window.
	view1 = gsl_matrix_submatrix(data->seq[seq1], pos1, 0, L, dim);
	view2 = gsl_matrix_submatrix(data->seq[seq2], pos2, 0, L, dim);

	//  Some error checking here would be nice!
	//  Did we get the matrices we wanted?

	//  This just makes it easier to handle the views
	mat1 = &view1.matrix;
	mat2 = &view2.matrix;
	cum = 1.0;

	//  Loop over each position
	for ( i=0 ; i<mat1->size1 ; i++){
		xsum=0.0;
		ysum=0.0;

		// First take the first dimension for elution time
		if (maxElut >= 0) {
			if(fabs(gsl_matrix_get(mat1,i,0)-
					gsl_matrix_get(mat2,i,0)) > maxElut) {
				cum = 0;
				break;
			}
		}
//		printf("\n");
//
		//  Loop over each subsequent dimension at each position
		for ( j=1 ; j<dim ; j++){
//			printf("mat1val=%lf,mat2val=%lf\n",gsl_matrix_get(mat1,i,j),
//					gsl_matrix_get(mat2,i,j));
			numerator += pow(j,weight) * sqrt(gsl_matrix_get(mat1, i, j)
						* gsl_matrix_get(mat2,i,j));
			ldenom += pow(j,weight) * gsl_matrix_get(mat1, i, j);
			rdenom += pow(j,weight) * gsl_matrix_get(mat2, i, j);
//			printf("numer=%lf,ldenom=%lf,rdenom=%lf\n",numerator,
//					ldenom,rdenom);
		}

		cum *= pow(numerator, 2.0)/(ldenom*rdenom);
	}

	return pow(cum,1.0/L);
}

double (*getCompFunc(int compFunc))(rdh_t*,int,int,int,double*) {
	double (*comparisonFunc) (rdh_t*,int,int,int,double*) = &rmsdCompare;
	
	switch (compFunc) {
		case 0:
			comparisonFunc = &rmsdCompare;
			break;
		case 1:
			comparisonFunc = &generalMatchFactor;
			break;
		case 2:
			comparisonFunc = &massSpecCompareWElut;
			break;
		default:
			comparisonFunc = &rmsdCompare;
			break;
	}

	return(comparisonFunc);
}

bitGraph_t *
realComparison(rdh_t *data, int L, double g, int compFunc, double* extraParams){
	int i,j;
	int seq1,pos1;
	int seq2,pos2;
	bitGraph_t *bg = NULL;
	double score;
	double (*comparisonFunc) (rdh_t*,int,int,int,double*) = &rmsdCompare;

	//  Initialize the rdh's index
	initRdhIndex(data, L, 1);

	//  Allocate a new bit graph
	bg = newBitGraph(data->indexSize);

	// Choose the comparison function, pass a reference to it
	
	comparisonFunc = getCompFunc(compFunc);
		
	for( i=0 ; i<data->indexSize ; i++){

		//  Skip seperators
		getRdhIndexSeqPos(data, i, &seq1, &pos1);
		if(seq1 == -1 || pos1 == -1){
			continue;
		}

		for( j=i ; j<data->indexSize ; j++){
			getRdhIndexSeqPos(data, j, &seq2, &pos2);
			if(seq2 == -1 || pos2 == -1){
				continue;
			}

			//  This is the comparison function
			score = comparisonFunc(data, i, j, L, extraParams);
//			printf("score (%2d,%2d) vs. (%2d, %2d) =\t%lf\n",seq1, pos1, seq2, pos2, score);
			if (compFunc == 0) {
				if(score <= g) {
					bitGraphSetTrueSym(bg,i,j);
				}
			} else if ((compFunc == 1) || (compFunc == 2)) {
				if(score >=g){
					/*printf("**");*/
					bitGraphSetTrueSym(bg,i,j);
				}
			} else {
				fprintf(stderr,"Comparison function undefined in "
							"realComparison function,\n located in "
							"realCompare.c.  Exiting.\n\n");
				fflush(stderr);
				exit(0);
			}
			/*printf("\n");*/
		}
	}

	return bg;
}

