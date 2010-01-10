#include "protAlign.h"
/************************************************************************
 * Note: the theory behind the optimal alignment calculations was taken fromi
 * the following paper:
 * Arun, KS, TS Huang, & SD Blostein "Least-Squares Fitting of Two 3-D 
 *  Point Sets." IEEE Transactions on Pattern Analysis and Machine Intelligence.
 *  VOL. PAMI-9, No 5, Sept 1987.
 *  ***********************************************************************/
/*************************************************************************
 * Master function for aligning the Pdbentries.  Note that the array of 
 *  alpha-carbon points has already been created, and that's what we'll use
 *  to do our manipulations.  This takes all of the alpha-carbon points and
 *  runs a window of length <L> along them, uses our centroid and SVD method
 *  to calculate the optimal translation and rotation, calculates the 
 *  RMSD for the two windows, and assigns a true bit in the master bitgraph
 *  if the RMSD is less than the threshold rmsd.
 * Input: an array of alpha-carbon points, an array of the Pdbentry files that
 *  contain the coordinates of the trace, the window length, the RMSD 
 *  threshold, and the number of windows as returned by getCaPoint.
 * Output: a bitGraph containing the similarity matrix of each window
 *  compared to each other window based on minimum RMSD and the given
 *  threshold.
 * ***********************************************************************/

bitGraph_t *
alignPdbs(caPoint_t *cap, Pdbentry_ **pdbs, int L, double g, int numWindows){
	int i,j;
	bitGraph_t *bg = newBitGraph(numWindows);
	gsl_matrix *a = gsl_matrix_alloc(L, 3);
	gsl_matrix *b = gsl_matrix_alloc(L, 3);

	for(i=0 ; i<numWindows ; i++){
		//  skip boundary entries
		if(cap[i].pdbNo == -1){
			continue;
		}
		pdb2mat(pdbs[cap[i].pdbNo], cap[i].chainNo, cap[i].caNo, cap[i].caNo+L, a);
		moveToCentroid(a);
		for(j=i+1 ; j<numWindows ; j++){
			if(cap[j].pdbNo == -1){
				continue;
			}
			pdb2mat(pdbs[cap[j].pdbNo], cap[j].chainNo, cap[j].caNo, cap[j].caNo+L, b);
			moveToCentroid(b);
			rotateMats(a,b);
			/*printf ("rmsd =\t%f (g = %lf)\n", gsl_matrix_rmsd(a,b), g);*/
			if(gsl_matrix_rmsd(a,b) <= g){
				bitGraphSetTrueSym(bg,i,j);
			}
		}
	}


	gsl_matrix_free(a);
	gsl_matrix_free(b);
	return bg;
}

/*************************************************************************
 * Function to get pointers to translate from the (pdb file number, chain
 *  number, residue number) triplet to a single array index.
 * Input: pointer to an array of Pdbentry pointers, the number of pdbfiles 
 *  given as input, the window length L, a pointer to the numnber of windows,
 *  and a pointer to the number of chains.
 * Output: pointer to an array of alpha-carbon points (also altered values
 *  for the number of windows and number of chains in the main function).
 * ************************************************************************/
// THIS ASSUMES pdbs HAS **ONLY** Ca ATOMS!!!
caPoint_t *getCaPoint(Pdbentry_ **pdbs, int num, int window, int *numWindows, int *numChains){
	int i,j,k,l;
	int totNumCa=0;
	caPoint_t *cap = NULL;
	*numChains = 0;

	// for each pdb file
	for( i=0 ; i<num ; i++){
		// for each chain in each pdb file
		for( j=0 ; j<pdbs[i]->Chainno ; j++){
			// Get the number of windows that we'll have.
			totNumCa += pdbs[i]->Chains[j].Aano - window + 1;
			// could change Aano to Atomno if we're
			// just reading in the Ca atoms...
			*numChains +=1;
			// Add one null ca between chains! --- this will
			//  make us able to see that a chain has ended and
			//  so we can't perform manipulations on part of a
			//  chain
			totNumCa++;
		}
	}

	// Now we'll make an array of caPoint structures.
	cap = (caPoint_t *) malloc (totNumCa * sizeof(caPoint_t));
	if (cap == NULL){
		fprintf(stderr, "Couldn't allocate mem!; %s\n", strerror(errno));
		exit(EXIT_FAILURE);
	}
	// Let the parent function know the number of windows.
	*numWindows = totNumCa;
	//  fill in cap
	l=0;
	for( i=0 ; i<num ; i++){
		for( j=0 ; j<pdbs[i]->Chainno ; j++){
			for( k=0 ; k<pdbs[i]->Chains[j].Aano - window +1 ; k++){
				cap[l].pdbNo = i;
				cap[l].chainNo = j;
				cap[l].caNo = k;
				l++;
			}
			// Don't forget -- add one null ca between chains!
			cap[l].pdbNo = -1;
			cap[l].chainNo = -1;
			cap[l].caNo = -1;
			l++;
		}
	}
	printf("L = %d ::: numWindows = %d\n", l, *numWindows);
	fflush(stdout);
	return cap;
}

/************************************************************************
 * Calculates the centroid of the gsl_matrix and moves the structure to
 *  have its centroid be at the origin.
 * Input: a gsl_matrix with the coordinates of two protein windows.
 * Output: Integer EXIT_SUCCESS value.
 * **********************************************************************/

int moveToCentroid(gsl_matrix *m){
	int i;
	gsl_vector *v = gsl_vector_calloc (m->size2);
	gsl_vector_view t;

	for(i=0; i<m->size1 ; i++){
		//  t is a row slice of m
		t = gsl_matrix_row(m, i);

		//  v is the sum of all rows
		gsl_vector_add(v, &t.vector);
	}
	//  v is now the average x,y,z
	gsl_vector_scale(v,  1/(double)m->size1);

	//  subtract off v from each row of m
	for(i=0; i<m->size1 ; i++){
		t = gsl_matrix_row(m, i);
		gsl_vector_sub(&t.vector, v);
	}
	gsl_vector_free(v);
	return EXIT_SUCCESS;
}

/*************************************************************************
 * Function to take a chain from a reduced-PDB file and create a gsl_matrix
 *  from it.
 * Input: a pdbentry structure that has the chain which is to be made into
 *  a matrix, the number of the desired chain, the number of the first
 *  residue to be put into the matrix, the number of the last residue to
 *  be put into the matrix, and the gsl_matrix into which it will be added.
 * Output: integer EXIT_SUCCESS value (and an altered gsl_matrix).
 * **********************************************************************/
/*gsl_matrix **/
int
pdb2mat(Pdbentry_ *s, int chainNum, int resStart, int resStop, gsl_matrix *m){
	Atom_ *a = NULL;
	int i,j;

	//  check for reasonable inputs
	if(		s == NULL ||
			s->Chains == NULL ||
			chainNum >= s->Chainno ||
			resStart >= resStop ||
			resStart < 0 ||
			resStop > s->Chains[chainNum].Aano){
		fprintf(stderr, "Bad pdb chain request! %d (%d->%d)", chainNum, resStart, resStop);
		fflush(stderr);
		exit(0);
	}
	if(m->size1 != resStop-resStart || m->size2 != 3){
		fprintf(stderr, "Bad matrix for pdb2mat!");
		fflush(stderr);
		exit(0);
	}

	//  a is a pointer to the array of atoms in this chain
	a = s->Chains[chainNum].Atoms;

	// Now pull the values out of the Pdbentry file.
	for( i=resStart, j=0 ; i<resStop; i++, j++){
		gsl_matrix_set(m, j, 0, (double)a[i].X);
		gsl_matrix_set(m, j, 1, (double)a[i].Y);
		gsl_matrix_set(m, j, 2, (double)a[i].Z);
	}
	/*return m;*/
	return EXIT_SUCCESS;
}

/***********************************************************************
 * Calculates the outer product of the coordinates of each location in the
 *  two windows and sums them together for the correlation matrix.
 * Input: two gsl_matrices with the coordinates of <windowlength> alpha-
 *  carbons.
 * Output: a 3x3 gsl_matrix with the sum of the outer products.
 * **********************************************************************/
gsl_matrix *
sumOuterProd(gsl_matrix *a, gsl_matrix *b){
	int i,j,k;
	double newVal;
	gsl_matrix *h;
	// Make sure the matrices are the same size
	if(a->size1 != b->size1 || a->size2 != b->size2){
		fprintf(stderr, "Bad matrices for correlation matrix!");
		fflush(stderr);
		exit(0);
	}
	// Create a matrix to store the result in.
	h = gsl_matrix_calloc(a->size2, a->size2);

	// Now sum each outer product into the result matrix.
	for( i=0 ; i<a->size1 ; i++){
		for( j=0 ; j<a->size2 ; j++){
			for( k=0 ; k<a->size2 ; k++){
				newVal = gsl_matrix_get(h,j,k);
				newVal += gsl_matrix_get(a,i,j) * gsl_matrix_get(b,i,k);
				gsl_matrix_set(h, j, k, newVal);
			}
		}
	}
	return h;
}

/**********************************************************************
 * An attractive way to print out the contents of a gsl_matrix.
 * Input: a file pointer for output, a gsl_matrix to be printed, and a string
 *  indicating the format of the printing.
 * Output: Integer success value (and the printing of the gsl_matrix to the
 *  file pointed to by fp).
 *  8*********************************************************************/

int
gsl_matrix_pretty_fprintf(FILE *fp, const gsl_matrix *a, const char *format){
	int i,j;
	for( i=0 ; i<a->size1 ; i++){
		for( j=0 ; j<a->size2 ; j++){
			fprintf(fp, format, gsl_matrix_get(a,i,j));
		}
		fprintf(fp,"\n");
	}
	return EXIT_SUCCESS;
}

/**************************************************************************
 * Finds the matrix rotation to produce the optimal RMSD for two protein 
 *  windows that have their centroids at the origin.
 * Input: two gsl_matrices with alpha-carbon windows that we are
 *  looking to rotate.
 * Output: integer EXIT_SUCCESS (and a changed matrix b to have the optimal
 *  rotation).
 * ***********************************************************************/
int
rotateMats(gsl_matrix *a, gsl_matrix *b){
	gsl_matrix *h = NULL;
	gsl_vector *S = NULL;
	gsl_matrix *V = NULL;
	gsl_matrix *R = NULL;
	gsl_vector *work = NULL;
	gsl_matrix *result = NULL;
	CBLAS_TRANSPOSE_t noT = CblasNoTrans;
	CBLAS_TRANSPOSE_t yesT = CblasTrans;

	// find the correlation matrix h
	h = sumOuterProd(a,b);

	S = gsl_vector_alloc(h->size2);
	V = gsl_matrix_alloc(h->size2, h->size2);
	work = gsl_vector_alloc(h->size2);

	//  svd		
	//  h = U S V^T
	//
	//  h is replaced with U on exit
	//
	gsl_linalg_SV_decomp(h, V, S, work);

	R = gsl_matrix_calloc(h->size1, h->size1);
	gsl_blas_dgemm(noT,yesT, 1.0, V, h, 1.0, R);


	//  the final result matrix
	result = gsl_matrix_calloc(b->size1, h->size2);

	gsl_blas_dgemm(noT, noT, 1.0, b, R, 0.0, result);

	gsl_matrix_swap(b, result);

	gsl_vector_free(S);
	gsl_vector_free(work);
	gsl_matrix_free(V);
	gsl_matrix_free(h);
	gsl_matrix_free(R);
	gsl_matrix_free(result);
	return EXIT_SUCCESS;
}

/***********************************************************************
 * Computes the rmsd between two matrices.
 * Input: two gsl_matrices with alpha-carbon traces.
 * Output: a double with the value of the RMSD.
 * **********************************************************************/
double gsl_matrix_rmsd(gsl_matrix *a, gsl_matrix *b){
	int i,j;
	double total_rmsd=0;
	gsl_matrix *c = gsl_matrix_alloc(a->size1, a->size2);

	//  Check for reasonable inputs
	if(a->size1 != b->size1 || a->size2 != b->size2){
		fprintf(stderr, "Bad matrices for correlation matrix!");
		fflush(stderr);
		exit(0);
	}
	gsl_matrix_memcpy(c, a);

	gsl_matrix_sub(c, b);

	gsl_matrix_mul_elements(c, c);

	// Find the sum of the square deviations.
	for (i=0 ; i<c->size1 ; i++){
		for (j=0 ; j<c->size2 ; j++){
			total_rmsd += gsl_matrix_get(c, i, j);
		}
	}

	// Now get the rmsd.
	total_rmsd = sqrt( total_rmsd / c->size1 );
	gsl_matrix_free(c);
	return total_rmsd;
}

/***********************************************************************
 * Outputs the protein motifs in any of a number of different formats,
 *  including pdb and an xml format useful with the gemoda-p pattern viewer
 *  (gpv).
 * Input: an array of alpha-carbon points, the array of initial Pdbentries
 *  containing all of the chains and coordinates, a linked list of all of the
 *  motifs to be output, the length of the smallest possible motif, an integer
 *  indicating how many unique structures a motif must be in to be returned, 
 *  a filehandle indicating where the output should be directed, and an 
 *  integer indicating the desired format of the output.
 * Output: integer success value (and output to the filehandle FH).
 * ********************************************************************/
int
outputProtPats(caPoint_t *cap, Pdbentry_ **pdbs, cll_t *allPats, int L, int us, FILE *FH, int format){
	int i,j, k, l;
	gsl_matrix *a;
	gsl_matrix *b;
	cll_t *curPat;
	int len, sup;
	int offset1, offset2;
	int pdbNo, chainNo, caNo;
	int check_us;
	int last;

	curPat = allPats;

	if(format == 0){
		fprintf(FH, "HEADER    gemoda-p output, many patterns in one file\n");
	}else{
		fprintf(FH, "<gemodap>\n");
	}
	i=0;
	// As long as there's another motif on the stack...
	while(curPat != NULL){
		len = curPat->length+L;
		sup = curPat->set->size;

		//  check for unique structure support
		if(us >= 1){
			check_us = 0;
			last = -100;
			for(j=0; j<sup && check_us<us ; j++){
				offset1 = curPat->set->members[j];
				pdbNo = cap[offset1].pdbNo;
				caNo = cap[offset1].caNo;
				if(pdbNo != last){
					check_us++;
					last = pdbNo;
				}
			}
			if(check_us<us){
				curPat = curPat->next;
				fflush(stdout);
				continue;
			}
		}
		
		// So it has enough unique structure support, so start the
		//   output.
		if(format == 0){
			fprintf(FH, "HEADER    gemoda-p output, many patterns in one file\n");
			fprintf(FH, "REMARK    pattern id=%d length=%d support=%d\n", i, len, sup);
		}else{
			fprintf(FH, "\t<pattern id=\"%d\" length=\"%d\" support=\"%d\">\n", i, len, sup);
		}

		// Make a couple of matrices to store the data in.
		a = gsl_matrix_alloc(len, 3);
		b = gsl_matrix_alloc(len, 3);

		//  we assume there's at least one member!
		//  we'll align all the rest to this when we print out!
		offset1 = curPat->set->members[0];

		pdbNo = cap[offset1].pdbNo;
		chainNo = cap[offset1].chainNo;
		caNo = cap[offset1].caNo;

		// Make the window into a matrix, and zero-centroid it.
		pdb2mat(pdbs[pdbNo], chainNo, caNo, caNo+len, a);
		moveToCentroid(a);

		// Now output it.
		if(format == 0){
			fprintf(FH, "REMARK    pattern id=%d length=%d support=%d\n", i, len, sup);
			fprintf(FH, "REMARK    instance struct=%s chain=%d pos=%d\n",
						pdbs[pdbNo]->Pdbcode, chainNo, caNo);
			fprintf(FH, "REMARK    sequence=");
			for(k=caNo ; k<caNo+len ; k++){
				fprintf(FH, "%c", pdbs[pdbNo]->Chains[chainNo].Atoms[k].Aa);
			}
			fprintf(FH, "\n");
		}else{
			fprintf(FH, "\t\t<instance struct=\"%s\" chain=\"%d\" pos=\"%d\">\n",
						pdbs[pdbNo]->Pdbcode, chainNo, caNo);
			fprintf(FH, "\t\t\t<seq>");
			for(k=caNo ; k<caNo+len ; k++){
				fprintf(FH, "%c", pdbs[pdbNo]->Chains[chainNo].Atoms[k].Aa);
			}
			fprintf(FH, "</seq>\n");
			fprintf(FH, "\t\t\t<coords>");
		}
		for(k=0 ; k<a->size1 ; k++){
			if(k!=0){
				fprintf(FH, " ");
			}
			for(l=0 ; l<a->size2 ; l++){
				if(l!=0){
					fprintf(FH, " ");
				}
				fprintf(FH, "%.3lf", gsl_matrix_get(a, k, l));
			}
		}
		fprintf(FH, "</coords>\n");
		fprintf(FH, "\t\t</instance>\n");

		//  now we have the first offset, let's align the
		//  other offsets to it using our svd algorithm
		//  and output them all so that they're aligned with the first
		for(j=1 ; j<sup ; j++){
			offset2 = curPat->set->members[j];
			pdbNo = cap[offset2].pdbNo;
			chainNo = cap[offset2].chainNo;
			caNo = cap[offset2].caNo;

			pdb2mat(pdbs[pdbNo], chainNo, caNo, caNo+len, b);

			moveToCentroid(b);
			rotateMats(a,b);

			// print out the 'b' instance
			fprintf(FH, "\t\t<instance struct=\"%s\" chain=\"%d\" pos=\"%d\">\n",
						pdbs[pdbNo]->Pdbcode, chainNo, caNo);

			fprintf(FH, "\t\t\t<seq>");
			for(k=caNo ; k<caNo+len ; k++){
				fprintf(FH, "%c", pdbs[pdbNo]->Chains[chainNo].Atoms[k].Aa);
			}
			fprintf(FH, "</seq>\n");

			fprintf(FH, "\t\t\t<coords>");
			for(k=0 ; k<b->size1 ; k++){
				if(k!=0){
					fprintf(FH, " ");
				}
				for(l=0 ; l<b->size2 ; l++){
					if(l!=0){
						fprintf(FH, " ");
					}
					fprintf(FH, "%.3f", gsl_matrix_get(b, k, l));
				}
			}
			fprintf(FH, "</coords>\n");
			fprintf(FH, "\t\t</instance>\n");

		}
		fprintf(FH, "\t</pattern>\n");
		curPat = curPat->next;
		i++;
		gsl_matrix_free(a);
		gsl_matrix_free(b);
	}
	fprintf(FH, "</gemodap>\n");
	return 0;
}
