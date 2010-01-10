#include "pdb.h"

/*  All of the functions in this file are prototyped in pdb.h
 */


//  Read a line from a file and make sure it is non-null.  If
//  the line is null (ie, we readed EOF) return EXIT_FAILURE
int
pdb_get_nonnull_line(char *line,int bufferSize, FILE *FH){
	if(fgets(line, bufferSize, FH) == NULL){
		fprintf(stderr, "Error, file is poorly formatted!\n");
		exit(EXIT_FAILURE);
	}
	return 0;
}

//  Read a pdb file and return an array of Pdbentry objects,
//  For most 
Pdbentry_ * get_pdb(char *filename){
	Pdbentry_ *myPdb = NULL;
	int bufferSize = 500;
	char line[bufferSize];
	FILE *FH = NULL;
	int i,j;
	char c, lc;
	Atom_ *a;

	//  Try to open the file, or exit with failure.
	if ((FH = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Couldn't open file %s; %s\n", filename, strerror(errno));
		exit(EXIT_FAILURE);
	}

	//  Allocate enough space for a Pdbentry_ object or
	//  exit with failure
	myPdb = (Pdbentry_ *) malloc (sizeof(Pdbentry_));
	if(myPdb == NULL){
		fprintf(stderr, "Cannot allocate memory!\n");
		exit(EXIT_FAILURE);
	}

	// scan the pdb id
	// this line must be non-null
	pdb_get_nonnull_line(line, bufferSize, FH);
	sscanf(line, "%4c", myPdb->Pdbcode);
	myPdb->Pdbcode[4] = '\0';

	// scan the number of chains
	// this line must be non-null
	pdb_get_nonnull_line(line, bufferSize, FH);
	sscanf(line, "%d", &(myPdb->Chainno));

	//  Allocate enough space for this chain, based on the
	//  size we read in, or die
	myPdb->Chains = (Chain_ *)malloc(myPdb->Chainno * sizeof(Chain_));
	if(myPdb->Chains == NULL){
		fprintf(stderr, "Cannot allocate memory!\n");
		exit(EXIT_FAILURE);
	}

	//  Read each x,y,z coordinate of the chain
	for( i=0 ; i<myPdb->Chainno ; i++){

		//  For each Chain,scan in the number of 
		//  amino acids that are in this chain
		pdb_get_nonnull_line(line, bufferSize, FH);
		sscanf(line, "%d", &(myPdb->Chains[i].Aano));

		//  Allocate enough space for this atom or die
		myPdb->Chains[i].Atoms = (Atom_ *)malloc(myPdb->Chains[i].Aano * sizeof(Atom_));
		if(myPdb->Chains == NULL){
			fprintf(stderr, "Cannot allocate memory!\n");
			exit(EXIT_FAILURE);
		}
	}
	//  Loop over each chain
	for( i=0 ; i<myPdb->Chainno ; i++){

		// And loop over each atom in this chain
		for( j=0 ; j<myPdb->Chains[i].Aano ; j++){

			//  Read in this atom's coordinates
			pdb_get_nonnull_line(line, bufferSize, FH);

			a = &(myPdb->Chains[i].Atoms[j]);
			sscanf(line, "%c %c %f %f %f", &(c), &(a->Aa), &(a->X), &(a->Y), &(a->Z));

			if(c != lc && j!=0){
				fprintf(stderr, "Error, poorly formatted file!\n");
				exit(EXIT_FAILURE);
			}
			lc = c;
		}
	}

	//  The following lines can be printed out for debugging
	/*
	for( i=0 ; i<myPdb->Chainno ; i++){
		printf("Chain %d\n", i);
		for( j=0 ; j<myPdb->Chains[i].Aano ; j++){
			a = &(myPdb->Chains[i].Atoms[j]);
			printf("(%c) %0.4f %0.4f %0.4f\n", a->Aa, a->X, a->Y, a->Z);
		}
	}

	*/

	return myPdb;
}

//  This frees up the memory for a pointer to a
//  Pdbentry_ object and sets the pointer to NULL
int
free_pdb(Pdbentry_ *e){
	int i;

	for( i=0 ; i<e->Chainno ; i++){
		if(e->Chains[i].Atoms != NULL){
			//  Free each atom in this chain
			free(e->Chains[i].Atoms);
			e->Chains[i].Atoms = NULL;
		}
	}
	//  Free each chain
	if(e->Chains != NULL){
		free(e->Chains);
		e->Chains = NULL;
	}
	return(EXIT_SUCCESS);
}
