#ifndef PDB_H
#define PDB_H

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

/* This code is based on pdbprot.h, a set
 * of Protein Data Bank I/O routines written in ANSI C,
 * on IRIX 4.0.5, on 9 Nov. 1995 by Andris Aszodi 
 *
 * http://mathbio.nimr.mrc.ac.uk/ftp/wtaylor/util/aa/incl/pdbprot.h
 *
 * It is a vastly simplified and hacked version of
 * his code.
 *
 */


/* 4-letter word for holding PDB codes */
typedef char Str4_[5]; 


typedef struct    /* entry for an atom */
{
	//  A one letter character for the amino acid this 
	//  alpha-carbon is taken from
	char Aa; 
	//  The residue number in the chain
	int Resno;  
	//  The x, y, z coordinates of the alpha carbon
	float X,Y,Z; 
} Atom_ ;

/* entry for a chain */
typedef struct 
{

	//  An array of atom objects
	Atom_ *Atoms; 
	//  The number of amino acids in this array
	int Aano;
	//  The chain ID for multichain protein structures
	char Chid;
} Chain_ ;

/* entry for a PDB record */
typedef struct	
{
	//  The 4-letter pdb code for this entry
    Str4_ Pdbcode;	/* standard PDB code "0XXX" */
	//  An array of chains in this 
    Chain_ *Chains;	
	//  Number of chains in the array of chains
    int Chainno;
} Pdbentry_ ;


//  A function to read a pdf structure
//  Takes a pointer to the name of a file,
Pdbentry_ *
get_pdb(char *filename);

int
free_pdb(Pdbentry_ *e);

#endif	/* PDB_H */
