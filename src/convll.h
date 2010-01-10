#ifndef CONVLL_H
#define CONVLL_H


#include <stdio.h>
#include <stdlib.h>
#include "bitSet.h"


typedef struct{
	int size;
	int *members;
} cSet_t;

typedef struct cnode{
	cSet_t *set;
	int id;
	//  have to use the struct as we haven't
	//  completed the typedef yet!
	int length;
	struct cnode *next;
	double stat; // Used to store the statistical store of a motif
} cll_t;

typedef struct mnode{
	int cliqueMembership;
	// must use struct because typedef is not complete yet
	struct mnode *next;
} mll_t;

// returns a pointer to a new linked list
// which has an added member
cll_t *
pushCll(cll_t * head);

// returns a pointer to a new linked list
// which has one less member (popped off top)
cll_t *
popCll(cll_t *head);


cll_t *
popAllCll(cll_t *head);

int
printCll(cll_t *head);

// Internal routine to initialize an empty linked list head
// Returns pointer to newly-modified linked list head
cll_t *
initheadCll(cll_t *head, cSet_t *newset);

// User function to push a cSet onto a linked list
// Returns pointer to a new linked list with complete head
cll_t *
pushcSet(cll_t *head, cSet_t *newset);

cll_t *
pushClique(bitSet_t *clique, cll_t *head, int* indexToSeq, int p);

mll_t * 
pushMemStack(mll_t *head, int cliqueNum);

mll_t *
popMemStack(mll_t *head);

mll_t *
popWholeMemStack(mll_t *head);

mll_t **
addToStacks(cll_t *node, mll_t **memberStacks);

mll_t **
fillMemberStacks(cll_t *head, mll_t **memberStacks);

mll_t **
emptyMemberStacks(mll_t **memberStacks, int size);

void
printMemberStacks(mll_t **memberStacks, int size);

bitSet_t *
searchMemsWithList(int *list, int listsize, mll_t **memList, int numOffsets,
			bitSet_t *queue);

bitSet_t *
setStackTrue(mll_t **memList, int i, bitSet_t *queue);

cll_t *
singleCliqueConv(cll_t *head, int firstClique, cll_t **firstGuess,int secondClique,
	cll_t **secondGuess, cll_t *nextPhase, bitSet_t *printStatus, int support);

mll_t *
mergeIntersect(cll_t *first, cll_t *second, mll_t *intersection,
		bitSet_t *printStatus, int *newSupport);

cll_t *
pushConvClique(mll_t *clique, cll_t *head);

cSet_t *
mllToCSet(mll_t *clique);

cSet_t *
bitSetToCSet(bitSet_t *clique);

cll_t *
wholeCliqueConv(cll_t *head, cll_t *node, cll_t **firstGuess, mll_t **memList, 
	int numOffsets, cll_t *nextPhase, bitSet_t *printStatus, int support);

cll_t *
wholeRoundConv(cll_t **head,mll_t **memList,int numOffsets,int support,int length,
		cll_t **allCliques);

cll_t *
completeConv(cll_t **head, int support, int numOffsets, int minLength,
		int *indexToSeq, int p);

int
printCllPattern(cll_t *node, int length);

int
uniqClique(cSet_t *clique, cll_t *head);

cll_t *
swapNodecSet(cll_t *head, int node, cSet_t *newClique);

int
yankCll(cll_t **head, cll_t *prev, cll_t **curr, cll_t **allCliques, int length);

cll_t *
removeSupers(cll_t *head, int node, cSet_t *newClique);

#endif

