#include <errno.h>
#include <string.h>
#include "convll.h"
#include "bitSet.h"
/*#include "efence.h"*/

extern cll_t * pruneCll(cll_t *head, int *indexToSeq, int p);

/************************************************************
 Pushes a new, empty head onto a linked list of cliques
Note: this should always be followed by a call to initheadCll, as the
  head pushed on here is empty and will be meaningless without any members.
  This function should NOT be used by the user; see pushcSet.
 Input: head of a linked list
 Output: head of a linked list
*************************************************************/
cll_t *
pushCll(cll_t * head)
{
	//Make a pointer, verify memory
	cll_t *a = NULL;
	a = (cll_t *) malloc(sizeof(cll_t));
	if (a == NULL) {
		fprintf(stderr, "\nMemory Error - pushCll\n%s\n", strerror(errno));
		fflush(stderr);
		exit(0);
	}
	// Initialize id (sequential) and pointer to next item, but not
	//  the cSet with the clique members
	if(head == NULL){
		a->id = 0;
		a->next = NULL;
	}else{
		a->next = head;
		a->id = head->id + 1;
	}
	a->set = NULL;
	a->length = -1;
	a->stat = -1;
	return a;
}



/***********************************************************************
 Removes the head of the clique linked list, returns the new head of the
     clique linked list, and frees the memory occupied by the old head.
Input: head of a linked list
Output: head of a linked list
*************************************************************************/

cll_t *
popCll(cll_t *head){
	// by default the new head is NULL...is important later
	cll_t *newHead = NULL;
	if (head == NULL) {
		fprintf(stderr, "\nCan't pop a null linked list\n");
		fflush(stderr);
		exit(0);
	}

	// unless this is the end of the linked list, set the new head
	//  to the next member of the list.  Otherwise, since by default the
	//  new head is NULL, it will properly return an empty list
	if( head->next != NULL){
		newHead = head->next;
	}

	// Check to see if there is a set.  If there is, and there are members,
	//  then first free the members.  And if there is a set, then free it.
	if(head->set != NULL) {
		if (head->set->members != NULL) {
			free(head->set->members);
			head->set->members = NULL;
		}
		free(head->set);
		head->set = NULL;
	}

	// Both the members and set have been freed, so now can free the cll_t
	//  without leaking anything.

	free(head);
	head = NULL;
	return newHead;
}

/********************************************************************
 * Shortcut function to pop all of the members of a linked list
 * Input: head of a linked list
 * Output: head of a now-empty linked list
********************************************************************/
cll_t *
popAllCll(cll_t *head){
	while(head != NULL){
		head = popCll(head);
	}
	return head;
}

/*******************************************************************
 * Prints the members (cliques) of a linked list in the format:
 * 	id = <unique id number of clique within linked list>
 * 	Length = <number of members of clique, if available>
 * 	Size = <length of each member of clique>
 * 	Members = <newline-separated list of members of the clique>
 * Input: head of a linked list
 * Output: Gives text output, returns (meaningless) exit value
 *************************************************************************/
int
printCll(cll_t *head){
	int i = 0;	
	cll_t *curr = head;
	while(curr != NULL){
		printf("id = %d\n", curr->id);
		// Make sure the clique is nonzero in size before attempting
		// to print it
		if ((curr->set != NULL) && (curr->set->size > 0)) {
			if (curr->length >= 0) {
				printf("Length = %d\n", curr->length);
			}
			printf("Size = %d\n", curr->set->size);
			printf("Members = \n");
			for (i = 0; i < curr->set->size; i++) {
				printf("\t%d\n", curr->set->members[i]);
			}
			printf("***********************************************\n");
		}
		else {
			fprintf(stderr, "\nClique has no members! -- printCll\n");
			fflush(stderr);
			exit(0);
		}
		curr = curr->next;
	}
	return EXIT_SUCCESS;
}

/*******************************************************************
 * Initializes the empty head of a linked list by adding a set to that head.
 *   Note: this is only called immediately after pushing onto a cll, because
 *   the push always creates a new empty head.
 *   This function should NOT be called by the user; see pushcSet.
 * Input: head of a linked list, pointer to a cSet_t list of clique members.
 * Output: head of a linked list
 * ******************************************************************/

cll_t *
initheadCll(cll_t *head, cSet_t *newset){
	// Check to make sure that the head is not already initialized.
	if (head->set != NULL) {
		printf("Stack head already initialized!");
		exit(0);
	}
	// Make the head's set pointer point to the new set.
	head->set = newset;
	return head;
}

/******************************************************************
 * Function that pushes the contents of a cSet (set of members
 * of a clique) onto a linked list of cliques.
 * Input: head of a linked list, new clique in the form of a cSet
 * Output: head of a linked list
 * ******************************************************************/

cll_t *
pushcSet(cll_t *head, cSet_t *newset){
	head = pushCll(head);
	head = initheadCll(head, newset);
	return head;
}


/***************************************************************
 * Converts a bitSet_t to a cSet_t for the purposes of pushing it onto
 * a linked list of cliques.  The bitSet_t data structure is used for
 * massive comparisons during clique-finding but is unwieldy/inefficient
 * when it is known that the structure is sparse.  The cSet_t allows for
 * efficient comparison of sparse bitSet_t's.
 * Use this just before pushing a newly-discovered clique onto a clique
 * linked list.
 * Input: a new clique in the form of a bitSet_t
 * Output: the same clique in the form of a cSet_t
 * ***************************************************************/

cSet_t *
bitSetToCSet(bitSet_t *clique){
	int cliqueSize = countSet(clique);
	int i = 0, start = 0;
	cSet_t *holder = (cSet_t *) malloc(sizeof(cSet_t));

	// Memory error checking
        if (holder == NULL) {
                fprintf(stderr, "\nMemory Error - bitSetToCSet - [1]\n%s\n", 
			strerror(errno));
                fflush(stderr);
                exit(0);
        }

	// More memory checking
	holder->members = (int *) malloc(cliqueSize*sizeof(int));
        if (holder->members == NULL) {
                fprintf(stderr,"\nMemory Error - bitSetToCSet - [2]\n%s\n", 
			strerror(errno));
                fflush(stderr);
                exit(0);
        }


	// For each member of the clique in the bitSet,
	for(i = 0; i < cliqueSize; i++){
		// Find the next one, add its location to the members array
		holder->members[i] = nextBitBitSet(clique, start);
		// (But check for errors... if we get to the end of the
		// bitSet, then something is wrong)
		if (holder->members[i] == -1) {
			fprintf(stderr,"\nClique error - not enough members\n");
			fflush(stderr);
			exit(0);
		}
		// Increment to move on in the nextBitBitSet search
		start = holder->members[i] + 1;
	}

	holder->size = cliqueSize;
	return holder;
}


/******************************************************************
 * Checks to enforce the -p flag (minimum number of unique input sequences
 *   in which the motif occurs).
 * Input: a clique in the form of a cSet_t, pointer to the index/sequence
 *   number data structure, the -p flag value
 * Output: An integer: 1 for success, 0 for failure
 * ******************************************************************/
int
checkCliquecSet(cSet_t *cliquecSet, int* indexToSeq, int p){
	int *seqNums = NULL;
	int thisSeq = 0, i = 0, j = 0;
	seqNums = (int *) malloc(p * sizeof(int));

	if (seqNums == NULL) {
		fprintf(stderr,"Memory error - checkCliquecSet\n%s\n",
				strerror(errno));
		fflush(stderr);
		exit(0);
	}

	// Initialize an array of integers of size p to sentinel values of -1
	for (i = 0; i < p; i++) {
		seqNums[i] = -1;
	}
	j = 0;

	if (cliquecSet->size < 1) {
		fprintf(stderr, "\nClique of zero size! - checkCliquecSet\n");
		fflush(stderr);
		exit(0);
	}
	// Find the first sequence number.
	seqNums[0] = indexToSeq[cliquecSet->members[0]];
	// Iterate over the remaining size of the clique
	for (i = 1; i < cliquecSet->size; i++) {
		// Find the next sequence number.
		thisSeq = indexToSeq[cliquecSet->members[i]];
		// The member list is in monotonic order, so we only need
		//  to compare the current member to the previous member to 
		//  find out if it comes from the same sequence.
		// If it's not from the same sequence, increment the unique
		//  sequence counter (j), store the next sequence number
		//  in the array.
		// Also check to see if we've already reached the p threshold,
		//  and if so, then bail out.
		if (thisSeq != seqNums[j]) {
			j++;
			seqNums[j] = thisSeq;
			if (j == p - 1) {
				break;
			}
		}
	}

	// Now just see what the value of the last number in the array is;
	//  if it's the sentinel, then we didn't find instances in p 
	//  unique sequences.  If it's not the sentinel, then we've met
	//  the -p criterion.
	if (seqNums[p - 1] == -1) {
		free(seqNums);
		return(0);
	} else {
		free(seqNums);
		return(1);
	}
}

/***********************************************************************
 * Pushes a bitSet onto a clique linked list, performing all necessary 
 *  manipulations in order to do so.
 * Input: new clique in the form of a bitSet_t, head of a linked list,
 *  pointer to the index/sequence number data structure, integer value of
 *  the -p flag
 * Output: head of an updated clique linked list
 * *********************************************************************/
// Yes, this variable is meant to be global; it's used for debugging.
int cliquecounter = 0;

cll_t *
pushClique(bitSet_t *clique, cll_t *head, int* indexToSeq, int p){
	cSet_t *cliquecSet = NULL;

	// Change the bitSet_t to a cSet_t
	cliquecSet = bitSetToCSet(clique);
	// If the -p flag has been assigned a value, then check the clique
	//  and only proceed if that criterion is met.  Otherwise, free the
	//  memory that we had allocated up to this point.
	if (p > 1) {
		if (checkCliquecSet(cliquecSet, indexToSeq, p)) {
			cliquecounter++;
			/*printf("%d\n",cliquecounter);*/
			/*fflush(stdout);*/
			head = pushcSet(head, cliquecSet);
		} else {
			free(cliquecSet->members);
			free(cliquecSet);
		}
	// If the -p flag wasn't set, then just push the cSet onto the linked
	//  list.
	} else {
		cliquecounter++;
		/*printf("%d\n",cliquecounter);*/
		/*fflush(stdout);*/
		head = pushcSet(head,cliquecSet);
	}
        return head;
}

/************************************************************************
 * This begins code for the member linked lists.  A single one of these
 *  linked lists functions somewhat similarly to the clique linked lists, 
 *  though with less information stored.
 * Functionally, an array of member linked lists is used to access the 
 *  "inverse" of what is contained in the clique linked lists.  That is,
 *  we would like to be able to look up the cliques that a given node is a
 *  member of, so we have an array of member linked lists of size equal
 *  to the number of nodes.
 * ***********************************************************************/

/*************************************************************************
 * Pushes a single clique membership onto a node's member stack.
 * Input: the head of a single member linked list, a clique number to be added
 * Output: the head of a single member linked list
 * *********************************************************************/

mll_t*
pushMemStack(mll_t *head, int cliqueNum){
        mll_t *a = NULL;
        a = (mll_t *) malloc(sizeof(mll_t));
	// Memory error checking
        if (a == NULL) {
                fprintf(stderr,"\nMemory Error - pushMemStack: %s\n", strerror(errno));
                fflush(stderr);
                exit(0);
        }
        if(head == NULL){
                a->next = NULL;
        }else{
                a->next = head;
        }
	// Store the number of the clique of which the node is a member.
	//  Note that we assume no duplication, which is guaranteed
	//  by our method of filling the member stacks, which is quite simple:
	//  go through all members of a clique (which have no duplicates
	//  because they are constructed from merge-intersections or from
	//  bitSet_t's) and add that clique to each node's membership list.
	a->cliqueMembership = cliqueNum;
        return a;
}

/**********************************************************************
 * Pops the head off of a single member linked list.
 * Input: head of a member linked list
 * Output: the new head of a member linked list after popping one item
 * *********************************************************************/

mll_t*
popMemStack(mll_t *head){
        // by default the new head is NULL...is important later
        mll_t *newHead = NULL;
        if (head == NULL) {
                fprintf(stderr, "\nCan't pop a null linked list - popMemStack\n");
                fflush(stderr);
                exit(0);
        }
        if(head->next != NULL){
                newHead = head->next;
        }
        free(head);
        head = NULL;
	return newHead;
}

/********************************************************************
 * Pops all items off of a member linked list.
 * Input: head of a member linked list
 * Output: empty head of a member linked list
 * ******************************************************************/
mll_t *
popWholeMemStack(mll_t *head){
        while(head != NULL){
                head = popMemStack(head);
        }
        return head;
}
                                                                                

/*****************************************************************
 * For one clique, it adds membership for that clique to all of its members'
 *  member stacks.
 * Input: a specific clique in a clique linked list, an array of member stacks
 * Output: the array of updated member stacks
 * ************************************************************************/

mll_t **
addToStacks(cll_t *node, mll_t **memberStacks){
	int i = 0;
	int cliqueNum = 0;

	// Make sure that we don't reference NULL values
	if (node->set != NULL) {
		// Go through each member of the clique's set
		for (i = 0; i < node->set->size; i++) {
			// Get the member's number
			cliqueNum = node->set->members[i];
			// Go to that member's linked list and push
			//  on the number of the current clique
			memberStacks[cliqueNum] = 
				pushMemStack(memberStacks[cliqueNum],node->id);
		}
	} else {
		fprintf(stderr, "\nNULL set for clique! - addToStacks\n");
		fflush(stderr);
		exit(0);
	}
	return memberStacks;
}

/***********************************************************************
 * Fills the entire memberStacks data structure by calling addToStacks
 *  for each clique in the clique linked list.
 * Input: head of a clique linked list, array of member linked lists
 * Output: the array of updated member linked lists
 * *********************************************************************/

mll_t **
fillMemberStacks(cll_t *head, mll_t **memberStacks){
	cll_t *curr = head;
	// Just go down the linked list calling addToStacks
	while (curr != NULL) {
		memberStacks = addToStacks(curr, memberStacks);
		curr = curr->next;
	}

	return memberStacks;
}

/***********************************************************************
 * After we have performed a round of convolution, this "empties" the 
 *  member stacks by popping all nodes off each member linked list.
 * Input: array of member linked lists, the size of that array (total number
 *  of offsets)
 * Output: the array of now-empty member linked lists
 * **********************************************************************/

mll_t **
emptyMemberStacks(mll_t **memberStacks, int size){
	int i = 0;

	for (i = 0; i < size; i++){
		memberStacks[i] = popWholeMemStack(memberStacks[i]);
	}

	return memberStacks;
}

/**********************************************************************
 * Prints the contents of the member stacks in the format:
 *     Offset <offset number>: <comma-separated list of cliques>
 * Input: array of member linked lists, size of that array (total number of
 *  offsets)
 * Output: only text output/no return value
 * **********************************************************************/

void
printMemberStacks(mll_t **memberStacks, int size){
	int i = 0;
	mll_t *curr = NULL;

	for (i = 0; i < size; i++){
		curr = memberStacks[i];
		printf("Offset %d: ", i);
		while (curr != NULL) {
			printf("%d,",curr->cliqueMembership);
			curr = curr->next;
		}
		printf("\n");
	}
}

/*********************************************************************
 * Adds all of the members of a given stack to a "queue" in the form of 
 *  a bitSet_t data structure.  That is, for each clique in the member linked
 *  list, it sets the corresponding bit in the bitSet_t true.
 * Input: array of member linked lists, an integer indicating a specific
 *  member linked list, and a bitSet_t of length >= the number of cliques
 *  in the current clique linked list
 * Ouput: the updated bitSet_t
 * *********************************************************************/

bitSet_t *
setStackTrue(mll_t **memList, int i, bitSet_t *queue) {
        mll_t *curr = memList[i];

	// Traverse down the member linked list
        while (curr != NULL) {
		// Set the bit in queue corresponding to the current clique
		//  membership true
                setTrue(queue, curr->cliqueMembership);
                curr = curr->next;
        }

        return queue;
}

/**************************************************************************
 * Creates one large queue by calling "setStackTrue" for each member of a
 *  list of offsets.  This then creates the union of clique membership
 *  for all offsets in the list being searched.
 * Input: an array of offset numbers, the length of that array, an array
 *  of member linked lists, the length of that array (the total number of
 *  offsets), and a bitSet_t to store the union/queue.
 * Output: the union/queue in a bitSet_t structure
 * ************************************************************************/

bitSet_t *
searchMemsWithList(int *list, int listsize, mll_t **memList, int numOffsets,
			bitSet_t *queue) {
	int i = 0;
	emptySet(queue);

	// Go through each offset in the list
	for (i = 0; i < listsize; i++) {
		// Check to make sure that's a valid offset number, and if so
		//  then set its stack true in the queue.
		if (list[i] + 1 < numOffsets) {
			queue = setStackTrue(memList, list[i] + 1, queue);
		} else {
			fprintf(stderr, "\nInvalid offset number! - searchMemsWithList\n");
			fprintf(stderr, "\nlist[i]+1 (%d) >= numOffsets (%d)\n", list[i] + 1, numOffsets);
			fflush(stderr);
			exit(0);
		}
	}
	
	return queue;
}

/***********************************************************************
 * Convolves one single clique against one other single clique.  Note that
 *  this is non-commutative, so exchanging firstClique and secondClique
 *  will not give the same results.  The "guess" pointers keep the location
 *  of the previous clique in the linked list so that we don't have to search
 *  the linked list from the beginning/end every time.  We exploit our 
 *  earlier tidiness in that we can reasonably guess that we will monotonically
 *  traverse down cliques.
 * Input: head of the current clique linked list, the id number of the 
 *  first clique, a pointer to a guess at the first clique, the id number 
 *  of the second clique, a pointer to a guess at the second clique, the 
 *  head of the clique linked list for the next round of convolution, a 
 *  bitSet indicating which cliques should be output as maximal, and the
 *  minimum support flag
 * Output: the head of clique linked list for the next round of convolution
 *  (which may have changed if the two cliques could be convolved)
 * *************************************************************************/

cll_t *
singleCliqueConv(cll_t *head, int firstClique, cll_t **firstGuess, int secondClique,
	cll_t **secondGuess, cll_t *nextPhase, bitSet_t *printStatus, int support){
	cll_t *first = NULL, *second = NULL;
	mll_t *survivingMems = NULL;
//	int flag = 0;
	int newSupport = 0;
//	cll_t *checker = head;		

	// Check to make sure we're looking for legitimate cliques.
	if ((firstClique > head->id) || (secondClique > head->id)) {
		fprintf(stderr, "\nNonexistent clique! - singleCliqueConv\n");
		fflush(stderr);
		exit(0);
	}

	// Our guesses depend on monotonic traversal.  If we don't find
	//  the first clique, then bail out.
	while((*firstGuess)->id != firstClique) {
		if ((*firstGuess)->next != NULL) {
			*firstGuess = (*firstGuess)->next;
		} else {
			fprintf(stderr, "\nFirst clique not found! - singleCliqueConv\n");
			fflush(stderr);
			exit(0);
		}
	}
	first = *firstGuess;

	// Our guesses depend on monotonic traversal.  If we don't find
	//  the second clique, then bail out.
	while((*secondGuess)->id != secondClique) {
		if ((*secondGuess)->next != NULL) {
			(*secondGuess) = (*secondGuess)->next;
		} else {
			fprintf(stderr, "\nSecond clique not found! - singleCliqueConv\n");
			fflush(stderr);
			exit(0);
		}
	}
	second = *secondGuess;
	// Find out what the surviving members are when the first clique
	//  is convolved with the second clique
	survivingMems = mergeIntersect(first,second,survivingMems,printStatus,
					&newSupport);

	// If the first clique is subsumed by the second, then it is not
	//  maximal, so don't print it.
        // printStatus true means print it!
        if (newSupport == first->set->size) {
                setFalse(printStatus,first->id);
	}
	// If the second clique is subsumed by the first, then it is not
	//  maximal, so don't print it.
	if (newSupport == second->set->size) {
                setFalse(printStatus,second->id);
        }


	// If the support of the clique just formed by convolution meets the
	//  support criterion, then push it on to the linked list for
	//  the next phase of convolution.
	if (newSupport >= support) {
//		printf("Push %d and %d\n",first->id,second->id);
		nextPhase = pushConvClique(survivingMems, nextPhase);
//		printf("---------\n");
//		printCll(nextPhase);
//		printf("---------\n");
	}

	// Pop the surviving members; they are no longer needed, as they
	//  either didn't meet the support criterion or have been pushed on
	//  already
	survivingMems = popWholeMemStack(survivingMems);

	return nextPhase;
}

/*************************************************************************
 * Convolves two cliques in a non-commutative manner.  It finds which members
 *  of the first clique are immediately followed by a member in the second
 *  clique.
 * Input: pointer to the location in the linked list of the first clique to
 *  be convolved, pointer to the location in the linked list of the second
 *  clique to be convolved, a member linked list used to store the intersection
 *  of the two cliques, the printstatus bitSet, and a pointer to an integer
 *  with the support of the clique formed by convolution.
 * Output: a member linked list with the intersection of the two cliques,
 *  plus the side effect of that intersection's cardinality being stored
 *  in the integer pointed to by newSupport.
 * ************************************************************************/

mll_t *
mergeIntersect(cll_t *first, cll_t *second, mll_t *intersection,
		bitSet_t *printstatus, int *newSupport) {

	int i = 0, j = 0, status = 0;
	
	// Make sure we are still in-bounds, otherwise we bail out
	// We'll refer to the offset currently being analyzed from the 
	//  first clique as the 'first offset' and the offset currently
	//  being analyzed from the second clique as the 'second offset'
	while((i < first->set->size) && (j < second->set->size)) {
		// If the second offset is earlier than the first offset plus
		//  one, then we move on to the next possible second offset
		if ((first->set->members[i] + 1) > second->set->members[j]) {
			j++;
		}
		// If the second offset is later than the first offset plus
		//  one, then we move on the next possible first offset
		else if ((first->set->members[i] + 1) < second->set->members[j]) {
			i++;
		}
		// Otherwise, the second offset is equal to the first offset
		//  plus one, so we have an extendable node.  Push that on
		//  to the intersection stack, move both the first and second
		//  offsets to their respective next possible offsets, and 
		//  increment the support counter for the new clique (status)
		else {
			intersection = pushMemStack(intersection,
							first->set->members[i]);
			i++;
			j++;
			status++;
		}
	}
	
	// Send the value of the clique's new support out of this function
	*newSupport = status;
	return intersection;
}

/***********************************************************************
 * Before we push a convolved clique onto the stack for the next level,
 *  this function ensures that it is not subsumed by and does not subsume any 
 *  other clique currently on that stack.
 * Input: a candidate clique for the next level in cSet_t form, and the head
 *  of the clique linked list for the next level
 * Output: an integer indicating the status of the proposed clique with 
 *  respect to the next level: -1 if the clique is unique, -2 if the clique
 *  is a subset/duplicate of an existing clique, or a clique id in the range
 *  [0,numcliques) representing the first clique of which the proposed one
 *  is a superset.
 * Note that by executing this each time a clique is added to the next level,
 *  we ensure that if the new clique is not unique, it can only be a superset
 *  or a subset of some other clique; it CANNOT be both a strictly superset of 
 *  one and a strictly subset of another.  One of those other two cliques 
 *  would have been identified in previous steps as being super- or sub-sets,
 *  so it is impossible for one clique now to be both a super and a subset.
 *  ***********************************************************************/
int
uniqClique(cSet_t *cliquecSet, cll_t *head) {
	int i = 0, j = 0;
	int asubbflag = 1, bsubaflag = 1;
    
	// Descend through all members of the next level's linked list
	while(head != NULL) {
		asubbflag = 1;
		bsubaflag = 1;
		i = 0;
		j = 0;
		// The proposed convolved clique will be referred to as the
		//  "first" clique, and the current clique being analyzed
		//  in the next level is the "second" clique.
		// Continue if we have more members in both cliques AND if it
		//  is still possible for one clique to be a subset of 
		//  the other.
		while((i < cliquecSet->size) && (j < head->set->size) &&
			((asubbflag == 1) || (bsubaflag == 1))) {
			// If the current member of the first clique is less
			//  than the current member of the second clique,
			//  it is impossible for the first clique to be a 
			//  subset of the second (since the members are
			//  traversed in ascending order.
			if (cliquecSet->members[i] < head->set->members[j]) {
				i++;
				asubbflag = 0;
			}
			// Similarly, if the current member of the second
			//  clique is less than the current member of the 
			//  second clique, the second can't be a subset
			//  of the first.
			else if (cliquecSet->members[i] > head->set->members[j]) {
				j++;
				bsubaflag = 0;
			}
			// Otherwise, they matched this time, so move them
			//  both on.
			else {
				i++;
				j++;
			}
		}

		// If the proposed clique is a subset of some other clique
		//  in the next level, then return -2, and it won't be added.
		//  (Note, this also is how exact duplicates are handled.)
		if ((asubbflag == 1) && (i == cliquecSet->size)) {
			return(-2);
		}
		// If the proposed clique is a superset of some other clique(s)
		//  in the next level, then return the id of the first clique
		//  of which it is a superset.  
		if ((bsubaflag == 1) && (j == head->set->size)) {
			return(head->id);
		}
		// If the proposed clique has not been found to be a superset
		//  or a subset yet, then move on to the next clique in
		//  the next level.
		head = head->next;
	}
	// If we've gotten here, we've checked all cliques in the previous
	//  level and haven't found the proposed clique to be a superset or
	//  a subset... if so, then we're all good, so return a -1.
	return(-1);
}
	
/*******************************************************************
 * Swaps out a node in a linked list that has been found to be a subset
 *  of a node that is not yet in the list.
 * Input: the head of a clique linked list, a specific node within that
 *  linked list that is to be removed, and the new clique that is the superset
 *  of the node to be removed (in cSet_t form)
 * Output: the head of the altered clique linked list
 * *********************************************************************/
cll_t *
swapNodecSet(cll_t *head, int node, cSet_t *newClique) {
	int foundflag = 0;	
	cll_t *curr = head;

	// First we find the node that needs to be swapped out
	while(curr != NULL) {
		if (curr->id == node) {
			foundflag = 1;
			break;
		}
		curr = curr->next;
	}

	// If we can't find it, then we get upset and exit.
	if (foundflag == 0) {
		fprintf(stderr, "\nClique not found! (in swapNode)\n");
		fflush(stderr);
		exit(0);
	}
	
	// Then we free the useless clique's members and its set data structure
	//  before pointing its set to the new clique.
	free(curr->set->members);
	free(curr->set);
	curr->set = newClique;
	return head;

}

/************************************************************************
 * This function finds all cliques in a linked list of which the proposed
 *  clique is a superset.  It starts looking AFTER the first clique which
 *  has already been found to be a subset.  In some senses, it is just a 
 *  continuation of the uniqclique function in order to take advantage
 *  of the fact that though a proposed clique can only be a subset of one
 *  existing next-level clique, it can be a superset of many existing next-
 *  level cliques.
 * Input: head of a clique linked list, the id of the first node found
 *  to be a subset of the proposed clique, and the proposed clique (in cSet_t
 *  form)
 * Output: the head of the clique linked list with all but the first subset
 *  (which was passed as an argument) removed.  This function is now ready
 *  for swapNode to be called.
 *  **********************************************************************/

cll_t *
removeSupers(cll_t *head, int node, cSet_t *newClique) {
	int foundStatus = 0;
	cll_t * curr = head;
	cll_t * prev = NULL;
	int i = 0, j = 0, breakFlag = 0;

	while(curr != NULL) {
		if (curr->id == node) {
			foundStatus = 1;
			break;
		}
		curr = curr->next;
	}

	if (foundStatus == 0) {
		fprintf(stderr, "\nFirst clique not found! (removeSupers)\n");
		fflush(stderr);
		exit(0);
	}
	
	// Now this is trickier, to remove nodes from the middle of a linked
	//  list; this means that we need to remember which node we were just
	//  at so that we can connect it to the node after the one we are 
	//  about to delete.
	prev = curr;
	curr = curr->next;

	// This code is similar to that in uniqClique.
	// Descend through all members of the next level's linked list.
	while(curr != NULL) {
		i = 0;
		j = 0;
		breakFlag = 0;
		// The proposed convolved clique will be referred to as the
		//  'first' clique, and the current clique being analyzed
		//  in the next level is the 'second' clique.
		// Continue if we have more members in both cliques.  We will
		//  have already broken out if it is not possible for this
		//  second clique to be a subset of the first.
		while((i < newClique->size) && (j < curr->set->size)) {
			// If the current member of the first clique is
			//  less than the current member of the second clique
			//  then it is still possible that the first is a 
			//  superset of the second, so move on to the next
			//  member.
			if (newClique->members[i] < curr->set->members[j]) {
				i++;
			}
			// If the current member of the first clique is greater
			//  than the current member of the second clique, then
			//  the proposed second clique cannot be a subset since
			//  its members are all in ascending order.  We also
			//  know that since the first clique already has
			//  a subset in this linked list, the current node
			//  cannot possibly be a superset of the proposed
			//  clique, so we can just disregard that.  Thus,
			//  we make a flag signifying this and break out.
			else if (newClique->members[i] > curr->set->members[j]) {
				breakFlag = 1;
				break;
			}
			else {
				i++;
				j++;
			}
		}
		// If the breakflag is 1, then we know
		//  that there is a member of the second clique not in the
		//  first, and so the second is not a subset.  If the breakflag
		//  is 0 but j is less than the second clique's size, then 
		//  we must have broken because we ran out of members in the
		//  first clique... thus, there is a member of the second 
		//  clique not in the first.  Thus, only if the breakflag is
		//  0 and j is equal to the size of the second clique do we
		//  know that every member of the second clique is in the first
		//  and that the second clique can thus be removed.
		if((breakFlag == 0) && (j == curr->set->size)) {
			// Make the previous clique point to the next one
			//  instead of the current one.
			prev->next = curr->next;
			// Free all of the memory used by the current clique.
			free(curr->set->members);
			free(curr->set);
			free(curr);
			curr = prev->next;
		} else {
			// Otherwise, the current second clique is not a 
			//  subset of the first, and we advance the prev and
			//  curr pointers.
			prev = curr;
			curr = curr->next;
		}
	}
	return head;
}

/*********************************************************************
 * Prints out the contents of a cSet_t in the following format:
 * 	Support = <number of nodes in clique>
 * 	Members = <newline-separated list of nodes in clique>
 * Input: a clique in the form of a cSet_t
 * Output: in text, the contents of the cSet.  An integer is returned as well,
 *  with 1 indicating success.
 *  *********************************************************************/
int
printCSet(cSet_t *node) {
        int i = 0;
	if (node->size == 0) {
		fprintf(stderr,"cSet has no members! - printCSet\n");
		fflush(stderr);
		exit(0);
	} else {
        	printf("\nSupport = %d\n", node->size);
        	printf("Members = \n");
        	for (i = 0; i < node->size; i++) {
                	printf("\t%d\n", node->members[i]);
        	}
        	return 1;
	}
}

/*********************************************************************
 * Pushes a freshly-convolved clique, currently in mll_t form, onto the 
 *  clique linked list for the next level.  Also checks to make sure that
 *  the convolved clique is unique, and if it isn't, it takes appropriate 
 *  action.
 * Input: a convolved clique in mll_t form, the head of a clique linked list
 *  for the next level
 * Output: (potentially new) head of the clique linked list for the next level
 * *************************************************************************/
cll_t *
pushConvClique(mll_t *clique, cll_t *head) {
        int status = 0;
	cSet_t *cliquecSet = NULL;

	// First change the clique to something we can used more easily
        cliquecSet = mllToCSet(clique);
	// Then check to make sure it's unique by finding out its status
	status = uniqClique(cliquecSet, head);
	
//	printf("Candidate:\n");
//	printCSet(cliquecSet);
	
	// If we get -2, then this clique is a subset, so just free 
	//  the cSet we just made and move on.
	if (status == -2) {
		free(cliquecSet->members);
		free(cliquecSet);
		cliquecSet = NULL;
	}
	// If we get -1, then this is a unique clique, so push it on.
	else if (status == -1) {
		head = pushcSet(head,cliquecSet);
	}
	// Otherwise, this clique is a superset, so we'll first remove
	//  all of the other cliques of which this is a superset.  Then
	//  we'll swap out the first clique of which this is a superset
	//  with this current clique.  The clique being removed is free'd 
	//  within the swapNode function.
	else {
		head = removeSupers(head,status,cliquecSet);
		head = swapNodecSet(head,status,cliquecSet);
	}
        return head;
}

/***********************************************************************
 * Turns a member linked list used to store the intersection of two cliques
 *  into something more useful: a cSet_t structure.
 * Input: a clique in mll_t form
 * Output: a clique in cSet_t form
 * **********************************************************************/
cSet_t *
mllToCSet(mll_t *clique) {
	int sizecount = 0, i = 0;
	cSet_t *cliqueCset = malloc(sizeof(cSet_t));
	mll_t *head = clique;
	if (cliqueCset == NULL) {
		fprintf(stderr,"Memory error - mllToCSet cSet\n%s\n",
				strerror(errno));
		fflush(stderr);
		exit(0);
	}

	// First count up how many members there are in the member linked list
	while (head != NULL) {
		sizecount++;
		head = head->next;
	}
	
	head = clique;
	cliqueCset->size = sizecount;
	cliqueCset->members = (int *) malloc(sizecount * sizeof(int));

	if (cliqueCset->members == NULL) {
		fprintf(stderr,"Memory error - mllTlCSet cliquemembers\n%s\n",
				strerror(errno));
		fflush(stderr);
		exit(0);
	}
	// In order to stay in the same format as with bitSet translation to
	//  cSet, we ensure that the ids of the members are ascending with
	//  ascending index number in the cSet.  This is accomplished by noting
	//  that since the intersection members are pushed onto the stack,
	//  a LIFO operation, that the first intersected nodes off the stack
	//  will have the highest ids, so we will put them at the end of
	//  the members array with the higher index values.
	for (i = sizecount - 1; i >= 0; i--) {
		cliqueCset->members[i] = head->cliqueMembership;
		head = head->next;
	}
	
	return cliqueCset;
}

/***************************************************************************
 * Convolves one single clique against all possible cliques that could 
 *  possibly be convolved.  It does not attempt to convolve all other cliques,
 *  but prunes that set by first looking at the offsets that are in the clique,
 *  then collecting all of the cliques who have members that are one greater
 *  than the offsets in this clique, and then convolving those cliques in
 *  a sort of "queue" using the bitSet data structure.
 * Input: the head of the clique linked list for the current level, the current
 *  node being convolved against in the linked list, the location of the
 *  previous node in the form of a pointer to a "guess", an array of member
 *  linked lists, the length of that array, the head of the clique linked list
 *  for the next level, a bitSet for the printStatus of maximality, and 
 *  the support criterion.
 * Output: the head of the (possibly modified) clique linked list for 
 *  the next level
 *  *************************************************************************/

cll_t *
wholeCliqueConv(cll_t *head, cll_t *node, cll_t **firstGuess, mll_t **memList, 
	int numOffsets, cll_t *nextPhase, bitSet_t *printStatus, int support){
	bitSet_t *queue = NULL;
	cSet_t *cliquesToSearch = NULL;
	int i = 0;
	cll_t **secondGuess = NULL;

	// This bitSet will be used to create a "queue" of the different
	//  cliques that must be convolved against the current primary clique.
	//  A bitset is used to make it easy to deal with duplicates, where
	//  multiple clique members' next offsets
	//  are all members of some other specific clique.
	queue = newBitSet(head->id + 1);
	queue = searchMemsWithList(node->set->members, node->set->size,
					memList, numOffsets, queue);
	// We'll use this "secondGuess" to store where the previous clique
	//  being convolved was... since we will progressing monotonically
	//  in descending order, this will save us some time in traversing the
	//  linked list looking for the clique that we want.
	secondGuess = (cll_t **) malloc(sizeof(cll_t *));
	if (secondGuess == NULL) {
		fprintf(stderr,"Memory error - wholeCliqueConv\n%s\n",
				strerror(errno));
		fflush(stderr);
		exit(0);
	}

	// If the offsets that we are looking for are in no other cliques,
	//  we can just bail out now.
	if(countSet(queue) == 0) {
		deleteBitSet(queue);
		return nextPhase;
	}
	
	// Otherwise, we start our secondGuess at the head and get going.
	*secondGuess = head;

	// We change the bitSet to something more useful.
	cliquesToSearch = bitSetToCSet(queue);
	
	// Note that we start from the end of the cSet member list so that
	//  we can convolve the highest-id cliques first, which are at the 
	//  beginning of our stack of cliques.
	for (i = cliquesToSearch->size - 1; i >= 0; i--) {
		nextPhase = singleCliqueConv(head,node->id,firstGuess,
			cliquesToSearch->members[i],secondGuess,
			nextPhase,printStatus,support);
	}

	// And then we free everything that we created
	deleteBitSet(queue);
	free(cliquesToSearch->members);
	free(cliquesToSearch);
	free(secondGuess);
	return nextPhase;
}

/************************************************************************
 * Performs convolution on all cliques in a linked list by repeatedly calling
 *  wholeCliqueConv.
 * Input: pointer to the head of a clique linked list for the current level,
 *  array of member linked lists, length of that array, minimum support
 *  threshold, the current length of motifs, and a pointer to a linked list
 *  containing all cliques that will be printed out.
 * Output: the head of the clique linked list for the next level of convolution
 * *************************************************************************/

cll_t *
wholeRoundConv(cll_t **head,mll_t **memList,int numOffsets,int support,int length,
		cll_t **allCliques){
	bitSet_t *printStatus = NULL;
	cll_t *curr = *head;
	cll_t *prev = NULL;
	cll_t *nextPhase = NULL;
	cll_t **firstGuess = NULL;

	// Create a bitset to keep track of print status for this level.
	//  It starts off all true, and gets changed to false if the patterns
	//  are not maximal.
	printStatus = newBitSet((*head)->id + 1);
	fillSet(printStatus);
	firstGuess = (cll_t **) malloc(sizeof(cll_t *));
	if (firstGuess == NULL) {
		fprintf(stderr,"Memory error - wholeRoundConv\n%s\n",
				strerror(errno));
		fflush(stderr);
		exit(0);
	}

	// Start off at the head.
	*firstGuess = *head;	
	// Convolve a whole clique at a time, traversing the linked list.
	//  Note that firstGuess gets altered within the function.
	while(curr != NULL) {
		nextPhase = wholeCliqueConv(*head,curr,firstGuess,memList,numOffsets,
						nextPhase,printStatus,support);
		curr = curr->next;
	}
	
	// Now go back to the head for printing output
	curr = *head;

//	printf("\n****************************************************\n");
//	printf("Length = %d", length);
//	printf("\n****************************************************\n");

	// For each clique that is still 'true' in printStatus and is thus
	//  maximal, perform some sort of output.  Yankcll will pull out the
	//  clique and save it for printing at a later time.
	while(curr != NULL) {
		if (checkBit(printStatus, curr->id)) {
// 			This is the line that makes the allCliques output.
//				Can either printcll, or add to allCliques.
//			printCllPattern(curr, length);
			yankCll(head,prev,&curr,allCliques,length);
		} else {
			prev = curr;
			curr = curr->next;
		}
	}

	// And clean up.
	deleteBitSet(printStatus);
	free(firstGuess);
	return nextPhase;
}

/***************************************************************************
 * Removes a clique from within a linked list in order to save it for later
 *  printing.  This is done so that the cliques are not printed as they are
 *  convolved, but rather after all rounds of convolution are complete.
 * Input: a pointer to the head of the current linked list, the clique prior
 *  to the one that is to be yanked (NULL if the clique to be yanked is
 *  the head), the clique that is to be yanked, a pointer to the head of the
 *  list with all cliques that are to be printed, and the length of the
 *  current motif.
 * Output: Nothing is returned beyond a success integer, but it alters
 *  the current level cll, the value of curr, and the linked list of all
 *  cliques that are to be printed.  
 *  ***********************************************************************/
int
yankCll(cll_t **head, cll_t *prev, cll_t **curr, cll_t **allCliques, int length){
	if (*curr == NULL) {
		fprintf(stderr,"\nCan't yank from end of cll!\n");
		fflush(stderr);
		exit(0);
	}
// If we're not on the head, change the previous node's "next".
// If we are on the head, make the new head be our current node's "next".
	if (prev != NULL) {
		prev->next = (*curr)->next;
	} else {
		*head = (*curr)->next;
	}

// Change next in curr, then change id and length information in curr
	(*curr)->next = *allCliques;

	if (*allCliques != NULL) {
		(*curr)->id = (*allCliques)->id + 1;
	} else {
		(*curr)->id = 0;
	}

	(*curr)->length = length;

	*allCliques = *curr;
	
	if (prev != NULL) {
		*curr = prev->next;
	} else {
		*curr = *head;
	}
	return (1);
}

/**********************************************************************
 * Performs complete convolution given the starting list of cliques.
 * Input: a pointer to the head of the initial clique linked list, the
 *  minimum support criterion value, the number of offsets in the sequence
 *  set, the minimum length of motifs (which is the length of motifs in
 *  the initial clique linked list), the index/Sequence data structure, and
 *  the value of the -p flag to prune based on unique sequence occurrences.
 * Output: a linked list of all maximal cliques based on the initial clique
 *  linked list.
 *  ********************************************************************/

cll_t *
completeConv(cll_t **head, int support, int numOffsets, int minLength, 
	int *indexToSeq, int p) {
	int i = 0;
	mll_t **memList = NULL;
	cll_t *nextPhase = NULL;
	cll_t *allCliques = NULL;
	int length = minLength;
	memList = (mll_t **) malloc(numOffsets * sizeof(mll_t *));
	if (memList == NULL) {
		fprintf(stderr,"Memory error - completeConv\n%s\n",
				strerror(errno));
		fflush(stderr);
		exit(0);
	}

	// The number of offsets will never change, so this can be defined
	//  now, though we will have to change what is in these arrays later.
	for (i = 0; i < numOffsets; i++) {
		memList[i] = NULL;
	}

// NOTE: This assumes that the elemPats all meet the support criterion

	// So we'll do this as long as the head is non-null.. that means that 
	//  the initial set of cliques must be non-null.  Those are then 
	//  convolved and the linked list for the next round is set to head,
	//  so this continues until the linked list for the "next round" at
	//  the end of some round is null.
	while(*head != NULL) {
		// First we get the inverse information for this round: find
		//  out which cliques each offset is a member of.
		memList = fillMemberStacks(*head, memList);
//		printf("numOffsets.bak = %d\n",numOffsets);
//		// Then we convolve a whole round.
		nextPhase = wholeRoundConv(head,memList,numOffsets,support,length,
						&allCliques);
		// Do some housekeeping.
		memList = emptyMemberStacks(memList, numOffsets);
		popAllCll(*head);
		// Enforce the -p flag for subsequent rounds.
		if (p > 1) {
			nextPhase = pruneCll(nextPhase,indexToSeq,p);
		}
		// And move on to the next round of convolution.
		*head = nextPhase;
		length++;
	}

	free(memList);

	return allCliques;
}

/***********************************************************************
 * Prints out the contents of a clique linked list node in this format:
 * 	Support = <number of motif occurrences> (id = <some id number>)
 * 	Members = <newline-separated list of offsets>
 * Input: a specific node to be output, the length of the motif inside it
 * Output: text per above, and an integer success value.
 * **********************************************************************/

int
printCllPattern(cll_t *node, int length) {
	int i = 0;

	printf("\nSupport = %d\t(id = %d)\n", node->set->size, node->id);
	printf("Members = \n");
	for (i = 0; i < node->set->size; i++) {
		printf("\t%d\n", node->set->members[i]);
	}
	return 1;
}
