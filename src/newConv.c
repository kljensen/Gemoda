#include "bitSet.h"
#include <errno.h>
#include "convll.h"
/*#include "efence.h"*/

/**********************************************************************
 * Recursive algorithm to exhaustively enumerate all of the maximal cliques
 *  that exist in the data.  This is one of the main workhorses of Gemoda
 *  when used in its exhaustive form.  This algorithm was originally 
 *  published by Etsuji Tomita, Akira Tanaka, and Haruhisa Takahasi
 *  as a Technical Report of IPSJ (Information Processing Society of
 *  Japan):
 *   Tomita, E, A Tanaka, & H Takahasi (1989).  "An optimal algorithm for 
 *    finding all of the cliques".  SIG Algorithms 12, pp 91-98.
 *  Input: a bitset with the nodes currently in the clique, a bitset with
 *   the candidates for expanding the clique, a bitset inidcating the current
 *   subgraph being searched, the bitGraph to be searched for cliques,
 *   the minimum support parameter, a counter variable for keeping track
 *   of how many nodes are in the current clique, a linked list of cliques
 *   that have been discovered so far, and a pointer to the data structure
 *   that dereferences offset indexes into sequence numbers, and the minimum
 *   number of unique sequences that must contain the motif.
 *  Output: integer success value of 0 (but more importantly, the elemPats
 *   clique linked list is expanded to contain all elementary (minimum-length)
 *   motif cliques.
 ************************************************************************/
int
findCliques(bitSet_t *Q, bitSet_t *cand, bitSet_t *mask, bitGraph_t *oG, int support, int qCount, cll_t **elemPats, int *indexToSeq, int p){
	bitSet_t **gammaOG = NULL;
	bitSet_t *candQ = newBitSet(oG->size);
	bitSet_t *newMask = newBitSet(oG->size);
	int i,q;
	int graphSize;
	int max = -1;
	int numBits;
	int u=0;
	int newMaskCount;
	int candQCount;
	graphSize = oG->size;

	//
	//  Find which vertex in subg maximizes |cand intersect gamma(u) |
	gammaOG = oG->graph;
	for ( i=0 ; i<graphSize ; i++){

		//  Don't check this vertex if it's masked
		if(!(checkBit(mask, i))){
			continue;
		}
		//  cand is always a subset of mask, so intersecting
		//  with mask is redundant
		bitSetIntersection(gammaOG[i], cand, candQ);
		numBits = countSet(candQ);
		if(numBits > max){
			u = i;
			max = numBits;
		}
	}

	//  Then do the extension of the q's
	qCount++;

	// This loop iterates over all possible values of cand - gamma() by
	//  iterating over all possible values of cand but immediately
	//  "continue"ing if the node is also in gamma(u)
	q = nextBitBitSet(cand, 0);
	while(q != -1){
		if (checkBit(gammaOG[u], q)) {
			q = nextBitBitSet(cand, q+1);
			continue;
		}

		//  SUBGq = SUBG i Gamma
		bitSetIntersection(mask, gammaOG[q], newMask);
		newMaskCount = countSet(newMask);

		setTrue(Q, q);
		// Only recurse if there are more candidates to be included,
		// and they will allow us to reach the minimum support. 
		if(newMaskCount > 0 && qCount + newMaskCount >= support){
			//  CANDq = CAND i Gamma
			bitSetIntersection(gammaOG[q], cand, candQ);
			candQCount = countSet(candQ);
			// only recurse if we can possibly get to a clique
			// of size with minimum support
			if(candQCount > 0 && qCount + candQCount >= support){
				// recursion with
				// new candidates, new mask, and original graph
				findCliques(Q, candQ, newMask, oG, support, qCount, elemPats, indexToSeq, p);
			}
		}else if(qCount >= support){
			// This should be done when:
			// 1. countSet(newMask) == 0 [connected subgraph is maximal]
			// 2. Qcount >= minCount [connected subgraph has enough nodes]
			*elemPats = pushClique(Q, *elemPats, indexToSeq, p);
		}

		// Remove q from Q, and remove q from cand
		setFalse(Q,q);
		setFalse(cand, q);

		q = nextBitBitSet(cand, q+1);
	}
	qCount--;
	deleteBitSet(candQ);
	deleteBitSet(newMask);
	return 0;
}

/*********************************************************************
 * A recursive routine for single linkage clustering.  This clustering is 
 *  much faster than exhaustively enumerating all cliques, but it puts 
 *  each node in only one cluster and is not guaranteed to give all 
 *  possible motifs.
 * Input: a bitSet containing the current motif, a bitSet containing
 *  candidates to be added to the current motif, a bitSet containing the
 *  current subgraph to be clustered, the original bitGraph to be clustered,
 *  the minimum support necessary for a motif to be returned, the current
 *  number of nodes in the motif, a linked list of elementary motifs (length
 *  is the same as the window size), pointer to a structure to derference 
 *  index values to sequence numbers, and the minimum number of unique
 *  sequences that a motif must be in to be returned.
 * Output: integer success value of 0 (but more importantly, the linked list
 *  elemPats is updated to contain all of the motifs of length = window size.
 *  *******************************************************************/
int
singleLinkage(bitSet_t *Q, bitSet_t *cand, bitSet_t *mask, bitGraph_t *oG, int support, int qCount, cll_t **elemPats, int *indexToSeq, int p){
	int i=0;
	int j=0;

	// go to the first vertex that has not been clustered yet
	i = nextBitBitSet(cand, 0);
	if(i != -1){

		// this vertex has been clustered
		setFalse(cand, i);


		// start a new cluster, Q
		copySet(oG->graph[i], Q);

		//  go over each vertex in the cluster
		j = nextBitBitSet(Q, 0);
		while( j != -1 ){

			//  if this vertex has been clustered already, skip it and go
			//  to the next one
			if(!checkBit(cand,j)){
				j = nextBitBitSet(Q, j+1);
				continue;
			}
			//  Add this vertex's neighbors to the current cluster
			bitSetUnion(Q, oG->graph[j], Q);

			//  This vertex has now been clustered
			setFalse(cand, j);
			//  go over each vertex in the cluster
			j = nextBitBitSet(Q, 0);
		}
		//  Did we make a cluster that was large enough?
		if(countSet(Q) >= support){
			*elemPats = pushClique(Q, *elemPats, indexToSeq, p);
		}
		// recurse
		singleLinkage(Q, cand, mask, oG, support, 0, elemPats, indexToSeq, p);

	}else{
		return 0;
	}
	return 0;
}

/**********************************************************************
 * The iterator used to "filter" the graph.  It takes information in the 
 *  bitset telling which nodes' rows have changed and only checks them...
 *  this should make it pretty efficient time-wise at only a small memory
 *  cost.
 * Note the convention that the first time this is called, the changed
 *  bitSet is empty... and that the master function is responsible for
 *  catching the signal that no changes were made in the last iteration.
 * Input: the bitGraph to be filtered, the minimum support required for
 *  a motif to be returned, a bitSet with nodes changed from the previous
 *  iteration, and a bitSet to export the nodes changed in this iteration.
 * Output: integer success value of 0 (and also a filtered bitGraph and
 *  a bitSet with the nodes changed in this iteration).
 *  ********************************************************************/

int
filterIter(bitGraph_t* graph, int support, bitSet_t* changed, bitSet_t* work) {
	int i = 0, j = 0;
	int lastBit = 0, nextBit = 0, lastRow = 0, nextRow = 0;
	int numNodes = 0;
	int changedSize = countSet(changed);
	emptySet(work);

	// Note the convention that the first time the function is called,
	//  it is done with an empty "changed" bitSet as a sentinel.  It is
	//  the responsibility of the master function calling the iterator
	//  to catch future empty changed sets to know that convergence has
	//  been achieved.
	//
	//  So, if it's your first time through, go through each node and make
	//   sure that each is connected to at least <support> - 1 others.
	if (changedSize == 0) {
		for (i = 0; i < graph->size; i++) {
			numNodes = countSet(graph->graph[i]);
			if (numNodes >= support - 1) {
				continue;
			} else {
				// Otherwise, zero it out, but going one by
				//  one so that you can also zero out the 
				//  symmetric bit.
				lastBit = 0;
				for (j = 0; j < numNodes; j++) {
					nextBit = nextBitBitSet(graph->graph[i],
						lastBit);
					if (nextBit == -1) {
						fprintf(stderr,
					"\nEnd of bitSet reached! - initial\n");
						fflush(stderr);
						exit(0);
					}
					setFalse(graph->graph[i],nextBit);
					setFalse(graph->graph[nextBit],i);
					// And set that corresponding bit true
					// in the work bitSet so that we
					// know we changed it for the next
					// round.
					setTrue(work,nextBit);
					lastBit = nextBit + 1;
				}
				
			}
		}
	} else {
		// Otherwise, we've been here before, so just follow what
		//  the changed bitSet says to do... only those bitSets that
		//  were changed could possibly have gone under the minimum
		//  support requirement.
		lastRow = 0;
		for (i = 0; i < changedSize; i++) {
			nextRow = nextBitBitSet(changed,lastRow);
			if (nextRow == -1) {
				fprintf(stderr, 
					"\nEnd of bitSet reached! - iter,row\n");
				fflush(stderr);
				exit(0);
			}
			// So now we've found the row that needs to be checked.
			//  We do the same thing we did above... either move
			//  on if it has enough possible support, or zero
			//  it out (with its symmetric locations) one by one.
			numNodes = countSet(graph->graph[nextRow]);
			if (numNodes >= support - 1) {
				lastRow = nextRow + 1;
				continue;
			} else {
				lastBit = 0;
				for (j = 0; j < numNodes; j++) {
					nextBit = nextBitBitSet(
							graph->graph[nextRow],
							lastBit);
					if (nextBit == -1) {
						fprintf(stderr,
					"\nEnd of BitSet reached! = iter,Bit\n");
						fflush(stderr);
						exit(0);
					}
					setFalse(graph->graph[nextRow],nextBit);
					setFalse(graph->graph[nextBit],nextRow);
					setTrue(work,nextBit);
					lastBit = nextBit + 1;
				}
				lastRow = nextRow + 1;
			}
		}
	}
	return 1;
}

/************************************************************************
 * Function to "filter" the initial bitGraph that is being clustered.  
 *  "Filtering" is the process of removing all nodes from the graph that
 *  cannot possibly be in motifs because they are not connected to enough
 *  other nodes.  This can be done once (if R != 1), or it can be done
 *  recursively (if R == 1).  When done recursively, it takes the just-filtered
 *  graph and checks all of the nodes that the recently removed node used
 *  to be connected to; since they have changed in connectivity, they may
 *  no longer be connected to enough nodes to be a member of a motif.  This
 *  is iterated until convergence.
 * Note that the default is to have recursive filtering on, as it ought to
 *  decrease the computational complexity of the clustering step and ought
 *  not have much of a computational footprint... in cases where it takes a
 *  while, it is probably having a good impact in the clustering step, 
 *  whereas if it is not effective, it probably won't take that long anyway.
 * Input: a bitGraph to be filtered, the minimum support that a motif must
 *  have, and the flag indicating recursive filtering or not.
 * Output: Integer success value of 0 (and an altered bitGraph so that all
 *  nodes with connections have at least <min support> connections).
 *  **********************************************************************/
int
filterGraph(bitGraph_t* graph, int support, int R) {
        bitSet_t * changed = newBitSet(graph->size);
        bitSet_t * work = newBitSet(graph->size);
        emptySet(changed);
        emptySet(work);
 
	// Iteratively call the filtering by copying the previous "work" into
	//  "changed" after each iteration step.
        if (R == 1) {
                do {
                        filterIter(graph,support,changed,work);
                        copySet(work,changed);
                }
                while (countSet(changed) > 0);
        } else {
		// Otherwise, just do it once.
                filterIter(graph,support,changed,work);
        }
                                                                                   
        deleteBitSet(changed);
        deleteBitSet(work);
                                                                                   
		return 0;
}

/**********************************************************************
 * Simple function (non-recursive) to prune off the first level of motifs
 *  that will not meet the "minimum number of unique sequences" criterion.
 *  This could have been implemented as above, but it may have gotten a little
 *  expensive with less yield, so only the first run through is done here.
 * Input: a bit graph to be pruned, a pointer to the structure that 
 *  dereferences offset indices to sequence numbers, a pointer to the 
 *  structure that dereferences seq/position to offsets, the number of 
 *  unique sequences in the input set, and the minimum number of unique
 *  sequences that must contain the motif.
 * Output: a pruned bitGraph.
 * **********************************************************************/
bitGraph_t *
pruneBitGraph(bitGraph_t *bg, int *indexToSeq, int **offsetToIndex, 
		int numOfSeqs, int p) {
	int i = 0, j = 0, nextBit = 0;
	int *seqNums = NULL;
	
	// Since we don't immediately know which node is in which source 
	//  sequence, we can't just count them up regularly.  Instead, we'll
	//  need to keep track of which sequences they come from and 
	//  increment _something_.  What we chose to do here is just make
	//  an array of integers of length = <p>.  Then, we try to put the
	//  source sequence number of each neighbor (including itself, since
	//  the main diagonal is still true at this time) into the next slot
	//  Since we will monotonically search the bitSet, we can just 
	//  move on to the first bit in the next sequence using the 
	//  offsetToIndex structure so that we know the next sequence number
	//  to be put in is always unique.
	seqNums = (int *) malloc(p * sizeof(int));
	if (seqNums == NULL) {
		fprintf(stderr,"Memory error - pruneBitGraph\n%s\n",strerror(errno));
		fflush(stderr);
		exit(0);
	}

	// So, for each row in the bitgraph...
	for (i = 0; i < bg->size; i++) {
		// Make sure the whole array is -1 sentinels.
		for (j = 0; j < p; j++) {
			seqNums[j] = -1;
		}
		j = 0;
		// Find the first neighbor of this bit.
		nextBit = nextBitBitSet(bg->graph[i],0);
		if (nextBit == -1) {
			continue;
		} else {
			// and put its sequence number in the array of ints.
			seqNums[0] = indexToSeq[nextBit];
		}
		// If it's the last sequence, then bail out so that we don't
		//  segfault in the next step.
		if (seqNums[0] >= numOfSeqs - 1) {
			emptySet(bg->graph[i]);
			continue;
		}
		// Find the next neighbor of this bit, STARTING AT the first
		//  bit in the next sequence.
		nextBit = nextBitBitSet(bg->graph[i],offsetToIndex[indexToSeq[nextBit] + 1][0]);
		// And iterate this until we run out of neighbors.
		while(nextBit >= 0) {
			j++;
			seqNums[j] = indexToSeq[nextBit];
			// Or until this new neighbor will fill up the array
			if (j == p - 1) {
				break;
			}
			// Or until this new neighbor is in the last sequence.
			if (seqNums[j] >= numOfSeqs - 1) {
				break;
			}
			// Get the next neighbor!
			nextBit = nextBitBitSet(bg->graph[i],
					offsetToIndex[indexToSeq[nextBit]+1][0]);
		}

		// If we didn't have enough unique sequences, and either a) we
		//  were in the nth-to-last sequence and there were no 
		//  neighbors after it, or b) we were in the last sequence,
		//  then the last number will still be our sentinel, -1.  If
		//  the last number is not a sentinel, then we have at least
		//  p distinct sequence occurrences, so we're OK.
		if (seqNums[p - 1] == -1) {
			emptySet(bg->graph[i]);
		}
	}
	free(seqNums);
	return(bg);
}

/***********************************************************************
 * Prunes a motif linked list of all motifs without support in at least <p>
 *  unique source sequences.
 * Input: head of a motif linked list, pointer to a structure that 
 *  dereferences offset indices to sequence numbers, minimum number of unique
 *  source sequences in which a motif must occur.
 * Output: head of a (potentially altered) motif linked list.
 * ********************************************************************/
cll_t * 
pruneCll(cll_t *head, int *indexToSeq, int p) {
	int i = 0, j = 0, thisSeq = 0;
	int *seqNums = NULL;
	cll_t *curr = head;
	cll_t *prev = NULL;
	cll_t *storage = NULL;

	// We'll do this similar to the pruneBitGraph function... we will
	//  keep track of which source sequence each motif occurrence was in.
	//  Again, since the occurrences are listed monotonically, we only
	//  need to compare the last non-sentinel index to the current
	//  sequence number.
	seqNums = (int *) malloc(p * sizeof(int));
	if (seqNums == NULL) {
		fprintf(stderr,"Memory error - pruneCll\n%s\n",strerror(errno));
		fflush(stderr);
		exit(0);
	}

	while(curr != NULL) {
		// First make sure the set size is at least p.
		// This is redundant, but extremely simple and not expensive,
		//  so we'll leave it in just as a check.
		if (curr->set->size < p) {
			if(prev != NULL) {
				prev->next = curr->next;
			} else {
				head = curr->next;
			}
			storage = curr->next;
			free(curr->set->members);
			free(curr->set);
			free(curr);
			curr = storage;
			continue;
		}

		for (i = 0; i < p; i++) {
			seqNums[i] = -1;
		}
		j = 0;

		seqNums[0] = indexToSeq[curr->set->members[0]];
		// Note, we've checked to make sure size > p, and we know
		// p must be 2 or greater, so we can start at 1 without
		// worrying about segfaulting
		for (i = 1; i < curr->set->size; i++) {
			thisSeq = indexToSeq[curr->set->members[i]];
			if (thisSeq != seqNums[j]) {
				j++;
				seqNums[j] = thisSeq;
				if (j == p - 1) {
					break;
				}
			}
		}

		// Same story as before... if the last number is -1,
		//  then we didn't have enough to fill up the <p> different
		//  slots, so this doesn't meet our criterion.
		if (seqNums[p - 1] == -1) {
			if (prev != NULL) {
				prev->next = curr->next;
			} else {
				head = curr->next;
			}
			storage = curr->next;
			free(curr->set->members);
			free(curr->set);
			free(curr);
			curr = storage;
		} else {
			prev = curr;
			curr = curr->next;
		}
	}
	free(seqNums);
	return(head);
}
		


/**********************************************************************
 * Our outer convolution function.  This function will call preliminary
 *  functions, cluster the data, and then call the main convolution function.
 *  This is the interface between the main gemoda-<x> code and the generic
 *  code that gets all of the work done.
 * Input: the bitGraph to be clustered and convolved, the minimum support
 *  necessary for a motif to be returned, a flag indicating whether 
 *  recursive filtering should be used, a pointer to the data structure that
 *  dereferences offset indices to sequence numbers, the number of unique
 *  source sequences that a motif must be present in, and a number indicating
 *  the clustering method that is to be used.
 * Output: the final motif linked list with all motifs that are to be 
 *  given as output to the user.
 *  *********************************************************************/

cll_t *
convolve(bitGraph_t *bg, int support, int R, int *indexToSeq, int p, 
	int clusterMethod, int **offsetToIndex,int numberOfSequences,
	int noConvolve, FILE *OUTPUT_FILE){
	bitSet_t *cand = NULL;
	bitSet_t *mask = NULL;
	bitSet_t *Q = NULL;
	int size = bg->size;
	cll_t *elemPats = NULL;
	cll_t *allCliques = NULL;
	cll_t *curr = NULL;

	// contains indices (rows) containing the threshold value.
	cand = newBitSet(size);
	mask = newBitSet(size);
	Q = newBitSet(size);
	fillSet(cand);
	fillSet(mask);

	// Note that we prune based on p before setting the diagonal false.
	if (p > 1) {
		bg = pruneBitGraph(bg,indexToSeq,offsetToIndex,numberOfSequences,p);
	}
	// Now we set the main diagonal false for clustering and filtering.
	bitGraphSetFalseDiagonal(bg);
	filterGraph(bg,support,R);

	fprintf(OUTPUT_FILE,"Graph filtered!  Now clustering...\n");
	fflush(NULL);

	if(clusterMethod==0){
		findCliques(Q, cand, mask, bg, support, 0, &elemPats, indexToSeq, p);
	}else{
		singleLinkage(Q, cand, mask, bg, support, 0, &elemPats, indexToSeq, p);
	}
	fprintf(OUTPUT_FILE,"Clusters found!  Now filtering clusters (if option set)...\n");
	fflush(NULL);
	if (p > 1) {
		elemPats = pruneCll(elemPats,indexToSeq,p);
	}

	deleteBitSet(cand);
	deleteBitSet(mask);
	deleteBitSet(Q);

	// Now let's convolve what we made.
	if (noConvolve == 0) {
		fprintf(OUTPUT_FILE,"Now convolving...\n");
		fflush(NULL);
		allCliques = completeConv(&elemPats,support,size,0,indexToSeq,p);
	}
	else {
		curr = elemPats;
		while(curr != NULL) {
			yankCll(&elemPats,NULL,&curr,&allCliques,0);
		}
	}


	return allCliques;
}
