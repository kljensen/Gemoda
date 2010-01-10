#include "errno.h"
#include "bitSet.h"
/*#include "efence.h"*/

/****************************************************************
 * Creates a bit array for use in high-throughput intersections/unions.
 * Input: desired size of bit array in bytes
 * Output: a new bit array in bit_t format
 * Note: this should not be called directly; see newBitSet
 * **************************************************************/
bit_t *
newBitArray(int bytes)
{
	bit_t *b = (bit_t *) malloc(bytes);
	if (b == NULL) {
		fprintf(stderr,"\nMemory error --- couldn't allocate bitArray!"
				" - newBitArray\n%s\n",strerror(errno));
		fflush(stderr);
		exit(0);
	}
	// Set them all false
	memset(b, 0, bytes);
	return b;
}

/****************************************************************
 * Creates a bitSet data structure that contains a bit array and
 *  information about that bit array that is necessary for quick and
 *  efficient access of the array.
 * Input: the desired length of the bit array
 * Output: a bitSet data structure
 * ***************************************************************/
bitSet_t *
newBitSet(int size)
{
	bitSet_t *s1 = (bitSet_t *) malloc(sizeof(bitSet_t));
	if (s1 == NULL) {
		fprintf(stderr,"\nMemory error --- couldn't allocate biSet!"
				" - newBitSet\n%s\n",strerror(errno));
		fflush(stderr);
		exit(0);
	}
	// Fill in details about the bitSet, allocate bitSet
	s1->max = size;
	s1->slots = BSNUMSLOTS(size);
	s1->bytes = s1->slots * sizeof(bit_t);
	s1->tf = newBitArray(s1->bytes);
	return s1;
}

/******************************************************************
 * Sets a specific bit in a bitSet as true.
 * Input: a bitSet, the number of the bit to be set as true
 * Output: integer success value of 0
 * ***************************************************************/
int
setTrue(bitSet_t * s1, int x)
{
/*	if (BSNUMSLOTS(x) > s1->slots) {
  Conditional changed, 5/25, by MPS: check x against s1->max, should
    be safer
*/
	if (x >= s1->max) {
		fprintf(stderr, "Set isn't large enough! - setTrue\n");
		fflush(stderr);
		exit(0);
	}
	BSSET(s1->tf, x);
	return 0;
}

/******************************************************************
 * Sets a specific bit in a bitSet as false.
 * Input: a bitSet, the number of the bit to be set as false
 * Output: integer success value of 0
 * ***************************************************************/

int
setFalse(bitSet_t * s1, int x)
{
/*	if (BSNUMSLOTS(x) > s1->slots) {
  Conditional changed, 5/25, by MPS: check x against s1->max, should
    be safer
*/
	if (x >= s1->max) {
		fprintf(stderr, "Set isn't large enough! - setFalse\n");
		fflush(stderr);
		exit(0);
	}
	BSCLEAR(s1->tf, x);
	return 0;
}

/***********************************************************************
 * Inverts all values in a bitSet, making all trues false and all 
 *  falses true.
 * Input: a bitSet
 * Output: integer success value of 0
 * ********************************************************************/
int flipBits(bitSet_t *s1){
	int i;
	for (i = 0; i < s1->slots; i++) {
		s1->tf[i] = ~s1->tf[i];
	}
	return 0;
}

/***********************************************************************
 * Sets all values in a bitSet to true.
 * Input: a bitSet
 * Output: integer success value of 0
 * *********************************************************************/
int
fillSet(bitSet_t * s1)
{
	memset(s1->tf, ~0, s1->bytes);
	return 0;
}

/***********************************************************************
 * Sets all values in a bitSet to false.
 * Input: a bitSet
 * Output: integer success value of 0
 * *********************************************************************/
int
emptySet(bitSet_t * s1)
{
	memset(s1->tf, 0, s1->bytes);
	return 0;
}

/**********************************************************************
 * Finds the value of a specific bit in a bitSet.
 * Input: a bitSet, the number of the bit being queried
 * Output: the value of the bit being queried (1 or 0)
 * ********************************************************************/
int
checkBit(bitSet_t * s1, int x)
{
	return BSTEST(s1->tf, x);
}

/*********************************************************************
 * Performs memory management for the deletion of a bitSet_t structure.
 * Input: a bitSet
 * Output: integer success value of 0
 * ********************************************************************/
int
deleteBitSet(bitSet_t * s1)
{
	if (s1->tf != NULL) {
		free(s1->tf);
		s1->tf = NULL;
	}
	if (s1 != NULL) {
		free(s1);
		s1 = NULL;
	}
	return 0;
}

/*****************************************************************
 * Finds the union of two bitSets
 * Input: first bit set for the union, second bit set for the union,
 *  a bit set in which to store the results
 * Output: an integer success value of 0 (and an altered third bit set
 *  with the results of the union
 *  *******************************************************************/
int
bitSetUnion(bitSet_t * s1, bitSet_t * s2, bitSet_t * s3)
{
	int i;
	if ((s1->slots != s2->slots) || (s1->slots != s3->slots)) {
		fprintf(stderr, "Sets aren't same size!\n");
		fflush(stderr);
		exit(0);
	}
	for (i = 0; i < s1->slots; i++) {
		s3->tf[i] = BSUNION(s1->tf[i], s2->tf[i]);
	}
	return 0;
}

/**********************************************************************
 * Copies the true/false contents of one bit set into an existing bit set.
 *  Both bit sets must be the same size.
 * Input: source bit set, destination bit set
 * Output: integer success value of 0 (and an altered destination bitset)
 * *******************************************************************/
int copySet(bitSet_t *s1, bitSet_t *s2){
	int i;
	if(s1->slots != s2->slots){
		fprintf(stderr, "Sets are not the same size!");
		fflush(stderr);
		exit(0);
	}
	for ( i=0 ; i<s1->slots ; i++){
		s2->tf[i] = s1->tf[i];
	}
	return 0;
}

/**********************************************************************
 * Copies the true/false contents of one bit graph into an existing bit
 *  graph.  Both bit graphs must be the same size, and each corresponding
 *  bit set between the two bit graphs must be the same size.
 * Input: source bit graph, destination bit graph
 * Output: integer success value of 0 (and an altered destination bit graph)
 * **********************************************************************/
int copyBitGraph(bitGraph_t *bg1, bitGraph_t *bg2){
	int i;
	if(bg1->size != bg2->size){
		fprintf(stderr, "Graphs are not the same size!");
		fflush(stderr);
		exit(0);
	}
	for ( i=0 ; i<bg1->size ; i++ ){
		copySet(bg1->graph[i], bg2->graph[i]);
	}
	return 0;
}

/**********************************************************************
 * Locates all differences between two bitSets.  The result bitSet contains
 *  a true at a given bit if the two source bitSets differ at that bit. 
 * Input: first bit set to be compared, second bit set to be compared,
 *  third bit set to store the results
 * Output: integer success value of 0 (and an altered destination bit set
 *  with a true where the two source bit sets differed)
 *  ********************************************************************/
int
bitSetDifference(bitSet_t * s1, bitSet_t * s2, bitSet_t * s3)
{
	int i;
	if ((s1->slots != s2->slots) || (s1->slots != s3->slots)) {
		fprintf(stderr, "Sets aren't same size!\n");
		fflush(stderr);
		exit(0);
	}
	for (i = 0; i < s1->slots; i++) {
			s3->tf[i] = (s1->tf[i] & (~s2->tf[i]));
	}
	return 0;
}

/*************************************************************************
 * "Adds" two bit sets together... Currently unknown functionality, not
 *  used in existing code.
 *  **********************************************************************/

int
bitSetSum(bitSet_t * s1, bitSet_t * s2, bitSet_t * s3)
{
	int i;
	if ((s1->slots != s2->slots) || (s1->slots != s3->slots)) {
		fprintf(stderr, "Sets aren't same size!\n");
		fflush(stderr);
		exit(0);
	}
	for (i = 0; i < s1->slots; i++) {
		s3->tf[i] = (s1->tf[i] + s2->tf[i]);
	}
	return 0;
}

/***********************************************************************
 * Finds the intersection of two bitsets.
 * Input: First bitSet to be intersected, second bitSet to be intersected,
 *  a bitSet to store the result of the intersection.
 * Output: Integer success value of 0 (and an altered destination bitSet
 *  with a true where both source bitSets had a true).
 *  *******************************************************************/

int
bitSetIntersection(bitSet_t * s1, bitSet_t * s2, bitSet_t * s3)
{
	int i;
	if ((s1->slots != s2->slots) || (s1->slots != s3->slots)) {
		fprintf(stderr, "Sets aren't same size!\n");
		fprintf(stderr, "set 1 slots = %d\n", s1->slots);
		fprintf(stderr, "set 2 slots = %d\n", s2->slots);
		fprintf(stderr, "set 3 slots = %d\n", s3->slots);
		fflush(stderr);
		exit(0);
	}
	for (i = 0; i < s1->slots; i++) {
		s3->tf[i] = BSINTERSECTION(s1->tf[i], s2->tf[i]);
	}
	return 0;
}

/**********************************************************************
 * Finds the intersection of 3 bitSets.
 * Input: First bitSet to be intersected, second bitset to be intersected,
 *  third bitSet to be intersected, a bitSet to store the result of the
 *  intersection.
 * Output: Integer success value of 0 (and an altered destination bitSet
 *  with a true where all three source bitSets had a true.)
 *  ********************************************************************/
int
bitSet3WayIntersection(bitSet_t * s1, bitSet_t * s2, bitSet_t * s3, bitSet_t * s4)
{
	int i;
	if ((s1->slots != s2->slots) || (s1->slots != s3->slots) || (s1->slots != s4->slots)) {
		fprintf(stderr, "Sets aren't same size!\n");
		fflush(stderr);
		exit(0);
	}
	for (i = 0; i < s1->slots; i++) {
		s4->tf[i] = BSINTERSECTION(s1->tf[i], s2->tf[i]);
		s4->tf[i] = BSINTERSECTION(s3->tf[i], s4->tf[i]);
	}
	return 0;
}

/*********************************************************************
 * Attempt at a fast way of counting how many true values are in a given 
 *  bitSet.  Currently deprecated, using precompiled version instead.
 *  *****************************************************************/

int
bitcount32(unsigned int n)
{
	/*
	   works for 32-bit numbers only 
	 */
	/*
	   fix last line for 64-bit numbers 
	 */

	register unsigned int tmp;

	tmp = n - ((n >> 1) & 033333333333)
		- ((n >> 2) & 011111111111);
	return ((tmp + (tmp >> 3)) & 030707070707) % 63;
}

/*********************************************************************
 * Data structure for storing the number of true bits in a given char (char
 *  assumed to be 8 bits).
 *  *******************************************************************/

static int bits_in_char[256] = {
	0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,	/* 0- 15 */
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,	/* 16 - 31 */
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,	/* 32 - 47 */
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,	/* 48 - 63 */
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,	/* 64 - 79 */
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,	/* 80 - 95 */
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,	/* 96 - 111 */
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,	/* 112 - 127 */
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,	/* 128 - 143 */
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,	/* 144 - 159 */
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,	/* 160 - 175 */
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,	/* 176 - 191 */
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,	/* 192 - 207 */
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,	/* 208 - 223 */
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,	/* 224 - 239 */
	4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8	/* 240 - 255 */
};

/**********************************************************************
 * Uses bits_in_char data structure to determine the number of true bits
 *  in a 32-bit int in an efficient manner.
 * Input: 32-bit int (equal to one slot in the bitSet).
 * Output: number of true bits in the input integer.
 * ********************************************************************/

int
bitcount32_precomp(unsigned int n)
{
	// works only for 32-bit ints

	return bits_in_char[n & 0xffu]
		+ bits_in_char[(n >> 8) & 0xffu]
		+ bits_in_char[(n >> 16) & 0xffu]
		+ bits_in_char[(n >> 24) & 0xffu];
}



/***********************************************************************
 * Currently there is no support for 64-bit architectures.
 * *********************************************************************/

/*
//  Count bits for 64bit integers
#define PCTWO(c)     (0x1u << (c))
#define PCMASK(c)    (((unsigned int)(-1)) / (PCTWO(PCTWO(c)) + 1u))
#define PCCOUNT(x,c) ((x) & PCMASK(c)) + (((x) >> (PCTWO(c))) & PCMASK(c))

int
bitcount64(unsigned int n)
{
	n = PCCOUNT(n, 0);
	n = PCCOUNT(n, 1);
	n = PCCOUNT(n, 2);
	n = PCCOUNT(n, 3);
	n = PCCOUNT(n, 4);
	n = PCCOUNT(n, 5) ;// for 64-bit integers 
	return n;
}
*/

/*************************************************************************
 * Counts the number of true values in a bitSet.
 * Input: bitSet_t
 * Output: number of true values in that bitSet.
 * **********************************************************************/

int
countSet(bitSet_t * s1)
{
	int i;
	int sum = 0;
	int (*bitCounter) () = &bitcount32_precomp;
// Currently there is no support for 64-bit architectures.
	/*
		if(sizeof(bit_t) * 8 == 32){
			bitCounter = &bitcount32_precomp;
		}else if(sizeof(bit_t) * 8 == 64){
			bitCounter = &bitcount64;
		}else{
			exit(0);
		}
	*/

	if(sizeof(bit_t)*8 != 32) {
		fprintf(stderr,"\nSorry, no support for 64-bit architectures just yet! - countSet\n");
		fflush(stderr);
		exit(0);
	}
	

	// Just count the number of true bits in each char, and do this for
	//  (num of chars per int) chars.
	for (i = 0; i < s1->slots; i++) {
		sum += bitCounter(s1->tf[i]);
	} 
	return sum;
}


/************************************************************************
 * Finds the index of the first non-zero bit at-or-after start.
 * Input: a bitSet_t to be searched, the index of the start bit.
 * Output: the index of the first non-zero bit at-or-after start.
 * *********************************************************************/

int
nextBitBitSet(bitSet_t *s1, int start)
{
	//  slot is our starting slot, the
	//  slot containing bit 'start'
	int slot = BITSLOT(start);
	int i;
	// stop is the bit to stop it --- it is equal to max, and it is
	//   the index of a bit that does NOT belong to the bitset
	int stop;
	bit_t bitFalse; 
	memset(&bitFalse,0,sizeof(bit_t));


	//  s1->max is the number of bits in s1
	//  test to see if we're looking too high
	if (start >= s1->max) {
		return -1;
	}

	//  s1->slots is the number of available slots
	//  skip over empty slots
	while(slot < s1->slots){
		/*printf("w");*/
		if(s1->tf[slot] != bitFalse){
			// this slot is not empty

			//  if each slot is, say 32 bits and 
			//  we asked for nextBitBitSet(s1, 5),
			//  then slot 0 will be non-zero.  but,
			//  instead of starting at 0, start at 5!
			if(BSBITSIZE * slot > start){
				//  set start to index of first
				//  bit in this slot
				start = BSBITSIZE * slot;
			}
			//  set the stop, with a a check against the 'max'
			//   element of the bitSet_t object
			if(BSBITSIZE*(slot+1) > s1->max){
				stop = s1->max;
			} else {
				stop = BSBITSIZE*(slot+1);
			}
			for( i=start ; i<stop ; i++){
				if(checkBit(s1, i)){
					return i;
				}
			}
		}
		slot++;
	}
	return -1;
}

/***********************************************************************
 * Counts the number of true (non-zero) values in a bitgraph.
 * Input: a bitGraph_t
 * Output: the integer number of true (non-zero) values in the bitGraph
 * **********************************************************************/

int
countBitGraphNonZero(bitGraph_t *bg){
	int i;
	int sum=0;
	// Iterate over all bitSets in the bitGraph
	for ( i=0 ; i<bg->size ; i++){
		sum += countSet(bg->graph[i]);
	}
	return sum;
}

/*****************************************************************
 * Prints a representation of a bitSet_t data structure in the following form:
 *  bitSet (addr = <address>; <# of members> members)
 *   max = <number of bits in the bitSet>
 *   slots = <number of slots (ints) in the bitSet>
 *   bytes = <number of bytes (chars) in the bitSet>
 *   members = <space-delimited list of true value indices>
 * Input: a bitSet_t to be displayed
 * Output: integer succses value of 0 (and the stdout text described above)
 * *********************************************************************/

int
printBitSet(bitSet_t * s1)
{
	int i;
	printf("bitSet (addr = %d; %d members)\n", (int) s1, countSet(s1));
	printf("\tmax = %d\n", s1->max);
	printf("\tslots = %d\n", s1->slots);
	printf("\tbytes = %d\n", s1->bytes);
	printf("\tmembers =");

/*	for (i = 0; i <= s1->max; i++) {
    Changed 5/25, by MPS: Rest of code consistently says that max is 
      the size of the array, and not an index.  So, if we want to print
      the bitSet, we want to go from index 0 to index (max - 1). 
      This change makes this function consistent with the 
      printBinaryBitSet function.
*/

	for (i = 0; i < s1->max; i++) {
		if (BSTEST(s1->tf, i)) {
			printf(" %d", i);
		}
	}
	printf("\n");
	return 0;
}

/************************************************************************
 * Finds the union of two rows (bitSets) within a bitGraph
 * Input: a bitGraph, first row to be compared, second row to be compared,
 *  a bitSet_t to store the union results
 * Output: integer success value of 0 (and an altered destination bitSet
 *  with a true value wherever one or both source bitSets had a true value)
 *  *********************************************************************/

int bitGraphRowUnion(bitGraph_t *bg, int row1, int row2, bitSet_t *s1){
	bitSetUnion(bg->graph[row1], bg->graph[row2], s1);
	return 0;
}

/**************************************************************************
 * Finds the intersection of two rows (bitSets) within a bitGraph.
 * Input: a bitGraph, first row to be compared, second row to be compared,
 *  a bitSet_t to store the intersection results
 * Output: integer success value of 0 (and an altered destination bitSet
 *  with a true value wherever both source bitSets had a true value)
 *  *********************************************************************/

int bitGraphRowIntersection(bitGraph_t *bg, int row1, int row2, bitSet_t *s1){
	bitSetIntersection(bg->graph[row1], bg->graph[row2], s1);
	return 0;
}


/*************************************************************************
 * Prints a representation of a bitSet_t structure as a string of 1's and 0's
 * Input: a bitSet to be printed
 * Output: integer success value of 0 (and the stdout text described above)
 * ***********************************************************************/
int
printBinaryBitSet(bitSet_t *s1){
	int i;
	for (i = 0; i < s1->max; i++) {
		printf("%d",  (BSTEST(s1->tf, i) ? 1: 0));
	}
	return 0;
}

/*************************************************************************
 * Checks the value of a bit in a bitGraph
 * Input: a bitGraph, the index of the row of the bitGraph with the bit to 
 * be checked, the index of the bit in that row that is to be checked
 * Output: the value of the bit in the bitGraph being checked.
 * **********************************************************************/

int
bitGraphCheckBit(bitGraph_t *bg, int x, int y){
	return checkBit(bg->graph[x], y);
}


/*************************************************************************
 * Sets a specific bit in a bitGraph true.
 * Input: a bitGraph, the index of the row of the bitGraph with the bit to
 *  be set, the index of the bit in that row that is to be set
 * Output: integer success value of 0 (and an altered bitGraph)
 * ***********************************************************************/

int
bitGraphSetTrue(bitGraph_t *bg, int x, int y){
	setTrue(bg->graph[x], y);
	return 0;
}

/*************************************************************************
 * Sets a specific bit in a bitGraph false.
 * Input: a bitGraph, the index of the row of the bitGraph with the bit to
 *  be set, the index of the bit in that row that is to be set
 * Output: integer success value of 0 (and an altered bitGraph)
 *************************************************************************/

int
bitGraphSetFalse(bitGraph_t *bg, int x, int y){
	setFalse(bg->graph[x], y);
	return 0;
}

/*************************************************************************
 * Sets a specific bit and its symmetric opposite in a bitGraph false.  For
 *  instance, given that we wanted to set the 3rd bit in the 5th row false,
 *  this would also set the 5th bit in the 3rd row.
 * Input: a bitGraph, the index of the row of the bitGraph with the bit to
 *  be set, the index of the bit in that row that is to be set
 * Output: integer success value of 0 (and an altered bitGraph)
 * **********************************************************************/
int
bitGraphSetFalseSym(bitGraph_t *bg, int x, int y){
	setFalse(bg->graph[x], y);
	setFalse(bg->graph[y], x);
	return 0;
}

/*************************************************************************
 * Sets a specific bit and its symmetric opposite in a bitGraph true.  For
 *  instance, given that we wanted to set the 3rd bit in the 5th row true,
 *  this would also set the 5th bit in the 3rd row.
 * Input: a bitGraph, the index of the row of the bitGraph with the bit to
 *  be set, the index of the bit in that row that is to be set
 * Output: integer success value of 0 (and an altered bitGraph)
 * **********************************************************************/

int
bitGraphSetTrueSym(bitGraph_t *bg, int x, int y){
	setTrue(bg->graph[x], y);
	setTrue(bg->graph[y], x);
	return 0;
}

/***********************************************************************
 * Sets the main diagonal of a bitGraph true.
 * Input: a bitGraph
 * Output: integer success value of 0 (and an altered bitGraph)
 * ********************************************************************/

int
bitGraphSetTrueDiagonal(bitGraph_t *bg){
	int i;
	for ( i=0 ; i<bg->size ; i++){
		setTrue(bg->graph[i], i);
	}
	return 0;
}

/***********************************************************************
 * Sets the main diagonal of a bitGraph false.
 * Input: a bitGraph
 * Output: integer success value of 0 (and an altered bitGraph)
 * ********************************************************************/

int
bitGraphSetFalseDiagonal(bitGraph_t *bg){
	int i;
	for ( i=0 ; i<bg->size ; i++){
		setFalse(bg->graph[i], i);
	}
	return 0;
}

/**********************************************************************
 * Prints a representation of a bitGraph using printBinaryBitSet.
 * Input: a bitGraph
 * Output: integer success value of 0 (and stdout text as described above)
 * *********************************************************************/

int
printBitGraph(bitGraph_t *bg){
	int i;
	for ( i=0 ; i<bg->size ; i++){
		printBinaryBitSet(bg->graph[i]);
		printf("\n");
	}
	return 0;
}

/***********************************************************************
 * Makes a bitGraph contain only true bits according to the bitmask given.
 *  Only locations with the row and column both true in the bitmask can be true
 *  if they were initially true.  If they were false, they remain false.
 *  If the location does not have both the row and the column in the
 *  bitmask, it is made false.
 * NOTE: not currently used in Gemoda.
 * Input: a bitGraph, a mask in the form of a bitSet
 * Output: integer success value of 0 (and an altered bitGraph)
 * *********************************************************************/

int
maskBitGraph(bitGraph_t *bg1, bitSet_t *bs){
	int i;
	for ( i=0 ; i<bg1->size ; i++ ){
		if(checkBit(bs, i)){
			bitSetIntersection(bg1->graph[i], bs, bg1->graph[i]);
		}else{
			emptySet(bg1->graph[i]);
		}
	}
	return 0;
}

/************************************************************************
 * Sets all bits in the bitGraph true.
 * Input: a bitGraph
 * Output: integer success value of 0 (and a bitGraph with all true bits)
 * **********************************************************************/

int
fillBitGraph(bitGraph_t *bg1){
	int i;
	for ( i=0 ; i<bg1->size ; i++ ){
		fillSet(bg1->graph[i]);
	}
	return 0;
}

/************************************************************************
 * Sets all bits in the bitGraph false.
 * Input: a bitGraph
 * Output: integer success value of 0 (and a bitGraph with all false bits)
 * **********************************************************************/

int
emptyBitGraph(bitGraph_t *bg1){
	int i;
	for ( i=0 ; i<bg1->size ; i++ ){
		emptySet(bg1->graph[i]);
	}
	return 0;
}

/*************************************************************************
 * Creates a bitGraph data structure.
 * Input: the size of the (square) bitGraph
 * Output: a new bitGraph data structure.
 * ***********************************************************************/

bitGraph_t *
newBitGraph(int size){
	bitGraph_t *bg = NULL;
	int i;
	bg = (bitGraph_t *) malloc (sizeof(bitGraph_t));
	if(bg == NULL){
		fprintf(stderr,"Memory error - Cannot allocate bitGraph - "
				"newBitGraph\n%s\n",strerror(errno));
		fflush(stderr);
		exit(0);
	}
	bg->size = size;
	bg->graph = (bitSet_t **) malloc (size * sizeof(bitSet_t *));
	if(bg->graph == NULL){
		fprintf(stderr,"Memory error - Cannot allocate bitGraphGraph - "
				"newBitGraph\n%s\n",strerror(errno));
		fflush(stderr);
		exit(0);
	}
	for ( i=0 ; i<size ; i++ ){
		bg->graph[i] = newBitSet(size);
	}
	return bg;
}

/***********************************************************************
 * Sets all bits in a bitGraph row (a bitSet) galse.
 * Input: a bitGraph, a row in the bitGraph to be emptied
 * Output: integer success value of 0 (and an altered bitGraph).
 * *********************************************************************/

int
emptyBitGraphRow(bitGraph_t * bg, int row){
	emptySet(bg->graph[row]);
	return 0;
}

/**************************************************************************
 * Deletes a bitGraph from memory.
 * Input: a bitGraph to be deleted.
 * Output: integer success value from 0 (and deletion of a bitGraph).
 * **********************************************************************/

int
deleteBitGraph(bitGraph_t *bg){
	int i;
	if(bg != NULL){
		if(bg->graph != NULL){
			for ( i=0 ; i<bg->size ; i++ ){
				deleteBitSet(bg->graph[i]);
			}
			free(bg->graph);
			bg->graph = NULL;
		}
		free(bg);
		bg = NULL;
	}
	return 0;
}
