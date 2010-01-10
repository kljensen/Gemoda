#ifndef BIT_SET_H
#define BIT_SET_H


#include "stdio.h"
#include "stdlib.h"
#include "string.h"

typedef unsigned int bit_t;

typedef struct {
	int max;
	int slots;
	int bytes;
	bit_t *tf;
} bitSet_t;

typedef struct{
	int size;
	bitSet_t **graph;
}bitGraph_t;

/************************************************************************
 * Some #define statments to make short functions that will hopefully be
 * optimized by the compiler.
 * ********************************************************************/

// Get the size of a bit_t, which is an unsigned int.
#define BSBITSIZE (sizeof(bit_t) * 8)

// Uses bit operations to make a mask the size of a bit_t with the y'th 
//  bit true and all other bits false.  Used for testing a single bit.
#define BSMASK(y) ( ((bit_t) 1) << y % BSBITSIZE )

// Finds which bit_t, or "slot", in a bitset_t the y'th bit belongs to.  
//  Also used for testing single bits.
#define BITSLOT(y) ( y / BSBITSIZE )

// Sets the y'th bit in x (a bitset_t) to be true using bitwise operators.
#define BSSET( x, y) ( x[BITSLOT(y)] |= BSMASK(y) )

// Sets the y'th bit in x (a bitset_t) to be false using bitwise operators.
#define BSCLEAR(x, y) ( x[BITSLOT(y)] &= ~BSMASK(y) )

// Tests whether the y'th bit in x (a bitset_t) is true.
#define BSTEST(x, y) ( x[BITSLOT(y)] & BSMASK(y) )

// Finds the total number of bit_t's ("slot"s) that are necessary for a
//  bitset of length n.  Uses integer division and supplements by 
//  BSBITSIZE - 1 to make sure that for slots that are less than full, a slot
//  is still allocated, and for slots that are full, no extra slot is
//  allocated.
#define BSNUMSLOTS(n) ((n + BSBITSIZE - 1) / BSBITSIZE)

// Performs a union operation on two bit_t's with bitwise operators.
#define BSUNION(x,y) ((x)|(y))

// Performs an intersection operation on two bit_t's with bitwise operators.
#define BSINTERSECTION(x,y) ((x)&(y))


bit_t *
newBitArray(int bytes);

bitSet_t *
newBitSet(int size);

int
setTrue(bitSet_t * s1, int x);

int
setFalse(bitSet_t * s1, int x);

int 
bitSetDifference(bitSet_t * s1, bitSet_t * s2, bitSet_t * s3);

int 
bitSet3WayDifference(bitSet_t * s1, bitSet_t * s2, bitSet_t * s3, bitSet_t * s4);

int 
bitSetSum(bitSet_t * s1, bitSet_t * s2, bitSet_t * s3);

int
flipBits(bitSet_t *s1);

int
fillSet(bitSet_t * s1);

int
emptySet(bitSet_t * s1);

int
checkBit(bitSet_t * s1, int x);

int
copyBitGraph(bitGraph_t *bg1, bitGraph_t *bg2);

int
copySet(bitSet_t *s1, bitSet_t *s2);

int
deleteBitSet(bitSet_t * s1);

int
bitSetUnion(bitSet_t * s1, bitSet_t * s2, bitSet_t * s3);

int
bitSetIntersection(bitSet_t * s1, bitSet_t * s2, bitSet_t * s3);

int
bitSet3WayIntersection(bitSet_t * s1, bitSet_t * s2, bitSet_t * s3, bitSet_t * s4);

int
countSet(bitSet_t * s1);

int
fillBitGraph(bitGraph_t *bg1);

int
emptyBitGraph(bitGraph_t *bg1);

int
printBitSet(bitSet_t * s1);

int
nextBitBitSet(bitSet_t *s1, int start);

int
countBitGraphNonZero(bitGraph_t *bg);

int
printBinaryBitSet(bitSet_t *s1);

int
bitGraphSetTrue(bitGraph_t *bg, int x, int y);

int
bitGraphSetTrueSym(bitGraph_t *bg, int x, int y);

int
bitGraphSetTrueDiagonal(bitGraph_t *bg);

int
bitGraphSetFalseDiagonal(bitGraph_t *bg);

int
printBitGraph(bitGraph_t *bg);

bitGraph_t *
newBitGraph(int size);

int
deleteBitGraph(bitGraph_t *bg);

int
bitGraphRowUnion(bitGraph_t *bg, int row1, int row2, bitSet_t *s1);

int
bitGraphRowIntersection(bitGraph_t *bg, int row1, int row2, bitSet_t *s1);

int
bitGraphCheckBit(bitGraph_t *bg, int x, int y);

int 
bitGraphSetFalseSym(bitGraph_t *bg, int x, int y);

int
bitGraphSetFalse(bitGraph_t *bg, int x, int y);

int
emptyBitGraphRow(bitGraph_t * bg, int row);

int
maskBitGraph(bitGraph_t *bg1, bitSet_t *bs);

#endif

