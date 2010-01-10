
#ifndef __MATRIXMAP_H_
#define __MATRIXMAP_H_

#include "matdata.h"
#include "matrices.h"

struct {
	char *name;
	const int (*mat)[MATRIX_SIZE];
} matrix_map[] = {
	{"dna_idmat",	dna_idmat},
	{"dna_negidmat",	dna_negidmat},
	{"identity_aa",	identity_aa},
	{"idmat",	idmat},
	{"blosum100",	blosum100},
	{"blosum30",	blosum30},
	{"blosum35",	blosum35},
	{"blosum40",	blosum40},
	{"blosum45",	blosum45},
	{"blosum50",	blosum50},
	{"blosum55",	blosum55},
	{"blosum60",	blosum60},
	{"blosum62",	blosum62},
	{"blosum65",	blosum65},
	{"blosum70",	blosum70},
	{"blosum75",	blosum75},
	{"blosum80",	blosum80},
	{"blosum85",	blosum85},
	{"blosum90",	blosum90},
	{"blosumn",	blosumn},
	{"dayhoff",	dayhoff},
	{"pam100",	pam100},
	{"pam110",	pam110},
	{"pam120",	pam120},
	{"pam130",	pam130},
	{"pam140",	pam140},
	{"pam150",	pam150},
	{"pam160",	pam160},
	{"pam190",	pam190},
	{"pam200",	pam200},
	{"pam210",	pam210},
	{"pam220",	pam220},
	{"pam230",	pam230},
	{"pam240",	pam240},
	{"pam250",	pam250},
	{"pam260",	pam260},
	{"pam280",	pam280},
	{"pam290",	pam290},
	{"pam300",	pam300},
	{"pam310",	pam310},
	{"pam320",	pam320},
	{"pam330",	pam330},
	{"pam340",	pam340},
	{"pam360",	pam360},
	{"pam370",	pam370},
	{"pam380",	pam380},
	{"pam390",	pam390},
	{"pam400",	pam400},
	{"pam430",	pam430},
	{"pam440",	pam440},
	{"pam450",	pam450},
	{"pam460",	pam460},
	{"pam490",	pam490},
	{"pam500",	pam500},
	{"phat_t75_b73",	phat_t75_b73},
	{"phat_t80_b78",	phat_t80_b78},
	{"phat_t85_b82",	phat_t85_b82},
	{"alpha_mat",	alpha_mat},
	{"beta_mat",	beta_mat},
	{"coil_mat",	coil_mat},

	{0, }
};

#endif	/*__MATRIXMAP_H_*/
