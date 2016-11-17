/*
 * init_nterms.c
 *
 *  Created on: 15 Nov 2016
 *      Author: jt273
 */

#include "mpfa.h"
#include <malloc.h>

void mpfa_init_nterms (mpfa_ptr x, unsigned n) {
	unsigned xTerm;
	x->symbols = malloc(n * sizeof(unsigned));
	x->deviations = malloc(n * sizeof(mpfr_t));
	for (xTerm = 0; xTerm < x->nTerms; xTerm++) {
		mpfr_init(&(x->deviations[xTerm]));
	}
	mpfr_init(&(x->centre));
	mpfr_init(&(x->radius));
}
