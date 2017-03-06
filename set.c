/*
 * set.c
 *
 *  Created on: 21 Oct 2016
 *      Author: jt273
 */

#include "mpfa.h"
#include <malloc.h>

void mpfa_set (mpfa_ptr z, mpfa_srcptr x) {
	unsigned zTerm;
	mpfr_prec_t prec;

	if (z == x) return;

	if (z->nTerms < x->nTerms) {
		// need to grow z then initialise extra terms
		if (z->nTerms == 0) {
			z->symbols = malloc(x->nTerms * sizeof(unsigned));
			z->deviations = malloc(x->nTerms * sizeof(mpfa_t));
		}
		else {
			z->symbols = realloc(z->symbols, x->nTerms * sizeof(unsigned));
			z->deviations = realloc(z->deviations, x->nTerms * sizeof(mpfa_t));
		}
		prec = mpfr_get_prec(&(z->centre));
		for (zTerm = z->nTerms; zTerm < x->nTerms; zTerm++) {
			mpfr_init2(&(z->deviations[zTerm]), prec);
		}
	}
	else if (z->nTerms > x->nTerms) {
		// need to clear extra terms then shrink z
		for (zTerm = x->nTerms; zTerm < z->nTerms; zTerm++) {
			mpfr_clear(&(z->deviations[zTerm]));
		}
		if (x->nTerms == 0) {
			free(z->symbols);
			free(z->deviations);
		}
		else {
			z->symbols = realloc(z->symbols, x->nTerms * sizeof(unsigned));
			z->deviations = realloc(z->deviations, x->nTerms * sizeof(mpfa_t));
		}
	}

	z->nTerms = x->nTerms;
	for (zTerm = 0; zTerm < z->nTerms; zTerm++) {
		z->symbols[zTerm] = x->symbols[zTerm];
		mpfr_set(&(z->deviations[zTerm]), &(x->deviations[zTerm]), MPFR_RNDN);
	}
	mpfr_set(&(z->centre), &(x->centre), MPFR_RNDN);
	mpfr_set(&(z->radius), &(x->radius), MPFR_RNDU);
}
