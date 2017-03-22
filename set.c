/*
 * set.c
 *
 *  Created on: 21 Oct 2016
 *      Author: jt273
 */

#include "mpfa.h"
#include <malloc.h>
#include <assert.h>

void mpfa_set (mpfa_ptr z, mpfa_srcptr x) {
	unsigned xTerm, zTerm;
	mpfr_prec_t prec;
	mpfr_t temp, error;

	if (z == x) return;

	prec = mpfr_get_prec(&(z->centre));
	mpfr_inits2(prec, temp, error, (mpfr_ptr) NULL);
	mpfr_set_si(error, 0, MPFR_RNDN);
	mpfr_set_si(&(z->radius), 0, MPFR_RNDN);

	if (mpfr_set(&(z->centre), &(x->centre), MPFR_RNDN)) {
		assert(mpfr_set_si(temp, (-prec + mpfr_get_exp(&(z->centre))), MPFR_RNDN) == 0);
		assert(mpfr_exp2(temp, temp, MPFR_RNDN) == 0);
		mpfr_add(error, error, temp, MPFR_RNDU);
	}

	if (z->nTerms < x->nTerms) {
		// need to grow z then initialise extra terms
		if (z->nTerms == 0) {
			z->symbols = malloc((x->nTerms + 1) * sizeof(unsigned));
			z->deviations = malloc((x->nTerms + 1) * sizeof(mpfa_t));
		}
		else {
			z->symbols = realloc(z->symbols, (x->nTerms + 1) * sizeof(unsigned));
			z->deviations = realloc(z->deviations, (x->nTerms + 1) * sizeof(mpfa_t));
		}
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
			z->symbols = realloc(z->symbols, (x->nTerms + 1) * sizeof(unsigned));
			z->deviations = realloc(z->deviations, (x->nTerms + 1) * sizeof(mpfa_t));
		}
	}

	for (xTerm = 0, zTerm = 0; xTerm < x->nTerms; xTerm++) {
		z->symbols[zTerm] = x->symbols[zTerm];

		if (mpfr_set(&(z->deviations[zTerm]), &(x->deviations[xTerm]), MPFR_RNDN)) {
			assert(mpfr_set_si(temp, (-prec + mpfr_get_exp(&(z->deviations[zTerm]))), MPFR_RNDN) == 0);
			assert(mpfr_exp2(temp, temp, MPFR_RNDN) == 0);
			mpfr_add(error, error, temp, MPFR_RNDU);
		}

		if (mpfr_zero_p(&(z->deviations[zTerm]))) {
			mpfr_clear(&(z->deviations[zTerm]));
		}
		else {
			mpfr_abs(temp, &(z->deviations[zTerm]), MPFR_RNDN);
			mpfr_add(&(z->radius), &(z->radius), temp, MPFR_RNDU);
			zTerm++;
		}
	}

	if (!mpfr_zero_p(error)) {
		z->symbols[zTerm] = mpfa_next_sym();
		mpfr_init2(&(z->deviations[zTerm]), prec);
		mpfr_set(&(z->deviations[zTerm]), error, MPFR_RNDN);
		mpfr_add(&(z->radius), &(z->radius), error, MPFR_RNDU);
		zTerm++;
	}

	z->nTerms = zTerm;
	if (z->nTerms == 0) {
		free(z->symbols);
		free(z->deviations);
	}
	else {
		z->symbols = realloc(z->symbols, z->nTerms * sizeof(unsigned));
		z->deviations = realloc(z->deviations, z->nTerms * sizeof(mpfr_t));
	}

	mpfr_clears(temp, error, (mpfr_ptr) NULL);
}
