/*
 * set_mpfr.c
 *
 *  Created on: 14 Mar 2017
 *      Author: james
 */

#include "mpfa.h"
#include <malloc.h>
#include <assert.h>

void mpfa_set_mpfr (mpfa_ptr z, mpfr_srcptr centre, mpfr_srcptr radius) {
	unsigned zTerm;
	mpfr_prec_t prec;

	prec = mpfr_get_prec(&(z->centre));

	if (mpfr_set(&(z->centre), centre, MPFR_RNDN)) {
		assert(mpfr_set_si(&(z->radius), (-prec + mpfr_get_exp(&(z->centre))), MPFR_RNDN) == 0);
		assert(mpfr_exp2(&(z->radius), &(z->radius), MPFR_RNDN) == 0);
		mpfr_add(&(z->radius), &(z->radius), radius, MPFR_RNDU);
	}
	else {
		mpfr_set(&(z->radius), radius, MPFR_RNDU);
	}

	if (mpfr_zero_p(&(z->radius))) {
		if (z->nTerms > 0) {
			for (zTerm = 0; zTerm < z->nTerms; zTerm++) {
				mpfr_clear(&(z->deviations[zTerm]));
			}
			z->nTerms = 0;
			free(z->symbols);
			free(z->deviations);
		}
	}
	else {
		if (z->nTerms == 0) {
			z->symbols = malloc(sizeof(unsigned));
			z->deviations = malloc(sizeof(mpfa_t));
			mpfr_init2(&(z->deviations[0]), prec);
		}
		else if (z->nTerms >= 2) {
			for (zTerm = 1; zTerm < z->nTerms; zTerm++) {
				mpfr_clear(&(z->deviations[zTerm]));
			}
			z->symbols = realloc(z->symbols, sizeof(unsigned));
			z->deviations = realloc(z->deviations, sizeof(mpfa_t));
		}
		z->nTerms = 1;
		z->symbols[0] = mpfa_next_sym();
		mpfr_set(&(z->deviations[0]), &(z->radius), MPFR_RNDN);
	}
}
