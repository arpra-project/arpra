/*
 * affine_1.c
 *
 *  Created on: 24 Oct 2016
 *      Author: jt273
 */

#include "mpfa.h"
#include <malloc.h>
#include <assert.h>

void mpfa_affine_1 (mpfa_ptr z, mpfr_srcptr alpha, mpfa_srcptr x, mpfr_ptr gamma, mpfr_ptr delta) {
	unsigned zTerm;
	int inexact;
	mpfr_prec_t prec;
	mpfr_t u, temp, error;

	prec = mpfr_get_prec(&(z->centre));
	mpfr_init2(u, prec);
	mpfr_init2(temp, prec);
	mpfr_init2(error, prec);
	mpfr_set_d(&(z->radius), 0.0, MPFR_RNDN);

	assert(!mpfr_set_si(u, -prec, MPFR_RNDN)); // fails if emax <= log2(prec)
	assert(!mpfr_exp2(u, u, MPFR_RNDN)); // fails if emin > 1-prec

	inexact = mpfr_mul(temp, alpha, &(x->centre), MPFR_RNDN);
	if (inexact) {
		mpfr_mul(error, u, temp, MPFR_RNDU);
		mpfr_add(delta, delta, error, MPFR_RNDU);
	}

	inexact = mpfr_add(&(z->centre), gamma, temp, MPFR_RNDN);
	if (inexact) {
		mpfr_mul(error, u, &(z->centre), MPFR_RNDU);
		mpfr_add(delta, delta, error, MPFR_RNDU);
	}

	for (zTerm = 0; zTerm < z->nTerms; zTerm++) {
		mpfr_clear(&(z->deviations[zTerm]));
	}
	z->nTerms = x->nTerms + 1;
	z->symbols = realloc(z->symbols, z->nTerms * sizeof(unsigned));
	z->deviations = realloc(z->symbols, z->nTerms * sizeof(mpfr_t));

	for (zTerm = 0; zTerm < (z->nTerms - 1); zTerm++) {
		z->symbols[zTerm] = x->symbols[zTerm];
		mpfr_init(&(z->deviations[zTerm]));

		inexact = mpfr_mul(&(z->deviations[zTerm]), alpha, &(x->deviations[zTerm]), MPFR_RNDN);
		if (inexact) {
			mpfr_mul(error, u, &(z->deviations[zTerm]), MPFR_RNDU);
			mpfr_add(delta, delta, error, MPFR_RNDU);
		}

		mpfr_abs(temp, &(z->deviations[zTerm]), MPFR_RNDN);
		mpfr_add(&(z->radius), &(z->radius), temp, MPFR_RNDU);
	}

	z->symbols[z->nTerms - 1] = mpfa_next_sym();
	mpfr_init_set(&(z->deviations[z->nTerms - 1]), delta, MPFR_RNDU);
	mpfr_add(&(z->radius), &(z->radius), delta, MPFR_RNDU);

	mpfr_clear(u);
	mpfr_clear(temp);
	mpfr_clear(error);
}
