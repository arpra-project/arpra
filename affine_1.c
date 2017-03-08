/*
 * affine_1.c
 *
 *  Created on: 24 Oct 2016
 *      Author: jt273
 */

#include "mpfa.h"
#include <malloc.h>
#include <assert.h>

void mpfa_affine_1 (mpfa_ptr z, mpfa_srcptr x, mpfr_ptr alpha, mpfr_ptr gamma, mpfr_ptr delta) {
	unsigned zTerm;
	mpfr_t u, temp, error;
	mpfr_prec_t prec;
	mpfa_t zNew;

	prec = mpfr_get_prec(&(z->centre));
	mpfr_init2(u, prec);
	mpfr_init2(temp, prec);
	mpfr_init2(error, prec);
	mpfa_init2(zNew, prec);
	mpfr_set_si(&(zNew->radius), 0, MPFR_RNDN);

	assert(mpfr_set_si(u, -prec, MPFR_RNDN) == 0); // fails if emax <= log2(prec)
	assert(mpfr_exp2(u, u, MPFR_RNDN) == 0); // fails if emin > 1-prec

	if (mpfr_mul(temp, alpha, &(x->centre), MPFR_RNDN)) {
		mpfr_mul(error, u, temp, MPFR_RNDU);
		mpfr_abs(error, error, MPFR_RNDN);
		mpfr_add(delta, delta, error, MPFR_RNDU);
	}

	if (mpfr_add(&(zNew->centre), gamma, temp, MPFR_RNDN)) {
		mpfr_mul(error, u, &(zNew->centre), MPFR_RNDU);
		mpfr_abs(error, error, MPFR_RNDN);
		mpfr_add(delta, delta, error, MPFR_RNDU);
	}

	zNew->nTerms = x->nTerms + 1;
	zNew->symbols = malloc(zNew->nTerms * sizeof(unsigned));
	zNew->deviations = malloc(zNew->nTerms * sizeof(mpfr_t));

	for (zTerm = 0; zTerm < (zNew->nTerms - 1); zTerm++) {
		zNew->symbols[zTerm] = x->symbols[zTerm];
		mpfr_init2(&(zNew->deviations[zTerm]), prec);

		if (mpfr_mul(&(zNew->deviations[zTerm]), alpha, &(x->deviations[zTerm]), MPFR_RNDN)) {
			mpfr_mul(error, u, &(zNew->deviations[zTerm]), MPFR_RNDU);
			mpfr_abs(error, error, MPFR_RNDN);
			mpfr_add(delta, delta, error, MPFR_RNDU);
		}

		mpfr_abs(temp, &(zNew->deviations[zTerm]), MPFR_RNDN);
		mpfr_add(&(zNew->radius), &(zNew->radius), temp, MPFR_RNDU);
	}

	zNew->symbols[zTerm] = mpfa_next_sym();
	mpfr_init2(&(zNew->deviations[zTerm]), prec);
	mpfr_set(&(zNew->deviations[zTerm]), delta, MPFR_RNDN);
	mpfr_add(&(zNew->radius), &(zNew->radius), delta, MPFR_RNDU);

	mpfr_clear(u);
	mpfr_clear(temp);
	mpfr_clear(error);
	mpfa_set(z, zNew);
	mpfa_clear(zNew);
}
