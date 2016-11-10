/*
 * affine_2.c
 *
 *  Created on: 2 Nov 2016
 *      Author: jt273
 */

#include "mpfa.h"

#include <malloc.h>

void mpfa_affine_2 (mpfa_ptr z, mpfr_srcptr alpha, mpfa_srcptr x, mpfr_srcptr beta, mpfa_srcptr y, mpfr_ptr gamma, mpfr_ptr delta) {
	mpfr_t temp;
	mpfr_init(temp);

	mpfa_affine_term_2(&(z->centre), alpha, &(x->centre), beta, &(y->centre), gamma, delta);
	mpfr_set_d(&(z->radius), 0.0, MPFR_RNDN);
	z->nTerms = x->nTerms + 1;
	z->symbols = malloc(sizeof(unsigned) * z->nTerms);
	z->deviations = malloc(sizeof(mpfr_t) * z->nTerms);

	unsigned xTerm, yTerm, zTerm;
	for (zTerm = 0; zTerm < (z->nTerms - 1); zTerm++) {
		z->symbols[zTerm] = x->symbols[zTerm];
		mpfr_init(&(z->deviations[zTerm]));
		mpfa_affine_term_2(&(z->deviations[zTerm]), alpha, &(x->deviations[zTerm]), beta, &(y->centre), NULL, delta);
		mpfr_abs(temp, &(z->deviations[zTerm]), MPFR_RNDN);
		mpfr_add(&(z->radius), &(z->radius), temp, MPFR_RNDU);
	}

	mpfr_add(&(z->radius), &(z->radius), delta, MPFR_RNDU);
	z->symbols[z->nTerms - 1] = mpfa_next_sym();
	mpfr_init_set(&(z->deviations[z->nTerms - 1]), delta, MPFR_RNDU);
	mpfr_clear(temp);
}
