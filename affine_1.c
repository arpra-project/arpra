/*
 * affine_1.c
 *
 *  Created on: 24 Oct 2016
 *      Author: jt273
 */

#include "mpfa.h"
#include <malloc.h>

void mpfa_affine_1 (mpfa_ptr z, mpfr_srcptr alpha, mpfa_srcptr x, mpfr_ptr gamma, mpfr_ptr delta) {
	unsigned zTerm;
	mpfa_t z_new;
	mpfr_t temp;

	mpfr_init(temp);
	z_new->nTerms = x->nTerms + 1;
	z_new->symbols = malloc(z_new->nTerms * sizeof(unsigned));
	z_new->deviations = malloc(z_new->nTerms * sizeof(mpfr_t));
	mpfr_init(&(z_new->centre));
	mpfa_term_linear_1(&(z_new->centre), alpha, &(x->centre), gamma, delta);
	mpfr_init_set_d(&(z_new->radius), 0.0, MPFR_RNDN);

	for (zTerm = 0; zTerm < (z_new->nTerms - 1); zTerm++) {
		z_new->symbols[zTerm] = x->symbols[zTerm];
		mpfr_init(&(z_new->deviations[zTerm]));
		mpfa_term_linear_1(&(z_new->deviations[zTerm]), alpha, &(x->deviations[zTerm]), NULL, delta);
		mpfr_abs(temp, &(z_new->deviations[zTerm]), MPFR_RNDN);
		mpfr_add(&(z_new->radius), &(z_new->radius), temp, MPFR_RNDU);
	}

	z_new->symbols[z_new->nTerms - 1] = mpfa_next_sym();
	mpfr_init_set(&(z_new->deviations[z_new->nTerms - 1]), delta, MPFR_RNDU);
	mpfr_add(&(z_new->radius), &(z_new->radius), delta, MPFR_RNDU);

	mpfa_clear(z);
	*z = *z_new;
	mpfr_clear(temp);
}
