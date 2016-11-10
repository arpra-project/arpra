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

	z->nTerms = x->nTerms + y->nTerms + 1;
	z->symbols = malloc(sizeof(unsigned) * z->nTerms);
	z->deviations = malloc(sizeof(mpfr_t) * z->nTerms);
	mpfa_affine_term_2(&(z->centre), alpha, &(x->centre), beta, &(y->centre), gamma, delta);
	mpfr_set_d(&(z->radius), 0.0, MPFR_RNDN);

	unsigned xTerm = 0, yTerm = 0, zTerm = 0;
	while ((xTerm < x->nTerms) && (yTerm < y->nTerms)) {
		if ((yTerm == y->nTerms) || (x->symbols[xTerm] < y->symbols[yTerm])) {
			z->symbols[zTerm] = x->symbols[xTerm];
			mpfr_init(&(z->deviations[zTerm]));
			mpfa_affine_term_1(&(z->deviations[zTerm]), alpha, &(x->deviations[xTerm]), NULL, delta);
			xTerm++;
		}
		else if ((xTerm == x->nTerms) || (y->symbols[yTerm] < x->symbols[xTerm])) {
			z->symbols[zTerm] = y->symbols[yTerm];
			mpfr_init(&(z->deviations[zTerm]));
			mpfa_affine_term_1(&(z->deviations[zTerm]), beta, &(y->deviations[yTerm]), NULL, delta);
			yTerm++;
		}
		else {
			z->symbols[zTerm] = x->symbols[xTerm];
			mpfr_init(&(z->deviations[zTerm]));
			mpfa_affine_term_2(&(z->deviations[zTerm]), alpha, &(x->deviations[xTerm]), beta, &(y->deviations[yTerm]), NULL, delta);
			xTerm++;
			yTerm++;
		}
		mpfr_abs(temp, &(z->deviations[zTerm]), MPFR_RNDN);
		mpfr_add(&(z->radius), &(z->radius), temp, MPFR_RNDU);
		zTerm++;
	}

	z->nTerms = zTerm + 1;
	z->symbols = realloc(z->symbols, sizeof(unsigned) * z->nTerms);
	z->symbols[zTerm] = mpfa_next_sym();
	z->deviations = realloc(z->symbols, sizeof(mpfr_t) * z->nTerms);
	mpfr_init_set(&(z->deviations[zTerm]), delta, MPFR_RNDU);
	mpfr_add(&(z->radius), &(z->radius), delta, MPFR_RNDU);

	mpfr_clear(temp);
}
