/*
 * mul.c
 *
 *  Created on: 11 Nov 2016
 *      Author: jt273
 */

#include "mpfa.h"
#include <malloc.h>

void mpfa_mul (mpfa_ptr z, mpfa_srcptr x, mpfa_srcptr y) {
	unsigned xTerm, yTerm, zTerm;
	mpfa_t z_new;
	mpfr_t temp, delta;

	mpfr_init(temp);
	mpfr_init_set_d(delta, 0.0, MPFR_RNDN);
	z_new->nTerms = x->nTerms + y->nTerms + 1;
	z_new->symbols = malloc(z_new->nTerms * sizeof(unsigned));
	z_new->deviations = malloc(z_new->nTerms * sizeof(mpfr_t));
	mpfr_init(&(z_new->centre));
	mpfa_term_linear_1 (&(z_new->centre), &(x->centre), &(y->centre), NULL, delta);
	mpfr_init_set_d(&(z_new->radius), 0.0, MPFR_RNDN);

	for (xTerm = 0, yTerm = 0, zTerm = 0; ((xTerm < x->nTerms) && (yTerm < y->nTerms)); zTerm++) {
		if ((yTerm == y->nTerms) || (x->symbols[xTerm] < y->symbols[yTerm])) {
			z_new->symbols[zTerm] = x->symbols[xTerm];
			mpfr_init(&(z_new->deviations[zTerm]));
			mpfa_term_linear_1(&(z_new->deviations[zTerm]), &(y->centre), &(x->deviations[xTerm]), NULL, delta);
			xTerm++;
		}
		else if ((xTerm == x->nTerms) || (y->symbols[yTerm] < x->symbols[xTerm])) {
			z_new->symbols[zTerm] = y->symbols[yTerm];
			mpfr_init(&(z_new->deviations[zTerm]));
			mpfa_term_linear_1(&(z_new->deviations[zTerm]), &(x->centre), &(y->deviations[yTerm]), NULL, delta);
			yTerm++;
		}
		else {
			z_new->symbols[zTerm] = x->symbols[xTerm];
			mpfr_init(&(z_new->deviations[zTerm]));
			mpfa_term_linear_2(&(z_new->deviations[zTerm]), &(y->centre), &(x->deviations[xTerm]), &(x->centre), &(y->deviations[yTerm]), NULL, delta);
			xTerm++;
			yTerm++;
		}
		mpfr_abs(temp, &(z_new->deviations[zTerm]), MPFR_RNDN);
		mpfr_add(&(z_new->radius), &(z_new->radius), temp, MPFR_RNDU);
	}

	// Trivial estimate of linearisation error
	mpfr_mul(temp, &(x->radius), &(y->radius), MPFR_RNDU);
	mpfr_add(delta, delta, temp, MPFR_RNDU);

	z_new->nTerms = zTerm + 1;
	z_new->symbols = realloc(z_new->symbols, z_new->nTerms * sizeof(unsigned));
	z_new->symbols[zTerm] = mpfa_next_sym();
	z_new->deviations = realloc(z_new->symbols, z_new->nTerms * sizeof(mpfr_t));
	mpfr_init_set(&(z_new->deviations[zTerm]), delta, MPFR_RNDU);
	mpfr_add(&(z_new->radius), &(z_new->radius), delta, MPFR_RNDU);

	mpfa_clear(z);
	*z = *z_new;
	mpfr_clear(temp);
	mpfr_clear(delta);
}
