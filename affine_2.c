/*
 * affine_2.c
 *
 *  Created on: 2 Nov 2016
 *      Author: jt273
 */

#include "mpfa.h"
#include <malloc.h>

void mpfa_affine_2 (mpfa_ptr z, mpfr_srcptr alpha, mpfa_srcptr x, mpfr_srcptr beta, mpfa_srcptr y, mpfr_ptr gamma, mpfr_ptr delta) {
	unsigned prec, xTerm, yTerm, zTerm;
	mpfr_t u, temp, error;
	mpfa_t z_new;

	prec = mpfr_get_prec(&(z->centre));
	mpfr_init2(u, prec);
	mpfr_init2(temp, prec);
	mpfr_init2(error, prec);
	mpfr_init2(&(z_new->centre), prec);
	mpfr_init2(&(z_new->radius), prec);
	mpfr_set_d(&(z_new->radius), 0.0, MPFR_RNDN);

	mpfr_set_si(u, -prec, MPFR_RNDN);
	mpfr_exp2(u, u, MPFR_RNDN);

	mpfr_mul(temp, alpha, &(x->centre), MPFR_RNDN);
	mpfr_mul(error, u, temp, MPFR_RNDU);
	mpfr_add(delta, delta, error, MPFR_RNDU);

	mpfr_add(&(z_new->centre), gamma, temp, MPFR_RNDN);
	mpfr_mul(error, u, &(z_new->centre), MPFR_RNDU);
	mpfr_add(delta, delta, error, MPFR_RNDU);

	mpfr_mul(temp, beta, &(y->centre), MPFR_RNDN);
	mpfr_mul(error, u, temp, MPFR_RNDU);
	mpfr_add(delta, delta, error, MPFR_RNDU);

	mpfr_add(&(z_new->centre), &(z_new->centre), temp, MPFR_RNDN);
	mpfr_mul(error, u, &(z_new->centre), MPFR_RNDU);
	mpfr_add(delta, delta, error, MPFR_RNDU);

	z_new->nTerms = x->nTerms + y->nTerms + 1;
	z_new->symbols = malloc(z_new->nTerms * sizeof(unsigned));
	z_new->deviations = malloc(z_new->nTerms * sizeof(mpfr_t));

	for (xTerm = 0, yTerm = 0, zTerm = 0; ((xTerm < x->nTerms) && (yTerm < y->nTerms)); zTerm++) {
		if ((yTerm == y->nTerms) || (x->symbols[xTerm] < y->symbols[yTerm])) {
			z_new->symbols[zTerm] = x->symbols[xTerm];
			mpfr_init(&(z_new->deviations[zTerm]));

			mpfr_mul(&(z_new->deviations[zTerm]), alpha, &(x->deviations[xTerm]), MPFR_RNDN);
			mpfr_mul(error, u, &(z_new->deviations[zTerm]), MPFR_RNDU);
			mpfr_add(delta, delta, error, MPFR_RNDU);

			xTerm++;
		}
		else if ((xTerm == x->nTerms) || (y->symbols[yTerm] < x->symbols[xTerm])) {
			z_new->symbols[zTerm] = y->symbols[yTerm];
			mpfr_init(&(z_new->deviations[zTerm]));

			mpfr_mul(&(z_new->deviations[zTerm]), beta, &(y->deviations[yTerm]), MPFR_RNDN);
			mpfr_mul(error, u, &(z_new->deviations[zTerm]), MPFR_RNDU);
			mpfr_add(delta, delta, error, MPFR_RNDU);

			yTerm++;
		}
		else {
			z_new->symbols[zTerm] = x->symbols[xTerm];
			mpfr_init(&(z_new->deviations[zTerm]));

			mpfr_mul(&(z_new->deviations[zTerm]), alpha, &(x->deviations[xTerm]), MPFR_RNDN);
			mpfr_mul(error, u, &(z_new->deviations[zTerm]), MPFR_RNDU);
			mpfr_add(delta, delta, error, MPFR_RNDU);

			mpfr_mul(temp, beta, &(y->deviations[yTerm]), MPFR_RNDN);
			mpfr_mul(error, u, temp, MPFR_RNDU);
			mpfr_add(delta, delta, error, MPFR_RNDU);

			mpfr_add(&(z_new->deviations[zTerm]), &(z_new->deviations[zTerm]), temp, MPFR_RNDN);
			mpfr_mul(error, u, &(z_new->deviations[zTerm]), MPFR_RNDU);
			mpfr_add(delta, delta, error, MPFR_RNDU);

			xTerm++;
			yTerm++;
		}
		mpfr_abs(temp, &(z_new->deviations[zTerm]), MPFR_RNDN);
		mpfr_add(&(z_new->radius), &(z_new->radius), temp, MPFR_RNDU);
	}

	z_new->nTerms = zTerm + 1;
	z_new->symbols = realloc(z_new->symbols, z_new->nTerms * sizeof(unsigned));
	z_new->symbols[zTerm] = mpfa_next_sym();
	z_new->deviations = realloc(z_new->symbols, z_new->nTerms * sizeof(mpfr_t));
	mpfr_init_set(&(z_new->deviations[zTerm]), delta, MPFR_RNDU);
	mpfr_add(&(z_new->radius), &(z_new->radius), delta, MPFR_RNDU);

	mpfr_clear(u);
	mpfr_clear(temp);
	mpfr_clear(error);
	mpfa_clear(z);
	*z = *z_new;
}
